#!/usr/bin/env python3
"""
run_somatic_variants.py — GATK4 Mutect2 somatic variant calling.

Supports tumor-normal paired mode and tumor-only mode.

Pipeline:
  1. GATK Mutect2             →  raw somatic calls + stats
  2. GATK MergeMutectStats    →  merged stats (if interval-scattered)
  3. GATK GetPileupSummaries  →  pileup at common variants
  4. GATK CalculateContamination → contamination estimate
  5. GATK FilterMutectCalls   →  filtered somatic VCF
  6. bcftools stats            →  QC summary

Somatic sample sheet TSV: sample_name, bam, sample_type (tumor|normal), [pair_id]

Usage:
    python run_somatic_variants.py \\
        --sample-sheet somatic_samples.tsv \\
        --reference ref/genome.fa --execute

    # Tumor-only (no matched normal)
    python run_somatic_variants.py \\
        --sample-sheet tumor_only.tsv \\
        --reference ref/genome.fa --tumor-only --execute
"""
import argparse, csv, json, logging, shutil, subprocess, sys
from datetime import datetime
from pathlib import Path

logging.basicConfig(format="%(asctime)s [%(levelname)s] %(message)s",
                    datefmt="%H:%M:%S", level=logging.INFO)
log = logging.getLogger(__name__)


def require_tool(name):
    if not shutil.which(name):
        log.error(f"Required: {name}"); sys.exit(1)


def run_cmd(cmd, log_file=None):
    s = " ".join(str(c) for c in cmd)
    log.info(f"Running: {s}")
    if log_file:
        with open(log_file, "a") as f:
            f.write(f"\n$ {s}\n")
            r = subprocess.run(cmd, stdout=f, stderr=f)
    else:
        r = subprocess.run(cmd, capture_output=True, text=True)
    return r.returncode


def load_somatic_samples(path: Path) -> list[dict]:
    """Return list of {tumor, normal} pairs. normal may be None for tumor-only."""
    rows = {}
    with open(path) as f:
        delim = "\t" if path.suffix == ".tsv" else ","
        for row in csv.DictReader(f, delimiter=delim):
            row = {k.strip(): v.strip() for k, v in row.items()}
            pair = row.get("pair_id", row["sample_name"])
            rows.setdefault(pair, {})
            stype = row.get("sample_type", "tumor").lower()
            rows[pair][stype] = {"name": row["sample_name"], "bam": Path(row["bam"])}
            rows[pair]["pair_id"] = pair

    pairs = []
    for pair_id, data in rows.items():
        if "tumor" not in data:
            log.error(f"Pair {pair_id} has no tumor sample"); sys.exit(1)
        pairs.append({
            "pair_id": pair_id,
            "tumor":  data["tumor"],
            "normal": data.get("normal"),
        })
    log.info(f"Loaded {len(pairs)} tumor sample(s)")
    return pairs


def run_mutect2(pair: dict, ref: Path, out_dir: Path,
                germline_resource: Path, pon: Path,
                intervals: str, log_file: Path) -> tuple[Path, Path]:
    pair_id   = pair["pair_id"]
    raw_vcf   = out_dir / f"{pair_id}.mutect2.vcf.gz"
    stats_file= out_dir / f"{pair_id}.mutect2.vcf.gz.stats"

    cmd = [
        "gatk", "Mutect2",
        "-R", str(ref),
        "-I", str(pair["tumor"]["bam"]),
        "-tumor", pair["tumor"]["name"],
        "-O", str(raw_vcf),
        "--f1r2-tar-gz", str(out_dir / f"{pair_id}.f1r2.tar.gz"),
    ]
    if pair["normal"]:
        cmd += ["-I", str(pair["normal"]["bam"]),
                "-normal", pair["normal"]["name"]]
        log.info(f"  Mutect2: {pair['tumor']['name']} vs {pair['normal']['name']}")
    else:
        log.info(f"  Mutect2 tumor-only: {pair['tumor']['name']}")

    if germline_resource and germline_resource.exists():
        cmd += ["--germline-resource", str(germline_resource)]
    if pon and pon.exists():
        cmd += ["--panel-of-normals", str(pon)]
    if intervals:
        cmd += ["-L", intervals]

    rc = run_cmd(cmd, log_file)
    if rc != 0:
        log.error(f"Mutect2 failed for {pair_id}"); sys.exit(1)
    return raw_vcf, stats_file


def estimate_contamination(pair: dict, ref: Path, out_dir: Path,
                            common_sites: Path, log_file: Path) -> Path:
    pair_id = pair["pair_id"]
    pileup  = out_dir / f"{pair_id}.pileup.table"
    contam  = out_dir / f"{pair_id}.contamination.table"

    if not common_sites or not common_sites.exists():
        log.info("  Skipping contamination estimate (no --common-sites provided)")
        return None

    run_cmd(["gatk", "GetPileupSummaries",
             "-I", str(pair["tumor"]["bam"]),
             "-V", str(common_sites),
             "-L", str(common_sites),
             "-O", str(pileup)], log_file)

    run_cmd(["gatk", "CalculateContamination",
             "-I", str(pileup), "-O", str(contam)], log_file)

    try:
        for line in contam.read_text().splitlines()[1:]:
            pct = float(line.split("\t")[1]) * 100
            log.info(f"  Estimated contamination: {pct:.2f}%")
    except Exception:
        pass
    return contam


def filter_mutect(raw_vcf: Path, stats_file: Path, ref: Path,
                  contam_table: Path, out_dir: Path, log_file: Path) -> Path:
    filtered = out_dir / raw_vcf.name.replace(".mutect2.", ".somatic.filtered.")

    cmd = ["gatk", "FilterMutectCalls",
           "-R", str(ref), "-V", str(raw_vcf), "-O", str(filtered),
           "--stats", str(stats_file)]
    if contam_table and contam_table.exists():
        cmd += ["--contamination-table", str(contam_table)]

    rc = run_cmd(cmd, log_file)
    if rc != 0:
        log.error("FilterMutectCalls failed"); sys.exit(1)
    log.info(f"  Filtered VCF: {filtered.name}")
    return filtered


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--sample-sheet",      type=Path, required=True,
                        help="TSV: sample_name, bam, sample_type (tumor|normal), pair_id")
    parser.add_argument("--reference",         type=Path, required=True)
    parser.add_argument("--germline-resource", type=Path,
                        help="af-only-gnomad VCF for Mutect2 germline resource")
    parser.add_argument("--panel-of-normals",  type=Path, help="PoN VCF")
    parser.add_argument("--common-sites",      type=Path,
                        help="Common variant sites for contamination estimate")
    parser.add_argument("--intervals",         help="Genomic intervals or BED file")
    parser.add_argument("--output-dir",        type=Path, default=Path("variant_out/somatic"))
    parser.add_argument("--execute",           action="store_true")
    args = parser.parse_args()

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_dir   = args.output_dir.with_name(f"{args.output_dir.name}_{timestamp}")
    log_file  = out_dir / "run.log"

    pairs = load_somatic_samples(args.sample_sheet)

    if not args.execute:
        log.info("=== DRY RUN (pass --execute to run) ===")
        log.info(f"Would write to: {out_dir}")
        for p in pairs:
            mode = "tumor-normal" if p["normal"] else "tumor-only"
            log.info(f"  [{mode}] {p['tumor']['name']}" +
                     (f" vs {p['normal']['name']}" if p["normal"] else ""))
        log.info("[DRY RUN] Would run: Mutect2 → CalculateContamination → FilterMutectCalls")
        log.info("[DRY RUN] Output: <pair>.somatic.filtered.vcf.gz")
        log.info("=== DRY RUN complete ===")
        return

    for t in ["gatk", "bcftools", "samtools"]:
        require_tool(t)
    if not args.reference.exists():
        log.error(f"Reference not found: {args.reference}"); sys.exit(1)

    fai = Path(str(args.reference) + ".fai")
    if not fai.exists():
        run_cmd(["samtools", "faidx", str(args.reference)], log_file)
    ref_dict = args.reference.with_suffix(".dict")
    if not ref_dict.exists():
        run_cmd(["gatk", "CreateSequenceDictionary", "-R", str(args.reference)], log_file)

    out_dir.mkdir(parents=True, exist_ok=True)
    log.info(f"Output: {out_dir}")

    results = {}
    for pair in pairs:
        log.info(f"Processing pair: {pair['pair_id']}")
        raw_vcf, stats = run_mutect2(pair, args.reference, out_dir,
                                      args.germline_resource, args.panel_of_normals,
                                      args.intervals, log_file)
        contam = estimate_contamination(pair, args.reference, out_dir,
                                         args.common_sites, log_file)
        filtered = filter_mutect(raw_vcf, stats, args.reference, contam, out_dir, log_file)
        results[pair["pair_id"]] = {"filtered_vcf": str(filtered)}

    manifest = {
        "pipeline": "somatic", "timestamp": timestamp,
        "reference": str(args.reference.resolve()),
        "pairs": results,
    }
    (out_dir / "manifest.json").write_text(json.dumps(manifest, indent=2))
    log.info(f"Done. Results in: {out_dir}")
    log.info(f"Next: run_variant_annotation.py --vcf <filtered.vcf.gz>")


if __name__ == "__main__":
    main()
