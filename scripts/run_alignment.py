#!/usr/bin/env python3
"""
run_alignment.py — Align DNA reads with BWA-MEM2, mark duplicates with GATK, and index.

Pipeline per sample:
  1. bwa-mem2 index  (build once, reuse)
  2. bwa-mem2 mem    (align paired-end reads → SAM)
  3. samtools sort   (coordinate-sorted BAM)
  4. GATK MarkDuplicatesSpark (mark/remove PCR duplicates)
  5. samtools index  (BAI index)
  6. samtools flagstat (alignment QC)

Sample sheet TSV columns: sample_name, fastq_r1, fastq_r2, [sample_type: tumor|normal]

Usage:
    python run_alignment.py --sample-sheet samplesheet.tsv \\
        --reference ref/genome.fa --execute
"""
import argparse, csv, hashlib, json, logging, shutil, subprocess, sys
from datetime import datetime
from pathlib import Path

logging.basicConfig(format="%(asctime)s [%(levelname)s] %(message)s",
                    datefmt="%H:%M:%S", level=logging.INFO)
log = logging.getLogger(__name__)


def require_tool(name):
    if not shutil.which(name):
        log.error(f"Required tool not found: {name}  →  python variant_preflight.py --pipeline alignment --emit-install-plan")
        sys.exit(1)


def run_cmd(cmd, log_file=None, env=None):
    s = " ".join(str(c) for c in cmd)
    log.info(f"Running: {s}")
    if log_file:
        log_file.parent.mkdir(parents=True, exist_ok=True)
        with open(log_file, "a") as f:
            f.write(f"\n$ {s}\n")
            r = subprocess.run(cmd, stdout=f, stderr=f, env=env)
    else:
        r = subprocess.run(cmd, capture_output=True, text=True)
    return r.returncode


def load_samples(path):
    samples = []
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t" if path.suffix == ".tsv" else ",")
        for row in reader:
            row = {k.strip(): v.strip() for k, v in row.items()}
            r1 = Path(row["fastq_r1"])
            r2 = Path(row["fastq_r2"])
            if not r1.exists():
                log.error(f"fastq_r1 not found: {r1}"); sys.exit(1)
            if not r2.exists():
                log.error(f"fastq_r2 not found: {r2}"); sys.exit(1)
            samples.append({
                "sample_name": row["sample_name"],
                "r1": r1, "r2": r2,
                "sample_type": row.get("sample_type", "germline"),
            })
    log.info(f"Loaded {len(samples)} sample(s)")
    return samples


def build_bwa_index(ref: Path, log_file: Path):
    log.info(f"Building BWA-MEM2 index for: {ref.name}")
    rc = run_cmd(["bwa-mem2", "index", str(ref)], log_file)
    if rc != 0:
        log.error("BWA-MEM2 index failed"); sys.exit(1)
    log.info("BWA-MEM2 index built")


def align_sample(sample: dict, ref: Path, out_dir: Path,
                 threads: int, log_file: Path) -> Path:
    name     = sample["sample_name"]
    sam_path = out_dir / f"{name}.sam"
    bam_path = out_dir / f"{name}.sorted.bam"
    out_dir.mkdir(parents=True, exist_ok=True)

    # Read group tag (required by GATK)
    rg = f"@RG\\tID:{name}\\tSM:{name}\\tPL:ILLUMINA\\tLB:{name}_lib1\\tPU:{name}_unit1"

    # BWA-MEM2 align
    align_cmd = [
        "bwa-mem2", "mem",
        "-t", str(threads),
        "-R", rg,
        str(ref), str(sample["r1"]), str(sample["r2"]),
    ]
    log.info(f"  Aligning {name}...")
    with open(log_file, "a") as lf, open(sam_path, "w") as sam_out:
        lf.write(f"\n$ {' '.join(str(c) for c in align_cmd)}\n")
        r = subprocess.run(align_cmd, stdout=sam_out, stderr=lf)
    if r.returncode != 0:
        log.error(f"BWA-MEM2 align failed for {name}"); sys.exit(1)

    # samtools sort
    log.info(f"  Sorting {name}...")
    rc = run_cmd(["samtools", "sort", "-@", str(threads),
                  "-o", str(bam_path), str(sam_path)], log_file)
    if rc != 0:
        log.error(f"samtools sort failed for {name}"); sys.exit(1)
    sam_path.unlink(missing_ok=True)

    return bam_path


def mark_duplicates(bam: Path, out_dir: Path, log_file: Path) -> Path:
    name      = bam.stem.replace(".sorted", "")
    dedup_bam = out_dir / f"{name}.markdup.bam"
    metrics   = out_dir / f"{name}.markdup_metrics.txt"

    log.info(f"  Marking duplicates: {name}")
    rc = run_cmd([
        "gatk", "MarkDuplicatesSpark",
        "-I", str(bam),
        "-O", str(dedup_bam),
        "-M", str(metrics),
        "--tmp-dir", str(out_dir),
    ], log_file)
    if rc != 0:
        log.error(f"MarkDuplicates failed for {name}"); sys.exit(1)
    bam.unlink(missing_ok=True)

    # Index
    run_cmd(["samtools", "index", str(dedup_bam)], log_file)
    return dedup_bam


def flagstat(bam: Path, out_dir: Path, log_file: Path) -> dict:
    out = out_dir / f"{bam.stem}_flagstat.txt"
    with open(out, "w") as f:
        r = subprocess.run(["samtools", "flagstat", str(bam)],
                           stdout=f, stderr=subprocess.PIPE, text=True)
    # Parse mapped %
    stats = {}
    try:
        text = out.read_text()
        for line in text.splitlines():
            if "mapped (" in line:
                pct = line.split("(")[1].split("%")[0]
                stats["mapped_pct"] = float(pct)
                break
        log.info(f"  flagstat: {stats.get('mapped_pct', '?')}% reads mapped → {out.name}")
    except Exception:
        pass
    return stats


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--sample-sheet", type=Path, required=True)
    parser.add_argument("--reference",    type=Path, required=True,
                        help="Reference genome FASTA (will be indexed if no .bwt file found)")
    parser.add_argument("--output-dir",   type=Path, default=Path("variant_out/alignment"))
    parser.add_argument("--threads",      type=int, default=4)
    parser.add_argument("--skip-markdup", action="store_true",
                        help="Skip duplicate marking (for amplicon / UMI data)")
    parser.add_argument("--execute",      action="store_true")
    args = parser.parse_args()

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_dir  = args.output_dir.with_name(f"{args.output_dir.name}_{timestamp}")
    log_file = out_dir / "run.log"

    if not args.execute:
        log.info("=== DRY RUN (pass --execute to run) ===")
        log.info(f"Would write outputs to: {out_dir}")
        samples = load_samples(args.sample_sheet)
        missing = [t for t in ["bwa-mem2", "samtools", "gatk"]
                   if not shutil.which(t)]
        if missing:
            log.info(f"[DRY RUN] Note: not installed: {', '.join(missing)}")
        log.info(f"[DRY RUN] Would align {len(samples)} sample(s) to {args.reference.name}")
        log.info("[DRY RUN] Would produce: <sample>.markdup.bam + .bai + flagstat.txt")
        log.info("[DRY RUN] Next step: run_germline_variants.py or run_somatic_variants.py")
        log.info("=== DRY RUN complete ===")
        return

    for t in ["bwa-mem2", "samtools", "gatk"]:
        require_tool(t)

    if not args.reference.exists():
        log.error(f"Reference not found: {args.reference}"); sys.exit(1)

    out_dir.mkdir(parents=True, exist_ok=True)
    log.info(f"Output directory: {out_dir}")

    samples = load_samples(args.sample_sheet)

    # Build BWA index if needed
    index_flag = Path(str(args.reference) + ".bwt.2bit.64")
    if not index_flag.exists():
        build_bwa_index(args.reference, log_file)

    bam_map = {}
    t0 = datetime.now()
    for s in samples:
        bam = align_sample(s, args.reference, out_dir, args.threads, log_file)
        if not args.skip_markdup:
            bam = mark_duplicates(bam, out_dir, log_file)
        stats = flagstat(bam, out_dir, log_file)
        bam_map[s["sample_name"]] = {
            "bam": str(bam), "bai": str(bam.with_suffix(".bam.bai")),
            "sample_type": s["sample_type"], **stats
        }
        log.info(f"  Done: {bam.name}")

    elapsed = (datetime.now() - t0).total_seconds()
    log.info(f"Alignment complete in {elapsed:.1f}s")

    manifest = {
        "pipeline": "alignment", "timestamp": timestamp,
        "reference": str(args.reference.resolve()),
        "threads": args.threads, "skip_markdup": args.skip_markdup,
        "samples": bam_map,
    }
    (out_dir / "manifest.json").write_text(json.dumps(manifest, indent=2))
    log.info(f"Manifest: {out_dir}/manifest.json")
    log.info(f"Done. BAMs in: {out_dir}")
    log.info("Next: run_germline_variants.py --bam-manifest variant_out/alignment_<ts>/manifest.json")


if __name__ == "__main__":
    main()
