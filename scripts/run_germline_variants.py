#!/usr/bin/env python3
"""
run_germline_variants.py — GATK4 germline variant calling (HaplotypeCaller → GenotypeGVCFs).

Pipeline:
  1. GATK HaplotypeCaller  →  per-sample GVCF
  2. GATK CombineGVCFs     →  multi-sample GVCF  (if > 1 sample)
  3. GATK GenotypeGVCFs    →  raw genotyped VCF
  4. GATK VariantFiltration →  hard-filtered VCF (SNPs + indels separately)
  5. bcftools stats         →  variant QC summary

Accepts either a BAM manifest (from run_alignment.py) or a direct --bam-list TSV.

Usage:
    python run_germline_variants.py \\
        --bam-manifest variant_out/alignment_<ts>/manifest.json \\
        --reference ref/genome.fa --execute

    python run_germline_variants.py \\
        --bam-list bam_list.tsv --reference ref/genome.fa \\
        --known-sites dbsnp.vcf.gz --execute
"""
import argparse, json, logging, shutil, subprocess, sys
from datetime import datetime
from pathlib import Path

logging.basicConfig(format="%(asctime)s [%(levelname)s] %(message)s",
                    datefmt="%H:%M:%S", level=logging.INFO)
log = logging.getLogger(__name__)


def require_tool(name):
    if not shutil.which(name):
        log.error(f"Required: {name}  →  conda install -c bioconda gatk4")
        sys.exit(1)


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


def load_bams(manifest_path=None, bam_list_path=None) -> list[dict]:
    if manifest_path:
        data = json.loads(manifest_path.read_text())
        return [{"name": n, "bam": Path(v["bam"])}
                for n, v in data["samples"].items()]
    # TSV: sample_name, bam
    samples = []
    import csv
    with open(bam_list_path) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            row = {k.strip(): v.strip() for k, v in row.items()}
            samples.append({"name": row["sample_name"], "bam": Path(row["bam"])})
    return samples


def haplotype_caller(sample: dict, ref: Path, out_dir: Path,
                     known_sites: list, intervals: str, log_file: Path) -> Path:
    gvcf = out_dir / f"{sample['name']}.g.vcf.gz"
    cmd = [
        "gatk", "HaplotypeCaller",
        "-R", str(ref),
        "-I", str(sample["bam"]),
        "-O", str(gvcf),
        "-ERC", "GVCF",
        "--sample-name", sample["name"],
    ]
    for ks in known_sites:
        cmd += ["--dbsnp", str(ks)]
    if intervals:
        cmd += ["-L", intervals]

    rc = run_cmd(cmd, log_file)
    if rc != 0:
        log.error(f"HaplotypeCaller failed for {sample['name']}"); sys.exit(1)
    log.info(f"  GVCF: {gvcf.name}")
    return gvcf


def combine_and_genotype(gvcfs: list[Path], ref: Path, out_dir: Path,
                          log_file: Path) -> Path:
    combined = out_dir / "combined.g.vcf.gz"
    raw_vcf  = out_dir / "raw_variants.vcf.gz"

    if len(gvcfs) > 1:
        log.info("Combining GVCFs...")
        cmd = ["gatk", "CombineGVCFs", "-R", str(ref), "-O", str(combined)]
        for g in gvcfs:
            cmd += ["-V", str(g)]
        rc = run_cmd(cmd, log_file)
        if rc != 0:
            log.error("CombineGVCFs failed"); sys.exit(1)
        joint_input = combined
    else:
        joint_input = gvcfs[0]

    log.info("Genotyping GVCFs...")
    rc = run_cmd([
        "gatk", "GenotypeGVCFs",
        "-R", str(ref),
        "-V", str(joint_input),
        "-O", str(raw_vcf),
    ], log_file)
    if rc != 0:
        log.error("GenotypeGVCFs failed"); sys.exit(1)
    log.info(f"Raw VCF: {raw_vcf.name}")
    return raw_vcf


def hard_filter(raw_vcf: Path, ref: Path, out_dir: Path, log_file: Path) -> Path:
    snp_vcf   = out_dir / "snps.filtered.vcf.gz"
    indel_vcf = out_dir / "indels.filtered.vcf.gz"
    final_vcf = out_dir / "germline.filtered.vcf.gz"

    # Select + filter SNPs
    run_cmd(["gatk", "SelectVariants", "-R", str(ref),
             "-V", str(raw_vcf), "--select-type-to-include", "SNP",
             "-O", str(out_dir / "snps.vcf.gz")], log_file)
    run_cmd([
        "gatk", "VariantFiltration",
        "-R", str(ref), "-V", str(out_dir / "snps.vcf.gz"),
        "--filter-expression", "QD < 2.0",     "--filter-name", "QD2",
        "--filter-expression", "FS > 60.0",    "--filter-name", "FS60",
        "--filter-expression", "MQ < 40.0",    "--filter-name", "MQ40",
        "--filter-expression", "SOR > 3.0",    "--filter-name", "SOR3",
        "-O", str(snp_vcf),
    ], log_file)

    # Select + filter Indels
    run_cmd(["gatk", "SelectVariants", "-R", str(ref),
             "-V", str(raw_vcf), "--select-type-to-include", "INDEL",
             "-O", str(out_dir / "indels.vcf.gz")], log_file)
    run_cmd([
        "gatk", "VariantFiltration",
        "-R", str(ref), "-V", str(out_dir / "indels.vcf.gz"),
        "--filter-expression", "QD < 2.0",   "--filter-name", "QD2",
        "--filter-expression", "FS > 200.0", "--filter-name", "FS200",
        "--filter-expression", "SOR > 10.0", "--filter-name", "SOR10",
        "-O", str(indel_vcf),
    ], log_file)

    # Merge back
    run_cmd(["gatk", "MergeVcfs",
             "-I", str(snp_vcf), "-I", str(indel_vcf),
             "-O", str(final_vcf)], log_file)
    log.info(f"Filtered VCF: {final_vcf.name}")
    return final_vcf


def vcf_stats(vcf: Path, out_dir: Path, log_file: Path) -> dict:
    stats_file = out_dir / f"{vcf.stem}_stats.txt"
    rc = run_cmd(["bcftools", "stats", str(vcf)], log_file)
    subprocess.run(["bcftools", "stats", str(vcf)],
                   stdout=open(stats_file, "w"), stderr=subprocess.DEVNULL)
    stats = {}
    try:
        for line in stats_file.read_text().splitlines():
            if line.startswith("SN") and "number of SNPs" in line:
                stats["n_snps"] = int(line.split("\t")[-1])
            if line.startswith("SN") and "number of indels" in line:
                stats["n_indels"] = int(line.split("\t")[-1])
        log.info(f"  Stats: {stats.get('n_snps','?')} SNPs, {stats.get('n_indels','?')} indels")
    except Exception:
        pass
    return stats


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    grp = parser.add_mutually_exclusive_group(required=True)
    grp.add_argument("--bam-manifest", type=Path, help="manifest.json from run_alignment.py")
    grp.add_argument("--bam-list",     type=Path, help="TSV: sample_name, bam")
    parser.add_argument("--reference",    type=Path, required=True)
    parser.add_argument("--known-sites",  type=Path, nargs="*", default=[],
                        help="Known variant VCFs for BQSR (optional but recommended)")
    parser.add_argument("--intervals",    help="Genomic intervals, e.g. chr22 or targets.bed")
    parser.add_argument("--output-dir",   type=Path, default=Path("variant_out/germline"))
    parser.add_argument("--execute",      action="store_true")
    args = parser.parse_args()

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_dir   = args.output_dir.with_name(f"{args.output_dir.name}_{timestamp}")
    log_file  = out_dir / "run.log"

    samples = load_bams(args.bam_manifest, args.bam_list)

    if not args.execute:
        log.info("=== DRY RUN (pass --execute to run) ===")
        log.info(f"Would write to: {out_dir}")
        log.info(f"Samples: {[s['name'] for s in samples]}")
        log.info("[DRY RUN] Would run: HaplotypeCaller → CombineGVCFs → GenotypeGVCFs → VariantFiltration")
        log.info("[DRY RUN] Output: germline.filtered.vcf.gz")
        log.info("=== DRY RUN complete ===")
        return

    for t in ["gatk", "bcftools"]:
        require_tool(t)
    if not args.reference.exists():
        log.error(f"Reference not found: {args.reference}"); sys.exit(1)

    out_dir.mkdir(parents=True, exist_ok=True)

    # samtools faidx if needed
    fai = Path(str(args.reference) + ".fai")
    if not fai.exists():
        log.info("Indexing reference with samtools faidx...")
        run_cmd(["samtools", "faidx", str(args.reference)], log_file)

    # GATK dict if needed
    ref_dict = args.reference.with_suffix(".dict")
    if not ref_dict.exists():
        log.info("Creating sequence dictionary...")
        run_cmd(["gatk", "CreateSequenceDictionary", "-R", str(args.reference)], log_file)
    log.info(f"Output: {out_dir}")

    log.info(f"HaplotypeCaller on {len(samples)} sample(s)...")
    gvcfs = []
    for s in samples:
        if not s["bam"].exists():
            log.error(f"BAM not found: {s['bam']}"); sys.exit(1)
        gvcfs.append(haplotype_caller(s, args.reference, out_dir,
                                       args.known_sites, args.intervals, log_file))

    final_vcf = combine_and_genotype(gvcfs, args.reference, out_dir, log_file)
    filtered  = hard_filter(final_vcf, args.reference, out_dir, log_file)
    stats     = vcf_stats(filtered, out_dir, log_file)

    manifest = {
        "pipeline": "germline", "timestamp": timestamp,
        "reference": str(args.reference.resolve()),
        "samples": [s["name"] for s in samples],
        "output_vcf": str(filtered),
        "stats": stats,
    }
    (out_dir / "manifest.json").write_text(json.dumps(manifest, indent=2))
    log.info(f"Done. VCF: {filtered}")
    log.info(f"Next: run_variant_annotation.py --vcf {filtered}")


if __name__ == "__main__":
    main()
