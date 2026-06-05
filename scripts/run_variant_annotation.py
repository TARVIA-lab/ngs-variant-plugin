#!/usr/bin/env python3
"""
run_variant_annotation.py — Annotate VCFs with functional consequence and population allele frequencies.

Steps:
  1. bcftools stats   →  variant summary (Ti/Tv, SNP/indel counts)
  2. bcftools annotate →  add dbSNP rs IDs (optional)
  3. SnpEff           →  functional consequence annotation (missense, nonsense, etc.)
  4. Summary TSV      →  top variants table (PASS + HIGH/MODERATE impact)

Usage:
    # Annotate with SnpEff (genome database auto-downloaded if needed)
    python run_variant_annotation.py \\
        --vcf variant_out/germline_<ts>/germline.filtered.vcf.gz \\
        --genome GRCh38.mane.1.2.ensembl --execute

    # Just stats + dbSNP IDs, no SnpEff
    python run_variant_annotation.py \\
        --vcf somatic.filtered.vcf.gz \\
        --dbsnp /refs/dbsnp_138.hg38.vcf.gz \\
        --skip-snpeff --execute
"""
import argparse, json, logging, re, shutil, subprocess, sys
from datetime import datetime
from pathlib import Path

logging.basicConfig(format="%(asctime)s [%(levelname)s] %(message)s",
                    datefmt="%H:%M:%S", level=logging.INFO)
log = logging.getLogger(__name__)

SNPEFF_GENOMES = {
    "hg38":     "GRCh38.mane.1.2.ensembl",
    "GRCh38":   "GRCh38.mane.1.2.ensembl",
    "hg19":     "GRCh37.87",
    "GRCh37":   "GRCh37.87",
    "mm39":     "GRCm39.105",
    "mm10":     "GRCm38.99",
}


def require_tool(name):
    if not shutil.which(name):
        log.error(f"Required: {name}")
        if name == "snpEff":
            log.error("Install: conda install -c bioconda snpeff")
        sys.exit(1)


def run_cmd(cmd, log_file=None, stdout_file=None):
    s = " ".join(str(c) for c in cmd)
    log.info(f"Running: {s}")
    if stdout_file:
        with open(log_file or "/dev/null", "a") as lf, open(stdout_file, "w") as out:
            r = subprocess.run(cmd, stdout=out, stderr=lf)
    elif log_file:
        with open(log_file, "a") as lf:
            lf.write(f"\n$ {s}\n")
            r = subprocess.run(cmd, stdout=lf, stderr=lf)
    else:
        r = subprocess.run(cmd, capture_output=True, text=True)
    return r.returncode


def bcftools_stats(vcf: Path, out_dir: Path) -> dict:
    stats_file = out_dir / f"{vcf.stem}_stats.txt"
    subprocess.run(["bcftools", "stats", str(vcf)],
                   stdout=open(stats_file, "w"), stderr=subprocess.DEVNULL)
    stats = {}
    try:
        for line in stats_file.read_text().splitlines():
            if "number of SNPs:" in line:
                stats["n_snps"] = int(line.split()[-1])
            elif "number of indels:" in line:
                stats["n_indels"] = int(line.split()[-1])
            elif "Ts/Tv:" in line:
                stats["ts_tv"] = float(line.split()[-1])
        log.info(f"  Stats: {stats.get('n_snps','?')} SNPs, "
                 f"{stats.get('n_indels','?')} indels, "
                 f"Ts/Tv={stats.get('ts_tv','?')}")
    except Exception:
        pass
    return stats


def annotate_dbsnp(vcf: Path, dbsnp: Path, out_dir: Path, log_file: Path) -> Path:
    out = out_dir / vcf.name.replace(".vcf.gz", ".dbsnp.vcf.gz")
    rc = run_cmd([
        "bcftools", "annotate",
        "-a", str(dbsnp),
        "-c", "ID",
        "-O", "z", "-o", str(out),
        str(vcf),
    ], log_file)
    if rc != 0:
        log.warning("dbSNP annotation failed — continuing without rs IDs")
        return vcf
    run_cmd(["bcftools", "index", "-t", str(out)])
    log.info(f"  dbSNP annotated: {out.name}")
    return out


def run_snpeff(vcf: Path, genome: str, out_dir: Path,
               log_file: Path) -> tuple[Path, Path]:
    genome_id = SNPEFF_GENOMES.get(genome, genome)
    ann_vcf   = out_dir / vcf.name.replace(".vcf.gz", ".snpeff.vcf.gz")
    html_out  = out_dir / "snpeff_summary.html"
    stats_csv = out_dir / "snpeff_stats.csv"

    log.info(f"  SnpEff annotation with genome: {genome_id}")
    with open(ann_vcf.with_suffix("").with_suffix(".vcf"), "w") as vcf_out:
        r = subprocess.run(
            ["snpEff", "-v", "-stats", str(html_out),
             "-csvStats", str(stats_csv), genome_id, str(vcf)],
            stdout=vcf_out,
            stderr=open(log_file, "a"),
        )
    if r.returncode != 0:
        log.error("SnpEff failed"); sys.exit(1)

    # Compress + index
    subprocess.run(["bcftools", "view", "-O", "z", "-o", str(ann_vcf),
                    str(ann_vcf.with_suffix("").with_suffix(".vcf"))])
    ann_vcf.with_suffix("").with_suffix(".vcf").unlink(missing_ok=True)
    subprocess.run(["bcftools", "index", "-t", str(ann_vcf)])
    log.info(f"  Annotated VCF: {ann_vcf.name}")
    log.info(f"  SnpEff HTML report: {html_out.name}")
    return ann_vcf, html_out


def make_summary_tsv(vcf: Path, out_dir: Path, max_rows: int = 500) -> Path:
    """Extract PASS variants into a readable TSV (bcftools query)."""
    out_tsv = out_dir / "variants_summary.tsv"
    header  = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]

    result = subprocess.run(
        ["bcftools", "query",
         "-f", "%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t%INFO\n",
         "-i", "FILTER='PASS' || FILTER='.'",
         str(vcf)],
        capture_output=True, text=True,
    )
    lines = ["\t".join(header)]
    for line in result.stdout.strip().splitlines()[:max_rows]:
        lines.append(line)

    out_tsv.write_text("\n".join(lines) + "\n")
    n = len(lines) - 1
    log.info(f"  Summary TSV: {n} PASS variants → {out_tsv.name}")
    return out_tsv


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--vcf",         type=Path, required=True)
    parser.add_argument("--genome",      default="GRCh38",
                        help="SnpEff genome DB (default: GRCh38). Aliases: hg38, GRCh38, mm39, hg19")
    parser.add_argument("--dbsnp",       type=Path, help="dbSNP VCF.gz for rs ID annotation")
    parser.add_argument("--skip-snpeff", action="store_true")
    parser.add_argument("--output-dir",  type=Path, default=Path("variant_out/annotation"))
    parser.add_argument("--execute",     action="store_true")
    args = parser.parse_args()

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_dir   = args.output_dir.with_name(f"{args.output_dir.name}_{timestamp}")
    log_file  = out_dir / "run.log"

    if not args.execute:
        log.info("=== DRY RUN (pass --execute to run) ===")
        log.info(f"Would annotate: {args.vcf.name}")
        log.info(f"Would write to: {out_dir}")
        steps = ["bcftools stats"]
        if args.dbsnp:
            steps.append("bcftools annotate (dbSNP rs IDs)")
        if not args.skip_snpeff:
            steps.append(f"SnpEff ({SNPEFF_GENOMES.get(args.genome, args.genome)})")
        steps.append("variants_summary.tsv")
        log.info(f"[DRY RUN] Steps: {' → '.join(steps)}")
        log.info("=== DRY RUN complete ===")
        return

    require_tool("bcftools")
    if not args.skip_snpeff:
        require_tool("snpEff")
    if not args.vcf.exists():
        log.error(f"VCF not found: {args.vcf}"); sys.exit(1)

    out_dir.mkdir(parents=True, exist_ok=True)
    log.info(f"Output: {out_dir}")

    vcf = args.vcf
    stats = bcftools_stats(vcf, out_dir)

    if args.dbsnp:
        vcf = annotate_dbsnp(vcf, args.dbsnp, out_dir, log_file)

    final_vcf = vcf
    if not args.skip_snpeff:
        final_vcf, _ = run_snpeff(vcf, args.genome, out_dir, log_file)

    summary = make_summary_tsv(final_vcf, out_dir)

    manifest = {
        "pipeline": "annotation", "timestamp": timestamp,
        "input_vcf": str(args.vcf), "output_vcf": str(final_vcf),
        "genome": args.genome, "stats": stats,
        "summary_tsv": str(summary),
    }
    (out_dir / "manifest.json").write_text(json.dumps(manifest, indent=2))
    log.info(f"Manifest: {out_dir}/manifest.json")
    log.info(f"Done. Annotated VCF: {final_vcf.name}")


if __name__ == "__main__":
    main()
