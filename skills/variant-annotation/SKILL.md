---
name: variant-annotation
description: Use when the user wants to annotate variants, add functional consequence to a VCF, run SnpEff, add dbSNP rs IDs, calculate Ti/Tv ratio, get variant statistics, or produce a readable variant table from a VCF file.
version: 1.0.0
---

# Variant Annotation

## Two Annotation Layers

1. **bcftools stats** — always runs: Ti/Tv ratio, SNP/indel counts, per-sample stats
2. **SnpEff** — functional consequence: missense, nonsense, splice site, synonymous, intergenic

## SnpEff Genome Aliases

| Alias | SnpEff database |
|-------|----------------|
| `hg38` / `GRCh38` | `GRCh38.mane.1.2.ensembl` |
| `hg19` / `GRCh37` | `GRCh37.87` |
| `mm39` | `GRCm39.105` |
| `mm10` | `GRCm38.99` |

## Usage

```bash
python scripts/run_variant_annotation.py \
  --vcf variant_out/germline_<ts>/germline.filtered.vcf.gz \
  --genome GRCh38 \
  --execute
```

## Interpreting Ti/Tv Ratio

- **WGS germline**: expected ~2.0–2.1
- **WES germline**: expected ~2.8–3.0 (exons are more constrained)
- **Somatic (cancer)**: can vary widely; COSMIC signatures may help
- Ti/Tv < 1.8 often indicates a quality problem

## `variants_summary.tsv`

A PASS-filtered variant table for quick inspection — up to 500 rows, no INFO parsing needed.
