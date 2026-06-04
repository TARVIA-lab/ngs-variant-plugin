---
name: variant-germline
description: Use when the user wants to call germline variants, run HaplotypeCaller, do joint genotyping, find inherited mutations, process constitutional DNA, or call SNPs and indels from normal tissue.
version: 1.0.0
---

# GATK Germline Variant Calling

## Pipeline

1. **HaplotypeCaller** → per-sample GVCF (joint calling mode)
2. **CombineGVCFs** → merged GVCF (multi-sample cohorts)
3. **GenotypeGVCFs** → raw genotyped VCF
4. **VariantFiltration** → hard-filtered VCF (SNPs and indels filtered separately)

## Hard Filter Thresholds (GATK Best Practices)

**SNPs**: QD < 2.0, FS > 60.0, MQ < 40.0, SOR > 3.0
**Indels**: QD < 2.0, FS > 200.0, SOR > 10.0

For cohorts > 30 samples, consider VQSR instead of hard filtering.

## Usage

```bash
python scripts/run_germline_variants.py \
  --bam-manifest variant_out/alignment_<ts>/manifest.json \
  --reference /refs/GRCh38/genome.fa \
  --known-sites /refs/dbsnp_138.hg38.vcf.gz \
  --intervals chr1-chr22 \
  --execute
```

## BQSR Note

Base Quality Score Recalibration (BQSR) is not included in this script to keep the pipeline lightweight. For production runs, add BQSR before HaplotypeCaller using `gatk BaseRecalibrator` + `gatk ApplyBQSR`.
