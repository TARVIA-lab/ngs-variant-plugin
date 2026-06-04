---
name: variant-somatic
description: Use when the user wants to call somatic mutations, run Mutect2, find tumor-specific variants, calculate tumor mutational burden (TMB), run tumor-normal variant calling, detect cancer driver mutations, or process oncology sequencing data.
version: 1.0.0
---

# GATK Mutect2 Somatic Variant Calling

## Tumor-Normal vs Tumor-Only

| Mode | Sensitivity | Specificity | When to use |
|------|------------|-------------|-------------|
| Tumor-normal | High | Highest | Matched normal available (preferred) |
| Tumor-only | Higher | Lower | No matched normal; use PoN to reduce noise |

## Somatic Sample Sheet Format

```
sample_name	bam	sample_type	pair_id
TUMOR_1	/path/TUMOR_1.markdup.bam	tumor	PAIR_1
NORMAL_1	/path/NORMAL_1.markdup.bam	normal	PAIR_1
```

`pair_id` links tumor to its matched normal. Multiple pairs can be in one sheet.

## Recommended Resources

- `--germline-resource` af-only-gnomad: helps distinguish somatic from germline
- `--panel-of-normals` 1000g_pon: removes systematic artifacts
- `--common-sites` small_exac_common: needed for contamination estimate

## Key Output Files

- `<pair>.somatic.filtered.vcf.gz` — PASS variants are high-confidence somatic calls
- `<pair>.contamination.table` — estimated cross-sample contamination %
- `<pair>.f1r2.tar.gz` — orientation bias model (can be used with LearnReadOrientationModel)
