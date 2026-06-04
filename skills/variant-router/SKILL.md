---
name: variant-router
description: Use when the user mentions "variant calling", "SNPs", "mutations", "somatic", "germline", "GATK", "Mutect2", "HaplotypeCaller", "BWA", "alignment", "VCF", "tumor-normal", "tumor mutational burden", "TMB", "WES", "whole exome", "whole genome sequencing", "DNA sequencing", or wants to find mutations in DNA data.
version: 1.0.0
---

# Variant Calling Router

## Two Routing Decisions

**Q1: Germline or somatic?**
- **Germline** → inherited variants, constitutional DNA, matched cohort → `run_germline_variants.py`
- **Somatic** → tumor mutations, cancer-specific changes, tumor-normal pairs → `run_somatic_variants.py`
- **Unsure?** Ask: "Is this cancer/tumor data or normal tissue?"

**Q2: Where in the pipeline?**

| User has | Route to |
|----------|----------|
| Raw FASTQs | Step 1: `run_alignment.py` |
| BAM files | Step 2: germline or somatic calling |
| VCF file | Step 3: `run_variant_annotation.py` |

## Always Dry-Run First

```bash
python scripts/run_alignment.py [args]           # dry run
python scripts/run_alignment.py [args] --execute # execute
```

## Pipeline Overview

```
FASTQ → run_alignment.py → BAMs
                             │
              ┌──────────────┴──────────────┐
              ▼                             ▼
  run_germline_variants.py    run_somatic_variants.py
  (HaplotypeCaller → GVCF)   (Mutect2 → contamination)
              │                             │
              └──────────────┬──────────────┘
                             ▼
              run_variant_annotation.py
              (SnpEff + bcftools stats)
```
