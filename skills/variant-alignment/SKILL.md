---
name: variant-alignment
description: Use when the user wants to align DNA reads to a reference genome, run BWA-MEM2, create BAM files, mark duplicates, or prepare BAMs for variant calling.
version: 1.0.0
---

# BWA-MEM2 Alignment

## Key Points

- **Read groups are required** by GATK — set automatically from sample name
- **Duplicate marking** is done with GATK MarkDuplicatesSpark — skip with `--skip-markdup` for amplicon/UMI data
- **Reference must be indexed** — done automatically if `.bwt.2bit.64` file is missing

## Usage

```bash
python scripts/run_alignment.py \
  --sample-sheet examples/germline_samplesheet.tsv \
  --reference /refs/GRCh38/genome.fa \
  --threads 8 --execute
```

## Checking Alignment Quality

Open `<sample>_flagstat.txt` — key metrics:
- **Mapped %** should be > 95% for WGS/WES
- **Properly paired %** should be > 90%
- **Duplicate %** for WES is typically 15–30%; > 50% suggests low library complexity

## Output Used by Next Steps

The `manifest.json` from alignment is passed directly to the variant callers:
```bash
python scripts/run_germline_variants.py \
  --bam-manifest variant_out/alignment_<ts>/manifest.json ...
```
