<div align="center">

# ngs-variant-plugin

**A Claude Code plugin for DNA variant calling — germline and somatic**

[![Python 3.10+](https://img.shields.io/badge/python-3.10%2B-blue?logo=python&logoColor=white)](https://www.python.org/)
[![BWA-MEM2](https://img.shields.io/badge/BWA--MEM2-2.2.1-green)](https://github.com/bwa-mem2/bwa-mem2)
[![GATK4](https://img.shields.io/badge/GATK4-4.6.2-orange)](https://gatk.broadinstitute.org/)
[![bcftools](https://img.shields.io/badge/bcftools-1.22-teal)](https://samtools.github.io/bcftools/)
[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](LICENSE)
[![Claude Code Plugin](https://img.shields.io/badge/Claude%20Code-Plugin-blueviolet?logo=anthropic)](https://claude.ai/code)
[![TARVIA-lab](https://img.shields.io/badge/TARVIA--lab-GitHub-black?logo=github)](https://github.com/TARVIA-lab)

Routes DNA sequencing data from FASTQ through BWA-MEM2 alignment, GATK4 germline or somatic variant calling, and bcftools/SnpEff annotation — locally, with open-source tools, via natural language in Claude Code.

[Quick Start](#quick-start) · [Pipelines](#pipeline) · [Claude Code Plugin](#claude-code-plugin-installation) · [Standalone Usage](#standalone-usage) · [Output Format](#output-format)

</div>

---

## Overview

`ngs-variant-plugin` is a TARVIA-lab tool for AI-assisted DNA variant discovery. Designed for oncology workflows — tumor mutational burden, somatic driver mutations, germline cancer predisposition — with full support for both tumor-normal paired and tumor-only Mutect2 modes.

```
"I have WES FASTQs from 5 tumor-normal pairs. Call somatic mutations."
```

Claude will:
1. Align reads with BWA-MEM2, mark duplicates with GATK
2. Run GATK Mutect2 per pair → contamination estimate → FilterMutectCalls
3. Annotate with SnpEff (functional consequence) and bcftools (Ti/Tv, stats)
4. Produce a clean `variants_summary.tsv` of PASS variants

Every step is a dry-run by default. Nothing executes until you confirm.

---

## Pipeline

```
FASTQ (R1 + R2)
      │
      ▼
┌────────────┐   ┌──────────────────────────────────────┐   ┌─────────────┐
│ Alignment  │──▶│        Variant Calling                │──▶│ Annotation  │
│ BWA-MEM2   │   │  ┌─────────────────────────────────┐  │   │ SnpEff +    │
│ + samtools │   │  │ Germline: HaplotypeCaller →      │  │   │ bcftools    │
│ sort       │   │  │ GenotypeGVCFs → hard filter      │  │   │ stats       │
│ + GATK     │   │  ├─────────────────────────────────┤  │   └─────────────┘
│ MarkDup    │   │  │ Somatic:  Mutect2 → Contamination│  │         │
└────────────┘   │  │ estimate → FilterMutectCalls     │  │   variants_summary.tsv
      │          │  └─────────────────────────────────┘  │   snpeff_summary.html
 *.markdup.bam   └──────────────────────────────────────┘
 flagstat.txt               │
                  germline.filtered.vcf.gz
                  <pair>.somatic.filtered.vcf.gz
```

| Step | Script | Tools | Output |
|------|--------|-------|--------|
| **1. Alignment** | `run_alignment.py` | BWA-MEM2, samtools, GATK | `*.markdup.bam` + flagstat |
| **2a. Germline** | `run_germline_variants.py` | GATK HaplotypeCaller | `germline.filtered.vcf.gz` |
| **2b. Somatic** | `run_somatic_variants.py` | GATK Mutect2 | `<pair>.somatic.filtered.vcf.gz` |
| **3. Annotation** | `run_variant_annotation.py` | SnpEff, bcftools | Annotated VCF + summary TSV |

---

## Claude Code Plugin Installation

### 1. Clone

```bash
git clone https://github.com/TARVIA-lab/ngs-variant-plugin.git
```

### 2. Register local marketplace

```
/plugin add-marketplace /path/to/ngs-variant-plugin/.agents/plugins/marketplace.json
```

### 3. Install

```
/plugin install ngs-variant
```

### Talking to Claude

```
"Align my WES samples to GRCh38: samplesheet.tsv"

"Call germline variants on my 3 BAMs — output a filtered VCF"

"I have tumor-normal pairs. Run Mutect2 and filter somatic calls."

"Annotate this VCF with SnpEff for human GRCh38"

"/variant-preflight --pipeline somatic"
```

---

## Standalone Usage

### Prerequisites

```bash
conda install -c bioconda -c conda-forge bwa-mem2 samtools gatk4 bcftools
pip install multiqc
conda install -c bioconda snpeff   # optional, for functional annotation
```

---

### Step 1 — Alignment

```bash
python scripts/run_alignment.py \
  --sample-sheet examples/germline_samplesheet.tsv \
  --reference /refs/GRCh38/genome.fa \
  --threads 8 \
  --execute
```

Produces `*.markdup.bam` + per-sample `flagstat.txt`. Pass `--skip-markdup` for amplicon data.

---

### Step 2a — Germline Variants

```bash
python scripts/run_germline_variants.py \
  --bam-manifest variant_out/alignment_<ts>/manifest.json \
  --reference /refs/GRCh38/genome.fa \
  --known-sites /refs/dbsnp_138.hg38.vcf.gz \
  --execute
```

Pass `--intervals chr22` or `--intervals targets.bed` to restrict to regions.

---

### Step 2b — Somatic Variants (tumor-normal)

```bash
python scripts/run_somatic_variants.py \
  --sample-sheet examples/somatic_samplesheet.tsv \
  --reference /refs/GRCh38/genome.fa \
  --germline-resource /refs/af-only-gnomad.hg38.vcf.gz \
  --panel-of-normals /refs/1000g_pon.hg38.vcf.gz \
  --common-sites /refs/small_exac_common_3.hg38.vcf.gz \
  --execute
```

Tumor-only (no matched normal):

```bash
python scripts/run_somatic_variants.py \
  --sample-sheet tumor_only.tsv \
  --reference /refs/GRCh38/genome.fa \
  --execute
```

---

### Step 3 — Annotation

```bash
# Full annotation (SnpEff functional consequence + stats)
python scripts/run_variant_annotation.py \
  --vcf variant_out/germline_<ts>/germline.filtered.vcf.gz \
  --genome GRCh38 \
  --execute

# Stats + dbSNP rs IDs only (no SnpEff)
python scripts/run_variant_annotation.py \
  --vcf somatic.filtered.vcf.gz \
  --dbsnp /refs/dbsnp_138.hg38.vcf.gz \
  --skip-snpeff \
  --execute
```

Supported `--genome` aliases: `hg38` / `GRCh38`, `hg19` / `GRCh37`, `mm39`, `mm10`.

---

## Input File Formats

### `germline_samplesheet.tsv`

```
sample_name	fastq_r1	fastq_r2	sample_type
SAMPLE_1	/data/SAMPLE_1_R1.fastq.gz	/data/SAMPLE_1_R2.fastq.gz	germline
SAMPLE_2	/data/SAMPLE_2_R1.fastq.gz	/data/SAMPLE_2_R2.fastq.gz	germline
```

### `somatic_samplesheet.tsv`

```
sample_name	bam	sample_type	pair_id
TUMOR_1	variant_out/alignment_<ts>/TUMOR_1.markdup.bam	tumor	PAIR_1
NORMAL_1	variant_out/alignment_<ts>/NORMAL_1.markdup.bam	normal	PAIR_1
```

---

## Output Format

```
variant_out/
├── alignment_<timestamp>/
│   ├── manifest.json
│   ├── SAMPLE_1.markdup.bam      ◀ coordinate-sorted, duplicate-marked
│   ├── SAMPLE_1.markdup.bam.bai
│   ├── SAMPLE_1.markdup_flagstat.txt
│   └── run.log
│
├── germline_<timestamp>/
│   ├── manifest.json
│   ├── SAMPLE_1.g.vcf.gz         # per-sample GVCF
│   ├── combined.g.vcf.gz         # merged GVCF
│   ├── raw_variants.vcf.gz       # pre-filter
│   ├── germline.filtered.vcf.gz  ◀ final hard-filtered VCF
│   └── run.log
│
├── somatic_<timestamp>/
│   ├── manifest.json
│   ├── PAIR_1.mutect2.vcf.gz     # raw Mutect2 calls
│   ├── PAIR_1.contamination.table
│   ├── PAIR_1.somatic.filtered.vcf.gz  ◀ final filtered somatic VCF
│   └── run.log
│
└── annotation_<timestamp>/
    ├── manifest.json
    ├── germline.filtered.snpeff.vcf.gz  ◀ annotated VCF
    ├── variants_summary.tsv             ◀ PASS variants as readable table
    ├── snpeff_summary.html              ◀ SnpEff HTML report
    └── run.log
```

### `variants_summary.tsv` columns

| Column | Description |
|--------|-------------|
| CHROM | Chromosome |
| POS | 1-based position |
| ID | rs ID (if dbSNP annotated) |
| REF | Reference allele |
| ALT | Alternate allele(s) |
| QUAL | Phred-scaled quality |
| FILTER | PASS or filter tag(s) |

---

## Running the Test Suite

`test_data/` contains a synthetic 80 KB reference genome and 3 samples with known SNPs:

```bash
# Step 1: Alignment
python scripts/run_alignment.py \
  --sample-sheet test_data/germline_samplesheet.tsv \
  --reference test_data/ref/test_genome.fa \
  --threads 4 --execute

# Step 2: Germline variant calling
python scripts/run_germline_variants.py \
  --bam-manifest variant_out/alignment_<ts>/manifest.json \
  --reference test_data/ref/test_genome.fa \
  --execute

# Step 3: Annotation (stats only — SnpEff requires a real genome DB)
python scripts/run_variant_annotation.py \
  --vcf variant_out/germline_<ts>/germline.filtered.vcf.gz \
  --skip-snpeff --execute
```

Expected output: **~584 SNPs** called across 3 samples, **100% mapping rate** per sample.

---

## Project Structure

```
ngs-variant-plugin/
├── .agents/plugins/marketplace.json
├── .claude-plugin/plugin.json
│
├── skills/
│   ├── variant-router/SKILL.md       # Routes to germline vs somatic
│   ├── variant-alignment/SKILL.md    # BWA-MEM2 + GATK MarkDup guidance
│   ├── variant-germline/SKILL.md     # HaplotypeCaller guidance
│   ├── variant-somatic/SKILL.md      # Mutect2 + contamination guidance
│   └── variant-annotation/SKILL.md  # SnpEff + bcftools guidance
│
├── scripts/
│   ├── variant_preflight.py          # Tool checker
│   ├── run_alignment.py              # BWA-MEM2 + samtools + MarkDup
│   ├── run_germline_variants.py      # HaplotypeCaller → GenotypeGVCFs → filter
│   ├── run_somatic_variants.py       # Mutect2 → contamination → filter
│   └── run_variant_annotation.py    # SnpEff + bcftools stats + summary TSV
│
├── references/
│   ├── pipeline-registry.json
│   └── tool-registry.json
│
├── examples/
│   ├── germline_samplesheet.tsv
│   └── somatic_samplesheet.tsv
│
└── test_data/
    ├── ref/test_genome.fa            # Synthetic 80 KB reference
    ├── fastq/                        # 3 samples × 8,000 read pairs
    └── germline_samplesheet.tsv
```

---

## Reference Downloads (GRCh38)

```bash
# Genome FASTA
wget https://ftp.ensembl.org/pub/release-111/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

# dbSNP (for known sites + rs ID annotation)
wget https://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/common_all_20180418.vcf.gz

# gnomAD AF-only (for Mutect2 germline resource)
wget https://storage.googleapis.com/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz

# 1000 Genomes Panel of Normals (for Mutect2 PoN)
wget https://storage.googleapis.com/gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz
```

---

## Tool Versions Tested

| Tool | Version | Install |
|------|---------|---------|
| BWA-MEM2 | 2.2.1 | `conda install -c bioconda bwa-mem2` |
| samtools | 1.22.1 | `conda install -c bioconda samtools` |
| GATK4 | 4.6.2.0 | `conda install -c bioconda gatk4` |
| bcftools | 1.22 | `conda install -c bioconda bcftools` |
| SnpEff | 5.2 | `conda install -c bioconda snpeff` |
| Python | 3.10 – 3.14 | — |

---

## Related Work

| Repo | Description |
|------|-------------|
| [ngs-rnaseq-plugin](https://github.com/TARVIA-lab/ngs-rnaseq-plugin) | Bulk RNA-seq: FastQC → Salmon → DESeq2 |
| [ngs-scrna-plugin](https://github.com/TARVIA-lab/ngs-scrna-plugin) | scRNA-seq: STARsolo → Scanpy → UMAP/Leiden |
| [Benchmarking-LLM-Scientific-Reasoning-in-Oncology](https://github.com/TARVIA-lab/Benchmarking-Large-Language-Model-Scientific-Reasoning-in-Oncology) | LLM benchmarks in oncology |
| **ngs-variant-plugin** | This repo — DNA variant calling |

---

## License

[Apache License 2.0](LICENSE)

---

## Acknowledgments

Built on [GATK4](https://gatk.broadinstitute.org/), [BWA-MEM2](https://github.com/bwa-mem2/bwa-mem2), and [SnpEff](https://pcingola.github.io/SnpEff/). Part of the TARVIA-lab suite of AI-assisted genomics tools.
