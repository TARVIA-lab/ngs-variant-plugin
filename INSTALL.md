# Installation & Setup Guide

## Python Dependencies (2 minutes)

```bash
git clone https://github.com/TARVIA-lab/ngs-variant-plugin.git
cd ngs-variant-plugin
pip install -r requirements.txt
python verify_setup.py
```

## Bioinformatics Tools (10-15 minutes)

This plugin requires external bioinformatics tools. Install via conda:

```bash
# Recommended: use conda for all tools
conda create -n ngs-variant python=3.10 -y
conda activate ngs-variant

# Install bioinformatics tools
conda install -c bioconda bwa samtools gatk4 bcftools snpEff -y

# Verify installation
python verify_bioinformatics.py
```

**Tools Required:**
- **BWA-MEM2** — Fast DNA alignment
- **GATK4** — Variant calling (germline + somatic)
- **bcftools** — VCF manipulation and stats
- **SnpEff** — Variant annotation

## Claude Code Plugin Setup

```bash
# Register as local marketplace
/plugin add-marketplace /path/to/ngs-variant-plugin/.agents/plugins/marketplace.json

# Install plugin
/plugin install ngs-variant
```

## Verify Everything

```bash
python verify_setup.py          # Python deps
python verify_bioinformatics.py # Bioinformatics tools
```

See [README.md](README.md) for usage examples.
