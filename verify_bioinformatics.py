#!/usr/bin/env python3
"""Verify bioinformatics tools for ngs-variant-plugin"""
import subprocess
import shutil
import sys

GREEN, YELLOW, RED, RESET = '\033[92m', '\033[93m', '\033[91m', '\033[0m'

tools = {
    'bwa': ['bwa', 'Version'],
    'samtools': ['samtools', '--version'],
    'gatk': ['gatk', '--version'],
    'bcftools': ['bcftools', '--version'],
    'snpEff': ['snpEff/snpEff.jar', '-version'],  # Special case
}

print("\nBioinformatics Tools\n" + "-"*40)

ok = True
for name, cmd in tools.items():
    if shutil.which(cmd[0]):
        print(f"{GREEN}✓{RESET} {name}")
    else:
        print(f"{RED}✗{RESET} {name} - conda install -c bioconda {name}")
        ok = False

if ok:
    print(f"\n{GREEN}✓ All tools available!{RESET}\n")
else:
    print(f"\n{YELLOW}Missing tools:{{RESET}}")
    print("conda install -c bioconda bwa samtools gatk4 bcftools snpEff\n")
    sys.exit(1)
