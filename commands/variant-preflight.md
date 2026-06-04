---
description: Check tool availability for the variant calling pipeline
argument-hint: [--pipeline alignment|germline|somatic|annotation|all] [--emit-install-plan]
allowed-tools: [Bash]
---

# Variant Calling Preflight Check

## Arguments

$ARGUMENTS

## Instructions

Run:
```bash
python scripts/variant_preflight.py $ARGUMENTS
```

If tools are missing and `--emit-install-plan` is passed, write `install_plan.json` and show the install commands. Offer to run them.
