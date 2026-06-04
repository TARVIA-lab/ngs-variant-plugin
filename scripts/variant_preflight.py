#!/usr/bin/env python3
"""
variant_preflight.py — Check tool availability for the variant calling pipeline.

Usage:
    python variant_preflight.py --list
    python variant_preflight.py --pipeline germline --emit-install-plan
    python variant_preflight.py --tool gatk
"""
import argparse, importlib, json, shutil, subprocess, sys
from datetime import datetime
from pathlib import Path

TOOLS = {
    "fastqc":   {"install": "conda install -c bioconda fastqc",         "pipelines": ["fastq_qc"]},
    "multiqc":  {"install": "pip install multiqc",                      "pipelines": ["fastq_qc"]},
    "bwa-mem2": {"install": "conda install -c bioconda bwa-mem2",       "pipelines": ["alignment"]},
    "samtools": {"install": "conda install -c bioconda samtools",       "pipelines": ["alignment", "germline", "somatic"]},
    "gatk":     {"install": "conda install -c bioconda gatk4",          "pipelines": ["germline", "somatic"]},
    "bcftools": {"install": "conda install -c bioconda bcftools",       "pipelines": ["annotation"]},
    "snpEff":   {"install": "conda install -c bioconda snpeff",         "pipelines": ["annotation"]},
}

PIPELINE_TOOLS = {
    "fastq_qc":  ["fastqc", "multiqc"],
    "alignment": ["bwa-mem2", "samtools", "gatk"],
    "germline":  ["bwa-mem2", "samtools", "gatk"],
    "somatic":   ["bwa-mem2", "samtools", "gatk"],
    "annotation":["bcftools", "snpEff"],
    "all":       list(TOOLS.keys()),
}


def check_tool(name: str) -> dict:
    path = shutil.which(name)
    version = None
    if path:
        try:
            r = subprocess.run([name, "--version"], capture_output=True, text=True, timeout=10)
            out = (r.stdout + r.stderr).strip().splitlines()
            version = out[0] if out else "unknown"
        except Exception:
            version = "unknown"
    return {"name": name, "found": bool(path), "path": path,
            "version": version, "install": TOOLS[name]["install"]}


def emit_install_plan(missing: list, outfile: str = "install_plan.json"):
    conda = [t for t in missing if "conda" in TOOLS[t]["install"]]
    pip   = [t for t in missing if "pip"   in TOOLS[t]["install"]]
    plan  = {
        "generated": datetime.now().isoformat(),
        "shell_commands": (
            (["conda install -c bioconda -c conda-forge " +
              " ".join(TOOLS[t]["install"].split()[-1] for t in conda)] if conda else []) +
            (["pip install " + " ".join(TOOLS[t]["install"].split()[-1] for t in pip)] if pip else [])
        ),
    }
    Path(outfile).write_text(json.dumps(plan, indent=2))
    print(f"\nInstall plan → {outfile}")
    for cmd in plan["shell_commands"]:
        print(f"  {cmd}")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--list", action="store_true")
    parser.add_argument("--pipeline", choices=list(PIPELINE_TOOLS.keys()))
    parser.add_argument("--tool")
    parser.add_argument("--emit-install-plan", action="store_true")
    parser.add_argument("--output", default="install_plan.json")
    args = parser.parse_args()

    if args.tool:
        r = check_tool(args.tool)
        s = "✓" if r["found"] else "✗"
        print(f"{s} {r['name']}: {'found — ' + r['version'] if r['found'] else 'NOT FOUND → ' + r['install']}")
        return

    pipeline = args.pipeline or "all"
    tool_names = PIPELINE_TOOLS[pipeline]
    missing = []

    print(f"\n=== Pipeline: {pipeline} ===")
    for t in tool_names:
        r = check_tool(t)
        s = "✓" if r["found"] else "✗"
        print(f"  {s} {r['name']:12s} {'v' + r['version'] if r['found'] else 'NOT FOUND  →  ' + r['install']}")
        if not r["found"]:
            missing.append(t)

    if missing:
        print(f"\n  {len(missing)} tool(s) missing.")
        if args.emit_install_plan:
            emit_install_plan(missing, args.output)
    else:
        print("\n  All requirements satisfied.")


if __name__ == "__main__":
    main()
