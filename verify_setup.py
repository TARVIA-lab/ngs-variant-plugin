#!/usr/bin/env python3
"""Verify Python dependencies for ngs-variant-plugin"""
import sys

GREEN, RED, RESET = '\033[92m', '\033[91m', '\033[0m'

def check_python():
    v = sys.version_info
    ok = v.major >= 3 and v.minor >= 10
    print(f"{'✓' if ok else '✗'} Python {v.major}.{v.minor}")
    return ok

def check_package(name, imp=None):
    try:
        __import__(imp or name)
        print(f"{GREEN}✓{RESET} {name}")
        return True
    except:
        print(f"{RED}✗{RESET} {name}")
        return False

print("\nPython Dependencies\n" + "-"*40)
ok = check_python()
for pkg in ['pyyaml', 'pandas', 'pydantic']:
    ok &= check_package(pkg)

if ok:
    print(f"\n{GREEN}✓ Python dependencies ready!{RESET}")
    print("Run: python verify_bioinformatics.py (for tools setup)\n")
else:
    print(f"\n{RED}✗ Missing dependencies{RESET}")
    print("pip install -r requirements.txt\n")
    sys.exit(1)
