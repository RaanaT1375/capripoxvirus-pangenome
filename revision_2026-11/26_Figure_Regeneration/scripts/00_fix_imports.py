#!/usr/bin/env python3
"""Stage 26, step 0 — make the other scripts find _common.py from anywhere.

The four scripts in this directory assumed _common.py sits exactly two
levels up. That holds only when the stage directory is a direct child of
revision_2026-11; if the zip was unpacked anywhere else, the import fails
with "No module named '_common'".

This replaces that fixed path with a search that walks up the directory tree
until it finds _common.py, so the scripts work wherever they are unpacked
and whatever directory you run them from. It also honours $COMMON_DIR if you
want to point at _common.py directly.

Run it once, from anywhere:

    python3 scripts/00_fix_imports.py

It is safe to re-run: scripts already patched are skipped. It prints where
_common.py was found, or says plainly that it could not find one.
"""
import os
import sys

OLD = 'sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))'

NEW = '''def _locate_common():
    env = os.environ.get("COMMON_DIR")
    if env and os.path.exists(os.path.join(env, "_common.py")):
        return env
    d = os.path.dirname(os.path.abspath(__file__))
    for _ in range(8):
        if os.path.exists(os.path.join(d, "_common.py")):
            return d
        parent = os.path.dirname(d)
        if parent == d:
            break
        d = parent
    sys.exit(
        "could not find _common.py in any parent of\\n  "
        + os.path.abspath(__file__)
        + "\\nIt belongs in the revision_2026-11 directory alongside the "
          "numbered stages.\\nIf it lives somewhere else, point at it:\\n"
          "  export COMMON_DIR=/path/to/the/directory/holding/_common.py")


sys.path.insert(0, _locate_common())'''

TARGETS = ["01_nonviral_screen.py", "02_refit_heaps.py",
           "03_figure1_data.py", "04_scoary_table.py"]

here = os.path.dirname(os.path.abspath(__file__))

# where is _common.py, from here?
found = None
d = here
for _ in range(8):
    if os.path.exists(os.path.join(d, "_common.py")):
        found = d
        break
    p = os.path.dirname(d)
    if p == d:
        break
    d = p

print(f"scripts directory: {here}")
if found:
    print(f"_common.py found:  {os.path.join(found, '_common.py')}")
else:
    print("_common.py NOT found in any parent of this directory.")
    print()
    print("Locate it and set COMMON_DIR before running the other scripts:")
    print("  find \"$PANG\" -maxdepth 3 -name _common.py")
    print("  export COMMON_DIR=/the/directory/it/is/in")
    print()
    print("Patching the scripts anyway, so they will pick up COMMON_DIR.")

print()
changed = skipped = missing = 0
for name in TARGETS:
    path = os.path.join(here, name)
    if not os.path.exists(path):
        print(f"  {name:<26} not present")
        missing += 1
        continue
    src = open(path).read()
    if "_locate_common" in src:
        print(f"  {name:<26} already patched")
        skipped += 1
        continue
    if OLD not in src:
        print(f"  {name:<26} import line not recognised — left alone")
        missing += 1
        continue
    open(path, "w").write(src.replace(OLD, NEW, 1))
    print(f"  {name:<26} patched")
    changed += 1

print(f"\n{changed} patched, {skipped} already done, {missing} skipped")
if found:
    print("\nNow run them in order:")
    print("  python3 scripts/01_nonviral_screen.py")
