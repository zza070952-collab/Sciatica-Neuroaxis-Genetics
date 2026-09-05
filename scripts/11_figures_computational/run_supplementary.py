"""Rebuild S1-S9 computational panels from locked tabular outputs."""
from pathlib import Path
import os, subprocess, sys
import supplementary_base as base
import supplementary_updates as updated

if __name__=='__main__':
    os.environ.setdefault('SCIATICA_PROJECT_ROOT',str(Path.cwd()))
    base.style()
    for fn in (base.s1,base.s2,base.s3,base.s4): fn()
    # Original all-section data renderer; final publication crop/assembly is documented.
    subprocess.run([sys.executable,str(Path(__file__).with_name('all_section_atlas.py'))],check=True)
    base.s7_s8_s9(7)
    updated.style()
    for fn in (updated.s6,updated.s8,updated.s9): fn()
