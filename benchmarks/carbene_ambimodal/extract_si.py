"""Pull every structure out of the SI of Wang & Xue, JACS 2026, 148, 8723 (doi 10.1021/jacs.5c13906).

Input: the pdftotext -layout dump of the SI (si/ja5c13906_si.txt). The SI lists each structure as
an xyz block: atom count, then "NAME G=<Hartree>", then coordinates; page numbers can fall inside
a block and are skipped. DFT level of every geometry: M06-2X/def2-TZVP-SMD(PhCl or toluene).

Writes structures/<section>/<name>.xyz (extended xyz: info carries name, section, G_hartree,
charge 0, spin 1 or 3) and structures/index.json.

Usage: python benchmarks/carbene_ambimodal/extract_si.py
"""

import json
import re
from pathlib import Path

from ase import Atoms
from ase.io import write

HERE = Path(__file__).resolve().parent
TEXT = HERE / "si" / "ja5c13906_si.txt"
OUT = HERE / "structures"

ATOM = re.compile(r"^\s*([A-Z][a-z]?)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s*$")
# "NAME G=<Hartree>" for optimized structures; "NAME Nu_i=<cm-1>" for the def2-SVP TSs used in QCMD
TITLE = re.compile(r"^\s*(\S.*?)\s+(G|Nu_i)\s*=\s*(-?\d+\.\d+)\s*$")
SECTIONS = {
    "Carbenes": "carbenes",
    "Olefins": "olefins",
    "Transition states": "ts",
    "Transition states (quasi-classical simulation)": "ts_qcmd",
    "Products": "products",
}


def is_count(line):
    return line.lstrip("\f").isdigit()


def is_page_number(line):
    return line.strip().isdigit() and line[:1] == " "


def slug(name):
    return re.sub(r"[^A-Za-z0-9+-]+", "_", name).strip("_")


def main():
    lines = TEXT.read_text().split("\n")  # not splitlines(): it also breaks at the page-feed \f
    start = next(i for i, line in enumerate(lines) if line.strip() == "Structural information")
    section, index, i = None, [], start
    while i < len(lines):
        stripped = lines[i].strip()
        if stripped in SECTIONS:
            section = SECTIONS[stripped]
        # The atom count starts the line; page numbers are centred (indented), and one can sit
        # between the count and the title (TS-6J) or just before the count (TS-CI2).
        k = i + 2 if i + 2 < len(lines) and is_page_number(lines[i + 1]) else i + 1
        if is_count(lines[i]) and k < len(lines) and TITLE.match(lines[k]):
            n = int(stripped)
            name, kind, value = TITLE.match(lines[k]).groups()
            symbols, positions, j = [], [], k + 1
            while len(symbols) < n:
                m = ATOM.match(lines[j])
                if m:
                    symbols.append(m.group(1))
                    positions.append([float(x) for x in m.group(2, 3, 4)])
                elif lines[j].strip() and not lines[j].strip().isdigit():
                    raise ValueError(f"unexpected line {j + 1} inside {name}: {lines[j]!r}")
                j += 1
            spin = 3 if "triplet" in name.lower() else 1
            atoms = Atoms(symbols, positions=positions)
            label = {"G": "G_hartree", "Nu_i": "imaginary_cm"}[kind]
            atoms.info.update(
                name=name, section=section, charge=0, spin=spin, **{label: float(value)}
            )
            stem, dup = slug(name), 1
            while any(row["file"] == f"structures/{section}/{stem}.xyz" for row in index):
                dup += 1  # same name listed twice (two Dimethoxycarbene entries)
                stem = f"{slug(name)}_{dup}"
            path = OUT / section / f"{stem}.xyz"
            path.parent.mkdir(parents=True, exist_ok=True)
            write(path, atoms, format="extxyz")
            index.append(
                dict(
                    name=name,
                    section=section,
                    file=str(path.relative_to(HERE)),
                    n_atoms=n,
                    formula=atoms.get_chemical_formula(),
                    spin=spin,
                    **{label: float(value)},
                )
            )
            i = j
            continue
        i += 1
    (OUT / "index.json").write_text(json.dumps(index, indent=1))
    for row in index:
        print(f"{row['section']:9s} {row['name']:28s} {row['formula']:12s} spin={row['spin']}")
    print(len(index), "structures")


if __name__ == "__main__":
    main()
