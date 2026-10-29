"""Tab-separated tables, shared by the workflow scripts: every result table is written and read here.

    write_tsv("summary.tsv", rows, ["sample", "length"])   # rows: dicts; missing keys -> "", extra keys ignored
    rows = read_tsv("summary.tsv")                          # [dict]
"""
import csv
from pathlib import Path


def write_tsv(path, rows, cols, append=False):
    """rows (dicts) as a TSV with the columns `cols` (header written unless appending to a non-empty file)"""
    path = Path(path)
    header = not (append and path.exists() and path.stat().st_size)
    with open(path, "a" if append else "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(cols), delimiter="\t", extrasaction="ignore", restval="")
        if header:
            w.writeheader()
        w.writerows(rows)
    return path


def read_tsv(path):
    """[dict] of a TSV with a header ([] when the file is missing or empty)"""
    path = Path(path)
    if not (path.exists() and path.stat().st_size):
        return []
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))
