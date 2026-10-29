"""Find paired-end FASTQ files and sample names in a reads directory.

A file is read 1 when its name contains a read-number token such as `_R1_001`, `_R1`, `.R1`, `-R1`, `_1` or `.1`
right before the extension or before a suffix like `.raw` / `.clean` / `_trimmed`; read 2 is the same name with 2.
The sample name is the part before the token. Files directly in the directory and one level of sub-directories
(one folder per sample) are both scanned.
"""
import csv
import re
from pathlib import Path

FQ_EXT = r"\.(?:fastq|fq)(?:\.gz)?"
# token, then an optional tag kept by many providers (Novogene .raw/.clean, Macrogen _trimmed, Illumina _001)
R1_RE = re.compile(rf"^(?P<sample>.+?)(?P<sep>[._-])(?P<r>R?)1(?P<tail>_001)?(?P<tag>[._](?:raw|clean|trimmed|sanfastq))?(?P<ext>{FQ_EXT})$")
LANE_RE = re.compile(r"_S\d+_L\d{3}$")   # Illumina bcl2fastq: SAMPLE_S1_L001_R1_001.fastq.gz


def mate(m):
    return f"{m['sample']}{m['sep']}{m['r']}2{m['tail'] or ''}{m['tag'] or ''}{m['ext']}"


def find_pairs(reads_dir):
    """Return {sample: (r1, r2)} sorted by sample. Raises ValueError on duplicate names or unpaired R1 files."""
    d = Path(reads_dir).expanduser().resolve()
    if not d.is_dir():
        raise ValueError(f"reads directory not found: {d}")
    files = [p for p in sorted(d.iterdir()) if p.is_file()] + \
            [p for p in sorted(d.glob("*/*")) if p.is_file()]
    pairs, unpaired = {}, []
    for f in files:
        m = R1_RE.match(f.name)
        if not m:
            continue
        r2 = f.with_name(mate(m))
        if not r2.exists():
            unpaired.append(f.name)
            continue
        sample = LANE_RE.sub("", m["sample"])
        if sample in pairs:
            raise ValueError(f"sample {sample!r} found twice: {pairs[sample][0]} and {f} "
                             "(merge lanes first, or use a samples.tsv)")
        pairs[sample] = (str(f), str(r2))
    if not pairs:
        raise ValueError(f"no paired FASTQ files (*_R1/_R2, *_1/_2, *.R1/.R2 ...) found in {d}")
    if unpaired:
        print(f"[wgs_skim] WARNING: R1 files without an R2 mate were skipped: {', '.join(unpaired)}")
    return dict(sorted(pairs.items()))


def read_sheet(path):
    """samples.tsv / .csv with columns sample, fastq_1, fastq_2 (nf-core style); relative paths are relative to the
    sheet. Further optional columns: see sheet_columns (scaffolds)."""
    p = Path(path).expanduser().resolve()
    delim = "," if p.suffix == ".csv" else "\t"
    pairs = {}
    with open(p, newline="", encoding="utf-8-sig") as fh:
        for row in csv.DictReader(fh, delimiter=delim):
            s = row["sample"].strip()
            r1, r2 = (Path(row[k].strip()).expanduser() for k in ("fastq_1", "fastq_2"))
            r1, r2 = (x if x.is_absolute() else p.parent / x for x in (r1, r2))
            if s in pairs:
                raise ValueError(f"sample {s!r} listed twice in {p}")
            pairs[s] = (str(r1.resolve()), str(r2.resolve()))
    return pairs


def sheet_columns(path, columns):
    """Optional per-sample files from extra sheet columns, e.g. `scaffolds` (an existing SPAdes assembly, used instead
    of running SPAdes): {sample: {column: absolute path}}, only for cells that are filled in."""
    p = Path(path).expanduser().resolve()
    delim = "," if p.suffix == ".csv" else "\t"
    out = {}
    with open(p, newline="", encoding="utf-8-sig") as fh:
        for row in csv.DictReader(fh, delimiter=delim):
            for c in columns:
                v = (row.get(c) or "").strip()
                if v:
                    f = Path(v).expanduser()
                    f = f if f.is_absolute() else p.parent / f
                    out.setdefault(row["sample"].strip(), {})[c] = str(f.resolve())
    return out


if __name__ == "__main__":
    import sys
    for s, (r1, r2) in find_pairs(sys.argv[1]).items():
        print(s, r1, r2, sep="\t")
