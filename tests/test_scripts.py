"""Unit tests of the pure-Python parts of the workflow scripts (no external tools, synthetic data, seconds).

Run: python -m pytest -q tests/        (python with Biopython for the nuclear_markers comparison tests)
The end-to-end check is the small test sample (README, "Testing").
"""
import gzip
import random
import sys
try:
    import defusedxml.minidom as minidom
except ImportError:                       # the SVG is our own output, not untrusted input
    import xml.dom.minidom as minidom
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "workflow" / "scripts"))
import assembly_check  # noqa: E402
import circularize  # noqa: E402
import collect_fasta  # noqa: E402
import compare_mt  # noqa: E402
import nuclear_markers as nm  # noqa: E402
import plot_compare  # noqa: E402
import samples  # noqa: E402
import seqtools  # noqa: E402
import summary  # noqa: E402
import topology  # noqa: E402

rnd = random.Random(3)


def dna(n):
    return "".join(rnd.choice("ACGT") for _ in range(n))


# ── samples.py ─────────────────────────────────────────────────────────────
def test_find_pairs_recognises_provider_names(tmp_path):
    for n in ("A_R1_001.fastq.gz", "A_R2_001.fastq.gz", "B_1.fq.gz", "B_2.fq.gz", "C.R1.raw.fastq.gz",
              "C.R2.raw.fastq.gz", "D_S3_L001_R1_001.fastq.gz", "D_S3_L001_R2_001.fastq.gz", "E_R1_trimmed.fq.gz",
              "E_R2_trimmed.fq.gz", "lonely_R1.fq.gz", "notes.txt"):
        (tmp_path / n).touch()
    (tmp_path / "F").mkdir()
    for n in ("F_1.clean.fq.gz", "F_2.clean.fq.gz"):
        (tmp_path / "F" / n).touch()
    pairs = samples.find_pairs(tmp_path)
    assert list(pairs) == ["A", "B", "C", "D", "E", "F"]
    assert pairs["C"][1].endswith("C.R2.raw.fastq.gz") and pairs["F"][0].endswith("F/F_1.clean.fq.gz")


def test_find_pairs_duplicate_sample_raises(tmp_path):
    for n in ("A_R1.fq.gz", "A_R2.fq.gz", "A_1.fq.gz", "A_2.fq.gz"):
        (tmp_path / n).touch()
    with pytest.raises(ValueError, match="twice"):
        samples.find_pairs(tmp_path)


def test_sheet_relative_paths_and_optional_columns(tmp_path):
    (tmp_path / "x_1.fq").touch()
    (tmp_path / "x_2.fq").touch()
    sheet = tmp_path / "s.tsv"
    sheet.write_text("sample\tfastq_1\tfastq_2\tscaffolds\nx\tx_1.fq\tx_2.fq\t\ny\tx_1.fq\tx_2.fq\tasm/y.fa\n")
    pairs = samples.read_sheet(sheet)
    assert pairs["x"] == (str(tmp_path / "x_1.fq"), str(tmp_path / "x_2.fq"))
    given = samples.sheet_columns(sheet, ["scaffolds"])
    assert given == {"y": {"scaffolds": str(tmp_path / "asm" / "y.fa")}}      # empty cells are not entries


# ── seqtools / collect_fasta ───────────────────────────────────────────────
def test_fasta_roundtrip_and_circle_detection(tmp_path):
    circle = dna(3000)
    linear = circle + circle[:400]                   # a circle assembled as a line: 400 bp end overlap
    f = tmp_path / "a.fa"
    with open(f, "w") as fh:
        seqtools.write_fasta([("c1 circular", linear.lower())], fh)
    (name, seq), = seqtools.read_fasta(f)
    assert name == "c1 circular" and seq == linear
    c = seqtools.find_circle(seq)
    assert c and c["P"] == 3000 and c["overlap"] == 400
    assert seqtools.find_circle(dna(3000)) is None


def test_collect_fasta_prefers_complete_graph(tmp_path):
    (tmp_path / "x.K105.scaffolds.graph1.1.path_sequence.fasta").write_text(">s\nACGT\n")
    (tmp_path / "x.K105.complete.graph1.1.path_sequence.fasta").write_text(">c (circular)\nACGTACGT\n")
    recs, circ = collect_fasta.best(tmp_path)
    assert circ and recs[0][1] == "ACGTACGT"
    assert collect_fasta.renamed(recs, "S1", "mt")[0][0] == "S1|mt|1 circular length=8"
    assert collect_fasta.best(tmp_path / "none") == ([], False)


# ── summary.py ─────────────────────────────────────────────────────────────
def test_summary_helpers(tmp_path):
    assert summary.n50([1, 2, 3, 4, 10]) == 10 and summary.n50([5, 5, 5, 5]) == 5 and summary.n50([]) == ""
    f = tmp_path / "a.fa"
    f.write_text(">a\nACGT\nAC\n>b\n\n>c\nA\n")
    assert summary.fasta_lengths(f) == [6, 1]
    b = tmp_path / "busco.txt"
    b.write_text("\tC:93.2%[S:90.1%,D:3.1%],F:2.0%,M:4.8%,n:954\t\n")
    assert summary.busco_line(b) == "C:93.2%[S:90.1%,D:3.1%],F:2.0%,M:4.8%,n:954"
    g = tmp_path / "gs.txt"
    g.write_text("property   min   max\nHeterozygous (ab)   1.23%   1.31%\n"
                 "Genome Haploid Length   1,000,000 bp   1,234,567 bp\n")
    assert summary.genomescope(g) == ("1.2", "1.31")
    assert summary.by_sample(["r/S1/S1.x", "r/S2/S2.x"]) == {"S1": Path("r/S1/S1.x"), "S2": Path("r/S2/S2.x")}


# ── topology.py ────────────────────────────────────────────────────────────
def test_split_point_follows_rotation_and_strand():
    s = dna(5000)
    nad5 = (1000, 1600)
    rot = s[700:] + s[:700]                          # the annotated record starts 700 bp later
    ann = [(rot, [("nad5", nad5[0] - 700, nad5[1] - 700)])]
    g, m = topology.split_point(s, ann)
    assert g == "nad5" and m == (nad5[0] + nad5[1]) // 2
    rev = rot.translate(topology.COMP)[::-1]         # the annotated record on the other strand: mirrored coordinates
    L = len(s)
    g, m2 = topology.split_point(s, [(rev, [("nad5", L - (nad5[1] - 700), L - (nad5[0] - 700))])])
    assert g == "nad5" and abs(m2 - m) <= 1
    assert topology.split_point(s, [(s[:100], [("nad5", 1, 50)])]) is None   # different length: not this sequence


def test_end_diagnosis_classes():
    assert "clean end" in topology.end_diagnosis("R", 40, [], "ACGTACGT", 40)
    tails = ["ATATATATATATATATAT" + dna(5) for _ in range(20)]
    d = topology.end_diagnosis("L", 40, tails, "ACGT", 40)
    assert "20/40 reads continue" in d and "dinucleotide repeat" in d and "consistently" in d
    assert "homopolymer" in topology.end_diagnosis("R", 10, [], "ACGTTTTTTTTT", 40)
    assert "coverage drops" in topology.end_diagnosis("R", 10, [], "ACGT", 40)


# ── circularize.py ─────────────────────────────────────────────────────────
def test_circularize_overlap_and_join():
    assert circularize.rc("AACGT") == "ACGTT"
    assert circularize.overlap_len("XXXXACGTAC" + "G" * 30, "G" * 30 + "TTTT") == 30
    circle = dna(4000)
    rows, out = circularize.circularize([("c1", circle + circle[:300])])
    assert rows[0]["method"] == "overlap" and len(out[0][1]) == 4000
    # three pieces of one circle, one of them reverse-complemented as MitoFinder writes it
    a, b, c = circle[:1500], circle[1400:3000], circle[2900:] + circle[:100]
    rows, out = circularize.circularize([("a", a), ("b (reverse)", circularize.rc(b)), ("c", c)])
    assert rows[-1]["method"] == "join" and len(out[0][1]) == 4000
    rows, _ = circularize.circularize([("lin", dna(3000))])
    assert rows[0]["method"] == "not_circular"


# ── compare_mt.py ──────────────────────────────────────────────────────────
def test_diff_blocks_and_description():
    ref = dna(2000)
    other = ref[:500] + "T" + ref[500:1200] + ref[1203:]          # 1-bp insertion, 3-bp deletion
    aln_a = ref[:500] + "-" + ref[500:]                           # alignment columns by hand
    aln_b = other[:501] + other[501:1201] + "---" + other[1201:]
    blocks = compare_mt.diff_blocks(aln_a, aln_b)
    assert [(s, e, len(x), len(y)) for s, e, x, y in blocks] == [(500, 500, 0, 1), (1200, 1203, 3, 0)]
    assert compare_mt.kind("", "T") == "insertion" and compare_mt.kind("ACG", "") == "deletion"
    assert compare_mt.kind("A", "C") == "mismatch" and compare_mt.kind("ACGT", "TGCA") == "complex"
    feats = [(400, 900, "cox1", 1, "CDS"), (1500, 1800, "nad2", 1, "CDS")]
    rows, explained = compare_mt.describe(ref, True, [(1190, 1210, "repeat segment 7")], feats, "mf", blocks)
    assert [r["gene"] for r in rows] == ["cox1", "intergenic cox1..nad2"]
    assert rows[1]["graph_region"] == "repeat segment 7" and explained == 1
    assert compare_mt.verdict(rows, explained) == "differ" and compare_mt.verdict([], 0) == "identical"


def test_orient_turns_and_rotates():
    ref = dna(3000)
    seq, strand, rot = compare_mt.orient(ref, circularize.rc(ref[1000:] + ref[:1000]), circular=True)
    assert seq == ref and strand == "-" and rot == 2000


def test_homopolymer_difference_is_explained():
    ref = dna(1000) + "T" * 12 + dna(1000)
    blocks = [(1005, 1005, "", "T")]                             # one T more inside the run
    rows, explained = compare_mt.describe(ref, True, [], [], "x", blocks)
    assert explained == 1 and "homopolymer Tx12" in rows[0]["graph_region"]


# ── assembly_check.py / plot_compare.py ────────────────────────────────────
def test_gff_cleaning_as_modify_gff2(tmp_path):
    g = tmp_path / "a.gff"
    g.write_text("##gff-version 3\nA\tGeneious\tregion\t1\t15817\t.\t+\t0\tIs_circular=true\n"
                 "A\tGeneious\tgene\t1\t1534\t.\t+\t.\tName=COX1 gene\n"
                 "A\tGeneious\tCDS\t1\t1534\t.\t+\t0\tName=COX1 gene\n"
                 "A\tGeneious\ttRNA\t1535\t1598\t.\t-\t.\tName=trnN%28guu%29\n")
    assert assembly_check.read_gff(g) == [(0, 1534, "COX1", 1, "CDS"), (1534, 1598, "trnN", -1, "tRNA")]
    assert assembly_check.gene_type("rrnL", "gene") == "rRNA" and assembly_check.gene_type("ND1", "gene") == "CDS"


def test_svg_figure_is_well_formed():
    rlen = 5000
    contigs = [("c1", "A" * 3000, False), ("c2", "C" * 1800, False)]
    blocks = [dict(contig="c1", qlen=3000, qs=0, qe=3000, strand="+", rs=0, re=3000, ops=[]),
              dict(contig="c2", qlen=1800, qs=0, qe=1800, strand="-", rs=3100, re=4900, ops=[])]
    genes = [(100, 1000, "cox1", 1, "CDS"), (3200, 3300, "trnA", -1, "tRNA")]
    events = [dict(contig="", event="not in assembly", ref_pos="3001-3100", contig_pos="", size=-100,
                   graph_region="")]
    diffs = [dict(ref_start="500", ref_end="500", type="mismatch", ref_bp="1", other_bp="1", gene="cox1",
                  graph_region="")]
    svg = plot_compare.draw(rlen, contigs, blocks, genes, [(4000, 4200, "repeat")], diffs, events, "mf", "t")
    minidom.parseString(svg)
    assert "cox1" in svg and "inverted" in svg


# ── nuclear_markers.py ─────────────────────────────────────────────────────
def test_reads_cap_and_completeness(tmp_path):
    fq = tmp_path / "r1.fq.gz"
    with gzip.open(fq, "wt") as fh:
        for i in range(100):
            fh.write(f"@r{i}\n{'A' * 150}\n+\n{'I' * 150}\n")
    assert nm.reads_per_file(fq, 3.0) == int(3e9 / 2 / 150)
    f = tmp_path / "x.path_sequence.fasta"
    f.write_text(">1(circular)\nACGT\n")
    assert nm.is_circular(f) and nm.is_complete(f) and not nm.is_complete(None)
    f.write_text(">scaffold_1\nACGT\n")
    assert not nm.is_circular(f)
    assert nm.is_complete(f, "rdna", cov={"18S": 100, "5.8S": 99, "28S": 96})      # linear, all genes whole
    assert not nm.is_complete(f, "rdna", cov={"18S": 100, "5.8S": 99, "28S": 80})
    assert nm.is_complete(f, "histone", cov={"H2A": 100, "H2B": 100, "H3": 100, "H4": 100})   # H1 not needed
    assert nm.label_gene("18SsmallsubunitribosomalRNA gene - X", "rdna") == "18S"
    assert nm.label_gene("5.8SribosomalRNA gene - X", "rdna") == "5.8S"
    assert nm.label_gene("H2A|CHY006|cov3377|NODE_1", "histone") == "H2A"
    assert nm.label_gene("MW622063 Hesionidae sp. A histone H3 (H3) gene, partial cds", "histone") == "H3"
    assert nm.label_gene("X53330.1 P.dumerilii histone gene cluster for core histones H2A, H2B, H3 and H4",
                         "histone") is None


def test_compare_and_status_from_consensus_dirs(tmp_path):
    pytest.importorskip("Bio")
    gene = dna(400)

    def side(d, seqs, classes=None):
        (d / "consensus").mkdir(parents=True)
        rows = ["gene\tlength\tpartial_5\tpartial_3\tclass"]
        for g, s in seqs.items():
            (d / "consensus" / f"{g}.fasta").write_text(f">{g}\n{s}\n")
            rows.append(f"{g}\t{len(s)}\t0\t0\t{(classes or {}).get(g, '')}")
        (d / "consensus.tsv").write_text("\n".join(rows) + "\n")
    go, sp = tmp_path / "histone" / "annotation" / "S1", tmp_path / "histone" / "spades_mitom" / "S1"
    side(go, {"H3": gene, "H4": gene[:300], "H2A": gene}, {"H3": "replication-dependent"})
    other = "T" if gene[200] != "T" else "A"
    side(sp, {"H3": gene, "H4": gene, "H2A": gene[:200] + other + gene[201:], "H1": gene})
    rows = {r["gene"]: r for r in nm.compare("S1", "histone", go, sp)}
    assert rows["H3"]["verdict"] == "identical" and rows["H4"]["verdict"] == "contained"
    assert rows["H2A"]["verdict"] == "differ" and rows["H2A"]["mismatches"] == 1
    assert rows["H1"]["verdict"] == "only SPAdes-MitoM" and rows["H3"]["getorganelle_class"] == "replication-dependent"
    for d in (go, sp):
        (d / "mitom.version").write_text("abc\n")
    nm.write_rows(tmp_path / "summary.tsv", list(rows.values()))
    (line,) = nm.status(tmp_path, ("histone",))
    assert line.startswith("S1") and "histone compared +spades" in line and "differ" in line
