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
import mitofinder_check  # noqa: E402
import mitofinder_summary  # noqa: E402
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
    assert summary.genomescope(g) == ("1.2", "1.31", None)
    # a genome skim: GenomeScope converges on the organelle peaks (e2e2 cwsc010: 0.28 Mb, 11.4 % read errors)
    assert summary.genome_size_check("0.3", 11.4, "0.30").startswith("unreliable: read error 11.4%")
    assert "coverage 3x < 15x" in summary.genome_size_check("1200", 0.4, "4.0")
    assert summary.genome_size_check("1200", 0.4, "40") == ""
    assert summary.by_sample(["r/S1/S1.x", "r/S2/S2.x"]) == {"S1": Path("r/S1/S1.x"), "S2": Path("r/S2/S2.x")}


# ── topology.py ────────────────────────────────────────────────────────────
def test_split_point_follows_rotation_and_strand():
    s = dna(5000)
    nad5 = (1000, 1600)
    rot = s[700:] + s[:700]                          # the annotated record starts 700 bp later
    ann = [(rot, [("nad5", nad5[0] - 700, nad5[1] - 700)])]
    g, m = topology.split_point(s, ann)
    assert g == "nad5" and m == (nad5[0] + nad5[1]) // 2
    rev = seqtools.rc(rot)         # the annotated record on the other strand: mirrored coordinates
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


def test_smk_shim_runs_snakemake_script(tmp_path):
    """nextflow/bin/smk_shim.py: a `script:` script (nuc_copies.py) reads the JSON as its `snakemake` object."""
    import json
    import subprocess
    busco = tmp_path / "busco" / "run_x"
    busco.mkdir(parents=True)
    (busco / "full_table.tsv").write_text("# head\n" + "".join(
        f"b{i}\tComplete\tNODE_{i}_length_900_cov_{c}\n" for i, c in enumerate((9.0, 10.0, 11.0))))
    (tmp_path / "f.tsv").write_text("gene\tcontig\n18S\tNODE_7_length_8000_cov_300.0\n")
    (tmp_path / "c.tsv").write_text("gene\tlength\tpartial_5\tpartial_3\n18S\t1800\t0\t0\n")
    spec = {"input": {"features": "f.tsv", "consensus": "c.tsv", "busco": "b.txt"}, "output": {"tsv": "out.tsv"},
            "params": {"busco_dir": "busco"}, "log": ["log"], "wildcards": {"sample": "S"},
            "scriptdir": str(HERE.parent / "workflow" / "scripts")}
    (tmp_path / "smk.json").write_text(json.dumps(spec))
    subprocess.run([sys.executable, str(HERE.parent / "nextflow" / "bin" / "smk_shim.py"), "smk.json",
                    str(HERE.parent / "workflow" / "scripts" / "nuc_copies.py")], cwd=tmp_path, check=True)
    rows = (tmp_path / "out.tsv").read_text().splitlines()
    assert rows[1].split("\t") == ["18S", "1800", "", "1", "300.0", "10.00", "30"]


def test_fasta_samples_dir_and_merged_split(tmp_path):
    """-i: one sample per FASTA file; a merged FASTA split by the first header word, repeated names all kept"""
    d = tmp_path / "in"
    d.mkdir()
    (d / "A.fasta").write_text(">c1\nACGT\n")
    (d / "B.fa.gz").write_bytes(gzip.compress(b">c1\nGG\n"))
    assert sorted(samples.fasta_samples([d])) == ["A", "B"]
    merged = tmp_path / "all.fasta"
    merged.write_text(">B2 topology=circular\nAAA\n>B17\nCC\n>B2 topology=linear\nTT\n")
    got = samples.fasta_samples([merged], tmp_path / "split", split_by_id=True)
    assert sorted(got) == ["B17", "B2"]
    assert Path(got["B2"]).read_text() == ">B2_contig1 topology=circular\nAAA\n>B2_contig2 topology=linear\nTT\n"
    assert Path(got["B17"]).read_text() == ">B17\nCC\n"
    man = (tmp_path / "split" / "split_manifest.tsv").read_text().splitlines()
    assert "B2\tB2_contig2\tB2 topology=linear" in man


def test_read_sheet_without_reads(tmp_path):
    """a sheet row with only an assembly (FASTA tools) has empty reads"""
    sheet = tmp_path / "s.tsv"
    sheet.write_text("sample\tfastq_1\tfastq_2\tscaffolds\nA\t\t\tA.fasta\n")
    assert samples.read_sheet(sheet) == {"A": ("", "")}
    assert samples.sheet_columns(sheet, ["scaffolds"])["A"]["scaffolds"] == str(tmp_path / "A.fasta")


def test_mitofinder_summary_gene_table(tmp_path):
    """gene lengths per sample (longest copy) and one FASTA per gene across samples"""
    for s, cox1 in (("A", "ATGAAA"), ("B", "ATG")):
        res = tmp_path / s / "Final_Results"
        res.mkdir(parents=True)
        (tmp_path / s / f"{s}.mitofinder.fasta").write_text(">c\n" + "A" * 100 + "\n")
        extra = f">{s}@ND6\nAAATTT\n>{s}@ND6\nAAA\n" if s == "B" else f">{s}@ND5\nCC\n"
        (res / f"{s}_final_genes_NT.fasta").write_text(f">{s}@COX1\n{cox1}\n>{s}@rrnL\nGGGG\n" + extra)
        (res / f"{s}.infos").write_text("Circularization: Yes\n" if s == "A" else "Circularization: Not found\n")
    out = tmp_path / "sum.tsv"
    sys.argv = ["x", "-o", str(out), str(tmp_path / "A"), str(tmp_path / "B")]
    mitofinder_summary.main()
    rows = [r.split("\t") for r in out.read_text().splitlines()]
    assert rows[0] == ["sample", "contigs", "length", "circular", "protein_genes", "rRNAs", "tRNAs", "COX1", "ND5",
                       "ND6", "rrnL"]
    assert rows[1] == ["A", "1", "100", "yes", "2", "1", "0", "6", "2", "N", "4"]
    assert rows[2] == ["B", "1", "100", "no", "2", "1", "0", "3", "N", "6", "4"]
    assert (tmp_path / "genes" / "COX1.fasta").read_text() == ">A\nATGAAA\n>B\nATG\n"
    assert (tmp_path / "genes" / "ND6.fasta").read_text() == ">B\nAAATTT\n>B_copy2\nAAA\n"   # unique names
    assert (tmp_path / "contigs" / "A_all_mt.fasta").exists()


def test_mitofinder_check_contamination_vs_no_close_reference():
    """another family is contamination only when close; distant off-family hits = no close relative in nt
    (MJ20260820041 CHY006: COX1 82 % Opheliidae, COX2 77 % Drosophila over 69 %, ATP6 83 % over 15 %)"""
    lin = {"Oph": "Metazoa; Annelida; Opheliidae; Travisia", "Dro": "Metazoa; Arthropoda; Drosophilidae; Drosophila",
           "Nep": "Metazoa; Annelida; Nephtyidae; Nephtys"}

    def hit(q, fam, pid, cov):
        return {"qseqid": q, "sacc": f"X{fam}", "organism": fam, "taxonomy": lin[fam], "pident": str(pid),
                "qcovs": str(cov)}
    rows = [hit("A|COX1", "Oph", 82.0, 99), hit("A|CYTB", "Oph", 78.1, 92), hit("A|COX2", "Dro", 77.5, 69),
            hit("A|ATP6", "Dro", 83.5, 15),
            hit("M|COX1", "Nep", 93.6, 100), hit("M|CYTB", "Oph", 78.1, 92), hit("M|ND5", "Oph", 78.5, 90)]
    qids = ["A|COX1", "A|CYTB", "A|COX2", "A|ATP6", "A|ND6", "M|COX1", "M|CYTB", "M|ND5"]
    genes, samples = mitofinder_check.check(rows, qids, 85)
    s = {x["sample"]: x for x in samples}
    assert s["A"]["family"] == "Opheliidae" and s["A"]["verdict"].startswith("consistent family, no close reference")
    assert s["M"]["verdict"] == "possible contamination: COX1"
    v = {g["query"]: g["verdict"] for g in genes}
    assert v["A|COX2"].startswith("no close reference (best Drosophilidae")
    assert v["A|ATP6"].startswith("partial hit") and v["A|ND6"] == "no hit"
    assert mitofinder_check.family("Metazoa; Annelida; Nephtyidae; Nephtys") == "Nephtyidae"


# ── shared modules: seqtools, tables, collect_fasta, runtools, mitofinder_run / prefilter ─────────────────
def test_seqtools_shared_helpers(tmp_path):
    assert seqtools.rc("ACGTNacgtRY") == "RYacgtNACGT"
    assert seqtools.parse_fasta(">a x\nac\ngt\n\n>b\nAA\n") == [("a x", "ACGT"), ("b", "AA")]
    seqtools.write_fasta([("a", "A" * 100), ("e", "")], tmp_path / "w.fa", width=0)
    assert (tmp_path / "w.fa").read_text() == ">a\n" + "A" * 100 + "\n>e\n"
    seqtools.write_fasta([("a", "A" * 100)], tmp_path / "w80.fa")
    assert [len(x) for x in (tmp_path / "w80.fa").read_text().split()] == [2, 80, 20]
    assert seqtools.is_circular("x|mt|1 circular length=5") and not seqtools.is_circular("x linear")
    gb = tmp_path / "r.gb"
    gb.write_text("LOCUS       R1   12 bp DNA circular\nFEATURES\nORIGIN\n        1 acgtacgtac gt\n//\n"
                  "LOCUS       R2   4 bp\nORIGIN\n        1 aaaa\n//\n")
    assert seqtools.genbank_seqs(gb) == [("R1", "ACGTACGTACGT"), ("R2", "AAAA")]


def test_tables_roundtrip(tmp_path):
    import tables
    f = tables.write_tsv(tmp_path / "t.tsv", [{"a": 1, "b": "x", "extra": 9}, {"a": 2}], ["a", "b"])
    assert f.read_text() == "a\tb\n1\tx\n2\t\n"
    tables.write_tsv(f, [{"a": 3, "b": "y"}], ["a", "b"], append=True)
    assert tables.read_tsv(f) == [{"a": "1", "b": "x"}, {"a": "2", "b": ""}, {"a": "3", "b": "y"}]
    assert tables.read_tsv(tmp_path / "missing.tsv") == []


def test_getorganelle_result_and_command(tmp_path):
    d = tmp_path / "go"
    d.mkdir()
    for n in ("animal_mt.K115.scaffolds.graph1.1.path_sequence.fasta", "animal_mt.K115.complete.graph1.1.path_sequence.fasta",
              "animal_mt.K115.complete.graph1.2.path_sequence.fasta"):
        (d / n).write_text(">p (circular)\nACGT\n")
    (d / "animal_mt.K115.complete.graph1.selected_graph.gfa").write_text("H\n")
    fasta, gfa, n = collect_fasta.result(d)
    assert fasta.name.startswith("animal_mt.K115.complete.graph1.1") and gfa.name.endswith("selected_graph.gfa")
    assert n == 2 and collect_fasta.result(tmp_path / "none") == (None, None, 0)
    assert collect_fasta.best(d) == ([("p (circular)", "ACGT")], True)
    cmd = collect_fasta.getorganelle_cmd("r1", "r2", "out", "anonym", 4, seed="s.fa", rounds=10, max_reads=None,
                                         genes="l.fa", extra=["-P", "0"])
    assert cmd[:10] == ["get_organelle_from_reads.py", "-1", "r1", "-2", "r2", "-o", "out", "-F", "anonym", "-t"]
    assert cmd[cmd.index("--max-reads") + 1] == "inf" and cmd[-2:] == ["-P", "0"]
    assert "-k" not in cmd and cmd[cmd.index("-s") + 1] == "s.fa" and cmd[cmd.index("--genes") + 1] == "l.fa"


def test_runtools_aligned_bp_and_run(tmp_path):
    import runtools
    hits = [{"sseqid": "a", "pident": "99.0", "length": "100"}, {"sseqid": "a", "pident": "90", "length": "50"},
            {"sseqid": "b", "pident": "100", "length": "10"}]
    assert runtools.aligned_bp(hits, min_pident=95) == {"a": 100, "b": 10}
    assert runtools.run(["sh", "-c", "echo hi; exit 3"], tmp_path / "x.log") == 3
    assert (tmp_path / "x.log").read_text() == "hi\n"
    with pytest.raises(RuntimeError):
        runtools.run(["false"], tmp_path / "y.log", check=True)


def test_mitofinder_collect_no_contig(tmp_path):
    """MitoFinder writes Final_Results even when it finds no mitochondrial contig (tools2: cwsc010) -> empty
    outputs and 0, not a failed job"""
    import mitofinder_run
    final = tmp_path / "S_MitoFinder_mitfi_Final_Results"
    final.mkdir()
    assert mitofinder_run.collect(final, tmp_path, "S") == 0
    assert (tmp_path / "S.mitofinder.fasta").read_text() == "" and (tmp_path / "S.mitofinder.gb").exists()
    (final / "S_mtDNA_contig.fasta").write_text(">c\nACGT\n")
    (final / "S_mtDNA_contig_genes_NT.fasta").write_text(">g\nAC\n")
    (final / "S_mtDNA_contig.gb").write_text("LOCUS c\n//\n")
    assert mitofinder_run.collect(final, tmp_path, "S") == 1
    assert (tmp_path / "S.mitofinder.fasta").read_text() == ">c\nACGT\n"
    assert mitofinder_run.mitofinder_exe("", []) == ["mitofinder"]
    exe = mitofinder_run.mitofinder_exe("/x/mf.sif", [tmp_path, tmp_path])
    assert exe[:3] == ["singularity", "exec", "--cleanenv"] and exe.count("-B") == 1 and exe[-1] == "mitofinder"


def test_split_by_id_sample_before_bar(tmp_path):
    f = tmp_path / "merged.fasta"
    f.write_text(">A|mt|1 circular\nAC\n>A|mt|2 linear\nGT\n>B2 COI\nTT\n")
    out = samples.split_fasta([f], tmp_path / "split")
    assert sorted(out) == ["A", "B2"]
    assert Path(out["A"]).read_text() == ">A_contig1 circular\nAC\n>A_contig2 linear\nGT\n"   # topology kept


def test_mitofinder_prefilter_coverage():
    import mitofinder_prefilter as mp
    assert mp.weighted_median([(100, 1.0), (100, 2.0), (1000, 50.0)]) == 50.0
    assert mp.coverage("NODE_3_length_1500_cov_12.5") == 12.5 and mp.coverage("contig1") is None


def test_every_script_imports_without_side_effects(tmp_path, monkeypatch):
    """docs/DEVELOPING.md rule 1: importing a module runs nothing (no files written, no exit)"""
    import importlib
    monkeypatch.chdir(tmp_path)
    scripts = HERE.parent / "workflow" / "scripts"
    for p in sorted(scripts.glob("*.py")):
        if p.stem == "build_marker_db" or p.stem.startswith("_"):
            continue
        importlib.import_module(p.stem)
    assert list(tmp_path.iterdir()) == []


def test_foreign_gene_same_species_only():
    """sim5: S5 COX1 100 % to the Nephtys mitogenome EU293739 that no other gene matches -> foreign; S4 rrnS
    88.7 % to Micronephthys (a related genus, other genes 76 %) -> not foreign"""
    def hit(q, acc, fam, pid, cov=99, slen=17000):
        return {"qseqid": q, "sacc": acc, "organism": fam, "taxonomy": f"Metazoa; Annelida; {fam}; x",
                "pident": str(pid), "qcovs": str(cov), "slen": str(slen)}
    rows = [hit("S5|COX1", "EU293739", "Nephtyidae", 100), hit("S5|rrnS", "Pis1", "Pisionidae", 82.4),
            hit("S5|COX3", "Are1", "Arenicolidae", 75.9),
            hit("S4|rrnS", "NC_087810", "Nephtyidae", 88.7), hit("S4|ND6", "NC_087810", "Nephtyidae", 76.3),
            hit("S4|COX1", "Eul1", "Eulepethidae", 77.2)]
    qids = ["S5|COX1", "S5|rrnS", "S5|COX3", "S4|rrnS", "S4|ND6", "S4|COX1"]
    genes, samples = mitofinder_check.check(rows, qids, 85)
    s = {x["sample"]: x for x in samples}
    assert s["S5"]["verdict"] == "possible contamination: COX1" and s["S5"]["family"] != "Nephtyidae"
    assert not s["S4"]["verdict"].startswith("possible contamination")


def test_extra_copies_from_joined_transcripts(tmp_path):
    """MitoFinder annotates COX1 on t1 only: t2 (also COX1 by blastx, reverse strand) becomes an extra copy"""
    import mitofinder_run
    out = tmp_path
    (out / "Final_Results").mkdir()
    t1, t2 = "ATG" + "A" * 97, "C" * 40 + "TTTGGG" + "C" * 54
    joined = t1 + "N" * 200 + t2
    (tmp_path / "joined.fasta").write_text(f">j\n{joined}\n")
    (out / "Final_Results" / "S_mtDNA_contig.gff").write_text(
        "j\tmitofinder\tCDS\t1\t100\t.\t+\t0\tName=COX1 CDS\n")
    places = [{"transcript": "t1", "start": "1", "end": "100", "gene": "cox1", "gene_start": "1", "gene_end": "100"},
              {"transcript": "t2", "start": "301", "end": "400", "gene": "COI", "gene_start": "46",
               "gene_end": "41"}]
    assert mitofinder_run.back_to_transcripts(out, "S", tmp_path / "joined.fasta", places) == 1
    extra = seqtools.read_fasta(out / "S.mitofinder.extra_copies.fasta")
    assert extra == [("COX1 t2 46-41", "CCCAAA")]
    assert seqtools.read_fasta(out / "S.mitofinder.fasta") == [("t1", t1)]


def test_hits_without_genbank_metadata_still_checked():
    """sim5 36 (NCBI efetch down): hits without organism / taxonomy are reported as such (not 'no hit'), the
    placeholder is no organism match, and a near-identical foreign copy is still flagged by its accession"""
    ph = "(no GenBank metadata)"

    def hit(q, acc, org, tax, pid, slen=17000):
        return {"qseqid": q, "sacc": acc, "organism": org, "taxonomy": tax, "pident": str(pid), "qcovs": "99",
                "slen": str(slen)}
    rows = [hit("S|COX1", "KY1", ph, "", 77.2), hit("S_copy2|COX1", "EU293739", ph, "", 100),
            hit("S|rrnS", "KY2", "Pisione sp.", "Metazoa; Pisionidae; Pisione", 82.4)]
    genes, samples = mitofinder_check.check(rows, ["S|COX1", "S_copy2|COX1", "S|rrnS"], 85)
    v = {g["query"]: g["verdict"] for g in genes}
    assert v["S|COX1"].startswith("hit KY1 77.2%, no GenBank metadata")
    assert v["S_copy2|COX1"].startswith("possible contamination: EU293739")
    assert samples[0]["verdict"].startswith("possible contamination: COX1") and "incomplete: 2 gene" in samples[0]["verdict"]


def test_genbank_meta_uses_the_cache(tmp_path):
    import blast_taxonomy
    cache = tmp_path / "meta.tsv"
    cache.write_text("accession\torganism\ttaxonomy\tdefinition\tseq_length\ttopology\treference\n"
                     "EU293739\tNephtys sp.\tMetazoa; Nephtyidae\tNephtys mt\t\tcircular\t\n")
    meta = blast_taxonomy.genbank_meta(["EU293739.1"], cache=cache)   # in the cache: no network
    assert meta["EU293739"]["organism"] == "Nephtys sp."


def test_two_copies_of_a_gene_are_two_organisms():
    """sim5 S5 (final): own COX1 / CYTB copies 75-80 % to scattered families, contaminant copies 100 % Nephtys -
    the two contaminant genes match the same mitogenome but must not vouch for each other"""
    def hit(q, acc, org, fam, pid, slen=17000):
        return {"qseqid": q, "sacc": acc, "organism": org, "taxonomy": f"Metazoa; Annelida; {fam}; x",
                "pident": str(pid), "qcovs": "99", "slen": str(slen)}
    rows = [hit("S|COX1", "KY1", "Eulepethus", "Eulepethidae", 77.2),
            hit("S_copy2|COX1", "EU293739", "Nephtys sp.", "Nephtyidae", 100),
            hit("S|CYTB", "PV1", "Terebellides", "Trichobranchidae", 74.9),
            hit("S_copy2|CYTB", "KY972462", "Nephtys sp. AB", "Nephtyidae", 100, slen=400),
            hit("S_copy2|CYTB", "EU293739", "Nephtys sp.", "Nephtyidae", 99.6),
            hit("S|rrnS", "KY2", "Pisione sp.", "Pisionidae", 82.4)]
    qids = ["S|COX1", "S_copy2|COX1", "S|CYTB", "S_copy2|CYTB", "S|rrnS"]
    _, samples = mitofinder_check.check(rows, qids, 85)
    assert samples[0]["verdict"] == "possible contamination: COX1, CYTB"
    assert samples[0]["family"] != "Nephtyidae"
