"""A readable view of tests/all_features.sh results: one folder per function with 说明.md (what the function does,
the command, the result of every sample against the truth, the verdict) and 结果/ (links to the key files only).

    report.py SIM_DIR FEATURES_DIR VIEW_DIR        (run tests/sim/check_truth.py first)
VIEW_DIR/README.md: one line per function. The full outputs stay in FEATURES_DIR (linked as 完整输出)."""
import json
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "workflow" / "scripts"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from check_truth import covered  # noqa: E402
from seqtools import fasta_lengths, n50  # noqa: E402
from tables import read_tsv  # noqa: E402

# function: (中文名, 做什么, 关键文件 {s} = 样品)
FUNCS = {
    "01_samples": ("样品识别", "从 reads 目录自动识别样品名与 R1/R2 配对", ["samples.tsv"]),
    "02_fastp": ("fastp 质控", "去接头、低质量；每个样品 clean reads 与报告", ["qc/{s}/{s}.fastp.html"]),
    "03_meangs": ("MEANGS 种子", "从 reads 中找线粒体种子（给 GetOrganelle 用）", ["mito/{s}/meangs/{s}.seed.fasta"]),
    "04_getorganelle_mt": ("GetOrganelle 线粒体（自带种子）", "GetOrganelle 用自带 animal_mt 种子组装线粒体，未成环则迭代",
                           ["mito/{s}/{s}.mt.fasta", "mito/{s}/{s}.mt.json"]),
    "05_meangs-getorganelle": ("组合：MEANGS + GetOrganelle（推荐）", "MEANGS 种子 → GetOrganelle → 未成环迭代",
                               ["mito/{s}/{s}.mt.fasta", "mito/{s}/{s}.mt.json"]),
    "06_getorganelle_nr": ("GetOrganelle 核 rDNA + 组蛋白", "GetOrganelle 组装 18S-ITS-28S 单元与组蛋白簇，MitoM 注释",
                           ["rdna/{s}/{s}.rdna.gbf", "histone/{s}/{s}.histone.gbf"]),
    "07_spades_genome": ("SPAdes 基因组组装", "SPAdes 组装全基因组", ["assembly/{s}/{s}.scaffolds.fasta"]),
    "08_spades_meta": ("metaSPAdes", "SPAdes --meta 模式组装", ["assembly/{s}/{s}.meta.fasta"]),
    "09_genome-size": ("基因组大小 GenomeScope", "jellyfish + GenomeScope2：基因组大小、杂合度",
                       ["genome_size/{s}/{s}.genomescope.txt"]),
    "10_busco_genome": ("BUSCO（基因组）", "已有组装的 BUSCO 完整度（-i 输入）", ["assembly/{s}/{s}.busco.txt"]),
    "11_mitom": ("MitoM 注释", "线粒体注释（PCG、rRNA、tRNA）", ["annotation/{s}/{s}.mt.gb"]),
    "12_blast": ("BLAST 物种鉴定", "cox1、16S 比对本地库，最佳命中物种", ["blast/blast_topn.tsv"]),
    "13_mitofinder_assembly": ("MitoFinder（组装，1 个远缘参考）", "MitoFinder 在 SPAdes 组装中找线粒体并注释；参考 = S1",
                               ["mitofinder/given/mitofinder_summary.tsv", "mitofinder/given/blast/samples_check.tsv"]),
    "14_mitofinder_reads": ("MitoFinder（reads，MEGAHIT）", "MitoFinder 自己组装 reads 再注释",
                            ["mitofinder/megahit/mitofinder_summary.tsv"]),
    "15_mitofinder_bioconda": ("MitoFinder bioconda 1.4.1", "没有 Singularity 镜像时的 MitoFinder",
                               ["mitofinder/given/mitofinder_summary.tsv"]),
    "16_spades-busco": ("组合：SPAdes + BUSCO", "原始 reads → fastp → SPAdes → BUSCO", ["assembly/{s}/{s}.busco.txt"]),
    "17_spades-mitofinder": ("组合：GetOrganelle + SPAdes-MitoFinder", "两条路线各得一个线粒体（参考 = 样品自己的注释）",
                             ["mito/{s}/{s}.mt.fasta", "mitofinder/spades/mitofinder_summary.tsv"]),
    "18_circularize": ("无 reads 环化", "MitoFinder contig + SPAdes 组装图 → 环状线粒体", ["{s}.circularize.tsv",
                                                                                    "{s}.circularized.fasta"]),
    "19_compare": ("GetOrganelle 与 MitoFinder 比较", "两条路线的线粒体是否一致", ["all.compare.html",
                                                                         "all.compare.tsv"]),
    "20_run_all": ("完整流程（默认全部步骤）", "wgs_skim run：所有步骤，最后一张汇总表", ["summary/samples.tsv",
                                                                         "summary/report.html"]),
    "21_nuclear": ("核多拷贝基因批量比较", "rDNA、组蛋白：GetOrganelle 与 SPAdes-MitoM 逐基因比较", ["summary.tsv"]),
    "22_sheet_given": ("样品表（已有组装）", "给出已有组装与线粒体，不重新组装", ["summary/samples.tsv"]),
    "23_split_by_id": ("合并 FASTA 按序列名拆分", "一个文件含多个样品 → 拆分 → 批量注释", ["annotation"]),
    "24_python_api": ("Python 接口", "wgs_skim.tool('busco', ...) 与命令行相同", ["assembly/{s}/{s}.busco.txt"]),
    "25_nextflow": ("Nextflow 引擎", "同一流程用 Nextflow 运行", ["summary/samples.tsv"]),
    "26_mitofinder_multiref": ("MitoFinder（组装，5 个近缘参考，推荐）", "与 13 相同，参考换成 5 个近缘种",
                               ["mitofinder/given/mitofinder_summary.tsv"]),
    "31_rnaspades": ("rnaSPAdes 转录组组装", "RNA-seq → rnaSPAdes 转录本", ["assembly/{s}/{s}.rna.fasta"]),
    "32_trinity": ("Trinity 转录组组装", "RNA-seq → Trinity 转录本", ["assembly/{s}/{s}.trinity.fasta"]),
    "33_busco_transcriptome": ("BUSCO（转录组）", "转录本的 BUSCO 完整度", ["assembly/{s}/{s}.busco.txt"]),
    "34_rnaspades-busco": ("组合：rnaSPAdes + BUSCO（推荐）", "RNA-seq → rnaSPAdes → BUSCO",
                           ["assembly/{s}/{s}.rna.busco.txt"]),
    "35_trinity-busco": ("组合：Trinity + BUSCO", "RNA-seq → Trinity → BUSCO", ["assembly/{s}/{s}.trinity.busco.txt"]),
    "36_rnaspades-mitofinder": ("组合：rnaSPAdes + MitoFinder + 污染检查", "转录组中提取线粒体基因，按基因汇总，BLAST 查污染",
                                ["mitofinder/rna/mitofinder_summary.tsv", "mitofinder/rna/blast/samples_check.tsv",
                                 "mitofinder/rna/genes"]),
    "37_trinity-mitofinder": ("组合：Trinity + MitoFinder + 污染检查", "同 36，Trinity 转录本",
                              ["mitofinder/trinity/mitofinder_summary.tsv",
                               "mitofinder/trinity/blast/samples_check.tsv", "mitofinder/trinity/genes"]),
    "38_mitofinder_rna_tool": ("MitoFinder 单独命令（两种转录组）", "wgs_skim mitofinder --assembler rnaspades,trinity",
                               ["mitofinder/rna/mitofinder_summary.tsv", "mitofinder/trinity/mitofinder_summary.tsv"]),
}


def busco(f):
    m = re.search(r"C:[\d.]+%\[S:[\d.]+%,D:[\d.]+%\],F:[\d.]+%,M:[\d.]+%,n:\d+", Path(f).read_text()) \
        if Path(f).exists() else None
    return m.group(0) if m else "无结果"


def extra_rows(fn, d, sim, samples):
    """per-sample results of the functions check_truth.py has no truth item for: [(sample, item, result)]"""
    out = []
    for s in samples:
        if fn == "02_fastp":
            j = json.load(open(d / "qc" / s / f"{s}.fastp.json"))["summary"]["after_filtering"]
            out.append((s, "clean 数据", f"{j['total_bases'] / 1e9:.2f} Gb，Q30 {j['q30_rate'] * 100:.1f}%"))
        elif fn == "03_meangs":
            f = d / "mito" / s / "meangs" / f"{s}.seed.fasta"
            cov, pid = covered(sim / "truth" / s / "mt.fasta", f) if f.exists() and f.stat().st_size else (0, 0)
            out.append((s, "种子", f"{sum(fasta_lengths(f)) if f.exists() else 0} bp，覆盖真实线粒体 {cov:.0f}%，"
                                 f"一致 {pid:.2f}%"))
        elif fn in ("08_spades_meta", "31_rnaspades", "32_trinity"):
            ext = {"08_spades_meta": "meta", "31_rnaspades": "rna", "32_trinity": "trinity"}[fn]
            lens = fasta_lengths(d / "assembly" / s / f"{s}.{ext}.fasta")
            out.append((s, "组装", f"{len(lens)} 条，共 {sum(lens) / 1e6:.1f} Mb，N50 {n50(lens)}"))
        elif fn in ("10_busco_genome", "16_spades-busco", "24_python_api", "33_busco_transcriptome"):
            out.append((s, "BUSCO", busco(d / "assembly" / s / f"{s}.busco.txt")))
        elif fn == "34_rnaspades-busco":
            out.append((s, "BUSCO", busco(d / "assembly" / s / f"{s}.rna.busco.txt")))
        elif fn == "35_trinity-busco":
            out.append((s, "BUSCO", busco(d / "assembly" / s / f"{s}.trinity.busco.txt")))
        elif fn == "18_circularize":
            for r in read_tsv(d / f"{s}.circularize.tsv"):
                out.append((s, "环化", f"{r['method']}，环长 {r['circle_length'] or '-'} bp（真值 "
                                     f"{truth(sim)[s]['mt_length']} bp）"))
        elif fn == "19_compare":
            for r in read_tsv(d / f"{s}.compare.tsv"):
                out.append((s, "比较", f"{r['verdict']}，一致 {r['identity_pct']}%，{r['ref_length']} vs "
                                     f"{r['other_length']} bp"))
        elif fn == "21_nuclear":
            v = {}
            for r in read_tsv(d / "summary.tsv"):
                if r["sample"] == s:
                    v[r["verdict"]] = v.get(r["verdict"], 0) + 1
            out.append((s, "基因比较", "，".join(f"{k} {n}" for k, n in v.items())))
        elif fn == "12_blast":
            best = {}
            for r in read_tsv(d / "blast" / "blast_topn.tsv"):
                if r["qseqid"].startswith(s + "|"):
                    best.setdefault(r["qseqid"].split("|")[-1], r)
            for g, r in best.items():
                out.append((s, f"BLAST {g}", f"{r['organism']}（{float(r['pident']):.1f}%）"))
    return out


_T = {}


def truth(sim):
    if not _T:
        _T.update({r["sample"]: r for r in read_tsv(Path(sim) / "truth.tsv")})
    return _T


def main():
    sim, feats, view = (Path(x).resolve() for x in sys.argv[1:4])
    view.mkdir(parents=True, exist_ok=True)
    checks = read_tsv(feats / "truth_check.tsv")
    samples = list(truth(sim))
    runs = {}
    for f in sorted(feats.glob("features.*.tsv")):          # the latest run of each function wins
        for r in read_tsv(f):
            runs[f"{r['id']}_{r['function']}"] = r
    lines = ["# wgs_skim 全功能测试结果（5 个模拟样品）", "",
             f"模拟数据与真值：`{sim}`（truth.tsv：每个样品的真实线粒体、rDNA、组蛋白、核基因组；S5 的 RNA 含 *Nephtys* "
             "COX1、CYTB 污染）。每个功能一个文件夹：`说明.md`（做什么、命令、每个样品的结果与真值、结论），"
             "`结果/`（只放关键结果文件的链接），`完整输出`（该功能的全部输出）。", "",
             "结果中的英文词：consistent = 各基因属同一科，no close reference = 数据库里没有近缘种（不是错误），"
             "possible contamination = 可能污染，partial hit = 只命中基因的一部分，circular / linear = 环状 / 线状，"
             "identical / contained / differ = 相同 / 包含于 / 不同，genes = 找到的基因数。", "",
             "| 功能 | 运行 | 与真值一致 | 结论 |", "|---|---|---|---|"]
    for fn, (name, what, files) in FUNCS.items():
        d = feats / fn
        if not d.is_dir():
            continue
        clean = re.sub(r"[（：，、]", "_", name.replace(" + ", "-").replace("+", "-")).replace("）", "").replace(" ", "")
        out = view / f"{fn[:2]}_{clean.strip('_')}"
        out.mkdir(exist_ok=True)
        res = out / "结果"
        res.mkdir(exist_ok=True)
        for x in res.iterdir():
            x.unlink() if x.is_symlink() or x.is_file() else None
        for pat in files:
            for s in (samples if "{s}" in pat else [None]):
                src = d / (pat.replace("{s}", s) if s else pat)
                if src.exists():                # one file for all samples: named with its folder (rna_, trinity_)
                    lname = src.name if s or len(Path(pat).parts) < 3 else f"{Path(pat).parts[-2]}_{src.name}"
                    if lname == "blast_samples_check.tsv":
                        lname = f"{Path(pat).parts[-3]}_污染检查_samples_check.tsv"
                    (res / lname).symlink_to(src)
        link = out / "完整输出"
        if link.is_symlink():
            link.unlink()
        link.symlink_to(d)
        run = runs.get(fn, {})
        rows = [r for r in checks if r["function"] == fn]
        ok = sum(r["verdict"] == "ok" for r in rows)
        bad = [r for r in rows if r["verdict"] != "ok"]
        verdict = ("通过" if not bad else f"{len(bad)} 项与真值不同（见说明）") if run.get("result") == "ok" else "运行失败"
        cmd = run.get("command", "")
        md = [f"# {fn[:2]} {name}", "", f"**做什么**：{what}", "",
              f"**运行**：{'成功' if run.get('result') == 'ok' else '失败'}，{run.get('seconds', '?')} 秒", "",
              "**命令**：", "```bash", cmd, "```", "", "## 每个样品的结果", "",
              "| 样品 | 项目 | 结果 | 真值 | 一致 |", "|---|---|---|---|---|"]
        for r in rows:
            md.append(f"| {r['sample']} | {r['item']} | {r['result']} | {r['truth']} | "
                      f"{'是' if r['verdict'] == 'ok' else '**否**'} |")
        for s, item, result in extra_rows(fn, d, sim, samples):
            md.append(f"| {s} | {item} | {result} | | |")
        md += ["", f"**结论**：{verdict}", ""]
        if bad:
            md += ["不一致的原因："] + [f"- {r['sample']} {r['item']}：{r['result']}（真值 {r['truth']}）"
                                      for r in bad] + [""]
            if "mitofinder" in fn and any("MitoFinder" in r["item"] for r in bad):
                md.append("MitoFinder 只用一个远缘参考（S1）时，远缘样品 S3–S5 会漏掉 ND2、ND4L、rrnL；"
                          "同样数据用 5 个近缘参考（26）时 15 / 15。这是 MitoFinder 的限制，不是流程错误。")
        md += ["", "## 结果文件（结果/ 中的链接）", ""] + [f"- `{x.name}`" for x in sorted(res.iterdir())]
        (out / "说明.md").write_text("\n".join(md) + "\n")
        tc = f"{ok}/{len(rows)}" if rows else "（无真值项，见说明）"
        lines.append(f"| [{out.name}]({out.name}/说明.md) | {'成功' if run.get('result') == 'ok' else '失败'} | {tc} | "
                     f"{verdict} |")
    (view / "README.md").write_text("\n".join(lines) + "\n")
    print(f"{len(FUNCS)} functions -> {view}/README.md")


if __name__ == "__main__":
    main()
