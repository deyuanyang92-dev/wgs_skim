"""说明.md in every function directory of tests/all_features.sh, and README.md above them.

    report.py SIM_DIR OUT_DIR          (after tests/sim/check_truth.py)
Each function directory is wgs_skim's own output (<tool>/<sample>/, <tool>/<tool>_summary.tsv, pipeline_info/);
说明.md says what the function does, the command, its directory layout, every sample's result against the truth
and the verdict."""
import json
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "workflow" / "scripts"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from check_truth import covered  # noqa: E402
from seqtools import fasta_lengths, n50  # noqa: E402
from tables import read_tsv  # noqa: E402

# what each tool directory holds (the 目录 section of 说明.md)
DIRS = {
    "fastp": "clean reads 与质控报告；fastp_summary.tsv = 每个样品的数据量、Q30",
    "meangs": "MEANGS 线粒体种子；meangs_summary.tsv",
    "getorganelle": "GetOrganelle 各轮组装（round0、round1…）；getorganelle_summary.tsv = 长度、环化、轮次",
    "mitogenome": "最终线粒体（GetOrganelle 或 SPAdes 中更好的一个）；mitogenome_summary.tsv",
    "spades": "SPAdes 组装（<样品>.scaffolds.fasta / .rna.fasta / .meta.fasta）；spades_summary.tsv = 条数、总长、N50",
    "trinity": "Trinity 转录本；trinity_summary.tsv",
    "busco": "BUSCO 完整度；busco_summary.tsv = 每个样品（每个组装）的 C / S / D / F / M",
    "genomescope": "k-mer 直方图与 GenomeScope；genomescope_summary.tsv = 基因组大小、杂合度、错误率",
    "mitom": "MitoM 注释（GenBank .gb、基因表）；mitom_summary.tsv = 蛋白基因 / rRNA / tRNA 数",
    "rdna": "rDNA 单元（GetOrganelle）与注释 .gbf；rdna_summary.tsv",
    "histone": "组蛋白簇（GetOrganelle）与注释 .gbf；histone_summary.tsv",
    "nuclear_markers": "SPAdes 组装中的 rDNA / 组蛋白拷贝（MitoM）、拷贝数、与 GetOrganelle 的比较",
    "topology": "reads 验证环状 / 线状",
    "blast": "blast_topn.tsv = 每条序列的最佳命中与物种分类",
    "summary": "samples.tsv（所有工具，每样品一行）、report.html",
    "pipeline_info": "运行记录：日志 logs/、用时 benchmarks/、配置 config_*.yaml、命令与版本 run_*.txt",
}


def dir_line(p):
    n = p.name
    if n.startswith("mitofinder_"):
        return (f"`{n}/` MitoFinder（{n.removeprefix('mitofinder_')} 的组装）：<样品>/；mitofinder_summary.tsv = 基因表"
                "（N = 未找到）；genes/<基因>.fasta；contamination/ = 污染检查")
    return f"`{n}/` {DIRS.get(n, '')}" if p.is_dir() else f"`{n}`"


def busco(f):
    m = re.search(r"C:[\d.]+%\[S:[\d.]+%,D:[\d.]+%\],F:[\d.]+%,M:[\d.]+%,n:\d+", Path(f).read_text()) \
        if Path(f).exists() else None
    return m.group(0) if m else "无结果"


def extra_rows(fn, d, sim, truth):
    """per-sample results of the functions check_truth.py has no truth item for: [(sample, item, result)]"""
    out = []
    for s in truth:
        if fn.endswith(("_fastp", "_fastp_rna")):
            j = json.load(open(d / "fastp" / s / f"{s}.fastp.json"))["summary"]["after_filtering"]
            out.append((s, "clean 数据", f"{j['total_bases'] / 1e9:.2f} Gb，Q30 {j['q30_rate'] * 100:.1f}%"))
        elif fn.endswith("_meangs"):
            f = d / "meangs" / s / f"{s}.seed.fasta"
            cov, pid = covered(sim / "truth" / s / "mt.fasta", f) if f.exists() and f.stat().st_size else (0, 0)
            out.append((s, "种子", f"{sum(fasta_lengths(f)) if f.exists() else 0} bp，覆盖真实线粒体 {cov:.0f}%，"
                                 f"一致 {pid:.2f}%"))
        elif fn.endswith(("spades_meta", "_rnaspades", "_trinity")):
            f = d / "trinity" / s / f"{s}.trinity.fasta" if (d / "trinity").is_dir() else \
                next(iter(sorted((d / "spades" / s).glob(f"{s}.*.fasta"))), None)
            lens = fasta_lengths(f) if f and f.exists() else []
            out.append((s, "组装", f"{len(lens)} 条，共 {sum(lens) / 1e6:.1f} Mb，N50 {n50(lens)}"))
        elif "busco" in fn or fn.endswith("python_api"):
            for f in sorted((d / "busco" / s).glob("*.busco.txt")):
                what = f.name.removeprefix(s + ".").removesuffix("busco.txt").strip(".") or "基因组"
                out.append((s, f"BUSCO {what}", busco(f)))
        elif fn.endswith("circularize"):
            for r in read_tsv(d / s / f"{s}.circularize.tsv"):
                out.append((s, "环化", f"{r['method']}，环长 {r['circle_length'] or '-'} bp（真值 "
                                     f"{truth[s]['mt_length']} bp）"))
        elif fn.endswith("_compare"):
            for r in read_tsv(d / s / f"{s}.compare.tsv"):
                out.append((s, "比较", f"{r['verdict']}，一致 {r['identity_pct']}%，{r['ref_length']} vs "
                                     f"{r['other_length']} bp"))
        elif fn.endswith("_nuclear"):
            v = {}
            for r in read_tsv(d / "summary.tsv"):
                if r["sample"] == s:
                    v[r["verdict"]] = v.get(r["verdict"], 0) + 1
            out.append((s, "基因比较", "，".join(f"{k} {n}" for k, n in v.items())))
        elif fn.endswith("_blast"):
            best = {}
            for r in read_tsv(d / "blast" / "blast_topn.tsv"):
                if r["qseqid"].startswith(s + "|"):
                    best.setdefault(r["qseqid"].split("|")[-1], r)
            for g, r in best.items():
                out.append((s, f"BLAST {g}", f"{r['organism']}（{float(r['pident']):.1f}%）"))
    return out


def main():
    sim, out = Path(sys.argv[1]).resolve(), Path(sys.argv[2]).resolve()
    checks = read_tsv(out / "truth_check.tsv")
    truth = {r["sample"]: r for r in read_tsv(sim / "truth.tsv")}
    runs = {}
    for f in sorted(out.glob("features.*.tsv")):            # the latest run of each function wins
        for r in read_tsv(f):
            runs[f"{r['id']}_{r['function']}"] = r
    index = ["# wgs_skim 全功能测试（5 个模拟样品）", "",
             f"模拟数据与真值：`{sim}`（truth.tsv）。每个功能一个目录，目录里就是 wgs_skim 的输出："
             "`<工具>/<样品>/`、`<工具>/<工具>_summary.tsv`（汇总）、`pipeline_info/`（日志、配置）。"
             "每个目录的 `说明.md`：做什么、命令、目录、每个样品与真值的对比、结论。clean reads 只做一次"
             "（02 DNA、30 RNA），其后的功能直接读取它们。", "",
             "结果中的英文词：consistent = 各基因属同一科，no close reference = 数据库里没有近缘种（不是错误），"
             "possible contamination = 可能污染，partial hit = 只命中基因的一部分，circular / linear = 环状 / 线状，"
             "identical / contained / differ = 相同 / 包含于 / 不同。", "",
             "| 功能 | 运行 | 与真值一致 | 结论 |", "|---|---|---|---|"]
    for d in sorted(p for p in out.iterdir() if p.is_dir() and re.match(r"\d\d_", p.name)):
        fn, run = d.name, runs.get(d.name, {})
        rows = [r for r in checks if r["function"] == fn]
        bad = [r for r in rows if r["verdict"] != "ok"]
        ran = run.get("result") == "ok"
        verdict = ("通过" if not bad else f"{len(bad)} 项与真值不同（见说明）") if ran else "运行失败"
        md = [f"# {fn}", "", f"**做什么 / 看什么**：{run.get('what', '')}", "",
              f"**运行**：{'成功' if ran else '失败 ' + run.get('missing', '')}，{run.get('seconds', '?')} 秒", "",
              "**命令**：", "```bash", run.get("command", ""), "```", "", "## 目录", ""]
        md += [f"- {dir_line(p)}" for p in sorted(d.iterdir()) if p.name != "说明.md"]
        md += ["", "## 每个样品的结果", "", "| 样品 | 项目 | 结果 | 真值 | 一致 |", "|---|---|---|---|---|"]
        md += [f"| {r['sample']} | {r['item']} | {r['result']} | {r['truth']} | "
               f"{'是' if r['verdict'] == 'ok' else '**否**'} |" for r in rows]
        md += [f"| {s} | {item} | {res} | | |" for s, item, res in extra_rows(fn, d, sim, truth)]
        md += ["", f"**结论**：{verdict}", ""]
        if bad:
            md += ["不一致：", *[f"- {r['sample']} {r['item']}：{r['result']}（真值 {r['truth']}）" for r in bad], ""]
            if "mitofinder" in fn and any("MitoFinder" in r["item"] for r in bad):
                md.append("MitoFinder 只用一个远缘参考时，远缘样品 S3–S5 会漏掉 ND2、ND4L、rrnL；同样数据用 5 个近缘参考"
                          "（16）时 15 / 15。这是 MitoFinder 的限制，不是流程错误。")
        (d / "说明.md").write_text("\n".join(md) + "\n")
        tc = f"{len(rows) - len(bad)}/{len(rows)}" if rows else "（无真值项，见说明）"
        index.append(f"| [{fn}]({fn}/说明.md) | {'成功' if ran else '失败'} | {tc} | {verdict} |")
    (out / "README.md").write_text("\n".join(index) + "\n")
    print(f"说明.md in every function directory; index -> {out / 'README.md'}")


if __name__ == "__main__":
    main()
