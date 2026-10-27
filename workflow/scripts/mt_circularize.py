#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
mt_circularize.py
判断线粒体组装结果是环形还是线性：能闭合的切成环形序列，并用 reads 验证首尾接合处；
不能闭合的，再判断是"真正的线性末端"还是"组装没拼完"。适用于 MitoFinder / SPAdes / MEGAHIT 的 contig，
也适用于 GetOrganelle / NOVOPlasty 已标为 circular 的结果（组装软件的 circular 标记不一定可靠）。

为什么需要：
  1. MitoFinder 1.4 的判断（circularizationCheck.py）是 contig 自身 BLAST：命中须从前 220 bp 内开始、
     落在末端 220 bp 内、长度 ≥ 40–50 bp。从头组装的线粒体 contig 两端常停在控制区的串联重复里，最外侧
     几十 bp 来自不同的重复拷贝（"毛边"），真正一致的重叠只有 30 多 bp，于是报 "Circularization: Not found"。
     例：MJ20260820041 CHY006_0729 的 SPAdes contig NODE_1248（16,828 bp）两端 58 bp 重叠中前 24 bp 不一致、
     后 34 bp 完全一致。
  2. 组装软件的 "circular" 也会出错：Ahuja et al. 2024（GBE 16:evae048，dunnlab/siph_skimming）发现 NOVOPlasty
     把管水母 Chelophyes appendiculata 的线粒体判为环形，用 reads 检查首尾接合处后确认是线性；GenBank 上几个
     标为 circular 的钵水母线粒体也没有 reads 跨过首尾接合处。

方法：
  A. 找首尾重叠（改进 MITObim circules.py，作者 2024 ZSE Perinereis 论文所用）：开头区段的 k-mer 在末端区段再次
     出现 → 候选环长 P（按支持的 k-mer 数投票，只统计"开头 ↔ 末端"）；在重叠区找 s[i]==s[i+P] 的最长一致片段
     （容忍两端毛边），≥ --min_run 才接受；环 = s[a:a+P]，a 为一致片段中点。找不到重叠时，按原样首尾相接（适用于
     GetOrganelle 等已去掉重叠的环形输出）。
  B. 首尾接合处检验（Ahuja et al. 2024 的做法）：把环从中间切开、首尾对接，使原来的接合处位于参考中部；
     reads 端到端比上并跨过接合处（两侧各 ≥ --flank bp、软剪切 < 10 bp、NM ≤ 2%）→ 环形；没有 → 不是环。
     与接合处两侧 150/300/600 bp 的局部对照点比较（控制区常 AT 富集、覆盖度低，用全环对照会低估）：
     ≥ 对照中位数 × 0.5 判为通过，0.2–0.5 为 weak，> 2 倍为 repeat_junction（接合处在重复区，长度不确定）。
     对照点 reads 中位数 < 30 时自动加大 reads（×4，最多 4000 万对）。
     跨越按 read 级比对范围计（NM ≤ 5%），不按比对块：接合处落在长同聚物里时 indel 会切断比对块。
     （实测：CHY006 真环 0.91–1.02；CHY002 接合处在 20 个 T 的 poly-T 中 0.56；控制区重复塌缩的 SPAdes 环 0.02；
       人为截断的 12 kb 片段 0；rDNA 片段 0）
  C. 逐位点跨越数（dunnlab count_spanning_reads.py 的思路）：每个相邻碱基对 i→i+1 有多少条 reads 的比对块同时
     覆盖两者；内部出现 0 的位置 = 没有 reads 支持的拼接点（可能拼错），列在 zero_span。
  D. 线性末端 vs 没拼完：reads 比到原始线性序列，看两端最外 30 bp 的 reads 有多少在外侧被软剪切（≥10 bp）：
     多数被剪切 = 序列还在往外延伸（组装不完整，或首尾实际相连但重叠没找到）；
     几乎不剪切、reads 整齐地终止在末端 = 真实的分子末端（线性染色体，如水母线粒体）。
  E. 可选 --compare：与另一组装（如 GetOrganelle）BLAST 比较长度与一致度，提示塌缩的串联重复等差异。

reads 用前 --validate_pairs 对（默认 1000 万）：先整体比一次，只取比上的 read 对，再在这批小数据上做 B–D。

结论（verdict）：
  circular_validated     首尾接合处有 reads 支持（≥ 对照点 50%）
  circular_repeat_junction 首尾相接，但接合处跨越 reads 是对照的 2 倍以上：接合处在（可能塌缩的）重复区，环长不确定
  circular_weak          首尾能接上但支持偏少（20–50%），需人工检查
  circular_rejected      找到首尾重叠，但接合处没有 reads 支持（如重复区塌缩，实际分子更长）
  linear_gap_circular    首尾不直接相连，但有 read 对一端在末尾朝外、另一端在开头朝外 → 分子是环，组装缺一段（给出缺口估计）
  linear_unjoined        首尾既无跨越 reads 也无跨首尾 read 对（Ahuja et al. 2024 的"线性"标准）。短 reads 只能到此：
                         可能是线性分子，也可能是缺口长于插入片段的环（如长重复控制区未拼出）。end_state 给出两端情况
                         （clean/extends/mixed，末端 reads 数为对照的倍数），仅供参考——已证实线性的管水母 mt 末端也常 extends。
  （--legacy_linear_classes：旧的 linear_terminus / linear_incomplete / linear_undetermined，按末端软剪切划分）
  *_unvalidated          没给 reads；*_low_reads 对照点跨越 reads < 5，无法判断

输出（-o DIR）：circularize.tsv；<prefix>.circular.fasta（verdict 为 circular_validated / circular_unvalidated）；
  <prefix>.other.fasta（其余，原样）；<prefix>.spanning.tsv.gz（逐位点跨越数，原始坐标）

用法：
  python3 mt_circularize.py -i mt_contigs.fasta -o circ/ --prefix CHY006_0729 \\
      -1 R1.clean.fastq.gz -2 R2.clean.fastq.gz -t 8 --compare getorganelle_mt.fasta
  # 批量：汇总 fasta（header 为 <样本>|locus=...，如 02.organelle/summary_fasta_all.fasta），reads 按样本自动找
  python3 mt_circularize.py --summary summary_fasta_all.fasta --fastp_dir 01.fastp -o circ/ --max_parallel 4 -t 8
  -> circ/<样本>/... 与 circ/topology_summary.tsv
"""
import argparse
import csv
import gzip
import logging
import os
import re
import shutil
import statistics
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path


def read_fasta(path):
    recs = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                recs.append([line[1:], []])
            elif recs:
                recs[-1][1].append(line)
    return [(n, "".join(q).upper()) for n, q in recs]


def write_fasta(path, recs):
    with open(path, "w") as fh:
        for n, q in recs:
            fh.write(f">{n}\n" + "\n".join(q[i:i + 80] for i in range(0, len(q), 80)) + "\n")


# ── 注释（Ahuja et al. 2024 的顺序：先注释，再按注释切开、检验接合处与基因顺序）──────────────

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(s):
    return s.translate(COMP)[::-1]


def norm_gene(n):
    """ND5/nad5 → nad5，CYTB → cob；tRNA 保留大小写（trnL1、trnS2）。"""
    n = n.strip().replace("-", "")
    if n.lower().startswith("trn"):
        return "trn" + n[3:]
    return re.sub(r"^nd(\d|4l)", r"nad\1", n.lower()).replace("cytb", "cob")


def read_genbank(path):
    """不依赖 Biopython 的最小 GenBank 解析：[(seq, [(type, gene, [(start0, end0_exclusive), ...]), ...])]。"""
    recs, feats, seq, cur, in_seq = [], [], [], None, False
    for line in open(path):
        if line.startswith("LOCUS"):
            feats, seq, cur, in_seq = [], [], None, False
        elif line.startswith("ORIGIN"):
            in_seq = True
        elif line.startswith("//"):
            recs.append(("".join(seq).upper(), feats))
            in_seq = False
        elif in_seq:
            seq.append(re.sub(r"[^A-Za-z]", "", line))
        elif re.match(r"^ {5}(CDS|rRNA|tRNA) +", line):
            t, loc = line.split()[:2]
            iv = [(int(a) - 1, int(b)) for a, b in re.findall(r"<?(\d+)\.\.>?(\d+)", loc)]
            cur = [t, "", iv]
            feats.append(cur)
        elif cur is not None and line.startswith(" " * 21) and "/gene=" in line and not cur[1]:
            cur[1] = norm_gene(line.split('"')[1])
        elif not line.startswith(" " * 21):
            cur = None
    return [(q, [tuple(f) for f in fs if f[1] and f[2]]) for q, fs in recs]


def load_annotations(paths):
    recs = []
    for p in paths:
        p = Path(p)
        if not p.exists():
            logging.warning(f"[CIRC] 注释不存在：{p}，按序列中点切开")
            continue
        for f in (sorted(p.glob("*.gb")) + sorted(p.glob("*.gbk")) if p.is_dir() else [p]):
            recs += read_genbank(f)
    return recs


def map_annotation(s, ann):
    """在注释记录里找与 s 相同的序列（注释软件常把序列旋转到 cox1 起点或反向互补），把基因坐标换回 s 的坐标。
    返回 [(type, gene, start0, end0_exclusive, crosses_ends)] 或 None；crosses_ends = 基因跨过 s 的首尾（环状注释才会出现）。"""
    L = len(s)
    for q, feats in ann:
        if len(q) != L:
            continue
        for strand, t in ((1, s), (-1, revcomp(s))):
            o = (t + t).find(q[:60])
            while 0 <= o < L and (t[o:] + t[:o]) != q:
                o = (t + t).find(q[:60], o + 1)
            if not 0 <= o < L:
                continue
            out = []
            for typ, g, iv in feats:
                # 注释坐标的首、末碱基（join 取外包）。b < a = 跨注释序列原点的 join(16001..16500,1..300)：
                # 不交换，保持从 a 顺时针到 b 的弧（交换会得到它的补集）
                a, b = iv[0][0], iv[-1][1] - 1
                pa, pb = (a + o) % L, (b + o) % L
                if strand < 0:
                    pa, pb = L - 1 - pb, L - 1 - pa
                cross = pa > pb
                out.append((typ, g, pa, pb + 1, cross))
            return out
    return None


def gene_split_point(genes, L):
    """dunnlab 做法：在最长的蛋白基因 nad5 中间切开（不跨首尾）；没有 nad5 用最长 CDS；都没有返回 None。"""
    cds = [g for g in genes if g[0] == "CDS" and not g[4]]
    pick = [g for g in cds if g[1] == "nad5"] or sorted(cds, key=lambda g: g[3] - g[2], reverse=True)
    if not pick:
        return None
    g = pick[0]
    return g[1], (g[2] + g[3]) // 2


def end_context(genes, L):
    """两端最近的基因及距离；跨过首尾的基因（序列在基因里断开）。"""
    inner = [g for g in genes if not g[4]]
    left = min(inner, key=lambda g: g[2], default=None)
    right = max(inner, key=lambda g: g[3], default=None)
    cross = [g[1] for g in genes if g[4]]
    txt = (f"L:{left[1]}+{left[2]}bp" if left else "L:-") + " | " + (f"R:{right[1]}+{L - right[3]}bp" if right else "R:-")
    return txt, cross


def adjacency_support(genes, prof):
    """相邻基因之间的拼接点有没有 reads 支持（dunnlab count_spanning_reads 用来确认重排的基因顺序）：
    每对相邻基因（按坐标）从前一基因末 50 bp 到后一基因首 50 bp 的最小跨越数。"""
    inner = sorted((g for g in genes if not g[4]), key=lambda g: g[2])
    res = []
    for g1, g2 in zip(inner, inner[1:]):
        lo, hi = max(0, g1[3] - 50), min(len(prof) - 1, g2[2] + 50)
        if hi > lo:
            res.append((f"{g1[1]}-{g2[1]}", min(prof[lo:hi])))
    return res


# ── A. 首尾重叠 → 环长 ──────────────────────────────────────────────────

def period_candidates(s, k, span, end_window):
    """开头 span bp 内的 k-mer 在末端 end_window bp 内再次出现 → 周期 P = 末端位置 - 开头位置。返回 {P: 票数}。"""
    L = len(s)
    if L < 2 * k:
        return {}
    tail = defaultdict(list)
    for j in range(max(span, L - end_window), L - k + 1):
        tail[s[j:j + k]].append(j)
    votes = defaultdict(int)
    for i in range(0, min(span, L - k + 1)):
        for j in tail.get(s[i:i + k], ()):
            votes[j - i] += 1
    return dict(votes)


def best_run(s, P):
    """重叠区（i < L-P）里 s[i]==s[i+P] 的最长连续片段 (start, length)。"""
    best, cur_start, cur = (0, 0), 0, 0
    for i in range(len(s) - P):
        if s[i] == s[i + P]:
            if cur == 0:
                cur_start = i
            cur += 1
            if cur > best[1]:
                best = (cur_start, cur)
        else:
            cur = 0
    return best


def find_circle(s, k=21, span=2000, end_window=5000, min_run=25, min_len=1000):
    """返回 dict(P, overlap, run_len, votes, anchor) 或 None。"""
    if len(s) < min_len:
        return None
    best = None
    for P, v in sorted(period_candidates(s, k, span, end_window).items(), key=lambda x: -x[1])[:10]:
        if P < min_len:
            continue
        rs, rl = best_run(s, P)
        if rl < min_run:
            continue
        cand = dict(P=P, overlap=len(s) - P, run_len=rl, votes=v, anchor=rs + rl // 2)
        if best is None or (rl, v) > (best["run_len"], best["votes"]):
            best = cand
    return best


# ── 比对与解析 ──────────────────────────────────────────────────────────

def sh(cmd):
    # pipefail：管道中间的 bwa / samtools 出错也要报错（否则只看最后一环，失败会变成"没有 reads"）
    subprocess.run(f"set -o pipefail; {cmd}", shell=True, check=True, executable="/bin/bash")


CIG = re.compile(r"(\d+)([MIDNSHP=X])")


def parse_sam(path):
    """逐条产出 (ref, start0, cigar_ops, nm, flag)；只含主比对。"""
    with open(path) as fh:
        for line in fh:
            if line.startswith("@"):
                continue
            f = line.split("\t", 12)
            if f[5] == "*":
                continue
            nm = next((int(t[5:]) for t in f[11:] if t.startswith("NM:i:")), 0)
            yield f[2], int(f[3]) - 1, [(int(n), op) for n, op in CIG.findall(f[5])], nm, int(f[1])


def bridging_pairs(sam, L, window=600):
    """线性序列上"跨首尾"的 read 对：一端在末尾 window bp 内、正向（朝右、朝外），另一端在开头 window bp 内、反向
    （朝左、朝外）。真正的线性分子不可能有这样的 read 对；有 → 首尾在分子上相邻（环形，中间可能还缺一段）。
    缺口估计：跨首尾 read 对的"已组装跨度"（末端 read 起点→序列末尾 + 序列开头→mate 比对终点，mate 跨度取 MC 标签）
    与序列内部若干点上用同样窗口选出的 read 对的片段长度比较，两者中位数之差即缺口（同样的窗口选择偏差相互抵消）。
    返回 (对数, 缺口中位估计 bp 或 None)。"""
    pts = [int(L * f) for f in (0.25, 0.5, 0.75)] if L > 4 * window else []
    bridge, ctrl = [], []
    with open(sam) as fh:
        for line in fh:
            if line.startswith("@"):
                continue
            f = line.split("\t", 11)
            flag = int(f[1])
            if not flag & 1 or flag & 0xC or f[6] != "=" or flag & 0x10 or not flag & 0x20:
                continue                                    # 只看正向 read、mate 反向的配对
            pos, mpos = int(f[3]) - 1, int(f[7]) - 1
            mc = next((t[5:].strip() for t in f[11].split("\t") if t.startswith("MC:Z:")), None) if len(f) > 11 else None
            mspan = sum(int(x) for x, op in CIG.findall(mc) if op in "MDN=X") if mc else 0
            if not mspan:
                continue
            if pos >= L - window and mpos < window:
                bridge.append((L - pos) + (mpos + mspan))
            for p in pts:                                   # 内部对照：同样的窗口跨过 p
                if p - window <= pos < p <= mpos < p + window:
                    ctrl.append(mpos + mspan - pos)
    return len(bridge), estimate_gap(bridge, ctrl)


def estimate_gap(bridge, ctrl):
    """缺口 g 的估计：内部对照的片段长度分布 F 减去 g、只保留仍能被观察到的（跨度 ≥ 观察到的最小跨度），
    取预测中位数与跨首尾 read 对实际跨度中位数最接近的 g。直接相减会低估：缺口越大，只有越长的片段能跨过。"""
    if len(bridge) < 3 or len(ctrl) < 20:
        return None
    obs, lo = statistics.median(bridge), min(bridge)
    best = None
    for g in range(-300, 3001, 5):
        pred = [x - g for x in ctrl if x - g >= lo]
        if len(pred) < 5:
            break
        d = abs(statistics.median(pred) - obs)
        if best is None or d < best[0]:
            best = (d, g)
    return best[1] if best else None


def aligned_blocks(start, ops):
    """比对块（M/=/X 连续覆盖参考的区段），D/N 断开块（与 count_spanning_reads.py 一致）。"""
    blocks, pos = [], start
    for n, op in ops:
        if op in "M=X":
            blocks.append((pos, pos + n))
        if op in "MDN=X":
            pos += n
    return blocks


def clips(ops):
    left = ops[0][0] if ops and ops[0][1] in "SH" else 0
    right = ops[-1][0] if ops and ops[-1][1] in "SH" else 0
    return left, right


def good_read(ops, nm):
    span = sum(n for n, op in ops if op in "MDN=X")
    return max(clips(ops)) < 10 and nm <= 0.02 * max(span, 1)


def spanning_profile(sam, lengths):
    """{ref: [跨越 i→i+1 的 reads 数]}，只用 good_read。差分数组累加。"""
    diff = {r: [0] * (L + 1) for r, L in lengths.items()}
    for ref, start, ops, nm, _ in parse_sam(sam):
        if ref not in diff or not good_read(ops, nm):
            continue
        for b0, b1 in aligned_blocks(start, ops):
            if b1 - b0 >= 2:
                diff[ref][b0] += 1
                diff[ref][b1 - 1] -= 1
    prof = {}
    for ref, d in diff.items():
        acc, out = 0, []
        for i in range(lengths[ref] - 1):
            acc += d[i]
            out.append(acc)
        prof[ref] = out
    return prof


def span_counts(sam, ref, points, flank, max_nm=0.05):
    """每个点 p：比对范围（read 级，不按比对块）覆盖 [p-flank, p+flank] 且 NM ≤ max_nm×比对长度的 reads 数。
    接合处检验用 read 级跨越：接合处若落在长同聚物（如线粒体控制区的 poly-T）里，测序 indel 会把比对块切断，
    按块计数会低估（CHY002：块计数 73 vs 对照 195.5）；比对块级的 0 值检查仍由 spanning_profile 做。"""
    cnt = {p: 0 for p in points}
    for r, start, ops, nm, _ in parse_sam(sam):
        if r != ref:
            continue
        end = start + sum(n for n, op in ops if op in "MDN=X")
        if nm > max_nm * max(end - start, 1):
            continue
        for p in points:
            if start <= p - flank and end >= p + flank:
                cnt[p] += 1
    return cnt


def end_clip_stats(sam, lengths, window=30):
    """原始线性序列两端：覆盖最外 window bp 的 reads 中，外侧软剪切 ≥10 bp 的比例。{ref: (左比例, 左n, 右比例, 右n)}"""
    st = {r: [0, 0, 0, 0] for r in lengths}
    for ref, start, ops, nm, _ in parse_sam(sam):
        if ref not in st:
            continue
        end = start + sum(n for n, op in ops if op in "MDN=X")
        lc, rc = clips(ops)
        if start < window:
            st[ref][1] += 1
            st[ref][0] += lc >= 10
        if end > lengths[ref] - window:
            st[ref][3] += 1
            st[ref][2] += rc >= 10
    return {r: (v[0] / v[1] if v[1] else None, v[1], v[2] / v[3] if v[3] else None, v[3]) for r, v in st.items()}


def zero_runs(prof, lo, hi):
    """prof[lo:hi] 中为 0 的连续区段 [(a,b)]（0-based，含 a 不含 b）。"""
    runs, a = [], None
    for i in range(lo, hi):
        if prof[i] == 0 and a is None:
            a = i
        elif prof[i] != 0 and a is not None:
            runs.append((a, i))
            a = None
    if a is not None:
        runs.append((a, hi))
    return runs


def analyse_reads(items, a, tools, wd, rot_fa, pairs, prof_rot, prof_lin, clipst):
    """前 pairs 对 reads：先整体比一次挑出比上的 read 对，再逐条序列（各自建索引）比到接合处居中的环和原始线性序列。
    不能把多条序列放进同一个索引：相同/相近的序列（如同一样本的不同组装）会让 reads 随机分给其中一条，
    跨越数被摊薄，线性片段末端的 reads 也会改比到别的完整序列上，看不到软剪切。"""
    head = f"| head -n {pairs * 4}" if pairs > 0 else ""
    sh(f"{tools['bwa']} mem -t {a.threads} {rot_fa} <(zcat -f {a.r1} {head}) <(zcat -f {a.r2} {head}) "
       f"2> {wd}/mem1.log | {tools['samtools']} view -u -G 12 - | {tools['samtools']} collate -O -u - {wd}/col "
       f"| {tools['samtools']} fastq -1 {wd}/s1.fq -2 {wd}/s2.fq -0 /dev/null -s /dev/null - 2>/dev/null")
    # reads 在 <(zcat ...) 里，读不出来（文件损坏、路径是目录）时 bwa 只收到空输入、照常退出 0：
    # 按 bwa 日志里实际读入的条数检查，0 条就报错，不把它当成"reads 太少"的线性判定
    n_read = sum(int(m) for m in re.findall(r"\[M::process\] read (\d+) sequences", Path(f"{wd}/mem1.log").read_text()))
    if n_read == 0:
        sys.exit(f"[ERROR] 从 {a.r1} / {a.r2} 读不到任何 reads（见 {wd}/mem1.log）")
    for it in items:
        k = it["i"]
        for tag, seq in (("rot", it["rot"]), ("lin", it["s"])):
            ref = f"{wd}/{tag}{k}.fasta"
            write_fasta(ref, [(f"{tag[0]}{k}", seq)])
            sh(f"{tools['bwa']} index {ref} 2>/dev/null && {tools['bwa']} mem -t {a.threads} {ref} "
               f"{wd}/s1.fq {wd}/s2.fq 2>/dev/null | {tools['samtools']} view -F 0x904 - > {wd}/{tag}{k}.sam")
        prof_rot.update(spanning_profile(f"{wd}/rot{k}.sam", {f"r{k}": len(it["rot"])}))
        Lr, J = len(it["rot"]), it["J"]
        # 对照点：接合处两侧 150/300/600 bp（局部，排除控制区 AT 富集、覆盖度低造成的偏差）+ 环上 4 个远处点（参考）
        local = [J + d for d in (-600, -300, -150, 150, 300, 600) if a.flank < J + d < Lr - a.flank]
        glob = [int((J + Lr * f) % Lr) for f in (0.2, 0.4, 0.6, 0.8)]
        glob = [p for p in glob if a.flank < p < Lr - a.flank]
        cnt = span_counts(f"{wd}/rot{k}.sam", f"r{k}", [J] + local + glob, a.flank)
        it["junction"] = cnt[J]
        it["local_med"] = statistics.median([cnt[p] for p in local]) if local else 0
        it["global_med"] = statistics.median([cnt[p] for p in glob]) if glob else 0
        prof_lin.update(spanning_profile(f"{wd}/lin{k}.sam", {f"l{k}": len(it["s"])}))
        clipst.update(end_clip_stats(f"{wd}/lin{k}.sam", {f"l{k}": len(it["s"])}))
        it["bridge"], it["gap_est"] = bridging_pairs(f"{wd}/lin{k}.sam", len(it["s"]))


def compare(seq, other_fa, blastn, workdir):
    q = os.path.join(workdir, "q.fasta")
    write_fasta(q, [("q", seq)])
    out = subprocess.run([blastn, "-query", q, "-subject", other_fa, "-outfmt",
                          "6 sseqid pident length qstart qend sstart send slen", "-evalue", "1e-20"],
                         capture_output=True, text=True).stdout.splitlines()
    if not out:
        return "no_hit"
    rows = [r.split("\t") for r in out]
    best = max(rows, key=lambda r: int(r[2]))
    cov = sum(int(r[2]) for r in rows if r[0] == best[0])
    return (f"{best[0].split()[0]}: {float(best[1]):.2f}% over {int(best[2])} bp; aligned {cov} bp of "
            f"{len(seq)} (this) vs {best[7]} (other)")


def tool(name, given):
    return given or shutil.which(name)   # in wgs_skim the rule's conda env provides bwa / samtools / blastn


# ── 主流程 ──────────────────────────────────────────────────────────────

def run_batch(a):
    """--summary 模式：按 header 的样本名拆分，每个样本单独调用本脚本（并行），最后合并 circularize.tsv。"""
    from concurrent.futures import ThreadPoolExecutor
    out = Path(a.outdir)
    out.mkdir(parents=True, exist_ok=True)
    groups = defaultdict(list)
    for n, q in read_fasta(a.summary):
        groups[n.split("|")[0].split()[0]].append((n, q))
    jobs = []
    for smp, recs in sorted(groups.items()):
        d = out / smp
        d.mkdir(exist_ok=True)
        old = (d / "input.fasta").read_text() if (d / "input.fasta").exists() else None
        write_fasta(d / "input.fasta", recs)
        if old is not None and old != (d / "input.fasta").read_text():   # 序列变了（如 mt_iter 换成迭代结果）→ 重算
            (d / "circularize.tsv").unlink(missing_ok=True)
        ann = " ".join(f"{f}:{f.stat().st_mtime_ns}" for x in (a.annotation or []) for f in
                       (sorted(Path(x).glob("*.gb")) if Path(x).is_dir() else [Path(x)]) if f.exists())
        sig = d / "annotation.sig"      # 注释变了（或第一次加注释）→ 重算
        if (sig.read_text() if sig.exists() else "") != ann:
            (d / "circularize.tsv").unlink(missing_ok=True)
            sig.write_text(ann)
        r1 = Path(a.fastp_dir) / smp / f"{smp}{a.suffix_fq[0]}"
        r2 = Path(a.fastp_dir) / smp / f"{smp}{a.suffix_fq[1]}"
        cmd = [sys.executable, __file__, "-i", str(d / "input.fasta"), "-o", str(d), "--prefix", smp,
               "-t", str(a.threads), "--validate_pairs", str(int(a.validate_pairs)),
               "--min_span_ratio", str(a.min_span_ratio), "--weak_span_ratio", str(a.weak_span_ratio)]
        if a.annotation:
            cmd += ["--annotation", *a.annotation]
        if r1.exists() and r2.exists():
            cmd += ["-1", str(r1), "-2", str(r2)]
        else:
            logging.warning(f"[CIRC] {smp}: 找不到 reads {r1}，只做序列层面的判断")
        jobs.append((smp, d, cmd))

    def one(job):
        smp, d, cmd = job
        if (d / "circularize.tsv").exists() and not a.force:
            return smp, 0
        (d / "circularize.tsv").unlink(missing_ok=True)   # 失败时不留旧结果冒充本次
        with open(d / "run.log", "w") as log:
            return smp, subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT).returncode

    failed = []
    with ThreadPoolExecutor(max(1, a.max_parallel)) as ex:
        for smp, rc in ex.map(one, jobs):
            logging.info(f"[CIRC] {smp}: rc={rc}")
            if rc != 0 or not (out / smp / "circularize.tsv").exists():
                failed.append(smp)
    rows, fields = [], None
    for smp, d, _ in jobs:
        f = d / "circularize.tsv"
        if f.exists():
            with open(f) as fh:
                r = csv.DictReader(fh, delimiter="\t")
                fields = fields or r.fieldnames
                rows += list(r)
    if fields:
        with open(out / "topology_summary.tsv", "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
            w.writeheader()
            w.writerows(rows)
        logging.info(f"[CIRC] {out / 'topology_summary.tsv'}: {len(rows)} 条序列")
    else:   # 没有任何结果：不留上次的汇总表
        (out / "topology_summary.tsv").unlink(missing_ok=True)
    if failed:   # 返回非 0：流程不写 topology.done，下次续跑只重算这些样本（其余已有 circularize.tsv）
        logging.error(f"[CIRC] 失败的样本（见 <样本>/run.log）：{', '.join(failed)}")
        return 1
    return 0


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-i", "--input", help="线粒体序列 fasta（一条或多条）")
    ap.add_argument("--summary", help="批量模式：汇总 fasta（header <样本>|...），配合 --fastp_dir")
    ap.add_argument("--fastp_dir", help="批量模式：clean reads 目录 <dir>/<样本>/<样本>.R1.clean.fastq.gz")
    ap.add_argument("--suffix_fq", nargs=2, default=[".R1.clean.fastq.gz", ".R2.clean.fastq.gz"])
    ap.add_argument("--max_parallel", type=int, default=4, help="批量模式同时检查几个样本")
    ap.add_argument("--force", action="store_true", help="批量模式：已有结果也重跑")
    ap.add_argument("-o", "--outdir", required=True)
    ap.add_argument("--prefix", default="sample")
    ap.add_argument("-1", "--r1", help="clean R1（给了才做 reads 检验）")
    ap.add_argument("-2", "--r2")
    ap.add_argument("-t", "--threads", type=int, default=8)
    ap.add_argument("-k", "--kmer", type=int, default=21, help="找首尾重叠用的 k-mer（默认 21）")
    ap.add_argument("--span", type=int, default=2000, help="开头取多少 bp 的 k-mer（默认 2000）")
    ap.add_argument("--end_window", type=int, default=5000, help="在末端多少 bp 内找（默认 5000）")
    ap.add_argument("--min_run", type=int, default=25, help="首尾一致片段最短长度（默认 25 bp）")
    ap.add_argument("--min_len", type=int, default=1000, help="环长下限（默认 1000）")
    ap.add_argument("--validate_pairs", type=float, default=1e7, help="检验用前多少对 reads（默认 1000 万；0 = 全部 reads）")
    ap.add_argument("--max_validate_pairs", type=float, default=4e7,
                    help="对照点跨越 reads 中位数 < --min_control 时自动加大 reads（×4），最多到此值（默认 4000 万）")
    ap.add_argument("--min_control", type=float, default=30, help="对照点跨越 reads 中位数的目标下限（默认 30）")
    ap.add_argument("--repeat_ratio", type=float, default=2.0,
                    help="接合处跨越 > 对照 × 此值判为 circular_repeat_junction：首尾相接但接合处在重复里，长度不确定（默认 2）")
    ap.add_argument("--min_bridge", type=int, default=5,
                    help="跨首尾 read 对（末端朝外 ↔ 开头朝外）≥ 此值（且 ≥ 对照 5%%）判为 linear_gap_circular（默认 5）")
    ap.add_argument("--min_end_reads", type=int, default=10, help="判'真实线性末端'时两端各至少多少条 reads（默认 10）")
    ap.add_argument("--flank", type=int, default=20, help="跨接合处两侧各至少多少 bp（默认 20）")
    ap.add_argument("--min_span_ratio", type=float, default=0.5,
                    help="接合处跨越 reads ≥ 对照点中位数 × 此值判为环形（默认 0.5）")
    ap.add_argument("--weak_span_ratio", type=float, default=0.2,
                    help="介于此值与 --min_span_ratio 之间判为 circular_weak（默认 0.2）。"
                         "实测：真环 0.56–1.02，重复区塌缩的 SPAdes 环 0.02，截断片段 0")
    ap.add_argument("--clip_incomplete", type=float, default=0.5,
                    help="末端 reads 外侧软剪切比例 ≥ 此值判为'没拼完'（默认 0.5）")
    ap.add_argument("--clip_terminus", type=float, default=0.2,
                    help="两端软剪切比例都 ≤ 此值判为'真实线性末端'（默认 0.2）")
    ap.add_argument("--compare", help="另一组装的 fasta（如 GetOrganelle 结果），比较长度与一致度")
    ap.add_argument("--legacy_linear_classes", action="store_true",
                    help="旧判定：按末端软剪切把首尾不相连的序列分成 linear_terminus / incomplete / undetermined"
                         "（线性对照校准后发现不可靠，默认改为 linear_unjoined + end_state）")
    ap.add_argument("--annotation", nargs="+",
                    help="注释 GenBank 文件或目录（如 06.mito_annotation/final_annotation）：按序列自动匹配；"
                         "有注释时在 nad5 中间切开检验首尾，并检查两端基因、相邻基因拼接点的 reads 支持")
    for t in ("bwa", "samtools", "blastn"):
        ap.add_argument(f"--{t}")
    a = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
    if a.summary:
        return run_batch(a)
    if not a.input:
        ap.error("需要 -i 或 --summary")
    tools = {t: tool(t, getattr(a, t)) for t in ("bwa", "samtools", "blastn")}
    out = Path(a.outdir)
    out.mkdir(parents=True, exist_ok=True)

    recs = read_fasta(a.input)
    ann = load_annotations(a.annotation) if a.annotation else []
    items = []
    for i, (name, s) in enumerate(recs):
        c = find_circle(s, a.kmer, a.span, a.end_window, a.min_run, a.min_len)
        circle = s[c["anchor"]:c["anchor"] + c["P"]] if c else s     # 无重叠：原样首尾相接
        genes = map_annotation(s, ann) if ann else None
        sp = gene_split_point(genes, len(s)) if genes else None
        if sp:      # 注释坐标（s）→ 环坐标；切点落在 c 去掉的重叠里时退回中点
            h = (sp[1] - c["anchor"]) % c["P"] if c else sp[1]
            split = f"{sp[0]}:{sp[1] + 1}"
        else:
            h, split = len(circle) // 2, "midpoint"
        h = h if 0 < h < len(circle) else len(circle) // 2
        items.append(dict(i=i, name=name, cid=name.split()[0], s=s, c=c, circle=circle, genes=genes, split=split,
                          rot=circle[h:] + circle[:h], J=len(circle) - h))

    have_reads = bool(a.r1 and a.r2 and tools["bwa"] and tools["samtools"])
    prof_rot, prof_lin, clipst = {}, {}, {}
    if have_reads:
        with tempfile.TemporaryDirectory(dir=out) as wd:
            rot_fa = f"{wd}/rot.fasta"      # 第 1 步只用来挑出线粒体 reads：所有序列放一起即可
            write_fasta(rot_fa, [(f"r{it['i']}", it["rot"]) for it in items] +
                        [(f"l{it['i']}", it["s"]) for it in items])
            sh(f"{tools['bwa']} index {rot_fa} 2>/dev/null")
            pairs = int(a.validate_pairs)
            while True:
                analyse_reads(items, a, tools, wd, rot_fa, pairs, prof_rot, prof_lin, clipst)
                low = min((it["local_med"] for it in items), default=0)
                if pairs <= 0 or low >= a.min_control or pairs >= a.max_validate_pairs:   # 0 = 已用全部 reads
                    break
                pairs = int(min(pairs * 4, a.max_validate_pairs))   # 对照点 reads 太少：加大 reads 再做一次
                logging.info(f"[CIRC] {a.prefix}: 对照点跨越 reads 中位数 {low:.0f} < {a.min_control}，改用前 {pairs:,} 对 reads")
            for it in items:
                it["pairs_used"] = pairs

    rows, circ_out, other_out, span_rows = [], [], [], []
    for it in items:
        s, c, L = it["s"], it["c"], len(it["circle"])
        row = dict(contig=it["cid"], length=len(s), method="end_overlap" if c else "as_is",
                   circle_length=L, overlap_bp=c["overlap"] if c else 0, agree_run_bp=c["run_len"] if c else "",
                   junction_reads="", control_reads_median="", global_control_median="", left_end_clip="",
                   right_end_clip="", bridging_pairs="", gap_estimate_bp="", zero_span="", pairs_used="",
                   split_at=it["split"], genes="", end_genes="", gene_across_ends="", gene_order_min_span="",
                   unsupported_adjacency="", end_state="", verdict="", compare="")
        genes = it["genes"]
        if genes:
            names = {g[1] for g in genes}
            ctx, cross = end_context(genes, len(s))
            row.update(genes=f"{len(names)} ({sum(g[0] == 'CDS' for g in genes)} CDS)", end_genes=ctx,
                       gene_across_ends=",".join(cross) or "none")
        if have_reads:
            j, med = it["junction"], it["local_med"]
            lc, ln, rcl, rn = clipst[f"l{it['i']}"]
            pl = prof_lin[f"l{it['i']}"]
            zr = zero_runs(pl, 150, max(150, len(pl) - 150))    # 内部（离两端 >150 bp）没有 reads 支持的拼接点
            row.update(bridging_pairs=it["bridge"],
                       gap_estimate_bp=f"{it['gap_est']:.0f}" if it["gap_est"] is not None else "",
                       junction_reads=j, control_reads_median=med, global_control_median=it["global_med"],
                       pairs_used=it["pairs_used"],
                       left_end_clip=f"{lc:.2f} (n={ln})" if lc is not None else "n=0",
                       right_end_clip=f"{rcl:.2f} (n={rn})" if rcl is not None else "n=0",
                       zero_span=";".join(f"{x + 1}-{y}" for x, y in zr[:10]) or "none")
            for p, v in enumerate(pl):
                span_rows.append((it["cid"], p + 1, v))
            if genes:
                adj = adjacency_support(genes, pl)
                if adj:
                    row.update(gene_order_min_span=min(v for _, v in adj),
                               unsupported_adjacency=";".join(k for k, v in adj if v == 0) or "none")
            if med < 5:
                verdict = ("circular" if c else "linear") + "_low_reads"
            elif j > a.repeat_ratio * med:
                verdict = "circular_repeat_junction"
            elif j >= a.min_span_ratio * med:
                verdict = "circular_validated"
            elif j >= a.weak_span_ratio * med:
                verdict = "circular_weak"
            elif c:
                verdict = "circular_rejected"
            elif it["bridge"] >= max(a.min_bridge, 0.05 * med):
                verdict = "linear_gap_circular"      # 首尾在分子上相邻：环形，组装中间缺约 gap_estimate_bp
            elif a.legacy_linear_classes:
                ends = [x for x in (lc, rcl) if x is not None]
                if any(x >= a.clip_incomplete for x in ends):
                    verdict = "linear_incomplete"
                elif ln >= a.min_end_reads and rn >= a.min_end_reads and all(x <= a.clip_terminus for x in ends):
                    # 环状注释能在首尾之间注释出完整基因 → 首尾在序列上是相连的，不当作线性末端
                    verdict = "linear_undetermined" if genes and row["gene_across_ends"] != "none" else "linear_terminus"
                else:
                    verdict = "linear_undetermined"
            else:
                # 首尾之间既无跨越 reads 也无跨首尾的 read 对 = Ahuja et al. 2024 判"线性"的标准。
                # 末端软剪切不能再细分"真末端 / 没拼完"：已证实线性的管水母 mt（dunnlab）末端外侧同样有大量
                # 未拼出的重复（Agalma clausi 右端 clip 0.86、reads 为对照 134 倍），见 end_state。
                # 短 reads 无法排除"环形但缺口长于插入片段（如长串联重复控制区）"，需长读长或 PCR 跨缺口确认。
                verdict = "linear_unjoined"
            if verdict.startswith("linear"):
                row["end_state"] = " | ".join(
                    f"{e}:{'clean' if x is not None and x <= a.clip_terminus else 'extends' if x is not None and x >= a.clip_incomplete else 'mixed' if x is not None else 'no_reads'}"
                    f",{n / max(med, 1):.1f}x" for e, x, n in (("L", lc, ln), ("R", rcl, rn)))
        else:
            verdict = "circular_unvalidated" if c else "linear_unvalidated"
        row["verdict"] = verdict
        if a.compare and tools["blastn"]:
            with tempfile.TemporaryDirectory(dir=out) as wd:
                row["compare"] = compare(it["circle"] if verdict.startswith("circular_v") or
                                         verdict == "circular_unvalidated" else s, a.compare, tools["blastn"], wd)
        if verdict in ("circular_validated", "circular_repeat_junction", "circular_unvalidated"):
            circ_out.append((f"{it['cid']} topology=circular length={L} method={row['method']} verdict={verdict}",
                             it["circle"]))
        else:
            other_out.append((f"{it['name']} verdict={verdict}", s))
        logging.info(f"[CIRC] {a.prefix} {it['cid']}: {verdict} | length {len(s)}, circle {L} ({row['method']}, "
                     f"overlap {row['overlap_bp']}, agree {row['agree_run_bp']}) | junction {row['junction_reads']} vs "
                     f"control {row['control_reads_median']} | end clip L {row['left_end_clip']} R {row['right_end_clip']}"
                     f" | zero-span {row['zero_span']}")
        rows.append(row)

    fields = ["sample"] + list(rows[0].keys()) if rows else ["sample"]
    with open(out / "circularize.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow({"sample": a.prefix, **r})
    write_fasta(out / f"{a.prefix}.circular.fasta", circ_out)
    write_fasta(out / f"{a.prefix}.other.fasta", other_out)
    if span_rows:
        with gzip.open(out / f"{a.prefix}.spanning.tsv.gz", "wt") as fh:
            fh.write("contig\tpos\treads_spanning_pos_to_next\n")
            for r in span_rows:
                fh.write("\t".join(map(str, r)) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
