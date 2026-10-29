# MitoM：线粒体注释

作者自己的注释软件（<https://github.com/deyuanyang92-dev/Mitom>）：蛋白编码基因、rRNA、tRNA，以 cox1 为起点，
带 NCBI table2asn 检查，可直接用于提交。

```bash
wgs_skim mitom -i mitogenomes/ -o mitom_out/                 # 每个 FASTA 一个样品（可多条 contig）
wgs_skim mitom -i mitogenomes/ -o mitom_out/ -t fish         # 脊椎动物：遗传密码 2，Chordata 蛋白集
```

- 类群 `--taxon` 决定遗传密码与蛋白参考集（`clade`）；参考数据 MITOS2 refseq89m（`wgs_skim install` 下载一次）。
- 输出：`annotation/<s>/<s>.mt.gb`（GenBank）、`<s>.mitom.tsv`（基因数）、`annotation/<s>/mitom/out/`（全部结果）。
- 实测：一条线粒体 15 s – 2 min、0.3 GB。

`mitom nuc`（rDNA 18S-ITS-28S、组蛋白 H1–H4）用于核多拷贝基因，见 [核多拷贝基因](Combo-nuclear.md)。

深入资料：[docs/tools/MitoM.md](../tools/MitoM.md)。
