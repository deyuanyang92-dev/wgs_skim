# Jellyfish + GenomeScope2：基因组大小

```bash
wgs_skim genome-size -r clean_reads/ -o gs_out/
wgs_skim genome-size -r clean_reads/ -o gs_out/ --config kmer_size=25 ploidy=2
```

- jellyfish 统计 k-mer（`kmer_size` = 21），GenomeScope2 拟合：基因组大小、杂合度、重复比例。
- 输出：`genomescope/<s>/<s>.k21.histo`、`<s>.genomescope.txt`、`genomescope/`（图）。
- 数据量要求：至少约 15× 覆盖度才可靠；浅测序（< 10×）样品的结果仅供参考，GenomeScope 可能拟合失败。
- 实测（0.3 Gb）：30 s – 2 min、6.8 GB（结果不可信：覆盖度太低）。

深入资料：[docs/tools/GenomeScope.md](../tools/GenomeScope.md)。
