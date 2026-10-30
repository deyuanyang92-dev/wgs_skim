# MEANGS：线粒体种子

从少量 reads（默认 200 万条）中快速组装出一段线粒体序列，作为 GetOrganelle 的种子。种子来自样品本身，
比 GetOrganelle 自带的通用种子更容易招募到 reads（尤其是与自带种子相距较远的类群）。

```bash
wgs_skim meangs -r clean_reads/ -o meangs_out/
wgs_skim meangs -r raw/ -o meangs_out/ --qc -t mollusc
```

- 类群（`--taxon` → `meangs_class`：A-worms、Mollusca、Chordata、Arthropoda …）决定 MEANGS 用的 HMM 蛋白集。
- 输出：`meangs/<s>/<s>.seed.fasta`；找不到时为空文件（GetOrganelle 将用自带种子）。
- 实测：0.3 Gb 样品 4–5 min、5.9 GB。配置：`meangs_reads`、`insert_size`。

一般不单独使用，而是 [MEANGS-GetOrganelle 组合](Combo-mitogenome.md)：`wgs_skim getorganelle --meangs`。

深入资料：[docs/tools/MEANGS.md](../tools/MEANGS.md)。
