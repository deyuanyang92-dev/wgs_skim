# BLAST+ 与 NCBI Entrez（wgs_skim 排错参考）

## 1. 简介与用途

物种核对（rule `blast`，`workflow/rules/blast.smk` → `workflow/scripts/blast_taxonomy.py`）：把每个样品选出的线粒体
（默认 cox1）和 rDNA（默认 18S、ITS、28S）序列用 `blastn` 比对 nt（本地库或 `-remote`），每条查询取 TopN 个不同 accession，
再用 Biopython `Entrez.efetch` 取 GenBank 注释头（organism、taxonomy、definition、reference），写 `blast/blast_topn.tsv`。

`wgs_skim build-blastdb`（`src/wgs_skim/cli.py` → `workflow/scripts/build_marker_db.py`）：从本地 nt 抽出一个分类群
（默认 Metazoa 33208）、长度 ≤ 50 kb 的记录，用 `blastdbcmd` + `makeblastdb -parse_seqids -taxid_map` 建一个小的标记库。

链接：
- 手册：https://www.ncbi.nlm.nih.gov/books/NBK279690/ ；发布说明：https://www.ncbi.nlm.nih.gov/books/NBK131777/
- 下载：https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/ ；数据库：https://ftp.ncbi.nlm.nih.gov/blast/db/
- 文档仓库：https://github.com/ncbi/blast_plus_docs
- 论文：Camacho et al. 2009 BMC Bioinformatics 10:421, doi:10.1186/1471-2105-10-421
- URL API 使用规范：https://blast.ncbi.nlm.nih.gov/doc/blast-help/developerinfo.html
- Entrez 使用规范：https://www.ncbi.nlm.nih.gov/books/NBK25497/

## 2. 版本

| | 版本 |
|---|---|
| wgs_skim 固定 | `blast =2.16.0`、`biopython =1.84`、`python =3.11`（`workflow/envs/blast.yaml`） |
| 最新 | 2.17.0（2025-07-22，FTP `LATEST/` 核对 2026-10-29） |

- 2.16.0（2024-07-22）：blastp-fast / blastx-fast、PSI-BLAST 多线程、SRA 过滤选项；对 blastn 无关键变化。
- 2.17.0：makeblastdb 可直接读压缩 FASTA；新增带表头的 CSV 输出格式；源码编译需 zlib / zstd / bzip2。
- 是否升级：**不需要**。2.17.0 对 wgs_skim 唯一可能的好处是 `build_marker_db.py` 可直接喂 `.gz`，但当前是管道喂
  makeblastdb，无收益。

数据库现状（`*-nucl-metadata.json`，2026-10-24/25）：

| 库 | 碱基 | 序列数 | 压缩包 | 解压后（bytes-total） |
|---|---|---|---|---|
| nt | 4.82 T | 132.5 M | 404 卷 | 1.24 TB |
| core_nt | 1.09 T | 132.3 M | 92 卷 | 0.31 TB |
| 本机 metazoa_marker（2026-10-28 建） | 204 G | 71.3 M | 30 卷左右 | ~68 GB |

**core_nt = nt 去掉大部分真核染色体级记录**，序列数几乎相同、体积约 1/4。自 2024-09 起 NCBI 把非交互式（URL API）对
`nt` 的搜索**自动改为 core_nt**（2024 BLAST News）；推测 `blastn -remote -db nt` 也属于此类（未逐字核实）。

## 3. wgs_skim 中的调用

```bash
blastn -query blast/query.fasta \
  -outfmt "6 qseqid sacc stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs slen" \
  -evalue 1e-10 -max_target_seqs <max(50, 10×topn)> -out blast/blast.out \
  -db <blast_db> -num_threads 8          # 本地
  # 或  -db nt -remote                      # blast_db 为空时
```

之后每个 (qseqid, sacc) 只留 bitscore 最高的 HSP，按 bitscore 取前 `blast_topn` 个；`Entrez.efetch(db=nucleotide,
rettype=gb, seq_start=1, seq_stop=1)` 每 100 个 accession 一批，失败重试 3 次（间隔 10/20/30 s），仍失败写 `(efetch failed)`。

可配置键：`blast`（开关）、`blast_db`（本地库前缀；空 = 在线）、`blast_mt`（full / cox1 / 16s / 12s / cytb / nad1 …）、
`blast_nr`（full / 18S / ITS / ITS1 / 5.8S / ITS2 / 28S）、`blast_topn`（默认 5）、`email`（Entrez 必填）。
CLI：`--blast-db`、`--email`。资源：`threads: 8`、`mem_mb: 48000`。

`build-blastdb`：`blastdbcmd -entry all -outfmt "%o %a %l %T"`（列表）→ 筛 taxid（`taxonomy4blast.sqlite3` 的 TaxidInfo 表求
后代）和长度 → `blastdbcmd -entry_batch ... -outfmt %f` → 按 OID 去重 → `makeblastdb -dbtype nucl -parse_seqids -taxid_map`。

## 4. 重要参数

| 参数 | 含义 | 说明 |
|---|---|---|
| `-task` | 默认 `megablast`（word 28） | wgs_skim 未设，即 megablast。跨种（同属 / 同科，COI 相似度 80–90 %）时 megablast 可能漏掉远缘命中；`dc-megablast` 专为跨种比较设计 |
| `-max_target_seqs` | 保留的目标序列数 | 不是"最好的 N 个"：它影响搜索过程（Shah et al. 2019，doi:10.1093/bioinformatics/bty833）。wgs_skim 取 ≥50 再自己排序，做法正确；不要设 < 5 |
| `-evalue 1e-10` | 阈值 | 18S / 28S 很保守，够用 |
| `qcovs` | 查询覆盖度（所有 HSP） | 判断"只命中一段"时有用 |
| `staxids sscinames` | 命中的 taxid / 学名 | 需要 `taxdb.btd/.bti` 在库目录或 `$BLASTDB`；可替代一部分 Entrez 查询（wgs_skim 未用） |
| `-taxids` / `-taxidlist` | 只搜某分类群 | 需库带 taxid（v5 库）；实测不减少本地 nt 时间（见下） |
| `-num_threads` | 线程 | `-remote` 时不能用（wgs_skim 只在本地给） |
| `-remote` | NCBI 服务器上搜 | 有 CPU 时间上限和速率限制 |

## 5. 已知问题与 issues

| 现象 / 报错 | 原因 | 解决 | 链接 |
|---|---|---|---|
| 本地 nt 一次搜索 40–50 分钟，与线程数、`-taxids` 无关 | nt 每次搜索都整库读，受磁盘限制 | `wgs_skim build-blastdb` 建标记库（22–121 s，命中相同） | 本项目 MEMORY / build_marker_db.py |
| `makeblastdb: Duplicate seq_ids` | nt v5 中相同序列只存一次、带多个 accession（NR_003027.2 有 920 个） | 按 OID 去重（已修） | build_marker_db.py 文档串 |
| `CPU usage limit was exceeded`（remote） | 每次远程搜索有固定 CPU 配额；查询多 / 长 / 低复杂度 | 少量查询、本地库 | https://support.nlm.nih.gov/knowledgebase/article/KA-03285/en-us ；https://biostars.org/p/9488143 |
| remote 部分查询被跳过、结果少 | 服务器限流 / 超时 | 本地库；或分批 | https://www.biostars.org/p/463990 、https://biostars.org/p/106819 |
| remote 很慢或被降速 | >100 次搜索 / 24 h 进入慢队列；轮询每 RID ≥ 1 min；提交间隔 ≥ 10 s；>50 次请在美东夜间 / 周末 | wgs_skim 一个批次只提交一次（所有查询一个文件），符合规范 | developerinfo.html |
| remote 结果与本地 nt 不同 | 远程 nt 实为 core_nt（2024-09 起） | 报告中注明数据库 | 2024 BLAST News |
| "最佳命中"不是最好的 | `-max_target_seqs` 很小时的提前截断 | ≥ 50 再按 bitscore 排序 | Shah et al. 2019；https://biostars.org/p/340129/ |
| `sscinames` 为 `N/A` | 找不到 taxdb | 把 taxdb.btd/bti 放进库目录或设 `BLASTDB` | NBK279690 |
| Entrez `HTTP 429` / 超时 | 无 API key 限 3 次/秒，有 key 10 次/秒；染色体级记录整条下载 | 只取注释头（seq_start=seq_stop=1，已修）；批量 100 个；可加 `Entrez.api_key` | NBK25497 |
| BLAST 进程 RSS ≈ 库大小（47 GB） | 库文件 mmap，页缓存计入 RSS | 属文件页，可回收；见第 6 节 | — |

## 6. 本项目中遇到的问题

- **TESTING.md #7**：efetch 下载整条序列（染色体级 rDNA 命中）超时 → organism 静默为空。修复：`seq_start=1, seq_stop=1` 只取注释头；失败时写 `(efetch failed)`。
- **TESTING.md #8**：样品名互为前缀（X 与 X.1）时用错样品序列 → 按样品目录对应。
- 本地 nt（本机 3.4 Tbases 版本）一次 COX1 + rDNA 搜索 42–52 min；改用 metazoa_marker（68 GB）后 22–121 s，同样的命中
  （MJ20260820041）。e2e2 实测：153 s，峰值 RSS 46.5 GB、`io_in` 30 GB。
- `mem_mb: 48000` 的依据是该 RSS；但 mmap 的库页是可回收的页缓存（推断，基于 max_uss ≈ max_rss 且无大块匿名内存的工作方式），
  真实"必须"内存远小于此。Snakemake 8 会把超出 `--resources mem_mb` 的声明**静默降到总量**（本机 8.18.2 实测：声明 48000、
  总量 20000 → 以 20000 运行），所以小内存机器上不会卡死，只是调度时独占内存。
- 早期（2026-10-27/28，run1/run5/run6）用旧的 `blast2metadata` 封装时 rule blast 失败数次，后改为现在的 `blast_taxonomy.py`。
- 本机环境变量 `BLASTDB=/media/cody/c9ef9bb9-.../Blastdb/blast` **指向不存在的目录**（2026-10-29 核对）；taxdb.btd/bti 实际在
  `/media/cody/18T_R2/Common_databases/Blastdb/blast/`。若以后用 `sscinames`，需在规则中显式设置 `BLASTDB`。

## 7. 排错清单

1. `logs/blast.log` 第一行是实际 `blastn` 命令（含 -db 路径）；最后一行 `N hits for M sequences`。
2. `blast/query.fasta`：查询是否为空、是否用了整条序列（日志中 `... not in the mt annotation: whole sequence used`）。
3. `blast/blast.out` 为空：数据库路径错（`BLAST Database error: No alias or index file found`）或 evalue 太严。
4. 表中 `(efetch failed)`：日志搜 `efetch failed`；网络恢复后 `-- --forcerun blast`。
5. 慢：`benchmarks/blast.tsv` 的 `io_in`；本地 nt → 换标记库；remote → 看 NCBI 状态与提交次数。
6. 库检查：`blastdbcmd -db <prefix> -info`（序列数、版本 v5、是否带 taxid）；`blastdbcheck -db <prefix>`。

## 8. 资料来源（查阅日期 2026-10-29）

- https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ ；https://ftp.ncbi.nlm.nih.gov/blast/db/nt-nucl-metadata.json 、core_nt-nucl-metadata.json
- https://blast.ncbi.nlm.nih.gov/doc/blast-news/2024-BLAST-News.html （2.16.0、core_nt 默认）
- BLAST+ 2.17.0 发布说明（NBK131777，经搜索摘要；页面本身需人机验证，未直接打开）
- https://blast.ncbi.nlm.nih.gov/doc/blast-help/developerinfo.html
- Shah N et al. 2019 Bioinformatics 35(9):1613, doi:10.1093/bioinformatics/bty833；Biostars 340129、392843
- 本项目：`workflow/rules/blast.smk`、`workflow/scripts/blast_taxonomy.py`、`build_marker_db.py`、`docs/TESTING.md` §4、
  `/media/cody/18T_R2/Common_databases/Blastdb/metazoa_marker/build.log`、`../wgs_skim_test/e2e2/benchmarks/blast.tsv`
