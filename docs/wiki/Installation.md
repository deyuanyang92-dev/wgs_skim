# 安装

需要：Linux、conda / mamba（推荐 miniforge）、git。所有生物信息软件由 wgs_skim 自动装进各自独立的 conda 环境，不需要手动安装。

```bash
git clone https://github.com/deyuanyang92-dev/wgs_skim.git
cd wgs_skim
conda env create -f environment.yml      # Snakemake + wgs_skim 命令本身
conda activate wgs_skim
pip install .
wgs_skim install                         # 一次：创建全部软件环境，并下载 MitoM、MITOS2 参考、GetOrganelle 数据库
```

- `wgs_skim install` 首次约 20–60 分钟（视网速）。环境放在 `--conda-prefix`（默认 `~/.wgs_skim/conda`），数据库放在
  `~/.wgs_skim/db`，所有项目共用，以后不再重装。
- 中断后重新运行 `wgs_skim install` 即可，已完成的环境保留。
- 环境目录要有 > 20 GB 空间；可以放在大盘上：`wgs_skim install --conda-prefix /data/wgs_skim_conda`，或在
  `~/.wgs_skim/config.yaml` 写一行 `conda_prefix: /data/wgs_skim_conda`。

## 个人默认设置 `~/.wgs_skim/config.yaml`

每次运行都会读取，适合写本机固定的东西：

```yaml
conda_prefix: /data/wgs_skim_conda          # 软件环境位置
db_dir: /data/wgs_skim_db                   # 数据库位置（默认 ~/.wgs_skim/db）
blast_db: /data/db/metazoa_marker/metazoa_marker   # 本地 BLAST 库（见 BLAST 页）
email: you@example.com                      # NCBI 要求
busco_lineage: /data/busco_downloads/lineages/metazoa_odb12   # 已下载的 BUSCO 谱系
```

命令行参数优先于这里的设置；这里的设置优先于软件默认值。

## 检查安装

```bash
wgs_skim --version
wgs_skim samples -r /path/to/rawdata         # 能识别出样品名即可
python -m pytest -q tests/                   # 单元测试，几秒钟
```
