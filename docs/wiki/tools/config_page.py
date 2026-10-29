"""Configuration.md from config/config.yaml (keeps the wiki table in step with the defaults):
    python docs/wiki/tools/config_page.py > docs/wiki/Configuration.md"""
import re
from pathlib import Path

HEAD = """# 配置

所有参数都在 `config/config.yaml`（下表由它生成：`python docs/wiki/tools/config_page.py`）。三种方式修改，后者优先：

1. `~/.wgs_skim/config.yaml`：个人长期设置（数据库目录、`blast_db`、`email`、`mitom_dir`），每次运行自动读取；
2. 命令行选项：`--taxon`、`--steps`、`--spades-modes`、`--blast-db`、`--mem-gb` ……（`wgs_skim run -h`）；
3. `--config 键=值 键=值`：任何一个键，例如 `--config mt_iter_max_data_gb=0 nuclear_max_gb=5,10`。

线程数和内存不在这里设置：`-j`（总核数）、`--mem-gb`（总内存），每个任务的用量见 [资源与续跑](Resources-and-resume)。"""

ABOVE = {   # keys whose explanation is the comment block above them
    "taxon": "annelid / mollusc / fish / crustacean：同时设定 MEANGS 类群、线粒体遗传密码、MitoM 类群、BUSCO 谱系（下面 4 个键可单独覆盖）",
    "steps": "只运行部分步骤（逗号分隔），所需的上游步骤自动加入；单个步骤 fastp meangs mito getorganelle nr spades rnaspades busco trinity rna genome_size annotate topology mitofinder blast nuclear summary，或组合 meangs-getorganelle spades-mitofinder spades-busco fastp-spades-busco rnaspades-busco trinity-busco（见 [组合流程](Combinations.md)）；空 = 下面打开的全部",
}
lines = (Path(__file__).resolve().parents[3] / "config" / "config.yaml").read_text().splitlines()
out, in_section = [HEAD], False
for i, line in enumerate(lines):
    m = re.match(r"# ── (.+?) ─", line)
    if m:
        out += ["", f"## {m.group(1)}", "", "| 键 | 默认 | 说明 |", "|---|---|---|"]
        in_section = True
        continue
    m = re.match(r'^([a-z_]+):\s*("[^"]*"|\S*)\s*(?:#\s*(.*))?$', line)
    if m and in_section:
        key, value, note = m.group(1), m.group(2) or '""', m.group(3) or ABOVE.get(m.group(1), "")
        j = i + 1
        while j < len(lines) and re.match(r"^\s{10,}#\s", lines[j]):      # continuation comment lines
            note += " " + lines[j].split("#", 1)[1].strip()
            j += 1
        out.append(f"| `{key}` | `{value}` | {note.replace('|', '/')} |")
print("\n".join(out))
