<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="assets/geno2flux-logo-dark.svg">
    <img alt="Geno2Flux" src="assets/geno2flux-logo-light.svg" width="420">
  </picture>
</p>

<p align="center"><a href="README.md">English</a> · <strong>简体中文</strong></p>

Geno2Flux 用来帮真菌代谢工程快速找到改造靶点：给它一个基因组，或者给它你手头已有的基因组尺度代谢模型（GEM），它就能跑出按优先级排好序的**过表达候选基因**、**单基因敲除候选**和**同源家族联合敲除候选**。

内置开源 [HiGHS](https://highs.dev/) 求解器（同时负责 LP 与 MILP），全程不需要购买商业 Gurobi 许可。

---

## 核心流程

1. **基因组到蛋白质组** — 调用 `ANNEVO` / `Helixer` / `Tiberius` / `funannotate` 从头预测基因，或直接导入已有的 GFF3 注释。
2. **模型重建与补洞** — 基于 AlphaGEM + Yeast9 构建草图模型，可选用 DeepLoc 细化亚细胞分室（没装会自动沿用默认分室，不会中断），再用 BLAST+ 与 HiGHS 做有基因组证据支撑的补洞。
3. **指定目标产物** — 直接选中模型里已有的反应，或挂载你自己的外源合成途径。
4. **靶点优先级排序** — 10 步 FSEOF 过表达扫描 → DIAMOND 同源家族聚类 → 单基因与同源家族联合敲除筛选 → 过滤掉影响生长的必需基因。

---

## 快速开始

先安装 [`uv`](https://docs.astral.sh/uv/)（会自动管理 Python 3.11–3.12），Linux 和 macOS 均可运行：

```sh
bash run.sh setup
make test
```

### 场景 A：手头已经有 GEM 模型，直接筛靶点

如果你已经有整理好的模型文件（SBML/YAML/JSON）和对应的蛋白质组 FASTA，可以跳过重建与补洞，直接做靶点筛选。

1. 安装用于同源家族聚类的 DIAMOND（也可以在配置里直接提供预计算的 `diamond_self_hits` 表）：
   ```sh
   bash run.sh fetch-resources --profile screening
   ```
2. 复制 [`configs/unified_pipeline.example.json`](configs/unified_pipeline.example.json) 为 `/work/run.local.json`，填好模型、蛋白质组、生长反应和目标反应：
   ```json
   {
     "localization_order": "skip",
     "gapfilled_model": "/work/model.xml",
     "proteome_fasta": "/work/proteome.fasta",
     "id_map": null,
     "biomass": "<你的生长反应 ID>",
     "target_reaction": "<你的目标反应 ID>",
     "compute_growth_essentiality": true,
     "allow_unverified_reuse": true
   }
   ```
   如果要挂载外源途径而不是选已有反应，在配置里写 `target_module` 即可（格式见 [配置与命令参考（英文）](docs/configuration.md)）。
3. 开始筛选：
   ```sh
   bash run.sh metabolic --config /work/run.local.json --output /work/run-001
   ```

### 场景 B：只有原始基因组，从头一路跑到靶点

从头执行 AlphaGEM 草图重建需要 **Linux x86_64** 以及外部数据库与独立环境，动手前先看 [环境搭建与外部资源（英文）](docs/setup.md)。

1. 由基因组生成蛋白质组与坐标输入（支持已有 GFF3，或直接调用 `annevo` / `helixer` / `tiberius` / `funannotate`）：
   ```sh
   bash run.sh predict-proteome \
     --tool gff \
     --genome /work/raw/genome.fasta \
     --gff /work/raw/prediction.gff3 \
     --case-id sample_strain \
     --prefix PROT \
     --output /work/pred
   ```
2. 在 `/work/run.local.json` 里填入 `/work/pred/inputs/` 和 AlphaGEM 环境路径，然后运行：
   ```sh
   bash run.sh metabolic --config /work/run.local.json --output /work/run-001
   ```

想确认当前环境缺什么，或查看全部子命令，随时运行 `bash run.sh check-env` 和 `bash run.sh help`。

---

## 跑完会拿到什么

结果统一输出在 `/work/run-001/<order>/` 下：

| 文件 | 业务用途 |
| --- | --- |
| `fseof/all_genes_fseof_essential_filtered.tsv` | 去掉生长必需基因后的过表达候选排序（按 10 步 FSEOF 斜率） |
| `paralog/paralog_family_rankings_essential_filtered.tsv` | DIAMOND 聚类后的同源基因家族排序（按家族聚合 FSEOF 得分） |
| `knockout/single_gene_ko_essential_filtered.tsv` | 单基因敲除候选清单（已过滤必需基因） |
| `knockout/family_ko_essential_filtered.tsv` | 同源基因家族联合敲除清单（已过滤必需基因） |
| `target_model_<order_slug>.xml` | 已挂载目标产物或外源途径的最终 SBML 模型 |

> 以上四张 `*_essential_filtered.tsv` 筛后表需要开启 `"compute_growth_essentiality": true`（上面最小配置和模板均已开启）或提供 `essential_genes_tsv` 才会生成；若两项都未配置，只输出未过滤的原始排序表。

---

## 延伸文档

- **[环境搭建与外部资源（英文）](docs/setup.md)** — 各阶段平台支持、一键下载 `fetch-resources`、DeepLoc 选配与降级、Linux AlphaGEM / DeepEC 环境说明。
- **[配置与命令参考（英文）](docs/configuration.md)** — 外源途径配置（`target_module`）、SHA-256 校验规则、逐阶段 `stage` 命令与仓库目录结构。
- **[第三方软件与数据声明](THIRD_PARTY.md)** · **[许可证（GPL-3.0）](LICENSE)**
