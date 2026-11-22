<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="assets/geno2flux-logo-dark.svg">
    <img alt="Geno2Flux" src="assets/geno2flux-logo-light.svg" width="420">
  </picture>
</p>

<p align="center"><a href="README.md">English</a> · <strong>简体中文</strong></p>

Geno2Flux 基于你自己的基因组/蛋白质组数据，为真核（真菌）菌株重建基因组尺度代谢模型（GEM），并通过双指标体系优先排序代谢工程改造靶点：

1. **基因组到 GEM 的重建与细化**：基因组→蛋白质组预测与编码 GFF3 标准化（`ANNEVO` / `Helixer` / `Tiberius` / `funannotate` + RefSeq/UniProt 蛋白证据，或直接导入已有 GFF3）→ AlphaGEM + Yeast9 草图重建 → 亚细胞定位分室细化（自动检测 DeepLoc 2.0/2.1、预计算表格，或自动降级保留默认分室）→ 基因组证据支撑的最小集合补洞（NCBI BLAST+ + HiGHS）。
2. **双指标靶点优先级排序**：可配置的产物途径或内源反应靶点 → 10 步 RAVEN 兼容 FSEOF 过表达扫描 → DIAMOND 旁系同源家族聚类 → 单反应、单基因及同源家族联合敲除筛选 → 生长必需性过滤。

---

## 1. 内置能力、外部依赖与平台边界

| 阶段 | 克隆 + `bash run.sh setup` 后即可运行 | 需要的外部依赖 | 平台边界 |
| --- | --- | --- | --- |
| 环境同步与测试（`make test`） | 是：uv 锁定环境 + 内置开源 HiGHS（`highspy`，覆盖 LP 与 MILP，无需 Gurobi） | 需先安装 [`uv`](https://docs.astral.sh/uv/)：它会自动管理 `bash run.sh setup` 所需的 Python 3.11–3.12 | Linux / macOS |
| GFF3 标准化与蛋白质组准备（`predict-proteome --tool gff`、`prepare-inputs`） | 是：内置纯 Python CDS 拼接器，`gffread` 缺失时自动回退 | 从头注释时才需要 ANNEVO / Helixer / Tiberius / funannotate | Linux / macOS |
| 基于已有 GEM 的靶点优先级排序（`target`/`fseof`/`paralogs`/`knockout`/`filter`） | 部分：靶点挂载、10 步 FSEOF、敲除筛选与生长必需性过滤内置（HiGHS LP） | 旁系同源聚类需 DIAMOND 2.1.10（或提供预计算的 `diamond_self_hits`） | Linux / macOS |
| 分室细化（`localize`） | 部分：未提供 `localization_tsv` 且 DeepLoc 缺失时自动降级，保留模板默认分室 | DeepLoc 2.0/2.1（DTU 个人学术许可，见第 4 节） | Linux / macOS |
| 最小集合补洞（`gapfill`） | 求解器内置（HiGHS MILP） | NCBI BLAST+ 2.16.0（`makeblastdb`、`tblastn`；`fetch-resources --profile gapfill` 校验下载）。`--profile gapfill` 因包含 `screening` 也会顺带安装 DIAMOND，但 `gapfill` 本身不依赖 DIAMOND | 官方二进制覆盖 linux-x86_64、linux-aarch64、darwin-aarch64、darwin-x86_64 |
| AlphaGEM 草图重建（从头重建步骤 1–3，`--profile paper`） | 否：不可裸跑 | 独立 AlphaGEM + DeepEC Python 环境、~13 GB Figshare 数据、ESM-1b 权重、eggNOG 数据库与契约匹配的 `fungi.dmnd`（见第 5 节） | 仅 Linux x86_64 |

**一句话边界**：[`uv`](https://docs.astral.sh/uv/)（自动管理 Python 3.11–3.12）是 `bash run.sh setup` 唯一的先决条件；在 Linux 与 macOS 上，克隆仓库并执行 `bash run.sh setup` 后即可直接运行测试套件、GFF3 标准化以及基于 HiGHS LP 的靶点阶段（`target`/`fseof`/`knockout`/`filter`；`localize` 默认降级保留原分室）；在已有 GEM 上跑完含旁系同源聚类的完整靶点筛选（`localization_order: "skip"` + `gapfilled_model`）还需 DIAMOND 2.1.10 或预计算的 `diamond_self_hits`，`gapfill` 还需 NCBI BLAST+ 2.16.0（`makeblastdb`、`tblastn`；`fetch-resources --profile gapfill` 会安装 BLAST+，并因该 profile 包含 `screening` 而一并安装 DIAMOND）；**从原始基因组出发的 AlphaGEM 从头重建（步骤 1–3）无法从裸克隆直接运行**，严格依赖第 5 节所述外部资源与 Linux x86_64 独立环境。

---

## 2. 安装与验证

使用 [`uv`](https://docs.astral.sh/uv/) 同步锁定的 Python 环境（Python 3.11–3.12，`>=3.11,<3.13`；`setup` 内部执行 `uv sync --frozen --extra test`），并通过自动化测试套件验证安装：

```sh
bash run.sh setup
make test
```

`make test` 执行 `uv run --frozen --extra test python -m pytest -q`。同步后的环境已内置开源 [HiGHS](https://highs.dev/) 求解器（`highspy`），同时覆盖 **LP**（FSEOF 扫描、敲除筛选、生长必需性过滤）与 **MILP**（最小集合补洞），无需商业 Gurobi 许可。

随时检查各阶段就绪状态（Python 解释器、外部二进制、已下载资源）：

```sh
bash run.sh check-env
```

查看全部 CLI 入口：

```sh
bash run.sh help
```

---

## 3. 靶点配置（已有 GEM 或自定义途径）

Geno2Flux 不预设产物、途径或菌株。必须在 JSON 配置（或 `scripts/run_stage.py target`）中为第 6 步选择两种方式之一：

1. **模型中已有的反应**（不新增反应）：
   ```json
   {
     "target_reaction": "r_2056",
     "target_module": null
   }
   ```

   `r_2056` 是捆绑 Yeast9 模型中琥珀酸交换反应的真实 ID（BiGG 注释为 `EX_succ_e`）。反应 ID 因模型而异，BiGG/注释别名**不能**作为反应 ID 使用。检查其他 SBML 模型的候选交换反应：
   ```sh
   uv run --frozen python - <<'PY'
   import cobra
   model = cobra.io.read_sbml_model("path/to/model.xml")
   for reaction in model.exchanges:
       print(reaction.id, reaction.name, reaction.bounds, reaction.annotation)
   PY
   ```
2. **自定义化学计量模块**（`"target_module"` 内联对象或 `.json` 文件路径；内部反应在筛选前会通过元素与电荷平衡校验）：
   ```json
   {
     "target_module": {
       "mode": "custom",
       "target_reaction": "DM_target_c",
       "expected_metabolites": {
         "chebi30031_c": { "formula": "C4H4O4", "charge": -2, "compartment": "c" }
       },
       "new_metabolites": [
         { "id": "target_product_c", "formula": "C4H4O4", "charge": -2, "compartment": "c" }
       ],
       "reactions": [
         { "id": "TARGET_SYN", "stoichiometry": { "chebi30031_c": -1.0, "target_product_c": 1.0 } },
         { "id": "DM_target_c", "stoichiometry": { "target_product_c": -1.0 } }
       ]
     }
   }
   ```

**已有 GEM 的快捷路径**：如果你已有整理好的 SBML/YAML/JSON GEM 和蛋白质组 FASTA，希望跳过 AlphaGEM 重建、定位与补洞直接做靶点排序，设置 `"localization_order": "skip"`、`"gapfilled_model": "path/to/model.xml"`、`"id_map": null`（基于示例配置修改时避免引用缺失的 `inputs/id_map.tsv`），并把 `"biomass"` 与 `"target_reaction"`（或 `"target_module"`）改为你模型中的实际反应 ID（示例配置默认为 Yeast9 的 `r_4041` / `r_2056`），同时配上对应的 `*_sha256` 或 `"allow_unverified_reuse": true`。此时旁系同源聚类仅需 DIAMOND（或提供预计算的 `diamond_self_hits` 以完全绕过 DIAMOND）。

**复用哈希校验门**：复用 `proteome_fasta`、`withgaps_model`、`gapfilled_model`、`localization_tsv`、`template_model` 或 `diamond_self_hits` 时，必须提供匹配的 `*_sha256` 字段（如 `proteome_sha256`、`gapfilled_model_sha256`、`diamond_self_hits_sha256` 等），除非显式设置 `"allow_unverified_reuse": true`。

---

## 4. DeepLoc 即插即用槽位与自动降级

DeepLoc 2.0 / 2.1 采用 DTU Health Tech 的个人学术许可，不能随公开仓库再分发，请直接向 DTU 申请：[DeepLoc 2.1 在线服务](https://services.healthtech.dtu.dk/services/DeepLoc-2.1/) · [2.1 下载申请表单](https://services.healthtech.dtu.dk/cgi-bin/sw_request?software=deeploc&version=2.1&packageversion=2.1&platform=All)，或 [DeepLoc 2.0 在线服务](https://services.healthtech.dtu.dk/services/DeepLoc-2.0/) · [2.0 下载申请表单](https://services.healthtech.dtu.dk/cgi-bin/sw_request?software=deeploc&version=2.0&packageversion=2.0&platform=All)。两种情况 Geno2Flux 都能无缝处理：

- **有 DeepLoc（零配置自动检测）**：拿到 DTU tarball 后执行 `bash run.sh fetch-resources --profile deeploc --deeploc-tarball /path/to/DeepLoc-2.1.tar.gz`，一键安装到 [`tools/deeploc/`](tools/deeploc/README.md)；也可以手动把 `deeploc2` 放/软链到 `tools/deeploc/` 下（如 `tools/deeploc/bin/deeploc2` 或 `tools/deeploc/deeploc2`）、设置 `G2F_DEEPLOC_BIN=/path/to/deeploc2`、将 `deeploc2` 加入 `PATH`，或直接提供预计算的 `localization_tsv`。
- **无 DeepLoc（自动降级）**：未提供 `localization_tsv` 且未安装 DeepLoc 时，Geno2Flux **不会崩溃**——记录一条提示、保留草图/模板模型的默认分室（`localization_refinement_audit.json` 中标记 `"deeploc_status": "fallback_default_compartments"`），其余流程照常运行。如需严格失败，在配置中设置 `"require_deeploc": true`（或在 `run_stage.py localize` 上加 `--require-deeploc`）。

---

## 5. 外部资源与 AlphaGEM 环境搭建

### 一键下载与就绪预检

`bash run.sh fetch-resources` 仅从各上游官方主机下载外部工具与数据库（`--dry-run` 可先预览全部动作）。校验强度按资产类别区分：

- **固定 SHA-256 / MD5 摘要**：Figshare `data_source.tar.xz`（SHA-256 + MD5 + 字节数）、NCBI BLAST+ 2.16.0 归档（MD5）、ESM-1b 接触回归文件（MD5 + 字节数），以及 eggNOG 系列——`eggnog.db`（解压后校验 SHA-256 + 字节数）、`e5.proteomes.faa` 与 `e5.taxid_info.tsv`（SHA-256）、过滤后的真菌 FASTA（对照固定契约 SHA-256）、构建出的 `fungi.dmnd`（SHA-256 + 字节数）。
- **仅字节数校验**（上游未发布摘要）：ESM-1b 主权重（`esm1b_t33_650M_UR50S.pt`，约 7.8 GB）与 `ANNEVO_Fungi.pt`。
- **上游未发布校验和**（以下载/解压成功为凭）：DIAMOND v2.1.10 Linux 压缩包、ANNEVO GitHub `main` 源码归档、`eggnog.db.gz` 压缩包（其**解压产物** `eggnog.db` 另行 SHA-256 校验）、可选 Helixer 权重、HuggingFace `prot_t5_xl_uniref50`。

模块化 `--profile` 取值——`screening`、`gapfill`、`alphagem`、`annotation`、`deeploc`、`all`——层层蕴含：`gapfill` 包含 `screening`，`alphagem` 包含 `gapfill`：

```sh
# 仅筛选：DIAMOND 装入 tools/bin/
bash run.sh fetch-resources --profile screening

# 补洞：DIAMOND + NCBI BLAST+ 2.16.0（tblastn、makeblastdb、blastp）装入 tools/bin/
bash run.sh fetch-resources --profile gapfill

# AlphaGEM 数据（从头重建）：暂存 AlphaGEM 源码快照 + 下载并校验
# ~13 GB Figshare data_source.tar.xz（data_available/、rhea/rhea_protembedding1~6.pkl、
# rhea2uniprot_sprot.tsv、rhea.rdf、plmsearchtools/、tools/）+ ESM-1b 权重 + eggNOG
#（--with-eggnog：~41 GB eggnog.db + 契约匹配的 fungi.dmnd；仅已暂存过才可省略）
# 仅在 alphagem_mode=structure_alignment 时追加 --with-prott5
bash run.sh fetch-resources --profile alphagem --runtime-root /work/alphagem-runtime --with-eggnog

# 从头基因组注释：ANNEVO 源码 + ANNEVO_Fungi.pt（可追加 --with-helixer）
bash run.sh fetch-resources --profile annotation --tool-root /work/annotation-tools
```

各 profile 明细：

- **`screening`**：下载 DIAMOND v2.1.10 至 `tools/bin/`（仅 `linux-x86_64` 有官方发布二进制）。macOS 与 `linux-aarch64` 无官方 DIAMOND 二进制，`fetch-resources` 会打印手动安装指引（`conda install -c conda-forge -c bioconda diamond=2.1.10` 或 `brew install diamond`）。
- **`gapfill`**：DIAMOND + NCBI BLAST+ 2.16.0（`tblastn`、`makeblastdb`、`blastp`）至 `tools/bin/`，官方归档覆盖 `linux-x86_64`、`linux-aarch64`、`darwin-aarch64`、`darwin-x86_64`。
- **`alphagem`**（`--runtime-root /work/alphagem-runtime --with-eggnog`，可选 `--with-prott5`）：
  - 通过 `prepare_alphagem_runtime.py` 暂存已校验的 `third_party/alphagem/` 源码快照，并初始化为干净的本地纯源码 git 提交；
  - 下载并校验约 13 GB 的 Figshare `data_source.tar.xz`（<https://doi.org/10.6084/m9.figshare.30630524>），解出 `data_available/`（`yeast.fasta`、`yeast.xlsx`）、`rhea/rhea_protembedding1~6.pkl`、`rhea/rhea2uniprot_sprot.tsv`、`rhea/rhea.rdf`、`plmsearchtools/` 与 `tools/`；
  - 下载 ESM-1b 权重（`esm1b_t33_650M_UR50S.pt` 约 7.8 GB + `esm1b_t33_650M_UR50S-contact-regression.pt`）；
  - `--with-eggnog`：下载解压 `eggnog.db`（约 41 GB，SHA-256 校验）并在 `tools/eggnog_mapper/data/` 下构建 `fungi.dmnd`（`879,076,733` 字节，SHA-256 校验）。中间产物（`eggnog.db.gz`、`e5.proteomes.faa`、`e5.taxid_info.tsv`、`e5.proteomes.fungi.faa`）隔离在 `<runtime>/var/downloads/eggnog/` 中，保证 `tools/eggnog_mapper/data/` 符合 `--profile paper` 的整洁要求；
  - `--with-prott5`：通过 `huggingface-cli` / `hf` 下载 `Rostlab/prot_t5_xl_uniref50`（仅 `alphagem_mode: "structure_alignment"` 需要）。
- **`annotation`**（`--tool-root /work/annotation-tools`，可选 `--with-helixer`）：获取 ANNEVO 源码 + `ANNEVO_Fungi.pt`（需要单独的 PyTorch 环境）。

### 从头重建（步骤 1–3，`--profile paper`）的严格边界

`bash run.sh setup`（uv）只安装轻量的 Geno2Flux 核心 + HiGHS 求解器。从头运行 AlphaGEM 步骤 1–3 **严格要求 Linux x86_64** 以及两个独立 Python 环境：

1. **AlphaGEM 环境**（配置 `"alphagem_python"` 或 `G2F_ALPHAGEM_PYTHON`）：由 [`third_party/alphagem/AlphaGEM.yml`](third_party/alphagem/AlphaGEM.yml) 创建（Linux x86_64；锁定 Python 3.10.15、COBRApy 0.29.1、BLAST+ 2.16.0、DIAMOND 2.1.10、Foldseek `9.427df8a`），**另需**安装 PyTorch、`fair-esm==2.0.0`，并按 [`third_party/alphagem/setup.sh`](third_party/alphagem/setup.sh) 从隔离副本安装 CLEAN——不要在 `tools/CLEAN/app` 内原地安装，也不要在 `tools/` 下残留 `__pycache__`/`.pyc`。
2. **DeepEC 环境**（配置 `"deepec_python"`）：由 [`third_party/alphagem/deepecenvironment.yml`](third_party/alphagem/deepecenvironment.yml) 创建（Linux x86_64；锁定 Python 3.7.7、`torch==1.7.0`、`transformers==3.5.1`）。
3. **`fungi.dmnd` 契约门**：统一流程调用 `AlphaGEM.py --profile paper`，要求 `fungi.dmnd` 精确字节长度 `879,076,733` 且 SHA-256 与 Linux x86_64 的 DIAMOND 2.1.10 构建结果匹配。在非 Linux 平台或用其他 DIAMOND 版本构建会导致 `fetch-resources --with-eggnog` 以 `contract_mismatch`（退出码 2）失败，`check-env` 报告 `BLOCKED`。`check-env --profile alphagem` 还会校验 `eggnog.db`（字节数 + SHA-256）以及 `tools/eggnog_mapper/data/` 下 9 个文件的 `eggnog_static_data` 契约——该目录必须恰好包含这 9 个文件加 `eggnog.db` 与 `fungi.dmnd`，中间产物（`eggnog.db.gz`、`e5.proteomes.faa`、`e5.taxid_info.tsv`、`e5.proteomes.fungi.faa`）会被判定为违规文件。
4. **`alphagem_mode: "structure_alignment"`**：还需在配置中提供 `prottrans_model_path`、`structure_maplist` 与 `structure_dir`（含 `*.pdb`），并依赖 Foldseek `9.427df8a` 与 `data_source.tar.xz` 内附带的 Linux x86_64 `USalign` 二进制。`structure_maplist` 为 Excel 工作簿（`.xlsx`），须包含 `Entry`、`Entry Name` 以及一个受支持的结构列（如 `Structure Data` / `AlphaFoldDB` / `PDB`），用于把每个 UniProt Entry 映射到 `structure_dir` 中的结构文件。（默认 `plmsearch` 模式使用 [PLMSearch](https://github.com/maovshao/PLMSearch)。）
5. 用以下命令验证环境交接：

   ```sh
   bash run.sh check-env --profile alphagem --runtime-root /work/alphagem-runtime \
     --alphagem-python /path/to/AlphaGEM/bin/python \
     --deepec-python /path/to/deepectransformer/bin/python
   ```

6. **运行目录归属规则**：流水线托管的运行目录（`"alphagem_runtime_root": null`）下，将 `reference_fasta` 与 `reference_xlsx` 指向你的 `<refname>.fasta` / `<refname>.xlsx` 文件，流水线会自行复制进 `data_available/`；若是另行准备的运行目录，请先把两个文件放入 `<runtime>/data_available/`（该目录按只读对待）。

### 上游官方来源

- 基因组注释：[ANNEVO](https://github.com/xjtu-omics/ANNEVO)、[Helixer](https://github.com/usadellab/Helixer)、[Tiberius](https://github.com/Gaius-Augustus/Tiberius)、[funannotate](https://github.com/nextgenusfs/funannotate)
- AlphaGEM 重建：[AlphaGEM Figshare 资源](https://doi.org/10.6084/m9.figshare.30630524)、[ESM-1b](https://github.com/facebookresearch/esm)、[PLMSearch](https://github.com/maovshao/PLMSearch)、[CLEAN](https://github.com/tttianhao/CLEAN)、[DeepECtransformer](https://github.com/kaistsystemsbiology/DeepProZyme)、[eggNOG-mapper](https://github.com/eggnogdb/eggnog-mapper)
- 比对与求解：[NCBI BLAST+ 2.16.0](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.16.0/)、[DIAMOND](https://github.com/bbuchfink/diamond)、[HiGHS](https://highs.dev/)、[uv](https://docs.astral.sh/uv/)

仓库自带 AlphaGEM 源码快照（GPL-3.0）、Yeast9 模板（CC BY 4.0）与参考反应池；**不包含**重建数据库与训练权重——由 `bash run.sh fetch-resources` 按需拉取，`prepare_alphagem_runtime.py` 只暂存已校验的源码文件，不下载这些资源。

---

## 6. 运行流程

### 第 1 步 —— 由基因组预测蛋白质组（或标准化已有编码 GFF3）

可选：为同源引导的预测工具（`funannotate` / `miniprot`）下载 NCBI RefSeq 或 UniProt 参考蛋白质组：

```sh
bash run.sh fetch-reference-proteomes \
  --uniprot-proteome YEAST=UP000002311 \
  --output /work/ref_proteomes
```

随后用 `predict-proteome` 一条命令完成基因预测（`annevo` / `helixer` / `tiberius` / `funannotate` / `custom`，或直接读入已有 `gff`）、最长转录本蛋白质组提取与质控，产出可直接使用的 `/work/pred/inputs/` AlphaGEM 靶点与坐标契约（也可用 `bash run.sh prepare-inputs` 单独由基因组 + GFF3 + FASTA 构建坐标契约）：

```sh
# 从已有编码 GFF3 出发（即使未安装 gffread 也能工作）：
bash run.sh predict-proteome \
  --tool gff \
  --genome /work/raw/genome.fasta \
  --gff /work/raw/prediction.gff3 \
  --case-id sample_strain \
  --prefix PROT \
  --output /work/pred

# 或直接调用 ANNEVO / Helixer / Tiberius / funannotate：
bash run.sh predict-proteome \
  --tool annevo \
  --genome /work/raw/genome.fasta \
  --annevo-dir /path/to/ANNEVO \
  --model-weights /path/to/ANNEVO_Fungi.pt \
  --tool-python /path/to/annevo-env/bin/python \
  --case-id sample_strain \
  --output /work/pred
```

### 第 2 步 —— 运行统一流水线（或逐阶段运行）

复制 [`configs/unified_pipeline.example.json`](configs/unified_pipeline.example.json) 为 `/work/run.local.json`，将其中的蛋白质组/坐标字段指向 `/work/pred/inputs/`。若进行 AlphaGEM 从头重建，把 `fetch-resources --profile alphagem` 打印的 `config_snippet` 合并进 `/work/run.local.json`——它会填好 `alphagem_runtime_root`、`reference_fasta`、`reference_xlsx` 与 `esm_model_path`，并给出 `alphagem_python` 与 `deepec_python` 的占位文本，合并后必须把这两个字段替换为你新建的 AlphaGEM 与 DeepEC 环境的实际 Python 路径（默认 `plmsearch` 模式保持 `prottrans_model_path: null`）。探索性运行时请记录 SHA-256 哈希或设置 `"allow_unverified_reuse": true`，然后执行：

```sh
bash run.sh metabolic --config /work/run.local.json --output /work/run-001
```

`/work/run-001/<order>/` 下的关键产物：

- `localization_refinement_audit.json` —— 被拉动的分室反应记录（或默认分室降级回执）。
- `target_model_<order_slug>.xml`（例如 `pre-gapfill/target_model_pre_gapfill.xml`，或 `skip` 模式下的 `screening/target_model_screening.xml`）与 `target_module_receipt.json` —— 已挂载靶点模块的目标模型（新增内部反应默认通过元素与电荷平衡校验）。
- `fseof/all_genes_fseof_essential_filtered.tsv` —— 按 10 步 FSEOF 通量斜率排序的非必需过表达候选。
- `paralog/paralog_family_rankings_essential_filtered.tsv` —— 经 DIAMOND 自比对聚类、按聚合 FSEOF 得分排序的同源基因家族。
- `knockout/single_gene_ko_essential_filtered.tsv` 与 `knockout/family_ko_essential_filtered.tsv` —— 排除生长必需基因后的单基因与同源家族联合敲除排序。

也可以先准备独立的 AlphaGEM 运行目录（`uv run --frozen python scripts/prepare_alphagem_runtime.py --prepare --runtime-root /work/alphagem-runtime`），再用 `bash run.sh stage <localize|gapfill|target|fseof|paralogs|knockout|filter>` 逐阶段执行：

```sh
export G2F_RUNTIME=/work/alphagem-runtime
export G2F_STAGE=/work/stages-001

# 在 AlphaGEM 草图重建（--name sample_strain --stop-after reaction_pool）之后：
bash run.sh stage localize --name sample_strain --inputs /work/pred/inputs
bash run.sh stage gapfill  --name sample_strain --inputs /work/pred/inputs
bash run.sh stage target   --inputs /work/pred/inputs --target-reaction r_2056
bash run.sh stage fseof    --inputs /work/pred/inputs --growth-reaction r_4041
bash run.sh stage paralogs --inputs /work/pred/inputs
bash run.sh stage knockout --inputs /work/pred/inputs
bash run.sh stage filter   --inputs /work/pred/inputs --growth-reaction r_4041
```

---

## 7. 仓库结构

- `src/geno2flux/` —— 统一流水线、亚细胞分室细化（DeepLoc 自动检测与默认降级）、可配置靶点模块、RAVEN 兼容 FSEOF 移植、DIAMOND 旁系同源分组、敲除筛选、生长必需性过滤与原生 HiGHS 后端。
- `scripts/` —— CLI 工具：基因组→蛋白质组预测（`predict_proteome.py`）、RefSeq/UniProt 参考蛋白质组拉取（`fetch_reference_proteomes.py`）、外部资源下载与就绪预检（`fetch_external_resources.py`）、GFF3 标准化（`standardize_annotation.py`）、坐标契约构建（`build_coordinate_contract.py`）、AlphaGEM 运行目录暂存（`prepare_alphagem_runtime.py`）与逐阶段执行（`run_stage.py`）。
- `tools/deeploc/` —— DeepLoc 2.0/2.1 的即插即用槽位（gitignore 排除，`README.md` 除外）。
- `tools/bin/` —— 下载的外部二进制（DIAMOND、BLAST+）投放目录，`run.sh` 自动前置到 `PATH`（gitignore 排除，`README.md` 除外）。
- `configs/unified_pipeline.example.json` —— 完整基因组→靶点运行的带注释模板配置。
- `third_party/alphagem/` —— 捆绑的 GPL-3.0 AlphaGEM 源码快照，含 HiGHS 求解器适配（`optlang_highs_interface.py`）、坐标契约与 CC BY 4.0 `yeast-GEM.xml` 模板及 `models/reaction_pool_final.json`。

---

## 8. 许可证与第三方声明

- Geno2Flux 以 GPL-3.0 许可证发布——见 [LICENSE](LICENSE)。
- `third_party/alphagem/` 捆绑 GPL-3.0 的 AlphaGEM 源码快照；Yeast9 模板采用 CC BY 4.0 许可；DeepLoc 2.0/2.1 仍属 DTU 个人学术许可，仓库不随附再分发。
- 所有内置组件与外部拉取组件的完整署名见 [THIRD_PARTY.md](THIRD_PARTY.md)；Python 依赖清单见 [requirements.txt](requirements.txt)。
