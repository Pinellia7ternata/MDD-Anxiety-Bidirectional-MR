# MDD-Anxiety MR分析 - 补充分析脚本说明

## 概述
本目录包含5个补充分析脚本，用于提升MDD-Anxiety双向MR分析的文章档次。

## 脚本清单

### 🔴 高优先级

#### 1. `10_run_cause_bayesian.R` - 完整贝叶斯CAUSE建模
**用途**：区分遗传相关与真实因果效应，提取causal/sharing posterior

**前置依赖**：
- R包：`causalMendelian`（脚本会自动安装）
- 数据：MDD/Anxiety formatted GWAS

**运行命令**：
```bash
Rscript scripts/10_run_cause_bayesian.R
```

**输出文件**：
- `10_CAUSE_bayesian/CAUSE_posterior_MDD_to_Anxiety.csv`
- `10_CAUSE_bayesian/CAUSE_posterior_Anxiety_to_MDD.csv`
- `10_CAUSE_bayesian/Figure5_CAUSE_MDD_to_Anxiety.png`
- `10_CAUSE_bayesian/Figure5_CAUSE_Anxiety_to_MDD.png`

**预计运行时间**：2-4小时（取决于SNP数量）

---

#### 2. `11_run_magma_gene_set.R` - MAGMA基因集富集分析
**用途**：功能注释，识别MDD/Anxiety共享的生物学通路

**前置依赖**：
- MAGMA软件（需单独下载）
  - 下载地址：https://ctg.cncr.nl/software/magma
  - 解压到：`D:/tools/magma/`
- 1000G EUR LD参考面板（g1000_eur.zip，约300MB）
  - 从MAGMA官网下载
- 基因位置文件（NCBI37.3.gene.loc）
  - 从MAGMA官网下载

**运行命令**：
```bash
Rscript scripts/11_run_magma_gene_set.R
```

**输出文件**：
- `11_MAGMA/MDD_gene_results.csv`
- `11_MAGMA/Anxiety_gene_results.csv`
- `11_MAGMA/shared_significant_genes.csv`
- `11_MAGMA/Figure6_MAGMA_MDD_manhattan.png`

**预计运行时间**：1-2小时

**⚠️ 重要提示**：如果不想安装MAGMA，可使用FUMA在线平台（见脚本13）

---

### 🟡 中优先级

#### 3. `12_run_gtex_eqtl_coloc.R` - GTEx brain eQTL共定位
**用途**：验证共享位点的脑组织表达调控机制

**前置依赖**：
- GTEx v8 eQTL数据（约50GB）
  - 下载地址：https://gtexportal.org/home/dataset/v8
  - 或使用Globus高速下载
  - 解压到：`D:/data/GTEx_v8_eQTL/`
- R包：`coloc`（已安装）

**运行命令**：
```bash
Rscript scripts/12_run_gtex_eqtl_coloc.R
```

**输出文件**：
- `12_GTEx_brain_eqtl_coloc/eQTL_coloc_summary.csv`
- 每个locus的详细coloc结果

**预计运行时间**：4-8小时（取决于locus数量和组织数量）

**⚠️ 备选方案**：可使用FUMA平台的GTEx eQTL注释功能，无需下载数据

---

#### 4. `13_prepare_fuma_input.R` - FUMA在线平台输入文件准备
**用途**：生成FUMA平台的输入文件，实现一站式功能注释

**优势**：
- 无需本地安装软件
- 自动完成gene mapping、eQTL、染色质交互、MAGMA、通路富集
- 可视化结果直接可用

**运行命令**：
```bash
Rscript scripts/13_prepare_fuma_input.R
```

**输出文件**：
- `13_FUMA_input/FUMA_input_MDD.txt`
- `13_FUMA_input/FUMA_input_Anxiety.txt`
- `13_FUMA_input/lead_SNPs_for_FUMA.txt`

**后续步骤**：
1. 访问 https://fuma.ctglab.nl/
2. 上传生成的文件
3. 等待2-24小时获取结果
4. 下载完整注释报告

**推荐指数**：⭐⭐⭐⭐⭐（最省事，功能最全）

---

### 🟢 低优先级

#### 5. `14_run_susie_finemapping.R` - SuSiE精细定位
**用途**：识别每个locus的95% credible set，精确定位因果变异

**前置依赖**：
- R包：`susieR`（已安装）
- 理想情况：真实LD矩阵（可从1000G计算）
- 简化情况：脚本使用Z-score近似

**运行命令**：
```bash
Rscript scripts/14_run_susie_finemapping.R
```

**输出文件**：
- `14_SuSiE_finemapping/SuSiE_MDD_locus*.csv`
- `14_SuSiE_finemapping/SuSiE_Anxiety_locus*.csv`
- `14_SuSiE_finemapping/SuSiE_summary_all_loci.csv`
- `14_SuSiE_finemapping/Figure7_SuSiE_PIP_comparison.png`

**预计运行时间**：30分钟-1小时

---

## 推荐运行顺序

### 方案A：本地完整分析
1. **FUMA输入文件准备**（脚本13）→ 上传FUMA → 等待结果
2. **CAUSE建模**（脚本10）→ 提取posterior
3. **SuSiE精细定位**（脚本14）→ credible set
4. 如FUMA无MAGMA → **本地MAGMA**（脚本11）
5. 如需组织特异性 → **GTEx eQTL**（脚本12）

### 方案B：快速补充（推荐）
1. **FUMA输入文件准备**（脚本13）→ 上传FUMA
2. **CAUSE建模**（脚本10）
3. **SuSiE精细定位**（脚本14）
4. FUMA结果返回后即可投稿

### 方案C：最小补充（时间紧迫时）
1. **FUMA输入文件准备**（脚本13）→ 上传FUMA
2. 使用FUMA结果 + 当前V9.1文档投稿

---

## 注意事项

1. **CAUSE运行时间较长**：全基因组约百万SNP，建议使用`subsample`参数
2. **MAGMA需要LD参考面板**：约300MB，从官网下载
3. **GTEx数据量较大**：完整eQTL约50GB，可只下载脑组织相关
4. **FUMA最省事**：无需本地安装，功能最全，强烈推荐
5. **真实LD矩阵**：所有基于LD的分析（MAGMA、SuSiE）使用真实LD更准确

---

## 脚本位置

所有脚本位于：`D:\2026年\文章\因果推断\scripts\`

---

## 联系

如有问题，请检查：
1. 文件路径是否正确
2. R包是否已安装
3. 外部软件（MAGMA）是否已下载
4. 数据文件是否存在

生成时间：2026-05-29 08:52
版本：V9.1补充脚本
