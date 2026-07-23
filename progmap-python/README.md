# ProgMap Python Tool（Linux/服务器版）

ProgMap 将正常组织、Stage I 和 Stage II/III 的表达谱与甲基化谱整合为 MECor 特征，并使用固定的 2048–128 跳连神经网络完成三分类和进展相关特征识别。该版本面向无图形界面的 Linux 工作站、计算服务器和 Slurm 集群，可通过一条 `progmap` 命令运行。

English documentation: [README_EN.md](README_EN.md)

## 1. 参考运行环境

- 操作系统：Ubuntu 22.04/24.04 x86_64 或兼容 Linux 发行版
- Python：3.9–3.11，推荐 3.11
- TensorFlow/Keras：2.14.0
- NumPy：1.26.4
- pandas：2.2.3
- SciPy：1.13.1
- scikit-learn：1.5.2
- statsmodels：0.14.4
- joblib：1.4.2

旧项目模型文件记录了 TensorFlow/Keras 2.14.0，旧缓存文件记录了 Python 3.9。为兼顾旧模型兼容性和服务器可重建性，本发布版固定 TensorFlow 2.14.0，并正式支持 Python 3.9–3.11。建议至少使用 16 GB 内存；约 13,000 个输入基因时，2048 单元的第一层参数量较大，32 GB 内存更稳妥。

## 2. Linux安装

### venv + pip（推荐）

```bash
git clone https://github.com/MengyanZhang-bioinfo/ProgMap.git
cd ProgMap/progmap-python
python3.11 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements-linux-cpu.txt
python -m pip install --no-deps .
progmap --help
```

`--no-deps` 是因为依赖已由锁定文件安装。如果希望由包管理器自动解析依赖，也可直接执行 `python -m pip install .`。

### Conda

```bash
conda env create -f environment-linux.yml
conda activate progmap-linux
progmap --help
```

## 3. 输入目录

`--data-root` 下每个癌种一个文件夹；名称为 `GEO` 的文件夹自动排除。每个癌种必须包含六个“基因 × 样本”CSV矩阵，第一列是基因标识，第一行是样本标识。

| 分组 | 表达谱 | 甲基化谱 |
|---|---|---|
| Normal | `en.csv` | `mn.csv` |
| Stage I | `exp1.csv` 或 `e1.csv` | `me1.csv` 或 `m1.csv` |
| Stage II/III | `exp2.csv` 或 `e2.csv` | `me2.csv` 或 `m2.csv` |

Linux文件名区分大小写，请保持上述小写文件名。表达与甲基化样本按样本名配对，六个矩阵按共有基因对齐。

示例结构：

```text
/data/PANCANCER/
├── BRCA/
│   ├── en.csv
│   ├── e1.csv
│   ├── e2.csv
│   ├── mn.csv
│   ├── m1.csv
│   └── m2.csv
├── CESC/
│   └── ...
└── GEO/                 # 自动排除
```

## 4. 一条命令运行

先检查输入，不启动TensorFlow训练：

```bash
progmap \
  --data-root /data/PANCANCER \
  --output /results/progmap_input_check \
  --cancers all \
  --dry-run \
  --device cpu
```

对全部癌种运行原始单侧Welch t检验并输出全部基因：

```bash
progmap \
  --data-root /data/PANCANCER \
  --output /results/progmap_ttest \
  --cancers all \
  --test ttest \
  --top-n all \
  --device auto \
  --threads 8
```

使用Wilcoxon rank-sum / Mann–Whitney U检验并输出前100个基因：

```bash
progmap \
  --data-root /data/PANCANCER \
  --output /results/progmap_wilcoxon \
  --cancers all \
  --test wilcoxon \
  --top-n 100 \
  --device auto \
  --threads 8
```

使用置换检验和特征效应量bootstrap区间：

```bash
progmap \
  --data-root /data/PANCANCER \
  --output /results/progmap_permutation \
  --cancers CESC \
  --test permutation \
  --permutations 10000 \
  --bootstrap-iterations 2000 \
  --top-n 100 \
  --device auto \
  --threads 8
```

查看全部参数：

```bash
progmap --help
```

## 5. 后台与集群运行

普通Linux服务器可用：

```bash
nohup progmap \
  --data-root /data/PANCANCER \
  --output /results/progmap \
  --cancers all \
  --test ttest \
  --top-n all \
  --device auto \
  --threads 8 \
  > progmap.log 2>&1 &
```

记录进程号并查看日志：

```bash
echo $!
tail -f progmap.log
```

Shell封装脚本：

```bash
chmod +x examples/run_all.sh
PROGMAP_THREADS=8 PROGMAP_DEVICE=auto \
  ./examples/run_all.sh /data/PANCANCER /results/progmap
```

Slurm示例：

```bash
sbatch server/slurm_run.sh /data/PANCANCER /results/progmap
```

提交前根据服务器策略调整 `server/slurm_run.sh` 中的内存、时间和CPU数量。如果需要GPU，请在集群脚本中增加相应的GPU资源请求，并使用 `--device gpu`；程序在TensorFlow无法识别GPU时会立即报错，不会悄悄退回CPU。

## 6. Docker

CPU镜像：

```bash
docker build -t progmap:0.2.0 .
docker run --rm \
  -v /absolute/path/PANCANCER:/data/PANCANCER:ro \
  -v /absolute/path/results:/results \
  progmap:0.2.0 \
  --data-root /data/PANCANCER \
  --output /results/run1 \
  --cancers all \
  --device cpu \
  --threads 8
```

GPU镜像使用TensorFlow 2.14.0官方GPU基础镜像；宿主机需安装兼容的NVIDIA驱动和NVIDIA Container Toolkit：

```bash
docker build -f Dockerfile.gpu -t progmap:0.2.0-gpu .
docker run --rm --gpus all \
  -v /absolute/path/PANCANCER:/data/PANCANCER:ro \
  -v /absolute/path/results:/results \
  progmap:0.2.0-gpu \
  --data-root /data/PANCANCER \
  --output /results/run1 \
  --cancers all \
  --device gpu
```

## 7. 固定模型和默认参数

- MECor输入：`expression_scaled × methylation_scaled × training_fold_Pearson_correlation`。
- 网络：Dense(2048, ReLU) → BatchNorm → Dropout(0.1) → Dense(128, ReLU) → BatchNorm → Dropout(0.1) → 与原MECor输入拼接 → Dense(3, Softmax)。
- 两个隐藏层均使用max-norm 3约束。
- 优化器：Adam；初始学习率 `1e-3`，每10,000步按0.9指数衰减。
- batch size：16；最大epoch：1,000。
- 早停：监测外层训练折内部验证集的 `val_auc`，patience 200，并恢复最佳权重。
- 默认类别权重：Normal/Stage I/Stage II–III = 0.25/0.50/0.25。
- 解释方法：多训练样本中位基线的增强积分梯度，默认64步、3个基线。
- 默认统计检验：单侧Welch t检验；默认Benjamini–Hochberg FDR，`alpha=0.01`。
- `--top-n all` 输出全部基因；例如 `--top-n 100` 输出最显著的100个不重复基因。

## 8. 数据泄漏控制

每个外层交叉验证折独立执行：

1. 划分外层训练集与外层测试集。
2. 仅在外层训练集内部划分拟合集和早停验证集。
3. 只用拟合集计算插补值、表达和甲基化0–1缩放参数以及逐基因Pearson相关系数。
4. 固定上述参数后再变换内部验证集和外层测试集。
5. 外层测试集仅用于该折最终预测和解释，不参与预处理拟合、相关计算、模型训练或早停。

每折相关系数写入 `fold_N/training_fold_correlations.csv`，每个样本的角色写入 `fold_N/sample_roles.csv`。`run_config.json` 同时记录Python、平台、依赖版本、线程设置和TensorFlow实际识别的设备。

## 9. 输出

- `cross_validated_predictions.csv`：每个外层测试样本的OOF预测概率。
- `fold_metrics.csv`：每折accuracy、balanced accuracy和多分类AUC。
- `features_by_class_all.csv`：三个目标类别的全部检验结果。
- `features_ranked_genes.csv`：每个基因保留最显著类别后的排名。
- `features_selected.csv`：由 `--top-n` 控制的最终特征表。
- `summary.json`、`run_config.json`：配置、样本量、版本、设备和汇总指标。
- `fold_N/preprocessor.joblib`、`model.keras`：每折预处理器和模型。

## 10. 测试

```bash
python -m pytest
python -m build
```

生成最小合成数据并检查输入：

```bash
python examples/create_synthetic_data.py --output /tmp/progmap/PANCANCER
progmap \
  --data-root /tmp/progmap/PANCANCER \
  --output /tmp/progmap/input_check \
  --dry-run \
  --device cpu
```

GitHub Actions配置位于 `.github/workflows/linux-ci.yml`，在Ubuntu 22.04上执行安装、单元测试、端到端训练测试、命令行冒烟测试和wheel/sdist构建。

## 11. 可重复性说明

相同随机种子提高可重复性，但不同CPU指令集、GPU型号、CUDA/cuDNN和TensorFlow算子仍可能造成微小浮点差异。正式分析应保存完整输出目录、`run_config.json`、输入文件校验值和Git提交号。数据不随软件包分发。

