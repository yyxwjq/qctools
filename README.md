# qctools# QCTools - Quantum Chemistry Analysis Toolkit

[English](#english) | [中文](#中文)

---

## English

### Overview

QCTools is a comprehensive Python toolkit for quantum chemistry calculations and analysis, designed for efficient processing of atomic structures and molecular dynamics trajectories. The package provides powerful tools for structural analysis, machine learning potential evaluation, and automated error detection in computational chemistry workflows.

### Features

- **📊 Radial Distribution Function (RDF)**: Calculate and visualize radial distribution functions with parallel processing support
- **📐 Angular Distribution Function (ADF)**: Compute angular distribution functions for structural analysis
- **🔢 Coordination Number**: Calculate coordination numbers between atom groups with efficient parallel processing
- **🔬 ML Potential Error Analysis**: Comprehensive error analysis for machine learning potentials (NEP, n2p2)
- **⚛️ Atomic Structure Editing**: Tools for manipulating atomic structures and removing specific elements
- **📏 Radius Analysis**: Advanced algorithms for minimum enclosing ball calculations
- **🔧 Logging System**: Integrated logging system for all analysis workflows

### Installation

#### Development Installation (Recommended)

```bash
# Clone the repository
git clone <repository-url>
cd qctools

# Activate conda environment (if using conda)
conda activate wizard

# Install in development mode
pip install -e .
```

#### Regular Installation

```bash
pip install qctools
```

### Quick Start

```python
import qctools
from qctools import rdf, adf, coord, ml
from ase.io import read

# Set up logging
qctools.qctools_logging()

# Load trajectory
images = read('trajectory.xyz', ':')

# Calculate RDF
rdf.get_rdf(images, cutoff=5.0, bin_size=0.1, cores=4)

# Calculate coordination numbers
coord_nums = coord.group_coordnum(images, group1=[0,1,2], group2=[3,4,5])

# ML potential error analysis
ml.error_img.main(
    trajname='trajectory.xyz',
    apps='nep',
    resource='software',
    fontsize=12,
    data={'energy': 'energy.txt', 'force': 'force.txt'}
)
```

### Module Documentation

#### 🔍 RDF Analysis (`qctools.rdf`)

Calculate radial distribution functions with support for first-neighbor analysis:

```python
from qctools.rdf import get_rdf

# Basic RDF calculation
get_rdf(images, cutoff=5.0, bin_size=0.1, first_neighbor=False, cores=4)

# First neighbor only
get_rdf(images, cutoff=5.0, bin_size=0.1, first_neighbor=True, cores=4)
```

#### 📐 ADF Analysis (`qctools.adf`)

Compute angular distribution functions:

```python
from qctools.adf import get_adf

get_adf(images, rcut=4.0, bin_size=5.0, cores=4)
```

#### 🔢 Coordination Number (`qctools.coord`)

Efficient coordination number calculations:

```python
from qctools.coord import group_coordnum

# Calculate coordination numbers between two groups
coord_nums = group_coordnum(
    traj=images,
    group1=[0, 1, 2],  # Indices of first group
    group2=[3, 4, 5],  # Indices of second group
    r0=4.0,            # Cutoff distance
    cores=4            # Number of cores
)
```

#### 🤖 ML Error Analysis (`qctools.ml.error_img`)

Comprehensive error analysis for machine learning potentials:

```python
from qctools.ml.error_img import main

# Analyze NEP potential errors
main(
    trajname='trajectory.xyz',
    apps='nep',                    # or 'n2p2'
    resource='software',           # or 'images'
    fontsize=12,
    data={'energy': 'energy.txt', 'force': 'force.txt'},
    er_bar=1.5,                   # Error threshold multiplier
    show_marginals=True           # Show marginal distributions
)
```

#### ⚛️ Structure Editing (`qctools.edit_atoms`)

Remove specific elements from structures:

```python
from qctools.edit_atoms import remove

# Remove hydrogen atoms
remove('structure.vasp', ['H'])
```

### Dependencies

- Python ≥ 3.7
- ASE (Atomic Simulation Environment)
- NumPy
- Matplotlib
- SciPy
- PyNEP (for NEP potential analysis)

### Contributing

We welcome contributions! Please feel free to submit pull requests or open issues for bug reports and feature requests.

### License

This project is licensed under the MIT License.

---

## 中文

### 概述

QCTools 是一个全面的量子化学计算分析Python工具包，专为高效处理原子结构和分子动力学轨迹而设计。该软件包为结构分析、机器学习势能评估和计算化学工作流程中的自动错误检测提供了强大的工具。

### 功能特性

- **📊 径向分布函数 (RDF)**: 计算和可视化径向分布函数，支持并行处理
- **📐 角度分布函数 (ADF)**: 计算用于结构分析的角度分布函数
- **🔢 配位数计算**: 高效并行计算原子团簇间的配位数
- **🔬 机器学习势能误差分析**: 全面的机器学习势能误差分析 (NEP, n2p2)
- **⚛️ 原子结构编辑**: 操作原子结构和移除特定元素的工具
- **📏 半径分析**: 最小包围球计算的高级算法
- **🔧 日志系统**: 集成的日志系统，用于所有分析工作流程

### 安装方法

#### 开发模式安装（推荐）

```bash
# 克隆仓库
git clone <repository-url>
cd qctools

# 激活conda环境（如果使用conda）
conda activate wizard

# 开发模式安装
pip install -e .
```

#### 常规安装

```bash
pip install qctools
```

### 快速开始

```python
import qctools
from qctools import rdf, adf, coord, ml
from ase.io import read

# 设置日志系统
qctools.qctools_logging()

# 加载轨迹文件
images = read('trajectory.xyz', ':')

# 计算径向分布函数
rdf.get_rdf(images, cutoff=5.0, bin_size=0.1, cores=4)

# 计算配位数
coord_nums = coord.group_coordnum(images, group1=[0,1,2], group2=[3,4,5])

# 机器学习势能误差分析
ml.error_img.main(
    trajname='trajectory.xyz',
    apps='nep',
    resource='software',
    fontsize=12,
    data={'energy': 'energy.txt', 'force': 'force.txt'}
)
```

### 模块文档

#### 🔍 RDF分析 (`qctools.rdf`)

计算径向分布函数，支持第一近邻分析：

```python
from qctools.rdf import get_rdf

# 基本RDF计算
get_rdf(images, cutoff=5.0, bin_size=0.1, first_neighbor=False, cores=4)

# 仅第一近邻
get_rdf(images, cutoff=5.0, bin_size=0.1, first_neighbor=True, cores=4)
```

#### 📐 ADF分析 (`qctools.adf`)

计算角度分布函数：

```python
from qctools.adf import get_adf

get_adf(images, rcut=4.0, bin_size=5.0, cores=4)
```

#### 🔢 配位数计算 (`qctools.coord`)

高效的配位数计算：

```python
from qctools.coord import group_coordnum

# 计算两组原子间的配位数
coord_nums = group_coordnum(
    traj=images,
    group1=[0, 1, 2],  # 第一组原子索引
    group2=[3, 4, 5],  # 第二组原子索引
    r0=4.0,            # 截止距离
    cores=4            # 核心数
)
```

#### 🤖 机器学习误差分析 (`qctools.ml.error_img`)

机器学习势能的全面误差分析：

```python
from qctools.ml.error_img import main

# 分析NEP势能误差
main(
    trajname='trajectory.xyz',
    apps='nep',                    # 或 'n2p2'
    resource='software',           # 或 'images'
    fontsize=12,
    data={'energy': 'energy.txt', 'force': 'force.txt'},
    er_bar=1.5,                   # 误差阈值乘数
    show_marginals=True           # 显示边际分布
)
```

#### ⚛️ 结构编辑 (`qctools.edit_atoms`)

从结构中移除特定元素：

```python
from qctools.edit_atoms import remove

# 移除氢原子
remove('structure.vasp', ['H'])
```

### 依赖项

- Python ≥ 3.7
- ASE (原子模拟环境)
- NumPy
- Matplotlib
- SciPy
- PyNEP (用于NEP势能分析)

### 贡献

我们欢迎贡献！请随时提交pull request或为错误报告和功能请求创建issue。

### 许可证

本项目采用MIT许可证。

---

### Examples / 示例

More examples can be found in the `example/` directory:
- `example/rdf/test.py` - RDF calculation example
- `example/adf/test.py` - ADF calculation example  
- `example/coord/test_group_coordnum.py` - Coordination number example
- `example/ml/error_img/` - ML error analysis examples

更多示例可在 `example/` 目录中找到：
- `example/rdf/test.py` - RDF计算示例
- `example/adf/test.py` - ADF计算示例  
- `example/coord/test_group_coordnum.py` - 配位数计算示例
- `example/ml/error_img/` - 机器学习误差分析示例

### Contact / 联系方式

For questions and support, please open an issue on GitHub.

如有问题和支持需求，请在GitHub上创建issue。
