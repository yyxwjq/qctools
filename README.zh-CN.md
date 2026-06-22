# QCTools - 量子化学分析工具包

[English](README.md)

QCTools 是一个面向日常量子化学和原子模拟分析的小型 Python 工具包。当前功能集中在轨迹结构分析、机器学习势误差检查、简单结构编辑，以及 VASP 任务目录巡检。

当前版本：`0.1.0`

## 安装

从 GitHub 以开发模式安装：

```bash
git clone https://github.com/yyxwjq/qctools.git
cd qctools
pip install -e .
```

也可以直接从 GitHub 安装：

```bash
pip install git+https://github.com/yyxwjq/qctools.git
```

NEP 相关功能需要额外安装 `pynep`：

```bash
pip install -e .[nep]
```

## Python 工具

### 日志

```python
import qctools

qctools.qctools_logging(filename="qctools.log", overwrite=False)
```

`overwrite=True` 会在配置日志前删除已有日志文件。

### RDF

```python
from ase.io import read
from qctools.rdf import get_rdf

images = read("trajectory.xyz", ":")
get_rdf(images, cutoff=5.0, bin_size=0.1, first_neighbor=False, cores=4)
```

RDF 输出到 `RDF/` 或 `RDF_first/`。当前实现使用按帧平均、按密度归一化的 `g(r)`。

### ADF

```python
from ase.io import read
from qctools.adf import get_adf

images = read("trajectory.xyz", ":")
get_adf(images, rcut=4.0, bin_size=5.0, cores=4)
```

ADF 输出到 `ADF/`。`rcut` 按中心原子到邻居原子的真实距离筛选。

### 配位数

```python
from qctools.coord import group_coordnum

coord_nums = group_coordnum(
    traj=images,
    group1=[0, 1, 2],
    group2=[3, 4, 5],
    r0=4.0,
    cores=4,
)
```

### 机器学习势误差分析

```python
from qctools.ml.error_img import main

main(
    trajname="trajectory.xyz",
    apps="nep",          # "nep" 或 "n2p2"
    resource="software", # "software" 或 "images"
    fontsize=12,
    data={"energy": "energy.txt", "force": "force.txt"},
    er_bar=1.5,
    show_marginals=True,
)
```

主要输出文件包括：

- `energy.data`, `force.data`
- `energy_error.txt`, `force_error.txt`
- `Err-energy.xyz`
- `Err-force-ini.xyz`, `Err-force-replaced.xyz`
- `leave-E-img.xyz`, `leave-F-img.xyz`
- `energy_error_analysis.png`, `force_error_analysis.png`

### 结构编辑

```python
from qctools.edit_atoms import remove

remove("structure.vasp", ["H"])
```

命令行入口：

```bash
qctools-edit remove structure.vasp H
```

## VASP 任务巡检

`qctools/vaspjob` 是一个 Bash 工具，用于扫描 VASP 计算目录。它会查找 `INCAR`，读取相关的 `POSCAR`、`OSZICAR`、`OUTCAR`，并把汇总表写入 `vaspjob.log`。

在包含 VASP 任务目录的路径下运行：

```bash
qctools/vaspjob
qctools/vaspjob -v
```

参数：

- `-v`：同时把报告打印到终端。
- `-m`：生成 `wavecar_list.txt`，列出所有 `WAVECAR` 文件及大小。
- `-r`：删除 `wavecar_list.txt` 中列出的 `WAVECAR` 文件。

`-r` 是破坏性操作。运行前请检查 `wavecar_list.txt`。

该脚本可识别常见任务类型，包括单点、离子优化、晶胞优化、频率、DOS、dimer、由 `lobsterin`/`LOBSTERIN` 标记的 COHP，以及 NEB 目录。

## 依赖

- Python >= 3.7
- ASE
- NumPy
- SciPy
- Matplotlib
- PyNEP，NEP calculator 工作流可选依赖

## 测试

回归测试使用标准库 `unittest`：

```bash
python -m unittest discover -s tests -v
python -m compileall -q qctools tests
```

## 许可证

MIT License.
