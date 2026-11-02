# Zero-Supercon

**从 α²F 到临界温度的一条零参数链。** 材料只提供声子谱与结构；库仑赝势 μ\*、X（ε_F/ω_D）、T_c 的其余部分全部由式子给出。

```
α²F(ω) ──→ (λ, ω_log) ──→ μ* ──→ T_c
            §11.13          导出自 (a, M) 的几何量
```

## 安装

```bash
pip install -e .
```

依赖只有 `numpy`。

## 用法

```python
from zero_supercon import ChainInput, chain

inp = ChainInput(a=3.21, M_amu=15.3, lam=0.70, wlog_K=846.0)
out = chain(inp)
print(out.mu_star, out.Tc_K)      # 0.0989  29.68
```

从 α²F 算：

```python
import numpy as np
from zero_supercon import from_a2f, chain

omega = np.linspace(0.1, 100.0, 400)          # meV
a2f   = np.exp(-0.5 * ((omega - 40.0) / 12.0) ** 2)
out = chain(from_a2f(omega, a2f, a=3.0, M_amu=50.0))
```

批量（JARVIS-SuperconDB）：

```bash
python3 examples/batch_jarvis.py
# n=845  median=0.8844  [0.5,2]=0.920  corr(ln)=0.9943
```

## 输入的四个数

| 量 | 来源 | 说明 |
|:--|:--|:--|
| $\lambda$ | 电声耦合 | 由 α²F 算（式 11.13），或直接给 |
| $\omega_{\log}$ | 声子谱 | 同上 |
| $a$ | 结构 | 晶格常数（Å），通常取 $abc$ 平均 |
| $M$ | 成分 | 平均离子质量（u） |

**其余全部导出**：$\mu$、$\mu^*$、$\varepsilon_F/\omega_D$、$T_c$。

## 精度（实测，非估计）

| 基准 | $n$ | 中位比 | $[0.5,2]$ 内 | $\ln$ 相关 |
|:--|--:|--:|--:|--:|
| 对 JARVIS 计算值 | 845 | 0.884 | 92.0% | 0.9943 |
| 对 SuperCon 实测 | 74 | 0.992 | 52.7% | 0.243 |

**单条不确定度 ±2 倍（IQR [0.58, 1.88]）。** 这是筛选级精度，不是预测级。

## 覆盖范围

**能算的材料 = 有 α²F 的材料。** 当前公开数据里只有 1,058 条（本仓库的 `data/jarvis_supercon.json`）。

| 档 | 材料数 |
|:--|--:|
| 能算 | **1,051** |
| 只有声子、缺 λ | 26,413 |
| 只有结构 | $10^5$–$10^6$ |

## 公式与出处

每个式子都对应 Zero 库的一篇（`docs/FORMULA_MAP.md` 有完整对照）：

| 量 | 式子 | 出处 |
|:--|:--|:--|
| $\lambda,\ \omega_{\log}$ | $2\int\alpha^2F/\omega\,d\omega$ 等 | Zb072 式 (11.13) |
| $\mu$ | $N(0)\langle V_C\rangle_{\rm FS}$，$V_C=4\pi e^2/(q^2+\kappa^2)$ | Zb072 命题 11.9 |
| $\varepsilon_F/\omega_D$ | $(ak_F)\sqrt{M/2m_e}$ | Zb072 命题 11.11 |
| $\hat e^2=e^2/(at)$ | $2(a/a_B)/(ak_F)^2$ | Zb072 命题 11.12 |
| $\mu^*$ | $\mu/(1+\mu\ln(\varepsilon_F/\omega_D))$ | Zb072 命题 11.13 |
| $T_c$ | $(\omega_{\log}/1.2)\exp[-1.04(1+\lambda)/(\lambda-\mu^*(1+0.62\lambda))]$ | Zb072 式 (11.11) |

## 数据

`data/jarvis_supercon.json` — JARVIS-SuperconDB（figshare [21370572](https://figshare.com/articles/dataset/JARVIS-SuperconDB/21370572)），1,058 条，含 α²F、λ、ω_log、结构。
`data/batch_tc_table.csv` — 本工具在全部 1,058 条上的输出（16 列）。

## 测试

```bash
python3 tests/test_regression.py
# T1 μ*=0.098930  Tc=29.680818 K
# T2 n=845 median=0.8844 [0.5,2]=0.920 corr=0.9943
# 全部通过
```

## 已知局限

1. **平均场**：$T_c$ 是平均场/交叉温度，有限谱上没有真正的相变。
2. **强耦合**：$\lambda>1.5$ 时 Allen–Dynes 形式有系统偏差（$\lambda\in[3,20)$ 的中位比降到 0.98）。
3. **相的身份**：同一化学式不同结构可差 74 倍，表里的每一行是**一个结构**，不是「这个材料」。
4. **不能从化学式直接算**：缺 α²F 就没有 λ，这一层本仓库不提供。

## 引用

Zero 库 Zb072《规范场的动力学与耗散：横向模式数、逐点守恒与电阻为零的判据》§11。

## 许可

MIT