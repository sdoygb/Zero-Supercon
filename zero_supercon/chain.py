# -*- coding: utf-8 -*-
"""
chain.py —— Zb072 §11 的核心链：从 (λ, ω_log, a, M) 到 T_c

公式索引（对应 Zero 库 Zb072）：

  式 (11.13)  λ = 2∫ α²F/ω dω
              ω_log = exp[(2/λ)∫ (α²F/ω) ln ω dω]
  命题 11.11  ε_F/ω_D = (a k_F) √(M/2m_e)
  命题 11.12  ê² = e²/(a t) = 2(a/a_B)/(a k_F)²        （取裸带 t = ε_F）
  命题 11.9   μ = N(0)·4πê²/(⟨q²⟩ + 4πê²N(0))
  命题 11.13  μ* = μ / (1 + μ ln(ε_F/ω_D))
  式 (11.11)  T_c = (ω_log/1.2) exp[−1.04(1+λ)/(λ − μ*(1+0.62λ))]

所有常数都是无量纲或由 Zb005/Zb036/Zb065 给出的纯数；
材料只提供 α²F 与结构。
"""
from __future__ import annotations

from dataclasses import dataclass, asdict
from typing import Optional

import numpy as np

# ---------------------------------------------------------------- 常数
# 编号即出处：见 docs/FORMULA_MAP.md
A_BOHR = 0.529177210903          # Å，玻尔半径（单位换算，不是锚）
AK_F = (3.0 * np.pi ** 2) ** (1.0 / 3.0)   # (3π²)^{1/3} = 3.0937，半满（Zb005 的谱）
N0_LATTICE = 0.3713              # 3D 立方格半满的 N(0)（Zb072 §11.7）
Q2_LATTICE = 18.97               # ⟨q²⟩_FS（同上，半满，k_F 归一）
MEV_TO_K = 11.604518             # 1 meV = 11.6045 K
M_E_TO_AMU = 1822.888486         # 1 u = 1822.888 m_e


# ---------------------------------------------------------------- 数据类
@dataclass
class ChainInput:
    """一条材料的输入。λ、ω_log 可由 α²F 算，也可直接给。"""
    a: float                     # 格距（Å），通常取 abc 的平均
    M_amu: float                 # 平均离子质量（原子质量单位）
    lam: float                   # λ（电声耦合）
    wlog_K: float                # ω_log（K）
    e2hat: Optional[float] = None   # ê² = e²/(at)；不给则由 a 与裸带算

    # 可选：来自结构的更精确信息
    v_atom: Optional[float] = None  # 每原子体积（Å³），用于 N(0) 若给
    n_electrons: Optional[float] = None  # 每胞价电子数


@dataclass
class ChainOutput:
    lam: float
    wlog_K: float
    mu: float
    mu_star: float
    e2hat: float
    ln_EF_over_wD: float
    Tc_K: float
    parts: dict

    def as_dict(self) -> dict:
        return asdict(self)


# ---------------------------------------------------------------- 各部
def moments_from_a2f(omega, a2f, unit_to_K: float = MEV_TO_K):
    """式 (11.13)：从 α²F 算 (λ, ω_log)。omega 与 a2f 等长，omega 单位与 unit_to_K 对应。"""
    omega = np.asarray(omega, dtype=float)
    a2f = np.asarray(a2f, dtype=float)
    m = (omega > 0) & (a2f >= 0)
    if m.sum() < 3:
        return float("nan"), float("nan")
    w = a2f[m] / omega[m]
    w = w[m] if False else w
    om = omega[m]
    den = float(np.trapezoid(w, om))
    if den <= 0:
        return float("nan"), float("nan")
    lam = 2.0 * den
    wlog = float(np.exp(np.trapezoid(w * np.log(om), om) / den)) * unit_to_K
    return lam, wlog


def mu_hat(e2hat: float, N0: float = N0_LATTICE, q2: float = Q2_LATTICE) -> float:
    """命题 11.9：μ = N(0)·⟨V_C⟩_FS，Thomas–Fermi 屏蔽。"""
    kappa2 = 4.0 * np.pi * e2hat * N0
    return N0 * 4.0 * np.pi * e2hat / (q2 + kappa2)


def mu_star(mu: float, ln_EF_over_wD: float) -> float:
    """命题 11.13：库仑赝势。"""
    return mu / (1.0 + mu * ln_EF_over_wD)


def allen_dynes(lam: float, wlog_K: float, ms: float) -> float:
    """式 (11.11)：Allen–Dynes 形式（μ* 由上面导出，不是输入）。"""
    den = lam - ms * (1.0 + 0.62 * lam)
    if den <= 0:
        return float("nan")
    return wlog_K / 1.2 * np.exp(-1.04 * (1.0 + lam) / den)


def chain(inp: ChainInput, bare_band: bool = True) -> ChainOutput:
    """完整链：结构 → μ* → T_c。

    bare_band=True 时取 t = ε_F（命题 11.12 的几何形式）；
    否则需要外部给 e2hat。
    """
    a, M = inp.a, inp.M_amu * M_E_TO_AMU

    if bare_band:
        e2hat = 2.0 * (a / A_BOHR) / AK_F ** 2
    else:
        if inp.e2hat is None:
            raise ValueError("bare_band=False 时必须给 e2hat")
        e2hat = inp.e2hat

    mu = mu_hat(e2hat)
    lnL = float(np.log(AK_F * np.sqrt(M / 2.0)))     # ln(ε_F/ω_D)
    ms = mu_star(mu, lnL)
    tc = allen_dynes(inp.lam, inp.wlog_K, ms)

    return ChainOutput(
        lam=inp.lam, wlog_K=inp.wlog_K, mu=mu, mu_star=ms, e2hat=e2hat,
        ln_EF_over_wD=lnL, Tc_K=tc,
        parts={"N0": N0_LATTICE, "q2": Q2_LATTICE, "akF": AK_F,
               "kappa2": 4.0 * np.pi * e2hat * N0_LATTICE})


# ---------------------------------------------------------------- 批量
def chain_batch(inputs, bare_band: bool = True):
    """向量化批量链：返回 (Tc 数组, μ* 数组, 明细 dict)。"""
    a = np.array([x.a for x in inputs], dtype=float)
    M = np.array([x.M_amu for x in inputs], dtype=float) * M_E_TO_AMU
    lam = np.array([x.lam for x in inputs], dtype=float)
    wl = np.array([x.wlog_K for x in inputs], dtype=float)
    if bare_band:
        e2h = 2.0 * (a / A_BOHR) / AK_F ** 2
    else:
        e2h = np.array([x.e2hat for x in inputs], dtype=float)
    k2 = 4.0 * np.pi * e2h * N0_LATTICE
    mu = N0_LATTICE * 4.0 * np.pi * e2h / (Q2_LATTICE + k2)
    lnL = np.log(AK_F * np.sqrt(M / 2.0))
    ms = mu / (1.0 + mu * lnL)
    den = lam - ms * (1.0 + 0.62 * lam)
    ok = den > 0
    tc = np.where(ok, wl / 1.2 * np.exp(-1.04 * (1.0 + lam) / np.where(ok, den, 1.0)), np.nan)
    return tc, ms, {"e2hat": e2h, "mu": mu, "ln_EF_over_wD": lnL}