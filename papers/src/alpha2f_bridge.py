#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""alpha2f_bridge.py —— Zd020：结构 → (λ, ω_log) 的三段桥

三段：
  (a) 键 → 力常数 → 声子频率：Morse 键 E_b=D_b(1-e^{-a(R-R0)})^2，
      Φ=E_b''(R0)=2D_b a^2，ω=sqrt(Φ/μ)=a sqrt(2D_b/μ)。
  (b) 力常数 → 声子谱：1D 双原子链的动力学矩阵，声学/光学两支与 ω_log。
  (c) 形变势 → 电声矩阵元：D=βR0t；g=D q u0，u0=sqrt(hbar/(2Mω))；
      q 的两种口径：k_F=κ/a（体相声子）与 1/a（界面局域模，保守）。
最后给出 λ=N(0)g²/(Mω²) 与 Zd009 的 C(D/t)² 的闭合关系：
      q=1/a 给 λ=c_N(D/t)²；q=k_F=κ/a 给 λ=C(D/t)²=C_N κ² (D/t)²。

输出：data/alpha2f_bridge.json
"""
from __future__ import annotations
import json
import math
import os

import numpy as np

from constants import AMU, C_N, HBAR, KAPPA2

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "data", "alpha2f_bridge.json")

# ---- (a) C–C 的 Morse 键 -----------------------------------------------------
DB_C = 3.6          # eV，键深
R0_C = 1.54e-10     # m
MU_C = 6 * AMU      # kg，C–C 伸缩的约化质量（m_C/2）
NU_C = 1332.0       # cm^-1，金刚石光学模
C_LIGHT = 2.99792458e10  # cm/s
OMEGA_C = 2 * math.pi * NU_C * C_LIGHT   # rad/s


def morse():
    # 由 ω 与 D_b 反解 a，再由 a 给 Φ
    a = OMEGA_C / math.sqrt(2 * DB_C * 1.602176634e-19 / MU_C)   # m^-1
    phi = 2 * DB_C * 1.602176634e-19 * a ** 2                    # N/m
    return {"a_per_A": a * 1e-10, "Phi_N_per_m": phi,
            "omega_cm-1": NU_C, "D_b_eV": DB_C}


# ---- (b) 1D 双原子链 ---------------------------------------------------------
def diatomic(m1=1.0, m2=2.0, K=1.0, nq=4000):
    """ω² = K(1/m1+1/m2) ± K sqrt((1/m1+1/m2)² - 4 sin²(qa/2)/(m1 m2))，q∈(0,π/a]。"""
    q = np.linspace(0, math.pi, nq, endpoint=False) + math.pi / nq
    s = np.sin(q / 2) ** 2
    A = K * (1 / m1 + 1 / m2)
    disc = np.sqrt(np.clip(A ** 2 - 4 * K ** 2 * s / (m1 * m2), 0, None))
    w_ac = np.sqrt(np.clip(A - disc, 0, None))
    w_op = np.sqrt(np.clip(A + disc, 0, None))
    ws = np.concatenate([w_ac, w_op])
    ws = ws[ws > 1e-12]
    return {"acoustic_min": float(w_ac.min()), "acoustic_max": float(w_ac.max()),
            "optical_min": float(w_op.min()), "optical_max": float(w_op.max()),
            "gap": float(w_op.min() - w_ac.max()),
            "wlog": float(math.exp(np.mean(np.log(ws))))}


# ---- (c) 形变势与 g ---------------------------------------------------------
def deform(beta=1.5, R0=1.54, t=2.5, M_amu=12.0, omega_cm=1332.0, a_ang=3.0):
    D = beta * R0 * t                       # eV（应变导数 ∂E/∂(δR/R0)）
    om = 2 * math.pi * omega_cm * C_LIGHT
    u0 = math.sqrt(HBAR / (2 * M_amu * AMU * om))   # m
    q_bond = 1 / (a_ang * 1e-10)                    # m^-1
    q_kf = (3 * math.pi ** 2) ** (1 / 3) / (a_ang * 1e-10)
    g_at_R0 = D * (1 / (R0 * 1e-10)) * u0
    return {"D_eV": D, "u0_A": u0 * 1e10,
            "g_bond_eV": D * q_bond * u0,
            "g_at_R0_eV": g_at_R0,
            "g_kF_eV": D * q_kf * u0,
            "g_no_length_factor": D * u0 * 1e10}   # 省略 1/a：单位 eV·Å，不是能量


def main():
    m = morse()
    d = diatomic()
    g = deform()
    # λ 的两种口径：N(0)=c_N/t、Mω²=t/a²
    t = 2.5
    lam_bond = C_N * (g["D_eV"] / t) ** 2
    lam_kf = C_N * KAPPA2 * (g["D_eV"] / t) ** 2
    out = {"morse": {k: round(v, 4) for k, v in m.items()},
           "diatomic": {k: round(v, 4) for k, v in d.items()},
           "deform": {k: round(v, 4) for k, v in g.items()},
           "lambda": {"q=1/a": round(lam_bond, 4), "q=kF": round(lam_kf, 4),
                      "C_gL1_times_(D/t)^2": round(C_N * KAPPA2 * (g["D_eV"] / t) ** 2, 4)},
           "ratios": {"kappa2": round(KAPPA2, 4), "g_kF_over_g_bond": round(g["g_kF_eV"] / g["g_bond_eV"], 3)}}
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    json.dump(out, open(OUT, "w"), ensure_ascii=False, indent=1)
    print(json.dumps(out, ensure_ascii=False, indent=1))
    print("saved", OUT)


if __name__ == "__main__":
    main()
