#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""bond_depth_omega.py —— 用 Zc036 定理 1.1 定 C–Ag 键深与界面模

链（全部可查）：
  Zc035 §5.7 / Zd024：C–Ag 跳跃 t_I = 1.72 eV，n_I = -dln|β|/dlnR = 5.03
  Zc036 定理 1.1：D_ij = 2 p_ij |β_ij|
  Zc034 §4.5：Φ = 2 D a²（Morse），ω = sqrt(Φ/μ)
  Zc038 命题 6.4：振子质量份额 m_C/ν（价闭合）
输入分档：p_ij（键级，C）、a（键程，C）、ν（价数，I）、t_I（Zd024 的导出）
"""
from __future__ import annotations
import json, math, os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "data", "bond_depth_omega.json")

EV_A2 = 1.602176634e-19 / 1e-20
AMU = 1.66054e-27
HB = 6.582119569e-16        # eV·s
KB = 8.617333e-5
CM = 1.4387769


def omega_kg(K_Nm, mu_amu):
    return HB * math.sqrt(K_Nm / (mu_amu * AMU)) / KB / CM     # cm^-1


def ad(lam, wK, mu=0.12):
    den = lam - mu * (1 + 0.62 * lam)
    return wK / 1.2 * math.exp(-1.04 * (1 + lam) / den)


def main():
    mC, mAg = 12.0, 107.87
    t_I, n_I = 1.72, 5.03     # Zd024 / Zc035 §5.7
    rows = []
    print("C–Ag：t_I=%.2f eV，n_I=%.2f（Zd024）" % (t_I, n_I))
    for p in (0.3, 0.5, 1.0):
        D = 2 * p * t_I
        for aname, a in (("C–C 反解 2.33", 2.33), ("Morse 2.0", 2.0)):
            K = 2 * D * a * a * EV_A2
            for mname, mu in (("整原子 μ=10.80", mC * mAg / (mC + mAg)),
                              ("价加权 μ=3.86", (mC / 3) * mAg / ((mC / 3) + mAg))):
                wcm = omega_kg(K, mu)
                wK = wcm * CM
                rows.append({"p": p, "D_eV": round(D, 2), "a_A": a, "K_Nm": round(K),
                             "mass": mname.split()[0], "mu_amu": round(mu, 2),
                             "omega_cm": round(wcm), "omega_K": round(wK),
                             "Tc_AD_134": round(ad(1.338, wK)), "Tc_AD_170": round(ad(1.703, wK))})
                print("  p=%.1f D=%.2feV a=%.2f  K=%4.0f N/m %-6s ω=%5.0f cm^-1 (%4.0f K)  Tc=%.0f/%.0f"
                      % (p, D, a, K, mname.split()[0], wcm, wK,
                         ad(1.338, wK), ad(1.703, wK)))
    json.dump({"t_I_eV": t_I, "n_I": n_I, "rows": rows}, open(OUT, "w"), ensure_ascii=False, indent=1)
    print("\nsaved", OUT)


if __name__ == "__main__":
    main()
