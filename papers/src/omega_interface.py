#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""omega_interface.py —— 界面模质量的三条取法与 Tc 的带

问题：界面键 C–M 的振子质量取什么？
  ①整原子：μ_I = m_C m_M/(m_C+m_M)（与 Zd029 的旧值一致）
  ②界面键份额：界面碳参与 1 个界面键 + (ν−1) 个面内键，界面键份额 m_C/ν_eff
  ③价加权：μ_I = (m_C/ν) m_M/((m_C/ν)+m_M)

核验：用 C–C 光学模的 Raman 实测反解 μ_eff（见 omega_fit.py），
界面的"质量份额"应有类似的量级；取①②③给 ω_I 的带。
"""
from __future__ import annotations
import json, math, os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "data", "omega_interface.json")

EV_A2 = 1.602176634e-19 / 1e-20
AMU = 1.66054e-27
HB = 6.582119569e-16
KB = 8.617333e-5


def omega_K(KI_Nm, mu_amu):
    return HB * math.sqrt(KI_Nm / (mu_amu * AMU)) / KB


def ad(lam, w, mu=0.12):
    den = lam - mu * (1 + 0.62 * lam)
    return w / 1.2 * math.exp(-1.04 * (1 + lam) / den)


def main():
    mC, mAg = 12.0, 107.87
    # sp2 石墨烯：K_I=651 N/m（D_b=4.2, a=2.2）；ν=3
    KI, nu = 651.0, 3
    parts = [
        ("整原子", mC),
        ("界面键份额 3/4", mC * 0.75),
        ("价加权 m_C/ν", mC / nu),
    ]
    rows = []
    print("sp2 石墨烯（K_I=651 N/m）的三条取法：")
    for tag, part in parts:
        muI = part * mAg / (part + mAg)
        w = omega_K(KI, muI)
        r = {"mass_model": tag, "mass_partial_amu": round(part, 2),
             "mu_I_amu": round(muI, 2), "omega_I_K": round(w)}
        for lam in (1.338, 1.703):
            r[f"Tc_AD_lam{lam}"] = round(ad(lam, w))
        rows.append(r)
        print("  %-16s 份额 %5.2f amu → μ_I=%5.2f → ω_I=%4.0f K；Tc(λ=1.34/1.70)=%3d/%3d K"
              % (tag, part, muI, w, r["Tc_AD_lam1.338"], r["Tc_AD_lam1.703"]))
    print("\n旧值（Zd029）：μ_I=10.80 → ω_I=1021 K")
    print("Raman 反解（omega_fit.py）：C–C 的 μ_eff=4.41–4.43 amu，落在 [1.5,6.0]")
    json.dump({"K_I_Nm": KI, "nu": nu, "rows": rows}, open(OUT, "w"), ensure_ascii=False, indent=1)
    print("saved", OUT)


if __name__ == "__main__":
    main()
