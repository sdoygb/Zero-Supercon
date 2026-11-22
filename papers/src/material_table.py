#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""material_table.py —— Zd020 §5.5：把 Φ／β／D／g 的 C 格用一组真实材料填上

口径（与 Zd020 §1–§5 一致）：
  t      裸带（Zd002 命题 2.2 的舞台约定）：t = hbar²k_F²/(2m_e)，k_F=κ/a，κ=(3π²)^{1/3}
  n_H    结构估计的 D/t：Harrison 的键长指数（sp: t∝a^{-2}；d: t∝a^{-5}）
  n_eff  数据反解的 D/t：由库表的 λ 反解（n=sqrt(λ/C)，C=c_N κ²，见 constants.py）
  D      = n·t；两条口径各给一个，比值就是"填充/几何"因子
  Φ      = μω²（由声子频率读出的力常数；μ=约化质量，键的伸缩模取 m/2）
  u0     = sqrt(hbar/(2Mω))；g = D·u0/a

输入分档：
  F  库表（`_staging/MP_声子_分档_全表.csv` 的 ω_log 与 λ；`_staging/批量Tc_全表.csv` 的 λ）
  C  结构/带参数（a、R₀、M、Harrison 指数；界面行按 Zd018 的 D_I 反解 t_I）

输出：data/material_table.json
"""
from __future__ import annotations
import csv
import json
import math
import os

import numpy as np

from constants import AMU, C_LAM_G1, C_N, HBAR, KAPPA2

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "data", "material_table.json")
MP = "/Users/oygb/Downloads/Zero/_staging/MP_声子_分档_全表.csv"
JARVIS = "/Users/oygb/Downloads/Zero/_staging/批量Tc_全表.csv"

KAPPA = (3 * math.pi ** 2) ** (1 / 3)
HB2_2ME_EVA2 = 3.80998  # hbar²/(2 m_e) = 3.81 eV·Å²
KB_EV = 8.617333e-5

# (名称, 结构, a[Å], M[amu], R0[Å], n_H, ω_log 取哪个 MP 条目)
BULK = [
    ("金刚石 C", "diamond", 3.567, 12.011, 3.567 * math.sqrt(3) / 4, 2, 1336.5),
    ("Si", "diamond", 5.431, 28.086, 5.431 * math.sqrt(3) / 4, 2, 403.1),
    ("W", "bcc", 3.165, 183.84, 3.165 * math.sqrt(3) / 2, 5, 171.1),
    ("Ta", "bcc", 3.306, 180.95, 3.306 * math.sqrt(3) / 2, 5, 129.6),
    ("Nb", "bcc", 3.301, 92.906, 3.301 * math.sqrt(3) / 2, 5, 175.4),
    ("Mo", "bcc", 3.147, 95.95, 3.147 * math.sqrt(3) / 2, 5, 230.9),
    ("Cu", "fcc", 3.615, 63.546, 3.615 / math.sqrt(2), 2, 223.3),
]

INTERFACE = {"name": "C–W 界面", "R0": 2.07, "M": 12.011, "omega_K": 1500.0,
             "a": 3.0, "D_I": 5.0, "beta": 1.5}


def load_mp():
    out = {}
    for r in csv.DictReader(open(MP)):
        try:
            out.setdefault(r["formula"], []).append(
                (float(r["wlog_K"]), float(r["lam"])))
        except Exception:
            pass
    return out


def load_jarvis():
    out = {}
    for r in csv.DictReader(open(JARVIS)):
        try:
            out.setdefault(r["formula"], []).append(
                (float(r["wlog_K"]), float(r["lam"])))
        except Exception:
            pass
    return out


def pick(mp, formula, wlog):
    """取与给定 ω_log 最接近的 MP 条目。"""
    cand = mp.get(formula) or []
    if not cand:
        return None
    return min(cand, key=lambda t: abs(t[0] - wlog))


def u0_A(M_amu, omega_K):
    om = KB_EV * omega_K * 1.602176634e-19 / HBAR   # rad/s
    return math.sqrt(HBAR / (2 * M_amu * AMU * om)) * 1e10


def main():
    mp = load_mp()
    jv = load_jarvis()
    rows = []
    for name, struct, a, M, R0, nH, wtarget in BULK:
        got = pick(mp, name.split()[-1] if name.startswith("金刚石") else name, wtarget)
        if got is None:
            raise SystemExit(f"MP 表里没有 {name}")
        wlog, lam = got
        t = HB2_2ME_EVA2 * (KAPPA / a) ** 2
        n_eff = math.sqrt(lam / C_LAM_G1)
        D_eff = n_eff * t
        D_H = nH * t
        u = u0_A(M, wlog)
        g_eff = D_eff * u / a
        g_H = D_H * u / a
        om = KB_EV * wlog * 1.602176634e-19 / HBAR
        phi = (M / 2) * AMU * om ** 2   # Φ=μω², μ=m/2（N/m）
        rows.append({
            "material": name, "structure": struct, "a_A": a, "M_amu": M, "R0_A": round(R0, 4),
            "wlog_K": round(wlog, 1), "lam_table": round(lam, 4),
            "t_bare_eV": round(t, 3), "n_H": nH, "n_eff": round(n_eff, 3),
            "D_H_eV": round(D_H, 2), "D_eff_eV": round(D_eff, 2),
            "u0_A": round(u, 4), "Phi_N_m": round(phi, 0),
            "g_H_eV": round(g_H, 4), "g_eff_eV": round(g_eff, 4),
            "lam_model_H": round(C_LAM_G1 * nH ** 2, 2),
            "lam_check": round(C_LAM_G1 * n_eff ** 2, 4),
            "jarvis_lam": [round(x[1], 3) for x in jv.get(name.split()[-1], [])][:3],
        })
    # 界面：t_I 由 Zd018 的 D_I=5 eV 与 β=1.5 反解
    it = INTERFACE
    tI = it["D_I"] / (it["beta"] * it["R0"])
    nI = it["beta"] * it["R0"]
    u = u0_A(it["M"], it["omega_K"])
    gI = it["D_I"] * u / it["a"]
    require = [(0.43 / gI) ** 2, (0.63 / gI) ** 2]
    inter = {"material": it["name"], "R0_A": it["R0"], "M_amu": it["M"],
             "omega_K": it["omega_K"], "a_A": it["a"],
             "t_I_eV": round(tI, 2), "n_I": round(nI, 2), "D_I_eV": it["D_I"],
             "u0_A": round(u, 4), "g_I_eV": round(gI, 4),
             "M_eff_req_for_g43": round(require[0], 0),
             "M_eff_req_for_g63": round(require[1], 0)}

    out = {"bulk": rows, "interface": inter,
           "notes": {
               "t": "bare band: t=3.81*(kappa/a)^2 eV (Zd002 命题 2.2)",
               "n_H": "Harrison exponent (sp 2, d 5) — C",
               "n_eff": "data-inverted D/t: sqrt(lambda/C_LAM_G1); lambda from MP table",
               "Phi": "mu*omega^2 with mu=M/2 (bond stretch)",
               "g": "D*u0/a; u0=sqrt(hbar/(2 M omega))",
           }}
    json.dump(out, open(OUT, "w"), ensure_ascii=False, indent=1)
    hdr = ("材料      a     M     ω_log   λ表    t裸带  n_H  n_eff  D_H   D_eff  u0    Φ     g_H   g_eff  λ_H")
    print(hdr)
    for r in rows:
        print(f"{r['material']:8s} {r['a_A']:5.3f} {r['M_amu']:6.1f} {r['wlog_K']:6.1f} "
              f"{r['lam_table']:5.3f} {r['t_bare_eV']:5.2f} {r['n_H']:3.0f} {r['n_eff']:5.2f} "
              f"{r['D_H_eV']:5.2f} {r['D_eff_eV']:5.2f} {r['u0_A']:5.4f} {r['Phi_N_m']:5.0f} "
              f"{r['g_H_eV']:6.4f} {r['g_eff_eV']:6.4f} {r['lam_model_H']:5.1f}")
    print()
    print("界面:", json.dumps(inter, ensure_ascii=False))
    print("saved", OUT)


if __name__ == "__main__":
    main()
