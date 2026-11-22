#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""rt_window.py —— 常温窗口：Tc 公式的必要条件与材料图谱

口径（Zb072 §11.6–§11.11）：
  akF   = (3π²)^{1/3}
  e2hat = 2(a/a_B)/(akF)²
  μ     = N0·4πe2hat/(⟨q²⟩+4πe2hat·N0),  N0=2·C_N=0.2854（总态密度）, ⟨q²⟩=18.97
  L     = akF·√(M/2m_e)          （M 以 m_e 计）
  μ*    = μ/(1+μ·ln L)
  Tc    = (ω_log/1.2)·exp[−1.04(1+λ)/(λ−μ*(1+0.62λ))]

输出：data/rt_window.json
"""
from __future__ import annotations
import csv, json, math, os, statistics as st

from constants import N_TOT

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = "/Users/oygb/Downloads/Zero/_staging/批量Tc_全表.csv"
OUT = os.path.join(ROOT, "data", "rt_window.json")

AKF = (3 * math.pi ** 2) ** (1 / 3)
AB = 0.529177
N0 = N_TOT  # 2·C_N：TF 屏蔽的 dn/dμ 计两自旋（见 constants.py）
Q2 = 18.97
ME_U = 1 / 1822.888


def mu_of(a):
    e2 = 2 * (a / AB) / AKF ** 2
    return N0 * 4 * math.pi * e2 / (Q2 + 4 * math.pi * e2 * N0)


def mustar_of(a, M_amu):
    M = M_amu / ME_U          # in m_e
    L = AKF * math.sqrt(M / 2.0)
    mu = mu_of(a)
    return mu / (1 + mu * math.log(L)), mu, math.log(L)


def ad(lam, wlog, mustar):
    den = lam - mustar * (1 + 0.62 * lam)
    if den <= 0:
        return float("nan")
    return wlog / 1.2 * math.exp(-1.04 * (1 + lam) / den)


def wlog_req(lam, mustar, tc=300.0):
    den = lam - mustar * (1 + 0.62 * lam)
    if den <= 1e-6:
        return float("inf")
    ex = 1.04 * (1 + lam) / den
    if ex > 300:
        return float("inf")
    return 1.2 * tc * math.exp(ex)


def main():
    rows = list(csv.DictReader(open(SRC)))
    for r in rows:
        r["lam"] = float(r["lam"])
        r["wlog"] = float(r["wlog_K"])
        r["a"] = float(r["a_A"])
        r["M"] = float(r["M_amu"])
        r["ms"], r["mu"], r["lnL"] = mustar_of(r["a"], r["M"])
        r["tc"] = ad(r["lam"], r["wlog"], r["ms"]) if r["lam"] > 0 else 0.0

    # ---- 1. formula check ----
    dmu = [abs(r["ms"] - float(r["mustar"])) for r in rows]
    dtc = [abs(r["tc"] - float(r["Tc_K"])) / max(float(r["Tc_K"]), 1e-9)
           for r in rows if r["Tc_K"] and float(r["Tc_K"]) > 1e-6]

    # ---- 2. ceiling ----
    ceiling = []
    for ms in (0.085, 0.10, 0.12):
        for lam in (1.0, 2.0, 3.0, 5.0, 10.0, 1e6):
            ceiling.append({"mustar": ms, "lam": lam,
                            "tc_over_wlog": ad(lam, 1.0, ms)})

    # ---- 3. room-temperature window ----
    window = []
    for ms in (0.085, 0.10, 0.12):
        wmin, lmin = float("inf"), None
        for i in range(1, 20001):
            lam = 0.02 * (1.0004 ** i)
            w = wlog_req(lam, ms)
            window.append({"mustar": ms, "lam": round(lam, 4), "wlog_req": w})
            if w < wmin:
                wmin, lmin = w, lam
        asym = 1.2 * 300.0 * math.exp(1.04 / (1 - 0.62 * ms))
        print(f"mu*={ms:.3f}: min wlog_req = {wmin:.1f} K at lam = {lmin:.3f} "
              f"(asymptote {asym:.1f} K)")

    # ---- 4. material map ----
    top = sorted(rows, key=lambda r: -r["tc"])[:15]
    top_out = [{"formula": r["formula"], "a": r["a"], "M": r["M"],
                "lam": round(r["lam"], 4), "wlog": round(r["wlog"], 1),
                "mustar": round(r["ms"], 4), "Tc": round(r["tc"], 2),
                "wlog_req300": round(wlog_req(r["lam"], r["ms"]), 1),
                "factor": round(wlog_req(r["lam"], r["ms"]) / r["wlog"], 2)}
               for r in top]

    # Pareto frontier in (lam, wlog): max wlog among lam >= x
    lam_max = max(r["lam"] for r in rows)
    bins = [0.2, 0.4, 0.6, 0.8, 1.0, 1.5, 2.0, 3.0, 5.0]
    pareto = []
    for x in bins:
        sub = [r["wlog"] for r in rows if r["lam"] >= x]
        if sub:
            pareto.append({"lam_ge": x, "n": len(sub),
                           "wlog_max": round(max(sub), 1),
                           "wlog_med": round(st.median(sub), 1)})

    # how far each material is: needed multiplier on wlog at its own lam
    fac = [(wlog_req(r["lam"], r["ms"]) / r["wlog"], r["formula"], r["lam"],
            r["wlog"], r["ms"], r["tc"]) for r in rows
           if r["lam"] > 0.02 and r["wlog"] > 1.0]
    fac.sort()
    closest = [{"formula": f, "lam": round(l, 4), "wlog": round(w, 1),
                "mustar": round(m, 4), "Tc": round(t, 2), "factor": round(c, 2)}
               for c, f, l, w, m, t in fac[:12]]

    # ---- 5. adiabatic/strong-coupling sanity: lam vs wlog correlation ----
    xs = [math.log(r["lam"]) for r in rows if r["lam"] > 0.02 and r["wlog"] > 1]
    ys = [math.log(r["wlog"]) for r in rows if r["lam"] > 0.02 and r["wlog"] > 1]
    mx, my = st.mean(xs), st.mean(ys)
    cov = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    corr = cov / math.sqrt(sum((x - mx) ** 2 for x in xs) * sum((y - my) ** 2 for y in ys))

    out = {
        "n": len(rows),
        "formula_check": {"max_dmustar": max(dmu), "median_rel_dTc": st.median(dtc)},
        "ceiling": ceiling,
        "window_min": {str(ms): round(min(w["wlog_req"] for w in window
                                         if w["mustar"] == ms), 1)
                       for ms in (0.085, 0.10, 0.12)},
        "window_asym": {str(ms): round(1.2 * 300.0 * math.exp(1.04 / (1 - 0.62 * ms)), 1)
                        for ms in (0.085, 0.10, 0.12)},
        "top": top_out,
        "pareto": pareto,
        "closest": closest,
        "corr_lnlam_lnwlog": round(corr, 4),
        "db_max_Tc": round(max(r["tc"] for r in rows), 3),
    }
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    json.dump(out, open(OUT, "w"), ensure_ascii=False, indent=1)

    print("check max|dmu*| =", max(dmu), " median rel dTc =", st.median(dtc))
    print("db max Tc =", out["db_max_Tc"], " corr(ln lam, ln wlog) =", round(corr, 3))
    print("\nTOP:")
    for t in top_out[:8]:
        print(f"  {t['formula']:8s} lam={t['lam']:.3f} wlog={t['wlog']:6.1f} "
              f"Tc={t['Tc']:6.2f} need_wlog={t['wlog_req300']:8.1f} x{t['factor']}")
    print("\nPARETO:")
    for p in pareto:
        print(f"  lam>={p['lam_ge']:.1f}: n={p['n']:4d} wlog_max={p['wlog_max']:6.1f}")
    print("\nCLOSEST:")
    for c in closest[:8]:
        print(f"  {c['formula']:8s} lam={c['lam']:.3f} wlog={c['wlog']:6.1f} "
              f"Tc={c['Tc']:6.2f} factor={c['factor']}")
    print("saved", OUT)


if __name__ == "__main__":
    main()
