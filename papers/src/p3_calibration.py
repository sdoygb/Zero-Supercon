#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""p3_calibration.py —— P3：μ* 的经验标定 与 二维 ρ_s 的量级

(a) μ* 的经验分布：对 74 条有实测 Tc 的材料，用 AD 形式反解
      μ* = [λ - 1.04(1+λ)/ln(ω_log/(1.2Tc))] / (1+0.62λ)
    与 Zd002 的公式值（0.077-0.113）比。
    对照：用库的计算值 Tc 反解（应按构造集中在 0.085）。
(b) 二维：T_BKT = (π/2)ρ_s，ρ_s = ħ²n/(4m*) -> T_BKT = πħ²n/(8m*k_B)。

输出：data/p3_calibration.json
"""
from __future__ import annotations
import json
import math
import os
import statistics as st

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "data", "p3_calibration.json")
EXP74 = "/Users/oygb/Downloads/Zero/_staging/Eliashberg_74条_试算_20261106.json"
JARVIS = "/Users/oygb/Downloads/Zero/_staging/批量Tc_全表.csv"

HBAR_SI = 1.054571817e-34
KB_SI = 1.380649e-23
ME_SI = 9.1093837e-31


def mustar(lam, wlog_K, Tc_K):
    """由 AD 形式反解 μ*；不可解返回 None。"""
    if Tc_K <= 0 or wlog_K <= 0:
        return None
    arg = wlog_K / (1.2 * Tc_K)
    if arg <= 1:
        return None
    ln = math.log(arg)
    if ln <= 0:
        return None
    mu = (lam - 1.04 * (1 + lam) / ln) / (1 + 0.62 * lam)
    return mu


def t_bkt(n_2d_cm2, m_star):
    n = n_2d_cm2 * 1e4                      # m^-2
    rho = HBAR_SI ** 2 * n / (4 * m_star * ME_SI)   # J（=K·k_B）
    return math.pi * HBAR_SI ** 2 * n / (8 * m_star * ME_SI * KB_SI)


def main():
    rows = {r["jid"]: r for r in open_csv()}
    d = json.load(open(EXP74))
    exp = []
    db = []
    for r in d:
        j = r["jid"]
        t = rows.get(j)
        if not t:
            continue
        lam, wlog = t["lam"], t["wlog"]
        me = mustar(lam, wlog, r.get("Tc_exp"))
        md = mustar(lam, wlog, r.get("Tc_db"))
        if me is not None:
            exp.append({"jid": j, "formula": t["formula"], "lam": lam, "wlog": wlog,
                        "Tc_exp": r["Tc_exp"], "mustar": round(me, 4)})
        if md is not None:
            db.append(md)

    def stats(vs):
        vs = sorted(vs)
        q = lambda p: vs[min(len(vs) - 1, int(p * len(vs)))]
        return {"n": len(vs), "median": round(q(0.5), 4), "q25": round(q(0.25), 4),
                "q75": round(q(0.75), 4), "min": round(vs[0], 4), "max": round(vs[-1], 4),
                "in_0.077_0.113": round(sum(1 for v in vs if 0.077 <= v <= 0.113) / len(vs), 3)}

    rho = [{"n_cm2": n, "m_star": m, "T_BKT_K": round(t_bkt(n, m), 1)}
           for n in (3e13, 1e14, 3e14) for m in (1.0, 3.0, 6.0)]
    out = {"exp74": stats([e["mustar"] for e in exp]),
           "db_control": stats(db),
           "exp_rows": exp, "bkt": rho}
    json.dump(out, open(OUT, "w"), ensure_ascii=False, indent=1)
    print("μ* 经验分布（74 条实测）:", out["exp74"])
    print("μ* 由库值反解（对照）:", out["db_control"])
    print("T_BKT:")
    for r in rho:
        print("  n=%.0e cm^-2  m*=%.0f  -> %.0f K" % (r["n_cm2"], r["m_star"], r["T_BKT_K"]))
    print("saved", OUT)


def open_csv():
    import csv
    for r in csv.DictReader(open(JARVIS)):
        try:
            r["lam"] = float(r["lam"])
            r["wlog"] = float(r["wlog_K"])
        except Exception:
            continue
        yield r


if __name__ == "__main__":
    main()
