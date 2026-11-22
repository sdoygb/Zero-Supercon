#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""anharmonic.py —— Zd014：边权对数导数 n 的分布

lambda = C n^2，C = c_N kappa^2 g_L，g_L = 1（纵声子）或 1/3（三维键向平均）。
输出：data/anharmonic.json
"""
from __future__ import annotations
import csv, json, math, os, statistics as st

from constants import C_N

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "data", "anharmonic.json")
CN = C_N  # 单自旋态密度（λ 链，见 constants.py）
KAPPA2 = (3 * math.pi ** 2) ** (2 / 3)


def load():
    out = []
    for r in csv.DictReader(open("/Users/oygb/Downloads/Zero/_staging/批量Tc_全表.csv")):
        try:
            l = float(r["lam"])
        except Exception:
            continue
        if l > 0:
            out.append(l)
    return sorted(out)


def main():
    lams = load(); n = len(lams)
    out = {}
    for gL in (1.0, 1 / 3):
        C = CN * KAPPA2 * gL
        ne = sorted(math.sqrt(l / C) for l in lams)
        out[str(round(gL, 2))] = {
            "C": round(C, 3), "n": n,
            "median": round(st.median(ne), 3), "q25": round(ne[n // 4], 3),
            "q75": round(ne[3 * n // 4], 3), "max": round(ne[-1], 2),
            "n_at_lam": {str(l): round(math.sqrt(l / C), 2) for l in (2, 3, 5)},
        }
    json.dump(out, open(OUT, "w"), ensure_ascii=False, indent=1)
    print(json.dumps(out, ensure_ascii=False, indent=1))
    print("saved", OUT)


if __name__ == "__main__":
    main()
