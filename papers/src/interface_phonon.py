#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""interface_phonon.py —— Zd012：候选对的界面声子

从两张表里取候选 A/B 相的 omega_log：
  A 端：C（金刚石）、CN、CHN、BN、BC
  B 端：TiC、ZrC、MoC、CSc、VN、NbN、TiN、HfN、MoN
判据：界面模的宿主是 B 端的碳/氮化物，其 omega_log 是否 > 1091 K。
输出：data/interface_phonon.json
"""
from __future__ import annotations
import csv, json, os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "data", "interface_phonon.json")
WMIN = 1091.0
WANT = {"C", "CN", "CHN", "BN", "BC", "CTi", "CZr", "CMo", "CSc", "NV", "NNb",
        "NTi", "HfN", "MoN", "B2CN", "B4C"}


def load(path, tag):
    out = {}
    for r in csv.DictReader(open(path)):
        f = r.get("formula", "")
        if f in WANT:
            try:
                out.setdefault(f, []).append({"w": float(r["wlog_K"]),
                                             "lam": float(r["lam"]), "src": tag})
            except Exception:
                pass
    return out


def main():
    j = load("/Users/oygb/Downloads/Zero/_staging/批量Tc_全表.csv", "JARVIS")
    m = load("/Users/oygb/Downloads/Zero/_staging/MP_声子_分档_全表.csv", "MP")
    allrows = {}
    for d in (m, j):
        for k, vs in d.items():
            allrows.setdefault(k, []).extend(vs)
    rows = {}
    for k, vs in allrows.items():
        hi = max(vs, key=lambda v: v["w"]); lo = min(vs, key=lambda v: v["w"])
        rows[k] = {"w": hi["w"], "w_min": lo["w"], "lam": hi["lam"], "src": hi["src"]}
    out = {"wmin": WMIN,
           "rows": [{"f": k, **v, "pass": bool(v["w"] > WMIN)}
                    for k, v in sorted(rows.items(), key=lambda kv: -kv[1]["w"])]}
    json.dump(out, open(OUT, "w"), ensure_ascii=False, indent=1)
    print(f"{'formula':8s} {'wlog':>7s} {'lam':>6s} {'pass':>5s}  src")
    for r in out["rows"]:
        print(f"{r['f']:8s} {r['w']:7.1f} [{r['w_min']:7.1f}] {r['lam']:6.3f} {str(r['pass']):>5s}  {r['src']}")
    print("saved", OUT)


if __name__ == "__main__":
    main()
