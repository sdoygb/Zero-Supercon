#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""natural_lambda.py —— Zd009：λ 的分布与"自然大小"

只做一件事：把两张表的 λ 分布统计出来（中位、四分位、尾部比例），
用来对照"λ 是 (D/t)² 的纯数"这一标度结论。
输出：data/natural_lambda.json
"""
from __future__ import annotations
import csv, json, os, statistics as st

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "data", "natural_lambda.json")


def load(path):
    out = []
    for r in csv.DictReader(open(path)):
        try:
            l = float(r["lam"])
        except Exception:
            continue
        if l > 0:
            out.append(l)
    return sorted(out)


def stats(xs):
    n = len(xs)
    return {"n": n, "median": round(st.median(xs), 3),
            "q25": round(xs[n // 4], 3), "q75": round(xs[3 * n // 4], 3),
            "max": round(max(xs), 3),
            "frac_gt": {str(t): round(sum(1 for x in xs if x > t) / n, 3)
                        for t in (0.5, 1, 2, 3, 5)}}


def main():
    a = load("/Users/oygb/Downloads/Zero/_staging/批量Tc_全表.csv")
    b = load("/Users/oygb/Downloads/Zero/_staging/MP_声子_分档_全表.csv")
    out = {"JARVIS": stats(a), "MP": stats(b)}
    json.dump(out, open(OUT, "w"), ensure_ascii=False, indent=1)
    print(json.dumps(out, ensure_ascii=False, indent=1))
    print("saved", OUT)


if __name__ == "__main__":
    main()
