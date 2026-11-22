#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""stability.py —— Zd015：稳定性标记与 lambda 的联合分布

数据：JARVIS 原始 JSON（含 stability 字段：stable / unstable）。
输出：data/stability.json
"""
from __future__ import annotations
import json, os, statistics as st

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = "/Users/oygb/Downloads/Zero/_staging/材料库/jarvis_epc_data_figshare_1058.json"
OUT = os.path.join(ROOT, "data", "stability.json")
WMIN = 1091.0


def stats(xs):
    n = len(xs)
    if n == 0:
        return None
    xs = sorted(xs)
    return {"n": n, "median": round(st.median(xs), 3),
            "q75": round(xs[3 * n // 4], 3), "max": round(max(xs), 2),
            "frac_gt": {str(t): round(sum(1 for x in xs if x > t) / n, 3)
                        for t in (1, 2, 3, 5)}}


def main():
    d = json.load(open(SRC))
    out = {}
    for flag in ("stable", "unstable"):
        sub = [x for x in d if x.get("stability") == flag and x["lamb"] > 0]
        out[flag] = {"lam": stats([x["lamb"] for x in sub]),
                     "wlog": stats([x["wlog"] for x in sub]),
                     "joint_lam2_w1091": sum(1 for x in sub if x["lamb"] > 2 and x["wlog"] > WMIN)}
    # 稳定集合里 lam>2 的名单
    hi = sorted([x for x in d if x.get("stability") == "stable" and x["lamb"] > 2],
                key=lambda x: -x["lamb"])
    out["stable_hi"] = [{"jid": x["jid"], "lam": round(x["lamb"], 2),
                         "wlog": round(x["wlog"], 1)} for x in hi]
    json.dump(out, open(OUT, "w"), ensure_ascii=False, indent=1)
    for flag in ("stable", "unstable"):
        print(flag, json.dumps(out[flag], ensure_ascii=False))
    print("stable lam>2 list:")
    for x in out["stable_hi"]:
        print("  ", x)
    print("saved", OUT)


if __name__ == "__main__":
    main()
