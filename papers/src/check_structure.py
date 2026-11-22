#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""check_structure.py —— Zd 正文的结构机检（AGENTS.md §0 第 3 步）

两项：
  1. `$$` 的个数为偶数（显示公式成对）；
  2. 每张 Markdown 表的列数与表头一致（按 `|` 切分，去掉首尾空列）。

用法：python3 src/check_structure.py [文件名...]（缺省检查根目录 *.md）
"""
from __future__ import annotations
import glob
import os
import sys


def split_row(line):
    cells = line.strip().strip("|").split("|")
    return [c.strip() for c in cells]


def check(path):
    lines = open(path, encoding="utf-8").read().splitlines()
    errs = []
    n_dd = sum(l.count("$$") for l in lines)
    if n_dd % 2:
        errs.append(f"$$ 个数为奇数（{n_dd}）")
    i = 0
    while i < len(lines):
        if lines[i].strip().startswith("|") and i + 1 < len(lines) and set(lines[i + 1].replace("|", "").replace(":", "").strip()) <= {"-", " "} and "-" in lines[i + 1]:
            ncol = len(split_row(lines[i]))
            j = i + 2
            while j < len(lines) and lines[j].strip().startswith("|"):
                n = len(split_row(lines[j]))
                if n != ncol:
                    errs.append(f"第 {j + 1} 行表格列数 {n} ≠ 表头 {ncol}")
                j += 1
            i = j
        else:
            i += 1
    return errs


def main():
    files = sys.argv[1:] or sorted(glob.glob("*.md"))
    bad = 0
    for f in files:
        errs = check(f)
        if errs:
            bad += 1
            print(f"FAIL {os.path.basename(f)}")
            for e in errs:
                print("   ", e)
    print(f"{len(files) - bad}/{len(files)} 通过")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
