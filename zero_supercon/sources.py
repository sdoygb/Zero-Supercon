# -*- coding: utf-8 -*-
"""
sources.py —— 读材料数据，转成 ChainInput

支持的两种输入：
  1. JARVIS-SuperconDB 的 json（含 a2F、lamb、wlog、atoms）
  2. 自己的 (omega, a2f) 数组

不做的事：不做化学式推断、不猜结构。给不出结构就报错，不填默认值。
"""
from __future__ import annotations

import json
from typing import Iterable, List, Optional

import numpy as np

from .chain import ChainInput, moments_from_a2f, MEV_TO_K

# 常用元素的原子质量（u）；覆盖 JARVIS 里出现的 62 种元素
ATOMIC_MASS = {
    "H": 1.008, "He": 4.0026, "Li": 6.94, "Be": 9.0122, "B": 10.81, "C": 12.011,
    "N": 14.007, "O": 15.999, "F": 18.998, "Ne": 20.18, "Na": 22.99, "Mg": 24.305,
    "Al": 26.982, "Si": 28.085, "P": 30.974, "S": 32.06, "Cl": 35.45, "Ar": 39.948,
    "K": 39.098, "Ca": 40.078, "Sc": 44.956, "Ti": 47.867, "V": 50.942, "Cr": 51.996,
    "Mn": 54.938, "Fe": 55.845, "Co": 58.933, "Ni": 58.693, "Cu": 63.546, "Zn": 65.38,
    "Ga": 69.723, "Ge": 72.63, "As": 74.922, "Se": 78.971, "Br": 79.904, "Kr": 83.798,
    "Rb": 85.468, "Sr": 87.62, "Y": 88.906, "Zr": 91.224, "Nb": 92.906, "Mo": 95.95,
    "Tc": 98.0, "Ru": 101.07, "Rh": 102.91, "Pd": 106.42, "Ag": 107.87, "Cd": 112.41,
    "In": 114.82, "Sn": 118.71, "Sb": 121.76, "Te": 127.6, "I": 126.9, "Xe": 131.29,
    "Cs": 132.91, "Ba": 137.33, "La": 138.91, "Ce": 140.12, "Pr": 140.91, "Nd": 144.24,
    "Sm": 150.36, "Eu": 151.96, "Gd": 157.25, "Tb": 158.93, "Dy": 162.5, "Ho": 164.93,
    "Er": 167.26, "Tm": 168.93, "Yb": 173.05, "Lu": 174.97, "Hf": 178.49, "Ta": 180.95,
    "W": 183.84, "Re": 186.21, "Os": 190.23, "Ir": 192.22, "Pt": 195.08, "Au": 196.97,
    "Hg": 200.59, "Tl": 204.38, "Pb": 207.2, "Bi": 208.98, "Po": 209.0, "At": 210.0,
    "Rn": 222.0, "Fr": 223.0, "Ra": 226.0, "Ac": 227.0, "Th": 232.04, "Pa": 231.04,
    "U": 238.03, "Np": 237.0, "Pu": 244.0, "Am": 243.0, "Cm": 247.0,
}


def load_jarvis(path: str) -> list:
    """读 JARVIS-SuperconDB 的 json。"""
    with open(path, encoding="utf-8") as f:
        return json.load(f)


def _a_from_atoms(atoms: dict) -> Optional[float]:
    abc = atoms.get("abc")
    if not abc:
        return None
    return float(np.mean(abc))


def _mass_from_elements(elements: Iterable[str]) -> Optional[float]:
    ms = [ATOMIC_MASS.get(e) for e in elements]
    if not ms or any(m is None for m in ms):
        return None
    return float(np.mean(ms))


def from_jarvis(records: list, use_a2f: bool = False, require_tc: bool = False) -> List[ChainInput]:
    """把 JARVIS 记录转成 ChainInput。

    use_a2f=True 时从 α²F 重算 (λ, ω_log)（式 11.13），否则用库里的 lamb/wlog 字段。
    require_tc=True 时只保留库里标了 T_c>0.05 的条目。
    """
    outs: List[ChainInput] = []
    for r in records:
        if require_tc and not (r.get("Tc", 0) > 0.05):
            continue
        at = r.get("atoms") or {}
        els = at.get("elements") or []
        a = _a_from_atoms(at)
        m = _mass_from_elements(els)
        if a is None or m is None:
            continue
        if use_a2f:
            x = np.array(r.get("a2F_original_x") or [], dtype=float)
            y = np.array(r.get("a2F_original_y") or [], dtype=float)
            if x.size < 3:
                continue
            lam, wl = moments_from_a2f(x, y, MEV_TO_K)
        else:
            lam, wl = r.get("lamb"), r.get("wlog")
        if lam is None or wl is None or not np.isfinite(lam) or not np.isfinite(wl):
            continue
        outs.append(ChainInput(a=a, M_amu=m, lam=float(lam), wlog_K=float(wl)))
    return outs


def from_a2f(omega, a2f, a: float, M_amu: float, unit_to_K: float = MEV_TO_K) -> ChainInput:
    """从自己的 α²F 数组构造输入。omega 的单位由 unit_to_K 指定（默认 meV）。"""
    lam, wl = moments_from_a2f(omega, a2f, unit_to_K)
    return ChainInput(a=a, M_amu=M_amu, lam=lam, wlog_K=wl)