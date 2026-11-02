# -*- coding: utf-8 -*-
"""
zero_supercon —— 从 α²F 到临界温度的零参数链
================================================
本模块实现 Zero 库 Zb072 §11 的整条链：

    α²F(ω) ──→ (λ, ω_log) ──→ μ* ──→ T_c

其中 λ、ω_log、μ* 全部按 Zb072 §11 的式子计算，材料只需提供：
  · α²F(ω)（或等价的 λ、ω_log）
  · 结构（晶格矩阵、原子数、元素）

公开接口：
    from zero_supercon import from_jarvis, from_a2f, chain, ChainInput

数值范围（对照 JARVIS-SuperconDB 845 条）：中位比 0.884，[0.5,2] 内 92.0%。
对实测（SuperCon 74 条）：中位比 0.992，单条 ±2 倍。
"""
from .chain import (ChainInput, ChainOutput, chain, mu_hat, mu_star,
                    moments_from_a2f, allen_dynes)
from .sources import from_jarvis, from_a2f, load_jarvis

__version__ = "1.0.0"
__all__ = [
    "ChainInput", "ChainOutput", "chain", "mu_hat", "mu_star",
    "moments_from_a2f", "allen_dynes",
    "from_jarvis", "from_a2f", "load_jarvis",
]