#!/usr/bin/env python3
"""最小例子：给定 α²F 或 (λ, ω_log) ＋ 结构，算 T_c。"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from zero_supercon import from_a2f, chain

# 1) 直接给 (λ, ω_log) 与结构（MgB2 的量级，仅示例）
from zero_supercon import ChainInput
inp = ChainInput(a=3.21, M_amu=15.3, lam=0.70, wlog_K=846.0)
out = chain(inp)
print("MgB2-like:", f"μ*={out.mu_star:.4f}", f"Tc={out.Tc_K:.2f} K")

# 2) 从 α²F 算（omega 单位 meV）
omega = np.linspace(0.1, 100.0, 400)     # meV
a2f = np.exp(-0.5 * ((omega - 40.0) / 12.0) ** 2)
inp2 = from_a2f(omega, a2f, a=3.0, M_amu=50.0)
out2 = chain(inp2)
print("from a2f:", f"λ={out2.lam:.4f}", f"ω_log={out2.wlog_K:.1f} K",
      f"μ*={out2.mu_star:.4f}", f"Tc={out2.Tc_K:.2f} K")
