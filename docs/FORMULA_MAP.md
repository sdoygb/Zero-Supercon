# 公式与出处对照

本工具实现的每个式子都来自 Zero 库的一篇。下表是完整索引。

## 链条

```
α²F(ω)
  │  式 (11.13)：λ = 2∫α²F/ω dω，ω_log = exp[(2/λ)∫(α²F/ω)lnω dω]
  ▼
(λ, ω_log)
  │  命题 11.11：ε_F/ω_D = (a k_F)√(M/2m_e)
  │  命题 11.12：ê² = e²/(a t) = 2(a/a_B)/(a k_F)²
  │  命题 11.9 ：μ = N(0)·4πê²/(⟨q²⟩ + 4πê²N(0))
  │  命题 11.13：μ* = μ/(1 + μ ln(ε_F/ω_D))
  ▼
μ*
  │  式 (11.11)：T_c = (ω_log/1.2)exp[−1.04(1+λ)/(λ−μ*(1+0.62λ))]
  ▼
T_c
```

## 逐项

| 式子 | Zero 库位置 | 本工具位置 | 常数 |
|:--|:--|:--|:--|
| $\lambda=2\int\alpha^2F/\omega\,d\omega$ | Zb072 式 (11.13) | `chain.moments_from_a2f` | — |
| $\omega_{\log}=\exp[(2/\lambda)\int(\alpha^2F/\omega)\ln\omega\,d\omega]$ | Zb072 式 (11.13) | 同上 | — |
| $k_F=(3\pi^2)^{1/3}$（半满） | Zb005 谱；Zb072 §11.10 | `chain.AK_F` | — |
| $N(0)=0.3713$（3D 立方格半满） | Zb072 §11.7 | `chain.N0_LATTICE` | — |
| $\langle q^2\rangle_{\rm FS}=18.97$ | Zb072 §11.7 | `chain.Q2_LATTICE` | — |
| $\hat e^2=2(a/a_B)/(ak_F)^2$ | Zb072 命题 11.12 | `chain.chain` | $a_B$ |
| $\mu$ | Zb072 命题 11.9 | `chain.mu_hat` | — |
| $\ln(\varepsilon_F/\omega_D)$ | Zb072 命题 11.11 | `chain.chain` | $m_e$ |
| $\mu^*$ | Zb072 命题 11.13 | `chain.mu_star` | — |
| $T_c$ | Zb072 式 (11.11) | `chain.allen_dynes` | — |
| $e^2=4\pi\alpha$ | Zb036 定理 4.1/5.3 | （只经 $\hat e^2$ 引入） | $\alpha^{-1}=137.036$ |
| $m_p/m_e=1836.153$ | Zb065 定理 1.1 | （只经 $M$ 引入） | — |

## 哪些是导出、哪些是输入

| 量 | 档位 |
|:--|:--|
| $\alpha$、$m_p/m_e$、$k_F$、$N(0)$、$\langle q^2\rangle$ | 导出（库内纯数） |
| $\mu$、$\mu^*$、$\varepsilon_F/\omega_D$ | 导出 |
| $\lambda$、$\omega_{\log}$ | 由 α²F 导出；α²F 本身是外部输入 |
| $a$、$M$ | 结构/成分输入 |
| 1.134、1.04、1.2、0.62 | Allen–Dynes 形式的标准常数（不是本库的常数） |

## 与标准公式的关系

标准 McMillan/Allen–Dynes 把 $\mu^*$ 当**自由参数**（常用 0.1–0.13）。
本工具把它**算出来**：$\mu^*=\mu(a,M)/(1+\mu\ln(\varepsilon_F/\omega_D))$，其中 $\mu$ 由 Thomas–Fermi 屏蔽与费米面平均给出。

在 $a=3$ Å、$M=m_p$ 时本工具给 $\mu^*=0.112$；中等金属给 0.087–0.092——都落在标准区间内，但**不是拟合来的**。
