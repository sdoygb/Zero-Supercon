# Zd021 · 形变势的通道分解：为什么 $D/t$ 不是键长指数

> 篇号：Zd021
> 状态：正式篇（推导篇）
> 前置：Zd009（$\lambda=C(D/t)^2$）、Zd014（$D/t=n$）、Zd020（材料表：$n_{\rm eff}$ 与 $n_H$）、Zb005（格点谱）、Zb072 §11（$N(0)$ 与 $\lambda$）
> 关联脚本：`src/deform_channels.py`、`src/material_table.py`
> 数据：`data/deform_channels.json`、`data/material_table.json`
> 登记日期：2026-11-22

## 摘要

Zd020 §5.5 把材料表填上后暴露一个差：结构给的 $D/t=n_H=2$／$5$，数据反解的 $n_{\rm eff}=0.34$–$0.79$——差 $4.6$–$10$ 倍。本篇说明这个差不该用"系数没算准"解释：**形变势有两个通道，进 $\lambda$ 的那个在费米面上取的是"对角＋非对角×填充"，不是键长指数本身。**

**三条命题**：

1. **矩阵元恒等式（命题 1.1）**：舞台的 hopping 随应变变时，带内矩阵元 $\langle k|\partial H/\partial\delta|k\rangle=-n\,\varepsilon_k$；加上对角项 $\Xi_d$，总形变势
$$\Xi(k)=\Xi_d-n\,\varepsilon_k .$$
2. **半满为零（命题 1.2）**：半满时 $\varepsilon_{k_F}=0$，hopping 通道在费米面上恒为零——**取 $c_N$ 的那个舞台，恰好是 hopping 通道消失的地方**。
3. **分解式（命题 2.1）**：$\lambda=C n_{\rm eff}^2$ 里的
$$n_{\rm eff}=\frac{\Xi_d}{t}+n\,\frac{\varepsilon_F}{t},$$
即"对角项＋键长指数×填充"。数据反解的 $n_{\rm eff}$ 测的是这个和，不是 $n$。

**一条判定**：$D/t=n$ 只在 $\lvert\varepsilon_F\rvert/t\simeq1$ 且 $\Xi_d$ 可忽略时成立；七个材料的 $\lambda$ 若全归 hopping 通道，要求 $\varepsilon_F/t=0.10$–$0.22$（填充 $1.03$–$1.06$，即全部贴着半满）——这个巧合不可信。故**对角通道必须存在**；它的量级由 $n_{\rm eff}$ 给（$\lvert\Xi_d\rvert/t$ 与 $n_{\rm eff}$ 同量级，当材料明显偏离半满时），但符号可与 hopping 项相反，从表里不能单独定出。

---

## §1 两个通道

舞台（Zb005 定理 1.2）：
$$\varepsilon_k=t\bigl(\lambda_k-2d\bigr),\qquad \lambda_k=4\sum_i\sin^2\frac{k_i}{2};\qquad \text{3D: } \varepsilon_k=-2t\sum_i\cos k_i .$$
应变 $\delta\equiv\delta a/a$ 进两个地方：格点能（对角）与 hopping（非对角，Zd014 的 $n=-\partial\ln t/\partial\ln a$）。

> **命题 1.1（矩阵元恒等式）.** 设 $t(\delta)=t(1-n\delta)$、格点能 $\varepsilon_i(\delta)=\varepsilon_0-\Xi_d\delta$。则带内矩阵元
> $$\langle k|\frac{\partial H}{\partial\delta}|k\rangle=\Xi_d-n\,\varepsilon_k .$$

*证明.* hopping 部分：$t(\delta)(\lambda_k-2d)=(1-n\delta)\varepsilon_k$，对 $\delta$ 求导得 $-n\varepsilon_k$；对角部分对每个 $k$ 都是 $-\Xi_d$（常数）。两项相加。$\blacksquare$

**读法.** 非对角项的矩阵元**正比于该 $k$ 处的带能**，不是常数；对角项才是 $k$ 无关的常数。

> **命题 1.2（半满为零）.** 半满（$\varepsilon_F=0$）时，hopping 通道在费米面上恒为零：$\langle k_F|\partial H_{\rm od}/\partial\delta|k_F\rangle=0$；此时 $\Xi(k_F)=\Xi_d$。

*证明.* 半满时 $\varepsilon_{k_F}=\varepsilon_F=0$，代入命题 1.1。$\blacksquare$

**读法（这条最要紧）.** Zb072 命题 11.2 取 $N(0)$ 的口径就是半满（$c_N=0.1427$）。在同一个舞台上，**hopping 通道的费米面矩阵元为零**——故 Zd009／Zd014 里那个"$D=-\beta R_0t$"不能是 $\lambda$ 的 $D$（至少不能是半满处的 $D$）。

---

## §2 分解式与数据

> **命题 2.1（$n_{\rm eff}$ 的分解）.** 取 $\lambda=N(0)\langle\Xi^2\rangle/(M\langle\omega^2\rangle)$、$N(0)=c_N/t$、$q\sim k_F=\kappa/a$，则
> $$\lambda=C\,n_{\rm eff}^2,\qquad n_{\rm eff}:=\frac{\Xi(k_F)}{t}=\frac{\Xi_d}{t}+n\,\frac{\varepsilon_F}{t},\qquad C=c_N\kappa^2g_L .$$

*证明.* $\lambda=(c_N/t)(\Xi q)^2/(t/a^2)=c_N\kappa^2(\Xi/t)^2$（Zd009 命题 1.1 的同一步），再用命题 1.1 在 $k_F$ 处取值。$\blacksquare$

**核验（填充因子）.** 三维简立方舞台上，$\varepsilon_F/t$ 与 $\lambda$ 的 hopping 份额 $(\varepsilon_F/6t)^2$：

| 填充 $n_e$ | $0.2$ | $0.4$ | $0.6$ | $0.8$ | $1.0$ | $1.2$ | $1.4$ | $1.6$ | $1.8$ |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| $\varepsilon_F/t$ | $-3.24$ | $-2.11$ | $-1.40$ | $-0.70$ | $0.00$ | $0.70$ | $1.40$ | $2.11$ | $3.24$ |
| $(\varepsilon_F/6t)^2$ | $0.292$ | $0.123$ | $0.054$ | $0.014$ | $0.000$ | $0.014$ | $0.054$ | $0.123$ | $0.292$ |

**核验（七个材料）.** 若 $\lambda$ 全来自 hopping 通道，则 $\varepsilon_F/t=n_{\rm eff}/n_H$；反推的填充 $n_e=1+2c_N(\varepsilon_F/t)$：

| 材料 | $\lambda$（表） | $n_H$ | $n_{\rm eff}$ | 若纯 hopping 的 $\varepsilon_F/t$ | 反推 $n_e$ |
|:--|--:|--:|--:|--:|--:|
| 金刚石 C | $0.259$ | $2$ | $0.43$ | $0.218$ | $1.062$ |
| Si | $0.209$ | $2$ | $0.39$ | $0.196$ | $1.056$ |
| W | $0.337$ | $5$ | $0.50$ | $0.099$ | $1.028$ |
| Ta | $0.778$ | $5$ | $0.76$ | $0.151$ | $1.043$ |
| Nb | $0.860$ | $5$ | $0.79$ | $0.159$ | $1.045$ |
| Mo | $0.695$ | $5$ | $0.71$ | $0.143$ | $1.041$ |
| Cu | $0.158$ | $2$ | $0.34$ | $0.170$ | $1.049$ |

**读法（两条）.**

1. **纯 hopping 图景要求所有材料贴在半满**（$n_e=1.03$–$1.06$）。七个材料横跨 sp 与 d、电子数从 $1$ 到 $6$，没有理由全部落在半满的 $5\%$ 以内。故**对角通道必须存在**；若材料偏离半满一个量级（$\lvert n_e-1\rvert\gtrsim0.1$），则 $\Xi_d/t$ 与 $n_{\rm eff}$ 同量级（$0.3$–$0.8$）。
2. **Cu 是最干净的反例**：s 带半满，$\varepsilon_F=0$，hopping 通道为零；而表里 $\lambda=0.158>0$（MP 表的模型值）。反解的对角形变势 $\Xi_d=n_{\rm eff}t=0.95\ \mathrm{eV}$（比 hopping 通道的 $D_{\rm hop}=5.58\ \mathrm{eV}$ 小 $5.9$ 倍，但**不为零**）。

---

## §3 对 $D/t=n$ 的读法修正

| 情形 | 用的 $D/t$ | 依据 |
|:--|:--|:--|
| 半满、且 $\Xi_d$ 可忽略 | $0$（hopping 通道为零） | 命题 1.2 |
| 偏离半满、$\Xi_d$ 可忽略 | $n\,\varepsilon_F/t$ | 命题 2.1 |
| 一般 | $\Xi_d/t+n\,\varepsilon_F/t$ | 命题 2.1 |

故 Zd014 的 $D/t=n$ 是**键长指数**，不是 $\lambda$ 直接用的那个比；两者只在"$\lvert\varepsilon_F\rvert/t\sim1$ 且 $\Xi_d$ 可忽略"时重合。七个材料的读数说明它们不重合：纯 hopping 要求 $n_e=1.03$–$1.06$，与材料实际填充不符。**$\Xi_d$ 与 $\varepsilon_F$ 不能从表里分离**——分解式给的是"两者怎么进 $\lambda$"，不是"各占多少"。

**对 Zd009 的意思**：$\lambda=C(D/t)^2$ 仍成立，但 $D$ 是**费米面处的形变势**（$\Xi(k_F)$），不是带边的 $\beta R_0t$。Zd020 §5.5 的 $D_{\rm eff}=n_{\rm eff}t$ 就是 $\Xi(k_F)$ 的读数。

---

### §3.5 $\Xi_d$ 的量级：格点能的 Coulomb 标度

> **命题 3.1（对角项的标度）.** 设格点能 $\varepsilon_0$ 的标度是 Coulomb 能 $e^2/a$（Zb042 的每格能量），则
> $$\Xi_d=-\frac{\partial\varepsilon_0}{\partial\delta}\sim c\,\frac{e^2}{a},\qquad c=O(0.1\text{–}0.6),$$
> 于是 $\Xi_d/t=c\,(e^2/a)/t$。

*证明.* $\varepsilon_0$ 是价电子在格点上的束缚能，其能量标度是 $e^2/a$（Zb042：每格能量正比于 $e^2/a$）；对 $\delta=\ln a$ 求导只差一个无量纲因子 $c$（有效价与局域场）。$\blacksquare$

**核验.** $e^2/a=14.4/a[\mathrm{\AA}]\ \mathrm{eV}$：七个材料 $3.98$–$4.58\ \mathrm{eV}$，$e^2/(at)=1.24$–$2.14$；取数据反解的 $n_{\rm eff}=0.34$–$0.79$，得 $c=0.18$–$0.61$：

| 材料 | $e^2/a$ (eV) | $e^2/(at)$ | $c=n_{\rm eff}/(e^2/at)$ |
|:--|--:|--:|--:|
| 金刚石 C | $4.04$ | $1.41$ | $0.31$ |
| Si | $2.65$ | $2.14$ | $0.18$ |
| W | $4.55$ | $1.25$ | $0.40$ |
| Ta | $4.36$ | $1.31$ | $0.58$ |
| Nb | $4.36$ | $1.30$ | $0.61$ |
| Mo | $4.58$ | $1.24$ | $0.57$ |
| Cu | $3.98$ | $1.43$ | $0.24$ |

**读法.** 对角项的标度是 **Coulomb 能**，不是带宽：$e^2/(at)\approx1.2$–$2.1$，乘 $c\sim0.2$–$0.6$ 给 $n_{\rm eff}$ 的 $0.34$–$0.79$。这正是"为什么 $n_{\rm eff}$ 是 $O(0.5)$ 而不是 $O(2$–$5)$"的答案——后者是 $t$ 的标度（$n_H$）。

**边界.** $c$ 是 C（有效价、局域场、屏蔽的残余）；$e^2/a$ 用裸 Coulomb：长程部分被电子气屏蔽，短程（on-site）部分留下。

---

## §4 边界

1. **单带、只含 hopping 与 on-site**：多带时 $\Xi$ 是矩阵（Zd019 的 $\lambda_{ab}$），本命题逐带成立、带间要另写。
2. **$\Xi_d$ 是 C**：对角形变势是材料量（格点能对应变的导数），本篇给它的**位置**（$n_{\rm eff}$ 的常数项）与量级（$0.3$–$0.8\,t$），不给第一性值。
3. **$\varepsilon_F/t$ 是 C／F**：七个材料的 $\varepsilon_F$ 未从材料算出；表里的反推值只在"纯 hopping"假设下成立。
4. **符号**：$\Xi_d$ 与 $\Xi_{\rm od}$ 可以异号（成键/反键两支，见 Zd020 §4 的界面电荷转移）；本篇只用绝对值定 $\lambda$。
5. **半满为零是模型内定理**：真实材料的 $\lambda\neq0$ 不否证它——真实材料有对角项，正是本篇的结论。

---

## 复跑

```bash
cd /Users/oygb/Downloads/Zero261117
python3 src/deform_channels.py
```

**核对数**：恒等式误差 $2.3\times10^{-10}$；半满 $\varepsilon_F/t=0.00$；纯 hopping 要求 $n_e=1.03$–$1.06$；Cu 的 $\Xi_d=0.95\ \mathrm{eV}$（$D_{\rm hop}=5.58\ \mathrm{eV}$）。

---

## 自检

| 项 | 内容 | 档位 | 依据 |
|:--|:--|:--|:--|
| 命题 1.1 | $\Xi(k)=\Xi_d-n\varepsilon_k$ | 定理（模型内） | §1 |
| 命题 1.2 | 半满时 hopping 通道为零 | 定理（模型内） | §1 |
| 命题 2.1 | $n_{\rm eff}=\Xi_d/t+n\varepsilon_F/t$ | 定理（模型内） | §2 |
| §2 表 | 填充因子与七材料反推 | 数值 | 脚本 |
| §3 表 | $D/t$ 的三种情形 | 推论 | §1–§2 |
| 命题 3.1 | $\Xi_d$ 的 Coulomb 标度 | 定理（标度）＋数值 | §3.5 |
| §4 边界 | 五条 | 边界 | — |
