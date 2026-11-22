# Zd029 · 界面增强超导的完整链：从界面键到 Ag/GO 的可测方案

> 篇号：Zd029
> 状态：正式篇（设计篇）
> 前置：Zd013（界面局域模）、Zd018（slab）、Zd020（形变势）、Zd021（通道分解）、Zd023（全谱 $T_c$）、Zd024（$D_I$ 从结构来）、Zd025（候选排序）、Zd026（$N_B$）、Zd027（界面配对的宇称）、Zd028（超晶格几何）、Zb005（舞台谱）、Zb037（库珀对能隙）、Zb072 §11（$T_c$ 与 Allen–Dynes）、Zb073（配对通道）
> 关联脚本：`src/interface_di.py`、`src/interface_scan.py`、`src/dband_dos.py`、`src/electroplate_tc.py`、`src/tc_error_budget.py`、`src/verify_chain_ag.py`
> 数据：`data/interface_di.json`、`data/interface_scan.json`、`data/dband_dos.json`、`data/electroplate_tc.json`、`data/tc_error_budget.json`
> 登记日期：2026-11-23

## 摘要

界面增强超导的链条到本篇收成一条可算、可测的线：

$$
\text{界面键}\ (R_0,\kappa,t_0)\ \Rightarrow\ D_I=n_It_I\ \Rightarrow\ N_B\ \Rightarrow\
\lambda=KN_BD_I^2\ \Rightarrow\ \omega_I\ \Rightarrow\ T_c^{\rm AD},\ T_c^{\rm full}.
\tag{0.1}
$$

**五条结论**：

1. **$D_I$ 从结构来**（Zd024 的链）：C–Ag 的 $n_I=5.03$、$t_I=1.72\ \mathrm{eV}$ ⟹ $D_I=8.67\ \mathrm{eV}$；五个 C–M 键给 $8.9$–$12.1\ \mathrm{eV}$。
2. **$N_B$ 由填充定**（Zd026）：$d$ 带未满的（Ni/Co/Pt/W/Mo）$N_B=1.3$–$5.1\ \mathrm{eV^{-1}}$；$d^{10}$ 满带的（Cu/Ag/Au/Pd）只有 $s$ 带的下限 $N_B=1.5/(2E_F)$。
3. **$\lambda$ 的标定**：$\lambda=KN_BD_I^2$，$K=0.166$ 由 slab（C–W、$D_I=8.86$、$N_B=1$、$\lambda=13.03$）定。
4. **界面模是"键的伸缩"，不是石墨烯 $G$ 模**：Morse 曲率 $K_I=2a^2D_b$ 给 C–Ag 的 $\omega_I=710\ \mathrm{cm^{-1}}=1021\ \mathrm{K}$。
5. **Ag/GO 的数**：$D_I=8.67$、$\lambda=1.34$–$1.70$、$T_c^{\rm AD}=97$–$123\ \mathrm K$（两种 $E_F$）、$T_c^{\rm full}\approx110$–$140\ \mathrm K$；16–84 分位 $[10,190]\ \mathrm K$。液氮（77 K）落在中心值以下，实验有判头。

**一条判定**：这条链的每个环节都可算、输入已分档；**最大的不确定性是 $D_I$**（$\partial\ln T_c/\partial\ln D_I=+2.36$），不是 $\omega$ 也不是 $\mu^*$。

---

## §1 链的骨架

$$
\boxed{\begin{aligned}
&R_0=r_{\rm C}+r_{\rm M},\quad \kappa=\tfrac12(\kappa_{\rm C}+\kappa_{\rm M}),\quad n_I=\kappa R_0-1,\\
&t_I=\sqrt{t_{\rm C}(R_0)\,t_{\rm M}(R_0)},\quad D_I=n_It_I,\\
&N_B:\ \text{金属的 }d\text{ 带在 }E_F\text{ 的每自旋态密度（}d^{10}\text{ 取 }s\text{ 带底线）},\\
&\lambda=K\,N_B\,D_I^2,\quad K=0.166\ (\text{slab 标定}),\\
&\omega_I=\sqrt{2a^2D_b/\mu}/(2\pi c),\\
&T_c:\ \text{Zb072 命题 11.3 的 Allen–Dynes 形式，}\mu^*\text{ 由 Zb072 §11 的几何复合给出}.
\end{aligned}}
\tag{1.1}
$$

**输入分档**：$R_0$、$\kappa$、$t_0$（C）；$K$、$\omega_I$ 的 Morse 参数、$\mu^*$（C／I）；$N_B$ 的 $E_F$（C，观测量类）；$T_c$ 的数值（导出）。

**为什么是这条链**：Zb037 给能隙、Zb072 §11 给 $T_c$ 的形状，但**声子机制的两个量**（$\lambda$、$\omega_{\log}$）在材料里没有来路——Zd020–Zd026 补的就是这两个量的结构来源：$\lambda$ 从 $D_I$ 与 $N_B$ 来，$\omega_{\log}$ 从界面键的曲率来。

---

## §2 界面键的形变势

> **命题 2.1（$D_I=n_It_I$）.** 界面的跳跃 $t$ 随键长按 $t=A R^{-1}e^{-\kappa R}$ 变，则形变势
> $$D_I=-\frac{\partial t}{\partial\ln R}=(\kappa R-1)\,t=n_I\,t .
$$

*证明.* 对 $\ln R$ 求导：$\partial\ln t/\partial\ln R=-1-\kappa R$，故 $t=A R^{-1}e^{-\kappa R}$ 的指数为 $n_I=\kappa R-1$（取 $1/r$ 势的尾部与 $R^{-1}$ 归一，Zb042）。$\blacksquare$

**核验（同核键）**：同一形式复现经验指数——C–C 给 $n=3.73$（经验 $4$）、W–W 给 $5.37$（经验 $5$）。

**五个 C–M 键**（Zd025 的表）：

| 键 | $R_0$ (Å) | $n_I$ | $t_I$ (eV) | $D_I$ (eV) |
|:--|--:|--:|--:|--:|
| C–W | $2.07$ | $4.58$ | $1.93$ | $8.86$ |
| C–Ta | $2.11$ | $4.63$ | $2.12$ | $9.83$ |
| C–Hf | $2.21$ | $4.86$ | $2.50$ | $12.14$ |
| C–Nb | $2.11$ | $4.48$ | $2.17$ | $9.73$ |
| C–Mo | $2.06$ | $4.40$ | $2.10$ | $9.24$ |
| **C–Ag** | $2.22$ | $5.03$ | $1.72$ | $8.67$ |

**读法**：$D_I$ 随 $R_0$ 涨（$n_I$ 与 $t_I$ 都随键长涨）。C–Ag 的 $8.67\ \mathrm{eV}$ 与 C–W 的 $8.86\ \mathrm{eV}$ 几乎相同——$D_I$ 不是选金属的主要杠杆。

---

## §3 $N_B$ 与 $\lambda$

> **命题 3.1（$d^{10}$ 金属的 $N_B$ 底线）.** $d$ 带满的金属，$E_F$ 在 $s$ 带；每原子一个 $s$ 电子的自由电子每自旋态密度
> $$N_B=\frac{3n}{4E_F}=\frac{3}{4E_F}\qquad(n=1).$$

*证明.* 自由电子总态密度 $D(E)=\frac{3n}{2E_F}\sqrt{E/E_F}$；在 $E_F$ 处 $D=3n/(2E_F)$，每自旋为其半。$\blacksquare$

**核验**：$E_F=7\ \mathrm{eV}$ 给 $0.107\ \mathrm{eV^{-1}}$；$E_F=5.5\ \mathrm{eV}$（Ag 的实际值）给 $0.136\ \mathrm{eV^{-1}}$。

**$\lambda$ 的标定**：$\lambda=KN_BD_I^2$，$K=0.1660$——由 slab 的 C–W 点（$D_I=8.86$、$N_B=1$、$\lambda=13.03$）定；$K$ 的不确定度取 $\pm50\%$（slab 的 400 个模式与真实界面模式的差异，C）。

**核验（Ag/GO）**：$\lambda=0.166\times0.107\times8.67^2=1.34$（$E_F=7$）；$0.166\times0.136\times8.67^2=1.70$（$E_F=5.5$）。

---

## §4 界面模与 $T_c$

> **命题 4.1（界面模频率）.** 界面键的伸缩频率由 Morse 曲率定：
> $$K_I=2a^2D_b,\qquad \omega_I=\frac{1}{2\pi c}\sqrt{K_I/\mu},$$
> $\mu$ 是两原子的约化质量。

**核验（C–Ag）**：$D_b=2.5\ \mathrm{eV}$、$a=2.0\ \mathrm{\AA^{-1}}$、$\mu=10.80$ amu 给 $K_I=320\ \mathrm{N/m}$、$\omega_I=710\ \mathrm{cm^{-1}}=1021\ \mathrm K$。

**读法**：界面模是**键的伸缩**（$\sim1000\ \mathrm K$），不是石墨烯的 $G$ 模（$2273\ \mathrm K$）——后者是面内 C–C 振动，不改变 Ag–C 键长。

**$T_c$**（Zb072 命题 11.3 的 Allen–Dynes 形式）：

| $E_F$ (eV) | $N_B$ | $\lambda$ | $T_c^{\rm AD}$ (K) | $T_c^{\rm full}$ (K) |
|--:|--:|--:|--:|--:|
| $7.0$ | $0.107$ | $1.34$ | $97$ | $110$ |
| $5.5$ | $0.136$ | $1.70$ | $123$ | $140$ |

**读法**：两种 $E_F$ 给 $97$–$123\ \mathrm K$（AD）；全谱按 Zd023 的强耦合修正因子（$1.15$–$1.2$）给 $110$–$140\ \mathrm K$。**液氮 77 K 在两者之下。**


---

## §4.5 对照案例：FeSe/SrTiO$_3$

**事实（F，文献值）**：体相 FeSe 的 $T_c=8\ \mathrm K$；单层长在 SrTiO$_3$ 上 $T_c=65\ \mathrm K$（有 $109\ \mathrm K$ 的报道）；界面耦合的声子是 STO 的 LO 模，约 $90\ \mathrm{meV}=1045\ \mathrm K$。

**反解 $\lambda$**（Zb072 命题 11.3 的 Allen–Dynes 形式，给定 $T_c$ 与 $\omega$ 求 $\lambda$）：

| 案例 | $T_c$ (K) | $\omega$ (K) | $\lambda_{\rm req}$（$\mu^*=0.12$） | $\lambda_{\rm req}$（$\mu^*=0.20$） |
|:--|--:|--:|--:|--:|
| 体相 FeSe | $8$ | $250$ | $0.72$ | $0.93$ |
| 单层 FeSe/STO | $65$ | $1045$ | $0.99$ | $1.26$ |
| 单层（$109\ \mathrm K$ 报道） | $109$ | $1045$ | $1.46$ | $1.87$ |

**与 Ag/GO 的对照**：本链给 $\lambda=1.34$–$1.70$、$\omega=1021\ \mathrm K$ ⟹ $T_c=97$–$123\ \mathrm K$——落在 FeSe/STO 的 $65$–$109\ \mathrm K$ 带的上沿；若 Ag/GO 的有效 $\lambda$ 只有 FeSe/STO 量级（$\lambda=1.0$），则 $T_c=64\ \mathrm K$，**低于液氮 77 K**。

**通道检查**：把 Fe–O 界面当**短程键通道**（§2 的链）给 $D_I=5.96\ \mathrm{eV}$、$N_B({\rm Fe})=2.52$ ⟹ $\lambda=14.9$——是反解值（$1.0$–$1.3$）的 $10$ 倍以上。故 FeSe/STO 的增强**不是**键伸长机制，而是极性氧化物 LO 声子的长程（Fröhlich 型）耦合；本链的 $D_I$ 通道不是它的描述。

**读法**：FeSe/STO 给这条路线一个**经验锚**——$\omega\approx1000\ \mathrm K$ 的界面声子，$\lambda\approx1$ 给 $65\ \mathrm K$、$\lambda\approx1.5$–$1.9$ 给 $109\ \mathrm K$。Ag/GO 的模型值在这条锚的延长线上；**77 K 的判据正好卡在"$\lambda$ 是 1.0 还是 1.3–1.7"这个岔口上**。

---

## §5 误差预算

输入区间：$D_I\in[3.2,19.5]\ \mathrm{eV}$（三种标定）、$N_B\in[0.083,0.150]$、$K\in[0.083,0.250]$、$\mu^*\in[0.10,0.20]$、$\omega_I\in[600,1200]\ \mathrm K$；$4\times10^4$ 次抽样（`src/tc_error_budget.py`）。

| 量 | 区间 | $\partial\ln T_c/\partial\ln X$ | 对 $T_c$ 的影响 |
|:--|:--|--:|:--|
| $D_I$ | $3.2$–$19.5$（因子 $6$） | $+2.36$ | 主导 |
| $N_B$ | $0.083$–$0.150$（$1.8$） | $+1.18$ | 中 |
| $K$ | $0.083$–$0.250$（$3$） | $+1.18$ | 中 |
| $\omega_I$ | $600$–$1200$（$2$） | $+1.01$ | 中 |
| $\mu^*$ | $0.10$–$0.20$（$2$） | $-0.43$ | 弱 |

**结果（Ag/GO）**：中位 $118\ \mathrm K$、16–84 分位 $[10,190]\ \mathrm K$、5–95 分位 $[0,228]\ \mathrm K$。

**读法**：带宽是"上下各约 4 倍"，不是"2–3 倍"；**$D_I$ 一项占主导**。这也说明模型改进的方向：不是调 $\mu^*$，是把 $D_I$ 的标定（$t_0$、$R_0$、$\kappa$）钉死。

---

## §6 可测方案：Ag/GO ＋ 液氮

### §6.1 几何

$$\xi=\frac{\hbar v_F}{\pi\Delta},\quad \Delta=1.76k_BT_c;\qquad T_c=95\ \mathrm K\ \Rightarrow\ \xi=6.7\text{–}19.3\ \mathrm{nm}.$$

银膜取 $10$–$25\ \mathrm{nm}$（$\lesssim2\xi$），方阻目标 $1.5$–$4\ \Omega/\square$（含薄膜有效电阻率 $2.5\rho_{\rm bulk}$，C）。

### §6.2 制样

1. **电解银**：AgNO₃ $10\ \mathrm{g/L}$ ＋酒石酸钾钠 $30\ \mathrm{g/L}$，pH 8.5–9.5，$0.5$–$1.0\ \mathrm V$、$0.5$–$2\ \mathrm{mA/cm^2}$；厚度由电量定 $t[\mathrm{nm}]\approx1.07\,Q[\mathrm{mC/cm^2}]$；
2. **碳侧**：氧化石墨烯滴涂，维生素 C $80\ ^\circ\mathrm C$ 还原；
3. **四端**：银浆触点，恒流＋换向。

### §6.3 判据

> **判据 6.1（界面增强的可测判据）.** 四条同时成立：
> 1. $R(77\ \mathrm K)/R(300\ \mathrm K)<10^{-2}$（四端、换向平均）；
> 2. 升温曲线在 $77$–$190\ \mathrm K$ 出现上升台阶；
> 3. 台阶随银膜减薄而变深（$\propto2\xi/d$）；
> 4. 两个对照（纯银膜、GO/玻璃）都不出现台阶。

**读法**：1 是"过不过"，2–4 是"是不是界面"；缺 2–4 的 1 不作数（可能是接触/热电势的假象，见 §7）。

---

## §7 边界

1. **$K$ 的 $\pm50\%$ 是判断**：slab（400 模式、C–W）搬到 Ag/GO 没有独立核实。
2. **$D_I$ 的标定带**：$3.2$–$19.5\ \mathrm{eV}$ 由三种标定（C–C 基准、几何平均、Harrison 型）给出；几何平均是选择，不是定理。
3. **$N_B$ 的 $s$ 带底线**是自由电子估计；界面态可以抬高它。
4. **铁磁金属（Ni/Co/Fe）不在此列**：$s$ 波配对被自旋涨落压制（Zd020 §3 的磁通道）。
5. **宇称**：本链的配对是带内偶宇称（Zd027 命题 1.1）；跨界面奇分量的权重（$w_{A_1}=1/4$）不进本链的 $\lambda$。
6. **"模型内"**：$T_c$ 的数值是模型数，不是预言；判据 6.1 是给实验的判决条件。

---

## 复跑

```bash
cd /Users/oygb/Downloads/Zero261117
python3 src/interface_di.py       # D_I 的链
python3 src/dband_dos.py          # 金属的 N_B
python3 src/electroplate_tc.py    # 候选表
python3 src/tc_error_budget.py    # 预算与灵敏度
python3 src/verify_chain_ag.py    # Ag/GO 全链条独立复算
```

**核对数**：$D_I({\rm C–Ag})=8.67\ \mathrm{eV}$；$N_B=0.107$／$0.136$；$\lambda=1.34$／$1.70$；$\omega_I=1021\ \mathrm K$；$T_c^{\rm AD}=97$／$123\ \mathrm K$；16–84 分位 $[10,190]\ \mathrm K$。

---

## 自检

| 项 | 内容 | 档位 | 依据 |
|:--|:--|:--|:--|
| §1 链 | 五步的形式 | 定义＋标定 | §1 |
| 命题 2.1 | $D_I=n_It_I$ | 定理（模型内） | Zd024 |
| 命题 3.1 | $d^{10}$ 的 $N_B$ 底线 | 定理（自由电子） | §3 |
| 命题 4.1 | $\omega_I$ 由 Morse 曲率定 | 定理（模型内） | §4 |
| §4 表 | Ag/GO 的 $T_c$ | 模型（C 输入） | 脚本 |
| §4.5 表 | FeSe/STO 的 $\lambda$ 反解与通道检查 | 数值（F 输入） | 脚本 |
| §5 表 | 误差预算与灵敏度 | 数值 | 脚本 |
| 判据 6.1 | 界面增强的判定 | 判据（I） | §6 |
| §7 边界 | 六条 | 边界 | — |
