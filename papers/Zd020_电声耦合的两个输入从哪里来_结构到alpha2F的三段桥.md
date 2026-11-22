# Zd020 · 电声耦合的两个输入从哪里来：结构 → $\alpha^2F$ 的三段桥

> 篇号：Zd020
> 状态：正式篇（推导篇）
> 前置：Zd002（四个材料数）、Zd009（$\lambda=C(D/t)^2$）、Zd014（$n$ 与 $g_L$）、Zd018（界面 $g$ 的标定）、Zb005（格点谱）、Zb012（图 Laplacian）、Zb020（Pauli 与交换）、Zb042（$1/r$）、Zb072 §11（$N(0)$ 与 $T_c$）
> 关联脚本：`src/alpha2f_bridge.py`、`src/dos_constant.py`、`src/material_table.py`
> 数据：`data/alpha2f_bridge.json`、`data/dos_constant.json`、`data/material_table.json`
> 登记日期：2026-11-22

## 摘要

Zd002／Zd003 把 $T_c$ 写成 $(\lambda,\omega_{\log})$ 的函数；Zd007 §2 把每个材料的 $\alpha^2F$ 记为"算不出（F）"。本篇把 $\alpha^2F$ 的**两个矩**各接一段结构来路：

$$
\boxed{\;\begin{aligned}
&\omega_{\log}\ \leftarrow\ \text{力常数 }\Phi\ \leftarrow\ \text{键能 }E_b(R)\ \leftarrow\ \text{键};\\
&N(0)\ \leftarrow\ \text{能带}\ \leftarrow\ \text{舞台谱（图／轨道）};\\
&\lambda\ \leftarrow\ g=\frac{D}{a}\sqrt{\frac{\hbar}{2M\omega}}\ \leftarrow\ \text{形变势 }D=\beta R_0t\ \leftarrow\ \text{键}.
\end{aligned}\;}
\tag{0.1}
$$

**四条结论**：

1. **力常数是键能的二阶导**：Morse 键给 $\Phi=2D_ba^2$、$\omega=a\sqrt{2D_b/\mu}$。C–C 的核验：$D_b=3.6\ \mathrm{eV}$、$\omega=1332\ \mathrm{cm^{-1}}$ 给 $a=2.33\ \mathrm{\AA^{-1}}$、$\Phi=627\ \mathrm{N/m}$。
2. **声子谱是动力学矩阵的本征值**：$D_{ij}(\mathbf q)=\frac1{\sqrt{m_im_j}}\sum\Phi_{ij}e^{i\mathbf q\cdot\mathbf R}$；两支进 $\omega_{\log}$，**声学支拉低它**。1D 双原子链核验：声学 $0$–$1.0$、光学 $1.414$–$1.732$、$\omega_{\log}=0.841$。
3. **$N(0)$ 由能带定，不是自由参数**：多带 $H_{\alpha\beta}(\mathbf k)$ 对角化给 partial DOS；单带舞台把它压成 $N(0)=c_N/t$，三维半满 $c_N=0.1427$（单自旋；两法复算）。
4. **$g$ 必须带长度**：$g=(D/a)\sqrt{\hbar/(2M\omega)}=Iu_0$（$I=Dq$）。$q=k_F$ 给 Zd009 的 $C(D/t)^2$；$q\sim1/a$ 给界面局域模的保守值。两种口径的 $g$ 差 $\kappa=(3\pi^2)^{1/3}=3.09$，进 $\lambda$ 差 $\kappa^2=9.57$。

**一条判定**：$\alpha^2F$ 从"算不出"变成"形式有了、值要数据"：$\Phi$、$g$、$D$ 的具体值要键能数据／第一性，本库给的是它们与结构量的**标度与形式**。

---

## §1 键 → 力常数 → 声子频率

> **命题 1.1（Morse 键的力常数）.** 取两中心键 $E_b(R)=D_b\bigl(1-e^{-a(R-R_0)}\bigr)^2$，则
> $$\Phi:=E_b''(R_0)=2D_ba^2,\qquad \omega=\sqrt{\frac{\Phi}{\mu}}=a\sqrt{\frac{2D_b}{\mu}} .$$

*证明.* $E_b'=2D_ba(1-e^{-a(R-R_0)})e^{-a(R-R_0)}$，$E_b''=2D_ba^2\bigl(2e^{-2a(R-R_0)}-e^{-a(R-R_0)}\bigr)$，在 $R_0$ 处 $E_b''=2D_ba^2$。$\omega=\sqrt{\Phi/\mu}$ 是单键谐振子。$\blacksquare$

**核验（C–C）.** $D_b=3.6\ \mathrm{eV}$、$\mu=m_C/2=6\,m_p$、$\nu=1332\ \mathrm{cm^{-1}}$ 给 $a=2.33\ \mathrm{\AA^{-1}}$、$\Phi=627\ \mathrm{N/m}$（经典共价键 $400$–$600\ \mathrm{N/m}$ 的量级）。

**读法.** 键深 $D_b$ 与键程尺度 $1/a$ 两个数定声子频率；力常数不是独立的第三个材料数。

---

## §2 力常数 → 声子谱 → $\omega_{\log}$

> **命题 2.1（动力学矩阵）.** 晶格的声子谱是对角化
> $$D_{ij}(\mathbf q)=\frac{1}{\sqrt{m_im_j}}\sum_{\boldsymbol\delta}\Phi_{ij}(\boldsymbol\delta)\,e^{i\mathbf q\cdot\boldsymbol\delta}$$
> 的本征值；声学支 $\omega(\mathbf q\to0)\to0$，光学支 $\omega(\mathbf q\to0)\to\omega_0=\sqrt{\Phi/\mu}$。

**核验（1D 双原子链，$m_1=1$、$m_2=2$、$K=1$）.** 声学支 $0$–$1.0$、光学支 $1.414$–$1.732$、两支之间隙 $0.414$；按 $\omega_{\log}=\exp\langle\ln\omega\rangle$ 对两支取几何平均给 $0.841$。

**读法.** $\omega_{\log}$ 不是光学支的 $\omega_0$：**声学支（$\omega\to0$）把它拉低**。这就是 Zd016 的候选类要付的第二笔代价（第一笔是 $\lambda=C(D/t)^2$ 里 $t$ 的压制）：轻元素只把光学支抬高，$\omega_{\log}$ 仍由低端定。

---

## §3 能带 → $N(0)$

> **命题 3.1（多带态密度）.** 每个顶点带一组轨道 $\alpha$，
> $$H_{\alpha\beta}(\mathbf k)=\varepsilon_\alpha\delta_{\alpha\beta}+\sum_{\boldsymbol\delta}t_{\alpha\beta}\,e^{i\mathbf k\cdot\boldsymbol\delta},$$
> 对角化给各带 $E_n(\mathbf k)$；$N(0)=\sum_n\int_{\rm FS}\frac{ds}{\lvert\nabla E_n\rvert}$，按轨道权重分解（partial DOS）。

**核验（单带、三维半满）.** 舞台谱 $\varepsilon_k=t\lambda_k-2dt$（Zb005 定理 1.2）给 $N(0)=c_N/t$，$c_N=0.1427$（单自旋，归一到 1；Monte Carlo $5\times10^9$ 样给 $0.14266$，$512^3$ 网格给 $0.1439$——网格有 $+0.9\%$ 离散偏置）。

**读法.** $N(0)$ 不是第四个材料数：它由舞台谱定（Zb005），材料只选谱。单带值是 $c_N/t$；多带只按轨道权重重新分配，总量仍由谱的带宽定（Zb072 命题 11.2）。

---

## §4 形变势 → $g$ → $\lambda$

> **命题 4.1（$g$ 的结构形式）.** 取形变势 $D=\beta R_0t$（$\beta=-\partial\ln t/\partial R$）、矩阵元 $I=Dq$ 与声子零点位移 $u_0=\sqrt{\hbar/(2M\omega)}$。在 $\hbar=1$ 的约定下
> $$g=I\,u_0=D\,q\,u_0=\frac{D}{a}\sqrt{\frac{\hbar}{2M\omega}}\quad(q\sim1/a),$$
> 两种写法给同一个耦合：
> $$\lambda=N(0)\frac{\langle I^2\rangle}{M\langle\omega^2\rangle}=2N(0)\frac{\langle g^2\rangle}{\omega}=C\Bigl(\frac{D}{t}\Bigr)^2,\qquad C=c_N\kappa^2g_L .$$

*证明.* 把 $I=Dq$、$N(0)=c_N/t$、$q=\kappa/a$、$M\omega^2=t/a^2$ 代入第一式即得 $C(D/t)^2$（Zd009 命题 1.1／Zd014 命题 2.1）。第二式用 $u_0^2=\hbar/(2M\omega)$ 与 $\hbar=1$：$2I^2u_0^2/\omega=I^2/(M\omega^2)$。$g$ 的表达式里 $D/a$ 是"单位位移的能量"：$D$ 是**应变**导数（能量），除以键长才是每单位长度的力。$\blacksquare$

**核验（C–C）.** $\beta=1.5\ \mathrm{\AA^{-1}}$、$R_0=1.54\ \mathrm{\AA}$、$t=2.5\ \mathrm{eV}$ 给 $D=5.78\ \mathrm{eV}$；$M=12\,m_p$、$\omega=1332\ \mathrm{cm^{-1}}$ 给 $u_0=0.0325\ \mathrm{\AA}$。于是

| $q$ 的口径 | $g$ (eV) | 用在哪 |
|:--|--:|:--|
| $q=1/a$，$a=3\ \mathrm{\AA}$ | $0.063$ | 界面局域模（Zd018 的保守口径） |
| $q=1/R_0$，$a=1.54\ \mathrm{\AA}$ | $0.122$ | 键长的上端 |
| $q=k_F=\kappa/a$ | $0.193$ | 体相声子（进 $\lambda$ 的口径） |

**读法（口径决定量级）.** 同一个 $D$ 下，三种 $q$ 给 $0.06$–$0.19\ \mathrm{eV}$；$q=k_F$ 与 $q=1/a$ 差 $\kappa=3.09$，进 $\lambda$ 差 $\kappa^2=9.57$。**"$g$ 是多少"必须先说 $q$ 是哪一支**；Zd018 的 $0.023$–$0.114\ \mathrm{eV}$ 正是 $q\sim1/a$ 一支，与本篇的 $0.063$–$0.122$ 同一量级。

**闭合核验.** $D/t=2.31$ 时：$q=1/a$ 给 $\lambda=c_N(D/t)^2=0.76$；$q=k_F$ 给 $\lambda=C(D/t)^2=7.29$。两者之比恰是 $\kappa^2=9.57$。

---

## §5 三段合起来：Zd007 第 1 件缺口的档位

| 段 | 输入 | 输出 | 档位 |
|:--|:--|:--|:--|
| 键 → 力常数 | $D_b$、$a$、$\mu$ | $\Phi$、$\omega_0$ | **C（键能数据）** |
| 力常数 → 声子谱 | $\Phi_{ij}$、$m_i$ | 两支、$\omega_{\log}$ | 形式（导出）；$\Phi$＝C |
| 能带 → $N(0)$ | 轨道、跳跃 $t_{\alpha\beta}$ | $N(0)$、partial DOS | 形式；$t_{\alpha\beta}$＝C |
| 形变势 → $g$ | $\beta$、$R_0$、$t$、$q$ | $g$、$\lambda$ | 标度（导出）；值＝C |

**读法.** 四段的**形式**都在本库（Zb005 的谱、Zb012 的图 Laplacian、Zb020 的 Pauli、Zb042 的 $1/r$）；每段的**值**都停在 C。故 Zd007 §2 的第 1 件（$\alpha^2F$）从"本库没有电声顶点"变成"形式已给、值要数据"：短板从"写不下"移到"拿不到数"。档位仍是 F（观测量），但不再是无路可走。

### §5.5 填表：七个材料与一个界面

把 §5 表里的 C 格用真实材料填上（脚本 `src/material_table.py`、数据 `data/material_table.json`）。输入分档：$a$、$M$、$R_0$、Harrison 指数是 C；$\omega_{\log}$、$\lambda$ 取库表（F，MP 分档）。

| 材料 | $a$ (Å) | $M$ ($m_p$) | $\omega_{\log}$ (K) | $\lambda$（表） | $t_{\rm bare}$ (eV) | $n_H$ | $n_{\rm eff}$ | $D_H$ (eV) | $D_{\rm eff}$ (eV) | $u_0$ (Å) | $\Phi$ (N/m) | $g_H$ (eV) | $g_{\rm eff}$ (eV) |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| 金刚石 C | $3.567$ | $12.0$ | $1336.5$ | $0.259$ | $2.87$ | $2$ | $0.43$ | $5.73$ | $1.25$ | $0.0389$ | $305$ | $0.063$ | $0.014$ |
| Si | $5.431$ | $28.1$ | $403.1$ | $0.209$ | $1.24$ | $2$ | $0.39$ | $2.47$ | $0.48$ | $0.0463$ | $65$ | $0.021$ | $0.004$ |
| W | $3.165$ | $183.8$ | $171.1$ | $0.337$ | $3.64$ | $5$ | $0.50$ | $18.20$ | $1.81$ | $0.0278$ | $77$ | $0.160$ | $0.016$ |
| Ta | $3.306$ | $180.9$ | $129.6$ | $0.778$ | $3.34$ | $5$ | $0.76$ | $16.68$ | $2.52$ | $0.0322$ | $43$ | $0.162$ | $0.025$ |
| Nb | $3.301$ | $92.9$ | $175.4$ | $0.860$ | $3.35$ | $5$ | $0.79$ | $16.73$ | $2.66$ | $0.0386$ | $41$ | $0.196$ | $0.031$ |
| Mo | $3.147$ | $96.0$ | $230.9$ | $0.695$ | $3.68$ | $5$ | $0.71$ | $18.41$ | $2.63$ | $0.0331$ | $73$ | $0.194$ | $0.028$ |
| Cu | $3.615$ | $63.5$ | $223.3$ | $0.158$ | $2.79$ | $2$ | $0.34$ | $5.58$ | $0.95$ | $0.0413$ | $45$ | $0.064$ | $0.011$ |

（$t_{\rm bare}=\hbar^2k_F^2/2m_e$，$k_F=\kappa/a$，Zd002 命题 2.2 的裸带；$n_H$＝Harrison 键长指数（sp 2、d 5，C）；$n_{\rm eff}=\sqrt{\lambda/C}$ 是数据反解的 $D/t$；$D_H=n_Ht$、$D_{\rm eff}=n_{\rm eff}t$；$\Phi=\mu\omega_{\log}^2$，$\mu=M/2$。）

**界面（C–W）**：$R_0=2.07\ \mathrm{\AA}$、$M=12$、$\omega_I=1500\ \mathrm K$、$a=3\ \mathrm{\AA}$；由 Zd018 的 $D_I=5\ \mathrm{eV}$ 与 $\beta=1.5\ \mathrm{\AA^{-1}}$ 反解 $t_I=D_I/(\beta R_0)=1.61\ \mathrm{eV}$、$n_I=3.1$。于是 $u_0=0.0367\ \mathrm{\AA}$、$g_I=0.061\ \mathrm{eV}$；到 $g_{\rm req}=0.43$／$0.63$ 要 $M_{\rm eff}=49$／$106$。

**读法（三条）.**

1. **结构估计与数据反解差 $4.6$–$10$ 倍**：$n_H=2$／$5$ 是键长指数，数据反解的 $n_{\rm eff}=0.34$–$0.79$。故进 $\lambda$ 的形变势**不是带边位移 $\beta R_0t$，而是费米面处的有效值**——差的这个因子就是 §5 表里"值＝C"的具体内容（填充与几何）。同一件事在 $\lambda$ 上放大成 $21$–$101$ 倍：$\lambda_H=5.5$–$34$ 而表值是 $0.16$–$0.86$。
2. **$\Phi$ 的两个口径**：由 $\omega_{\log}$ 读出 $41$–$305\ \mathrm{N/m}$；由金刚石光学模（$1332\ \mathrm{cm^{-1}}$）读出 $627\ \mathrm{N/m}$。$\omega_{\log}$ 是两支的几何平均，不是键的伸缩频率——这是 §2 的声学支代价在数上的样子。
3. **$g$ 的自然值是 $0.004$–$0.031\ \mathrm{eV}$**，比 $g_{\rm req}$ 小 $14$–$150$ 倍；界面键把 $D$ 抬到 $5\ \mathrm{eV}$ 才给 $g_I=0.061\ \mathrm{eV}$，多模需求 $M_{\rm eff}=49$–$106$——与 Zd018 的"$30$–$150$"一致。填表后这条判定更紧：体相的 $g_{\rm eff}$ 比界面的 $0.061\ \mathrm{eV}$ 还小 $2$–$15$ 倍。

**边界.** $t$ 用裸带是**上限**：真实 d 带 $t$ 更小（W 的 d 带 $\sim1\ \mathrm{eV}$），$g$ 会更小。$\lambda$ 取 MP 分档模型值（F／模型）；Nb、Ta 与 JARVIS 的 $\alpha^2F$ 派生值差 $3\%$–$11\%$。$n_H$ 是 Harrison 指数（C），$n_{\rm eff}$ 是反解（把 $\lambda$ 全归到 $n$）。界面行的 $D_I$ 取 Zd018 的 $5\ \mathrm{eV}$（C），不是从结构算的。

---

## §6 边界

1. **$\Phi$、$g$、$D$ 的具体值是 C**：本篇给的是与结构量的标度与形式；Morse 参数、$\beta$、轨道跳跃都要键能／第一性数据。
2. **$q$ 的口径必须写明**：$q=k_F$ 进 $\lambda$（体相），$q\sim1/a$ 进界面局域模的 $g$（Zd018）；进 $\lambda$ 差 $\kappa^2$，不能互换。
3. **Migdal 条件**：$\Gamma=\omega_D/\varepsilon_F\ll1$ 才单圈；适用区与失效区见 Zd016 §3.5（轻网络 $\Gamma\ll1$，窄带／平带才失效）。
4. **电子关联不在本篇**：多带 $H_{\alpha\beta}$ 是单粒子图像；关联泛函是开放（背靠 Zb012 的 Laplacian 与 Zb020 的交换）。
5. **两支的 $\omega_{\log}$**：声学支的 $\omega\to0$ 使 $\ln\omega$ 发散，实际的 $\omega_{\log}$ 依赖低频截断（温度、边界、无序），截断是材料量（F）。

---

## 复跑

```bash
cd /Users/oygb/Downloads/Zero261117
python3 src/alpha2f_bridge.py     # 力常数/双原子链/形变势与 g
python3 src/dos_constant.py       # c_N 的两法复算（约 40 s）
```

**核对数**：$a=2.33\ \mathrm{\AA^{-1}}$、$\Phi=627\ \mathrm{N/m}$；声学 $0$–$1.0$、光学 $1.414$–$1.732$、$\omega_{\log}=0.841$；$D=5.78\ \mathrm{eV}$、$u_0=0.0325\ \mathrm{\AA}$；$g=0.063/0.122/0.193\ \mathrm{eV}$；$\lambda=0.76$ 与 $7.29$（比 $9.57$）。

---

## 自检

| 项 | 内容 | 档位 | 依据 |
|:--|:--|:--|:--|
| 命题 1.1 | $\Phi=2D_ba^2$、$\omega=a\sqrt{2D_b/\mu}$ | 定理（模型内） | §1 |
| 命题 2.1 | 动力学矩阵与两支 | 定理（模型内） | §2 |
| 命题 3.1 | 多带 $N(0)$ | 定理（模型内） | §3 |
| 命题 4.1 | $g=(D/a)\sqrt{\hbar/2M\omega}$、$\lambda=C(D/t)^2$ | 定理（模型内） | §4 |
| §4 表 | 三种 $q$ 的 $g$ | 数值 | 脚本 |
| §5 表 | 四段的档位 | 汇总 | §1–§4 |
| §5.5 表 | 七个材料＋界面的 $\Phi, D, g$ | 数值（值＝C） | 脚本 |
| §6 边界 | 五条 | 边界 | — |
