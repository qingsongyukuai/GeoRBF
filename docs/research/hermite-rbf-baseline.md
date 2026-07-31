# GeoRBF 广义 Hermite RBF 数学与数值基线

状态：研究结论，供后续 spec、tickets 与 v1.0.0 路线使用  
研究票据：[GeoRBF #3](https://github.com/qingsongyukuai/GeoRBF/issues/3)  
Surfe 参照点：[`290dbe0`](https://github.com/MichaelHillier/surfe/tree/290dbe0ab344f4258a4935f05cad0f153f0f69a4)  
范围：三维、小中规模、稠密直接法优先；不选择 Rust 数值 crate，不实现产品代码，不覆盖 VTK、Qt 或可视化

## 结论摘要

GeoRBF 不应重写 Surfe 的五套模型内核。v1 的共同数学核心应是：

1. 用一个有限线性泛函表示所有观测原语：

   \[
   Lf=\sum_{a=1}^{s}\left(\alpha_a f(x_a)+\beta_a^\mathsf T\nabla f(x_a)\right).
   \]

   evaluation、gradient component、directional derivative、同层位点之间的 evaluation difference，都只是 \((x_a,\alpha_a,\beta_a)\) 的不同组合。组合后的系统统一由

   \[
   K_{ij}=L_i^xL_j^y k(x,y),\qquad
   P_{i\gamma}=L_i p_\gamma
   \]

   装配，而不是按“界面模型”“层序模型”“属性模型”等分支手写成对 block。

2. 对严格正定（SPD）核，泛函必须线性独立且核必须能承受所请求的导数。对 \(m\) 阶条件正定（CPD）核，必须同时声明符号、CPD 阶数、对应多项式空间 \(\Pi_{m-1}\)、side condition、unisolvence 与 gauge。核名称或 Surfe 的枚举值不能替代这些证明义务。

3. 纯等式插值的自然系统是对称鞍点系统

   \[
   \begin{bmatrix}K&P\\P^\mathsf T&0\end{bmatrix}
   \begin{bmatrix}c\\d\end{bmatrix}
   =
   \begin{bmatrix}y\\0\end{bmatrix}.
   \]

   对 CPD 核，\(K\) 只需在 \(\ker(P^\mathsf T)\) 上正定，通常并非全空间 PSD。因此不能直接把原始 \(K\) 当作通用 QP Hessian。可行的数学路线是：等式问题使用稳定的对称不定分解；凸优化问题先以 null-space \(c=Zu\) 约化，使 \(Z^\mathsf TKZ\succ0\)，或使用另有证明的投影核。Surfe 的 modified kernel 和坐标 epsilon 扰动不能当作一般证明。

4. 硬/软约束必须是逐观测的属性，而不是整个模型的 mode：

   - \(Lf=y\)、\(\ell\le Lf\le u\) 是仿射硬约束；
   - 带权平方残差是凸二次软约束；
   - 带非负 slack 的二次或 \(L_1\) 罚项仍可落为 QP/锥问题；
   - 有方向且极性已知的 normal 角锥可落为 SOCP；
   - 未知极性、精确物理厚度、\(\|\nabla f\|=1\)、仅以 tangent 定义的双侧角带，一般是非凸或不可识别的，不能伪装成线性 box。

5. v1 可以承诺“小中规模、稠密、可诊断”，但接口应把装配、问题形式和后端能力分开。求解结果必须报告 rank/gauge、缩放、残差、约束违背、迭代/终止状态与不可行性；仅检查 `isfinite` 不够。

6. Surfe 可保留的思想是 generalized interpolation、层位 difference、normal/tangent 的导数语义、restricted-range 的变分方向，以及把各向异性写成坐标度量。必须拒绝的是五套重复装配/求值、按 mode 切 solver、未经 unisolvence/gauge 验证的 modified-kernel 实现、依赖 \(\|\nabla f\|\approx1\) 或 \(\approx2\) 的角度/厚度换算、经验 eigenvalue 修补、自制无完整状态的 QP，以及把手工 VTK 程序当测试。

下面以三个标签明确区分结论：

- **数学上必须**：不满足就不能声称问题有定义、可解或凸。
- **Surfe 历史选择**：固定提交中的事实，只作为参考或失败边界。
- **GeoRBF 待决定**：需要后续 spec 作产品选择，但本研究没有擅自定案。

## 1. 证据与适用边界

### 1.1 一手数学资料

核心依据如下：

- Hillier 等人的地质建模论文明确把多种地质数据表示为线性泛函，并讨论 inequality、数据驱动 anisotropy 与额外方向泛函：[Hillier et al., 2014, DOI 10.1007/s11004-014-9540-3](https://doi.org/10.1007/s11004-014-9540-3)。其 [erratum, DOI 10.1007/s11004-014-9554-x](https://doi.org/10.1007/s11004-014-9554-x) 修正了 contact、gradient 与 tangent 泛函公式，并重申泛函集合须线性独立；应连同正文一起读。
- Wu 给出一泛函/一 center 参数化下、由有限阶导数 evaluations 构成的 Hermite–Birkhoff RBF 表示、广义 block 系统、rank 条件与多项式再现：[Wu, 1992, author-accessible PDF](https://pages.mtu.edu/~struther/Courses/OLD/5630/Refs/RBFs/HermiteBirkhoofInter_RadBasisFun.pdf)，[DOI 10.1007/BF02836101](https://doi.org/10.1007/BF02836101)。
- Micchelli 给出 CPD、辅助多项式与 side condition 的经典可解性基础，并证明 Gaussian、inverse multiquadric、幂函数/多二次型相关结果：[Micchelli, 1986, PDF](https://pages.stat.wisc.edu/~wahba/stat860public/pdf1/micchelli.interpolation.86.pdf)，[DOI 10.1007/BF01893414](https://doi.org/10.1007/BF01893414)。
- Madych–Nelson 提供 CPD 插值的变分/半 Hilbert 空间与误差框架：[Madych & Nelson, 1990, DOI 10.1090/S0025-5718-1990-0993931-7](https://doi.org/10.1090/S0025-5718-1990-0993931-7)。
- Narcowich–Ward 把 generalized Hermite 数据推广到紧支撑分布与线性 side conditions，并给出 well-poised/invertibility 与逆矩阵范数结果：[Narcowich & Ward, 1994, DOI 10.1090/S0025-5718-1994-1254147-6](https://doi.org/10.1090/S0025-5718-1994-1254147-6)，[AMS PDF](https://www.ams.org/mcom/1994-63-208/S0025-5718-1994-1254147-6/S0025-5718-1994-1254147-6.pdf)，[作者出版目录](https://people.tamu.edu/~f-narcowich/pubs.html)。
- 紧支撑、给定空间维数下的正定 Wendland 函数依据原始论文：[Wendland, 1995, DOI 10.1007/BF02123482](https://doi.org/10.1007/BF02123482)。
- restricted-range 的原始工作将问题描述为：在满足有限约束的函数中最小化平滑度二次泛函：[Beatson et al., 2004, publication record and abstract](https://www.researchgate.net/publication/265264671_Surface_Reconstruction_via_Smoothest_Restricted_Range_Approximation)。可访问材料没有给出足够细节来证明 Surfe 的具体 modified-kernel 实现。以 unisolvent anchors 为 CPD seminorm 增加 gauge/norm、从而得到 SPD reproducing kernel 的构造，可由 [Beatson, Light & Billings, 2000, DOI 10.1137/S1064827599361771](https://doi.org/10.1137/S1064827599361771) 核对；这仍不为 Surfe 的 anchor 选择和 epsilon 扰动背书。
- SOCP 的标准形式与凸性依据作者页及原始论文：[Lobo et al., 1998](https://web.stanford.edu/~boyd/papers/socp.html)。凸优化的标准定义、KKT 与缩放背景见作者提供的 [Boyd & Vandenberghe, 2004](https://www.seas.ucla.edu/~vandenbe/cvxbook.html)。
- 对称不定系统需要专门的对称 pivoting，而非把一般 LU 当作数值契约：[Bunch & Kaufman, 1977, DOI 10.1090/S0025-5718-1977-0428694-0](https://doi.org/10.1090/S0025-5718-1977-0428694-0)。

### 1.2 证据限制

本研究不根据二手 RBF 表或 Surfe 类名补全 CPD 表。尤其是符号约定很关键：某些常见 radial formula 本身是 conditionally negative definite，乘以 \(-1\) 后才是 CPD；同一插值空间在纯线性求解中可以通过系数符号吸收，但在“把矩阵当凸目标 Hessian”时不能忽略。

因此，核表只把由上述一手资料与公式直接支持的结论标为“已准入候选”。其余项明确保留为证明缺口。v1 spec 若要增加这些核，必须附带维数、符号、阶数和原点 jet 的独立证明或权威出处。

## 2. 一个核心：有限线性泛函

### 2.1 规范表示

在三维域 \(\Omega\subset\mathbb R^3\) 上，v1 需要的观测泛函可限制为

\[
L_i f
=
\sum_{a=1}^{s_i}
\left(
\alpha_{ia} f(x_{ia})
+
\beta_{ia}^\mathsf T\nabla f(x_{ia})
\right),
\]

其中 \(s_i\) 有限，\(\alpha_{ia}\in\mathbb R\)，\(\beta_{ia}\in\mathbb R^3\)。这是“值与一阶导数的有限线性组合”，足以覆盖本次范围，不需要公开二阶导数观测。

各语义的 lowering 如下：

| 上层语义 | 线性泛函 |
|---|---|
| 点值 | \(\delta_x f=f(x)\) |
| 第 \(q\) 个 gradient component | \(\partial_q\delta_x f=e_q^\mathsf T\nabla f(x)\) |
| directional derivative | \(D_v\delta_x f=v^\mathsf T\nabla f(x)\) |
| 同界面/同层位 | \((\delta_x-\delta_{x_0})f=0\) |
| 层位增量 | \((\delta_{x_a}-\delta_{x_b})f=\Delta\) 或带上下界 |
| 多个观测的线性关系 | 对上述原子的有限加权和 |

Narcowich–Ward 以紧支撑 distributions 处理最一般的同点、多点和导数组合；Wu 给出一泛函/一 center 参数化下有限光滑 radial kernel 的具体 Hermite–Birkhoff block、rank 与 regularity 条件。二者共同支持本抽象。不能把 Wu 的参数化无条件外推为所有同点多泛函情形；也不能用 Narcowich–Ward 的 \(C^\infty\) kernel 假设替代有限光滑幂核的原点检查。

**数学上必须**

- lowering 后的每个 \(L_i\) 必须保留来源 ID、单位、缩放与方向约定；否则无法做诊断或把结果映射回输入。
- 方向向量是否要求归一化必须由语义层声明。数学核心只计算给定 \(\beta\)，不会偷偷归一化。
- 必须在装配前检测重复、零泛函和线性依赖的明显情形；“观测记录不同”不等于“泛函线性独立”。

**GeoRBF 待决定**

- 上层 `Interface`、`Normal`、`Tangent`、`Order` 等公开类型的最终命名与序列化格式；
- 是否允许任意多个 support points 的公开泛函，还是 v1 公开少量安全构造器、内部使用通用表示；
- observation covariance 是只支持对角权重，还是从一开始允许小块相关协方差。

### 2.2 广义 representer 与装配

令 \(k(x,y)\) 为对称核，\(\{p_\gamma\}_{\gamma=1}^M\) 是所需多项式空间的基。候选解写为

\[
f(x)
=
\sum_{j=1}^{N} c_j L_j^y k(x,y)
+
\sum_{\gamma=1}^{M}d_\gamma p_\gamma(x).
\]

应用 \(L_i^x\) 得

\[
K_{ij}=L_i^xL_j^y k(x,y),\qquad
P_{i\gamma}=L_i p_\gamma.
\]

因此 evaluation–evaluation、evaluation–gradient、gradient–gradient、difference–difference 等 block 都从同一个双线性 pairing 生成。一个 kernel jet 只需提供：

\[
k(x,y),\quad
\nabla_x k,\quad
\nabla_y k,\quad
\nabla_x\nabla_y^\mathsf T k,
\]

上层系数完成 contraction。比如

\[
D_u^xD_v^y k=u^\mathsf T
(\nabla_x\nabla_y^\mathsf T k)v,
\]

而 difference–difference 自动展开成四项，不需要专用模型代码。

这也统一了求值。任意查询泛函 \(Q\) 的预测为

\[
Qf
=
\sum_j c_j Q^xL_j^y k(x,y)
+
\sum_\gamma d_\gamma Qp_\gamma.
\]

scalar evaluation 与 gradient evaluation 只是不同的 \(Q\)，不能再各模型复制一套循环。

### 2.3 对称性与符号

当核对称、泛函为实线性且 mixed derivatives 存在时，

\[
K_{ij}=K_{ji}.
\]

实现应由一个 canonical pairing 同时写入对称项，或独立计算后用误差阈值验证；不应依赖五套手写 block 恰好符号一致。对平移不变径向核，\(\nabla_x k=-\nabla_y k\)，但这只是可测试的解析恒等式，不应通过散落的 `dx_p1` / `dx_p2` 分支维护。

## 3. 可解性：光滑性、PD/CPD、多项式与 gauge

### 3.1 导数 regularity

因为系统同时在行、列应用一阶泛函，collocated derivative–derivative 元素需要

\[
\partial_{x_a}\partial_{y_b}k(x,y)
\]

在 \(x=y\) 有定义且极限一致。对 v1，一条简单、偏保守但可执行的准入规则是：核 \(k\) 作为 \((x,y)\) 的函数至少有连续 mixed derivatives 到 \((1,1)\)，并为 \(r=0\) 和 compact-support cutoff 提供显式解析极限。

**数学上必须**

- “函数值可算”不代表 Hermite Gram 可算；
- 原点不能用任意 epsilon 替代解析极限；
- compact-support 核还必须在 cutoff 两侧满足声明的连续性；
- anisotropic coordinate transform 后必须对物理坐标正确应用 chain rule。

Surfe 的 `R` 核正是反例：它的 value 返回 \(r\)，但所有一阶与二阶方法都直接 `throw -666`，[`basis.cpp` L1300–L1369](https://github.com/MichaelHillier/surfe/blob/290dbe0ab344f4258a4935f05cad0f153f0f69a4/surfe_lib/basis.cpp#L1300-L1369)。它不能进入涵盖 normal/tangent 的 GeoRBF v1 kernel registry。

### 3.2 SPD 情形

若 \(k\) 对所考虑的泛函集合严格正定，则 \(K\succ0\)，不需要数学上强制的多项式块。Wu 对 Gaussian 与 inverse multiquadric 指出它们是 \(C^\infty\)，Fourier transform 为正且衰减，在 rank 条件下 generalized interpolation 唯一存在。

可选 trend 与 CPD 所需多项式必须分开：

- **required polynomial**：由 CPD 阶数决定，不允许用户随意关闭；
- **optional trend**：产品建模选择，需要另外说明可识别性与 regularization。

把二者合并成 `poly_term: bool` 会制造无效配置。

### 3.3 CPD 情形

采用如下明确约定：

> \(k\) 是严格 \(m\) 阶 CPD，若对任何非零 \(c\) 满足
> \(\sum_i c_i L_i p=0\)（所有 \(p\in\Pi_{m-1}\)），都有
> \(c^\mathsf TKc>0\)。

对 generalized functionals，点值 Vandermonde 被

\[
P_{i\gamma}=L_i p_\gamma
\]

取代。应求解

\[
\begin{bmatrix}
K&P\\
P^\mathsf T&0
\end{bmatrix}
\begin{bmatrix}c\\d\end{bmatrix}
=
\begin{bmatrix}y\\0\end{bmatrix},
\qquad
P^\mathsf Tc=0.
\]

在一般 functional IR 中，Narcowich–Ward 的 well-poised 条件要求数据 distributions 在 admissible space 上独立，并正确消去 polynomial modes；其有限维表现包括
\(\operatorname{rank}(P)=\dim\Pi_{m-1}\) 且 \(K\) 在 side-condition 子空间上严格正定。Wu 给出其一泛函/一 center 特例中的 rank 条件和 finite-smoothness 条件，Micchelli 给出点值 CPD 的经典基础。它们共同导出 v1 的直接 preflight：

1. kernel metadata 决定 \(m\) 与多项式基；
2. 在归一化坐标中构造 \(P\)；
3. 用 rank-revealing 分解估计 rank 与 null space；
4. 若 \(P\) 不满列秩，报告缺失的 polynomial modes；
5. 若满秩，再求解鞍点系统或约化问题。

在三维，\(\dim\Pi_{m-1}=\binom{m+2}{3}\)。\(P\) 满列秩仍不能取代 functional independence 检查；重复或依赖的同点泛函必须另外诊断。

### 3.4 gauge 不是“加一点 epsilon”

仅含 derivative 或 evaluation difference 的泛函会消灭常数：

\[
L_i1=0.
\]

因此常数 offset 不可识别。更高阶的组合也可能消灭线性或更高多项式。此时有三种不同情况：

- 核的 CPD side condition 所需 \(P\) 满秩：正常唯一；
- 目标本来只定义到某个 gauge：应显式声明 quotient，并要求一个 gauge condition 才输出绝对 field value；
- 数据本应识别该 mode 但没有识别：输入退化，必须报错。

允许的 gauge 例子包括固定一个 field value、固定多项式系数，或显式零均值约束；选择哪种会影响输出语义，必须由 spec 决定。修改一个输入坐标使四点“不共面”不是 gauge。

Surfe 的 difference 模型知道常数项消失，因此手工把多项式项数减一，[`lajaunie.cpp` L194–L208](https://github.com/MichaelHillier/surfe/blob/290dbe0ab344f4258a4935f05cad0f153f0f69a4/surfe_lib/lajaunie.cpp#L194-L208)。这个观察是对的，但 GeoRBF 应从 \(P\) 的 rank/null space 得出，而不是由模型名触发。

### 3.5 CPD 与凸目标

CPD 只保证

\[
c^\mathsf TKc>0\quad\text{for }c\ne0,\ P^\mathsf Tc=0,
\]

不保证 \(K\succeq0\) 于全空间。标准 convex QP 要求 Hessian 在自由变量空间上 PSD；SOCP 也不能接受任意不定二次目标。

一个严格可行的约化是取 \(Z\) 的列为 \(\ker(P^\mathsf T)\) 的正交基，令 \(c=Zu\)。则

\[
H=Z^\mathsf TKZ\succ0.
\]

对所有观测与查询建立线性 map 后，硬等式、box、slack 和角锥都可在 \(u\) 及其他辅助变量上表达。另一条路线是使用经过证明的 projected/modified kernel。两条路线都必须处理 polynomial/gauge 变量，且证明与具体 functional set 一致。

**GeoRBF 待决定**

- 约化采用显式 \(Z\)、range-space 还是一个已证明的投影核；
- 多项式系数如何与优化变量共同消元；
- rank tolerance、坐标缩放与 null-space 更新策略；
- SPD 与 CPD 是否共享同一个 problem-form 层，后端只看到标准化后的 linear/QP/SOCP。

这些是 solver spec 的问题，不是 Rust crate 选择。

## 4. Kernel registry：证明随实现走

### 4.1 必需 metadata

每个核定义必须是一个不可拆开的契约：

- 精确公式及所有参数的合法域；
- 支持的空间维数；
- `SPD` 或“带明确符号的 `CPD(m)`”；
- CPD 对应的 \(\Pi_{m-1}\)；
- 支持的 functional derivative order；
- \(r=0\) 的 value/gradient/mixed-Hessian 极限；
- 若紧支撑，cutoff 处的连续性和左右极限；
- shape/length scale 的单位与预缩放语义；
- global linear anisotropy 下的 chain rule；
- 引用的定理或针对实现公式的证明测试。

只有公式相同、参数化相同、维数条件相同，文献结论才可转移。

### 4.2 v1 研究准入表

下表不是承诺的产品列表，而是下一步 spec 的证据门槛。

| 公式族 | 本研究可确认 | 对 v1 的结论 |
|---|---|---|
| Gaussian \(e^{-(\varepsilon r)^2}\), \(\varepsilon>0\) | 严格 PD，\(C^\infty\)；Wu 与 Micchelli 均覆盖 | **已准入候选**；仍需定义 shape 单位、缩放和稳定参数范围 |
| inverse multiquadric \((c^2+r^2)^{-\beta}\), \(c,\beta>0\) | 严格 PD，\(C^\infty\)；Micchelli §4、Wu Example 4.1 | **已准入候选**；Surfe 的参数写法是 `shape + r²`，必须重命名或明确定义 |
| cubic \(r^3\) | 由 Micchelli 的幂函数结果，按本文符号是 `CPD(2)`，需 \(\Pi_1\)；在原点有本范围所需 mixed second derivative 极限 | **已准入候选**；须以三维实现 jet 与 null-space positivity 测试锁定 |
| Wendland \(C^2\) 的三维实例 | Wendland 原论文给出按维数的紧支撑严格 PD 构造；Surfe 公式为 \((1-r/\rho)^4_+(1+4r/\rho)\) | **可准入候选**，前提是 spec 固定为 \(d\le3\)、\(\rho>0\)，并验证原点/cutoff jet |
| distance \(r\) | raw 正号是 conditionally negative，\(-r\) 为 `CPD(1)`；但 \(r=0\) 不支持 collocated first-derivative Hermite | **拒绝进入本范围**；改变符号也不能补足所需光滑性，Surfe 的 `throw -666` 与此一致 |
| multiquadric \((c^2+r^2)^{1/2}\) | 按本文约定，raw 正号是 conditionally negative；\(-\sqrt{c^2+r^2}\) 才是 `CPD(1)`，对应 \(\Pi_0\)。\(c>0\) 时 \(C^\infty\) | **正号拒绝作凸 Hessian**；负号版本完成 metadata/约化测试后才可准入 |
| power multiquadric \((c^2+r^2)^{3/2}\) | Micchelli 的一般幂函数结果给出 raw 正号 `CPD(2)`，对应 \(\Pi_1\)；\(c>0\) 时 \(C^\infty\) | **可准入候选**；不能继承 Surfe 任意 polynomial setting |
| \(r^4\log r\) | 按本文约定，\(-r^4\log r\) 是 `CPD(3)`，对应 \(\Pi_2\)；raw 正号为相反符号。它足以支持本范围的一阶泛函 | **Surfe 正号拒绝作凸 Hessian**；负号/P2 版本完成证明测试后才可准入 |
| Surfe `MaternC4` 公式 | 源码公式类似缩放后的 Matérn \(5/2\)，但类名、smoothness 命名与参数语义不充分 | **未准入**；需要以精确公式重新命名并给 PD/jet 证明 |

这样保留了“系统化核与 CPD 阶数”，同时没有伪造完整表。一个核若没有 metadata，就不是“暂时少一个优化”，而是无效配置。

### 4.3 Surfe kernel 事实

固定提交中：

- cubic 实现为 \(r^3\)，[`basis.cpp` L261–L348](https://github.com/MichaelHillier/surfe/blob/290dbe0ab344f4258a4935f05cad0f153f0f69a4/surfe_lib/basis.cpp#L261-L348)；
- MQ 实现为 \(\sqrt{\text{shape}+r^2}\)，[`basis.cpp` L636–L710](https://github.com/MichaelHillier/surfe/blob/290dbe0ab344f4258a4935f05cad0f153f0f69a4/surfe_lib/basis.cpp#L636-L710)；
- 所谓 TPS 实现为 \(r^4\log r\)，[`basis.cpp` L825–L936](https://github.com/MichaelHillier/surfe/blob/290dbe0ab344f4258a4935f05cad0f153f0f69a4/surfe_lib/basis.cpp#L825-L936)；
- Wendland C2 公式和 cutoff 分支位于 [`basis.cpp` L1442–L1570](https://github.com/MichaelHillier/surfe/blob/290dbe0ab344f4258a4935f05cad0f153f0f69a4/surfe_lib/basis.cpp#L1442-L1570)；
- kernel 接口枚举了 point/planar/tangent 的所有成对方法，[`basis.h` L49–L81](https://github.com/MichaelHillier/surfe/blob/290dbe0ab344f4258a4935f05cad0f153f0f69a4/surfe_lib/basis.h#L49-L81)，而通用 contraction 已经隐藏在部分 forwarding 中，[`basis.cpp` L155–L258](https://github.com/MichaelHillier/surfe/blob/290dbe0ab344f4258a4935f05cad0f153f0f69a4/surfe_lib/basis.cpp#L155-L258)。

GeoRBF 应把后者提升为 first-class functional/jet，而不是继续扩展 pairwise virtual-method 矩阵。

## 5. 约束 lowering 与凸性边界

### 5.1 统一变量视图

一旦 \(f\) 由有限系数 \(\theta\) 表示，每个线性泛函都有

\[
Lf=a^\mathsf T\theta.
\]

约束与损失由此分类。分类必须发生在数学问题层，不应由模型类决定 solver。

| 语义 | 规范形式 | 问题类别 | 备注 |
|---|---|---|---|
| 硬点值/导数/差值 | \(a^\mathsf T\theta=y\) | linear equality | 可进入鞍点或凸问题 |
| 硬范围/顺序 | \(\ell\le a^\mathsf T\theta\le u\) | LP/QP/SOCP 中的 affine constraint | 可诊断 infeasible |
| 软等式 | \(\frac{w}{2}(a^\mathsf T\theta-y)^2\) | convex QP，\(w\ge0\) | \(w\) 有单位 |
| 软 box，二次 slack | \(\ell-s^-\le a^\mathsf T\theta\le u+s^+\), \(s^\pm\ge0\) | convex QP | 罚项需明示 |
| 软 box，\(L_1\) slack | 同上，加 \(\lambda(s^++s^-)\) | LP/QP/SOCP | 便于 sparse violation |
| normal 角锥，极性已知 | SOC 约束 | SOCP | 需排除零梯度 |
| 未知极性 normal | 两个相反锥的并 | 非凸/离散 | 不能当一个 SOC |
| 精确 \(\|\nabla f\|=s\) | norm equality | 非凸 | 不能当 QP |

### 5.2 逐点硬软混合

同一观测批次中可以有：

- 界面 A 的两个点为 hard difference；
- 第三个点为 soft difference；
- 某 normal 的方向为 hard、模长为 soft；
- 某顺序为 hard lower bound；
- 另一个顺序允许有带罚 slack。

核心不需要“linear mode / restricted-range mode”。每条 lowered constraint 携带：

- relation：equality、lower、upper、interval、cone；
- enforcement：hard 或 soft；
- loss/penalty 与非负参数；
- measurement units / covariance；
- source ID 与 group ID。

**数学上必须**

- 权重不得为负；
- hard 与无限权重不是同一个数值实现；
- 标准化必须同步变换目标、bounds、报告残差与用户单位；
- 相互矛盾的 hard constraints 必须返回 infeasible 或 rank conflict，不得落成 NaN 或“最接近”而不告知。

### 5.3 界面与层位关系

同一 horizon 的最稳健线性语义是选 anchor \(x_0\)，对其他点施加

\[
f(x_i)-f(x_0)=0.
\]

这与给所有点指定未知共同 level 等价，但差值表示直接暴露 gauge。层位顺序可以写成

\[
f(x_\text{above})-f(x_\text{below})\ge\Delta_f,
\]

其中 \(\Delta_f\) 是 **field-value 单位** 的间隔。若 \(\Delta_f=0\)，只表达顺序；若非零，它仍是 affine constraint，但不是自动等于物理距离。

Surfe 的 Lajaunie 路径正是用 anchor pairs，[`lajaunie.cpp` L138–L150](https://github.com/MichaelHillier/surfe/blob/290dbe0ab344f4258a4935f05cad0f153f0f69a4/surfe_lib/lajaunie.cpp#L138-L150)，pair–pair kernel 手工展开四项，[`lajaunie.cpp` L688–L730](https://github.com/MichaelHillier/surfe/blob/290dbe0ab344f4258a4935f05cad0f153f0f69a4/surfe_lib/lajaunie.cpp#L688-L730)。GeoRBF 保留 difference 语义，但让通用 \(L_i^xL_j^yk\) 自动展开。

### 5.4 厚度为何不可直接等于 field increment

沿曲线 \(\gamma\) 的 field difference 是

\[
f(\gamma(1))-f(\gamma(0))
=
\int_0^1\nabla f(\gamma(t))^\mathsf T\gamma'(t)\,dt.
\]

只有在方向与积分路径明确、且 gradient scale 已校准时，field increment 才能转换为长度。尤其：

- 仅有 orientation 通常只约束方向，不确定 \(\|\nabla f\|\)；
- field 的仿射重标度会改变 increment，而不改变零水平面几何；
- “最小厚度 = 两点 field difference 下界”实际是 field-unit 间隔，不是无条件的米数。

精确几何厚度往往涉及未知交点、积分路径或 \(1/\|\nabla f\|\)，会成为 nonlinear/nonconvex 问题。v1 若坚持保持 convex，应把产品能力命名为“stratigraphic field separation”，除非另有 scale calibration spec。

Surfe 自己在代码注释中承认 thickness bound 依赖 \(\|\nabla f\|\approx1\)，且“现实中并不成立”，[`stratigraphic_surfaces.cpp` L420–L445](https://github.com/MichaelHillier/surfe/blob/290dbe0ab344f4258a4935f05cad0f153f0f69a4/surfe_lib/stratigraphic_surfaces.cpp#L420-L445)。这应成为 GeoRBF 明确拒绝的失败边界。

**GeoRBF 待决定**

- v1 是否只提供 field separation；
- 是否提供有单位的 gradient magnitude calibration；
- 物理厚度是否明确推迟到 nonlinear backend；
- horizon polarity 与“above/below”的坐标/地质约定。

### 5.5 Normal

需要区分四种产品语义：

1. **完整 gradient observation**

   \[
   \nabla f(x)=g_0.
   \]

   是三个 affine equalities，同时约束方向和局部 scale。

2. **仅方向、极性已知、精确**

   可用两个独立 tangent 向量 \(t_1,t_2\perp n\)：

   \[
   t_1^\mathsf T\nabla f=0,\quad
   t_2^\mathsf T\nabla f=0,
   \]

   再以 \(n^\mathsf T\nabla f\ge s_{\min}>0\) 排除零梯度。这仍是 affine。

3. **角度容差、极性已知**

   令 \(g=\nabla f(x)\)，单位 normal 为 \(n\)，\(0\le\theta<\pi/2\)。锥

   \[
   \|(I-nn^\mathsf T)g\|_2
   \le
   \tan\theta\,(n^\mathsf Tg),
   \qquad n^\mathsf Tg\ge0
   \]

   是 convex SOC；但 \(g=0\) 也满足，所以若“normal observation”要求非零，还需线性下界 \(n^\mathsf Tg\ge s_{\min}>0\) 或其他 scale 信息。

4. **极性未知**

   可行集是以 \(n\) 与 \(-n\) 为轴的两个锥之并，通常非凸。v1 不应无声选择一个极性。

Surfe 把 normal 拆为三个 gradient components，这一基本语义可以保留：装配 block 见 [`single_surface.cpp` L1007–L1027](https://github.com/MichaelHillier/surfe/blob/290dbe0ab344f4258a4935f05cad0f153f0f69a4/surfe_lib/single_surface.cpp#L1007-L1027)。但其 uncertainty 是把各分量独立转为 lower/upper box，[`modelling_input.cpp` L201–L269](https://github.com/MichaelHillier/surfe/blob/290dbe0ab344f4258a4935f05cad0f153f0f69a4/surfe_lib/modelling_input.cpp#L201-L269)。分量 box 不是旋转不变的角锥，也无法严谨表达 polar angle。

### 5.6 Tangent

精确 tangent observation 的自然语义是

\[
t^\mathsf T\nabla f(x)=0,
\]

即一个 directional derivative equality。Surfe 的 tangent kernel forwarding 也做这个 contraction，[`basis.cpp` L169–L175](https://github.com/MichaelHillier/surfe/blob/290dbe0ab344f4258a4935f05cad0f153f0f69a4/surfe_lib/basis.cpp#L169-L175)。

“gradient 与 tangent 近似正交，角误差为 \(\theta\)”若写成

\[
|t^\mathsf Tg|\le\sin\theta\,\|g\|_2,
\]

一般不是 convex constraint：它允许一个双锥/互补形状，右侧还是变量 norm。若改成

\[
|t^\mathsf Tg|\le\tau
\]

则为 affine interval，但 \(\tau\) 是 directional-derivative 单位，不是无需 scale 的角度。

Surfe 的 `setAngleBounds` 把角度转换为 inner-product 上下界，并显式假设 gradient norm 约为 2，[`modelling_input.h` L238–L252](https://github.com/MichaelHillier/surfe/blob/290dbe0ab344f4258a4935f05cad0f153f0f69a4/surfe_lib/modelling_input.h#L238-L252)。GeoRBF 不应保留这一隐含比例。

**GeoRBF 待决定**

- v1 tangent uncertainty 是 field-unit interval，还是推迟角度语义；
- normal 是否默认带 polarity；
- normal vector 的 magnitude 是数据还是输入时归一化；
- \(s_{\min}\) 的来源、单位与默认行为；
- normal angle SOCP 是否属于 v1 必做后端能力。

### 5.7 各向异性

最安全的 v1 形式是一个用户明确给出的全局可逆线性变换 \(A\)，定义

\[
k_A(x,y)=\phi(\|A(x-y)\|_2).
\]

若 \(M=A^\mathsf TA\succ0\)，这就是全局 SPD metric。因为它等价于在变换坐标 \(z=Ax\) 中使用原核，原核的 PD/CPD 性质在可逆坐标变换下保留，多项式空间也通过坐标变换对应。

物理坐标导数必须遵循

\[
\nabla_x k_A=A^\mathsf\nabla_z k,\qquad
\nabla_x\nabla_y^\mathsf T k_A
=
A^\mathsf H_{zw} A.
\]

Surfe 的 global anisotropy 从 normal 外积和 eigenvalues 构造 transform，并把小 eigenvalue 截到 \(10^{-4}\)，[`basis.cpp` L75–L152](https://github.com/MichaelHillier/surfe/blob/290dbe0ab344f4258a4935f05cad0f153f0f69a4/surfe_lib/basis.cpp#L75-L152)。这是经验估计器，不是 metric 输入契约，也没有 uncertainty/identifiability 诊断。

实验分支只能作为未来方向：

- [`TensorInterpolation` commit `d0beb9f`](https://github.com/MichaelHillier/surfe/tree/d0beb9fc697d6c0bcd1e38fba31d4615a5267674) 定义局部 orientation tensors、邻域 eigen-analysis 与 SPD matrix 插值；
- 其中 local anisotropy 调用在执行路径被注释，[`modeling_methods.cpp` L437–L451](https://github.com/MichaelHillier/surfe/blob/d0beb9fc697d6c0bcd1e38fba31d4615a5267674/surfe_lib/modeling_methods.cpp#L437-L451)；
- 邻域 eigenvalues 还包含绝对值、floor 与经验增量修补，[`anisotropy_input.cpp` L227–L320](https://github.com/MichaelHillier/surfe/blob/d0beb9fc697d6c0bcd1e38fba31d4615a5267674/surfe_lib/anisotropy_input.cpp#L227-L320)。

这说明“局部 SPD metric field”值得预留接口，但不证明把任意 \(A(x)\) 塞入 \(\phi(\|A(x)(x-y)\|)\) 后仍对称或正定。v1 只应承诺全局 \(A\)；局部 nonstationary kernel 必须有独立 PD 构造。

## 6. 数值问题与求解契约

### 6.1 三类问题，不按模型分 solver

装配后应按数学形式选择：

1. **等式插值/回归的对称系统**
   - SPD：Cholesky 类分解；
   - CPD augmented system：带 pivoting 的对称不定分解；
   - rank-deficient：rank-revealing QR/SVD 或明确 gauge 处理。

2. **带 affine bounds 的凸 QP**
   - PSD/PD 的约化 Hessian；
   - hard/soft equalities、bounds、slacks；
   - 可靠 primal/dual status 与 infeasibility 诊断。

3. **带 normal angle cone 的 SOCP**
   - affine equalities/bounds 加二阶锥；
   - 可靠 primal/dual residual、gap、infeasibility certificate。

未来 sparse/matrix-free backend 应实现同一标准化 problem form；不应要求 domain model 知道矩阵存储。

### 6.2 前处理与缩放

RBF 系统对坐标尺度、shape parameter、不同 observation units 和 polynomial columns 敏感。v1 spec 至少要规定：

- 坐标平移与 characteristic length 缩放；
- functional row scaling：value 与 derivative 的单位不同；
- polynomial basis 在 normalized coordinates 中计算；
- objective/constraint scaling 可逆，报告恢复到用户单位；
- kernel 参数是在物理坐标还是 normalized coordinates 中解释；
- condition estimate 或 factorization diagnostics；
- tolerance 相对问题 scale 定义，不写死绝对常数。

### 6.3 必须返回的诊断

成功结果至少包括：

- problem class 与后端 capability；
- \(N\)、polynomial dimension、估计 rank/nullity；
- coordinate、row、objective scaling；
- primal residual、hard constraint 最大违背；
- soft residual 分组统计；
- 若优化，dual residual、duality gap、iteration count、termination reason；
- conditioning/pivot 警告；
- gauge 选择；
- 每个失败或高违背项对应的 source ID。

失败应是结构化状态：invalid kernel/functional pairing、rank deficient、infeasible、unbounded、numerical failure、iteration limit。不能把所有失败折叠成 `false`。

### 6.4 Surfe 求解边界

Surfe 的线性求解只调用 `partialPivLu()`，随后只检查权重是否 finite，[`matrix_solver.cpp` L47–L58](https://github.com/MichaelHillier/surfe/blob/290dbe0ab344f4258a4935f05cad0f153f0f69a4/surfe_lib/matrix_solver.cpp#L47-L58)。这忽略了系统对称结构、rank/gauge、残差和 conditioning。

其 QP 是自制 predictor-corrector/LOQO 变体；代码中可见显式对角逆与由一般 LU 解 KKT 的路径，例如 [`math_methods.cpp` L94–L448](https://github.com/MichaelHillier/surfe/blob/290dbe0ab344f4258a4935f05cad0f153f0f69a4/math_lib/math_methods.cpp#L94-L448) 和 [`math_methods.cpp` L539–L975](https://github.com/MichaelHillier/surfe/blob/290dbe0ab344f4258a4935f05cad0f153f0f69a4/math_lib/math_methods.cpp#L539-L975)。这不能成为 GeoRBF 的稳定性或终止状态基线。

GeoRBF 本研究只规定 solver contract 和数学 form，不选择库。

## 7. Surfe 固定提交审计

### 7.1 五套重复内核

固定提交有五个独立 `get_interpolation_matrix`：

- [`single_surface.cpp` L924](https://github.com/MichaelHillier/surfe/blob/290dbe0ab344f4258a4935f05cad0f153f0f69a4/surfe_lib/single_surface.cpp#L924)
- [`lajaunie.cpp` L672](https://github.com/MichaelHillier/surfe/blob/290dbe0ab344f4258a4935f05cad0f153f0f69a4/surfe_lib/lajaunie.cpp#L672)
- [`stratigraphic_surfaces.cpp` L482](https://github.com/MichaelHillier/surfe/blob/290dbe0ab344f4258a4935f05cad0f153f0f69a4/surfe_lib/stratigraphic_surfaces.cpp#L482)
- [`continuous_property.cpp` L475](https://github.com/MichaelHillier/surfe/blob/290dbe0ab344f4258a4935f05cad0f153f0f69a4/surfe_lib/continuous_property.cpp#L475)
- [`vector_field.cpp` L105](https://github.com/MichaelHillier/surfe/blob/290dbe0ab344f4258a4935f05cad0f153f0f69a4/surfe_lib/vector_field.cpp#L105)

同样的 scalar/vector evaluation 循环也分别存在：

- [`single_surface.cpp` L577–L718](https://github.com/MichaelHillier/surfe/blob/290dbe0ab344f4258a4935f05cad0f153f0f69a4/surfe_lib/single_surface.cpp#L577-L718)
- [`lajaunie.cpp` L849–L970](https://github.com/MichaelHillier/surfe/blob/290dbe0ab344f4258a4935f05cad0f153f0f69a4/surfe_lib/lajaunie.cpp#L849-L970)
- [`stratigraphic_surfaces.cpp` L767–L892](https://github.com/MichaelHillier/surfe/blob/290dbe0ab344f4258a4935f05cad0f153f0f69a4/surfe_lib/stratigraphic_surfaces.cpp#L767-L892)
- [`continuous_property.cpp` L344–L447](https://github.com/MichaelHillier/surfe/blob/290dbe0ab344f4258a4935f05cad0f153f0f69a4/surfe_lib/continuous_property.cpp#L344-L447)
- [`vector_field.cpp` L163–L214](https://github.com/MichaelHillier/surfe/blob/290dbe0ab344f4258a4935f05cad0f153f0f69a4/surfe_lib/vector_field.cpp#L163-L214)

solve setup 也嵌在各 model method 中，例如 [`single_surface.cpp` L197–L272](https://github.com/MichaelHillier/surfe/blob/290dbe0ab344f4258a4935f05cad0f153f0f69a4/surfe_lib/single_surface.cpp#L197-L272)、[`lajaunie.cpp` L235–L272](https://github.com/MichaelHillier/surfe/blob/290dbe0ab344f4258a4935f05cad0f153f0f69a4/surfe_lib/lajaunie.cpp#L235-L272) 与 [`stratigraphic_surfaces.cpp` L648–L704](https://github.com/MichaelHillier/surfe/blob/290dbe0ab344f4258a4935f05cad0f153f0f69a4/surfe_lib/stratigraphic_surfaces.cpp#L648-L704)。重复不只在 kernel block，也贯穿 problem selection、solve 与 evaluation。

这不是 GeoRBF 的兼容目标。五种上层 workflow 可以 lowering 成不同的 functional/constraint collections，但只能共享一个：

- kernel jet；
- functional pairing；
- polynomial block；
- problem builder；
- solve dispatch；
- query evaluator；
- diagnostics。

### 7.2 modified kernel 与 restricted range

Surfe 在 restricted-range 分支关闭显式 polynomial term、打开 modified basis 并切换到 quadratic problem，[`lajaunie.cpp` L189–L208](https://github.com/MichaelHillier/surfe/blob/290dbe0ab344f4258a4935f05cad0f153f0f69a4/surfe_lib/lajaunie.cpp#L189-L208)。modified kernel：

- 从某个 horizon 中启发式选择四个点；
- 只支持一阶多项式；
- 只检测部分共面特例；
- 若找不到非共面点，直接给一个坐标加 `Epilson`，[`basis.cpp` L1572–L1823](https://github.com/MichaelHillier/surfe/blob/290dbe0ab344f4258a4935f05cad0f153f0f69a4/surfe_lib/basis.cpp#L1572-L1823)；
- 再手写 point/gradient/tangent 的 modified pairings，[`basis.cpp` L1963–L2833](https://github.com/MichaelHillier/surfe/blob/290dbe0ab344f4258a4935f05cad0f153f0f69a4/surfe_lib/basis.cpp#L1963-L2833)。

restricted-range 的“最小 seminorm 且满足有限不等式”思想有 Beatson 等人的变分依据。Beatson–Light–Billings 还说明：若 anchors 对 \(\Pi_{m-1}\) 真正 unisolvent，并以对应 Lagrange polynomials 给 CPD seminorm 增加 anchor-value norm，就能得到显式 SPD reproducing kernel。Surfe 的公式体现了这项思想；失败在于只支持四个一阶多项式 anchors、启发式从单个 horizon 取点，并在失败时改动用户坐标。GeoRBF 可保留“显式 gauge 后得到 SPD kernel”的数学选项，但必须拒绝其 anchor 与 epsilon 实现，并从 \(P\)、null space 与 CPD proof 推导凸问题。

### 7.3 可保留与必须拒绝

| Surfe 观察 | 处理 |
|---|---|
| 广义插值可组合值、gradient、tangent、inequality | **保留思想**，改为 first-class functional + relation |
| 同界面用 evaluation difference | **保留语义**，通用装配 |
| normal 为 gradient 信息 | **保留但拆清方向/模长/极性** |
| tangent 为 directional derivative | **保留精确线性语义** |
| restricted range 是 seminorm 下的有限约束优化 | **保留变分方向** |
| 全局 linear coordinate anisotropy | **保留显式 metric 形式** |
| 五套 matrix/evaluation/solver 流程 | **拒绝** |
| kernel pairwise virtual-method 笛卡尔积 | **拒绝** |
| 按 modelling mode 选择 polynomial/modified/QP | **拒绝** |
| 四点 modified kernel + epsilon 改坐标 | **拒绝** |
| normal component boxes 代表 angular uncertainty | **拒绝** |
| tangent angle 假定 gradient norm \(\approx2\) | **拒绝** |
| thickness 假定 gradient norm \(\approx1\) | **拒绝** |
| normal covariance eigenvalue floor 推断 metric | **不作默认；可作为未来 estimator 研究输入** |
| `partialPivLu` + finite check | **拒绝作为求解契约** |
| 自制 QP 的 bool 状态 | **拒绝作为优化契约** |
| 手工 GUI/VTK demo 当 `test` | **拒绝作为验证基线** |

### 7.4 测试事实

Surfe 的 CMake 只构造一个 `test` executable，没有注册 CTest 或断言式 suite，[`CMakeLists.txt` L82–L93](https://github.com/MichaelHillier/surfe/blob/290dbe0ab344f4258a4935f05cad0f153f0f69a4/CMakeLists.txt#L82-L93)。`test/main.cpp` 是带本机路径、VTK 输出与可视化的手工程序，[`test/main.cpp` L1–L95](https://github.com/MichaelHillier/surfe/blob/290dbe0ab344f4258a4935f05cad0f153f0f69a4/test/main.cpp#L1-L95)。

因此 Surfe 的数值输出只能作为迁移样例或差异调查输入，不能作为 GeoRBF 的正确性 oracle。固定源码里还存在 `SetTangentConstraints` 实际调用 `AddPlanarConstraintwNormal` 的明显 API 错配，[`surfe_api.cpp` L162–L180](https://github.com/MichaelHillier/surfe/blob/290dbe0ab344f4258a4935f05cad0f153f0f69a4/surfe_lib/surfe_api.cpp#L162-L180)，进一步说明历史 API 不应成为兼容目标。

## 8. GeoRBF v1 数学基线

以下条目可以直接进入 spec 的 normative requirements。

### 8.1 Functional 与 query

- MUST 支持有限个 value/first-derivative atoms 的线性组合。
- MUST 用同一个 bilinear pairing 装配 \(K_{ij}\)。
- MUST 用同一个 query functional 路径求 value、gradient 与 directional derivative。
- MUST 追踪 source ID、units、normalization 与 support points。
- MUST 在装配前拒绝 non-finite、零方向、非法 kernel 参数和不支持的 derivative order。
- SHOULD 允许内部 functional representation 扩展，而不承诺任意公开回调 ABI。

### 8.2 Kernel 与 polynomial

- MUST 在 registry 中携带精确 formula、dimension、SPD/带符号 CPD order、required polynomial、jet 与 parameter domain。
- MUST 从 kernel metadata 自动选择 required polynomial。
- MUST 对 \(P\) 做 rank/unisolvence 检查。
- MUST 显式处理 gauge；不得通过 coordinate jitter 或无报告 diagonal jitter 伪造唯一性。
- MUST 对 origin/cutoff 用解析值，并验证 Gram symmetry。
- MUST 拒绝缺少证明 metadata 的 kernel/functional 组合。

### 8.3 Constraint

- MUST 允许逐观测 hard/soft。
- MUST 支持 affine equality、lower、upper、interval。
- MUST 明确有限 support 上的 pointwise constraints 不等于整个连续域上的全局保证。
- MUST 支持带非负 slack 的显式 soft constraint。
- MUST 把同层位、层位差和层位顺序 lowering 成 evaluation differences。
- MUST 把完整 gradient observation 与仅 normal direction 分成不同语义。
- MUST 把精确 tangent lowering 成 directional derivative equality。
- MUST 不把 field separation 标为 physical thickness，除非存在 scale calibration。
- MUST 不把 component boxes 标为 rotation-invariant angular uncertainty。
- MAY 在 v1 支持已知 polarity normal SOC；若支持，必须处理零梯度退化。

### 8.4 Problem form 与 solver

- MUST 保留 symmetric structure。
- MUST 对 CPD convex problems 使用有证明的 null-space/projected formulation。
- MUST 以 problem capability 选择 linear/QP/SOCP backend，而非以模型名选择。
- MUST 返回结构化 termination、rank、residual、violation、scaling 与 source diagnostics。
- MUST 区分 invalid、rank deficient、infeasible、unbounded、iteration limit 和 numerical failure。
- SHOULD 首先实现小中规模稠密路径。
- SHOULD 让 matrix storage/operator 与 domain lowering 解耦，以便未来 sparse/matrix-free backend。
- MUST NOT 在本 spec 中绑定某个 Rust crate。

### 8.5 Anisotropy

- MUST 以对称正定 metric 或其可逆 factor \(A\) 表示全局 anisotropy。
- MUST 验证 eigenvalues、condition 和单位。
- MUST 对 gradient/Hessian pairing 使用物理坐标 chain rule。
- MUST NOT 把数据驱动 eigen-analysis 当默认真值。
- MUST NOT 在没有 nonstationary PD proof 时承诺局部 \(A(x)\) kernel。

## 9. 验证 oracle

### 9.1 Kernel jet

对每个准入核：

- 在随机非奇异点比较解析 gradient/mixed Hessian 与高精度 finite difference 或独立自动微分；
- 在 \(r\to0\) 验证解析极限，不通过把 \(r\) 改成 epsilon；
- 对紧支撑核验证 cutoff 左右的 value/derivative 连续阶；
- 验证 \(k(x,y)=k(y,x)\)、\(\nabla_xk=-\nabla_yk\)（平移不变时）及 mixed-Hessian transpose；
- 在 anisotropic transform 下比较 chain rule 与 transformed-coordinate oracle。

### 9.2 Functional algebra

- 随机生成 finite atoms，比较通用 contraction 与显式展开；
- difference–difference 必须等于四个 value pairing；
- tangent–normal 必须等于向量对 mixed Hessian 的 contraction；
- 输入重排只对矩阵做相同 permutation，不改变预测；
- 合并重复 atoms 前后结果一致；
- linearity：\((aL+bQ)f=aLf+bQf\)。

### 9.3 CPD/unisolvence

- 用 QR/SVD 构造 \(Z\approx\ker(P^\mathsf T)\)，验证
  \(Z^\mathsf TKZ\) 对准入 CPD 核为对称正定；
- manufactured polynomial \(p\in\Pi_{m-1}\) 必须被精确再现到容差；
- 删除识别某个 polynomial mode 的观测后，必须报告具体 rank loss；
- 仅 difference/derivative 的数据必须触发预期 constant gauge；
- 添加显式 gauge 后绝对 field 唯一，几何不变量保持不变。

### 9.4 Manufactured fields

至少使用：

- affine \(f(x)=a+b^\mathsf Tx\)：精确 value、gradient、tangent、difference；
- quadratic field：验证 value/gradient query 与 mixed pairings；
- 平行 horizons：验证 difference/order 与 affine coordinate transform；
- 已知 global anisotropy 的坐标变换对；
- 带已知噪声的 soft/hard 混合集。

验收不依赖 Surfe 数值。Surfe fixture 可以验证“我们理解了输入语义”，不能覆盖解析 oracle。

### 9.5 优化与诊断

对很小的问题构造可手算或可枚举的 oracle：

- hard equality + box 的 active set；
- soft equality 的 closed-form ridge 解；
- 一个 active、一个 inactive inequality；
- quadratic 与 \(L_1\) slack；
- normal SOC 的边界、内部、零梯度退化；
- 两条相互矛盾 hard constraints，必须报告 infeasible；
- 对约化 QP 检查 primal feasibility、dual feasibility、stationarity 与 complementarity；
- 对缩放前后问题验证物理单位中的结果一致。

### 9.6 数值与性能

- 扫描坐标尺度、shape 参数、相近 support points 与权重动态范围；
- 验证残差随 condition 恶化时诊断会升级，而非静默成功；
- benchmark 分开记录 lowering、dense assembly、factorization、solve 与 batch query；
- v1 性能门槛只针对小中规模，不能以没有 memory/operator seam 为代价；
- sparse 与 matrix-free 在 v1 只需接口/数据流可替换，不需要伪 benchmark。

## 10. 仍需由后续 spec 决定的问题

这些问题已经被缩小到可作产品选择的范围，不应在实现中隐式决定：

1. v1 最小 kernel 集是否只取 Gaussian、IMQ、cubic，还是同时纳入经验证的三维 Wendland C2。
2. kernel shape/length scale 的用户单位与自动缩放交互。
3. gauge 默认策略：要求用户提供、自动 anchor，还是 quotient-only 输出。
4. soft constraint 的损失族、权重/方差语义与是否支持相关 covariance。
5. normal 输入默认是完整 gradient、unit direction，还是两者显式不同类型。
6. normal polarity 的必填/未知行为。
7. normal angle SOCP 是否属于 v1；若属于，\(s_{\min}\) 如何规定。
8. tangent uncertainty 是否只提供 field-unit interval；角度语义是否推迟。
9. “thickness”在 v1 是否改名为 field separation；是否另立 scale-calibrated thickness 研究。
10. CPD QP 的 null-space/projected formulation 及 rank tolerance。
11. dense direct solver 所需的精确诊断和后端 capability interface。
12. global anisotropy 是只接收 \(A\)，还是接收 SPD \(M\) 并规范化 factor。
13. inequality/order 只承诺有限观测点约束，还是要求连续域上的保证；后者需要另一个证明与离散化/认证问题。

以下事项不再是开放方向：

- 不兼容 Surfe 历史 API；
- 不保留五套重复模型内核；
- 不准入没有 PD/CPD/jet 证据的核；
- 不以 \(\|\nabla f\|\approx1\) 或 \(\approx2\) 解释厚度/角度；
- 不用 epsilon 改输入几何来处理 unisolvence；
- 不把 `partialPivLu + isfinite` 或无完整状态的自制 QP 当 v1 稳定性标准；
- 不把 VTK/GUI 手工程序作为测试 oracle。

## 11. 可直接转为后续工作的边界

本研究支持把后续工作分成彼此可验收、不会重新制造五套内核的 slices：

1. functional algebra 与 canonical kernel jet spec；
2. kernel registry、PD/CPD metadata 与首批核证明；
3. polynomial/unisolvence/gauge preflight；
4. unified dense assembly 与 query evaluator；
5. hard/soft affine constraint lowering；
6. symmetric equality solver contract；
7. CPD null-space convex problem form；
8. QP 与可选 normal-angle SOCP capability；
9. geology semantics：horizon/order/field separation；
10. orientation semantics：gradient/normal/tangent/polarity；
11. global SPD anisotropy；
12. diagnostics、manufactured/property tests 与分阶段 benchmark；
13. sparse/matrix-free backend seam（仅接口，不承诺 v1 大规模实现）。

每个 slice 都应引用本文相应的“数学上必须”，并把第 10 节的产品选择写成显式 acceptance criteria。这样 tickets 不需要从 Surfe 类名出发，也不会为尚未消除的不确定性提前扩张范围。
