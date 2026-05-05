# JSPACE 后处理重构：数学流程与代码说明

这套补丁把 post-EM 后处理重写成一条更自洽、也更接近 HIGGS 的路径：

\[
\Gamma_{\text{biased}}
\;\xrightarrow{\text{debias}}\;
\Gamma_{\text{db}}
\;\xrightarrow{\text{Rayleigh search}}\;
M^\star
\;\xrightarrow{\text{optional rescue}}\;
M_{\text{final}}
\;\xrightarrow{\text{HIGGS-like refit}}\;
\Gamma_{\text{refit}}
\;\xrightarrow{\text{recolor}}\;
\Omega_{\text{final}}.
\]

---

## 1. Debias：只做 one-step correction，不再塞固定 Rayleigh

对每个频点 \(f\)：

\[
\Gamma_{\text{db},f}
=
2\Gamma_{\text{hat},f}
-
\Gamma_{\text{hat},f}\,\Sigma^e_f\,\Gamma_{\text{hat},f}.
\]

这里 \(\Sigma^e_f\) 是最终 EM 结束后的 whitened covariance。

为了给后面的 Rayleigh search 提供标准化尺度，同时返回方差 proxy：

\[
V_{f,ij}
=
|G_{ii}|\,|G_{jj}| + |G_{ij}|^2,
\]

其中 \(G\) 可以选 biased estimator（默认）或 debiased estimator。

> 关键改动：`module_debias.m` 不再做任何固定 `r=3` 的 Rayleigh 裁剪。

---

## 2. Search：用 debiased statistic 出 mask，但把 mask 施加到 ridge/Frobenius 模板上

### 2.1 Rayleigh mask

对任意候选 \(r\)，定义

\[
M_f(r)_{ij}
=
\mathbf 1\!\left(
|\Gamma_{\text{db},f,ij}|
\ge
\frac{r}{\sqrt{M_{\text{eff}}}}
\sqrt{\max(V_{f,ij},0)}
\right),
\]

并强制对角线保留。

这一步只决定 support，不直接把 \(\Gamma_{\text{db}}\) 当作最终 precision。

### 2.2 HIGGS-like ridge/Frobenius template

令

\[
\Sigma^e_f = U_f \operatorname{diag}(s_f) U_f^H.
\]

构造单频 ridge 模板：

\[
\Gamma_{\text{ridge},f}
=
U_f \operatorname{diag}\!\left(
\frac{\sqrt{s_f^2 + 4\rho}-s_f}{2\rho}
\right) U_f^H.
\]

这正是你给的 HIGGS 代码中的 spectral ridge/Frobenius 路径。

于是每个 \(r\) 的候选 precision 写成

\[
\Gamma_f(r)
=
\Pi_{H_{++}}
\!\left(
\mathcal H\big(M_f(r)\odot \Gamma_{\text{ridge},f}\big)
\right).
\]

也就是说：
- debiased matrix 负责“统计筛边”；
- ridge template 负责“给一个可 refit 的 SPD 基底”。

### 2.3 Score

补丁里提供三种 score：

1. `higgs`
\[
\sum_f \Big[\log\det \Gamma_f(r) - \operatorname{tr}(\Sigma^e_f\Gamma_f(r))\Big]
- \lambda_2 \sum_f \|W\odot \Gamma_f(r)\|_{1,\text{off}}.
\]

2. `jspace`
\[
\sum_f \Big[\log\det \Gamma_f(r) - \operatorname{tr}(\Sigma^e_f\Gamma_f(r))\Big]
- \lambda_1 R_{\text{freq}}
- \lambda_3 R_{\text{quad}}.
\]

3. `hybrid`（默认）
\[
\sum_f \Big[\log\det \Gamma_f(r) - \operatorname{tr}(\Sigma^e_f\Gamma_f(r))\Big]
- \lambda_1 R_{\text{freq}}
- \lambda_2 \|W\odot \Gamma\|_{1,\text{off}}
- \lambda_3 R_{\text{quad}}.
\]

再叠加与当前 JSPACE 一致的 density penalty。

---

## 3. Rescue：基于最终 whitened covariance 重新算阈值

Rescue 不再复用最开始的 `thr0`。  
补丁改为对最终的 `target_S = \Sigma^{e,\star}` 重新调用：

```matlab
[~, thr_end] = compute_scales_thresholds_(target_S, K_dev, W_dev, Nr, cfg);
```

得到：
- `thr_end.t_active`
- `thr_end.t_rescue`

默认 rescue mask：

\[
M^{\text{rescue}}_f
=
\mathbf 1\!\left(|\Sigma^{e,\star}_{f,ij}| > t^{\star}_{\text{rescue}}\right).
\]

### rescue 模式
- `post_rescue_mode = 'off'`
- `post_rescue_mode = 'auto'`
- `post_rescue_mode = 'always'`

在 `auto` 下，如果 Rayleigh search 得到的中位密度超出允许范围，就触发 rescue。

---

## 4. Refit：尽量靠近 HIGGS，同时保留跨频近似解析路径

### 4.1 单频 HIGGS-like refit（解析）
直接把最终 mask 施加到单频 ridge 模板：

\[
\Gamma^{\text{single}}_f
=
\Pi_{H_{++}}
\!\left(
\mathcal H\big(M^{\text{final}}_f \odot \Gamma_{\text{ridge},f}\big)
\right).
\]

这一步就是最接近 HIGGS 的“mask + eigendecomposition / SPD cleanup”。

### 4.2 跨频解析 surrogate refit（新增）
对每条边 \((i,j)\)，收集它在所有频点上的单频 ridge 初值：

\[
b_{ij} = [b_{ij}(1),\dots,b_{ij}(F)]^T.
\]

在 active 频点子集上，解

\[
\min_x
\frac12\|x-b_{ij}\|_2^2
+
\lambda_1 w_{ij} x^H L_K x
+
\lambda_3 w_{ij}\|x\|_2^2.
\]

闭式解为

\[
x_{ij}
=
\big(I + 2w_{ij}(\lambda_1 L_K + \lambda_3 I)\big)^{-1} b_{ij},
\]

如果不同频点 mask 不同，就只在 active 子集对应的主子矩阵上求解，再把 inactive 频点设成 0。

> 这一层是**解析的**，但它是针对 Frobenius / graph-Tikhonov surrogate 的解析解，  
> 不是完整 logdet + SPD + support-coupled JSPACE 原目标的一步闭式解。

---

## 5. 和现有文献的关系（只说最关键的）

- 单频 ridge 这条谱分解公式，和 HIGGS 代码完全同宗。
- 多图 / 多类的 ridge 型 precision estimation，Bilgrau et al. 2020 给出了显式单类 closed form，以及把 fused ridge 更新改写成非 fused 形式后再套 closed form 的递推方案。
- 但对于“带 logdet、带 support mask、带一般图拉普拉斯跨频耦合”的完整问题，我没有找到一个文献里已经给出的**一-shot 同时闭式解**。更常见的是：
  - joint graphical lasso：数值算法；
  - LASICH（Laplacian shrinkage）：ADMM；
  - targeted fused ridge：显式 special cases + 递推更新。

因此，这个补丁里给出的跨频解析 refit，最合适的数学定位是：

> **HIGGS 单频解析 ridge 路径 + JSPACE 跨频图拉普拉斯的 analytic surrogate smoothing。**

---

## 6. 新增 / 修改文件

- `module_debias.m`
- `module_higgs_ridge_template.m`
- `module_rayleigh_search.m`
- `module_higgs_refit.m`
- `solver_jspace_3d_opt_stoch.m`

---

## 7. 推荐默认配置

```matlab
cfg.post_run_search    = true;
cfg.post_rescue_mode   = 'auto';   % off | auto | always
cfg.post_rescue_union  = false;
cfg.post_refit_mode    = 'fused_ridge_surrogate';  % debias_only | single_ridge | fused_ridge_surrogate | stoch_mstep
cfg.post_search_score  = 'hybrid'; % higgs | jspace | hybrid
cfg.post_ridge_penalty = [];       % [] -> defaults to max(lambda2_final^2, 1e-3)
cfg.post_min_eig       = 1e-8;
cfg.post_store_aux     = false;
```

---

## 8. 你后面如果还想继续推，我建议优先比较三条后处理路径

1. `debias_only`
2. `single_ridge`
3. `fused_ridge_surrogate`

比较它们在：
- log-likelihood
- edge density
- SPD stability
- recolored precision 的 downstream 生物学一致性
- 频率连续性

上的差异，再决定是否真的需要保留 `stoch_mstep` refit。
