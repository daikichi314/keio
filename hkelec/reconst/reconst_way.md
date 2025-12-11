承知いたしました。これまでの対話内容、および文献（*NIM A 433, PTEP 2019* 等）に基づき、**APFit**と**fiTQun**、そして**あなたの手法**の違いを、数式を中心として体系的にまとめます。

このまとめは、そのまま指導教員への報告や、論文・修論の「既存手法との比較」の章の骨子として使えるレベルの完全性を目指します。

---

# スーパーカミオカンデ再構成アルゴリズムの完全な数理的比較

## 1. APFit (Atmospheric Pulse Fitter)
**特徴:** 「ロバストな頂点決定」と「パターン認識による方向決定」を段階的に行う手法。
**統計的立場:** M-推定（M-estimator）の一種。外れ値（散乱光）の影響を極小化するスコア関数（Goodness）を最大化する。

### 1-1. 頂点再構成 (Vertex Reconstruction)
APFitは、まず時間を手がかりに頂点 $\boldsymbol{x}$ を決定します。

**【最大化するスコア関数（Goodness）】**
$$G(\boldsymbol{x}, t_0) = \frac{1}{\sum_{i=1}^{N_{hit}} \frac{1}{\sigma_i^2}} \sum_{i=1}^{N_{hit}} \frac{1}{\sigma_i^2} \exp\left( -\frac{(t_{res, i}(\boldsymbol{x}, t_0))^2}{2(1.5 \cdot \sigma_i)^2} \right)$$

**各項の定義:**
* $N_{hit}$: ヒットしたPMTの総数（**Unhit PMTは含まない**）。
* $\sigma_i$: PMT $i$ の固有の時間分解能（観測電荷 $Q_i$ に依存する関数）。
* $1.5$: ガウス分布の裾野を少し広げ、多少の遅れを許容するための経験的係数。
* **$t_{res, i}$ (時間残差 Time Residual):**
    $$t_{res, i} = \underbrace{(t_{raw, i} - t_{cable, i} - t_{walk}(Q_i))}_{t_{obs, i} \text{: 補正済み観測時刻}} - \underbrace{\left( t_0 + \frac{|\boldsymbol{x}_{PMT, i} - \boldsymbol{x}|}{v_{group}} \right)}_{t_{exp, i} \text{: 予想時刻}}$$
    * $v_{group}$: 水中の光の群速度（約 $21.7$ cm/ns）。

**数学的意味:**
「残差がゼロに近い（＝ガウス分布のピークに乗っている）」PMTだけがスコアに寄与し、残差が大きい（＝散乱光やノイズ）PMTは $\exp(-\text{Large}) \approx 0$ となり無視されます。これにより**外れ値に対して極めて頑健**になります。

### 1-2. 方向再構成 (Direction Reconstruction)
頂点 $\boldsymbol{x}$ を固定した後、粒子の方向 $\boldsymbol{d}$ を決定します。

**【最大化する評価関数】**
$$Q_{dir}(\boldsymbol{d}) = \frac{\int_{0}^{\Theta_{max}} Q(\Theta) \cdot W(\Theta) d\Theta}{\int_{0}^{\Theta_{max}} Q(\Theta) d\Theta}$$
※ 実際は離散和で計算されます。

**各項の定義:**
* $\Theta$: 試行方向 $\boldsymbol{d}$ と、頂点から各PMTへのベクトルがなす角（Opening Angle）。
* $Q(\Theta)$: 角度 $\Theta$ の方向にあるPMTで観測された電荷の合計。
* $W(\Theta)$: **重み関数（チェレンコフリングの期待分布）**。
    * $$W(\Theta) \propto \exp\left( -\frac{(\Theta - \theta_C)^2}{2\sigma_\theta^2} \right)$$
    * $\theta_C$: チェレンコフ角（水中では約 $42^\circ$）。

**数学的意味:**
観測された電荷分布 $Q(\Theta)$ と、理想的なリング分布 $W(\Theta)$ の相関（マッチング度）を最大化しています。

---

## 2. fiTQun (Fast Inverse Technique ...)
**特徴:** 全てのパラメータ（頂点、方向、エネルギー、時間）を同時に決定する手法。
**統計的立場:** 最尤推定法（Maximum Likelihood Estimation）。物理現象の確率密度関数（PDF）を厳密に構築し、その積を最大化する。

### 2-1. 尤度関数 (Likelihood Function)
フィッティングパラメータ $\mathbf{x} = (x, y, z, t, \theta, \phi, p)$ に対する尤度 $L(\mathbf{x})$ は、全PMTの確率の積です。

$$L(\mathbf{x}) = \prod_{j \in Unhit} P_j(\text{unhit}|\mathbf{x}) \times \prod_{i \in Hit} P_i(hit|\mathbf{x})$$

解析では、この**負の対数尤度（Negative Log-Likelihood: NLL）**を最小化します。

$$-\ln L(\mathbf{x}) = \sum_{j \in Unhit} \underbrace{(-\ln P_j(unhit))}_{\text{Unhit項}} + \sum_{i \in Hit} \underbrace{(-\ln P_i(q_i, t_i))}_{\text{Hit項}}$$

### 2-2. Unhit項（光らなかった情報の活用）
PMTの応答がポアソン分布に従うと仮定します。
予想光量（期待値）を $\mu_j(\mathbf{x})$ とすると、観測値が0（Unhit）である確率は：
$$P(0|\mu_j) = e^{-\mu_j}$$
したがって、Unhit項はシンプルになります：
$$-\ln P_j(unhit) = -\ln(e^{-\mu_j}) = \mu_j(\mathbf{x})$$
**意味:** 「光るはず（$\mu$が大）なのに光らなかった」場合、その予想光量 $\mu$ がそのままペナルティとして加算されます。**壁際解析における最大の強みです。**

### 2-3. Hit項（電荷と時間の確率）
HitしたPMT $i$ の確率は、電荷と時間の同時確率です。
$$P_i(q_i, t_i|\mathbf{x}) \approx P_q(q_i | \mu_i) \times P_t(t_i | \mu_i, \mathbf{x})$$

**A. 電荷項 $P_q$ (Charge PDF)**
単純なポアソン分布ではなく、PMTの増幅揺らぎなども考慮した関数ですが、基本骨格はポアソンです。
$$- \ln P_q \approx \mu_i - q_i \ln \mu_i + \text{const.}$$
（これはあなたのコードに追加したBaker-Cousins $\chi^2$ の主要部と同じです）

**B. 時間項 $P_t$ (Time PDF)**
ガウス分布ではなく、**散乱光を含んだ混合分布**を使用します。
$$P_t(t) = (1 - w_{scat}) \cdot G(t - t_{exp}, \sigma) + w_{scat} \cdot S(t - t_{exp})$$
* $G(t, \sigma)$: 直接光成分（ガウス分布）。
* $S(t)$: 散乱・反射光成分（指数関数的減衰を持つ広い分布、またはテーブル参照）。
* $w_{scat}$: 散乱光が混入する確率（壁からの距離や水の透明度に依存）。

### 2-4. 予想光量モデル (Predicted Charge Model)
全ての基礎となる予想光量 $\mu_i(\mathbf{x})$ は、以下の物理モデルで計算されます。

$$\mu_i(\mathbf{x}) = \Phi(p) \times \int_{0}^{L(p)} \underbrace{g(\theta_{dir}(s))}_{\text{Cherenkov Profile}} \times \underbrace{\frac{e^{-d(s)/L_{att}}}{d(s)^2}}_{\text{Geometry \& Water}} \times \underbrace{\Omega(s)}_{\text{Solid Angle}} \times \underbrace{\epsilon(\eta(s))}_{\text{Angular Eff.}} ds$$

**各項の定義（完全版）:**
* **$\Phi(p)$ (Total Light Yield):** 運動量 $p$ の粒子が発生させる全光子数に比例する量。シミュレーション等で作成された**ルックアップテーブル**を参照します（関数として近似も可）。
* **$L(p)$ (Track Length):** 粒子の飛程。運動量 $p$ に依存するテーブル値。
* **$g(\theta_{dir})$:** チェレンコフ放射の角度分布。理想的には $\delta(\theta - 42^\circ)$ ですが、実際は散乱で広がったピークを持ちます。
* **$L_{att}$:** 水の光減衰長（Attenuation Length）。
* **$\epsilon(\eta)$:** PMTの入射角依存性（Angular Acceptance）。SKでは多項式 $\sum C_k \cos^k \eta$ で定義されます。

---

## 3. あなたの手法、APFit、fiTQunの比較表

| 項目 | あなたの現在の手法 (Reconstructor) | APFit (従来手法) | fiTQun (最新・壁際推奨) |
| :--- | :--- | :--- | :--- |
| **数学的原理** | **最小二乗法 ($\chi^2$)** | **適合度最大化 (Goodness)** | **最尤推定法 (Likelihood)** |
| **目的関数** | $\sum (t_{obs}-t_{exp})^2/\sigma^2 + \dots$ | $\sum \exp\left(-\frac{\Delta t^2}{2\sigma^2}\right)$ | $-\sum \ln P_{total}$ |
| **Unhit PMT** | **使用しない** (Hitのみ) | **使用しない** (Hitのみ) | **使用する** (Unhit項 $\sum \mu$ を加算) |
| **時間の扱い** | ガウス分布 ($\Delta t^2$) | ガウス分布 (ただし外れ値は無視) | **混合分布** (直接光 + 散乱光テール) |
| **電荷の扱い** | ガウス近似 (またはPoisson $\chi^2$) | 方向決定のみに使用 (Vertexには不使用) | ポアソン確率として頂点決定にもフル活用 |
| **光源モデル** | 点光源 (等方性または簡易指向性) | 点光源 (Vertex時) / リング (Dir時) | **有限長トラック + チェレンコフリング** |
| **外れ値耐性** | **弱い** (二乗で効くため引っ張られる) | **強い** (無視されるため頑健) | **強い** (散乱光PDFで確率的に処理) |
| **壁際での性能** | PMTが少ないと解が不定になりやすい | リングが欠けると精度・ロバスト性が崩壊 | Unhit情報が制約となり、精度を維持可能 |

---

## 4. 指導教員への報告・実装に向けたアドバイス

### あなたのコードで目指すべき修正（優先度順）

1.  **統計モデルの変更（ガウス $\to$ ポアソン）**
    * $\chi^2$ の項を Baker-Cousins 型 ($2(\mu - n + n \ln(n/\mu))$) に変更済みであればOKです。

2.  **Unhit情報の導入（必須）**
    * コード内のループを「HitしたPMT」だけでなく、「全PMT」に拡張してください。
    * 観測電荷 $n_i=0$ のPMTに対して、ペナルティ項 $+2\mu_i$ を $\chi^2$ に加算してください。これが壁際解析の決定打になります。

3.  **時間項のロバスト化**
    * 単純な $\Delta t^2$ の代わりに、**EMG (Exponentially Modified Gaussian)** 関数を用いた対数尤度項 $-\ln P_{EMG}(t)$ を導入することを強く推奨します。これにより、散乱光による「遅れ」がフィットを狂わせるのを防げます。

これらの情報を整理して伝えれば、物理的な理解も実装の方針も完璧であることを示せるはずです。頑張ってください！