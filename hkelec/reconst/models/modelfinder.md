# コード使い方ガイド（詳細解説）

## 1. 全体の流れ

このコードは **PMT半径 `r_pmt` の最適値を自動探索** するシステムです。以下の3段階で動作します：

````
段階1: データ準備
  ↓
段階2: 全r_pmt値を試して各モデルのフィット品質を計算
  ↓
段階3: 最適なr_pmtを決定し、その値で詳細解析を実行
````

---

## 2. 各関数の役割と使い方

### **関数① `analyze_pmt_center_height(r_pmt, df, show_plots=True)`**

**役割**: 固定されたr_pmt値で、f関数とg関数のフィットを実行

**入力パラメータ**:
- `r_pmt`: PMT半径 [cm]（例：20, 25, 30）
- `df`: 入力データフレーム
- `show_plots`: グラフを表示するか（True/False）

**返り値**:
- `results`: フィット結果のリスト（各PMTごと）
- `fit_functions`: フィット済み関数オブジェクト
- `chi2_totals`: χ²値（フィット品質）

**使用例**:
```python
# 例1: r_pmt = 20 cm でテスト実行（グラフ表示）
results, fit_funcs, chi2 = analyze_pmt_center_height(20, df, show_plots=True)

# 例2: r_pmt = 25 cm で実行（グラフなし）
results, fit_funcs, chi2 = analyze_pmt_center_height(25, df, show_plots=False)

# 結果の確認
print(results[0])  # PMT1の結果
# 出力例:
# {
#   'PMT_num': 1,
#   'n_data': 45,
#   'r_pmt': 25,
#   'z_pmt_center': 55.5,
#   'f_c0': 1234.56,
#   'f_c0_err': 12.34,
#   'f_chi2_ndf': 1.23,
#   'g_c0': 5678.90,
#   'g_c0_err': 56.78,
#   'g_chi2_ndf': 1.45
# }
```

**内部処理**:
1. PMT中心の z 座標を計算: `z_pmt_center = 80.5 - r_pmt`
2. 各PMTから光源までの距離を計算
3. f関数 (`y = c0 * (1 - sqrt(1 - (r_pmt/r)^2))`) でフィット
4. g関数 (`y = c0 / r^2`) でフィット
5. χ²/ndf（カイ二乗値）を計算

---

### **関数② `find_optimal_r_pmt(df, r_pmt_range=(10, 45), step=0.5)`**

**役割**: r_pmt を 10~45 cm で 0.5 cm 刻みに試し、最適値を探索

**入力パラメータ**:
- `df`: 入力データフレーム
- `r_pmt_range`: 探索範囲 (最小, 最大) [cm]
- `step`: 刻み幅 [cm]

**返り値**: 6個の値
```python
optimal_r_pmt_f,      # f関数での最適r_pmt
results_f,            # f関数での詳細結果
fit_functions_f,      # f関数での関数オブジェクト
optimal_r_pmt_g,      # g関数での最適r_pmt
results_g,            # g関数での詳細結果
fit_functions_g       # g関数での関数オブジェクト
```

**実行例**（コード内のセル#VSC-cd8d9a0c find_models-2.ipynb）:

```python
# 最適なr_pmtを探索（時間がかかります ≈ 5~10分）
optimal_r_pmt_f, results_f, fit_functions_f, \
optimal_r_pmt_g, results_g, fit_functions_g = find_optimal_r_pmt(df)

# 結果を表示
print(f"\n{'='*80}")
print(f"探索結果サマリー")
print(f"{'='*80}")
print(f"f関数での最適r_pmt: {optimal_r_pmt_f:.1f} cm")
print(f"g関数での最適r_pmt: {optimal_r_pmt_g:.1f} cm")
print(f"{'='*80}")

# 出力例:
# ================================================================================
# 探索結果サマリー
# ================================================================================
# f関数での最適r_pmt: 23.5 cm
# g関数での最適r_pmt: 24.0 cm
# ================================================================================
```

---

## 3. 具体的な使用シナリオ

### **シナリオA: クイックテスト（1つのr_pmt値を試す）**

```python
# Step 1: 単一のr_pmt値でテスト
print("=== r_pmt = 20 cm でのテスト ===")
results_test, fit_func_test, chi2_test = analyze_pmt_center_height(
    r_pmt=20, 
    df=df, 
    show_plots=True  # グラフを表示
)

# Step 2: 結果をPandasデータフレームで確認
import pandas as pd
results_df = pd.DataFrame(results_test)
print("\nフィット結果:")
print(results_df[['PMT_num', 'f_c0', 'f_chi2_ndf', 'g_c0', 'g_chi2_ndf']])

# 出力例:
#   PMT_num       f_c0  f_chi2_ndf       g_c0  g_chi2_ndf
# 0        1  1234.56        1.23  5678.90        1.45
# 1        2  1245.67        1.31  5689.01        1.52
# 2        3  1256.78        1.28  5700.12        1.48
# 3        4  1267.89        1.35  5711.23        1.55

# Step 3: χ²値で比較（f vs g）
for res in results_test:
    pmt = res['PMT_num']
    f_chi2 = res['f_chi2_ndf']
    g_chi2 = res['g_chi2_ndf']
    better = 'f' if f_chi2 < g_chi2 else 'g'
    print(f"PMT {pmt}: f({f_chi2:.3f}) vs g({g_chi2:.3f}) → {better} 関数が良い")
```

---

### **シナリオB: 複数のr_pmt値を比較する**

```python
# Step 1: 複数のr_pmt値を試す
r_pmt_candidates = [15, 20, 25, 30, 35]
comparison_results = {}

for r_pmt in r_pmt_candidates:
    print(f"\n探索中: r_pmt = {r_pmt} cm...")
    results, _, chi2_totals = analyze_pmt_center_height(
        r_pmt=r_pmt,
        df=df,
        show_plots=False  # グラフは表示しない
    )
    
    # chi2の合計を計算
    chi2_sum_f = sum([chi2_dict['f_chi2_ndf'] for chi2_dict in chi2_totals.values() 
                     if not np.isnan(chi2_dict['f_chi2_ndf'])])
    chi2_sum_g = sum([chi2_dict['g_chi2_ndf'] for chi2_dict in chi2_totals.values() 
                     if not np.isnan(chi2_dict['g_chi2_ndf'])])
    
    comparison_results[r_pmt] = {
        'chi2_f': chi2_sum_f,
        'chi2_g': chi2_sum_g,
        'results': results
    }

# Step 2: 結果を表にまとめる
print("\n" + "="*60)
print("r_pmt 比較表")
print("="*60)
print(f"{'r_pmt [cm]':>12} {'χ²(f)':>15} {'χ²(g)':>15}")
print("-"*60)
for r_pmt in r_pmt_candidates:
    f = comparison_results[r_pmt]['chi2_f']
    g = comparison_results[r_pmt]['chi2_g']
    print(f"{r_pmt:>12.1f} {f:>15.4f} {g:>15.4f}")
print("="*60)

# Step 3: 最適値を決定
best_r_pmt_f = min(comparison_results.items(), 
                   key=lambda x: x[1]['chi2_f'])[0]
best_r_pmt_g = min(comparison_results.items(), 
                   key=lambda x: x[1]['chi2_g'])[0]
print(f"\n最適なr_pmt:")
print(f"  f関数: {best_r_pmt_f:.1f} cm")
print(f"  g関数: {best_r_pmt_g:.1f} cm")
```

**出力例**:
```
============================================================
r_pmt 比較表
============================================================
r_pmt [cm]      χ²(f)      χ²(g)
------------------------------------------------------------
       15.0        5.234        6.123
       20.0        2.456        3.234
       25.0        1.234        2.123  ← f関数で最小
       30.0        1.789        1.234  ← g関数で最小
       35.0        3.456        2.345
============================================================

最適なr_pmt:
  f関数: 25.0 cm
  g関数: 30.0 cm
```

---

### **シナリオC: 全自動探索（推奨方法）**

```python
# Step 1: 最適なr_pmtを自動探索（5~10分かかります）
optimal_r_pmt_f, results_f, fit_functions_f, \
optimal_r_pmt_g, results_g, fit_functions_g = find_optimal_r_pmt(
    df=df,
    r_pmt_range=(10, 45),  # 10~45 cm を探索
    step=0.5               # 0.5 cm 刻み（計71点）
)

# Step 2: 結果の詳細確認
print("\n" + "="*80)
print("最適化結果")
print("="*80)

# f関数での最適値の詳細
print(f"\n【f関数での最適解】")
print(f"  最適r_pmt: {optimal_r_pmt_f:.1f} cm")
print(f"  Z_pmt_center: {80.5 - optimal_r_pmt_f:.1f} cm")
print(f"\n  PMTごとのフィット結果:")
for res in results_f:
    pmt = res['PMT_num']
    c0 = res['f_c0']
    c0_err = res['f_c0_err']
    chi2 = res['f_chi2_ndf']
    print(f"    PMT{pmt}: c0 = {c0:8.1f}±{c0_err:6.1f}  χ²/ndf = {chi2:.4f}")

# g関数での最適値の詳細
print(f"\n【g関数での最適解】")
print(f"  最適r_pmt: {optimal_r_pmt_g:.1f} cm")
print(f"  Z_pmt_center: {80.5 - optimal_r_pmt_g:.1f} cm")
print(f"\n  PMTごとのフィット結果:")
for res in results_g:
    pmt = res['PMT_num']
    c0 = res['g_c0']
    c0_err = res['g_c0_err']
    chi2 = res['g_chi2_ndf']
    print(f"    PMT{pmt}: c0 = {c0:8.1f}±{c0_err:6.1f}  χ²/ndf = {chi2:.4f}")

print("="*80)

# Step 3: 得られた関数を使って予測する
print("\n【予測例】")
test_r = 50  # r = 50 cm での予測
print(f"r = {test_r} cm での応答:")

# f関数での予測
if 1 in fit_functions_f:
    y_f = fit_functions_f[1]['f'](test_r)
    print(f"  f関数 (PMT1): {y_f:.2f}")

# g関数での予測
if 1 in fit_functions_g:
    y_g = fit_functions_g[1]['g'](test_r)
    print(f"  g関数 (PMT1): {y_g:.2f}")
```

**出力例**:
```
================================================================================
最適化結果
================================================================================

【f関数での最適解】
  最適r_pmt: 23.5 cm
  Z_pmt_center: 57.0 cm

  PMTごとのフィット結果:
    PMT1: c0 =  1234.5± 12.3  χ²/ndf = 1.2345
    PMT2: c0 =  1245.6± 13.4  χ²/ndf = 1.3456
    PMT3: c0 =  1256.7± 14.5  χ²/ndf = 1.2567
    PMT4: c0 =  1267.8± 15.6  χ²/ndf = 1.4678

【g関数での最適解】
  最適r_pmt: 24.0 cm
  Z_pmt_center: 56.5 cm

  PMTごとのフィット結果:
    PMT1: c0 =  5678.9± 56.7  χ²/ndf = 1.3456
    PMT2: c0 =  5689.0± 57.8  χ²/ndf = 1.4567
    PMT3: c0 =  5700.1± 58.9  χ²/ndf = 1.3678
    PMT4: c0 =  5711.2± 60.0  χ²/ndf = 1.5789

================================================================================

【予測例】
r = 50 cm での応答:
  f関数 (PMT1): 142.35
  g関数 (PMT1): 5.68
```

---

## 4. よくある質問

### **Q1: グラフが表示されない場合は？**

```python
# show_plots=True でもグラフが出ない場合
import matplotlib.pyplot as plt
plt.ion()  # インタラクティブモードをON

# その後で実行
analyze_pmt_center_height(25, df, show_plots=True)
```

### **Q2: χ²値が小さいほど良いモデル？**

```
はい。χ² が小さいほどデータにフィットしている。
一般に χ²/ndf ≈ 1 が良いとされている。
```

### **Q3: f関数とg関数どちらを使うべき？**

```python
# 結果を比較して判定
if optimal_r_pmt_f と optimal_r_pmt_g が近い:
    → どちらでもOK。物理的意味で選択
else:
    → χ²値が小さい方を選ぶ
```

---

## 5. コード実行の推奨順序

```python
# 1️⃣ 準備
df = pd.read_csv('~/lab/data/dataset/dataset_all.csv')

# 2️⃣ テスト（1分程度）
print("テスト実行...")
analyze_pmt_center_height(20, df, show_plots=True)

# 3️⃣ 本探索（5~10分）
print("最適値探索中...")
optimal_r_pmt_f, results_f, fit_functions_f, \
optimal_r_pmt_g, results_g, fit_functions_g = find_optimal_r_pmt(df)

# 4️⃣ 結果確認
print(f"f関数: {optimal_r_pmt_f:.1f} cm")
print(f"g関数: {optimal_r_pmt_g:.1f} cm")
```