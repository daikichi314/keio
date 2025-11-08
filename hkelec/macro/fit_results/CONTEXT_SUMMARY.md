# 解析ワークフロー — 状況サマリ

作成日: 2025-11-04

このファイルは、これまでの対話で行った変更・実装・残作業をまとめたものです。将来の作業の参照としてここに保存します。

## 目的
複数手法（gausfit / peakfinder / meanfinder）で得られる電荷情報（pC）と時間情報（time_diff, ns）を用いて
- HV (V) vs Charge (pC)
- Charge (pC) vs Time (ns)
のグラフ（エラーバー含む）を作成するための一連の解析パイプラインを整備する。

## 解析パターン（再掲）
1. gausfit: 電荷ヒストグラムをガウスでフィット -> peakを取得。time_diffをEMGでフィット。
2. peakfinder: 電荷ヒストグラムの最大ビン中心を電荷値とする（フィッティングなし）。time_diffはEMGでフィット。
3. meanfinder: 電荷ヒストグラムの平均値（GetMean）を電荷値とする。time_diffはEMGでフィット。

## これまでに行った編集・追加（概要）
- 新規追加
  - `plot_ct.C` — Charge vs Time のエラーバー付きプロットを作成する小さなROOT実行ファイルを追加。
  - `run_ct_plotter.sh` — `plot_ct` をディレクトリ内のファイルに対してまとめて実行するスクリプトを追加。

- 変更（編集）
  - `Makefile` — `plot_ct` ターゲットを追加し `all` に含めるよう修正。
  - `select_gain.C` — 以下の改良を実施:
    - gausfit系とmeanfinder系の両方のサマリーファイル形式を受け取れるよう、行を柔軟に解析するように変更。
    - サチュレーション判定を強化（ROOTヒストグラム名パターンを複数試行するようにし、最後のビンと右方向に2番目（空ビンをスキップしたもの）を比較して5倍ルールで判定）。
    - 出力に電荷誤差を追加（pedestal誤差とpeak誤差の伝播）。
  - `meanfinder.C` — サチュレーション判定ロジックを内部で行わないように変更（各ヒストの mean, mean_err, rms, root_file を出力する形式に変更）。
  - `create_ct_plot.C` — （注: その後ユーザがこのファイルの編集を元に戻しました。最新のワークスペースでは編集が取り消されている可能性があります。必ず実際のファイル内容を確認してください。）

## ユーザによって元に戻されたファイル
- `/home/daiki/keio/hkelec/macro/fit_results/create_ct_plot.C` — このファイルへの編集はユーザによって取り消されました。現在のワークツリーで編集が反映されているか必ず確認してください。

## 追加 / 変更したファイル一覧（作業履歴）
- 新規: `fit_results/plot_ct.C`
- 新規: `run_ct_plotter.sh`
- 変更: `fit_results/select_gain.C`
- 変更: `fit_results/meanfinder.C`
- 変更: `fit_results/Makefile`
- 変更（編集後に元に戻された可能性あり）: `fit_results/create_ct_plot.C`

## 出力ファイル・命名規則（現状の方針）
- gausfit: `_gausfit.txt`、HV vs Charge: `HV_vs_Charge_ch<ch>.txt`（gaus系用）
- peakfinder: `_peak.txt`、HV vs Charge: `HV_vs_Charge_peak_hgain_ch<ch>.txt` （peak専用命名）
- meanfinder: `_mean.txt`（各ヒストの mean を出力）、HV vs Charge: `HV_vs_Charge_mean_ch<ch>.txt`
- timefit（共通）: `_timefit.txt`（全手法同一フォーマット）
- Charge vs Time（最終的に作るファイル）: `Charge_vs_Time_ch<ch>.txt`（中身: Charge, Charge_err, Time_peak, Time_peak_err）

## 動作確認・実行方法（現時点）
1. ビルド

```bash
cd ~/hkelec/macro/fit_results
make
```

2. Charge vs Time をプロット

```bash
cd ~/hkelec/macro
./run_ct_plotter.sh <対象ディレクトリ> <gaus|peak|mean>
```

3. meanfinder バッチ（注: まだ更新が必要）
- `run_meanfinder_batch.sh` は meanfinder を起動して `_mean.txt` を集約しますが、現状 `select_gain_mean` 呼び出しの確認・修正が残っています。

## 未完了（残作業）
- run_meanfinder_batch.sh の更新: meanfinder の出力を取りまとめ、`select_gain_mean` を呼んで HV vs Charge(mean) ファイルを作成するように修正（未完）。
- `select_gain_mean.C` の最終チェック（既存 `select_gain_mean.C` がある場合は内容と整合させる必要あり）。
- 全体のビルド & ランによる統合テスト（`make` → バッチ実行 → 出力ファイル確認）。
- `create_ct_plot.C` のリストア or 再適用（ユーザの意向に沿って、必要なら再実装）。

## 現在の todo（短期）
1. `run_meanfinder_batch.sh` を `select_gain_mean` 呼び出しに対応させる（優先）。
2. `select_gain_mean.C` の存在確認と調整（gaus/mean/peak 共通で使えるか確認）。
3. `create_ct_plot.C` の現在版を確認し、必要なら再編集。
4. `make` → 簡易テスト（sampleファイルがあれば1チャネルで動かす）。

## メモ・注意点
- ROOT のインクルードパスは開発環境に依存します（`root-config --cflags` / `--glibs` を Makefile で使用）。
- サチュレーション判定は「最後のビンのエントリ数が、右から非空のビンのエントリ数の5倍以上」というルールに統一しました。
- gaus / mean の summary フォーマットは微妙に異なるため、`select_gain.C` を両方に対応するように拡張しています。出力の列順・区切りの厳密なフォーマットが異なるとパースに失敗する可能性があるため、実際のsummaryファイル（`summary_mean_all.txt` / `summary_peak_all.txt` / `summary_*`）を一度確認してください。

---

必要ならこの要約を別名で保存（例: `WORKLOG.md` や issue テンプレート）します。次に何を保存/編集しますか？

