# hkelec/macro ワークフロー README

このドキュメントは、gausfit / peakfinder / meanfinder の各ワークフローを実行し、HV vs Charge および Charge vs Time のプロットを得る手順と出力フォーマットをまとめたものです。

対象ディレクトリ例
- データ/出力ディレクトリ: /home/daiki/lab/data/20251027/100pe

主要スクリプト
- `run_gausfit_batch.sh <target_dir> [--fit-all] [--make-ct-plot]`
  - ガウスフィット（既存の `gausfit` マクロ）をバッチ実行します。
  - `--fit-all` を付けると hgain/lgain/tot 等を全て処理します。
  - `--make-ct-plot` を付けると最終的に Charge vs Time 用ファイルを生成し、`plot_ct` を呼んで PDF を作ります。

- `run_peakfinder_batch.sh <target_dir> [--fit-all] [--make-ct-plot]`
  - ピーク検出（peakfinder）と時間フィット（TTS 系）をバッチで実行します。

- `run_meanfinder_batch.sh <target_dir> [--fit-all] [--make-ct-plot]`
  - meanfinder をバッチ実行します（平均値を用いる解析）。

- `run_hv_fitter.sh <target_dir> <gaus|mean|peak>`
  - 各チャンネルの HV vs Charge ファイルに対して `fit_hv_gain` を回し、フィット PDF を生成します。

- `run_ct_plotter.sh <target_dir> <mode>`
  - 既に生成された `Charge_vs_Time_ch<N>.txt` を `plot_ct` で PDF に変換します。

ビルド（必要なバイナリは `fit_results/` にあります）
- `cd fit_results && make create_ct_plot fit_hv_gain plot_ct`

出力ファイルとフォーマット
- Per-scan 出力（scan 単位）
  - `LDhkelec_HVScan-XXX-..._gausfit.txt` : ガウスフィットのサマリ
  - `LDhkelec_HVScan-XXX-..._timefit.txt` : 時間フィットの出力（各チャンネルごと）
  - `LDhkelec_HVScan-XXX-..._h_hgain_ch0_fit.pdf` 等: 各チャネルのフィット図

- HV vs Charge 系（チャンネル毎）
  - `HV_vs_Charge_ch<N>.txt` : 生の HV vs Charge（選択前）
  - `HV_vs_ChargeSelected_mean_ch<N>.txt` : 選択後（mean/peak により命名規則あり）。フォーマット:
    - 3 カラム: HV Charge Charge_err
    - 4 カラム (一部): HV Charge Charge_err extra
  - `HV_vs_ChargeSelected_mean_ch<N>_fit.pdf` : `fit_hv_gain` によるフィット図

- Charge vs Time（各チャンネル）
  - `Charge_vs_Time_ch<N>.txt` : 4 カラムで出力されます (ヘッダ付き)
    - charge (pC)
    - charge_err (pC) — 利用可能な場合
    - time_peak (ns)
    - time_err (ns) — 現状、time_err が 0 の場合があります（どの time-fit 列を使うか未決定）

注意・現在の挙動
- `create_ct_plot` は 2 通りの入力を受け取れるようにしました。
  - (A) mean-summary CSV（例: `summary_HV_vs_Charge_mean.txt`）
  - (B) per-channel `HV_vs_ChargeSelected_*_ch<N>.txt` ファイル（優先して使用）
- `plot_ct` はカンマ区切り/空白区切りの両方を受け取れるように柔軟化しています。
- `fit_hv_gain` は y エラーがあれば `TGraphErrors` を用いて加重フィットを行います。

サチュレーション / 選択ロジック
- 選択ロジックは既存の `select_gain` / `select_gain_mean` により行われます。要求に従い、飽和判定は「最右ビンとその左隣で比較し、右端 / 非ゼロ-2nd_right > 5 の場合は飽和」といった基準を採用しています（実装済み）。

トラブルシューティング
- Pedestal ROOT ファイルが見つからない場合、警告が出ますがテキストベースの代替ファイル (hkelec_pedestal_hithist_fits.txt) を作成して進める処理を組み込んでいます。
- `create_ct_plot` が空の `Charge_vs_Time_ch<N>.txt` を作る場合、まずディレクトリに `HV_vs_ChargeSelected_*_ch<N>.txt` が存在するか、`summary_timefit_all.txt` ができているか確認してください。

連絡先
- 自動化パッチとドキュメントは `hkelec/macro` 配下にあります。不明点があればこのリポジトリの編集履歴を確認してください。
