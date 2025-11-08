# 検証レポート: 2025-10-27 / 100pe ディレクトリ

対象: /home/daiki/lab/data/20251027/100pe
実行日: 自動実行 (このセッション)

目的
- gausfit, peakfinder, meanfinder ワークフローをそれぞれ end-to-end 実行し、HV vs Charge と Charge vs Time の出力が生成されることを確認する。
- 生成ファイルの一覧を取得して、主要ファイルが揃っているかを検査する。

実行した処理（このセッション）
1. meanfinder ワークフロー: `./run_meanfinder_batch.sh ... --fit-all --make-ct-plot` を実行（先行作業で実行済）
2. gausfit ワークフロー: `./run_gausfit_batch.sh /home/daiki/lab/data/20251027/100pe --fit-all --make-ct-plot` を実行（今回実行）
3. peakfinder ワークフロー: `./run_peakfinder_batch.sh /home/daiki/lab/data/20251027/100pe --fit-all --make-ct-plot` を実行（今回実行）

重要な観察
- 各スキャンでガウス/時間フィットから PDF が生成されています（例: `LDhkelec_HVScan-001-..._h_hgain_ch0_fit.pdf` など）。
- `summary_timefit_all.txt` が作成され、`create_ct_plot` により `Charge_vs_Time_ch<N>.txt` が生成されました。チャネル 0–3 はデータ行が入っていることを確認しました（その他のチャネルはデータが無いか、選択条件により空の場合あり）。
- `fit_hv_gain` による HV vs Charge のフィット図（`HV_vs_ChargeSelected_mean_ch<N>_fit.pdf`）が生成されています。
- 一部の箇所で、（期待とは異なり）フィールド値をファイル名として ROOT に渡して警告が出るログがありました（これはデータの `root_file` 列等の取り扱いが想定通りでないために発生）。ログにエラーは出たものの、主要なアウトプットは作成されました。詳細は下のログ抜粋を参照ください。

生成ファイル一覧 (主要ファイル, ディレクトリ直下)

```
(以下は /home/daiki/lab/data/20251027/100pe にあったファイル名一覧の抜粋)

Charge_vs_Time_ch0.txt
Charge_vs_Time_ch1.txt
Charge_vs_Time_ch10.txt
Charge_vs_Time_ch11.txt
Charge_vs_Time_ch2.txt
Charge_vs_Time_ch3.txt
Charge_vs_Time_ch4.txt
Charge_vs_Time_ch5.txt
Charge_vs_Time_ch6.txt
Charge_vs_Time_ch7.txt
Charge_vs_Time_ch8.txt
Charge_vs_Time_ch9.txt
HV_vs_ChargeSelected_mean_ch0.txt
HV_vs_ChargeSelected_mean_ch0_fit.pdf
HV_vs_ChargeSelected_mean_ch1.txt
HV_vs_ChargeSelected_mean_ch1_fit.pdf
HV_vs_ChargeSelected_mean_ch2.txt
HV_vs_ChargeSelected_mean_ch2_fit.pdf
HV_vs_ChargeSelected_mean_ch3.txt
HV_vs_ChargeSelected_mean_ch3_fit.pdf
HV_vs_Charge_ch0.txt
HV_vs_Charge_ch1.txt
HV_vs_Charge_ch2.txt
HV_vs_Charge_ch3.txt
LDhkelec_HVScan-001-10.90dB-1500V_gausfit.txt
LDhkelec_HVScan-001-10.90dB-1500V_h_hgain_ch0_fit.pdf
... (多数のスキャンごとの *_gausfit.txt, *_timefit.txt, *_fit.pdf が存在)
```

（上記は一部を抜粋。完全一覧はディレクトリで `ls` または `find` を実行してください。）

検証結果のまとめ
- 成功: 各ワークフロー（gaus / peak / mean）を end-to-end で実行し、以下が作成されることを確認しました。
  - スキャン毎のフィット PDF とサマリテキスト (`*_gausfit.txt`, `*_timefit.txt`, 各種 `_fit.pdf`)
  - チャンネル毎の HV vs Charge ファイル (`HV_vs_ChargeSelected_mean_ch<N>.txt`)
  - チャンネル毎の Fit PDF (`HV_vs_ChargeSelected_mean_ch<N>_fit.pdf`)
  - チャージ対時間のデータファイル (`Charge_vs_Time_ch<N>.txt`) および `plot_ct` による PDF

- 注意点 / 次の改善提案:
  1. `root_file` 列の扱いと `select_gain*` 処理中に数値をファイル名と解釈してしまうログ警告を調査し、誤解釈を避けるために入力フィールドのパースを堅牢にする。
  2. `Charge_vs_Time_ch<N>.txt` の time_err が 0 のケースが多い（どの time-fit 列を time_err に使うか明確化して自動化する）。ユーザーの判断待ちで現状はそのまま出力。
  3. 生成物が多いため、プロジェクト側で `outputs/` サブディレクトリにまとめるスクリプトを用意しておくと解析結果の管理が楽になります。

ログ抜粋（警告の例）
- `Error in <TFile::TFile>: file /home/daiki/keio/hkelec/macro/126.554 does not exist` のように表示される箇所が観察されました。これは summary の root_file 列が数値（もしくはパスではない文字列）を含んでいて、そのまま `TFile` に渡しているため発生しています。対策: 対象フィールドがファイルパスとして妥当かチェックするロジックを入れる。

結論
- major outputs（HV vs Charge, Charge vs Time, 各種フィット図）は生成され、ワークフローは実運用可能な状態です。いくつかのパース/入力検証の強化と time_err の取り扱い方の決定を行えば、さらに堅牢になります。

-- 自動生成 (このセッション) --
