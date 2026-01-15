#!/home/daiki/keio/.venv/bin/python3
# -*- coding: utf-8 -*-

import os
import sys
import glob
import re
import pandas as pd
import numpy as np
import argparse

# ------------------------------------------------------------------
# 1. 使い方・仕様を表示する関数
# ------------------------------------------------------------------
def print_usage(prog_name):
    print("=" * 70)
    print(" [概要] ")
    print(" 解析結果テキスト(*_fit_results.txt)を集計し、外部PDデータと結合するスクリプト")
    print(" ファイル名から真の位置、再構成設定、減衰量を抽出し、")
    print(" 真の光量(PDデータより計算)と比較したサマリーCSVを作成します。")
    print(" ")
    print(" [使い方] ")
    print(f" {prog_name} <ResultDir> <PDSummaryCSV>")
    print(" ")
    print(" [引数] ")
    print(" <ResultDir>    : batch_analysis.py の出力先ディレクトリ (reconst_images)")
    print(" <PDSummaryCSV> : 真の光量計算に用いるPDデータのCSVファイルパス")
    print("                  (RunId, PDvolt_mean 等が含まれるファイル)")
    print(" ")
    print(" [出力] ")
    print(" summary_results.csv : ResultDir内に保存されます")
    print(" ")
    print(" [真の光量計算式] ")
    print(" A_true = (PDvolt / 5.0) * 10^((15.0 - Attenuation_dB) / 10.0)")
    print("======================================================================")

# ------------------------------------------------------------------
# 2. ファイル名解析関数
# ------------------------------------------------------------------
def parse_filename_info(filename):
    """
    ファイル名からパラメータを抽出する。
    
    期待形式:
      [BaseName]_reconst_[Hits]_[Qdist]_[Qmodel]_[Tdist]_fit_results.txt
      例: LDhkelec_x0_y0_z162-001-15.00dB_reconst_3hits_bc_func_f_gausT_fit_results.txt
    
    戻り値:
      (base_run_key, x_true, y_true, z_true, db, hits, qdist, q_mu_model, tdist)
    """
    # 末尾の _fit_results.txt を除去
    name_body = filename.replace("_fit_results.txt", "")
    
    # _reconst_ で分割
    if "_reconst_" not in name_body:
        return None
    
    parts = name_body.split("_reconst_")
    base_run_key = parts[0] # PDデータとの結合キー (例: LDhkelec_...-15.00dB)
    suffixes = parts[1]     # オプション部分 (例: 3hits_bc_func_f_gausT)
    
    # --- 真の位置とdBの抽出 ---
    # パターン: x(数値)_y(数値)_z(数値)...(数値)dB
    xyz_pat = r"x([+-]?\d+\.?\d*)_y([+-]?\d+\.?\d*)_z([+-]?\d+\.?\d*)"
    db_pat  = r"([\d\.]+)dB"
    
    x_true = y_true = z_true = np.nan
    db = np.nan
    
    m_xyz = re.search(xyz_pat, base_run_key)
    if m_xyz:
        x_true, y_true, z_true = map(float, m_xyz.groups())
        
    m_db = re.search(db_pat, base_run_key)
    if m_db:
        db = float(m_db.group(1))
        
    # --- 再構成設定の抽出 ---
    # suffixes を "_" で分割して解析
    opts = suffixes.split("_")
    
    hits = np.nan
    qdist = "unknown"
    q_mu_model = "unknown"
    tdist = "unknown"
    
    # Hits
    if "3hits" in opts: hits = 3
    elif "4hits" in opts: hits = 4
    
    # Q Dist
    if "bc" in opts: qdist = "bc"
    elif "gausQ" in opts: qdist = "gausQ"
    elif "noQ" in opts: qdist = "noQ"
    
    # Q Model (noQの場合はnone)
    if qdist == "noQ":
        q_mu_model = "none"
    else:
        if "func" in opts:
            idx = opts.index("func")
            if idx + 1 < len(opts):
                q_mu_model = "func_" + opts[idx+1] # func_f or func_g
    
    # T Dist
    if "gausT" in opts: tdist = "gausT"
    elif "goodness" in opts: tdist = "goodness"
    elif "emg" in opts: tdist = "emg"
    elif "noT" in opts: tdist = "noT"
    
    return base_run_key, x_true, y_true, z_true, db, hits, qdist, q_mu_model, tdist

# ------------------------------------------------------------------
# 3. テキスト結果読み込み関数
# ------------------------------------------------------------------
def parse_result_txt(filepath):
    """
    fit_results.txt の中身を読み取り、辞書を返す。
    値が見つからない場合は -9999 (NaN相当) を入れる。
    """
    data = {
        "x_res": -9999, "x_sigma": -9999,
        "y_res": -9999, "y_sigma": -9999,
        "z_res": -9999, "z_sigma": -9999,
        "t": -9999,     "t_sigma": -9999,
        "A_res": -9999, "A_sigma": -9999,
        "r_sigma": -9999
    }
    
    current_section = None
    
    try:
        with open(filepath, "r") as f:
            for line in f:
                line = line.strip()
                if not line: continue
                
                # セクション判定
                if "[Position X]" in line: current_section = "x"
                elif "[Position Y]" in line: current_section = "y"
                elif "[Position Z]" in line: current_section = "z"
                elif "[Time]" in line: current_section = "t"
                elif "[Amplitude A]" in line: current_section = "A"
                elif "[Radial Distance r" in line: current_section = "r"
                
                # 値の抽出
                if current_section in ["x", "y", "z", "t", "A"]:
                    if line.startswith("Mean"):
                        parts = line.split(":")
                        if len(parts) > 1:
                            val_str = parts[1].split("+/-")[0].strip()
                            key = "t" if current_section == "t" else f"{current_section}_res"
                            try: data[key] = float(val_str)
                            except: pass
                            
                    elif line.startswith("Sigma"):
                        parts = line.split(":")
                        if len(parts) > 1:
                            val_str = parts[1].split("+/-")[0].strip()
                            key = f"{current_section}_sigma"
                            try: data[key] = float(val_str)
                            except: pass

                elif current_section == "r":
                    if line.startswith("Sigma"):
                        parts = line.split(":")
                        if len(parts) > 1:
                            val_str = parts[1].split("+/-")[0].strip()
                            try: data["r_sigma"] = float(val_str)
                            except: pass
                            
    except Exception as e:
        print(f"Error reading {filepath}: {e}")
        
    return data

# ------------------------------------------------------------------
# 4. 真の光量計算
# ------------------------------------------------------------------
def calculate_true_A(pd_volt, db):
    # A_true = (PDvolt / 5.0) * 10^((15.0 - dB) / 10.0)
    if pd_volt is None or np.isnan(pd_volt) or np.isnan(db):
        return -9999
    return (pd_volt / 5.0) * (10 ** ((15.0 - db) / 10.0))

# ------------------------------------------------------------------
# メイン処理
# ------------------------------------------------------------------
if __name__ == "__main__":
    if len(sys.argv) < 3:
        print_usage(sys.argv[0])
        sys.exit(1)
        
    result_dir = sys.argv[1]
    pd_csv_path = sys.argv[2]
    
    if not os.path.exists(result_dir):
        print(f"Error: Result directory not found: {result_dir}")
        sys.exit(1)
    if not os.path.exists(pd_csv_path):
        print(f"Error: PD CSV file not found: {pd_csv_path}")
        sys.exit(1)
        
    # PDデータの読み込みとマップ化
    print("Loading PD Summary Data...")
    try:
        df_pd = pd.read_csv(pd_csv_path)
        # RunIdのパスからファイル名(キー)を抽出
        # キー: "LDhkelec_x...-15.00dB", 値: PDvolt_mean
        pd_map = {}
        for _, row in df_pd.iterrows():
            run_path = row['RunId']
            key = os.path.basename(run_path)
            pd_map[key] = row['PDvolt_mean']
    except Exception as e:
        print(f"Error loading PD CSV: {e}")
        sys.exit(1)
        
    # 結果ファイルの処理
    files = glob.glob(os.path.join(result_dir, "*_fit_results.txt"))
    if not files:
        print("No result files found.")
        sys.exit(0)
        
    print(f"Processing {len(files)} result files...")
    
    all_data = []
    
    for f in files:
        fname = os.path.basename(f)
        
        # 1. ファイル名解析
        f_info = parse_filename_info(fname)
        if not f_info:
            print(f"Skipping unknown format: {fname}")
            continue
            
        base_key, x_t, y_t, z_t, db, hits, qdist, q_mu, tdist = f_info
        
        # 2. テキスト解析
        res_vals = parse_result_txt(f)
        
        # 3. 真の光量計算
        pd_val = pd_map.get(base_key, np.nan)
        A_true = calculate_true_A(pd_val, db)
        
        # 4. 行データの作成
        row = {
            # True Pos
            "x_true": x_true, "y_true": y_true, "z_true": z_t,
            # Reconst Settings
            "hits": hits, "qdist": qdist, "q_mu_model": q_mu, "tdist": tdist,
            "attenuation_dB": db,
            # Parsed Results
            **res_vals,
            # True A
            "A_true": A_true
        }
        
        # 5. 計算列 (Diff = Res - True)
        # 値が-9999(欠損)の場合は計算結果もNaNまたは無効値にする処理が必要だが、
        # ここでは単純計算し、後でフィルタリングする方針とする
        
        if row["x_res"] != -9999:
            row["d"] = np.sqrt((row["x_res"] - x_t)**2 + (row["y_res"] - y_t)**2 + (row["z_res"] - z_t)**2)
            row["x_diff"] = row["x_res"] - x_t
            row["y_diff"] = row["y_res"] - y_t
            row["z_diff"] = row["z_res"] - z_t
        else:
            row["d"] = -9999
            row["x_diff"] = -9999
            row["y_diff"] = -9999
            row["z_diff"] = -9999
            
        if row["A_res"] != -9999 and row["A_true"] != -9999:
            row["A_diff"] = row["A_res"] - row["A_true"]
        else:
            row["A_diff"] = -9999
            
        all_data.append(row)
        
    # DataFrame化と列の並び替え
    df_out = pd.DataFrame(all_data)
    
    # 指定された順序
    cols_order = [
        "x_res", "y_res", "z_res",          # 再構成位置
        "x_sigma", "y_sigma", "z_sigma",    # 位置不確かさ
        "r_sigma",                          # 位置分解能
        "x_true", "y_true", "z_true",       # 真の位置
        "d",                                # 距離
        "x_diff", "y_diff", "z_diff",       # 差分 (Res - True)
        "t", "t_sigma",                     # 再構成時間
        "A_res", "A_sigma",                 # 光量パラメータ
        "A_true",                           # 真の光量
        "A_diff",                           # A差分
        "attenuation_dB",                   # dB
        "hits", "qdist", "q_mu_model", "tdist" # 再構成設定
    ]
    
    # 存在しない列を除外して並べ替え
    final_cols = [c for c in cols_order if c in df_out.columns]
    df_out = df_out[final_cols]
    
    # 保存
    out_path = os.path.join(result_dir, "summary_results.csv")
    df_out.to_csv(out_path, index=False)
    
    print(f"Done. Summary saved to: {out_path}")
    print(df_out.head())