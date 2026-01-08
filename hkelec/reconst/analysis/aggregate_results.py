#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import glob
import re
import pandas as pd
import numpy as np

# ------------------------------------------------------------------
# 1. 使い方を表示する関数
# ------------------------------------------------------------------
def print_usage(prog_name):
    print("=" * 70)
    print(" [概要] ")
    print(" 解析結果のテキストファイル(*_fit_results.txt)を集計するスクリプト")
    print(" ファイル名から真の位置(x,y,z)を読み取り、解析結果と結合してCSVを出力します。")
    print(" ")
    print(" [使い方] ")
    print(f" {prog_name} <InputDirectory>")
    print(" ")
    print(" [引数] ")
    print(" <InputDirectory> : 解析結果txtファイルが入っているディレクトリパス")
    print("                    (例: ./data/20251120/reconst_images/)")
    print(" ")
    print(" [出力] ")
    print(" summary_results.csv : 集計結果CSV (入力ディレクトリに保存されます)")
    print("======================================================================")

# ------------------------------------------------------------------
# 2. ファイル名から真の位置を抽出する関数
# ------------------------------------------------------------------
def parse_filename_for_true_pos(filename):
    """
    ファイル名から xX_yY_zZ のパターンを探して数値を抽出する。
    マッチしない場合は (NaN, NaN, NaN) を返す。
    例: "run_x0_y-50_z100_fit_results.txt" -> (0.0, -50.0, 100.0)
    """
    # 正規表現: x(数字)_y(数字)_z(数字)
    # [+-]? : プラスマイナス符号があってもよい
    # \d+\.?\d* : 整数または小数
    pattern = r"x([+-]?\d+\.?\d*)_y([+-]?\d+\.?\d*)_z([+-]?\d+\.?\d*)"
    
    match = re.search(pattern, filename)
    if match:
        try:
            x_true = float(match.group(1))
            y_true = float(match.group(2))
            z_true = float(match.group(3))
            return x_true, y_true, z_true
        except ValueError:
            return np.nan, np.nan, np.nan
    else:
        return np.nan, np.nan, np.nan

# ------------------------------------------------------------------
# 3. テキストファイルの中身を解析する関数
# ------------------------------------------------------------------
def parse_result_file(filepath):
    """
    fit_results.txt の中身を読んで辞書形式で返す
    """
    data = {
        "x_res": np.nan, "x_sigma": np.nan,
        "y_res": np.nan, "y_sigma": np.nan,
        "z_res": np.nan, "z_sigma": np.nan,
        "t": np.nan,     "t_sigma": np.nan,
        "r_sigma": np.nan
    }
    
    current_section = None
    
    try:
        with open(filepath, "r") as f:
            for line in f:
                line = line.strip()
                if not line: continue
                
                # セクションの特定
                if "[Position X]" in line: current_section = "x"
                elif "[Position Y]" in line: current_section = "y"
                elif "[Position Z]" in line: current_section = "z"
                elif "[Time]" in line: current_section = "t"
                elif "[Radial Distance r" in line: current_section = "r"
                
                # 値の読み取り (行フォーマット: "Mean  : 1.23 +/- 0.1 cm" 等)
                if current_section in ["x", "y", "z", "t"]:
                    if line.startswith("Mean"):
                        # "Mean : " の後ろの数値を取得
                        parts = line.split(":")
                        if len(parts) > 1:
                            # "+/-" で分割して前の部分をとる
                            val_str = parts[1].split("+/-")[0].strip()
                            data[f"{current_section}_res" if current_section != "t" else "t"] = float(val_str)
                            
                    elif line.startswith("Sigma"):
                        parts = line.split(":")
                        if len(parts) > 1:
                            val_str = parts[1].split("+/-")[0].strip()
                            data[f"{current_section}_sigma"] = float(val_str)

                elif current_section == "r":
                    # rの場合は "Sigma (Resolution): ..." という行を探す
                    if line.startswith("Sigma"):
                        parts = line.split(":")
                        if len(parts) > 1:
                            val_str = parts[1].split("+/-")[0].strip()
                            data["r_sigma"] = float(val_str)
                            
    except Exception as e:
        print(f"Error reading {filepath}: {e}")
        
    return data

# ------------------------------------------------------------------
# メイン処理
# ------------------------------------------------------------------
if __name__ == "__main__":
    if len(sys.argv) < 2:
        print_usage(sys.argv[0])
        sys.exit(1)
        
    input_dir = sys.argv[1]
    
    if not os.path.isdir(input_dir):
        print(f"Error: Directory not found: {input_dir}")
        sys.exit(1)
        
    # 結果ファイルを検索
    # batch_analysis.py の出力は *_fit_results.txt
    search_pattern = os.path.join(input_dir, "*_fit_results.txt")
    files = glob.glob(search_pattern)
    
    if not files:
        print(f"No result files found in {input_dir}")
        sys.exit(0)
        
    print(f"Found {len(files)} files. Processing...")
    
    all_data = []
    
    for filepath in files:
        filename = os.path.basename(filepath)
        
        # 1. ファイル名から真値を抽出
        x_true, y_true, z_true = parse_filename_for_true_pos(filename)
        
        # 2. ファイルの中身から解析結果を抽出
        res_data = parse_result_file(filepath)
        
        # 3. データの結合
        row = {
            "x_true": x_true, "y_true": y_true, "z_true": z_true,
            **res_data # 辞書の展開
        }
        
        # 4. 計算列の追加 (d, diff)
        # 再構成位置が読み取れていない場合はNaNになるので計算もNaNになる
        row["x_diff"] = row["x_true"] - row["x_res"]
        row["y_diff"] = row["y_true"] - row["y_res"]
        row["z_diff"] = row["z_true"] - row["z_res"]
        
        # 距離 d = sqrt( (x_true - x_res)^2 + ... )
        row["d"] = np.sqrt(
            (row["x_diff"])**2 + 
            (row["y_diff"])**2 + 
            (row["z_diff"])**2
        )
        
        all_data.append(row)
        
    # DataFrame作成
    df = pd.DataFrame(all_data)
    
    # カラムの並び順を指定通りに
    cols_order = [
        "x_true", "y_true", "z_true", 
        "x_res", "x_sigma", "y_res", "y_sigma", "z_res", "z_sigma", 
        "t", "t_sigma", 
        "r_sigma", "d", 
        "x_diff", "y_diff", "z_diff"
    ]
    
    # 存在しないカラムがあれば除外して並べ替え
    final_cols = [c for c in cols_order if c in df.columns]
    df = df[final_cols]
    
    # CSV保存
    output_csv = os.path.join(input_dir, "summary_results.csv")
    df.to_csv(output_csv, index=False)
    
    print(f"Done! Summary saved to: {output_csv}")
    print(df.head())