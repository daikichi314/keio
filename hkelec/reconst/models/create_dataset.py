#!/home/daiki/keio/.venv/bin/python3
# -*- coding: utf-8 -*-

"""
ファイル名: create_dataset.py
概要:
    測定データ(*mean.txt)とPDモニタリングデータ(CSV)を統合し、
    電荷期待値モデルの構築・検証用データセット(CSV)を作成するスクリプト。

詳細:
    1. 指定されたディレクトリ内の電荷測定データ (*mean.txt) を読み込みます。
    2. ファイル名から光源位置(x, y, z)と減衰量(dB)情報を抽出します。
    3. PDモニタリングデータ (CSV) を読み込み、RunIDに基づいてリファレンス光量を計算します。
    4. 幾何学的パラメータ（距離 r, 入射角 cos(alpha)）を計算します。
       - PMT表面中心と、半球中心の2パターンの座標で計算を行います。
    5. 電荷データ (hgain または lgain) を抽出し、これらを統合して1つのCSVに出力します。

    ※ pc_by_h (High Gain) と pc_by_l (Low Gain) は同じチャンネルで重複しない前提で、
       存在する方のデータを採用します。
    ※ PMT番号は、データ上の ch(0~3) を 1~4 に変換して出力します。

使い方:
    python3 create_dataset.py [オプション]

オプション:
    -h, --help      このヘルプメッセージを表示して終了します。
    -i, --input_dir 測定データ(*mean.txt)が格納されているディレクトリ (必須)
    -p, --pd_csv    PDモニタリングデータのCSVファイルパス (必須)
    -o, --output    出力するCSVファイルパス。
                    省略時は入力ディレクトリ(-i)直下に 'dataset_for_modeling.csv' を作成します。

作成者: Gemini (Modified based on user request)
作成日: 2025-01-07
"""

import pandas as pd
import numpy as np
import glob
import os
import re
import sys
import argparse

# ==========================================
# 定数・設定定義
# ==========================================

# PMTの配置座標 (x, y) [cm]
# CH0~CH3 の位置定義
PMT_XY = {
    0: (-35.0,  35.0), # CH0 (PMT1)
    1: ( 35.0,  35.0), # CH1 (PMT2)
    2: (-35.0, -35.0), # CH2 (PMT3)
    3: ( 35.0, -35.0)  # CH3 (PMT4)
}

# Z座標の設定 [cm]
Z_SURFACE = 80.5    # PMT表面の中心
Z_HEMISPHERE = 48.0 # PMTを半球とみなした時の中心

# PMTの向き（法線ベクトル）
# ここでは全PMTがZ軸正方向を向いていると仮定 (0, 0, 1)
PMT_DIR = np.array([0.0, 0.0, 1.0])

# 光量計算のリファレンス設定
REF_ATTENUATION = 15.0 # dB (基準減衰量)
REF_PD_VOLT = 5.0      # V  (基準PD電圧)

# ==========================================
# 関数定義
# ==========================================

def calculate_light_power(pd_volt, attenuation_db):
    """
    PD電圧と減衰量から相対的な実際の光量(light_power)を計算する関数。
    
    計算式:
      Power = (V_meas / V_ref) * 10^((Att_ref - Att_meas) / 10)
    
    Args:
        pd_volt (float): 測定されたPD電圧 [V]
        attenuation_db (float): 設定された減衰量 [dB]
    
    Returns:
        float: 基準状態(15dB, 5V)を1とした時の相対光量
    """
    if pd_volt is None or np.isnan(pd_volt):
        return None
    
    # 電圧による比率 (電圧に比例と仮定)
    volt_ratio = pd_volt / REF_PD_VOLT
    
    # 減衰量による補正
    # 減衰量が小さいほど光量は大きくなるため、(Ref - Meas) / 10 乗となる
    attenuation_factor = 10 ** ((REF_ATTENUATION - attenuation_db) / 10.0)
    
    return volt_ratio * attenuation_factor

def calculate_geometry(source_pos, pmt_pos):
    """
    光源とPMTの位置関係から、距離(r)と入射角の余弦(cos alpha)を計算する関数。
    
    Args:
        source_pos (tuple): 光源の座標 (x, y, z)
        pmt_pos (tuple): PMTの座標 (x, y, z)
        
    Returns:
        tuple: (dist, cos_alpha)
            - dist: 距離 [cm]
            - cos_alpha: 光源方向ベクトルとPMT向きベクトルの内積 (正対時=1.0)
    """
    src = np.array(source_pos)
    pmt = np.array(pmt_pos)
    
    # PMTから光源へのベクトル (Incoming Light Directionの逆ベクトル)
    # これとPMTの向き(0,0,1)との内積をとる
    vec_p_to_s = src - pmt
    dist = np.linalg.norm(vec_p_to_s)
    
    if dist == 0:
        return 0.0, 0.0
    
    # 単位ベクトル化
    unit_vec = vec_p_to_s / dist
    
    # cos(alpha) = 単位ベクトル同士の内積
    cos_alpha = np.dot(unit_vec, PMT_DIR)
    
    return dist, cos_alpha

def parse_filename(filename):
    """
    ファイル名から実験パラメータを抽出する関数。
    
    期待する形式:
        BASENAME_xX_yY_zZ-RUNNUMBER-XXdB_mean.txt
        例: LDhkelec_x-35_y-35_z147-003-15.00dB_mean.txt
        
    Args:
        filename (str): ファイルパス
        
    Returns:
        tuple: (run_key, x, y, z, db)
            - run_key: PDデータと紐づけるためのキー（ファイル名ベース）
            - x, y, z: 光源座標
            - db: 減衰量
            解析失敗時は (None, None, None, None, None) を返す。
    """
    # 拡張子と末尾の識別子を除去してベース名を取得
    base_name_with_ext = os.path.basename(filename)
    base = base_name_with_ext.replace('_mean.txt', '')
    
    # 正規表現で数値を抽出 (負の数、小数に対応)
    # パターン: x(数値)_y(数値)_z(数値)-(ラン番号)-(数値)dB
    pattern = r"x([\d\.-]+)_y([\d\.-]+)_z([\d\.-]+)-(\d+)-([\d\.]+)dB"
    match = re.search(pattern, base)
    
    if match:
        try:
            x = float(match.group(1))
            y = float(match.group(2))
            z = float(match.group(3))
            # run = match.group(4) # 今回は使用しないが抽出可能
            db = float(match.group(5))
            return base, x, y, z, db
        except ValueError:
            return None, None, None, None, None
            
    return None, None, None, None, None

def process_data(mean_files_dir, pd_summary_csv, output_csv):
    """
    データの読み込み、統合、計算、出力を行うメイン処理関数。
    """
    print("========================================================")
    print(" データセット作成処理を開始します")
    print(f"  入力ディレクトリ: {mean_files_dir}")
    print(f"  PDデータ        : {pd_summary_csv}")
    print(f"  出力ファイル    : {output_csv}")
    print("========================================================")

    # 1. PD Summaryデータの読み込み
    print(f"[1/4] PD Summaryファイルを読み込んでいます...")
    try:
        df_pd = pd.read_csv(pd_summary_csv)
    except Exception as e:
        print(f"エラー: PDファイルの読み込みに失敗しました。\n詳細: {e}")
        return

    # PDデータの辞書化 (高速検索のため)
    # キー: RunIdのパス末尾のディレクトリ名, 値: PDvolt_mean
    # RunId例: /home/.../LDhkelec_x-0_y-35_z127-001-15.00dB
    pd_map = {}
    for _, row in df_pd.iterrows():
        try:
            run_id_path = row['RunId']
            # パスから末尾のディレクトリ名を抽出
            run_key = os.path.basename(run_id_path)
            pd_map[run_key] = row['PDvolt_mean']
        except KeyError:
            print("警告: PDファイルに 'RunId' または 'PDvolt_mean' 列がありません。スキップします。")
            continue

    # 2. mean.txtファイルの探索
    print(f"[2/4] 電荷測定データ(*mean.txt)を探索しています...")
    search_path = os.path.join(mean_files_dir, "*mean.txt")
    files = glob.glob(search_path)
    
    if not files:
        print(f"エラー: ディレクトリ '{mean_files_dir}' に '*mean.txt' ファイルが見つかりません。")
        return

    print(f" -> {len(files)} 個のファイルを検出しました。")

    output_rows = []
    
    # 3. 各ファイルの解析と統合
    print(f"[3/4] データを解析・統合しています...")
    
    for i, f in enumerate(files):
        # 進捗表示
        if (i+1) % 100 == 0:
            print(f"  処理中... {i+1}/{len(files)}")

        # ファイル名からパラメータ取得
        run_key, src_x, src_y, src_z, db = parse_filename(f)
        
        if run_key is None:
            # パースできないファイルはスキップ（警告はうるさいので出さないか、必要なら出す）
            # print(f"警告: ファイル名解析スキップ - {os.path.basename(f)}")
            continue
            
        # PDデータとのマッチング
        if run_key in pd_map:
            pd_volt = pd_map[run_key]
            light_power = calculate_light_power(pd_volt, db)
        else:
            # PDデータが見つからない場合はスキップ（光量が計算できないため）
            # print(f"警告: PDデータ不一致 - {run_key}")
            continue

        # ファイル内容読み込み
        try:
            # ファイル形式: # ch,type,mean,mean_err,rms,root_file
            df_meas = pd.read_csv(f, comment='#', names=['ch', 'type', 'mean', 'mean_err', 'rms', 'root_file'])
            
            # 文字列の前後の空白を除去
            df_meas['type'] = df_meas['type'].str.strip()
            
            # 必要なゲインタイプ (pc_by_h または pc_by_l) を抽出
            # ※同じchについて両方が現れることはない前提
            target_df = df_meas[df_meas['type'].isin(['pc_by_h', 'pc_by_l'])]
            
            if target_df.empty:
                continue

            for _, row in target_df.iterrows():
                ch = int(row['ch'])
                charge = row['mean']
                charge_err = row['mean_err']
                
                # PMT番号が定義済みのものかチェック
                if ch not in PMT_XY:
                    continue

                # PMT位置の取得
                pmt_x, pmt_y = PMT_XY[ch]
                
                # --- ジオメトリ計算 1: 表面中心 (Surface) ---
                pos_surf = (pmt_x, pmt_y, Z_SURFACE)
                r_surf, cos_surf = calculate_geometry((src_x, src_y, src_z), pos_surf)
                
                # --- ジオメトリ計算 2: 半球中心 (Hemisphere Center) ---
                pos_hemi = (pmt_x, pmt_y, Z_HEMISPHERE)
                r_hemi, cos_hemi = calculate_geometry((src_x, src_y, src_z), pos_hemi)

                # 結果リストに追加
                # ch (0~3) を PMT_num (1~4) に変換
                pmt_num = ch + 1

                output_rows.append({
                    '#PMT_num': pmt_num,
                    'Charge(pC)': charge,
                    'Charge_err(pC)': charge_err,
                    'light_power': light_power,
                    'x': src_x,
                    'y': src_y,
                    'z': src_z,
                    'r': r_surf,
                    'cos(alpha)': cos_surf,
                    'r_from_center': r_hemi,
                    'cos(alpha_from_center)': cos_hemi
                })

        except Exception as e:
            print(f"エラー: ファイル処理中に例外が発生しました {os.path.basename(f)}: {e}")
            continue

    # 4. CSV出力
    print(f"[4/4] CSVファイルに出力しています...")
    if output_rows:
        df_out = pd.DataFrame(output_rows)
        
        # カラム順序の指定
        cols = ['#PMT_num', 'Charge(pC)', 'Charge_err(pC)', 'light_power', 
                'x', 'y', 'z', 'r', 'cos(alpha)', 
                'r_from_center', 'cos(alpha_from_center)']
        
        # 念のためカラムの存在を確認して並べ替え
        df_out = df_out[cols]
        
        try:
            df_out.to_csv(output_csv, index=False)
            print(f"完了: {output_csv} を作成しました。")
            print(f"データ件数: {len(df_out)} 行")
        except Exception as e:
            print(f"エラー: CSVファイルの書き込みに失敗しました。: {e}")
    else:
        print("警告: 出力すべきデータが見つかりませんでした。入力ディレクトリや条件を確認してください。")

# ==========================================
# メインエントリポイント
# ==========================================
if __name__ == "__main__":
    # コマンドライン引数の設定
    parser = argparse.ArgumentParser(
        description='測定データ(*mean.txt)とPDデータ(CSV)を統合して解析用データセットを作成します。',
        epilog='使用例: ./create_dataset.py -i ./results -p ./data/PD_summary.csv',
        add_help=False
    )
    
    parser.add_argument('-h', '--help', action='store_true', help='ヘルプメッセージを表示します')
    parser.add_argument('-i', '--input_dir', type=str, help='測定データ(*mean.txt)があるディレクトリパス')
    parser.add_argument('-p', '--pd_csv', type=str, help='PDモニタリングデータのCSVファイルパス')
    parser.add_argument('-o', '--output', type=str, default=None, 
                        help='出力ファイル名 (省略時は入力ディレクトリ直下の dataset_for_modeling.csv)')

    args = parser.parse_args()

    # ヘルプ表示または引数不足時の処理
    if args.help or not args.input_dir or not args.pd_csv:
        print(__doc__) # モジュールのdocstringを表示
        if not args.help:
            print("\n[エラー] 必須引数 (-i, -p) が指定されていません。\n")
        sys.exit(0)

    # パスチェック
    if not os.path.exists(args.input_dir):
        print(f"エラー: 入力ディレクトリが見つかりません: {args.input_dir}")
        sys.exit(1)
    
    if not os.path.exists(args.pd_csv):
        print(f"エラー: PDファイルが見つかりません: {args.pd_csv}")
        sys.exit(1)

    # 出力パスの決定
    # -o が指定されていない場合は -i のディレクトリ直下に保存
    if args.output is None:
        output_path = os.path.join(args.input_dir, 'dataset_for_modeling.csv')
    else:
        output_path = args.output

    # 処理実行
    process_data(args.input_dir, args.pd_csv, output_path)