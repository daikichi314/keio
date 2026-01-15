#!/home/daiki/keio/.venv/bin/python3
# -*- coding: utf-8 -*-

import sys
import os
import glob
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from mpl_toolkits.mplot3d import Axes3D

# プロットのフォントサイズ設定
plt.rcParams['font.size'] = 12

# ------------------------------------------------------------------
# 1. 使い方・仕様を表示する関数
# ------------------------------------------------------------------
def print_usage(prog_name):
    print("=" * 70)
    print(" [概要] ")
    print(" 再構成結果CSVファイルから詳細な解析画像とフィット結果を出力するスクリプト")
    print(" ")
    print(" [処理内容] ")
    print(" 1. 1次元分布 (X, Y, Z, t, A) のヒストグラムとガウスフィット")
    print("    ※ Aが計算されていない(-9999)場合はスキップします")
    print(" 2. 3次元距離 r のヒストグラム (X,Y,Zのフィット中心からの距離)")
    print(" 3. カイ二乗 (chi2) 分布のヒストグラム")
    print(" 4. 2次元分布 (XY, YZ, XZ) のヒストグラム (広域・拡大)")
    print(" 5. 3次元散布図 (カラー: 時間)")
    print(" ")
    print(" [使い方] ")
    print(f" {prog_name} <InputPath>")
    print(" ")
    print(" [引数] ")
    print(" <InputPath> : CSVファイルパス または ディレクトリパス")
    print(" ")
    print(" [出力] ")
    print(" 入力ディレクトリ直下の 'reconst_images/' に画像(.pdf)と結果(.txt)を保存")
    print(" 具体的には、各CSVファイル 'basename.csv' に対し、'basename_hist_X.pdf' などの画像と 'basename_results.txt' などの結果ファイルを作成します。")
    print("======================================================================")

# ------------------------------------------------------------------
# ガウス関数定義
# ------------------------------------------------------------------
def gaussian(x, a, mu, sigma):
    return a * np.exp(-(x - mu)**2 / (2 * sigma**2))

# ------------------------------------------------------------------
# 1次元ヒストグラム + ガウスフィットを行う関数
# ------------------------------------------------------------------
def process_1d_fit(df, col, label, unit, output_dir, base_name, text_lines):
    """
    戻り値: (center_val, fit_success_flag)
    ※ center_val はフィット成功時はフィットMean、失敗時は算術平均を返す
    """
    # データが存在しない、または全て-9999(計算除外)の場合はスキップ
    if col not in df.columns:
        text_lines.append(f"[{label}]")
        text_lines.append("  Not Found in CSV\n")
        return 0, False

    data = df[col]
    # -9990以下は計算されていない値とみなす
    valid_data = data[data > -9990]

    if len(valid_data) == 0:
        text_lines.append(f"[{label}]")
        text_lines.append("  Not Calculated (-9999)\n")
        return 0, False

    plt.figure(figsize=(8, 6))
    
    counts, bin_edges = np.histogram(valid_data, bins=50)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    
    # 初期パラメータ
    arithmetic_mean = np.mean(valid_data)
    p0 = [max(counts), arithmetic_mean, np.std(valid_data)]
    
    fit_success = False
    popt = [0, 0, 0]
    perr = [0, 0, 0]

    try:
        popt, pcov = curve_fit(gaussian, bin_centers, counts, p0=p0, maxfev=5000)
        perr = np.sqrt(np.diag(pcov))
        fit_success = True
    except:
        pass 

    plt.hist(valid_data, bins=50, color='skyblue', edgecolor='white', label='Data')
    
    if fit_success:
        x_smooth = np.linspace(min(valid_data), max(valid_data), 200)
        plt.plot(x_smooth, gaussian(x_smooth, *popt), 'r-', lw=2, label='Fit')
        
        info = (f"Mean: {popt[1]:.2f} {unit}\n"
                f"Sigma: {abs(popt[2]):.2f} {unit}")
        plt.gca().text(0.05, 0.95, info, transform=plt.gca().transAxes,
                       verticalalignment='top', 
                       bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        
        text_lines.append(f"[{label}]")
        text_lines.append(f"  Mean  : {popt[1]:.4f} +/- {perr[1]:.4f} {unit}")
        text_lines.append(f"  Sigma : {abs(popt[2]):.4f} +/- {perr[2]:.4f} {unit}\n")
    else:
        text_lines.append(f"[{label}]")
        text_lines.append("  Fit Failed\n")

    plt.title(f'{label} Distribution')
    plt.xlabel(f'{label} [{unit}]')
    plt.ylabel('Events')
    plt.legend()
    plt.tight_layout()
    
    save_path = os.path.join(output_dir, f"{base_name}_hist_{col}.pdf")
    plt.savefig(save_path)
    plt.close()
    
    # 次の処理のために中心値を返す
    center_val = popt[1] if fit_success else arithmetic_mean
    return center_val, fit_success

# ------------------------------------------------------------------
# 3次元距離 r の分布 (ミラーリングフィット)
# ------------------------------------------------------------------
def process_radial_fit(df, center_pos, output_dir, base_name, text_lines):
    """
    center_pos: (mu_x, mu_y, mu_z) のタプル。ガウスフィットの中心を使う。
    """
    label = "Radial Distance r"
    unit = "cm"
    
    mu_x, mu_y, mu_z = center_pos
    
    # -9999を含まないデータのみ使用（念のためフィルタ）
    df_valid = df[(df['fit_x'] > -9990) & (df['fit_y'] > -9990) & (df['fit_z'] > -9990)]
    
    if len(df_valid) == 0:
        text_lines.append(f"[{label} (Mirrored)]")
        text_lines.append("  Not Calculated (No valid XYZ)\n")
        return

    # 1. 指定された中心(ガウスフィット結果)からの距離 r を計算
    r = np.sqrt( (df_valid['fit_x'] - mu_x)**2 + 
                 (df_valid['fit_y'] - mu_y)**2 + 
                 (df_valid['fit_z'] - mu_z)**2 )
    
    # 2. データをミラーリング (r と -r を結合)
    data_mirror = np.concatenate([-r, r])
    
    plt.figure(figsize=(8, 6))
    
    counts, bin_edges = np.histogram(data_mirror, bins=100)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    
    # 初期パラメータ
    p0 = [max(counts), 0, np.std(r)]
    
    fit_success = False
    popt = [0, 0, 0]
    perr = [0, 0, 0]

    try:
        popt, pcov = curve_fit(gaussian, bin_centers, counts, p0=p0, maxfev=5000)
        perr = np.sqrt(np.diag(pcov))
        fit_success = True
    except:
        pass

    plt.hist(data_mirror, bins=100, color='lightgreen', edgecolor='white', label='Mirrored Data (-r, r)')
    
    if fit_success:
        x_smooth = np.linspace(min(data_mirror), max(data_mirror), 200)
        plt.plot(x_smooth, gaussian(x_smooth, *popt), 'r-', lw=2, label='Fit')
        
        info = (f"Center: {popt[1]:.2f} {unit}\n"
                f"Sigma: {abs(popt[2]):.2f} {unit}")
        plt.gca().text(0.05, 0.95, info, transform=plt.gca().transAxes,
                       verticalalignment='top', 
                       bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        
        text_lines.append(f"[{label} (Mirrored)]")
        text_lines.append(f"  Center Used: ({mu_x:.2f}, {mu_y:.2f}, {mu_z:.2f})")
        text_lines.append(f"  Sigma (Resolution): {abs(popt[2]):.4f} +/- {perr[2]:.4f} {unit}\n")
    else:
        text_lines.append(f"[{label}]")
        text_lines.append("  Fit Failed\n")

    plt.title(f'{label} Distribution (Mirrored)')
    plt.xlabel(f'Distance from Fit Center [{unit}]')
    plt.ylabel('Events')
    plt.legend()
    plt.grid(True, linestyle=':', alpha=0.6)
    plt.tight_layout()
    
    save_path = os.path.join(output_dir, f"{base_name}_hist_radius_mirror.pdf")
    plt.savefig(save_path)
    plt.close()

# ------------------------------------------------------------------
# メイン処理: 1つのCSVファイルを解析
# ------------------------------------------------------------------
def process_csv_file(file_path):
    print(f"Processing: {file_path}")

    target_dir = os.path.dirname(file_path)
    base_name = os.path.splitext(os.path.basename(file_path))[0]
    output_dir = os.path.join(target_dir, "reconst_images")
    
    if not os.path.exists(output_dir):
        try:
            os.makedirs(output_dir)
        except OSError:
            return

    try:
        df = pd.read_csv(file_path)
        
        # 最低限必要なカラムチェック
        required = ['fit_x', 'fit_y', 'fit_z', 't_light', 'chi2']
        if not all(c in df.columns for c in required):
            print(f"Skipped {file_path}: Missing columns")
            return

        fit_results_text = [f"=== Fit Results for {base_name} ===\n"]
        
        # 1. 1次元ヒストグラム (戻り値としてフィット中心を受け取る)
        mu_x, ok_x = process_1d_fit(df, 'fit_x', 'Position X', 'cm', output_dir, base_name, fit_results_text)
        mu_y, ok_y = process_1d_fit(df, 'fit_y', 'Position Y', 'cm', output_dir, base_name, fit_results_text)
        mu_z, ok_z = process_1d_fit(df, 'fit_z', 'Position Z', 'cm', output_dir, base_name, fit_results_text)
        mu_t, ok_t = process_1d_fit(df, 't_light', 'Time', 'ns', output_dir, base_name, fit_results_text)
        
        # [追加] 振幅A (カラムが存在し、かつ有効な値がある場合のみフィット)
        # 'A' がCSVにあるか確認
        if 'A' in df.columns:
            mu_A, ok_A = process_1d_fit(df, 'A', 'Amplitude A', 'arb.', output_dir, base_name, fit_results_text)
        else:
            fit_results_text.append("[Amplitude A]")
            fit_results_text.append("  Column 'A' not found\n")

        # 2. 距離rのヒストグラム (X,Y,Zのフィット中心を渡す)
        process_radial_fit(df, (mu_x, mu_y, mu_z), output_dir, base_name, fit_results_text)

        # 3. カイ二乗分布
        plt.figure(figsize=(8, 6))
        plt.hist(df['chi2'], bins=50, color='orange', edgecolor='black')
        plt.title('Chi-square Distribution')
        plt.xlabel(r'$\chi^2$')
        plt.ylabel('Events')
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f"{base_name}_chi2.pdf"))
        plt.close()

        # 4. 2次元ヒストグラム
        pairs = [
            ('fit_x', 'fit_y', 'X', 'Y', mu_x, mu_y, ok_x, ok_y),
            ('fit_y', 'fit_z', 'Y', 'Z', mu_y, mu_z, ok_y, ok_z),
            ('fit_x', 'fit_z', 'X', 'Z', mu_x, mu_z, ok_x, ok_z)
        ]

        for col1, col2, lab1, lab2, m1, m2, s1, s2 in pairs:
            # -9999のデータを除外してプロット
            df_2d = df[(df[col1] > -9990) & (df[col2] > -9990)]
            if len(df_2d) == 0: continue

            # A. 広域
            plt.figure(figsize=(8, 6))
            h = plt.hist2d(df_2d[col1], df_2d[col2], bins=50, 
                           range=[[-400, 400], [-400, 400]], cmap='viridis', cmin=1)
            plt.title(f'{lab1}-{lab2} Distribution (Wide)')
            plt.xlabel(f'{lab1} [cm]')
            plt.ylabel(f'{lab2} [cm]')
            plt.colorbar(h[3], label='Counts')
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, f"{base_name}_2d_{lab1}{lab2}_wide.pdf"))
            plt.close()

            # B. 拡大 (フィット成功時のみ)
            if s1 and s2:
                plt.figure(figsize=(8, 6))
                range_zoom = [[m1 - 20, m1 + 20], [m2 - 20, m2 + 20]]
                h = plt.hist2d(df_2d[col1], df_2d[col2], bins=50, 
                               range=range_zoom, cmap='viridis', cmin=1)
                plt.title(f'{lab1}-{lab2} Distribution (Zoom: Mean +/- 20)')
                plt.xlabel(f'{lab1} [cm]')
                plt.ylabel(f'{lab2} [cm]')
                plt.colorbar(h[3], label='Counts')
                plt.tight_layout()
                plt.savefig(os.path.join(output_dir, f"{base_name}_2d_{lab1}{lab2}_zoom.pdf"))
                plt.close()

        # 5. 3次元散布図
        df_3d = df[(df['fit_x'] > -9990) & (df['fit_y'] > -9990) & (df['fit_z'] > -9990)]
        if len(df_3d) > 0:
            fig = plt.figure(figsize=(10, 8))
            ax = fig.add_subplot(111, projection='3d')
            sc = ax.scatter(df_3d['fit_x'], df_3d['fit_y'], df_3d['fit_z'], 
                            c=df_3d['t_light'], cmap='plasma', s=10, alpha=0.6)
            ax.set_title('3D Reconstruction Position')
            ax.set_xlabel('X [cm]')
            ax.set_ylabel('Y [cm]')
            ax.set_zlabel('Z [cm]')
            cbar = fig.colorbar(sc, ax=ax, shrink=0.6)
            cbar.set_label('Time [ns]')
            plt.savefig(os.path.join(output_dir, f"{base_name}_3d_scatter.pdf"))
            plt.close()

        # 結果テキスト保存
        txt_path = os.path.join(output_dir, f"{base_name}_fit_results.txt")
        with open(txt_path, "w") as f:
            f.writelines("\n".join(fit_results_text))
        
        print(f"Saved results to: {output_dir}")

    except Exception as e:
        print(f"Error processing {file_path}: {e}")

# ------------------------------------------------------------------
# メイン処理
# ------------------------------------------------------------------
if __name__ == "__main__":
    if len(sys.argv) < 2:
        print_usage(sys.argv[0])
        sys.exit(1)

    input_path = sys.argv[1]

    if os.path.isdir(input_path):
        search_pattern = os.path.join(input_path, "*.csv")
        csv_files = glob.glob(search_pattern)
        if not csv_files:
            print(f"No CSV files found in {input_path}")
        else:
            print(f"Found {len(csv_files)} CSV files.")
            for f in csv_files:
                process_csv_file(f)

    elif os.path.isfile(input_path):
        if input_path.endswith(".csv"):
            process_csv_file(input_path)
        else:
            print("Error: Not a CSV file.")
    else:
        print(f"Error: Invalid path: {input_path}")