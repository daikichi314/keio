#!/bin/bash

# 出力ファイルヘッダーの作成
echo "Voltage Gain_ch00 Gain_ch01 Gain_ch02 Gain_ch03" > /home/daiki/lab/data/20250809/HV_vs_gain.txt

# 電圧の初期値と増分
voltage=1500
increment=50

# ファイル名の連番と電圧を対応させてループ
for i in {6..18}
do
    # XXの値を2桁にフォーマット
    XX=$(printf "%02d" $i)

    # ファイル名を定義（0を追加）
    input_file="test0807-0${XX}_hithist.root"
    output_fit_file="test0807-0${XX}_fit01.root"
    output_txt_file="test0807-0${XX}_hithist.txt"

    echo "Processing ${input_file}..."

    # C++プログラムを実行してtxtファイルを生成 
    # USAGE : ./simplefit (Input Root filename) (Output Root filename) (optional: num of peak finding) (optional, 1:Add in existing root, 0:Recreate new (default) )
    # n_peak_finding = 2 (default)
    # Add_in_existing_root = 0 (default)
    ./simplefit ${input_file} ${output_fit_file} 2 0

    # 生成されたtxtファイルからゲインの値を抽出
    # BBB RESULT: の後に続く6番目の列を抽出
    gain_ch00=$(awk 'NR==1' /home/daiki/lab/data/20250809/${output_txt_file} | awk '{print $6}')
    gain_ch01=$(awk 'NR==2' /home/daiki/lab/data/20250809/${output_txt_file} | awk '{print $6}')
    gain_ch02=$(awk 'NR==3' /home/daiki/lab/data/20250809/${output_txt_file} | awk '{print $6}')
    gain_ch03=$(awk 'NR==4' /home/daiki/lab/data/20250809/${output_txt_file} | awk '{print $6}')

    # 電圧とゲインの値をHV_vs_gain.txtに追記
    echo "${voltage} ${gain_ch00} ${gain_ch01} ${gain_ch02} ${gain_ch03}" | sed 's/,//g' >> /home/daiki/lab/data/20250809/HV_vs_gain.txt

    # 次のループのために電圧を更新
    voltage=$((voltage + increment))
done

echo "Processing complete. The result is saved in /home/daiki/lab/data/20250809/HV_vs_gain.txt"
