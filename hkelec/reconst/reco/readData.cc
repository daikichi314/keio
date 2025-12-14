/*
 * id: readData.cc
 * Place: /home/daiki/keio/hkelec/reconst/reco/
 * Author: Gemini 3 Pro
 * Last Edit: 2025-12-12
 */
#include "readData.hh"

// ... (定数定義やreadPedestalsは変更なし) ...

// データ変換ロジック
PMTData DataReader::createPMTData() {
    PMTData data;
    data.eventID = b_eventID;
    data.ch = b_ch;
    data.time = b_time_diff;
    data.isHit = true; // 10. ファイルから読み込まれたデータは基本的にHitしたもの

    // ... (ペデスタル取得、電荷計算ロジックは変更なし) ...
    // charge計算部分
    // ...

    // 座標セット
    if (b_ch >= 0 && b_ch < 4) {
        data.x = PMT_POSITIONS[b_ch][0];
        data.y = PMT_POSITIONS[b_ch][1];
        data.z = PMT_POSITIONS[b_ch][2];
        data.dir_x = PMT_DIR[0];
        data.dir_y = PMT_DIR[1];
        data.dir_z = PMT_DIR[2];
    } else {
        data.x = 0; data.y = 0; data.z = 0;
    }

    return data;
}

// nextEvent等は変更なし