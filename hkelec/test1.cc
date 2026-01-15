#include "readData.hh"

// ADCから電荷(pC)への変換係数
const double K_HGAIN = 0.073; // High Gain用 pC/ADC
const double K_LGAIN = 0.599; // Low Gain用 pC/ADC
const double SATURATION_THRESHOLD = 4000.0; // High Gainが飽和したとみなす閾値

// ペデスタル読み込み関数の実装
int readPedestals(const std::string &filename, std::map<int, PedestalData> &pedestalMap) {
    std::ifstream infile(filename); // ファイルを開く
    if (!infile) { // 開けなかった場合
        std::cerr << "Error: Cannot open pedestal file " << filename << std::endl;
        return 1;
    }

    std::string line;
    // ファイルを1行ずつ読み込む
    while (std::getline(infile, line)) {
        if (line.empty() || line[0] == '#') continue; // 空行やコメント行(#で開始)はスキップ
        
        // カンマ区切りCSVをスペース区切りに変換して扱いやすくする
        for (char &c : line) if (c == ',') c = ' ';
        
        std::istringstream iss(line); // 文字列ストリームを作成
        int ch;
        std::string type;
        double mean, err;
        
        // チャンネル、タイプ、平均値、誤差を読み取る
        if (!(iss >> ch >> type >> mean >> err)) continue; // 読み取り失敗ならスキップ

        // タイプに応じてマップに値をセット
        if (type == "hgain") {
            pedestalMap[ch].hgain_mean = mean;
        } else if (type == "lgain") {
            pedestalMap[ch].lgain_mean = mean;
        }
    }
    return 0; // 成功
}

// DataReaderコンストラクタの実装
DataReader::DataReader(const std::string &filename, const std::map<int, PedestalData> &pedMap) 
    : pedestalMap(pedMap), currentEntry(0), hasBufferedHit(false) { // 初期化子リスト
    
    file = TFile::Open(filename.c_str()); // ROOTファイルを開く
    if (!file || file->IsZombie()) { // 開けなかった、または壊れている場合
        std::cerr << "Error: Cannot open ROOT file " << filename << std::endl;
        tree = nullptr;
        return;
    }

    tree = dynamic_cast<TTree*>(file->Get("processed_hits")); // TTreeを取得
    if (!tree) {
        std::cerr << "Error: Cannot find 'processed_hits' tree" << std::endl;
        return;
    }

    // ブランチのアドレスをメンバ変数にセット
    tree->SetBranchAddress("eventID", &b_eventID);
    tree->SetBranchAddress("ch", &b_ch);
    tree->SetBranchAddress("hgain", &b_hgain);
    tree->SetBranchAddress("lgain", &b_lgain);
    tree->SetBranchAddress("tot", &b_tot);
    tree->SetBranchAddress("time_diff", &b_time_diff);

    nEntries = tree->GetEntries(); // 総エントリー数を取得
}

// DataReaderデストラクタの実装
DataReader::~DataReader() {
    if (file) file->Close(); // ファイルが開いていれば閉じる
}

// 1ヒット分のデータを作成する内部関数
PMTData DataReader::createPMTData() {
    PMTData data;
    data.eventID = b_eventID; // イベントIDをセット
    data.ch = b_ch; // チャンネルをセット
    
    // processed_hitsのtime_diffブランチの値を時刻としてセット
    data.time = b_time_diff; 

    // ペデスタル値を取得
    double ped_h = 0, ped_l = 0;
    auto it = pedestalMap.find(b_ch); // マップからチャンネルを検索
    if (it != pedestalMap.end()) { // 見つかった場合
        ped_h = it->second.hgain_mean;
        ped_l = it->second.lgain_mean;
    }

    // High Gainが飽和しているかチェック
    if (b_hgain >= SATURATION_THRESHOLD) {
        // 飽和時はLow Gainを使用: (ADC - ペデスタル) * 係数
        data.charge = (b_lgain - ped_l) * K_LGAIN;
    } else {
        // 通常時はHigh Gainを使用
        data.charge = (b_hgain - ped_h) * K_HGAIN;
    }

    // 負の電荷は物理的にあり得ないので0にする（ノイズ対策）
    if (data.charge < 0) data.charge = 0;

    // PMTの位置情報を定数配列から取得してセット
    if (b_ch >= 0 && b_ch < 4) {
        data.x = PMT_POSITIONS[b_ch][0];
        data.y = PMT_POSITIONS[b_ch][1];
        data.z = PMT_POSITIONS[b_ch][2];
        data.dir_x = PMT_DIR[0];
        data.dir_y = PMT_DIR[1];
        data.dir_z = PMT_DIR[2];
    } else {
        // 想定外のチャンネルの場合 (0,0,0) にする
        data.x = 0; data.y = 0; data.z = 0;
    }

    return data;
}

// 次のイベントを取得する関数の実装
bool DataReader::nextEvent(std::vector<PMTData> &eventHits) {
    eventHits.clear(); // ベクトルをクリア
    
    if (!tree) return false; // ツリーが無ければ終了

    // 前回のループで「次のイベントの先頭」を読み込んでバッファしていた場合
    if (hasBufferedHit) {
        // バッファの値を現在のブランチ変数に戻す
        b_eventID = buf_eventID;
        b_ch = buf_ch;
        b_hgain = buf_hgain;
        b_lgain = buf_lgain;
        b_tot = buf_tot;
        b_time_diff = buf_time_diff;
        
        // データを作成して追加
        eventHits.push_back(createPMTData());
        hasBufferedHit = false; // バッファ使用済み
    }

    int currentEventID = -1;
    if (!eventHits.empty()) {
        currentEventID = eventHits[0].eventID; // 現在処理中のイベントIDを記録
    }

    // TTreeを読み進めるループ
    while (currentEntry < nEntries) {
        tree->GetEntry(currentEntry); // エントリーを読み込む
        currentEntry++; // カウンタを進める

        // まだヒットがない場合、このヒットがイベントの開始
        if (eventHits.empty()) {
            currentEventID = b_eventID;
        }

        // イベントIDが変わった場合（＝次のイベントのデータが来た）
        if (b_eventID != currentEventID) {
            hasBufferedHit = true; // このデータをバッファに保存
            buf_eventID = b_eventID;
            buf_ch = b_ch;
            buf_hgain = b_hgain;
            buf_lgain = b_lgain;
            buf_tot = b_tot;
            buf_time_diff = b_time_diff;
            return true; // 現在のイベント収集完了としてリターン
        }

        // 現在のイベントと同じIDなら、リストに追加して次へ
        eventHits.push_back(createPMTData());
    }

    // ファイルの最後まで読み終わった場合、リストにデータがあればtrueを返す
    return !eventHits.empty();
}