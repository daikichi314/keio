/*id: eventtree2hist.C*/
/*Place: hkpd@hkdaq:~/hkelec/analysis_macro/ */
/*Last Edit: 2025-10-13 Gemini (memory-efficient version)*/
/*eventtree.root to triggered data TTree and optional histograms*/
/*コンパイル可能*/

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TString.h>
#include <TMath.h> // TMath::MinElement, MaxElement を使うために追加
#include <iostream>
#include <vector>
#include <set>

// 外部で定義されたHitクラスとMetaDataクラスのヘッダを読み込む
#include "RootInterface/Hit.h"
#include "RootInterface/MetaData.h"

void read_event_tree(TString input_file, TString output_file) {

    // --- 1. ファイルの準備 ---
    // 入力ROOTファイルを開く
    auto ifile = TFile::Open(input_file, "READ");
    // 出力ROOTファイルを新規作成する（同名ファイルがあれば上書き）
    auto ofile = TFile::Open(output_file, "RECREATE");

    // --- 2. 出力用TTreeの定義 ---
    // "processed_hits"という名前の新しいTTreeオブジェクトを作成
    TTree* new_tree = new TTree("processed_hits", "Processed Hit Data per Channel");

    // TTreeに格納するための変数を宣言
    int eventID;
    int ch;
    double hgain, lgain, tot;
    double tdc_diff, time_diff;

    // TTreeに「ブランチ」を作成し、上で宣言した変数のアドレスを紐付ける
    // これにより、new_tree->Fill()が呼ばれるたびに変数の現在の値が書き込まれる
    new_tree->Branch("eventID", &eventID, "eventID/I");
    new_tree->Branch("ch", &ch, "ch/I");
    new_tree->Branch("hgain", &hgain, "hgain/D");
    new_tree->Branch("lgain", &lgain, "lgain/D");
    new_tree->Branch("tot", &tot, "tot/D");
    new_tree->Branch("tdc_diff", &tdc_diff, "tdc_diff/D");
    new_tree->Branch("time_diff", &time_diff, "time_diff/D");

    // --- 3. 入力TTreeの読み込み ---
    // 入力ファイルから"event"という名前のTTreeを取得
    auto tree = ifile->Get<TTree>("event");
    // Hitオブジェクトのベクトルへのポインタを準備
    std::vector<Hit>* v_hit = nullptr;
    std::vector<Hit>* v_trig = nullptr;
    // TTreeのブランチとポインタを紐付ける
    tree->SetBranchAddress("NormalHits", &v_hit);
    tree->SetBranchAddress("TriggerHits", &v_trig);

    // メタデータを読み込んで表示する
    auto metadata = ifile->Get<MetaData>("metadata");
    metadata->Print();

    // --- 4. イベントループによるデータ処理 ---
    // ユニークなチャンネル番号を効率的に保存するためのsetコンテナ
    std::set<int> unique_channels;

    // 全イベント数でループを実行
    long nEvents = tree->GetEntries();
    for (long iEvent = 0; iEvent < nEvents; ++iEvent) {
        // iEvent番目のデータをメモリに読み込む
        tree->GetEntry(iEvent);
        // 1000イベントごとに進捗状況をターミナルに表示
        if (iEvent % 1000 == 0) {
            std::cout << "Processing Event: " << iEvent << " / " << nEvents << std::endl;
        }

        // トリガーヒットがなければ、そのイベントは処理しない
        if (v_trig->empty()) {
            continue;
        }
        // 最初のトリガーヒットの情報を基準値とする
        double trigger_tdc = v_trig->at(0).tdc;
        double trigger_time = v_trig->at(0).time;
        
        // 現在のイベント番号を保存
        eventID = iEvent;

        // イベント内の全ノーマルヒットでループ
        for (const Hit& hit : *v_hit) {
            // TTreeに書き込むための各変数に値を代入
            ch = hit.channel;
            hgain = hit.hgain;
            lgain = hit.lgain;
            tot = hit.tot;
            tdc_diff = hit.tdc - trigger_tdc;
            time_diff = hit.time - trigger_time;
            // ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
            // ★ 以下を必要に応じて編集: ヒットを選択する条件を追加                                  ★
            // ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
            if (time_diff >= 2.5e-7||time_diff <= 1.5e-7/*ケーブルの長さに依存*/) {
                continue; // 条件を満たさないヒットはここで処理を中断し、次のヒットへ移る
            }
            
            // このヒットの情報をTTreeに1エントリーとして書き込む
            new_tree->Fill();

            // このヒットのチャンネル番号をsetに追加（重複は自動で無視される）
            unique_channels.insert(ch);
        }
    }

    // ==============================================================================
    // (オプション) TTreeからチャンネルごとのヒストグラムを作成する
    // ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
    // ★ このブロックのコメントアウトを必要に応じて外して使用してください。                   ★
    // ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
    // /*
    if(!unique_channels.empty()){
        std::cout << "Found " << unique_channels.size() << " unique channels. Creating histograms..." << std::endl;
        ofile->cd(); // ヒストグラムが正しいファイルに作成されるようにディレクトリを移動

        // 見つかった全てのユニークなチャンネル番号でループ
        for (int ch_num : unique_channels) {
            
            TString selection = Form("ch==%d", ch_num);
            Long64_t n_selected = 0;

            // --- hgain ---
            new_tree->Draw("hgain", selection, "goff");
            n_selected = new_tree->GetSelectedRows();
            if (n_selected > 0) {
                double min_val = TMath::MinElement(n_selected, new_tree->GetV1());
                double max_val = TMath::MaxElement(n_selected, new_tree->GetV1());
                double margin = (max_val - min_val) * 0.05;
                if (margin == 0) margin = 1;
                
                TH1D* h_hgain = new TH1D(Form("h_hgain_ch%d", ch_num), Form("High Gain ADC Ch %d", ch_num), 400, min_val - margin, max_val + margin);
                new_tree->Draw(Form("hgain>>%s", h_hgain->GetName()), selection, "goff");
            }

            // --- lgain ---
            new_tree->Draw("lgain", selection, "goff");
            n_selected = new_tree->GetSelectedRows();
            if (n_selected > 0) {
                double min_val = TMath::MinElement(n_selected, new_tree->GetV1());
                double max_val = TMath::MaxElement(n_selected, new_tree->GetV1());
                double margin = (max_val - min_val) * 0.05;
                if (margin == 0) margin = 1;

                TH1D* h_lgain = new TH1D(Form("h_lgain_ch%d", ch_num), Form("Low Gain ADC Ch %d", ch_num), 400, min_val - margin, max_val + margin);
                new_tree->Draw(Form("lgain>>%s", h_lgain->GetName()), selection, "goff");
            }

            // --- tot ---
            // new_tree->Draw("tot", selection, "goff");
            // n_selected = new_tree->GetSelectedRows();
            // if (n_selected > 0) {
            //     double min_val = TMath::MinElement(n_selected, new_tree->GetV1());
            //     double max_val = TMath::MaxElement(n_selected, new_tree->GetV1());
            //     double margin = (max_val - min_val) * 0.05;
            //     if (margin == 0) margin = 1;

            //     TH1D* h_tot = new TH1D(Form("h_tot_ch%d", ch_num), Form("Time over Threshold Ch %d", ch_num), 400, min_val - margin, max_val + margin);
            //     new_tree->Draw(Form("tot>>%s", h_tot->GetName()), selection, "goff");
            // }

            // --- tdc_diff ---
            new_tree->Draw("tdc_diff", selection, "goff");
            n_selected = new_tree->GetSelectedRows();
            if (n_selected > 0) {
                double min_val = TMath::MinElement(n_selected, new_tree->GetV1());
                double max_val = TMath::MaxElement(n_selected, new_tree->GetV1());
                double margin = (max_val - min_val) * 0.05;
                if (margin == 0) margin = 1;
                
                TH1D* h_tdc_diff = new TH1D(Form("h_tdc_diff_ch%d", ch_num), Form("TDC - Trigger TDC Ch %d", ch_num), 400, min_val - margin, max_val + margin);
                new_tree->Draw(Form("tdc_diff>>%s", h_tdc_diff->GetName()), selection, "goff");
            }

            // --- time_diff ---
            new_tree->Draw("time_diff", selection, "goff");
            n_selected = new_tree->GetSelectedRows();
            if (n_selected > 0) {
                double min_val = TMath::MinElement(n_selected, new_tree->GetV1());
                double max_val = TMath::MaxElement(n_selected, new_tree->GetV1());
                double margin = (max_val - min_val) * 0.05;
                if (margin == 0) margin = 1.0e-9; // time_diffは値が小さいのでマージンも小さくする

                TH1D* h_time_diff = new TH1D(Form("h_time_diff_ch%d", ch_num), Form("Time - Trigger Time (s) Ch %d", ch_num), 400, min_val - margin, max_val + margin);
                new_tree->Draw(Form("time_diff>>%s", h_time_diff->GetName()), selection, "goff");
            }
        }
    }
    // */
    // ==============================================================================


    // --- 5. ファイルの書き込みとクローズ ---
    std::cout << "Writing TTree and histograms to " << output_file << std::endl;
    // メモリ上にあるTTreeやヒストグラムを全て出力ファイルに書き込む
    ofile->Write();
    // ファイルを閉じる
    ofile->Close();
    ifile->Close();
}

// --- スタンドアロン実行のためのmain関数 ---
// シェルから ./eventtree2hist <input> <output> のように実行された場合に呼ばれる
int main(int argc, char* argv[]) {
    // 引数の数が2つでない場合は使い方を表示して終了
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <input_file.root> <output_file.root>" << std::endl;
        return 1;
    }
    // メインの解析関数を呼び出す
    read_event_tree(argv[1], argv[2]);
    return 0;
}