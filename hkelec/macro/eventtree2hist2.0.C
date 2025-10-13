/*id: eventtree2hist.C*/
/*Place: hkpd@hkdaq:~/hkelec/analysis_macro/ */
/*Last Edit: 2025-10-10 Gemini*/
/*eventtree.root to triggered data*/
/*コンパイル可能*/

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TString.h>
#include <iostream>
#include <vector>

// 外部のヘッダファイルをインクルード
#include "RootInterface/Hit.h"
#include "RootInterface/MetaData.h"

void read_event_tree(TString input_file, TString output_file) {

    auto ifile = TFile::Open(input_file, "READ");
    auto ofile = TFile::Open(output_file, "RECREATE");

    /* 出力用の新しいTTreeを作成 */
    TTree* new_tree = new TTree("processed_hits", "Processed Hit Data per Channel");

    /* TTreeに保存するための変数を定義 */
    int eventID;
    int ch;
    double hgain, lgain, tot;
    double tdc_diff, time_diff;

    /* Branchを作成し、変数と紐付ける */
    new_tree->Branch("eventID", &eventID, "eventID/I");
    new_tree->Branch("ch", &ch, "ch/I");
    new_tree->Branch("hgain", &hgain, "hgain/D");
    new_tree->Branch("lgain", &lgain, "lgain/D");
    new_tree->Branch("tot", &tot, "tot/D");
    new_tree->Branch("tdc_diff", &tdc_diff, "tdc_diff/D");
    new_tree->Branch("time_diff", &time_diff, "time_diff/D");


    /* Load original tree */
    auto tree = ifile->Get<TTree>("event");
    std::vector<Hit>* v_hit = nullptr;
    std::vector<Hit>* v_trig = nullptr;
    tree->SetBranchAddress("NormalHits", &v_hit);
    tree->SetBranchAddress("TriggerHits", &v_trig);

    /* Print metadata */
    auto metadata = ifile->Get<MetaData>("metadata");
    metadata->Print();

    /* Process events */
    long nEvents = tree->GetEntries();
    for (long iEvent = 0; iEvent < nEvents; ++iEvent) {
        tree->GetEntry(iEvent);
        std::cout << "Processing Event: " << iEvent << std::endl;

        // イベントにTrigger Hitがある場合のみ、その時間・TDCを基準とする
        if (v_trig->empty()) {
            continue; // Trigger Hitがないイベントはスキップ
        }
        double trigger_tdc = v_trig->at(0).tdc;
        double trigger_time = v_trig->at(0).time;
        
        eventID = iEvent; // イベント番号を保存

        for (const Hit& hit : *v_hit) {
            // 変数に値をセット
            ch = hit.channel;
            hgain = hit.hgain;
            lgain = hit.lgain;
            tot = hit.tot;
            tdc_diff = hit.tdc - trigger_tdc;
            time_diff = hit.time - trigger_time;

            // TTreeに1エントリー分を書き込む
            new_tree->Fill();
        }
    }

    ofile->cd(); // 出力ファイルに移動
    new_tree->Write(); // TTreeを書き込む
    ofile->Close();
    ifile->Close();
}

// スタンドアロン実行のためのメイン関数
int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <input_file.root> <output_file.root>" << std::endl;
        return 1;
    }
    read_event_tree(argv[1], argv[2]);
    return 0;
}