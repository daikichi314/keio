/*id: eventtree2hist.C*/
/*Place: hkpd@hkdaq:~/hkelec/analysis_macro/ */
/*Last Edit: 2025-10-14 Gemini (two-pass analysis, full histogram block)*/
/* (修正: 2025-10-27 Gemini (ns単位に変更) ) */
/* (修正: 2025-10-30 Gemini (ヒストグラムのビン幅を最小単位に指定) ) */
/* (修正: 2025-11-21 Gemini (トリガーのhgain >= 850 の条件を追加) ) */
/* (修正: 2026-01-15 Gemini (チャンネルごとのtime_cutを適用するように変更) ) */
/*eventtree.root to triggered data TTree and optional histograms*/
/*コンパイル可能*/

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TString.h>
#include <TMath.h>
#include <iostream>
#include <vector>
#include <set>
#include <map> // [追加] チャンネルごとの管理に使用

// 外部で定義されたHitクラスとMetaDataクラスのヘッダを読み込む
#include "RootInterface/Hit.h"
#include "RootInterface/MetaData.h"

void read_event_tree(TString input_file, TString output_file) {

    // --- 1. ファイルの準備 ---
    auto ifile = TFile::Open(input_file, "READ");
    if (!ifile || ifile->IsZombie()) {
        std::cerr << "Error: Could not open input file " << input_file << std::endl;
        return;
    }
    auto ofile = TFile::Open(output_file, "RECREATE");
    if (!ofile || ofile->IsZombie()) {
        std::cerr << "Error: Could not create output file " << output_file << std::endl;
        ifile->Close();
        return;
    }

    // --- 2. 入力TTreeの準備 ---
    auto tree = ifile->Get<TTree>("event");
    if (!tree) {
        std::cerr << "Error: Could not find TTree 'event' in " << input_file << std::endl;
        ifile->Close();
        ofile->Close();
        return;
    }
    std::vector<Hit>* v_hit = nullptr;
    std::vector<Hit>* v_trig = nullptr;
    tree->SetBranchAddress("NormalHits", &v_hit);
    tree->SetBranchAddress("TriggerHits", &v_trig);

    // メタデータを読み込んで表示する
    auto metadata = ifile->Get<MetaData>("metadata");
    if (metadata) {
        metadata->Print();
    } else {
        std::cout << "Warning: Could not find 'metadata' object in " << input_file << std::endl;
    }

    long nEvents = tree->GetEntries();

    // [追加設定] トリガー選別の閾値
    const double TRIGGER_HGAIN_THRESHOLD = 800.0;

    // ==============================================================================
    // --- 3. 1回目のスキャン： チャンネルごとに time_diff のピーク位置を特定する ---
    // ==============================================================================
    std::cout << "\n--- Pass 1: Finding time_diff peak per Channel (Trigger hgain >= " << TRIGGER_HGAIN_THRESHOLD << ") ---" << std::endl;

    // [変更] チャンネルごとのヒストグラムを保持するマップ
    std::map<int, TH1D*> prescan_hists;

    for (long iEvent = 0; iEvent < nEvents; ++iEvent) {
        tree->GetEntry(iEvent);
        if (iEvent % 10000 == 0) {
            std::cout << "Scanning event: " << iEvent << " / " << nEvents << std::endl;
        }

        if (v_trig->empty()) continue;

        // [追加箇所] 1. トリガーのhgainによる選別 (Pass 1)
        if (v_trig->at(0).hgain < TRIGGER_HGAIN_THRESHOLD) {
            continue; // 条件を満たさないイベントはスキップ
        }

        double trigger_time = v_trig->at(0).time;

        for (const Hit& hit : *v_hit) {
            double time_diff = hit.time - trigger_time;
            int ch = hit.channel;

            // 0 ns < time_diff < 1000 ns のヒットのみ対象
            if (time_diff > 0 && time_diff < 1000e-9) {
                // まだこのチャンネルのヒストグラムがなければ作成する
                if (prescan_hists.find(ch) == prescan_hists.end()) {
                    TString hname = Form("h_time_prescan_ch%d", ch);
                    TString htitle = Form("Pre-scan time peak Ch %d; Time Diff (s); Counts", ch);
                    // 0 ns から 1000 ns までを 1000 分割 => ビン幅 1 ns
                    prescan_hists[ch] = new TH1D(hname, htitle, 1000, 0, 1000e-9);
                }
                prescan_hists[ch]->Fill(time_diff);
            }
        }
    }

    // --- 各チャンネルのピーク位置とカット範囲を決定して保持する ---
    // キー: チャンネル番号, 値: pair<low, high>
    std::map<int, std::pair<double, double>> channel_time_cuts;

    std::cout << "\n--- Calculating Time Cuts per Channel ---" << std::endl;
    for (auto const& [ch, hist] : prescan_hists) {
        // データが少なすぎるチャンネル (10未満) は破棄する仕様
        if (hist->GetEntries() < 10) {
            std::cout << "Channel " << ch << ": Not enough entries (" << hist->GetEntries() << "). Skipping." << std::endl;
            continue; 
        }

        int peak_bin = hist->GetMaximumBin();
        double time_peak = hist->GetXaxis()->GetBinCenter(peak_bin);

        // ピーク値からカット範囲を決定 (ピーク ± 50 ns)
        // double time_cut_low  = time_peak - 50.0e-9;
        // double time_cut_high = time_peak + 50.0e-9;
        double time_cut_low  = time_peak - 8.0e-9;
        double time_cut_high = time_peak + 8.0e-9;

        channel_time_cuts[ch] = std::make_pair(time_cut_low, time_cut_high);

        std::cout << "Channel " << ch << ": Peak=" << time_peak * 1e9 
                  << " ns, Window=[" << time_cut_low * 1e9 << ", " << time_cut_high * 1e9 << "] ns" << std::endl;
    }


    // --- 4. 出力用TTreeの定義 ---
    TTree* new_tree = new TTree("processed_hits", "Processed Hit Data per Channel");
    int eventID;
    int ch;
    double hgain, lgain, tot;
    double tdc_diff, time_diff; // time_diff は ns 単位で保存されるように変更
    new_tree->Branch("eventID", &eventID, "eventID/I");
    new_tree->Branch("ch", &ch, "ch/I");
    new_tree->Branch("hgain", &hgain, "hgain/D");
    new_tree->Branch("lgain", &lgain, "lgain/D");
    new_tree->Branch("tot", &tot, "tot/D");
    new_tree->Branch("tdc_diff", &tdc_diff, "tdc_diff/D");
    new_tree->Branch("time_diff", &time_diff, "time_diff/D"); 

    // ==============================================================================
    // --- 5. 2回目のスキャン：決定した時間範囲でヒットを選択し、TTreeを作成 ---
    // ==============================================================================
    std::cout << "\n--- Pass 2: Processing hits with PER-CHANNEL time cut (Trigger hgain >= " << TRIGGER_HGAIN_THRESHOLD << ") ---" << std::endl;
    std::set<int> unique_channels; // 実際にTTreeに保存されたチャンネルを記録

    for (long iEvent = 0; iEvent < nEvents; ++iEvent) {
        tree->GetEntry(iEvent);
        if (iEvent % 10000 == 0) {
            std::cout << "Processing event: " << iEvent << " / " << nEvents << std::endl;
        }

        if (v_trig->empty()) continue;

        // [追加箇所] 2. トリガーのhgainによる選別 (Pass 2)
        if (v_trig->at(0).hgain < TRIGGER_HGAIN_THRESHOLD) {
            continue;
        }

        double trigger_tdc = v_trig->at(0).tdc;
        double trigger_time = v_trig->at(0).time;
        eventID = iEvent;

        for (const Hit& hit : *v_hit) {
            ch = hit.channel;

            // [変更] チャンネルごとのカット設定が存在するか確認
            // (Pass 1でデータ不足で破棄されたチャンネルはここで除外される)
            if (channel_time_cuts.find(ch) == channel_time_cuts.end()) {
                continue;
            }

            double local_cut_low = channel_time_cuts[ch].first;
            double local_cut_high = channel_time_cuts[ch].second;

            time_diff = hit.time - trigger_time; // この時点では (s) 単位

            // [変更] チャンネル個別の時間範囲でヒットを選択 (s 単位で比較)
            if (time_diff < local_cut_low || time_diff > local_cut_high) {
                continue;
            }

            hgain = hit.hgain;
            lgain = hit.lgain;
            tot = hit.tot;
            tdc_diff = hit.tdc - trigger_tdc;

            // 1. TTreeに保存する直前に (s) から (ns) へ単位を変換
            time_diff *= 1e9;

            new_tree->Fill();
            unique_channels.insert(ch);
        }
    }

    // ==============================================================================
    // --- 6. (オプション) TTreeからチャンネルごとのヒストグラムを作成 ---
    // ==============================================================================
    // このブロックのコメントアウトを外して使用してください。
    // /*
    if(!unique_channels.empty()){
        std::cout << "\nFound " << unique_channels.size() << " unique channels. Creating histograms..." << std::endl;
        ofile->cd();
        for (int ch_num : unique_channels) {
            TString selection = Form("ch==%d", ch_num);
            Long64_t n_selected = 0;

            // [変更箇所] --- hgain (ビン幅 1.0 で自動範囲設定) ---
            new_tree->Draw("hgain", selection, "goff");
            n_selected = new_tree->GetSelectedRows();
            if (n_selected > 0) {
                double min_val = TMath::MinElement(n_selected, new_tree->GetV1());
                double max_val = TMath::MaxElement(n_selected, new_tree->GetV1());
                double margin = (max_val - min_val) * 0.05;
                if (margin == 0) margin = 1;

                const double bin_width = 1.0;
                double xlow = TMath::Floor(min_val - margin);
                double xup = TMath::Ceil(max_val + margin);
                int nbins = TMath::Nint((xup - xlow) / bin_width);
                if (nbins <= 0) nbins = 1;

                TH1D* h_hgain = new TH1D(Form("h_hgain_ch%d", ch_num), Form("High Gain ADC Ch %d", ch_num), nbins, xlow + bin_width/2, xup + bin_width/2);
                new_tree->Draw(Form("hgain>>%s", h_hgain->GetName()), selection, "goff");
            }

            // [変更箇所] --- lgain (ビン幅 1.0 で自動範囲設定) ---
            new_tree->Draw("lgain", selection, "goff");
            n_selected = new_tree->GetSelectedRows();
            if (n_selected > 0) {
                double min_val = TMath::MinElement(n_selected, new_tree->GetV1());
                double max_val = TMath::MaxElement(n_selected, new_tree->GetV1());
                double margin = (max_val - min_val) * 0.05;
                if (margin == 0) margin = 1;

                const double bin_width = 1.0;
                double xlow = TMath::Floor(min_val - margin);
                double xup = TMath::Ceil(max_val + margin);
                int nbins = TMath::Nint((xup - xlow) / bin_width);
                if (nbins <= 0) nbins = 1;

                TH1D* h_lgain = new TH1D(Form("h_lgain_ch%d", ch_num), Form("Low Gain ADC Ch %d", ch_num), nbins, xlow + bin_width/2, xup + bin_width/2);
                new_tree->Draw(Form("lgain>>%s", h_lgain->GetName()), selection, "goff");
            }

            // [変更箇所] --- tot (ビン幅 1.0 で自動範囲設定) ---
            new_tree->Draw("tot", selection, "goff");
            n_selected = new_tree->GetSelectedRows();
            if (n_selected > 0) {
                double min_val = TMath::MinElement(n_selected, new_tree->GetV1());
                double max_val = TMath::MaxElement(n_selected, new_tree->GetV1());
                double margin = (max_val - min_val) * 0.05;
                if (margin == 0) margin = 1;

                const double bin_width = 1.0;
                double xlow = TMath::Floor(min_val - margin);
                double xup = TMath::Ceil(max_val + margin);
                int nbins = TMath::Nint((xup - xlow) / bin_width);
                if (nbins <= 0) nbins = 1;

                TH1D* h_tot = new TH1D(Form("h_tot_ch%d", ch_num), Form("Time over Threshold Ch %d", ch_num), nbins, xlow + bin_width/2, xup + bin_width/2);
                new_tree->Draw(Form("tot>>%s", h_tot->GetName()), selection, "goff");
            }

            // [変更箇所] --- tdc_diff (ビン幅 1.0 で自動範囲設定) ---
            new_tree->Draw("tdc_diff", selection, "goff");
            n_selected = new_tree->GetSelectedRows();
            if (n_selected > 0) {
                double min_val = TMath::MinElement(n_selected, new_tree->GetV1());
                double max_val = TMath::MaxElement(n_selected, new_tree->GetV1());
                double margin = (max_val - min_val) * 0.05;
                if (margin == 0) margin = 1;

                const double bin_width = 1.0;
                double xlow = TMath::Floor(min_val - margin);
                double xup = TMath::Ceil(max_val + margin);
                int nbins = TMath::Nint((xup - xlow) / bin_width);
                if (nbins <= 0) nbins = 1;

                TH1D* h_tdc_diff = new TH1D(Form("h_tdc_diff_ch%d", ch_num), Form("TDC - Trigger TDC Ch %d", ch_num), nbins, xlow+bin_width/2, xup+bin_width/2);
                new_tree->Draw(Form("tdc_diff>>%s", h_tdc_diff->GetName()), selection, "goff");
            }

            // [変更箇所] --- time_diff (ビン幅 0.25 ns で自動範囲設定) ---
            new_tree->Draw("time_diff", selection, "goff");
            n_selected = new_tree->GetSelectedRows();
            if (n_selected > 0) {
                double min_val = TMath::MinElement(n_selected, new_tree->GetV1());
                double max_val = TMath::MaxElement(n_selected, new_tree->GetV1());
                double margin = (max_val - min_val) * 0.05;

                if (margin == 0) margin = 1.0;

                const double bin_width = 0.25;
                double xlow_raw = min_val - margin;
                double xup_raw = max_val + margin;

                double xlow = TMath::Floor(xlow_raw / bin_width) * bin_width;
                double xup = TMath::Ceil(xup_raw / bin_width) * bin_width;

                int nbins = TMath::Nint((xup - xlow) / bin_width);
                if (nbins <= 0) nbins = 1;

                TH1D* h_time_diff = new TH1D(Form("h_time_diff_ch%d", ch_num), Form("Time - Trigger Time (ns) Ch %d", ch_num), nbins, xlow+bin_width/2, xup+bin_width/2);
                new_tree->Draw(Form("time_diff>>%s", h_time_diff->GetName()), selection, "goff");
            }
        }
    }
    // */
    // ==============================================================================


    // --- 7. ファイルの書き込みとクローズ ---
    std::cout << "\nWriting TTree and histograms to " << output_file << std::endl;
    
    // [変更] チャンネルごとのプレキャン用ヒストグラムも保存する
    for (auto const& [ch, hist] : prescan_hists) {
        hist->Write();
    }
    
    ofile->Write(); // new_treeと、もし作成されていれば他のヒストグラムを保存
    ofile->Close();
    ifile->Close();
}

// --- スタンドアロン実行のためのmain関数 ---
int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <input_file.root> <output_file.root>" << std::endl;
        return 1;
    }
    read_event_tree(argv[1], argv[2]);
    return 0;
}