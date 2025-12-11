/*id: eventtree2hist.C*/
    /*Place: hkpd@hkdaq:~/hkelec/analysis_macro/ */
    /*Last Edit: 2025-10-14 Gemini (two-pass analysis, full histogram block)*/
    /* (修正: 2025-10-27 Gemini (ns単位に変更) ) */
    /* (修正: 2025-10-30 Gemini (ヒストグラムのビン幅を最小単位に指定) ) */
    /* (修正: 2025-11-21 Gemini (トリガーのhgain >= 850 の条件を追加) ) */
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

        // [追加設定] トリガー選別の閾値////////////////////////////////////////////////////////////
        const double TRIGGER_HGAIN_THRESHOLD = 800.0;

        // ==============================================================================
        // --- 3. 1回目のスキャン： time_diff のピーク位置を特定する ---
        // ==============================================================================
        std::cout << "\n--- Pass 1: Finding time_diff peak (Trigger hgain >= " << TRIGGER_HGAIN_THRESHOLD << ") ---" << std::endl;
        
        // ピーク特定用の仮ヒストグラム (0 ns から 1000 ns までを 1000 分割 => ビン幅 1 ns)
        TH1D* h_time_prescan = new TH1D("h_time_prescan", "Pre-scan for time peak; Time Difference (s); Counts", 1000, 0, 1000e-9);

        for (long iEvent = 0; iEvent < nEvents; ++iEvent) {
            tree->GetEntry(iEvent);
            if (iEvent % 10000 == 0) {
                std::cout << "Scanning event: " << iEvent << " / " << nEvents << std::endl;
            }

            if (v_trig->empty()) continue;

            // [追加箇所] 1. トリガーのhgainによる選別 (Pass 1)
            // トリガーヒットはTriggerHitsの0番目にあると仮定
            if (v_trig->at(0).hgain < TRIGGER_HGAIN_THRESHOLD) {
                continue; // 条件を満たさないイベントはスキップ
            }

            double trigger_time = v_trig->at(0).time;

            for (const Hit& hit : *v_hit) {
                double time_diff = hit.time - trigger_time;
                // 0 ns < time_diff < 1000 ns のヒットのみでヒストグラムを作成
                if (time_diff > 0 && time_diff < 1000e-9) {
                    h_time_prescan->Fill(time_diff);
                }
            }
        }

        // ピークのビンを探し、その中心値を time_peak として取得
        int peak_bin = h_time_prescan->GetMaximumBin();
        double time_peak = h_time_prescan->GetXaxis()->GetBinCenter(peak_bin);

        // ピーク値からカット範囲を決定 (ピーク ± 50 ns)
        double time_cut_low  = time_peak - 50.0e-9;
        double time_cut_high = time_peak + 50.0e-9;
        
        std::cout << "Time peak found at: " << time_peak * 1e9 << " ns" << std::endl;
        std::cout << "Applying new time window: " << time_cut_low * 1e9 << " ns to " << time_cut_high * 1e9 << " ns" << std::endl;

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
        new_tree->Branch("time_diff", &time_diff, "time_diff/D"); // ブランチ名は gausfit.C のために "time_diff" のままにする

        // ==============================================================================
        // --- 5. 2回目のスキャン：決定した時間範囲でヒットを選択し、TTreeを作成 ---
        // ==============================================================================
        std::cout << "\n--- Pass 2: Processing hits with dynamic time cut (Trigger hgain >= " << TRIGGER_HGAIN_THRESHOLD << ") ---" << std::endl;
        std::set<int> unique_channels;

        for (long iEvent = 0; iEvent < nEvents; ++iEvent) {
            tree->GetEntry(iEvent);
            if (iEvent % 10000 == 0) {
                std::cout << "Processing event: " << iEvent << " / " << nEvents << std::endl;
            }

            if (v_trig->empty()) continue;

            // [追加箇所] 2. トリガーのhgainによる選別 (Pass 2)
            // Pass 1と同じ条件でフィルタリングしないと、ピーク位置の整合性が取れない可能性があるためここでも適用
            if (v_trig->at(0).hgain < TRIGGER_HGAIN_THRESHOLD) {
                continue; 
            }

            double trigger_tdc = v_trig->at(0).tdc;
            double trigger_time = v_trig->at(0).time;
            eventID = iEvent;

            for (const Hit& hit : *v_hit) {
                time_diff = hit.time - trigger_time; // この時点では (s) 単位

                // 1回目のスキャンで決定した動的な時間範囲でヒットを選択 (s 単位で比較)
                if (time_diff < time_cut_low || time_diff > time_cut_high) {
                    continue; 
                }
                
                ch = hit.channel;
                hgain = hit.hgain;
                lgain = hit.lgain;
                to
                t = hit.tot;
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

                    // 1. ビン幅を 1.0 (adc) に設定
                    const double bin_width = 1.0;
                    // 2. ビンの境界が整数になるように調整
                    double xlow = TMath::Floor(min_val - margin);
                    double xup = TMath::Ceil(max_val + margin);
                    // 3. ビン数を計算
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

                    // 1. ビン幅を 1.0 (adc) に設定
                    const double bin_width = 1.0;
                    // 2. ビンの境界が整数になるように調整
                    double xlow = TMath::Floor(min_val - margin);
                    double xup = TMath::Ceil(max_val + margin);
                    // 3. ビン数を計算
                    int nbins = TMath::Nint((xup - xlow) / bin_width);
                    if (nbins <= 0) nbins = 1;

                    TH1D* h_lgain = new TH1D(Form("h_lgain_ch%d", ch_num), Form("Low Gain ADC Ch %d", ch_num), nbins, xlow + bin_width/2, xup + bin_width/2);
                    new_tree->Draw(Form("lgain>>%s", h_lgain->GetName()), selection, "goff");
                }

                // [変更箇所] --- tot (ビン幅 1.0 で自動範囲設定) ---
                // (totもTDC値の一種である可能性が高いため、tdc_diffと同様にビン幅1.0に設定します)
                new_tree->Draw("tot", selection, "goff");
                n_selected = new_tree->GetSelectedRows();
                if (n_selected > 0) {
                    double min_val = TMath::MinElement(n_selected, new_tree->GetV1());
                    double max_val = TMath::MaxElement(n_selected, new_tree->GetV1());
                    double margin = (max_val - min_val) * 0.05;
                    if (margin == 0) margin = 1;

                    // 1. ビン幅を 1.0 (tdc) に設定
                    const double bin_width = 1.0;
                    // 2. ビンの境界が整数になるように調整
                    double xlow = TMath::Floor(min_val - margin);
                    double xup = TMath::Ceil(max_val + margin);
                    // 3. ビン数を計算
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
                    
                    // 1. ビン幅を 1.0 (tdc) に設定
                    const double bin_width = 1.0;
                    // 2. ビンの境界が整数になるように調整
                    double xlow = TMath::Floor(min_val - margin);
                    double xup = TMath::Ceil(max_val + margin);
                    // 3. ビン数を計算
                    int nbins = TMath::Nint((xup - xlow) / bin_width);
                    if (nbins <= 0) nbins = 1;

                    TH1D* h_tdc_diff = new TH1D(Form("h_tdc_diff_ch%d", ch_num), Form("TDC - Trigger TDC Ch %d", ch_num), nbins, xlow+bin_width/2, xup+bin_width/2);
                    new_tree->Draw(Form("tdc_diff>>%s", h_tdc_diff->GetName()), selection, "goff");
                }

                // [変更箇所] --- time_diff (ビン幅 0.25 ns で自動範囲設定) ---
                // TTreeの "time_diff" ブランチには (ns) 単位のデータが入っている
                new_tree->Draw("time_diff", selection, "goff");
                n_selected = new_tree->GetSelectedRows();
                if (n_selected > 0) {
                    double min_val = TMath::MinElement(n_selected, new_tree->GetV1());
                    double max_val = TMath::MaxElement(n_selected, new_tree->GetV1());
                    double margin = (max_val - min_val) * 0.05;
                    
                    // 2. マージンのデフォルト値を (ns) 単位に修正
                    if (margin == 0) margin = 1.0; 

                    // 1. ビン幅を 0.25 (ns) に設定
                    const double bin_width = 0.25;
                    double xlow_raw = min_val - margin;
                    double xup_raw = max_val + margin;

                    // 2. ビンの境界が 0.25 の倍数になるように調整
                    double xlow = TMath::Floor(xlow_raw / bin_width) * bin_width;
                    double xup = TMath::Ceil(xup_raw / bin_width) * bin_width;

                    // 3. ビン数を計算
                    int nbins = TMath::Nint((xup - xlow) / bin_width);
                    if (nbins <= 0) nbins = 1;

                    // 4. ヒストグラムのタイトルを (s) から (ns) に修正 (これは既に修正済み)
                    TH1D* h_time_diff = new TH1D(Form("h_time_diff_ch%d", ch_num), Form("Time - Trigger Time (ns) Ch %d", ch_num), nbins, xlow+bin_width/2, xup+bin_width/2);
                    new_tree->Draw(Form("time_diff>>%s", h_time_diff->GetName()), selection, "goff");
                }
            }
        }
        // */
        // ==============================================================================


        // --- 7. ファイルの書き込みとクローズ ---
        std::cout << "\nWriting TTree and histograms to " << output_file << std::endl;
        h_time_prescan->Write(); // ピーク特定に使ったヒストグラムも参考として保存
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

// /*id: eventtree2hist.C*/
// /*Place: hkpd@hkdaq:~/hkelec/analysis_macro/ */
// /*Last Edit: 2025-10-14 Gemini (two-pass analysis, full histogram block)*/
// /* (修正: 2025-10-27 Gemini (ns単位に変更) ) */
// /* (修正: 2025-10-30 Gemini (ヒストグラムのビン幅を最小単位に指定) ) */
// /*eventtree.root to triggered data TTree and optional histograms*/
// /*コンパイル可能*/

// #include <TFile.h>
// #include <TTree.h>
// #include <TH1D.h>
// #include <TString.h>
// #include <TMath.h>
// #include <iostream>
// #include <vector>
// #include <set>

// // 外部で定義されたHitクラスとMetaDataクラスのヘッダを読み込む
// #include "RootInterface/Hit.h"
// #include "RootInterface/MetaData.h"

// void read_event_tree(TString input_file, TString output_file) {

//     // --- 1. ファイルの準備 ---
//     auto ifile = TFile::Open(input_file, "READ");
//     if (!ifile || ifile->IsZombie()) {
//         std::cerr << "Error: Could not open input file " << input_file << std::endl;
//         return;
//     }
//     auto ofile = TFile::Open(output_file, "RECREATE");
//     if (!ofile || ofile->IsZombie()) {
//         std::cerr << "Error: Could not create output file " << output_file << std::endl;
//         ifile->Close();
//         return;
//     }

//     // --- 2. 入力TTreeの準備 ---
//     auto tree = ifile->Get<TTree>("event");
//     if (!tree) {
//         std::cerr << "Error: Could not find TTree 'event' in " << input_file << std::endl;
//         ifile->Close();
//         ofile->Close();
//         return;
//     }
//     std::vector<Hit>* v_hit = nullptr;
//     std::vector<Hit>* v_trig = nullptr;
//     tree->SetBranchAddress("NormalHits", &v_hit);
//     tree->SetBranchAddress("TriggerHits", &v_trig);

//     // メタデータを読み込んで表示する
//     auto metadata = ifile->Get<MetaData>("metadata");
//     if (metadata) {
//         metadata->Print();
//     } else {
//         std::cout << "Warning: Could not find 'metadata' object in " << input_file << std::endl;
//     }
    
//     long nEvents = tree->GetEntries();

//     // ==============================================================================
//     // --- 3. 1回目のスキャン： time_diff のピーク位置を特定する ---
//     // ==============================================================================
//     std::cout << "\n--- Pass 1: Finding time_diff peak ---" << std::endl;
    
//     // ピーク特定用の仮ヒストグラム (0 ns から 1000 ns までを 1000 分割 => ビン幅 1 ns)
//     TH1D* h_time_prescan = new TH1D("h_time_prescan", "Pre-scan for time peak; Time Difference (s); Counts", 1000, 0, 1000e-9);

//     for (long iEvent = 0; iEvent < nEvents; ++iEvent) {
//         tree->GetEntry(iEvent);
//         if (iEvent % 10000 == 0) {
//             std::cout << "Scanning event: " << iEvent << " / " << nEvents << std::endl;
//         }

//         if (v_trig->empty()) continue;
//         double trigger_time = v_trig->at(0).time;

//         for (const Hit& hit : *v_hit) {
//             double time_diff = hit.time - trigger_time;
//             // 0 ns < time_diff < 1000 ns のヒットのみでヒストグラムを作成
//             if (time_diff > 0 && time_diff < 1000e-9) {
//                 h_time_prescan->Fill(time_diff);
//             }
//         }
//     }

//     // ピークのビンを探し、その中心値を time_peak として取得
//     int peak_bin = h_time_prescan->GetMaximumBin();
//     double time_peak = h_time_prescan->GetXaxis()->GetBinCenter(peak_bin);

//     // ピーク値からカット範囲を決定 (ピーク ± 50 ns)
//     double time_cut_low  = time_peak - 50.0e-9;
//     double time_cut_high = time_peak + 50.0e-9;
    
//     std::cout << "Time peak found at: " << time_peak * 1e9 << " ns" << std::endl;
//     std::cout << "Applying new time window: " << time_cut_low * 1e9 << " ns to " << time_cut_high * 1e9 << " ns" << std::endl;

//     // --- 4. 出力用TTreeの定義 ---
//     TTree* new_tree = new TTree("processed_hits", "Processed Hit Data per Channel");
//     int eventID;
//     int ch;
//     double hgain, lgain, tot;
//     double tdc_diff, time_diff; // time_diff は ns 単位で保存されるように変更
//     new_tree->Branch("eventID", &eventID, "eventID/I");
//     new_tree->Branch("ch", &ch, "ch/I");
//     new_tree->Branch("hgain", &hgain, "hgain/D");
//     new_tree->Branch("lgain", &lgain, "lgain/D");
//     new_tree->Branch("tot", &tot, "tot/D");
//     new_tree->Branch("tdc_diff", &tdc_diff, "tdc_diff/D");
//     new_tree->Branch("time_diff", &time_diff, "time_diff/D"); // ブランチ名は gausfit.C のために "time_diff" のままにする

//     // ==============================================================================
//     // --- 5. 2回目のスキャン：決定した時間範囲でヒットを選択し、TTreeを作成 ---
//     // ==============================================================================
//     std::cout << "\n--- Pass 2: Processing hits with dynamic time cut ---" << std::endl;
//     std::set<int> unique_channels;

//     for (long iEvent = 0; iEvent < nEvents; ++iEvent) {
//         tree->GetEntry(iEvent);
//         if (iEvent % 10000 == 0) {
//             std::cout << "Processing event: " << iEvent << " / " << nEvents << std::endl;
//         }

//         if (v_trig->empty()) continue;
//         double trigger_tdc = v_trig->at(0).tdc;
//         double trigger_time = v_trig->at(0).time;
//         eventID = iEvent;

//         for (const Hit& hit : *v_hit) {
//             time_diff = hit.time - trigger_time; // この時点では (s) 単位

//             // 1回目のスキャンで決定した動的な時間範囲でヒットを選択 (s 単位で比較)
//             if (time_diff < time_cut_low || time_diff > time_cut_high) {
//                 continue; 
//             }
            
//             ch = hit.channel;
//             hgain = hit.hgain;
//             lgain = hit.lgain;
//             tot = hit.tot;
//             tdc_diff = hit.tdc - trigger_tdc;

//             // 1. TTreeに保存する直前に (s) から (ns) へ単位を変換
//             time_diff *= 1e9; 

//             new_tree->Fill();
//             unique_channels.insert(ch);
//         }
//     }

//     // ==============================================================================
//     // --- 6. (オプション) TTreeからチャンネルごとのヒストグラムを作成 ---
//     // ==============================================================================
//     // このブロックのコメントアウトを外して使用してください。
//     // /*
//     if(!unique_channels.empty()){
//         std::cout << "\nFound " << unique_channels.size() << " unique channels. Creating histograms..." << std::endl;
//         ofile->cd(); 
//         for (int ch_num : unique_channels) {
//             TString selection = Form("ch==%d", ch_num);
//             Long64_t n_selected = 0;

//             // [変更箇所] --- hgain (ビン幅 1.0 で自動範囲設定) ---
//             new_tree->Draw("hgain", selection, "goff");
//             n_selected = new_tree->GetSelectedRows();
//             if (n_selected > 0) {
//                 double min_val = TMath::MinElement(n_selected, new_tree->GetV1());
//                 double max_val = TMath::MaxElement(n_selected, new_tree->GetV1());
//                 double margin = (max_val - min_val) * 0.05;
//                 if (margin == 0) margin = 1;

//                 // 1. ビン幅を 1.0 (adc) に設定
//                 const double bin_width = 1.0;
//                 // 2. ビンの境界が整数になるように調整
//                 double xlow = TMath::Floor(min_val - margin);
//                 double xup = TMath::Ceil(max_val + margin);
//                 // 3. ビン数を計算
//                 int nbins = TMath::Nint((xup - xlow) / bin_width);
//                 if (nbins <= 0) nbins = 1;
                
//                 TH1D* h_hgain = new TH1D(Form("h_hgain_ch%d", ch_num), Form("High Gain ADC Ch %d", ch_num), nbins, xlow + bin_width/2, xup + bin_width/2);
//                 new_tree->Draw(Form("hgain>>%s", h_hgain->GetName()), selection, "goff");
//             }

//             // [変更箇所] --- lgain (ビン幅 1.0 で自動範囲設定) ---
//             new_tree->Draw("lgain", selection, "goff");
//             n_selected = new_tree->GetSelectedRows();
//             if (n_selected > 0) {
//                 double min_val = TMath::MinElement(n_selected, new_tree->GetV1());
//                 double max_val = TMath::MaxElement(n_selected, new_tree->GetV1());
//                 double margin = (max_val - min_val) * 0.05;
//                 if (margin == 0) margin = 1;

//                 // 1. ビン幅を 1.0 (adc) に設定
//                 const double bin_width = 1.0;
//                 // 2. ビンの境界が整数になるように調整
//                 double xlow = TMath::Floor(min_val - margin);
//                 double xup = TMath::Ceil(max_val + margin);
//                 // 3. ビン数を計算
//                 int nbins = TMath::Nint((xup - xlow) / bin_width);
//                 if (nbins <= 0) nbins = 1;

//                 TH1D* h_lgain = new TH1D(Form("h_lgain_ch%d", ch_num), Form("Low Gain ADC Ch %d", ch_num), nbins, xlow + bin_width/2, xup + bin_width/2);
//                 new_tree->Draw(Form("lgain>>%s", h_lgain->GetName()), selection, "goff");
//             }

//             // [変更箇所] --- tot (ビン幅 1.0 で自動範囲設定) ---
//             // (totもTDC値の一種である可能性が高いため、tdc_diffと同様にビン幅1.0に設定します)
//             new_tree->Draw("tot", selection, "goff");
//             n_selected = new_tree->GetSelectedRows();
//             if (n_selected > 0) {
//                 double min_val = TMath::MinElement(n_selected, new_tree->GetV1());
//                 double max_val = TMath::MaxElement(n_selected, new_tree->GetV1());
//                 double margin = (max_val - min_val) * 0.05;
//                 if (margin == 0) margin = 1;

//                 // 1. ビン幅を 1.0 (tdc) に設定
//                 const double bin_width = 1.0;
//                 // 2. ビンの境界が整数になるように調整
//                 double xlow = TMath::Floor(min_val - margin);
//                 double xup = TMath::Ceil(max_val + margin);
//                 // 3. ビン数を計算
//                 int nbins = TMath::Nint((xup - xlow) / bin_width);
//                 if (nbins <= 0) nbins = 1;

//                 TH1D* h_tot = new TH1D(Form("h_tot_ch%d", ch_num), Form("Time over Threshold Ch %d", ch_num), nbins, xlow + bin_width/2, xup + bin_width/2);
//                 new_tree->Draw(Form("tot>>%s", h_tot->GetName()), selection, "goff");
//             }

//             // [変更箇所] --- tdc_diff (ビン幅 1.0 で自動範囲設定) ---
//             new_tree->Draw("tdc_diff", selection, "goff");
//             n_selected = new_tree->GetSelectedRows();
//             if (n_selected > 0) {
//                 double min_val = TMath::MinElement(n_selected, new_tree->GetV1());
//                 double max_val = TMath::MaxElement(n_selected, new_tree->GetV1());
//                 double margin = (max_val - min_val) * 0.05;
//                 if (margin == 0) margin = 1;
                
//                 // 1. ビン幅を 1.0 (tdc) に設定
//                 const double bin_width = 1.0;
//                 // 2. ビンの境界が整数になるように調整
//                 double xlow = TMath::Floor(min_val - margin);
//                 double xup = TMath::Ceil(max_val + margin);
//                 // 3. ビン数を計算
//                 int nbins = TMath::Nint((xup - xlow) / bin_width);
//                 if (nbins <= 0) nbins = 1;

//                 TH1D* h_tdc_diff = new TH1D(Form("h_tdc_diff_ch%d", ch_num), Form("TDC - Trigger TDC Ch %d", ch_num), nbins, xlow+bin_width/2, xup+bin_width/2);
//                 new_tree->Draw(Form("tdc_diff>>%s", h_tdc_diff->GetName()), selection, "goff");
//             }

//             // [変更箇所] --- time_diff (ビン幅 0.25 ns で自動範囲設定) ---
//             // TTreeの "time_diff" ブランチには (ns) 単位のデータが入っている
//             new_tree->Draw("time_diff", selection, "goff");
//             n_selected = new_tree->GetSelectedRows();
//             if (n_selected > 0) {
//                 double min_val = TMath::MinElement(n_selected, new_tree->GetV1());
//                 double max_val = TMath::MaxElement(n_selected, new_tree->GetV1());
//                 double margin = (max_val - min_val) * 0.05;
                
//                 // 2. マージンのデフォルト値を (ns) 単位に修正
//                 if (margin == 0) margin = 1.0; 

//                 // 1. ビン幅を 0.25 (ns) に設定
//                 const double bin_width = 0.25;
//                 double xlow_raw = min_val - margin;
//                 double xup_raw = max_val + margin;

//                 // 2. ビンの境界が 0.25 の倍数になるように調整
//                 double xlow = TMath::Floor(xlow_raw / bin_width) * bin_width;
//                 double xup = TMath::Ceil(xup_raw / bin_width) * bin_width;

//                 // 3. ビン数を計算
//                 int nbins = TMath::Nint((xup - xlow) / bin_width);
//                 if (nbins <= 0) nbins = 1;

//                 // 4. ヒストグラムのタイトルを (s) から (ns) に修正 (これは既に修正済み)
//                 TH1D* h_time_diff = new TH1D(Form("h_time_diff_ch%d", ch_num), Form("Time - Trigger Time (ns) Ch %d", ch_num), nbins, xlow+bin_width/2, xup+bin_width/2);
//                 new_tree->Draw(Form("time_diff>>%s", h_time_diff->GetName()), selection, "goff");
//             }
//         }
//     }
//     // */
//     // ==============================================================================


//     // --- 7. ファイルの書き込みとクローズ ---
//     std::cout << "\nWriting TTree and histograms to " << output_file << std::endl;
//     h_time_prescan->Write(); // ピーク特定に使ったヒストグラムも参考として保存
//     ofile->Write(); // new_treeと、もし作成されていれば他のヒストグラムを保存
//     ofile->Close();
//     ifile->Close();
// }

// // --- スタンドアロン実行のためのmain関数 ---
// int main(int argc, char* argv[]) {
//     if (argc != 3) {
//         std::cerr << "Usage: " << argv[0] << " <input_file.root> <output_file.root>" << std::endl;
//         return 1;
//     }
//     read_event_tree(argv[1], argv[2]);
//     return 0;
// }
