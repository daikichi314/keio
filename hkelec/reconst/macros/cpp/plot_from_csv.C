/*
 * id: plot_from_csv.C
 * Place: /home/daiki/keio/hkelec/reconst/macros/cpp/
 * Last Edit: 2026-01-08 Gemini
 *
 * 概要: 指定ディレクトリ内のCSVファイルを読み込み、
 * Charge vs 各種パラメータのグラフを作成し、フィッティングを行う。
 * plot_summary.Cで作成されたCSVファイルから外れ値を手動で除去した後の処理用。
 *
 * コンパイル:
 * g++ plot_from_csv.C -o plot_from_csv $(root-config --cflags --glibs)
 */

#include <TSystem.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TAxis.h>
#include <TString.h>
#include <TStyle.h>
#include <TMath.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TMatrixDSym.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <map>
#include <algorithm>
#include <functional>
#include <cmath>

// 誤差伝播: 派生量 (例: 最小値) の不確かさを数値微分で推定
Double_t GetDerivedError(TF1 *f, const TMatrixDSym &cov, std::function<Double_t(TF1*)> func)
{
    int nPar = f->GetNpar();
    std::vector<Double_t> params(nPar);
    for(int i=0; i<nPar; ++i) params[i] = f->GetParameter(i);

    Double_t epsilon = 1e-4;
    std::vector<Double_t> grad(nPar);

    for (int i = 0; i < nPar; ++i) {
        Double_t original_val = params[i];
        f->SetParameter(i, original_val + epsilon);
        Double_t val_plus = func(f);
        f->SetParameter(i, original_val - epsilon);
        Double_t val_minus = func(f);
        grad[i] = (val_plus - val_minus) / (2.0 * epsilon);
        f->SetParameter(i, original_val);
    }

    Double_t variance = 0.0;
    for (int i = 0; i < nPar; ++i) {
        for (int j = 0; j < nPar; ++j) {
            variance += grad[i] * grad[j] * cov(i, j);
        }
    }
    return (variance > 0) ? TMath::Sqrt(variance) : 0.0;
}

// CSVファイルからデータを読み込むための構造体
struct GraphData {
    std::vector<double> x, ex, y, ey;
    std::vector<int> include_in_fit;
};

// CSV読み込み関数
bool ReadCSV(const std::string& filename, GraphData& data) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open " << filename << std::endl;
        return false;
    }

    std::string line;
    // ヘッダー行をスキップ
    std::getline(file, line);

    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string token;
        std::vector<double> vals;

        while (std::getline(ss, token, ',')) {
            try {
                vals.push_back(std::stod(token));
            } catch (...) {
                continue;
            }
        }

        if (vals.size() >= 4) {
            data.x.push_back(vals[0]);
            data.ex.push_back(vals[1]);
            data.y.push_back(vals[2]);
            data.ey.push_back(vals[3]);
            // 5番目の列がinclude_in_fit（デフォルト=1）
            int flag = (vals.size() >= 5) ? static_cast<int>(vals[4]) : 1;
            data.include_in_fit.push_back(flag);
        }
    }
    file.close();
    return !data.x.empty();
}

// ファイル名からチャンネル番号とタイプを抽出
bool ParseFilename(const std::string& filename, int& ch, std::string& type) {
    // 期待形式: Charge_vs_<type>_ch<XX>.csv
    size_t pos_vs = filename.find("Charge_vs_");
    if (pos_vs == std::string::npos) return false;

    size_t pos_ch = filename.find("_ch");
    if (pos_ch == std::string::npos) return false;

    size_t pos_csv = filename.find(".csv");
    if (pos_csv == std::string::npos) return false;

    // タイプ抽出
    type = filename.substr(pos_vs + 10, pos_ch - pos_vs - 10);

    // チャンネル番号抽出
    std::string ch_str = filename.substr(pos_ch + 3, pos_csv - pos_ch - 3);
    try {
        ch = std::stoi(ch_str);
    } catch (...) {
        return false;
    }

    return true;
}

// メイン処理関数
void ProcessCSVDirectory(const std::string& csv_dir, bool save_pdf = true) {
    // 出力ファイル（CSV）
    TString out_csv_path = TString(csv_dir) + "/fit_results_from_csv.csv";
    std::ofstream outfile(out_csv_path.Data());
    outfile << "ch,graph_type,p0,p0_err,p1,p1_err,p2,p2_err,p3,p3_err,chi2,ndf,min_val,min_err,at_charge" << std::endl;

    // ディレクトリ内のCSVファイルを取得
    void* dirp = gSystem->OpenDirectory(csv_dir.c_str());
    if (!dirp) {
        std::cerr << "Error: Cannot open directory " << csv_dir << std::endl;
        return;
    }

    std::map<int, std::map<std::string, GraphData>> data_map;
    const char* entry;
    while ((entry = gSystem->GetDirEntry(dirp))) {
        std::string filename(entry);
        if (filename.find("Charge_vs_") != std::string::npos && filename.find(".csv") != std::string::npos) {
            int ch;
            std::string type;
            if (ParseFilename(filename, ch, type)) {
                std::string fullpath = csv_dir + "/" + filename;
                GraphData gdata;
                if (ReadCSV(fullpath, gdata)) {
                    data_map[ch][type] = gdata;
                    std::cout << "Loaded: " << filename << " (ch=" << ch << ", type=" << type << ", points=" << gdata.x.size() << ")" << std::endl;
                }
            }
        }
    }
    gSystem->FreeDirectory(dirp);

    // 各チャンネル・タイプについて処理
    for (auto& [ch, type_map] : data_map) {
        for (auto& [type, gdata] : type_map) {
            if (gdata.x.empty()) continue;

            std::cout << "\nProcessing Ch" << ch << " " << type << " (" << gdata.x.size() << " points)" << std::endl;

            // データ範囲取得
            double x_min_data = *std::min_element(gdata.x.begin(), gdata.x.end());
            double x_max_data = *std::max_element(gdata.x.begin(), gdata.x.end());
            double y_min_data = *std::min_element(gdata.y.begin(), gdata.y.end());
            double y_max_data = *std::max_element(gdata.y.begin(), gdata.y.end());

            // x方向の描画範囲（全データの範囲）
            double display_min = (x_min_data < 0) ? x_min_data * 1.1 : 0.0;
            double display_max = (x_max_data > 0) ? x_max_data * 1.1 : 100.0;
            if (display_min <= 0) display_min = 1e-6;

            // y方向の描画範囲（外れ値も含めてマージンを確保）
            double y_span = y_max_data - y_min_data;
            if (y_span <= 0) y_span = std::max(std::abs(y_max_data), 1.0);
            double y_min_plot = y_min_data - 0.1 * y_span;
            double y_max_plot = y_max_data + 0.1 * y_span;

            // include_in_fit=1 のデータのみでフィッティング用TGraphErrors作成
            std::vector<double> x_fit, ex_fit, y_fit, ey_fit;
            for (size_t i = 0; i < gdata.x.size(); ++i) {
                if (gdata.include_in_fit[i] == 1) {
                    x_fit.push_back(gdata.x[i]);
                    ex_fit.push_back(gdata.ex[i]);
                    y_fit.push_back(gdata.y[i]);
                    ey_fit.push_back(gdata.ey[i]);
                }
            }

            // フィット範囲決定（フィットに使う点だけで計算）
            double x_min_fit = *std::min_element(x_fit.begin(), x_fit.end());
            double x_max_fit = *std::max_element(x_fit.begin(), x_fit.end());
            double range_min = (x_min_fit < 0) ? x_min_fit * 1.1 : x_min_fit * 0.9;
            double range_max = x_max_fit * 1.1;
            if (range_min <= 0) range_min = 1e-6;

            // フィッティング用グラフ
            TGraphErrors* gr_fit = new TGraphErrors(x_fit.size(), x_fit.data(), y_fit.data(), 
                                                     ex_fit.data(), ey_fit.data());

            // 描写用グラフ：全データポイント
            TGraphErrors* gr_all = new TGraphErrors(gdata.x.size(), gdata.x.data(), gdata.y.data(), 
                                                     gdata.ex.data(), gdata.ey.data());

            std::string y_unit = "[ns]";
            if (type == "Gamma" || type.find("Amp") != std::string::npos) y_unit = "[arb. units]";
            if (type == "Mean" || type == "Peak" || type == "GausMu") y_unit = "[ns (abs)]";
            gr_all->SetTitle(Form("Ch%d %s;Charge [pC];%s %s", ch, type.c_str(), type.c_str(), y_unit.c_str()));
            
            // フィットモデル
            TF1* f_model = new TF1("f_model", "[0]*pow(x,-0.5) + [1] + [2]*x + [3]*x*x", range_min, range_max);
            f_model->SetLineColor(kRed);

            // 初期パラメータを設定（1段階目のフィット前）
            // y_min_data に基づいて初期パラメータを切り替え
            if (y_min_data >= 100.0) {
                f_model->SetParameter(0, 18.66);   // p0の初期値
                f_model->SetParameter(1, 247.0);  // p1の初期値
                f_model->SetParameter(2, -0.005);    // p2の初期値
                f_model->SetParameter(3, 0.0);    // p3の初期値
            } else {
                f_model->SetParameter(0, 3.5);   // p0の初期値
                f_model->SetParameter(1, 0.0);   // p1の初期値
                f_model->SetParameter(2, 0.0);    // p2の初期値
                f_model->SetParameter(3, 0.0);    // p3の初期値
            }

            // フィットに使える点が足りない場合はスキップ
            if (x_fit.size() < 4) {
                std::cerr << "Warning: Not enough include_in_fit points for fitting (Ch" << ch << ", " << type << "). Skipping fit." << std::endl;
                // PDFだけ出力して次へ
                if (save_pdf) {
                    TCanvas* c = new TCanvas("c", "c", 800, 600);
                    c->SetGrid();

                    TGraphErrors* gr_included = new TGraphErrors();
                    TGraphErrors* gr_excluded = new TGraphErrors();
                    for (size_t i = 0; i < gdata.x.size(); ++i) {
                        if (gdata.include_in_fit[i] == 1) {
                            gr_included->AddPoint(gdata.x[i], gdata.y[i]);
                            gr_included->SetPointError(gr_included->GetN()-1, gdata.ex[i], gdata.ey[i]);
                        } else {
                            gr_excluded->AddPoint(gdata.x[i], gdata.y[i]);
                            gr_excluded->SetPointError(gr_excluded->GetN()-1, gdata.ex[i], gdata.ey[i]);
                        }
                    }

                    gr_included->SetTitle(Form("Ch%d %s;Charge [pC];%s %s", ch, type.c_str(), type.c_str(), y_unit.c_str()));
                    gr_included->SetMarkerStyle(20);
                    gr_included->SetMarkerColor(kBlack);
                    gr_included->SetMarkerSize(0.8);
                    gr_included->GetXaxis()->SetLimits(display_min, display_max);
                    gr_included->GetYaxis()->SetRangeUser(y_min_plot, y_max_plot);

                    gr_excluded->SetMarkerStyle(5);
                    gr_excluded->SetMarkerColor(kGray);
                    gr_excluded->SetMarkerSize(0.8);

                    gr_included->Draw("APE");
                    gr_excluded->Draw("PE");

                    TString pdf_name = Form("%s/Refitted_Charge_vs_%s_ch%02d.pdf", csv_dir.c_str(), type.c_str(), ch);
                    c->SaveAs(pdf_name);

                    delete gr_included;
                    delete gr_excluded;
                    delete c;
                }

                delete f_model;
                delete gr_fit;
                delete gr_all;
                continue;
            }

            // 2段階フィット（フィッティング用データのみ）
            TFitResultPtr r1 = gr_fit->Fit(f_model, "QS", "", range_min, range_max);
            double p0_init = f_model->GetParameter(0);
            double p1_init = f_model->GetParameter(1);
            double p2_init = f_model->GetParameter(2);
            double p3_init = f_model->GetParameter(3);
            f_model->SetParameters(p0_init, p1_init, p2_init, p3_init);
            TFitResultPtr r2 = gr_fit->Fit(f_model, "S", "", range_min, range_max);

            // 最小値計算
            double min_val = f_model->GetMinimum(range_min, range_max);
            double at_charge = f_model->GetMinimumX(range_min, range_max);

            // 最小値の誤差
            TMatrixDSym cov = r2 ? r2->GetCovarianceMatrix() : TMatrixDSym();
            auto min_value_func = [&](TF1* func) -> Double_t {
                double xm = func->GetMinimumX(range_min, range_max);
                return func->Eval(xm);
            };
            double min_err = 0.0;
            if (cov.GetNrows() == f_model->GetNpar()) {
                min_err = GetDerivedError(f_model, cov, min_value_func);
            } else {
                std::cerr << "Warning: covariance matrix unavailable for Ch" << ch << ", " << type << " (min_err set to 0)" << std::endl;
            }

            // 結果出力
            outfile << ch << "," << type;
            for(int i=0; i<4; ++i) {
                outfile << "," << f_model->GetParameter(i) << "," << f_model->GetParError(i);
            }
            outfile << "," << f_model->GetChisquare() << "," << f_model->GetNDF() 
                    << "," << min_val << "," << min_err << "," << at_charge << std::endl;

            // PDF出力
            if (save_pdf) {
                TCanvas* c = new TCanvas("c", "c", 800, 600);
                c->SetGrid();
                gr_all->GetXaxis()->SetLimits(display_min, display_max);
                
                // グラフを別々に描画
                // include_in_fit=1 のデータ（黒、マーカー●）
                TGraphErrors* gr_included = new TGraphErrors();
                // include_in_fit=0 のデータ（グレー、マーカー×）
                TGraphErrors* gr_excluded = new TGraphErrors();
                
                for (size_t i = 0; i < gdata.x.size(); ++i) {
                    if (gdata.include_in_fit[i] == 1) {
                        gr_included->AddPoint(gdata.x[i], gdata.y[i]);
                        gr_included->SetPointError(gr_included->GetN()-1, gdata.ex[i], gdata.ey[i]);
                    } else {
                        gr_excluded->AddPoint(gdata.x[i], gdata.y[i]);
                        gr_excluded->SetPointError(gr_excluded->GetN()-1, gdata.ex[i], gdata.ey[i]);
                    }
                }
                
                // included グラフ設定
                gr_included->SetTitle(Form("Ch%d %s;Charge [pC];%s %s", ch, type.c_str(), type.c_str(), y_unit.c_str()));
                gr_included->SetMarkerStyle(20);  // ●
                gr_included->SetMarkerColor(kBlack);
                gr_included->SetMarkerSize(0.8);
                gr_included->GetXaxis()->SetLimits(display_min, display_max);
                gr_included->GetYaxis()->SetRangeUser(y_min_plot, y_max_plot);
                
                // excluded グラフ設定
                gr_excluded->SetMarkerStyle(5);   // ×
                gr_excluded->SetMarkerColor(kBlue);
                gr_excluded->SetMarkerSize(2.0);
                
                // 描画
                gr_included->Draw("APE");
                gr_excluded->Draw("PE");
                f_model->Draw("same");

                TString pdf_name = Form("%s/Refitted_Charge_vs_%s_ch%02d.pdf", csv_dir.c_str(), type.c_str(), ch);
                c->SaveAs(pdf_name);
                delete c;
                delete gr_included;
                delete gr_excluded;
            }

            delete f_model;
            delete gr_fit;
            delete gr_all;
        }
    }

    outfile.close();
    std::cout << "\nProcessing completed." << std::endl;
    std::cout << " - Results: " << out_csv_path << std::endl;
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <csv_directory> [--no-pdf]" << std::endl;
        std::cerr << "  csv_directory: Directory containing Charge_vs_*.csv files" << std::endl;
        std::cerr << "  --no-pdf: Skip PDF generation (optional)" << std::endl;
        return 1;
    }

    std::string csv_dir = argv[1];
    bool save_pdf = true;

    if (argc >= 3 && std::string(argv[2]) == "--no-pdf") {
        save_pdf = false;
    }

    ProcessCSVDirectory(csv_dir, save_pdf);
    return 0;
}
