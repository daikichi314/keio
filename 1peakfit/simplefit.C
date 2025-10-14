#include <iostream>
#include <fstream>
#include <TString.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TF1.h>
#include <TStyle.h>
#include <TMath.h>
#include <TFile.h>
#include <cstdlib>
#include <TGraph.h>
#include <TSpectrum.h>
#include <TLine.h>
#include <TList.h>
#include <TPaveStats.h>
#include <TLatex.h>

// adcresult_fill関数のための構造体定義
struct AdcResult {
  Int_t ch;
  Double_t mean;
  Double_t variance;
  Double_t peakx;
  Double_t peaky;
  Double_t FWHMlow;
  Double_t FWHMhigh;
  Double_t FWHM;
  Double_t sigmalow;
  Double_t sigmahigh;
  Double_t sigma;
  Double_t peak1pex;
  Double_t peak1pesigma;
  Double_t peak1pexerr;
  Double_t peak1pesigmaerr;
  Double_t valley;
  Double_t valleyx;
  Double_t peak;
  Double_t pv;
  Double_t areasigmalow;
  Double_t areasigmahigh;
  Double_t areasigma;
  Double_t gainpeakx;
  Double_t gainpeak1pex;
  Double_t gainpeak1pexerr;
  Double_t chi2;
  Double_t ndf;
};

// FitPedPe関数からアクセスするグローバル変数を宣言
Double_t gainpeak, respeak, respv, resfwhm, ressigma, bsratio, chi2ndf, peakval;

// adcresult_fill関数の宣言（元のコードに定義がないため仮に宣言）
void adcresult_fill(TFile*, AdcResult&);

// FitPedPe関数のプロトタイプ宣言
void FitPedPe(TH1* h, const char* name = "", Int_t nPeaks = 1, Int_t rebin = 4, Bool_t DrawEff = kFALSE);

#ifndef rcomp
void simplefit(const char* rootfilename = "input.root", const char* outfilename = "test.root", bool Update = kFALSE, Int_t npeak = 1) {
#else
#define inputpath "/home/daiki/lab/data/20250809"
#define outputpdfpath "/home/daiki/lab/data/20250809"
#define outputtxtpath "/home/daiki/lab/data/20250809"
int main(int argc, char* argv[]){
    if (argc < 3) {
      std::cout << "USAGE : ./simplefit (Input Root filename) (Output Root filename) (optional: num of peak finding) (optional, 1:Add in existing root, 0:Recreate new (default) )" << std::endl;
      std::cout << "例 : ./simplefit input.root output.root 3 0" << std::endl;
      std::cout << "入出力先を要確認" << std::endl;
      return 0;
    }
    const char* rootfilename = argv[1];
    const char* outfilename = argv[2];
    int npeak = 2;
    if (argc > 3) {
      npeak = atoi(argv[3]);
    }
    bool Update = kFALSE;
    if (argc > 4) {
        Update = atoi(argv[4]);
    }
#endif
    
    // AdcResult構造体を宣言
    AdcResult result;

    bool OutputTxt = kTRUE;
    TString rootfilename_base = TString(rootfilename).ReplaceAll(".root", "");
    TString pdffilename = Form("%s/%s.pdf", outputpdfpath, rootfilename_base.Data());
    TString textfilename = Form("%s/%s.txt", outputtxtpath, rootfilename_base.Data());
    //出力ファイル名を調整
    gErrorIgnoreLevel = 5000;
    
    TCanvas *c = new TCanvas("c","c",800,600);
    
    TString full_input_path = Form("%s/%s", inputpath, rootfilename);
    TFile* file = new TFile(full_input_path.Data(), "READ");
    // 入力元をカレントディレクトリから変更

    TH1 *h = (TH1*)0;
    
    c->Print(Form("%s[",pdffilename.Data()));

    // TFile *rootf;
    // if (!Update) {
    // TString full_output_root_path = Form("%s/%s", outputpdfpath, outfilename);
    // rootf = TFile::Open(full_output_root_path.Data(), "RECREATE", outfilename);
    // rootf->Close();
    // }
    // // 出力先をカレントディレクトリから変更
    // ルートファイルを出力しない場合はここをコメントアウト(1)

    std::ofstream resultfile;
    if (OutputTxt) 
      resultfile.open(textfilename.Data());
    
    for (int ich = 0; ich < 4; ich++) {
      h = (TH1*)(file->Get(Form("all_ch%02d_hgain",ich))); // high gain  何のグラフを読むかでここを変える
      
      TF1* fit_func = nullptr; //0829追加
      int fit_success_flag = 0; // 0: 成功, 1: 失敗

      // ヒストグラムが見つからない場合、またはエントリ数が少ない場合の処理
      if (h == nullptr || h->GetEntries() <= 100) {
          fit_success_flag = 1; // 失敗とマーク
          
          if (h == nullptr) {
              std::cout << Form("Histogram all_ch%02d_hgain not found. Skipping.", ich) << std::endl;
          } else {
              std::cout << "Channel " << ich << " has too few entries (" << h->GetEntries() << "). Skipping fit." << std::endl;
          }

          // フィット失敗時は結果変数を無効な値に設定
          gainpeak = -1;
          respeak = -1;
          respv = -1;
          resfwhm = -1;
          ressigma = -1;
          bsratio = -1;
          chi2ndf = -1;
      } else {
          // ヒストグラムが存在し、エントリ数が十分な場合はフィットを実行
          std::cout << Form("all_ch%02d_hgain",ich) << std::endl;
          std::cout << h << std::endl;
          std::cout << "ENTRIES " << h->GetEntries() << " " << h->GetTitle() << std::endl;

          h->SetTitle(Form("%s", rootfilename));
          //       h->GetXaxis()->SetTitle("Charge [pC]");
          h->GetXaxis()->SetTitle("High Gain [ADC]");
          h->GetXaxis()->SetRangeUser(400,800);
          //h->GetXaxis()->SetRangeUser(-0.5,15);

          //c->SetLogy(1);
          
          c->cd();
          /////////////////////////////////////////////////////
          FitPedPe(h, Form("BLPMT"), 1/*npeak変更*/, 0, kFALSE); // 1ピークフィットに変更
          /////////////////////////////////////////////////////
          h->Draw();
          
          // フィット関数を描画
          TF1 *fit_func = (TF1*)h->GetFunction("gaus");
          if (fit_func) {
              fit_func->Draw("same");
          }
      }
      
      // コンソールとテキストファイルに結果を出力
      std::cout << Form("BBB RESULT : %02d, %d, ", ich, fit_success_flag) << gainpeak << ", " << respeak << ", " << respv << ", " << resfwhm << ", " << ressigma << ", " << bsratio << ", " << chi2ndf << ", " << std::endl;
      //(respeak!=0?TMath::Sqrt(resfwhm*resfwhm*respeak*respeak*1e-4-adcpedestalsigma[ich]*50.*1.e-3*0.1*adcpedestalsigma[ich]*50.*1.e-3*0.1)/respeak*100.:-1) << std::endl;
      if (OutputTxt)
          resultfile << Form("BBB RESULT : %02d, %d, ", ich, fit_success_flag) << gainpeak << ", " << respeak << ", " << respv << ", " << resfwhm << ", " << ressigma << ", " << bsratio << ", " << chi2ndf << ", " << std::endl;
      //(respeak!=0?TMath::Sqrt(resfwhm*resfwhm*respeak*respeak*1e-4-adcpedestalsigma[ich]*50.*1.e-3*0.1*adcpedestalsigma[ich]*50.*1.e-3*0.1*2.355*2.355)/respeak*100.:-1) << std::endl;
      
      result.ch = ich;

////////////////////////////////////////////////
      // rootf = TFile::Open(outfilename, "UPDATE", outfilename);
      // adcresult_fill(rootf, result);
      // rootf->Close();
      // ルートファイルを出力しない場合はここをコメントアウト(2)
////////////////////////////////////////////////

      gStyle->SetOptFit(1111);
      gStyle->SetOptTitle(1);
      
      //c->SaveAs(Form("%srun%06d_1pechargelog_p%d_%s.pdf",path,run,npeak,serial.Data()));
      //c->Print(pdffilename);
      
      // c->SetLogy(0);
      // if (peakval < 10) peakval = 10;
      // else if (peakval > 1e4) peakval = 1e4;   //0829変更
      if (h != nullptr) {
        double ymax = h->GetMaximum();
        h->GetYaxis()->SetRangeUser(0, ymax*1.2);
        h->Draw("e1");
        if (fit_func) fit_func->Draw("same");
      }

      
      // hが存在する場合にのみGetYaxis()を呼び出す
      if (h != nullptr) {
          h->GetYaxis()->SetRangeUser(0, peakval*1.6);
      }
      
      //c->SaveAs(Form("%srun%06d_1pecharge_p%d_%s.root",path,run,npeak,serial.Data()));
      //c->SaveAs(Form("%srun%06d_1pecharge_p%d_%s.pdf",path,run,npeak,serial.Data()));
      
      c->Print(pdffilename);
    }
    
    c->Print(Form("%s]",pdffilename.Data()));
    
    file->Close();
    if (OutputTxt) resultfile.close();

    return 0;
    
}

// FitPedPe関数
void FitPedPe(TH1* h, const char* name, Int_t nPeaks, Int_t rebin, Bool_t DrawEff) {

    const Int_t npx = 1000;
    //Int_t rebin = 4; // After pedestal fit
    //Int_t nPeaks = 3; // 変更: 1に固定
    //Bool_t DrawEff = kFALSE;

    gStyle->SetOptFit(1);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    //////////////////////////////
    ////   Define Functions   //////
    //////////////////////////////
    TF1 *funcped      = new TF1(Form("funcped%s",name),"gaus", h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());
    // 1ピークフィットのため、1pe以上の関数はコメントアウト
    // TF1 *func1pe      = new TF1(Form("func1pe%s",name),"gaus", h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());
    // TF1 *func2pe      = new TF1(Form("func2pe%s",name),"gaus", h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());
    funcped   ->SetLineColor(kMagenta);
    // func1pe   ->SetLineColor(kGreen+1);
    // func2pe   ->SetLineColor(kGreen);
    funcped ->SetLineStyle(3);
    // func1pe ->SetLineStyle(2);
    // func2pe ->SetLineStyle(2);
    // func1pe ->SetLineWidth(1);
    // func2pe ->SetLineWidth(1);
    funcped ->SetNpx(npx);
    // func1pe ->SetNpx(npx);
    // func2pe ->SetNpx(npx);

    // 2ピーク以上の複雑な関数をコメントアウト
    // TF1 *func2peak          = new TF1(Form("func2peak%s",name),"gaus(0)+gaus(3)+(0.5*(([6]*[3]*((TMath::Erf(((x-[1])/[2])))+(1-TMath::Erf((x-[4])/[5]))-1))))");
    // TF1 *func2peak_1peall = new TF1(Form("func1peall%s",name),"gaus(3)+(0.5*(([6]*[3]*((TMath::Erf(((x-[1])/[2])))+(1-TMath::Erf((x-[4])/[5]))-1))))");
    // TF1 *func2peak_bs     = new TF1(Form("func1bs%s",name),"(0.5*(([6]*[3]*((TMath::Erf(((x-[1])/[2])))+(1-TMath::Erf((x-[4])/[5]))-1))))");
    // func2peak       ->SetNpx(npx);
    // func2peak_1peall->SetNpx(npx);
    // func2peak_bs    ->SetNpx(npx);

    // TF1 *func3peak       = new TF1(Form("func3peak%s",name),"gaus(0)+gaus(3)+gaus(6)+(0.5*(([9]*[3]*((TMath::Erf(((x-[1])/[2])))+(1-TMath::Erf((x-[4])/[5]))-1))"
    //                                         "+((([4]-[1])/([7]-[1]))*[6]*[9]*[9]*((TMath::Erf((x-[1])/[2]))+((1-TMath::Erf((x-[7])/[8]))-1)))"
    //                                         "+((([4]-[1])/([7]-[4]))*[6]*2*(1-[9])*[9]*((TMath::Erf((x-[4])/[5]))+(1-TMath::Erf((x-[7])/[8]))-1)))");
    // //TF1 *func3peak_1peall = func2peak_1peall;//new TF1(Form("func1peall2%s",name),"gaus(3)+(0.5*(([6]*[3]*((TMath::Erf(((x-[1])/[2])))+(1-TMath::Erf((x-[4])/[5]))-1))))");//
    // TF1 *func3peak_1peall = new TF1(Form("func3peak_1peall%s",name),"gaus(3)+(0.5*(([9]*[3]*((TMath::Erf(((x-[1])/[2])))+(1-TMath::Erf((x-[4])/[5]))-1))))");
    // //TF1 *func3peak_bs     = func2peak_bs;     //new TF1(Form("func1bs2%s",name),"(0.5*(([6]*[3]*((TMath::Erf(((x-[1])/[2])))+(1-TMath::Erf((x-[4])/[5]))-1))))");          //
    // TF1 *func3peak_bs     = new TF1(Form("func3peak_1bs%s",name),"(0.5*(([9]*[3]*((TMath::Erf(((x-[1])/[2])))+(1-TMath::Erf((x-[4])/[5]))-1))))");
    // TF1 *func3peak_2bs    = new TF1(Form("func2bs%s",name),"0.5*(([9]*[3]*((TMath::Erf(((x-[1])/[2])))+(1-TMath::Erf((x-[4])/[5]))-1))"
    //                                         "+((([4]-[1])/([7]-[1]))*[6]*[9]*[9]*((TMath::Erf((x-[1])/[2]))+((1-TMath::Erf((x-[7])/[8]))-1)))"
    //                                         "+((([4]-[1])/([7]-[4]))*[6]*2*(1-[9])*[9]*((TMath::Erf((x-[4])/[5]))+(1-TMath::Erf((x-[7])/[8]))-1)))");
    // func3peak->SetNpx(npx);
    // func3peak->SetLineColor(kBlue);
    // func3peak_1peall->SetNpx(npx);
    // func3peak_2bs->SetNpx(npx);
    // func3peak_2bs->SetLineColor(kCyan+1);
    // func3peak_bs->SetNpx(npx);
    // func3peak_bs->SetLineColor(kCyan+1);
    // func3peak_bs->SetLineStyle(2);

    // 1ピークフィット用の関数を定義 (funcpedを流用)
    TF1 *func_1peak = funcped;
    func_1peak->SetLineColor(kRed); // 色を変更
    func_1peak->SetLineStyle(1);

    TF1 *func, *func1peall, *funcbs;

    // switch文を1ピーク用に修正
    switch(nPeaks) {
        case 1:
             func       = func_1peak;
             func1peall = func_1peak; // 1ピークのため、全体も同じ関数
             funcbs     = nullptr;   // バックグラウンドなし
             break;
        case 3:
        case 2:
        default:
             // 既存の2ピーク・3ピーク処理はコメントアウト
             // func       = func2peak;
             // func1peall = func2peak_1peall;
             // funcbs     = func2peak_bs;
             return; // 1ピーク以外は処理を中断
    }

    // 1ピーク用のスタイル設定
    func->SetLineColor(kBlue);
    func1peall->SetLineColor(kRed);
    // funcbsはnullptrなので設定しない
    func->SetLineStyle(1);
    func1peall->SetLineStyle(2);
    // funcbs->SetLineStyle(2); // nullptrなので設定しない

    ///////////////////////
    //// PEAK FINDING /////
    ///////////////////////
    
    /* should modify these initial values... */
    Double_t peakx_pre[3] = {440,470,500};//{1500, 1600, 1700};//{0, 2.4, 4.8};//[pC]
    Double_t peaky_pre[3] = {h->GetBinContent(h->FindBin(peakx_pre[0])), h->GetBinContent(h->FindBin(peakx_pre[1])), h->GetBinContent(h->FindBin(peakx_pre[1]))} ;

#ifdef PEAKFINDING
    TSpectrum s(nPeaks);
    s.Search(h,4,"goff",0.001);
    Int_t FoundPeaks = s.GetNPeaks();

    for(int iPeak=0;iPeak<TMath::Min(nPeaks,FoundPeaks);iPeak++){
     peakx_pre[iPeak]=s.GetPositionX()[iPeak];
     peaky_pre[iPeak]=s.GetPositionY()[iPeak];
     std::cout << "Peak Position " << iPeak << " is X: "<<peakx_pre[iPeak] <<" Y:"<< peaky_pre[iPeak] << std::endl;
    }
    if (nPeaks == 1) { // 1ピークの場合の処理を追加
      if (FoundPeaks > 0) {
        peakx_pre[0] = s.GetPositionX()[0];
        peaky_pre[0] = s.GetPositionY()[0];
      }
    } else 
    if (peakx_pre[0] > peakx_pre[1]) {
      peakx_pre[2] = peakx_pre[0];
      peakx_pre[0] = peakx_pre[1];
      peakx_pre[1] = peakx_pre[2];
      peaky_pre[2] = peaky_pre[0];
      peaky_pre[0] = peaky_pre[1];
      peaky_pre[1] = peaky_pre[2];
    }
#endif

    // 1ピークフィットのためのフィット範囲を設定
    h->Fit(func_1peak,"NQ0","",peakx_pre[0]*0.8,peakx_pre[0]*1.2);
    
    // フィット後のパラメータをセット
    Double_t fit_scale = func_1peak->GetParameter(0);
    Double_t fit_peak  = func_1peak->GetParameter(1);
    Double_t fit_sigma = func_1peak->GetParameter(2);
    
    // 不要なRebin処理をコメントアウト
    // if (rebin>1) {
    //  h->Rebin(rebin);
    //  funcped->SetParameter(0, funcped->GetParameter(0)*(Double_t)(rebin));
    //  func1pe->SetParameter(0, func1pe->GetParameter(0)*(Double_t)(rebin));
    //  funcped->SetParError(0, funcped->GetParError(0)*(Double_t)(rebin));
    //  func1pe->SetParError(0, func1pe->GetParError(0)*(Double_t)(rebin));
    // }

    // パラメータ名と制限の設定を1ピーク用に修正
    // ... 修正が複雑なため、ここは元のコードをコメントアウトし、簡潔な処理に置き換え ...
    // func_1peak->SetParName(0, "Scale");
    // func_1peak->SetParName(1, "Peak");
    // func_1peak->SetParName(2, "Sigma");
    // func_1peak->SetParLimits(2, 0, 100); // Sigmaを正の値に制限
    
    // h->Fit(func_1peak, "B", "", h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());

    // 1ピーク用のResultパラメータ計算
    // 不要な変数をコメントアウト
    // Double_t areaped      = 0;
    // Double_t area1pe      = 0;
    // Double_t areabs       = 0;
    // Double_t area1peall   = 0;

    // 1ピーク用の物理量計算
    Double_t peakx = func_1peak->GetMaximumX();
    Double_t peaky = func_1peak->GetMaximum();
    Double_t FWHMlow = peakx - func_1peak->GetX(peaky*0.5, h->GetXaxis()->GetXmin(), peakx);
    Double_t FWHMhigh = func_1peak->GetX(peaky*0.5, peakx, h->GetXaxis()->GetXmax()) - peakx;
    Double_t FWHM = FWHMhigh + FWHMlow;
    Double_t sigma = TMath::Abs(func_1peak->GetParameter(2));
    Double_t gainpeakx = peakx/1.60217662e-19*1.e-12;
    Double_t chi2 = func_1peak->GetChisquare();
    Double_t ndf = (Double_t)func_1peak->GetNDF();

    // 1ピークフィットに不要な計算をコメントアウト
    // Double_t pv = -1;
    // Double_t areabs = 0;
    // Double_t areasigma = 0;
    
    // Result parameters
    // Double_t mean = func_1peak->Mean(h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());
    // Double_t variance = func_1peak->Variance(h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());
    
    // 不要な変数にデフォルト値を設定
    Double_t mean = -1;
    Double_t variance = -1;
    Double_t pv = -1;
    Double_t areaped = -1;
    Double_t area1pe = -1;
    Double_t area1peall = -1;
    Double_t areabs = -1;
    Double_t areasigmalow = -1;
    Double_t areasigmahigh = -1;
    Double_t areasigma = -1;
    Double_t peak1pex = -1;
    Double_t peak1pesigma = -1;
    Double_t peak1pexerr = -1;
    Double_t peak1pesigmaerr = -1;
    Double_t gainpeak1pex = -1;
    
    peakval = peakx; // peakvalを1ピークのピーク位置に設定

    std::cout << std::endl << std::endl;
    std::cout << "=========== RESULT ==========" << std::endl;
    std::cout << " peakx            = " << " " <<  peakx         << " "      << std::endl;
    std::cout << " FWHM             = " << " " <<  FWHM          << " " << " (" <<  (peakx!=0?FWHM/peakx*100.:-1)  << " [%])"  << std::endl;
    std::cout << " sigma            = " << " " <<  sigma         << " " << " (" <<  (peakx!=0?sigma/peakx*100.:-1) << " [%])" << std::endl;
    std::cout << " Gain             = " << " " <<  gainpeakx     << " "      << std::endl;
    std::cout << " Chi2/NDF         = " << " " <<  chi2/ndf      << " "      << std::endl;
    std::cout << "=============================" << std::endl;

    gainpeak = gainpeakx;
    respeak = peakx;
    respv = pv;
    resfwhm = (peakx!=0?FWHM/peakx*100.:-1);
    ressigma = (peakx!=0?sigma/peakx*100.:-1);
    bsratio = 0; // 1ピークフィットのため0に固定
    chi2ndf = ndf!=0?chi2/ndf:0;

    // 1ピークフィットのため、不要な描画をコメントアウト
    // if (!DrawEff) {
    //  gPad->SetLogy(1);
    //  h->Draw("e1");
    //  funcped ->Draw("same");
    //  func1pe ->Draw("same");
    //  if (nPeaks>2) { func2pe->Draw("same"); }
    //  funcbs  ->Draw("same");
    //  if (nPeaks>2) { func3peak_bs->Draw("same"); }
    //  func1peall->Draw("same");
    //  func->Draw("same");
    //  gPad->Update();
    // } else if (DrawEff) {
    //   // ... 効率グラフの描画コード ...
    // }
    
    // // 1ピークの関数だけを描画
    // gPad->SetLogy(1);
    // h->Draw("e1");
    // func->Draw("same");//0829変更
}

// #include <iostream>
// #include <fstream>
// #include <TString.h>
// #include <TCanvas.h>
// #include <TH1.h>
// #include <TF1.h>
// #include <TStyle.h>
// #include <TMath.h>
// #include <TFile.h>
// #include <cstdlib>
// #include <TGraph.h>
// #include <TSpectrum.h>
// #include <TLine.h>
// #include <TList.h>
// #include <TPaveStats.h>
// #include <TLatex.h>

// // adcresult_fill関数のための構造体定義
// struct AdcResult {
//   Int_t ch;
//   Double_t mean;
//   Double_t variance;
//   Double_t peakx;
//   Double_t peaky;
//   Double_t FWHMlow;
//   Double_t FWHMhigh;
//   Double_t FWHM;
//   Double_t sigmalow;
//   Double_t sigmahigh;
//   Double_t sigma;
//   Double_t peak1pex;
//   Double_t peak1pesigma;
//   Double_t peak1pexerr;
//   Double_t peak1pesigmaerr;
//   Double_t valley;
//   Double_t valleyx;
//   Double_t peak;
//   Double_t pv;
//   Double_t areasigmalow;
//   Double_t areasigmahigh;
//   Double_t areasigma;
//   Double_t gainpeakx;
//   Double_t gainpeak1pex;
//   Double_t gainpeak1pexerr;
//   Double_t chi2;
//   Double_t ndf;
// };

// // FitPedPe関数からアクセスするグローバル変数を宣言
// Double_t gainpeak, respeak, respv, resfwhm, ressigma, bsratio, chi2ndf, peakval;

// // adcresult_fill関数の宣言（元のコードに定義がないため仮に宣言）
// void adcresult_fill(TFile*, AdcResult&);

// // FitPedPe関数のプロトタイプ宣言
// void FitPedPe(TH1* h, const char* name = "", Int_t nPeaks = 1, Int_t rebin = 4, Bool_t DrawEff = kFALSE);

// #ifndef rcomp
// void simplefit(const char* rootfilename = "input.root", const char* outfilename = "test.root", bool Update = kFALSE, Int_t npeak = 1) {
// #else
// #define inputpath "/home/daiki/lab/data/20250809"
// #define outputpdfpath "/home/daiki/lab/data/20250809"
// #define outputtxtpath "/home/daiki/lab/data/20250809"
// int main(int argc, char* argv[]){
//     if (argc < 3) {
//       std::cout << "USAGE : ./simplefit (Input Root filename) (Output Root filename) (optional: num of peak finding) (optional, 1:Add in existing root, 0:Recreate new (default) )" << std::endl;
//       std::cout << "例 : ./simplefit input.root output.root 3 0" << std::endl;
//       std::cout << "入出力先を要確認" << std::endl;
//       return 0;
//     }
//     const char* rootfilename = argv[1];
//     const char* outfilename = argv[2];
//     int npeak = 2;
//     if (argc > 3) {
//       npeak = atoi(argv[3]);
//     }
//     bool Update = kFALSE;
//     if (argc > 4) {
//         Update = atoi(argv[4]);
//     }
// #endif
    
//     // AdcResult構造体を宣言
//     AdcResult result;

//     bool OutputTxt = kTRUE;
//     TString rootfilename_base = TString(rootfilename).ReplaceAll(".root", "");
//     TString pdffilename = Form("%s/%s.pdf", outputpdfpath, rootfilename_base.Data());
//     TString textfilename = Form("%s/%s.txt", outputtxtpath, rootfilename_base.Data());
//     //出力ファイル名を調整
//     gErrorIgnoreLevel = 5000;
    
//     TCanvas *c = new TCanvas("c","c",800,600);
    
//     TString full_input_path = Form("%s/%s", inputpath, rootfilename);
//     TFile* file = new TFile(full_input_path.Data(), "READ");
//     // 入力元をカレントディレクトリから変更

//     TH1 *h = (TH1*)0;
    
//     c->Print(Form("%s[",pdffilename.Data()));

//     // TFile *rootf;
//     // if (!Update) {
//     // TString full_output_root_path = Form("%s/%s", outputpdfpath, outfilename);
//     // rootf = TFile::Open(full_output_root_path.Data(), "RECREATE", outfilename);
//     // rootf->Close();
//     // }
//     // // 出力先をカレントディレクトリから変更
//     // ルートファイルを出力しない場合はここをコメントアウト(1)

//     std::ofstream resultfile;
//     if (OutputTxt) 
//       resultfile.open(textfilename.Data());
    
//     for (int ich = 0; ich < 4; ich++) {
//       h = (TH1*)(file->Get(Form("all_ch%02d_hgain",ich))); // high gain  何のグラフを読むかでここを変える
      
//       int fit_success_flag = 0; // 0: 成功, 1: 失敗

//       // ヒストグラムが見つからない場合、またはエントリ数が少ない場合の処理
//       if (h == nullptr || h->GetEntries() <= 100) {
//           fit_success_flag = 1; // 失敗とマーク
          
//           if (h == nullptr) {
//               std::cout << Form("Histogram all_ch%02d_hgain not found. Skipping.", ich) << std::endl;
//           } else {
//               std::cout << "Channel " << ich << " has too few entries (" << h->GetEntries() << "). Skipping fit." << std::endl;
//           }

//           // フィット失敗時は結果変数を無効な値に設定
//           gainpeak = -1;
//           respeak = -1;
//           respv = -1;
//           resfwhm = -1;
//           ressigma = -1;
//           bsratio = -1;
//           chi2ndf = -1;
//       } else {
//           // ヒストグラムが存在し、エントリ数が十分な場合はフィットを実行
//           std::cout << Form("all_ch%02d_hgain",ich) << std::endl;
//           std::cout << h << std::endl;
//           std::cout << "ENTRIES " << h->GetEntries() << " " << h->GetTitle() << std::endl;

//           h->SetTitle(Form("%s", rootfilename));
//           //       h->GetXaxis()->SetTitle("Charge [pC]");
//           h->GetXaxis()->SetTitle("High Gain [ADC]");
//           h->GetXaxis()->SetRangeUser(400,800);
//           //h->GetXaxis()->SetRangeUser(-0.5,15);

//           //c->SetLogy(1);
          
//           c->cd();
//           /////////////////////////////////////////////////////
//           FitPedPe(h, Form("BLPMT"), 1/*npeak変更*/, 0, kFALSE); // 1ピークフィットに変更
//           /////////////////////////////////////////////////////
//           h->Draw();
          
//           // フィット関数を描画
//           TF1 *fit_func = (TF1*)h->GetFunction("gaus");
//           if (fit_func) {
//               fit_func->Draw("same");
//           }
//       }
      
//       // コンソールとテキストファイルに結果を出力
//       std::cout << Form("BBB RESULT : %02d, %d, ", ich, fit_success_flag) << gainpeak << ", " << respeak << ", " << respv << ", " << resfwhm << ", " << ressigma << ", " << bsratio << ", " << chi2ndf << ", " << std::endl;
//       //(respeak!=0?TMath::Sqrt(resfwhm*resfwhm*respeak*respeak*1e-4-adcpedestalsigma[ich]*50.*1.e-3*0.1*adcpedestalsigma[ich]*50.*1.e-3*0.1)/respeak*100.:-1) << std::endl;
//       if (OutputTxt)
//           resultfile << Form("BBB RESULT : %02d, %d, ", ich, fit_success_flag) << gainpeak << ", " << respeak << ", " << respv << ", " << resfwhm << ", " << ressigma << ", " << bsratio << ", " << chi2ndf << ", " << std::endl;
//       //(respeak!=0?TMath::Sqrt(resfwhm*resfwhm*respeak*respeak*1e-4-adcpedestalsigma[ich]*50.*1.e-3*0.1*adcpedestalsigma[ich]*50.*1.e-3*0.1*2.355*2.355)/respeak*100.:-1) << std::endl;
      
//       result.ch = ich;

// ////////////////////////////////////////////////
//       // rootf = TFile::Open(outfilename, "UPDATE", outfilename);
//       // adcresult_fill(rootf, result);
//       // rootf->Close();
//       // ルートファイルを出力しない場合はここをコメントアウト(2)
// ////////////////////////////////////////////////

//       gStyle->SetOptFit(1111);
//       gStyle->SetOptTitle(1);
      
//       //c->SaveAs(Form("%srun%06d_1pechargelog_p%d_%s.pdf",path,run,npeak,serial.Data()));
//       //c->Print(pdffilename);
      
//       c->SetLogy(0);
//       if (peakval < 10) peakval = 10;
//       else if (peakval > 1e4) peakval = 1e4;
      
//       // hが存在する場合にのみGetYaxis()を呼び出す
//       if (h != nullptr) {
//           h->GetYaxis()->SetRangeUser(0, peakval*1.6);
//       }
      
//       //c->SaveAs(Form("%srun%06d_1pecharge_p%d_%s.root",path,run,npeak,serial.Data()));
//       //c->SaveAs(Form("%srun%06d_1pecharge_p%d_%s.pdf",path,run,npeak,serial.Data()));
      
//       c->Print(pdffilename);
//     }
    
//     c->Print(Form("%s]",pdffilename.Data()));
    
//     file->Close();
//     if (OutputTxt) resultfile.close();

//     return 0;
    
// }

// // FitPedPe関数
// void FitPedPe(TH1* h, const char* name, Int_t nPeaks, Int_t rebin, Bool_t DrawEff) {

//     const Int_t npx = 1000;
//     //Int_t rebin = 4; // After pedestal fit
//     //Int_t nPeaks = 3; // 変更: 1に固定
//     //Bool_t DrawEff = kFALSE;

//     gStyle->SetOptFit(1);
//     gStyle->SetOptStat(0);
//     gStyle->SetOptTitle(0);

//     //////////////////////////////
//     ////   Define Functions   //////
//     //////////////////////////////
//     TF1 *funcped      = new TF1(Form("funcped%s",name),"gaus", h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());
//     // 1ピークフィットのため、1pe以上の関数はコメントアウト
//     // TF1 *func1pe      = new TF1(Form("func1pe%s",name),"gaus", h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());
//     // TF1 *func2pe      = new TF1(Form("func2pe%s",name),"gaus", h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());
//     funcped   ->SetLineColor(kMagenta);
//     // func1pe   ->SetLineColor(kGreen+1);
//     // func2pe   ->SetLineColor(kGreen);
//     funcped ->SetLineStyle(3);
//     // func1pe ->SetLineStyle(2);
//     // func2pe ->SetLineStyle(2);
//     // func1pe ->SetLineWidth(1);
//     // func2pe ->SetLineWidth(1);
//     funcped ->SetNpx(npx);
//     // func1pe ->SetNpx(npx);
//     // func2pe ->SetNpx(npx);

//     // 2ピーク以上の複雑な関数をコメントアウト
//     // TF1 *func2peak          = new TF1(Form("func2peak%s",name),"gaus(0)+gaus(3)+(0.5*(([6]*[3]*((TMath::Erf(((x-[1])/[2])))+(1-TMath::Erf((x-[4])/[5]))-1))))");
//     // TF1 *func2peak_1peall = new TF1(Form("func1peall%s",name),"gaus(3)+(0.5*(([6]*[3]*((TMath::Erf(((x-[1])/[2])))+(1-TMath::Erf((x-[4])/[5]))-1))))");
//     // TF1 *func2peak_bs     = new TF1(Form("func1bs%s",name),"(0.5*(([6]*[3]*((TMath::Erf(((x-[1])/[2])))+(1-TMath::Erf((x-[4])/[5]))-1))))");
//     // func2peak       ->SetNpx(npx);
//     // func2peak_1peall->SetNpx(npx);
//     // func2peak_bs    ->SetNpx(npx);

//     // TF1 *func3peak       = new TF1(Form("func3peak%s",name),"gaus(0)+gaus(3)+gaus(6)+(0.5*(([9]*[3]*((TMath::Erf(((x-[1])/[2])))+(1-TMath::Erf((x-[4])/[5]))-1))"
//     //                                         "+((([4]-[1])/([7]-[1]))*[6]*[9]*[9]*((TMath::Erf((x-[1])/[2]))+((1-TMath::Erf((x-[7])/[8]))-1)))"
//     //                                         "+((([4]-[1])/([7]-[4]))*[6]*2*(1-[9])*[9]*((TMath::Erf((x-[4])/[5]))+(1-TMath::Erf((x-[7])/[8]))-1)))");
//     // //TF1 *func3peak_1peall = func2peak_1peall;//new TF1(Form("func1peall2%s",name),"gaus(3)+(0.5*(([6]*[3]*((TMath::Erf(((x-[1])/[2])))+(1-TMath::Erf((x-[4])/[5]))-1))))");//
//     // TF1 *func3peak_1peall = new TF1(Form("func3peak_1peall%s",name),"gaus(3)+(0.5*(([9]*[3]*((TMath::Erf(((x-[1])/[2])))+(1-TMath::Erf((x-[4])/[5]))-1))))");
//     // //TF1 *func3peak_bs     = func2peak_bs;     //new TF1(Form("func1bs2%s",name),"(0.5*(([6]*[3]*((TMath::Erf(((x-[1])/[2])))+(1-TMath::Erf((x-[4])/[5]))-1))))");          //
//     // TF1 *func3peak_bs     = new TF1(Form("func3peak_1bs%s",name),"(0.5*(([9]*[3]*((TMath::Erf(((x-[1])/[2])))+(1-TMath::Erf((x-[4])/[5]))-1))))");
//     // TF1 *func3peak_2bs    = new TF1(Form("func2bs%s",name),"0.5*(([9]*[3]*((TMath::Erf(((x-[1])/[2])))+(1-TMath::Erf((x-[4])/[5]))-1))"
//     //                                         "+((([4]-[1])/([7]-[1]))*[6]*[9]*[9]*((TMath::Erf((x-[1])/[2]))+((1-TMath::Erf((x-[7])/[8]))-1)))"
//     //                                         "+((([4]-[1])/([7]-[4]))*[6]*2*(1-[9])*[9]*((TMath::Erf((x-[4])/[5]))+(1-TMath::Erf((x-[7])/[8]))-1)))");
//     // func3peak->SetNpx(npx);
//     // func3peak->SetLineColor(kBlue);
//     // func3peak_1peall->SetNpx(npx);
//     // func3peak_2bs->SetNpx(npx);
//     // func3peak_2bs->SetLineColor(kCyan+1);
//     // func3peak_bs->SetNpx(npx);
//     // func3peak_bs->SetLineColor(kCyan+1);
//     // func3peak_bs->SetLineStyle(2);

//     // 1ピークフィット用の関数を定義 (funcpedを流用)
//     TF1 *func_1peak = funcped;
//     func_1peak->SetLineColor(kRed); // 色を変更
//     func_1peak->SetLineStyle(1);

//     TF1 *func, *func1peall, *funcbs;

//     // switch文を1ピーク用に修正
//     switch(nPeaks) {
//         case 1:
//              func       = func_1peak;
//              func1peall = func_1peak; // 1ピークのため、全体も同じ関数
//              funcbs     = nullptr;   // バックグラウンドなし
//              break;
//         case 3:
//         case 2:
//         default:
//              // 既存の2ピーク・3ピーク処理はコメントアウト
//              // func       = func2peak;
//              // func1peall = func2peak_1peall;
//              // funcbs     = func2peak_bs;
//              return; // 1ピーク以外は処理を中断
//     }

//     // 1ピーク用のスタイル設定
//     func->SetLineColor(kBlue);
//     func1peall->SetLineColor(kRed);
//     // funcbsはnullptrなので設定しない
//     func->SetLineStyle(1);
//     func1peall->SetLineStyle(2);
//     // funcbs->SetLineStyle(2); // nullptrなので設定しない

//     ///////////////////////
//     //// PEAK FINDING /////
//     ///////////////////////
    
//     /* should modify these initial values... */
//     Double_t peakx_pre[3] = {440,470,500};//{1500, 1600, 1700};//{0, 2.4, 4.8};//[pC]
//     Double_t peaky_pre[3] = {h->GetBinContent(h->FindBin(peakx_pre[0])), h->GetBinContent(h->FindBin(peakx_pre[1])), h->GetBinContent(h->FindBin(peakx_pre[1]))} ;

// #ifdef PEAKFINDING
//     TSpectrum s(nPeaks);
//     s.Search(h,4,"goff",0.001);
//     Int_t FoundPeaks = s.GetNPeaks();

//     for(int iPeak=0;iPeak<TMath::Min(nPeaks,FoundPeaks);iPeak++){
//      peakx_pre[iPeak]=s.GetPositionX()[iPeak];
//      peaky_pre[iPeak]=s.GetPositionY()[iPeak];
//      std::cout << "Peak Position " << iPeak << " is X: "<<peakx_pre[iPeak] <<" Y:"<< peaky_pre[iPeak] << std::endl;
//     }
//     if (nPeaks == 1) { // 1ピークの場合の処理を追加
//       if (FoundPeaks > 0) {
//         peakx_pre[0] = s.GetPositionX()[0];
//         peaky_pre[0] = s.GetPositionY()[0];
//       }
//     } else 
//     if (peakx_pre[0] > peakx_pre[1]) {
//       peakx_pre[2] = peakx_pre[0];
//       peakx_pre[0] = peakx_pre[1];
//       peakx_pre[1] = peakx_pre[2];
//       peaky_pre[2] = peaky_pre[0];
//       peaky_pre[0] = peaky_pre[1];
//       peaky_pre[1] = peaky_pre[2];
//     }
// #endif

//     // 1ピークフィットのためのフィット範囲を設定
//     h->Fit(func_1peak,"NQ0","",peakx_pre[0]*0.8,peakx_pre[0]*1.2);
    
//     // フィット後のパラメータをセット
//     Double_t fit_scale = func_1peak->GetParameter(0);
//     Double_t fit_peak  = func_1peak->GetParameter(1);
//     Double_t fit_sigma = func_1peak->GetParameter(2);
    
//     // 不要なRebin処理をコメントアウト
//     // if (rebin>1) {
//     //  h->Rebin(rebin);
//     //  funcped->SetParameter(0, funcped->GetParameter(0)*(Double_t)(rebin));
//     //  func1pe->SetParameter(0, func1pe->GetParameter(0)*(Double_t)(rebin));
//     //  funcped->SetParError(0, funcped->GetParError(0)*(Double_t)(rebin));
//     //  func1pe->SetParError(0, func1pe->GetParError(0)*(Double_t)(rebin));
//     // }

//     // パラメータ名と制限の設定を1ピーク用に修正
//     // ... 修正が複雑なため、ここは元のコードをコメントアウトし、簡潔な処理に置き換え ...
//     // func_1peak->SetParName(0, "Scale");
//     // func_1peak->SetParName(1, "Peak");
//     // func_1peak->SetParName(2, "Sigma");
//     // func_1peak->SetParLimits(2, 0, 100); // Sigmaを正の値に制限
    
//     // h->Fit(func_1peak, "B", "", h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());

//     // 1ピーク用のResultパラメータ計算
//     // 不要な変数をコメントアウト
//     // Double_t areaped      = 0;
//     // Double_t area1pe      = 0;
//     // Double_t areabs       = 0;
//     // Double_t area1peall   = 0;

//     // 1ピーク用の物理量計算
//     Double_t peakx = func_1peak->GetMaximumX();
//     Double_t peaky = func_1peak->GetMaximum();
//     Double_t FWHMlow = peakx - func_1peak->GetX(peaky*0.5, h->GetXaxis()->GetXmin(), peakx);
//     Double_t FWHMhigh = func_1peak->GetX(peaky*0.5, peakx, h->GetXaxis()->GetXmax()) - peakx;
//     Double_t FWHM = FWHMhigh + FWHMlow;
//     Double_t sigma = TMath::Abs(func_1peak->GetParameter(2));
//     Double_t gainpeakx = peakx/1.60217662e-19*1.e-12;
//     Double_t chi2 = func_1peak->GetChisquare();
//     Double_t ndf = (Double_t)func_1peak->GetNDF();

//     // 1ピークフィットに不要な計算をコメントアウト
//     // Double_t pv = -1;
//     // Double_t areabs = 0;
//     // Double_t areasigma = 0;
    
//     // Result parameters
//     // Double_t mean = func_1peak->Mean(h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());
//     // Double_t variance = func_1peak->Variance(h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());
    
//     // 不要な変数にデフォルト値を設定
//     Double_t mean = -1;
//     Double_t variance = -1;
//     Double_t pv = -1;
//     Double_t areaped = -1;
//     Double_t area1pe = -1;
//     Double_t area1peall = -1;
//     Double_t areabs = -1;
//     Double_t areasigmalow = -1;
//     Double_t areasigmahigh = -1;
//     Double_t areasigma = -1;
//     Double_t peak1pex = -1;
//     Double_t peak1pesigma = -1;
//     Double_t peak1pexerr = -1;
//     Double_t peak1pesigmaerr = -1;
//     Double_t gainpeak1pex = -1;
    
//     peakval = peakx; // peakvalを1ピークのピーク位置に設定

//     std::cout << std::endl << std::endl;
//     std::cout << "=========== RESULT ==========" << std::endl;
//     std::cout << " peakx            = " << " " <<  peakx         << " "      << std::endl;
//     std::cout << " FWHM             = " << " " <<  FWHM          << " " << " (" <<  (peakx!=0?FWHM/peakx*100.:-1)  << " [%])"  << std::endl;
//     std::cout << " sigma            = " << " " <<  sigma         << " " << " (" <<  (peakx!=0?sigma/peakx*100.:-1) << " [%])" << std::endl;
//     std::cout << " Gain             = " << " " <<  gainpeakx     << " "      << std::endl;
//     std::cout << " Chi2/NDF         = " << " " <<  chi2/ndf      << " "      << std::endl;
//     std::cout << "=============================" << std::endl;

//     gainpeak = gainpeakx;
//     respeak = peakx;
//     respv = pv;
//     resfwhm = (peakx!=0?FWHM/peakx*100.:-1);
//     ressigma = (peakx!=0?sigma/peakx*100.:-1);
//     bsratio = 0; // 1ピークフィットのため0に固定
//     chi2ndf = ndf!=0?chi2/ndf:0;

//     // 1ピークフィットのため、不要な描画をコメントアウト
//     // if (!DrawEff) {
//     //  gPad->SetLogy(1);
//     //  h->Draw("e1");
//     //  funcped ->Draw("same");
//     //  func1pe ->Draw("same");
//     //  if (nPeaks>2) { func2pe->Draw("same"); }
//     //  funcbs  ->Draw("same");
//     //  if (nPeaks>2) { func3peak_bs->Draw("same"); }
//     //  func1peall->Draw("same");
//     //  func->Draw("same");
//     //  gPad->Update();
//     // } else if (DrawEff) {
//     //   // ... 効率グラフの描画コード ...
//     // }
    
//     // 1ピークの関数だけを描画
//     gPad->SetLogy(1);
//     h->Draw("e1");
//     func->Draw("same");
// }