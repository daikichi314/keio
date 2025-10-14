//ttsfit.C
//By Y.Nishimura

#define header_cxx
#include <vector>
#include <string>
//#include <iostream>
//#include <fstream>

#include "./header.h"
#define nhv 8
#define inputpath "./"
#define outputpdfpath "./"

#include "TROOT.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TTree.h"
#include "TSpectrum.h"
#include "TMath.h"
#include "TH1.h"
#include "TF1.h"
#include "TLine.h"
#include "TPaveText.h"
#include "TPaveStats.h"
#include "TLatex.h"

double FWHM;
double chi2ndf;

//fit histogram and return the result
Double_t ExpGaus(Double_t *x, Double_t *par)
{
   //http://pibeta.phys.virginia.edu/~pibeta/docs/publications/penny_diss/node35.html
   //par[0] : height of Gaussian
   //par[1] : peak position
   //par[2] : sigma of Gaussian
   //par[3] : transition point between gaussian and exponential

   Double_t fitval;
   //if (par[1] > par[3]) par[3] = par[1];

   //if(x[0] > par[1] + par[3]) { // left tail
   if(x[0] < par[1] + par[3]) { // right tail
      if(par[2] != 0) {
         fitval = par[0] * TMath::Exp(-1 * (x[0] - par[1]) * (x[0] - par[1]) / 2 / par[2] / par[2]);
      } else {
         fitval = 0;
      }
   } else {
      if(par[2] != 0) {
         fitval = par[0] * TMath::Exp(par[3] / par[2] / par[2] * (par[3] / 2 - (x[0] - par[1])));
      } else {
         fitval = 0;
      }
   }

   return fitval;
}

Double_t EMG(Double_t *x, Double_t *par)
{
   return 0.5*par[3]*exp(0.5*par[3]*(2*par[0]+par[3]*par[2]*par[2]-2.*x[0]))
          *TMath::Erfc((par[0]+par[3]*par[2]*par[2]-x[0])/(sqrt(2.)*par[2]))*par[1];
}

Double_t GetFWHM(TF1 *f)
{

  return
  f->GetX(f->GetMaximum()*0.5,f->GetMaximumX(),f->GetXaxis()->GetXmax())
  -
  f->GetX(f->GetMaximum()*0.5,f->GetXaxis()->GetXmin(),f->GetMaximumX());

}

Double_t GetPeak(TF1 *f)
{

  return
  f->GetMaximumX();

}


void ttshistofit(TH1 *h)
{

  if (h->GetEntries() < 10) { return; }
  Bool_t IsAsymGaus = kFALSE;
  Bool_t IsEMG = kTRUE;//FALSE;
 // Bool_t IsEMG = kFALSE;
  Bool_t IsExpGaus = kFALSE;

  //h->GetXaxis()->SetRangeUser(-50,100);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);

  Double_t var[10], varerr[10];

  //Step 0 : Gaussian
  std::cout << " ============ Gaussian ============ " << std::endl;
  TF1 *fgaus = new TF1("fgaus", "gaus", -500, 500);
     fgaus->SetNpx(9000);
     fgaus->SetLineColor(kCyan);
     fgaus->SetLineWidth(1);
   
     //fgaus->SetParameter(1, h->GetBinCenter(h->GetMaximumBin()));
     fgaus->FixParameter(1, h->GetBinCenter(h->GetMaximumBin()));
     fgaus->SetParameter(2, 1.5);

     h->Fit(fgaus,"QBN","", h->GetBinCenter(h->GetMaximumBin())-12, h->GetBinCenter(h->GetMaximumBin())+12);
     //h->Fit(fgaus,"","", h->GetBinCenter(h->GetMaximumBin())-12, h->GetBinCenter(h->GetMaximumBin())+5);
     fgaus->ReleaseParameter(1);
     h->Fit(fgaus,"N","", h->GetBinCenter(h->GetMaximumBin())-12, h->GetBinCenter(h->GetMaximumBin())+25);
     for (Int_t ipar = 0; ipar < 3; ipar++) {
        var[ipar]    = fgaus->GetParameter(ipar);
        varerr[ipar] = fgaus->GetParError(ipar);
     }

  TF1 *fasymgaus = 0;
  if (IsAsymGaus) {
  //Step 2 : Asymmetric Gaussian
  std::cout << " ============ Asymmetric Gaussian ============ " << std::endl;
     fasymgaus = new TF1("fasymgaus", "[0]*((x<=[1])*TMath::Gaus(x,[1],[2],0)+(x>[1])*TMath::Gaus(x,[1],[3],0))", -500, 500);
     fasymgaus->SetNpx(9000);
     fasymgaus->SetLineColor(kMagenta);
     fasymgaus->SetLineWidth(1);
     fasymgaus->SetParameter(0, fgaus->GetParameter(0));
     fasymgaus->SetParameter(1, fgaus->GetParameter(1));
     fasymgaus->SetParameter(2, fgaus->GetParameter(2));
     fasymgaus->SetParameter(3, fgaus->GetParameter(2));
     fasymgaus->SetParName(0, "Scale");
     fasymgaus->SetParName(1, "Peak");
     fasymgaus->SetParName(2, "#sigma_{Left}");
     fasymgaus->SetParName(3, "#sigma_{Right}");
     //h->Fit(fasymgaus,"","", h->GetBinCenter(h->GetMaximumBin())-20, h->GetBinCenter(h->GetMaximumBin())+30);
     //h->Fit(fasymgaus,"","", h->GetBinCenter(h->GetMaximumBin())-20, h->GetBinCenter(h->GetMaximumBin())+5);
     h->Fit(fasymgaus,"N+","", h->GetBinCenter(h->GetMaximumBin())-20, h->GetBinCenter(h->GetMaximumBin())+25);
     for (Int_t ipar = 0; ipar < 4; ipar++) {
        var[ipar]    = fasymgaus->GetParameter(ipar);
        varerr[ipar] = fasymgaus->GetParError(ipar);
     }
     std::cout << "Scale     " <<  var[0]    << "[ns]" << std::endl;
     std::cout << "Peak      " <<  var[1]    << "[ns]" << std::endl;
     std::cout << "Sigma     " <<  var[2]    << "[ns]" << std::endl;
     std::cout << "Sigma     " <<  var[3]    << "[ns]" << std::endl;
     std::cout << "PeakErr   " <<  varerr[1] << "[ns]" << std::endl;
     std::cout << "SigmaErr  " <<  varerr[2] << "[ns]" << std::endl;
     std::cout << "SigmaErr  " <<  varerr[3] << "[ns]" << std::endl;
   }

  //TF1 *f2 = new TF1("f2","expo",230,250);
  //h->Fit("f2","","",-15,-2);

  TF1 *emg = 0;
  if (IsEMG) {
  std::cout << " ============ Exponential Modified Gaussian ============ " << std::endl;
  emg = new TF1("emg", EMG, -100, 500, 4);
  emg->SetLineColor(kRed);
  emg->SetLineStyle(2);
  emg->SetNpx(9000);
  emg->SetParName(0, "#mu");
  emg->SetParName(1, "#lambda");
  emg->SetParName(2, "#sigma");
  emg->SetParName(3, "#gamma");
  emg->SetParameter(0,var[1]);
  emg->SetParameter(1,var[0]*10.);
  //emg->SetParameter(1,var[0]);
  emg->SetParameter(2,var[2]*0.7);
  emg->SetParameter(3,var[2]*0.1);
  //emg->SetParameter(3,var[2]);
  emg->SetParLimits(0, h->GetBinCenter(h->GetMaximumBin())-3, h->GetBinCenter(h->GetMaximumBin())+3);
  emg->SetParLimits(1, 1, 1000000);
  emg->SetParLimits(2, 0.3, 5);
  emg->SetParLimits(3, 0.1, 5);

  //h->Fit(emg,"","", h->GetBinCenter(h->GetMaximumBin())-20, h->GetBinCenter(h->GetMaximumBin())+30);
  //h->Fit(emg,"N","", h->GetBinCenter(h->GetMaximumBin())-20, h->GetBinCenter(h->GetMaximumBin())+5);
  //h->Fit(emg,"N","", h->GetBinCenter(h->GetMaximumBin())-20, h->GetBinCenter(h->GetMaximumBin())+5);
  emg->FixParameter(2,var[2]*0.7);
  h->Fit(emg,"BN0","", h->GetBinCenter(h->GetMaximumBin())-20, h->GetBinCenter(h->GetMaximumBin())+5);
  emg->ReleaseParameter(1);
  emg->ReleaseParameter(2);
  emg->ReleaseParameter(3);
  h->Fit(emg,"B+","+", h->GetBinCenter(h->GetMaximumBin())-20, h->GetBinCenter(h->GetMaximumBin())+25);
  }

  TF1 *expgaus = 0;
  if (IsExpGaus) {
  std::cout << " ============ Exponential + Gaussian ============ " << std::endl;
  expgaus = new TF1("expgaus", ExpGaus, -100, 500,4);
  expgaus->SetLineColor(kBlue);
  expgaus->SetNpx(9000);
  expgaus->SetParName(0, "Scale");
  expgaus->SetParName(1, "Peak");
  expgaus->SetParName(2, "#sigma");
  expgaus->SetParName(3, "TransPoint");
  expgaus->SetParameter(0,var[0]);
  expgaus->SetParLimits(0, var[0]*0.1, var[0]*1000.);
  expgaus->SetParameter(1,var[1]);
  expgaus->SetParLimits(1, var[1]-var[2]*1., var[1]+var[3]*1.);
  expgaus->SetParameter(2,var[3]);
  expgaus->SetParLimits(2, var[2]*0.5, var[3]*1.5);
  expgaus->SetParameter(3,var[1]+var[3]*2.);
  expgaus->SetParLimits(3,var[1], var[1]+var[3]*5.);
  h->Fit(expgaus,"BNQ","N", h->GetBinCenter(h->GetMaximumBin())-20, h->GetBinCenter(h->GetMaximumBin())+25);
  expgaus->FixParameter(0,expgaus->GetParameter(0));
  expgaus->FixParameter(1,expgaus->GetParameter(1));
  expgaus->FixParameter(2,expgaus->GetParameter(2));
  h->Fit(expgaus,"BNQ","N", h->GetBinCenter(h->GetMaximumBin())-20, h->GetBinCenter(h->GetMaximumBin())+25);
  expgaus->ReleaseParameter(0);
  expgaus->ReleaseParameter(1);
  expgaus->ReleaseParameter(2);
  //expgaus->SetParLimits(3,emg->GetMaximumX(), var[1]+var[3]*5.);
  //expgaus->SetParLimits(3, h->GetBinCenter(h->GetMaximumBin()), h->GetBinCenter(h->GetMaximumBin())*var[2]*4.);
  //expgaus->SetParameter(3,2.1);
  //expgaus->FixParameter(3,f2->GetParameter(1));
  //h->Fit(expgaus,"","", h->GetBinCenter(h->GetMaximumBin())-20, h->GetBinCenter(h->GetMaximumBin())+30);
  //h->Fit(expgaus,"N","", h->GetBinCenter(h->GetMaximumBin())-20, h->GetBinCenter(h->GetMaximumBin())+5);
  h->Fit(expgaus,"N+","", h->GetBinCenter(h->GetMaximumBin())-20, h->GetBinCenter(h->GetMaximumBin())+25);


  h->Draw();
  // BLPMT
  h->GetXaxis()->SetRangeUser(h->GetBinCenter(h->GetMaximumBin())-10, h->GetBinCenter(h->GetMaximumBin())+20);
  // MCPPMTPMT
  //h->GetXaxis()->SetRangeUser(-30, 70);

/*
  fgaus->Draw("same");
  fasymgaus->Draw("same");
  expgaus->Draw("same");
*/
  emg->Draw("same");
  }
  
  std::cout << " ======================================================= " << std::endl;
  std::cout << " ============          Result               ============ " << std::endl;
  std::cout << " ============ sigma, FWHM, peak, chi2, ndf  ============ " << std::endl;
  std::cout << "FitRes " << fgaus->GetName()     << " " << fgaus->GetParameter(2)      << "[ns] " << GetFWHM(fgaus)   << "[ns] " << GetPeak(fgaus)   << "[ns] "<< " CHI2:" << fgaus->GetChisquare() << " NDF:" << fgaus->GetNDF() << std::endl;
  if (IsAsymGaus)
  std::cout << "FitRes " << fasymgaus->GetName() << " (" << fasymgaus->GetParameter(2) << "[ns] " << ", " << fasymgaus->GetParameter(3) << "[ns]) " << GetFWHM(fasymgaus) << "[ns] " << GetPeak(fasymgaus) << "[ns] "<< " CHI2:" << fasymgaus->GetChisquare() << " NDF:" << fasymgaus->GetNDF() << std::endl;
  if (IsExpGaus)
  std::cout << "FitRes " << expgaus->GetName()   << " " << expgaus->GetParameter(2)    << "[ns] " << GetFWHM(expgaus) << "[ns] " << GetPeak(expgaus) << "[ns] "<< " CHI2:" << expgaus->GetChisquare() << " NDF:" << expgaus->GetNDF() << std::endl;
  if (IsEMG)
  std::cout << "FitRes " << emg->GetName()       << " " << emg->GetParameter(2)        << "[ns] " << GetFWHM(emg)     << "[ns] " << GetPeak(emg)     << "[ns] "<< " CHI2:" << emg->GetChisquare() << " NDF:" << emg->GetNDF() << std::endl;
  std::cout << " ======================================================= " << std::endl;

  TF1 *fout = emg;//expgaus;
  std::cout << "TTSOUT " << ", " << fout->GetParameter(2) << ", " << GetFWHM(fout) << ", " << GetPeak(fout) << std::endl;
  if (emg!=0)
  std::cout << "TTSERROUT " << ", " << ((emg->GetParameter(2)!=0)?emg->GetParError(2)*GetFWHM(fout)/emg->GetParameter(2):0) << ", " << emg->GetParError(0) << std::endl;

     gPad->Update();

     Float_t st0_lower_left_x = 0.6;
     Float_t st0_lower_left_y = 0.35;
     Float_t st_Width = 0.4;
     Float_t st_Height = 0.65;

     TPaveStats *st0 = (TPaveStats*)gPad->GetPrimitive("stats");
     //TPaveStats *st0 = (TPaveStats*)h->FindObject("stats");
     if (st0) {
       st0->SetName("result");
       TList *listOfLines = st0->GetListOfLines();
       
       gStyle->SetTitleXOffset(0.);
       gStyle->SetTitleX(0.);
       gStyle->SetTitleAlign(13);//left a=10*h+v, h=1,2,3 for left, center, top,  v=1,2,3 for bottom, middle, top
       st0->SetX1NDC(st0_lower_left_x);
       st0->SetY1NDC(st0_lower_left_y);
       st0->SetX2NDC(st0_lower_left_x + st_Width);
       st0->SetY2NDC(st0_lower_left_y + st_Height);
       TLatex *myt[4];
       myt[0] = new TLatex(0,0,Form("Single PE RESULT"));
       myt[1] = new TLatex(0,0,Form("Relative Peak  = %2.2f [ns]",GetPeak(fout)));
       myt[2] = new TLatex(0,0,Form("Upper Sigma = %2.2f [ns]",fout->GetParameter(2)));
       myt[3] = new TLatex(0,0,Form("FWHM = %3.2f [ns]", GetFWHM(fout)));
       FWHM = GetFWHM(fout);
       chi2ndf = fout->GetChisquare()/(double)fout->GetNDF();
       for (Int_t i = 0; i < 4; i++) {
          myt[i] ->SetTextFont(42);
          myt[i] ->SetTextSize(0.04);
          myt[i] ->SetTextColor(kBlue);
          listOfLines->Add(myt[i]);
       }
          myt[0] ->SetTextSize(0.02);
       h->SetStats(0);
       gPad->Modified();
     }

  //h->GetXaxis()->SetRangeUser(-15,20);
  h->GetXaxis()->SetRangeUser(h->GetBinCenter(h->GetMaximumBin())-15, h->GetBinCenter(h->GetMaximumBin())+20);

}

std::vector<std::string>* GetSerial(header hd) {

      hd.LoadTree(0);
      hd.fChain->GetEntry(0);
      for (int i = 0; i < hd.serial->size(); i++) {
         std::cout << hd.run << " " << hd.serial->at(i).c_str() << " " << hd.HVS[i] << std::endl;
      }
   return hd.serial;
}

#ifdef __CINT__
int ttsfit(int run = 123, const char* outfilename = "test.root", bool Update = kFALSE)
{
#else
int main(int argc, char* argv[]){
   if (argc < 2) {
      std::cout << "USAGE : ./ttsfit (Run) (Output Root filename) (1:Add in existing root, 0:Recreate new)" << std::endl;
      return 0;
   }
   int run = atoi(argv[1]);
   const char* outfilename = argv[2];
   bool Update = kFALSE;
   if (argc > 3) {
       Update = atoi(argv[4]);
   }
#endif
  int ch = 0; 
  double offset = 72;

   const Int_t oldErrorIgnoreLevel = gErrorIgnoreLevel;
   gErrorIgnoreLevel = 5000;

  TCanvas *c = new TCanvas("c","c",800,600);

  TFile *file = TFile::Open(Form("./run%06d.root",run));

  TString serial[nhv];
  Double_t HVS[nhv];
  Double_t THR[nhv];
  Int_t LDA, LDW;
   header hd((TTree*)(file->Get("header")));
//   std::vector<std::string>* serial = GetSerial(hd);
//   //header *hd = new header((TTree*)(file->Get("header")));
         hd.LoadTree(0);
         hd.fChain->GetEntry(0);
      for (int i = 0; i < nhv; i++) {
         std::cout << hd.run << " " << hd.serial->at(i).c_str() << " " << hd.HVS[i] << std::endl;
         if (i<hd.serial->size()) {
            serial[i] = hd.serial->at(i);
         } else {
            serial[i] = "";
         }
         HVS[i] = hd.HVS[i];
         THR[i] = hd.THR[i];
      }
      //result.run = run;
      //result.unixtime = hd.start;
      LDA = hd.LSamp;
      LDW = hd.LSwidth;
  
  TTree *tdc = (TTree*)(file->Get("tdc"));
  c->Divide(3,3);
  c->Print(Form("%s/run%06d_tdc.pdf[", outputpdfpath, run));


  for (ch = 0; ch < 7; ch++) {
     c->cd(ch+1);
     tdc->Draw(Form("(tdcL%d-tdcL%d[0])-(%f)>>h%02d(%d,%d,%d)",ch+1,8,offset,ch,4800,-300,300));
 
     TH1* h = (TH1*)(gDirectory->Get(Form("h%02d",ch)));
     h->SetTitle(Form("Ch%01d %s %4.1fV Run%d", ch+1, serial[ch].Data(), HVS[ch], run));

     FWHM = 0;
     chi2ndf = 0;
     ttshistofit(h);
     std::cout << "BBB RESULT : " << hd.run << ", " << ch << ", " << hd.serial->at(ch).c_str() << ", " << hd.HVS[ch] << ", " << FWHM << ", " << chi2ndf << std::endl;

     //c->SaveAs(Form("1pe_tts_Run%06d_ch%02d.root",run,ch));
     //c->SaveAs(Form("1pe_tts_Run%06d.root",run));
  }

     c->Print(Form("%s/run%06d_tdc.pdf", outputpdfpath, run));
  c->Print(Form("%s/run%06d_tdc.pdf]", outputpdfpath, run));

  return 0;

}
