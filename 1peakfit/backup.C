//$Id$
//By Y.Nishimura

//#define __CINT__
#define PEAKFINDING

#define header_cxx
#include <vector>
#include <string>
//#include <iostream>
#include <fstream>


#include "TROOT.h"
#include "TFile.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TTree.h"
#include "TList.h"
#include "TSpectrum.h"
#include "TMath.h"
#include "TH1.h"
#include "TF1.h"
#include "TLine.h"
#include "TPaveText.h"
#include "TPaveStats.h"
#include "TLatex.h"

#include "adcresult.h"

//#define inputpath "."
#define inputpath "/home/daiki/lab/data/20250809"
// //#define outputpdfpath "."
#define outputpdfpath "/home/daiki/lab/data/20250809"
// //#define outputtxtpath "."
#define outputtxtpath "/home/daiki/lab/data/20250809"
// #define inputpath (TString(getenv("HOME")) + "/lab/data/20250809").Data()
// #define outputpdfpath "."
// #define outputtxtpath (TString(getenv("HOME")) + "/lab/data/20250809").Data()
   
static FitResultSPE result;
Double_t peakval;

static double gainpeak, respeak, respeakerr, respv, resfwhm, ressigma, chi2ndf, bsratio;

void FitPedPe(TH1* h, const char* name = "", Int_t nPeaks = 3, Int_t rebin = 4, Bool_t DrawEff = kFALSE) {

   const Int_t npx = 1000;
   //Int_t rebin = 4; // After pedestal fit
   //Int_t nPeaks = 3;
   //Bool_t DrawEff = kFALSE;

   gStyle->SetOptFit(1);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);

   //////////////////////////////
   ////  Define Functions  //////
   //////////////////////////////
   TF1 *funcped    = new TF1(Form("funcped%s",name),"gaus", h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());
   TF1 *func1pe    = new TF1(Form("func1pe%s",name),"gaus", h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());
   TF1 *func2pe    = new TF1(Form("func2pe%s",name),"gaus", h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());
   funcped  ->SetLineColor(kMagenta);
   func1pe  ->SetLineColor(kGreen+1);
   func2pe  ->SetLineColor(kGreen);
   funcped ->SetLineStyle(3);
   func1pe ->SetLineStyle(2);
   func2pe ->SetLineStyle(2);
   func1pe ->SetLineWidth(1);
   func2pe ->SetLineWidth(1);
   funcped ->SetNpx(npx);
   func1pe ->SetNpx(npx);
   func2pe ->SetNpx(npx);

   TF1 *func2peak        = new TF1(Form("func2peak%s",name),"gaus(0)+gaus(3)+(0.5*(([6]*[3]*((TMath::Erf(((x-[1])/[2])))+(1-TMath::Erf((x-[4])/[5]))-1))))");
   TF1 *func2peak_1peall = new TF1(Form("func1peall%s",name),"gaus(3)+(0.5*(([6]*[3]*((TMath::Erf(((x-[1])/[2])))+(1-TMath::Erf((x-[4])/[5]))-1))))");
   TF1 *func2peak_bs     = new TF1(Form("func1bs%s",name),"(0.5*(([6]*[3]*((TMath::Erf(((x-[1])/[2])))+(1-TMath::Erf((x-[4])/[5]))-1))))");
   func2peak     ->SetNpx(npx);
   func2peak_1peall->SetNpx(npx);
   func2peak_bs  ->SetNpx(npx);

   TF1 *func3peak       = new TF1(Form("func3peak%s",name),"gaus(0)+gaus(3)+gaus(6)+(0.5*(([9]*[3]*((TMath::Erf(((x-[1])/[2])))+(1-TMath::Erf((x-[4])/[5]))-1))"
                                            "+((([4]-[1])/([7]-[1]))*[6]*[9]*[9]*((TMath::Erf((x-[1])/[2]))+((1-TMath::Erf((x-[7])/[8]))-1)))" 
                                            "+((([4]-[1])/([7]-[4]))*[6]*2*(1-[9])*[9]*((TMath::Erf((x-[4])/[5]))+(1-TMath::Erf((x-[7])/[8]))-1))))");
   //TF1 *func3peak_1peall = func2peak_1peall;//new TF1(Form("func1peall2%s",name),"gaus(3)+(0.5*(([6]*[3]*((TMath::Erf(((x-[1])/[2])))+(1-TMath::Erf((x-[4])/[5]))-1))))");//
   TF1 *func3peak_1peall = new TF1(Form("func3peak_1peall%s",name),"gaus(3)+(0.5*(([9]*[3]*((TMath::Erf(((x-[1])/[2])))+(1-TMath::Erf((x-[4])/[5]))-1))))");
   //TF1 *func3peak_bs     = func2peak_bs;    //new TF1(Form("func1bs2%s",name),"(0.5*(([6]*[3]*((TMath::Erf(((x-[1])/[2])))+(1-TMath::Erf((x-[4])/[5]))-1))))");           //
   TF1 *func3peak_bs     = new TF1(Form("func3peak_1bs%s",name),"(0.5*(([9]*[3]*((TMath::Erf(((x-[1])/[2])))+(1-TMath::Erf((x-[4])/[5]))-1))))");
   TF1 *func3peak_2bs    = new TF1(Form("func2bs%s",name),"0.5*(([9]*[3]*((TMath::Erf(((x-[1])/[2])))+(1-TMath::Erf((x-[4])/[5]))-1))"
                                            "+((([4]-[1])/([7]-[1]))*[6]*[9]*[9]*((TMath::Erf((x-[1])/[2]))+((1-TMath::Erf((x-[7])/[8]))-1)))" 
                                            "+((([4]-[1])/([7]-[4]))*[6]*2*(1-[9])*[9]*((TMath::Erf((x-[4])/[5]))+(1-TMath::Erf((x-[7])/[8]))-1)))");
   func3peak->SetNpx(npx);
   func3peak->SetLineColor(kBlue);
   func3peak_1peall->SetNpx(npx);
   func3peak_2bs->SetNpx(npx);
   func3peak_2bs->SetLineColor(kCyan+1);
   func3peak_bs->SetNpx(npx);
   func3peak_bs->SetLineColor(kCyan+1);
   func3peak_bs->SetLineStyle(2);

   TF1 *func, *func1peall, *funcbs;

  switch(nPeaks) {
     case 3:
            func       = func3peak;
            func1peall = func3peak_1peall;
            funcbs     = func3peak_2bs;
            break;
     default:
     case 2:
            func       = func2peak;
            func1peall = func2peak_1peall;
            funcbs     = func2peak_bs;
            break;
  }

   func      ->SetLineColor(kBlue);
   func1peall->SetLineColor(kRed);
   funcbs    ->SetLineColor(kCyan+1);
   func      ->SetLineStyle(1);
   func1peall->SetLineStyle(2);
   funcbs    ->SetLineStyle(2);




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

  //for(int iini=0; iini<2; iini++){
  //  peakx_pre[iini]=0;
  //  peaky_pre[iini]=0;
  //}
  for(int iPeak=0;iPeak<TMath::Min(nPeaks,FoundPeaks);iPeak++){
    peakx_pre[iPeak]=s.GetPositionX()[iPeak];
    peaky_pre[iPeak]=s.GetPositionY()[iPeak];
    std::cout << "Peak Position " << iPeak << " is X: "<<peakx_pre[iPeak] <<" Y:"<< peaky_pre[iPeak] << std::endl;
  }
  if (nPeaks == 1) {
     peakx_pre[1] = 2.4;
     peaky_pre[1] = peaky_pre[0]*0.005;
  } else 
  if (peakx_pre[0] > peakx_pre[1]) {
     peakx_pre[2] = peakx_pre[0];
     peakx_pre[0] = peakx_pre[1];
     peakx_pre[1] = peakx_pre[2];
     peaky_pre[2] = peaky_pre[0];
     peaky_pre[0] = peaky_pre[1];
     peaky_pre[1] = peaky_pre[2];

  }
  //if (FoundPeaks < nPeaks) {
     peakx_pre[2] = peakx_pre[0] + (peakx_pre[1]-peakx_pre[0])*2;
     peaky_pre[2] = h->GetBinContent(h->FindBin(peakx_pre[2]));
  //}
#endif

  h->Fit(funcped,"NQ0","",peakx_pre[0]*0.8,peakx_pre[0]+(peakx_pre[1]-peakx_pre[0])*0.2);
  //h->Fit(func1pe,"NQ0","", 0.5,3);//h->GetXaxis()->GetXmax());
  h->Fit(func1pe,"NQ0","", peakx_pre[1]-(peakx_pre[1]-peakx_pre[0])*0.2,peakx_pre[1]+(peakx_pre[1]-peakx_pre[0])*0.2);//h->GetXaxis()->GetXmax());
  for(int iSet=0; iSet<3;iSet++){
      func2peak->SetParameter(iSet  , funcped ->GetParameter(iSet));
      func2peak->SetParameter(iSet+3, func1pe ->GetParameter(iSet));
  }
  if (rebin>1) {
     h->Rebin(rebin);
     funcped->SetParameter(0, funcped->GetParameter(0)*(Double_t)(rebin));
     func1pe->SetParameter(0, func1pe->GetParameter(0)*(Double_t)(rebin));
     funcped->SetParError(0, funcped->GetParError(0)*(Double_t)(rebin));
     func1pe->SetParError(0, func1pe->GetParError(0)*(Double_t)(rebin));
  }

  for(int iSet=0; iSet<3/*nPeaks*/;iSet++){

    func3peak->SetParName(iSet*3,Form("Scale_{%d}",iSet));
    func3peak->SetParName(iSet*3+1,Form("Peak_{%d}",iSet));
    func3peak->SetParName(iSet*3+2,Form("#sigma_{%d}",iSet));
    func3peak_1peall->SetParName(iSet*3,Form("Scale_{%d}",iSet));
    func3peak_1peall->SetParName(iSet*3+1,Form("Peak_{%d}",iSet));
    func3peak_1peall->SetParName(iSet*3+2,Form("#sigma_{%d}",iSet));
    if (iSet>0) {
    func3peak_bs->SetParName(iSet*3,Form("Scale_{%d}",iSet));
    }
    func3peak_bs->SetParName(iSet*3+1,Form("Peak_{%d}",iSet));
    func3peak_bs->SetParName(iSet*3+2,Form("#sigma_{%d}",iSet));
    func3peak_2bs->SetParName(iSet*3,Form("Scale_{%d}",iSet));
    func3peak_2bs->SetParName(iSet*3+1,Form("Peak_{%d}",iSet));
    func3peak_2bs->SetParName(iSet*3+2,Form("#sigma_{%d}",iSet));
    if (iSet>1) continue;

    func2peak->SetParName(iSet*3,Form("Scale_{%d}",iSet));
    func2peak->SetParName(iSet*3+1,Form("Peak_{%d}",iSet));
    func2peak->SetParName(iSet*3+2,Form("#sigma_{%d}",iSet));
    func2peak_1peall->SetParName(iSet*3,Form("Scale_{%d}",iSet));
    func2peak_1peall->SetParName(iSet*3+1,Form("Peak_{%d}",iSet));
    func2peak_1peall->SetParName(iSet*3+2,Form("#sigma_{%d}",iSet));
    if (iSet>0) {
    func2peak_bs->SetParName(iSet*3,Form("Scale_{%d}",iSet));
    }
    func2peak_bs->SetParName(iSet*3+1,Form("Peak_{%d}",iSet));
    func2peak_bs->SetParName(iSet*3+2,Form("#sigma_{%d}",iSet));

    //func2peak->SetParameter(iSet*3,  peaky[iSet]);
    //func2peak->SetParameter(iSet*3+1,peakx[iSet]);
    //func2peak->SetParameter(iSet*3+2,peaky[iSet]*0.3);//unit pC
  }

  func2peak->SetParName(6, "1pe BS ratio");
  func2peak_1peall->SetParName(6, "1pe BS ratio");
  func2peak_bs->SetParName(6, "1pe BS ratio");
  func2peak->SetParLimits(6,0,1);
  func3peak->SetParName(9, "1pe BS ratio");
  func3peak->SetParLimits(9,0,1);

  //h->Fit(func2peak,(nPeaks==2?"":"NQ0"),"",h->GetXaxis()->GetXmin(), 2.5);//h->GetXaxis()->GetXmax());
  h->Fit(func2peak,(nPeaks==2?"":"NQ0"),"",h->GetXaxis()->GetXmin(), peakx_pre[1]*1.5);//h->GetXaxis()->GetXmax());

  for(int iSet=0; iSet<3;iSet++){
      funcped ->SetParameter(iSet, func2peak->GetParameter(iSet));
      func1pe ->SetParameter(iSet, func2peak->GetParameter(iSet+3));
      func2peak_bs  ->SetParameter(iSet, func2peak->GetParameter(iSet));
      func2peak_bs  ->SetParameter(iSet+3, func2peak->GetParameter(iSet+3));
      func2peak_1peall->SetParameter(iSet, func2peak->GetParameter(iSet));
      func2peak_1peall->SetParameter(iSet+3, func2peak->GetParameter(iSet+3));
      func3peak->SetParameter(iSet, func2peak->GetParameter(iSet));
      func3peak->SetParameter(iSet+3, func2peak->GetParameter(iSet+3));
  }
  func2peak_bs   ->SetParameter(6, func2peak->GetParameter(6));
  func2peak_1peall->SetParameter(6, func2peak->GetParameter(6));
  func3peak->SetParameter(9, func2peak->GetParameter(6));

  if (nPeaks == 3) {
     //Histograms should be rebinned if fitting failed
     func3peak->SetParameter(6, func2peak->GetParameter(3)*0.1);
     func3peak->SetParameter(7, func2peak->GetParameter(4)*2-func2peak->GetParameter(1));
     func3peak->SetParameter(8, TMath::Abs(func2peak->GetParameter(5))*TMath::Sqrt(2));
     func3peak->SetParLimits(8, TMath::Abs(func3peak->GetParameter(5))*0.5,TMath::Abs(func3peak->GetParameter(5))*TMath::Sqrt(2)*1.2);
     for(int iSet=0; iSet<6;iSet++){
         func3peak->FixParameter(iSet, func3peak->GetParameter(iSet));
     }
         //func3peak->FixParameter(9, func3peak->GetParameter(9));
         func3peak->FixParameter(9, 0.25); //201807
         func3peak->SetParLimits(7,func3peak->GetParameter(4)+func3peak->GetParameter(5),func3peak->GetParameter(4)*2-func3peak->GetParameter(1)+func3peak->GetParameter(5));
         func3peak->FixParameter(7,func3peak->GetParameter(4)*2-func3peak->GetParameter(1));
     h->Fit(func3peak,"B","",h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());
     
     // Set as 3 peaks
     for(int iSet=0; iSet<10;iSet++){
         func3peak->ReleaseParameter(iSet);
     }
     func3peak->SetParLimits(7,func3peak->GetParameter(7)-func3peak->GetParameter(8)*0.3,func3peak->GetParameter(7)+func3peak->GetParameter(8)*0.3);
     func3peak->SetParLimits(8,TMath::Abs(func3peak->GetParameter(8))*0.7,TMath::Abs(func3peak->GetParameter(8))*1.3);
     func3peak->SetParLimits(8,TMath::Abs(func3peak->GetParameter(5)),TMath::Abs(func3peak->GetParameter(5))*1.);//tmp
     func3peak->SetParLimits(9,0.15,0.5);
     //func3peak->FixParameter(9,0.3147);
     h->Fit(func3peak,"B","",h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());
     func3peak->FixParameter(7, func3peak->GetParameter(4)*2-func3peak->GetParameter(1));//tmp
     func3peak->SetParLimits(8,TMath::Abs(func3peak->GetParameter(8))*0.1,TMath::Abs(func3peak->GetParameter(8))*1.3);//tmp
     //func3peak->FixParameter(8, func3peak->GetParameter(5)*TMath::Sqrt(2));//tmp
     func3peak->FixParameter(8, func3peak->GetParameter(5));//tmp
     h->Fit(func3peak,"B","",-0.3, h->GetXaxis()->GetXmax());//tmp
     func3peak->ReleaseParameter(7);
     func3peak->ReleaseParameter(8);
     func3peak->SetParLimits(8,TMath::Abs(func3peak->GetParameter(8))*0.1,TMath::Abs(func3peak->GetParameter(8))*1.7);//tmp
     if (func3peak->GetParameter(3) < func3peak->GetParameter(6) || func3peak->GetParameter(3) < 0) {
        Double_t tmp = func3peak->GetParameter(6); 
        //func3peak->SetParameter(6, func3peak->GetParameter(3));
        func3peak->SetParameter(3, tmp); 
        func3peak->SetParameter(4, func3peak->GetParameter(7));
        //func3peak->SetParameter(7, func3peak->GetParameter(7)*2.);
        func3peak->FixParameter(6, 0); 
        func3peak->FixParameter(7, 0); 
        func3peak->FixParameter(8, 0); 
        std::cout << "1. Failed ... fitting again." << std::endl;
        h->Fit(func3peak,"B","",-0.3, h->GetXaxis()->GetXmax());//tmp
     }
     //h->Fit(func3peak,"B","",-0.3, 5);//h->GetXaxis()->GetXmax());//tmp
     //h->Fit(func3peak,"B","",-0.3, 6);//h->GetXaxis()->GetXmax());//tmp
     //h->Fit(func3peak,"B","",-0.3, h->GetXaxis()->GetXmax());//tmp
     //for(int iSet=0; iSet<10;iSet++){
     //    func3peak->ReleaseParameter(iSet);
     //}
     //h->Fit(func3peak,"","",h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());
     //h->Fit(func3peak,"","",h->GetXaxis()->GetXmin(), func3peak->GetParameter(7)+func3peak->GetParameter(8)*2.);
     //
     //BS fit
     func3peak->FixParameter(0, func3peak->GetParameter(0));
     func3peak->FixParameter(1, func3peak->GetParameter(1));
     func3peak->FixParameter(2, func3peak->GetParameter(2));
     if (func3peak->GetParameter(3) < func3peak->GetParameter(6) || func3peak->GetParameter(3) < 0) {
        func3peak->FixParameter(6, 0); 
        func3peak->FixParameter(7, 0); 
        func3peak->FixParameter(8, 0); 
     } else {
     func3peak->FixParameter(6, func3peak->GetParameter(6)); 
     func3peak->FixParameter(7, func3peak->GetParameter(7)); 
     func3peak->FixParameter(8, func3peak->GetParameter(8)); 
     }
     func3peak->FixParameter(3, func3peak->GetParameter(3));
     func3peak->FixParameter(4, func3peak->GetParameter(4));
     func3peak->FixParameter(5, func3peak->GetParameter(5));
     func3peak->SetParLimits(9,0.01,1.);
     //h->GetXaxis()->SetRangeUser(func3peak->GetMinimum(func3peak->GetParameter(1),func3peak->GetParameter(4))*0.9, (func3peak->GetParameter(4)-func3peak->GetParameter(5))*1.2);
     std::cout << "Fit BS Range DDD : " << func3peak->GetMinimumX(func3peak->GetParameter(1),func3peak->GetParameter(4))*0.9 << " - " << (func3peak->GetParameter(4)-func3peak->GetParameter(5))*1.2 << std::endl;
     h->Fit(func3peak,"B","",func3peak->GetMinimumX(func3peak->GetParameter(1),func3peak->GetParameter(4))*0.9, (func3peak->GetParameter(4)-func3peak->GetParameter(5))*1.2);
     //Final fit
     for(int iSet=3; iSet<9;iSet++){
         func3peak->ReleaseParameter(iSet);
         func3peak->SetParLimits(iSet, TMath::Min(func3peak->GetParameter(iSet)*0.8, func3peak->GetParameter(iSet)*1.2), TMath::Max(func3peak->GetParameter(iSet)*0.8, func3peak->GetParameter(iSet)*1.2));
     }
     //func3peak->FixParameter(0, func3peak->GetParameter(0));
     //func3peak->FixParameter(1, func3peak->GetParameter(1));
     //func3peak->FixParameter(2, func3peak->GetParameter(2));
     ////if (func3peak->GetParameter(3) < func3peak->GetParameter(6) || func3peak->GetParameter(3) < 0) {
     ////   func3peak->FixParameter(6, 0); 
     ////   func3peak->FixParameter(7, 0); 
     ////   func3peak->FixParameter(8, 0); 
     ////} else {
     //func3peak->FixParameter(6, func3peak->GetParameter(6)); 
     //func3peak->FixParameter(7, func3peak->GetParameter(7)); 
     //func3peak->FixParameter(8, func3peak->GetParameter(8)); 
     ////}
     func3peak->FixParameter(9, func3peak->GetParameter(9)); 
     //h->GetXaxis()->SetRangeUser(0.7, func3peak->GetParameter(7)+func3peak->GetParameter(8)*2.);
     h->Fit(func3peak,"B","",0.7, func3peak->GetParameter(7)+func3peak->GetParameter(8)*2.);
     if (func3peak->GetParameter(3) < func3peak->GetParameter(6) || func3peak->GetParameter(3) < 0) {
        func3peak->FixParameter(6, 0); 
        func3peak->FixParameter(7, 0); 
        func3peak->FixParameter(8, 0); 
        func3peak->FixParameter(9, func3peak->GetParameter(9)); 
        std::cout << "Failed ... fitting again." << std::endl;
        h->Fit(func3peak,"B","",0.7, func3peak->GetParameter(7)+func3peak->GetParameter(8)*2.);
     }
     for(int iSet=0; iSet<10;iSet++){
         func3peak->ReleaseParameter(iSet);
     }
     //h->GetXaxis()->SetRangeUser(-0.5,15);

     for(int iSet=0; iSet<3;iSet++){
         funcped         ->SetParameter(iSet,   func3peak->GetParameter(iSet));
         func1pe         ->SetParameter(iSet,   func3peak->GetParameter(iSet+3));
         func2pe         ->SetParameter(iSet,   func3peak->GetParameter(iSet+6));
         func3peak_bs    ->SetParameter(iSet,   func3peak->GetParameter(iSet));
         func3peak_bs    ->SetParameter(iSet+3, func3peak->GetParameter(iSet+3));
         func3peak_bs    ->SetParameter(iSet+6, func3peak->GetParameter(iSet+6));
         func3peak_1peall->SetParameter(iSet,   func3peak->GetParameter(iSet));
         func3peak_1peall->SetParameter(iSet+3, func3peak->GetParameter(iSet+3));
         func3peak_1peall->SetParameter(iSet+6, func3peak->GetParameter(iSet+6));
         func3peak_2bs   ->SetParameter(iSet,   func3peak->GetParameter(iSet));
         func3peak_2bs   ->SetParameter(iSet+3, func3peak->GetParameter(iSet+3));
         func3peak_2bs   ->SetParameter(iSet+6, func3peak->GetParameter(iSet+6));
     }
     func3peak_2bs  ->SetParameter(9, func3peak->GetParameter(9));
     func3peak_bs   ->SetParameter(9, func3peak->GetParameter(9));
     func3peak_1peall->SetParameter(9, func3peak->GetParameter(9));
     //func1peak->SetParameter(6, func3peak->GetParameter(9));
     func1pe->SetParameter(6, func3peak->GetParameter(9));
  }

  //func3peak->SetRange(h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());
  //funcped ->SetRange(h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());
  //func1pe ->SetRange(h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());
  //func1peall->SetRange(h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());
  //funcbs  ->SetRange(h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());

  std::cout << " Integral " << h->GetXaxis()->GetXmin() << " - " <<  h->GetXaxis()->GetXmax() << std::endl;
  //Double_t areaped    = funcped   ->Integral(h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());
  Double_t areaped    = funcped   ->Integral(funcped->GetParameter(1)-funcped->GetParameter(2)*5, funcped->GetParameter(1)+funcped->GetParameter(2)*5);
  //std::cout << " Integral " << funcped->GetParameter(0) << " - " <<  funcped->GetParameter(1) << " - " <<  funcped->GetParameter(2) << std::endl;
  Double_t area1pe    = func1pe   ->Integral(h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());
  Double_t area1peall = func1peall->Integral(h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());
  Double_t areabs     = funcbs    ->Integral(h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());
  //Double_t area1pe    = func1pe   ->Integral(funcped->GetParameter(1)-funcped->GetParameter(4)*4, h->GetXaxis()->GetXmax());
  //Double_t area1peall = func1peall->Integral(funcped->GetParameter(1)-funcped->GetParameter(4)*4, h->GetXaxis()->GetXmax());
  //Double_t areabs     = funcbs    ->Integral(funcped->GetParameter(1)-funcped->GetParameter(4)*4, h->GetXaxis()->GetXmax());

  // Result parameters
  Double_t mean      = func1peall->Mean(h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());
  Double_t variance  = func1peall->Variance(h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());
  Double_t peakx = func1peall->GetMaximumX();
  Double_t peaky = func1peall->GetMaximum();
  Double_t FWHMlow  = peakx - func1peall->GetX(peaky*0.5, h->GetXaxis()->GetXmin(), peakx);
  Double_t FWHMhigh = func1peall->GetX(peaky*0.5, peakx, h->GetXaxis()->GetXmax()) - peakx;
  Double_t FWHM = FWHMhigh + FWHMlow;
  Double_t sigmalow  = FWHMlow/(2.*TMath::Sqrt(TMath::Log(2)*2));
  Double_t sigmahigh = FWHMhigh/(2.*TMath::Sqrt(TMath::Log(2)*2));
  Double_t sigma = (sigmahigh + sigmalow);
  Double_t peak1pex     = func1pe->GetParameter(1);
  Double_t peak1pesigma = func1pe->GetParameter(2);
  Double_t peak1pexerr     = func1pe->GetParError(1);
  Double_t peak1pesigmaerr = func1pe->GetParError(2);
  Double_t valley = func->GetMinimum(func->GetParameter(1),func->GetParameter(4));
  Double_t valleyx = func->GetMinimumX(func->GetParameter(1),func->GetParameter(4));
  Double_t peak = func->GetMaximum(func->GetParameter(4)-func->GetParameter(5),func->GetParameter(4)+func->GetParameter(5));
  Double_t pv   = (valley!=0?peak/valley:-1);
  Double_t areasigmalow;
  Double_t areasigmahigh;
  Double_t areasigma;
  Double_t gainpeakx        = peakx      /1.60217662e-19*1.e-12;
  Double_t gainpeak1pex     = peak1pex   /1.60217662e-19*1.e-12;
  Double_t gainpeak1pexerr  = peak1pexerr/1.60217662e-19*1.e-12;
  Double_t chi2 = func->GetChisquare();
  Double_t ndf  = (Double_t)func->GetNDF();

  peakval = peak;

     std::cout << std::endl << std::endl;
     std::cout << "=========== RESULT ==========" << std::endl;
     //std::cout << " mean               = " << " " <<  mean          << " (Mean/Peak = " <<  mean/peakx*100.  << " [%])"    << std::endl;
     //std::cout << " variance           = " << " " <<  variance      << " "    << std::endl;
     //std::cout << " var/mean           = " << " " <<  mean/variance*100. << " [%]"    << std::endl;
     //std::cout << std::endl;
     std::cout << " peakx              = " << " " <<  peakx         << " (Peak/Mean = " <<  peakx/mean*100. << " [%])"  << std::endl;
     std::cout << " peaky              = " << " " <<  peaky         << " "    << std::endl;
     std::cout << " FWHM               = " << " " <<  FWHM          << " " << " (" <<  FWHM    /peakx*100.  << " [%])"  << std::endl;
     std::cout << " FWHMlow            = " << " " <<  FWHMlow       << " " << " (" <<  FWHMlow /peakx*100.  << " [%])"  << std::endl;
     std::cout << " FWHMhigh           = " << " " <<  FWHMhigh      << " " << " (" <<  FWHMhigh/peakx*100.  << " [%])"  << std::endl;
     std::cout << " sigma              = " << " " <<  sigma         << " " << " (" <<  sigma    /peakx*100. << " [%])" << std::endl;
     std::cout << " sigmalow           = " << " " <<  sigmalow      << " " << " (" <<  sigmalow /peakx*100. << " [%])" << std::endl;
     std::cout << " sigmahigh          = " << " " <<  sigmahigh     << " " << " (" <<  sigmahigh/peakx*100. << " [%])" << std::endl;
     std::cout << std::endl;
     std::cout << " 1pe peak in Gauss  = " << " " <<  peak1pex      << " +/-" <<   peak1pexerr      << " "    << std::endl;
     std::cout << " 1pe sigma in Gauss = " << " " <<  peak1pesigma  << " +/-" << peak1pesigmaerr << " (" <<  ((peak1pex!=0)?peak1pesigma/peak1pex*100.:0)  << " [%])" << std::endl;
     std::cout << std::endl;
     std::cout << " Valley Height      = " << " " <<  pv            << " "    << std::endl;
     std::cout << " 1 Peak Height      = " << " " <<  peak          << " "    << std::endl;
     std::cout << " P/V ratio by func  = " << " " <<  (valley!=0?peak/valley:-1)   << " "    << std::endl;
     std::cout << " P/V ratio by histo = " << " " <<  (h->GetBinContent(h->FindBin(valleyx))!=0?h->GetBinContent(h->FindBin(peak1pex))/h->GetBinContent(h->FindBin(valleyx)):-1)   << " "    << std::endl;
     std::cout << std::endl;
     std::cout << " Area of pedestal   = " << " " <<  areaped      << " "    << std::endl;
     std::cout << " Area of 1pe peak   = " << " " <<  area1pe      << " "    << std::endl;
     std::cout << " Area of BS         = " << " " <<  areabs       << " "    << std::endl;
     std::cout << " Area of 1pe total  = " << " " <<  area1peall   << " "    << std::endl;
     std::cout << " Occupancy of 1pe   = " << " " <<  ((area1peall+areaped)!=0?area1peall/(area1peall+areaped)*100.:-1)   << " [%]"    << std::endl;
     if (areaped>0) std::cout << " 1PE Occupancy      = " << " " <<  (areaped!=0?area1pe/areaped*100.:-1) << " [%]"    << std::endl;
     std::cout << " Gain               = " << " " <<  peakx/1.60217662e-19*1.e-12 << " "    << std::endl;
     std::cout << "=============================" << std::endl;

     std::cout << std::endl << "=============================" << std::endl;
     std::cout << "AAAAA  Gain = " << " " <<  peakx/1.60217662e-19*1.e-12 << " "    << std::endl;
     std::cout << "AAAAA  P/V  = " << " " <<  pv                          << " "    << std::endl;
     std::cout << "AAAAA  FWHM = " << " " <<  FWHM          << " " << " (" <<  (peakx!=0?FWHM/peakx*100.:-1)  << " [%])"  << std::endl;

     std::cout << "=============================" << std::endl;
     std::cout << std::endl;

     gainpeak = peakx/1.60217662e-19*1.e-12;
     respeak = peakx;
     respeakerr = peak1pexerr;
     resfwhm = (peakx!=0?FWHM/peakx*100.:-1);
     ressigma = (peak1pex!=0)?peak1pesigma/peak1pex*100.:0;
     chi2ndf = ndf!=0?chi2/ndf:0;
     bsratio = area1peall>0 ? areabs / area1peall : 0;

//  if (DrawEff) {
    TGraph *geff = (TGraph*)(func1peall->DrawIntegral("al"));
    geff->SetName(Form("eff1peall%s",name));
    TGraph *geffinv = new TGraph();
    geffinv->SetName(Form("effinv1peall%s",name));
    Double_t  yoffset = 1.-0.93;
    if (kTRUE) {//Normalization
      for (Int_t ig = 0; ig < geff->GetN(); ig++) {
         geff->GetY()[ig] = geff->GetY()[ig] / geff->GetY()[geff->GetN()-1]*100.;
         geffinv->SetPoint(ig, geff->GetY()[ig], geff->GetX()[ig]);
      }
      if (yoffset != 0) {
      for (Int_t ig = 0; ig < geff->GetN(); ig++) {
         geff->GetY()[ig] = (geff->GetY()[ig] + yoffset) / (geff->GetY()[geff->GetN()-1] + yoffset)*100.;
         //std::cout << ig << " " << geff->GetY()[ig] << std::endl;
         geffinv->SetPoint(ig, geff->GetY()[ig], geff->GetX()[ig]);
      }
      }
    }
    geff->GetXaxis()->SetTitle(h->GetXaxis()->GetTitle());
    geff->GetYaxis()->SetTitle("Efficiency [%]");
//  }

    Double_t eff_center = geff->Eval(peakx);
     areasigmalow  = peakx - geffinv->Eval(eff_center - 34.1);
     areasigmahigh = geffinv->Eval(eff_center + 34.1) - peakx;
     areasigma = (areasigmalow+areasigmahigh)*0.5;
     std::cout << "========== Sigma in area ===================" << std::endl;
     std::cout << " Gaussian sigma = " <<  peak1pesigma      << " " << " (" <<  peak1pesigma/peak1pex*100. << " [%])" << std::endl;
     std::cout << " sigma        = " << " " <<  areasigma         << " " << " (" <<  areasigma    /peakx*100. << " [%])" << std::endl;
     std::cout << " sigmalow     = " << " " <<  areasigmalow      << " " << " (" <<  areasigmalow /peakx*100. << " [%])" << std::endl;
     std::cout << " sigmahigh    = " << " " <<  areasigmahigh     << " " << " (" <<  areasigmahigh/peakx*100. << " [%])" << std::endl;
     std::cout << "=============================" << std::endl;

  //for(int iSet=0; iSet<nPeaks*3+1;iSet++){
  //   std::cout << func->GetParameter(iSet) << " ";
  //}
  //std::cout << std::endl;

     result.mean            = mean           ;
     result.variance        = variance       ;
     result.peakx           = peakx          ;
     result.peaky           = peaky          ;
     result.FWHMlow         = FWHMlow        ;
     result.FWHMhigh        = FWHMhigh       ;
     result.FWHM            = FWHM           ;
     result.sigmalow        = sigmalow       ;
     result.sigmahigh       = sigmahigh      ;
     result.sigma           = sigma          ;
     result.peak1pex        = peak1pex       ;
     result.peak1pesigma    = peak1pesigma   ;
     result.peak1pexerr     = peak1pexerr    ;
     result.peak1pesigmaerr = peak1pesigmaerr;
     result.valley          = valley         ;
     result.valleyx         = valleyx        ;
     result.peak            = peak           ;
     result.pv              = pv             ;
     result.areasigmalow    = areasigmalow   ;
     result.areasigmahigh   = areasigmahigh  ;
     result.areasigma       = areasigma      ;
     result.gainpeakx       = gainpeakx      ;
     result.gainpeak1pex    = gainpeak1pex   ;
     result.gainpeak1pexerr = gainpeak1pexerr;
     result.chi2            = chi2           ;
     result.ndf             = ndf            ;


  if (!DrawEff) {
     gPad->SetLogy(1);
     //h->SetMaximum(peaky*1.4);
     h->Draw("e1");
     funcped ->Draw("same");
     func1pe ->Draw("same");
     if (nPeaks>2) { func2pe->Draw("same"); }
     funcbs  ->Draw("same");
     if (nPeaks>2) { func3peak_bs->Draw("same"); }
     func1peall->Draw("same");
     //func3peak->Draw("same");
     func->Draw("same");

     gPad->Update();

     gStyle->SetTitleXOffset(0.);
     gStyle->SetTitleX(0.);
     gStyle->SetTitleAlign(13);//left a=10*h+v, h=1,2,3 for left, center, top,  v=1,2,3 for bottom, middle, top

     /// Stat Box ///
     Float_t st0_lower_left_x = 0.6;
     Float_t st0_lower_left_y = 0.35;
     Float_t st_Width = 0.4;
     Float_t st_Height = 0.65;

     /* //     TPaveStats *st0 = (TPaveStats*)gPad->GetPrimitive("stats"); */
     /* //TPaveStats *st0 = (TPaveStats*)h->FindObject("stats"); */
     /* TPaveStats *st0 = new TPaveStats(); */
     /* st0->SetName("result"); */
     /* TList *listOfLines = st0->GetListOfLines(); */
     /* if (st0) { */
     /*   /\* st0->SetX1NDC(st0_lower_left_x); *\/ */
     /*   /\* st0->SetY1NDC(st0_lower_left_y); *\/ */
     /*   /\* st0->SetX2NDC(st0_lower_left_x + st_Width); *\/ */
     /*   /\* st0->SetY2NDC(st0_lower_left_y + st_Height); *\/ */
     /*   TLatex *myt[9]; */
     /*   /\* myt[0] = new TLatex(0,0,Form("Single PE RESULT")); *\/ */
     /*   /\* myt[1] = new TLatex(0,0,Form("Gain  = %1.2e",result.gainpeakx)); *\/ */
     /*   /\* myt[2] = new TLatex(0,0,Form("Peak  = %f ",result.peakx)); *\/ */
     /*   /\* myt[3] = new TLatex(0,0,Form("Sigma = %2.1f [%%]",((result.peakx!=0)?result.sigma/result.peakx*100.:0))); *\/ */
     /*   /\* myt[4] = new TLatex(0,0,Form("P/V ratio = %2.2f",result.pv)); *\/ */
     /*   /\* myt[5] = new TLatex(0,0,Form("FWHM = %3.1f [%%]", (result.peakx!=0?result.FWHM/result.peakx*100.:-1))); *\/ */
     /*   //myt[6] = new TLatex(0,0,Form("Peak (1pe only) = %2.2f [pC]",    result.peak1pex)); */
     /*   //myt[7] = new TLatex(0,0,Form("Sigma (1pe only) = %2.1f [%%]",     result.peak1pesigma)); */
     /*   //myt[8] = new TLatex(0,0,Form("Gain (1pe only) = %1.3e",  result.gainpeak1pex)); */
     /*   for (Int_t i = 0; i < 6; i++) { */
     /* 	 //	 myt[i] ->SetTextFont(42); */
     /* 	 //	 myt[i] ->SetTextSize(0.03); */
     /* 	 //          myt[i] ->SetTextColor(kBlue); */
     /* 	 listOfLines->Add(myt[i]); */
     /*   } */
     /*   h->SetStats(0); */

     /*   gPad->Modified(); */
     /* } */


  } else if (DrawEff) {
     geff->Draw("AL");
     TLine *lmean = new TLine(mean, 0, mean, geff->Eval(mean));
     TLine *lpeak = new TLine(peakx, 0, peakx, geff->Eval(peakx));
     lmean->SetLineStyle(2);
     lpeak->SetLineStyle(1);
     lmean->SetLineColor(kRed);
     lpeak->SetLineColor(kRed);
     lmean->Draw("same");
     lpeak->Draw("same");
  }

}

//#ifdef __CINT__
#ifndef rcomp
void simplefit(const char* rootfilename = "input.root", const char* outfilename = "test.root", bool Update = kFALSE, Int_t npeak = 3) {
#else
int main(int argc, char* argv[]){
   if (argc < 3) {
      //std::cout << "USAGE : ./simplefit (Input Root filename) (Output Root filename) (optional: num of peak finding) (optional, 1:Add in existing root, 0:Recreate new (default) )" << std::endl;
      std::cout << "USAGE : ./simplefit (Input Root filename) (Output Text filename) (optional: num of peak finding) (optional, 1:Add in existing root, 0:Recreate new (default) )" << std::endl;
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

   bool OutputTxt = kTRUE;
   TString pdffilename = Form("%s/%s.pdf",outputpdfpath,rootfilename);
   TString textfilename = Form("%s/%s.txt", outputtxtpath, rootfilename);
   gErrorIgnoreLevel = 5000;
   
   TCanvas *c = new TCanvas("c","c",800,600);
   
   TFile* file = new TFile(rootfilename, "READ");

   TH1 *h = (TH1*)0;
   
   //c->Print(Form("%s[",pdffilename.Data()));

   //TFile *rootf;
   //if (!Update) {
   //  rootf = TFile::Open(outfilename, "RECREATE", outfilename);
   //   rootf->Close();
   // }
   std::ofstream resultfile;
   if (OutputTxt) 
     //resultfile.open(textfilename.Data());
     resultfile.open(outfilename);
   //   for (int ich = 0; ich < 1; ich++) {
   for (int ich : {0,1,2,3}) {
   //for (int ich : {1}) {
     
      h = (TH1*)(file->Get(Form("all_ch%02d_hgain",ich)));
   // ここから変更箇所
   if (h == nullptr || h->GetMaximum() <= 10) {
   std::cout << "Channel " << ich << " has too small peak height (" << (h ? h->GetMaximum() : 0) << "). Skipping." << std::endl;
   continue;
   }
   // ここまで変更箇所
      std::cout << Form("all_ch%02d_hgain",ich) << std::endl;
      std::cout << h << std::endl;
      std::cout << "ENTRIES " << h->GetEntries() << " " << h->GetTitle() << std::endl;
      
      h->SetTitle(Form("%s", rootfilename));
      //      h->GetXaxis()->SetTitle("Charge [pC]");
      h->GetXaxis()->SetTitle("High Gain [ADC]");
      h->GetXaxis()->SetRangeUser(400,800);
      //h->GetXaxis()->SetRangeUser(-0.5,15);
      
      
      //c->SetLogy(1);
      
      c->cd();
      /////////////////////////////////////////////////////
      FitPedPe(h, Form("BLPMT"), npeak/*Peak*/, 0, kFALSE);
      /////////////////////////////////////////////////////
      h->Draw();
      std::cout << "BBB RESULT : " << gainpeak << ", " << respeak << ", " << respv << ", " << resfwhm << ", " << ressigma << ", " << bsratio << ", " << chi2ndf << ", " << std::endl;
      //(respeak!=0?TMath::Sqrt(resfwhm*resfwhm*respeak*respeak*1e-4-adcpedestalsigma[ich]*50.*1.e-3*0.1*adcpedestalsigma[ich]*50.*1.e-3*0.1)/respeak*100.:-1) << std::endl;
      if (OutputTxt) 
	resultfile << ich << ", " << "BBB RESULT : " << gainpeak << ", " << respeak << ", " << respv << ", " << resfwhm << ", " << ressigma << ", " << bsratio << ", " << chi2ndf << ", " << std::endl;
	//resultfile << "BBB RESULT : " << gainpeak << ", " << respeak << ", " << respv << ", " << resfwhm << ", " << ressigma << ", " << bsratio << ", " << chi2ndf << ", " << std::endl;
   //(respeak!=0?TMath::Sqrt(resfwhm*resfwhm*respeak*respeak*1e-4-adcpedestalsigma[ich]*50.*1.e-3*0.1*adcpedestalsigma[ich]*50.*1.e-3*0.1*2.355*2.355)/respeak*100.:-1) << std::endl;
      result.ch = ich;

////////////////////////////////////////////////
      // rootf = TFile::Open(outfilename, "UPDATE", outfilename);
      // adcresult_fill(rootf, result);
      // rootf->Close();
////////////////////////////////////////////////

      gStyle->SetOptFit(1111);
      gStyle->SetOptTitle(1);
      
      //c->SaveAs(Form("%srun%06d_1pechargelog_p%d_%s.pdf",path,run,npeak,serial.Data()));
      //c->Print(pdffilename);
      
      c->SetLogy(0);
      if (peakval < 10) peakval = 10;
      else if (peakval > 1e4) peakval = 1e4;
      h->GetYaxis()->SetRangeUser(0, peakval*1.6);
      //c->SaveAs(Form("%srun%06d_1pecharge_p%d_%s.root",path,run,npeak,serial.Data()));
      //c->SaveAs(Form("%srun%06d_1pecharge_p%d_%s.pdf",path,run,npeak,serial.Data()));
      
      //c->Print(pdffilename);
   
   }
   
   //c->Print(Form("%s]",pdffilename.Data()));
   
   file->Close();
   if (OutputTxt) resultfile.close();

   return 0;
   
}
