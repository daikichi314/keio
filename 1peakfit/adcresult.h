//#include "TRoot.h"
#include <string.h>

class FitResultSPE
{
   public:
//Setting
  Int_t    run;
  Int_t    ch;
  //std::string   serial;
  char serial[16];
  Float_t  HVS;
  Double_t unixtime;
  Int_t    LDA;
  Int_t    LDW;

//Result
  Double_t mean      ;
  Double_t variance  ;
  Double_t peakx ;
  Double_t peaky ;
  Double_t FWHMlow  ;
  Double_t FWHMhigh ;
  Double_t FWHM ;
  Double_t sigmalow  ;
  Double_t sigmahigh ;
  Double_t sigma ;
  Double_t peak1pex     ;
  Double_t peak1pesigma ;
  Double_t peak1pexerr     ;
  Double_t peak1pesigmaerr ;
  Double_t valley ;
  Double_t valleyx ;
  Double_t peak ;
  Double_t pv ;
  Double_t areasigmalow;
  Double_t areasigmahigh;
  Double_t areasigma;
  Double_t gainpeakx      ;
  Double_t gainpeak1pex   ;
  Double_t gainpeak1pexerr;

  Double_t chi2;
  Double_t ndf;

};


int adcresult_fill(TFile *rootf, FitResultSPE res)
{
    const char *adctrname = "adcr";
    //std::vector<std::string> serialstrings;
    //    serialstrings.push_back(str);

    TTree *adctr = 0;
    TBranch        *b_run            ;
    TBranch        *b_ch             ;
    //TBranch        *b_nserial        ;
    TBranch        *b_serial         ;
    TBranch        *b_HVS            ;
    TBranch        *b_unixtime       ;
    TBranch        *b_LDA            ;
    TBranch        *b_LDW            ;
    TBranch        *b_mean           ;
    TBranch        *b_variance       ;
    TBranch        *b_peakx          ;
    TBranch        *b_peaky          ;
    TBranch        *b_FWHMlow        ;
    TBranch        *b_FWHMhigh       ;
    TBranch        *b_FWHM           ;
    TBranch        *b_sigmalow       ;
    TBranch        *b_sigmahigh      ;
    TBranch        *b_sigma          ;
    TBranch        *b_peak1pex       ;
    TBranch        *b_peak1pesigma   ;
    TBranch        *b_peak1pexerr    ;
    TBranch        *b_peak1pesigmaerr;
    TBranch        *b_valley         ;
    TBranch        *b_valleyx        ;
    TBranch        *b_peak           ;
    TBranch        *b_pv             ;
    TBranch        *b_areasigmalow   ;
    TBranch        *b_areasigmahigh  ;
    TBranch        *b_areasigma      ;
    TBranch        *b_gainpeakx      ;
    TBranch        *b_gainpeak1pex   ;
    TBranch        *b_gainpeak1pexerr;
    TBranch        *b_chi2;
    TBranch        *b_ndf;
    //TBranch        *b_length;
    //TBranch        *b_description;
    //TBranch        *b_start;
    //TBranch        *b_serial;
    //TBranch        *b_nhv;
    //TBranch        *b_hv;
    //TBranch        *b_hvi;
    //TBranch        *b_hvs;
    //TBranch        *b_nthr;
    //TBranch        *b_thr;
    //TBranch        *b_lsfactor;
    //TBranch        *b_lsamp;
    //TBranch        *b_lswidth;
    //TBranch        *b_nadc;
    //TBranch        *b_gaincorr;
    //TBranch        *b_pedestal;
    //TBranch        *b_ntdc;
    //TBranch        *b_delay;

    if (rootf->Get(adctrname)) {
       adctr = (TTree*)(rootf->Get(adctrname));
       adctr->SetBranchAddress("run"            , &(res.run            ),   &b_run            );
       adctr->SetBranchAddress("ch"             , &(res.ch             ),   &b_ch             );
       adctr->SetBranchAddress("serial"         ,   res.serial          ,   &b_serial         );
       adctr->SetBranchAddress("HVS"            , &(res.HVS            ),   &b_HVS            );
       adctr->SetBranchAddress("unixtime"       , &(res.unixtime       ),   &b_unixtime       );
       adctr->SetBranchAddress("LDA"            , &(res.LDA            ),   &b_LDA            );
       adctr->SetBranchAddress("LDW"            , &(res.LDW            ),   &b_LDW            );
       adctr->SetBranchAddress("mean"           , &(res.mean           ),   &b_mean           );
       adctr->SetBranchAddress("variance"       , &(res.variance       ),   &b_variance       );
       adctr->SetBranchAddress("peakx"          , &(res.peakx          ),   &b_peakx          );
       adctr->SetBranchAddress("peaky"          , &(res.peaky          ),   &b_peaky          );
       adctr->SetBranchAddress("FWHMlow"        , &(res.FWHMlow        ),   &b_FWHMlow        );
       adctr->SetBranchAddress("FWHMhigh"       , &(res.FWHMhigh       ),   &b_FWHMhigh       );
       adctr->SetBranchAddress("FWHM"           , &(res.FWHM           ),   &b_FWHM           );
       adctr->SetBranchAddress("sigmalow"       , &(res.sigmalow       ),   &b_sigmalow       );
       adctr->SetBranchAddress("sigmahigh"      , &(res.sigmahigh      ),   &b_sigmahigh      );
       adctr->SetBranchAddress("sigma"          , &(res.sigma          ),   &b_sigma          );
       adctr->SetBranchAddress("peak1pex"       , &(res.peak1pex       ),   &b_peak1pex       );
       adctr->SetBranchAddress("peak1pesigma"   , &(res.peak1pesigma   ),   &b_peak1pesigma   );
       adctr->SetBranchAddress("peak1pexerr"    , &(res.peak1pexerr    ),   &b_peak1pexerr    );
       adctr->SetBranchAddress("peak1pesigmaerr", &(res.peak1pesigmaerr),   &b_peak1pesigmaerr);
       adctr->SetBranchAddress("valley"         , &(res.valley         ),   &b_valley         );
       adctr->SetBranchAddress("valleyx"        , &(res.valleyx        ),   &b_valleyx        );
       adctr->SetBranchAddress("peak"           , &(res.peak           ),   &b_peak           );
       adctr->SetBranchAddress("pv"             , &(res.pv             ),   &b_pv             );
       adctr->SetBranchAddress("areasigmalow"   , &(res.areasigmalow   ),   &b_areasigmalow   );
       adctr->SetBranchAddress("areasigmahigh"  , &(res.areasigmahigh  ),   &b_areasigmahigh  );
       adctr->SetBranchAddress("areasigma"      , &(res.areasigma      ),   &b_areasigma      );
       adctr->SetBranchAddress("gainpeakx"      , &(res.gainpeakx      ),   &b_gainpeakx      );
       adctr->SetBranchAddress("gainpeak1pex"   , &(res.gainpeak1pex   ),   &b_gainpeak1pex   );
       adctr->SetBranchAddress("gainpeak1pexerr", &(res.gainpeak1pexerr),   &b_gainpeak1pexerr);
       adctr->SetBranchAddress("chi2"           , &(res.chi2           ),   &b_chi2           );
       adctr->SetBranchAddress("ndf"            , &(res.ndf            ),   &b_ndf            );
    } else {
       adctr = new TTree(adctrname,"ADC SPE Result");
       adctr->Branch("run"            , &(res.run            ),"run/I"            );
       adctr->Branch("ch"             , &(res.ch             ),"ch/I"             );
       adctr->Branch("serial"         , &(res.serial         ),"serial[16]/C"     );
       adctr->Branch("HVS"            , &(res.HVS            ),"HVS/F"            );
       adctr->Branch("unixtime"       , &(res.unixtime       ),"unixtime/D"       );
       adctr->Branch("LDA"            , &(res.LDA            ),"LDA/I"            );
       adctr->Branch("LDW"            , &(res.LDW            ),"LDW/I"            );
       adctr->Branch("mean"           , &(res.mean           ),"mean/D"           );
       adctr->Branch("variance"       , &(res.variance       ),"variance/D"       );
       adctr->Branch("peakx"          , &(res.peakx          ),"peakx/D"          );
       adctr->Branch("peaky"          , &(res.peaky          ),"peaky/D"          );
       adctr->Branch("FWHMlow"        , &(res.FWHMlow        ),"FWHMlow/D"        );
       adctr->Branch("FWHMhigh"       , &(res.FWHMhigh       ),"FWHMhigh/D"       );
       adctr->Branch("FWHM"           , &(res.FWHM           ),"FWHM/D"           );
       adctr->Branch("sigmalow"       , &(res.sigmalow       ),"sigmalow/D"       );
       adctr->Branch("sigmahigh"      , &(res.sigmahigh      ),"sigmahigh/D"      );
       adctr->Branch("sigma"          , &(res.sigma          ),"sigma/D"          );
       adctr->Branch("peak1pex"       , &(res.peak1pex       ),"peak1pex/D"       );
       adctr->Branch("peak1pesigma"   , &(res.peak1pesigma   ),"peak1pesigma/D"   );
       adctr->Branch("peak1pexerr"    , &(res.peak1pexerr    ),"peak1pexerr/D"    );
       adctr->Branch("peak1pesigmaerr", &(res.peak1pesigmaerr),"peak1pesigmaerr/D");
       adctr->Branch("valley"         , &(res.valley         ),"valley/D"         );
       adctr->Branch("valleyx"        , &(res.valleyx        ),"valleyx/D"        );
       adctr->Branch("peak"           , &(res.peak           ),"peak/D"           );
       adctr->Branch("pv"             , &(res.pv             ),"pv/D"             );
       adctr->Branch("areasigmalow"   , &(res.areasigmalow   ),"areasigmalow/D"   );
       adctr->Branch("areasigmahigh"  , &(res.areasigmahigh  ),"areasigmahigh/D"  );
       adctr->Branch("areasigma"      , &(res.areasigma      ),"areasigma/D"      );
       adctr->Branch("gainpeakx"      , &(res.gainpeakx      ),"gainpeakx/D"      );
       adctr->Branch("gainpeak1pex"   , &(res.gainpeak1pex   ),"gainpeak1pex/D"   );
       adctr->Branch("gainpeak1pexerr", &(res.gainpeak1pexerr),"gainpeak1pexerr/D");
       adctr->Branch("chi2"           , &(res.chi2           ),"chi2/D"           );
       adctr->Branch("ndf"            , &(res.ndf            ),"ndf/D"            );
    }
    adctr->Fill();
    adctr->Write(adctrname,TObject::kOverwrite);
    rootf->Close();

    return 0;
}

