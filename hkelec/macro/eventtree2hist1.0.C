/*id: eventtree2hist.C*/
/*Place: hkpd@hkdaq:~/hkelec/analysis_macro/ */
/*Last Edit: 2025-10-10 Daisuke Kishita*/
/*eventtree.root to triggered data*/
/*コンパイル可能*/

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TString.h>
#include <iostream>
#include <vector>

// ⚠️ 外部のヘッダファイルをインクルード（パスはMakefileで指定します）
#include "RootInterface/Hit.h"
#include "RootInterface/MetaData.h" 
// (注: ヘッダファイルが RootInterface/ のサブディレクトリにあるため、このように記述するのが適切です)

/* Template macro for reading Event Tree
 * Tutorial: https://root.cern/doc/master/hvector_8C.html */

void read_event_tree (TString input_file, TString output_file) {

  /* Input and output files */
  auto ifile = TFile::Open(input_file, "READ");
  auto ofile = TFile::Open(output_file, "RECREATE");

  /* Histograms */
  auto h_hgain = new TH1D("h_hgain", "High Gain ADC", 4096, 0., 4096.);

  /* Load tree */
  auto tree = ifile->Get<TTree>("event");
  std::vector<Hit> *v_hit = nullptr;
  std::vector<Hit> *v_trig = nullptr;
  tree->SetBranchAddress("NormalHits", &v_hit);
  tree->SetBranchAddress("TriggerHits", &v_trig);

  /* Print metadata */
  auto metadata = ifile->Get<MetaData>("metadata");
  metadata->Print();

  /* Process events */
  long nEvents = tree->GetEntries();
  for (long iEvent=0; iEvent<nEvents; ++iEvent) {
    tree->GetEntry(iEvent);
    std::cout << "-----------------------------------------------------------------------------" << std::endl
              << "                          Event: " << iEvent                                   << std::endl
              << "-----------------------------------------------------------------------------" << std::endl;
    for (const Hit &trig : *v_trig) {
      std::cout << "Trigger hit:"
                << " type=" << trig.type
                << " channel=" << trig.channel
                << " hgain="   << trig.hgain
                << " lgain="   << trig.lgain
                << " tdc="     << trig.tdc
                << " time="    << trig.time
                << " tot="     << trig.tot
                << std::endl;

    }
    for (const Hit &hit : *v_hit) {
      std::cout << "Normal hit:"
                << " type=" << hit.type
                << " channel=" << hit.channel
                << " hgain="   << hit.hgain
                << " lgain="   << hit.lgain
                << " tdc="     << hit.tdc
                << " time="    << hit.time
                << " tot="     << hit.tot
                << std::endl;
      h_hgain->Fill(hit.hgain);
    }
    std::cout << std::endl;

  }

  ofile->Write();
  ofile->Close();
  ifile->Close();

}

// スタンドアロン実行のためのメイン関数
int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <input_file.root> <output_file.root>" << std::endl;
        return 1;
    }
    
    // 既存関数を呼び出す
    read_event_tree(argv[1], argv[2]); 
    
    return 0;
}


