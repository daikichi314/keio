//header.h
//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Apr  3 11:48:40 2018 by ROOT version 5.34/36
// from TTree header/Run information
// found on file: /data/run/unified/run000040.root
//////////////////////////////////////////////////////////

#ifndef header_h
#define header_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include <vector>
#include <string>

// Fixed size dimensions of array or collections stored in the TTree if any.

class header {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           run;
   Int_t           length;
   Double_t        start;
   Char_t          description[256];   //[length]
   Int_t           nhv;
   Float_t         HV[16];   //[nhv]
   Float_t         HVI[16];   //[nhv]
   Float_t         HVS[16];   //[nhv]
   Int_t           nthr;
   Float_t         THR[8];   //[nthr]
   Float_t         LSfactor;
   Float_t         LSamp;
   Float_t         LSwidth;
   Int_t           nadc;
   Float_t         gaincorr[8];   //[nadc]
   Float_t         pedestal[8];   //[nadc]
   Int_t           ntdc;
   Float_t         delay[16];   //[ntdc]
   ////vector<string>  *serial;
   std::vector<std::string>  *serial;

   //std::string          *serial;
   //std::string          serial;
   //std::vector<std::string> **serial;
   //std::vector<std::string> *serial;
   //std::vector<std::string> serial;

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_length;   //!
   TBranch        *b_start;   //!
   TBranch        *b_description;   //!
   TBranch        *b_nhv;   //!
   TBranch        *b_HV;   //!
   TBranch        *b_HVI;   //!
   TBranch        *b_HVS;   //!
   TBranch        *b_nthr;   //!
   TBranch        *b_THR;   //!
   TBranch        *b_LSfactor;   //!
   TBranch        *b_LSamp;   //!
   TBranch        *b_LSwidth;   //!
   TBranch        *b_nadc;   //!
   TBranch        *b_gaincorr;   //!
   TBranch        *b_pedestal;   //!
   TBranch        *b_ntdc;   //!
   TBranch        *b_delay;   //!
   TBranch        *b_serial;   //!

   header(TTree *tree=0);
   virtual ~header();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   //virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef header_cxx
header::header(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/data/run/unified/run000040.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/data/run/unified/run000040.root");
      }
      f->GetObject("header",tree);

   }
   Init(tree);
}

header::~header()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t header::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t header::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void header::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   serial = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   //fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("length", &length, &b_length);
   fChain->SetBranchAddress("start", &start, &b_start);
   fChain->SetBranchAddress("description", description, &b_description);
   fChain->SetBranchAddress("nhv", &nhv, &b_nhv);
   fChain->SetBranchAddress("HV", HV, &b_HV);
   fChain->SetBranchAddress("HVI", HVI, &b_HVI);
   fChain->SetBranchAddress("HVS", HVS, &b_HVS);
   fChain->SetBranchAddress("nthr", &nthr, &b_nthr);
   fChain->SetBranchAddress("THR", THR, &b_THR);
   fChain->SetBranchAddress("LSfactor", &LSfactor, &b_LSfactor);
   fChain->SetBranchAddress("LSamp", &LSamp, &b_LSamp);
   fChain->SetBranchAddress("LSwidth", &LSwidth, &b_LSwidth);
   fChain->SetBranchAddress("nadc", &nadc, &b_nadc);
   fChain->SetBranchAddress("gaincorr", gaincorr, &b_gaincorr);
   fChain->SetBranchAddress("pedestal", pedestal, &b_pedestal);
   fChain->SetBranchAddress("ntdc", &ntdc, &b_ntdc);
   fChain->SetBranchAddress("delay", delay, &b_delay);
   fChain->SetBranchAddress("serial", &serial, &b_serial);
   Notify();
}

Bool_t header::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void header::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t header::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef header_cxx
