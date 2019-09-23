//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Nov  9 01:59:17 2016 by ROOT version 5.34/34
// from TTree T/USBWC Data Tree
// found on file: Run_1_Data_9_15_2016_Binary.bin.root
//////////////////////////////////////////////////////////

#ifndef run1vis_h
#define run1vis_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class run1vis {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           EventNumber_usbwc;
   Float_t         SamplingPeriod_usbwc;
   Double_t        EpochTime_usbwc;
   Int_t           Year_usbwc;
   Int_t           Month_usbwc;
   Int_t           Day_usbwc;
   Int_t           Hour_usbwc;
   Int_t           Minute_usbwc;
   Int_t           Second_usbwc;
   Int_t           Millisecond_usbwc;
   Int_t           channel_usbwc;
   Int_t           EventIDsamIndex_usbwc;
   Int_t           FirstCellToPlotsamIndex_usbwc;
   Float_t         MeasuredBaseline_usbwc[8];
   Float_t         AmplitudeValue_usbwc[8];
   Float_t         ComputedCharge_usbwc[8];
   Float_t         RiseTimeInstant_usbwc[8];
   Float_t         FallTimeInstant_usbwc[8];
   Float_t         RawTriggerRate_usbwc[8];
   Short_t         amplValues_usbwc[8][1024];

   // List of branches
   TBranch        *b_EventNumber_usbwc;   //!
   TBranch        *b_SamplingPeriod_usbwc;   //!
   TBranch        *b_EpochTime_usbwc;   //!
   TBranch        *b_Year_usbwc;   //!
   TBranch        *b_Month_usbwc;   //!
   TBranch        *b_Day_usbwc;   //!
   TBranch        *b_Hour_usbwc;   //!
   TBranch        *b_Minute_usbwc;   //!
   TBranch        *b_Second_usbwc;   //!
   TBranch        *b_Millisecond_usbwc;   //!
   TBranch        *b_channel_usbwc;   //!
   TBranch        *b_EventIDsamIndex_usbwc;   //!
   TBranch        *b_FirstCellToPlotsamIndex_usbwc;   //!
   TBranch        *b_MeasuredBaseline_usbwc;   //!
   TBranch        *b_AmplitudeValue_usbwc;   //!
   TBranch        *b_ComputedCharge_usbwc;   //!
   TBranch        *b_RiseTimeInstant_usbwc;   //!
   TBranch        *b_RawTriggerRate_usbwc;   //!
   TBranch        *b_RawTriggerRate_usbwc;   //!
   TBranch        *b_amplValues_usbwc;   //!

   run1vis(TTree *tree=0);
   void DrawOneChannel_waveformSimple(Int_t evID, Int_t chID, TString key_lineOn, TString key_printInfoTrue);
   virtual ~run1vis();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef run1vis_cxx
run1vis::run1vis(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("./wfData/Run_1_Data_9_15_2016_Binary.bin.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("./wfData/Run_1_Data_9_15_2016_Binary.bin.root");
      }
      f->GetObject("T",tree);

   }
   Init(tree);
}

run1vis::~run1vis()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t run1vis::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t run1vis::LoadTree(Long64_t entry)
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

void run1vis::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("EventNumber_usbwc", &EventNumber_usbwc, &b_EventNumber_usbwc);
   fChain->SetBranchAddress("SamplingPeriod_usbwc", &SamplingPeriod_usbwc, &b_SamplingPeriod_usbwc);
   fChain->SetBranchAddress("EpochTime_usbwc", &EpochTime_usbwc, &b_EpochTime_usbwc);
   fChain->SetBranchAddress("Year_usbwc", &Year_usbwc, &b_Year_usbwc);
   fChain->SetBranchAddress("Month_usbwc", &Month_usbwc, &b_Month_usbwc);
   fChain->SetBranchAddress("Day_usbwc", &Day_usbwc, &b_Day_usbwc);
   fChain->SetBranchAddress("Hour_usbwc", &Hour_usbwc, &b_Hour_usbwc);
   fChain->SetBranchAddress("Minute_usbwc", &Minute_usbwc, &b_Minute_usbwc);
   fChain->SetBranchAddress("Second_usbwc", &Second_usbwc, &b_Second_usbwc);
   fChain->SetBranchAddress("Millisecond_usbwc", &Millisecond_usbwc, &b_Millisecond_usbwc);
   fChain->SetBranchAddress("channel_usbwc", &channel_usbwc, &b_channel_usbwc);
   fChain->SetBranchAddress("EventIDsamIndex_usbwc", &EventIDsamIndex_usbwc, &b_EventIDsamIndex_usbwc);
   fChain->SetBranchAddress("FirstCellToPlotsamIndex_usbwc", &FirstCellToPlotsamIndex_usbwc, &b_FirstCellToPlotsamIndex_usbwc);
   fChain->SetBranchAddress("MeasuredBaseline_usbwc", MeasuredBaseline_usbwc, &b_MeasuredBaseline_usbwc);
   fChain->SetBranchAddress("AmplitudeValue_usbwc", AmplitudeValue_usbwc, &b_AmplitudeValue_usbwc);
   fChain->SetBranchAddress("ComputedCharge_usbwc", ComputedCharge_usbwc, &b_ComputedCharge_usbwc);
   fChain->SetBranchAddress("RiseTimeInstant_usbwc", RiseTimeInstant_usbwc, &b_RiseTimeInstant_usbwc);
   fChain->SetBranchAddress("FallTimeInstant_usbwc", FallTimeInstant_usbwc, &b_RawTriggerRate_usbwc);
   fChain->SetBranchAddress("RawTriggerRate_usbwc", RawTriggerRate_usbwc, &b_RawTriggerRate_usbwc);
   fChain->SetBranchAddress("amplValues_usbwc", amplValues_usbwc, &b_amplValues_usbwc);
   Notify();
}

Bool_t run1vis::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void run1vis::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t run1vis::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef run1vis_cxx
