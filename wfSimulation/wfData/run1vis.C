#define run1vis_cxx
#include "run1vis.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "../src/waveformSimple.hh";

void run1vis::DrawOneChannel_waveformSimple(Int_t evID, Int_t chID, TString key_lineOn, TString key_printInfoTrue){
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  Long64_t jentry = evID;
  Long64_t ientry = LoadTree(jentry);
  nb = fChain->GetEntry(jentry);   nbytes += nb;
  Int_t i = 0;
  const Int_t np = 1024;
  const Int_t nCh = 8;
  Int_t npBaseLine = 20;
  Double_t tLeft  = 5;
  Double_t tRight = 10;
  Double_t vVal = 0.06;
  Int_t nSplinePoints = 0;
  //Int_t nSplinePoints = 1;
  Double_t* wfT = new Double_t[np];
  Double_t* wfA = new Double_t[np];
  for(i = 0;i<np;i++){
    wfA[i] = amplValues_usbwc[chID][i]/10.0*0.00061;
    wfT[i] = i*SamplingPeriod_usbwc/1000.0;
  }
  waveformSimple *wf = new waveformSimple( wfT, wfA, np, npBaseLine);
  //waveformSimple *wf = new waveformSimple( wfT, wfA, np, npBaseLine, nSplinePoints);
  wf->findPar(tLeft,tRight,vVal);
  //wf->normalizeWF(wf->GetAmax());
  //wf->alignInTimeOfMaximum(150);
  wf->Draw(key_lineOn);
  //wf->Draw("no");
  if(key_printInfoTrue == "printInfoTrue"){
    wf->printPar();
  }
}

void run1vis::Loop(){
  //   In a ROOT session, you can do:
  //      Root > .L run1vis.C
  //      Root > run1vis t
  //      Root > t.GetEntry(12); // Fill t data members with entry number 12
  //      Root > t.Show();       // Show values of entry 12
  //      Root > t.Show(16);     // Read and show values of entry 16
  //      Root > t.Loop();       // Loop on all entries
  //
  //     This is the loop skeleton where:
  //    jentry is the global entry number in the chain
  //    ientry is the entry number in the current Tree
  //  Note that the argument to GetEntry must be:
  //    jentry for TChain::GetEntry
  //    ientry for TTree::GetEntry and TBranch::GetEntry
  //
  //       To read only selected branches, Insert statements like:
  // METHOD1:
  //    fChain->SetBranchStatus("*",0);  // disable all branches
  //    fChain->SetBranchStatus("branchname",1);  // activate branchname
  // METHOD2: replace line
  //    fChain->GetEntry(jentry);       //read all branches
  //by  b_branchname->GetEntry(ientry); //read only this branch
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();
  cout<<nentries<<endl;
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;
  }
}
