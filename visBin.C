#define visBin_cxx
#include "visBin.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "./wfSimulation/src/waveformSimple.hh";

void visBin::Loop(){
  //  In a ROOT session, you can do:
  //  Root > .L visBin.C
  //  Root > visBin t
  //  Root > t.GetEntry(12); // Fill t data members with entry number 12
  //  Root > t.Show();       // Show values of entry 12
  //  Root > t.Show(16);     // Read and show values of entry 16
  //  Root > t.Loop();       // Loop on all entries
  //
  // This is the loop skeleton where:
  // jentry is the global entry number in the chain
  // ientry is the entry number in the current Tree
  // Note that the argument to GetEntry must be:
  // jentry for TChain::GetEntry
  // ientry for TTree::GetEntry and TBranch::GetEntry
  //
  // To read only selected branches, Insert statements like:
  // METHOD1:
  //    fChain->SetBranchStatus("*",0);  // disable all branches
  //    fChain->SetBranchStatus("branchname",1);  // activate branchname
  // METHOD2: replace line
  //    fChain->GetEntry(jentry);       //read all branches
  //by  b_branchname->GetEntry(ientry); //read only this branch
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();
  cout<<"nentries = "<<nentries<<endl;
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;
  }
}

void visBin::DrawOneChannel(Int_t evID, Int_t chID){
   if (fChain == 0) return;
   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nbytes = 0, nb = 0;
   Long64_t jentry = evID;
   Long64_t ientry = LoadTree(jentry);
   nb = fChain->GetEntry(jentry);   nbytes += nb;
   Int_t i = 0;
   const Int_t np = 1024;
   const Int_t nCh = 8;
   Float_t t[np];
   Float_t A[np];
   for(i = 0;i<np;i++){
     A[i] = amplValues_usbwc[chID][i]/10.0*0.00061;
     t[i] = i*SamplingPeriod_usbwc/1000.0;
   }
   TString grtitle;
   grtitle  = " EventID = "; grtitle += evID;
   grtitle += " chID = "; grtitle += chID;
   TGraph *gr = new TGraph(np,t,A); 
   gr->SetTitle(grtitle.Data());
   gr->GetXaxis()->SetTitle("Time, ns");
   gr->GetYaxis()->SetTitle("Amplitude, V");
   TCanvas *c1 = new TCanvas("c1","canva",10,10,1200,800);
   gr->Draw("APL");
}

void visBin::DrawOneChannel_waveformSimple(Int_t evID, Int_t chID, TString key_lineOn, TString key_printInfoTrue,TString fileAnaSetups){
   readAnaSetup(fileAnaSetups);
   if (fChain == 0) return;
   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nbytes = 0, nb = 0;
   Long64_t jentry = evID;
   Long64_t ientry = LoadTree(jentry);
   nb = fChain->GetEntry(jentry);   nbytes += nb;
   Int_t i = 0;
   const Int_t np = 1024;
   const Int_t nCh = 8;
   Int_t npBaseLine;
   Double_t tLeft;
   Double_t tRight;
   Double_t vVal;
   npBaseLine = _nPointBaseLine[chID];
   tLeft  = _chargeWinL[chID];
   tRight = _chargeWinR[chID];
   vVal = _sigConstTh[chID];
   Double_t* wfT = new Double_t[np];
   Double_t* wfA = new Double_t[np];
   for(i = 0;i<np;i++){
     wfA[i] = amplValues_usbwc[chID][i]/10.0*0.00061;
     wfT[i] = i*SamplingPeriod_usbwc/1000.0;
   }
   waveformSimple *wf = new waveformSimple( wfT, wfA, np, npBaseLine);
   wf->findPar(tLeft,tRight,vVal);
   wf->Draw(key_lineOn);
   if(key_printInfoTrue == "printInfoTrue"){
     wf->printPar();
   }
}

void visBin::readAnaSetup(TString fileAnaSetups){
  /////////WF analiz////////////
  ifstream indataSetup;
  TString mot;
  indataSetup.open(fileAnaSetups.Data()); 
  assert(indataSetup.is_open());  
  indataSetup>>mot; indataSetup>>mot; indataSetup>>mot;
  indataSetup>>mot; indataSetup>>mot;
  for(Int_t i = 0;i<_USBWCnChannels;i++){
    indataSetup>>mot;
    indataSetup>>_nPointBaseLine[i];
    indataSetup>>_sigConstTh[i];
    indataSetup>>_chargeWinL[i];
    indataSetup>>_chargeWinR[i];
  }
  /////////////////////////////////////////
  cout<<setw(20)<<"channel"<<setw(20)<<"nPointBaseLine"<<setw(20)<<"sigConstTh"<<setw(20)<<"chargeWinL"<<setw(20)<<"chargeWinR"<<endl;
  for(i = 0;i<_USBWCnChannels;i++){
    cout<<setw(20)<<i
	<<setw(20)<<_nPointBaseLine[i]
	<<setw(20)<<_sigConstTh[i]
	<<setw(20)<<_chargeWinL[i]
	<<setw(20)<<_chargeWinR[i]<<endl;
  }
  //////////////////////////////
}
