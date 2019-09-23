//my
#include "src/wfSim.hh"
#include "src/waveform.hh"

//root
#include "TROOT.h"
#include "TString.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TRandom3.h"
#include "TGraph.h"

//C, C++
#include <iostream>
#include <stdlib.h>
#include <assert.h>
#include <fstream>
#include <iomanip>
#include <time.h>

using namespace std;

const Int_t nMcpPMmax = 100000;

int main(int argc, char *argv[]){

  clock_t start, finish;
  start = clock();

  if(argc == 5 && atoi(argv[1])==0){
    /////DTIMETEST/////
    TString histF = argv[2];
    Int_t Nwf = atoi(argv[3]);
    Int_t MySeed = atoi(argv[4]);
    Int_t i = 0;
    Int_t j = 0;
    cout<<"--> wavefrom dTime test <--"<<endl
	<<" runID  = "<<atoi(argv[1])<<endl
	<<" histF  = "<<histF<<endl
	<<" Nwf    = "<<Nwf<<endl
	<<" MySeed = "<<MySeed<<endl;

    TRandom3 *rnd = new TRandom3(MySeed); 
    double t0;
    double dTime;
    double ctdt;
    Double_t time1;
    Double_t time1ConstTh;
    Double_t time1PosAtLevelRiseEdge50mV;
    Double_t time1PosAtLevelFallEdge50mV;
    Double_t time1PosAtLevelRiseEdge100mV;
    Double_t time1PosAtLevelFallEdge100mV;
    Double_t time1PosAtLevelRiseEdge150mV;
    Double_t time1PosAtLevelFallEdge150mV;



    /// WF simulation ///
    Double_t noiseRMS = 0.0015;
    /////////////////////

    //// WF analisys ////
    const Int_t nSpline = 0;
    double cfdL = 0.5;
    double constTh = 0.15;
    double constTh_50mV  = 0.05;//V
    double constTh_100mV = 0.10;//V
    double constTh_150mV = 0.15;//V
    Int_t nPointsNoise = 15;
    Double_t signalThreshold = 30.0/1000.0;
    Double_t crossTalkThreshold = -10.0/1000.0;
    const Int_t nP = 256 + 255*nSpline;
    /////////////////////
      
    ////////HISTOS///////////
    TH1D *h1TimeFirstPos = new TH1D("h1TimeFirstPos","Time First pos",10000,-10.0,100);
    /////////////////////////

    t0 = rnd->Uniform(0.0,5.0);
    wfSim *wfS1 = new wfSim( rnd, 256, 0.3125, t0);
    
    for( i=0; i<Nwf; i++){
      if(i%1000==0)
	cout<<i<<endl;

      t0 = rnd->Uniform(0.0,5.0);
      wfS1->SetDigitTime0(t0);
      
      //gaus
      wfS1->genGaussWF( 20.0, 0.2);
      
      if(noiseRMS != 0.0){
	wfS1->generateNoiseGauss(noiseRMS);
      }

      waveform *wf1  = new waveform( wfS1->getWFT(), wfS1->getWFA(), 256, nSpline);
      wf1->SetSignalThreshold(signalThreshold);
      wf1->SetCrossTalkThreshold(crossTalkThreshold);
      wf1->findParametersOftheWaveform();
      time1 = wf1->findFirstTimePosCFD(cfdL);
      time1ConstTh = wf1->findFirstTimePosConstThreas(constTh);
      
      time1PosAtLevelRiseEdge50mV  = wf1->findFirstTimePosAtLevelRiseEdge(constTh_50mV);
      time1PosAtLevelRiseEdge100mV = wf1->findFirstTimePosAtLevelRiseEdge(constTh_100mV);
      time1PosAtLevelRiseEdge150mV = wf1->findFirstTimePosAtLevelRiseEdge(constTh_150mV);

      time1PosAtLevelFallEdge50mV  = wf1->findFirstTimePosAtLevelFallEdge(constTh_50mV);
      time1PosAtLevelFallEdge100mV = wf1->findFirstTimePosAtLevelFallEdge(constTh_100mV);
      time1PosAtLevelFallEdge150mV = wf1->findFirstTimePosAtLevelFallEdge(constTh_150mV);


      cout<<"time1                        = "<<time1<<endl
	  <<"time1ConstTh                 = "<<time1ConstTh<<endl
	  <<"time1PosAtLevelRiseEdge50mV  = "<<time1PosAtLevelRiseEdge50mV<<endl
	  <<"time1PosAtLevelRiseEdge100mV = "<<time1PosAtLevelRiseEdge100mV<<endl
	  <<"time1PosAtLevelRiseEdge150mV = "<<time1PosAtLevelRiseEdge150mV<<endl
	  <<"time1PosAtLevelFallEdge50mV  = "<<time1PosAtLevelFallEdge50mV<<endl
	  <<"time1PosAtLevelFallEdge100mV = "<<time1PosAtLevelFallEdge100mV<<endl
	  <<"time1PosAtLevelFallEdge150mV = "<<time1PosAtLevelFallEdge150mV<<endl;

      wfS1->CleanWf();

      wf1->PrintWaveFormInfo();      

      delete wf1;
    }

    delete wfS1;

    TFile* rootFile = new TFile(histF.Data(), "RECREATE", " Histograms from wfsim  ", 1);
    rootFile->cd();
    if (rootFile->IsZombie()){
      cout<<"  ERROR ---> file "<<histF<<" is zombi"<<endl;
      assert(0);
    }
    else
      cout<<"  Output Histos file ---> "<<histF<<endl;


  }
  else{
    cout<<"  wavefrom (width ampl)           "<<endl
	<<"  runID [1] = 0                   "<<endl
	<<"        [2] - file with histos    "<<endl
	<<"        [3] - number of waveforms "<<endl
	<<"        [4] - seed                "<<endl;
  }


  finish = clock();

  cout<<" time of work "<<((finish - start)/CLOCKS_PER_SEC )<<endl;


  return 0;
}
