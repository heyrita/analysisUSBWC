//my
#include "waveformSimple.hh"

//root
#include <TGraph.h>
#include <TString.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TPaveLabel.h>
#include <TString.h>
#include <TSpline.h>
#include <TMultiGraph.h>
#include <TFile.h>
#include <TH2D.h>
#include <TSpline.h>

//C, C++
#include <iostream>
#include <fstream>
#include <assert.h>
#include <iomanip>
#include <stdio.h>

using namespace std;

waveformSimple::waveformSimple(const Double_t *wfT, const Double_t *wfA, Int_t nPoints, Int_t nPointsBaseLine){
  _nSplinePoints = 0;
  _nPoints = nPoints;
  _wfA = new double[_nPoints];
  _wfT = new double[_nPoints];
  for(int i = 0;i<_nPoints;i++){
      _wfA[i] = wfA[i];
  }
  double abasel = findBaseLineAmpl(nPointsBaseLine);
  _abasel = abasel;
  for(int i = 0;i<_nPoints;i++)
    _wfA[i] = _wfA[i] - abasel;
  double amin = _wfA[0];
  double amax = _wfA[0];
  for(int i = 0;i<_nPoints;i++){
    _wfT[i] = wfT[i];
    if(_wfA[i]<amin)
      amin = _wfA[i];
    if(_wfA[i]>amax)
      amax = _wfA[i];
  }
  double aminabs = TMath::Abs(amin);
  double amaxabs = TMath::Abs(amax);
  if(amaxabs<aminabs)
    for(int i = 0;i<_nPoints;i++)
      _wfA[i] = -_wfA[i];
}

waveformSimple::waveformSimple(const Double_t *wfT, const Double_t *wfA, Int_t nPoints, Int_t nPointsBaseLine, Int_t nSplinePoints){
  if(nSplinePoints<=0){
    cout<<"  ERROR --> nSplinePoints <= 0"<<endl
	<<"            nSplinePoints  = "<<nSplinePoints<<endl;
    assert(0);
  }
  _nSplinePoints = nSplinePoints;
  int i = 0;
  int j = 0;
  int k = 0;
  double *wfAtmp = new double[nPoints];
  double *wfTtmp = new double[nPoints];
  for(i = 0;i<nPoints;i++){
    wfAtmp[i] = wfA[i];
    wfTtmp[i] = wfT[i];
  }
  /////////////////////////////////////////////////////////////
  TSpline3 *spline3=0;
  spline3 = new TSpline3("wfspline", wfTtmp, wfAtmp, nPoints);
  _nPoints = nPoints + (nPoints-1)*nSplinePoints;
  _wfA = new double[_nPoints];
  _wfT = new double[_nPoints];
  Double_t dTime;
  for(i = 0;i<(nPoints-1);i++){
    dTime = wfT[i+1] - wfT[i];
    for(j=0;j<=nSplinePoints;j++){
      _wfT[k] = wfT[i] + dTime/(nSplinePoints+1)*j;
      _wfA[k] = spline3->Eval(_wfT[k]);
      //cout<<"_wfT["<<k<<"] "<< _wfT[k]<<endl;
      k++;
    }
  }
  _wfT[_nPoints-1] = wfT[nPoints-1];
  _wfA[_nPoints-1] = spline3->Eval(_wfT[_nPoints-1]);
  /////////////////////////////////////////////////////////////
  double abasel = findBaseLineAmpl(nPointsBaseLine,wfA);
  _abasel = abasel;
  for(int i = 0;i<_nPoints;i++)
    _wfA[i] = _wfA[i] - abasel;
  double amin = _wfA[0];
  double amax = _wfA[0];
  for(int i = 0;i<_nPoints;i++){
    if(_wfA[i]<amin)
      amin = _wfA[i];
    if(_wfA[i]>amax)
      amax = _wfA[i];
  }
  double aminabs = TMath::Abs(amin);
  double amaxabs = TMath::Abs(amax);
  if(amaxabs<aminabs)
    for(int i = 0;i<_nPoints;i++)
      _wfA[i] = -_wfA[i];
  delete spline3;
  delete wfAtmp;
  delete wfTtmp;
}

waveformSimple::~waveformSimple(){
  delete _wfT;
  delete _wfA;
}

void waveformSimple::normalizeWF(double normVal){
  if(normVal>0.0){
    for(int i = 0;i<_nPoints;i++)
      _wfA[i] = _wfA[i]/normVal;
  }
  else{
    cout<<" ERROR ---> normVal<= 0.0"<<endl
	<<"            normVal = "<<normVal<<endl;
    assert(0);
  }
}

Double_t waveformSimple::Get_wfAi(int i){
  if(i>=0 && i<_nPoints)
    return _wfA[i];
  assert(0);
  return -999.0;
}

Double_t waveformSimple::Get_wfTi(int i){
  if(i>=0 && i<_nPoints)
    return _wfT[i];
  assert(0);
  return -999.0;
}

void waveformSimple::alignInTimeOfMaximum(double tVal){
  double dT = _wfT[_i_amax] - tVal;
  for(Int_t i = 0;i<_nPoints;i++)
    _wfT[i] = _wfT[i] - dT;
}

Double_t waveformSimple::findBaseLineAmpl(Int_t np){
  Double_t abasel = 0;
  for(int i = 0;i<np;i++)
    abasel += _wfA[i];
  return abasel/np;
}

Double_t waveformSimple::findBaseLineAmpl(Int_t np, const Double_t *wfA){
  Double_t abasel = 0;
  for(int i = 0;i<np;i++)
    abasel += wfA[i];
  return abasel/np;
}

Double_t waveformSimple::findChargeInWindow(Double_t tLeft,Double_t tRight){
  Double_t chargeInWindow = 0.0;
  Double_t dT = 0.0;
  for(Int_t i = 0;i<(_nPoints - 1);i++){
    dT = (_wfT[i+1] - _wfT[i]);
    if(dT<0.0){
      cout<<" ERROR -> dT<0.0, dT = (_wfT[i+1] - _wfT[i])"<<endl
	  <<"dT        = "<<dT<<endl
	  <<"_wfT[i+1] = "<<_wfT[i+1]<<endl
	  <<"_wfT[i]   = "<<_wfT[i]<<endl
	  <<"i         = "<<i<<endl;
      assert(0);
    }
    if(_wfT[i]>=tLeft && _wfT[i]<=tRight)
      chargeInWindow = chargeInWindow + _wfA[i]*dT;
  }
  return chargeInWindow;
}

void waveformSimple::findPar(Double_t tLeft,Double_t tRight,Double_t  vVal){
  _chargeWindow_tLeft = tLeft;
  _chargeWindow_tRight = tRight;
  _vVal = vVal;
  double timeCFD;        //cfdVal = 0.5;
  double timeCFD_fall;   //cfdVal = 0.5;
  double timeCFD01;      //cfdVal = 0.1;
  double timeCFD01_fall; //cfdVal = 0.1;
  double timeCFD09 ;     //cfdVal = 0.9;
  double timeCFD09_fall; //cfdVal = 0.9;
  double sigWidth = 0;
  double sigRise = 0;
  double sigFall = 0;
  double chargeInW = 0;
  double timeConstThreas = 0;
  double amin = _wfA[0];
  double amax = _wfA[0];
  double tmin = _wfT[0];
  double tmax = _wfT[0];
  Int_t i_amin = 0;
  Int_t i_amax = 0;
  for(int i = 0;i<_nPoints;i++){
    if(_wfA[i]<amin){
      amin = _wfA[i];
      tmin = _wfT[i];
      i_amin = i;
    }
    if(_wfA[i]>amax){
      amax = _wfA[i];
      tmax = _wfT[i];
      i_amax = i;
    }
  }
  _i_amax = i_amax;
  Double_t cfdVal = 0.5;
  for(Int_t i = i_amax; i>0 ; i--){
    if(_wfA[i]>=amax*cfdVal && _wfA[i-1]<amax*cfdVal){
      timeCFD = linInterpol(amax*cfdVal,
			 _wfA[i],_wfA[i-1],
			 _wfT[i],_wfT[i-1]);
      break;
    }
  }
  for(Int_t i = i_amax; i<(_nPoints-1); i++){
    if(_wfA[i]>=amax*cfdVal && _wfA[i+1]<amax*cfdVal){
      timeCFD_fall = linInterpol(amax*cfdVal,
				 _wfA[i+1],_wfA[i],
				 _wfT[i+1],_wfT[i]);
      break;
    }
  }

  cfdVal = 0.1;
  for(Int_t i = i_amax; i>0 ; i--){
    if(_wfA[i]>=amax*cfdVal && _wfA[i-1]<amax*cfdVal){
      timeCFD01 = linInterpol(amax*cfdVal,
			    _wfA[i],_wfA[i-1],
			    _wfT[i],_wfT[i-1]);
      break;
    }
  }
  for(Int_t i = i_amax; i<(_nPoints-1); i++){
    if(_wfA[i]>=amax*cfdVal && _wfA[i+1]<amax*cfdVal){
      timeCFD01_fall = linInterpol(amax*cfdVal,
				 _wfA[i+1],_wfA[i],
				 _wfT[i+1],_wfT[i]);
      break;
    }
  }

  cfdVal = 0.9;
  for(Int_t i = i_amax; i>0 ; i--){
    if(_wfA[i]>=amax*cfdVal && _wfA[i-1]<amax*cfdVal){
      timeCFD09 = linInterpol(amax*cfdVal,
			    _wfA[i],_wfA[i-1],
			    _wfT[i],_wfT[i-1]);
      break;
    }
  }
  for(Int_t i = i_amax; i<(_nPoints-1); i++){
    if(_wfA[i]>=amax*cfdVal && _wfA[i+1]<amax*cfdVal){
      timeCFD09_fall = linInterpol(amax*cfdVal,
				 _wfA[i+1],_wfA[i],
				 _wfT[i+1],_wfT[i]);
      break;
    }
  }
  
  chargeInW = findChargeInWindow(tmax - tLeft,tmax + tRight);

  for(Int_t i = i_amax; i>0 ; i--){
    if(_wfA[i]>=vVal && _wfA[i-1]<=vVal){
      timeConstThreas = linInterpol(vVal,
				    _wfA[i],_wfA[i-1],
				    _wfT[i],_wfT[i-1]);
      //cout<<"timeConstThreas = "<<timeConstThreas<<endl;
      break;
    }
  }
  
  sigWidth = timeCFD_fall -  timeCFD;
  sigRise  = timeCFD09 - timeCFD01;
  sigFall  = timeCFD01_fall - timeCFD09_fall;
  
  _sigWidth = sigWidth;
  _sigRise = sigRise;
  _sigFall = sigFall;
  _timeCFD = timeCFD;  
  _chargeInW = chargeInW;
  _timeConstThreas = timeConstThreas;
  _amin = amin;
  _amax = amax;
  _tmin = tmin;
  _tmax = tmax;  
  _timeCFD_fall = timeCFD_fall;
  _timeCFD01 = timeCFD01;
  _timeCFD01_fall = timeCFD01_fall;
  _timeCFD09 =  timeCFD09;
  _timeCFD09_fall =  timeCFD09_fall; 
}

Double_t waveformSimple::linInterpol(Double_t y,Double_t y2,Double_t y1,Double_t x2,Double_t x1){
  if(x2<x1 || (y2 - y1)==0.0){
    cout<<endl<<" ERROR--> x2<x1 || (y2 - y1)==0.0"<<endl
	<<" x1 = "<<x1<<endl
	<<" x2 = "<<x2<<endl
	<<" y1 = "<<y1<<endl
	<<" y2 = "<<y2<<endl;
    assert(0);
  }
  return ((y-y2)*(x2-x1)/(y2-y1)+x2);
}

void waveformSimple::printPar(){
  cout<<"_nPoints         = "<<_nPoints<<endl
      <<"_nSplinePoints   = "<<_nSplinePoints<<endl
      <<"_abasel          = "<<_abasel<<endl
      <<"_sigWidth        = "<<_sigWidth<<endl
      <<"_sigRise         = "<<_sigRise<<endl
      <<"_sigFall         = "<<_sigFall<<endl
      <<"_timeCFD         = "<<_timeCFD<<endl
      <<"_chargeInW       = "<<_chargeInW<<endl;
  cout<<"_timeConstThreas = "<<_timeConstThreas<<endl
      <<"_amin            = "<<_amin<<endl
      <<"_amax            = "<<_amax<<endl
      <<"_tmin            = "<<_tmin<<endl
      <<"_tmax            = "<<_tmax<<endl;
  cout<<"_timeCFD01       = "<<_timeCFD01<<endl
      <<"_timeCFD01_fall  = "<<_timeCFD01_fall<<endl
      <<"_timeCFD_fall    = "<<_timeCFD_fall<<endl
      <<"_timeCFD09       = "<<_timeCFD09<<endl
      <<"_timeCFD09_fall  = "<<_timeCFD09_fall<<endl;
  cout<<"_chargeWindow_tLeft  = "<<_chargeWindow_tLeft<<endl;
  cout<<"_chargeWindow_tRight = "<<_chargeWindow_tRight<<endl;
  cout<<"_chargeWindow_tLeft  val = "<<_tmax - _chargeWindow_tLeft<<endl;
  cout<<"_chargeWindow_tRight val = "<<_tmax + _chargeWindow_tRight<<endl;
}

void waveformSimple::Draw(TString key_lineOn){

  Double_t amplMax;
  Double_t amplMin;

  if(key_lineOn == "yes"){  
    amplMax = _amax + 0.002;
    amplMin = _amin - 0.002;
  }
  else{
    amplMax =  1.2;
    amplMin = -0.2;
  }

  TCanvas *c1 = new TCanvas("c1","wavefor",10,10,800,800);
  c1->SetFillColor(kWhite);

  TGraph *gr1 = new TGraph( _nPoints, _wfT, _wfA); 
  gr1->GetXaxis()->SetRangeUser(_wfT[0],_wfT[_nPoints-1]);
  gr1->SetMaximum(amplMax);
  gr1->SetMinimum(amplMin);

  TLine *ln_abasel = new TLine( _wfT[0], 0.0, _wfT[_nPoints-1], 0.0);
  TLine *ln_amin = new TLine( _wfT[0], _amin, _wfT[_nPoints-1], _amin);
  TLine *ln_amax = new TLine( _wfT[0], _amax, _wfT[_nPoints-1], _amax);
  TLine *ln_tmin = new TLine( _tmin, amplMax, _tmin, amplMin);
  TLine *ln_tmax = new TLine( _tmax, amplMax, _tmax, amplMin);
  TLine *ln_chargeWindow_tLeft  = new TLine( _tmax -  _chargeWindow_tLeft, amplMax, _tmax -  _chargeWindow_tLeft, amplMin);
  TLine *ln_chargeWindow_tRight = new TLine( _tmax + _chargeWindow_tRight, amplMax, _tmax + _chargeWindow_tRight, amplMin);
  TLine *ln_RT = new TLine( _timeCFD01, _amax*0.1, _timeCFD09, _amax*0.9);
  TLine *ln_FT = new TLine( _timeCFD09_fall, _amax*0.9, _timeCFD01_fall, _amax*0.1);
  TLine *ln_WT = new TLine( _timeCFD, _amax*0.5, _timeCFD_fall, _amax*0.5);
  TLine *ln_FixT = new TLine( _wfT[0], _vVal, _timeConstThreas, _vVal);

  ln_abasel->SetLineColor(kRed);
  ln_amin->SetLineColor(kMagenta);
  ln_amax->SetLineColor(kBlue);
  ln_tmin->SetLineColor(kMagenta);
  ln_tmax->SetLineColor(kBlue);
  ln_chargeWindow_tLeft->SetLineColor(kRed);
  ln_chargeWindow_tRight->SetLineColor(kRed);
  ln_RT->SetLineColor(kRed+2);
  ln_FT->SetLineColor(kRed+2);
  ln_WT->SetLineColor(kRed);
  ln_FixT->SetLineColor(kBlack);

  ln_abasel->SetLineWidth(2);
  ln_amin->SetLineWidth(2);
  ln_amax->SetLineWidth(2);
  ln_tmin->SetLineWidth(2);
  ln_tmax->SetLineWidth(2);
  ln_chargeWindow_tLeft->SetLineWidth(2);
  ln_chargeWindow_tRight->SetLineWidth(2);
  ln_RT->SetLineWidth(2);
  ln_FT->SetLineWidth(2);
  ln_WT->SetLineWidth(2);
  ln_FixT->SetLineWidth(2);

  ln_abasel->SetLineStyle(kDashed);
  ln_FixT->SetLineStyle(kDashed);

  /*
  lnAmplPos->SetLineColor(kRed);
  lnAmplNeg->SetLineColor(kRed);
  lnTimeCFDPos->SetLineColor(kGreen);
  lnChanDerAftFirPosAmplAmpl->SetLineStyle(kDashed);
  lnChanDerAftFirPosAmplTime->SetLineStyle(kDashed);
  lnChanDerAftFirPosAmplAmpl->SetLineColor(kMagenta);
  lnChanDerAftFirPosAmplTime->SetLineColor(kMagenta);
  lnChanDerBefFirPosAmplAmpl->SetLineStyle(kDashed);
  lnChanDerBefFirPosAmplTime->SetLineStyle(kDashed);
  lnChanDerBefFirPosAmplAmpl->SetLineColor(kYellow);
  lnChanDerBefFirPosAmplTime->SetLineColor(kYellow);
  lnFirstRiseTime->SetLineColor(kBlue);
  lnFirstFallTime->SetLineColor(kBlue);
  lnFirstWidthTime->SetLineColor(kRed);
  lnSignalTh->SetLineStyle(kDashDotted);
  lnSignalTh->SetLineColor(kRed);
  lnCrossTalkTh->SetLineStyle(kDashDotted);
  lnCrossTalkTh->SetLineColor(kRed);
  */

  /*
  srtChargeTOT = "ChargeTOT ";
  srtChargeTOT += floor(_chargeTOT*100.0)/100.0;
  TPaveLabel *lchargeTOT = new TPaveLabel(_wfT[_nPoints-1] - _wfT[_nPoints-1]*0.2
					  ,amplMax-0.1, _wfT[_nPoints-1],
					  amplMax, srtChargeTOT.Data());
  lchargeTOT->SetFillColor(kWhite);
  */

  
  gr1->SetMarkerStyle(20);
  gr1->SetMarkerSize(1);

  gr1->SetTitle("");
  gr1->Draw("APL");

  gr1->GetXaxis()->SetTitle("Time, ns");
  gr1->GetYaxis()->SetTitle("Amplitude, V");

  if(key_lineOn == "yes"){
    ln_abasel->Draw("same");
    ln_amin->Draw("same");
    ln_amax->Draw("same");
    ln_tmin->Draw("same");
    ln_tmax->Draw("same");
    ln_chargeWindow_tLeft->Draw("same");
    ln_chargeWindow_tRight->Draw("same");
    ln_RT->Draw("same");
    ln_FT->Draw("same");
    ln_WT->Draw("same");
    ln_FixT->Draw("same");
  }

  //lnTimeCFDNeg->Draw("same");
  //lnAmplCFDNeg->Draw("same");

  //lchargeTOT->Draw();
  //lchargeTOT_p->Draw();
}
