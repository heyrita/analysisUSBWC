//my
#include "waveform.hh"

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

//C, C++
#include <iostream>
#include <fstream>
#include <assert.h>
#include <iomanip>
#include <stdio.h>


using namespace std;

waveform::waveform(const Double_t *wfT, const Double_t *wfA, Int_t nPoints){
  _nPoints = nPoints;
  _wfA = new double[_nPoints];
  _wfT = new double[_nPoints];
  _wfApublic = new double[_nPoints];
  _wfTpublic = new double[_nPoints];
  for(int i = 0;i<_nPoints;i++){
    _wfA[i] = wfA[i];
    _wfT[i] = wfT[i];    
  }
  InitParam();
}

waveform::waveform(const Double_t *wfT, const Double_t *wfA, Int_t nPoints, TString autoInv){
  if(autoInv == "yes"){
    _nPoints = nPoints;
    _wfA = new double[_nPoints];
    _wfT = new double[_nPoints];
    _wfApublic = new double[_nPoints];
    _wfTpublic = new double[_nPoints];
    double amin = wfA[0];
    double amax = wfA[0];
    for(int i = 0;i<_nPoints;i++){
      if(wfA[i]<amin)
	amin = wfA[i];
      if(wfA[i]>amax)
	amax = wfA[i];
    }
    double aminabs = TMath::Abs(amin);
    double amaxabs = TMath::Abs(amax);
    if(amaxabs>=aminabs){
      for(int i = 0;i<_nPoints;i++){
	_wfA[i] = wfA[i];
	_wfT[i] = wfT[i];    
      }
    }
    else{
      for(int i = 0;i<_nPoints;i++){
	_wfA[i] = -wfA[i];
	_wfT[i] = wfT[i];    
      }
    }
    InitParam();
  }
  else{
    assert(0);
  }
}

waveform::waveform(const Double_t *wfT,const Double_t *wfA, Int_t nPoints, Int_t nSplinePoints){
  if(nSplinePoints == 0){
    _nPoints = nPoints;
    _wfA = new double[_nPoints];
    _wfT = new double[_nPoints];
    _wfApublic = new double[_nPoints];
    _wfTpublic = new double[_nPoints];
    for(int i = 0;i<_nPoints;i++){
      _wfA[i] = wfA[i];
      _wfT[i] = wfT[i];    
    }
    InitParam();
  }
  else{
    if(nSplinePoints<0){
      cout<<"  ERROR --> nSplinePoints < 1"<<endl
	  <<"            nSplinePoints  =  "<<nSplinePoints<<endl;
      assert(0);
    }
    int i = 0;
    int j = 0;
    int k = 0;
    double *wfAtmp = new double[nPoints];
    double *wfTtmp = new double[nPoints];
    for(i = 0;i<nPoints;i++){
      wfAtmp[i] = wfA[i];
      wfTtmp[i] = wfT[i];
    }
    TSpline3 *spline3=0;
    spline3 = new TSpline3("wfspline", wfTtmp, wfAtmp, nPoints);
    _nPoints = nPoints + (nPoints-1)*nSplinePoints;
    //cout<<"_nPoints = "<<_nPoints<<endl;
    _wfA = new double[_nPoints];
    _wfT = new double[_nPoints];
    _wfApublic = new double[_nPoints];
    _wfTpublic = new double[_nPoints];
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
    //cout<<"_wfT["<<_nPoints-1<<"] "<< _wfT[_nPoints-1]<<endl
    //<<"_wfA["<<_nPoints-1<<"] "<< _wfA[_nPoints-1]<<endl;
    delete spline3;
    delete wfAtmp;
    delete wfTtmp;
    InitParam();
    SetNpointBaseLine((9*nSplinePoints + 10));
  }
}

waveform::~waveform(){  
  delete _wfT;
  delete _wfA;
  delete _wfApublic;
  delete _wfTpublic;
}

void waveform::makeSmooth(Int_t ponitsBeforeAndAfter){
  Double_t amplaverage = 0.0;
  Double_t *wfA = new double[_nPoints];
  Int_t iStart = 0;
  Int_t iStop = 0;
  Int_t i = 0;
  for(i = 0;i<(_nPoints-1);i++){
    iStart = i - ponitsBeforeAndAfter;
    iStop = i + ponitsBeforeAndAfter;
    if(iStart<0)
      iStart = 0;
    if(iStop>(_nPoints-1))
      iStop = _nPoints-1;
    wfA[i] = findAmplAverage(iStart,iStop);
  }
  for(i = 0;i<_nPoints;i++){
    _wfA[i] = wfA[i];
  }
  delete wfA;
}

Double_t waveform::findAmplAverage(Int_t iStart,Int_t iStop){
  Double_t amplaverage = 0.0;
  Double_t np = 0.0;
  for(Int_t i = iStart;i<=iStop;i++){
    amplaverage = amplaverage + _wfA[i];
  }
  np = iStop - iStart + 1;
  //cout<<"iStart         = "<<iStart<<endl
  //  <<"iStop          = "<<iStop<<endl
  //  <<"amplaverage/np = "<<amplaverage/np<<endl
  //  <<"_wfA[i]        = "<<_wfA[iStop]<<endl;

  return amplaverage/np;
}

void waveform::findParametersOftheWaveforShortList(){
  findAndSubtractBaseLineAmpl();
  findTotalCharge();
  findFirstAmplitudePos();
  findTotAmplitudePos();
  findTotAmplitudeNeg();
  findFirstRiseTimePos();
  findFirstFallTimePos();
  findFirstWidthTimePos();  
}

void waveform::findParametersOftheWaveform(){
  //Int_t i = 0;
  findAndSubtractBaseLineAmpl();
  findTotalCharge();
  findFirstAmplitudePos();
  findFirstAmplitudeNeg();
  findTotAmplitudePos();
  findTotAmplitudeNeg();
  findFirstRiseTimePos();
  findFirstFallTimePos();
  findFirstWidthTimePos();  
  _ChangDerID = findFirstRightChangDer(_FirstAmplitudePosID);
  if(_ChangDerID != -999){
    _ChangDerAmpl = _wfA[_ChangDerID];
    _ChangDerTime = _wfT[_ChangDerID];
  }
  _FirstChangDerID = findFirstLeftChangDer(_FirstAmplitudePosID);
  if(_FirstChangDerID != -999){
    _FirstChangDerAmpl = _wfA[_FirstChangDerID];
    _FirstChangDerTime = _wfT[_FirstChangDerID];
  }
  findMinAndMaxAmplBeforeFirstChangDerAmpl();
  findFirstTimeBasis();
  findChargeFirst();
  findwfID();
  finddTimeFirtsAmplPosChangLeftDer();
  finddTimeFirtsAmplPosFirstTimeAmplNeg();
  //cout<<"_FirstAmplitudePosID  "<<_FirstAmplitudePosID<<endl
  //<<"_FirstAmplitudeNegID  "<<_FirstAmplitudeNegID<<endl;
}

Double_t waveform::findArea(Double_t timeMin, Double_t timeMax){

  Int_t leftID  = findClosestTimeID(timeMin);
  Int_t rightID = findClosestTimeID(timeMax);
  
  //if(leftID>rightID){
  //cout<<" ----> ERROR leftID>rightID"<<endl
  //<<" leftID   "<<leftID<<endl
  //<<" rightID  "<<rightID<<endl;
  //assert(0);
  //}

  if(leftID<0 && leftID>=_nPoints)
    return -999.0;

  if(rightID<0 && (rightID+1)>=_nPoints)
    return -999.0;

  Double_t area = 0.0;  
  for(Int_t i = leftID;i<rightID;i++)
    area = area+_wfA[i]*(_wfT[i+1]-_wfT[i]);

  return area;
}

Int_t waveform::findNumberOfPosIntersectionAtLevel(Double_t level){
  Int_t nINT = 0;
  for(Int_t i = 0;i<(_nPoints-1);i++){
    if(_wfA[i]<=level && _wfA[i+1]>level){
      nINT++;
    }
  }
  return nINT;
}

Int_t waveform::findNumberOfNegIntersectionAtLevel(Double_t level){
  Int_t nINT = 0;
  for(Int_t i = 0;i<(_nPoints-1);i++){
    if(_wfA[i]>=level && _wfA[i+1]<level){
      nINT++;
    }
  }
  return nINT;
}

void waveform::findMinAndMaxAmplBeforeFirstChangDerAmpl(){
  if(_FirstChangDerID !=-999){
    _MaxAmplBeforeFirstChangDerAmpl = _wfA[0];
    _MinAmplBeforeFirstChangDerAmpl = _wfA[0];
    for(Int_t i = 0;i<_FirstChangDerID;i++){
      if(_MaxAmplBeforeFirstChangDerAmpl<_wfA[i])
	_MaxAmplBeforeFirstChangDerAmpl = _wfA[i];
      if(_MinAmplBeforeFirstChangDerAmpl>_wfA[i])
	_MinAmplBeforeFirstChangDerAmpl = _wfA[i];
    }
  }
  else{
    _MaxAmplBeforeFirstChangDerAmpl = -999.0;
    _MinAmplBeforeFirstChangDerAmpl = -999.0;    
  }
}

void  waveform::finddTimeFirtsAmplPosChangLeftDer(){
  if(_FirstTimeAmplitudePos != -999.0 && _FirstChangDerTime != -999.0){
    _dTimeFirtsAmplPosChangLeftDer = _FirstTimeAmplitudePos - _FirstChangDerTime;
  }
  else{
    _dTimeFirtsAmplPosChangLeftDer = -999.0;
  }
}

void  waveform::finddTimeFirtsAmplPosFirstTimeAmplNeg(){
  if(_FirstTimeAmplitudePos != -999.0 && _FirstTimeAmplitudeNeg != -999.0){
    _dTimeFirtsAmplPosFirstTimeAmplNeg = _FirstTimeAmplitudePos - _FirstTimeAmplitudeNeg;
  }
  else{
    _dTimeFirtsAmplPosFirstTimeAmplNeg = -999.0;
  }
}

Int_t waveform::findClosestTimeID(Double_t t0){
  Double_t dT = 999;
  Double_t dTMin = 999;
  Int_t iMin = -999;
  for(int i = 0;i<_nPoints;i++){
    dT = TMath::Abs(_wfT[i] - t0);
    if(dT<dTMin){
      dTMin = dT;
      iMin = i;
    }
  }
  return iMin;
}

bool waveform::makeAlignmentAndNormalisation( Double_t *wfNorAligT, Double_t *wfNorAligA,  
					      Int_t pointID, Double_t t0, Double_t CFDRation){
  Double_t norm = 0.0;
  if(pointID !=-999){
    norm = _wfA[pointID];
  }
  Int_t t0Align      = findClosestTimeID(t0);
  Int_t pointIDalign = findClosestTimeID(findTimeLeftCFD( pointID, CFDRation));
  //t0Align      --> Where 
  //pointIDalign --> Which
  Int_t i;
  //if(pointIDalign != pointID){
  //cout<<"pointIDalign  = "<<pointIDalign<<endl
  //<<"pointID       = "<<pointID<<endl
  //<<"_wfA[pointID] = "<<_wfA[pointID]<<endl;
  //}
  if(norm == 0.0){
    cout<<" ERROR --> norm == 0"<<endl;
    assert(0);
  }  
  if(t0Align == -999 || pointIDalign == -999){
    return false;
  }
  Int_t anmID = 0;
  for(i = 0;i<_nPoints;i++){
    anmID = pointIDalign + (i - t0Align);
    if(anmID>=0 && anmID<_nPoints)
      wfNorAligA[i] = _wfA[anmID]/norm;
    else
      wfNorAligA[i] = 0.0;
    wfNorAligT[i] = _wfT[i];
  }
  return true;
  //if(wfNorAligA[pointID]!=1.0){
  //cout<<wfNorAligA[pointID]<<endl;
  //}
  
  //TGraph *gr1 = new TGraph( _nPoints, wfNorAligT, wfNorAligA); 
  //gr1->SetMaximum(3.0);
  //gr1->SetMinimum(-3.0);
  //gr1->SaveAs("makeAlignmentAndNormalisation.C");
  //tAlign, Int_t pointIDAlign;
}


void waveform::findwfID(){
  //-4 - no left  end and no right end of the signal
  //-3 - no left  end of the signal
  //-2 - no right end of the signal
  //-1 - no time measurements
  // 0 - CT
  // 1 - SP
  // 2 - MP
  if(_FirstAmplitudePosID == -999){
    _wfID = -1;
  }
  else if(_FirstAmplitudePosID != -999 && _ChangDerID == -999 && _FirstChangDerID != -999){
    _wfID = -2;
  }
  else if(_FirstAmplitudePosID != -999 && _ChangDerID != -999 && _FirstChangDerID == -999){
    _wfID = -3;
  }
  else if(_FirstAmplitudePosID != -999 && _ChangDerID == -999 && _FirstChangDerID == -999){
    _wfID = -4;
  }
  else if (_FirstAmplitudePosID != -999 && _ChangDerID != -999 && _FirstChangDerID != -999){
    if(_FirstTimeAmplitudeNeg == -999.0 && (_FirstAmplitudePos - _ChangDerAmpl)>_FirstAmplitudePos*_mpFraction )
      _wfID = 1;
    else if(_FirstTimeAmplitudePos<_FirstTimeAmplitudeNeg && (_FirstAmplitudePos - _ChangDerAmpl)>_FirstAmplitudePos*_mpFraction)
      _wfID = 1;
    else if(_FirstTimeAmplitudePos>_FirstTimeAmplitudeNeg && _FirstTimeAmplitudeNeg !=-999.0)
      _wfID = 0;
    else if(_FirstTimeAmplitudePos>_FirstTimeAmplitudeNeg && _FirstTimeAmplitudeNeg == -999.0 && ( (_FirstAmplitudePos - _ChangDerAmpl)>_FirstAmplitudePos*_mpFraction) )
      _wfID = 1;
    else if((_FirstAmplitudePos - _ChangDerAmpl)<=_FirstAmplitudePos*_mpFraction)
      _wfID = 2;
    else
      _wfID = -1;
  }
  else{
    cout<<"_FirstAmplitudePosID = "<<_FirstAmplitudePosID<<endl
	<<"_ChangDerID          = "<<_ChangDerID<<endl
	<<"_FirstChangDerID     = "<<_FirstChangDerID<<endl;
  }
}

void waveform::findFirstTimeBasis(){
  if(_FirstTimeAmplitude00Pos != -999.0 && _FirstTimeAmplitudeFall00Pos != -999.0)
    _FirstTimeBasis = _FirstTimeAmplitudeFall00Pos - _FirstTimeAmplitude00Pos;
  else
    _FirstTimeBasis = -999.0;
}

void waveform::findChargeFirst(){
  Double_t dT = 0.0;
  if(_FirstChangDerID != -999 && _ChangDerID != -999){
    if (_ChangDerID == (_nPoints-1)){
      _ChangDerID = _ChangDerID - 1;
    }
    for(Int_t i = _FirstChangDerID;i<=_ChangDerID;i++){
      dT = (_wfT[i+1] - _wfT[i]);
      _chargeFirst = _chargeFirst + _wfA[i]*dT;
      //cout<<"   i    = "<<i<<"   "
      //  <<" _wf[i] = "<<_wf[i]<<endl;
    }
    _chargeFirst = _chargeFirst*1000.0/50.0;//pC
    Double_t factor = 1.0/1.6*(1.0e7)/(TMath::Power(10.0,36.0/20.0))/7.0/1.0e5;
    _chargeFirst = _chargeFirst*factor;//p.e.
  }
}

Int_t waveform::findFirstRightChangDer(Int_t FirstAmplitudePosID){
  Double_t dAmpl = -999.0;
  Int_t ChangDerID = -999;
  if(FirstAmplitudePosID != -999){
    for(Int_t i = FirstAmplitudePosID; i<(_nPoints-1); i++){
      dAmpl = _wfA[i] - _wfA[i+1];
      if(dAmpl<0.0){
	ChangDerID = i;
	//cout<<_ChangDerID<<endl
	//  <<_ChangDerAmpl<<endl
	//  <<_ChangDerTime<<endl;
	break;
      }
    }
  }
  else{
    return -999;
  }
  return ChangDerID;
}

Int_t waveform::findFirstLeftChangDer(Int_t FirstAmplitudePosID){
  //cout<<"FirstAmplitudePosID "<<FirstAmplitudePosID<<endl;
  Double_t dAmpl = -999.0;
  Int_t FirstChangDerID = -999;
  if(FirstAmplitudePosID != -999 ){
    for(Int_t i = FirstAmplitudePosID;i>0;i--){
      dAmpl = _wfA[i] - _wfA[i-1];
      //cout<<"i   "<<i<<endl;
      if(dAmpl<0.0){
	FirstChangDerID = i;
	//cout<<FirstChangDerID<<endl;
	//cout<<_ChangDerID<<endl
	//  <<_ChangDerAmpl<<endl
	//  <<_ChangDerTime<<endl;
	break;
      }
    }
  }
  else{
    return -999;
  }
  //cout<<"--> "<<FirstChangDerID<<endl;
  return FirstChangDerID;
}

void waveform::findFirstWidthTimePos(){
  if( _FirstTimeAmplitudeFall05Pos != -999.0 && _FirstTimeAmplitude05Pos != -999.0)
    _FirstWidthTimePos = (_FirstTimeAmplitudeFall05Pos - _FirstTimeAmplitude05Pos);
  else
    _FirstWidthTimePos = -999.0;
}

void waveform::findFirstRiseTimePos(){
  Int_t i;
  if(_FirstAmplitudePosID != -999){
    for(i = _FirstAmplitudePosID; i>0 ; i--){
      if(_wfA[i]>=_FirstAmplitudePos*0.1 && _wfA[i-1]<_FirstAmplitudePos*0.1){
	_FirstTimeAmplitude01Pos = linInterpol(_FirstAmplitudePos*0.1,
					       _wfA[i],_wfA[i-1],
					       _wfT[i],_wfT[i-1]);
	break;
      }
    }
    for(i = _FirstAmplitudePosID; i>0 ; i--){
      if(_wfA[i]>=_FirstAmplitudePos*0.5 && _wfA[i-1]<_FirstAmplitudePos*0.5){
	_FirstTimeAmplitude05Pos = linInterpol(_FirstAmplitudePos*0.5,
					       _wfA[i],_wfA[i-1],
					       _wfT[i],_wfT[i-1]);
	break;
      }
    }    
    for(i = _FirstAmplitudePosID; i>0 ; i--){
      if(_wfA[i]>=_FirstAmplitudePos*0.9 && _wfA[i-1]<_FirstAmplitudePos*0.9){
	_FirstTimeAmplitude09Pos = linInterpol(_FirstAmplitudePos*0.9,
					       _wfA[i],_wfA[i-1],
					       _wfT[i],_wfT[i-1]);
	break;
      }
    }    
    for(i = _FirstAmplitudePosID; i>0 ; i--){
      if(_wfA[i]>=0.0 && _wfA[i-1]<0.0){
	_FirstTimeAmplitude00Pos = linInterpol(0.0,
					       _wfA[i],_wfA[i-1],
					       _wfT[i],_wfT[i-1]);
	break;
      }
    }  
    //cout<<"---------"<<endl
    //<<"_FirstTimeAmplitude01Pos = "<<_FirstTimeAmplitude01Pos<<endl
    //<<"_FirstTimeAmplitude05Pos = "<<_FirstTimeAmplitude05Pos<<endl
    //<<"_FirstTimeAmplitude09Pos = "<<_FirstTimeAmplitude09Pos<<endl
    //<<"_FirstAmplitudePos       = "<<_FirstAmplitudePos<<endl
    //<<"_FirstAmplitudePosID     = "<<_FirstAmplitudePosID<<endl
    //<<"_wfT[i]                  = "<<_wfT[_FirstAmplitudePosID]<<endl;
    if(_FirstTimeAmplitude09Pos != -999.0 && _FirstTimeAmplitude01Pos != -999.0){
      _FirstRiseTimePos = _FirstTimeAmplitude09Pos - _FirstTimeAmplitude01Pos;
    }
  }
  else{
    _FirstTimeAmplitude01Pos = -999.0;
    _FirstTimeAmplitude05Pos = -999.0;
    _FirstTimeAmplitude09Pos = -999.0;
    _FirstTimeAmplitude00Pos = -999.0;
    _FirstRiseTimePos = -999.0;
  }
}

void waveform::findFirstFallTimePos(){
  Int_t i;
  if(_FirstAmplitudePosID != -999){
    for(i = _FirstAmplitudePosID; i<(_nPoints-1) ; i++){
      //cout<<"_FirstAmplitudePos*0.1 "<<_FirstAmplitudePos*0.1<<endl
      //<<"_wfA[i]                "<<_wfA[i]<<endl
      //<<"_wfA[i+1]              "<<_wfA[i+1]<<endl;
      if(_wfA[i]>=_FirstAmplitudePos*0.1 && _wfA[i+1]<_FirstAmplitudePos*0.1){
	_FirstTimeAmplitudeFall01Pos = linInterpol(_FirstAmplitudePos*0.1,
						   _wfA[i+1],_wfA[i],
						   _wfT[i+1],_wfT[i]);
	break;
      }
    }     
    for(i = _FirstAmplitudePosID; i<(_nPoints-1) ; i++){
      if(_wfA[i]>=_FirstAmplitudePos*0.5 && _wfA[i+1]<_FirstAmplitudePos*0.5){
	_FirstTimeAmplitudeFall05Pos = linInterpol(_FirstAmplitudePos*0.5,
						   _wfA[i+1],_wfA[i],
						   _wfT[i+1],_wfT[i]);
	break;
      }
    }
    for(i = _FirstAmplitudePosID; i<(_nPoints-1) ; i++){
      if(_wfA[i]>=_FirstAmplitudePos*0.9 && _wfA[i+1]<_FirstAmplitudePos*0.9){
	_FirstTimeAmplitudeFall09Pos = linInterpol(_FirstAmplitudePos*0.9,
						   _wfA[i+1],_wfA[i],
						   _wfT[i+1],_wfT[i]);
	break;
      }
    }
    for(i = _FirstAmplitudePosID; i<(_nPoints-1) ; i++){
      if(_wfA[i]>=0.0 && _wfA[i+1]<0.0){
	_FirstTimeAmplitudeFall00Pos = linInterpol(0.0,
						   _wfA[i+1],_wfA[i],
						   _wfT[i+1],_wfT[i]);
	break;
      }
    }
    if(_FirstTimeAmplitudeFall01Pos !=-999.0 && _FirstTimeAmplitudeFall09Pos != -999.0){
      _FirstFallTimePos = _FirstTimeAmplitudeFall01Pos - _FirstTimeAmplitudeFall09Pos;
    }
  }
  else{
    _FirstTimeAmplitudeFall01Pos = -999.0;
    _FirstTimeAmplitudeFall05Pos = -999.0;
    _FirstTimeAmplitudeFall09Pos = -999.0;
    _FirstTimeAmplitudeFall00Pos = -999.0;
    _FirstFallTimePos = -999.0;
  }
}

Double_t waveform::findFirstTimeNegCFD(Double_t cfdVal){
  if(_FirstAmplitudeNegID != -999){
    for(Int_t i = _FirstAmplitudeNegID; i>0 ; i--){
      if(_wfA[i]<=_FirstAmplitudeNeg*cfdVal && _wfA[i-1]>_FirstAmplitudeNeg*cfdVal){
	return linInterpol(_FirstAmplitudeNeg*cfdVal,
			   _wfA[i],_wfA[i-1],
			   _wfT[i],_wfT[i-1]);
      }
    }
  }
  return -999.0;
}

Double_t waveform::findFirstTimePosCFD(Double_t cfdVal){
  if(_FirstAmplitudePosID !=-999){
    for(Int_t i = _FirstAmplitudePosID; i>0 ; i--){
      if(_wfA[i]>=_FirstAmplitudePos*cfdVal && _wfA[i-1]<_FirstAmplitudePos*cfdVal){
	return linInterpol(_FirstAmplitudePos*cfdVal,
			   _wfA[i],_wfA[i-1],
			   _wfT[i],_wfT[i-1]);
	//cout<<"i = "<<i<<endl;
      }
    }
  }
  return -999.0;
}

Double_t waveform::findFirstTimePosConstThreas(Double_t vVal){
  if(_FirstAmplitudePosID !=-999 && _wfA[_FirstAmplitudePosID]>vVal){
    for(Int_t i = _FirstAmplitudePosID; i>0 ; i--){
      if(_wfA[i]>=vVal && _wfA[i-1]<vVal){
	return linInterpol(vVal,
			   _wfA[i],_wfA[i-1],
			   _wfT[i],_wfT[i-1]);
	//cout<<"i = "<<i<<endl;
      }
    }
  }
  return -999.0;
}

Double_t waveform::findFirstTimePosAtLevelRiseEdge(Double_t level){
  for(Int_t i = 0;i<(_nPoints-1);i++){
    if(_wfA[i]<=level && _wfA[i+1]>level){
      return linInterpol(level,
			 _wfA[i+1],_wfA[i],
			 _wfT[i+1],_wfT[i]);
    }
  }
  return -999.0;
}

Double_t waveform::findFirstTimePosAtLevelFallEdge(Double_t level){
  for(Int_t i = 0;i<(_nPoints-1);i++){
    if(_wfA[i]>=level && _wfA[i+1]<level){
      return linInterpol(level,
			 _wfA[i+1],_wfA[i],
			 _wfT[i+1],_wfT[i]);
    }
  }
  return -999.0;
}

Double_t waveform::findTimeLeftCFD( Int_t pointID, Double_t cfdVal){
  if(pointID !=-999){
    for(Int_t i = pointID; i>0 ; i--){
      if(_wfA[i]>=_wfA[pointID]*cfdVal && _wfA[i-1]<_wfA[pointID]*cfdVal){
	return linInterpol(_wfA[pointID]*cfdVal,
			   _wfA[i],_wfA[i-1],
			   _wfT[i],_wfT[i-1]);
	//cout<<"i = "<<i<<endl;
      }
    }
  }
  return -999.0;
}

void waveform::findTotAmplitudePos(){
  Double_t AmplMax = -999.0;
  for(Int_t i = 0;i<_nPoints;i++){
    if(_wfA[i]>AmplMax){
      AmplMax = _wfA[i];
      _TotAmplitudePosID = i;
    }
  }
  if(_TotAmplitudePosID>=0 && _TotAmplitudePosID<_nPoints){
    _TotAmplitudePos = _wfA[_TotAmplitudePosID];
    _TotTimeAmplitudePos = _wfT[_TotAmplitudePosID];
  }
}

void waveform::findTotAmplitudeNeg(){
  Double_t AmplMin = _wfA[0];
  for(Int_t i = 0;i<_nPoints;i++){
    if(_wfA[i]<AmplMin){
      AmplMin = _wfA[i];
      _TotAmplitudeNegID = i;
    }
  }
  if(_TotAmplitudeNegID>=0 && _TotAmplitudeNegID<_nPoints){
    _TotAmplitudeNeg = _wfA[_TotAmplitudeNegID];
    _TotTimeAmplitudeNeg = _wfT[_TotAmplitudeNegID];
  }
}

void waveform::findFirstAmplitudeNeg(){
  Double_t dAmpl = -999.0;
  for(Int_t i = 0;i<_nPoints;i++){
    if(_wfA[i]<_crossTalkThreshold){
      if((i+1)<_nPoints){
	dAmpl = _wfA[i] - _wfA[i+1];
	if(dAmpl<0.0){
	  _FirstAmplitudeNegID = i;
	  i = _nPoints;
	}
      }
      else{
	_FirstAmplitudeNegID = -999;
      }
    }
  }
  if(_FirstAmplitudeNegID != -999){
    _FirstAmplitudeNeg = _wfA[_FirstAmplitudeNegID];
    _FirstTimeAmplitudeNeg = _wfT[_FirstAmplitudeNegID]; 
  }
  else{
    _FirstAmplitudeNeg = -999.0;
    _FirstTimeAmplitudeNeg = -999.0;
  }  
}

void waveform::findFirstAmplitudePos(){
  //_FirstAmplitudePosID = -999;
  Double_t dAmpl = -999.0;
  for(Int_t i = 0;i<_nPoints;i++){
    if(_wfA[i]>_SignalThreshold){
      if((i+1)<_nPoints){
	dAmpl = _wfA[i] - _wfA[i+1];
	if(dAmpl>0.0){
	  _FirstAmplitudePosID = i;
	  i = _nPoints;
	}
      }
      else{
	_FirstAmplitudePosID = -999;
      }
    }
  }
  if(_FirstAmplitudePosID != -999){
    _FirstAmplitudePos = _wfA[_FirstAmplitudePosID];
    _FirstTimeAmplitudePos = _wfT[_FirstAmplitudePosID];
  }
  else{
    _FirstAmplitudePos = -999.0;
    _FirstTimeAmplitudePos =  -999.0;
  }
}

Double_t waveform::findChargeInWindow(Double_t tLeft,Double_t tRight){
  _chargeInWindow = 0.0;
  _chargeOutWindow = 0.0;
  Double_t dT = 0.0;
  Double_t _chargeAll = 0.0;
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
      _chargeInWindow = _chargeInWindow + _wfA[i]*dT;
    else
      _chargeOutWindow = _chargeOutWindow + _wfA[i]*dT;
    _chargeAll = _chargeAll + + _wfA[i]*dT;
  }
  return _chargeAll;
}

void waveform::findTotalCharge(){
  _chargeTOT   = 0.0;
  _chargeTOTVT = 0.0;
  _chargeTOT_p = 0.0;
  _chargeTOT_m = 0.0;
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
    _chargeTOT = _chargeTOT + _wfA[i]*dT;
    _chargeTOTVT = _chargeTOTVT + _wfA[i]*dT;
    if(_wfA[i]>0.0)
      _chargeTOT_p = _chargeTOT_p + _wfA[i]*dT;
    else
      _chargeTOT_m = _chargeTOT_m + _wfA[i]*dT;
  }
  //pC
  _chargeTOT   = _chargeTOT*1000.0/50.0;//pC
  _chargeTOT_p = _chargeTOT_p*1000.0/50.0;//pC
  _chargeTOT_m = _chargeTOT_m*1000.0/50.0;//pC
  //number off p.e.
  //number of e charges taking into account amplification of the signal
  //with MCPPMT itself ~7*10^5 and 36 dB amplifier
  Double_t factor = 1.0/1.6*(1.0e7)/(TMath::Power(10.0,36.0/20.0))/7.0/1.0e5;
  _chargeTOT = _chargeTOT*factor;
  _chargeTOT_p = _chargeTOT_p*factor;
  _chargeTOT_m = _chargeTOT_m*factor;
}

void waveform::findSamplingPeriod(){
  if(_nPoints<=2){
    cout<<" ERROR ---> _nPoints < = 2 "<<endl
	<<"            _nPoints   =   "<<_nPoints<<endl;
    assert(0);
  }
  Double_t dT = 0.0;
  Int_t i = 0;
  for( i = 0; i<(_nPoints - 1); i++){
    dT = dT + (_wfT[i+1] - _wfT[i]);
  }
  _SamplingPeriod = dT/(_nPoints - 1);
  dT = 0.0;
  _SamplingPeriodStanDev = 0.0;
  for( i = 0; i<(_nPoints - 1); i++){
    dT = (_wfT[i+1] - _wfT[i]);
    _SamplingPeriodStanDev = _SamplingPeriodStanDev +  (_SamplingPeriod - dT)*(_SamplingPeriod - dT);
  }
  _SamplingPeriodStanDev = TMath::Sqrt(_SamplingPeriodStanDev/(_nPoints - 2));
}

void waveform::findAndSubtractBaseLineAmpl(){
  if(_nPointBaisLine<=1){
    cout<<" ERROR ---> _nPointBaisLine <= 1"<<endl
	<<"            _nPointBaisLine  =  "<<_nPointBaisLine<<endl;
    assert(0);
  }
  Int_t i;
  _BaseLineAmpl = 0.0;
  for(i = 0;i<_nPointBaisLine;i++){
    _BaseLineAmpl = _BaseLineAmpl + _wfA[i];
  }
  _BaseLineAmpl = _BaseLineAmpl/_nPointBaisLine;
  _BaseLineAmplStanDev = 0.0;
  for( i = 0; i<_nPointBaisLine; i++){
    _BaseLineAmplStanDev = _BaseLineAmplStanDev + (_BaseLineAmpl - _wfA[i])*(_BaseLineAmpl - _wfA[i]);
  }
  _BaseLineAmplStanDev = TMath::Sqrt(_BaseLineAmplStanDev/(_nPointBaisLine - 1));
  ///substruct baseline amplitude///
  for( i = 0; i<_nPoints; i++){
    _wfA[i] = _wfA[i] - _BaseLineAmpl;
    _wfApublic[i] = _wfA[i];
    _wfTpublic[i] = _wfT[i];
 }
  _MaxBaseLineAmpl = _wfA[0];
  _MinBaseLineAmpl = _wfA[0];
  for(i = 0;i<_nPointBaisLine;i++){
    if(_MaxBaseLineAmpl<_wfA[i])
      _MaxBaseLineAmpl = _wfA[i];
    if(_MinBaseLineAmpl>_wfA[i])
      _MinBaseLineAmpl = _wfA[i];
  }
}

void waveform::InitParam(){

  _nPointBaisLine = 6;
  _SignalThreshold = 30.0/1000.0;     //V
  _crossTalkThreshold = -10.0/1000.0; //V
  _mpFraction = 0.8;

  _BaseLineAmpl        = -999.0;
  _BaseLineAmplStanDev = -999.0;
  _MaxBaseLineAmpl     = -999.0;
  _MinBaseLineAmpl     = -999.0;

  _SamplingPeriod        = -999.0;
  _SamplingPeriodStanDev = -999.0;

  _TotAmplitudePosID   = -999;
  _TotAmplitudePos     = -999.0;
  _TotTimeAmplitudePos = -999.0;

  _TotAmplitudeNegID   = -999;
  _TotAmplitudeNeg     = -999.0;
  _TotTimeAmplitudeNeg = -999.0;

  _FirstAmplitudePosID   = -999;
  _FirstAmplitudePos     = -999.0;
  _FirstTimeAmplitudePos = -999.0;

  _FirstAmplitudeNegID   = -999;
  _FirstAmplitudeNeg     = -999.0;
  _FirstTimeAmplitudeNeg = -999.0;

  _ChangDerID = -999;
  _ChangDerAmpl = -999.0;
  _ChangDerTime = -999.0;

  _FirstChangDerID = -999;
  _FirstChangDerAmpl = -999.0;
  _FirstChangDerTime = -999.0;
  _MaxAmplBeforeFirstChangDerAmpl = -999.0;
  _MinAmplBeforeFirstChangDerAmpl = -999.0;

  _FirstTimeBasis = -999.0;

  _FirstTimeAmplitude01Pos = -999.0;
  _FirstTimeAmplitude05Pos = -999.0;
  _FirstTimeAmplitude09Pos = -999.0;
  _FirstTimeAmplitude00Pos = -999.0;
  _FirstRiseTimePos = -999.0;

  _FirstTimeAmplitudeFall01Pos = -999.0;
  _FirstTimeAmplitudeFall05Pos = -999.0;
  _FirstTimeAmplitudeFall09Pos = -999.0;
  _FirstTimeAmplitudeFall00Pos = -999.0;
  _FirstFallTimePos = -999.0;

  _FirstWidthTimePos = -999.0;

  _dTimeFirtsAmplPosChangLeftDer = -999.0;
  _dTimeFirtsAmplPosFirstTimeAmplNeg = -999.0;
  _dTimeChangRightDerChangLeftDer = -999.0;

  _chargeTOT = -999.0;
  _chargeTOT_p = -999.0;
  _chargeTOT_m = -999.0;
  _chargeFirst = -999.0;

  _wfID = -999;
}

Double_t waveform::linInterpol(Double_t y,Double_t y2,Double_t y1,Double_t x2,Double_t x1){
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

void waveform::PrintWaveFormInfo(){
  cout<<"----------------"<<endl
      <<"_nPoints             "<<_nPoints<<endl
      <<"_nPointBaisLine      "<<_nPointBaisLine<<endl
      <<"_SignalThreshold     "<<_SignalThreshold<<" V"<<endl
      <<"_crossTalkThreshold  "<<_crossTalkThreshold<<" V"<<endl
      <<"_mpFraction          "<<_mpFraction<<endl;
  cout<<endl<<"----------------"<<endl
      <<"_BaseLineAmpl        "<<_BaseLineAmpl<<" V"<<endl     
      <<"_BaseLineAmplStanDev "<<_BaseLineAmplStanDev<<" V"<<endl
      <<"_MaxBaseLineAmpl     "<<_MaxBaseLineAmpl<<" V"<<endl
      <<"_MinBaseLineAmpl     "<<_MinBaseLineAmpl<<" V"<<endl;
  cout<<endl<<"----------------"<<endl
      <<"_SamplingPeriod        "<<_SamplingPeriod<<" ns"<<endl
      <<"_SamplingPeriodStanDev "<<_SamplingPeriodStanDev<<" ns"<<endl;
  cout<<endl<<"----------------"<<endl
      <<"_TotAmplitudePosID   "<<_TotAmplitudePosID<<endl
      <<"_TotAmplitudePos     "<<_TotAmplitudePos<<" V"<<endl
      <<"_TotTimeAmplitudePos "<<_TotTimeAmplitudePos<<" ns"<<endl;
  cout<<endl<<"----------------"<<endl
      <<"_TotAmplitudeNegID   "<<_TotAmplitudeNegID<<endl
      <<"_TotAmplitudeNeg     "<<_TotAmplitudeNeg<<" V"<<endl
      <<"_TotTimeAmplitudeNeg "<<_TotTimeAmplitudeNeg<<" ns"<<endl;
  cout<<endl<<"----------------"<<endl
      <<"_FirstAmplitudePosID   "<<_FirstAmplitudePosID<<endl
      <<"_FirstAmplitudePos     "<<_FirstAmplitudePos<<" V"<<endl
      <<"_FirstTimeAmplitudePos "<<_FirstTimeAmplitudePos<<" ns"<<endl;
  cout<<endl<<"----------------"<<endl
      <<"_FirstAmplitudeNegID   "<<_FirstAmplitudeNegID<<endl 
      <<"_FirstAmplitudeNeg     "<<_FirstAmplitudeNeg<<" V"<<endl
      <<"_FirstTimeAmplitudeNeg "<<_FirstTimeAmplitudeNeg<<" ns"<<endl;
  cout<<endl<<"----------------"<<endl
      <<"_ChangDerID   "<<_ChangDerID<<endl      
      <<"_ChangDerAmpl "<<_ChangDerAmpl<<" V"<<endl
      <<"_ChangDerTime "<<_ChangDerTime<<" ns"<<endl; 
  cout<<endl<<"----------------"<<endl
      <<"_FirstChangDerID                "<<_FirstChangDerID<<endl
      <<"_FirstChangDerAmpl              "<<_FirstChangDerAmpl<<" V"<<endl
      <<"_FirstChangDerTime              "<<_FirstChangDerTime<<" ns"<<endl
      <<"_MaxAmplBeforeFirstChangDerAmpl "<<_MaxAmplBeforeFirstChangDerAmpl<<" V"<<endl
      <<"_MinAmplBeforeFirstChangDerAmpl "<<_MinAmplBeforeFirstChangDerAmpl<<" V"<<endl;
  cout<<endl<<"----------------"<<endl
      <<"_FirstTimeBasis    "<<_FirstTimeBasis<<" ns"<<endl;
  cout<<endl<<"----------------"<<endl
      <<"_FirstTimeAmplitude01Pos "<<_FirstTimeAmplitude01Pos<<" ns"<<endl 
      <<"_FirstTimeAmplitude05Pos "<<_FirstTimeAmplitude05Pos<<" ns"<<endl
      <<"_FirstTimeAmplitude09Pos "<<_FirstTimeAmplitude09Pos<<" ns"<<endl
      <<"_FirstRiseTimePos        "<<_FirstRiseTimePos<<" ns"<<endl;
  cout<<endl<<"----------------"<<endl
      <<"_FirstTimeAmplitudeFall01Pos "<<_FirstTimeAmplitudeFall01Pos<<" ns"<<endl
      <<"_FirstTimeAmplitudeFall05Pos "<<_FirstTimeAmplitudeFall05Pos<<" ns"<<endl
      <<"_FirstTimeAmplitudeFall09Pos "<<_FirstTimeAmplitudeFall09Pos<<" ns"<<endl
      <<"_FirstFallTimePos            "<<_FirstFallTimePos<<" ns"<<endl;
  cout<<endl<<"----------------"<<endl
      <<"_FirstWidthTimePos "<<_FirstWidthTimePos<<" ns"<<endl;
  cout<<endl<<"----------------"<<endl
      <<"_chargeTOT   "<<_chargeTOT<<" p.e."<<endl
      <<"_chargeTOT_p "<<_chargeTOT_p<<" p.e."<<endl
      <<"_chargeTOT_m "<<_chargeTOT_m<<" p.e."<<endl
      <<"_chargeFirst "<<_chargeFirst<<" p.e."<<endl;
  cout<<endl<<"----------------"<<endl
      <<"_wfID "<<_wfID<<endl;
}

void waveform::Draw(Double_t cfdVal){

  Double_t amplMax =  1.2;
  Double_t amplMin = -0.2;

  TString srtChargeTOT;
  TString srtChargeFirst;
  TString srtChargeTOT_p;
  TString srtChargeTOT_m;
  TString srtWfID;

  TCanvas *c1 = new TCanvas("c1","wavefor",10,10,1000,800);
  c1->SetFillColor(kWhite);

  TGraph *gr1 = new TGraph( _nPoints, _wfT, _wfA); 

  gr1->SetMaximum(amplMax);
  gr1->SetMinimum(amplMin);

  Double_t firstTimePosCFD = findFirstTimePosCFD(cfdVal);

  TLine *lnBase    = new TLine(0.0, 0.0, _wfT[_nPoints-1], 0.0);
  TLine *lnAmplPos = new TLine(0.0,_FirstAmplitudePos, _wfT[_nPoints-1],_FirstAmplitudePos);
  TLine *lnAmplNeg = new TLine(0.0,_FirstAmplitudeNeg, _wfT[_nPoints-1],_FirstAmplitudeNeg);
  TLine *lnTimeCFDPos = new TLine( firstTimePosCFD, amplMax, firstTimePosCFD, amplMin);
  TLine *lnAmplCFDPos = new TLine(0.0,_FirstAmplitudePos*cfdVal,_wfT[_nPoints-1],_FirstAmplitudePos*cfdVal);
  //TLine *lnTimeCFDNeg = new TLine(_FirstTimeNegCFD,amplMax,_FirstTimeNegCFD, amplMin);
  //TLine *lnAmplCFDNeg = new TLine(0.0,_FirstAmplitudeNeg*cfdVal,_wfT[_nPoints-1],_FirstAmplitudeNeg*cfdVal);

  TLine *lnChanDerAftFirPosAmplAmpl = new TLine(0.0,_ChangDerAmpl,_wfT[_nPoints-1],_ChangDerAmpl);
  TLine *lnChanDerAftFirPosAmplTime = new TLine(_ChangDerTime,amplMin,_ChangDerTime,amplMax);

  TLine *lnChanDerBefFirPosAmplAmpl = new TLine(0.0,_FirstChangDerAmpl,_wfT[_nPoints-1],_FirstChangDerAmpl);
  TLine *lnChanDerBefFirPosAmplTime = new TLine(_FirstChangDerTime,amplMin,_FirstChangDerTime,amplMax);
  
  TLine *lnFirstRiseTime = new TLine(_FirstTimeAmplitude01Pos, _FirstAmplitudePos*0.1,
				     _FirstTimeAmplitude09Pos, _FirstAmplitudePos*0.9);
  TLine *lnFirstFallTime = new TLine(_FirstTimeAmplitudeFall01Pos, _FirstAmplitudePos*0.1,
				     _FirstTimeAmplitudeFall09Pos, _FirstAmplitudePos*0.9);

  TLine *lnFirstWidthTime = new TLine(_FirstTimeAmplitude05Pos, _FirstAmplitudePos*0.5,
				      _FirstTimeAmplitudeFall05Pos, _FirstAmplitudePos*0.5);

  TLine *lnSignalTh = new TLine(0.0, _SignalThreshold,
				_wfT[_nPoints-1], _SignalThreshold);
  TLine *lnCrossTalkTh = new TLine(0.0, _crossTalkThreshold,
				   _wfT[_nPoints-1], _crossTalkThreshold);
  
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

  srtChargeTOT = "ChargeTOT ";
  srtChargeTOT += floor(_chargeTOT*100.0)/100.0;
  TPaveLabel *lchargeTOT = new TPaveLabel(_wfT[_nPoints-1] - _wfT[_nPoints-1]*0.2
					  ,amplMax-0.1, _wfT[_nPoints-1],
					  amplMax, srtChargeTOT.Data());
  lchargeTOT->SetFillColor(kWhite);

  srtChargeTOT_p = "ChargeTOT_p ";
  srtChargeTOT_p += floor(_chargeTOT_p*100.0)/100.0;
  TPaveLabel *lchargeTOT_p = new TPaveLabel(_wfT[_nPoints-1] - _wfT[_nPoints-1]*0.2
					    ,amplMax-0.2, _wfT[_nPoints-1],
					    amplMax-0.1, srtChargeTOT_p.Data());
  lchargeTOT_p->SetFillColor(kWhite);

  srtChargeTOT_m = "ChargeTOT_m ";
  srtChargeTOT_m += floor(_chargeTOT_m*100.0)/100.0;
  TPaveLabel *lchargeTOT_m = new TPaveLabel(_wfT[_nPoints-1] - _wfT[_nPoints-1]*0.2
					    ,amplMax-0.3, _wfT[_nPoints-1],
					    amplMax-0.2, srtChargeTOT_m.Data());
  lchargeTOT_m->SetFillColor(kWhite);


  srtChargeFirst = "ChargeFirst ";
  srtChargeFirst += floor(_chargeFirst*100.0)/100.0;
  TPaveLabel *lchargeFirst = new TPaveLabel(_wfT[_nPoints-1] - _wfT[_nPoints-1]*0.2
					    ,amplMax-0.4, _wfT[_nPoints-1],
					    amplMax-0.3, srtChargeFirst.Data());
  lchargeFirst->SetFillColor(kWhite);

  srtWfID = "";
  if(_wfID == -1)
    srtWfID += "-1";
  else if(_wfID == 0)
    srtWfID += "CT";
  else if(_wfID == 1)
    srtWfID += "SP";
  else if(_wfID == 2)
    srtWfID += "MP";
  else
    srtWfID += "-999";
  TPaveLabel *lwfID = new TPaveLabel(_wfT[_nPoints-1] - _wfT[_nPoints-1]*0.2
				     ,amplMax-0.5, _wfT[_nPoints-1],
				     amplMax-0.4, srtWfID.Data());
  lwfID->SetFillColor(kWhite);
  
  gr1->SetMarkerStyle(20);
  gr1->SetMarkerSize(1);

  gr1->Draw("APL");

  lnBase->Draw("same");
  lnAmplPos->Draw("same");
  lnAmplNeg->Draw("same");
  lnTimeCFDPos->Draw("same");
  lnAmplCFDPos->Draw("same");
  //lnTimeCFDNeg->Draw("same");
  //lnAmplCFDNeg->Draw("same");
  lnChanDerAftFirPosAmplAmpl->Draw("same");
  lnChanDerAftFirPosAmplTime->Draw("same");
  lnFirstRiseTime->Draw("same");
  lnFirstFallTime->Draw("same");
  lnFirstWidthTime->Draw("same");
  lnSignalTh->Draw("same");
  lnCrossTalkTh->Draw("same");
  lnChanDerBefFirPosAmplAmpl->Draw("same");
  lnChanDerBefFirPosAmplTime->Draw("same");

  lchargeTOT->Draw();
  lchargeTOT_p->Draw();
  lchargeTOT_m->Draw();
  lchargeFirst->Draw();
  lwfID->Draw();
}

void waveform::DrawCharge( Double_t cfdVal, Double_t dletaLeft, Double_t dletaRight){

  Double_t amplMax =  1.2;
  Double_t amplMin = -0.2;

  TString srtChargeTOT;
  TString srtChargeInWindow;
  TString srtChargeOutWindow;

  TCanvas *c1 = new TCanvas("c1","wavefor",10,10,1000,800);
  c1->SetFillColor(kWhite);

  TGraph *gr1 = new TGraph( _nPoints, _wfT, _wfA); 

  gr1->SetMaximum(amplMax);
  gr1->SetMinimum(amplMin);

  Double_t firstTimePosCFD = findFirstTimePosCFD(cfdVal);
  Double_t sigChargeAllVT = GetChargeVTTOT();
  Double_t t_left = GetTotTimeAmplitudePos() - dletaLeft;
  Double_t t_right = GetTotTimeAmplitudePos() + dletaRight;
  Double_t sigChargeAll = findChargeInWindow(t_left,t_right);
  Double_t sigChargeInWindowVT = GetChargeInWindow();
  Double_t sigChargeOutWindowVT = GetChargeOutWindow();

  TLine *lnBase    = new TLine(0.0, 0.0, _wfT[_nPoints-1], 0.0);
  TLine *lnAmplPos = new TLine(0.0,_FirstAmplitudePos, _wfT[_nPoints-1],_FirstAmplitudePos);
  TLine *lnAmplNeg = new TLine(0.0,_FirstAmplitudeNeg, _wfT[_nPoints-1],_FirstAmplitudeNeg);
  TLine *lnTimeCFDPos = new TLine( firstTimePosCFD, amplMax, firstTimePosCFD, amplMin);
  TLine *lnAmplCFDPos = new TLine(0.0,_FirstAmplitudePos*cfdVal,_wfT[_nPoints-1],_FirstAmplitudePos*cfdVal);
  //TLine *lnTimeCFDNeg = new TLine(_FirstTimeNegCFD,amplMax,_FirstTimeNegCFD, amplMin);
  //TLine *lnAmplCFDNeg = new TLine(0.0,_FirstAmplitudeNeg*cfdVal,_wfT[_nPoints-1],_FirstAmplitudeNeg*cfdVal);

  TLine *lnChanDerAftFirPosAmplAmpl = new TLine(0.0,_ChangDerAmpl,_wfT[_nPoints-1],_ChangDerAmpl);
  TLine *lnChanDerAftFirPosAmplTime = new TLine(_ChangDerTime,amplMin,_ChangDerTime,amplMax);

  TLine *lnChanDerBefFirPosAmplAmpl = new TLine(0.0,_FirstChangDerAmpl,_wfT[_nPoints-1],_FirstChangDerAmpl);
  TLine *lnChanDerBefFirPosAmplTime = new TLine(_FirstChangDerTime,amplMin,_FirstChangDerTime,amplMax);
  
  TLine *lnFirstRiseTime = new TLine(_FirstTimeAmplitude01Pos, _FirstAmplitudePos*0.1,
				     _FirstTimeAmplitude09Pos, _FirstAmplitudePos*0.9);
  TLine *lnFirstFallTime = new TLine(_FirstTimeAmplitudeFall01Pos, _FirstAmplitudePos*0.1,
				     _FirstTimeAmplitudeFall09Pos, _FirstAmplitudePos*0.9);

  TLine *lnFirstWidthTime = new TLine(_FirstTimeAmplitude05Pos, _FirstAmplitudePos*0.5,
				      _FirstTimeAmplitudeFall05Pos, _FirstAmplitudePos*0.5);

  TLine *lnSignalTh = new TLine(0.0, _SignalThreshold,
				_wfT[_nPoints-1], _SignalThreshold);
  TLine *lnCrossTalkTh = new TLine(0.0, _crossTalkThreshold,
				   _wfT[_nPoints-1], _crossTalkThreshold);
  
  TLine *lnChargeWindowL = new TLine( t_left, amplMax, t_left, amplMin);
  TLine *lnChargeWindowR = new TLine( t_right, amplMax, t_right, amplMin);

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

  lnChargeWindowL->SetLineColor(kBlack);
  lnChargeWindowR->SetLineColor(kBlack);
  lnChargeWindowL->SetLineWidth(2);
  lnChargeWindowR->SetLineWidth(2);


  char buffer [50];
  srtChargeTOT = "ChargeTOT ";
  //srtChargeTOT += floor(sigChargeAllVT*100000.0)/100000.0;
  //srtChargeTOT += sigChargeAllVT;
  sprintf (buffer, "%10.5f",sigChargeAllVT);
  srtChargeTOT += buffer;
  TPaveLabel *lchargeTOT = new TPaveLabel(_wfT[_nPoints-1] - _wfT[_nPoints-1]*0.2
					  ,amplMax-0.1, _wfT[_nPoints-1],
					  amplMax, srtChargeTOT.Data());
  lchargeTOT->SetFillColor(kWhite);

  srtChargeInWindow = "Charge in window ";
  //srtChargeInWindow += floor(sigChargeInWindowVT*100000.0)/100000.0;
  //srtChargeInWindow += sigChargeInWindowVT;
  sprintf (buffer, "%10.5f",sigChargeInWindowVT);
  srtChargeInWindow += buffer;
  TPaveLabel *lchargeInWindow = new TPaveLabel(_wfT[_nPoints-1] - _wfT[_nPoints-1]*0.2
					       ,amplMax-0.2, _wfT[_nPoints-1],
					       amplMax-0.1, srtChargeInWindow.Data());
  lchargeInWindow->SetFillColor(kWhite);

  srtChargeOutWindow = "Charge out window ";
  //srtChargeOutWindow += floor(sigChargeOutWindowVT*100000.0)/100000.0;
  //srtChargeOutWindow += sigChargeOutWindowVT;
  sprintf (buffer, "%10.5f",sigChargeOutWindowVT);
  srtChargeOutWindow += buffer;
  TPaveLabel *lchargeOutWindow = new TPaveLabel(_wfT[_nPoints-1] - _wfT[_nPoints-1]*0.2
					    ,amplMax-0.3, _wfT[_nPoints-1],
					    amplMax-0.2, srtChargeOutWindow.Data());
  lchargeOutWindow->SetFillColor(kWhite);

  
  gr1->SetMarkerStyle(20);
  gr1->SetMarkerSize(1);

  gr1->Draw("APL");

  lnBase->Draw("same");
  lnAmplPos->Draw("same");
  lnAmplNeg->Draw("same");
  lnTimeCFDPos->Draw("same");
  lnAmplCFDPos->Draw("same");
  //lnTimeCFDNeg->Draw("same");
  //lnAmplCFDNeg->Draw("same");
  lnChanDerAftFirPosAmplAmpl->Draw("same");
  lnChanDerAftFirPosAmplTime->Draw("same");
  lnFirstRiseTime->Draw("same");
  lnFirstFallTime->Draw("same");
  lnFirstWidthTime->Draw("same");
  lnSignalTh->Draw("same");
  lnCrossTalkTh->Draw("same");
  lnChanDerBefFirPosAmplAmpl->Draw("same");
  lnChanDerBefFirPosAmplTime->Draw("same");
  lnChargeWindowL->Draw("same");
  lnChargeWindowR->Draw("same");

  lchargeTOT->Draw();
  lchargeInWindow->Draw();
  lchargeOutWindow->Draw();

  cout<<"sigChargeAllVT       = "<<sigChargeAllVT<<endl
      <<"sigChargeAll         = "<<sigChargeAll<<endl
      <<"sigChargeInWindowVT  = "<<sigChargeInWindowVT<<endl
      <<"sigChargeOutWindowVT = "<<sigChargeOutWindowVT<<endl;


}

void waveform::DrawSmooth(Int_t ponitsBeforeAndAfter,Int_t numberOfSmoothIterations){

  Double_t amplMax =  1.2;
  Double_t amplMin = -0.2;

  TCanvas *c1 = new TCanvas("c1","wavefor",10,10,1000,800);
  c1->SetFillColor(kWhite);

  Double_t *wfA  = new double[_nPoints];
  Double_t *wfAs = new double[_nPoints];
  Int_t i = 0;
  for(i=0;i<_nPoints;i++){
    wfA[i] = _wfA[i];
  }
  for(i = 0;i<numberOfSmoothIterations;i++)
    makeSmooth(ponitsBeforeAndAfter);
  for(i=0;i<_nPoints;i++){
    wfAs[i] = _wfA[i];
  }

  TGraph *gr1 = new TGraph( _nPoints, _wfT, wfA); 
  TGraph *gr2 = new TGraph( _nPoints, _wfT, wfAs); 

  gr1->SetMaximum(amplMax);
  gr1->SetMinimum(amplMin);
  gr2->SetMaximum(amplMax);
  gr2->SetMinimum(amplMin);

  gr1->SetMarkerStyle(20);
  gr2->SetMarkerStyle(20);
  gr2->SetMarkerColor(632);
  gr2->SetLineColor(632);
  //cout<<"kRed = "<<kRed<<endl;

  TMultiGraph *mg = new TMultiGraph();
  mg->Add(gr1,"pl");
  mg->Add(gr2,"pl");
  mg->SetMaximum(amplMax);
  mg->SetMinimum(amplMin);
  mg->Draw("APL");
}


void waveform::Draw15(Double_t level){

  Double_t amplMax =  1.2;
  Double_t amplMin = -1.2;

  TCanvas *c1 = new TCanvas("c1","wavefor",10,10,1000,800);
  c1->SetFillColor(kWhite);

  TGraph *gr1 = new TGraph( _nPoints, _wfT, _wfA); 

  gr1->SetMaximum(amplMax);
  gr1->SetMinimum(amplMin);

  Double_t timeFirstPosAtLevelFallEdge = findFirstTimePosAtLevelFallEdge(level);

  TLine *lnlevel = new TLine(0.0, level, _wfT[_nPoints-1],level);
  TLine *lnTime = new TLine( timeFirstPosAtLevelFallEdge, amplMax, timeFirstPosAtLevelFallEdge, amplMin);
  
  lnlevel->SetLineColor(kRed);
  lnTime->SetLineColor(kGreen);

  gr1->SetMarkerStyle(20);
  gr1->SetMarkerSize(1);
  gr1->Draw("APL");
  lnlevel->Draw("same");
  lnTime->Draw("same");
}


/*
void waveform::SaveWaveForm2Cfile(TString nameff){
  TGraph *gr1 = new TGraph( usbConst::usbwcNsamplingPoint, _wfTime, _wf); 
  gr1->SetMarkerStyle(20);
  gr1->SetMarkerSize(1);
  gr1->Draw("APL");
  gr1->SaveAs(nameff.Data());
  //PtintWaveFormInfo();
  //graph->SetMarkerStyle(20);
  //graph->SetMarkerSize(1);
  //graph->Draw("APL");
}
*/

/*
void waveform::DrawSpline(){

  Double_t amplMax =  1.2;
  Double_t amplMin = -0.2;

  TCanvas *c1 = new TCanvas("c1","wavefor",10,10,1000,800);
  c1->SetFillColor(kWhite);

  TGraph *gr1 = new TGraph( usbConst::usbwcNPoint, _wfTimeSp, _wfSp); 
  TGraph *gr2 = new TGraph( usbConst::usbwcNsamplingPoint, _wfTime, _wf); 

  gr1->SetMaximum(amplMax);
  gr1->SetMinimum(amplMin);  
  gr1->SetMarkerStyle(20);
  gr1->SetMarkerSize(1);

  gr2->SetMaximum(amplMax);
  gr2->SetMinimum(amplMin);  
  gr2->SetMarkerStyle(20);
  gr2->SetMarkerSize(1);
  gr2->SetMarkerColor(kRed);


  ////////////////
  TLine *lnTimeCFDPos = new TLine(_FirstTimePosCFD,amplMax,_FirstTimePosCFD, amplMin);
  TLine *lnAmplCFDPos = new TLine(0.0,_FirstAmplitudePos*_CFD,_SamplingPeriod*(usbConst::usbwcNsamplingPoint-1),_FirstAmplitudePos*_CFD);

  TLine *lnTimeCFDPos_Spline = new TLine(_FirstTimePosCFD_Spline,amplMax,_FirstTimePosCFD_Spline, amplMin);
  TLine *lnAmplCFDPos_Spline = new TLine(0.0,_FirstAmplitudePos_Spline*_CFD,_SamplingPeriod*(usbConst::usbwcNsamplingPoint-1),_FirstAmplitudePos_Spline*_CFD);
  ////////////////

  TMultiGraph *mg = new TMultiGraph("mg"," wf vs wfSp");
  mg->Add(gr1);
  mg->Add(gr2);
  mg->Draw("APL");

  lnTimeCFDPos->Draw();
  lnAmplCFDPos->Draw();
  lnTimeCFDPos_Spline->Draw();
  lnAmplCFDPos_Spline->Draw();

  //gr1->Draw("APL");

}
*/

/*
void waveform::SetIdealCosmicWFnorm(TString fileN){
  TFile *ff = new TFile(fileN.Data());
  TString histN;
  TString txtFileN;
  Int_t i = 0;
  Int_t j = 0;
  for(i = 0; i<usbConst::usbwcNchannels; i++){
    //for(i = 3; i<4; i++){
    histN = "h2waveFormNorm_ch";
    histN += i;
    cout<<"histN = "<<histN<<endl;
    TH2D *h2 = (TH2D*)ff->Get(histN.Data());
    //for(j = 0;j<usbConst::usbwcNPoint;j++){
    //_wfSpNormCosmic[j][i] = FindMaxProbAmpl(h2,_wfTimeSpAligment[j]);
    //}
    txtFileN = "../../idealWF/cosmic_ch";
    txtFileN += i;
    txtFileN += ".txt";
    SaveIdealWFintoFile(h2,txtFileN);
  }
  ff->Close();
}

void waveform::SetIdealLaserWFnorm(TString fileN){
  TFile *ff = new TFile(fileN.Data());
  //ff->ls();
  TString histN;
  Int_t i = 0;
  Int_t j = 0;
  //for(i = 0; i<usbConst::usbwcNchannels; i++){
  for(i = 3; i<4; i++){
    histN = "h2waveFormNorm_ch";
    histN += i;
    cout<<"histN = "<<histN<<endl;
    TH2D *h2 = (TH2D*)ff->Get(histN.Data());
    for(j = 0;j<usbConst::usbwcNPoint;j++){
      _wfSpNormLaser[j][i] = FindMaxProbAmpl(h2,_wfTimeSpAligment[j]);
      //cout<<_wfSpNormLaser[j][i]<<endl;
    }
  }
  ff->Close();
}


Double_t waveform::FindMaxProbAmpl(TH2D *h2, Double_t t){
  Int_t i = 0;
  Int_t j = 0;
  Double_t tMin;
  Double_t tMax;
  Int_t nBinX = h2->GetNbinsX();
  Int_t nBinY = h2->GetNbinsY();
  Double_t ampl = -999.0;
  Double_t prob = -999.0;
  for(i = 0; i<nBinX; i++){
    tMin = h2->GetXaxis()->GetBinLowEdge(i+1);
    tMax = h2->GetXaxis()->GetBinUpEdge(i+1);
    if(t>=tMin && t<tMax){
      //cout<<"tMin  =  "<<tMin<<endl
      //  <<"tMax  =  "<<tMax<<endl
      //  <<"t     =  "<<t<<endl;
      for(j = 0;j<nBinY;j++){
	if(prob<h2->GetBinContent(i+1,j+1)){
	  prob = h2->GetBinContent(i+1,j+1);
	  ampl = h2->GetYaxis()->GetBinCenter(j+1);
	}
      }
      break;
    }
  }
  return ampl;
}


void waveform::DrawIdealWFnorm(Int_t chID){

  Double_t amplMax =  1.2;
  Double_t amplMin = -0.2;

  Double_t wfSpCos[usbConst::usbwcNPoint];
  Double_t wfSpLas[usbConst::usbwcNPoint];
  
  Int_t i;

  for( i = 0;i<usbConst::usbwcNPoint;i++){
    if(chID>=0 && chID<usbConst::usbwcNchannels){
      wfSpCos[i] = _wfSpNormCosmic[i][chID];
    }
    else{
      cout<<endl<<"  ERROR --> chID is wrong"<<endl
	  <<"                  chID  = "<<chID<<endl;
      assert(0);
    }
  }

  for( i = 0;i<usbConst::usbwcNPoint;i++){
    if(chID>=0 && chID<usbConst::usbwcNchannels){
      wfSpLas[i] = _wfSpNormLaser[i][chID];
      //cout<<wfSpLas[i]<<endl;
    }
    else{
      cout<<endl<<"  ERROR --> chID is wrong"<<endl
	  <<"                  chID  = "<<chID<<endl;
      assert(0);
    }
  }
  

  TGraph *gr1 = new TGraph( usbConst::usbwcNPoint, _wfTimeSpAligment, wfSpCos); 
  gr1->SetMaximum(amplMax);
  gr1->SetMinimum(amplMin);  
  gr1->SetMarkerStyle(20);
  gr1->SetMarkerSize(1);


  TGraph *gr2 = new TGraph( usbConst::usbwcNPoint, _wfTimeSpAligment, wfSpLas); 
  gr2->SetMaximum(amplMax);
  gr2->SetMinimum(amplMin);  
  gr2->SetMarkerStyle(20);
  gr2->SetMarkerSize(1);
  gr2->SetMarkerColor(kRed);

  TMultiGraph *mg = new TMultiGraph("mg"," wf vs wfSp");
  mg->Add(gr1);
  mg->Add(gr2);
  mg->SetMaximum(amplMax);
  mg->SetMinimum(amplMin);  
  mg->Draw("APL");

  //gr2->Draw("APL");
}

void waveform::SaveIdealWFintoFile(TH2D* h2, TString txtFileN){
  ofstream dataTxtFile;
  dataTxtFile.open(txtFileN.Data()); 
  assert(dataTxtFile.is_open());

  TString grName;


  Int_t i = 0;
  Int_t j = 0;
  Int_t k = 0;

  const Int_t maxNP = 1000;
  Double_t xx[maxNP];
  Double_t yy[maxNP];

  for(i = 0;i<maxNP;i++){
    xx[i] = -999.0;
    yy[i] = -999.0;  
  }


  Double_t tMin = 13.0;
  Double_t tMax = 14.5;
  Double_t t = -999.0;
  Int_t nBinX = h2->GetNbinsX();
  Int_t nBinY = h2->GetNbinsY();
  Double_t ampl = -999.0;
  Double_t prob = -999.0;
  for(i = 0; i<nBinX; i++){
    t = h2->GetXaxis()->GetBinCenter(i+1);
    if(t>=tMin && t<tMax){
      //cout<<"tMin  =  "<<tMin<<endl
      //  <<"tMax  =  "<<tMax<<endl
      //  <<"t     =  "<<t<<endl;
      for(j = 0;j<nBinY;j++){
	if(prob<h2->GetBinContent(i+1,j+1)){
	  prob = h2->GetBinContent(i+1,j+1);
	  ampl = h2->GetYaxis()->GetBinCenter(j+1);
	}
      }
      dataTxtFile<<setw(5)<<k
		 <<setw(15)<<t
		 <<setw(15)<<ampl<<endl;
      xx[k] = t;
      yy[k] = ampl;
      k++;

      prob = -999.0;
    }
  }
  dataTxtFile.close();


  TGraph *gr1 = new TGraph(k,xx,yy);
  grName = txtFileN;
  grName += ".C";
  gr1->SaveAs(grName.Data());
}

*/
