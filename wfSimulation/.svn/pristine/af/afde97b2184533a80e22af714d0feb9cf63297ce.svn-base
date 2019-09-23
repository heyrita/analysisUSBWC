#ifndef waveform_h
#define waveform_h

//root
#include <TROOT.h>
#include <TMath.h>

class TStirng;
class TCanvas;
class TH2D;
class TString;

using namespace std;

class waveform{
public :
  waveform(const Double_t *wfT,const Double_t *wfA, Int_t nPoints);
  waveform(const Double_t *wfT, const Double_t *wfA, Int_t nPoints, Int_t autoInv);
  waveform(const Double_t *wfT,const Double_t *wfA, Int_t nPoints, Int_t nSplinePoints);
  ~waveform();

  Double_t *_wfApublic;//Amplitude of the wavform with substructed baseline
  Double_t *_wfTpublic;//time of the points

  void makeSmooth(Int_t ponitsBeforeAndAfter);
  Double_t findAmplAverage(Int_t iStart,Int_t iStop);

  void findParametersOftheWaveform();
  void findParametersOftheWaveforShortList();
  void PrintWaveFormInfo();
  void Draw(Double_t cfdVal);
  void Draw15(Double_t level);
  void DrawCharge( Double_t cfdVal, Double_t dletaLeft, Double_t dletaRight);
  void DrawSmooth(Int_t ponitsBeforeAndAfter,Int_t numberOfSmoothIterations);

  void SetNpointBaseLine(Int_t nPointBaisLine){
    _nPointBaisLine = nPointBaisLine;
  }

  void SetSignalThreshold(Double_t SignalThreshold){_SignalThreshold = SignalThreshold;}
  void SetCrossTalkThreshold(Double_t crossTalkThreshold){_crossTalkThreshold = crossTalkThreshold;}
  void SetMpFraction(Double_t mpFraction){_mpFraction = mpFraction;}
  
  Int_t GetnPoints(){return _nPoints;}

  //LB 05.05.2014  
  //void findChargeInWindow(Double_t tLeft,Double_t tRight);
  //LB 10.05.2014
  Double_t findChargeInWindow(Double_t tLeft,Double_t tRight);

  Double_t findFirstTimePosCFD(Double_t cfdVal);
  Double_t findFirstTimePosConstThreas(Double_t vVal);
  //LB 11.02.2011
  //Fall edge was used for the first time 
  //all most all other function will use first positive edge
  Double_t findFirstTimePosAtLevelFallEdge(Double_t level);
  Double_t findFirstTimePosAtLevelRiseEdge(Double_t level);
  Double_t findFirstTimeNegCFD(Double_t cfdVal);
  Double_t findTimeLeftCFD( Int_t pointID, Double_t cfdVal);
  Double_t findArea(Double_t timeMin, Double_t timeMax);
  bool makeAlignmentAndNormalisation( Double_t *wfNorAligA, Double_t *wfNorAligT, Int_t pointID, Double_t t0, Double_t CFDRation);
  void findSamplingPeriod();
  Int_t findNumberOfPosIntersectionAtLevel(Double_t level);
  Int_t findNumberOfNegIntersectionAtLevel(Double_t level);

  inline Double_t GetMaxBaseLineAmpl(){return _MaxBaseLineAmpl;}
  inline Double_t GetMinBaseLineAmpl(){return _MinBaseLineAmpl;}
  inline Double_t GetSamplingPeriod(){return _SamplingPeriod;}
  inline Double_t GetSamplingPeriodStanDev(){return _SamplingPeriodStanDev;}
  inline Double_t GetBaseLineAmpl(){return _BaseLineAmpl;}      //baseline amplitude
  inline Int_t GetFirstAmplitudePosID(){return _FirstAmplitudePosID;}
  inline Double_t GetFirstAmplitudePos(){return _FirstAmplitudePos;} //point ID with first positive Amplitude of the waveform 
  inline Double_t GetFirstAmplitudeNeg(){return _FirstAmplitudeNeg;} //point ID with sirst negative Amplitude of the waveform
  inline Double_t GetChargeTOT(){return _chargeTOT;} //total charge of the waveform
  inline Double_t GetChargeVTTOT(){return _chargeTOTVT;} //total charge (VT vns) of the waveform
  inline Double_t GetChargeFirst(){return _chargeFirst;} // charge of first peak of the waveform
  inline Double_t GetChargeTOT_p(){return _chargeTOT_p;} //total charge of the waveform takein into account point bigger then baseline
  inline Double_t GetChargeTOT_m(){return _chargeTOT_m;} //total charge of the waveform takein into account point less then baseline
  inline Double_t GetChargeTOT_pm(){return (_chargeTOT_p - _chargeTOT_m);}

  inline Double_t GetChargeInWindow(){return _chargeInWindow;}
  inline Double_t GetChargeOutWindow(){return _chargeOutWindow;}

 
  inline Int_t GetWfID(){return _wfID;}
  inline Double_t GetFirstRiseTimePos(){return _FirstRiseTimePos;}
  inline Double_t GetFirstFallTimePos(){return _FirstFallTimePos;}
  inline Int_t GetTotAmplitudePosID(){return _TotAmplitudePosID;}  
  inline Double_t GetTotAmplitudePos(){return _TotAmplitudePos;}
  inline Double_t GetTotTimeAmplitudePos(){return _TotTimeAmplitudePos;}
  inline Int_t GetTotAmplitudeNegID(){return _TotAmplitudeNegID;}  
  inline Double_t GetTotAmplitudeNeg(){return _TotAmplitudeNeg;}
  inline Double_t GetTotTimeAmplitudeNeg(){return _TotTimeAmplitudeNeg;}
  inline Double_t GetFirstWidthTimePos(){return _FirstWidthTimePos;}
  inline Double_t GetFirstTimeBasis(){return _FirstTimeBasis;}
  inline Double_t GetdTimeFirtsAmplPosChangLeftDer() {return _dTimeFirtsAmplPosChangLeftDer;}
  inline Double_t GetdTimeFirtsAmplPosFirstTimeAmplNeg() {return _dTimeFirtsAmplPosFirstTimeAmplNeg;}
  inline Double_t GetMaxAmplBeforeFirstChangDerAmpl(){return _MaxAmplBeforeFirstChangDerAmpl;}
  inline Double_t GetMinAmplBeforeFirstChangDerAmpl(){return _MinAmplBeforeFirstChangDerAmpl;}

  inline Double_t GetSignalThreshold(){return _SignalThreshold;}
  inline Double_t GetCrossTalkThreshold(){return _crossTalkThreshold;}
  inline Double_t GetMpFraction(){return _mpFraction;}
  
  //public :
private : 
  Int_t _nPoints;

  Double_t *_wfA;    //Amplitude of the wavform with substructed baseline
  Double_t *_wfT;    //time of the points

  Int_t    _nPointBaisLine;    //number of points used to calculate baseline. 
  Double_t _SignalThreshold;   //signal threashold, mV (when this treashold Xing we consider it as a signal)
  Double_t _crossTalkThreshold;//crossTalk threasold, mV (when this treashold Xing we consider it as a crostalck)
  Double_t _mpFraction;        //fraction of the _FirstAmplitudePos to remove MP signals 

  Double_t _BaseLineAmpl;          //baseline amplitude
  Double_t _BaseLineAmplStanDev;   //Standard deviation of the baseline amplitude  
  Double_t _MaxBaseLineAmpl;       //Max baseline amplitude
  Double_t _MinBaseLineAmpl;       //Min baseline amplitude

  Double_t _SamplingPeriod;        //average time between two point (ns)
  Double_t _SamplingPeriodStanDev; //Standard deviation time between two point (ns)

  Int_t    _TotAmplitudePosID;   //point ID with maximum amplitude
  Double_t _TotAmplitudePos;     //maximum amplitude
  Double_t _TotTimeAmplitudePos; //time if the point with maximum amplitude

  Int_t    _TotAmplitudeNegID;   //point ID with min amplitude
  Double_t _TotAmplitudeNeg;     //min amplitude
  Double_t _TotTimeAmplitudeNeg; //time if the point with min amplitude

  Int_t    _FirstAmplitudePosID;   //point ID with first positive Amplitude of the waveform 
  Double_t _FirstAmplitudePos;     //first positive Amplitude of the waveform 
  Double_t _FirstTimeAmplitudePos; //time of the first positive Amplitude of the waveform 

  Int_t    _FirstAmplitudeNegID;   //point ID with first negative Amplitude of the waveform
  Double_t _FirstAmplitudeNeg;     //first negative Amplitude of the waveform
  Double_t _FirstTimeAmplitudeNeg; //time of the first negative Amplitude of the waveform

  Int_t    _ChangDerID;        //point ID change derivative after  First Pos amplitude (to the right from first ampl)
  Double_t _ChangDerAmpl;      //change derivative after First Pos amplitude (amplitude)
  Double_t _ChangDerTime;      //change derivative after First Pos amplitude (time)

  Int_t    _FirstChangDerID;   //point ID change derivative before First Pos amplitude (to the left from first ampl)
  Double_t _FirstChangDerAmpl; //amplitude change derivative before First Pos amplitude
  Double_t _FirstChangDerTime; //time change derivative before First Pos amplitude
  Double_t _MaxAmplBeforeFirstChangDerAmpl; //max amplitude before _FirstChangDerAmpl 
  Double_t _MinAmplBeforeFirstChangDerAmpl; //max amplitude before _FirstChangDerAmpl 

  Double_t _dTimeFirtsAmplPosChangLeftDer;    //Time difference between first positive amplitude and chanfe der left
  Double_t _dTimeFirtsAmplPosFirstTimeAmplNeg;//Time difference between first positive amplitude and first neg amplitude
  Double_t _dTimeChangRightDerChangLeftDer;   //time difference between left and right der change

  //rise time of the signal
  Double_t _FirstTimeAmplitude01Pos; //time of the first positive Amplitude of the waveform takent ot 0.1 of amplitude 
  Double_t _FirstTimeAmplitude05Pos; //time of the first positive Amplitude of the waveform takent ot 0.5 of amplitude 
  Double_t _FirstTimeAmplitude09Pos; //time of the first positive Amplitude of the waveform takent ot 0.9 of amplitude 
  Double_t _FirstTimeAmplitude00Pos; //time of the first positive Amplitude of the waveform takent ot 0.0 of amplitude 
  Double_t _FirstRiseTimePos;        //rise time of the first pos signal

  //fall time of the signal
  Double_t _FirstTimeAmplitudeFall01Pos; //time of the first positive Amplitude of the waveform takent ot 0.1 of amplitude in the fallinf age 
  Double_t _FirstTimeAmplitudeFall05Pos; //time of the first positive Amplitude of the waveform takent ot 0.5 of amplitude in the fallinf age 
  Double_t _FirstTimeAmplitudeFall09Pos; //time of the first positive Amplitude of the waveform takent ot 0.9 of amplitude in the fallinf age 
  Double_t _FirstTimeAmplitudeFall00Pos; //time of the first positive Amplitude of the waveform takent ot 0.0 of amplitude in the fallinf age 
  Double_t _FirstFallTimePos;            //fall time of the first pos signal

  //width of the signal
  Double_t _FirstWidthTimePos;

  Double_t _FirstTimeBasis;    //time difference between xing 0 from the left and right

  Double_t _chargeTOT;         //total charge of the waveform
  Double_t _chargeTOTVT;         //total charge of the waveform
  Double_t _chargeTOT_p;       //total charge of the waveform takein into account point bigger then baseline
  Double_t _chargeTOT_m;       //total charge of the waveform takein into account point less then baseline

  Double_t _chargeFirst;       //total charge of the waveform

  Double_t _chargeInWindow;// Charge in window
  Double_t _chargeOutWindow;// Charge out window

  Int_t _wfID;
  //-1 - no time measurements
  // 0 - CT
  // 1 - SP
  // 2 - MP

  //function to find time
  //make linear interpolation between points P1{x1,y1} and P2{x2,y2} and then find intersection x with y 
  Double_t linInterpol(Double_t y,Double_t y2,Double_t y1,Double_t x2,Double_t x1);

  void  InitParam();
  void  findAndSubtractBaseLineAmpl();
  void  findTotalCharge();
  void  findChargeFirst();
  void  findFirstAmplitudePos();
  void  findFirstAmplitudeNeg();
  void  findFirstRiseTimePos();
  void  findFirstFallTimePos();
  void  findFirstWidthTimePos();
  Int_t findFirstRightChangDer(Int_t FirstAmplitudePosID);
  Int_t findFirstLeftChangDer(Int_t FirstAmplitudePosID);
  void  findMinAndMaxAmplBeforeFirstChangDerAmpl();
  void  findFirstTimeBasis();
  void  findwfID();
  Int_t findClosestTimeID(Double_t t0);
  void  finddTimeFirtsAmplPosChangLeftDer();
  void  finddTimeFirtsAmplPosFirstTimeAmplNeg();
  void  findTotAmplitudePos();
  void  findTotAmplitudeNeg();

  //void waveformSpline( Int_t nSpPoints, Double_t idealWF_t[100], Double_t  idealWF_A[100]);
  //void SaveWaveForm2Cfile(TString nameff);
  //void SetIdealCosmicWFnorm(TString fileN);
  //void SetIdealLaserWFnorm(TString fileN);
  //Double_t FindMaxProbAmpl(TH2D *h2, Double_t t);
  //void DrawIdealWFnorm(Int_t chID);
  //void SaveIdealWFintoFile(TH2D* h2, TString txtFileN);
  //separate function we do not need it al the time
  //ci zminni povunni butu zminnumu funkcii
  //Double_t *_wfNorm; //Mormalized waveform 
  //Double_t *_wfTime;
  //Double_t *_wfTimeAligment;
  ////////////////////////////////
  // in the special function
  //////fill _wfTimeAligment//////
  //for(i = 0;i<usbConst::usbwcNsamplingPoint;i++){
  //_wfTimeAligment[i] = _wfTime[i] - _FirstTimeAmplitudePos + _TimeAligment;
  //_wfNorm[i] = _wf[i]/_FirstAmplitudePos;
  //}
  ////////////////////////////////
  //povunna butu zminnoju funkcii jaka bude znahodutu 4as
  //Double_t _CFD;               //Constant fraction threshold 
  //povunna butu zminnoju funkcii jaka bude znahodutu 4as
  //Double_t _TimeAligment;      //ns
  /////////////////////Splines///////////////////////////
  //Double_t _wfSp[usbConst::usbwcNPoint];//container of the wavform with substructed baseline amplitude with additional points taken from splines
  //Double_t _wfSpNorm[usbConst::usbwcNPoint];//container of the wavform with substructed baseline amplitude with additional points taken from splines
  //Double_t _wfTimeSp[usbConst::usbwcNPoint];//container of the wavform with substructed baseline amplitude with additional points taken from splines
  //Double_t _wfTimeSpAligment[usbConst::usbwcNPoint];//container of the wavform with substructed baseline amplitude with additional points taken from splines
  //Double_t _wfSpNormCosmic[usbConst::usbwcNPoint][usbConst::usbwcNchannels];
  //Double_t _wfSpNormLaser[usbConst::usbwcNPoint][usbConst::usbwcNchannels];
  //Double_t _chi2;
  //Int_t _FirstAmplitudePosID_Spline;
  //Double_t _FirstAmplitudePos_Spline;
  //Double_t _FirstTimeAmplitudePos_Spline;
  //Double_t _FirstTimePosCFD_Spline;
  //for(i = 0;i<usbConst::usbwcNPoint;i++){
  //_wfSp[i] = -999.0;
  //_wfSpNorm[i] = -999.0;
  //_wfTimeSp[i] = -999.0;
  //_wfTimeSpAligment[i] = -999.0;
  //}
  //_chi2 = -999.0;
  //_FirstAmplitudePosID_Spline = -999;
  //_FirstAmplitudePos_Spline = -999.0;
  //_FirstTimeAmplitudePos_Spline = -999.0;
  //_FirstTimePosCFD_Spline = -999.0;
};

#endif
