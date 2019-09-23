#ifndef waveformSimple_h
#define waveformSimple_h

//root
#include <TROOT.h>
#include <TMath.h>

class TStirng;
class TCanvas;
class TH2D;
class TString;
class TSpline3;

using namespace std;

class waveformSimple{
public :
  waveformSimple(const Double_t *wfT, const Double_t *wfA, Int_t nPoints, Int_t nPointsBaseLine);
  waveformSimple(const Double_t *wfT, const Double_t *wfA, Int_t nPoints, Int_t nPointsBaseLine, Int_t nSplinePoints);
  ~waveformSimple();

  Double_t findBaseLineAmpl(Int_t np);
  Double_t findBaseLineAmpl(Int_t np, const Double_t *wfA);
  Double_t findChargeInWindow(Double_t tLeft,Double_t tRight);
  void findPar(Double_t tLeft,Double_t tRight,Double_t  vVal);
  Double_t linInterpol(Double_t y,Double_t y2,Double_t y1,Double_t x2,Double_t x1);
  void printPar();
  void Draw(TString key_lineOn);

  void normalizeWF(double normVal);
  void alignInTimeOfMaximum(double tVal);

  double GetAbasel(){return _abasel;}
  double GetSigWidth(){return _sigWidth;}
  double GetSigRise(){return _sigRise;}
  double GetSigFall(){return _sigFall;}
  double GetTimeCFD(){return _timeCFD;}
  double GetChargeInW(){return _chargeInW;}
  double GetTimeConstThreas(){return _timeConstThreas;}
  double GetAmin(){return _amin;}
  double GetAmax(){return _amax;}
  double GetTmin(){return _tmin;}
  double GetTmax(){return _tmax;}
  Int_t Get_nPoints(){return _nPoints;}
  Double_t Get_wfAi(int i);
  Double_t Get_wfTi(int i);

private : 
  Int_t _nPoints;
  Int_t _nPointBaisLine;
  Int_t _nSplinePoints;

  Double_t *_wfA;    //Amplitude of the wavform with substructed baseline
  Double_t *_wfT;    //time of the points

  double _abasel;
  double _sigWidth;
  double _sigRise;
  double _sigFall;
  double _timeCFD;
  double _chargeInW;
  double _timeConstThreas;
  double _amin;
  double _amax;
  double _tmin;
  double _tmax;
  double _timeCFD_fall;
  double _timeCFD01;
  double _timeCFD01_fall;
  double _timeCFD09;
  double _timeCFD09_fall;
  double _chargeWindow_tLeft;
  double _chargeWindow_tRight;
  int _i_amax;
  double _vVal;
};

#endif
