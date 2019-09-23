#ifndef wfSim_h
#define wfSim_h

//my

//C, C++
#include <string>
#include <iostream>
#include <math.h>

//root
#include <TString.h>
#include <TMath.h>

using namespace std;

class TRandom3;

class wfSim {

public :
  wfSim( TRandom3 *rnd, int nDigitPoint, double dTimeDigit, double digitTime0);
  ~wfSim();

public :
  void genTriangleWF( double timeTrue, double basis);
  void genGaussWF( double timeTrue, double sigma);
  void genNinoWF( double timeTrue, double amplitude, double width);
  void genMCPPMT_SP_WF( double timeTrue);
  void genMCPPMT_CT_WF( double timeTrue);
  void genMCPPMT_CS_WF( double timeTrue);
  void generateNoiseGauss(double noiseRMS);
  void Draw();
  void CleanWf();
  void SetIdealShape_SP(int nPiontsMcpPM, double *xxMcpPM, double *yyMcpPM, double timeMax);
  void SetIdealShape_CT(int nPiontsMcpPM, double *xxMcpPM, double *yyMcpPM, double timeMax);
  void SetIdealShape_CS(int nPiontsMcpPM, double *xxMcpPM, double *yyMcpPM, double timeMax);

  const double *getWFA()const {return _wfA;};
  const double *getWFT()const {return _wfT;};

  void SetDigitTime0(double digitTime0){_digitTime0 = digitTime0;}
  double GetDigitTime0(){return _digitTime0;}  
  
private:
  //const int _maxDigitPoint = ;
  int _nDigitPoint;
  double _dTimeDigit;

  TRandom3 *_rnd;

  double *_wfT;
  double *_wfA;

  double _digitTime0;

  double _constDelay;//ns

  double triangleShape( double x, double timeTrue, double basis);
  //double genAmplitudeDistPMT();
  double genAmplitudeDistGauss(double meanA, double rmsA);
  double genAmplitudeDistMCPMPT_SP();
  double genAmplitudeDistMCPMPT_CT();
  double genAmplitudeDistMCPMPT_CS();

  //MCP pmt
  double mcpPMShapeSP( double x, double timeTrue);
  double mcpPMShapeCT( double x, double timeTrue);
  double mcpPMShapeCS( double x, double timeTrue);

  double *_xxMcpPM_SP;
  double *_yyMcpPM_SP;
  int _nPiontsMcpPM_SP;
  double _timeMax_SP;

  double *_xxMcpPM_CT;
  double *_yyMcpPM_CT;
  int _nPiontsMcpPM_CT;
  double _timeMax_CT;

  double *_xxMcpPM_CS;
  double *_yyMcpPM_CS;
  int _nPiontsMcpPM_CS;
  double _timeMax_CS;

  double linInterpol( double x, double y2, double y1, double x2, double x1);
  int getLinearPointID( double x, double xxMin, double xxMax, int nPoints);
};

#endif
