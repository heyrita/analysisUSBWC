////////////////////////////////////////////////////////////////////
//                Date: Tue Sep  1 09:13:39 CEST 2015             //
//               Autor: Leonid Burmistrov                         //
// Program description: Simple analisys program of the waweforms. //
//                      Functions will be integrated into usbwc   //
//                      library.                                  //
////////////////////////////////////////////////////////////////////

#ifndef usbwclib_h
#define usbwclib_h

using namespace std;

enum TimeMethod {
  FixedThreshold,
  CFD
};

enum PolarityType {
  Postive,
  Negative
};

struct SoftwareMeasurementsParameters {
  bool ForceBaseline;
  double EnforcedBaseline;     // in Volts
  int nbOPointsForBaseline;
  TimeMethod FixedThresholdOrCFD;
  double fixedThresholdValue;  // in Volts,
  double CFDratio;
  int RefCellForCharge;        // in samples
  int ChargeLength;            // in samples
  SoftwareMeasurementsParameters() :
    ForceBaseline(0),
    EnforcedBaseline(0),    // in Volts
    nbOPointsForBaseline(20),
    FixedThresholdOrCFD(CFD),
    fixedThresholdValue(0), // in Volts,
    CFDratio(0.5),
    RefCellForCharge(0),    // in samples
    ChargeLength(1024)      // in samples
  {;}
  ~SoftwareMeasurementsParameters() 
  {;}
};

struct SoftwareMeasurements {
  double ComputedBaseline;
  PolarityType PeakPolarity;
  double PeakAmplitude;      //in Volts
  double PeakTime;           // in ns
  double LeadingEdgeTime;    // in ns
  double TrailingEdgeTime;   // in ns
  double Charge;             // in pC
  SoftwareMeasurements() :
    ComputedBaseline(-999.0),
    PeakPolarity(Postive),
    PeakAmplitude(-999.0),    //in Volts
    PeakTime(-999.0),         // in ns
    LeadingEdgeTime(-999.0),  // in ns
    TrailingEdgeTime(-999.0), // in ns
    Charge(-999.0)            // in pC
  {;}
  ~SoftwareMeasurements()
  {;}
};

class usbwclib{
public :
  //            volts                                 ps
  usbwclib(double *waveform, int nbOfPoints, double samplingPeriod);
  ~usbwclib();
  
  void Make_WaveformMeasurements(double *waveform, int nbOfPoints, double samplingPeriod, SoftwareMeasurementsParameters softwareParameters, SoftwareMeasurements softwareMeasurements);
  double findBaseLineAmpl(int np, double *waveformA, int nbOfPoints);
  double findChargeInWindow(double tLeft, double tRight, double *waveformA, double *waveformT, int nbOfPoints);
  PolarityType findSignalPolarity(double *waveformA, int nbOfPoints);
  double findMaximumAmplitude(double *waveformA, int nbOfPoints, int &i_amax);
  void findTimeCFD(double *waveformA, double *waveformT, int nbOfPoints, int i_amax, double CFDratio, double &leadingEdgeTime, double &trailingEdgeTime);
  void findTimeFixThreshold(double *waveformA, double *waveformT, int nbOfPoints, int i_amax, double fixedThresholdValue, double &leadingEdgeTime, double &trailingEdgeTime);
  double findChargeInWindow(double *waveformA, double *waveformT, int nbOfPoints, int refCellForCharge, int chargeLength);
  double linInterpol(double y, double y2, double y1, double x2, double x1);

};

#endif
