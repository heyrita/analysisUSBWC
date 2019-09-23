////////////////////////////////////////////////////////////////////
//                Date: Tue Sep  1 09:13:39 CEST 2015             //
//               Autor: Leonid Burmistrov                         //
// Program description: Simple analisys program of the waweforms. //
//                      Functions will be integrated into usbwc   //
//                      library.                                  //
////////////////////////////////////////////////////////////////////

//my
#include "usbwclib.hh"

//C, C++
#include <iostream>
#include <fstream>
#include <assert.h>
#include <iomanip>
#include <stdio.h>
#include <cmath>  
#include <stdlib.h>

using namespace std;

usbwclib::usbwclib(double *waveform, int nbOfPoints, double samplingPeriod){
  SoftwareMeasurementsParameters softPar;
  SoftwareMeasurements softMeas;
  Make_WaveformMeasurements(waveform, nbOfPoints, samplingPeriod,softPar,softMeas);
}

usbwclib::~usbwclib(){
}

void usbwclib::Make_WaveformMeasurements(double *waveform, int nbOfPoints, double samplingPeriod, SoftwareMeasurementsParameters softwareParameters, SoftwareMeasurements softwareMeasurements){

  int i = 0;

  cout<<"Make_WaveformMeasurements()"<<endl;
  cout<<"nbOfPoints     "<<nbOfPoints<<endl
      <<"samplingPeriod "<<samplingPeriod<<" ps"<<endl;
  
  cout<<"Parameters for the waveform analysis -> "<<endl
      <<"softwareParameters.ForceBaseline        "<<softwareParameters.ForceBaseline<<endl
      <<"softwareParameters.EnforcedBaseline     "<<softwareParameters.EnforcedBaseline<<endl
      <<"softwareParameters.nbOPointsForBaseline "<<softwareParameters.nbOPointsForBaseline<<endl
      <<"softwareParameters.FixedThresholdOrCFD  "<<softwareParameters.FixedThresholdOrCFD<<endl
      <<"softwareParameters.fixedThresholdValue  "<<softwareParameters.fixedThresholdValue<<endl
      <<"softwareParameters.CFDratio             "<<softwareParameters.CFDratio<<endl
      <<"softwareParameters.RefCellForCharge     "<<softwareParameters.RefCellForCharge<<endl
      <<"softwareParameters.ChargeLength         "<<softwareParameters.ChargeLength<<endl;
  
  double *waveformA = new double [nbOfPoints];//Amplitude of the waveform with subtraction baseline
  double *waveformT = new double [nbOfPoints];//time of the points

  for(i = 0;i<nbOfPoints;i++){
    waveformA[i] = waveform[i];
    waveformT[i] = i*samplingPeriod;
  }
  
  //Baseline subtraction
  if(softwareParameters.ForceBaseline == true){
    for(i = 0;i<nbOfPoints;i++)
      waveformA[i] = waveformA[i] - softwareParameters.EnforcedBaseline;
    softwareMeasurements.ComputedBaseline = softwareParameters.EnforcedBaseline;
  }
  else{
    double abasel = findBaseLineAmpl(softwareParameters.nbOPointsForBaseline, waveformA, nbOfPoints);
    for(i = 0;i<nbOfPoints;i++)
      waveformA[i] = waveformA[i] - abasel;
    softwareMeasurements.ComputedBaseline = abasel;
  }

  //PolarityType
  int i_amax = 0;
  softwareMeasurements.PeakPolarity = findSignalPolarity(waveformA, nbOfPoints);
  if(softwareMeasurements.PeakPolarity == Postive){
    //PeakAmplitude
    softwareMeasurements.PeakAmplitude = findMaximumAmplitude(waveformA, nbOfPoints, i_amax); //in Volts
  }
  else if(softwareMeasurements.PeakPolarity == Negative){
    softwareMeasurements.PeakAmplitude = -findMaximumAmplitude(waveformA, nbOfPoints, i_amax); //in Volts
    //softwareMeasurements.PeakAmplitude = -softwareMeasurements.PeakAmplitude;
  }
  else{
    assert(0);
  }
  //PeakTime
  softwareMeasurements.PeakTime = waveformT[i_amax];  // in ns
  
  //LeadingEdgeTime
  //TrailingEdgeTime
  double leadingEdgeTime = -999.0; 
  double trailingEdgeTime =-999.0;  
  if(softwareParameters.FixedThresholdOrCFD == CFD)
    findTimeCFD(waveformA, waveformT, nbOfPoints, i_amax, softwareParameters.CFDratio, leadingEdgeTime, trailingEdgeTime);
  else if(softwareParameters.FixedThresholdOrCFD == FixedThreshold)
    findTimeFixThreshold(waveformA, waveformT, nbOfPoints, i_amax, softwareParameters.fixedThresholdValue, leadingEdgeTime, trailingEdgeTime);
  else
    assert(0);
  softwareMeasurements.LeadingEdgeTime  = leadingEdgeTime;    // in ns
  softwareMeasurements.TrailingEdgeTime = trailingEdgeTime;   // in ns

  softwareMeasurements.Charge = findChargeInWindow(waveformA, waveformT, nbOfPoints, softwareParameters.RefCellForCharge, softwareParameters.ChargeLength);
}

double usbwclib::findBaseLineAmpl(int np, double *waveformA, int nbOfPoints){
  double abasel = 0;
  if(np <= 0 || np > nbOfPoints){
    cout<<endl<<" ERROR--> np < = 0 || np > nbOfPoints"<<endl
	<<"        np = "<<np<<endl
	<<"nbOfPoints = "<<nbOfPoints<<endl;
    assert(0);
  }
  for(int i = 0;i<np;i++)
    abasel += waveformA[i];
  return abasel/np;
}

double usbwclib::findChargeInWindow(double tLeft, double tRight, double *waveformA, double *waveformT, int nbOfPoints){
 double chargeInWindow = 0.0;
 double dT = 0.0;
 for(int i = 0;i<(nbOfPoints - 1);i++){
   dT = (waveformT[i+1] - waveformT[i]);
   if(dT<0.0){
     cout<<" ERROR -> dT<0.0, dT = (waveformT[i+1] - waveformT[i])"<<endl
	 <<"dT        = "<<dT<<endl
	 <<"waveformT[i+1] = "<<waveformT[i+1]<<endl
	 <<"waveformT[i]   = "<<waveformT[i]<<endl
	 <<"i         = "<<i<<endl;
     assert(0);
   }
   if(waveformT[i]>=tLeft && waveformT[i]<=tRight)
     chargeInWindow = chargeInWindow + waveformA[i]*dT;
 }
 return chargeInWindow;
}

PolarityType usbwclib::findSignalPolarity(double *waveformA, int nbOfPoints){
  double amin = waveformA[0];
  double amax = waveformA[0];
  int i = 0;
  for(i = 0;i<nbOfPoints;i++){
    if(waveformA[i]<amin)
      amin = waveformA[i];
    if(waveformA[i]>amax)
      amax = waveformA[i];
  }
  double aminabs = abs(amin);
  double amaxabs = abs(amax);
  if(amaxabs<aminabs){
    for(i = 0;i<nbOfPoints;i++)
      waveformA[i] = -waveformA[i];
    return Negative;
  }
  else{
    return Postive;
  }
  return Negative;
}

double usbwclib::findMaximumAmplitude(double *waveformA, int nbOfPoints, int &i_amax){
  i_amax = 0;
  double amax = waveformA[0];
  for(int i = 0;i<nbOfPoints;i++){
    if(waveformA[i]>amax){
      amax = waveformA[i];
      i_amax = i;
    }
  }
  if(i_amax<0 || i_amax>=nbOfPoints){
    cout<<endl<<" ERROR--> i_amax<0 || i_amax>=nbOfPoints"<<endl
	<<"                i_amax = "<<i_amax<<endl
	<<"            nbOfPoints = "<<nbOfPoints<<endl;
    assert(0);
  }
  return amax;
}

void usbwclib::findTimeCFD(double *waveformA, double *waveformT, int nbOfPoints, int i_amax, double CFDratio, double &leadingEdgeTime, double &trailingEdgeTime){
  if(i_amax<0 || i_amax>=nbOfPoints){
    cout<<endl<<" ERROR--> i_amax<0 || i_amax>=nbOfPoints"<<endl
	<<"                i_amax = "<<i_amax<<endl
	<<"            nbOfPoints = "<<nbOfPoints<<endl;
    assert(0);
  }
  leadingEdgeTime  = -999.0;
  trailingEdgeTime = -999.0;
  double amax = waveformA[i_amax];
  double tmax = waveformT[i_amax];
  int i = 0;
  for(i = i_amax; i>0 ; i--){
    if(waveformA[i]>=amax*CFDratio && waveformA[i-1]<amax*CFDratio){
      leadingEdgeTime = linInterpol(amax*CFDratio,
				    waveformA[i],waveformA[i-1],
				    waveformT[i],waveformT[i-1]);
      break;
    }
  }
  for(i = i_amax; i<(nbOfPoints-1); i++){
    if(waveformA[i]>=amax*CFDratio && waveformA[i+1]<amax*CFDratio){
      trailingEdgeTime = linInterpol(amax*CFDratio,
				     waveformA[i+1],waveformA[i],
				     waveformT[i+1],waveformT[i]);
      break;
    }
  }
}

void usbwclib::findTimeFixThreshold(double *waveformA, double *waveformT, int nbOfPoints, int i_amax, double fixedThresholdValue, double &leadingEdgeTime, double &trailingEdgeTime){
  if(i_amax<0 || i_amax>=nbOfPoints){
    cout<<endl<<" ERROR--> i_amax<0 || i_amax>=nbOfPoints"<<endl
	<<"                i_amax = "<<i_amax<<endl
	<<"            nbOfPoints = "<<nbOfPoints<<endl;
    assert(0);
  }
  leadingEdgeTime  = -999.0;
  trailingEdgeTime = -999.0;
  int i = 0;
  for(int i = i_amax; i>0 ; i--){
    if(waveformA[i]>=fixedThresholdValue && waveformA[i-1]<fixedThresholdValue){
      leadingEdgeTime = linInterpol(fixedThresholdValue,
				    waveformA[i],waveformA[i-1],
				    waveformT[i],waveformT[i-1]);
      break;
    }
  }
  for(i = i_amax; i<(nbOfPoints-1); i++){
    if(waveformA[i]>=fixedThresholdValue && waveformA[i+1]<fixedThresholdValue){
      trailingEdgeTime = linInterpol(fixedThresholdValue,
				     waveformA[i+1],waveformA[i],
				     waveformT[i+1],waveformT[i]);
      break;
    }
  }
}

double usbwclib::findChargeInWindow(double *waveformA, double *waveformT, int nbOfPoints, int refCellForCharge, int chargeLength){
  double chargeInWindow = 0.0;
  double dT = 0.0;
  int i_cell = 0;
  for(int i = 0;i<chargeLength;i++){
    i_cell = refCellForCharge + i;
    if(i_cell<0 || i_cell>=(nbOfPoints-1)){
      cout<<" ERROR -> i_cell<0 || i_cell>=(nbOfPoints-1)"<<endl
	  <<"i_cell           = "<<i_cell<<endl
	  <<"nbOfPoints       = "<<nbOfPoints<<endl
	  <<"refCellForCharge = "<<refCellForCharge<<endl
	  <<"chargeLength     = "<<chargeLength<<endl;
      assert(0);
    }
    dT = (waveformT[i_cell+1] - waveformT[i_cell]);
    if(dT<0.0){
      cout<<" ERROR -> dT<0.0, dT = (waveformT[i+1] - waveformT[i])"<<endl
	  <<"dT             = "<<dT<<endl
	  <<"waveformT[i+1] = "<<waveformT[i+1]<<endl
	  <<"waveformT[i]   = "<<waveformT[i]<<endl
	  <<"i              = "<<i_cell<<endl;
      assert(0);
    }
    chargeInWindow = chargeInWindow + waveformA[i_cell]*dT;
  }
  double resistance50Om = 50.0;
  chargeInWindow = chargeInWindow/resistance50Om;
  return chargeInWindow;
}

double usbwclib::linInterpol(double y, double y2, double y1, double x2, double x1){
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
