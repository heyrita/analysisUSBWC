//my
#include "src/usbwclib.hh"

//C, C++
#include <iostream>
#include <stdlib.h>
#include <assert.h>
#include <fstream>
#include <iomanip>
#include <time.h>

using namespace std;

int main(int argc, char *argv[]){

  clock_t start, finish;
  start = clock();

  if(argc == 1){
    double *waveform;
    int nbOfPoints = 1024;
    double samplingPeriod = 312.5;
    int i = 0;
    waveform = new double [nbOfPoints];
    for(i = 0;i<nbOfPoints;i++)
      waveform[i] = 2.22;
    usbwclib *c = new usbwclib(waveform,nbOfPoints,samplingPeriod);
  }
  /*
    else{
    cout<<"  wavefrom (width ampl)           "<<endl
    <<"  runID [1] = 0                   "<<endl
    <<"        [2] - file with histos    "<<endl
    <<"        [3] - number of waveforms "<<endl
    <<"        [4] - seed                "<<endl;
  }
  */

  finish = clock();

  cout<<" time of work "<<((finish - start)/CLOCKS_PER_SEC)<<endl;
  //cout<<" time of work "<<(finish - start)<<endl;


  return 0;
}
