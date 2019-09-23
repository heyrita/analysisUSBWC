{
  TRandom3 *rnd = new TRandom3(12431231);
  gROOT->LoadMacro("./src/wfSim.cpp");
  gROOT->LoadMacro("./src/waveform.cpp");
  wfSim *wfS = new wfSim( rnd, 256, 0.3125, 1.2);
  //wfS->genGaussWF( 20.0, 2.0);
  wfS->genNinoWF( 50.0, 0.22, 16.0);
  //wfS->Draw();
  waveform *wf = new waveform( wfS->getWFT(), wfS->getWFA(), 256, 0);
  wf->findParametersOftheWaveform();
  wf->Draw(0.5);
  //wf->PrintWaveFormInfo();
}
