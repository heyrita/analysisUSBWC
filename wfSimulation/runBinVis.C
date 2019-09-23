{
  gROOT->LoadMacro("./wfData/run1vis.C");
  gROOT->LoadMacro("./src/waveformSimple.cpp");
  run1vis t;
  //t.Loop();
  Int_t evID;
  Int_t chID;
  TString key_lineOn;
  TString key_printInfoTrue;
  evID = 1;
  chID = 0;
  key_lineOn = "yes";
  key_printInfoTrue = "printInfoTrue";
  //t.DrawOneChannel(evID, chID);
  t.DrawOneChannel_waveformSimple(evID, chID, key_lineOn, key_printInfoTrue);
}
