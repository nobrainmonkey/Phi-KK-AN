// Run from the repository root after building with cons in the STAR environment.
void RunMCTruthTests(const char *mode="io", const char *input="", const char *output="mc_io.root", const char *legacy="", Long64_t expectedEntries=-1) {
  gROOT->Macro("loadMuDst.C");
  if(gSystem->Load("StUPCFilterMaker.so") < 0) { gSystem->Exit(1); return; }
  if(expectedEntries >= 0) {
    TFile checkFile(output);
    TTree *checkTree = dynamic_cast<TTree*>(checkFile.Get("mUPCTree"));
    if(checkFile.IsZombie() || !checkTree || checkTree->GetEntries()!=expectedEntries) {
      std::cerr << "Unexpected UPCDst event count; expected " << expectedEntries << std::endl;
      gSystem->Exit(1); return;
    }
    checkFile.Close();
  }
  TString star = gSystem->Getenv("STAR");
  gSystem->AddIncludePath(Form("-I./include -I%s/StRoot -I%s/StRoot/StarClassLibrary -I%s/.sl73_x8664_gcc485/include", star.Data(), star.Data(), star.Data()));
  TString compile = gSystem->GetMakeSharedLib();
  compile.ReplaceAll("g++", "g++ -std=c++11");
  gSystem->SetMakeSharedLib(compile);
  if(gROOT->LoadMacro("tests/TestMCTruth.C+") < 0) { gSystem->Exit(1); return; }
  Int_t error = 0;
  Long_t result = gROOT->ProcessLine(Form("TestMCTruth(\"%s\",\"%s\",\"%s\",\"%s\");", mode, input, output, legacy), &error);
  gSystem->Exit(error || result ? 1 : 0);
}
