#include "../events/VzFit.h"
#include "TFile.h"
#include "TH1D.h"
#include "TLine.h"
#include "TNamed.h"
#include "TParameter.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"
#include <cmath>
#include <cstdio>
#include <fstream>

// Run from plots/signal: root -l -b -q 'coherent_rap_range.C'
// Mass spectra in four signed-rapidity windows, before the final track-pT cut.
// All other nominal cuts are applied; the stored preselection is unchanged.
// Windows: [-0.5,-0.25), [-0.25,0.25], (0.25,0.5], (0.5,1].
void coherent_rap_range(TString input_file = "../../data/phi-kk-mudst.root") {
  double pairPtMax = 0.1; // GeV/c
  double chi2Max = 15;
  double vzMax = 100, vrMax = 5; // cm
  int nHitsMin = 20;
  double massMin = 1.00, massMax = 1.04; // GeV/c^2
  TString outputDir = "coherent_rap_range/";

  TFile *file = TFile::Open(input_file, "READ");
  if (!file || file->IsZombie()) {
    printf("Cannot open input file\n");
    return;
  }
  TTree *tree = (TTree *)file->Get("uls_data");
  if (!tree) {
    printf("Cannot find uls_data\n");
    file->Close();
    return;
  }
  TFile *calibration =
      TFile::Open("../track/PID/nsigmakaon_mean_vs_p.root", "READ");
  if (!calibration || calibration->IsZombie()) {
    printf("Cannot open kaon PID calibration\n");
    file->Close();
    return;
  }
  TH1D *shift = (TH1D *)calibration->Get("hMeanNSigmaKaon");
  if (!shift) {
    printf("Cannot find hMeanNSigmaKaon\n");
    calibration->Close();
    file->Close();
    return;
  }
  double correctionMin = shift->GetBinCenter(1);
  double correctionMax = shift->GetBinCenter(shift->GetNbinsX());

  UShort_t bbcEast = 0, bbcWest = 0, nhits1 = 0, nhits2 = 0;
  UInt_t nPrimary = 0;
  double vx = 0, vy = 0, vz = 0;
  double mass = 0, pairPt = 0, pairY = 0;
  double pt1 = 0, pt2 = 0, eta1 = 0, eta2 = 0, nsk1 = 0, nsk2 = 0;
  tree->SetBranchAddress("uls_pair_BBC_Veto_East", &bbcEast);
  tree->SetBranchAddress("uls_pair_BBC_Veto_West", &bbcWest);
  tree->SetBranchAddress("uls_pair_event_n_primary_tracks", &nPrimary);
  tree->SetBranchAddress("uls_pair_event_vx", &vx);
  tree->SetBranchAddress("uls_pair_event_vy", &vy);
  tree->SetBranchAddress("uls_pair_event_vz", &vz);
  tree->SetBranchAddress("uls_phi_mass", &mass);
  tree->SetBranchAddress("uls_phi_pt", &pairPt);
  tree->SetBranchAddress("uls_phi_rapidity", &pairY);
  tree->SetBranchAddress("uls_kaon1_pt", &pt1);
  tree->SetBranchAddress("uls_kaon2_pt", &pt2);
  tree->SetBranchAddress("uls_kaon1_eta", &eta1);
  tree->SetBranchAddress("uls_kaon2_eta", &eta2);
  tree->SetBranchAddress("uls_kaon1_nhits", &nhits1);
  tree->SetBranchAddress("uls_kaon2_nhits", &nhits2);
  tree->SetBranchAddress("uls_kaon1_nsigma_kaon", &nsk1);
  tree->SetBranchAddress("uls_kaon2_nsigma_kaon", &nsk2);

  gROOT->cd();
  TString labels[4] = {"-0.50 #leq y_{KK} < -0.25",
                       "-0.25 #leq y_{KK} #leq 0.25", "0.25 < y_{KK} #leq 0.50",
                       "0.50 < y_{KK} #leq 1.00"};
  TString intervals[4] = {"[-0.5,-0.25)", "[-0.25,0.25]", "(0.25,0.5]",
                          "(0.5,1]"};
  TH1D *h_mass[4];
  int total[4] = {0}, window[4] = {0};
  int outsideScan = 0, outsideScanWindow = 0;
  for (int i = 0; i < 4; ++i) {
    h_mass[i] =
        new TH1D(Form("mass_rapidity_bin%d", i),
                 ";M_{KK} [GeV/c^{2}];Pairs / 4 MeV/c^{2}", 55, 0.98, 1.20);
  }
  TH1D *h_y_all = new TH1D("pair_y_all_masses", ";y_{KK};Pairs", 24, -1.5, 1.5);
  TH1D *h_y_window =
      new TH1D("pair_y_mass_window", ";y_{KK};Pairs", 24, -1.5, 1.5);
  TH1D *h_y_outside_window =
      new TH1D("pair_y_outside_mass_window", ";y_{KK};Pairs", 24, -1.5, 1.5);

  for (Long64_t entry = 0; entry < tree->GetEntries(); ++entry) {
    if (tree->GetEntry(entry) <= 0) {
      printf("Failed to read entry %lld\n", entry);
      calibration->Close();
      file->Close();
      return;
    }
    if (bbcEast != 0 || bbcWest != 0)
      continue;
    if (nPrimary != 2)
      continue;
    if (std::abs(vz) > vzMax || std::hypot(vx, vy) > vrMax)
      continue;
    if (nhits1 < nHitsMin || nhits2 < nHitsMin)
      continue;
    if (pairPt > pairPtMax)
      continue;

    double correctionP1 = pt1 * std::cosh(eta1);
    double correctionP2 = pt2 * std::cosh(eta2);
    if (correctionP1 < correctionMin)
      correctionP1 = correctionMin;
    if (correctionP1 > correctionMax)
      correctionP1 = correctionMax;
    if (correctionP2 < correctionMin)
      correctionP2 = correctionMin;
    if (correctionP2 > correctionMax)
      correctionP2 = correctionMax;
    double nskCorr1 = nsk1 - shift->Interpolate(correctionP1);
    double nskCorr2 = nsk2 - shift->Interpolate(correctionP2);
    double chi2 = nskCorr1 * nskCorr1 + nskCorr2 * nskCorr2;
    if (chi2 > chi2Max)
      continue;

    // Exact signed-y boundaries requested for the four panels.
    bool inMassWindow = mass > massMin && mass < massMax;
    h_y_all->Fill(pairY);
    if (inMassWindow)
      h_y_window->Fill(pairY);
    else
      h_y_outside_window->Fill(pairY);

    int bin = -1;
    if (pairY >= -0.5 && pairY < -0.25)
      bin = 0;
    else if (pairY >= -0.25 && pairY <= 0.25)
      bin = 1;
    else if (pairY > 0.25 && pairY <= 0.5)
      bin = 2;
    else if (pairY > 0.5 && pairY <= 1.0)
      bin = 3;
    if (bin < 0) {
      outsideScan++;
      if (inMassWindow)
        outsideScanWindow++;
      continue;
    }
    h_mass[bin]->Fill(mass);
    total[bin]++;
    if (inMassWindow)
      window[bin]++;
  }
  calibration->Close();
  file->Close();

  gSystem->mkdir(outputDir, true);
  std::ofstream table((outputDir + "counts.csv").Data());
  table << "signed_y_interval,selected_pairs,mass_window_pairs,outside_window_"
           "pairs\n";
  for (int i = 0; i < 4; ++i) {
    printf("%s: %d selected, %d window, %d outside window\n",
           intervals[i].Data(), total[i], window[i], total[i] - window[i]);
    table << "\"" << intervals[i] << "\"," << total[i] << "," << window[i]
          << "," << total[i] - window[i] << "\n";
  }
  table.close();
  printf("Outside the four windows: %d selected, %d window\n", outsideScan,
         outsideScanWindow);

  // Keep cloned plot frames out of ROOT directories.
  bool addDirectory = TH1::AddDirectoryStatus();
  TH1::AddDirectory(false);

  // Same plot style and linear axes as the nominal extraction.
  gStyle->SetPaperSize(26, 19);
  PlotStyle::Options style;
  style.title = "Run-19 Au+Au";
  style.subtitle = "#sqrt{s_{NN}} = 200 GeV";
  style.note = "ZDCMB; no final track-p_{T} cut";
  style.legendX1 = 0.56;
  style.headroom =
      1.65; // Space for the rapidity and count labels above the errors.
  style.xlabel = "M_{KK} [GeV/c^{2}]";
  style.ylabel = "Pairs / 4 MeV/c^{2}";
  style.XRange(0.98, 1.20);

  TCanvas *summary = new TCanvas("rapidity_windows", "", 1500, 1100);
  summary->Divide(2, 2, 0.002, 0.002);
  for (int i = 0; i < 4; ++i) {
    PlotStyle::Options opt = style;
    TCanvas *c = PlotStyle::Draw({{h_mass[i], "ULS", "E1"}}, opt);
    DrawLines({labels[i], Form("Window: %d; outside: %d", window[i],
                               total[i] - window[i])},
              opt.annotationX, BelowNote(opt), 0.032);
    c->Update();
    TLine lower(massMin, 0, massMin, 0.45 * c->GetUymax());
    lower.SetLineStyle(2);
    lower.SetLineColor(kRed + 1);
    lower.DrawClone();
    TLine upper(massMax, 0, massMax, 0.45 * c->GetUymax());
    upper.SetLineStyle(2);
    upper.SetLineColor(kRed + 1);
    upper.DrawClone();
    c->Print(outputDir + Form("mass_rapidity_bin%d.pdf", i));
    summary->cd(i + 1);
    c->DrawClonePad();
    delete c;
  }
  summary->Print(outputDir + "mass_rapidity_windows.pdf");
  summary->Print(outputDir + "mass_rapidity_windows.png");
  {
    PlotStyle::Options opt = style;
    opt.XRange(-1.5, 1.5);
    opt.xlabel = "y_{KK}";
    opt.ylabel = "Pairs";
    TCanvas *c = PlotStyle::Draw({{h_y_window, "Mass window", "E1"},
                                  {h_y_outside_window, "Outside window", "E1"}},
                                 opt);
    DrawLines({"No pair-rapidity or final track-p_{T} cut"}, opt.annotationX,
              BelowNote(opt), 0.032);
    c->Update();
    TLine lower(-0.5, 0, -0.5, 0.45 * c->GetUymax());
    lower.SetLineStyle(2);
    lower.SetLineColor(kRed + 1);
    lower.DrawClone();
    TLine upper(0.5, 0, 0.5, 0.45 * c->GetUymax());
    upper.SetLineStyle(2);
    upper.SetLineColor(kRed + 1);
    upper.DrawClone();
    c->Print(outputDir + "pair_rapidity.pdf");
    c->Print(outputDir + "pair_rapidity.png");
    delete c;
  }

  TH1::AddDirectory(addDirectory);
  TFile *output =
      TFile::Open(outputDir + "coherent_rap_range.root", "RECREATE");
  if (!output || output->IsZombie()) {
    printf("Cannot create coherent_rap_range.root\n");
    return;
  }
  for (int i = 0; i < 4; ++i) {
    h_mass[i]->Write();
    TParameter<int>(Form("rapidity_bin%d_window_count", i), window[i]).Write();
    TNamed(Form("rapidity_bin%d_interval", i), intervals[i].Data()).Write();
  }
  h_y_all->Write();
  h_y_window->Write();
  h_y_outside_window->Write();
  TParameter<int>("outside_scan_selected", outsideScan).Write();
  TParameter<int>("outside_scan_window", outsideScanWindow).Write();
  TParameter<int>("final_track_pt_cut_applied", 0).Write();
  TParameter<double>("pair_pt_max_GeV", pairPtMax).Write();
  TParameter<double>("mass_window_low_GeV", massMin).Write();
  TParameter<double>("mass_window_high_GeV", massMax).Write();
  TNamed("selection",
         "BBC veto; two primary tracks; |Vz|<=100 cm; Vr<=5 cm; both "
         "nHits>=20; pair pT<=0.1 GeV/c; corrected chi2KK<=15; "
         "no final track pT cut; no pair rapidity cut; stored preselection "
         "unchanged.")
      .Write();
  TNamed("rapidity_bins",
         "Signed y: [-0.5,-0.25), [-0.25,0.25], (0.25,0.5], (0.5,1]. Raw "
         "counts, not divided by rapidity-bin width.")
      .Write();
  output->Close();
  delete summary;
}
