#include "../../events/VzFit.h" // Existing PlotStyle and 2D plotting style.
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
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

// Run from plots/signal/incoherent: root -l -b -q 'incoherent_yield.C'
// Nominal raw yield: ULS count in 1.00 < M_KK < 1.04 GeV/c^2.
// Saves the selection plots, track/pair plots, and histograms in this
// directory.
void incoherent_yield(TString input_file = "../../../data/phi-kk-mudst.root") {
  // Nominal cuts.
  double massMin = 1.00, massMax = 1.04;   // GeV/c^2; strict window edges
  double pairPtMin = 0.3, pairPtMax = 0.9; // GeV/c
  double pairYMax = 0.5;
  double chi2Max = 10; // corrected kaon-pair chi2
  double chi2PionMin = 15;
  // No nominal electron rejection; its effect is studied in pid_study/.
  int nHitsMin = 20;
  double vzMax = 100, vrMax = 5; // cm
  TString outputDir = "./";

  // Input tree and kaon PID calibration.
  TFile *file = TFile::Open(input_file, "READ");
  if (!file || file->IsZombie()) {
    printf("Cannot open input file: %s\n", input_file.Data());
    return;
  }
  TTree *tree = (TTree *)file->Get("uls_data");
  if (!tree) {
    printf("Cannot find uls_data\n");
    file->Close();
    return;
  }
  TFile *calibration =
      TFile::Open("../../track/PID/nsigmakaon_mean_vs_p.root", "READ");
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

  // Variables read from each tree entry.
  UShort_t bbcEast = 0;
  UShort_t bbcWest = 0;
  UInt_t nPrimary = 0;
  Double_t vx = 0;
  Double_t vy = 0;
  Double_t vz = 0;
  Double_t mass = 0;
  Double_t pairPt = 0;
  Double_t pairY = 0;
  Double_t pt1 = 0;
  Double_t pt2 = 0;
  Double_t eta1 = 0;
  Double_t eta2 = 0;
  UShort_t nhits1 = 0;
  UShort_t nhits2 = 0;
  Double_t nsk1 = 0;
  Double_t nsk2 = 0;
  Double_t pairPhi = 0;
  Double_t phi1 = 0;
  Double_t phi2 = 0;
  UShort_t nd1 = 0;
  UShort_t nd2 = 0;
  Double_t nsPi1 = 0;
  Double_t nsPi2 = 0;
  Double_t nsE1 = 0, nsE2 = 0;
  Double_t nsP1 = 0, nsP2 = 0;
  Double_t chi2Pion = 0;
  Double_t de1 = 0;
  Double_t de2 = 0;
  Double_t dx1 = 0;
  Double_t dy1 = 0;
  Double_t dz1 = 0;
  Double_t dx2 = 0;
  Double_t dy2 = 0;
  Double_t dz2 = 0;

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
  tree->SetBranchAddress("uls_phi_azimuthal_angle", &pairPhi);
  tree->SetBranchAddress("uls_kaon1_azimuthal_angle", &phi1);
  tree->SetBranchAddress("uls_kaon2_azimuthal_angle", &phi2);
  tree->SetBranchAddress("uls_kaon1_nhitsDEdx", &nd1);
  tree->SetBranchAddress("uls_kaon2_nhitsDEdx", &nd2);
  tree->SetBranchAddress("uls_kaon1_nsigma_pion", &nsPi1);
  tree->SetBranchAddress("uls_kaon2_nsigma_pion", &nsPi2);
  tree->SetBranchAddress("uls_kaon1_nsigma_electron", &nsE1);
  tree->SetBranchAddress("uls_kaon2_nsigma_electron", &nsE2);
  tree->SetBranchAddress("uls_kaon1_nsigma_proton", &nsP1);
  tree->SetBranchAddress("uls_kaon2_nsigma_proton", &nsP2);
  tree->SetBranchAddress("uls_phi_x2pipi", &chi2Pion);
  tree->SetBranchAddress("uls_kaon1_dEdx", &de1);
  tree->SetBranchAddress("uls_kaon2_dEdx", &de2);
  tree->SetBranchAddress("uls_kaon1_dca_global_x", &dx1);
  tree->SetBranchAddress("uls_kaon1_dca_global_y", &dy1);
  tree->SetBranchAddress("uls_kaon1_dca_global_z", &dz1);
  tree->SetBranchAddress("uls_kaon2_dca_global_x", &dx2);
  tree->SetBranchAddress("uls_kaon2_dca_global_y", &dy2);
  tree->SetBranchAddress("uls_kaon2_dca_global_z", &dz2);

  // Histograms belong to the ROOT session, not to either input file.
  gROOT->cd();
  TH1D *h_mass =
      new TH1D("uls_nom_mass", ";M_{KK} [GeV/c^{2}];Pairs / 4 MeV/c^{2}", 55,
               0.98, 1.20);

  // Selection checks: omit the cut on the plotted variable.
  TH1D *h_nm1_vz = new TH1D("nm1_vz", ";V_{z} [cm];Pairs", 40, -200, 200);
  TH1D *h_nm1_vr = new TH1D("nm1_vr", ";V_{r} [cm];Pairs", 30, 0, 15);
  TH1D *h_nm1_pair_pt =
      new TH1D("nm1_pair_pt", ";p_{T}^{KK} [GeV/c];Pairs", 40, 0, 1.0);
  TH1D *h_nm1_pair_y = new TH1D("nm1_pair_y", ";y_{KK};Pairs", 24, -1.5, 1.5);
  TH1D *h_nm1_pid =
      new TH1D("nm1_pid", ";#chi^{2}_{KK} (corrected);Pairs", 30, 0, 150);
  // Pair and track diagnostics inside the mass window.
  TH1D *h_pair_pt = new TH1D("pair_pt", ";p_{T}^{KK} [GeV/c];Pairs", 12,
                             pairPtMin, pairPtMax);
  TH1D *h_pair_rapidity =
      new TH1D("pair_rapidity", ";y_{KK};Pairs", 10, -0.5, 0.5);
  TH1D *h_pair_phi =
      new TH1D("pair_phi", ";#varphi_{KK} [rad];Pairs", 12, -M_PI, M_PI);
  TH1D *h_track_pt =
      new TH1D("track_pt", ";p_{T}^{track} [GeV/c];Tracks", 25, 0, 1.0);
  TH1D *h_track_p =
      new TH1D("track_p", ";p^{track} [GeV/c];Tracks", 25, 0, 2.0);
  TH1D *h_track_eta = new TH1D("track_eta", ";#eta^{track};Tracks", 24, -3, 3);
  TH1D *h_track_phi =
      new TH1D("track_phi", ";#varphi^{track} [rad];Tracks", 12, -M_PI, M_PI);
  TH1D *h_track_nhits = new TH1D("track_nhits", ";nHits;Tracks", 16, 0, 80);
  TH1D *h_track_nhitsdedx =
      new TH1D("track_nhitsdedx", ";nHitsDedx;Tracks", 16, 0, 80);
  TH1D *h_track_dca =
      new TH1D("track_dca", ";Global DCA [cm];Tracks", 20, 0, 5);
  TH1D *h_track_nsigma_kaon =
      new TH1D("track_nsigma_kaon", ";n#sigma_{K} (uncorrected);Tracks", 24, -6, 6);
  TH1D *h_track_nsigma_kaon_corrected = new TH1D(
      "track_nsigma_kaon_corrected", ";n#sigma_{K}^{corr};Tracks", 24, -6, 6);
  TH1D *h_track_nsigma_pion =
      new TH1D("track_nsigma_pion", ";n#sigma_{#pi};Tracks", 25, -5, 45);
  TH2D *h_track_p_dedx =
      new TH2D("track_p_dedx", ";p^{track} [GeV/c];dE/dx [keV/cm]", 40, 0, 2,
               60, 0, 120);

  TH1D *h_track_nsigma_electron =
      new TH1D("track_nsigma_electron", ";n#sigma_{e};Tracks", 55, -15, 40);
  TH1D *h_track_nsigma_proton =
      new TH1D("track_nsigma_proton", ";n#sigma_{p};Tracks", 40, -30, 10);

  // Also save the same diagnostics before the mass-window requirement.
  TH1D *h_all_pair_pt = (TH1D *)h_pair_pt->Clone("all_pair_pt");
  TH1D *h_all_pair_rapidity =
      (TH1D *)h_pair_rapidity->Clone("all_pair_rapidity");
  TH1D *h_all_pair_phi = (TH1D *)h_pair_phi->Clone("all_pair_phi");
  TH1D *h_all_track_pt = (TH1D *)h_track_pt->Clone("all_track_pt");
  TH1D *h_all_track_p = (TH1D *)h_track_p->Clone("all_track_p");
  TH1D *h_all_track_eta = (TH1D *)h_track_eta->Clone("all_track_eta");
  TH1D *h_all_track_phi = (TH1D *)h_track_phi->Clone("all_track_phi");
  TH1D *h_all_track_nhits = (TH1D *)h_track_nhits->Clone("all_track_nhits");
  TH1D *h_all_track_nhitsdedx =
      (TH1D *)h_track_nhitsdedx->Clone("all_track_nhitsdedx");
  TH1D *h_all_track_dca = (TH1D *)h_track_dca->Clone("all_track_dca");
  TH1D *h_all_track_nsigma_kaon =
      (TH1D *)h_track_nsigma_kaon->Clone("all_track_nsigma_kaon");
  TH1D *h_all_track_nsigma_kaon_corrected =
      (TH1D *)h_track_nsigma_kaon_corrected->Clone(
          "all_track_nsigma_kaon_corrected");
  TH1D *h_all_track_nsigma_pion =
      (TH1D *)h_track_nsigma_pion->Clone("all_track_nsigma_pion");
  TH2D *h_all_track_p_dedx = (TH2D *)h_track_p_dedx->Clone("all_track_p_dedx");

  TH1D *h_all_track_nsigma_electron =
      (TH1D *)h_track_nsigma_electron->Clone("all_track_nsigma_electron");
  TH1D *h_all_track_nsigma_proton =
      (TH1D *)h_track_nsigma_proton->Clone("all_track_nsigma_proton");

  Long64_t counts[9] = {0}, windows[9] = {0};
  int nSelected = 0, nWindow = 0;

  // Loop over unlike-sign pairs.
  for (Long64_t entry = 0; entry < tree->GetEntries(); ++entry) {
    if (tree->GetEntry(entry) <= 0) {
      printf("Failed to read entry %lld\n", entry);
      calibration->Close();
      file->Close();
      return;
    }
    double p1 = pt1 * std::cosh(eta1);
    double p2 = pt2 * std::cosh(eta2);

    // Interpolate the TOF calibration; use the endpoint outside its range.
    double correctionP1 = p1;
    double correctionP2 = p2;
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
    double vr = std::hypot(vx, vy);
    bool inMassWindow = mass > massMin && mass < massMax;

    // Named cuts for the plots with one requirement omitted.
    bool passBBC = bbcEast == 0 && bbcWest == 0;
    bool passMultiplicity = nPrimary == 2;
    bool passVz = std::abs(vz) <= vzMax;
    bool passVr = vr <= vrMax;
    bool passHits = nhits1 >= nHitsMin && nhits2 >= nHitsMin;
    bool passPairPt = pairPt > pairPtMin && pairPt < pairPtMax;
    bool passRapidity = std::abs(pairY) < pairYMax;
    bool passKaonPID = chi2 <= chi2Max;
    bool passPionRejection = chi2Pion >= chi2PionMin;
    bool passPID = passKaonPID && passPionRejection;

    // Each selection-check histogram applies all the other cuts, with no mass
    // cut.
    if (passBBC && passMultiplicity && passVr && passHits && passPairPt &&
        passRapidity && passPID) {
      h_nm1_vz->Fill(vz);
    }
    if (passBBC && passMultiplicity && passVz && passHits && passPairPt &&
        passRapidity && passPID) {
      h_nm1_vr->Fill(vr);
    }
    if (passBBC && passMultiplicity && passVz && passVr && passHits &&
        passRapidity && passPID) {
      h_nm1_pair_pt->Fill(pairPt);
    }
    if (passBBC && passMultiplicity && passVz && passVr && passHits &&
        passPairPt && passPID) {
      h_nm1_pair_y->Fill(pairY);
    }
    if (passBBC && passMultiplicity && passVz && passVr && passHits &&
        passPairPt && passRapidity && passPionRejection) {
      h_nm1_pid->Fill(chi2);
    }
    // Apply the nominal cuts in order and count pairs after each step.
    counts[0]++;
    if (inMassWindow)
      windows[0]++;
    if (!passBBC)
      continue;
    counts[1]++;
    if (inMassWindow)
      windows[1]++;

    if (!passMultiplicity)
      continue;
    counts[2]++;
    if (inMassWindow)
      windows[2]++;

    if (!(passVz && passVr))
      continue;
    counts[3]++;
    if (inMassWindow)
      windows[3]++;

    if (!passHits)
      continue;
    counts[4]++;
    if (inMassWindow)
      windows[4]++;

    if (!passPairPt)
      continue;
    counts[5]++;
    if (inMassWindow)
      windows[5]++;

    if (!passRapidity)
      continue;
    counts[6]++;
    if (inMassWindow)
      windows[6]++;

    if (!passKaonPID)
      continue;
    counts[7]++;
    if (inMassWindow)
      windows[7]++;

    if (!passPionRejection)
      continue;
    counts[8]++;
    if (inMassWindow)
      windows[8]++;

    nSelected++;
    h_mass->Fill(mass);

    // Wrap azimuths to [-pi, pi] for the plotted range.
    double pairAngle = std::atan2(std::sin(pairPhi), std::cos(pairPhi));
    double trackAngle1 = std::atan2(std::sin(phi1), std::cos(phi1));
    double trackAngle2 = std::atan2(std::sin(phi2), std::cos(phi2));
    double dca1 = std::sqrt(dx1 * dx1 + dy1 * dy1 + dz1 * dz1);
    double dca2 = std::sqrt(dx2 * dx2 + dy2 * dy2 + dz2 * dz2);

    // Full selected sample: one entry per pair, two entries per track
    // histogram.
    h_all_pair_pt->Fill(pairPt);
    h_all_pair_rapidity->Fill(pairY);
    h_all_pair_phi->Fill(pairAngle);
    h_all_track_pt->Fill(pt1);
    h_all_track_pt->Fill(pt2);
    h_all_track_p->Fill(p1);
    h_all_track_p->Fill(p2);
    h_all_track_eta->Fill(eta1);
    h_all_track_eta->Fill(eta2);
    h_all_track_phi->Fill(trackAngle1);
    h_all_track_phi->Fill(trackAngle2);
    h_all_track_nhits->Fill(nhits1);
    h_all_track_nhits->Fill(nhits2);
    h_all_track_nhitsdedx->Fill(nd1);
    h_all_track_nhitsdedx->Fill(nd2);
    h_all_track_dca->Fill(dca1);
    h_all_track_dca->Fill(dca2);
    h_all_track_nsigma_kaon->Fill(nsk1);
    h_all_track_nsigma_kaon->Fill(nsk2);
    h_all_track_nsigma_kaon_corrected->Fill(nskCorr1);
    h_all_track_nsigma_kaon_corrected->Fill(nskCorr2);
    h_all_track_nsigma_pion->Fill(nsPi1);
    h_all_track_nsigma_pion->Fill(nsPi2);
    h_all_track_nsigma_electron->Fill(nsE1);
    h_all_track_nsigma_electron->Fill(nsE2);
    h_all_track_nsigma_proton->Fill(nsP1);
    h_all_track_nsigma_proton->Fill(nsP2);
    h_all_track_p_dedx->Fill(p1, de1 * 1e6); // GeV/cm to keV/cm
    h_all_track_p_dedx->Fill(p2, de2 * 1e6);

    // Nominal raw yield and the diagnostics shown in the note.
    if (!inMassWindow)
      continue;
    nWindow++;
    h_pair_pt->Fill(pairPt);
    h_pair_rapidity->Fill(pairY);
    h_pair_phi->Fill(pairAngle);
    h_track_pt->Fill(pt1);
    h_track_pt->Fill(pt2);
    h_track_p->Fill(p1);
    h_track_p->Fill(p2);
    h_track_eta->Fill(eta1);
    h_track_eta->Fill(eta2);
    h_track_phi->Fill(trackAngle1);
    h_track_phi->Fill(trackAngle2);
    h_track_nhits->Fill(nhits1);
    h_track_nhits->Fill(nhits2);
    h_track_nhitsdedx->Fill(nd1);
    h_track_nhitsdedx->Fill(nd2);
    h_track_dca->Fill(dca1);
    h_track_dca->Fill(dca2);
    h_track_nsigma_kaon->Fill(nsk1);
    h_track_nsigma_kaon->Fill(nsk2);
    h_track_nsigma_kaon_corrected->Fill(nskCorr1);
    h_track_nsigma_kaon_corrected->Fill(nskCorr2);
    h_track_nsigma_pion->Fill(nsPi1);
    h_track_nsigma_pion->Fill(nsPi2);
    h_track_nsigma_electron->Fill(nsE1);
    h_track_nsigma_electron->Fill(nsE2);
    h_track_nsigma_proton->Fill(nsP1);
    h_track_nsigma_proton->Fill(nsP2);
    h_track_p_dedx->Fill(p1, de1 * 1e6); // GeV/cm to keV/cm
    h_track_p_dedx->Fill(p2, de2 * 1e6);
  }
  calibration->Close();
  file->Close();
  printf("Selected ULS pairs: %d; nominal raw yield: %d\n", nSelected, nWindow);

  // Save the cutflow table used by the analysis note.
  gSystem->mkdir(outputDir, true);
  TString stages[9] = {"Stored ULS pairs",
                       "BBC veto",
                       "Two primary tracks",
                       "Vertex",
                       "Both tracks nHits >= 20",
                       "0.3 < Pair pT < 0.9",
                       "Abs(y) < 0.5",
                       Form("Corrected kaon chi2 <= %.0f", chi2Max),
                       Form("Pion chi2 >= %.0f", chi2PionMin)};
  // Math-formatted labels for the Typst note.
  TString tableStages[9] = {
      "Preselected ULS pairs",
      "BBC veto",
      R"($N_"prim" = 2$)",
      Form(R"($abs(V_z) <= %.0f thin "cm"$; $V_r <= %.0f thin "cm"$)", vzMax,
           vrMax),
      Form(R"(Both tracks: $N_"hits" >= %d$)", nHitsMin),
      Form(R"($%.1f < p_T^(K K) < %.1f thin "GeV"\/c$)", pairPtMin, pairPtMax),
      Form(R"($abs(y_(K K)) < %.1f$)", pairYMax),
      Form(R"($chi^2_(K K) <= %.0f$)", chi2Max),
      Form(R"($chi^2_(pi pi) >= %.0f$)", chi2PionMin)};
  TH1D *h_cutflow = new TH1D("cutflow", "Cumulative ULS pair counts", 9, 0, 9);
  TH1D *h_cutflow_window =
      new TH1D("cutflow_window", "Cumulative ULS window counts", 9, 0, 9);
  std::ofstream table((outputDir + "cutflow.typ").Data());
  table << "#table(columns: (2fr, 1fr, 1fr), align: (left, right, right), "
           "stroke: none, table.hline(), table.header([Cumulative selection], "
           "[All masses], [Mass window]), table.hline(),\n";
  for (int i = 0; i < 9; ++i) {
    h_cutflow->SetBinContent(i + 1, counts[i]);
    h_cutflow_window->SetBinContent(i + 1, windows[i]);
    h_cutflow->GetXaxis()->SetBinLabel(i + 1, stages[i]);
    h_cutflow_window->GetXaxis()->SetBinLabel(i + 1, stages[i]);
    table << "[" << tableStages[i] << "], [" << counts[i] << "], ["
          << windows[i] << "],\n";
    printf("%-28s %9lld total %9lld window\n", stages[i].Data(), counts[i],
           windows[i]);
  }
  table << "table.hline())\n";
  table.close();

  // Use the existing plot style for every figure.
  gStyle->SetPaperSize(26, 19);
  PlotStyle::Options style;
  style.title = "Run-19 Au+Au";
  style.subtitle = "#sqrt{s_{NN}} = 200 GeV";
  style.note = "ZDCMB; incoherent region";
  style.legendX1 = 0.56;

  // Selection plots and their cut lines, with linear axes.
  {
    PlotStyle::Options opt = style;
    opt.xlabel = h_nm1_vz->GetXaxis()->GetTitle();
    opt.ylabel = "Pairs";
    TCanvas *c = PlotStyle::Draw({{h_nm1_vz, "ULS", "E1"}}, opt);
    DrawLines({"All other nominal cuts; no mass cut"}, opt.annotationX,
              BelowNote(opt), 0.032);
    c->Update();
    TLine line0(-vzMax, 0, -vzMax, 0.55 * c->GetUymax());
    line0.SetLineColor(kRed + 1);
    line0.SetLineStyle(2);
    line0.DrawClone();
    TLine line1(vzMax, 0, vzMax, 0.55 * c->GetUymax());
    line1.SetLineColor(kRed + 1);
    line1.SetLineStyle(2);
    line1.DrawClone();
    c->Print(outputDir + "nm1_vz.pdf");
    delete c;
  }
  {
    PlotStyle::Options opt = style;
    opt.xlabel = h_nm1_vr->GetXaxis()->GetTitle();
    opt.ylabel = "Pairs";
    TCanvas *c = PlotStyle::Draw({{h_nm1_vr, "ULS", "E1"}}, opt);
    DrawLines({"All other nominal cuts; no mass cut"}, opt.annotationX,
              BelowNote(opt), 0.032);
    c->Update();
    TLine line0(vrMax, 0, vrMax, 0.55 * c->GetUymax());
    line0.SetLineColor(kRed + 1);
    line0.SetLineStyle(2);
    line0.DrawClone();
    c->Print(outputDir + "nm1_vr.pdf");
    delete c;
  }
  {
    PlotStyle::Options opt = style;
    opt.xlabel = h_nm1_pair_pt->GetXaxis()->GetTitle();
    opt.ylabel = "Pairs";
    TCanvas *c = PlotStyle::Draw({{h_nm1_pair_pt, "ULS", "E1"}}, opt);
    DrawLines({"All other nominal cuts; no mass cut"}, opt.annotationX,
              BelowNote(opt), 0.032);
    c->Update();
    TLine lowerPt(pairPtMin, 0, pairPtMin, 0.55 * c->GetUymax());
    lowerPt.SetLineColor(kRed + 1);
    lowerPt.SetLineStyle(2);
    lowerPt.DrawClone();
    TLine line0(pairPtMax, 0, pairPtMax, 0.55 * c->GetUymax());
    line0.SetLineColor(kRed + 1);
    line0.SetLineStyle(2);
    line0.DrawClone();
    c->Print(outputDir + "nm1_pair_pt.pdf");
    delete c;
  }
  {
    PlotStyle::Options opt = style;
    opt.xlabel = h_nm1_pair_y->GetXaxis()->GetTitle();
    opt.ylabel = "Pairs";
    TCanvas *c = PlotStyle::Draw({{h_nm1_pair_y, "ULS", "E1"}}, opt);
    DrawLines({"All other nominal cuts; no mass cut"}, opt.annotationX,
              BelowNote(opt), 0.032);
    c->Update();
    TLine line0(-pairYMax, 0, -pairYMax, 0.55 * c->GetUymax());
    line0.SetLineColor(kRed + 1);
    line0.SetLineStyle(2);
    line0.DrawClone();
    TLine line1(pairYMax, 0, pairYMax, 0.55 * c->GetUymax());
    line1.SetLineColor(kRed + 1);
    line1.SetLineStyle(2);
    line1.DrawClone();
    c->Print(outputDir + "nm1_pair_y.pdf");
    delete c;
  }
  {
    PlotStyle::Options opt = style;
    opt.xlabel = h_nm1_pid->GetXaxis()->GetTitle();
    opt.ylabel = "Pairs";
    TCanvas *c = PlotStyle::Draw({{h_nm1_pid, "ULS", "E1"}}, opt);
    DrawLines({"All other nominal cuts; no mass cut"}, opt.annotationX,
              BelowNote(opt), 0.032);
    c->Update();
    TLine line0(chi2Max, 0, chi2Max, 0.55 * c->GetUymax());
    line0.SetLineColor(kRed + 1);
    line0.SetLineStyle(2);
    line0.DrawClone();
    c->Print(outputDir + "nm1_pid.pdf");
    delete c;
  }
  // The pair and track plots share the same drawing settings and error bars.
  TH1D *diagnostics[] = {h_pair_pt,
                         h_pair_rapidity,
                         h_pair_phi,
                         h_track_pt,
                         h_track_p,
                         h_track_eta,
                         h_track_phi,
                         h_track_nhits,
                         h_track_nhitsdedx,
                         h_track_dca,
                         h_track_nsigma_kaon,
                         h_track_nsigma_kaon_corrected,
                         h_track_nsigma_pion,
                         h_track_nsigma_electron,
                         h_track_nsigma_proton};
  for (int i = 0; i < 15; ++i) {
    TH1D *hist = diagnostics[i];
    PlotStyle::Options opt = style;
    opt.xlabel = hist->GetXaxis()->GetTitle();
    opt.ylabel = hist->GetYaxis()->GetTitle();
    if (hist == h_track_nsigma_electron || hist == h_track_nsigma_proton) {
      int first = hist->FindFirstBinAbove(0);
      int last = hist->FindLastBinAbove(0);
      if (first > 0 && last >= first)
        opt.XRange(hist->GetXaxis()->GetBinLowEdge(first) - 1,
                   hist->GetXaxis()->GetBinUpEdge(last) + 1);
    }
    if (hist == h_track_p)
      opt.XRange(0, 0.9);
    TCanvas *c = PlotStyle::Draw({{hist, "ULS", "E1"}}, opt);
    DrawLines({"Selected ULS, 1.00 < M_{KK} < 1.04 GeV/c^{2}"}, opt.annotationX,
              BelowNote(opt), 0.032);
    c->Print(outputDir + hist->GetName() + ".pdf");
    if (hist == h_track_p)
      c->Print(outputDir + "track_p.png");
    delete c;
  }
  {
    PlotStyle::Options opt = style;
    opt.annotationX = 0.48;
    opt.xlabel = "p^{track} [GeV/c]";
    opt.ylabel = "dE/dx [keV/cm]";
    TCanvas *c = Draw2D(h_track_p_dedx, opt);
    DrawLines({Form("%d tracks from %d window pairs", 2 * nWindow, nWindow)},
              opt.annotationX, BelowNote(opt) - 0.050, 0.032);
    c->Print(outputDir + "track_p_dedx.pdf");
    delete c;
  }
  {
    PlotStyle::Options opt = style;
    opt.xlabel = "M_{KK} [GeV/c^{2}]";
    opt.ylabel = "Pairs / 4 MeV/c^{2}";
    opt.XRange(0.98, 1.20);
    TCanvas *c = PlotStyle::Draw({{h_mass, "ULS", "E1"}}, opt);
    DrawLines({Form("Nominal raw yield: %d", nWindow)}, opt.annotationX,
              BelowNote(opt), 0.032);
    c->Update();
    TLine lower(massMin, 0, massMin, 0.45 * c->GetUymax());
    lower.SetLineStyle(2);
    lower.SetLineColor(kRed + 1);
    lower.DrawClone();
    TLine upper(massMax, 0, massMax, 0.45 * c->GetUymax());
    upper.SetLineStyle(2);
    upper.SetLineColor(kRed + 1);
    upper.DrawClone();
    c->Print(outputDir + "mass.pdf");
    c->Print(outputDir + "mass.png");
    delete c;
  }

  // Write every histogram to the ROOT file.
  TFile *output = TFile::Open(outputDir + "incoherent_yield.root", "RECREATE");
  if (!output || output->IsZombie()) {
    printf("Cannot create incoherent_yield.root\n");
    return;
  }
  h_cutflow->Write();
  h_cutflow_window->Write();
  h_mass->Write();
  h_mass->Write("pair_mass");
  h_nm1_vz->Write();
  h_nm1_vr->Write();
  h_nm1_pair_pt->Write();
  h_nm1_pair_y->Write();
  h_nm1_pid->Write();
  h_pair_pt->Write();
  h_pair_rapidity->Write();
  h_pair_phi->Write();
  h_track_pt->Write();
  h_track_p->Write();
  h_track_eta->Write();
  h_track_phi->Write();
  h_track_nhits->Write();
  h_track_nhitsdedx->Write();
  h_track_dca->Write();
  h_track_nsigma_kaon->Write();
  h_track_nsigma_kaon_corrected->Write();
  h_track_nsigma_pion->Write();
  h_track_nsigma_electron->Write();
  h_track_nsigma_proton->Write();
  h_track_p_dedx->Write();
  TParameter<int>("raw_uls_selected", nSelected).Write();
  TParameter<int>("raw_uls_window", nWindow).Write();
  TParameter<double>("corrected_chi2_kk_max", chi2Max).Write();
  TParameter<double>("chi2_pipi_min", chi2PionMin).Write();
  TParameter<int>("electron_rejection_applied", 0).Write();
  TParameter<int>("final_track_pt_cut_applied", 0).Write();
  TParameter<double>("pair_pt_min_GeV", pairPtMin).Write();
  TParameter<double>("pair_pt_max_GeV", pairPtMax).Write();
  TParameter<double>("pair_abs_y_max", pairYMax).Write();
  TParameter<double>("mass_window_low_GeV", massMin).Write();
  TParameter<double>("mass_window_high_GeV", massMax).Write();
  TNamed(
      "selection",
      "BBC veto; two primary tracks; |Vz|<=100 cm; Vr<=5 cm; both nHits>=20; "
      "0.3<pair pT<0.9 GeV/c; |pair y|<0.5; corrected chi2KK<=10; "
      "chi2pipi>=15; "
      "no electron or proton rejection; no "
      "additional track pT cut; stored preselection unchanged.")
      .Write();
  TNamed("histogram_samples",
         "Top-level pair_pt, pair_rapidity, pair_phi and track_*: strict mass "
         "window; pair_mass and uls_nom_mass: all nominal selected pairs; "
         "selected_all_masses/: all nominal selected pairs and their tracks. "
         "Raw unweighted counts.")
      .Write();

  output->mkdir("selected_all_masses");
  output->cd("selected_all_masses");
  h_all_pair_pt->Write("pair_pt");
  h_all_pair_rapidity->Write("pair_rapidity");
  h_all_pair_phi->Write("pair_phi");
  h_all_track_pt->Write("track_pt");
  h_all_track_p->Write("track_p");
  h_all_track_eta->Write("track_eta");
  h_all_track_phi->Write("track_phi");
  h_all_track_nhits->Write("track_nhits");
  h_all_track_nhitsdedx->Write("track_nhitsdedx");
  h_all_track_dca->Write("track_dca");
  h_all_track_nsigma_kaon->Write("track_nsigma_kaon");
  h_all_track_nsigma_kaon_corrected->Write("track_nsigma_kaon_corrected");
  h_all_track_nsigma_pion->Write("track_nsigma_pion");
  h_all_track_nsigma_electron->Write("track_nsigma_electron");
  h_all_track_nsigma_proton->Write("track_nsigma_proton");
  h_all_track_p_dedx->Write("track_p_dedx");
  h_mass->Write("pair_mass");
  output->Close();
}
