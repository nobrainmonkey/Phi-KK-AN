#include "../events/VzFit.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"
#include <cmath>
#include <iostream>
#include <memory>
#include <vector>

// QA of the TOF kaon sample (KaonTree from StPhiKKMaker) used for the nSigmaKaon calibration in
// kaon_tof.C. Production selection: triggers 700001/27/26, bad runs removed, |Vz| < 50 cm,
// Vr < 2 cm, centrality 70-80% (cent9 = 0; the producer code accepts
// cent9 = 0 or 1, but this file contains only cent9 = 0); primary tracks of vertex 0 with nHits >= 20,
// nHits/nHitsPoss > 0.52, TOF matched, 0 <= beta < 1, 0.22 < m^2 < 0.28 GeV^2/c^4 (the pasted
// producer code says 0.23, but this file is filled from 0.22).
// Histograms stored in the file (hM2, hBeta, ...) are filled before the m^2 window.
// Output: PDF files in PID/.

const TString kOut = "PID/";

PlotStyle::Options BaseOptions(const TString& xlabel, const TString& ylabel) {
  PlotStyle::Options opt;
  opt.title = "Run-19 Au+Au";
  opt.subtitle = "#sqrt{s_{NN}} = 200 GeV";
  opt.note = "70-80%, TOF kaon sample";
  opt.xlabel = xlabel;
  opt.ylabel = ylabel;
  opt.showLegend = false;
  return opt;
}

void Plot1D(TH1* h, const TString& xlabel, const TString& ylabel, const TString& name, bool logY = false,
            double xlo = 0, double xhi = 0) {
  PlotStyle::Options opt = BaseOptions(xlabel, ylabel);
  opt.logY = logY;
  if (xhi > xlo) opt.XRange(xlo, xhi);
  TCanvas* c = PlotStyle::Draw(h, opt);
  if (c) c->Print(kOut + name + ".pdf");
}

void Plot2D(TH2* h, const TString& xlabel, const TString& ylabel, const TString& name, double annotationX = 0.18) {
  PlotStyle::Options opt = BaseOptions(xlabel, ylabel);
  opt.logY = true; // log z
  opt.annotationX = annotationX;
  Draw2D(h, opt)->Print(kOut + name + ".pdf");
}

void kaon_tof_qa(TString input_file = "/Users/xihehan/phi-kk-results/data/kaon_mudst.root") {
  std::unique_ptr<TFile> file(TFile::Open(input_file, "READ"));
  if (!file || file->IsZombie()) return;
  TTree* tree = file->Get<TTree>("KaonTree");
  if (!tree) {
    std::cerr << "Missing KaonTree in " << input_file << "\n";
    return;
  }
  gSystem->mkdir(kOut, true);

  // ---- Histograms stored by the producer (before the m^2 window) --------------
  for (const char* n : {"hEventCounter", "hVertexZ", "hVertexR", "hCent9", "hRefMult", "hM2", "hBeta"}) {
    TH1* h = file->Get<TH1>(n);
    if (!h) continue;
    h->SetDirectory(nullptr);
    const TString name(n);
    if (name == "hEventCounter") {
      // Labelled bins: draw directly so the bin labels are kept.
      TCanvas* c = new TCanvas("c_event_counter", "", 1100, 800);
      c->SetLeftMargin(0.14);
      c->SetBottomMargin(0.22);
      c->SetLogy();
      h->SetStats(false);
      h->SetTitle("");
      h->SetLineColor(PlotStyle::Color(0));
      h->SetLineWidth(2);
      h->GetYaxis()->SetTitle("Events");
      h->SetMarkerSize(1.2);
      gStyle->SetPaintTextFormat(".3g");
      h->Draw("HIST TEXT0");
      c->Print(kOut + "event_counter.pdf");
    }
    if (name == "hVertexZ")      Plot1D(h, "V_{z} [cm] (before cut)", "Events", "event_vz_all");
    if (name == "hVertexR")      Plot1D(h, "V_{r} [cm] (before cut)", "Events", "event_vr_all", true);
    if (name == "hCent9")        Plot1D(h, "Centrality bin (cent9)", "Events", "event_cent9");
    if (name == "hRefMult")      Plot1D(h, "RefMult", "Events", "event_refmult", true);
    if (name == "hM2")           Plot1D(h, "m^{2} [GeV^{2}/c^{4}] (all TOF-matched tracks)", "Tracks", "track_m2_all", true);
    if (name == "hBeta")         Plot1D(h, "#beta (all TOF-matched tracks)", "Tracks", "track_beta_all", true, 0, 1.05);
  }

  // ---- Tree branches ------------------------------------------------------------
  std::vector<double>* track_pt = nullptr;
  std::vector<double>* track_eta = nullptr;
  std::vector<double>* track_phi = nullptr;
  std::vector<double>* track_dca = nullptr;
  std::vector<double>* track_nSigmaKaon = nullptr;
  std::vector<double>* track_nSigmaPion = nullptr;
  std::vector<double>* track_tof_beta = nullptr;
  std::vector<double>* track_m2 = nullptr;
  std::vector<UShort_t>* track_nHits = nullptr;
  std::vector<UShort_t>* track_nHitsFit = nullptr;
  std::vector<UShort_t>* track_nHitsPoss = nullptr;
  std::vector<Short_t>* track_charge = nullptr;
  Float_t event_vz, event_vr;
  Int_t event_nPrimaryTracks, event_refMult;

  tree->SetBranchStatus("*", false);
  auto bind = [&](const char* name, void* address) {
    tree->SetBranchStatus(name, true);
    tree->SetBranchAddress(name, address);
  };
  bind("track_pt", &track_pt);
  bind("track_eta", &track_eta);
  bind("track_phi", &track_phi);
  bind("track_dca", &track_dca);
  bind("track_nSigmaKaon", &track_nSigmaKaon);
  bind("track_nSigmaPion", &track_nSigmaPion);
  bind("track_tof_beta", &track_tof_beta);
  bind("track_m2", &track_m2);
  bind("track_nHits", &track_nHits);
  bind("track_nHitsFit", &track_nHitsFit);
  bind("track_nHitsPoss", &track_nHitsPoss);
  bind("track_charge", &track_charge);
  bind("event_vz", &event_vz);
  bind("event_vr", &event_vr);
  bind("event_nPrimaryTracks", &event_nPrimaryTracks);
  bind("event_refMult", &event_refMult);

  // Events (saved events only: at least one kaon candidate)
  TH1D* hVz = new TH1D("hVz", "", 100, -50, 50);
  TH1D* hVr = new TH1D("hVr", "", 100, 0, 2);
  TH1D* hNPrim = new TH1D("hNPrim", "", 200, 0, 400);
  TH1D* hNKaon = new TH1D("hNKaon", "", 10, 0.5, 10.5);

  // Kaon candidates
  TH1D* hPt = new TH1D("hPt", "", 125, 0, 2.5);
  TH1D* hP = new TH1D("hP", "", 150, 0, 3);
  TH1D* hEta = new TH1D("hEta", "", 100, -1.5, 1.5);
  TH1D* hPhi = new TH1D("hPhi", "", 90, -M_PI, M_PI);
  TH1D* hCharge = new TH1D("hCharge", "", 3, -1.5, 1.5);
  TH1D* hNHits = new TH1D("hNHits", "", 80, 0, 80);
  TH1D* hNHitsFit = new TH1D("hNHitsFit", "", 80, 0, 80);
  TH1D* hNHitsPoss = new TH1D("hNHitsPoss", "", 80, 0, 80);
  TH1D* hDca = new TH1D("hDca", "", 150, 0, 3);
  TH1D* hM2 = new TH1D("hM2", "", 100, 0.2, 0.3);
  TH1D* hNSigK = new TH1D("hNSigK", "", 200, -10, 10);
  TH1D* hNSigPi = new TH1D("hNSigPi", "", 200, -10, 20);
  TH2D* hPM2 = new TH2D("hPM2", "", 150, 0, 3, 100, 0.2, 0.3);
  TH2D* hPInvBeta = new TH2D("hPInvBeta", "", 150, 0, 3, 150, 0.9, 2.5);
  TH2D* hPNSigK = new TH2D("hPNSigK", "", 150, 0, 3, 224, -8, 20);
  TH2D* hPNSigPi = new TH2D("hPNSigPi", "", 150, 0, 3, 200, -5, 20);
  TH2D* hEtaPhi = new TH2D("hEtaPhi", "", 100, -1.5, 1.5, 90, -M_PI, M_PI);

  const Long64_t nEntries = tree->GetEntries();
  for (Long64_t i = 0; i < nEntries; ++i) {
    tree->GetEntry(i);
    hVz->Fill(event_vz);
    hVr->Fill(event_vr);
    hNPrim->Fill(event_nPrimaryTracks);
    hNKaon->Fill(track_pt->size());

    for (size_t j = 0; j < track_pt->size(); ++j) {
      const double pt = (*track_pt)[j], eta = (*track_eta)[j], p = pt * std::cosh(eta);
      const double nsk = (*track_nSigmaKaon)[j], nspi = (*track_nSigmaPion)[j];
      const double beta = (*track_tof_beta)[j], m2 = (*track_m2)[j];

      hPt->Fill(pt);
      hP->Fill(p);
      hEta->Fill(eta);
      hPhi->Fill((*track_phi)[j]);
      hCharge->Fill((*track_charge)[j]);
      hNHits->Fill((*track_nHits)[j]);
      hNHitsFit->Fill((*track_nHitsFit)[j]);
      hNHitsPoss->Fill((*track_nHitsPoss)[j]);
      hDca->Fill((*track_dca)[j]);
      hM2->Fill(m2);
      hNSigK->Fill(nsk);
      hNSigPi->Fill(nspi);
      hPM2->Fill(p, m2);
      if (beta > 0) hPInvBeta->Fill(p, 1 / beta);
      hPNSigK->Fill(p, nsk);
      hPNSigPi->Fill(p, nspi);
      hEtaPhi->Fill(eta, (*track_phi)[j]);
    }
  }

  // ---- Plots ------------------------------------------------------------------------
  Plot1D(hVz, "V_{z} [cm]", "Events", "event_vz");
  Plot1D(hVr, "V_{r} [cm]", "Events", "event_vr");
  Plot1D(hNPrim, "Primary tracks (vertex 0)", "Events", "event_nprimary", true);
  Plot1D(hNKaon, "Kaon candidates per event", "Events", "event_nkaon", true);

  Plot1D(hPt, "p_{T} [GeV/c]", "Tracks", "track_pt");
  Plot1D(hP, "p [GeV/c]", "Tracks", "track_p");
  Plot1D(hEta, "#eta", "Tracks", "track_eta");
  Plot1D(hPhi, "#varphi [rad]", "Tracks", "track_phi");
  Plot1D(hCharge, "Charge", "Tracks", "track_charge");
  Plot1D(hNHits, "nHits", "Tracks", "track_nhits");
  Plot1D(hNHitsFit, "nHitsFit", "Tracks", "track_nhitsfit");
  Plot1D(hNHitsPoss, "nHitsPoss", "Tracks", "track_nhitsposs");
  Plot1D(hDca, "Global DCA [cm]", "Tracks", "track_dca", true);
  Plot1D(hM2, "m^{2} [GeV^{2}/c^{4}]", "Tracks", "track_m2");
  Plot1D(hNSigK, "n#sigma_{K}", "Tracks", "track_nsigma_kaon");
  Plot1D(hNSigPi, "n#sigma_{#pi}", "Tracks", "track_nsigma_pion");

  Plot2D(hPM2, "p [GeV/c]", "m^{2} [GeV^{2}/c^{4}]", "track_p_m2", 0.45);
  Plot2D(hPInvBeta, "p [GeV/c]", "1/#beta", "track_p_invbeta", 0.45);
  Plot2D(hPNSigK, "p [GeV/c]", "n#sigma_{K}", "track_p_nsigma_kaon", 0.45);
  Plot2D(hPNSigPi, "p [GeV/c]", "n#sigma_{#pi}", "track_p_nsigma_pion", 0.45);
  Plot2D(hEtaPhi, "#eta", "#varphi [rad]", "track_eta_phi");
}
