#include "../events/PlotStyle.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1D.h"
#include "TLine.h"
#include "TMath.h"
#include "TParameter.h"
#include "TSystem.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include <cmath>
#include <iostream>
#include <memory>
#include <vector>

// Mean nSigmaKaon vs momentum for kaons from phi -> K+K- in the UPC data, for comparison with the
// TOF kaons of kaon_tof.C (run that first for the overlay).
//
// Selection (ULS pairs): BBC veto, exactly 2 primary tracks in the event, |Vz| <= 100 cm,
// Vr < 5 cm (event vertex), pair rapidity |y| < 0.5; both tracks |nSigmaKaon| <= 5.
// Tag and probe: a track is used (probe) when the other track (tag) passes the pion rejection
// |nSigmaPion| >= 3, so the pion rejection is never applied to the track whose nSigmaKaon is measured.
// Probes from pairs in the phi window 1.00 < M_KK < 1.04 GeV/c^2 are used. The non-phi background
// is described by the sideband 1.06 < M_KK < 1.10 GeV/c^2, scaled to the background in the window
// from a fit of the probe mass spectrum (Voigtian with the PDG width + pol1).
// In each 0.1 GeV/c momentum bin with at least 50 window probes, the window nSigmaKaon distribution
// is fitted (binned Poisson likelihood) with a Gaussian for the kaons plus the scaled sideband
// distribution as a fixed background template.
// Output: PID/phi/.

const double kGammaPhi = 0.004249;           // PDG phi width [GeV]
const double kWinLo = 1.00, kWinHi = 1.04;   // phi window [GeV/c^2]
const double kSideLo = 1.06, kSideHi = 1.10; // sideband [GeV/c^2]
const int kMinTracks = 50;                   // window tracks needed for a momentum bin
const TString kOut = "PID/phi/";

double MassModel(double* x, double* p) {
  return p[0] * TMath::Voigt(x[0] - p[1], p[2], kGammaPhi) + p[3] + p[4] * (x[0] - 1.02);
}

// Background part of MassModel (p[0] + p[1] (m - 1.02)), integrated for the background count.
double MassBkg(double* x, double* p) { return p[0] + p[1] * (x[0] - 1.02); }

// nSigmaKaon model in one momentum bin: N_K kaons in a normalised Gaussian plus the sideband
// histogram gSide scaled by gScale as a fixed background. gBinW is the nSigmaKaon bin width, so
// p[0] is the number of kaons; p[1] = mu, p[2] = sigma. Set gSide, gScale, gBinW before fitting.
TH1D* gSide = nullptr;
double gScale = 0, gBinW = 0;
double NSigmaModel(double* x, double* p) {
  return p[0] * gBinW * TMath::Gaus(x[0], p[1], p[2], true) + gScale * gSide->GetBinContent(gSide->FindBin(x[0]));
}

void kaon_phi(TString input_file = "../../data/phi-kk-mudst.root") {
  std::unique_ptr<TFile> file(TFile::Open(input_file, "READ"));
  if (!file || file->IsZombie()) return;

  TTreeReader reader("uls_data", file.get());
  TTreeReaderValue<UShort_t> bbc_east(reader, "uls_pair_BBC_Veto_East");
  TTreeReaderValue<UShort_t> bbc_west(reader, "uls_pair_BBC_Veto_West");
  TTreeReaderValue<UInt_t> n_primary(reader, "uls_pair_event_n_primary_tracks");
  TTreeReaderValue<Double_t> vx(reader, "uls_pair_event_vx"), vy(reader, "uls_pair_event_vy"), vz(reader, "uls_pair_event_vz");
  TTreeReaderValue<Double_t> mass(reader, "uls_phi_mass");
  TTreeReaderValue<Double_t> rapidity(reader, "uls_phi_rapidity");
  TTreeReaderValue<Double_t> pt1(reader, "uls_kaon1_pt"), pt2(reader, "uls_kaon2_pt");
  TTreeReaderValue<Double_t> eta1(reader, "uls_kaon1_eta"), eta2(reader, "uls_kaon2_eta");
  TTreeReaderValue<Double_t> nsk1(reader, "uls_kaon1_nsigma_kaon"), nsk2(reader, "uls_kaon2_nsigma_kaon");
  TTreeReaderValue<Double_t> nspi1(reader, "uls_kaon1_nsigma_pion"), nspi2(reader, "uls_kaon2_nsigma_pion");

  // Momentum bins: 0.1 - 1.0 GeV/c in steps of 0.1 (bins with too few tracks are skipped).
  const std::vector<double> pBins = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
  const int nP = pBins.size() - 1;

  TH1D* hMass = new TH1D("hMass", "", 60, 0.98, 1.10);
  std::vector<TH1D*> hWin(nP), hSide(nP);
  for (int b = 0; b < nP; ++b) {
    hWin[b] = new TH1D(Form("hWin_pbin_%d", b), Form("%.1f < p < %.1f;n#sigma_{K};Tracks", pBins[b], pBins[b + 1]), 50, -5, 5);
    hSide[b] = new TH1D(Form("hSide_pbin_%d", b), "", 50, -5, 5);
    hWin[b]->Sumw2();
    hSide[b]->Sumw2();
  }

  while (reader.Next()) {
    if (*bbc_east != 0 || *bbc_west != 0) continue;
    if (*n_primary != 2) continue;
    if (std::abs(*vz) > 100 || std::hypot(*vx, *vy) >= 5) continue;
    if (std::abs(*rapidity) >= 0.5) continue;
    if (std::abs(*nsk1) > 5 || std::abs(*nsk2) > 5) continue;
    const bool inWin = *mass > kWinLo && *mass < kWinHi;
    const bool inSide = *mass > kSideLo && *mass < kSideHi;

    const double p[2] = {*pt1 * std::cosh(*eta1), *pt2 * std::cosh(*eta2)};
    const double nsk[2] = {*nsk1, *nsk2};
    const double nspi[2] = {*nspi1, *nspi2};
    for (int k = 0; k < 2; ++k) {
      if (std::abs(nspi[1 - k]) < 3) continue; // tag (other track) pion rejection
      hMass->Fill(*mass);
      if (!inWin && !inSide) continue;
      for (int b = 0; b < nP; ++b) {
        if (p[k] >= pBins[b] && p[k] < pBins[b + 1]) {
          (inWin ? hWin : hSide)[b]->Fill(nsk[k]);
          break;
        }
      }
    }
  }

  // ---- Background in the window from the mass fit, and the sideband scale ----------
  TF1* fMass = new TF1("fMass", MassModel, 0.99, 1.10, 5);
  fMass->SetParameters(5, 1.019, 0.0015, 5, 0);
  fMass->SetParLimits(2, 0.0003, 0.005);
  hMass->Fit(fMass, "QLR0");
  TF1 fBkg("fBkg", MassBkg, 0.99, 1.10, 2);
  fBkg.SetParameters(fMass->GetParameter(3), fMass->GetParameter(4));
  const double nWin = hMass->Integral(hMass->FindBin(kWinLo + 1e-6), hMass->FindBin(kWinHi - 1e-6));
  const double nSide = hMass->Integral(hMass->FindBin(kSideLo + 1e-6), hMass->FindBin(kSideHi - 1e-6));
  const double bWin = fBkg.Integral(kWinLo, kWinHi) / hMass->GetBinWidth(1);
  const double scale = bWin / nSide;
  printf("phi window: %.0f probes, background %.0f (signal fraction %.2f); sideband %.0f probes, scale %.3f\n",
         nWin, bWin, 1 - bWin / nWin, nSide, scale);

  // ---- Fit of the window nSigmaKaon in each momentum bin ------------------------------
  // Binned Poisson likelihood fit of the window histogram with
  //   N_K * Gaussian(mu, sigma) + scale * (sideband histogram),
  // so that low-count bins are treated correctly (no fit of a subtracted histogram).
  // The sideband template is fixed; N_K, mu, sigma are free.
  TH1D* hMean = new TH1D("hMeanNSigmaKaonPhi", "", nP, pBins.data());
  TH1D* hNK = new TH1D("hNKaonPhi", "", nP, pBins.data());     // fitted number of kaons per bin
  std::vector<TH1D*> hSub(nP);
  std::vector<TF1*> fits(nP, nullptr);   // signal Gaussian only, for drawing on the subtracted histogram
  for (int b = 0; b < nP; ++b) {
    hSub[b] = static_cast<TH1D*>(hWin[b]->Clone(Form("hSub_pbin_%d", b)));
    hSub[b]->Add(hSide[b], -scale);
    if (hWin[b]->GetEntries() < kMinTracks) continue;

    gSide = hSide[b];
    gScale = scale;
    gBinW = hWin[b]->GetBinWidth(1);
    TF1 model(Form("model_%d", b), NSigmaModel, -5, 5, 3);
    model.SetParameters(std::max(1., hSub[b]->Integral()), hSub[b]->GetMean(), 1);
    model.SetParLimits(2, 0.3, 3);
    hWin[b]->Fit(&model, "QLR0");

    fits[b] = new TF1(Form("signal_%d", b), Form("[0] * %g * TMath::Gaus(x, [1], [2], true)", gBinW), -5, 5);
    fits[b]->SetParameters(model.GetParameter(0), model.GetParameter(1), model.GetParameter(2));
    fits[b]->SetParError(1, model.GetParError(1));
    hMean->SetBinContent(b + 1, model.GetParameter(1));
    hMean->SetBinError(b + 1, model.GetParError(1));
    hNK->SetBinContent(b + 1, model.GetParameter(0));
    printf("p [%.1f, %.1f): %4.0f window probes, N_K = %5.1f, <nSigmaKaon> = %6.3f +/- %.3f, sigma = %.2f\n", pBins[b],
           pBins[b + 1], hWin[b]->GetEntries(), model.GetParameter(0), model.GetParameter(1), model.GetParError(1),
           std::abs(model.GetParameter(2)));
  }

  // ---- Output ----------------------------------------------------------------
  gSystem->mkdir(kOut + "fits", true);
  TFile out(kOut + "nsigmakaon_phi_mean_vs_p.root", "RECREATE");
  hMass->Write();
  hMean->Write();
  hNK->Write();
  TParameter<double>("sidebandScale", scale).Write();
  for (int b = 0; b < nP; ++b) {
    hWin[b]->Write();
    hSide[b]->Write();
    hSub[b]->Write();
  }
  out.Close();

  // Mass spectrum with the fit, window, and sideband.
  PlotStyle::Options om;
  om.title = "Run-19 Au+Au";
  om.subtitle = "#sqrt{s_{NN}} = 200 GeV";
  om.note = "Tag |n#sigma_{#pi}| #geq 3, probe |n#sigma_{K}| #leq 5, |y| < 0.5";
  om.xlabel = "M_{KK} [GeV/c^{2}]";
  om.ylabel = "Probes";
  om.showLegend = false;
  TCanvas* cm = PlotStyle::Draw({{hMass, "", "E1", kBlack}}, om);
  if (cm) {
    cm->Update();
    fMass->SetLineColor(kRed + 1);
    fMass->SetNpx(500);
    fMass->Draw("SAME");
    for (double x : {kWinLo, kWinHi, kSideLo, kSideHi}) {
      TLine* l = new TLine(x, cm->GetUymin(), x, 0.55 * cm->GetUymax());
      l->SetLineStyle(2);
      l->SetLineColor(x < 1.05 ? kAzure + 2 : kGray + 2);
      l->Draw();
    }
    cm->Print(kOut + "mass.pdf");
  }

  // <nSigmaKaon> vs p, with the TOF kaons for comparison.
  TH1D* hTof = nullptr;
  std::unique_ptr<TFile> fTof(TFile::Open("PID/nsigmakaon_mean_vs_p.root", "READ"));
  if (fTof && !fTof->IsZombie() && (hTof = fTof->Get<TH1D>("hMeanNSigmaKaon"))) hTof->SetDirectory(nullptr);
  PlotStyle::Options op;
  op.title = "Run-19 Au+Au";
  op.subtitle = "#sqrt{s_{NN}} = 200 GeV";
  op.note = "#phi #rightarrow K^{+}K^{-} tag and probe";
  op.xlabel = "p [GeV/c]";
  op.ylabel = "#LTn#sigma_{K}#GT";
  op.showEntries = false;
  op.XRange(0.1, 1.0);
  op.YRange(-0.5, 3);
  op.legendX1 = 0.55;
  std::vector<PlotStyle::Series> series = {{hMean, "#phi kaons", "E1", kRed + 1}};
  if (hTof) series.push_back({hTof, "TOF kaons (70-80%)", "E1", kGray + 2});
  TCanvas* cp = PlotStyle::Draw(series, op);
  if (cp) {
    TLine zero(0.1, 0, 1.0, 0);
    zero.SetLineStyle(3);
    zero.SetLineColor(kGray + 1);
    zero.Draw();
    cp->Print(kOut + "mean_vs_p.pdf");
  }

  // Fit in each momentum bin: window, scaled sideband, subtracted, Gaussian.
  for (int b = 0; b < nP; ++b) {
    if (!fits[b]) continue;
    TH1D* side = static_cast<TH1D*>(hSide[b]->Clone(Form("hSideScaled_pbin_%d", b)));
    side->Scale(scale);
    PlotStyle::Options o;
    o.title = Form("%.1f < p < %.1f GeV/c", pBins[b], pBins[b + 1]);
    o.subtitle = Form("#mu = %.2f #pm %.2f, #sigma = %.2f", fits[b]->GetParameter(1), fits[b]->GetParError(1),
                      std::abs(fits[b]->GetParameter(2)));
    o.xlabel = "n#sigma_{K}";
    o.ylabel = "Probes";
    o.showEntries = false;
    o.legendX1 = 0.62;
    TCanvas* cb = PlotStyle::Draw({{hWin[b], "window", "HIST", kGray + 2},
                                   {side, "sideband (scaled)", "HIST", kAzure + 1},
                                   {hSub[b], "subtracted", "E1", kBlack}}, o);
    if (!cb) continue;
    fits[b]->SetLineColor(kRed + 1);
    fits[b]->SetLineWidth(2);
    fits[b]->SetNpx(400);
    fits[b]->Draw("SAME");
    cb->Print(Form("%sfits/pbin_%02d.pdf", kOut.Data(), b));
  }
}
