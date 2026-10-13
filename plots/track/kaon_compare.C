#include "../events/PlotStyle.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TMath.h"
#include "TParameter.h"
#include <cmath>
#include <iostream>
#include <memory>
#include <vector>

// nSigmaKaon of phi kaons (kaon_phi.C) and TOF kaons (kaon_tof.C) in one merged momentum bin,
// 0.3 < p < 0.6 GeV/c (bins 02-04 of both macros). Run after kaon_tof.C and kaon_phi.C.
// The phi kaons are fitted as in kaon_phi.C: binned Poisson likelihood fit of the summed window
// histogram with a Gaussian plus the scaled, summed sideband histogram as a fixed template.
// Because <nSigmaKaon> changes quickly with p inside the merged bin, the TOF mean is also given
// weighted to the momentum distribution of the phi kaons (fitted N_K in each 0.1 GeV/c bin).
// Output: PID/nsigmakaon_compare_merged.pdf.

const int kFirst = 2, kLast = 4; // momentum bins 0.3-0.4, 0.4-0.5, 0.5-0.6 GeV/c

// Same model as kaon_phi.C: N_K kaons in a normalised Gaussian + gScale * sideband histogram.
TH1D* gSide = nullptr;
double gScale = 0, gBinW = 0;
double NSigmaModel(double* x, double* p) {
  return p[0] * gBinW * TMath::Gaus(x[0], p[1], p[2], true) + gScale * gSide->GetBinContent(gSide->FindBin(x[0]));
}

void kaon_compare() {
  std::unique_ptr<TFile> fPhi(TFile::Open("PID/phi/nsigmakaon_phi_mean_vs_p.root", "READ"));
  std::unique_ptr<TFile> fTof(TFile::Open("PID/nsigmakaon_mean_vs_p.root", "READ"));
  if (!fPhi || fPhi->IsZombie() || !fTof || fTof->IsZombie()) {
    std::cerr << "Run kaon_tof.C and kaon_phi.C first.\n";
    return;
  }
  const double scale = fPhi->Get<TParameter<double>>("sidebandScale")->GetVal();
  TH1D* hNK = fPhi->Get<TH1D>("hNKaonPhi");
  TH1D* hTofMean = fTof->Get<TH1D>("hMeanNSigmaKaon");

  // Sum the three momentum bins; TOF mean weighted by the fitted phi kaons in each bin.
  TH1D *hWin = nullptr, *hSide = nullptr, *hTof = nullptr;
  double sw = 0, swmu = 0;
  for (int b = kFirst; b <= kLast; ++b) {
    TH1D* w = fPhi->Get<TH1D>(Form("hWin_pbin_%d", b));
    TH1D* s = fPhi->Get<TH1D>(Form("hSide_pbin_%d", b));
    TH1D* t = fTof->Get<TH1D>(Form("hNSigma_pbin_%d", b));
    if (!w || !s || !t) return;
    if (!hWin) {
      hWin = static_cast<TH1D*>(w->Clone("hWinMerged"));
      hSide = static_cast<TH1D*>(s->Clone("hSideMerged"));
      hTof = static_cast<TH1D*>(t->Clone("hTofMerged"));
      for (TH1D* h : {hWin, hSide, hTof}) h->SetDirectory(nullptr);
    } else {
      hWin->Add(w);
      hSide->Add(s);
      hTof->Add(t);
    }
    const double nk = hNK->GetBinContent(b + 1);
    sw += nk;
    swmu += nk * hTofMean->GetBinContent(b + 1);
  }
  const double muTofWeighted = swmu / sw;

  // phi kaons: likelihood fit of the merged window histogram.
  gSide = hSide;
  gScale = scale;
  gBinW = hWin->GetBinWidth(1);
  TF1 model("model", NSigmaModel, -5, 5, 3);
  TH1D* hSub = static_cast<TH1D*>(hWin->Clone("hSubMerged"));
  hSub->Add(hSide, -scale);
  model.SetParameters(hSub->Integral(), hSub->GetMean(), 1);
  model.SetParLimits(2, 0.3, 3);
  hWin->Fit(&model, "QLR0");
  const double nK = model.GetParameter(0), mu = model.GetParameter(1), muErr = model.GetParError(1), sig = model.GetParameter(2);

  // TOF kaons: Gaussian fit of the merged histogram (as in kaon_tof.C).
  TF1 fTofFit("fTofFit", "gaus", -5, 5);
  hTof->Fit(&fTofFit, "Q0");

  printf("0.3 < p < 0.6 GeV/c\n");
  printf("  phi kaons (likelihood fit): mu = %.3f +/- %.3f, sigma = %.3f, N_K = %.0f\n", mu, muErr, sig, nK);
  printf("  TOF kaons:                  mu = %.3f +/- %.3f, sigma = %.3f\n", fTofFit.GetParameter(1),
         fTofFit.GetParError(1), fTofFit.GetParameter(2));
  printf("  TOF weighted to phi p:      mu = %.3f\n", muTofWeighted);

  // Overlay, both normalised to unit area; the phi kaons are shown sideband subtracted.
  hTof->Rebin(2);  // 0.1 -> 0.2 in nSigmaKaon, the phi binning
  const double nTof = hTof->Integral();
  hTof->Scale(1 / nTof);
  hSub->Scale(1 / nK);
  TF1 fPhiDraw("fPhiDraw", Form("%g * TMath::Gaus(x, %g, %g, true)", gBinW, mu, sig), -5, 5);
  TF1 fTofDraw("fTofDraw", "gaus", -5, 5);
  fTofDraw.SetParameters(fTofFit.GetParameter(0) * 2 / nTof, fTofFit.GetParameter(1), fTofFit.GetParameter(2));

  PlotStyle::Options opt;
  opt.title = "Run-19 Au+Au";
  opt.subtitle = "#sqrt{s_{NN}} = 200 GeV";
  opt.note = "0.3 < p < 0.6 GeV/c";
  opt.xlabel = "n#sigma_{K}";
  opt.ylabel = "Normalised tracks";
  opt.showEntries = false;
  opt.XRange(-5, 5);
  opt.legendX1 = 0.55;
  TCanvas* c = PlotStyle::Draw({{hTof, "TOF kaons (70-80%)", "HIST", kBlack},
                                {hSub, "#phi kaons, sideband subtracted", "E1", kRed + 1}}, opt);
  if (!c) return;
  fTofDraw.SetLineColor(kBlack);
  fTofDraw.SetLineStyle(2);
  fPhiDraw.SetLineColor(kRed + 1);
  for (TF1* f : {&fTofDraw, &fPhiDraw}) {
    f->SetNpx(400);
    f->SetLineWidth(2);
    f->Draw("SAME");
  }
  TLatex t;
  t.SetNDC();
  t.SetTextFont(42);
  t.SetTextSize(0.032);
  t.SetTextAlign(13);
  t.DrawLatex(0.18, 0.66, Form("#phi: #mu = %.2f #pm %.2f", mu, muErr));
  t.DrawLatex(0.18, 0.61, Form("TOF: #mu = %.2f (%.2f weighted to #phi p)", fTofFit.GetParameter(1), muTofWeighted));
  c->Print("PID/nsigmakaon_compare_merged.pdf");
}
