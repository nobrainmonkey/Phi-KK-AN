#include "../events/PlotStyle.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TLine.h"
#include "TSystem.h"
#include "TTree.h"
#include <cmath>
#include <iostream>
#include <memory>
#include <vector>

// Mean nSigmaKaon vs momentum for TOF-identified kaons (reproduces phi-kk-results/kaon/NSigma.C).
//
// Input: KaonTree from StPhiKKMaker (peripheral 70-80% Au+Au, TOF-matched primary tracks with
// 0.22 < m^2 < 0.28 GeV^2/c^4 already applied in production). Here each track must also pass
//   |nSigmaKaon| < 5 and 0.22 < m^2 < 0.27 GeV^2/c^4.
// In each momentum bin the nSigmaKaon distribution is fitted with a Gaussian in [-5, 5]
// (if it has more than 30 entries; otherwise the histogram mean is used).
// Output: PID/nsigmakaon_mean_vs_p.root, PID/nsigmakaon_mean_vs_p.pdf, and the fit in each
// momentum bin in PID/fits/.

void kaon_tof(TString input_file = "/Users/xihehan/phi-kk-results/data/kaon_mudst.root") {
  std::unique_ptr<TFile> file(TFile::Open(input_file, "READ"));
  if (!file || file->IsZombie()) return;
  TTree* tree = file->Get<TTree>("KaonTree");
  if (!tree) {
    std::cerr << "Missing KaonTree in " << input_file << "\n";
    return;
  }

  // One entry per event; track branches are vectors.
  std::vector<double>* track_pt = nullptr;
  std::vector<double>* track_eta = nullptr;
  std::vector<double>* track_nSigmaKaon = nullptr;
  std::vector<double>* track_m2 = nullptr;
  tree->SetBranchStatus("*", false);
  for (const char* b : {"track_pt", "track_eta", "track_nSigmaKaon", "track_m2"}) tree->SetBranchStatus(b, true);
  tree->SetBranchAddress("track_pt", &track_pt);
  tree->SetBranchAddress("track_eta", &track_eta);
  tree->SetBranchAddress("track_nSigmaKaon", &track_nSigmaKaon);
  tree->SetBranchAddress("track_m2", &track_m2);

  // Momentum bins: 0.1 - 2.5 GeV/c in steps of 0.1 (edges written out, as in NSigma.C).
  const std::vector<double> pBins = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3,
                                     1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5};
  const int nP = pBins.size() - 1;
  const double pMin = pBins.front(), pMax = pBins.back();

  std::vector<TH1D*> hNSigma(nP);
  for (int b = 0; b < nP; ++b) {
    hNSigma[b] = new TH1D(Form("hNSigma_pbin_%d", b),
                          Form("nSigmaKaon, %.2f < p < %.2f;n#sigma_{K};Tracks", pBins[b], pBins[b + 1]),
                          120, -6, 6);
    hNSigma[b]->Sumw2();
  }

  const Long64_t nEntries = tree->GetEntries();
  for (Long64_t i = 0; i < nEntries; ++i) {
    tree->GetEntry(i);
    for (size_t j = 0; j < track_pt->size(); ++j) {
      const double nsk = (*track_nSigmaKaon)[j], m2 = (*track_m2)[j];
      if (std::abs(nsk) >= 5) continue;
      if (m2 <= 0.22 || m2 >= 0.27) continue;

      const double p = (*track_pt)[j] * std::cosh((*track_eta)[j]);
      for (int b = 0; b < nP; ++b) {
        if (p >= pBins[b] && p < pBins[b + 1]) {
          hNSigma[b]->Fill(nsk);
          break;
        }
      }
    }
  }

  // Mean nSigmaKaon per bin, as a graph (NSigma.C output) and as a histogram (for plotting).
  TGraphErrors* gMean = new TGraphErrors(nP);
  gMean->SetName("gMeanNSigmaKaon");
  gMean->SetTitle("Mean nSigmaKaon vs p;p (GeV/c);<nSigmaKaon>");
  TH1D* hMean = new TH1D("hMeanNSigmaKaon", "", nP, pBins.data());

  std::vector<TF1*> fits(nP, nullptr);
  for (int b = 0; b < nP; ++b) {
    double mean, err;
    if (hNSigma[b]->GetEntries() > 30) {
      fits[b] = new TF1(Form("gaus_%d", b), "gaus", -5, 5);
      hNSigma[b]->Fit(fits[b], "Q0");
      mean = fits[b]->GetParameter(1);
      err = fits[b]->GetParError(1);
    } else {
      mean = hNSigma[b]->GetMean();
      err = hNSigma[b]->GetMeanError();
    }
    const double pc = hMean->GetBinCenter(b + 1);
    gMean->SetPoint(b, pc, mean);
    gMean->SetPointError(b, 0.5 * hMean->GetBinWidth(b + 1), err);
    hMean->SetBinContent(b + 1, mean);
    hMean->SetBinError(b + 1, err);
    printf("p [%.1f, %.1f): <nSigmaKaon> = %6.3f +/- %.3f  (%.0f tracks)\n",
           pBins[b], pBins[b + 1], mean, err, hNSigma[b]->GetEntries());
  }

  // ---- Output ----------------------------------------------------------------
  gSystem->mkdir("PID", true);
  TFile out("PID/nsigmakaon_mean_vs_p.root", "RECREATE");
  gMean->Write();
  hMean->Write();
  for (TH1D* h : hNSigma) h->Write();
  out.Close();

  PlotStyle::Options opt;
  opt.title = "Run-19 Au+Au";
  opt.subtitle = "#sqrt{s_{NN}} = 200 GeV";
  opt.note = "70-80%, TOF kaons";
  opt.xlabel = "p [GeV/c]";
  opt.ylabel = "#LTn#sigma_{K}#GT";
  opt.showLegend = false;
  opt.showEntries = false;
  opt.YRange(-0.5, 2.5);
  TCanvas* c = PlotStyle::Draw({{hMean, "", "E1", kBlack}}, opt);
  if (!c) return;
  TLine zero(pMin, 0, pMax, 0);
  zero.SetLineStyle(2);
  zero.SetLineColor(kGray + 2);
  zero.Draw();
  TLatex t;
  t.SetNDC();
  t.SetTextFont(42);
  t.SetTextSize(0.032);
  t.SetTextAlign(31);
  t.DrawLatex(0.93, 0.20, "|n#sigma_{K}| < 5, 0.22 < m^{2} < 0.27 GeV^{2}/c^{4}");
  c->Print("PID/nsigmakaon_mean_vs_p.pdf");

  // ---- Gaussian fit in each momentum bin: PID/fits/nsigmakaon_fit_pbin_<b>.pdf ----
  gSystem->mkdir("PID/fits", true);
  for (int b = 0; b < nP; ++b) {
    PlotStyle::Options o;
    o.title = Form("%.1f < p < %.1f GeV/c", pBins[b], pBins[b + 1]);
    o.xlabel = "n#sigma_{K}";
    o.ylabel = "Tracks";
    o.showLegend = false;
    o.labelSize = 0.05;
    o.titleSize = 0.058;
    o.left = 0.17;
    o.yTitleOffset = 1.35;
    o.annotationX = 0.21;
    if (fits[b])
      o.subtitle = Form("#mu = %.3f #pm %.3f", fits[b]->GetParameter(1), fits[b]->GetParError(1)),
      o.note = Form("#sigma = %.3f", fits[b]->GetParameter(2));
    else
      o.subtitle = Form("mean = %.3f (no fit)", hNSigma[b]->GetMean());
    TCanvas* cb = PlotStyle::Draw({{hNSigma[b], "", "E1", kBlack}}, o);
    if (!cb) continue;
    if (fits[b]) {
      fits[b]->SetLineColor(kRed + 1);
      fits[b]->SetLineWidth(2);
      fits[b]->SetNpx(400);
      fits[b]->Draw("SAME");
    }
    cb->Print(Form("PID/fits/nsigmakaon_fit_pbin_%02d.pdf", b));
  }
}
