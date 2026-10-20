#include "../../../events/PlotStyle.h"
#include "TEllipse.h"
#include "TFile.h"
#include "TGraph.h"
#include "TLine.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"
#include "cuts.h"
#include <fstream>

// Run here: root -l -b -q pid_plots.C
// Histograms and the small selected-pair tree go into pid_plots.root.
void pid_plots() {
  TH1::AddDirectory(false);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(kBird);
  gSystem->mkdir("figures", true);
  TFile input(pairFile);
  TTree *tree = (TTree *)input.Get("uls_data");
  TFile correction(correctionFile);
  TH1D *shift = (TH1D *)correction.Get("hMeanNSigmaKaon");
  if (!tree || !shift) {
    printf("Cannot find uls_data or hMeanNSigmaKaon.\n");
    return;
  }
  double correctionMin = shift->GetBinCenter(1);
  double correctionMax = shift->GetBinCenter(shift->GetNbinsX());

  UShort_t bbcEast, bbcWest, hits1, hits2;
  UInt_t nPrimary;
  Short_t charge1, charge2;
  Double_t vx, vy, vz, mass, pairPt, y, pt1, pt2, eta1, eta2;
  Double_t k1, k2, e1, e2, pi1, pi2, pr1, pr2, ee, pipi, pp;
  tree->SetBranchStatus("*", 0);
  tree->SetBranchStatus("uls_pair_*", 1);
  tree->SetBranchStatus("uls_phi_*", 1);
  tree->SetBranchStatus("uls_kaon*_pt", 1);
  tree->SetBranchStatus("uls_kaon*_eta", 1);
  tree->SetBranchStatus("uls_kaon*_nhits", 1);
  tree->SetBranchStatus("uls_kaon*_charge", 1);
  tree->SetBranchStatus("uls_kaon*_nsigma_*", 1);
  tree->SetBranchAddress("uls_pair_BBC_Veto_East", &bbcEast);
  tree->SetBranchAddress("uls_pair_BBC_Veto_West", &bbcWest);
  tree->SetBranchAddress("uls_pair_event_n_primary_tracks", &nPrimary);
  tree->SetBranchAddress("uls_pair_event_vx", &vx);
  tree->SetBranchAddress("uls_pair_event_vy", &vy);
  tree->SetBranchAddress("uls_pair_event_vz", &vz);
  tree->SetBranchAddress("uls_phi_mass", &mass);
  tree->SetBranchAddress("uls_phi_pt", &pairPt);
  tree->SetBranchAddress("uls_phi_rapidity", &y);
  tree->SetBranchAddress("uls_kaon1_pt", &pt1);
  tree->SetBranchAddress("uls_kaon2_pt", &pt2);
  tree->SetBranchAddress("uls_kaon1_eta", &eta1);
  tree->SetBranchAddress("uls_kaon2_eta", &eta2);
  tree->SetBranchAddress("uls_kaon1_nhits", &hits1);
  tree->SetBranchAddress("uls_kaon2_nhits", &hits2);
  tree->SetBranchAddress("uls_kaon1_charge", &charge1);
  tree->SetBranchAddress("uls_kaon2_charge", &charge2);
  tree->SetBranchAddress("uls_kaon1_nsigma_kaon", &k1);
  tree->SetBranchAddress("uls_kaon2_nsigma_kaon", &k2);
  tree->SetBranchAddress("uls_kaon1_nsigma_electron", &e1);
  tree->SetBranchAddress("uls_kaon2_nsigma_electron", &e2);
  tree->SetBranchAddress("uls_kaon1_nsigma_pion", &pi1);
  tree->SetBranchAddress("uls_kaon2_nsigma_pion", &pi2);
  tree->SetBranchAddress("uls_kaon1_nsigma_proton", &pr1);
  tree->SetBranchAddress("uls_kaon2_nsigma_proton", &pr2);
  tree->SetBranchAddress("uls_phi_x2ee", &ee);
  tree->SetBranchAddress("uls_phi_x2pipi", &pipi);
  tree->SetBranchAddress("uls_phi_x2pp", &pp);

  TFile output("pid_plots.root", "RECREATE");
  double p1, p2, kc1, kc2, kk;
  TTree pairs("pairs", "ULS after event/kinematic cuts, before PID");
  pairs.Branch("mass", &mass);
  pairs.Branch("pt", &pairPt);
  pairs.Branch("y", &y);
  pairs.Branch("p1", &p1);
  pairs.Branch("p2", &p2);
  pairs.Branch("charge1", &charge1);
  pairs.Branch("charge2", &charge2);
  pairs.Branch("kk", &kk);
  pairs.Branch("ee", &ee);
  pairs.Branch("pipi", &pipi);
  pairs.Branch("pp", &pp);
  pairs.Branch("k1_corr", &kc1);
  pairs.Branch("k2_corr", &kc2);

  // First index: 0 = no PID cuts, 1 = KK and pion cuts, before electron cuts.
  // Second index: KK, ee, pipi, pp, in that order.
  const char *names[4] = {"kk", "ee", "pipi", "pp"};
  const char *trackNames[4] = {"kaon", "electron", "pion", "proton"};
  const char *labels[4] = {"#chi^{2}_{KK} (corrected)", "#chi^{2}_{ee}",
                           "#chi^{2}_{#pi#pi}", "#chi^{2}_{pp}"};
  const char *trackLabels[4] = {"n#sigma_{K} (corrected)", "n#sigma_{e}",
                                "n#sigma_{#pi}", "n#sigma_{p}"};
  double chiMax[4] = {150, 1200, 1200, 1200};
  double nsMin[4] = {-12, -12, -6, -30};
  double nsMax[4] = {12, 32, 34, 10};
  TH1D *chi[2][4], *nsigma[2][4];
  for (int selection = 0; selection < 2; ++selection) {
    for (int species = 0; species < 4; ++species) {
      chi[selection][species] =
          new TH1D(Form("chi2_%s_%d", names[species], selection), "",
                   species == 0 ? 300 : 600, 0, chiMax[species]);
      nsigma[selection][species] =
          new TH1D(Form("nsigma_%s_%d", trackNames[species], selection), "", 48,
                   nsMin[species], nsMax[species]);
      chi[selection][species]->Sumw2();
      nsigma[selection][species]->Sumw2();
    }
  }
  TH2D *maps[3];
  TGraph all[3], kaon[3];
  for (int species = 0; species < 3; ++species) {
    maps[species] = new TH2D(Form("kk_vs_%s", names[species + 1]), "", 300, 0,
                             150, 600, 0, 1200);
    maps[species]->Sumw2();
    all[species].SetName(Form("kk_vs_%s_all", names[species + 1]));
    kaon[species].SetName(Form("kk_vs_%s_kaon", names[species + 1]));
  }
  TH1D before("mass_without_ee", "", 55, .98, 1.20);
  TH1D after("mass_with_ee", "", 55, .98, 1.20);
  TH1D kaonOnly("mass_without_pipi", "", 55, .98, 1.20);
  TH1D electronRemoved("mass_rejected_ee", "", 55, .98, 1.20);
  TH2D kaonPlane("nsigma_kaon_pair", "", 120, -12, 12, 120, -12, 12);
  // Two mass regions, before pion/electron rejection.
  TH2D *pionMaps[2];
  for (int region = 0; region < 2; ++region) {
    pionMaps[region] = new TH2D(Form("pion_mass_region_%d", region), "", 300, 0,
                                150, 600, 0, 1200);
    pionMaps[region]->Sumw2();
  }
  int kaonCounts[3][2] = {}; // KK limits 8, 10, 15; all masses and mass window.
  int pionBefore[2] = {}, pionAfter[2] = {}; // Window and upper mass region.
  int nKaonOnlyWindow = 0;
  before.Sumw2();
  after.Sumw2();
  kaonOnly.Sumw2();
  electronRemoved.Sumw2();
  kaonPlane.Sumw2();
  int nWindow = 0, nBefore = 0, nAfter = 0;
  for (Long64_t entry = 0; entry < tree->GetEntries(); ++entry) {
    tree->GetEntry(entry);
    if (bbcEast != 0 || bbcWest != 0 || nPrimary != 2)
      continue;
    if (std::abs(vz) > vzMax || std::hypot(vx, vy) > vrMax)
      continue;
    if (hits1 < hitsMin || hits2 < hitsMin)
      continue;
    if (pairPt <= pairPtMin || pairPt >= pairPtMax || std::abs(y) >= pairYMax)
      continue;

    p1 = pt1 * std::cosh(eta1);
    p2 = pt2 * std::cosh(eta2);
    double cp1 = p1, cp2 = p2;
    if (cp1 < correctionMin)
      cp1 = correctionMin;
    if (cp1 > correctionMax)
      cp1 = correctionMax;
    if (cp2 < correctionMin)
      cp2 = correctionMin;
    if (cp2 > correctionMax)
      cp2 = correctionMax;
    kc1 = k1 - shift->Interpolate(cp1);
    kc2 = k2 - shift->Interpolate(cp2);
    kk = kc1 * kc1 + kc2 * kc2;
    pairs.Fill();

    bool passKaonPion = kk <= kkMax && pipi >= pipiMin;
    bool passElectron = ee >= eeMin && ee >= eeFactor * kk;
    bool inWindow = mass > massMin && mass < massMax;
    int region = inWindow ? 0 : (mass >= massMax ? 1 : -1);
    for (int choice = 0; choice < 3; ++choice) {
      double limit = choice == 0 ? 8 : (choice == 1 ? 10 : 15);
      if (kk <= limit) {
        ++kaonCounts[choice][0];
        if (inWindow)
          ++kaonCounts[choice][1];
      }
    }
    if (region >= 0) {
      pionMaps[region]->Fill(kk, pipi);
      if (kk <= kkMax) {
        ++pionBefore[region];
        if (pipi >= pipiMin)
          ++pionAfter[region];
      }
    }
    if (kk <= kkMax) {
      kaonOnly.Fill(mass);
      if (inWindow)
        ++nKaonOnlyWindow;
    }
    if (passKaonPion) {
      before.Fill(mass);
      if (passElectron)
        after.Fill(mass);
      else
        electronRemoved.Fill(mass);
    }
    if (mass <= massMin || mass >= massMax)
      continue;
    ++nWindow;
    kaonPlane.Fill(kc1, kc2);
    if (passKaonPion) {
      ++nBefore;
      if (passElectron)
        ++nAfter;
    }
    double x2[4] = {kk, ee, pipi, pp};
    double ns1[4] = {kc1, e1, pi1, pr1};
    double ns2[4] = {kc2, e2, pi2, pr2};
    for (int species = 0; species < 4; ++species) {
      chi[0][species]->Fill(x2[species]);
      // If both tracks enter the same bin, the pair contributes 2 counts:
      // variance 4, rather than the variance 2 of two independent tracks.
      int b1 = nsigma[0][species]->FindBin(ns1[species]);
      int b2 = nsigma[0][species]->FindBin(ns2[species]);
      nsigma[0][species]->Fill(ns1[species]);
      nsigma[0][species]->Fill(ns2[species]);
      if (b1 == b2) {
        double error = nsigma[0][species]->GetBinError(b1);
        nsigma[0][species]->SetBinError(b1, std::sqrt(error * error + 2));
      }
      if (!passKaonPion)
        continue;
      chi[1][species]->Fill(x2[species]);
      nsigma[1][species]->Fill(ns1[species]);
      nsigma[1][species]->Fill(ns2[species]);
      if (b1 == b2) {
        double error = nsigma[1][species]->GetBinError(b1);
        nsigma[1][species]->SetBinError(b1, std::sqrt(error * error + 2));
      }
    }
    for (int species = 0; species < 3; ++species) {
      maps[species]->Fill(kk, x2[species + 1]);
      all[species].SetPoint(all[species].GetN(), kk, x2[species + 1]);
      if (passKaonPion)
        kaon[species].SetPoint(kaon[species].GetN(), kk, x2[species + 1]);
    }
  }
  output.cd();
  pairs.Write();
  before.Write();
  after.Write();
  kaonOnly.Write();
  electronRemoved.Write();
  kaonPlane.Write();

  PlotStyle::Options opt;
  opt.title = "Run-19 Au+Au";
  opt.subtitle = "#sqrt{s_{NN}} = 200 GeV";
  opt.note = "ZDCMB; incoherent region";
  opt.showEntries = false;
  opt.headroom = 2.1;
  opt.legendX1 = .53;
  opt.legendTextSize = .027;
  TString baseLabel =
      Form("ULS, #chi^{2}_{KK} #leq %.0f, #chi^{2}_{#pi#pi} #geq %.0f", kkMax,
           pipiMin);
  TLatex selectionText;
  selectionText.SetNDC();
  selectionText.SetTextFont(42);
  selectionText.SetTextSize(.029);
  double chiZoom[4] = {20, 50, 100, 150};
  for (int species = 0; species < 4; ++species) {
    TCanvas comparison(Form("chi2_%s_canvas", names[species]), "", 1800, 760);
    comparison.Divide(2, 1, .002, .002);
    for (int panel = 0; panel < 2; ++panel) {
      TH1D *h0 = (TH1D *)chi[0][species]->Clone("display_all");
      TH1D *h1 = (TH1D *)chi[1][species]->Clone("display_kaon");
      if (panel == 0) {
        h0->Rebin(species == 0 ? 5 : 10);
        h1->Rebin(species == 0 ? 5 : 10);
      }
      opt.xlabel = labels[species];
      opt.ylabel = Form("Pairs / %.1f", h0->GetBinWidth(1));
      opt.XRange(0, panel == 0 ? chiMax[species] : chiZoom[species]);
      TCanvas *c = PlotStyle::Draw(
          {{h0, "ULS, before final PID", "E1"}, {h1, baseLabel, "E1"}}, opt);
      selectionText.DrawLatex(.18, .75,
                              "0.3 < p_{T}^{KK} < 0.9 GeV/c; |y_{KK}| < 0.5");
      selectionText.DrawLatex(
          .18, .69, "1.00 < M_{KK} < 1.04 GeV/c^{2}; no electron rejection");
      if (panel == 0)
        c->Print(Form("figures/chi2_%s_full.pdf", names[species]));
      comparison.cd(panel + 1);
      c->DrawClonePad();
      delete c;
      delete h0;
      delete h1;
    }
    comparison.Print(Form("figures/chi2_%s.png", names[species]));
    comparison.Print(Form("figures/chi2_%s.pdf", names[species]));
    opt.hasXRange = false;
    opt.xlabel = trackLabels[species];
    opt.ylabel = Form("Tracks / %.2f", nsigma[0][species]->GetBinWidth(1));
    TCanvas *c =
        PlotStyle::Draw({{nsigma[0][species], "ULS, before final PID", "E1"},
                         {nsigma[1][species], baseLabel, "E1"}},
                        opt);
    selectionText.DrawLatex(.18, .75,
                            "0.3 < p_{T}^{KK} < 0.9 GeV/c; |y_{KK}| < 0.5");
    selectionText.DrawLatex(
        .18, .69, "1.00 < M_{KK} < 1.04 GeV/c^{2}; no electron rejection");
    c->Print(Form("figures/nsigma_%s.png", trackNames[species]));
    c->Print(Form("figures/nsigma_%s.pdf", trackNames[species]));
    delete c;
    for (int selection = 0; selection < 2; ++selection) {
      chi[selection][species]->Write();
      nsigma[selection][species]->Write();
    }
  }
  // Corrected kaon response of the two tracks; the nominal pair cut is a
  // circle.
  {
    TCanvas c("kaon_pair_canvas", "", 1800, 820);
    c.Divide(2, 1, .002, .002);
    for (int panel = 0; panel < 2; ++panel) {
      c.cd(panel + 1);
      gPad->SetLeftMargin(.14);
      gPad->SetRightMargin(.15);
      gPad->SetBottomMargin(.13);
      gPad->SetTopMargin(.055);
      gPad->SetTicks();
      double extent = panel == 0 ? 12 : 4.5;
      kaonPlane.GetXaxis()->SetRangeUser(-extent, extent);
      kaonPlane.GetYaxis()->SetRangeUser(-extent, extent);
      PlotStyle::StyleAxis(kaonPlane.GetXaxis(), opt);
      PlotStyle::StyleAxis(kaonPlane.GetYaxis(), opt);
      PlotStyle::StyleAxis(kaonPlane.GetZaxis(), opt);
      kaonPlane.GetXaxis()->SetTitle("n#sigma_{K,1}^{corr}");
      kaonPlane.GetYaxis()->SetTitle("n#sigma_{K,2}^{corr}");
      kaonPlane.GetYaxis()->SetTitleOffset(1.35);
      kaonPlane.DrawCopy("COLZ");
      TObject *boundaries[2];
      for (int boundary = 0; boundary < 2; ++boundary) {
        double limit = boundary == 0 ? kkMax : 15;
        TEllipse circle(0, 0, std::sqrt(limit), std::sqrt(limit));
        circle.SetFillStyle(0);
        circle.SetLineColor(boundary == 0 ? kRed + 1 : kBlack);
        circle.SetLineStyle(boundary == 0 ? 1 : 2);
        circle.SetLineWidth(2);
        boundaries[boundary] = circle.DrawClone();
      }
      // The full-range panel carries the shared dataset and boundary legend.
      if (panel == 0) {
        TLatex text;
        text.SetNDC();
        text.SetTextFont(62);
        text.SetTextSize(.035);
        text.DrawLatex(.18, .90, "Run-19 Au+Au");
        text.SetTextFont(42);
        text.SetTextSize(.029);
        text.DrawLatex(.18, .85, "#sqrt{s_{NN}} = 200 GeV");
        text.DrawLatex(.18, .80, "ZDCMB; ULS");
        TLegend legend(.48, .67, .83, .925);
        legend.SetFillStyle(0);
        legend.SetBorderSize(0);
        legend.SetTextFont(42);
        legend.SetTextSize(.025);
        legend.AddEntry((TObject *)0, "1.00 < M_{KK} < 1.04 GeV/c^{2}", "");
        legend.AddEntry((TObject *)0, "0.3 < p_{T}^{KK} < 0.9 GeV/c", "");
        legend.AddEntry((TObject *)0,
                        "|y_{KK}| < 0.5; no extra track-p_{T} cut", "");
        legend.AddEntry((TObject *)0, "Before final PID", "");
        legend.AddEntry(boundaries[0],
                        Form("#chi^{2}_{KK} = %.0f (nominal)", kkMax), "l");
        legend.AddEntry(boundaries[1], "#chi^{2}_{KK} = 15 (comparison)", "l");
        legend.DrawClone();
      }
    }
    c.Print("figures/nsigma_kaon_pair.png");
    c.Print("figures/nsigma_kaon_pair.pdf");
    kaonPlane.GetXaxis()->SetRange(0, 0);
    kaonPlane.GetYaxis()->SetRange(0, 0);
  }
  // Show explicitly where the pion rejection acts in the two mass regions.
  {
    TCanvas c("pion_regions_canvas", "", 1800, 850);
    c.Divide(2, 1, .002, .002);
    for (int region = 0; region < 2; ++region) {
      c.cd(region + 1);
      gPad->SetLeftMargin(.14);
      gPad->SetRightMargin(.15);
      gPad->SetBottomMargin(.13);
      gPad->SetTopMargin(.055);
      gPad->SetTicks();
      pionMaps[region]->GetXaxis()->SetRangeUser(0, 20);
      pionMaps[region]->GetYaxis()->SetRangeUser(0, 80);
      PlotStyle::StyleAxis(pionMaps[region]->GetXaxis(), opt);
      PlotStyle::StyleAxis(pionMaps[region]->GetYaxis(), opt);
      PlotStyle::StyleAxis(pionMaps[region]->GetZaxis(), opt);
      pionMaps[region]->GetXaxis()->SetTitle(labels[0]);
      pionMaps[region]->GetYaxis()->SetTitle(labels[2]);
      pionMaps[region]->GetYaxis()->SetTitleOffset(1.4);
      pionMaps[region]->Draw("COLZ");
      TLine kcut(kkMax, 0, kkMax, 80), picut(0, pipiMin, 20, pipiMin);
      kcut.SetLineStyle(2);
      picut.SetLineStyle(2);
      picut.SetLineColor(kRed + 1);
      TObject *kline = kcut.DrawClone();
      TObject *pionLine = picut.DrawClone();
      TLatex text;
      text.SetNDC();
      text.SetTextFont(62);
      text.SetTextSize(.035);
      text.DrawLatex(.18, .90, "Run-19 Au+Au");
      text.SetTextFont(42);
      text.SetTextSize(.029);
      text.DrawLatex(.18, .85, "#sqrt{s_{NN}} = 200 GeV");
      text.DrawLatex(.18, .80, "ZDCMB; ULS");
      TLegend legend(.48, .68, .83, .925);
      legend.SetFillStyle(0);
      legend.SetBorderSize(0);
      legend.SetTextFont(42);
      legend.SetTextSize(.025);
      legend.AddEntry((TObject *)0,
                      region == 0 ? "1.00 < M_{KK} < 1.04 GeV/c^{2}"
                                  : "1.04 #leq M_{KK} < 1.20 GeV/c^{2}",
                      "");
      legend.AddEntry((TObject *)0, "0.3 < p_{T}^{KK} < 0.9 GeV/c", "");
      legend.AddEntry((TObject *)0, "|y_{KK}| < 0.5; no extra track-p_{T} cut",
                      "");
      legend.AddEntry(kline, Form("#chi^{2}_{KK} = %.0f", kkMax), "l");
      legend.AddEntry(pionLine, Form("#chi^{2}_{#pi#pi} = %.0f", pipiMin), "l");
      legend.DrawClone();
    }
    c.Print("figures/kk_vs_pipi_massregions.png");
    c.Print("figures/kk_vs_pipi_massregions.pdf");
    for (int region = 0; region < 2; ++region) {
      pionMaps[region]->GetXaxis()->SetRange(0, 0);
      pionMaps[region]->GetYaxis()->SetRange(0, 0);
      pionMaps[region]->Write();
    }
  }
  // Same two-dimensional plot twice: full range and the cut region.
  double zoomMax[3] = {50, 100, 150};
  for (int species = 0; species < 3; ++species) {
    TCanvas c(Form("plane_%d", species), "", 1800, 760);
    c.Divide(2, 1, .002, .002);
    for (int panel = 0; panel < 2; ++panel) {
      c.cd(panel + 1);
      gPad->SetLeftMargin(.14);
      gPad->SetRightMargin(.15);
      gPad->SetBottomMargin(.13);
      gPad->SetTopMargin(.055);
      gPad->SetTicks();
      double xmax = panel == 0 ? 150 : 15;
      double ymax = panel == 0 ? 1200 : zoomMax[species];
      TH2D *display = (TH2D *)maps[species]->Clone("plane_display");
      if (panel == 0)
        display->Rebin2D(5, 10);
      display->GetXaxis()->SetRangeUser(0, xmax);
      display->GetYaxis()->SetRangeUser(0, ymax);
      PlotStyle::StyleAxis(display->GetXaxis(), opt);
      PlotStyle::StyleAxis(display->GetYaxis(), opt);
      PlotStyle::StyleAxis(display->GetZaxis(), opt);
      display->GetXaxis()->SetTitle(labels[0]);
      display->GetYaxis()->SetTitle(labels[species + 1]);
      display->GetXaxis()->SetTitleOffset(1.15);
      display->GetYaxis()->SetTitleOffset(1.4);
      display->DrawCopy("COLZ");
      delete display;
      TLine cutK(kkMax, 0, kkMax, ymax);
      cutK.SetLineStyle(2);
      cutK.SetLineColor(kGray + 2);
      cutK.SetLineWidth(2);
      TObject *kline = cutK.DrawClone();
      TLegend legend(.46, .65, .83, .925);
      legend.SetFillStyle(0);
      legend.SetBorderSize(0);
      legend.SetTextFont(42);
      legend.SetTextSize(.025);
      legend.AddEntry((TObject *)0, "1.00 < M_{KK} < 1.04 GeV/c^{2}", "");
      legend.AddEntry((TObject *)0, "0.3 < p_{T}^{KK} < 0.9 GeV/c", "");
      legend.AddEntry((TObject *)0, "|y_{KK}| < 0.5; no extra track-p_{T} cut",
                      "");
      legend.AddEntry((TObject *)0, "Before final PID", "");
      legend.AddEntry(kline, Form("#chi^{2}_{KK} = %.0f", kkMax), "l");
      if (species == 0) {
        TLine cutE(0, eeMin, xmax, eeMin);
        TLine ratio(0, 0, xmax, eeFactor * xmax);
        cutE.SetLineStyle(2);
        ratio.SetLineStyle(3);
        cutE.SetLineColor(kRed + 1);
        ratio.SetLineColor(kRed + 1);
        cutE.SetLineWidth(2);
        ratio.SetLineWidth(2);
        TObject *eline = cutE.DrawClone();
        TObject *ratioLine = ratio.DrawClone();
        legend.AddEntry(eline, Form("#chi^{2}_{ee} = %.0f", eeMin), "l");
        legend.AddEntry(ratioLine,
                        Form("#chi^{2}_{ee} = %.0f#chi^{2}_{KK}", eeFactor),
                        "l");
      }
      if (species == 1) {
        TLine cutPi(0, pipiMin, xmax, pipiMin);
        cutPi.SetLineStyle(2);
        cutPi.SetLineColor(kRed + 1);
        TObject *pionLine = cutPi.DrawClone();
        legend.AddEntry(pionLine, Form("#chi^{2}_{#pi#pi} = %.0f", pipiMin),
                        "l");
      }
      if (panel == 0) {
        TLatex text;
        text.SetNDC();
        text.SetTextFont(62);
        text.SetTextSize(.035);
        text.DrawLatex(.18, .90, "Run-19 Au+Au");
        text.SetTextFont(42);
        text.SetTextSize(.029);
        text.DrawLatex(.18, .85, "#sqrt{s_{NN}} = 200 GeV");
        text.DrawLatex(.18, .80, "ZDCMB; ULS");
        legend.DrawClone();
      }
    }
    c.Print(Form("figures/kk_vs_%s.png", names[species + 1]));
    c.Print(Form("figures/kk_vs_%s.pdf", names[species + 1]));
    maps[species]->GetXaxis()->SetRange(0, 0);
    maps[species]->GetYaxis()->SetRange(0, 0);
    maps[species]->Write();
    all[species].Write();
    kaon[species].Write();
  }
  // Raw window yield, with and without BOTH electron requirements.
  int lost = nBefore - nAfter;
  double fraction = nBefore ? double(lost) / nBefore : 0;
  double fractionError =
      nBefore ? std::sqrt(fraction * (1 - fraction) / nBefore) : 0;
  opt.xlabel = "M_{KK} [GeV/c^{2}]";
  opt.ylabel = "Pairs / 4 MeV/c^{2}";
  opt.headroom = 2.2;
  TCanvas *c = PlotStyle::Draw(
      {{&before, "ULS, no electron rejection", "E1"},
       {&after,
        Form("ULS, #chi^{2}_{ee} #geq %.0f and #geq %.0f#chi^{2}_{KK}", eeMin,
             eeFactor),
        "E1"}},
      opt);
  TLatex text;
  text.SetNDC();
  text.SetTextFont(42);
  text.SetTextSize(.031);
  text.DrawLatex(.18, .75,
                 Form("0.3 < p_{T}^{KK} < 0.9 GeV/c; |y_{KK}| < 0.5; "
                      "#chi^{2}_{KK} #leq %.0f",
                      kkMax));
  text.DrawLatex(
      .18, .69,
      Form("#chi^{2}_{#pi#pi} #geq %.0f; no additional track-p_{T} cut",
           pipiMin));
  text.DrawLatex(.18, .63,
                 Form("%.2f < M_{KK} < %.2f: %d #rightarrow %d pairs", massMin,
                      massMax, nBefore, nAfter));
  c->Update();
  TLine left(massMin, 0, massMin, .50 * c->GetUymax());
  TLine right(massMax, 0, massMax, .50 * c->GetUymax());
  left.SetLineStyle(2);
  right.SetLineStyle(2);
  left.DrawClone();
  right.DrawClone();
  c->Print("figures/mass_electron_cuts.png");
  c->Print("figures/mass_electron_cuts.pdf");
  delete c;
  c = PlotStyle::Draw(
      {{&kaonOnly, "ULS, no pion rejection", "E1"},
       {&before, Form("ULS, #chi^{2}_{#pi#pi} #geq %.0f", pipiMin), "E1"}},
      opt);
  text.DrawLatex(.18, .75,
                 Form("#chi^{2}_{KK} #leq %.0f; no electron rejection", kkMax));
  text.DrawLatex(
      .18, .69,
      Form("Mass window: %d #rightarrow %d pairs", nKaonOnlyWindow, nBefore));
  text.DrawLatex(.18, .63,
                 Form("M_{KK} #geq %.2f: %d #rightarrow %d pairs", massMax,
                      pionBefore[1], pionAfter[1]));
  c->Update();
  left.SetY2(.50 * c->GetUymax());
  right.SetY2(.50 * c->GetUymax());
  left.DrawClone();
  right.DrawClone();
  c->Print("figures/mass_pion_cut.png");
  c->Print("figures/mass_pion_cut.pdf");
  delete c;
  c = PlotStyle::Draw({{&electronRemoved,
                        Form("ULS, #chi^{2}_{ee} < %.0f or < %.0f#chi^{2}_{KK}",
                             eeMin, eeFactor),
                        "E1"}},
                      opt);
  text.DrawLatex(.18, .75,
                 Form("0.3 < p_{T}^{KK} < 0.9 GeV/c; |y_{KK}| < 0.5; "
                      "#chi^{2}_{KK} #leq %.0f",
                      kkMax));
  text.DrawLatex(
      .18, .69,
      Form("#chi^{2}_{#pi#pi} #geq %.0f; no additional track-p_{T} cut",
           pipiMin));
  text.DrawLatex(
      .18, .63,
      Form("Rejected in 1.00 < M_{KK} < 1.04 GeV/c^{2}: %d pairs", lost));
  c->Update();
  left.SetY2(.50 * c->GetUymax());
  right.SetY2(.50 * c->GetUymax());
  left.DrawClone();
  right.DrawClone();
  c->Print("figures/mass_rejected_electron.png");
  c->Print("figures/mass_rejected_electron.pdf");
  delete c;
  TH1D counts("window_yields", "", 3, 0, 3);
  counts.GetXaxis()->SetBinLabel(1, "without ee");
  counts.GetXaxis()->SetBinLabel(2, "with ee");
  counts.GetXaxis()->SetBinLabel(3, "removed");
  counts.SetBinContent(1, nBefore);
  counts.SetBinContent(2, nAfter);
  counts.SetBinContent(3, lost);
  for (int b = 1; b <= 3; ++b)
    counts.SetBinError(b, std::sqrt(counts.GetBinContent(b)));
  counts.Write();
  output.Close();
  std::ofstream result("window_yields.txt");
  result << "Raw ULS counts, " << massMin << " < M_KK < " << massMax
         << " GeV/c^2\n"
         << "Before PID: " << nWindow << "\n"
         << "KK <= " << kkMax << ", pipi >= " << pipiMin << "\n"
         << "Without electron cuts: " << nBefore << "\n"
         << "With ee >= " << eeMin << " and ee >= " << eeFactor
         << "*KK: " << nAfter << "\nRemoved: " << lost << " +/- "
         << std::sqrt(double(lost)) << "\nLoss: " << 100 * fraction << " +/- "
         << 100 * fractionError << "%\n";
  printf("Before PID: %d; without ee: %d; with ee: %d; removed: %d\n", nWindow,
         nBefore, nAfter, lost);
  result << "\nKaon cut only (no pion/electron cuts):\n";
  for (int choice = 0; choice < 3; ++choice)
    result << "KK <= " << (choice == 0 ? 8 : (choice == 1 ? 10 : 15)) << ": "
           << kaonCounts[choice][0] << " all masses, " << kaonCounts[choice][1]
           << " window\n";
  result << "Pion cut at nominal KK: window " << pionBefore[0] << " -> "
         << pionAfter[0] << "; upper mass region " << pionBefore[1] << " -> "
         << pionAfter[1] << "\n";
}
