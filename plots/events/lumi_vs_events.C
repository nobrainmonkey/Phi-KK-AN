#include "TAxis.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TProfile.h"
#include "TVirtualPad.h"
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

// Integrated luminosity per run (lumi_by_run.txt) vs. number of trigger-700001
// events per run in the MuDST (hLumiBeforeBadRunEventsByRun: includes runs that
// the GoodRun list drops; those are drawn as open markers).
// Second figure: live ZDC coincidences per run, <R_ZDC> x (t_end - t_start) x
// livetime (t_run = t_end - t_start), from pLumiBeforeBadRunZDCCoincidenceRateVsRun (event-averaged
// runInfo rate) and the table's start/end times [s], vs. the table L_int.
namespace {
struct RunLumi { int run; double lumi; double lumiEvents; double duration; double livetime; };

void StyleAxes(TAxis* x, TAxis* y) {
  for (TAxis* a : {x, y}) {
    a->SetLabelFont(42); a->SetTitleFont(42);
    a->SetLabelSize(0.040); a->SetTitleSize(0.047);
    a->SetTickLength(0.025); a->SetNdivisions(505);
  }
  x->SetTitleOffset(1.15); y->SetTitleOffset(1.40);
}

void StylePad(TVirtualPad* p, double bottom) {
  p->SetMargin(0.14, 0.045, bottom, 0.055);
  p->SetTicks(1, 1);
}

void Annotate() {
  TLatex t; t.SetNDC(); t.SetTextFont(42);
  t.SetTextSize(0.040); t.DrawLatex(0.18, 0.88, "Run-19 Au+Au");
  t.SetTextSize(0.034); t.DrawLatex(0.18, 0.83, "#sqrt{s_{NN}} = 200 GeV");
  t.DrawLatex(0.18, 0.78, "ZDCMB");
}
}

void lumi_vs_events(TString input_file = "../../data/phi-kk-mudst.root",
                    TString lumi_file = "lumi_by_run.txt",
                    TString lumi_unit = "#mub^{-1}") {
  std::vector<RunLumi> table;
  std::ifstream in(lumi_file.Data());
  std::string line;
  while (std::getline(in, line)) {
    if (line.empty() || line[0] == '#') continue;
    std::istringstream ss(line);
    RunLumi r; double t0, t1, prescale, live2; int fill; std::string det;
    if (ss >> r.run >> t0 >> t1 >> fill >> r.lumi >> prescale >> r.livetime >> det >> live2 >> r.lumiEvents) {
      r.duration = t1 - t0;
      table.push_back(r);
    }
  }
  if (table.empty()) { std::cerr << "No runs read from " << lumi_file << std::endl; return; }

  std::unique_ptr<TFile> file(TFile::Open(input_file, "READ"));
  if (!file || file->IsZombie()) return;
  TH1* hAll = dynamic_cast<TH1*>(file->Get("hLumiBeforeBadRunEventsByRun"));
  TH1* hGood = dynamic_cast<TH1*>(file->Get("hLumiGoodRunEventsByRun"));
  TProfile* pZdc = dynamic_cast<TProfile*>(file->Get("pLumiBeforeBadRunZDCCoincidenceRateVsRun"));
  if (!hAll || !hGood || !pZdc) { std::cerr << "Missing per-run event histograms" << std::endl; return; }

  TGraph* gGood = new TGraph();
  TGraph* gDropped = new TGraph();
  TGraph* gRatio = new TGraph();
  TGraph* gRatioDropped = new TGraph();
  printf("%9s %12s %10s %10s %8s %12s\n", "run", "L_int", "N_mudst", "N_lumi", "good", "N/L");
  for (std::size_t i = 0; i < table.size(); ++i) {
    const RunLumi& r = table[i];
    const double n = hAll->GetBinContent(hAll->FindFixBin(r.run));
    const bool good = hGood->GetBinContent(hGood->FindFixBin(r.run)) > 0;
    const double ratio = r.lumi > 0 ? n / r.lumi : 0;
    printf("%9d %12.6g %10.0f %10.0f %8s %12.4g\n", r.run, r.lumi, n, r.lumiEvents,
           good ? "yes" : "no", ratio);
    TGraph* g = good ? gGood : gDropped;
    g->SetPoint(g->GetN(), r.lumi, n);
    TGraph* gr = good ? gRatio : gRatioDropped;
    gr->SetPoint(gr->GetN(), i, ratio);
  }

  const int blue = TColor::GetColor("#0072B2");
  const int orange = TColor::GetColor("#D55E00");
  for (TGraph* g : {gGood, gRatio}) { g->SetMarkerStyle(20); g->SetMarkerSize(1.4); g->SetMarkerColor(blue); }
  for (TGraph* g : {gDropped, gRatioDropped}) { g->SetMarkerStyle(24); g->SetMarkerSize(1.6); g->SetMarkerColor(orange); g->SetLineColor(orange); }

  // --- N_events vs L_int, with a proportional fit to GoodRun points ---
  TCanvas* c = new TCanvas("cLumiVsEvents", "cLumiVsEvents", 2200, 1000);
  c->Divide(2, 1, 0.001, 0.001);
  TVirtualPad* c1 = c->cd(1);
  StylePad(c1, 0.13);
  double lmax = 0, nmax = 0;
  for (const RunLumi& r : table) lmax = std::max(lmax, r.lumi);
  for (int i = 1; i <= hAll->GetNbinsX(); ++i) nmax = std::max(nmax, hAll->GetBinContent(i));
  TH1* frame = c1->DrawFrame(0, 0, 1.1 * lmax, 1.45 * nmax);
  frame->GetXaxis()->SetTitle(TString::Format("Integrated luminosity [%s]", lumi_unit.Data()));
  frame->GetYaxis()->SetTitle("Events per run (MuDST)");
  StyleAxes(frame->GetXaxis(), frame->GetYaxis());

  TF1* fit = new TF1("fLumiProp", "[0]*x", 0, 1.1 * lmax);
  fit->SetLineColor(kGray + 2); fit->SetLineStyle(2); fit->SetLineWidth(2);
  gGood->Fit(fit, "QN");
  fit->Draw("SAME");
  gGood->Draw("P SAME");
  if (gDropped->GetN() > 0) gDropped->Draw("P SAME");

  TLegend* leg = new TLegend(0.50, 0.74, 0.93, 0.925);
  leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextFont(42); leg->SetTextSize(0.033);
  leg->AddEntry(gGood, TString::Format("Good runs (%d)", gGood->GetN()), "P");
  if (gDropped->GetN() > 0)
    leg->AddEntry(gDropped, TString::Format("Not in GoodRun list (%d)", gDropped->GetN()), "P");
  leg->AddEntry(fit, TString::Format("N = %.4g #times L", fit->GetParameter(0)), "L");
  leg->Draw();
  Annotate();

  // --- Events per unit luminosity, run by run ---
  TVirtualPad* c2 = c->cd(2);
  StylePad(c2, 0.20); // Room for vertical run labels.
  double rmax = 0;
  for (TGraph* g : {gRatio, gRatioDropped})
    for (int i = 0; i < g->GetN(); ++i) rmax = std::max(rmax, g->GetY()[i]);
  TH1* frame2 = c2->DrawFrame(-0.5, 0, table.size() - 0.5, 1.45 * rmax);
  TString inverse_unit = lumi_unit.EndsWith("^{-1}")
    ? TString(lumi_unit(0, lumi_unit.Length() - 5)) : "(" + lumi_unit + ")^{-1}";
  frame2->GetYaxis()->SetTitle(TString::Format("Events / L_{int} [%s]", inverse_unit.Data()));
  StyleAxes(frame2->GetXaxis(), frame2->GetYaxis());
  frame2->GetXaxis()->SetNdivisions(table.size(), false);
  frame2->GetXaxis()->SetLabelSize(0);
  TLatex lab; lab.SetTextFont(42); lab.SetTextSize(0.022); lab.SetTextAngle(90); lab.SetTextAlign(32);
  for (std::size_t i = 0; i < table.size(); ++i)
    lab.DrawLatex(i, -0.02 * 1.45 * rmax, TString::Format("%d", table[i].run));
  TF1* mean = new TF1("fMean", "[0]", -0.5, table.size() - 0.5);
  mean->SetParameter(0, fit->GetParameter(0));
  mean->SetLineColor(kGray + 2); mean->SetLineStyle(2); mean->SetLineWidth(2);
  mean->Draw("SAME");
  gRatio->Draw("P SAME");
  if (gRatioDropped->GetN() > 0) gRatioDropped->Draw("P SAME");
  TLegend* leg2 = new TLegend(0.50, 0.74, 0.93, 0.925);
  leg2->SetBorderSize(0); leg2->SetFillStyle(0); leg2->SetTextFont(42); leg2->SetTextSize(0.033);
  leg2->AddEntry(gRatio, "Good runs", "P");
  if (gRatioDropped->GetN() > 0) leg2->AddEntry(gRatioDropped, "Not in GoodRun list", "P");
  leg2->AddEntry(mean, "Fit slope (good runs)", "L");
  leg2->Draw();
  Annotate();
  c->Print("lumi_vs_events.png");
  c->Print("lumi_vs_events.pdf");

  // --- ZDC coincidence rate x live time vs. table L_int ---
  // Runs below minLumi (20191015, 20193019: ~no events) have an ill-defined
  // event-averaged rate and are left out of the cross-section panel.
  const double minLumi = 0.01;
  TGraph* zGood = new TGraph();
  TGraph* zDropped = new TGraph();
  TGraph* sGood = new TGraph();
  TGraph* sDropped = new TGraph();
  int nSkipped = 0;
  printf("\n%9s %10s %8s %10s %12s %12s\n", "run", "<R_ZDC>", "dt[s]", "livetime", "N_ZDC,live", "sigma[b]");
  for (const RunLumi& r : table) {
    const int bin = pZdc->FindFixBin(r.run);
    const double rate = pZdc->GetBinContent(bin);
    const double nZdc = rate * r.duration * r.livetime;
    const double sigma = r.lumi > 0 ? nZdc / r.lumi * 1e-6 : 0; // ub -> b
    const bool good = hGood->GetBinContent(hGood->FindFixBin(r.run)) > 0;
    printf("%9d %10.1f %8.0f %10.4f %12.4g %12.4g\n", r.run, rate, r.duration, r.livetime, nZdc, sigma);
    TGraph* g = good ? zGood : zDropped;
    g->SetPoint(g->GetN(), r.lumi, nZdc);
    if (r.lumi < minLumi || pZdc->GetBinEntries(bin) <= 0) { ++nSkipped; continue; }
    TGraph* gs = good ? sGood : sDropped;
    gs->SetPoint(gs->GetN(), rate, sigma);
  }
  for (TGraph* g : {zGood, sGood}) { g->SetMarkerStyle(20); g->SetMarkerSize(1.4); g->SetMarkerColor(blue); }
  for (TGraph* g : {zDropped, sDropped}) { g->SetMarkerStyle(24); g->SetMarkerSize(1.6); g->SetMarkerColor(orange); }

  TCanvas* cz = new TCanvas("cZdcVsLumi", "cZdcVsLumi", 2200, 1000);
  cz->Divide(2, 1, 0.001, 0.001);
  TVirtualPad* z1 = cz->cd(1);
  StylePad(z1, 0.13);
  double zmax = 0;
  for (TGraph* g : {zGood, zDropped})
    for (int i = 0; i < g->GetN(); ++i) zmax = std::max(zmax, g->GetY()[i]);
  TH1* zframe = z1->DrawFrame(0, 0, 1.1 * lmax, 1.45 * zmax);
  zframe->GetXaxis()->SetTitle(TString::Format("Integrated luminosity [%s]", lumi_unit.Data()));
  zframe->GetYaxis()->SetTitle("#LT R_{ZDC} #GT #times t_{run} #times livetime");
  StyleAxes(zframe->GetXaxis(), zframe->GetYaxis());
  TF1* zfit = new TF1("fZdcProp", "[0]*x", 0, 1.1 * lmax);
  zfit->SetLineColor(kGray + 2); zfit->SetLineStyle(2); zfit->SetLineWidth(2);
  zGood->Fit(zfit, "QN");
  zfit->Draw("SAME");
  zGood->Draw("P SAME");
  if (zDropped->GetN() > 0) zDropped->Draw("P SAME");
  TLegend* zleg = new TLegend(0.50, 0.74, 0.93, 0.925);
  zleg->SetBorderSize(0); zleg->SetFillStyle(0); zleg->SetTextFont(42); zleg->SetTextSize(0.033);
  zleg->AddEntry(zGood, TString::Format("Good runs (%d)", zGood->GetN()), "P");
  if (zDropped->GetN() > 0)
    zleg->AddEntry(zDropped, TString::Format("Not in GoodRun list (%d)", zDropped->GetN()), "P");
  zleg->AddEntry(zfit, TString::Format("slope = %.3g b", zfit->GetParameter(0) * 1e-6), "L");
  zleg->Draw();
  Annotate();

  // Implied ZDC cross-section per run vs. its mean coincidence rate.
  TVirtualPad* z2 = cz->cd(2);
  StylePad(z2, 0.13);
  double smin = 1e30, smax = 0, rmaxZ = 0;
  for (TGraph* g : {sGood, sDropped})
    for (int i = 0; i < g->GetN(); ++i) {
      smin = std::min(smin, g->GetY()[i]); smax = std::max(smax, g->GetY()[i]);
      rmaxZ = std::max(rmaxZ, g->GetX()[i]);
    }
  TH1* sframe = z2->DrawFrame(0, 0, 1.1 * rmaxZ, 1.45 * smax); // y axis from zero
  sframe->GetXaxis()->SetTitle("#LT R_{ZDC} #GT [Hz]");
  sframe->GetYaxis()->SetTitle("#LT R_{ZDC} #GT t_{run} livetime / L_{int} [b]");
  StyleAxes(sframe->GetXaxis(), sframe->GetYaxis());
  TF1* sline = new TF1("fZdcSigma", "[0]", 0, 1.1 * rmaxZ);
  sline->SetParameter(0, zfit->GetParameter(0) * 1e-6);
  sline->SetLineColor(kGray + 2); sline->SetLineStyle(2); sline->SetLineWidth(2);
  sline->Draw("SAME");
  sGood->Draw("P SAME");
  if (sDropped->GetN() > 0) sDropped->Draw("P SAME");
  TLegend* sleg = new TLegend(0.50, 0.74, 0.93, 0.925);
  sleg->SetBorderSize(0); sleg->SetFillStyle(0); sleg->SetTextFont(42); sleg->SetTextSize(0.033);
  sleg->AddEntry(sGood, "Good runs", "P");
  if (sDropped->GetN() > 0) sleg->AddEntry(sDropped, "Not in GoodRun list", "P");
  sleg->AddEntry(sline, "Fit slope (left panel)", "L");
  sleg->Draw();
  Annotate();
  if (nSkipped > 0) {
    TLatex note; note.SetNDC(); note.SetTextFont(42); note.SetTextSize(0.028);
    note.DrawLatex(0.18, 0.17, TString::Format("%d run(s) with L_{int} < %g %s not shown",
                                              nSkipped, minLumi, lumi_unit.Data()));
  }
  cz->Print("zdc_rate_vs_lumi.png");
  cz->Print("zdc_rate_vs_lumi.pdf");
}
