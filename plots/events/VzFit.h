#ifndef VZ_FIT_H
#define VZ_FIT_H

// Shared helpers for the vz efficiency macros (vz.C, vz_pion.C).
#include "PlotStyle.h"
#include "TF1.h"
#include "TH1D.h"
#include "TH2.h"
#include "TLine.h"
#include "TMath.h"
#include <cstdio>
#include <vector>

// 2D counterpart of PlotStyle::Draw, sharing its axis style and annotation.
TCanvas* Draw2D(TH2* h, const PlotStyle::Options& opt) {
  TCanvas* canvas = new TCanvas(TString::Format("c_%s", h->GetName()), opt.title,
                                opt.width, opt.height);
  canvas->SetLeftMargin(opt.left);
  canvas->SetRightMargin(0.15); // room for the palette
  canvas->SetBottomMargin(opt.bottom);
  canvas->SetTopMargin(opt.top);
  canvas->SetTickx(1);
  canvas->SetTicky(1);
  canvas->SetLogz(opt.logY); // logY doubles as log z for 2D
  h->SetTitle("");
  h->SetStats(false);
  PlotStyle::StyleAxis(h->GetXaxis(), opt);
  PlotStyle::StyleAxis(h->GetYaxis(), opt);
  PlotStyle::StyleAxis(h->GetZaxis(), opt);
  h->GetXaxis()->SetTitle(opt.xlabel);
  h->GetYaxis()->SetTitle(opt.ylabel);
  h->GetXaxis()->SetTitleOffset(opt.xTitleOffset);
  h->GetYaxis()->SetTitleOffset(opt.yTitleOffset);
  h->DrawCopy("COLZ");

  TLatex text;
  text.SetNDC();
  text.SetTextAlign(13);
  double textY = opt.annotationY;
  text.SetTextFont(62);
  text.SetTextSize(0.041);
  text.DrawLatex(opt.annotationX, textY, opt.title);
  textY -= 0.056;
  text.SetTextFont(42);
  text.SetTextSize(0.036);
  text.DrawLatex(opt.annotationX, textY, opt.subtitle);
  textY -= 0.050;
  text.DrawLatex(opt.annotationX, textY, opt.note);
  if (opt.showEntries) {
    textY -= 0.050;
    text.DrawLatex(opt.annotationX, textY,
                   ("Entries = " + PlotStyle::ShortNumber(h->GetEntries())).c_str());
  }
  canvas->RedrawAxis();
  return canvas;
}

// Extra annotation lines (e.g. selection criteria), drawn in a column starting at (x, y) in NDC.
// align follows TLatex: 13 = left/top (default), 33 = right/top.
void DrawLines(const std::vector<TString>& lines, double x, double y, double size = 0.028,
               int align = 13) {
  TLatex text;
  text.SetNDC();
  text.SetTextFont(42);
  text.SetTextSize(size);
  text.SetTextAlign(align);
  for (const TString& line : lines) {
    text.DrawLatex(x, y, line);
    y -= 1.3 * size;
  }
}

// Right edge of the frame in Draw2D (left of the palette), for right-aligned text.
const double kDraw2DRight = 0.83;

// NDC height just below the title / subtitle / note block of PlotStyle and Draw2D.
inline double BelowNote(const PlotStyle::Options& opt) { return opt.annotationY - 0.056 - 0.050 - 0.050; }

// Fraction of a Gaussian inside |vz| < zcut.
double GausFraction(double mean, double sigma, double zcut = 100) {
  return 0.5 * (TMath::Erf((zcut - mean) / (std::sqrt(2.) * sigma)) -
                TMath::Erf((-zcut - mean) / (std::sqrt(2.) * sigma)));
}

// Double-Gaussian vz fit; draws data, total, narrow and wide components, and
// reports |vz| < 100 cm fractions for the narrow Gaussian, both Gaussians, and a direct count.
// extra: optional lines drawn below the annotation (e.g. selection criteria).
void FitVz(TH1D* h, PlotStyle::Options opt, TString label, TString name,
           const std::vector<TString>& extra = {}) {
  const double range = 190; // stay inside the +-200 cm vertex-finder edge
  TF1 fit("fit_" + name, "gaus(0)+gaus(3)", -range, range);
  fit.SetParameters(h->GetMaximum(), -5, 20, 0.05 * h->GetMaximum(), -5, 80);
  fit.SetParLimits(2, 5, 40);    // narrow
  fit.SetParLimits(5, 25, 2000); // wide
  h->Fit(&fit, "QRN0");

  const double m1 = fit.GetParameter(1), s1 = std::abs(fit.GetParameter(2));
  const double m2 = fit.GetParameter(4), s2 = std::abs(fit.GetParameter(5));
  const double a1 = fit.GetParameter(0) * s1, a2 = fit.GetParameter(3) * s2;
  const double w2 = a2 / (a1 + a2);
  const double f_narrow = GausFraction(m1, s1);
  const double f_both = (1 - w2) * f_narrow + w2 * GausFraction(m2, s2);
  const double f_count = h->Integral(h->FindBin(-99.999), h->FindBin(99.999)) / h->Integral();
  printf("%-28s mu1=%6.2f s1=%5.2f | mu2=%6.2f s2=%6.1f w2=%.3f | chi2/ndf=%6.1f/%d |"
         " f_narrow=%.4f f_both=%.4f f_count=%.4f\n",
         name.Data(), m1, s1, m2, s2, w2, fit.GetChisquare(), fit.GetNDF(),
         f_narrow, f_both, f_count);

  opt.xlabel = "v_{z} [cm]";
  opt.ylabel = "Events / 5 cm";
  opt.logY = false;
  // Compact legend in the middle so the fit summary sits beside it on the same rows.
  opt.legendX1 = 0.49;
  opt.legendX2 = 0.69;
  opt.legendTextSize = 0.027;
  TCanvas* canvas = PlotStyle::Draw({{h, label, "E1", kBlack}}, opt);
  if (!canvas) return;
  canvas->cd();

  TF1* total = new TF1("total_" + name, "gaus(0)+gaus(3)", -range, range);
  total->SetParameters(fit.GetParameters());
  TF1* narrow = new TF1("narrow_" + name, "gaus", -range, range);
  narrow->SetParameters(fit.GetParameter(0), m1, s1);
  TF1* wide = new TF1("wide_" + name, "gaus", -range, range);
  wide->SetParameters(fit.GetParameter(3), m2, s2);
  const int colors[3] = {PlotStyle::Color(1), PlotStyle::Color(0), PlotStyle::Color(2)};
  TF1* curves[3] = {total, narrow, wide};
  for (int i = 0; i < 3; ++i) {
    curves[i]->SetNpx(1000);
    curves[i]->SetLineColor(colors[i]);
    curves[i]->SetLineWidth(i == 0 ? 3 : 2);
    curves[i]->SetLineStyle(i == 0 ? 1 : 2);
    curves[i]->Draw("SAME");
  }

  // |vz| = 100 cm selection boundaries, stopping below the annotation and fit summary.
  canvas->Update();
  const double line_top = canvas->GetUymin() + 0.70 * (canvas->GetUymax() - canvas->GetUymin());
  for (double z : {-100., 100.}) {
    TLine* edge = new TLine(z, canvas->GetUymin(), z, line_top);
    edge->SetLineColor(kGray + 2);
    edge->SetLineStyle(2);
    edge->SetLineWidth(2);
    edge->Draw();
  }

  TLegend* legend = nullptr;
  for (TObject* obj : *canvas->GetListOfPrimitives())
    if ((legend = dynamic_cast<TLegend*>(obj))) break;
  if (legend) {
    legend->SetY1NDC(opt.legendY2 - 0.04 * 4);
    legend->AddEntry(total, "Double Gaussian", "l");
    legend->AddEntry(narrow, "Narrow", "l");
    legend->AddEntry(wide, "Wide", "l");
  }

  TLatex text;
  text.SetNDC();
  text.SetTextFont(42);
  text.SetTextSize(0.027);
  text.SetTextAlign(13);
  const double x = opt.legendX2 + 0.01;
  double y = opt.legendY2 - 0.005;
  const TString lines[5] = {
    TString::Format("#sigma_{1} = %.1f cm, #sigma_{2} = %.1f cm", s1, s2),
    TString::Format("w_{2} = %.3f, #chi^{2}/ndf = %.0f/%d", w2, fit.GetChisquare(), fit.GetNDF()),
    TString::Format("f_{narrow}(|v_{z}|<100) = %.4f", f_narrow),
    TString::Format("f_{both}(|v_{z}|<100) = %.4f", f_both),
    TString::Format("f_{count}(|v_{z}|<100) = %.4f", f_count)};
  for (const TString& line : lines) {
    text.DrawLatex(x, y, line);
    y -= 0.037;
  }
  DrawLines(extra, opt.annotationX, BelowNote(opt));
  canvas->Modified();
  canvas->Update();
  canvas->Print(name + ".png");
  canvas->Print(name + ".pdf");
}

// ULS and ULS - LS double-Gaussian fits on a shared y range; the top ~28% of the
// frame is kept free for the annotation, legend and fit summary.
void FitVzPair(TH1D* uls, TH1D* ls, PlotStyle::Options opt, TString name,
               const std::vector<TString>& extra = {}) {
  uls->Sumw2();
  ls->Sumw2();
  TH1D* sub = static_cast<TH1D*>(uls->Clone(name + "_sub_h"));
  sub->Add(ls, -1);

  double ylow = 0, yhigh = 0;
  for (TH1D* h : {uls, sub})
    for (int bin = 1; bin <= h->GetNbinsX(); ++bin) {
      ylow = std::min(ylow, h->GetBinContent(bin) - h->GetBinError(bin));
      yhigh = std::max(yhigh, h->GetBinContent(bin) + h->GetBinError(bin));
    }
  if (ylow < 0) ylow *= 1.2;

  opt.YRange(ylow, ylow + (yhigh - ylow) / 0.72);
  opt.showEntries = false; // ULS - LS entries would count both inputs
  FitVz(uls, opt, "ULS", name + "_uls", extra);
  FitVz(sub, opt, "ULS #minus LS", name + "_sub", extra);
}

#endif
