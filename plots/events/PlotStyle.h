#ifndef PLOT_STYLE_H
#define PLOT_STYLE_H

// Reusable ROOT 6 / C++11 style for one-dimensional histograms.
// Include this header beside your ROOT macro. Original histograms are unchanged.
#include "TAxis.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TError.h"
#include "TH1.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TString.h"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <memory>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace PlotStyle {

struct Options {
  TString xlabel = "";
  TString ylabel = "Counts";
  TString title = "";       // Upper-left annotation, e.g. Run-19 Au+Au
  TString subtitle = "";    // e.g. #sqrt{s_{NN}} = 200 GeV
  TString note = "";        // e.g. ZDCMB trigger

  bool logX = false;
  bool logY = false;
  bool showLegend = true;
  bool showEntries = true; // Total TH1::GetEntries(), independent of display range.
  bool shortEntries = true; // 850, 1.2K, 3.45M; false displays the full count.
  int width = 1100;
  int height = 800;

  // Shared layout: edit these defaults to change every later plot.
  double left = 0.14, right = 0.045, bottom = 0.13, top = 0.055;
  double labelSize = 0.040, titleSize = 0.047;
  double xTitleOffset = 1.15, yTitleOffset = 1.40;
  double headroom = 1.45;
  double legendX1 = 0.50, legendX2 = 0.93, legendY2 = 0.925;
  double legendTextSize = 0.033;
  double annotationX = 0.18, annotationY = 0.915;

  // Automatic ranges are used until XRange / YRange is called.
  bool hasXRange = false, hasYRange = false;
  double xmin = 0, xmax = 1, ymin = 0, ymax = 1;

  Options& XRange(double low, double high) {
    xmin = low; xmax = high; hasXRange = true; return *this;
  }
  Options& YRange(double low, double high) {
    ymin = low; ymax = high; hasYRange = true; return *this;
  }
};

struct Series {
  const TH1* hist;
  TString label;
  TString draw;
  int color;
  int lineStyle;

  // draw: "HIST" for steps, "E1" for points with errors.
  // color < 0 and lineStyle == 0 select the shared palette automatically.
  Series(const TH1* h, TString name = "", TString option = "HIST",
         int c = -1, int style = 0)
    : hist(h), label(name), draw(option), color(c), lineStyle(style) {}
};

// Below 1000, keep the number; otherwise use K/M with up to two decimals.
inline std::string ShortNumber(double value) {
  std::ostringstream formatted;
  const double magnitude = std::abs(value);
  if (!std::isfinite(value) || magnitude < 1000) {
    formatted << std::setprecision(15) << value;
    return formatted.str();
  }
  const bool millions = magnitude >= 1e6;
  const double scale = millions ? 1e6 : 1e3;
  std::string suffix = millions ? "M" : "K";
  formatted << std::fixed << std::setprecision(2) << value / scale;
  std::string number = formatted.str();
  // Promote a rounded 1000K to 1M.
  if (suffix == "K" && (number == "1000.00" || number == "-1000.00")) {
    number = value < 0 ? "-1.00" : "1.00";
    suffix = "M";
  }
  while (!number.empty() && number.back() == '0') number.pop_back();
  if (!number.empty() && number.back() == '.') number.pop_back();
  return number + suffix;
}

inline int Color(std::size_t index) {
  const char* colors[] = {
    "#0072B2", "#D55E00", "#009E73", "#CC79A7", "#56B4E9", "#303030"
  };
  return TColor::GetColor(colors[index % 6]);
}

inline void StyleAxis(TAxis* axis, const Options& opt) {
  axis->SetLabelFont(42);
  axis->SetTitleFont(42);
  axis->SetLabelSize(opt.labelSize);
  axis->SetTitleSize(opt.titleSize);
  axis->SetLabelOffset(0.008);
  axis->SetTickLength(0.025);
  axis->SetNdivisions(505);
  axis->SetAxisColor(kBlack);
  axis->SetLabelColor(kBlack);
  axis->SetTitleColor(kBlack);
}

inline TCanvas* Draw(const std::vector<Series>& series,
                     const Options& opt = Options()) {
  if (series.empty()) {
    Error("PlotStyle::Draw", "Supply at least one histogram.");
    return nullptr;
  }
  if ((opt.hasXRange && !(opt.xmax > opt.xmin)) ||
      (opt.hasYRange && !(opt.ymax > opt.ymin)) ||
      (opt.logX && opt.hasXRange && opt.xmin <= 0) ||
      (opt.logY && opt.hasYRange && opt.ymin <= 0) || opt.headroom <= 1) {
    Error("PlotStyle::Draw", "Invalid ranges or headroom; log ranges must be positive.");
    return nullptr;
  }

  static unsigned long sequence = 0;
  const unsigned long id = ++sequence;
  std::vector<std::unique_ptr<TH1> > copies;
  std::vector<TString> labels;
  double dataMin = 0, dataMax = 0, smallestPositive = 0;
  std::size_t legendEntries = 0;
  const int styles[] = {1, 2, 3, 7, 9, 5};

  for (std::size_t i = 0; i < series.size(); ++i) {
    if (!series[i].hist || series[i].hist->GetDimension() != 1) {
      Error("PlotStyle::Draw", "Series %lu is not a valid 1D histogram.",
            static_cast<unsigned long>(i));
      return nullptr;
    }
    const TString name = TString::Format("plot_%lu_hist_%lu", id,
                                       static_cast<unsigned long>(i));
    std::unique_ptr<TH1> h(static_cast<TH1*>(series[i].hist->Clone(name)));
    h->SetDirectory(nullptr);
    h->SetTitle("");
    h->SetStats(false);
    h->SetFillStyle(0);
    const int color = series[i].color < 0 ? Color(i) : series[i].color;
    h->SetLineColor(color);
    h->SetMarkerColor(color);
    h->SetLineWidth(3);
    h->SetLineStyle(series[i].lineStyle == 0 ? styles[i % 6] : series[i].lineStyle);
    h->SetMarkerStyle(20 + static_cast<int>(i % 4));
    h->SetMarkerSize(0.9);

    // Reset inherited display ranges before applying this plot's configuration.
    h->GetXaxis()->SetRange(0, 0);
    if (opt.hasXRange) h->GetXaxis()->SetRangeUser(opt.xmin, opt.xmax);
    StyleAxis(h->GetXaxis(), opt);
    StyleAxis(h->GetYaxis(), opt);
    h->GetXaxis()->SetTitle(opt.xlabel);
    h->GetYaxis()->SetTitle(opt.ylabel);
    h->GetXaxis()->SetTitleOffset(opt.xTitleOffset);
    h->GetYaxis()->SetTitleOffset(opt.yTitleOffset);

    TString drawOption = series[i].draw;
    drawOption.ToUpper();
    const bool withErrors = drawOption.Contains("E");
    // Find the shared y range using only the bins visible in this plot.
    const int first = std::max(1, h->GetXaxis()->GetFirst());
    const int last = std::min(h->GetNbinsX(), h->GetXaxis()->GetLast());
    for (int bin = first; bin <= last; ++bin) {
      const double value = h->GetBinContent(bin);
      const double error = withErrors ? h->GetBinError(bin) : 0;
      if (!std::isfinite(value) || !std::isfinite(error)) continue;
      dataMin = std::min(dataMin, value - error);
      dataMax = std::max(dataMax, value + error);
      const double positive = value - error > 0 ? value - error : value;
      if (positive > 0 && (smallestPositive == 0 || positive < smallestPositive))
        smallestPositive = positive;
    }
    TString label = series[i].label;
    if (opt.showEntries) {
      const double entries = series[i].hist->GetEntries();
      const TString number = opt.shortEntries ? TString(ShortNumber(entries).c_str())
                                              : TString::Format("%.0f", entries);
      const TString count = "Entries = " + number;
      label = label.IsNull() ? count : label + "  (" + count + ")";
    }
    if (!label.IsNull()) ++legendEntries;
    labels.push_back(label);
    copies.push_back(std::move(h));
  }

  // Without XRange, use the union of the histogram domains.
  double xlow = copies[0]->GetXaxis()->GetXmin();
  double xhigh = copies[0]->GetXaxis()->GetXmax();
  for (std::size_t i = 1; i < copies.size(); ++i) {
    xlow = std::min(xlow, copies[i]->GetXaxis()->GetXmin());
    xhigh = std::max(xhigh, copies[i]->GetXaxis()->GetXmax());
  }
  if (opt.hasXRange) { xlow = opt.xmin; xhigh = opt.xmax; }
  if (opt.logX && xlow <= 0) {
    Error("PlotStyle::Draw", "Set a positive XRange for a logarithmic x axis.");
    return nullptr;
  }

  double ylow = dataMin < 0 ? dataMin - 0.08 * (dataMax - dataMin) : 0;
  double yhigh = dataMax > ylow ? ylow + opt.headroom * (dataMax - ylow) : ylow + 1;
  if (opt.logY) {
    ylow = smallestPositive > 0 ? smallestPositive * 0.5 : 0.1;
    const double peak = std::max(dataMax, ylow * 10);
    yhigh = ylow * std::pow(peak / ylow, opt.headroom);
  }
  if (opt.hasYRange) { ylow = opt.ymin; yhigh = opt.ymax; }

  TCanvas* canvas = new TCanvas(TString::Format("plot_canvas_%lu", id),
                                opt.title, opt.width, opt.height);
  canvas->SetFillColor(kWhite);
  canvas->SetBorderMode(0);
  canvas->SetFrameFillColor(kWhite);
  canvas->SetFrameBorderMode(0);
  canvas->SetFrameLineWidth(1);
  canvas->SetLeftMargin(opt.left);
  canvas->SetRightMargin(opt.right);
  canvas->SetBottomMargin(opt.bottom);
  canvas->SetTopMargin(opt.top);
  canvas->SetTickx(1);
  canvas->SetTicky(1);
  canvas->SetLogx(opt.logX);
  canvas->SetLogy(opt.logY);
  canvas->cd();

  // A separate frame gives overlays the same axis range, even with different bins.
  TH1F* frame = canvas->DrawFrame(xlow, ylow, xhigh, yhigh);
  frame->SetTitle("");
  frame->SetStats(false);
  StyleAxis(frame->GetXaxis(), opt);
  StyleAxis(frame->GetYaxis(), opt);
  frame->GetXaxis()->SetTitle(opt.xlabel);
  frame->GetYaxis()->SetTitle(opt.ylabel);
  frame->GetXaxis()->SetTitleOffset(opt.xTitleOffset);
  frame->GetYaxis()->SetTitleOffset(opt.yTitleOffset);

  TLegend* legend = nullptr;
  if (opt.showLegend && legendEntries > 0) {
    const double ybottom = opt.legendY2 - 0.051 * legendEntries;
    legend = new TLegend(opt.legendX1, ybottom, opt.legendX2, opt.legendY2);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetTextFont(42);
    legend->SetTextSize(opt.legendTextSize);
    legend->SetMargin(0.22);
    legend->SetBit(TObject::kCanDelete);
  }
  for (std::size_t i = 0; i < copies.size(); ++i) {
    TString option = series[i].draw;
    option += " SAME";
    // DrawCopy makes the displayed histogram belong to the canvas.
    TH1* drawn = copies[i]->DrawCopy(option);
    if (legend && !labels[i].IsNull()) {
      TString upper = series[i].draw;
      upper.ToUpper();
      legend->AddEntry(drawn, labels[i],
                       upper.Contains("E") || upper.Contains("P") ? "lep" : "l");
    }
  }

  TLatex text;
  text.SetNDC();
  text.SetTextAlign(13);
  text.SetTextColor(kBlack);
  double textY = opt.annotationY;
  if (!opt.title.IsNull()) {
    text.SetTextFont(62);
    text.SetTextSize(0.041);
    text.DrawLatex(opt.annotationX, textY, opt.title);
    textY -= 0.056;
  }
  text.SetTextFont(42);
  text.SetTextSize(0.036);
  if (!opt.subtitle.IsNull()) {
    text.DrawLatex(opt.annotationX, textY, opt.subtitle);
    textY -= 0.050;
  }
  if (!opt.note.IsNull()) text.DrawLatex(opt.annotationX, textY, opt.note);
  // If the legend is hidden, keep the requested entry counters as plain text.
  if (!opt.showLegend && opt.showEntries) {
    text.SetTextAlign(33);
    text.SetTextSize(opt.legendTextSize);
    for (std::size_t i = 0; i < labels.size(); ++i) {
      text.SetTextColor(copies[i]->GetLineColor());
      text.DrawLatex(opt.legendX2, opt.legendY2 - 0.051 * i, labels[i]);
    }
  }
  if (legend) legend->Draw();
  canvas->RedrawAxis();
  canvas->Modified();
  canvas->Update();
  return canvas;
}

inline TCanvas* Draw(const TH1* histogram, const Options& opt = Options(),
                     TString drawOption = "HIST") {
  return Draw(std::vector<Series>(1, Series(histogram, "", drawOption)), opt);
}

} // namespace PlotStyle
#endif
