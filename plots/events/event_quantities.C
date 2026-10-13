#include "PlotStyle.h"
#include "TFile.h"
#include <memory>

void event_quantities(TString input_file = "../../data/phi-kk-mudst.root") {
  std::unique_ptr<TFile> file(TFile::Open(input_file, "READ"));

  TH1* zdc = dynamic_cast<TH1*>(file->Get("hAllEventZDCCoincidenceRate"));
  TH1* bbc = dynamic_cast<TH1*>(file->Get("hAllEventBBCCoincidenceRate"));

  PlotStyle::Options opt;
  opt.xlabel = "Coincidence rate [Hz]"; // Assumes stored x values are in Hz.
  opt.ylabel = "Counts";
  opt.showEntries = true;             // Show each histogram's total entry count.
  opt.XRange(0, 6000);                 // 0-6 kHz in the stored Hz coordinates
  // opt.YRange(0, 50000);            // Optional; otherwise auto-scale both curves.
  // opt.logY = true;                 // Automatic positive y range on a log axis.
  opt.title = "Run-19 Au+Au";
  opt.subtitle = "#sqrt{s_{NN}} = 200 GeV";
  opt.note = "ZDCMB";

  TCanvas* canvas = PlotStyle::Draw({
    {zdc, "ZDC"},
    {bbc, "BBC"}
  }, opt);
  if (!canvas) return;
  canvas->Print("all_bbc_zdc_coin.png");
  canvas->Print("all_bbc_zdc_coin.pdf"); // Vector output for analysis notes.


  TH1* all_zdc_ua_est = dynamic_cast<TH1*>(file->Get("hAllEventZDCUAEast"));
  TH1* all_zdc_ua_west = dynamic_cast<TH1*>(file->Get("hAllEventZDCUAWest"));
  PlotStyle::Options opt2_zdcua;
  opt2_zdcua.xlabel = "ZDC UA [ADC]";
  opt2_zdcua.logY = true;
    opt2_zdcua.ylabel = "Counts";
    opt2_zdcua.showEntries = true;
    opt2_zdcua.title = "Run-19 Au+Au";
    opt2_zdcua.subtitle = "#sqrt{s_{NN}} = 200 GeV";
    opt2_zdcua.note = "ZDCMB";
    TCanvas* canvas_zdcua = PlotStyle::Draw({
      {all_zdc_ua_est, "ZDC UA East"},
      {all_zdc_ua_west, "ZDC UA West"}
    }, opt2_zdcua);
    if (!canvas_zdcua) return;
    canvas_zdcua->Print("all_zdc_ua.png");
    canvas_zdcua->Print("all_zdc_ua.pdf"); // Vector output for



  

  // The input file closes here; canvas-owned histogram copies stay visible.
}
