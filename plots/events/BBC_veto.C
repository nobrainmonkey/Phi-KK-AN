#include "DiPion.h"
#include "TFile.h"
#include "TBox.h"
#include "TStyle.h"
#include "TTree.h"
#include <memory>
#include <vector>

// BBC veto efficiency from coherent rho0 -> pi+ pi-: the modified Soding fit to the ULS - LS
// mass spectrum without the BBC veto, and separately for events passing and rejected by it.
// The efficiency is Y_pass / (Y_pass + Y_rejected); the two samples are disjoint, so their fit
// uncertainties are independent. Y_pass / Y_all is kept as a cross-check.
// Input and tree preselection as in vz_pion.C (UPCDst rho -> pi pi trees).
//
// Di-pion selection without the BBC veto ("P" is replaced by uls / ls). No TofMult cut:
// TOF activity is correlated with the BBC, so a TofMult cut would bias the veto efficiency.
const TString kNoVeto =
  "P_pair_n_tof_matched>=1"                                 // >= 1 TOF-matched track
  "&&abs(P_kaon1_nsigma_pion)<3&&abs(P_kaon2_nsigma_pion)<3" // pion dE/dx
  "&&P_phi_pt<0.1";                                        // coherent region
const TString kVeto = "&&P_pair_BBC_Veto_East==0&&P_pair_BBC_Veto_West==0";
const TString kRejected = "&&!(P_pair_BBC_Veto_East==0&&P_pair_BBC_Veto_West==0)";

const std::vector<TString> kSelectionLines = {
  "2 primary tracks",
  "N_{TOF match} #geq 1, |n#sigma_{#pi}| < 3",
  "p_{T}^{#pi#pi} < 0.1 GeV/c"};

struct SodingYield { double total, total_err, rho, rho_err; };

// Yields from the fit: everything but the remnant background over the fit range, and the rho
// Breit-Wigner from 2 m_pi to M_rho + 5 Gamma_rho (as in the paper).
SodingYield Yields(TF1* fit, const TFitResultPtr& result, TString name) {
  const double lo = fit->GetXmin(), hi = fit->GetXmax();
  const TMatrixDSym cov = result->GetCovarianceMatrix();
  TF1 signal(name + "_signal", SodingSignal, lo, hi, 10);
  TF1 rho(name + "_rho", SodingRho, 2 * kMpi, 2.0, 10);
  signal.SetParameters(fit->GetParameters());
  rho.SetParameters(fit->GetParameters());
  const double rho_hi = fit->GetParameter(0) + 5 * fit->GetParameter(1);
  return {signal.Integral(lo, hi), signal.IntegralError(lo, hi, fit->GetParameters(), cov.GetMatrixArray()),
          rho.Integral(2 * kMpi + 1e-6, rho_hi),
          rho.IntegralError(2 * kMpi + 1e-6, rho_hi, fit->GetParameters(), cov.GetMatrixArray())};
}

void BBC_veto(TString input_file = "../../data/rho_pipi_upcdst.root") {
  std::unique_ptr<TFile> file(TFile::Open(input_file, "READ"));
  TTree* uls_tree = dynamic_cast<TTree*>(file->Get("uls_data"));
  TTree* ls_tree = dynamic_cast<TTree*>(file->Get("ls_data"));
  gStyle->SetPalette(kBird);

  PlotStyle::Options opt;
  opt.title = "Run-19 Au+Au";
  opt.subtitle = "#sqrt{s_{NN}} = 200 GeV";

  const double plot_lo = 0.5, lo = 0.6, hi = 1.2; // shown from plot_lo, fitted from lo
  const TString cut[3] = {kNoVeto, kNoVeto + kVeto, kNoVeto + kRejected};
  const TString note[3] = {"ZDCMB, no BBC veto", "ZDCMB, BBC veto", "ZDCMB, BBC rejected"};
  const TString name[3] = {"BBC_veto_pre", "BBC_veto_post", "BBC_veto_rejected"};
  TH1D* spectrum[3];
  SodingYield yield[3];
  printf("%-18s %-9s %-15s %-15s %-15s %-15s %-15s %-12s\n", "sample", "chi2/ndf", "M_rho [GeV]",
         "G_rho [GeV]", "B/A", "C/A", "G_w [MeV]", "phi_w");
  for (int i = 0; i < 3; ++i) {
    spectrum[i] = MassSpectrum(uls_tree, ls_tree, cut[i], name[i] + "_h", plot_lo, hi);
    TF1* fit = new TF1(name[i] + "_fit", SodingFull, lo, hi, 10);
    const TFitResultPtr result = FitSoding(spectrum[i], fit);
    yield[i] = Yields(fit, result, name[i]);
    const TMatrixDSym cov = result->GetCovarianceMatrix();
    const auto ba = AmplitudeRatio(fit, cov, 3), ca = AmplitudeRatio(fit, cov, 4);
    printf("%-18s %4.1f/%-4d %.4f+-%.4f  %.4f+-%.4f  %+.3f+-%.3f  %.3f+-%.3f  %4.1f+-%4.1f  %.2f+-%.2f\n",
           name[i].Data(), fit->GetChisquare(), fit->GetNDF(), fit->GetParameter(0), fit->GetParError(0),
           fit->GetParameter(1), fit->GetParError(1), ba.first, ba.second, ca.first, ca.second,
           1000 * fit->GetParameter(6), 1000 * fit->GetParError(6), fit->GetParameter(7), fit->GetParError(7));
    printf("%-18s yield(all - f_p) = %.0f +- %.0f, rho BW = %.0f +- %.0f\n", "",
           yield[i].total, yield[i].total_err, yield[i].rho, yield[i].rho_err);

    PlotStyle::Options opt_fit = opt;
    opt_fit.note = note[i];
    DrawSodingFit(spectrum[i], fit, opt_fit, {
      TString::Format("M_{#rho} = %.4f, #Gamma_{#rho} = %.3f GeV/c^{2}", fit->GetParameter(0), fit->GetParameter(1)),
      TString::Format("#Gamma_{#omega} = %.1f MeV/c^{2}, #chi^{2}/ndf = %.0f/%d", 1000 * fit->GetParameter(6),
                      fit->GetChisquare(), fit->GetNDF()),
      TString::Format("Yield (all #minus f_{p}) = %.0f #pm %.0f", yield[i].total, yield[i].total_err),
      TString::Format("#rho^{0} BW yield = %.0f #pm %.0f", yield[i].rho, yield[i].rho_err)},
      kSelectionLines, name[i]);
  }

  // Efficiency from the disjoint pass / rejected samples, with independent uncertainties.
  auto Split = [](double pass, double pass_err, double rej, double rej_err) {
    const double sum = pass + rej;
    return std::make_pair(pass / sum, std::hypot(rej * pass_err, pass * rej_err) / (sum * sum));
  };
  const auto eff_total = Split(yield[1].total, yield[1].total_err, yield[2].total, yield[2].total_err);
  const auto eff_rho = Split(yield[1].rho, yield[1].rho_err, yield[2].rho, yield[2].rho_err);
  // Direct ULS - LS count over the fit range, with its statistical uncertainty.
  auto Count = [&](TH1D* h, double& err) {
    return h->IntegralAndError(h->FindBin(lo + 1e-6), h->FindBin(hi - 1e-6), err, "width");
  };
  double count_err[3], count[3];
  for (int i = 0; i < 3; ++i) count[i] = Count(spectrum[i], count_err[i]);
  const auto eff_count = Split(count[1], count_err[1], count[2], count_err[2]);
  printf("BBC veto efficiency, pass / (pass + rejected):\n");
  printf("  yield(all - f_p) %.4f +- %.4f | rho BW %.4f +- %.4f | ULS-LS count %.4f +- %.4f\n",
         eff_total.first, eff_total.second, eff_rho.first, eff_rho.second, eff_count.first, eff_count.second);
  printf("cross-check, pass / all: yield(all - f_p) %.4f | rho BW %.4f | count %.4f; "
         "closure (pass + rejected) / all: %.4f\n",
         yield[1].total / yield[0].total, yield[1].rho / yield[0].rho, count[1] / count[0],
         (yield[1].total + yield[2].total) / yield[0].total);

  // Spectra without and with the veto on one canvas.
  PlotStyle::Options opt_cmp = opt;
  opt_cmp.note = "ZDCMB, ULS #minus LS";
  opt_cmp.xlabel = "M_{#pi#pi} [GeV/c^{2}]";
  opt_cmp.ylabel = "dN/dM_{#pi#pi} [(GeV/c^{2})^{-1}]";
  opt_cmp.showEntries = false;
  opt_cmp.headroom = 1.8; // room for the selection lines above the rho peak
  TCanvas* c_cmp = PlotStyle::Draw({
    {spectrum[0], "No BBC veto", "E1", -1, 1},
    {spectrum[1], "BBC veto", "E1", -1, 1}
  }, opt_cmp);
  if (!c_cmp) return;
  DrawLines(kSelectionLines, opt.annotationX, BelowNote(opt));
  DrawLines({TString::Format("#varepsilon_{BBC} (all #minus f_{p}) = %.3f #pm %.3f", eff_total.first, eff_total.second),
             TString::Format("#varepsilon_{BBC} (#rho^{0} BW) = %.3f #pm %.3f", eff_rho.first, eff_rho.second),
             TString::Format("#varepsilon_{BBC} (ULS #minus LS) = %.3f #pm %.3f", eff_count.first, eff_count.second)},
            opt.legendX1 + 0.05, opt.legendY2 - 0.051 * 2 - 0.02);
  c_cmp->Print("BBC_veto_compare.png");
  c_cmp->Print("BBC_veto_compare.pdf");

  // Efficiency vs M_pipi and rapidity in the rho window, from ULS - LS counts with binomial
  // uncertainties. The counts include remnant background, which the veto removes preferentially,
  // so they sit slightly below the fit-based nominal value; flatness is judged against the
  // integrated count (dashed line).
  const TString window = "&&P_phi_mass>0.6&&P_phi_mass<0.95";
  auto Pairs = [&](TString c) {
    TString u = c, l = c;
    u.ReplaceAll("P_", "uls_");
    l.ReplaceAll("P_", "ls_");
    return double(uls_tree->GetEntries(u)) - double(ls_tree->GetEntries(l));
  };
  const double window_eff = Pairs(kNoVeto + window + kVeto) / Pairs(kNoVeto + window);
  struct Scan { TString name, var, xlabel; std::vector<double> edges; };
  const Scan scans[2] = {
    {"BBC_veto_eff_mass", "P_phi_mass", "M_{#pi#pi} [GeV/c^{2}]", {0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95}},
    {"BBC_veto_eff_rapidity", "P_phi_rapidity", "y_{#pi#pi}", {-1.0, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1.0}}};
  for (const Scan& scan : scans) {
    TH1D* eff = new TH1D(scan.name + "_h", "", scan.edges.size() - 1, scan.edges.data());
    for (int bin = 1; bin <= eff->GetNbinsX(); ++bin) {
      const TString range = TString::Format("&&%s>=%g&&%s<%g", scan.var.Data(), eff->GetBinLowEdge(bin),
                                            scan.var.Data(), eff->GetBinLowEdge(bin + 1));
      const double all = Pairs(kNoVeto + window + range), pass = Pairs(kNoVeto + window + range + kVeto);
      const double e = pass / all;
      eff->SetBinContent(bin, e);
      eff->SetBinError(bin, std::sqrt(e * (1 - e) / all));
    }
    PlotStyle::Options opt_e = opt;
    opt_e.note = "ZDCMB, ULS #minus LS";
    opt_e.xlabel = scan.xlabel;
    opt_e.ylabel = "#varepsilon_{BBC}";
    opt_e.showEntries = false;
    opt_e.YRange(0.70, 0.95);
    TCanvas* c_e = PlotStyle::Draw({{eff, "ULS #minus LS count per bin", "E1", kBlack}}, opt_e);
    if (!c_e) return;
    c_e->cd();
    const double x0 = eff->GetXaxis()->GetXmin(), x1 = eff->GetXaxis()->GetXmax();
    TBox* band = new TBox(x0, eff_total.first - eff_total.second, x1, eff_total.first + eff_total.second);
    band->SetFillColorAlpha(PlotStyle::Color(1), 0.25);
    band->SetLineWidth(0);
    band->Draw();
    TLine* nominal = new TLine(x0, eff_total.first, x1, eff_total.first);
    nominal->SetLineColor(PlotStyle::Color(1));
    nominal->SetLineWidth(3);
    nominal->Draw();
    TLine* counted = new TLine(x0, window_eff, x1, window_eff);
    counted->SetLineColor(kGray + 2);
    counted->SetLineStyle(2);
    counted->SetLineWidth(2);
    counted->Draw();
    eff->SetLineWidth(3); // redraw on top of the band, in the legend style
    eff->SetLineColor(kBlack);
    eff->SetMarkerColor(kBlack);
    eff->SetMarkerStyle(20);
    eff->SetMarkerSize(0.9);
    eff->Draw("E1 SAME");
    TLegend* legend = nullptr;
    for (TObject* obj : *c_e->GetListOfPrimitives())
      if ((legend = dynamic_cast<TLegend*>(obj))) break;
    if (legend) {
      legend->SetY1NDC(opt.legendY2 - 0.051 * 3);
      legend->AddEntry(nominal, TString::Format("Nominal %.3f #pm %.3f", eff_total.first, eff_total.second), "l");
      legend->AddEntry(counted, TString::Format("Count, 0.6 < M_{#pi#pi} < 0.95: %.3f", window_eff), "l");
    }
    std::vector<TString> lines = kSelectionLines;
    lines.push_back("0.60 < M_{#pi#pi} < 0.95 GeV/c^{2}");
    DrawLines(lines, opt.annotationX, BelowNote(opt));
    c_e->Print(scan.name + ".png");
    c_e->Print(scan.name + ".pdf");
  }
}
