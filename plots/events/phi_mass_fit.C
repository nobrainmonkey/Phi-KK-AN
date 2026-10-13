#include "TCanvas.h"
#include "TH1.h"
#include "TLegend.h"
#include "TF1.h"
#include "TFile.h"
#include "TFitResult.h"
#include "TLatex.h"
#include "TLine.h"
#include "TGraph.h"
#include "TMath.h"
#include "TPad.h"
#include "TSystem.h"
#include "TTree.h"
#include <memory>

// phi -> K+K- ULS mass spectrum fitted with a Voigtian (Breit-Wigner (x) Gaussian) + pol1.
// The BW width is fixed to the PDG phi width, so sigma is the detector mass resolution.
// Binned Poisson likelihood fit; chi2/ndf is the Baker-Cousins likelihood-ratio chi2.
// Exactly 2 primary tracks in the whole event (the production cut is per vertex only). The cut
// vr <= vr_max is an argument of phi_mass_fit; vz and vr are both from the event vertex, and vr
// is about (0,0) as in vz.C.
const TString kCut =
  "uls_phi_pt>=0.3&&uls_phi_pt<=0.9&&abs(uls_phi_rapidity)<0.5"
  "&&uls_pair_BBC_Veto_East==0&&uls_pair_BBC_Veto_West==0&&abs(uls_pair_event_vz)<=100"
  "&&uls_pair_event_n_primary_tracks==2"
  "&&uls_phi_x2kk<=8&&uls_phi_x2pipi>=15&&uls_phi_x2ee>=10&&uls_phi_x2ee>2*uls_phi_x2kk";

const double kGammaPhi = 0.004249; // PDG phi width [GeV]
const double kXref = 1.02;         // pol1 pivot [GeV]
// 45 bins over 0.98-1.10 GeV (2.67 MeV), as in the reference plot. 1.00 and 1.04 are not bin
// edges, so window counts come from a histogram with kFine times finer bins.
const double kPlotLo = 0.98, kPlotHi = 1.10;
const int kNbins = 45, kFine = 8;
const double kBinWidth = (kPlotHi - kPlotLo) / kNbins; // [GeV]

// par: 0 N_sig (total), 1 mass, 2 sigma, 3 Gamma, 4 bkg density at kXref [1/GeV], 5 slope [1/GeV^2].
// Function values are counts per bin; integrals divided by kBinWidth are yields.
double Signal(double* x, double* p) { return kBinWidth * p[0] * TMath::Voigt(x[0] - p[1], p[2], p[3]); }
double Background(double* x, double* p) { return kBinWidth * (p[4] + p[5] * (x[0] - kXref)); }
double Total(double* x, double* p) { return Signal(x, p) + Background(x, p); }

// Yield of one component in [lo, hi] and its error from the fit covariance (central-difference
// gradient), integrating with a relative tolerance the Voigtian can actually reach.
std::pair<double, double> WindowYield(TF1& f, const double* par, const TMatrixDSym& cov, double lo, double hi) {
  auto yield = [&](const std::vector<double>& q) {
    f.SetParameters(q.data());
    return f.Integral(lo, hi, 1e-9) / kBinWidth;
  };
  const int n = f.GetNpar();
  const std::vector<double> p(par, par + n);
  std::vector<double> grad(n, 0);
  for (int i = 0; i < n; ++i) {
    if (cov(i, i) <= 0) continue;
    const double step = 1e-3 * std::sqrt(cov(i, i));
    std::vector<double> up = p, down = p;
    up[i] += step;
    down[i] -= step;
    grad[i] = (yield(up) - yield(down)) / (2 * step);
  }
  double var = 0;
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j) var += grad[i] * grad[j] * cov(i, j);
  return {yield(p), std::sqrt(var)};
}

struct Gof { double bc, pearson; int nbins; TGraph* pull; };

// Goodness of fit over the bins in [lo, hi], and the pulls (n - mu) / sqrt(mu) with mu the
// bin-integrated fit, so that their squares sum to the Pearson chi2.
Gof GoodnessOfFit(const TH1* h, TF1* f, double lo, double hi) {
  Gof g{0, 0, 0, new TGraph()};
  for (int b = h->GetXaxis()->FindBin(lo + 1e-9); b <= h->GetXaxis()->FindBin(hi - 1e-9); ++b) {
    const double n = h->GetBinContent(b);
    const double mu = f->Integral(h->GetBinLowEdge(b), h->GetBinLowEdge(b + 1)) / kBinWidth;
    g.bc += 2 * (mu - n + (n > 0 ? n * std::log(n / mu) : 0));
    g.pearson += (n - mu) * (n - mu) / mu;
    g.pull->SetPoint(g.nbins++, h->GetBinCenter(b), (n - mu) / std::sqrt(mu));
  }
  return g;
}

void phi_mass_fit(TString input_file = "../../data/phi-kk-mudst.root", double fit_lo = 1.00,
                  double fit_hi = 1.10, double win_lo = 1.00, double win_hi = 1.04, double vr_max = 5) {
  const TString cut = kCut + TString::Format("&&sqrt(uls_pair_event_vx**2+uls_pair_event_vy**2)<=%g", vr_max);
  const TString name = TString("phi_mass_fit") + (vr_max != 5 ? TString::Format("_vr%g", vr_max) : TString());
  const TString vertex_label = TString::Format("N_{trk}^{prim} = 2, v_{r} #leq %g cm", vr_max);
  std::unique_ptr<TFile> file(TFile::Open(input_file, "READ"));
  TTree* tree = dynamic_cast<TTree*>(file->Get("uls_data"));
  tree->SetBranchStatus("*", 0);
  for (const char* b : {"uls_phi_mass", "uls_phi_pt", "uls_phi_rapidity", "uls_pair_BBC_Veto_East",
                        "uls_pair_BBC_Veto_West", "uls_pair_event_vx", "uls_pair_event_vy",
                        "uls_pair_event_vz", "uls_pair_event_n_primary_tracks", "uls_phi_x2kk",
                        "uls_phi_x2pipi", "uls_phi_x2ee"})
    tree->SetBranchStatus(b, 1);

  const double plot_lo = kPlotLo, plot_hi = kPlotHi;
  TH1D* h_fine = new TH1D("phi_mass_fine", "", kFine * kNbins, plot_lo, plot_hi);
  h_fine->Sumw2(false);
  tree->Draw("uls_phi_mass>>phi_mass_fine", cut, "goff");
  TH1D* h = dynamic_cast<TH1D*>(h_fine->Rebin(kFine, "phi_mass"));
  // Fit whole bins: widen [fit_lo, fit_hi] to the edges of the bins that contain them.
  fit_lo = h->GetBinLowEdge(h->FindBin(fit_lo + 1e-9));
  fit_hi = h->GetBinLowEdge(h->FindBin(fit_hi - 1e-9) + 1);

  TF1* fit = new TF1("phi_fit", Total, fit_lo, fit_hi, 6);
  fit->SetParNames("N_{sig}", "M", "#sigma", "#Gamma", "B_{0}", "B_{1}");
  fit->SetParameters(250, 1.0195, 0.003, kGammaPhi, 2000, 0);
  fit->FixParameter(3, kGammaPhi);
  fit->SetParLimits(2, 1e-4, 0.02);
  const TFitResultPtr result = h->Fit(fit, "LSRQ0");
  const double* p = fit->GetParameters();
  const TMatrixDSym cov = result->GetCovarianceMatrix();

  TF1 sig("phi_sig", Signal, plot_lo, plot_hi, 6), bkg("phi_bkg", Background, plot_lo, plot_hi, 6);
  const auto [S, S_err] = WindowYield(sig, p, cov, win_lo, win_hi);
  const auto [B, B_err] = WindowYield(bkg, p, cov, win_lo, win_hi);
  const double signif = S / std::sqrt(S + B);
  const double observed = h_fine->Integral(h_fine->FindBin(win_lo + 1e-9), h_fine->FindBin(win_hi - 1e-9));
  const Gof gof = GoodnessOfFit(h, fit, fit_lo, fit_hi);
  const int ndf = gof.nbins - fit->GetNumberFreeParameters();

  printf("output: %s\nvertex cuts: %s\n", name.Data(), vertex_label.Data());
  printf("input: %s\nselected ULS pairs: %.0f (%.0f in %.2f-%.2f GeV)\n", input_file.Data(), h->GetEntries(),
         h->Integral(), plot_lo, plot_hi);
  printf("fit %.5f-%.5f GeV, %.3f MeV bins, status %d\n", fit_lo, fit_hi, 1000 * kBinWidth, (int)result);
  printf("M      = %.2f +- %.2f MeV\n", 1000 * p[1], 1000 * fit->GetParError(1));
  printf("sigma  = %.2f +- %.2f MeV   (Gamma fixed %.3f MeV)\n", 1000 * p[2], 1000 * fit->GetParError(2),
         1000 * kGammaPhi);
  printf("N_sig (total) = %.1f +- %.1f\n", p[0], fit->GetParError(0));
  printf("window %.3f-%.3f: S = %.1f +- %.1f, B = %.1f +- %.1f, S+B = %.1f (observed %.0f)\n", win_lo,
         win_hi, S, S_err, B, B_err, S + B, observed);
  printf("S/sqrt(S+B) = %.2f\n", signif);
  printf("chi2/ndf (Baker-Cousins) = %.1f/%d = %.2f, Pearson = %.1f/%d = %.2f, p = %.2f\n", gof.bc, ndf,
         gof.bc / ndf, gof.pearson, ndf, gof.pearson / ndf, TMath::Prob(gof.bc, ndf));

  // Drawing: Times bold, ULS black / LS red with blue error bars, total fit red dashed,
  // background blue dashed, pulls in a separate lower pad. Text sizes are in pixels (font
  // precision 3), so both pads share them.
  TH1D* h_ls_fine = new TH1D("phi_mass_ls_fine", "", kFine * kNbins, plot_lo, plot_hi);
  h_ls_fine->Sumw2(false);
  TTree* ls_tree = dynamic_cast<TTree*>(file->Get("ls_data"));
  TString ls_cut = cut;
  ls_cut.ReplaceAll("uls_", "ls_");
  ls_tree->Draw("ls_phi_mass>>phi_mass_ls_fine", ls_cut, "goff");
  TH1D* h_ls = dynamic_cast<TH1D*>(h_ls_fine->Rebin(kFine, "phi_mass_ls"));
  printf("LS pairs: %.0f (%.0f in window)\n", h_ls->GetEntries(),
         h_ls_fine->Integral(h_ls_fine->FindBin(win_lo + 1e-9), h_ls_fine->FindBin(win_hi - 1e-9)));

  // Not drawn: the threshold bin containing 0.987, where LS split-track pairs (one track found as
  // two fragments on disjoint pad rows, so M = 2 m_K) pile up. It is outside the fit range, so no
  // fit result changes. Empty bins are not drawn with option "E".
  const double skip_mass = 0.987;
  for (TH1D* hist : {h, h_ls}) {
    const int bin = hist->FindBin(skip_mass);
    printf("skipped bin %.4f-%.4f in %s: %.0f pairs\n", hist->GetBinLowEdge(bin), hist->GetBinLowEdge(bin + 1),
           hist->GetName(), hist->GetBinContent(bin));
    hist->SetBinContent(bin, 0);
    hist->SetBinError(bin, 0);
  }

  const int err_color = kBlue + 2;
  const double label_px = 22, title_px = 26, text_px = 24;
  TCanvas* canvas = new TCanvas("phi_mass_fit_canvas", "", 800, 780);
  TPad* top = new TPad("top_pad", "", 0, 0.32, 1, 1);
  TPad* bottom = new TPad("pull_pad", "", 0, 0, 1, 0.30);
  top->SetMargin(0.13, 0.075, 0.07, 0.04);
  bottom->SetMargin(0.13, 0.075, 0.32, 0.05);
  top->Draw();
  bottom->Draw();
  // Tick lengths are fractions of the pad height, so rescale them for the shorter pull pad.
  auto style_axis = [&](TAxis* axis, double tick) {
    axis->SetLabelFont(23);
    axis->SetLabelSize(label_px);
    axis->SetTitleFont(23);
    axis->SetTitleSize(title_px);
    axis->SetTickLength(tick);
  };

  top->cd();
  double ymax = 0;
  for (int b = 1; b <= h->GetNbinsX(); ++b) ymax = std::max(ymax, h->GetBinContent(b) + h->GetBinError(b));
  for (TH1D* hist : {h, h_ls}) {
    hist->SetStats(false);
    hist->SetLineColor(err_color);
    hist->SetLineWidth(2);
    hist->SetMarkerStyle(20);
    hist->SetMarkerSize(0.9);
  }
  h->SetMarkerColor(kBlack);
  h_ls->SetMarkerColor(kRed);
  h->SetMinimum(0);
  h->SetMaximum(1.35 * ymax);
  style_axis(h->GetXaxis(), 0.03);
  style_axis(h->GetYaxis(), 0.02);
  h->GetXaxis()->SetNdivisions(506);
  h->GetYaxis()->SetTitle("Counts");
  h->GetYaxis()->CenterTitle();
  h->GetYaxis()->SetTitleOffset(1.7);
  h->Draw("E");
  h_ls->Draw("E SAME");

  TF1* bkg_draw = new TF1("phi_bkg_draw", Background, fit_lo, fit_hi, 6);
  bkg_draw->SetParameters(p);
  for (TF1* f : {fit, bkg_draw}) {
    f->SetNpx(1000);
    f->SetLineWidth(3);
    f->SetLineStyle(2);
  }
  fit->SetLineColor(kRed);
  bkg_draw->SetLineColor(kBlue);
  bkg_draw->Draw("SAME");
  fit->Draw("SAME");
  h->Draw("E SAME");

  TLegend* legend = new TLegend(0.60, 0.76, 0.78, 0.87);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->SetTextFont(23);
  legend->SetTextSize(text_px);
  legend->AddEntry(h, "K^{#pm}K^{#mp}", "lep");
  legend->AddEntry(h_ls, "K^{#pm}K^{#pm}", "lep");
  legend->Draw();

  TLatex text;
  text.SetNDC();
  text.SetTextFont(33);
  text.SetTextSize(text_px);
  const TString lines[] = {
    TString::Format("#chi^{2}/NDF = %.4f / %d", gof.bc, ndf),
    TString::Format("p-value = %.4f", TMath::Prob(gof.bc, ndf)),
    TString::Format("N_{sig} = %.2f  #pm %.2f", S, S_err),
    TString::Format("M_{KK} = %.4f  #pm %.4f", p[1], fit->GetParError(1)),
    TString::Format("#sigma = %.4f  #pm %.4f", p[2], fit->GetParError(2)),
    TString::Format("#Gamma = %.4f  #pm %.4f", p[3], fit->GetParError(3)),
    TString::Format("S/#sqrt{S+B} = %.4f", signif)};
  const double x = 0.52;
  double y = 0.70;
  text.SetTextFont(23);
  text.DrawLatex(0.17, 0.89, vertex_label);
  text.DrawLatex(x, y, "All Runs");
  text.SetTextFont(33);
  for (const TString& line : lines) {
    y -= 0.058;
    text.DrawLatex(x, y, line);
  }

  bottom->cd();
  const double tick_scale = (1 - 0.32) / 0.30; // top pad height / pull pad height
  TH1F* pull_frame = bottom->DrawFrame(plot_lo, -3, plot_hi, 3);
  style_axis(pull_frame->GetXaxis(), 0.03 * tick_scale);
  style_axis(pull_frame->GetYaxis(), 0.02);
  pull_frame->GetXaxis()->SetNdivisions(506);
  pull_frame->GetYaxis()->SetNdivisions(506);
  pull_frame->GetXaxis()->SetTitle("M_{KK} (GeV/c^{2})");
  pull_frame->GetXaxis()->SetTitleOffset(1.1);
  pull_frame->GetYaxis()->SetTitle("Pull");
  pull_frame->GetYaxis()->CenterTitle();
  pull_frame->GetYaxis()->SetTitleOffset(1.7);
  for (double level : {-2., -1., 0., 1., 2.}) {
    TLine* line = new TLine(plot_lo, level, plot_hi, level);
    line->SetLineStyle(3);
    line->SetLineWidth(2);
    line->Draw();
  }
  gof.pull->SetMarkerStyle(20);
  gof.pull->SetMarkerSize(0.9);
  gof.pull->SetMarkerColor(kRed);
  gof.pull->Draw("P SAME");

  canvas->Modified();
  canvas->Update();
  canvas->Print(name + ".pdf");
  // ROOT's batch PNG writer restarts the dash pattern on every curve segment, so dashed fits come
  // out solid. Instead crop the PDF (ROOT puts it on a larger page) and render it with
  // Ghostscript, whose standard fonts keep Times bold; fall back to ROOT's PNG without them.
  if (gSystem->Exec(TString::Format("pdfcrop --margins 10 %s.pdf %s.pdf > /dev/null && "
                                    "gs -q -dSAFER -dBATCH -dNOPAUSE -sDEVICE=png16m -r250 -dTextAlphaBits=4 "
                                    "-dGraphicsAlphaBits=4 -sOutputFile=%s.png %s.pdf",
                                    name.Data(), name.Data(), name.Data(), name.Data())) != 0)
    canvas->Print(name + ".png");
}
