#ifndef DI_PION_H
#define DI_PION_H

// Shared di-pion helpers (vz_pion.C, BBC_veto.C): tree projection, the modified Soding
// rho0 mass model, the ULS - LS mass histogram, its fit, and the fit plot with components.
#include "VzFit.h"
#include "TColor.h"
#include "TFitResult.h"
#include "TMatrixDSym.h"
#include "TTree.h"
#include <complex>
#include <vector>

// Fill h from one charge tree, replacing the "P" prefix in the expression and the cut.
void Project(TTree* tree, TString prefix, TH1* h, TString expr, TString cut) {
  expr.ReplaceAll("P_", prefix + "_");
  cut.ReplaceAll("P_", prefix + "_");
  tree->Project(h->GetName(), expr, cut);
}

// Modified Soding parametrization (STAR, PRC 96 (2017) 054904, Eq. 2-5):
// |A rho_BW + B + C e^{i phi} omega_BW|^2 + p0 + p1 M, with A, B, C real.
// p: 0 M_rho, 1 Gamma_rho, 2 A, 3 B, 4 C, 5 M_omega, 6 Gamma_omega, 7 phi, 8 p0, 9 p1.
const double kMpi = 0.13957;
struct Amplitudes { std::complex<double> rho, omega; double B; };
Amplitudes SodingAmplitudes(double M, const double* p) {
  using C = std::complex<double>;
  if (M <= 2 * kMpi) return {0, 0, 0};
  const double Mr = p[0], Mw = p[5];
  const double Gr = p[1] * (Mr / M) * std::pow((M*M - 4*kMpi*kMpi) / (Mr*Mr - 4*kMpi*kMpi), 1.5);
  const double Gw = M > 3 * kMpi  // 3-pion threshold of the total omega width
                    ? p[6] * (Mw / M) * std::pow((M*M - 9*kMpi*kMpi) / (Mw*Mw - 9*kMpi*kMpi), 1.5) : 0;
  const double Gwpp = 0.0153 * p[6] * (Mw / M) *
                      std::pow((M*M - 4*kMpi*kMpi) / (Mw*Mw - 4*kMpi*kMpi), 1.5);
  return {p[2] * std::sqrt(M * Mr * Gr) / C(M*M - Mr*Mr, Mr * Gr),
          p[4] * std::polar(1.0, p[7]) * std::sqrt(M * Mw * Gwpp) / C(M*M - Mw*Mw, Mw * Gw),
          p[3]};
}
double SodingFull(double* x, double* p) {
  const Amplitudes a = SodingAmplitudes(x[0], p);
  return std::norm(a.rho + a.B + a.omega) + p[8] + p[9] * x[0];
}
// Everything but the remnant background, and the rho Breit-Wigner alone (same 10 parameters,
// so the fit covariance can be used for integral errors).
double SodingSignal(double* x, double* p) {
  const Amplitudes a = SodingAmplitudes(x[0], p);
  return std::norm(a.rho + a.B + a.omega);
}
double SodingRho(double* x, double* p) { return std::norm(SodingAmplitudes(x[0], p).rho); }
// Amplitude ratio p[j] / A_rho (B_pipi / A: j = 3, C_omega / A: j = 4) and its uncertainty from the
// fit covariance, including the correlation with A.
std::pair<double, double> AmplitudeRatio(const TF1* fit, const TMatrixDSym& cov, int j) {
  const double x = fit->GetParameter(j), a = fit->GetParameter(2);
  const double rel2 = cov(j, j) / (x * x) + cov(2, 2) / (a * a) - 2 * cov(j, 2) / (x * a);
  return {x / a, std::abs(x / a) * std::sqrt(std::max(rel2, 0.))};
}

// Terms of |rho + B + omega|^2, selected by p[10]; together they sum to the full amplitude squared.
enum SodingTerm { kRhoBW, kRhoPiPi, kOmegaBW, kOmegaRho, kOmegaPiPi, kNonRes, kNTerms };
double SodingPart(double* x, double* p) {
  const Amplitudes a = SodingAmplitudes(x[0], p);
  switch (static_cast<int>(p[10])) {
    case kRhoBW:     return std::norm(a.rho);
    case kRhoPiPi:   return 2 * a.B * a.rho.real();
    case kOmegaBW:   return std::norm(a.omega);
    case kOmegaRho:  return 2 * (a.rho * std::conj(a.omega)).real();
    case kOmegaPiPi: return 2 * a.B * a.omega.real();
    case kNonRes:    return a.B * a.B;
  }
  return 0;
}

// ULS - LS mass spectrum as dN/dM over [lo, hi]: 2.5 MeV bins across the rho-omega interference
// region, 10 MeV elsewhere. The fit range is set by the TF1, so the histogram may extend below it.
TH1D* MassSpectrum(TTree* uls_tree, TTree* ls_tree, TString cut, TString name,
                   double lo = 0.5, double hi = 1.2) {
  std::vector<double> edges;
  for (double m = lo; m < 0.72 - 1e-9; m += 0.01) edges.push_back(m);
  for (double m = 0.72; m < 0.84 - 1e-9; m += 0.0025) edges.push_back(m);
  for (double m = 0.84; m < hi + 1e-9; m += 0.01) edges.push_back(m);
  TH1D* uls = new TH1D(name + "_uls", "", edges.size() - 1, edges.data());
  TH1D* ls = new TH1D(name + "_ls", "", edges.size() - 1, edges.data());
  Project(uls_tree, "uls", uls, "P_phi_mass", cut);
  Project(ls_tree, "ls", ls, "P_phi_mass", cut);
  uls->Sumw2();
  ls->Sumw2();
  TH1D* sub = static_cast<TH1D*>(uls->Clone(name));
  sub->Add(ls, -1);
  sub->Scale(1, "width");
  return sub;
}

// Soding fit with omega mass and width floating; the omega phase has local minima, so keep the
// best of several starts and refit from it for consistent errors and chi2.
TFitResultPtr FitSoding(TH1D* sub, TF1* fit) {
  double best_chi2 = -1;
  std::vector<double> best(10);
  const double a0 = std::sqrt(sub->GetBinContent(sub->FindBin(0.76)) * 0.15);
  for (double phase : {-2.0, 0.0, 1.0, 2.5}) {
    fit->SetParameters(0.775, 0.15, a0, -0.64 * a0, 0.3 * a0, 0.78266, 0.00868, phase, 0, 0);
    fit->SetParLimits(0, 0.74, 0.81);
    fit->SetParLimits(1, 0.10, 0.22);
    fit->SetParLimits(4, 0, 10 * a0);
    fit->SetParLimits(5, 0.77, 0.795);
    fit->SetParLimits(6, 0.002, 0.06);
    fit->SetParLimits(7, -TMath::Pi(), TMath::Pi());
    sub->Fit(fit, "QRN0");
    if (best_chi2 < 0 || fit->GetChisquare() < best_chi2) {
      best_chi2 = fit->GetChisquare();
      for (int i = 0; i < 10; ++i) best[i] = fit->GetParameter(i);
    }
  }
  fit->SetParameters(best.data());
  return sub->Fit(fit, "QRNS0");
}

// Fit plot with all Soding components (as in Fig. 4 of the paper). The data are shown over the
// histogram range, the fit and its components over the fit range only. summary is printed beside
// the legend, selection below it.
void DrawSodingFit(TH1D* sub, TF1* fit, PlotStyle::Options opt, const std::vector<TString>& summary,
                   const std::vector<TString>& selection, TString name) {
  const double lo = fit->GetXmin(), hi = fit->GetXmax();
  struct Component { int term; TString label; int color; int style; };
  const Component components[kNTerms] = {
    {kRhoBW, "#rho^{0} BW", PlotStyle::Color(0), 1},
    {kRhoPiPi, "#rho^{0} interference #pi#pi", PlotStyle::Color(4), 2},
    {kOmegaBW, "#omega BW", PlotStyle::Color(2), 1},
    {kOmegaRho, "#omega interference #rho^{0}", TColor::GetColor("#E69F00"), 2},
    {kOmegaPiPi, "#omega interference #pi#pi", PlotStyle::Color(5), 3},
    {kNonRes, "Non-resonant #pi#pi", PlotStyle::Color(3), 2}};
  std::vector<TF1*> parts;
  double ylow = 0;
  for (const Component& c : components) {
    TF1* part = new TF1(TString::Format("%s_part_%d", name.Data(), c.term), SodingPart, lo, hi, 11);
    part->SetParameters(fit->GetParameters());
    part->SetParameter(10, c.term);
    part->SetNpx(1000);
    part->SetLineColor(c.color);
    part->SetLineStyle(c.style);
    part->SetLineWidth(3);
    ylow = std::min(ylow, part->GetMinimum(lo, hi));
    parts.push_back(part);
  }
  TF1* fp = new TF1(name + "_fp", "[0]+[1]*x", lo, hi);
  fp->SetParameters(fit->GetParameter(8), fit->GetParameter(9));
  fp->SetLineColor(kGray + 1);
  fp->SetLineStyle(9);
  fp->SetLineWidth(3);
  ylow = std::min(ylow, std::min(fp->Eval(lo), fp->Eval(hi)));
  ylow *= 1.15;

  opt.xlabel = "M_{#pi#pi} [GeV/c^{2}]";
  opt.ylabel = "dN/dM_{#pi#pi} [(GeV/c^{2})^{-1}]";
  opt.logY = false;
  opt.showEntries = false;
  opt.legendX1 = 0.49;
  opt.legendX2 = 0.69;
  opt.legendTextSize = 0.027;
  double ymax = 0;
  for (int bin = 1; bin <= sub->GetNbinsX(); ++bin)
    ymax = std::max(ymax, sub->GetBinContent(bin) + sub->GetBinError(bin));
  opt.YRange(ylow, ylow + (ymax - ylow) / 0.72);
  TCanvas* canvas = PlotStyle::Draw({{sub, "ULS #minus LS", "E1", kBlack}}, opt);
  if (!canvas) return;
  canvas->cd();

  TLine* zero = new TLine(sub->GetXaxis()->GetXmin(), 0, sub->GetXaxis()->GetXmax(), 0);
  zero->SetLineColor(kGray + 1);
  zero->Draw();
  for (TF1* part : parts) part->Draw("SAME");
  fp->Draw("SAME");
  fit->SetNpx(1000);
  fit->SetLineColor(PlotStyle::Color(1));
  fit->SetLineWidth(4);
  fit->Draw("SAME");

  TLegend* legend = nullptr;
  for (TObject* obj : *canvas->GetListOfPrimitives())
    if ((legend = dynamic_cast<TLegend*>(obj))) break;
  if (legend) {
    legend->SetY1NDC(opt.legendY2 - 0.04 * (kNTerms + 3));
    legend->AddEntry(fit, "S#ddot{o}ding fit", "l");
    for (int i = 0; i < kNTerms; ++i) legend->AddEntry(parts[i], components[i].label, "l");
    legend->AddEntry(fp, "Remnant bkg.", "l");
  }

  const double x = opt.legendX2 + 0.01;
  const double y = opt.legendY2 - 0.005;
  DrawLines(summary, x, y, 0.027);
  DrawLines(selection, x, y - 1.3 * 0.027 * summary.size() - 0.015);
  canvas->Modified();
  canvas->Update();
  canvas->Print(name + ".png");
  canvas->Print(name + ".pdf");
}

#endif
