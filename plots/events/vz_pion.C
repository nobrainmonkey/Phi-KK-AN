#include "DiPion.h"
#include "TFile.h"
#include "TH2D.h"
#include "TStyle.h"
#include "TTree.h"
#include <memory>
#include <vector>

// Di-pion version of vz.C, from the UPCDst rho -> pi pi trees (uls_data, ls_data).
// Tree preselection (analysis_definition in the file): trigger 700001, bad runs removed,
// exactly two primary tracks on one vertex, x2pipi < 50, pion pT < 1 GeV/c,
// 0.5 < M < 1.3 GeV/c^2, pair pT < 1 GeV/c. Branch names keep the phi/kaon labels.
// Coordinates are the pair vertex (phi_vx/vy/vz); vr is about (0,0).
//
// Di-pion selection added here ("P" is replaced by uls / ls):
const TString kSelection =
  "P_pair_BBC_Veto_East==0&&P_pair_BBC_Veto_West==0"      // BBC veto
  "&&P_pair_TofMult<=6&&P_pair_n_tof_matched>=1"           // >= 1 TOF-matched track
  "&&abs(P_kaon1_nsigma_pion)<3&&abs(P_kaon2_nsigma_pion)<3" // pion dE/dx
  "&&P_phi_pt<0.1";                                        // coherent region
const TString kMassWindow = "&&P_phi_mass>0.6&&P_phi_mass<0.95"; // around the rho pole
const TString kVr = "sqrt(P_phi_vx*P_phi_vx+P_phi_vy*P_phi_vy)";

// Selection as printed on the plots.
const std::vector<TString> kSelectionLines = {
  "2 primary tracks",
  "TofMult #leq 6, N_{TOF match} #geq 1",
  "|n#sigma_{#pi}| < 3, BBC veto",
  "p_{T}^{#pi#pi} < 0.1 GeV/c"};
std::vector<TString> WithMassWindow() {
  std::vector<TString> lines = kSelectionLines;
  lines.push_back("0.60 < M_{#pi#pi} < 0.95 GeV/c^{2}");
  return lines;
}

// Mass spectrum of the selected di-pions (no vertex or mass-window cut): ULS vs LS,
// and the Soding fit to ULS - LS, to show the sample is rho dominated.
void FitMass(TTree* uls_tree, TTree* ls_tree, PlotStyle::Options opt) {
  const double plot_lo = 0.5, lo = 0.6, hi = 1.2; // shown from plot_lo, fitted from lo
  const int nbins = std::lround((hi - plot_lo) / 0.01);
  TH1D* uls = new TH1D("mass_uls", "", nbins, plot_lo, hi);
  TH1D* ls = new TH1D("mass_ls", "", nbins, plot_lo, hi);
  Project(uls_tree, "uls", uls, "P_phi_mass", kSelection);
  Project(ls_tree, "ls", ls, "P_phi_mass", kSelection);
  uls->Sumw2();
  ls->Sumw2();

  opt.xlabel = "M_{#pi#pi} [GeV/c^{2}]";
  opt.ylabel = "Pairs / 10 MeV/c^{2}";
  opt.logY = true;
  TCanvas* c_ls = PlotStyle::Draw({{uls, "ULS"}, {ls, "LS"}}, opt);
  if (!c_ls) return;
  DrawLines(kSelectionLines, opt.legendX1 + 0.05, opt.legendY2 - 0.051 * 2 - 0.02);
  c_ls->Print("vz_pion_mass.png");
  c_ls->Print("vz_pion_mass.pdf");

  TH1D* sub = MassSpectrum(uls_tree, ls_tree, kSelection, "mass_sub", plot_lo, hi);
  TF1* fit = new TF1("mass_fit", SodingFull, lo, hi, 10);
  FitSoding(sub, fit);

  // Composition inside the mass window used for the vz study.
  const double wlo = 0.6, whi = 0.95;
  const int b1 = uls->FindBin(wlo + 1e-6), b2 = uls->FindBin(whi - 1e-6);
  const double n_uls = uls->Integral(b1, b2), n_ls = ls->Integral(b1, b2);
  const double n_fp = fit->GetParameter(8) * (whi - wlo) +
                      0.5 * fit->GetParameter(9) * (whi*whi - wlo*wlo); // dN/dM integral
  const double purity = 1 - (n_ls + n_fp) / n_uls;
  printf("mass fit: chi2/ndf=%.1f/%d M_rho=%.4f G_rho=%.4f B/A=%.3f C/A=%.3f M_w=%.4f G_w=%.4f phi=%.2f"
         " | window %.2f-%.2f: ULS=%.0f LS=%.0f f_p=%.0f purity=%.4f\n",
         fit->GetChisquare(), fit->GetNDF(), fit->GetParameter(0), fit->GetParameter(1),
         fit->GetParameter(3) / fit->GetParameter(2), fit->GetParameter(4) / fit->GetParameter(2),
         fit->GetParameter(5), fit->GetParameter(6), fit->GetParameter(7),
         wlo, whi, n_uls, n_ls, n_fp, purity);

  DrawSodingFit(sub, fit, opt, {
    TString::Format("M_{#rho} = %.4f, #Gamma_{#rho} = %.3f GeV/c^{2}", fit->GetParameter(0), fit->GetParameter(1)),
    TString::Format("#Gamma_{#omega} = %.1f MeV/c^{2}, #chi^{2}/ndf = %.0f/%d", 1000 * fit->GetParameter(6),
                    fit->GetChisquare(), fit->GetNDF()),
    TString::Format("%.2f < M_{#pi#pi} < %.2f GeV/c^{2}:", wlo, whi),
    TString::Format("N_{LS}/N_{ULS} = %.4f", n_ls / n_uls),
    TString::Format("1 #minus (N_{LS}+f_{p})/N_{ULS} = %.4f", purity)},
    kSelectionLines, "vz_pion_mass_fit");
}

void vz_pion(TString input_file = "../../data/rho_pipi_upcdst.root") {
  std::unique_ptr<TFile> file(TFile::Open(input_file, "READ"));
  TTree* uls_tree = dynamic_cast<TTree*>(file->Get("uls_data"));
  TTree* ls_tree = dynamic_cast<TTree*>(file->Get("ls_data"));
  gStyle->SetPalette(kBird);

  const TString selection = kSelection + kMassWindow;

  PlotStyle::Options opt;
  opt.title = "Run-19 Au+Au";
  opt.subtitle = "#sqrt{s_{NN}} = 200 GeV";
  opt.note = "ZDCMB";

  FitMass(uls_tree, ls_tree, opt);

  // vr of the pair vertex: All / ULS / LS.
  TH1D* vr_uls = new TH1D("vr_uls", "", 150, 0, 15);
  TH1D* vr_ls = new TH1D("vr_ls", "", 150, 0, 15);
  Project(uls_tree, "uls", vr_uls, kVr, selection);
  Project(ls_tree, "ls", vr_ls, kVr, selection);
  TH1D* vr_all = static_cast<TH1D*>(vr_uls->Clone("vr_all"));
  vr_all->Add(vr_ls);

  PlotStyle::Options opt_vr = opt;
  opt_vr.xlabel = "v_{r} [cm]";
  opt_vr.ylabel = "Events / 0.1 cm";
  opt_vr.logY = true;
  opt_vr.headroom = 1.9; // room for the selection lines
  TCanvas* c_vr = PlotStyle::Draw({
    {vr_all, "All"},
    {vr_uls, "ULS"},
    {vr_ls, "LS"}
  }, opt_vr);
  if (!c_vr) return;
  DrawLines(WithMassWindow(), opt.annotationX, BelowNote(opt));
  c_vr->Print("vz_pion_vr.png");
  c_vr->Print("vz_pion_vr.pdf");

  // vz before and after vr < 3 cm, All / ULS / LS on one canvas each.
  const TString vr_cut[2] = {"", "&&" + kVr + "<3"};
  const TString vz_note[2] = {"ZDCMB",
                              "ZDCMB, v_{r} < 3 cm"};
  const TString vz_name[2] = {"vz_pion_vz", "vz_pion_vz_vr3"};
  for (int i = 0; i < 2; ++i) {
    TH1D* vz_uls = new TH1D(vz_name[i] + "_uls", "", 100, -250, 250);
    TH1D* vz_ls = new TH1D(vz_name[i] + "_ls", "", 100, -250, 250);
    Project(uls_tree, "uls", vz_uls, "P_phi_vz", selection + vr_cut[i]);
    Project(ls_tree, "ls", vz_ls, "P_phi_vz", selection + vr_cut[i]);
    TH1D* vz_all = static_cast<TH1D*>(vz_uls->Clone(vz_name[i] + "_all"));
    vz_all->Add(vz_ls);

    PlotStyle::Options opt_vz = opt;
    opt_vz.xlabel = "v_{z} [cm]";
    opt_vz.ylabel = "Events / 5 cm";
    opt_vz.note = vz_note[i];
    TCanvas* c_vz = PlotStyle::Draw({
      {vz_all, "All"},
      {vz_uls, "ULS"},
      {vz_ls, "LS"}
    }, opt_vz);
    if (!c_vz) return;
    DrawLines(WithMassWindow(), opt.annotationX, BelowNote(opt));
    c_vz->Print(vz_name[i] + ".png");
    c_vz->Print(vz_name[i] + ".pdf");
  }

  // Double-Gaussian vz fits: ULS only and ULS - LS, with and without vr < 3 cm.
  const TString fit_note[2] = {"ZDCMB, no v_{r} cut",
                               "ZDCMB, v_{r} < 3 cm"};
  const TString fit_name[2] = {"vz_pion_fit", "vz_pion_fit_vr3"};
  for (int i = 0; i < 2; ++i) {
    TH1D* fit_uls = new TH1D(fit_name[i] + "_uls_h", "", 80, -200, 200);
    TH1D* fit_ls = new TH1D(fit_name[i] + "_ls_h", "", 80, -200, 200);
    Project(uls_tree, "uls", fit_uls, "P_phi_vz", selection + vr_cut[i]);
    Project(ls_tree, "ls", fit_ls, "P_phi_vz", selection + vr_cut[i]);

    PlotStyle::Options opt_fit = opt;
    opt_fit.note = fit_note[i];
    FitVzPair(fit_uls, fit_ls, opt_fit, fit_name[i], WithMassWindow());
  }

  PlotStyle::Options opt2 = opt;
  opt2.xlabel = "v_{x} [cm]";
  opt2.ylabel = "v_{y} [cm]";
  opt2.logY = true;
  opt2.showEntries = false; // entries are listed with the selection, top right

  // vx:vy panels (ULS + LS unless noted); x spans are 1.2x the y spans so circles stay round.
  struct Panel { TString name, cut, note; int nx, ny; double x, ylo, yhi; bool uls_only; };
  const Panel panels[4] = {
    {"vz_pion_vxvy", "", "ZDCMB", 360, 300, 18, -11, 19, false},
    {"vz_pion_vxvy_ring", "&&" + kVr + ">2&&" + kVr + "<6", "ZDCMB, 2 < v_{r} < 6 cm",
     228, 190, 11.4, -6, 13, false},
    {"vz_pion_vxvy_vr3", "&&" + kVr + "<3", "ZDCMB, v_{r} < 3 cm",
     234, 195, 7.8, -4, 9, false},
    {"vz_pion_vxvy_vr3_uls", "&&" + kVr + "<3", "ZDCMB, v_{r} < 3 cm, ULS",
     234, 195, 7.8, -4, 9, true}};
  for (const Panel& panel : panels) {
    TH2D* xy = new TH2D(panel.name + "_h", "", panel.nx, -panel.x, panel.x,
                        panel.ny, panel.ylo, panel.yhi);
    Project(uls_tree, "uls", xy, "P_phi_vy:P_phi_vx", selection + panel.cut);
    if (!panel.uls_only) {
      TH2D* xy_ls = static_cast<TH2D*>(xy->Clone(panel.name + "_ls_h"));
      xy_ls->Reset();
      Project(ls_tree, "ls", xy_ls, "P_phi_vy:P_phi_vx", selection + panel.cut);
      xy->Add(xy_ls);
    }
    opt2.note = panel.note;
    Draw2D(xy, opt2);
    std::vector<TString> lines = WithMassWindow();
    lines.push_back(("Entries = " + PlotStyle::ShortNumber(xy->GetEntries())).c_str());
    DrawLines(lines, kDraw2DRight, opt2.annotationY, 0.028, 33);
    gPad->Print(panel.name + ".png");
    gPad->Print(panel.name + ".pdf");
  }
}
