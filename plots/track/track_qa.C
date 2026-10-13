#include "../events/VzFit.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TLine.h"
#include "TProfile.h"
#include "TString.h"
#include "TSystem.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TTreeReaderValue.h"
#include <cmath>
#include <iostream>
#include <map>
#include <memory>

// Track, pair and event QA for the ULS and LS pair trees of phi-kk-mudst.root, built on
// track_skim.C. Selection: BBC veto (east and west), exactly 2 primary tracks in the event and
// |Vz| < 100 cm; nothing else beyond the production cuts. Both tracks of every pair are filled.
// Output: PDF files in qa/{track,track2d,pair,event}/.

const double kC = 29.9792458; // speed of light [cm/ns]

// Expected dE/dx of each hypothesis vs p, taken from the data itself: mean dE/dx of tracks whose
// nsigma_X is within kExpWindow of zero, i.e. tracks sitting on the reference curve that nsigma_X
// is computed from (Eq. nsigma of the note). Filled from ULS and LS tracks together.
const double kExpWindow = 0.3;
const char* kExpName[4] = {"e", "#pi", "K", "p"};
const int kExpColor[4] = {kMagenta + 1, kRed + 1, kBlack, kGray + 2};
TProfile* gExpected[4] = {nullptr, nullptr, nullptr, nullptr};

void FillExpected(double p, double de, const double nsigma[4]) {
  for (int i = 0; i < 4; ++i) {
    if (!gExpected[i]) gExpected[i] = new TProfile(TString::Format("expected_%d", i), "", 100, 0, 2.5);
    if (std::abs(nsigma[i]) < kExpWindow) gExpected[i]->Fill(p, de);
  }
}

// Expected-dE/dx curve as a graph; p bins with fewer than 5 tracks are skipped.
TGraph* ExpectedGraph(int i) {
  TGraph* g = new TGraph();
  for (int b = 1; b <= gExpected[i]->GetNbinsX(); ++b)
    if (gExpected[i]->GetBinEntries(b) >= 5)
      g->SetPoint(g->GetN(), gExpected[i]->GetBinCenter(b), gExpected[i]->GetBinContent(b));
  g->SetLineColor(kExpColor[i]);
  g->SetLineWidth(3);
  g->SetLineStyle(i == 2 ? 1 : 2);
  return g;
}

// One track of a pair.
struct Track {
  int nhits, nhitsfit, nhitsdedx, nhitsposs, charge;
  double pt, eta, phi, dedx; // dE/dx in GeV/cm as stored
  double nsk, nspi, nse, nsp;
  double dca[3], gdca[3];    // primary / global DCA to the vertex [cm]
  float tof_time, tof_path;  // -999 if not matched
  ULong64_t itpc_map;         // iTPC (inner) pad-row topology bits
  UInt_t map0, map1;          // main topology map, low and high 32 bits
};

// One histogram per name; ULS and LS get separate sets with the same binning.
struct Hists {
  TString tag;
  std::map<TString, TH1*> h;
  explicit Hists(const TString& t) : tag(t) {}
  TH1* H1(const TString& n, int nb, double lo, double hi) {
    auto& p = h[n];
    if (!p) { p = new TH1D(tag + "_" + n, "", nb, lo, hi); p->Sumw2(false); }
    return p;
  }
  TH1* H2(const TString& n, int nx, double xlo, double xhi, int ny, double ylo, double yhi) {
    auto& p = h[n];
    if (!p) p = new TH2D(tag + "_" + n, "", nx, xlo, xhi, ny, ylo, yhi);
    return p;
  }
  void F1(const TString& n, int nb, double lo, double hi, double x) { H1(n, nb, lo, hi)->Fill(x); }
  void F2(const TString& n, int nx, double xlo, double xhi, int ny, double ylo, double yhi,
          double x, double y) {
    static_cast<TH2*>(H2(n, nx, xlo, xhi, ny, ylo, yhi))->Fill(x, y);
  }

  void FillTrack(const Track& t) {
    const double p = t.pt * std::cosh(t.eta), de = t.dedx * 1e6; // keV/cm
    const double dcaxy = std::hypot(t.dca[0], t.dca[1]);
    const double dca3 = std::sqrt(dcaxy * dcaxy + t.dca[2] * t.dca[2]);
    const double gdcaxy = std::hypot(t.gdca[0], t.gdca[1]);
    const double gdca3 = std::sqrt(gdcaxy * gdcaxy + t.gdca[2] * t.gdca[2]);
    const bool tof = t.tof_time > 0 && t.tof_path > 0;

    F1("nhits", 80, 0, 80, t.nhits);
    F1("nhitsfit", 80, 0, 80, t.nhitsfit);
    F1("nhitsdedx", 80, 0, 80, t.nhitsdedx);
    F1("nhitsposs", 80, 0, 80, t.nhitsposs);
    // Pad rows in iTPC numbering: rows 1-40 (inner) are bits 1-40 of the iTPC word, rows 41-72
    // (outer) are bits 21-52 of the main map (StTrackTopologyMap::iTpcBit).
    const ULong64_t main = (ULong64_t(t.map1) << 32) | t.map0;
    int n_inner = 0, n_outer = 0;
    for (int row = 1; row <= 72; ++row) {
      const bool hit = row <= 40 ? (t.itpc_map >> row & 1ULL) : (main >> (row - 20) & 1ULL);
      if (!hit) continue;
      F1("row_occupancy", 72, 0.5, 72.5, row);
      (row <= 40 ? n_inner : n_outer)++;
    }
    F1("itpc_rows", 41, -0.5, 40.5, n_inner);
    F1("outer_rows", 33, -0.5, 32.5, n_outer);
    F2("inner_outer_rows", 41, -0.5, 40.5, 33, -0.5, 32.5, n_inner, n_outer);
    F1("nhitsfit_ratio", 100, 0, 1.2, t.nhitsposs > 0 ? double(t.nhitsfit) / t.nhitsposs : -1);
    F1("nhitsdedx_ratio", 100, 0, 1.2, t.nhitsfit > 0 ? double(t.nhitsdedx) / t.nhitsfit : -1);
    F1("charge", 3, -1.5, 1.5, t.charge);
    F1("pt", 100, 0, 1.0, t.pt);
    F1("p", 125, 0, 2.5, p);
    F1("eta", 100, -2.5, 2.5, t.eta);
    F1("phi", 90, -M_PI, M_PI, t.phi);
    F1("dedx", 200, 0, 40, de);
    F1("nsigma_kaon", 160, -20, 20, t.nsk);
    F1("nsigma_pion", 160, -20, 20, t.nspi);
    F1("nsigma_electron", 160, -20, 20, t.nse);
    F1("nsigma_proton", 160, -20, 20, t.nsp);
    F1("dca_x", 200, -0.5, 0.5, t.dca[0]);
    F1("dca_y", 200, -0.5, 0.5, t.dca[1]);
    F1("dca_z", 200, -0.5, 0.5, t.dca[2]);
    F1("dca_xy", 200, 0, 0.5, dcaxy);
    F1("dca_3d", 200, 0, 0.5, dca3);
    F1("gdca_x", 200, -4, 4, t.gdca[0]);
    F1("gdca_y", 200, -4, 4, t.gdca[1]);
    F1("gdca_z", 200, -4, 4, t.gdca[2]);
    F1("gdca_xy", 200, 0, 4, gdcaxy);
    F1("gdca_3d", 200, 0, 4, gdca3);
    F1("tof_matched", 2, -0.5, 1.5, tof);
    if (tof) {
      F1("tof_invbeta", 150, 0.5, 3.5, kC * t.tof_time / t.tof_path);
      F2("p_invbeta", 250, 0, 2.5, 150, 0.5, 3.5, p, kC * t.tof_time / t.tof_path);
    }

    F2("p_dedx", 250, 0, 2.5, 200, 0, 40, p, de);
    const double ns[4] = {t.nse, t.nspi, t.nsk, t.nsp};
    FillExpected(p, de, ns);
    F2("qp_dedx", 500, -2.5, 2.5, 200, 0, 40, t.charge * p, de);
    F2("p_nsigma_kaon", 250, 0, 2.5, 168, -12, 30, p, t.nsk);
    F2("p_nsigma_pion", 250, 0, 2.5, 160, -20, 20, p, t.nspi);
    F2("eta_phi", 100, -2.5, 2.5, 90, -M_PI, M_PI, t.eta, t.phi);
    F2("eta_pt", 100, -2.5, 2.5, 100, 0, 1.0, t.eta, t.pt);
    F2("eta_nhitsfit", 100, -2.5, 2.5, 80, 0, 80, t.eta, t.nhitsfit);
    F2("pt_nhitsfit", 100, 0, 1.0, 80, 0, 80, t.pt, t.nhitsfit);
    F2("nhitsfit_nhitsdedx", 80, 0, 80, 80, 0, 80, t.nhitsfit, t.nhitsdedx);
    F2("nhitsfit_gdca_3d", 80, 0, 80, 100, 0, 4, t.nhitsfit, gdca3);
    F2("eta_gdca_3d", 100, -2.5, 2.5, 100, 0, 4, t.eta, gdca3);
    F2("gdca_xy_z", 200, -4, 4, 200, -4, 4, t.gdca[2], gdcaxy);
    F2("dca_xy_z", 200, -0.5, 0.5, 200, 0, 0.5, t.dca[2], dcaxy);
  }

  void FillPair(double pt, double y, double mass, double x2kk, double x2pipi, double x2ee,
                double x2pp, double nsk1, double nsk2) {
    F1("pair_pt", 100, 0, 1.0, pt);
    F1("pair_pt2", 100, 0, 1.0, pt * pt);
    F1("pair_rapidity", 100, -2.5, 2.5, y);
    F1("pair_mass", 220, 0.98, 1.2, mass);
    F1("pair_x2kk", 100, 0, 100, x2kk);
    F1("pair_x2pipi", 200, 0, 400, x2pipi);
    F1("pair_x2ee", 200, 0, 400, x2ee);
    F1("pair_x2pp", 200, 0, 400, x2pp);
    F2("pair_nsigma_kaon_12", 160, -20, 20, 160, -20, 20, nsk1, nsk2);
    F2("pair_x2kk_x2pipi", 100, 0, 100, 200, 0, 400, x2kk, x2pipi);
    F2("pair_mass_pt", 220, 0.98, 1.2, 100, 0, 1.0, mass, pt);
  }

  void FillEvent(double vx, double vy, double vz, int nvtx, int nglobal, int tofmult,
                 int ntrk_used, double ranking, int zdc_e, int zdc_w) {
    F1("event_vz", 100, -100, 100, vz);
    F1("event_vr", 100, 0, 5, std::hypot(vx, vy));
    F2("event_vx_vy", 200, -5, 5, 200, -5, 5, vx, vy);
    F1("event_nvtx", 10, 0.5, 10.5, nvtx);
    F1("event_n_global_tracks", 100, 0, 200, nglobal);
    F1("event_TofMult", 30, -1.5, 28.5, tofmult);
    F1("event_vertex_n_tracks_used", 10, -0.5, 9.5, ntrk_used);
    F1("event_vertex_ranking", 100, -4, 3, ranking);
    F1("event_ZDC_UA_East", 200, 0, 2000, zdc_e);
    F1("event_ZDC_UA_West", 200, 0, 2000, zdc_w);
  }
};

// ---- Plotting -------------------------------------------------------------------
const std::vector<TString> kSelection = {"BBC veto, N_{prim} = 2", "|V_{z}| < 100 cm"};

PlotStyle::Options BaseOptions() {
  PlotStyle::Options opt;
  opt.title = "Run-19 Au+Au";
  opt.subtitle = "#sqrt{s_{NN}} = 200 GeV";
  opt.note = "ZDCMB";
  opt.legendX1 = 0.56;
  return opt;
}

void Overlay(Hists& uls, Hists& ls, const TString& name, const TString& xlabel,
             const TString& ylabel, const TString& out, bool logY = false) {
  PlotStyle::Options opt = BaseOptions();
  opt.xlabel = xlabel;
  opt.ylabel = ylabel;
  opt.logY = logY;
  TCanvas* c = PlotStyle::Draw({{uls.h.at(name), "ULS"}, {ls.h.at(name), "LS"}}, opt);
  if (!c) return;
  DrawLines(kSelection, opt.annotationX, BelowNote(opt), 0.032);
  c->Print(out + ".pdf");
}

// Pair pT with the nominal coherent (pT <= 0.1) and incoherent (0.3 < pT < 0.9 GeV/c) regions.
void PtRegions(Hists& uls, Hists& ls, const TString& out) {
  PlotStyle::Options opt = BaseOptions();
  opt.xlabel = "Pair p_{T} [GeV/c]";
  opt.ylabel = "Pairs";
  opt.logY = true;
  TCanvas* c = PlotStyle::Draw({{uls.h.at("pair_pt"), "ULS"}, {ls.h.at("pair_pt"), "LS"}}, opt);
  if (!c) return;
  c->Update();
  const double ylo = std::pow(10, c->GetUymin()), yhi = std::pow(10, c->GetUymax());
  const double ytop = std::pow(10, c->GetUymin() + 0.62 * (c->GetUymax() - c->GetUymin()));
  struct Region { double lo, hi; int color; const char* label; };
  for (const Region& r : {Region{0.0, 0.1, kAzure + 1, "Coherent"}, Region{0.3, 0.9, kOrange + 1, "Incoherent"}}) {
    for (double x : {r.lo, r.hi}) {
      if (x <= 0) continue;
      TLine* l = new TLine(x, ylo, x, ytop);
      l->SetLineStyle(2);
      l->SetLineColor(kGray + 2);
      l->Draw();
    }
    TLatex t;
    t.SetTextFont(42);
    t.SetTextSize(0.034);
    t.SetTextAlign(11);
    t.SetTextColor(r.color + 2);
    // Labels sit at the bottom, right of each region's lower edge, clear of the histograms.
    // Arrows point from the boundary line into the region.
    if (r.lo > 0) t.DrawLatex(r.lo + 0.01, ylo * 2.0, TString(r.label) + " #rightarrow");
    else          t.DrawLatex(r.hi + 0.01, ylo * 2.0, TString("#leftarrow ") + r.label);
  }
  (void)yhi;
  DrawLines(kSelection, opt.annotationX, BelowNote(opt), 0.032);
  c->RedrawAxis();
  c->Print(out + ".pdf");
}

void Plot2D(Hists& hs, const TString& name, const TString& xlabel, const TString& ylabel,
            const TString& out, double annotationX = 0.18, bool expected = false,
            double legX1 = 0.55, double legY1 = 0.40) {
  PlotStyle::Options opt = BaseOptions();
  opt.xlabel = xlabel;
  opt.ylabel = ylabel;
  opt.logY = true; // log z
  opt.note = hs.tag == "uls" ? "ZDCMB, ULS" : hs.tag == "ls" ? "ZDCMB, LS" : "ZDCMB, ULS + LS";
  opt.annotationX = annotationX;
  TCanvas* c = Draw2D(static_cast<TH2*>(hs.h.at(name)), opt);
  DrawLines(kSelection, annotationX, BelowNote(opt) - 0.050, 0.032);
  if (expected) {
    TLegend* leg = new TLegend(legX1, legY1, legX1 + 0.26, legY1 + 0.13);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.030);
    leg->SetNColumns(2);
    leg->SetHeader("Expected dE/dx");
    for (int i = 0; i < 4; ++i) {
      TGraph* g = ExpectedGraph(i);
      if (g->GetN() == 0) continue;
      g->Draw("L SAME");
      leg->AddEntry(g, kExpName[i], "l");
    }
    leg->Draw();
  }
  c->Print(out + "_" + hs.tag + ".pdf");
}

void track_qa(TString input_file = "../../data/phi-kk-mudst.root") {
  std::unique_ptr<TFile> file(TFile::Open(input_file, "READ"));
  if (!file || file->IsZombie()) return;

  // ---- ULS pairs -------------------------------------------------------------
  TTreeReader uls_reader("uls_data", file.get());
  TTreeReaderValue<Double_t> uls_phi_pt(uls_reader, "uls_phi_pt");
  TTreeReaderValue<Double_t> uls_phi_eta(uls_reader, "uls_phi_eta");
  TTreeReaderValue<Double_t> uls_phi_rapidity(uls_reader, "uls_phi_rapidity");
  TTreeReaderValue<Double_t> uls_phi_x2kk(uls_reader, "uls_phi_x2kk");
  TTreeReaderValue<Double_t> uls_phi_x2pipi(uls_reader, "uls_phi_x2pipi");
  TTreeReaderValue<Double_t> uls_phi_x2ee(uls_reader, "uls_phi_x2ee");
  TTreeReaderValue<Double_t> uls_phi_x2pp(uls_reader, "uls_phi_x2pp");
  TTreeReaderValue<Double_t> uls_phi_vz(uls_reader, "uls_phi_vz");
  TTreeReaderValue<Double_t> uls_phi_vx(uls_reader, "uls_phi_vx");
  TTreeReaderValue<Double_t> uls_phi_vy(uls_reader, "uls_phi_vy");
  TTreeReaderValue<Double_t> uls_phi_mass(uls_reader, "uls_phi_mass");
  TTreeReaderValue<Double_t> uls_phi_azimuthal_angle(uls_reader, "uls_phi_azimuthal_angle");
  TTreeReaderValue<Double_t> uls_kaon1_pt(uls_reader, "uls_kaon1_pt");
  TTreeReaderValue<Double_t> uls_kaon1_eta(uls_reader, "uls_kaon1_eta");
  TTreeReaderValue<Double_t> uls_kaon1_rapidity(uls_reader, "uls_kaon1_rapidity");
  TTreeReaderValue<Double_t> uls_kaon1_azimuthal_angle(uls_reader, "uls_kaon1_azimuthal_angle");
  TTreeReaderValue<Double_t> uls_kaon1_dca_x(uls_reader, "uls_kaon1_dca_x");
  TTreeReaderValue<Double_t> uls_kaon1_dca_y(uls_reader, "uls_kaon1_dca_y");
  TTreeReaderValue<Double_t> uls_kaon1_dca_z(uls_reader, "uls_kaon1_dca_z");
  TTreeReaderValue<Double_t> uls_kaon1_dca_global_x(uls_reader, "uls_kaon1_dca_global_x");
  TTreeReaderValue<Double_t> uls_kaon1_dca_global_y(uls_reader, "uls_kaon1_dca_global_y");
  TTreeReaderValue<Double_t> uls_kaon1_dca_global_z(uls_reader, "uls_kaon1_dca_global_z");
  TTreeReaderValue<Double_t> uls_kaon2_pt(uls_reader, "uls_kaon2_pt");
  TTreeReaderValue<Double_t> uls_kaon2_eta(uls_reader, "uls_kaon2_eta");
  TTreeReaderValue<Double_t> uls_kaon2_rapidity(uls_reader, "uls_kaon2_rapidity");
  TTreeReaderValue<Double_t> uls_kaon2_azimuthal_angle(uls_reader, "uls_kaon2_azimuthal_angle");
  TTreeReaderValue<Double_t> uls_kaon2_dca_x(uls_reader, "uls_kaon2_dca_x");
  TTreeReaderValue<Double_t> uls_kaon2_dca_y(uls_reader, "uls_kaon2_dca_y");
  TTreeReaderValue<Double_t> uls_kaon2_dca_z(uls_reader, "uls_kaon2_dca_z");
  TTreeReaderValue<Double_t> uls_kaon2_dca_global_x(uls_reader, "uls_kaon2_dca_global_x");
  TTreeReaderValue<Double_t> uls_kaon2_dca_global_y(uls_reader, "uls_kaon2_dca_global_y");
  TTreeReaderValue<Double_t> uls_kaon2_dca_global_z(uls_reader, "uls_kaon2_dca_global_z");
  TTreeReaderValue<UInt_t> uls_kaon1_topology_map0(uls_reader, "uls_kaon1_topology_map0");
  TTreeReaderValue<UInt_t> uls_kaon1_topology_map1(uls_reader, "uls_kaon1_topology_map1");
  TTreeReaderValue<ULong64_t> uls_kaon1_topology_map_itpc(uls_reader, "uls_kaon1_topology_map_itpc");
  TTreeReaderValue<UShort_t> uls_kaon1_nhits(uls_reader, "uls_kaon1_nhits");
  TTreeReaderValue<UInt_t> uls_kaon2_topology_map0(uls_reader, "uls_kaon2_topology_map0");
  TTreeReaderValue<UInt_t> uls_kaon2_topology_map1(uls_reader, "uls_kaon2_topology_map1");
  TTreeReaderValue<ULong64_t> uls_kaon2_topology_map_itpc(uls_reader, "uls_kaon2_topology_map_itpc");
  TTreeReaderValue<UShort_t> uls_kaon2_nhits(uls_reader, "uls_kaon2_nhits");
  TTreeReaderValue<UShort_t> uls_kaon1_nhitsDEdx(uls_reader, "uls_kaon1_nhitsDEdx");
  TTreeReaderValue<UShort_t> uls_kaon2_nhitsDEdx(uls_reader, "uls_kaon2_nhitsDEdx");
  TTreeReaderValue<UShort_t> uls_kaon1_nhitsposs(uls_reader, "uls_kaon1_nhitsposs");
  TTreeReaderValue<UShort_t> uls_kaon2_nhitsposs(uls_reader, "uls_kaon2_nhitsposs");
  TTreeReaderValue<UShort_t> uls_kaon1_nhitsfit(uls_reader, "uls_kaon1_nhitsfit");
  TTreeReaderValue<UShort_t> uls_kaon2_nhitsfit(uls_reader, "uls_kaon2_nhitsfit");
  TTreeReaderValue<Double_t> uls_kaon1_nsigma_kaon(uls_reader, "uls_kaon1_nsigma_kaon");
  TTreeReaderValue<Double_t> uls_kaon1_nsigma_pion(uls_reader, "uls_kaon1_nsigma_pion");
  TTreeReaderValue<Double_t> uls_kaon1_nsigma_proton(uls_reader, "uls_kaon1_nsigma_proton");
  TTreeReaderValue<Double_t> uls_kaon1_nsigma_electron(uls_reader, "uls_kaon1_nsigma_electron");
  TTreeReaderValue<Double_t> uls_kaon2_nsigma_kaon(uls_reader, "uls_kaon2_nsigma_kaon");
  TTreeReaderValue<Double_t> uls_kaon2_nsigma_pion(uls_reader, "uls_kaon2_nsigma_pion");
  TTreeReaderValue<Double_t> uls_kaon2_nsigma_proton(uls_reader, "uls_kaon2_nsigma_proton");
  TTreeReaderValue<Double_t> uls_kaon2_nsigma_electron(uls_reader, "uls_kaon2_nsigma_electron");
  TTreeReaderValue<UShort_t> uls_pair_ZDC_UA_East(uls_reader, "uls_pair_ZDC_UA_East");
  TTreeReaderValue<UShort_t> uls_pair_ZDC_UA_West(uls_reader, "uls_pair_ZDC_UA_West");
  TTreeReaderValue<UShort_t> uls_pair_BBC_Small_East(uls_reader, "uls_pair_BBC_Small_East");
  TTreeReaderValue<UShort_t> uls_pair_BBC_Small_West(uls_reader, "uls_pair_BBC_Small_West");
  TTreeReaderValue<UInt_t> uls_pair_BBC_Small_ADC_East_raw(uls_reader, "uls_pair_BBC_Small_ADC_East_raw");
  TTreeReaderValue<UInt_t> uls_pair_BBC_Small_ADC_West_raw(uls_reader, "uls_pair_BBC_Small_ADC_West_raw");
  TTreeReaderValue<UInt_t> uls_pair_BBC_Small_TDC_East_raw(uls_reader, "uls_pair_BBC_Small_TDC_East_raw");
  TTreeReaderValue<UInt_t> uls_pair_BBC_Small_TDC_West_raw(uls_reader, "uls_pair_BBC_Small_TDC_West_raw");
  TTreeReaderValue<UShort_t> uls_pair_BBC_Veto_East(uls_reader, "uls_pair_BBC_Veto_East");
  TTreeReaderValue<UShort_t> uls_pair_BBC_Veto_West(uls_reader, "uls_pair_BBC_Veto_West");
  TTreeReaderValue<UShort_t> uls_pair_TofMult(uls_reader, "uls_pair_TofMult");
  TTreeReaderValue<Float_t> uls_kaon1_tof_time(uls_reader, "uls_kaon1_tof_time");
  TTreeReaderValue<Float_t> uls_kaon1_tof_pathlength(uls_reader, "uls_kaon1_tof_pathlength");
  TTreeReaderValue<Float_t> uls_kaon2_tof_time(uls_reader, "uls_kaon2_tof_time");
  TTreeReaderValue<Float_t> uls_kaon2_tof_pathlength(uls_reader, "uls_kaon2_tof_pathlength");
  TTreeReaderValue<Int_t> uls_pair_event_nvtx(uls_reader, "uls_pair_event_nvtx");
  TTreeReaderValue<UInt_t> uls_pair_event_n_primary_tracks(uls_reader, "uls_pair_event_n_primary_tracks");
  TTreeReaderValue<UInt_t> uls_pair_event_n_global_tracks(uls_reader, "uls_pair_event_n_global_tracks");
  TTreeReaderValue<UShort_t> uls_pair_has_trigger_data(uls_reader, "uls_pair_has_trigger_data");
  TTreeReaderValue<Float_t> uls_pair_vertex_ranking(uls_reader, "uls_pair_vertex_ranking");
  TTreeReaderValue<Int_t> uls_pair_vertex_index(uls_reader, "uls_pair_vertex_index");
  TTreeReaderValue<UShort_t> uls_pair_is_highest_ranking(uls_reader, "uls_pair_is_highest_ranking");
  TTreeReaderValue<Double_t> uls_pair_event_vx(uls_reader, "uls_pair_event_vx");
  TTreeReaderValue<Double_t> uls_pair_event_vy(uls_reader, "uls_pair_event_vy");
  TTreeReaderValue<Double_t> uls_pair_event_vz(uls_reader, "uls_pair_event_vz");
  TTreeReaderValue<Int_t> uls_pair_run_number(uls_reader, "uls_pair_run_number");
  TTreeReaderValue<Int_t> uls_pair_event_number(uls_reader, "uls_pair_event_number");
  TTreeReaderValue<Double_t> uls_kaon1_dEdx(uls_reader, "uls_kaon1_dEdx");
  TTreeReaderValue<Double_t> uls_kaon2_dEdx(uls_reader, "uls_kaon2_dEdx");
  TTreeReaderValue<Short_t> uls_kaon1_trackid(uls_reader, "uls_kaon1_trackid");
  TTreeReaderValue<Short_t> uls_kaon2_trackid(uls_reader, "uls_kaon2_trackid");
  TTreeReaderValue<Short_t> uls_kaon1_charge(uls_reader, "uls_kaon1_charge");
  TTreeReaderValue<Short_t> uls_kaon2_charge(uls_reader, "uls_kaon2_charge");
  TTreeReaderValue<Float_t> uls_pair_vertex_error_x(uls_reader, "uls_pair_vertex_error_x");
  TTreeReaderValue<Float_t> uls_pair_vertex_error_y(uls_reader, "uls_pair_vertex_error_y");
  TTreeReaderValue<Float_t> uls_pair_vertex_error_z(uls_reader, "uls_pair_vertex_error_z");
  TTreeReaderValue<Float_t> uls_pair_vertex_chi2(uls_reader, "uls_pair_vertex_chi2");
  TTreeReaderValue<UShort_t> uls_pair_vertex_n_tracks_used(uls_reader, "uls_pair_vertex_n_tracks_used");
  TTreeReaderValue<Int_t> uls_pair_vertex_n_primary_tracks(uls_reader, "uls_pair_vertex_n_primary_tracks");
  TTreeReaderValue<Int_t> uls_pair_vertex_finder_id(uls_reader, "uls_pair_vertex_finder_id");
  TTreeReaderValue<Int_t> uls_pair_vertex_flag(uls_reader, "uls_pair_vertex_flag");
  TTreeReaderValue<UShort_t> uls_pair_has_vertex(uls_reader, "uls_pair_has_vertex");
  TTreeReaderValue<UShort_t> uls_pair_trigger_700001(uls_reader, "uls_pair_trigger_700001");
  TTreeReaderValue<UShort_t> uls_pair_detector_valid(uls_reader, "uls_pair_detector_valid");
  TTreeReaderValue<UInt_t> uls_pair_trigger_error_flag(uls_reader, "uls_pair_trigger_error_flag");
  TTreeReaderValue<UShort_t> uls_pair_has_zdc_data(uls_reader, "uls_pair_has_zdc_data");
  TTreeReaderValue<UShort_t> uls_pair_lastDSM1(uls_reader, "uls_pair_lastDSM1");
  TTreeReaderValue<Double_t> uls_pair_BBC_CoincidenceRate(uls_reader, "uls_pair_BBC_CoincidenceRate");
  TTreeReaderValue<Double_t> uls_pair_ZDC_CoincidenceRate(uls_reader, "uls_pair_ZDC_CoincidenceRate");
  TTreeReaderValue<UInt_t> uls_pair_n_uls_pairs(uls_reader, "uls_pair_n_uls_pairs");
  TTreeReaderValue<UInt_t> uls_pair_n_ls_pairs(uls_reader, "uls_pair_n_ls_pairs");
  TTreeReaderArray<UShort_t> uls_pair_BBC_Small_ADC_East(uls_reader, "uls_pair_BBC_Small_ADC_East");
  TTreeReaderArray<UShort_t> uls_pair_BBC_Small_ADC_West(uls_reader, "uls_pair_BBC_Small_ADC_West");
  TTreeReaderArray<UShort_t> uls_pair_BBC_Small_TDC_East(uls_reader, "uls_pair_BBC_Small_TDC_East");
  TTreeReaderArray<UShort_t> uls_pair_BBC_Small_TDC_West(uls_reader, "uls_pair_BBC_Small_TDC_West");

  Hists uls_h("uls");
  while (uls_reader.Next()) {
    if (*uls_pair_BBC_Veto_East != 0 || *uls_pair_BBC_Veto_West != 0) continue;
    if (*uls_pair_event_n_primary_tracks != 2) continue;
    if (std::abs(*uls_pair_event_vz) >= 100) continue;
    {
      Track t{*uls_kaon1_nhits, *uls_kaon1_nhitsfit, *uls_kaon1_nhitsDEdx, *uls_kaon1_nhitsposs,
              *uls_kaon1_charge, *uls_kaon1_pt, *uls_kaon1_eta, *uls_kaon1_azimuthal_angle,
              *uls_kaon1_dEdx, *uls_kaon1_nsigma_kaon, *uls_kaon1_nsigma_pion,
              *uls_kaon1_nsigma_electron, *uls_kaon1_nsigma_proton,
              {*uls_kaon1_dca_x, *uls_kaon1_dca_y, *uls_kaon1_dca_z},
              {*uls_kaon1_dca_global_x, *uls_kaon1_dca_global_y, *uls_kaon1_dca_global_z},
              *uls_kaon1_tof_time, *uls_kaon1_tof_pathlength, *uls_kaon1_topology_map_itpc,
              *uls_kaon1_topology_map0, *uls_kaon1_topology_map1};
      uls_h.FillTrack(t);
    }
    {
      Track t{*uls_kaon2_nhits, *uls_kaon2_nhitsfit, *uls_kaon2_nhitsDEdx, *uls_kaon2_nhitsposs,
              *uls_kaon2_charge, *uls_kaon2_pt, *uls_kaon2_eta, *uls_kaon2_azimuthal_angle,
              *uls_kaon2_dEdx, *uls_kaon2_nsigma_kaon, *uls_kaon2_nsigma_pion,
              *uls_kaon2_nsigma_electron, *uls_kaon2_nsigma_proton,
              {*uls_kaon2_dca_x, *uls_kaon2_dca_y, *uls_kaon2_dca_z},
              {*uls_kaon2_dca_global_x, *uls_kaon2_dca_global_y, *uls_kaon2_dca_global_z},
              *uls_kaon2_tof_time, *uls_kaon2_tof_pathlength, *uls_kaon2_topology_map_itpc,
              *uls_kaon2_topology_map0, *uls_kaon2_topology_map1};
      uls_h.FillTrack(t);
    }
    uls_h.FillPair(*uls_phi_pt, *uls_phi_rapidity, *uls_phi_mass, *uls_phi_x2kk, *uls_phi_x2pipi,
                *uls_phi_x2ee, *uls_phi_x2pp, *uls_kaon1_nsigma_kaon, *uls_kaon2_nsigma_kaon);
    uls_h.FillEvent(*uls_pair_event_vx, *uls_pair_event_vy, *uls_pair_event_vz, *uls_pair_event_nvtx,
                 *uls_pair_event_n_global_tracks, *uls_pair_TofMult, *uls_pair_vertex_n_tracks_used,
                 *uls_pair_vertex_ranking, *uls_pair_ZDC_UA_East, *uls_pair_ZDC_UA_West);
  }

  // ---- LS pairs --------------------------------------------------------------
  TTreeReader ls_reader("ls_data", file.get());
  TTreeReaderValue<Double_t> ls_phi_pt(ls_reader, "ls_phi_pt");
  TTreeReaderValue<Double_t> ls_phi_eta(ls_reader, "ls_phi_eta");
  TTreeReaderValue<Double_t> ls_phi_rapidity(ls_reader, "ls_phi_rapidity");
  TTreeReaderValue<Double_t> ls_phi_vz(ls_reader, "ls_phi_vz");
  TTreeReaderValue<Double_t> ls_phi_vx(ls_reader, "ls_phi_vx");
  TTreeReaderValue<Double_t> ls_phi_vy(ls_reader, "ls_phi_vy");
  TTreeReaderValue<Double_t> ls_phi_mass(ls_reader, "ls_phi_mass");
  TTreeReaderValue<Double_t> ls_phi_azimuthal_angle(ls_reader, "ls_phi_azimuthal_angle");
  TTreeReaderValue<Double_t> ls_phi_x2kk(ls_reader, "ls_phi_x2kk");
  TTreeReaderValue<Double_t> ls_phi_x2pipi(ls_reader, "ls_phi_x2pipi");
  TTreeReaderValue<Double_t> ls_phi_x2ee(ls_reader, "ls_phi_x2ee");
  TTreeReaderValue<Double_t> ls_phi_x2pp(ls_reader, "ls_phi_x2pp");
  TTreeReaderValue<Double_t> ls_kaon1_pt(ls_reader, "ls_kaon1_pt");
  TTreeReaderValue<Double_t> ls_kaon1_eta(ls_reader, "ls_kaon1_eta");
  TTreeReaderValue<Double_t> ls_kaon1_rapidity(ls_reader, "ls_kaon1_rapidity");
  TTreeReaderValue<Double_t> ls_kaon1_azimuthal_angle(ls_reader, "ls_kaon1_azimuthal_angle");
  TTreeReaderValue<Double_t> ls_kaon1_dca_x(ls_reader, "ls_kaon1_dca_x");
  TTreeReaderValue<Double_t> ls_kaon1_dca_y(ls_reader, "ls_kaon1_dca_y");
  TTreeReaderValue<Double_t> ls_kaon1_dca_z(ls_reader, "ls_kaon1_dca_z");
  TTreeReaderValue<Double_t> ls_kaon1_dca_global_x(ls_reader, "ls_kaon1_dca_global_x");
  TTreeReaderValue<Double_t> ls_kaon1_dca_global_y(ls_reader, "ls_kaon1_dca_global_y");
  TTreeReaderValue<Double_t> ls_kaon1_dca_global_z(ls_reader, "ls_kaon1_dca_global_z");
  TTreeReaderValue<Double_t> ls_kaon2_pt(ls_reader, "ls_kaon2_pt");
  TTreeReaderValue<Double_t> ls_kaon2_eta(ls_reader, "ls_kaon2_eta");
  TTreeReaderValue<Double_t> ls_kaon2_rapidity(ls_reader, "ls_kaon2_rapidity");
  TTreeReaderValue<Double_t> ls_kaon2_azimuthal_angle(ls_reader, "ls_kaon2_azimuthal_angle");
  TTreeReaderValue<Double_t> ls_kaon2_dca_x(ls_reader, "ls_kaon2_dca_x");
  TTreeReaderValue<Double_t> ls_kaon2_dca_y(ls_reader, "ls_kaon2_dca_y");
  TTreeReaderValue<Double_t> ls_kaon2_dca_z(ls_reader, "ls_kaon2_dca_z");
  TTreeReaderValue<Double_t> ls_kaon2_dca_global_x(ls_reader, "ls_kaon2_dca_global_x");
  TTreeReaderValue<Double_t> ls_kaon2_dca_global_y(ls_reader, "ls_kaon2_dca_global_y");
  TTreeReaderValue<Double_t> ls_kaon2_dca_global_z(ls_reader, "ls_kaon2_dca_global_z");
  TTreeReaderValue<UInt_t> ls_kaon1_topology_map0(ls_reader, "ls_kaon1_topology_map0");
  TTreeReaderValue<UInt_t> ls_kaon1_topology_map1(ls_reader, "ls_kaon1_topology_map1");
  TTreeReaderValue<ULong64_t> ls_kaon1_topology_map_itpc(ls_reader, "ls_kaon1_topology_map_itpc");
  TTreeReaderValue<UShort_t> ls_kaon1_nhits(ls_reader, "ls_kaon1_nhits");
  TTreeReaderValue<UInt_t> ls_kaon2_topology_map0(ls_reader, "ls_kaon2_topology_map0");
  TTreeReaderValue<UInt_t> ls_kaon2_topology_map1(ls_reader, "ls_kaon2_topology_map1");
  TTreeReaderValue<ULong64_t> ls_kaon2_topology_map_itpc(ls_reader, "ls_kaon2_topology_map_itpc");
  TTreeReaderValue<UShort_t> ls_kaon2_nhits(ls_reader, "ls_kaon2_nhits");
  TTreeReaderValue<UShort_t> ls_kaon1_nhitsDEdx(ls_reader, "ls_kaon1_nhitsDEdx");
  TTreeReaderValue<UShort_t> ls_kaon2_nhitsDEdx(ls_reader, "ls_kaon2_nhitsDEdx");
  TTreeReaderValue<UShort_t> ls_kaon1_nhitsposs(ls_reader, "ls_kaon1_nhitsposs");
  TTreeReaderValue<UShort_t> ls_kaon2_nhitsposs(ls_reader, "ls_kaon2_nhitsposs");
  TTreeReaderValue<UShort_t> ls_kaon1_nhitsfit(ls_reader, "ls_kaon1_nhitsfit");
  TTreeReaderValue<UShort_t> ls_kaon2_nhitsfit(ls_reader, "ls_kaon2_nhitsfit");
  TTreeReaderValue<Double_t> ls_kaon1_nsigma_kaon(ls_reader, "ls_kaon1_nsigma_kaon");
  TTreeReaderValue<Double_t> ls_kaon1_nsigma_pion(ls_reader, "ls_kaon1_nsigma_pion");
  TTreeReaderValue<Double_t> ls_kaon1_nsigma_proton(ls_reader, "ls_kaon1_nsigma_proton");
  TTreeReaderValue<Double_t> ls_kaon1_nsigma_electron(ls_reader, "ls_kaon1_nsigma_electron");
  TTreeReaderValue<Double_t> ls_kaon2_nsigma_kaon(ls_reader, "ls_kaon2_nsigma_kaon");
  TTreeReaderValue<Double_t> ls_kaon2_nsigma_pion(ls_reader, "ls_kaon2_nsigma_pion");
  TTreeReaderValue<Double_t> ls_kaon2_nsigma_proton(ls_reader, "ls_kaon2_nsigma_proton");
  TTreeReaderValue<Double_t> ls_kaon2_nsigma_electron(ls_reader, "ls_kaon2_nsigma_electron");
  TTreeReaderValue<UShort_t> ls_pair_ZDC_UA_East(ls_reader, "ls_pair_ZDC_UA_East");
  TTreeReaderValue<UShort_t> ls_pair_ZDC_UA_West(ls_reader, "ls_pair_ZDC_UA_West");
  TTreeReaderValue<UShort_t> ls_pair_BBC_Small_East(ls_reader, "ls_pair_BBC_Small_East");
  TTreeReaderValue<UShort_t> ls_pair_BBC_Small_West(ls_reader, "ls_pair_BBC_Small_West");
  TTreeReaderValue<UInt_t> ls_pair_BBC_Small_ADC_East_raw(ls_reader, "ls_pair_BBC_Small_ADC_East_raw");
  TTreeReaderValue<UInt_t> ls_pair_BBC_Small_ADC_West_raw(ls_reader, "ls_pair_BBC_Small_ADC_West_raw");
  TTreeReaderValue<UInt_t> ls_pair_BBC_Small_TDC_East_raw(ls_reader, "ls_pair_BBC_Small_TDC_East_raw");
  TTreeReaderValue<UInt_t> ls_pair_BBC_Small_TDC_West_raw(ls_reader, "ls_pair_BBC_Small_TDC_West_raw");
  TTreeReaderValue<UShort_t> ls_pair_BBC_Veto_East(ls_reader, "ls_pair_BBC_Veto_East");
  TTreeReaderValue<UShort_t> ls_pair_BBC_Veto_West(ls_reader, "ls_pair_BBC_Veto_West");
  TTreeReaderValue<UShort_t> ls_pair_TofMult(ls_reader, "ls_pair_TofMult");
  TTreeReaderValue<Float_t> ls_kaon1_tof_time(ls_reader, "ls_kaon1_tof_time");
  TTreeReaderValue<Float_t> ls_kaon1_tof_pathlength(ls_reader, "ls_kaon1_tof_pathlength");
  TTreeReaderValue<Float_t> ls_kaon2_tof_time(ls_reader, "ls_kaon2_tof_time");
  TTreeReaderValue<Float_t> ls_kaon2_tof_pathlength(ls_reader, "ls_kaon2_tof_pathlength");
  TTreeReaderValue<Int_t> ls_pair_event_nvtx(ls_reader, "ls_pair_event_nvtx");
  TTreeReaderValue<UInt_t> ls_pair_event_n_primary_tracks(ls_reader, "ls_pair_event_n_primary_tracks");
  TTreeReaderValue<UInt_t> ls_pair_event_n_global_tracks(ls_reader, "ls_pair_event_n_global_tracks");
  TTreeReaderValue<UShort_t> ls_pair_has_trigger_data(ls_reader, "ls_pair_has_trigger_data");
  TTreeReaderValue<Float_t> ls_pair_vertex_ranking(ls_reader, "ls_pair_vertex_ranking");
  TTreeReaderValue<Int_t> ls_pair_vertex_index(ls_reader, "ls_pair_vertex_index");
  TTreeReaderValue<UShort_t> ls_pair_is_highest_ranking(ls_reader, "ls_pair_is_highest_ranking");
  TTreeReaderValue<Double_t> ls_pair_event_vx(ls_reader, "ls_pair_event_vx");
  TTreeReaderValue<Double_t> ls_pair_event_vy(ls_reader, "ls_pair_event_vy");
  TTreeReaderValue<Double_t> ls_pair_event_vz(ls_reader, "ls_pair_event_vz");
  TTreeReaderValue<Int_t> ls_pair_run_number(ls_reader, "ls_pair_run_number");
  TTreeReaderValue<Int_t> ls_pair_event_number(ls_reader, "ls_pair_event_number");
  TTreeReaderValue<Double_t> ls_kaon1_dEdx(ls_reader, "ls_kaon1_dEdx");
  TTreeReaderValue<Double_t> ls_kaon2_dEdx(ls_reader, "ls_kaon2_dEdx");
  TTreeReaderValue<Short_t> ls_kaon1_trackid(ls_reader, "ls_kaon1_trackid");
  TTreeReaderValue<Short_t> ls_kaon2_trackid(ls_reader, "ls_kaon2_trackid");
  TTreeReaderValue<Short_t> ls_kaon1_charge(ls_reader, "ls_kaon1_charge");
  TTreeReaderValue<Short_t> ls_kaon2_charge(ls_reader, "ls_kaon2_charge");
  TTreeReaderValue<Float_t> ls_pair_vertex_error_x(ls_reader, "ls_pair_vertex_error_x");
  TTreeReaderValue<Float_t> ls_pair_vertex_error_y(ls_reader, "ls_pair_vertex_error_y");
  TTreeReaderValue<Float_t> ls_pair_vertex_error_z(ls_reader, "ls_pair_vertex_error_z");
  TTreeReaderValue<Float_t> ls_pair_vertex_chi2(ls_reader, "ls_pair_vertex_chi2");
  TTreeReaderValue<UShort_t> ls_pair_vertex_n_tracks_used(ls_reader, "ls_pair_vertex_n_tracks_used");
  TTreeReaderValue<Int_t> ls_pair_vertex_n_primary_tracks(ls_reader, "ls_pair_vertex_n_primary_tracks");
  TTreeReaderValue<Int_t> ls_pair_vertex_finder_id(ls_reader, "ls_pair_vertex_finder_id");
  TTreeReaderValue<Int_t> ls_pair_vertex_flag(ls_reader, "ls_pair_vertex_flag");
  TTreeReaderValue<UShort_t> ls_pair_has_vertex(ls_reader, "ls_pair_has_vertex");
  TTreeReaderValue<UShort_t> ls_pair_trigger_700001(ls_reader, "ls_pair_trigger_700001");
  TTreeReaderValue<UShort_t> ls_pair_detector_valid(ls_reader, "ls_pair_detector_valid");
  TTreeReaderValue<UInt_t> ls_pair_trigger_error_flag(ls_reader, "ls_pair_trigger_error_flag");
  TTreeReaderValue<UShort_t> ls_pair_has_zdc_data(ls_reader, "ls_pair_has_zdc_data");
  TTreeReaderValue<UShort_t> ls_pair_lastDSM1(ls_reader, "ls_pair_lastDSM1");
  TTreeReaderValue<Double_t> ls_pair_BBC_CoincidenceRate(ls_reader, "ls_pair_BBC_CoincidenceRate");
  TTreeReaderValue<Double_t> ls_pair_ZDC_CoincidenceRate(ls_reader, "ls_pair_ZDC_CoincidenceRate");
  TTreeReaderValue<UInt_t> ls_pair_n_uls_pairs(ls_reader, "ls_pair_n_uls_pairs");
  TTreeReaderValue<UInt_t> ls_pair_n_ls_pairs(ls_reader, "ls_pair_n_ls_pairs");
  TTreeReaderArray<UShort_t> ls_pair_BBC_Small_ADC_East(ls_reader, "ls_pair_BBC_Small_ADC_East");
  TTreeReaderArray<UShort_t> ls_pair_BBC_Small_ADC_West(ls_reader, "ls_pair_BBC_Small_ADC_West");
  TTreeReaderArray<UShort_t> ls_pair_BBC_Small_TDC_East(ls_reader, "ls_pair_BBC_Small_TDC_East");
  TTreeReaderArray<UShort_t> ls_pair_BBC_Small_TDC_West(ls_reader, "ls_pair_BBC_Small_TDC_West");

  Hists ls_h("ls");
  while (ls_reader.Next()) {
    if (*ls_pair_BBC_Veto_East != 0 || *ls_pair_BBC_Veto_West != 0) continue;
    if (*ls_pair_event_n_primary_tracks != 2) continue;
    if (std::abs(*ls_pair_event_vz) >= 100) continue;
    {
      Track t{*ls_kaon1_nhits, *ls_kaon1_nhitsfit, *ls_kaon1_nhitsDEdx, *ls_kaon1_nhitsposs,
              *ls_kaon1_charge, *ls_kaon1_pt, *ls_kaon1_eta, *ls_kaon1_azimuthal_angle,
              *ls_kaon1_dEdx, *ls_kaon1_nsigma_kaon, *ls_kaon1_nsigma_pion,
              *ls_kaon1_nsigma_electron, *ls_kaon1_nsigma_proton,
              {*ls_kaon1_dca_x, *ls_kaon1_dca_y, *ls_kaon1_dca_z},
              {*ls_kaon1_dca_global_x, *ls_kaon1_dca_global_y, *ls_kaon1_dca_global_z},
              *ls_kaon1_tof_time, *ls_kaon1_tof_pathlength, *ls_kaon1_topology_map_itpc,
              *ls_kaon1_topology_map0, *ls_kaon1_topology_map1};
      ls_h.FillTrack(t);
    }
    {
      Track t{*ls_kaon2_nhits, *ls_kaon2_nhitsfit, *ls_kaon2_nhitsDEdx, *ls_kaon2_nhitsposs,
              *ls_kaon2_charge, *ls_kaon2_pt, *ls_kaon2_eta, *ls_kaon2_azimuthal_angle,
              *ls_kaon2_dEdx, *ls_kaon2_nsigma_kaon, *ls_kaon2_nsigma_pion,
              *ls_kaon2_nsigma_electron, *ls_kaon2_nsigma_proton,
              {*ls_kaon2_dca_x, *ls_kaon2_dca_y, *ls_kaon2_dca_z},
              {*ls_kaon2_dca_global_x, *ls_kaon2_dca_global_y, *ls_kaon2_dca_global_z},
              *ls_kaon2_tof_time, *ls_kaon2_tof_pathlength, *ls_kaon2_topology_map_itpc,
              *ls_kaon2_topology_map0, *ls_kaon2_topology_map1};
      ls_h.FillTrack(t);
    }
    ls_h.FillPair(*ls_phi_pt, *ls_phi_rapidity, *ls_phi_mass, *ls_phi_x2kk, *ls_phi_x2pipi,
                *ls_phi_x2ee, *ls_phi_x2pp, *ls_kaon1_nsigma_kaon, *ls_kaon2_nsigma_kaon);
    ls_h.FillEvent(*ls_pair_event_vx, *ls_pair_event_vy, *ls_pair_event_vz, *ls_pair_event_nvtx,
                 *ls_pair_event_n_global_tracks, *ls_pair_TofMult, *ls_pair_vertex_n_tracks_used,
                 *ls_pair_vertex_ranking, *ls_pair_ZDC_UA_East, *ls_pair_ZDC_UA_West);
  }

  // ---- Plots -----------------------------------------------------------------
  for (const char* d : {"qa/track", "qa/track2d", "qa/pair", "qa/event"}) gSystem->mkdir(d, true);
  const TString T = "qa/track/track_", T2 = "qa/track2d/track_", P = "qa/pair/", E = "qa/event/";
  const TString kp = "p [GeV/c]", kde = "dE/dx [keV/cm]";

  // Track, 1D (ULS vs LS)
  Overlay(uls_h, ls_h, "nhits", "nHits", "Tracks", T + "nhits");
  Overlay(uls_h, ls_h, "nhitsfit", "nHitsFit", "Tracks", T + "nhitsfit");
  Overlay(uls_h, ls_h, "nhitsfit", "nHitsFit", "Tracks", T + "nhitsfit_log", true);
  Overlay(uls_h, ls_h, "nhitsdedx", "nHitsDedx", "Tracks", T + "nhitsdedx");
  Overlay(uls_h, ls_h, "nhitsposs", "nHitsPoss", "Tracks", T + "nhitsposs");
  Overlay(uls_h, ls_h, "itpc_rows", "Number of iTPC (inner) pad rows hit per track", "Tracks", T + "itpc_rows");
  Overlay(uls_h, ls_h, "outer_rows", "Number of outer pad rows hit per track", "Tracks", T + "outer_rows");
  Overlay(uls_h, ls_h, "row_occupancy", "Pad row (1-40 iTPC, 41-72 outer)", "Tracks with a hit in the row", T + "row_occupancy");
  Overlay(uls_h, ls_h, "nhitsfit_ratio", "nHitsFit / nHitsPoss", "Tracks", T + "nhitsfit_ratio");
  Overlay(uls_h, ls_h, "nhitsdedx_ratio", "nHitsDedx / nHitsFit", "Tracks", T + "nhitsdedx_ratio");
  Overlay(uls_h, ls_h, "charge", "Charge", "Tracks", T + "charge");
  Overlay(uls_h, ls_h, "pt", "p_{T} [GeV/c]", "Tracks", T + "pt");
  Overlay(uls_h, ls_h, "p", kp, "Tracks", T + "p");
  Overlay(uls_h, ls_h, "eta", "#eta", "Tracks", T + "eta");
  Overlay(uls_h, ls_h, "phi", "#varphi [rad]", "Tracks", T + "phi");
  Overlay(uls_h, ls_h, "dedx", kde, "Tracks", T + "dedx", true);
  Overlay(uls_h, ls_h, "nsigma_kaon", "n#sigma_{K}", "Tracks", T + "nsigma_kaon");
  Overlay(uls_h, ls_h, "nsigma_pion", "n#sigma_{#pi}", "Tracks", T + "nsigma_pion");
  Overlay(uls_h, ls_h, "nsigma_electron", "n#sigma_{e}", "Tracks", T + "nsigma_electron");
  Overlay(uls_h, ls_h, "nsigma_proton", "n#sigma_{p}", "Tracks", T + "nsigma_proton");
  Overlay(uls_h, ls_h, "dca_x", "Primary DCA_{x} [cm]", "Tracks", T + "dca_x", true);
  Overlay(uls_h, ls_h, "dca_y", "Primary DCA_{y} [cm]", "Tracks", T + "dca_y", true);
  Overlay(uls_h, ls_h, "dca_z", "Primary DCA_{z} [cm]", "Tracks", T + "dca_z", true);
  Overlay(uls_h, ls_h, "dca_xy", "Primary DCA_{xy} [cm]", "Tracks", T + "dca_xy", true);
  Overlay(uls_h, ls_h, "dca_3d", "Primary DCA [cm]", "Tracks", T + "dca_3d", true);
  Overlay(uls_h, ls_h, "gdca_x", "Global DCA_{x} [cm]", "Tracks", T + "gdca_x", true);
  Overlay(uls_h, ls_h, "gdca_y", "Global DCA_{y} [cm]", "Tracks", T + "gdca_y", true);
  Overlay(uls_h, ls_h, "gdca_z", "Global DCA_{z} [cm]", "Tracks", T + "gdca_z", true);
  Overlay(uls_h, ls_h, "gdca_xy", "Global DCA_{xy} [cm]", "Tracks", T + "gdca_xy", true);
  Overlay(uls_h, ls_h, "gdca_3d", "Global DCA [cm]", "Tracks", T + "gdca_3d", true);
  Overlay(uls_h, ls_h, "tof_matched", "TOF matched", "Tracks", T + "tof_matched");
  Overlay(uls_h, ls_h, "tof_invbeta", "1/#beta (TOF matched)", "Tracks", T + "tof_invbeta");

  // Track, 2D, all tracks (ULS + LS)
  Hists all_h("all");
  for (const char* n : {"p_dedx", "p_nsigma_kaon"}) {
    TH1* sum = static_cast<TH1*>(uls_h.h.at(n)->Clone(TString("all_") + n));
    sum->Add(ls_h.h.at(n));
    all_h.h[n] = sum;
  }
  all_h.h.at("p_dedx")->GetXaxis()->SetRangeUser(0, 1.6);
  all_h.h.at("p_dedx")->GetYaxis()->SetRangeUser(0, 36);
  all_h.h.at("p_nsigma_kaon")->GetXaxis()->SetRangeUser(0, 1.6);
  all_h.h.at("p_nsigma_kaon")->GetYaxis()->SetRangeUser(-12, 28);
  Plot2D(all_h, "p_dedx", kp, kde, T2 + "p_dedx_expected", 0.50, true, 0.55, 0.36);
  Plot2D(all_h, "p_nsigma_kaon", kp, "n#sigma_{K}", T2 + "p_nsigma_kaon", 0.50);

  // Track, 2D (ULS and LS separately)
  for (Hists* hs : {&uls_h, &ls_h}) {
    Plot2D(*hs, "p_dedx", kp, kde, T2 + "p_dedx", 0.45);
    Plot2D(*hs, "p_dedx", kp, kde, T2 + "p_dedx_expected", 0.45, true);
    Plot2D(*hs, "qp_dedx", "q #times p [GeV/c]", kde, T2 + "qp_dedx", 0.45);
    Plot2D(*hs, "p_nsigma_kaon", kp, "n#sigma_{K}", T2 + "p_nsigma_kaon", 0.45);
    Plot2D(*hs, "p_nsigma_pion", kp, "n#sigma_{#pi}", T2 + "p_nsigma_pion", 0.45);
    Plot2D(*hs, "p_invbeta", kp, "1/#beta (TOF matched)", T2 + "p_invbeta", 0.45);
    Plot2D(*hs, "eta_phi", "#eta", "#varphi [rad]", T2 + "eta_phi");
    Plot2D(*hs, "eta_pt", "#eta", "p_{T} [GeV/c]", T2 + "eta_pt");
    Plot2D(*hs, "eta_nhitsfit", "#eta", "nHitsFit", T2 + "eta_nhitsfit");
    Plot2D(*hs, "pt_nhitsfit", "p_{T} [GeV/c]", "nHitsFit", T2 + "pt_nhitsfit");
    Plot2D(*hs, "inner_outer_rows", "iTPC (inner) pad rows hit", "Outer pad rows hit", T2 + "inner_outer_rows", 0.45);
    Plot2D(*hs, "nhitsfit_nhitsdedx", "nHitsFit", "nHitsDedx", T2 + "nhitsfit_nhitsdedx");
    Plot2D(*hs, "nhitsfit_gdca_3d", "nHitsFit", "Global DCA [cm]", T2 + "nhitsfit_gdca_3d", 0.45);
    Plot2D(*hs, "eta_gdca_3d", "#eta", "Global DCA [cm]", T2 + "eta_gdca_3d");
    Plot2D(*hs, "gdca_xy_z", "Global DCA_{z} [cm]", "Global DCA_{xy} [cm]", T2 + "gdca_xy_z");
    Plot2D(*hs, "dca_xy_z", "Primary DCA_{z} [cm]", "Primary DCA_{xy} [cm]", T2 + "dca_xy_z");
  }

  // Pair
  Overlay(uls_h, ls_h, "pair_pt", "Pair p_{T} [GeV/c]", "Pairs", P + "pair_pt");
  Overlay(uls_h, ls_h, "pair_pt", "Pair p_{T} [GeV/c]", "Pairs", P + "pair_pt_log", true);
  PtRegions(uls_h, ls_h, P + "pair_pt_regions");
  Overlay(uls_h, ls_h, "pair_pt2", "Pair p_{T}^{2} [GeV^{2}/c^{2}]", "Pairs", P + "pair_pt2", true);
  Overlay(uls_h, ls_h, "pair_rapidity", "Pair rapidity (K^{+}K^{-} hypothesis)", "Pairs", P + "pair_rapidity");
  Overlay(uls_h, ls_h, "pair_mass", "M_{KK} [GeV/c^{2}]", "Pairs", P + "pair_mass");
  Overlay(uls_h, ls_h, "pair_x2kk", "#chi^{2}_{KK}", "Pairs", P + "pair_x2kk");
  Overlay(uls_h, ls_h, "pair_x2pipi", "#chi^{2}_{#pi#pi}", "Pairs", P + "pair_x2pipi", true);
  Overlay(uls_h, ls_h, "pair_x2ee", "#chi^{2}_{ee}", "Pairs", P + "pair_x2ee", true);
  Overlay(uls_h, ls_h, "pair_x2pp", "#chi^{2}_{pp}", "Pairs", P + "pair_x2pp", true);
  for (Hists* hs : {&uls_h, &ls_h}) {
    Plot2D(*hs, "pair_nsigma_kaon_12", "n#sigma_{K} (track 1)", "n#sigma_{K} (track 2)", P + "pair_nsigma_kaon_12");
    Plot2D(*hs, "pair_x2kk_x2pipi", "#chi^{2}_{KK}", "#chi^{2}_{#pi#pi}", P + "pair_x2kk_x2pipi", 0.45);
    Plot2D(*hs, "pair_mass_pt", "M_{KK} [GeV/c^{2}]", "Pair p_{T} [GeV/c]", P + "pair_mass_pt", 0.45);
  }

  // Event (one entry per pair)
  Overlay(uls_h, ls_h, "event_vz", "V_{z} [cm]", "Pairs", E + "event_vz");
  Overlay(uls_h, ls_h, "event_vr", "V_{r} [cm]", "Pairs", E + "event_vr", true);
  Overlay(uls_h, ls_h, "event_nvtx", "Number of vertices", "Pairs", E + "event_nvtx", true);
  Overlay(uls_h, ls_h, "event_n_global_tracks", "Number of global tracks", "Pairs", E + "event_n_global_tracks", true);
  Overlay(uls_h, ls_h, "event_TofMult", "TOF multiplicity", "Pairs", E + "event_TofMult", true);
  Overlay(uls_h, ls_h, "event_vertex_n_tracks_used", "Tracks used in vertex fit", "Pairs", E + "event_vertex_n_tracks_used", true);
  Overlay(uls_h, ls_h, "event_vertex_ranking", "Vertex ranking", "Pairs", E + "event_vertex_ranking");
  Overlay(uls_h, ls_h, "event_ZDC_UA_East", "ZDC East UA [ADC]", "Pairs", E + "event_ZDC_UA_East", true);
  Overlay(uls_h, ls_h, "event_ZDC_UA_West", "ZDC West UA [ADC]", "Pairs", E + "event_ZDC_UA_West", true);
  for (Hists* hs : {&uls_h, &ls_h}) Plot2D(*hs, "event_vx_vy", "V_{x} [cm]", "V_{y} [cm]", E + "event_vx_vy");
}
