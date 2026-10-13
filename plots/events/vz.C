#include "VzFit.h"
#include "TFile.h"
#include "TH2D.h"
#include "TStyle.h"
#include "TTree.h"
#include <memory>

// Tree preselection (stored in vz_data UserInfo): trigger 700001, bad runs removed,
// lastDSM(1) bits 1 and 2 both zero (BBC veto), exactly two primary tracks on one vertex.
// Coordinates are the common two-track vertex (vertex_*); vr is about (0,0).

void vz(TString input_file = "../../data/phi-kk-mudst-vz.root") {
  std::unique_ptr<TFile> file(TFile::Open(input_file, "READ"));
  TTree* tree = dynamic_cast<TTree*>(file->Get("vz_data"));
  gStyle->SetPalette(kBird);

  const TString valid = "vertex_position_valid"; // finite coordinates only, not a physics cut

  // vr before any cut: all two-track events, and split by pair charge.
  TH1D* vr_all = new TH1D("vr_all", "", 150, 0, 15);
  TH1D* vr_uls = new TH1D("vr_uls", "", 150, 0, 15);
  TH1D* vr_ls = new TH1D("vr_ls", "", 150, 0, 15);
  tree->Draw("vertex_vr>>vr_all", valid, "goff");
  tree->Draw("vertex_vr>>vr_uls", valid + "&&is_uls", "goff");
  tree->Draw("vertex_vr>>vr_ls", valid + "&&is_ls", "goff");

  PlotStyle::Options opt;
  opt.xlabel = "v_{r} [cm]";
  opt.ylabel = "Events / 0.1 cm";
  opt.logY = true;
  opt.title = "Run-19 Au+Au";
  opt.subtitle = "#sqrt{s_{NN}} = 200 GeV";
  opt.note = "ZDCMB, BBC veto, two-track vertex";
  TCanvas* c_vr = PlotStyle::Draw({
    {vr_all, "All"},
    {vr_uls, "ULS"},
    {vr_ls, "LS"}
  }, opt);
  if (!c_vr) return;
  c_vr->Print("vz_vr.png");
  c_vr->Print("vz_vr.pdf");

  // vz before and after vr < 3 cm, All / ULS / LS on one canvas each.
  const TString vr_cut[2] = {"", "&&vertex_vr<3"};
  const TString vz_note[2] = {"ZDCMB, BBC veto, two-track vertex", "ZDCMB, BBC veto, v_{r} < 3 cm"};
  const TString vz_name[2] = {"vz_vz", "vz_vz_vr3"};
  for (int i = 0; i < 2; ++i) {
    TH1D* vz_all = new TH1D(vz_name[i] + "_all", "", 100, -250, 250);
    TH1D* vz_uls = new TH1D(vz_name[i] + "_uls", "", 100, -250, 250);
    TH1D* vz_ls = new TH1D(vz_name[i] + "_ls", "", 100, -250, 250);
    tree->Draw("vertex_vz>>" + vz_name[i] + "_all", valid + vr_cut[i], "goff");
    tree->Draw("vertex_vz>>" + vz_name[i] + "_uls", valid + vr_cut[i] + "&&is_uls", "goff");
    tree->Draw("vertex_vz>>" + vz_name[i] + "_ls", valid + vr_cut[i] + "&&is_ls", "goff");

    PlotStyle::Options opt_vz = opt;
    opt_vz.xlabel = "v_{z} [cm]";
    opt_vz.ylabel = "Events / 5 cm";
    opt_vz.logY = false;
    opt_vz.note = vz_note[i];
    TCanvas* c_vz = PlotStyle::Draw({
      {vz_all, "All"},
      {vz_uls, "ULS"},
      {vz_ls, "LS"}
    }, opt_vz);
    if (!c_vz) return;
    c_vz->Print(vz_name[i] + ".png");
    c_vz->Print(vz_name[i] + ".pdf");
  }

  // Double-Gaussian vz fits: ULS only and ULS - LS, with and without vr < 3 cm.
  const TString fit_note[2] = {"ZDCMB, BBC veto, no v_{r} cut", "ZDCMB, BBC veto, v_{r} < 3 cm"};
  const TString fit_name[2] = {"vz_fit", "vz_fit_vr3"};
  for (int i = 0; i < 2; ++i) {
    TH1D* fit_uls = new TH1D(fit_name[i] + "_uls_h", "", 80, -200, 200);
    TH1D* fit_ls = new TH1D(fit_name[i] + "_ls_h", "", 80, -200, 200);
    tree->Draw("vertex_vz>>" + fit_name[i] + "_uls_h", valid + vr_cut[i] + "&&is_uls", "goff");
    tree->Draw("vertex_vz>>" + fit_name[i] + "_ls_h", valid + vr_cut[i] + "&&is_ls", "goff");

    PlotStyle::Options opt_fit = opt;
    opt_fit.note = fit_note[i];
    FitVzPair(fit_uls, fit_ls, opt_fit, fit_name[i]);
  }

  PlotStyle::Options opt2 = opt;
  opt2.xlabel = "v_{x} [cm]";
  opt2.ylabel = "v_{y} [cm]";
  opt2.showEntries = true;

  // x spans are 1.2x the y spans so circles stay round on the 1100x800 canvas.
  // vx:vy before any cut.
  // Selection and entries go top right, so the range only needs headroom above |v_y| = 10 cm.
  TH2D* xy_all = new TH2D("xy_all", "", 336, -16.8, 16.8, 280, -11, 17);
  tree->Draw("vertex_vy:vertex_vx>>xy_all", valid, "goff");
  PlotStyle::Options opt_all = opt2;
  opt_all.note = "UPC-like events";
  opt_all.showEntries = false;
  Draw2D(xy_all, opt_all);
  DrawLines({"ZDCMB, bad runs excluded", "BBC veto", "2 primary tracks on one vertex",
             "No PID or kinematic cut",
             ("Entries = " + PlotStyle::ShortNumber(xy_all->GetEntries())).c_str()},
            kDraw2DRight, opt.annotationY, 0.028, 33);
  gPad->Print("vz_vxvy.png");
  gPad->Print("vz_vxvy.pdf");

  // Ring region: 2 < vr < 6 cm, where the r ~ 3.8 cm ring sits.
  TH2D* xy_ring = new TH2D("xy_ring", "", 384, -9.6, 9.6, 320, -6, 10);
  tree->Draw("vertex_vy:vertex_vx>>xy_ring", valid + "&&vertex_vr>2&&vertex_vr<6", "goff");
  opt2.note = "ZDCMB, BBC veto, 2 < v_{r} < 6 cm";
  Draw2D(xy_ring, opt2)->Print("vz_vxvy_ring.png");
  gPad->Print("vz_vxvy_ring.pdf");

  // vx:vy after vr < 3 cm.
  TH2D* xy_vr3 = new TH2D("xy_vr3", "", 360, -6, 6, 300, -4, 6);
  tree->Draw("vertex_vy:vertex_vx>>xy_vr3", valid + "&&vertex_vr<3", "goff");
  opt2.note = "ZDCMB, BBC veto, v_{r} < 3 cm";
  Draw2D(xy_vr3, opt2)->Print("vz_vxvy_vr3.png");
  gPad->Print("vz_vxvy_vr3.pdf");

  // Same, ULS pairs only.
  TH2D* xy_vr3_uls = new TH2D("xy_vr3_uls", "", 360, -6, 6, 300, -4, 6);
  tree->Draw("vertex_vy:vertex_vx>>xy_vr3_uls", valid + "&&vertex_vr<3&&is_uls", "goff");
  opt2.note = "ZDCMB, BBC veto, v_{r} < 3 cm, ULS";
  Draw2D(xy_vr3_uls, opt2)->Print("vz_vxvy_vr3_uls.png");
  gPad->Print("vz_vxvy_vr3_uls.pdf");
}
