#include "TFile.h"
#include "TString.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TTreeReaderValue.h"
#include <iostream>
#include <memory>

// Pair-tree skimmer for phi-kk-mudst.root. Each tree entry is one K+K- candidate pair
// (uls_data: unlike-sign, ls_data: like-sign). Every branch of both trees is bound below.

void track_skim(TString input_file = "../../data/phi-kk-mudst.root") {
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

  while (uls_reader.Next()) {

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

  while (ls_reader.Next()) {

  }
}
