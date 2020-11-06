// Standard Root packages needed for the analysis.
#include <iostream>
#include <fstream>
#include <TFile.h>
#include <TChain.h>
#include <TTree.h>

using namespace std;

void Generating_Signal_And_Background_Trees_With_Background_Quantities() {

  TFile* output_file = new TFile("input_file_testing_numi_bkgd_overlay_with_all_weights_full_stats.root", "RECREATE");

  // Declare the variables that are output into the signal tree.
  int   run;
  int   subrun;
  int   event;
  float truth_neutrino_energy;
  int   NC_channel;
  int   num_muminus_tracks;
  int   num_muplus_tracks;
  int   num_piplus_tracks;
  int   num_piminus_tracks;
  int   num_pi0_tracks;
  int   num_proton_tracks;
  int   num_electron_tracks;
  int   num_positron_tracks;
  int   num_photon_tracks;
  float truth_reco_vtx_distance;
  float spline_fix_mcweight;
  float special_kdar_weight;
  float sampleweight;
  float treeweight;
  int   ext;
  float neutrino_energy;
  float neutrino_cos_angle;
  float reco_vtx_x;
  float reco_vtx_y;
  float reco_vtx_z;
  float muon_x_momentum_normalized;
  float muon_y_momentum_normalized;
  float muon_z_momentum_normalized;
  float x_component_of_muon_momentum;
  float y_component_of_muon_momentum;
  float z_component_of_muon_momentum;
  float num_beam_flashes;
  float flash_PEs;
  float flash_z;
  float flash_y;
  float flash_z_difference;
  float flash_y_difference;
  float truncated_dQdx_of_TPCObject_muon_candidate_u_plane;
  float truncated_dQdx_of_TPCObject_muon_candidate_v_plane;
  float truncated_dQdx_of_TPCObject_muon_candidate_y_plane;
  float truncated_dQdx_of_TPCObject_muon_candidate_sum;
  float median_dQdx_of_TPCObject_muon_candidate_u_plane;
  float median_dQdx_of_TPCObject_muon_candidate_v_plane;
  float median_dQdx_of_TPCObject_muon_candidate_y_plane;
  float median_dQdx_of_TPCObject_muon_candidate_sum;
  float u_plane_hit_sum;
  float v_plane_hit_sum;
  float y_plane_hit_sum;
  float hit_sum;
  float u_plane_vicinity_hit_sum;
  float v_plane_vicinity_hit_sum;
  float y_plane_vicinity_hit_sum;
  float vicinity_hit_sum;
  float u_plane_vicinity_hit_sum_difference;
  float v_plane_vicinity_hit_sum_difference;
  float y_plane_vicinity_hit_sum_difference;
  float vicinity_hit_sum_difference;
  float u_plane_number_of_slice_hits_in_vicinity_ADCs;
  float v_plane_number_of_slice_hits_in_vicinity_ADCs;
  float y_plane_number_of_slice_hits_in_vicinity_ADCs;
  float total_number_of_slice_hits_in_vicinity_ADCs;
  float num_of_top_pandora_crossings;
  float num_of_bottom_pandora_crossings;
  float num_of_front_pandora_crossings;
  float num_of_back_pandora_crossings;
  int   number_of_tracks_in_TPCObject;
  int   num_of_tracks_originating_from_vertex;
  float muon_candidate_length;
  float sum_of_TPCObject_track_lengths;
  float total_length_of_tracks_originating_from_vertex;
  float rootino_fix_mcweight;
  float central_value_mcweight;
  float other_universe_mcweights[100];
  float axialff_mcweights[2];
  float rpaccqe_mcweights[2];
  float xsecshape_mcweights[2];
 
  TTree* signal_tree = new TTree("signal_tree", "");

  signal_tree->Branch("run", &run, "run/I");
  signal_tree->Branch("subrun", &subrun, "subrun/I");
  signal_tree->Branch("event", &event, "event/I");
  signal_tree->Branch("truth_neutrino_energy", &truth_neutrino_energy, "truth_neutrino_energy/F");
  signal_tree->Branch("NC_channel", &NC_channel, "NC_channel/I");
  signal_tree->Branch("truth_reco_vtx_distance", &truth_reco_vtx_distance, "truth_reco_vtx_distance/F");
  signal_tree->Branch("num_muminus_tracks", &num_muminus_tracks, "num_muminus_tracks/I");
  signal_tree->Branch("num_muplus_tracks", &num_muplus_tracks, "num_muplus_tracks/I");
  signal_tree->Branch("num_piplus_tracks", &num_piplus_tracks, "num_piplus_tracks/I");
  signal_tree->Branch("num_piminus_tracks", &num_piminus_tracks, "num_piminus_tracks/I");
  signal_tree->Branch("num_pi0_tracks", &num_pi0_tracks, "num_pi0_tracks/I");
  signal_tree->Branch("num_proton_tracks", &num_proton_tracks, "num_proton_tracks/I");
  signal_tree->Branch("num_electron_tracks", &num_electron_tracks, "num_electron_tracks/I");
  signal_tree->Branch("num_positron_tracks", &num_positron_tracks, "num_positron_tracks/I");
  signal_tree->Branch("num_photon_tracks", &num_photon_tracks, "num_photon_tracks/I");
  signal_tree->Branch("spline_fix_mcweight", &spline_fix_mcweight, "spline_fix_mcweight/F");
  signal_tree->Branch("neutrino_energy", &neutrino_energy, "neutrino_energy/F");
  signal_tree->Branch("neutrino_cos_angle", &neutrino_cos_angle, "neutrino_cos_angle/F");
  signal_tree->Branch("reco_vtx_x", &reco_vtx_x, "reco_vtx_x/F");
  signal_tree->Branch("reco_vtx_y", &reco_vtx_y, "reco_vtx_y/F");
  signal_tree->Branch("reco_vtx_z", &reco_vtx_z, "reco_vtx_z/F");
  signal_tree->Branch("muon_x_momentum_normalized", &muon_x_momentum_normalized, "muon_x_momentum_normalized/F");
  signal_tree->Branch("muon_y_momentum_normalized", &muon_y_momentum_normalized, "muon_y_momentum_normalized/F");
  signal_tree->Branch("muon_z_momentum_normalized", &muon_z_momentum_normalized, "muon_z_momentum_normalized/F");
  signal_tree->Branch("x_component_of_muon_momentum", &x_component_of_muon_momentum, "x_component_of_muon_momentum/F");
  signal_tree->Branch("y_component_of_muon_momentum", &y_component_of_muon_momentum, "y_component_of_muon_momentum/F");
  signal_tree->Branch("z_component_of_muon_momentum", &z_component_of_muon_momentum, "z_component_of_muon_momentum/F");
  signal_tree->Branch("num_beam_flashes", &num_beam_flashes, "num_beam_flashes/F");
  signal_tree->Branch("flash_PEs", &flash_PEs, "flash_PEs/F");
  signal_tree->Branch("flash_z_difference", &flash_z_difference, "flash_z_difference/F");
  signal_tree->Branch("flash_y_difference", &flash_y_difference, "flash_y_difference/F");
  signal_tree->Branch("truncated_dQdx_of_TPCObject_muon_candidate_u_plane", &truncated_dQdx_of_TPCObject_muon_candidate_u_plane, "truncated_dQdx_of_TPCObject_muon_candidate_u_plane/F");
  signal_tree->Branch("truncated_dQdx_of_TPCObject_muon_candidate_v_plane", &truncated_dQdx_of_TPCObject_muon_candidate_v_plane, "truncated_dQdx_of_TPCObject_muon_candidate_v_plane/F");
  signal_tree->Branch("truncated_dQdx_of_TPCObject_muon_candidate_y_plane", &truncated_dQdx_of_TPCObject_muon_candidate_y_plane, "truncated_dQdx_of_TPCObject_muon_candidate_y_plane/F");
  signal_tree->Branch("truncated_dQdx_of_TPCObject_muon_candidate_sum", &truncated_dQdx_of_TPCObject_muon_candidate_sum, "truncated_dQdx_of_TPCObject_muon_candidate_sum/F");
  signal_tree->Branch("median_dQdx_of_TPCObject_muon_candidate_u_plane", &median_dQdx_of_TPCObject_muon_candidate_u_plane, "median_dQdx_of_TPCObject_muon_candidate_u_plane/F");
  signal_tree->Branch("median_dQdx_of_TPCObject_muon_candidate_v_plane", &median_dQdx_of_TPCObject_muon_candidate_v_plane, "median_dQdx_of_TPCObject_muon_candidate_v_plane/F");
  signal_tree->Branch("median_dQdx_of_TPCObject_muon_candidate_y_plane", &median_dQdx_of_TPCObject_muon_candidate_y_plane, "median_dQdx_of_TPCObject_muon_candidate_y_plane/F");
  signal_tree->Branch("median_dQdx_of_TPCObject_muon_candidate_sum", &median_dQdx_of_TPCObject_muon_candidate_sum, "median_dQdx_of_TPCObject_muon_candidate_sum/F");
  signal_tree->Branch("u_plane_hit_sum", &u_plane_hit_sum, "u_plane_hit_sum/F");
  signal_tree->Branch("v_plane_hit_sum", &v_plane_hit_sum, "v_plane_hit_sum/F");
  signal_tree->Branch("y_plane_hit_sum", &y_plane_hit_sum, "y_plane_hit_sum/F");
  signal_tree->Branch("hit_sum", &hit_sum, "hit_sum/F");
  signal_tree->Branch("u_plane_vicinity_hit_sum", &u_plane_vicinity_hit_sum, "u_plane_vicinity_hit_sum/F");
  signal_tree->Branch("v_plane_vicinity_hit_sum", &v_plane_vicinity_hit_sum, "v_plane_vicinity_hit_sum/F");
  signal_tree->Branch("y_plane_vicinity_hit_sum", &y_plane_vicinity_hit_sum, "y_plane_vicinity_hit_sum/F");
  signal_tree->Branch("vicinity_hit_sum", &vicinity_hit_sum, "vicinity_hit_sum/F");
  signal_tree->Branch("u_plane_vicinity_hit_sum_difference", &u_plane_vicinity_hit_sum_difference, "u_plane_vicinity_hit_sum_difference/F");
  signal_tree->Branch("v_plane_vicinity_hit_sum_difference", &v_plane_vicinity_hit_sum_difference, "v_plane_vicinity_hit_sum_difference/F");
  signal_tree->Branch("y_plane_vicinity_hit_sum_difference", &y_plane_vicinity_hit_sum_difference, "y_plane_vicinity_hit_sum_difference/F");
  signal_tree->Branch("vicinity_hit_sum_difference", &vicinity_hit_sum_difference, "vicinity_hit_sum_difference/F");
  signal_tree->Branch("u_plane_number_of_slice_hits_in_vicinity_ADCs", &u_plane_number_of_slice_hits_in_vicinity_ADCs, "u_plane_number_of_slice_hits_in_vicinity_ADCs/F");
  signal_tree->Branch("v_plane_number_of_slice_hits_in_vicinity_ADCs", &v_plane_number_of_slice_hits_in_vicinity_ADCs, "v_plane_number_of_slice_hits_in_vicinity_ADCs/F");
  signal_tree->Branch("y_plane_number_of_slice_hits_in_vicinity_ADCs", &y_plane_number_of_slice_hits_in_vicinity_ADCs, "y_plane_number_of_slice_hits_in_vicinity_ADCs/F");
  signal_tree->Branch("total_number_of_slice_hits_in_vicinity_ADCs", &total_number_of_slice_hits_in_vicinity_ADCs, "total_number_of_slice_hits_in_vicinity_ADCs/F");
  signal_tree->Branch("num_of_top_pandora_crossings", &num_of_top_pandora_crossings, "num_of_top_pandora_crossings/F");
  signal_tree->Branch("num_of_bottom_pandora_crossings", &num_of_bottom_pandora_crossings, "num_of_bottom_pandora_crossings/F");
  signal_tree->Branch("number_of_tracks_in_TPCObject", &number_of_tracks_in_TPCObject, "number_of_tracks_in_TPCObject/I");
  signal_tree->Branch("num_of_tracks_originating_from_vertex", &num_of_tracks_originating_from_vertex, "num_of_tracks_originating_from_vertex/I");
  signal_tree->Branch("muon_candidate_length", &muon_candidate_length, "muon_candidate_length/F");
  signal_tree->Branch("sum_of_TPCObject_track_lengths", &sum_of_TPCObject_track_lengths, "sum_of_TPCObject_track_lengths/F");
  signal_tree->Branch("total_length_of_tracks_originating_from_vertex", &total_length_of_tracks_originating_from_vertex, "total_length_of_tracks_originating_from_vertex/F");

  // Declare the variables that are input from the KDAR events.
  int    KDAR_input_run;
  int    KDAR_input_subrun;
  int    KDAR_input_event;
  int    KDAR_input_NC_Channel;
  float  KDAR_input_truth_neutrino_energy;
  int    KDAR_input_num_muminus_tracks;
  int    KDAR_input_num_muplus_tracks;
  int    KDAR_input_num_piplus_tracks;
  int    KDAR_input_num_piminus_tracks;
  int    KDAR_input_num_pi0_tracks;
  int    KDAR_input_num_proton_tracks;
  int    KDAR_input_num_electron_showers;
  int    KDAR_input_num_positron_showers;
  int    KDAR_input_num_photon_showers;
  float  KDAR_input_spline_fix_mcweight;
  float  KDAR_input_truth_reco_vtx_distance;
  float  KDAR_input_neutrino_energy;
  float  KDAR_input_neutrino_cos_angle;
  float  KDAR_input_reco_vtx_x;
  float  KDAR_input_reco_vtx_y;
  float  KDAR_input_reco_vtx_z;
  float  KDAR_input_muon_x_momentum_normalized;
  float  KDAR_input_muon_y_momentum_normalized;
  float  KDAR_input_muon_z_momentum_normalized;
  float  KDAR_input_x_component_of_muon_momentum;
  float  KDAR_input_y_component_of_muon_momentum;
  float  KDAR_input_z_component_of_muon_momentum;
  float  KDAR_input_num_beam_flashes;
  float  KDAR_input_flash_PEs;
  float  KDAR_input_flash_z;
  float  KDAR_input_flash_y;
  float  KDAR_input_truncated_dQdx_of_TPCObject_muon_candidate_u_plane;
  float  KDAR_input_truncated_dQdx_of_TPCObject_muon_candidate_v_plane;
  float  KDAR_input_truncated_dQdx_of_TPCObject_muon_candidate_y_plane;
  float  KDAR_input_truncated_dQdx_of_TPCObject_muon_candidate_sum;
  float  KDAR_input_median_dQdx_of_TPCObject_muon_candidate_u_plane;
  float  KDAR_input_median_dQdx_of_TPCObject_muon_candidate_v_plane;
  float  KDAR_input_median_dQdx_of_TPCObject_muon_candidate_y_plane;
  float  KDAR_input_median_dQdx_of_TPCObject_muon_candidate_sum;
  float  KDAR_input_u_plane_hit_sum;
  float  KDAR_input_v_plane_hit_sum;
  float  KDAR_input_y_plane_hit_sum;
  float  KDAR_input_hit_sum;
  float  KDAR_input_u_plane_vicinity_hit_sum;
  float  KDAR_input_v_plane_vicinity_hit_sum;
  float  KDAR_input_y_plane_vicinity_hit_sum;
  float  KDAR_input_vicinity_hit_sum;
  float  KDAR_input_u_plane_number_of_slice_hits_in_vicinity_ADCs;
  float  KDAR_input_v_plane_number_of_slice_hits_in_vicinity_ADCs;
  float  KDAR_input_y_plane_number_of_slice_hits_in_vicinity_ADCs;
  float  KDAR_input_total_number_of_slice_hits_in_vicinity_ADCs;
  float  KDAR_input_num_of_top_pandora_crossings;
  float  KDAR_input_num_of_bottom_pandora_crossings;
  float  KDAR_input_num_of_front_pandora_crossings;
  float  KDAR_input_num_of_back_pandora_crossings;
  int    KDAR_input_number_of_tracks_in_TPCObject;
  int    KDAR_input_num_of_tracks_originating_from_vertex;
  float  KDAR_input_muon_candidate_length;
  float  KDAR_input_sum_of_TPCObject_track_lengths;
  float  KDAR_input_total_length_of_tracks_originating_from_vertex;

  TChain* KDAR_input_tree = new TChain("t0ana/tree");
  KDAR_input_tree->Add("/home/barnchri/Preparing_To_Perform_BDT_Studies/NuMI_KDAR_Box_Opening_Testing_File_All_Cuts.root");
  KDAR_input_tree->SetBranchAddress("run", &KDAR_input_run);
  KDAR_input_tree->SetBranchAddress("subrun", &KDAR_input_subrun);
  KDAR_input_tree->SetBranchAddress("event", &KDAR_input_event);
  KDAR_input_tree->SetBranchAddress("NC_channel", &KDAR_input_NC_Channel);
  KDAR_input_tree->SetBranchAddress("truth_energy_in_tree", &KDAR_input_truth_neutrino_energy);
  KDAR_input_tree->SetBranchAddress("truth_reco_vtx_distance", &KDAR_input_truth_reco_vtx_distance);
  KDAR_input_tree->SetBranchAddress("num_muminus_tracks", &KDAR_input_num_muminus_tracks);
  KDAR_input_tree->SetBranchAddress("num_muplus_tracks", &KDAR_input_num_muplus_tracks);
  KDAR_input_tree->SetBranchAddress("num_piplus_tracks", &KDAR_input_num_piplus_tracks);
  KDAR_input_tree->SetBranchAddress("num_piminus_tracks", &KDAR_input_num_piminus_tracks);
  KDAR_input_tree->SetBranchAddress("num_pi0_tracks", &KDAR_input_num_pi0_tracks);
  KDAR_input_tree->SetBranchAddress("num_proton_tracks", &KDAR_input_num_proton_tracks);
  KDAR_input_tree->SetBranchAddress("num_electron_showers", &KDAR_input_num_electron_showers);
  KDAR_input_tree->SetBranchAddress("num_positron_showers", &KDAR_input_num_positron_showers);
  KDAR_input_tree->SetBranchAddress("num_photon_showers", &KDAR_input_num_photon_showers);
  KDAR_input_tree->SetBranchAddress("spline_fix_mcweight", &KDAR_input_spline_fix_mcweight);
  KDAR_input_tree->SetBranchAddress("reco_energy_in_tree", &KDAR_input_neutrino_energy);
  KDAR_input_tree->SetBranchAddress("cos_angle_between_truth_and_reco_direction_vectors_in_tree", &KDAR_input_neutrino_cos_angle);
  KDAR_input_tree->SetBranchAddress("reco_vtx_x", &KDAR_input_reco_vtx_x);
  KDAR_input_tree->SetBranchAddress("reco_vtx_y", &KDAR_input_reco_vtx_y);
  KDAR_input_tree->SetBranchAddress("reco_vtx_z", &KDAR_input_reco_vtx_z);
  KDAR_input_tree->SetBranchAddress("muon_x_momentum_normalized", &KDAR_input_muon_x_momentum_normalized);
  KDAR_input_tree->SetBranchAddress("muon_y_momentum_normalized", &KDAR_input_muon_y_momentum_normalized);
  KDAR_input_tree->SetBranchAddress("muon_z_momentum_normalized", &KDAR_input_muon_z_momentum_normalized);
  KDAR_input_tree->SetBranchAddress("x_component_of_muon_momentum", &KDAR_input_x_component_of_muon_momentum);
  KDAR_input_tree->SetBranchAddress("y_component_of_muon_momentum", &KDAR_input_y_component_of_muon_momentum);
  KDAR_input_tree->SetBranchAddress("z_component_of_muon_momentum", &KDAR_input_z_component_of_muon_momentum);
  KDAR_input_tree->SetBranchAddress("num_beam_flashes", &KDAR_input_num_beam_flashes);
  KDAR_input_tree->SetBranchAddress("flash_PEs", &KDAR_input_flash_PEs);
  KDAR_input_tree->SetBranchAddress("flash_z", &KDAR_input_flash_z);
  KDAR_input_tree->SetBranchAddress("flash_y", &KDAR_input_flash_y);
  KDAR_input_tree->SetBranchAddress("truncated_dQdx_of_TPCObject_muon_candidate_u_plane", &KDAR_input_truncated_dQdx_of_TPCObject_muon_candidate_u_plane);
  KDAR_input_tree->SetBranchAddress("truncated_dQdx_of_TPCObject_muon_candidate_v_plane", &KDAR_input_truncated_dQdx_of_TPCObject_muon_candidate_v_plane);
  KDAR_input_tree->SetBranchAddress("truncated_dQdx_of_TPCObject_muon_candidate_y_plane", &KDAR_input_truncated_dQdx_of_TPCObject_muon_candidate_y_plane);
  KDAR_input_tree->SetBranchAddress("truncated_dQdx_of_TPCObject_muon_candidate_sum", &KDAR_input_truncated_dQdx_of_TPCObject_muon_candidate_sum);
  KDAR_input_tree->SetBranchAddress("median_dQdx_of_TPCObject_muon_candidate_u_plane", &KDAR_input_median_dQdx_of_TPCObject_muon_candidate_u_plane);
  KDAR_input_tree->SetBranchAddress("median_dQdx_of_TPCObject_muon_candidate_v_plane", &KDAR_input_median_dQdx_of_TPCObject_muon_candidate_v_plane);
  KDAR_input_tree->SetBranchAddress("median_dQdx_of_TPCObject_muon_candidate_y_plane", &KDAR_input_median_dQdx_of_TPCObject_muon_candidate_y_plane);
  KDAR_input_tree->SetBranchAddress("median_dQdx_of_TPCObject_muon_candidate_sum", &KDAR_input_median_dQdx_of_TPCObject_muon_candidate_sum);
  KDAR_input_tree->SetBranchAddress("ADC_hit_sum_u_plane", &KDAR_input_u_plane_hit_sum);
  KDAR_input_tree->SetBranchAddress("ADC_hit_sum_v_plane", &KDAR_input_v_plane_hit_sum);
  KDAR_input_tree->SetBranchAddress("ADC_hit_sum_y_plane", &KDAR_input_y_plane_hit_sum);
  KDAR_input_tree->SetBranchAddress("ADC_hit_sum_total", &KDAR_input_hit_sum);
  KDAR_input_tree->SetBranchAddress("ADC_hits_vicinity_sum_u_plane", &KDAR_input_u_plane_vicinity_hit_sum);
  KDAR_input_tree->SetBranchAddress("ADC_hits_vicinity_sum_v_plane", &KDAR_input_v_plane_vicinity_hit_sum);
  KDAR_input_tree->SetBranchAddress("ADC_hits_vicinity_sum_y_plane", &KDAR_input_y_plane_vicinity_hit_sum);
  KDAR_input_tree->SetBranchAddress("ADC_hits_vicinity_sum_total", &KDAR_input_vicinity_hit_sum);
  KDAR_input_tree->SetBranchAddress("hits_vicinity_num_u_plane", &KDAR_input_u_plane_number_of_slice_hits_in_vicinity_ADCs);
  KDAR_input_tree->SetBranchAddress("hits_vicinity_num_v_plane", &KDAR_input_v_plane_number_of_slice_hits_in_vicinity_ADCs);
  KDAR_input_tree->SetBranchAddress("hits_vicinity_num_y_plane", &KDAR_input_y_plane_number_of_slice_hits_in_vicinity_ADCs);
  KDAR_input_tree->SetBranchAddress("hits_vicinity_num_total", &KDAR_input_total_number_of_slice_hits_in_vicinity_ADCs);
  KDAR_input_tree->SetBranchAddress("num_of_top_pandora_crossings", &KDAR_input_num_of_top_pandora_crossings);
  KDAR_input_tree->SetBranchAddress("num_of_bottom_pandora_crossings", &KDAR_input_num_of_bottom_pandora_crossings);
  KDAR_input_tree->SetBranchAddress("num_of_front_pandora_crossings", &KDAR_input_num_of_front_pandora_crossings);
  KDAR_input_tree->SetBranchAddress("num_of_back_pandora_crossings", &KDAR_input_num_of_back_pandora_crossings);
  KDAR_input_tree->SetBranchAddress("number_of_tracks_in_TPCObject", &KDAR_input_number_of_tracks_in_TPCObject);
  KDAR_input_tree->SetBranchAddress("num_of_tracks_originating_from_vertex", &KDAR_input_num_of_tracks_originating_from_vertex);
  KDAR_input_tree->SetBranchAddress("muon_candidate_length", &KDAR_input_muon_candidate_length);
  KDAR_input_tree->SetBranchAddress("sum_of_TPCObject_track_lengths", &KDAR_input_sum_of_TPCObject_track_lengths);
  KDAR_input_tree->SetBranchAddress("total_length_of_tracks_originating_from_vertex", &KDAR_input_total_length_of_tracks_originating_from_vertex);

  //int KDAR_input_entries = KDAR_input_tree->GetEntries();
  int KDAR_input_entries = 0;

  std::cout << "The number of entries in the KDAR input tree = " << KDAR_input_entries << "." << std::endl;

  for ( size_t iter = 0; iter < KDAR_input_entries; iter++ ) { // So it doesn't run over any events and we don't fill the final tree twice.

    // Set the tree to the current entry.
    KDAR_input_tree->GetEntry( iter );

    if ( KDAR_input_flash_PEs < 50.0 || KDAR_input_flash_PEs > 2000 || KDAR_input_sum_of_TPCObject_track_lengths > 65.0 ) continue;

    // Set all of the variables and fill the new tree.
    run                                                = KDAR_input_run;
    subrun                                             = KDAR_input_subrun;
    event                                              = KDAR_input_event;
    NC_channel                                         = KDAR_input_NC_Channel;
    truth_neutrino_energy                              = KDAR_input_truth_neutrino_energy;
    num_muminus_tracks                                 = KDAR_input_num_muminus_tracks;
    num_muplus_tracks                                  = KDAR_input_num_muplus_tracks;
    num_piplus_tracks                                  = KDAR_input_num_piplus_tracks;
    num_piminus_tracks                                 = KDAR_input_num_piminus_tracks;
    num_pi0_tracks                                     = KDAR_input_num_pi0_tracks;
    num_proton_tracks                                  = KDAR_input_num_proton_tracks;
    num_electron_tracks                                = KDAR_input_num_electron_showers;
    num_positron_tracks                                = KDAR_input_num_positron_showers;
    num_photon_tracks                                  = KDAR_input_num_photon_showers;
    spline_fix_mcweight                                = KDAR_input_spline_fix_mcweight;
    truth_reco_vtx_distance                            = KDAR_input_truth_reco_vtx_distance;
    neutrino_energy                                    = KDAR_input_neutrino_energy;
    neutrino_cos_angle                                 = KDAR_input_neutrino_cos_angle;
    reco_vtx_x                                         = KDAR_input_reco_vtx_x;
    reco_vtx_y                                         = KDAR_input_reco_vtx_y;
    reco_vtx_z                                         = KDAR_input_reco_vtx_z;
    muon_x_momentum_normalized                         = KDAR_input_muon_x_momentum_normalized;
    muon_y_momentum_normalized                         = KDAR_input_muon_y_momentum_normalized;
    muon_z_momentum_normalized                         = KDAR_input_muon_z_momentum_normalized;
    x_component_of_muon_momentum                       = KDAR_input_x_component_of_muon_momentum;
    y_component_of_muon_momentum                       = KDAR_input_y_component_of_muon_momentum;
    z_component_of_muon_momentum                       = KDAR_input_z_component_of_muon_momentum;
    num_beam_flashes                                   = KDAR_input_num_beam_flashes;
    flash_PEs                                          = KDAR_input_flash_PEs;
    flash_z                                            = KDAR_input_flash_z;
    flash_y                                            = KDAR_input_flash_y;
    truncated_dQdx_of_TPCObject_muon_candidate_u_plane = KDAR_input_truncated_dQdx_of_TPCObject_muon_candidate_u_plane;
    truncated_dQdx_of_TPCObject_muon_candidate_v_plane = KDAR_input_truncated_dQdx_of_TPCObject_muon_candidate_v_plane;
    truncated_dQdx_of_TPCObject_muon_candidate_y_plane = KDAR_input_truncated_dQdx_of_TPCObject_muon_candidate_y_plane;
    truncated_dQdx_of_TPCObject_muon_candidate_sum     = KDAR_input_truncated_dQdx_of_TPCObject_muon_candidate_sum;
    median_dQdx_of_TPCObject_muon_candidate_u_plane    = KDAR_input_median_dQdx_of_TPCObject_muon_candidate_u_plane;
    median_dQdx_of_TPCObject_muon_candidate_v_plane    = KDAR_input_median_dQdx_of_TPCObject_muon_candidate_v_plane;
    median_dQdx_of_TPCObject_muon_candidate_y_plane    = KDAR_input_median_dQdx_of_TPCObject_muon_candidate_y_plane;
    median_dQdx_of_TPCObject_muon_candidate_sum        = KDAR_input_median_dQdx_of_TPCObject_muon_candidate_sum;
    u_plane_hit_sum                                    = KDAR_input_u_plane_hit_sum;
    v_plane_hit_sum                                    = KDAR_input_v_plane_hit_sum;
    y_plane_hit_sum                                    = KDAR_input_y_plane_hit_sum;
    hit_sum                                            = KDAR_input_hit_sum;
    u_plane_vicinity_hit_sum                           = KDAR_input_u_plane_vicinity_hit_sum;
    v_plane_vicinity_hit_sum                           = KDAR_input_v_plane_vicinity_hit_sum;
    y_plane_vicinity_hit_sum                           = KDAR_input_y_plane_vicinity_hit_sum;
    vicinity_hit_sum                                   = KDAR_input_vicinity_hit_sum;
    u_plane_number_of_slice_hits_in_vicinity_ADCs      = KDAR_input_u_plane_number_of_slice_hits_in_vicinity_ADCs;
    v_plane_number_of_slice_hits_in_vicinity_ADCs      = KDAR_input_v_plane_number_of_slice_hits_in_vicinity_ADCs;
    y_plane_number_of_slice_hits_in_vicinity_ADCs      = KDAR_input_y_plane_number_of_slice_hits_in_vicinity_ADCs;
    total_number_of_slice_hits_in_vicinity_ADCs        = KDAR_input_total_number_of_slice_hits_in_vicinity_ADCs;
    num_of_top_pandora_crossings                       = KDAR_input_num_of_top_pandora_crossings;
    num_of_bottom_pandora_crossings                    = KDAR_input_num_of_bottom_pandora_crossings;
    num_of_front_pandora_crossings                     = KDAR_input_num_of_front_pandora_crossings;
    num_of_back_pandora_crossings                      = KDAR_input_num_of_back_pandora_crossings;
    number_of_tracks_in_TPCObject                      = KDAR_input_number_of_tracks_in_TPCObject;
    num_of_tracks_originating_from_vertex              = KDAR_input_num_of_tracks_originating_from_vertex;
    muon_candidate_length                              = KDAR_input_muon_candidate_length;
    sum_of_TPCObject_track_lengths                     = KDAR_input_sum_of_TPCObject_track_lengths;
    total_length_of_tracks_originating_from_vertex     = KDAR_input_total_length_of_tracks_originating_from_vertex;

    flash_z_difference                                 = ( reco_vtx_z - KDAR_input_flash_z );
    flash_y_difference                                 = ( reco_vtx_y - KDAR_input_flash_y );
    vicinity_hit_sum_difference                        = ( vicinity_hit_sum - hit_sum );
    u_plane_vicinity_hit_sum_difference                = ( u_plane_vicinity_hit_sum - u_plane_hit_sum );
    v_plane_vicinity_hit_sum_difference                = ( v_plane_vicinity_hit_sum - v_plane_hit_sum );
    y_plane_vicinity_hit_sum_difference                = ( y_plane_vicinity_hit_sum - y_plane_hit_sum );

    // Reset 'mcweight'.
    spline_fix_mcweight                                = 5.;

    signal_tree->Fill();

  }

  TTree* background_tree = new TTree("background_tree", "");

  background_tree->Branch("run", &run, "run/I");
  background_tree->Branch("subrun", &subrun, "subrun/I");
  background_tree->Branch("event", &event, "event/I");
  background_tree->Branch("truth_neutrino_energy", &truth_neutrino_energy, "truth_neutrino_energy/F");
  background_tree->Branch("truth_reco_vtx_distance", &truth_reco_vtx_distance, "truth_reco_vtx_distance/F");
  background_tree->Branch("NC_channel", &NC_channel, "NC_channel/I");
  background_tree->Branch("num_muminus_tracks", &num_muminus_tracks, "num_muminus_tracks/I");
  background_tree->Branch("num_muplus_tracks", &num_muplus_tracks, "num_muplus_tracks/I");
  background_tree->Branch("num_piplus_tracks", &num_piplus_tracks, "num_piplus_tracks/I");
  background_tree->Branch("num_piminus_tracks", &num_piminus_tracks, "num_piminus_tracks/I");
  background_tree->Branch("num_pi0_tracks", &num_pi0_tracks, "num_pi0_tracks/I");
  background_tree->Branch("num_proton_tracks", &num_proton_tracks, "num_proton_tracks/I");
  background_tree->Branch("num_electron_tracks", &num_electron_tracks, "num_electron_tracks/I");
  background_tree->Branch("num_positron_tracks", &num_positron_tracks, "num_positron_tracks/I");
  background_tree->Branch("num_photon_tracks", &num_photon_tracks, "num_photon_tracks/I");
  background_tree->Branch("special_kdar_weight", &special_kdar_weight, "special_kdar_weight/F");
  background_tree->Branch("spline_fix_mcweight", &spline_fix_mcweight, "spline_fix_mcweight/F");
  background_tree->Branch("sampleweight", &sampleweight, "sampleweight/F");
  background_tree->Branch("treeweight", &treeweight, "treeweight/F");
  background_tree->Branch("ext", &ext, "ext/I");
  background_tree->Branch("neutrino_energy", &neutrino_energy, "neutrino_energy/F");
  background_tree->Branch("neutrino_cos_angle", &neutrino_cos_angle, "neutrino_cos_angle/F");
  background_tree->Branch("reco_vtx_x", &reco_vtx_x, "reco_vtx_x/F");
  background_tree->Branch("reco_vtx_y", &reco_vtx_y, "reco_vtx_y/F");
  background_tree->Branch("reco_vtx_z", &reco_vtx_z, "reco_vtx_z/F");
  background_tree->Branch("muon_x_momentum_normalized", &muon_x_momentum_normalized, "muon_x_momentum_normalized/F");
  background_tree->Branch("muon_y_momentum_normalized", &muon_y_momentum_normalized, "muon_y_momentum_normalized/F");
  background_tree->Branch("muon_z_momentum_normalized", &muon_z_momentum_normalized, "muon_z_momentum_normalized/F");
  background_tree->Branch("x_component_of_muon_momentum", &x_component_of_muon_momentum, "x_component_of_muon_momentum/F");
  background_tree->Branch("y_component_of_muon_momentum", &y_component_of_muon_momentum, "y_component_of_muon_momentum/F");
  background_tree->Branch("z_component_of_muon_momentum", &z_component_of_muon_momentum, "z_component_of_muon_momentum/F");
  background_tree->Branch("num_beam_flashes", &num_beam_flashes, "num_beam_flashes/F");
  background_tree->Branch("flash_PEs", &flash_PEs, "flash_PEs/F");
  background_tree->Branch("flash_z_difference", &flash_z_difference, "flash_z_difference/F");
  background_tree->Branch("flash_y_difference", &flash_y_difference, "flash_y_difference/F");
  background_tree->Branch("truncated_dQdx_of_TPCObject_muon_candidate_u_plane", &truncated_dQdx_of_TPCObject_muon_candidate_u_plane, "truncated_dQdx_of_TPCObject_muon_candidate_u_plane/F");
  background_tree->Branch("truncated_dQdx_of_TPCObject_muon_candidate_v_plane", &truncated_dQdx_of_TPCObject_muon_candidate_v_plane, "truncated_dQdx_of_TPCObject_muon_candidate_v_plane/F");
  background_tree->Branch("truncated_dQdx_of_TPCObject_muon_candidate_y_plane", &truncated_dQdx_of_TPCObject_muon_candidate_y_plane, "truncated_dQdx_of_TPCObject_muon_candidate_y_plane/F");
  background_tree->Branch("truncated_dQdx_of_TPCObject_muon_candidate_sum", &truncated_dQdx_of_TPCObject_muon_candidate_sum, "truncated_dQdx_of_TPCObject_muon_candidate_sum/F");
  background_tree->Branch("median_dQdx_of_TPCObject_muon_candidate_u_plane", &median_dQdx_of_TPCObject_muon_candidate_u_plane, "median_dQdx_of_TPCObject_muon_candidate_u_plane/F");
  background_tree->Branch("median_dQdx_of_TPCObject_muon_candidate_v_plane", &median_dQdx_of_TPCObject_muon_candidate_v_plane, "median_dQdx_of_TPCObject_muon_candidate_v_plane/F");
  background_tree->Branch("median_dQdx_of_TPCObject_muon_candidate_y_plane", &median_dQdx_of_TPCObject_muon_candidate_y_plane, "median_dQdx_of_TPCObject_muon_candidate_y_plane/F");
  background_tree->Branch("median_dQdx_of_TPCObject_muon_candidate_sum", &median_dQdx_of_TPCObject_muon_candidate_sum, "median_dQdx_of_TPCObject_muon_candidate_sum/F");
  background_tree->Branch("u_plane_hit_sum", &u_plane_hit_sum, "u_plane_hit_sum/F");
  background_tree->Branch("v_plane_hit_sum", &v_plane_hit_sum, "v_plane_hit_sum/F");
  background_tree->Branch("y_plane_hit_sum", &y_plane_hit_sum, "y_plane_hit_sum/F");
  background_tree->Branch("hit_sum", &hit_sum, "hit_sum/F");
  background_tree->Branch("u_plane_vicinity_hit_sum", &u_plane_vicinity_hit_sum, "u_plane_vicinity_hit_sum/F");
  background_tree->Branch("v_plane_vicinity_hit_sum", &v_plane_vicinity_hit_sum, "v_plane_vicinity_hit_sum/F");
  background_tree->Branch("y_plane_vicinity_hit_sum", &y_plane_vicinity_hit_sum, "y_plane_vicinity_hit_sum/F");
  background_tree->Branch("vicinity_hit_sum", &vicinity_hit_sum, "vicinity_hit_sum/F");
  background_tree->Branch("u_plane_vicinity_hit_sum_difference", &u_plane_vicinity_hit_sum_difference, "u_plane_vicinity_hit_sum_difference/F");
  background_tree->Branch("v_plane_vicinity_hit_sum_difference", &v_plane_vicinity_hit_sum_difference, "v_plane_vicinity_hit_sum_difference/F");
  background_tree->Branch("y_plane_vicinity_hit_sum_difference", &y_plane_vicinity_hit_sum_difference, "y_plane_vicinity_hit_sum_difference/F");
  background_tree->Branch("vicinity_hit_sum_difference", &vicinity_hit_sum_difference, "vicinity_hit_sum_difference/F");
  background_tree->Branch("u_plane_number_of_slice_hits_in_vicinity_ADCs", &u_plane_number_of_slice_hits_in_vicinity_ADCs, "u_plane_number_of_slice_hits_in_vicinity_ADCs/F");
  background_tree->Branch("v_plane_number_of_slice_hits_in_vicinity_ADCs", &v_plane_number_of_slice_hits_in_vicinity_ADCs, "v_plane_number_of_slice_hits_in_vicinity_ADCs/F");
  background_tree->Branch("y_plane_number_of_slice_hits_in_vicinity_ADCs", &y_plane_number_of_slice_hits_in_vicinity_ADCs, "y_plane_number_of_slice_hits_in_vicinity_ADCs/F");
  background_tree->Branch("total_number_of_slice_hits_in_vicinity_ADCs", &total_number_of_slice_hits_in_vicinity_ADCs, "total_number_of_slice_hits_in_vicinity_ADCs/F");
  background_tree->Branch("num_of_top_pandora_crossings", &num_of_top_pandora_crossings, "num_of_top_pandora_crossings/F");
  background_tree->Branch("num_of_bottom_pandora_crossings", &num_of_bottom_pandora_crossings, "num_of_bottom_pandora_crossings/F");
  background_tree->Branch("num_of_front_pandora_crossings", &num_of_front_pandora_crossings, "num_of_front_pandora_crossings/F");
  background_tree->Branch("num_of_back_pandora_crossings", &num_of_back_pandora_crossings, "num_of_back_pandora_crossings/F");
  background_tree->Branch("number_of_tracks_in_TPCObject", &number_of_tracks_in_TPCObject, "number_of_tracks_in_TPCObject/I");
  background_tree->Branch("num_of_tracks_originating_from_vertex", &num_of_tracks_originating_from_vertex, "num_of_tracks_originating_from_vertex/I");
  background_tree->Branch("muon_candidate_length", &muon_candidate_length, "muon_candidate_length/F");
  background_tree->Branch("sum_of_TPCObject_track_lengths", &sum_of_TPCObject_track_lengths, "sum_of_TPCObject_track_lengths/F");
  background_tree->Branch("total_length_of_tracks_originating_from_vertex", &total_length_of_tracks_originating_from_vertex, "total_length_of_tracks_originating_from_vertex/F");
  background_tree->Branch("rootino_fix_mcweight", &rootino_fix_mcweight, "rootino_fix_mcweight/F");
  background_tree->Branch("central_value_mcweight", &central_value_mcweight, "central_value_mcweight/F");
  background_tree->Branch("other_universe_mcweights", &other_universe_mcweights, "other_universe_mcweights[100]/F");
  background_tree->Branch("axialff_mcweights", &axialff_mcweights, "axialff_mcweights[2]/F");
  background_tree->Branch("rpaccqe_mcweights", &rpaccqe_mcweights, "rpaccqe_mcweights[2]/F");
  background_tree->Branch("xsecshape_mcweights", &xsecshape_mcweights, "xsecshape_mcweights[2]/F");

  // Declare the variables that are input from the background events.
  int    Background_input_run;
  int    Background_input_subrun;
  int    Background_input_event;
  int    Background_input_NC_Channel;
  float  Background_input_truth_neutrino_energy;
  int    Background_input_num_muminus_tracks;
  int    Background_input_num_muplus_tracks;
  int    Background_input_num_piplus_tracks;
  int    Background_input_num_piminus_tracks;
  int    Background_input_num_pi0_tracks;
  int    Background_input_num_proton_tracks;
  int    Background_input_num_electron_showers;
  int    Background_input_num_positron_showers;
  int    Background_input_num_photon_showers;
  float  Background_input_spline_fix_mcweight;
  float  Background_input_truth_reco_vtx_distance;
  float  Background_input_neutrino_energy;
  float  Background_input_neutrino_cos_angle;
  float  Background_input_reco_vtx_x;
  float  Background_input_reco_vtx_y;
  float  Background_input_reco_vtx_z;
  float  Background_input_muon_x_momentum_normalized;
  float  Background_input_muon_y_momentum_normalized;
  float  Background_input_muon_z_momentum_normalized;
  float  Background_input_x_component_of_muon_momentum;
  float  Background_input_y_component_of_muon_momentum;
  float  Background_input_z_component_of_muon_momentum;
  float  Background_input_num_beam_flashes;
  float  Background_input_flash_PEs;
  float  Background_input_flash_z;
  float  Background_input_flash_y;
  float  Background_input_truncated_dQdx_of_TPCObject_muon_candidate_u_plane;
  float  Background_input_truncated_dQdx_of_TPCObject_muon_candidate_v_plane;
  float  Background_input_truncated_dQdx_of_TPCObject_muon_candidate_y_plane;
  float  Background_input_truncated_dQdx_of_TPCObject_muon_candidate_sum;
  float  Background_input_median_dQdx_of_TPCObject_muon_candidate_u_plane;
  float  Background_input_median_dQdx_of_TPCObject_muon_candidate_v_plane;
  float  Background_input_median_dQdx_of_TPCObject_muon_candidate_y_plane;
  float  Background_input_median_dQdx_of_TPCObject_muon_candidate_sum;
  float  Background_input_u_plane_hit_sum;
  float  Background_input_v_plane_hit_sum;
  float  Background_input_y_plane_hit_sum;
  float  Background_input_hit_sum;
  float  Background_input_u_plane_vicinity_hit_sum;
  float  Background_input_v_plane_vicinity_hit_sum;
  float  Background_input_y_plane_vicinity_hit_sum;
  float  Background_input_vicinity_hit_sum;
  float  Background_input_u_plane_number_of_slice_hits_in_vicinity_ADCs;
  float  Background_input_v_plane_number_of_slice_hits_in_vicinity_ADCs;
  float  Background_input_y_plane_number_of_slice_hits_in_vicinity_ADCs;
  float  Background_input_total_number_of_slice_hits_in_vicinity_ADCs;
  float  Background_input_num_of_top_pandora_crossings;
  float  Background_input_num_of_bottom_pandora_crossings;
  float  Background_input_num_of_front_pandora_crossings;
  float  Background_input_num_of_back_pandora_crossings;
  int    Background_input_number_of_tracks_in_TPCObject;
  int    Background_input_num_of_tracks_originating_from_vertex;
  float  Background_input_muon_candidate_length;
  float  Background_input_sum_of_TPCObject_track_lengths;
  float  Background_input_total_length_of_tracks_originating_from_vertex;
  float  Background_input_rootino_fix_mcweight;
  float  Background_input_central_value_mcweight;
  float  Background_input_other_universe_mcweights[100];
  float  Background_input_axialff_mcweights[2];
  float  Background_input_rpaccqe_mcweights[2];
  float  Background_input_xsecshape_mcweights[2];

  TChain* Background_input_tree = new TChain("t0ana/tree");
  Background_input_tree->Add("/home/barnchri/Performing_Analysis_With_MCWeights/NuMI_KDAR_Test_With_All_Weights.root");
  Background_input_tree->SetBranchAddress("run", &Background_input_run);
  Background_input_tree->SetBranchAddress("subrun", &Background_input_subrun);
  Background_input_tree->SetBranchAddress("event", &Background_input_event);
  Background_input_tree->SetBranchAddress("NC_channel", &Background_input_NC_Channel);
  Background_input_tree->SetBranchAddress("truth_energy_in_tree", &Background_input_truth_neutrino_energy);
  Background_input_tree->SetBranchAddress("NC_channel", &Background_input_NC_Channel);
  Background_input_tree->SetBranchAddress("num_muminus_tracks", &Background_input_num_muminus_tracks);
  Background_input_tree->SetBranchAddress("num_muplus_tracks", &Background_input_num_muplus_tracks);
  Background_input_tree->SetBranchAddress("num_piplus_tracks", &Background_input_num_piplus_tracks);
  Background_input_tree->SetBranchAddress("num_piminus_tracks", &Background_input_num_piminus_tracks);
  Background_input_tree->SetBranchAddress("num_pi0_tracks", &Background_input_num_pi0_tracks);
  Background_input_tree->SetBranchAddress("num_proton_tracks", &Background_input_num_proton_tracks);
  Background_input_tree->SetBranchAddress("num_electron_showers", &Background_input_num_electron_showers);
  Background_input_tree->SetBranchAddress("num_positron_showers", &Background_input_num_positron_showers);
  Background_input_tree->SetBranchAddress("num_photon_showers", &Background_input_num_photon_showers);
  Background_input_tree->SetBranchAddress("spline_fix_mcweight", &Background_input_spline_fix_mcweight);
  Background_input_tree->SetBranchAddress("truth_reco_vtx_distance", &Background_input_truth_reco_vtx_distance);
  Background_input_tree->SetBranchAddress("reco_energy_in_tree", &Background_input_neutrino_energy);
  Background_input_tree->SetBranchAddress("cos_angle_between_truth_and_reco_direction_vectors_in_tree", &Background_input_neutrino_cos_angle);
  Background_input_tree->SetBranchAddress("reco_vtx_x", &Background_input_reco_vtx_x);
  Background_input_tree->SetBranchAddress("reco_vtx_y", &Background_input_reco_vtx_y);
  Background_input_tree->SetBranchAddress("reco_vtx_z", &Background_input_reco_vtx_z);
  Background_input_tree->SetBranchAddress("muon_x_momentum_normalized", &Background_input_muon_x_momentum_normalized);
  Background_input_tree->SetBranchAddress("muon_y_momentum_normalized", &Background_input_muon_y_momentum_normalized);
  Background_input_tree->SetBranchAddress("muon_z_momentum_normalized", &Background_input_muon_z_momentum_normalized);
  Background_input_tree->SetBranchAddress("x_component_of_muon_momentum", &Background_input_x_component_of_muon_momentum);
  Background_input_tree->SetBranchAddress("y_component_of_muon_momentum", &Background_input_y_component_of_muon_momentum);
  Background_input_tree->SetBranchAddress("z_component_of_muon_momentum", &Background_input_z_component_of_muon_momentum);
  Background_input_tree->SetBranchAddress("num_beam_flashes", &Background_input_num_beam_flashes);
  Background_input_tree->SetBranchAddress("flash_PEs", &Background_input_flash_PEs);
  Background_input_tree->SetBranchAddress("flash_z", &Background_input_flash_z);
  Background_input_tree->SetBranchAddress("flash_y", &Background_input_flash_y);
  Background_input_tree->SetBranchAddress("truncated_dQdx_of_TPCObject_muon_candidate_u_plane", &Background_input_truncated_dQdx_of_TPCObject_muon_candidate_u_plane);
  Background_input_tree->SetBranchAddress("truncated_dQdx_of_TPCObject_muon_candidate_v_plane", &Background_input_truncated_dQdx_of_TPCObject_muon_candidate_v_plane);
  Background_input_tree->SetBranchAddress("truncated_dQdx_of_TPCObject_muon_candidate_y_plane", &Background_input_truncated_dQdx_of_TPCObject_muon_candidate_y_plane);
  Background_input_tree->SetBranchAddress("truncated_dQdx_of_TPCObject_muon_candidate_sum", &Background_input_truncated_dQdx_of_TPCObject_muon_candidate_sum);
  Background_input_tree->SetBranchAddress("median_dQdx_of_TPCObject_muon_candidate_u_plane", &Background_input_median_dQdx_of_TPCObject_muon_candidate_u_plane);
  Background_input_tree->SetBranchAddress("median_dQdx_of_TPCObject_muon_candidate_v_plane", &Background_input_median_dQdx_of_TPCObject_muon_candidate_v_plane);
  Background_input_tree->SetBranchAddress("median_dQdx_of_TPCObject_muon_candidate_y_plane", &Background_input_median_dQdx_of_TPCObject_muon_candidate_y_plane);
  Background_input_tree->SetBranchAddress("median_dQdx_of_TPCObject_muon_candidate_sum", &Background_input_median_dQdx_of_TPCObject_muon_candidate_sum);
  Background_input_tree->SetBranchAddress("ADC_hit_sum_u_plane", &Background_input_u_plane_hit_sum);
  Background_input_tree->SetBranchAddress("ADC_hit_sum_v_plane", &Background_input_v_plane_hit_sum);
  Background_input_tree->SetBranchAddress("ADC_hit_sum_y_plane", &Background_input_y_plane_hit_sum);
  Background_input_tree->SetBranchAddress("ADC_hit_sum_total", &Background_input_hit_sum);
  Background_input_tree->SetBranchAddress("ADC_hits_vicinity_sum_u_plane", &Background_input_u_plane_vicinity_hit_sum);
  Background_input_tree->SetBranchAddress("ADC_hits_vicinity_sum_v_plane", &Background_input_v_plane_vicinity_hit_sum);
  Background_input_tree->SetBranchAddress("ADC_hits_vicinity_sum_y_plane", &Background_input_y_plane_vicinity_hit_sum);
  Background_input_tree->SetBranchAddress("hits_vicinity_num_u_plane", &Background_input_u_plane_number_of_slice_hits_in_vicinity_ADCs);
  Background_input_tree->SetBranchAddress("hits_vicinity_num_v_plane", &Background_input_v_plane_number_of_slice_hits_in_vicinity_ADCs);
  Background_input_tree->SetBranchAddress("hits_vicinity_num_y_plane", &Background_input_y_plane_number_of_slice_hits_in_vicinity_ADCs);
  Background_input_tree->SetBranchAddress("hits_vicinity_num_total", &Background_input_total_number_of_slice_hits_in_vicinity_ADCs);
  Background_input_tree->SetBranchAddress("num_of_top_pandora_crossings", &Background_input_num_of_top_pandora_crossings);
  Background_input_tree->SetBranchAddress("num_of_bottom_pandora_crossings", &Background_input_num_of_bottom_pandora_crossings);
  Background_input_tree->SetBranchAddress("num_of_front_pandora_crossings", &Background_input_num_of_front_pandora_crossings);
  Background_input_tree->SetBranchAddress("num_of_back_pandora_crossings", &Background_input_num_of_back_pandora_crossings);
  Background_input_tree->SetBranchAddress("number_of_tracks_in_TPCObject", &Background_input_number_of_tracks_in_TPCObject);
  Background_input_tree->SetBranchAddress("num_of_tracks_originating_from_vertex", &Background_input_num_of_tracks_originating_from_vertex);
  Background_input_tree->SetBranchAddress("muon_candidate_length", &Background_input_muon_candidate_length);
  Background_input_tree->SetBranchAddress("sum_of_TPCObject_track_lengths", &Background_input_sum_of_TPCObject_track_lengths);
  Background_input_tree->SetBranchAddress("total_length_of_tracks_originating_from_vertex", &Background_input_total_length_of_tracks_originating_from_vertex);
  Background_input_tree->SetBranchAddress("rootino_fix_mcweight", &Background_input_rootino_fix_mcweight);
  Background_input_tree->SetBranchAddress("central_value_mcweight", &Background_input_central_value_mcweight);
  Background_input_tree->SetBranchAddress("other_universe_mcweights", &Background_input_other_universe_mcweights);
  Background_input_tree->SetBranchAddress("axialff_mcweights", &Background_input_axialff_mcweights);
  Background_input_tree->SetBranchAddress("rpaccqe_mcweights", &Background_input_rpaccqe_mcweights);
  Background_input_tree->SetBranchAddress("xsecshape_mcweights", &Background_input_xsecshape_mcweights);

  int Background_input_entries = Background_input_tree->GetEntries();

  special_kdar_weight     = 1.0;
  sampleweight            = 0.441;
  treeweight              = 3.682;
  ext                     = 1;
  
  for ( size_t iter = 0; iter < Background_input_entries; iter++ ) {

    // Set the tree to the current entry.
    Background_input_tree->GetEntry( iter );

    if ( Background_input_flash_PEs < 50.0 || Background_input_flash_PEs > 2000. || Background_input_sum_of_TPCObject_track_lengths > 65.0 ) continue;

    // Set all of the variables and fill the new tree.
    run                                                = Background_input_run;
    subrun                                             = Background_input_subrun;
    event                                              = Background_input_event;
    NC_channel                                         = Background_input_NC_Channel;
    truth_neutrino_energy                              = Background_input_truth_neutrino_energy;
    num_muminus_tracks                                 = Background_input_num_muminus_tracks;
    num_muplus_tracks                                  = Background_input_num_muplus_tracks;
    num_piplus_tracks                                  = Background_input_num_piplus_tracks;
    num_piminus_tracks                                 = Background_input_num_piminus_tracks;
    num_pi0_tracks                                     = Background_input_num_pi0_tracks;
    num_proton_tracks                                  = Background_input_num_proton_tracks;
    num_electron_tracks                                = Background_input_num_electron_showers;
    num_positron_tracks                                = Background_input_num_positron_showers;
    num_photon_tracks                                  = Background_input_num_photon_showers;
    spline_fix_mcweight                                = Background_input_spline_fix_mcweight;
    truth_reco_vtx_distance                            = Background_input_truth_reco_vtx_distance;
    neutrino_energy                                    = Background_input_neutrino_energy;
    neutrino_cos_angle                                 = Background_input_neutrino_cos_angle;
    reco_vtx_x                                         = Background_input_reco_vtx_x;
    reco_vtx_y                                         = Background_input_reco_vtx_y;
    reco_vtx_z                                         = Background_input_reco_vtx_z;
    muon_x_momentum_normalized                         = Background_input_muon_x_momentum_normalized;
    muon_y_momentum_normalized                         = Background_input_muon_y_momentum_normalized;
    muon_z_momentum_normalized                         = Background_input_muon_z_momentum_normalized;
    x_component_of_muon_momentum                       = Background_input_x_component_of_muon_momentum;
    y_component_of_muon_momentum                       = Background_input_y_component_of_muon_momentum;
    z_component_of_muon_momentum                       = Background_input_z_component_of_muon_momentum;
    num_beam_flashes                                   = Background_input_num_beam_flashes;
    flash_PEs                                          = Background_input_flash_PEs;
    flash_z                                            = Background_input_flash_z;
    flash_y                                            = Background_input_flash_y;
    truncated_dQdx_of_TPCObject_muon_candidate_u_plane = Background_input_truncated_dQdx_of_TPCObject_muon_candidate_u_plane;
    truncated_dQdx_of_TPCObject_muon_candidate_v_plane = Background_input_truncated_dQdx_of_TPCObject_muon_candidate_v_plane;
    truncated_dQdx_of_TPCObject_muon_candidate_y_plane = Background_input_truncated_dQdx_of_TPCObject_muon_candidate_y_plane;
    truncated_dQdx_of_TPCObject_muon_candidate_sum     = Background_input_truncated_dQdx_of_TPCObject_muon_candidate_sum;
    median_dQdx_of_TPCObject_muon_candidate_u_plane    = Background_input_median_dQdx_of_TPCObject_muon_candidate_u_plane;
    median_dQdx_of_TPCObject_muon_candidate_v_plane    = Background_input_median_dQdx_of_TPCObject_muon_candidate_v_plane;
    median_dQdx_of_TPCObject_muon_candidate_y_plane    = Background_input_median_dQdx_of_TPCObject_muon_candidate_y_plane;
    median_dQdx_of_TPCObject_muon_candidate_sum        = Background_input_median_dQdx_of_TPCObject_muon_candidate_sum;
    u_plane_hit_sum                                    = Background_input_u_plane_hit_sum;
    v_plane_hit_sum                                    = Background_input_v_plane_hit_sum;
    y_plane_hit_sum                                    = Background_input_y_plane_hit_sum;
    hit_sum                                            = Background_input_hit_sum;
    u_plane_vicinity_hit_sum                           = Background_input_u_plane_vicinity_hit_sum;
    v_plane_vicinity_hit_sum                           = Background_input_v_plane_vicinity_hit_sum;
    y_plane_vicinity_hit_sum                           = Background_input_y_plane_vicinity_hit_sum;
    vicinity_hit_sum                                   = Background_input_vicinity_hit_sum;
    u_plane_number_of_slice_hits_in_vicinity_ADCs      = Background_input_u_plane_number_of_slice_hits_in_vicinity_ADCs;
    v_plane_number_of_slice_hits_in_vicinity_ADCs      = Background_input_v_plane_number_of_slice_hits_in_vicinity_ADCs;
    y_plane_number_of_slice_hits_in_vicinity_ADCs      = Background_input_y_plane_number_of_slice_hits_in_vicinity_ADCs;
    total_number_of_slice_hits_in_vicinity_ADCs        = Background_input_total_number_of_slice_hits_in_vicinity_ADCs;
    num_of_top_pandora_crossings                       = Background_input_num_of_top_pandora_crossings;
    num_of_bottom_pandora_crossings                    = Background_input_num_of_bottom_pandora_crossings;
    num_of_front_pandora_crossings                     = Background_input_num_of_front_pandora_crossings;
    num_of_back_pandora_crossings                      = Background_input_num_of_back_pandora_crossings;
    number_of_tracks_in_TPCObject                      = Background_input_number_of_tracks_in_TPCObject;
    num_of_tracks_originating_from_vertex              = Background_input_num_of_tracks_originating_from_vertex;
    muon_candidate_length                              = Background_input_muon_candidate_length;
    sum_of_TPCObject_track_lengths                     = Background_input_sum_of_TPCObject_track_lengths;
    total_length_of_tracks_originating_from_vertex     = Background_input_total_length_of_tracks_originating_from_vertex;

    flash_z_difference                                 = ( reco_vtx_z - flash_z );
    flash_y_difference                                 = ( reco_vtx_y - flash_y );

    vicinity_hit_sum_difference                        = ( vicinity_hit_sum - hit_sum );
    u_plane_vicinity_hit_sum_difference                = ( u_plane_vicinity_hit_sum - u_plane_hit_sum );
    v_plane_vicinity_hit_sum_difference                = ( v_plane_vicinity_hit_sum - v_plane_hit_sum );
    y_plane_vicinity_hit_sum_difference                = ( y_plane_vicinity_hit_sum - y_plane_hit_sum );

    // Fill the other weighting variables.
    rootino_fix_mcweight                               = Background_input_rootino_fix_mcweight;
    central_value_mcweight                             = Background_input_central_value_mcweight;

    for ( size_t i = 0; i < 100; i++ ) {

      other_universe_mcweights[i] = Background_input_other_universe_mcweights[i];

    }

    for ( size_t i = 0; i < 2; i++ ) {

      axialff_mcweights[i] = Background_input_axialff_mcweights[i];

    }

    for ( size_t i = 0; i < 2; i++ ) {

      rpaccqe_mcweights[i] = Background_input_rpaccqe_mcweights[i];

    }

    for ( size_t i = 0; i < 2; i++ ) {

      xsecshape_mcweights[i] = Background_input_xsecshape_mcweights[i];

    }

    // Change the MCWeight Based On If It's a KDAR Event.
    if ( fabs( truth_neutrino_energy - 235.532 ) < 0.001 ) special_kdar_weight = 5.;

    background_tree->Fill();

  }

  std::cout << "Number of entries in signal tree = " << signal_tree->GetEntries() << "." << std::endl;
  std::cout << "Number of entries in background tree = " << background_tree->GetEntries() << "." << std::endl;

  output_file->Write();
  output_file->Close();

}
