#include <cstdlib>
#include <vector>
#include <iostream>
#include <map>
#include <string>
#include <fstream>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TStopwatch.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TMarker.h"

#if not defined(__CINT__) || defined(__MAKECINT__)
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"
#endif

using namespace TMVA;
using namespace std;

void TMVAClassificationApplication_MC_edited_background_only( TString myMethodList = "" ) {   

#ifdef __CINT__
  gROOT->ProcessLine( ".O0" ); // turn off optimization in CINT
#endif
  //---------------------------------------------------------------
  // This loads the library
  TMVA::Tools::Instance();
  
  // Default MVA methods to be trained + tested
  std::map<std::string,int> Use;
  
  // --- Boosted Decision Trees
  Use["BDTG"] = 1; // uses Gradient Boost
  std::cout << std::endl;
  std::cout << "==> Start TMVAClassificationApplication" << std::endl;
  
  // Select methods (don't look at this code - not of interest)
  if (myMethodList != "") {
  
    for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) it->second = 0;
    std::vector<TString> mlist = gTools().SplitString( myMethodList, ',' );
    
    for (UInt_t i=0; i<mlist.size(); i++) {
      std::string regMethod(mlist[i]);
      if (Use.find(regMethod) == Use.end()) {
	std::cout << "Method \"" << regMethod 
		  << "\" not known in TMVA under this name. Choose among the following:" << std::endl;

	for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) {
	  std::cout << it->first << " ";
	}

	std::cout << std::endl;
	return;
      }
      Use[regMethod] = 1;
    }
  }

  // --- Create the Reader object
  TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );    

  // Create a set of variables and declare them to the reader
  // - the variable names MUST corresponds in name and type to those given in the weight file(s) used

  int   run;
  int   subrun;
  int   event;
  int   NC_channel;
  float truth_neutrino_energy;
  float truth_reco_vtx_distance;
  double truth_theta_angle_for_weighting;
  int   num_muminus_tracks;
  int   num_muplus_tracks;
  int   num_piplus_tracks;
  int   num_piminus_tracks;
  int   num_pi0_tracks;
  int   num_proton_tracks;
  int   num_electron_tracks;
  int   num_positron_tracks;
  int   num_photon_tracks;
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
  float flash_y_difference;
  float flash_z_difference;
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
  
  //reader->AddVariable("neutrino_energy", &neutrino_energy);
  reader->AddVariable("neutrino_cos_angle", &neutrino_cos_angle);
  reader->AddVariable("reco_vtx_x", &reco_vtx_x);
  reader->AddVariable("reco_vtx_y", &reco_vtx_y);
  reader->AddVariable("reco_vtx_z", &reco_vtx_z);
  reader->AddVariable("x_component_of_muon_momentum", &x_component_of_muon_momentum);
  reader->AddVariable("y_component_of_muon_momentum", &y_component_of_muon_momentum);
  reader->AddVariable("z_component_of_muon_momentum", &z_component_of_muon_momentum);
  //reader->AddVariable("flash_PEs", &flash_PEs);
  reader->AddVariable("flash_z_difference", &flash_z_difference);
  reader->AddVariable("flash_y_difference", &flash_y_difference);
  reader->AddVariable("truncated_dQdx_of_TPCObject_muon_candidate_u_plane", &truncated_dQdx_of_TPCObject_muon_candidate_u_plane);
  reader->AddVariable("truncated_dQdx_of_TPCObject_muon_candidate_v_plane", &truncated_dQdx_of_TPCObject_muon_candidate_v_plane);
  reader->AddVariable("truncated_dQdx_of_TPCObject_muon_candidate_y_plane", &truncated_dQdx_of_TPCObject_muon_candidate_y_plane);
  reader->AddVariable("truncated_dQdx_of_TPCObject_muon_candidate_sum", &truncated_dQdx_of_TPCObject_muon_candidate_sum);
  reader->AddVariable("hit_sum", &hit_sum);
  reader->AddVariable("u_plane_hit_sum", &u_plane_hit_sum);
  reader->AddVariable("v_plane_hit_sum", &v_plane_hit_sum);
  reader->AddVariable("y_plane_hit_sum", &y_plane_hit_sum);
  reader->AddVariable("vicinity_hit_sum_difference", &vicinity_hit_sum_difference);
  reader->AddVariable("u_plane_vicinity_hit_sum_difference", &u_plane_vicinity_hit_sum_difference);
  reader->AddVariable("v_plane_vicinity_hit_sum_difference", &v_plane_vicinity_hit_sum_difference);
  reader->AddVariable("y_plane_vicinity_hit_sum_difference", &y_plane_vicinity_hit_sum_difference);
  reader->AddVariable("total_number_of_slice_hits_in_vicinity_ADCs", &total_number_of_slice_hits_in_vicinity_ADCs);
  reader->AddVariable("u_plane_number_of_slice_hits_in_vicinity_ADCs", &u_plane_number_of_slice_hits_in_vicinity_ADCs);
  reader->AddVariable("v_plane_number_of_slice_hits_in_vicinity_ADCs", &v_plane_number_of_slice_hits_in_vicinity_ADCs);
  reader->AddVariable("y_plane_number_of_slice_hits_in_vicinity_ADCs", &y_plane_number_of_slice_hits_in_vicinity_ADCs);
  
  // Book method(s)
  reader->BookMVA("BDTG method", "dataset_calibrated_all_stats/weights/TMVAClassification_minus_BDTG.weights.xml" ); // WEIGHT FILE IS GREEN FILE .XML
  
  // Prepare input tree (this must be replaced by your data source)
  // in this example, there is a toy tree with signal and one with background events
  // we'll later on use only the "signal" events for the test in this example.
  //

  int   output_run;
  int   output_subrun;
  int   output_event;
  float output_score;
  int   output_NC_channel;      
  int   output_num_muminus_tracks;
  int   output_num_muplus_tracks;
  int   output_num_piplus_tracks;
  int   output_num_piminus_tracks;
  int   output_num_pi0_tracks;
  int   output_num_proton_tracks;
  int   output_num_electron_tracks;
  int   output_num_positron_tracks;
  int   output_num_photon_tracks;
  int   output_ext;
  double output_truth_theta_angle_for_weighting;
  float output_truth_reco_vtx_distance;
  float output_spline_fix_mcweight;
  float output_special_kdar_weight;
  float output_sampleweight;
  float output_treeweight;
  float output_truth_neutrino_energy;
  float output_neutrino_energy;
  float output_neutrino_cos_angle;
  float output_reco_vtx_x;
  float output_reco_vtx_y;
  float output_reco_vtx_z;
  float output_muon_x_momentum_normalized;
  float output_muon_y_momentum_normalized;
  float output_muon_z_momentum_normalized;
  float output_x_component_of_muon_momentum;
  float output_y_component_of_muon_momentum;
  float output_z_component_of_muon_momentum;
  float output_num_beam_flashes;
  float output_flash_PEs;
  float output_flash_z;
  float output_flash_y;
  float output_flash_z_difference;
  float output_flash_y_difference;
  float output_truncated_dQdx_of_TPCObject_muon_candidate_u_plane;
  float output_truncated_dQdx_of_TPCObject_muon_candidate_v_plane;
  float output_truncated_dQdx_of_TPCObject_muon_candidate_y_plane;
  float output_truncated_dQdx_of_TPCObject_muon_candidate_sum;
  float output_median_dQdx_of_TPCObject_muon_candidate_u_plane;
  float output_median_dQdx_of_TPCObject_muon_candidate_v_plane;
  float output_median_dQdx_of_TPCObject_muon_candidate_y_plane;
  float output_median_dQdx_of_TPCObject_muon_candidate_sum;
  float output_u_plane_hit_sum;
  float output_v_plane_hit_sum;
  float output_y_plane_hit_sum;
  float output_hit_sum;
  float output_u_plane_vicinity_hit_sum_difference;
  float output_v_plane_vicinity_hit_sum_difference;
  float output_y_plane_vicinity_hit_sum_difference;
  float output_vicinity_hit_sum_difference;
  float output_u_plane_number_of_slice_hits_in_vicinity_ADCs;
  float output_v_plane_number_of_slice_hits_in_vicinity_ADCs;
  float output_y_plane_number_of_slice_hits_in_vicinity_ADCs;
  float output_total_number_of_slice_hits_in_vicinity_ADCs;
  float output_num_of_top_pandora_crossings;
  float output_num_of_bottom_pandora_crossings;
  float output_num_of_front_pandora_crossings;
  float output_num_of_back_pandora_crossings;
  int   output_number_of_tracks_in_TPCObject;
  int   output_num_of_tracks_originating_from_vertex;
  float output_muon_candidate_length;
  float output_sum_of_TPCObject_track_lengths;
  float output_total_length_of_tracks_originating_from_vertex;
  float output_rootino_fix_mcweight;
  float output_central_value_mcweight;
  float output_other_universe_mcweights[100];
  float output_axialff_mcweights[2];
  float output_rpaccqe_mcweights[2];
  float output_xsecshape_mcweights[2];
  
  TFile* output_file        = new TFile("output_file_cv_testing_sample.root", "RECREATE");
  
  TTree* output_bkgd_tree   = new TTree("background_tree", "A Tree With Information on Which Backgrounds Are Passing");

  output_bkgd_tree->Branch("run", &output_run, "run/I");
  output_bkgd_tree->Branch("subrun", &output_subrun, "subrun/I");
  output_bkgd_tree->Branch("event", &output_event, "event/I");
  output_bkgd_tree->Branch("score", &output_score, "score/F");
  output_bkgd_tree->Branch("truth_neutrino_energy", &output_truth_neutrino_energy, "truth_neutrino_energy/F");
  output_bkgd_tree->Branch("NC_channel", &output_NC_channel, "NC_channel/I");
  output_bkgd_tree->Branch("truth_reco_vtx_distance", &output_truth_reco_vtx_distance, "truth_reco_vtx_distance/F");
  output_bkgd_tree->Branch("truth_theta_angle_for_weighting", &output_truth_theta_angle_for_weighting, "truth_theta_angle_for_weighting/D");
  output_bkgd_tree->Branch("num_muminus_tracks", &output_num_muminus_tracks, "num_muminus_tracks/I");
  output_bkgd_tree->Branch("num_muplus_tracks",&output_num_muplus_tracks, "num_muplus_tracks/I");
  output_bkgd_tree->Branch("num_piplus_tracks", &output_num_piplus_tracks, "num_piplus_tracks/I");
  output_bkgd_tree->Branch("num_piminus_tracks",&output_num_piminus_tracks, "num_piminus_tracks/I");
  output_bkgd_tree->Branch("num_pi0_tracks", &output_num_pi0_tracks, "num_pi0_tracks/I");
  output_bkgd_tree->Branch("num_proton_tracks", &output_num_proton_tracks, "num_proton_tracks/I");
  output_bkgd_tree->Branch("num_electron_tracks", &output_num_electron_tracks, "num_electron_tracks/I");
  output_bkgd_tree->Branch("num_positron_tracks", &output_num_positron_tracks, "num_positron_tracks/I");
  output_bkgd_tree->Branch("num_photon_tracks", &output_num_photon_tracks, "num_photon_tracks/I");
  output_bkgd_tree->Branch("spline_fix_mcweight", &output_spline_fix_mcweight, "spline_fix_mcweight/F");
  output_bkgd_tree->Branch("special_kdar_weight", &output_special_kdar_weight, "special_kdar_weight/F");
  output_bkgd_tree->Branch("sampleweight", &output_sampleweight, "sampleweight/F");
  output_bkgd_tree->Branch("treeweight", &output_treeweight, "treeweight/F");
  output_bkgd_tree->Branch("ext", &output_ext, "ext/I");
  output_bkgd_tree->Branch("neutrino_energy", &output_neutrino_energy, "neutrino_energy/F");
  output_bkgd_tree->Branch("neutrino_cos_angle", &output_neutrino_cos_angle, "neutrino_cos_angle/F");
  output_bkgd_tree->Branch("reco_vtx_x", &output_reco_vtx_x, "reco_vtx_x/F");
  output_bkgd_tree->Branch("reco_vtx_y", &output_reco_vtx_y, "reco_vtx_y/F");
  output_bkgd_tree->Branch("reco_vtx_z", &output_reco_vtx_z, "reco_vtx_z/F");
  output_bkgd_tree->Branch("muon_x_momentum_normalized", &output_muon_x_momentum_normalized, "muon_x_momentum_normalized/F");
  output_bkgd_tree->Branch("muon_y_momentum_normalized", &output_muon_y_momentum_normalized, "muon_y_momentum_normalized/F");
  output_bkgd_tree->Branch("muon_z_momentum_normalized", &output_muon_z_momentum_normalized, "muon_z_momentum_normalized/F");
  output_bkgd_tree->Branch("x_component_of_muon_momentum", &output_x_component_of_muon_momentum, "x_component_of_muon_momentum/F");
  output_bkgd_tree->Branch("y_component_of_muon_momentum", &output_y_component_of_muon_momentum, "y_component_of_muon_momentum/F");
  output_bkgd_tree->Branch("z_component_of_muon_momentum", &output_z_component_of_muon_momentum, "z_component_of_muon_momentum/F");
  output_bkgd_tree->Branch("num_beam_flashes", &output_num_beam_flashes, "num_beam_flashes/F");
  output_bkgd_tree->Branch("flash_PEs", &output_flash_PEs, "flash_PEs/F");
  output_bkgd_tree->Branch("flash_z_difference", &output_flash_z_difference, "flash_z_difference/F");
  output_bkgd_tree->Branch("flash_y_difference", &output_flash_y_difference, "flash_y_difference/F");
  output_bkgd_tree->Branch("truncated_dQdx_of_TPCObject_muon_candidate_u_plane", &output_truncated_dQdx_of_TPCObject_muon_candidate_u_plane, "truncated_dQdx_of_TPCObject_muon_candidate_u_plane/F");
  output_bkgd_tree->Branch("truncated_dQdx_of_TPCObject_muon_candidate_v_plane", &output_truncated_dQdx_of_TPCObject_muon_candidate_v_plane, "truncated_dQdx_of_TPCObject_muon_candidate_v_plane/F");
  output_bkgd_tree->Branch("truncated_dQdx_of_TPCObject_muon_candidate_y_plane", &output_truncated_dQdx_of_TPCObject_muon_candidate_y_plane, "truncated_dQdx_of_TPCObject_muon_candidate_y_plane/F");
  output_bkgd_tree->Branch("truncated_dQdx_of_TPCObject_muon_candidate_sum", &output_truncated_dQdx_of_TPCObject_muon_candidate_sum, "truncated_dQdx_of_TPCObject_muon_candidate_sum/F");
  output_bkgd_tree->Branch("median_dQdx_of_TPCObject_muon_candidate_u_plane", &output_median_dQdx_of_TPCObject_muon_candidate_u_plane, "median_dQdx_of_TPCObject_muon_candidate_u_plane/F");
  output_bkgd_tree->Branch("median_dQdx_of_TPCObject_muon_candidate_v_plane", &output_median_dQdx_of_TPCObject_muon_candidate_v_plane, "median_dQdx_of_TPCObject_muon_candidate_v_plane/F");
  output_bkgd_tree->Branch("median_dQdx_of_TPCObject_muon_candidate_y_plane", &output_median_dQdx_of_TPCObject_muon_candidate_y_plane, "median_dQdx_of_TPCObject_muon_candidate_y_plane/F");
  output_bkgd_tree->Branch("median_dQdx_of_TPCObject_muon_candidate_sum", &output_median_dQdx_of_TPCObject_muon_candidate_sum, "median_dQdx_of_TPCObject_muon_candidate_sum/F");
  output_bkgd_tree->Branch("u_plane_hit_sum", &output_u_plane_hit_sum, "u_plane_hit_sum/F");
  output_bkgd_tree->Branch("v_plane_hit_sum", &output_v_plane_hit_sum, "v_plane_hit_sum/F");
  output_bkgd_tree->Branch("y_plane_hit_sum", &output_y_plane_hit_sum, "y_plane_hit_sum/F");
  output_bkgd_tree->Branch("hit_sum", &output_hit_sum, "hit_sum/F");
  output_bkgd_tree->Branch("u_plane_vicinity_hit_sum_difference", &output_u_plane_vicinity_hit_sum_difference, "u_plane_vicinity_hit_sum_difference/F");
  output_bkgd_tree->Branch("v_plane_vicinity_hit_sum_difference", &output_v_plane_vicinity_hit_sum_difference, "v_plane_vicinity_hit_sum_difference/F");
  output_bkgd_tree->Branch("y_plane_vicinity_hit_sum_difference", &output_y_plane_vicinity_hit_sum_difference, "y_plane_vicinity_hit_sum_difference/F");
  output_bkgd_tree->Branch("vicinity_hit_sum_difference", &output_vicinity_hit_sum_difference, "vicinity_hit_sum_difference/F");
  output_bkgd_tree->Branch("u_plane_number_of_slice_hits_in_vicinity_ADCs", &output_u_plane_number_of_slice_hits_in_vicinity_ADCs, "u_plane_number_of_slice_hits_in_vicinity_ADCs/F");
  output_bkgd_tree->Branch("v_plane_number_of_slice_hits_in_vicinity_ADCs", &output_v_plane_number_of_slice_hits_in_vicinity_ADCs, "v_plane_number_of_slice_hits_in_vicinity_ADCs/F");
  output_bkgd_tree->Branch("y_plane_number_of_slice_hits_in_vicinity_ADCs", &output_y_plane_number_of_slice_hits_in_vicinity_ADCs, "y_plane_number_of_slice_hits_in_vicinity_ADCs/F");
  output_bkgd_tree->Branch("total_number_of_slice_hits_in_vicinity_ADCs", &output_total_number_of_slice_hits_in_vicinity_ADCs, "total_number_of_slice_hits_in_vicinity_ADCs/F");
  output_bkgd_tree->Branch("num_of_top_pandora_crossings", &output_num_of_top_pandora_crossings, "num_of_top_pandora_crossings/F");
  output_bkgd_tree->Branch("num_of_bottom_pandora_crossings", &output_num_of_bottom_pandora_crossings, "num_of_bottom_pandora_crossings/F");
  output_bkgd_tree->Branch("num_of_front_pandora_crossings", &output_num_of_front_pandora_crossings, "num_of_front_pandora_crossings/F");
  output_bkgd_tree->Branch("num_of_back_pandora_crossings", &output_num_of_back_pandora_crossings, "num_of_back_pandora_crossings/F");
  output_bkgd_tree->Branch("number_of_tracks_in_TPCObject", &output_number_of_tracks_in_TPCObject, "number_of_tracks_in_TPCObject/I");
  output_bkgd_tree->Branch("num_of_tracks_originating_from_vertex", &output_num_of_tracks_originating_from_vertex, "num_of_tracks_originating_from_vertex/I");
  output_bkgd_tree->Branch("muon_candidate_length", &output_muon_candidate_length, "muon_candidate_length/F");
  output_bkgd_tree->Branch("sum_of_TPCObject_track_lengths", &output_sum_of_TPCObject_track_lengths, "sum_of_TPCObject_track_lengths/F");
  output_bkgd_tree->Branch("length_of_tracks_originating_from_vertex", &output_total_length_of_tracks_originating_from_vertex, "length_of_tracks_originating_from_vertex/F");
  output_bkgd_tree->Branch("rootino_fix_mcweight", &output_rootino_fix_mcweight, "rootino_fix_mcweight/F");
  output_bkgd_tree->Branch("central_value_mcweight", &output_central_value_mcweight, "central_value_mcweight/F");
  output_bkgd_tree->Branch("other_universe_mcweights", &output_other_universe_mcweights, "other_universe_mcweights[100]/F");
  output_bkgd_tree->Branch("axialff_mcweights", &output_axialff_mcweights, "axialff_mcweights[2]/F");
  output_bkgd_tree->Branch("rpaccqe_mcweights", &output_rpaccqe_mcweights, "rpaccqe_mcweights[2]/F");
  output_bkgd_tree->Branch("xsecshape_mcweights", &output_xsecshape_mcweights, "xsecshape_mcweights[2]/F");
  
  TFile *input         = TFile::Open("./input_file_cv_testing_sample.root");
  std::cout << "--- Select background sample" << std::endl;

  TTree* theBackgroundTree = (TTree*)input->Get("background_tree");
  std::cout << "--- TMVAClassificationApp    : Using input file: " << input->GetName() << std::endl; 

  theBackgroundTree->SetBranchAddress("run", &run);
  theBackgroundTree->SetBranchAddress("subrun", &subrun);
  theBackgroundTree->SetBranchAddress("event", &event);
  theBackgroundTree->SetBranchAddress("NC_channel", &NC_channel);
  theBackgroundTree->SetBranchAddress("truth_neutrino_energy", &truth_neutrino_energy);
  theBackgroundTree->SetBranchAddress("num_muminus_tracks", &num_muminus_tracks);
  theBackgroundTree->SetBranchAddress("num_muplus_tracks", &num_muplus_tracks);
  theBackgroundTree->SetBranchAddress("num_piplus_tracks", &num_piplus_tracks);
  theBackgroundTree->SetBranchAddress("num_piminus_tracks", &num_piminus_tracks);
  theBackgroundTree->SetBranchAddress("num_pi0_tracks", &num_pi0_tracks);
  theBackgroundTree->SetBranchAddress("num_proton_tracks", &num_proton_tracks);
  theBackgroundTree->SetBranchAddress("num_electron_tracks", &num_electron_tracks);
  theBackgroundTree->SetBranchAddress("num_positron_tracks", &num_positron_tracks);
  theBackgroundTree->SetBranchAddress("num_photon_tracks", &num_photon_tracks);
  theBackgroundTree->SetBranchAddress("spline_fix_mcweight", &spline_fix_mcweight);
  theBackgroundTree->SetBranchAddress("special_kdar_weight", &special_kdar_weight);
  theBackgroundTree->SetBranchAddress("sampleweight", &sampleweight);
  theBackgroundTree->SetBranchAddress("treeweight", &treeweight);
  theBackgroundTree->SetBranchAddress("ext", &ext);
  theBackgroundTree->SetBranchAddress("truth_reco_vtx_distance", &truth_reco_vtx_distance);
  theBackgroundTree->SetBranchAddress("truth_theta_angle_for_weighting", &truth_theta_angle_for_weighting);
  theBackgroundTree->SetBranchAddress("neutrino_energy", &neutrino_energy);
  theBackgroundTree->SetBranchAddress("neutrino_cos_angle", &neutrino_cos_angle);
  theBackgroundTree->SetBranchAddress("reco_vtx_x", &reco_vtx_x);
  theBackgroundTree->SetBranchAddress("reco_vtx_y", &reco_vtx_y);
  theBackgroundTree->SetBranchAddress("reco_vtx_z", &reco_vtx_z);
  theBackgroundTree->SetBranchAddress("muon_x_momentum_normalized", &muon_x_momentum_normalized);
  theBackgroundTree->SetBranchAddress("muon_y_momentum_normalized", &muon_y_momentum_normalized);
  theBackgroundTree->SetBranchAddress("muon_z_momentum_normalized", &muon_z_momentum_normalized);
  theBackgroundTree->SetBranchAddress("x_component_of_muon_momentum", &x_component_of_muon_momentum);
  theBackgroundTree->SetBranchAddress("y_component_of_muon_momentum", &y_component_of_muon_momentum);
  theBackgroundTree->SetBranchAddress("z_component_of_muon_momentum", &z_component_of_muon_momentum);
  theBackgroundTree->SetBranchAddress("num_beam_flashes", &num_beam_flashes);
  theBackgroundTree->SetBranchAddress("flash_PEs", &flash_PEs);
  theBackgroundTree->SetBranchAddress("flash_z_difference", &flash_z_difference);
  theBackgroundTree->SetBranchAddress("flash_y_difference", &flash_y_difference);
  theBackgroundTree->SetBranchAddress("truncated_dQdx_of_TPCObject_muon_candidate_u_plane", &truncated_dQdx_of_TPCObject_muon_candidate_u_plane);
  theBackgroundTree->SetBranchAddress("truncated_dQdx_of_TPCObject_muon_candidate_v_plane", &truncated_dQdx_of_TPCObject_muon_candidate_v_plane);
  theBackgroundTree->SetBranchAddress("truncated_dQdx_of_TPCObject_muon_candidate_y_plane", &truncated_dQdx_of_TPCObject_muon_candidate_y_plane);
  theBackgroundTree->SetBranchAddress("truncated_dQdx_of_TPCObject_muon_candidate_sum", &truncated_dQdx_of_TPCObject_muon_candidate_sum);
  theBackgroundTree->SetBranchAddress("median_dQdx_of_TPCObject_muon_candidate_u_plane", &median_dQdx_of_TPCObject_muon_candidate_u_plane);
  theBackgroundTree->SetBranchAddress("median_dQdx_of_TPCObject_muon_candidate_v_plane", &median_dQdx_of_TPCObject_muon_candidate_v_plane);
  theBackgroundTree->SetBranchAddress("median_dQdx_of_TPCObject_muon_candidate_y_plane", &median_dQdx_of_TPCObject_muon_candidate_y_plane);
  theBackgroundTree->SetBranchAddress("median_dQdx_of_TPCObject_muon_candidate_sum", &median_dQdx_of_TPCObject_muon_candidate_sum);
  theBackgroundTree->SetBranchAddress("u_plane_hit_sum", &u_plane_hit_sum);
  theBackgroundTree->SetBranchAddress("v_plane_hit_sum", &v_plane_hit_sum);
  theBackgroundTree->SetBranchAddress("y_plane_hit_sum", &y_plane_hit_sum);
  theBackgroundTree->SetBranchAddress("hit_sum", &hit_sum);
  theBackgroundTree->SetBranchAddress("vicinity_hit_sum_difference", &vicinity_hit_sum_difference);
  theBackgroundTree->SetBranchAddress("u_plane_vicinity_hit_sum_difference", &u_plane_vicinity_hit_sum_difference);
  theBackgroundTree->SetBranchAddress("v_plane_vicinity_hit_sum_difference", &v_plane_vicinity_hit_sum_difference);
  theBackgroundTree->SetBranchAddress("y_plane_vicinity_hit_sum_difference", &y_plane_vicinity_hit_sum_difference);
  theBackgroundTree->SetBranchAddress("total_number_of_slice_hits_in_vicinity_ADCs", &total_number_of_slice_hits_in_vicinity_ADCs);
  theBackgroundTree->SetBranchAddress("u_plane_number_of_slice_hits_in_vicinity_ADCs", &u_plane_number_of_slice_hits_in_vicinity_ADCs);
  theBackgroundTree->SetBranchAddress("v_plane_number_of_slice_hits_in_vicinity_ADCs", &v_plane_number_of_slice_hits_in_vicinity_ADCs);
  theBackgroundTree->SetBranchAddress("y_plane_number_of_slice_hits_in_vicinity_ADCs", &y_plane_number_of_slice_hits_in_vicinity_ADCs);
  theBackgroundTree->SetBranchAddress("num_of_top_pandora_crossings", &num_of_top_pandora_crossings);
  theBackgroundTree->SetBranchAddress("num_of_bottom_pandora_crossings", &num_of_bottom_pandora_crossings);
  theBackgroundTree->SetBranchAddress("num_of_front_pandora_crossings", &num_of_front_pandora_crossings);
  theBackgroundTree->SetBranchAddress("num_of_back_pandora_crossings", &num_of_back_pandora_crossings);
  theBackgroundTree->SetBranchAddress("number_of_tracks_in_TPCObject", &number_of_tracks_in_TPCObject);
  theBackgroundTree->SetBranchAddress("num_of_tracks_originating_from_vertex", &num_of_tracks_originating_from_vertex);
  theBackgroundTree->SetBranchAddress("muon_candidate_length", &muon_candidate_length);
  theBackgroundTree->SetBranchAddress("sum_of_TPCObject_track_lengths", &sum_of_TPCObject_track_lengths);
  theBackgroundTree->SetBranchAddress("total_length_of_tracks_originating_from_vertex", &total_length_of_tracks_originating_from_vertex);
  theBackgroundTree->SetBranchAddress("rootino_fix_mcweight", &rootino_fix_mcweight);
  theBackgroundTree->SetBranchAddress("central_value_mcweight", &central_value_mcweight);
  theBackgroundTree->SetBranchAddress("other_universe_mcweights", &other_universe_mcweights);
  theBackgroundTree->SetBranchAddress("axialff_mcweights", &axialff_mcweights);
  theBackgroundTree->SetBranchAddress("rpaccqe_mcweights", &rpaccqe_mcweights);
  theBackgroundTree->SetBranchAddress("xsecshape_mcweights", &xsecshape_mcweights);
  
  output_ext = 0;
  
  std::cout << "Number of entries in 'theBackgroundTree' = " << theBackgroundTree->GetEntries() << "." << std::endl;
  
  for (Long64_t ievt = 0; ievt < theBackgroundTree->GetEntries(); ievt++){
    theBackgroundTree->GetEntry(ievt);
    Float_t value = reader->EvaluateMVA("BDTG method");

    output_run                                                       = run;
    output_subrun                                                    = subrun;
    output_event                                                     = event;
    output_score                                                     = value;
    output_truth_neutrino_energy                                     = truth_neutrino_energy;
    output_num_muminus_tracks                                        = num_muminus_tracks;
    output_num_muplus_tracks                                         = num_muplus_tracks;
    output_num_piplus_tracks                                         = num_piplus_tracks;
    output_num_piminus_tracks                                        = num_piminus_tracks;
    output_num_pi0_tracks                                            = num_pi0_tracks;
    output_num_proton_tracks                                         = num_proton_tracks;
    output_num_electron_tracks                                       = num_electron_tracks;
    output_num_positron_tracks                                       = num_positron_tracks;
    output_num_photon_tracks                                         = num_photon_tracks;
    output_spline_fix_mcweight                                       = spline_fix_mcweight;
    output_special_kdar_weight                                       = special_kdar_weight;
    output_sampleweight                                              = sampleweight;
    output_treeweight                                                = treeweight;
    output_ext                                                       = ext;
    output_NC_channel                                                = NC_channel;
    output_truth_reco_vtx_distance                                   = truth_reco_vtx_distance;
    output_truth_theta_angle_for_weighting                           = truth_theta_angle_for_weighting;
    output_neutrino_energy                                           = neutrino_energy;
    output_neutrino_cos_angle                                        = neutrino_cos_angle;
    output_reco_vtx_x                                                = reco_vtx_x;
    output_reco_vtx_y                                                = reco_vtx_y;
    output_reco_vtx_z                                                = reco_vtx_z;
    output_muon_x_momentum_normalized                                = muon_x_momentum_normalized;
    output_muon_y_momentum_normalized                                = muon_y_momentum_normalized;
    output_muon_z_momentum_normalized                                = muon_z_momentum_normalized;
    output_x_component_of_muon_momentum                              = x_component_of_muon_momentum;
    output_y_component_of_muon_momentum                              = y_component_of_muon_momentum;
    output_z_component_of_muon_momentum                              = z_component_of_muon_momentum;
    output_num_beam_flashes                                          = num_beam_flashes;
    output_flash_PEs                                                 = flash_PEs;
    output_flash_z_difference                                        = flash_z_difference;
    output_flash_y_difference                                        = flash_y_difference;
    output_truncated_dQdx_of_TPCObject_muon_candidate_u_plane        = truncated_dQdx_of_TPCObject_muon_candidate_u_plane;
    output_truncated_dQdx_of_TPCObject_muon_candidate_v_plane        = truncated_dQdx_of_TPCObject_muon_candidate_v_plane;
    output_truncated_dQdx_of_TPCObject_muon_candidate_y_plane        = truncated_dQdx_of_TPCObject_muon_candidate_y_plane;
    output_truncated_dQdx_of_TPCObject_muon_candidate_sum            = truncated_dQdx_of_TPCObject_muon_candidate_sum;
    output_median_dQdx_of_TPCObject_muon_candidate_u_plane           = median_dQdx_of_TPCObject_muon_candidate_u_plane;
    output_median_dQdx_of_TPCObject_muon_candidate_v_plane           = median_dQdx_of_TPCObject_muon_candidate_v_plane;
    output_median_dQdx_of_TPCObject_muon_candidate_y_plane           = median_dQdx_of_TPCObject_muon_candidate_y_plane;
    output_median_dQdx_of_TPCObject_muon_candidate_sum               = median_dQdx_of_TPCObject_muon_candidate_sum;
    output_u_plane_hit_sum                                           = u_plane_hit_sum;
    output_v_plane_hit_sum                                           = v_plane_hit_sum;
    output_y_plane_hit_sum                                           = y_plane_hit_sum;
    output_hit_sum                                                   = hit_sum;
    output_u_plane_vicinity_hit_sum_difference                       = u_plane_vicinity_hit_sum_difference;
    output_v_plane_vicinity_hit_sum_difference                       = v_plane_vicinity_hit_sum_difference;
    output_y_plane_vicinity_hit_sum_difference                       = y_plane_vicinity_hit_sum_difference;
    output_vicinity_hit_sum_difference                               = vicinity_hit_sum_difference;
    output_u_plane_number_of_slice_hits_in_vicinity_ADCs             = u_plane_number_of_slice_hits_in_vicinity_ADCs;
    output_v_plane_number_of_slice_hits_in_vicinity_ADCs             = v_plane_number_of_slice_hits_in_vicinity_ADCs;
    output_y_plane_number_of_slice_hits_in_vicinity_ADCs             = y_plane_number_of_slice_hits_in_vicinity_ADCs;
    output_total_number_of_slice_hits_in_vicinity_ADCs               = total_number_of_slice_hits_in_vicinity_ADCs;
    output_num_of_top_pandora_crossings                              = num_of_top_pandora_crossings;
    output_num_of_bottom_pandora_crossings                           = num_of_bottom_pandora_crossings;
    output_num_of_front_pandora_crossings                            = num_of_front_pandora_crossings;
    output_num_of_back_pandora_crossings                             = num_of_back_pandora_crossings;
    output_number_of_tracks_in_TPCObject                             = number_of_tracks_in_TPCObject;
    output_num_of_tracks_originating_from_vertex                     = num_of_tracks_originating_from_vertex;
    output_muon_candidate_length                                     = muon_candidate_length;
    output_sum_of_TPCObject_track_lengths                            = sum_of_TPCObject_track_lengths;
    output_total_length_of_tracks_originating_from_vertex            = total_length_of_tracks_originating_from_vertex;

    output_rootino_fix_mcweight                                      = rootino_fix_mcweight;
    output_central_value_mcweight                                    = central_value_mcweight;

    for ( size_t i = 0; i < 100; i++ ) {
      
      output_other_universe_mcweights[i] = other_universe_mcweights[i];

    }

    for ( size_t i = 0; i < 2; i++ ) {

      output_axialff_mcweights[i] = axialff_mcweights[i];

    }

    for ( size_t i = 0; i < 2; i++ ) {

      output_rpaccqe_mcweights[i] = rpaccqe_mcweights[i];

    }


    for ( size_t i = 0; i < 2; i++ ) {

      output_xsecshape_mcweights[i] = xsecshape_mcweights[i];

    }
    
    if ( output_treeweight > 1.1 )
      output_ext = 1;
    
    output_bkgd_tree->Fill();
    
  }
    
  std::cout << "--- End of event loop: "; 
  
  delete reader;
    
  std::cout << "==> TMVAClassificationApplication is done!" << endl << std::endl;
  
  output_file->Write();
  output_file->Close();

} 
