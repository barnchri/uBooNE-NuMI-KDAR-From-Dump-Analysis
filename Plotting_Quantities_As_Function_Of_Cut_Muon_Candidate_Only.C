#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TCanvas.h"
#include "TH1.h"

void Plotting_Quantities_As_Function_Of_Cut_Muon_Candidate_Only() {
  
  TFile* file = new TFile("fooling_around_with_legend.root", "RECREATE");

  float score_cut_value = 0.7;

  // Put the new flux stuff up here.
  TFile* input_1D_flux_file  = new TFile("/home/barnchri/Experimenting_With_Different_Flux_Universes/output_ratio_file_1D.root");
  TH1D*  input_1D_flux_hist  = (TH1D*)input_1D_flux_file->Get("numu_CV_AV_TPC_5MeV_bin");

  TFile* input_2D_flux_file  = new TFile("/home/barnchri/Experimenting_With_Different_Flux_Universes/output_ratio_file_2D.root");
  
  int num_flux_universes = 100;
      
  float data_normalized_events;
  float data_score;
  float data_muon_candidate_length;

  float signal_spline_fix_mcweight;
  float signal_rootino_fix_mcweight;
  float signal_central_value_mcweight;
  float signal_special_kdar_weight;
  float signal_normalized_events;
  float signal_score;
  float signal_muon_candidate_length;
  
  float background_spline_fix_mcweight;
  float background_rootino_fix_mcweight;
  float background_central_value_mcweight;
  float background_special_kdar_weight;
  float background_flux_cv_weight; // Calculated using the input histograms.
  float background_flux_other_universe_weight;
  float background_truth_neutrino_energy; // In MeV!!!
  double background_truth_angle; 
  float background_other_universe_mcweights[100];
  float background_axialff_mcweights[2];
  float background_rpaccqe_mcweights[2];
  float	background_xsecshape_mcweights[2];
  int   background_ext;
  float background_numi_ext_normalized_events;
  float background_numi_nu_normalized_events;
  float background_muon_candidate_length;
  float background_score;

  float systematic_normalization_factor;

  // Detector systematic samples.
  float rayleigh_length_normalized_events;
  float rayleigh_length_spline_fix_mcweight;
  float rayleigh_length_rootino_fix_mcweight;
  float rayleigh_length_central_value_mcweight;
  float rayleigh_length_special_kdar_weight;
  int   rayleigh_length_ext;
  float rayleigh_length_numi_ext_normalized_events;
  float rayleigh_length_numi_nu_normalized_events;
  float rayleigh_length_muon_candidate_length;
  float rayleigh_length_score;

  float	light_yield_normalized_events;
  float light_yield_spline_fix_mcweight;
  float light_yield_rootino_fix_mcweight;
  float light_yield_central_value_mcweight;
  float light_yield_special_kdar_weight;
  int   light_yield_ext;
  float light_yield_numi_ext_normalized_events;
  float light_yield_numi_nu_normalized_events;
  float light_yield_muon_candidate_length;
  float light_yield_score;

  float	attenuation_normalized_events;
  float attenuation_spline_fix_mcweight;
  float attenuation_rootino_fix_mcweight;
  float attenuation_central_value_mcweight;
  float attenuation_special_kdar_weight;
  int   attenuation_ext;
  float attenuation_numi_ext_normalized_events;
  float attenuation_numi_nu_normalized_events;
  float attenuation_muon_candidate_length;
  float attenuation_score;

  float	sce_normalized_events;
  float sce_spline_fix_mcweight;
  float sce_rootino_fix_mcweight;
  float sce_central_value_mcweight;
  float sce_special_kdar_weight;
  int   sce_ext;
  float sce_numi_ext_normalized_events;
  float sce_numi_nu_normalized_events;
  float sce_muon_candidate_length;
  float sce_score;

  float	recombination_normalized_events;
  float recombination_spline_fix_mcweight;
  float recombination_rootino_fix_mcweight;
  float recombination_central_value_mcweight;
  float recombination_special_kdar_weight;
  int   recombination_ext;
  float recombination_numi_ext_normalized_events;
  float recombination_numi_nu_normalized_events;
  float recombination_muon_candidate_length;
  float recombination_score;

  float	scalex_normalized_events;
  float scalex_spline_fix_mcweight;
  float scalex_rootino_fix_mcweight;
  float scalex_central_value_mcweight;
  float scalex_special_kdar_weight;
  int   scalex_ext;
  float scalex_numi_ext_normalized_events;
  float scalex_numi_nu_normalized_events;
  float scalex_muon_candidate_length;
  float scalex_score;

  float	scaleyz_normalized_events;
  float scaleyz_spline_fix_mcweight;
  float scaleyz_rootino_fix_mcweight;
  float scaleyz_central_value_mcweight;
  float scaleyz_special_kdar_weight;
  int   scaleyz_ext;
  float scaleyz_numi_ext_normalized_events;
  float scaleyz_numi_nu_normalized_events;
  float scaleyz_muon_candidate_length;
  float scaleyz_score;

  float	scaleangleyz_normalized_events;
  float scaleangleyz_spline_fix_mcweight;
  float scaleangleyz_rootino_fix_mcweight;
  float scaleangleyz_central_value_mcweight;
  float scaleangleyz_special_kdar_weight;
  int   scaleangleyz_ext;
  float scaleangleyz_numi_ext_normalized_events;
  float scaleangleyz_numi_nu_normalized_events;
  float scaleangleyz_muon_candidate_length;
  float scaleangleyz_score;

  float	scaleanglexz_normalized_events;
  float scaleanglexz_spline_fix_mcweight;
  float scaleanglexz_rootino_fix_mcweight;
  float scaleanglexz_central_value_mcweight;
  float scaleanglexz_special_kdar_weight;
  int   scaleanglexz_ext;
  float scaleanglexz_numi_ext_normalized_events;
  float scaleanglexz_numi_nu_normalized_events;
  float scaleanglexz_muon_candidate_length;
  float scaleanglexz_score;

  float	scalededx_normalized_events;
  float scalededx_spline_fix_mcweight;
  float scalededx_rootino_fix_mcweight;
  float scalededx_central_value_mcweight;
  float scalededx_special_kdar_weight;
  int   scalededx_ext;
  float scalededx_numi_ext_normalized_events;
  float scalededx_numi_nu_normalized_events;
  float scalededx_muon_candidate_length;
  float scalededx_score;

  // XSec Systematic Samples.
  std::vector< float > other_xsec_universes_normalized_events;
  std::vector< float > axialff_universes_normalized_events;
  std::vector< float > rpaccqe_universes_normalized_events;
  std::vector< float > xsecshape_universes_normalized_events;
  
  other_xsec_universes_normalized_events.resize( 100, 0. );
  axialff_universes_normalized_events.resize( 2, 0. );
  rpaccqe_universes_normalized_events.resize( 2, 0. );
  xsecshape_universes_normalized_events.resize( 2, 0. );

    // Flux Systematic Samples.                                                                                                                                                                            
  float flux_cv_normalized_events;
  std::vector< float > flux_other_universes_normalized_events;

  flux_other_universes_normalized_events.resize( num_flux_universes, 0. );

  // Declare the two background muon candidate histograms.
  TH1F*    muon_candidate_length_ext_background                      = new TH1F("muon_candidate_length_ext_background", "NuMI EXT Muon Candidate Length", 20, 0., 40.);
  TH1F*    muon_candidate_length_ext_background_before_normalization = new TH1F("muon_candidate_length_ext_background_before_normalization", "NuMI EXT Muon Candidate Length", 20, 0., 40.);
  TH1F*    muon_candidate_length_nu_background                       = new TH1F("muon_candidate_length_nu_background", "NuMI Neutrino Background Muon Candidate Length", 20, 0., 40.);
  TH1F*    muon_candidate_length_nu_background_without_weights       = new TH1F("muon_candidate_length_nu_background_without_weights", "NuMI Neutrino Background Muon Candidate Length", 20, 0., 40.);
  
  // Declare the xsec systematic histograms.
  int num_other_xsec_universes = 100;
  
  std::vector< TH1F* > other_xsec_universes;
  other_xsec_universes.resize( num_other_xsec_universes );

  for ( int other_universe_iter = 0; other_universe_iter < num_other_xsec_universes; other_universe_iter++ ) {

    TH1F* muon_candidate_length_nu_background_other_xsec_universe = new TH1F(Form("other_xsec_universe_%d", other_universe_iter), Form("Other XSec Universe #%d", other_universe_iter), 20, 0., 40. );

    other_xsec_universes.at( other_universe_iter ) = muon_candidate_length_nu_background_other_xsec_universe;

    delete muon_candidate_length_nu_background_other_xsec_universe;

  }

  int num_axialff_corrections = 2;

  std::vector< TH1F* > axialff_universes;
  axialff_universes.resize( num_axialff_corrections );

  for ( int axialff_iter = 0; axialff_iter < num_axialff_corrections; axialff_iter++ ) {

    TH1F* muon_candidate_length_nu_background_axialff_universe = new TH1F(Form("axialff_universe_%d", axialff_iter), Form("Axialff Universe #%d", axialff_iter), 20, 0., 40. );

    axialff_universes.at( axialff_iter ) = muon_candidate_length_nu_background_axialff_universe;

    delete muon_candidate_length_nu_background_axialff_universe;

  }

  int num_rpaccqe_corrections = 2;

  std::vector< TH1F* > rpaccqe_universes;
  rpaccqe_universes.resize( num_rpaccqe_corrections );

  for ( int rpaccqe_iter = 0; rpaccqe_iter < num_rpaccqe_corrections; rpaccqe_iter++ ) {

    TH1F* muon_candidate_length_nu_background_rpaccqe_universe = new TH1F(Form("rpaccqe_universe_%d", rpaccqe_iter), Form("RPACCQE Universe #%d", rpaccqe_iter), 20, 0., 40. );

    rpaccqe_universes.at( rpaccqe_iter ) = muon_candidate_length_nu_background_rpaccqe_universe;

    delete muon_candidate_length_nu_background_rpaccqe_universe;

  }
  
  int num_xsecshape_corrections = 2;

  std::vector< TH1F* > xsecshape_universes;
  xsecshape_universes.resize( num_xsecshape_corrections );

  for ( int xsecshape_iter = 0; xsecshape_iter < num_xsecshape_corrections; xsecshape_iter++ ) {

    TH1F* muon_candidate_length_nu_background_xsecshape_universe = new TH1F(Form("xsecshape_universe_%d", xsecshape_iter), Form("XSecShape Universe #%d", xsecshape_iter), 20, 0., 40. );

    xsecshape_universes.at( xsecshape_iter ) = muon_candidate_length_nu_background_xsecshape_universe;

    delete muon_candidate_length_nu_background_xsecshape_universe;

  }

  // Flux Systematic samples.
  TH1F*    muon_candidate_length_nu_background_flux_cv         = new TH1F("muon_candidate_length_nu_background_flux_cv", "NuMI Neutrino Background Muon Candidate Length", 20, 0., 40.);
  
  TH1F*   muon_candidate_length_nu_background_CV_flux_sample   = new TH1F("muon_candidate_length_nu_background_CV_flux_sample", "NuMI Neutrino Background Muon Candidate Length", 20, 0., 40.);

  std::vector< TH1F* > flux_universes;
  flux_universes.resize( num_flux_universes );

  for ( int flux_iter = 0; flux_iter < num_flux_universes; flux_iter++ ) {

    TH1F* muon_candidate_length_nu_background_flux_universe = new TH1F(Form("flux_universe_%d", flux_iter), Form("Flux Universe #%d", flux_iter), 20, 0., 40.);

    flux_universes.at( flux_iter ) = muon_candidate_length_nu_background_flux_universe;

    delete muon_candidate_length_nu_background_flux_universe;

  }
  
  // Declare the detector systematic histograms.
  TH1F*    muon_candidate_length_nu_background_rayleigh_length = new TH1F("muon_candidate_length_nu_background_rayleigh_length", "NuMI Neutrino Background Muon Candidate Length", 20, 0., 40.);
  TH1F*    muon_candidate_length_nu_background_light_yield     = new TH1F("muon_candidate_length_nu_background_light_yield", "NuMI Neutrino Background Muon Candidate Length", 20, 0., 40.);
  TH1F*    muon_candidate_length_nu_background_attenuation     = new TH1F("muon_candidate_length_nu_background_attenuation", "NuMI Neutrino Background Muon Candidate Length", 20, 0., 40.);
  TH1F*    muon_candidate_length_nu_background_sce             = new TH1F("muon_candidate_length_nu_background_sce", "NuMI Neutrino Background Muon Candidate Length", 20, 0., 40.);
  TH1F*    muon_candidate_length_nu_background_recombination   = new TH1F("muon_candidate_length_nu_background_recombination", "NuMI Neutrino Background Muon Candidate Length", 20, 0., 40.);
  TH1F*    muon_candidate_length_nu_background_scalex          = new TH1F("muon_candidate_length_nu_background_scalex", "NuMI Neutrino Background Muon Candidate Length", 20, 0., 40.);
  TH1F*    muon_candidate_length_nu_background_scaleyz         = new TH1F("muon_candidate_length_nu_background_scaleyz", "NuMI Neutrino Background Muon Candidate Length", 20, 0., 40.);
  TH1F*    muon_candidate_length_nu_background_scaleangleyz    = new TH1F("muon_candidate_length_nu_background_scaleangleyz", "NuMI Neutrino Background Muon Candidate Length", 20, 0., 40.);
  TH1F*    muon_candidate_length_nu_background_scaleanglexz    = new TH1F("muon_candidate_length_nu_background_scaleanglexz", "NuMI Neutrino Background Muon Candidate Length", 20, 0., 40.);
  TH1F*    muon_candidate_length_nu_background_scalededx       = new TH1F("muon_candidate_length_nu_background_scalededx", "NuMI Neutrino Background Muon Candidate Length", 20, 0., 40.);
  TH1F*    muon_candidate_length_signal                        = new TH1F("muon_candidate_length_signal", "Signal Muon Candidate Length", 20, 0., 40.);
  TH1F*    muon_candidate_length_data                          = new TH1F("muon_candidate_length_data", "Data Muon Candidate Length", 20, 0., 40.);
  TH1F*    muon_candidate_length_data_before_subtraction       = new TH1F("muon_candidate_length_data_before_subtraction", "Data Muon Candidate Length", 20, 0., 40.);
  THStack* muon_candidate_length_stack                         = new THStack("muon_candidate_length_stack", "");
  
  TH1F*    score_ext_background                                = new TH1F("score_ext_background", "NuMI EXT Score", 100, -1.0, 1.0);
  TH1F*    score_nu_background                                 = new TH1F("score_nu_background", "NuMI Neutrino Background Score", 100, -1.0, 1.0);
  TH1F*    score_signal                                        = new TH1F("score_signal", "Signal Score", 100, -1.0, 1.0);
  TH1F*    score_data                                          = new TH1F("score_data", "Data Score", 100, -1.0, 1.0);
  THStack* score_stack                                         = new THStack("score_stack", "");

  // Load in the data tree.
  TChain* data  = new TChain("signal_tree");
  data->Add("/home/barnchri/Finishing_Off_Initial_Systematic_Studies/output_file_beamdata_all_stats.root");
  data->SetBranchAddress("score", &data_score);
  data->SetBranchAddress("muon_candidate_length", &data_muon_candidate_length);

  // Load in the signal tree.
  TChain* signal = new TChain("signal_tree");
  signal->Add("/home/barnchri/Finishing_Off_Initial_Systematic_Studies/output_file_cv_testing.root");
  signal->SetBranchAddress("score", &signal_score);
  signal->SetBranchAddress("spline_fix_mcweight", &signal_spline_fix_mcweight);
  signal->SetBranchAddress("rootino_fix_mcweight", &signal_rootino_fix_mcweight);
  signal->SetBranchAddress("central_value_mcweight", &signal_central_value_mcweight);
  signal->SetBranchAddress("special_kdar_weight", &signal_special_kdar_weight);
  signal->SetBranchAddress("muon_candidate_length", &signal_muon_candidate_length);

  // Load in the background tree.
  TChain* background = new TChain("background_tree");
  background->Add("/home/barnchri/Finishing_Off_Initial_Systematic_Studies/output_file_cv_testing.root");
  background->SetBranchAddress("score", &background_score);
  background->SetBranchAddress("spline_fix_mcweight", &background_spline_fix_mcweight);
  background->SetBranchAddress("rootino_fix_mcweight", &background_rootino_fix_mcweight);
  background->SetBranchAddress("central_value_mcweight", &background_central_value_mcweight);
  background->SetBranchAddress("special_kdar_weight", &background_special_kdar_weight);
  background->SetBranchAddress("truth_neutrino_energy", &background_truth_neutrino_energy);
  background->SetBranchAddress("truth_theta_angle_for_weighting", &background_truth_angle);
  background->SetBranchAddress("other_universe_mcweights", &background_other_universe_mcweights);
  background->SetBranchAddress("axialff_mcweights", &background_axialff_mcweights);
  background->SetBranchAddress("rpaccqe_mcweights", &background_rpaccqe_mcweights);
  background->SetBranchAddress("xsecshape_mcweights", &background_xsecshape_mcweights);
  background->SetBranchAddress("ext", &background_ext);
  background->SetBranchAddress("muon_candidate_length", &background_muon_candidate_length);

  // Include information for the detector variation samples.

  // rayleigh length
  TChain* rayleigh_length = new TChain("background_tree");
  rayleigh_length->Add("/home/barnchri/Finishing_Off_Initial_Systematic_Studies/output_file_rayleigh_length_testing.root");
  rayleigh_length->SetBranchAddress("score", &rayleigh_length_score);
  rayleigh_length->SetBranchAddress("spline_fix_mcweight", &rayleigh_length_spline_fix_mcweight);
  rayleigh_length->SetBranchAddress("rootino_fix_mcweight", &rayleigh_length_rootino_fix_mcweight);
  rayleigh_length->SetBranchAddress("central_value_mcweight", &rayleigh_length_central_value_mcweight);
  rayleigh_length->SetBranchAddress("special_kdar_weight", &rayleigh_length_special_kdar_weight);
  rayleigh_length->SetBranchAddress("ext", &rayleigh_length_ext);
  rayleigh_length->SetBranchAddress("muon_candidate_length", &rayleigh_length_muon_candidate_length);

  // light yield
  TChain* light_yield = new TChain("background_tree");
  light_yield->Add("/home/barnchri/Finishing_Off_Initial_Systematic_Studies/output_file_light_yield_testing.root");
  light_yield->SetBranchAddress("score", &light_yield_score);
  light_yield->SetBranchAddress("spline_fix_mcweight", &light_yield_spline_fix_mcweight);
  light_yield->SetBranchAddress("rootino_fix_mcweight", &light_yield_rootino_fix_mcweight);
  light_yield->SetBranchAddress("central_value_mcweight", &light_yield_central_value_mcweight);
  light_yield->SetBranchAddress("special_kdar_weight", &light_yield_special_kdar_weight);
  light_yield->SetBranchAddress("ext", &light_yield_ext);
  light_yield->SetBranchAddress("muon_candidate_length", &light_yield_muon_candidate_length);

  // attenuation
  TChain* attenuation = new TChain("background_tree");
  attenuation->Add("/home/barnchri/Finishing_Off_Initial_Systematic_Studies/output_file_attenuation_testing.root");
  attenuation->SetBranchAddress("score", &attenuation_score);
  attenuation->SetBranchAddress("spline_fix_mcweight", &attenuation_spline_fix_mcweight);
  attenuation->SetBranchAddress("rootino_fix_mcweight", &attenuation_rootino_fix_mcweight);
  attenuation->SetBranchAddress("central_value_mcweight", &attenuation_central_value_mcweight);
  attenuation->SetBranchAddress("special_kdar_weight", &attenuation_special_kdar_weight);
  attenuation->SetBranchAddress("ext", &attenuation_ext);
  attenuation->SetBranchAddress("muon_candidate_length", &attenuation_muon_candidate_length);

  // sce
  TChain* sce = new TChain("background_tree");
  sce->Add("/home/barnchri/Finishing_Off_Initial_Systematic_Studies/output_file_sce_testing.root");
  sce->SetBranchAddress("score", &sce_score);
  sce->SetBranchAddress("spline_fix_mcweight", &sce_spline_fix_mcweight);
  sce->SetBranchAddress("rootino_fix_mcweight", &sce_rootino_fix_mcweight);
  sce->SetBranchAddress("central_value_mcweight", &sce_central_value_mcweight);
  sce->SetBranchAddress("special_kdar_weight", &sce_special_kdar_weight);
  sce->SetBranchAddress("ext", &sce_ext);
  sce->SetBranchAddress("muon_candidate_length", &sce_muon_candidate_length);

  // recombination
  TChain* recombination = new TChain("background_tree");
  recombination->Add("/home/barnchri/Finishing_Off_Initial_Systematic_Studies/output_file_recombination_testing.root");
  recombination->SetBranchAddress("score", &recombination_score);
  recombination->SetBranchAddress("spline_fix_mcweight", &recombination_spline_fix_mcweight);
  recombination->SetBranchAddress("rootino_fix_mcweight", &recombination_rootino_fix_mcweight);
  recombination->SetBranchAddress("central_value_mcweight", &recombination_central_value_mcweight);
  recombination->SetBranchAddress("special_kdar_weight", &recombination_special_kdar_weight);
  recombination->SetBranchAddress("ext", &recombination_ext);
  recombination->SetBranchAddress("muon_candidate_length", &recombination_muon_candidate_length);

  // scalex
  TChain* scalex = new TChain("background_tree");
  scalex->Add("/home/barnchri/Finishing_Off_Initial_Systematic_Studies/output_file_scalex_testing.root");
  scalex->SetBranchAddress("score", &scalex_score);
  scalex->SetBranchAddress("spline_fix_mcweight", &scalex_spline_fix_mcweight);
  scalex->SetBranchAddress("rootino_fix_mcweight", &scalex_rootino_fix_mcweight);
  scalex->SetBranchAddress("central_value_mcweight", &scalex_central_value_mcweight);
  scalex->SetBranchAddress("special_kdar_weight", &scalex_special_kdar_weight);
  scalex->SetBranchAddress("ext", &scalex_ext);
  scalex->SetBranchAddress("muon_candidate_length", &scalex_muon_candidate_length);

  // scaleyz
  TChain* scaleyz = new TChain("background_tree");
  scaleyz->Add("/home/barnchri/Finishing_Off_Initial_Systematic_Studies/output_file_scaleyz_testing.root");
  scaleyz->SetBranchAddress("score", &scaleyz_score);
  scaleyz->SetBranchAddress("spline_fix_mcweight", &scaleyz_spline_fix_mcweight);
  scaleyz->SetBranchAddress("rootino_fix_mcweight", &scaleyz_rootino_fix_mcweight);
  scaleyz->SetBranchAddress("central_value_mcweight", &scaleyz_central_value_mcweight);
  scaleyz->SetBranchAddress("special_kdar_weight", &scaleyz_special_kdar_weight);
  scaleyz->SetBranchAddress("ext", &scaleyz_ext);
  scaleyz->SetBranchAddress("muon_candidate_length", &scaleyz_muon_candidate_length);

  // scaleangleyz
  TChain* scaleangleyz = new TChain("background_tree");
  scaleangleyz->Add("/home/barnchri/Finishing_Off_Initial_Systematic_Studies/output_file_scaleangleyz_testing.root");
  scaleangleyz->SetBranchAddress("score", &scaleangleyz_score);
  scaleangleyz->SetBranchAddress("spline_fix_mcweight", &scaleangleyz_spline_fix_mcweight);
  scaleangleyz->SetBranchAddress("rootino_fix_mcweight", &scaleangleyz_rootino_fix_mcweight);
  scaleangleyz->SetBranchAddress("central_value_mcweight", &scaleangleyz_central_value_mcweight);
  scaleangleyz->SetBranchAddress("special_kdar_weight", &scaleangleyz_special_kdar_weight);
  scaleangleyz->SetBranchAddress("ext", &scaleangleyz_ext);
  scaleangleyz->SetBranchAddress("muon_candidate_length", &scaleangleyz_muon_candidate_length);

  // scaleanglexz
  TChain* scaleanglexz = new TChain("background_tree");
  scaleanglexz->Add("/home/barnchri/Finishing_Off_Initial_Systematic_Studies/output_file_scaleanglexz_testing.root");
  scaleanglexz->SetBranchAddress("score", &scaleanglexz_score);
  scaleanglexz->SetBranchAddress("spline_fix_mcweight", &scaleanglexz_spline_fix_mcweight);
  scaleanglexz->SetBranchAddress("rootino_fix_mcweight", &scaleanglexz_rootino_fix_mcweight);
  scaleanglexz->SetBranchAddress("central_value_mcweight", &scaleanglexz_central_value_mcweight);
  scaleanglexz->SetBranchAddress("special_kdar_weight", &scaleanglexz_special_kdar_weight);
  scaleanglexz->SetBranchAddress("ext", &scaleanglexz_ext);
  scaleanglexz->SetBranchAddress("muon_candidate_length", &scaleanglexz_muon_candidate_length);

  // scalededx
  TChain* scalededx = new TChain("background_tree");
  scalededx->Add("/home/barnchri/Finishing_Off_Initial_Systematic_Studies/output_file_scalededx_testing.root");
  scalededx->SetBranchAddress("score", &scalededx_score);
  scalededx->SetBranchAddress("spline_fix_mcweight", &scalededx_spline_fix_mcweight);
  scalededx->SetBranchAddress("rootino_fix_mcweight", &scalededx_rootino_fix_mcweight);
  scalededx->SetBranchAddress("central_value_mcweight", &scalededx_central_value_mcweight);
  scalededx->SetBranchAddress("special_kdar_weight", &scalededx_special_kdar_weight);
  scalededx->SetBranchAddress("ext", &scalededx_ext);
  scalededx->SetBranchAddress("muon_candidate_length", &scalededx_muon_candidate_length);

  // Loop over the data entries to fill the histograms.
  int num_data_entries                        = data->GetEntries();
  float num_data_events_above_score_threshold = 0.;

  for ( int data_iter = 0; data_iter < num_data_entries; data_iter++ ) {

    data->GetEntry( data_iter );

    if ( data_score > score_cut_value ) {

      muon_candidate_length_data->Fill( data_muon_candidate_length);
      muon_candidate_length_data_before_subtraction->Fill( data_muon_candidate_length );
      score_data->Fill( data_score );

      num_data_events_above_score_threshold += 1.0;
      
    }

  }

  std::cout << "The number of data events with a score > " << score_cut_value << " = " << num_data_events_above_score_threshold << "." << std::endl;

  // Loop over the signal entries to fill the histograms.
  int num_signal_entries = signal->GetEntries();

  float   num_unweighted_signal_events_above_score_threshold = 0;
  float   num_weighted_signal_events_above_score_threshold   = 0;

  for ( int signal_iter = 0; signal_iter < num_signal_entries; signal_iter++ ) {

    signal->GetEntry( signal_iter );

    if ( signal_score > score_cut_value ) {

      muon_candidate_length_signal->Fill( signal_muon_candidate_length, signal_spline_fix_mcweight * signal_rootino_fix_mcweight * signal_central_value_mcweight * signal_special_kdar_weight );
      score_signal->Fill( signal_score, signal_spline_fix_mcweight * signal_rootino_fix_mcweight * signal_central_value_mcweight * signal_special_kdar_weight );
      
      num_unweighted_signal_events_above_score_threshold++;
      num_weighted_signal_events_above_score_threshold += signal_spline_fix_mcweight * signal_rootino_fix_mcweight * signal_central_value_mcweight * signal_special_kdar_weight;
      
    }

  }

  std::cout << "The number of unweighted KDAR signal events with a score > " << score_cut_value << " = " << num_unweighted_signal_events_above_score_threshold << "." << std::endl;

  // Loop over the background entries to fill the histograms..
  // Fill the flux and xsec histograms in the same loop.
  int    num_background_entries                                                 = background->GetEntries();
  float num_unweighted_nu_events_above_score_threshold                          = 0.;
  float num_weighted_nu_events_above_score_threshold                            = 0.;
  float num_unweighted_ext_events_above_score_threshold                         = 0.;

  // Number of events in the xsec systematic samples.
  std::vector < float > num_weighted_nu_events_above_threshold_other_universes;
  std::vector < float > num_weighted_nu_events_above_threshold_axialff;
  std::vector < float > num_weighted_nu_events_above_threshold_rpaccqe;
  std::vector < float > num_weighted_nu_events_above_threshold_xsecshape;
  
  num_weighted_nu_events_above_threshold_other_universes.resize( 100, 0. );
  num_weighted_nu_events_above_threshold_axialff.resize( 2, 0. );
  num_weighted_nu_events_above_threshold_rpaccqe.resize( 2, 0. );
  num_weighted_nu_events_above_threshold_xsecshape.resize( 2, 0. );

  // Number of events in flux systematic samples.
  float num_flux_cv_systematic_weighted_nu_events                               = 0.;
  std::vector< float > num_flux_other_universes_systematic_weighted_nu_events;
  num_flux_other_universes_systematic_weighted_nu_events.resize( 100, 0. );
  
  // Declare variables for the bin numbers that you will use for the flux systematics.
  int   cv_x_bin_num                                                            = 0;

  int   other_universe_x_bin_num                                                = 0;
  int   other_universe_y_bin_num                                                = 0;
  
  for ( int background_iter = 0; background_iter < num_background_entries; background_iter++ ) {

    background->GetEntry( background_iter );
    
    if ( background_score > score_cut_value ) {

      if ( background_ext == 1 ) {

	muon_candidate_length_ext_background->Fill( background_muon_candidate_length );
	muon_candidate_length_ext_background_before_normalization->Fill( background_muon_candidate_length );
	score_ext_background->Fill( background_score );
	
	num_unweighted_ext_events_above_score_threshold++;

      }
	
      else {

	muon_candidate_length_nu_background->Fill( background_muon_candidate_length, background_spline_fix_mcweight * background_rootino_fix_mcweight * background_central_value_mcweight * background_special_kdar_weight );
	muon_candidate_length_nu_background_without_weights->Fill( background_muon_candidate_length );
	score_nu_background->Fill( background_score, background_spline_fix_mcweight * background_rootino_fix_mcweight * background_central_value_mcweight * background_special_kdar_weight );

	// Other universe xsec systematic.
	for ( size_t other_universe_iter = 0; other_universe_iter < 100; other_universe_iter++ ) {

	  other_xsec_universes.at( other_universe_iter )->Fill( background_muon_candidate_length, background_other_universe_mcweights[ other_universe_iter ] * background_spline_fix_mcweight * background_rootino_fix_mcweight * background_special_kdar_weight );
	  num_weighted_nu_events_above_threshold_other_universes.at( other_universe_iter ) += ( background_other_universe_mcweights[ other_universe_iter ] * background_spline_fix_mcweight * background_rootino_fix_mcweight * background_special_kdar_weight );
	  
	}

	// Axialff xsec systematic.
	for ( size_t axialff_iter = 0; axialff_iter < 2; axialff_iter++ ) {

	  axialff_universes.at( axialff_iter )->Fill( background_muon_candidate_length, background_axialff_mcweights[ axialff_iter ] * background_spline_fix_mcweight * background_rootino_fix_mcweight * background_special_kdar_weight );
	  num_weighted_nu_events_above_threshold_axialff.at( axialff_iter ) += ( background_axialff_mcweights[ axialff_iter ] * background_spline_fix_mcweight * background_rootino_fix_mcweight * background_special_kdar_weight );
	    
	}

	// RPACCQE xsec systematic.
	for ( size_t rpaccqe_iter = 0; rpaccqe_iter < 2; rpaccqe_iter++ ) {

	  rpaccqe_universes.at( rpaccqe_iter )->Fill( background_muon_candidate_length, background_rpaccqe_mcweights[ rpaccqe_iter ] * background_spline_fix_mcweight * background_rootino_fix_mcweight * background_special_kdar_weight );
	  num_weighted_nu_events_above_threshold_rpaccqe.at( rpaccqe_iter ) += ( background_rpaccqe_mcweights[ rpaccqe_iter ] * background_spline_fix_mcweight * background_rootino_fix_mcweight * background_special_kdar_weight );
	    
	}

	// XSecShape xsec systematic.
	for ( size_t xsecshape_iter = 0; xsecshape_iter < 2; xsecshape_iter++ ) {
	  
	  xsecshape_universes.at( xsecshape_iter )->Fill( background_muon_candidate_length, background_xsecshape_mcweights[ xsecshape_iter ] * background_spline_fix_mcweight * background_rootino_fix_mcweight * background_special_kdar_weight );
	  num_weighted_nu_events_above_threshold_xsecshape.at( xsecshape_iter ) += ( background_xsecshape_mcweights[ xsecshape_iter ] * background_spline_fix_mcweight * background_rootino_fix_mcweight * background_special_kdar_weight );
	  
	}

	num_weighted_nu_events_above_score_threshold += ( background_spline_fix_mcweight * background_rootino_fix_mcweight * background_central_value_mcweight * background_special_kdar_weight );
	num_unweighted_nu_events_above_score_threshold++;

	// Find the flux cv correction weight for the background nu event.
	cv_x_bin_num = ( int( ( background_truth_neutrino_energy / 5.0 ) ) + 1 );
	
	background_flux_cv_weight = input_1D_flux_hist->GetBinContent( cv_x_bin_num );

	// Fill the correct histogram with this weight and the others.
	muon_candidate_length_nu_background_flux_cv->Fill( background_muon_candidate_length, background_flux_cv_weight * background_spline_fix_mcweight * background_rootino_fix_mcweight * background_central_value_mcweight * background_special_kdar_weight );

	num_flux_cv_systematic_weighted_nu_events += ( background_flux_cv_weight * background_spline_fix_mcweight * background_rootino_fix_mcweight * background_central_value_mcweight * background_special_kdar_weight );

	// Set the x and y bin numbers for the 2D histogram by hand.
	other_universe_x_bin_num = 0;
	other_universe_y_bin_num = 0;

	if ( ( background_truth_neutrino_energy / 1000. ) < 0.025 )
	  other_universe_x_bin_num = 1;

	else if ( ( background_truth_neutrino_energy / 1000. ) > 0.025 && ( background_truth_neutrino_energy / 1000. ) < 0.03 )
	  other_universe_x_bin_num = 2;

	else if ( ( background_truth_neutrino_energy / 1000. ) > 0.03 && ( background_truth_neutrino_energy / 1000. ) < 0.235 )
	  other_universe_x_bin_num = 3;

	else if ( ( background_truth_neutrino_energy / 1000. ) > 0.235 && ( background_truth_neutrino_energy / 1000. ) < 0.24 )
	  other_universe_x_bin_num = 4;

	else if ( ( background_truth_neutrino_energy / 1000. ) > 0.24 && ( background_truth_neutrino_energy / 1000. ) < 0.50 )
	  other_universe_x_bin_num = 5;

	else if ( ( background_truth_neutrino_energy / 1000. ) > 0.50 && ( background_truth_neutrino_energy / 1000. ) < 1.00 )
	  other_universe_x_bin_num = 6;

	else if ( ( background_truth_neutrino_energy / 1000. ) > 1.00 && ( background_truth_neutrino_energy / 1000. ) < 1.25 )
	  other_universe_x_bin_num = 7;

	else if ( ( background_truth_neutrino_energy / 1000. ) > 1.25 && ( background_truth_neutrino_energy / 1000. ) < 1.50 )
	  other_universe_x_bin_num = 8;

	else if	( ( background_truth_neutrino_energy / 1000. ) > 1.50 && ( background_truth_neutrino_energy / 1000. ) < 1.75 )
	  other_universe_x_bin_num = 9;

	else if ( ( background_truth_neutrino_energy / 1000. ) > 1.75 && ( background_truth_neutrino_energy / 1000. ) < 2.00 )
	  other_universe_x_bin_num = 10;

	else if ( ( background_truth_neutrino_energy / 1000. ) > 2.00 && ( background_truth_neutrino_energy / 1000. ) < 2.25 )
	  other_universe_x_bin_num = 11;

	else if ( ( background_truth_neutrino_energy / 1000. ) > 2.25 && ( background_truth_neutrino_energy / 1000. ) < 2.50 )
	  other_universe_x_bin_num = 12;

	else if	( ( background_truth_neutrino_energy / 1000. ) > 2.50 && ( background_truth_neutrino_energy / 1000. ) < 3.00 )
	  other_universe_x_bin_num = 13;

	else if ( ( background_truth_neutrino_energy / 1000. ) > 3.00 && ( background_truth_neutrino_energy / 1000. ) < 4.00 )
	  other_universe_x_bin_num = 14;

	else if ( ( background_truth_neutrino_energy / 1000. ) > 4.00 && ( background_truth_neutrino_energy / 1000. ) < 5.00 )
	  other_universe_x_bin_num = 15;

	else if ( ( background_truth_neutrino_energy / 1000. ) > 5.00 && ( background_truth_neutrino_energy / 1000. ) < 6.00 )
	  other_universe_x_bin_num = 16;

	else if	( ( background_truth_neutrino_energy / 1000. ) > 6.00 && ( background_truth_neutrino_energy / 1000. ) < 7.00 )
	  other_universe_x_bin_num = 17;

	else
	  other_universe_x_bin_num = 18;

	if ( background_truth_angle < 20.00 )
	  other_universe_y_bin_num = 1;

	else if ( background_truth_angle > 20.00 && background_truth_angle < 110.00 )
	  other_universe_y_bin_num = 2;

	else
	  other_universe_y_bin_num = 3;
	
	// Loop through the universes and find out what the weight is for each one.
	for ( int flux_universe_iter = 0; flux_universe_iter < num_flux_universes; flux_universe_iter++ ) {

	  TH2D* input_other_universe_2D_hist = (TH2D*)input_2D_flux_file->Get(Form("numu_PPFXMaster_Uni_%d_AV_TPC_2D", flux_universe_iter));
	  
	  background_flux_other_universe_weight = input_other_universe_2D_hist->GetBinContent( other_universe_x_bin_num, other_universe_y_bin_num );
	  
	  flux_universes.at( flux_universe_iter )->Fill( background_muon_candidate_length, background_flux_other_universe_weight * background_spline_fix_mcweight * background_rootino_fix_mcweight * background_central_value_mcweight * background_special_kdar_weight );
	  num_flux_other_universes_systematic_weighted_nu_events.at( flux_universe_iter ) += ( background_flux_other_universe_weight * background_spline_fix_mcweight * background_rootino_fix_mcweight * background_central_value_mcweight * background_special_kdar_weight );

	  delete input_other_universe_2D_hist;

	}
	  
      }
	
    }
    
  }

  // Loop through the systematic samples separately to fill those histograms.

  // rayleigh length
  int   num_rayleigh_length_entries                               = rayleigh_length->GetEntries();
  float num_weighted_rayleigh_length_events_above_score_threshold = 0;
  
  for ( int iter = 0; iter < num_rayleigh_length_entries; iter++ ) {

    rayleigh_length->GetEntry( iter );

    if ( rayleigh_length_score > score_cut_value ) {

      if ( rayleigh_length_ext == 0 ) {

	muon_candidate_length_nu_background_rayleigh_length->Fill( rayleigh_length_muon_candidate_length, rayleigh_length_spline_fix_mcweight * rayleigh_length_rootino_fix_mcweight * rayleigh_length_central_value_mcweight * rayleigh_length_special_kdar_weight );
	num_weighted_rayleigh_length_events_above_score_threshold += ( rayleigh_length_spline_fix_mcweight * rayleigh_length_rootino_fix_mcweight * rayleigh_length_central_value_mcweight * rayleigh_length_special_kdar_weight );

      }

    }

  }

  // light yield
  int   num_light_yield_entries                               = light_yield->GetEntries();
  float	num_weighted_light_yield_events_above_score_threshold = 0;

  for ( int iter = 0; iter < num_light_yield_entries; iter++ ) {

    light_yield->GetEntry( iter );

    if ( light_yield_score > score_cut_value ) {

      if ( light_yield_ext == 0 ) {

        muon_candidate_length_nu_background_light_yield->Fill( light_yield_muon_candidate_length, light_yield_spline_fix_mcweight * light_yield_rootino_fix_mcweight * light_yield_central_value_mcweight * light_yield_special_kdar_weight );
	num_weighted_light_yield_events_above_score_threshold += ( light_yield_spline_fix_mcweight * light_yield_rootino_fix_mcweight * light_yield_central_value_mcweight * light_yield_special_kdar_weight );

      }

    }

  }

  // attenuation
  int   num_attenuation_entries                               = attenuation->GetEntries();
  float	num_weighted_attenuation_events_above_score_threshold = 0;
  
  for ( int iter = 0; iter < num_attenuation_entries; iter++ ) {

    attenuation->GetEntry( iter );

    if ( attenuation_score > score_cut_value ) {

      if ( attenuation_ext == 0 ) {

        muon_candidate_length_nu_background_attenuation->Fill( attenuation_muon_candidate_length, attenuation_spline_fix_mcweight * attenuation_rootino_fix_mcweight * attenuation_central_value_mcweight * attenuation_special_kdar_weight );
	num_weighted_attenuation_events_above_score_threshold += ( attenuation_spline_fix_mcweight * attenuation_rootino_fix_mcweight * attenuation_central_value_mcweight * attenuation_special_kdar_weight );
	
      }

    }

  }

  // sce
  int   num_sce_entries                               = sce->GetEntries();
  float	num_weighted_sce_events_above_score_threshold = 0;
  
  for ( int iter = 0; iter < num_sce_entries; iter++ ) {

    sce->GetEntry( iter );

    if ( sce_score > score_cut_value ) {

      if ( sce_ext == 0 ) {

        muon_candidate_length_nu_background_sce->Fill( sce_muon_candidate_length, sce_spline_fix_mcweight * sce_rootino_fix_mcweight * sce_central_value_mcweight * sce_special_kdar_weight );
	num_weighted_sce_events_above_score_threshold += ( sce_spline_fix_mcweight * sce_rootino_fix_mcweight * sce_central_value_mcweight * sce_special_kdar_weight );
	
      }

    }

  }

  // recombination
  int   num_recombination_entries                               = recombination->GetEntries();
  float	num_weighted_recombination_events_above_score_threshold = 0;
  
  for ( int iter = 0; iter < num_recombination_entries; iter++ ) {

    recombination->GetEntry( iter );

    if ( recombination_score > score_cut_value ) {

      if ( recombination_ext == 0 ) {

        muon_candidate_length_nu_background_recombination->Fill( recombination_muon_candidate_length, recombination_spline_fix_mcweight * recombination_rootino_fix_mcweight * recombination_central_value_mcweight * recombination_special_kdar_weight );
	num_weighted_recombination_events_above_score_threshold += ( recombination_spline_fix_mcweight * recombination_rootino_fix_mcweight * recombination_central_value_mcweight * recombination_special_kdar_weight );
	
      }

    }

  }

  // scalex
  int   num_scalex_entries                               = scalex->GetEntries();
  float	num_weighted_scalex_events_above_score_threshold = 0;
  
  for ( int iter = 0; iter < num_scalex_entries; iter++ ) {

    scalex->GetEntry( iter );

    if ( scalex_score > score_cut_value ) {

      if ( scalex_ext == 0 ) {

        muon_candidate_length_nu_background_scalex->Fill( scalex_muon_candidate_length, scalex_spline_fix_mcweight * scalex_rootino_fix_mcweight * scalex_central_value_mcweight * scalex_special_kdar_weight );
	num_weighted_scalex_events_above_score_threshold += ( scalex_spline_fix_mcweight * scalex_rootino_fix_mcweight * scalex_central_value_mcweight * scalex_special_kdar_weight );

      }

    }

  }

  // scaleyz
  int   num_scaleyz_entries                               = scaleyz->GetEntries();
  float	num_weighted_scaleyz_events_above_score_threshold = 0;
  
  for ( int iter = 0; iter < num_scaleyz_entries; iter++ ) {

    scaleyz->GetEntry( iter );

    if ( scaleyz_score > score_cut_value ) {

      if ( scaleyz_ext == 0 ) {

        muon_candidate_length_nu_background_scaleyz->Fill( scaleyz_muon_candidate_length, scaleyz_spline_fix_mcweight * scaleyz_rootino_fix_mcweight * scaleyz_central_value_mcweight * scaleyz_special_kdar_weight );
	num_weighted_scaleyz_events_above_score_threshold += ( scaleyz_spline_fix_mcweight * scaleyz_rootino_fix_mcweight * scaleyz_central_value_mcweight * scaleyz_special_kdar_weight );

      }

    }

  }

  // scaleangleyz
  int   num_scaleangleyz_entries                               = scaleangleyz->GetEntries();
  float num_weighted_scaleangleyz_events_above_score_threshold = 0;
  
  for ( int iter = 0; iter < num_scaleangleyz_entries; iter++ ) {

    scaleangleyz->GetEntry( iter );

    if ( scaleangleyz_score > score_cut_value ) {

      if ( scaleangleyz_ext == 0 ) {

        muon_candidate_length_nu_background_scaleangleyz->Fill( scaleangleyz_muon_candidate_length, scaleangleyz_spline_fix_mcweight * scaleangleyz_rootino_fix_mcweight * scaleangleyz_central_value_mcweight * scaleangleyz_special_kdar_weight );
	num_weighted_scaleangleyz_events_above_score_threshold += ( scaleangleyz_spline_fix_mcweight * scaleangleyz_rootino_fix_mcweight * scaleangleyz_central_value_mcweight * scaleangleyz_special_kdar_weight );
	
      }

    }

  }

  // scaleanglexz
  int   num_scaleanglexz_entries                               = scaleanglexz->GetEntries();
  float num_weighted_scaleanglexz_events_above_score_threshold = 0;
  
  for ( int iter = 0; iter < num_scaleanglexz_entries; iter++ ) {

    scaleanglexz->GetEntry( iter );

    if ( scaleanglexz_score > score_cut_value ) {

      if ( scaleanglexz_ext == 0 ) {

        muon_candidate_length_nu_background_scaleanglexz->Fill( scaleanglexz_muon_candidate_length, scaleanglexz_spline_fix_mcweight * scaleanglexz_rootino_fix_mcweight * scaleanglexz_central_value_mcweight * scaleanglexz_special_kdar_weight );
	num_weighted_scaleanglexz_events_above_score_threshold += ( scaleanglexz_spline_fix_mcweight * scaleanglexz_rootino_fix_mcweight * scaleanglexz_central_value_mcweight * scaleanglexz_special_kdar_weight );
	
      }

    }

  }

  // scalededx
  int   num_scalededx_entries                               = scalededx->GetEntries();
  float num_weighted_scalededx_events_above_score_threshold = 0;
  
  for ( int iter = 0; iter < num_scalededx_entries; iter++ ) {

    scalededx->GetEntry( iter );

    if ( scalededx_score > score_cut_value ) {

      if ( scalededx_ext == 0 ) {

        muon_candidate_length_nu_background_scalededx->Fill( scalededx_muon_candidate_length, scalededx_spline_fix_mcweight * scalededx_rootino_fix_mcweight * scalededx_central_value_mcweight * scalededx_special_kdar_weight );
        num_weighted_scalededx_events_above_score_threshold += ( scalededx_spline_fix_mcweight * scalededx_rootino_fix_mcweight * scalededx_central_value_mcweight * scalededx_special_kdar_weight );
	
      }

    }

  }

  
  std::cout << "The number of unweighted EXT events with a score > " << score_cut_value << " = " << num_unweighted_ext_events_above_score_threshold << "." << std::endl;
  std::cout << "The number of unweighted NuMI background events with a score > " << score_cut_value << " = " << num_unweighted_nu_events_above_score_threshold << "." << std::endl;

  // Normalize each of the samples.
  signal_normalized_events              = ( ( num_weighted_signal_events_above_score_threshold /  51794.259 ) *  74.515 );
  background_numi_ext_normalized_events = ( ( num_unweighted_ext_events_above_score_threshold / 5472.0 ) *  875.856 );
  background_numi_nu_normalized_events  = ( ( num_weighted_nu_events_above_score_threshold / 28225.5 ) * 1335.183  );
  data_normalized_events                = num_data_events_above_score_threshold;

  // Introduce the factor that you will scale all of the systematic plots by.
  systematic_normalization_factor       = 0.0473;

  // Calculate the number of normalized events for each of the detector systematic samples.
  rayleigh_length_normalized_events     = ( systematic_normalization_factor * num_weighted_rayleigh_length_events_above_score_threshold );
  light_yield_normalized_events         = ( systematic_normalization_factor * num_weighted_light_yield_events_above_score_threshold );
  attenuation_normalized_events         = ( systematic_normalization_factor * num_weighted_attenuation_events_above_score_threshold );
  sce_normalized_events                 = ( systematic_normalization_factor * num_weighted_sce_events_above_score_threshold );
  recombination_normalized_events       = ( systematic_normalization_factor * num_weighted_recombination_events_above_score_threshold );
  scalex_normalized_events              = ( systematic_normalization_factor * num_weighted_scalex_events_above_score_threshold );
  scaleyz_normalized_events             = ( systematic_normalization_factor * num_weighted_scaleyz_events_above_score_threshold );
  scaleangleyz_normalized_events        = ( systematic_normalization_factor * num_weighted_scaleangleyz_events_above_score_threshold );
  scaleanglexz_normalized_events        = ( systematic_normalization_factor * num_weighted_scaleanglexz_events_above_score_threshold );
  scalededx_normalized_events           = ( systematic_normalization_factor * num_weighted_scalededx_events_above_score_threshold );
  
  // Calculate the number of normalized events for each of the xsec systematic samples.
  for ( size_t other_universe_iter = 0; other_universe_iter < 100; other_universe_iter++ ) {

    other_xsec_universes_normalized_events.at( other_universe_iter ) = ( num_weighted_nu_events_above_threshold_other_universes.at( other_universe_iter ) * systematic_normalization_factor );

  }

  for ( size_t axialff_iter = 0; axialff_iter < 2; axialff_iter++ ) {

    axialff_universes_normalized_events.at( axialff_iter ) = ( num_weighted_nu_events_above_threshold_axialff.at( axialff_iter ) * systematic_normalization_factor );

  }

  for ( size_t rpaccqe_iter = 0; rpaccqe_iter < 2; rpaccqe_iter++ ) {

    rpaccqe_universes_normalized_events.at( rpaccqe_iter ) = ( num_weighted_nu_events_above_threshold_rpaccqe.at( rpaccqe_iter ) * systematic_normalization_factor );

  }

  for ( size_t xsecshape_iter = 0; xsecshape_iter < 2; xsecshape_iter++ ) {

    xsecshape_universes_normalized_events.at( xsecshape_iter ) = ( num_weighted_nu_events_above_threshold_xsecshape.at( xsecshape_iter ) * systematic_normalization_factor );

  }

  // Calculate the number of normalized events for each of the flux systematic samples.
  flux_cv_normalized_events = ( systematic_normalization_factor * num_flux_cv_systematic_weighted_nu_events ); 

  for ( size_t other_flux_universe_iter = 0; other_flux_universe_iter < num_flux_universes; other_flux_universe_iter++ ) {

    flux_other_universes_normalized_events.at( other_flux_universe_iter ) = ( systematic_normalization_factor * num_flux_other_universes_systematic_weighted_nu_events.at( other_flux_universe_iter ) );

  }
    
  std::cout << "The number of normalized signal events = " << signal_normalized_events << "." << std::endl;
  std::cout << "The number of normalized NuMI EXT events = " << background_numi_ext_normalized_events << "." << std::endl;
  std::cout << "The number of normalized background nu events = " << background_numi_nu_normalized_events << "." << std::endl;
  std::cout << "The total number of normalized background events (NuMI EXT + background nu) = " << ( background_numi_ext_normalized_events + background_numi_nu_normalized_events ) << "." << std::endl;
  std::cout << "The number of data events = " << data_normalized_events << "." << std::endl;

  std::cout << "The significance ( Signal / Sqrt( Background + Signal ) ) = " << ( signal_normalized_events / TMath::Sqrt( ( background_numi_ext_normalized_events + background_numi_nu_normalized_events ) ) ) << "." << std::endl;  

  // Scale each of the histograms with the appropriate number of events.
  muon_candidate_length_signal->Scale( signal_normalized_events / muon_candidate_length_signal->Integral() );
  muon_candidate_length_ext_background->Scale(  background_numi_ext_normalized_events / muon_candidate_length_ext_background->Integral() );
  muon_candidate_length_nu_background->Scale( background_numi_nu_normalized_events / muon_candidate_length_nu_background->Integral() );
  muon_candidate_length_data->Scale( data_normalized_events / muon_candidate_length_data->Integral() );

  // Normalize the detector systematic samples.
  muon_candidate_length_nu_background_rayleigh_length->Scale( rayleigh_length_normalized_events / muon_candidate_length_nu_background_rayleigh_length->Integral() );
  muon_candidate_length_nu_background_light_yield->Scale(     light_yield_normalized_events     / muon_candidate_length_nu_background_light_yield->Integral() );
  muon_candidate_length_nu_background_attenuation->Scale(     attenuation_normalized_events     / muon_candidate_length_nu_background_attenuation->Integral() );
  muon_candidate_length_nu_background_sce->Scale(             sce_normalized_events             / muon_candidate_length_nu_background_sce->Integral() );
  muon_candidate_length_nu_background_recombination->Scale(   recombination_normalized_events   / muon_candidate_length_nu_background_recombination->Integral() );
  muon_candidate_length_nu_background_scalex->Scale(          scalex_normalized_events          / muon_candidate_length_nu_background_scalex->Integral() );
  muon_candidate_length_nu_background_scaleyz->Scale(         scaleyz_normalized_events         / muon_candidate_length_nu_background_scaleyz->Integral() );
  muon_candidate_length_nu_background_scaleangleyz->Scale(    scaleangleyz_normalized_events    / muon_candidate_length_nu_background_scaleangleyz->Integral() );
  muon_candidate_length_nu_background_scaleanglexz->Scale(    scaleanglexz_normalized_events    / muon_candidate_length_nu_background_scaleanglexz->Integral() );
  muon_candidate_length_nu_background_scalededx->Scale(       scalededx_normalized_events       / muon_candidate_length_nu_background_scalededx->Integral() );

  // Normalize the xsec systematic samples and calculate the rms value for each bin as well.
  std::vector< float > xsec_other_universe_rms_values;
  std::vector< float > axialff_rms_values;
  std::vector< float > rpaccqe_rms_values;
  std::vector< float > xsecshape_rms_values;

  xsec_other_universe_rms_values.resize( 20, 0. );
  axialff_rms_values.resize( 20, 0. );
  rpaccqe_rms_values.resize( 20, 0. );
  xsecshape_rms_values.resize( 20, 0. );

  // XSec Other Universe
  for ( size_t xsec_other_universe_hist_iter = 0; xsec_other_universe_hist_iter < 100; xsec_other_universe_hist_iter++ ) {

    other_xsec_universes.at( xsec_other_universe_hist_iter )->Scale( other_xsec_universes_normalized_events.at( xsec_other_universe_hist_iter ) / other_xsec_universes.at( xsec_other_universe_hist_iter )->Integral() );

    for ( size_t bin_iter = 1; bin_iter < 21; bin_iter++ ) {

      xsec_other_universe_rms_values.at( bin_iter - 1 ) += ( TMath::Power( ( other_xsec_universes.at( xsec_other_universe_hist_iter )->GetBinContent( bin_iter ) - muon_candidate_length_nu_background->GetBinContent( bin_iter ) ), 2 ) * 0.01 );

    }
    
    std::cout << "The number of events in bin #6 of universe #" << xsec_other_universe_hist_iter << " = " << other_xsec_universes.at( xsec_other_universe_hist_iter )->GetBinContent( 6 ) << "." << std::endl;
    
  }

  // Axialff
  for ( size_t axialff_hist_iter = 0; axialff_hist_iter < 2; axialff_hist_iter++ ) {

    axialff_universes.at( axialff_hist_iter )->Scale( axialff_universes_normalized_events.at( axialff_hist_iter ) / axialff_universes.at( axialff_hist_iter )->Integral() );

    for ( size_t bin_iter = 1; bin_iter < 21; bin_iter++ ) {

      axialff_rms_values.at( bin_iter - 1 ) += ( TMath::Power( ( axialff_universes.at( axialff_hist_iter )->GetBinContent( bin_iter ) - muon_candidate_length_nu_background->GetBinContent( bin_iter ) ), 2 ) * 0.5 );

    }
    
    
  }

  // RPACCQE
  for ( size_t rpaccqe_hist_iter = 0; rpaccqe_hist_iter < 2; rpaccqe_hist_iter++ ) {

    rpaccqe_universes.at( rpaccqe_hist_iter )->Scale( rpaccqe_universes_normalized_events.at( rpaccqe_hist_iter ) / rpaccqe_universes.at( rpaccqe_hist_iter )->Integral() );

    for ( size_t bin_iter = 1; bin_iter < 21; bin_iter++ ) {

      rpaccqe_rms_values.at( bin_iter - 1 ) += ( TMath::Power( ( rpaccqe_universes.at( rpaccqe_hist_iter )->GetBinContent( bin_iter ) - muon_candidate_length_nu_background->GetBinContent( bin_iter ) ), 2 ) * 0.5 );

    }
    
    
  }

  // XSecShape
  for ( size_t xsecshape_hist_iter = 0; xsecshape_hist_iter < 2; xsecshape_hist_iter++ ) {

    xsecshape_universes.at( xsecshape_hist_iter )->Scale( xsecshape_universes_normalized_events.at( xsecshape_hist_iter ) / xsecshape_universes.at( xsecshape_hist_iter )->Integral() );

    if ( xsecshape_hist_iter == 1 ) {

      for ( size_t bin_iter = 1; bin_iter < 21; bin_iter++ ) {

	xsecshape_rms_values.at( bin_iter - 1 ) = ( xsecshape_universes.at( xsecshape_hist_iter )->GetBinContent( bin_iter ) - muon_candidate_length_nu_background->GetBinContent( bin_iter ) );

      }

    }
    
  }
  
  // Declare the variables that will be used to
  float bin_rms_holder_value              = 0.;

  // Declare the histograms with the rms values.
  TH1F* muon_candidate_length_nu_background_other_xsec_universes_rms_values    = new TH1F("other_xsec_universes_rms_values", "Other XSec Universes RMS Values", 20, 0., 40. );
  TH1F* muon_candidate_length_nu_background_axialff_rms_values                 = new TH1F("axialff_rms_values", "Axialff Universes RMS Values", 20, 0., 40. );
  TH1F* muon_candidate_length_nu_background_rpaccqe_rms_values                 = new TH1F("rpaccqe_rms_values", "RPACCQE Universes RMS Values", 20, 0., 40. );
  TH1F* muon_candidate_length_nu_background_xsecshape_rms_values               = new TH1F("xsecshape_rms_values", "XSecShape Universes RMS Values", 20, 0., 40.);
  
  // Generate a histogram with the rms values for all xsec uncertainties but the 'xsecshape_universes' (we only take the difference between that one and the CV).

  // Other RMS Values
  for ( size_t other_universe_bin_num = 1; other_universe_bin_num < 21; other_universe_bin_num++ ) {

    muon_candidate_length_nu_background_other_xsec_universes_rms_values->SetBinContent( other_universe_bin_num, TMath::Sqrt( xsec_other_universe_rms_values.at( other_universe_bin_num - 1 ) ) );
    
  }

  // Axialff RMS Values
  for ( size_t axialff_bin_num = 1; axialff_bin_num < 21; axialff_bin_num++ ) {

    muon_candidate_length_nu_background_axialff_rms_values->SetBinContent( axialff_bin_num, TMath::Sqrt( axialff_rms_values.at( axialff_bin_num - 1 ) ) );

  }

  // RPACCQE RMS Values
  for ( size_t rpaccqe_bin_num = 1; rpaccqe_bin_num < 21; rpaccqe_bin_num++ ) {

    muon_candidate_length_nu_background_rpaccqe_rms_values->SetBinContent( rpaccqe_bin_num, TMath::Sqrt( rpaccqe_rms_values.at( rpaccqe_bin_num - 1 ) ) );

  }

  // XSecShape RMS Values
  for ( size_t xsecshape_bin_num = 1; xsecshape_bin_num < 21; xsecshape_bin_num++ ) {

    muon_candidate_length_nu_background_xsecshape_rms_values->SetBinContent( xsecshape_bin_num, xsecshape_rms_values.at( xsecshape_bin_num - 1 ) );

  }

  // Normalize the background flux CV sample.
  muon_candidate_length_nu_background_flux_cv->Scale( flux_cv_normalized_events / muon_candidate_length_nu_background_flux_cv->Integral() );

  // Make a histogram of the difference of the flux cv correction wrt the background cv sample for the bins in the muon candidate length plot.
  TH1F* muon_candidate_length_nu_background_flux_cv_difference = new TH1F("muon_candidate_length_nu_background_flux_cv_difference", "Flux CV Correction Differences", 20, 0., 40. );

  for ( size_t flux_cv_bin_iter = 1;  flux_cv_bin_iter < 21; flux_cv_bin_iter++ ) {

    muon_candidate_length_nu_background_flux_cv_difference->SetBinContent( flux_cv_bin_iter, fabs( muon_candidate_length_nu_background_flux_cv->GetBinContent( flux_cv_bin_iter ) - muon_candidate_length_nu_background->GetBinContent( flux_cv_bin_iter ) ) );

  }
  
  std::vector< float > flux_other_universe_rms_values;
  flux_other_universe_rms_values.resize( 20, 0. );
  
  // Normalize the background other flux universes samples.
  for ( int flux_universe_iter = 0; flux_universe_iter < num_flux_universes; flux_universe_iter++ ) {

    flux_universes.at( flux_universe_iter )->Scale( flux_other_universes_normalized_events.at( flux_universe_iter ) / flux_universes.at( flux_universe_iter )->Integral() );

    for ( size_t bin_iter = 1; bin_iter < 21; bin_iter++ ) {

      flux_other_universe_rms_values.at( bin_iter - 1 ) += ( TMath::Power( flux_universes.at( flux_universe_iter )->GetBinContent( bin_iter ) - muon_candidate_length_nu_background->GetBinContent( bin_iter ), 2 ) * 0.01 );

    }
    
  }

  // Make a histogram of the rms values of the flux universe distribution for the bins in the muon candidate length plot.
  TH1F* muon_candidate_length_nu_background_other_flux_universes_rms_values    = new TH1F("other_flux_universes_rms_values", "Other Flux Universes RMS Values", 20, 0., 40. );

  for ( int flux_universe_bin_num = 1; flux_universe_bin_num < 21; flux_universe_bin_num++ ) {

    muon_candidate_length_nu_background_other_flux_universes_rms_values->SetBinContent( flux_universe_bin_num, TMath::Sqrt( flux_other_universe_rms_values.at( flux_universe_bin_num - 1 ) ) );
    
  }
  
  // Create variables for the difference between the entries in these histograms and those in the background cv sample.                                                                                 
  float rayleigh_length_and_cv_difference      = 0.;
  float light_yield_and_cv_difference          = 0.;
  float attenuation_and_cv_difference          = 0.;
  float sce_and_cv_difference                  = 0.;
  float recombination_and_cv_difference        = 0.;
  float scalex_and_cv_difference               = 0.;
  float scaleyz_and_cv_difference              = 0.;
  float scaleangleyz_and_cv_difference         = 0.;
  float scaleanglexz_and_cv_difference         = 0.;
  float scalededx_and_cv_difference            = 0.;

  // Declare the difference for the xsec systematics and the central value.
  float other_universes_and_cv_difference      = 0.;
  float axialff_and_cv_differences             = 0.;
  float rpaccqe_and_cv_differences             = 0.;
  float xsecshape_and_cv_differences           = 0.;

  // Declare the difference for the flux systematics and the central value.
  float flux_cv_and_cv_difference              = 0.;
  float flux_universes_and_cv_difference       = 0.;

  float total_detector_syst_contribution       = 0.;
  float total_flux_syst_contribution           = 0.;
  float total_xsec_syst_contribution           = 0.;
  float total_bkgd_statistical_contribution    = 0.;
  float total_data_statistical_contribution    = 0.;
  
  // Declare new histograms that will still contain all of the information but can be manipulated properly.
  TH1D* muon_candidate_length_ext_background_copy                      = new TH1D("muon_candidate_length_ext_background_copy", "NuMI EXT Muon Candidate Length", 20, 0., 40.);
  TH1D* muon_candidate_length_nu_background_copy                       = new TH1D("muon_candidate_length_nu_background_copy", "NuMI Nu Background Muon Candidate Length", 20, 0., 40.);
  TH1D* muon_candidate_length_signal_copy                              = new TH1D("muon_candidate_length_signal_copy", "Signal Muon Candidate Length", 20, 0., 40.);
  
  for ( size_t bin_num = 1; bin_num < 21; bin_num++ ) {

    muon_candidate_length_ext_background_copy->SetBinContent( bin_num, muon_candidate_length_ext_background->GetBinContent( bin_num ) );
    muon_candidate_length_nu_background_copy->SetBinContent( bin_num, muon_candidate_length_nu_background->GetBinContent( bin_num ) );
    muon_candidate_length_signal_copy->SetBinContent( bin_num, muon_candidate_length_signal->GetBinContent( bin_num ) );

    // Calculate the differences between the detector systematic samples and the cv sample.
    rayleigh_length_and_cv_difference = ( muon_candidate_length_nu_background->GetBinContent( bin_num ) - muon_candidate_length_nu_background_rayleigh_length->GetBinContent( bin_num ) );
    light_yield_and_cv_difference     = ( muon_candidate_length_nu_background->GetBinContent( bin_num ) - muon_candidate_length_nu_background_light_yield->GetBinContent( bin_num ) );
    attenuation_and_cv_difference     = ( muon_candidate_length_nu_background->GetBinContent( bin_num ) - muon_candidate_length_nu_background_attenuation->GetBinContent( bin_num ) );
    sce_and_cv_difference             = ( muon_candidate_length_nu_background->GetBinContent( bin_num ) - muon_candidate_length_nu_background_sce->GetBinContent( bin_num ) );
    recombination_and_cv_difference   = ( muon_candidate_length_nu_background->GetBinContent( bin_num ) - muon_candidate_length_nu_background_recombination->GetBinContent( bin_num ) );
    scalex_and_cv_difference          = ( muon_candidate_length_nu_background->GetBinContent( bin_num ) - muon_candidate_length_nu_background_scalex->GetBinContent( bin_num ) );
    scaleyz_and_cv_difference         = ( muon_candidate_length_nu_background->GetBinContent( bin_num ) - muon_candidate_length_nu_background_scaleyz->GetBinContent( bin_num ) );
    scaleangleyz_and_cv_difference    = ( muon_candidate_length_nu_background->GetBinContent( bin_num ) - muon_candidate_length_nu_background_scaleangleyz->GetBinContent( bin_num ) );
    scaleanglexz_and_cv_difference    = ( muon_candidate_length_nu_background->GetBinContent( bin_num ) - muon_candidate_length_nu_background_scaleanglexz->GetBinContent( bin_num ) );
    scalededx_and_cv_difference       = ( muon_candidate_length_nu_background->GetBinContent( bin_num ) - muon_candidate_length_nu_background_scalededx->GetBinContent( bin_num ) );

    // Calculate the differences between the xsec systematic samples and the CV sample.
    other_universes_and_cv_difference = muon_candidate_length_nu_background_other_xsec_universes_rms_values->GetBinContent( bin_num );
    axialff_and_cv_differences        = muon_candidate_length_nu_background_axialff_rms_values->GetBinContent( bin_num );
    rpaccqe_and_cv_differences        = muon_candidate_length_nu_background_rpaccqe_rms_values->GetBinContent( bin_num );
    xsecshape_and_cv_differences      = muon_candidate_length_nu_background_xsecshape_rms_values->GetBinContent( bin_num );

    // Calculate the differences between the flux systematic samples and the CV sample.
    flux_cv_and_cv_difference         = muon_candidate_length_nu_background_flux_cv_difference->GetBinContent( bin_num );

    flux_universes_and_cv_difference  = muon_candidate_length_nu_background_other_flux_universes_rms_values->GetBinContent( bin_num );

    // Calculate the total of each type of uncertainty.
    total_detector_syst_contribution      = TMath::Power( rayleigh_length_and_cv_difference, 2 ) + TMath::Power( light_yield_and_cv_difference, 2 ) + TMath::Power( attenuation_and_cv_difference, 2 ) + TMath::Power( sce_and_cv_difference, 2 ) + TMath::Power( recombination_and_cv_difference, 2 ) + TMath::Power( scalex_and_cv_difference, 2 ) + TMath::Power( scaleyz_and_cv_difference, 2 ) + TMath::Power( scaleangleyz_and_cv_difference, 2 ) + TMath::Power( scaleanglexz_and_cv_difference, 2 ) + TMath::Power( scalededx_and_cv_difference, 2 );
    total_xsec_syst_contribution          = TMath::Power( other_universes_and_cv_difference, 2 ) + TMath::Power( axialff_and_cv_differences, 2 ) + TMath::Power( rpaccqe_and_cv_differences, 2 ) + TMath::Power( xsecshape_and_cv_differences, 2 );
    total_flux_syst_contribution          = TMath::Power( flux_cv_and_cv_difference, 2 ) + TMath::Power( flux_universes_and_cv_difference, 2 );
    total_data_statistical_contribution   = ( muon_candidate_length_data_before_subtraction->GetBinContent( bin_num ) );

    total_bkgd_statistical_contribution    = ( TMath::Power( ( 875.856 / 5472.0 ), 2 ) * muon_candidate_length_ext_background_before_normalization->GetBinContent( bin_num ) + TMath::Power( ( 1335.183 / 28225.5 ), 2 ) * muon_candidate_length_nu_background_without_weights->GetBinContent( bin_num ) );

    // Set the histograms that will be used to make the plot.
    muon_candidate_length_data->SetBinError( bin_num, TMath::Sqrt( total_data_statistical_contribution ) );
    muon_candidate_length_ext_background->SetBinContent( bin_num, muon_candidate_length_ext_background->GetBinContent( bin_num ) + muon_candidate_length_nu_background->GetBinContent( bin_num ) );
    muon_candidate_length_ext_background->SetBinError( bin_num, TMath::Sqrt( total_bkgd_statistical_contribution + total_detector_syst_contribution + total_xsec_syst_contribution + total_flux_syst_contribution ) );
						
    if ( bin_num > 1 && bin_num < 15 ) {

      std::cout << "Bin #" << bin_num << ":" << std::endl;
      std::cout << "CV = " << muon_candidate_length_nu_background->GetBinContent( bin_num ) << "." << std::endl;
      std::cout << "Total Detector Systematic Contribution = " << TMath::Sqrt( total_detector_syst_contribution ) << "." << std::endl;
      std::cout << "Flux CV Correction Systematic = " << flux_cv_and_cv_difference << "." << std::endl;
      std::cout << "Flux Universes Systematic = " << flux_universes_and_cv_difference << "." << std::endl;
      std::cout << "Total Flux Systematic Contribution = " << TMath::Sqrt( total_flux_syst_contribution ) << "." << std::endl;
      std::cout << "XSec Other Universes Systematic = " << other_universes_and_cv_difference << "." << std::endl;
      std::cout << "XSec Axialff Systematic = " << axialff_and_cv_differences << "." << std::endl;
      std::cout << "XSec RPACCQE Systematic = " << rpaccqe_and_cv_differences << "." << std::endl;
      std::cout << "XSec XSecShape Systematic = " << xsecshape_and_cv_differences << "." << std::endl;
      std::cout << "Total XSec Systematic Contribution = " << TMath::Sqrt( total_xsec_syst_contribution ) << "." << std::endl;
      std::cout << "Total Systematic Contribution = " << TMath::Sqrt( total_detector_syst_contribution + total_xsec_syst_contribution + total_flux_syst_contribution ) << "." << std::endl;
      std::cout << "Total Data Statistical Contribution = " <<  total_data_statistical_contribution << "." << std::endl;
      std::cout << "Total Bkgd Statistical Contribution = " << total_bkgd_statistical_contribution << "." << std::endl;
      std::cout << std::endl;
      
    }
    
  }

  // Print out the integral of the KDAR signal plot.
  std::cout << "The integral of the KDAR signal plot copy = " << muon_candidate_length_signal_copy->Integral() << "." << std::endl;

  // Fill the "stacked" plot.
  muon_candidate_length_nu_background_copy->SetFillColor(kOrange);
  muon_candidate_length_ext_background_copy->SetFillColor(kGreen);
  muon_candidate_length_stack->Add(muon_candidate_length_nu_background_copy);
  muon_candidate_length_stack->Add(muon_candidate_length_ext_background_copy);


  // Declare the canvas and fill it with each of the plots.
  TCanvas* muon_candidate_length_canvas = new TCanvas("muon_candidate_length_canvas", "muon_candidate_length_canvas", 600, 600);
  muon_candidate_length_canvas->cd();
  muon_candidate_length_stack->SetMaximum( 50.0 );
  gStyle->SetOptStat(0);
  muon_candidate_length_stack->Draw();
  muon_candidate_length_stack->GetXaxis()->SetTitle("Muon Candidate Length [cm]");
  muon_candidate_length_stack->GetYaxis()->SetTitle("Events Per Bin");
  muon_candidate_length_stack->GetXaxis()->SetTitleOffset(1.10);
  muon_candidate_length_stack->GetYaxis()->SetTitleOffset(1.40);
  muon_candidate_length_stack->GetXaxis()->SetTickLength(0.015);
  muon_candidate_length_stack->GetYaxis()->SetTickLength(0.015);
  muon_candidate_length_stack->GetXaxis()->SetNdivisions( 8 );
  muon_candidate_length_stack->GetYaxis()->SetNdivisions( 10 );
  muon_candidate_length_stack->SetTitle(Form("Reconstructed Muon Candidate Length W/ All Uncertainties: BDT Score Cut = %.01f", score_cut_value));
  gPad->Modified();
  muon_candidate_length_data->SetLineColor(kBlack);
  muon_candidate_length_data->Draw("Sames E1");
  muon_candidate_length_ext_background->SetMarkerSize(0.);
  muon_candidate_length_ext_background->SetFillStyle(3003);
  muon_candidate_length_ext_background->SetMarkerColor(kBlack);
  muon_candidate_length_ext_background->SetFillColor(kRed);
  muon_candidate_length_ext_background->Draw("Same E2 P");
  muon_candidate_length_ext_background_copy->SetMarkerColor(kWhite); // Test For the Legends.

  // Save the canvas to the output file.
  file->cd();
  muon_candidate_length_canvas->Write();
    
  file->Write();
  file->Close();

}
