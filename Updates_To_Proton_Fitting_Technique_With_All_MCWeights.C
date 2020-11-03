#include <TMath.h>
#include <TGraph2D.h>
#include <TRandom2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TF2.h>
#include <TH1.h>
#include <TChain.h>
#include <TTree.h>
#include <TVirtualFitter.h>
#include <TPolyLine3D.h>
#include <Math/Vector3D.h>
#include <iostream>

using namespace ROOT::Math;


// define the parameteric line equation 
void line(double t, double *p, double &x, double &y, double &z) { 
   // a parameteric line is define from 6 parameters but 4 are independent
   // x0,y0,z0,z1,y1,z1 which are the coordinates of two points on the line
   // can choose z0 = 0 if line not parallel to x-y plane and z1 = 1; 
   x = p[0] + p[1]*t; 
   y = p[2] + p[3]*t;
   z = t; 
} 

// calculate distance line-point 
double distance2(double x,double y,double z, double *p) { 
   // distance line point is D= | (xp-x0) cross  ux | 
   // where ux is direction of line and x0 is a point in the line (like t = 0) 
   XYZVector xp(x,y,z); 
   XYZVector x0(p[0], p[2], 0. ); 
   XYZVector x1(p[0] + p[1], p[2] + p[3], 1. ); 
   XYZVector u = (x1-x0).Unit(); 
   double d2 = ((xp-x0).Cross(u)) .Mag2(); 
   return d2; 
}
bool first = true; 


// function to be minimized 
void SumDistance2(int &, double *, double & sum, double * par, int ) { 
   // the TGraph must be a global variable
   TGraph2D * gr = dynamic_cast<TGraph2D*>( (TVirtualFitter::GetFitter())->GetObjectFit() );
   assert(gr != 0);
   double * x = gr->GetX();
   double * y = gr->GetY();
   double * z = gr->GetZ();
   int npoints = gr->GetN();
   sum = 0;
   for (int i  = 0; i < npoints; ++i) { 
      double d = distance2(x[i],y[i],z[i],par); 
      sum += d;
   }
   first = false;
}

void Updates_To_Proton_Fitting_Technique_With_All_MCWeights() 
{
   gStyle->SetOptStat(0);
   gStyle->SetOptFit();

   // Declare the neutrino direction vector up here.
   double truth_neutrino_direction[3];
   truth_neutrino_direction[0] = 0.08717;
   truth_neutrino_direction[1] = 0.178;
   truth_neutrino_direction[2] = -0.1025;
   
   int    run;
   int    subrun;
   int    event;
   double truth_neutrino_energy;
   double truth_total_proton_KE;
   double truth_leading_proton_KE;
   double truth_reco_vtx_distance;
   int    kdar_from_dump_event;
   int    NC_channel;
   int    num_muminus_tracks;
   int    num_muplus_tracks;
   int    num_piplus_tracks;
   int    num_piminus_tracks;
   int    num_pi0_tracks;
   int    num_proton_tracks;
   int    num_electron_showers;
   int    num_positron_showers;
   int    num_photon_showers;
   double spline_fix_mcweight;
   double rootino_fix_mcweight;
   double central_value_mcweight;
   double other_universe_mcweights[100];
   double axialff_mcweights[2];
   double rpaccqe_mcweights[2];
   double xsecshape_mcweights[2];
   float  reco_vtx_x;
   float  reco_vtx_y;
   float  reco_vtx_z;
   int    muon_track_is_correctly_oriented;
   double x_slope;
   double x_intercept;
   double y_slope;
   double y_intercept;
   double x_slope_err;
   double x_intercept_err;
   double y_slope_err;
   double y_intercept_err;
   double distance_between_furthest_calo_points;
   double distance_using_track_projection;
   float  muon_track_length;
   float  reco_muon_KE;
   float  muon_x_momentum;
   float  muon_y_momentum;
   float  muon_z_momentum;
   float  muon_phi_without_orientation;
   float  muon_cos_theta_without_orientation;
   float  proton_kinetic_energy_from_additional_tracks;
   float  x_momentum_component_sum_additional_tracks;
   float  y_momentum_component_sum_additional_tracks;
   float  z_momentum_component_sum_additional_tracks;
   float  proton_x_projection;
   float  proton_y_projection;
   float  proton_z_projection;
   float  x_component_of_muon_momentum;
   float  y_component_of_muon_momentum;
   float  z_component_of_muon_momentum;
   int    num_calo_points_used_in_fit;
   float  num_beam_flashes;
   float  flash_PEs;
   float  flash_z;
   float  flash_y;
   float  truncated_dQdx_of_TPCObject_muon_candidate_u_plane;
   float  truncated_dQdx_of_TPCObject_muon_candidate_v_plane;
   float  truncated_dQdx_of_TPCObject_muon_candidate_y_plane;
   float  truncated_dQdx_of_TPCObject_muon_candidate_sum;
   float  median_dQdx_of_TPCObject_muon_candidate_u_plane;
   float  median_dQdx_of_TPCObject_muon_candidate_v_plane;
   float  median_dQdx_of_TPCObject_muon_candidate_y_plane;
   float  median_dQdx_of_TPCObject_muon_candidate_sum;
   float  ADC_hits_sum_u_plane;
   float  ADC_hits_sum_v_plane;
   float  ADC_hits_sum_y_plane;
   float  ADC_hits_sum;
   int    hits_num_u_plane;
   int    hits_num_v_plane;
   int    hits_num_y_plane;
   int    hits_num;
   float  ADC_hits_vicinity_sum_u_plane;
   float  ADC_hits_vicinity_sum_v_plane;
   float  ADC_hits_vicinity_sum_y_plane;
   float  ADC_hits_vicinity_sum;
   int    hits_vicinity_num_u_plane;
   int    hits_vicinity_num_v_plane;
   int    hits_vicinity_num_y_plane;
   int    hits_vicinity_num_total;
   int    num_of_top_pandora_crossings;
   int    num_of_bottom_pandora_crossings;
   int    num_of_front_pandora_crossings;
   int    num_of_back_pandora_crossings;
   int    number_of_tracks_in_TPCObject;
   int    num_of_tracks_originating_from_vertex;
   float  muon_candidate_length;
   float  sum_of_TPCObject_track_lengths;
   float  total_length_of_tracks_originating_from_vertex;
   
   TFile* file  = new TFile("NuMI_Background_Recombination_Proton_Fitting_Output_W_All_MCWeights_No_Bad_Events.root", "RECREATE");
   TTree* ftree = new TTree("proton_energy_tree", "Tree With Proton Energy/Direction Parameters");
   ftree->Branch("run", &run, "run/I");
   ftree->Branch("subrun", &subrun, "subrun/I");
   ftree->Branch("event", &event, "event/I");
   ftree->Branch("truth_neutrino_energy", &truth_neutrino_energy, "truth_neutrino_energy/D");
   ftree->Branch("truth_total_proton_KE", &truth_total_proton_KE, "truth_total_proton_KE/D");
   ftree->Branch("truth_leading_proton_KE", &truth_leading_proton_KE, "truth_leading_proton_KE/D");
   ftree->Branch("truth_reco_vtx_distance", &truth_reco_vtx_distance, "truth_reco_vtx_distance/D");
   ftree->Branch("kdar_from_dump_event", &kdar_from_dump_event, "kdar_from_dump_event/I");
   ftree->Branch("NC_channel", &NC_channel, "NC_channel/I");
   ftree->Branch("num_muminus_tracks", &num_muminus_tracks, "num_muminus_tracks/I");
   ftree->Branch("num_muplus_tracks", &num_muplus_tracks, "num_muplus_tracks/I");
   ftree->Branch("num_piplus_tracks", &num_piplus_tracks, "num_piplus_tracks/I");
   ftree->Branch("num_piminus_tracks", &num_piminus_tracks, "num_piminus_tracks/I");
   ftree->Branch("num_pi0_tracks", &num_pi0_tracks, "num_pi0_tracks/I");
   ftree->Branch("num_proton_tracks", &num_proton_tracks, "num_proton_tracks/I");
   ftree->Branch("num_electron_showers", &num_electron_showers, "num_electron_showers/I");
   ftree->Branch("num_positron_showers", &num_positron_showers, "num_positron_showers/I");
   ftree->Branch("num_photon_showers", &num_photon_showers, "num_photon_showers/I");
   ftree->Branch("spline_fix_mcweight", &spline_fix_mcweight, "spline_fix_mcweight/D");
   ftree->Branch("rootino_fix_mcweight", &rootino_fix_mcweight, "rootino_fix_mcweight/D");
   ftree->Branch("central_value_mcweight", &central_value_mcweight, "central_value_mcweight/D");
   ftree->Branch("other_universe_mcweights", &other_universe_mcweights, "other_universe_mcweights[100]/D");
   ftree->Branch("axialff_mcweights", &axialff_mcweights, "axialff_mcweights[2]/D");
   ftree->Branch("rpaccqe_mcweights", &rpaccqe_mcweights, "rpaccqe_mcweights[2]/D");
   ftree->Branch("xsecshape_mcweights", &xsecshape_mcweights, "xsecshape_mcweights[2]/D");
   ftree->Branch("reco_vtx_x", &reco_vtx_x, "reco_vtx_x/F");
   ftree->Branch("reco_vtx_y", &reco_vtx_y, "reco_vtx_y/F");
   ftree->Branch("reco_vtx_z", &reco_vtx_z, "reco_vtx_z/F");
   ftree->Branch("muon_track_is_correctly_oriented", &muon_track_is_correctly_oriented, "muon_track_is_correctly_oriented/I");
   ftree->Branch("x_slope", &x_slope, "x_slope/D");
   ftree->Branch("x_intercept", &x_intercept, "x_intercept/D");
   ftree->Branch("y_slope", &y_slope, "y_slope/D");
   ftree->Branch("y_intercept", &y_intercept, "y_intercept/D");
   ftree->Branch("x_slope_err", &x_slope_err, "x_slope_err/D");
   ftree->Branch("x_intercept_err", &x_intercept_err, "x_intercept_err/D");
   ftree->Branch("y_slope_err", &y_slope_err, "y_slope_err/D");
   ftree->Branch("y_intercept_err", &y_intercept_err, "y_intercept_err/D");
   ftree->Branch("distance_between_furthest_calo_points", &distance_between_furthest_calo_points, "distance_between_furthest_calo_points/D");
   ftree->Branch("distance_using_track_projection", &distance_using_track_projection, "distance_using_track_projection/D");
   ftree->Branch("muon_track_length", &muon_track_length, "muon_track_length/F");
   ftree->Branch("reco_muon_KE", &reco_muon_KE, "reco_muon_KE/F");
   ftree->Branch("muon_phi_without_orientation", &muon_phi_without_orientation, "muon_phi_without_orientation/F");
   ftree->Branch("muon_cos_theta_without_orientation", &muon_cos_theta_without_orientation, "muon_cos_theta_without_orientation/F");
   ftree->Branch("muon_x_momentum", &muon_x_momentum, "muon_x_momentum/F");
   ftree->Branch("muon_y_momentum", &muon_y_momentum, "muon_y_momentum/F");
   ftree->Branch("muon_z_momentum", &muon_z_momentum, "muon_z_momentum/F");
   ftree->Branch("proton_kinetic_energy_from_additional_tracks", &proton_kinetic_energy_from_additional_tracks, "proton_kinetic_energy_from_additional_tracks/F");
   ftree->Branch("x_momentum_component_sum_additional_tracks", &x_momentum_component_sum_additional_tracks, "x_momentum_component_sum_additional_tracks/F");
   ftree->Branch("y_momentum_component_sum_additional_tracks", &y_momentum_component_sum_additional_tracks, "y_momentum_component_sum_additional_tracks/F");
   ftree->Branch("z_momentum_component_sum_additional_tracks", &z_momentum_component_sum_additional_tracks, "z_momentum_component_sum_additional_tracks/F");
   ftree->Branch("proton_x_projection", &proton_x_projection, "proton_x_projection/F");
   ftree->Branch("proton_y_projection", &proton_y_projection, "proton_y_projection/F");
   ftree->Branch("proton_z_projection", &proton_z_projection, "proton_z_projection/F");
   ftree->Branch("x_component_of_muon_momentum", &x_component_of_muon_momentum, "x_component_of_muon_momentum/F");
   ftree->Branch("y_component_of_muon_momentum", &y_component_of_muon_momentum, "y_component_of_muon_momentum/F");
   ftree->Branch("z_component_of_muon_momentum", &z_component_of_muon_momentum, "z_component_of_muon_momentum/F");
   ftree->Branch("num_calo_points_used_in_fit", &num_calo_points_used_in_fit, "num_calo_points_used_in_fit/I");
   ftree->Branch("num_beam_flashes", &num_beam_flashes, "num_beam_flashes/F");
   ftree->Branch("flash_PEs", &flash_PEs, "flash_PEs/F");
   ftree->Branch("flash_z", &flash_z, "flash_z/F");
   ftree->Branch("flash_y", &flash_y, "flash_y/F");
   ftree->Branch("truncated_dQdx_of_TPCObject_muon_candidate_u_plane", &truncated_dQdx_of_TPCObject_muon_candidate_u_plane, "truncated_dQdx_of_TPCObject_muon_candidate_u_plane/F");
   ftree->Branch("truncated_dQdx_of_TPCObject_muon_candidate_v_plane", &truncated_dQdx_of_TPCObject_muon_candidate_v_plane, "truncated_dQdx_of_TPCObject_muon_candidate_v_plane/F");
   ftree->Branch("truncated_dQdx_of_TPCObject_muon_candidate_y_plane", &truncated_dQdx_of_TPCObject_muon_candidate_y_plane, "truncated_dQdx_of_TPCObject_muon_candidate_y_plane/F");
   ftree->Branch("truncated_dQdx_of_TPCObject_muon_candidate_sum", &truncated_dQdx_of_TPCObject_muon_candidate_sum, "truncated_dQdx_of_TPCObject_muon_candidate_sum/F");
   ftree->Branch("median_dQdx_of_TPCObject_muon_candidate_u_plane", &median_dQdx_of_TPCObject_muon_candidate_u_plane, "median_dQdx_of_TPCObject_muon_candidate_u_plane/F");
   ftree->Branch("median_dQdx_of_TPCObject_muon_candidate_v_plane", &median_dQdx_of_TPCObject_muon_candidate_v_plane, "median_dQdx_of_TPCObject_muon_candidate_v_plane/F");
   ftree->Branch("median_dQdx_of_TPCObject_muon_candidate_y_plane", &median_dQdx_of_TPCObject_muon_candidate_y_plane, "median_dQdx_of_TPCObject_muon_candidate_y_plane/F");
   ftree->Branch("median_dQdx_of_TPCObject_muon_candidate_sum", &median_dQdx_of_TPCObject_muon_candidate_sum, "median_dQdx_of_TPCObject_muon_candidate_sum/F");
   ftree->Branch("ADC_hits_sum_u_plane", &ADC_hits_sum_u_plane, "ADC_hits_sum_u_plane/F");
   ftree->Branch("ADC_hits_sum_v_plane", &ADC_hits_sum_v_plane, "ADC_hits_sum_v_plane/F");
   ftree->Branch("ADC_hits_sum_y_plane", &ADC_hits_sum_y_plane, "ADC_hits_sum_y_plane/F");
   ftree->Branch("ADC_hits_sum", &ADC_hits_sum, "ADC_hits_sum/F");
   ftree->Branch("hits_num_u_plane", &hits_num_u_plane, "hits_num_u_plane/I");
   ftree->Branch("hits_num_v_plane", &hits_num_v_plane, "hits_num_v_plane/I");
   ftree->Branch("hits_num_y_plane", &hits_num_y_plane, "hits_num_y_plane/I");
   ftree->Branch("hits_num", &hits_num, "hits_num/I");
   ftree->Branch("ADC_hits_vicinity_sum_u_plane", &ADC_hits_vicinity_sum_u_plane, "ADC_hits_vicinity_sum_u_plane/F");
   ftree->Branch("ADC_hits_vicinity_sum_v_plane", &ADC_hits_vicinity_sum_v_plane, "ADC_hits_vicinity_sum_v_plane/F");
   ftree->Branch("ADC_hits_vicinity_sum_y_plane", &ADC_hits_vicinity_sum_y_plane, "ADC_hits_vicinity_sum_y_plane/F");
   ftree->Branch("ADC_hits_vicinity_sum", &ADC_hits_vicinity_sum, "ADC_hits_vicinity_sum/F");
   ftree->Branch("hits_vicinity_num_u_plane", &hits_vicinity_num_u_plane, "hits_vicinity_num_u_plane/I");
   ftree->Branch("hits_vicinity_num_v_plane", &hits_vicinity_num_v_plane, "hits_vicinity_num_v_plane/I");
   ftree->Branch("hits_vicinity_num_y_plane", &hits_vicinity_num_y_plane, "hits_vicinity_num_y_plane/I");
   ftree->Branch("hits_vicinity_num_total", &hits_vicinity_num_total, "hits_vicinity_num_total/I");
   ftree->Branch("num_of_top_pandora_crossings", &num_of_top_pandora_crossings, "num_of_top_pandora_crossings/I");
   ftree->Branch("num_of_bottom_pandora_crossings", &num_of_bottom_pandora_crossings, "num_of_bottom_pandora_crossings/I");
   ftree->Branch("num_of_front_pandora_crossings", &num_of_front_pandora_crossings, "num_of_front_pandora_crossings/I");
   ftree->Branch("num_of_back_pandora_crossings", &num_of_back_pandora_crossings, "num_of_back_pandora_crossings/I");
   ftree->Branch("number_of_tracks_in_TPCObject", &number_of_tracks_in_TPCObject, "number_of_tracks_in_TPCObject/I");
   ftree->Branch("num_of_tracks_originating_from_vertex", &num_of_tracks_originating_from_vertex, "num_of_tracks_originating_from_vertex/I");
   ftree->Branch("muon_candidate_length", &muon_candidate_length, "muon_candidate_length/F");
   ftree->Branch("sum_of_TPCObject_track_lengths", &sum_of_TPCObject_track_lengths, "sum_of_TPCObject_track_lengths/F");
   ftree->Branch("total_length_of_tracks_originating_from_vertex", &total_length_of_tracks_originating_from_vertex, "total_length_of_tracks_originating_from_vertex/F");   
   
   int    run_input;
   int    subrun_input;
   int    event_input;
   double truth_neutrino_energy_input;
   double truth_total_proton_KE_input;
   double truth_leading_proton_KE_input;
   double mcweight_input;
   int    num_calo_points_above_threshold;
   double x_points_above_threshold[10000];
   double y_points_above_threshold[10000];
   double z_points_above_threshold[10000];
   double reconstructed_muon_KE;
   double muon_track_first_point_x;
   double muon_track_first_point_y;
   double muon_track_first_point_z;
   double muon_track_last_point_x;
   double muon_track_last_point_y;
   double muon_track_last_point_z;
   double vertex_location_x;
   double vertex_location_y;
   double vertex_location_z;
   double other_end_location_x;
   double other_end_location_y;
   double other_end_location_z;
   double input_proton_kinetic_energy_from_additional_tracks;
   double input_x_momentum_component_sum_additional_tracks;
   double input_y_momentum_component_sum_additional_tracks;
   double input_z_momentum_component_sum_additional_tracks;
   int    proton_x_direction;
   int    proton_y_direction;
   int    proton_z_direction;
   int    muon_track_is_correctly_oriented_input;
   
   // Load in all of the necessary quantities from the input file.
   TChain* tree = new TChain("UBXSec/proton_hit_identification_tree");
   tree->Add("/home/barnchri/Performing_Analysis_With_MCWeights/Recombination_Variation_Filter_Output_No_Bad_Events.root");
   tree->SetBranchAddress("run", &run_input);
   tree->SetBranchAddress("subrun", &subrun_input);
   tree->SetBranchAddress("event", &event_input);
   tree->SetBranchAddress("neutrino_energy", &truth_neutrino_energy_input);
   tree->SetBranchAddress("total_proton_KE", &truth_total_proton_KE_input);
   tree->SetBranchAddress("leading_proton_KE", &truth_leading_proton_KE_input);
   tree->SetBranchAddress("num_points_above_threshold", &num_calo_points_above_threshold);
   tree->SetBranchAddress("x_points_above_threshold", &x_points_above_threshold);
   tree->SetBranchAddress("y_points_above_threshold", &y_points_above_threshold);
   tree->SetBranchAddress("z_points_above_threshold", &z_points_above_threshold);
   tree->SetBranchAddress("reconstructed_muon_KE", &reconstructed_muon_KE);
   tree->SetBranchAddress("muon_track_first_point_x", &muon_track_first_point_x);
   tree->SetBranchAddress("muon_track_first_point_y", &muon_track_first_point_y);
   tree->SetBranchAddress("muon_track_first_point_z", &muon_track_first_point_z);
   tree->SetBranchAddress("muon_track_last_point_x", &muon_track_last_point_x);
   tree->SetBranchAddress("muon_track_last_point_y", &muon_track_last_point_y);
   tree->SetBranchAddress("muon_track_last_point_z", &muon_track_last_point_z);
   tree->SetBranchAddress("vertex_location_x", &vertex_location_x);
   tree->SetBranchAddress("vertex_location_y", &vertex_location_y);
   tree->SetBranchAddress("vertex_location_z", &vertex_location_z);
   tree->SetBranchAddress("other_end_location_x", &other_end_location_x);
   tree->SetBranchAddress("other_end_location_y", &other_end_location_y);
   tree->SetBranchAddress("other_end_location_z", &other_end_location_z);
   tree->SetBranchAddress("proton_kinetic_energy_from_additional_tracks", &input_proton_kinetic_energy_from_additional_tracks);
   tree->SetBranchAddress("x_momentum_component_sum_additional_tracks", &input_x_momentum_component_sum_additional_tracks);
   tree->SetBranchAddress("y_momentum_component_sum_additional_tracks", &input_y_momentum_component_sum_additional_tracks);
   tree->SetBranchAddress("z_momentum_component_sum_additional_tracks", &input_z_momentum_component_sum_additional_tracks);
   tree->SetBranchAddress("proton_x_direction", &proton_x_direction);
   tree->SetBranchAddress("proton_y_direction", &proton_y_direction);
   tree->SetBranchAddress("proton_z_direction", &proton_z_direction);
   tree->SetBranchAddress("muon_track_is_correctly_oriented", &muon_track_is_correctly_oriented_input);
   
   int    passing_events_run;
   int    passing_events_subrun;
   int    passing_events_event;
   int    kdar_from_dump_event_input;
   int    NC_channel_input;
   double spline_fix_mcweight_input;
   double central_value_mcweight_input;
   double rootino_fix_mcweight_input;
   double other_universe_mcweights_input[100];
   double axialff_mcweights_input[2];
   double rpaccqe_mcweights_input[2];
   double xsecshape_mcweights_input[2];
   int    num_muminus_tracks_input;
   int    num_muplus_tracks_input;
   int    num_piplus_tracks_input;
   int    num_piminus_tracks_input;
   int    num_pi0_tracks_input;
   int    num_proton_tracks_input;
   int    num_electron_showers_input;
   int    num_positron_showers_input;
   int    num_photon_showers_input;
   double truth_reco_vtx_distance_input;
   double muon_candidate_track_length;
   double ubxsec_muon_phi;
   double ubxsec_muon_cos_theta;
   int    num_beam_flashes_input;
   double flash_PEs_input;
   double flash_z_input;
   double flash_y_input;
   double truncated_dQdx_of_TPCObject_muon_candidate_u_plane_input;
   double truncated_dQdx_of_TPCObject_muon_candidate_v_plane_input;
   double truncated_dQdx_of_TPCObject_muon_candidate_y_plane_input;
   double median_dQdx_of_TPCObject_muon_candidate_u_plane_input;
   double median_dQdx_of_TPCObject_muon_candidate_v_plane_input;
   double median_dQdx_of_TPCObject_muon_candidate_y_plane_input;
   double u_plane_sum_of_slice_associated_hits_ADCs;
   double v_plane_sum_of_slice_associated_hits_ADCs;
   double y_plane_sum_of_slice_associated_hits_ADCs;
   double total_sum_of_slice_associated_hits_ADCs;
   int    u_plane_num_of_slice_associated_hits;
   int    v_plane_num_of_slice_associated_hits;
   int    y_plane_num_of_slice_associated_hits;
   double total_sum_of_slice_hits_in_vicinity_ADCs;
   double u_plane_sum_of_slice_hits_in_vicinity_ADCs;
   double v_plane_sum_of_slice_hits_in_vicinity_ADCs;
   double y_plane_sum_of_slice_hits_in_vicinity_ADCs;
   int    total_num_of_slice_associated_hits;
   int    u_plane_num_of_slice_hits_in_vicinity;
   int    v_plane_num_of_slice_hits_in_vicinity;
   int    y_plane_num_of_slice_hits_in_vicinity;
   int    total_num_of_slice_hits_in_vicinity;
   int    num_of_top_pandora_crossings_input;
   int    num_of_bottom_pandora_crossings_input;
   int    num_of_front_pandora_crossings_input;
   int    num_of_back_pandora_crossings_input;
   int    number_of_tracks_in_TPCObject_input;
   int    num_of_tracks_originating_from_vertex_input;
   double muon_candidate_length_input;
   double sum_of_TPCObject_track_lengths_input;
   
   
   //     For the cuts that we have to apply afterwards.
   double total_length_of_tracks_originating_from_vertex_input;
   int    event_fails_new_two_track_selection_input;
   
   TChain* passing_events_tree = new TChain("UBXSec/_passing_events_tree");
   passing_events_tree->Add("/home/barnchri/Performing_Analysis_With_MCWeights/Recombination_Variation_Filter_Output_No_Bad_Events.root");
   passing_events_tree->SetBranchAddress("run", &passing_events_run);
   passing_events_tree->SetBranchAddress("subrun", &passing_events_subrun);
   passing_events_tree->SetBranchAddress("event", &passing_events_event);
   passing_events_tree->SetBranchAddress("kdar_from_dump_event", &kdar_from_dump_event_input);
   passing_events_tree->SetBranchAddress("NC_channel", &NC_channel_input);
   passing_events_tree->SetBranchAddress("spline_fix_mcweight", &spline_fix_mcweight_input);
   passing_events_tree->SetBranchAddress("rootino_fix_mcweight", &rootino_fix_mcweight_input);
   passing_events_tree->SetBranchAddress("central_value_mcweight", &central_value_mcweight_input);
   passing_events_tree->SetBranchAddress("other_universe_mcweights", &other_universe_mcweights_input);
   passing_events_tree->SetBranchAddress("axialff_mcweights", &axialff_mcweights_input);
   passing_events_tree->SetBranchAddress("rpaccqe_mcweights", &rpaccqe_mcweights_input);
   passing_events_tree->SetBranchAddress("xsecshape_mcweights", &xsecshape_mcweights_input);
   passing_events_tree->SetBranchAddress("num_muminus_tracks", &num_muminus_tracks_input);
   passing_events_tree->SetBranchAddress("num_muplus_tracks", &num_muplus_tracks_input);
   passing_events_tree->SetBranchAddress("num_piplus_tracks", &num_piplus_tracks_input);
   passing_events_tree->SetBranchAddress("num_piminus_tracks", &num_piminus_tracks_input);
   passing_events_tree->SetBranchAddress("num_pi0_tracks", &num_pi0_tracks_input);
   passing_events_tree->SetBranchAddress("num_proton_tracks", &num_proton_tracks_input);
   passing_events_tree->SetBranchAddress("num_electron_showers", &num_electron_showers_input);
   passing_events_tree->SetBranchAddress("num_positron_showers", &num_positron_showers_input);
   passing_events_tree->SetBranchAddress("num_photon_showers", &num_photon_showers_input);
   passing_events_tree->SetBranchAddress("truth_reco_vtx_distance", &truth_reco_vtx_distance_input );
   passing_events_tree->SetBranchAddress("muon_candidate_length", &muon_candidate_track_length);
   passing_events_tree->SetBranchAddress("ubxsec_muon_phi", &ubxsec_muon_phi);
   passing_events_tree->SetBranchAddress("ubxsec_muon_cos_theta", &ubxsec_muon_cos_theta);
   passing_events_tree->SetBranchAddress("num_beam_flashes", &num_beam_flashes_input);
   passing_events_tree->SetBranchAddress("flash_PEs", &flash_PEs_input);
   passing_events_tree->SetBranchAddress("flash_z", &flash_z_input);
   passing_events_tree->SetBranchAddress("flash_y", &flash_y_input);
   passing_events_tree->SetBranchAddress("truncated_dQdx_of_TPCObject_muon_candidate_u_plane", &truncated_dQdx_of_TPCObject_muon_candidate_u_plane_input);
   passing_events_tree->SetBranchAddress("truncated_dQdx_of_TPCObject_muon_candidate_v_plane", &truncated_dQdx_of_TPCObject_muon_candidate_v_plane_input);
   passing_events_tree->SetBranchAddress("truncated_dQdx_of_TPCObject_muon_candidate_y_plane", &truncated_dQdx_of_TPCObject_muon_candidate_y_plane_input);
   passing_events_tree->SetBranchAddress("median_dQdx_of_TPCObject_muon_candidate_u_plane", &median_dQdx_of_TPCObject_muon_candidate_u_plane_input);
   passing_events_tree->SetBranchAddress("median_dQdx_of_TPCObject_muon_candidate_v_plane", &median_dQdx_of_TPCObject_muon_candidate_v_plane_input);
   passing_events_tree->SetBranchAddress("median_dQdx_of_TPCObject_muon_candidate_y_plane", &median_dQdx_of_TPCObject_muon_candidate_y_plane_input);
   passing_events_tree->SetBranchAddress("u_plane_sum_of_slice_associated_hits_ADCs", &u_plane_sum_of_slice_associated_hits_ADCs);
   passing_events_tree->SetBranchAddress("v_plane_sum_of_slice_associated_hits_ADCs", &v_plane_sum_of_slice_associated_hits_ADCs);
   passing_events_tree->SetBranchAddress("y_plane_sum_of_slice_associated_hits_ADCs", &y_plane_sum_of_slice_associated_hits_ADCs);
   passing_events_tree->SetBranchAddress("total_sum_of_slice_associated_hits_ADCs", &total_sum_of_slice_associated_hits_ADCs);
   passing_events_tree->SetBranchAddress("u_plane_number_of_slice_associated_hits", &u_plane_num_of_slice_associated_hits);
   passing_events_tree->SetBranchAddress("v_plane_number_of_slice_associated_hits", &v_plane_num_of_slice_associated_hits);
   passing_events_tree->SetBranchAddress("y_plane_number_of_slice_associated_hits", &y_plane_num_of_slice_associated_hits);
   passing_events_tree->SetBranchAddress("total_number_of_slice_associated_hits", &total_num_of_slice_associated_hits);
   passing_events_tree->SetBranchAddress("u_plane_sum_of_slice_hits_in_vicinity_ADCs", &u_plane_sum_of_slice_hits_in_vicinity_ADCs);
   passing_events_tree->SetBranchAddress("v_plane_sum_of_slice_hits_in_vicinity_ADCs", &v_plane_sum_of_slice_hits_in_vicinity_ADCs);
   passing_events_tree->SetBranchAddress("y_plane_sum_of_slice_hits_in_vicinity_ADCs", &y_plane_sum_of_slice_hits_in_vicinity_ADCs);
   passing_events_tree->SetBranchAddress("total_sum_of_slice_hits_in_vicinity_ADCs", &total_sum_of_slice_hits_in_vicinity_ADCs);
   passing_events_tree->SetBranchAddress("u_plane_number_of_slice_hits_in_vicinity", &u_plane_num_of_slice_hits_in_vicinity);
   passing_events_tree->SetBranchAddress("v_plane_number_of_slice_hits_in_vicinity", &v_plane_num_of_slice_hits_in_vicinity);
   passing_events_tree->SetBranchAddress("y_plane_number_of_slice_hits_in_vicinity", &y_plane_num_of_slice_hits_in_vicinity);
   passing_events_tree->SetBranchAddress("total_number_of_slice_hits_in_vicinity", &total_num_of_slice_hits_in_vicinity);
   passing_events_tree->SetBranchAddress("num_of_top_pandora_crossings", &num_of_top_pandora_crossings_input);
   passing_events_tree->SetBranchAddress("num_of_bottom_pandora_crossings", &num_of_bottom_pandora_crossings_input);
   passing_events_tree->SetBranchAddress("num_of_front_pandora_crossings", &num_of_front_pandora_crossings_input);
   passing_events_tree->SetBranchAddress("num_of_back_pandora_crossings", &num_of_back_pandora_crossings_input);
   passing_events_tree->SetBranchAddress("total_length_of_tracks_originating_from_vertex", &total_length_of_tracks_originating_from_vertex_input);
   passing_events_tree->SetBranchAddress("event_fails_new_two_track_selection", &event_fails_new_two_track_selection_input);
   passing_events_tree->SetBranchAddress("number_of_tracks_in_TPCObject", &number_of_tracks_in_TPCObject_input);
   passing_events_tree->SetBranchAddress("num_of_tracks_originating_from_vertex", &num_of_tracks_originating_from_vertex_input);
   passing_events_tree->SetBranchAddress("muon_candidate_length", &muon_candidate_length_input);
   passing_events_tree->SetBranchAddress("sum_of_TPCObject_track_lengths", &sum_of_TPCObject_track_lengths_input);
   passing_events_tree->SetBranchAddress("total_length_of_tracks_originating_from_vertex", &total_length_of_tracks_originating_from_vertex_input);
   
   int num_entries                = tree->GetEntries();

   int passing_events_num_entries = passing_events_tree->GetEntries();

   bool event_fails_additional_cuts = false;
   int  num_entries_passing = 0;
   
   for ( size_t i = 0; i < num_entries; i++ ) {

     tree->GetEntry( i );

     event_fails_additional_cuts = false;
     
     for ( size_t j = 0; j < passing_events_num_entries; j++ ) {
       
       passing_events_tree->GetEntry( j );

       if ( passing_events_run == run_input && passing_events_subrun == subrun_input && passing_events_event == event_input ) {
       
	 break;
	   
       }

     }

     // Skip the event if it doesn't pass the selection requirements.
     if ( total_length_of_tracks_originating_from_vertex_input > 40.0 || event_fails_new_two_track_selection_input == 1 ) continue;

     num_entries_passing++;

     run                                          = run_input;
     subrun                                       = subrun_input;
     event                                        = event_input;
     truth_neutrino_energy                        = truth_neutrino_energy_input;
     truth_total_proton_KE                        = truth_total_proton_KE_input;
     truth_leading_proton_KE                      = truth_leading_proton_KE_input;
     truth_reco_vtx_distance                      = truth_reco_vtx_distance_input;
     kdar_from_dump_event                         = kdar_from_dump_event_input;
     spline_fix_mcweight                          = spline_fix_mcweight_input;
     rootino_fix_mcweight                         = rootino_fix_mcweight_input;
     central_value_mcweight                       = central_value_mcweight_input;
     NC_channel                                   = NC_channel_input;
     num_muminus_tracks                           = num_muminus_tracks_input;
     num_muplus_tracks                            = num_muplus_tracks_input;
     num_piplus_tracks                            = num_piplus_tracks_input;
     num_piminus_tracks                           = num_piminus_tracks_input;
     num_pi0_tracks                               = num_pi0_tracks_input;
     num_proton_tracks                            = num_proton_tracks_input;
     num_electron_showers                         = num_electron_showers_input;
     num_positron_showers                         = num_positron_showers_input;
     num_photon_showers                           = num_photon_showers_input;
     truth_reco_vtx_distance                      = truth_reco_vtx_distance_input;
     reco_vtx_x                                   = vertex_location_x;
     reco_vtx_y                                   = vertex_location_y;
     reco_vtx_z                                   = vertex_location_z;
     muon_track_is_correctly_oriented             = muon_track_is_correctly_oriented_input;
     num_calo_points_used_in_fit                  = num_calo_points_above_threshold;
     muon_track_length                            = muon_candidate_track_length;
     reco_muon_KE                                 = reconstructed_muon_KE;
     muon_phi_without_orientation                 = ubxsec_muon_phi;
     muon_cos_theta_without_orientation           = ubxsec_muon_cos_theta;
     proton_kinetic_energy_from_additional_tracks = input_proton_kinetic_energy_from_additional_tracks;
     x_momentum_component_sum_additional_tracks   = input_x_momentum_component_sum_additional_tracks;
     y_momentum_component_sum_additional_tracks   = input_y_momentum_component_sum_additional_tracks;
     z_momentum_component_sum_additional_tracks   = input_z_momentum_component_sum_additional_tracks;
     num_beam_flashes                             = num_beam_flashes_input;
     flash_PEs                                    = flash_PEs_input;
     flash_z                                      = flash_z_input;
     flash_y                                      = flash_y_input;
     ADC_hits_sum_u_plane                         = u_plane_sum_of_slice_associated_hits_ADCs;
     ADC_hits_sum_v_plane                         = v_plane_sum_of_slice_associated_hits_ADCs;
     ADC_hits_sum_y_plane                         = y_plane_sum_of_slice_associated_hits_ADCs;
     ADC_hits_sum                                 = total_sum_of_slice_associated_hits_ADCs;
     hits_num_u_plane                             = u_plane_num_of_slice_associated_hits;
     hits_num_v_plane                             = v_plane_num_of_slice_associated_hits;
     hits_num_y_plane                             = y_plane_num_of_slice_associated_hits;
     hits_num                                     = total_num_of_slice_associated_hits;
     ADC_hits_vicinity_sum_u_plane                = u_plane_sum_of_slice_hits_in_vicinity_ADCs;
     ADC_hits_vicinity_sum_v_plane                = v_plane_sum_of_slice_hits_in_vicinity_ADCs;
     ADC_hits_vicinity_sum_y_plane                = y_plane_sum_of_slice_hits_in_vicinity_ADCs;
     ADC_hits_vicinity_sum                        = total_sum_of_slice_hits_in_vicinity_ADCs;
     hits_vicinity_num_u_plane                    = u_plane_num_of_slice_hits_in_vicinity;
     hits_vicinity_num_v_plane                    = v_plane_num_of_slice_hits_in_vicinity;
     hits_vicinity_num_y_plane                    = y_plane_num_of_slice_hits_in_vicinity;
     hits_vicinity_num_total                      = total_num_of_slice_hits_in_vicinity;
     num_of_top_pandora_crossings                 = num_of_top_pandora_crossings_input;
     num_of_bottom_pandora_crossings              = num_of_bottom_pandora_crossings_input;
     num_of_front_pandora_crossings               = num_of_front_pandora_crossings_input;
     num_of_back_pandora_crossings                = num_of_back_pandora_crossings_input;
     number_of_tracks_in_TPCObject                = number_of_tracks_in_TPCObject_input;
     num_of_tracks_originating_from_vertex        = num_of_tracks_originating_from_vertex_input;
     muon_candidate_length                        = muon_candidate_length_input;
     sum_of_TPCObject_track_lengths               = sum_of_TPCObject_track_lengths_input;
     total_length_of_tracks_originating_from_vertex = total_length_of_tracks_originating_from_vertex_input;
     
     for ( size_t i = 0; i < 100; i++ ) {

       other_universe_mcweights[ i ] = other_universe_mcweights_input[ i ];

       std::cout << "Other Universe Input #" << i << " = " << other_universe_mcweights_input[ i ] << "." << std::endl;

     }

     for ( size_t i = 0; i < 2; i++ ) {

       axialff_mcweights[ i ] = axialff_mcweights_input[ i ];

       std::cout << "AxialFF Input #" << i << " = " << axialff_mcweights_input[ i ] << "." << std::endl;

     }

     for ( size_t i = 0; i < 2; i++ ) {

       rpaccqe_mcweights[ i ] = rpaccqe_mcweights_input[ i ];

       std::cout << "RPACCQE Input #" << i << " = " << rpaccqe_mcweights_input[ i ] << "." << std::endl;
       
     }

     for ( size_t i = 0; i < 2; i++ ) {

       xsecshape_mcweights[ i ] = xsecshape_mcweights_input[ i ];

       std::cout << "XSecShape Input #" << i << " = " << xsecshape_mcweights_input[ i ] << "." << std::endl;
       
     }
     
     truncated_dQdx_of_TPCObject_muon_candidate_u_plane = truncated_dQdx_of_TPCObject_muon_candidate_u_plane_input;
     truncated_dQdx_of_TPCObject_muon_candidate_v_plane = truncated_dQdx_of_TPCObject_muon_candidate_v_plane_input;
     truncated_dQdx_of_TPCObject_muon_candidate_y_plane = truncated_dQdx_of_TPCObject_muon_candidate_y_plane_input;

     if	( truncated_dQdx_of_TPCObject_muon_candidate_u_plane < 0.0 )
       truncated_dQdx_of_TPCObject_muon_candidate_u_plane = 0.0;
     
     if	( truncated_dQdx_of_TPCObject_muon_candidate_v_plane < 0.0 )
       truncated_dQdx_of_TPCObject_muon_candidate_v_plane = 0.0;
     
     if ( truncated_dQdx_of_TPCObject_muon_candidate_y_plane < 0.0 )
       truncated_dQdx_of_TPCObject_muon_candidate_y_plane = 0.0;

     truncated_dQdx_of_TPCObject_muon_candidate_sum = ( truncated_dQdx_of_TPCObject_muon_candidate_u_plane + truncated_dQdx_of_TPCObject_muon_candidate_v_plane + truncated_dQdx_of_TPCObject_muon_candidate_y_plane );

     median_dQdx_of_TPCObject_muon_candidate_u_plane = median_dQdx_of_TPCObject_muon_candidate_u_plane_input;
     median_dQdx_of_TPCObject_muon_candidate_v_plane = median_dQdx_of_TPCObject_muon_candidate_v_plane_input;
     median_dQdx_of_TPCObject_muon_candidate_y_plane = median_dQdx_of_TPCObject_muon_candidate_y_plane_input;

     if ( median_dQdx_of_TPCObject_muon_candidate_u_plane < 0.0 )
       median_dQdx_of_TPCObject_muon_candidate_u_plane  = 0.0;

     if ( median_dQdx_of_TPCObject_muon_candidate_v_plane < 0.0 )
       median_dQdx_of_TPCObject_muon_candidate_v_plane  = 0.0;

     if ( median_dQdx_of_TPCObject_muon_candidate_y_plane < 0.0 )
       median_dQdx_of_TPCObject_muon_candidate_y_plane  = 0.0;

     median_dQdx_of_TPCObject_muon_candidate_sum = ( median_dQdx_of_TPCObject_muon_candidate_u_plane + median_dQdx_of_TPCObject_muon_candidate_v_plane + median_dQdx_of_TPCObject_muon_candidate_y_plane );
     
     // Calculate the y-component of momentum.
     x_component_of_muon_momentum                      = fabs( muon_track_last_point_x - muon_track_first_point_x ) / TMath::Sqrt( ( muon_track_last_point_x - muon_track_first_point_x ) * ( muon_track_last_point_x - muon_track_first_point_x ) + ( muon_track_last_point_y - muon_track_first_point_y ) * ( muon_track_last_point_y - muon_track_first_point_y ) + ( muon_track_last_point_z - muon_track_first_point_z) * ( muon_track_last_point_z - muon_track_first_point_z ) );
     y_component_of_muon_momentum                      = fabs( muon_track_last_point_y - muon_track_first_point_y ) / TMath::Sqrt( ( muon_track_last_point_x - muon_track_first_point_x ) * ( muon_track_last_point_x - muon_track_first_point_x ) + ( muon_track_last_point_y - muon_track_first_point_y ) * ( muon_track_last_point_y - muon_track_first_point_y ) + ( muon_track_last_point_z - muon_track_first_point_z ) * ( muon_track_last_point_z - muon_track_first_point_z ) );
     z_component_of_muon_momentum                      = fabs( muon_track_last_point_z - muon_track_first_point_z ) / TMath::Sqrt( ( muon_track_last_point_x - muon_track_first_point_x ) * ( muon_track_last_point_x - muon_track_first_point_x ) + ( muon_track_last_point_y - muon_track_first_point_y ) * ( muon_track_last_point_y - muon_track_first_point_y ) + ( muon_track_last_point_z - muon_track_first_point_z) * ( muon_track_last_point_z - muon_track_first_point_z ) );
     
     std::cout << "run = " << run << " subrun = " << subrun << " event = " << event << "." << std::endl;

     double closest_distance_to_outer_point = 10000000.;
     double closest_distance_to_inner_point = 10000000.;

     // Clear out the 'idx_removed_points' vector at the very end.
     bool outliers_found_last_iteration       = false;
     int  num_outliers_removed                = 0;
     distance_between_furthest_calo_points    = -1.;
     std::vector< size_t > idx_removed_points;
     idx_removed_points.clear();

     size_t idx_of_outer_point                = -1;
     size_t idx_of_inner_point                = -1;

     do {

       // Reset 'outliers_found_last_iteration' to false.
       outliers_found_last_iteration = false;

       // Find the two points in 'x_points_above_threshold' located closest to the end of the track (ignoring outliers).
       idx_of_outer_point                     = -1;
       idx_of_inner_point                     = -1;
       distance_between_furthest_calo_points  = -1.;
       
       // Find the maximum distance between two points above threshold.
       for ( size_t outer_point_iter = 0; outer_point_iter < num_calo_points_above_threshold; outer_point_iter++ ) {

	 bool outer_point_is_an_outlier = false;

	 // Skip this point if it is one of the outliers on the track.
	 for ( size_t outer_outlier_point_iter = 0; outer_outlier_point_iter < idx_removed_points.size(); outer_outlier_point_iter++ ) {

	   if ( outer_point_iter == idx_removed_points.at( outer_outlier_point_iter ) ) outer_point_is_an_outlier = true;

	 }

	 if ( outer_point_is_an_outlier == true ) continue;
	 
	 for ( size_t inner_point_iter = ( outer_point_iter + 1 ); inner_point_iter < num_calo_points_above_threshold; inner_point_iter++ ) {

	   bool inner_point_is_an_outlier = false;

	   // Skip this point if it is one of the outliers on the track.
	   for ( size_t inner_outlier_point_iter = 0; inner_outlier_point_iter < idx_removed_points.size(); inner_outlier_point_iter++ ) {

	     if ( inner_point_iter == idx_removed_points.at( inner_outlier_point_iter ) ) inner_point_is_an_outlier = true;
	   
	   }

	   if ( inner_point_is_an_outlier == true ) continue;

	   double distance = TMath::Sqrt( ( x_points_above_threshold[ outer_point_iter ] - x_points_above_threshold[ inner_point_iter ] ) * ( x_points_above_threshold[ outer_point_iter ] - x_points_above_threshold[ inner_point_iter ] ) + ( y_points_above_threshold[ outer_point_iter ] - y_points_above_threshold[ inner_point_iter ] ) * ( y_points_above_threshold[ outer_point_iter ] - y_points_above_threshold[ inner_point_iter ] ) + ( z_points_above_threshold[ outer_point_iter ] - z_points_above_threshold[ inner_point_iter ] ) * ( z_points_above_threshold[ outer_point_iter ] - z_points_above_threshold[ inner_point_iter ] ) );

	   if ( distance > distance_between_furthest_calo_points ) {

	     distance_between_furthest_calo_points = distance;
	     idx_of_outer_point                    = outer_point_iter;
	     idx_of_inner_point                    = inner_point_iter;

	   } // End of replacing the furthest distance if you find the points furthest away.

	 } // End of the loop over the inner points.

       } // End of the loop over the outer points.

       std::cout << "Distance Between The Furthest Points = " << distance_between_furthest_calo_points << " cm." << std::endl;

       std::cout << "Outer Point Coordinates: x = " << x_points_above_threshold[idx_of_outer_point] << " cm y = " << y_points_above_threshold[idx_of_outer_point] << " cm z = " << z_points_above_threshold[idx_of_outer_point] << " cm." << std::endl;
       std::cout << "Inner Point Coordinates: x = " << x_points_above_threshold[idx_of_inner_point] << " cm y = " << y_points_above_threshold[idx_of_inner_point] << " cm z = " << z_points_above_threshold[idx_of_inner_point] << " cm." << std::endl;

       std::cout << "The number of points already tagged as outliers = " << idx_removed_points.size() << "." << std::endl;
       
       // Only look for more outliers if the distance between the outer points is > 3.0 cm.
       if ( distance_between_furthest_calo_points > 3.0 ) {

	 // Find the next closest point to the two outer points.
	 for ( size_t dist_to_last_point_iter = 0; dist_to_last_point_iter < num_calo_points_above_threshold; dist_to_last_point_iter++ ) {

	   // Skip this point if it is an outlier.
	   bool point_is_an_outlier = false;

	   // Skip this point if it is one of the outliers on the track.
	   for ( size_t point_iter = 0; point_iter < idx_removed_points.size(); point_iter++ ) {

	     if ( point_iter == idx_removed_points.at( point_iter ) ) point_is_an_outlier = true;

	   }

	   if ( point_is_an_outlier == true ) continue;

	   // Skip itself and the other point under consideration.
	   if ( dist_to_last_point_iter == idx_of_outer_point || dist_to_last_point_iter == idx_of_inner_point ) continue;

	   double distance = TMath::Sqrt( ( x_points_above_threshold[ dist_to_last_point_iter ] - x_points_above_threshold[ idx_of_outer_point ] ) * ( x_points_above_threshold[ dist_to_last_point_iter ] - x_points_above_threshold[ idx_of_outer_point ] ) + ( y_points_above_threshold[ dist_to_last_point_iter ] - y_points_above_threshold[ idx_of_outer_point ] ) * ( y_points_above_threshold[ dist_to_last_point_iter ] - y_points_above_threshold[ idx_of_outer_point ] ) + ( z_points_above_threshold[ dist_to_last_point_iter ] - z_points_above_threshold[ idx_of_outer_point ] ) * ( z_points_above_threshold[ dist_to_last_point_iter ] - z_points_above_threshold[ idx_of_outer_point ] ) );
	 
	   if ( distance < closest_distance_to_outer_point ) {
	     closest_distance_to_outer_point = distance;
	   }
	   
	 }
       
	 // Find the next closest point to the two inner points.
	 for ( size_t dist_to_first_point_iter = 0; dist_to_first_point_iter < num_calo_points_above_threshold; dist_to_first_point_iter++ ) {
	 
	   // Skip this point if it is an outlier.
	   bool point_is_an_outlier = false;

           // Skip this point if it is one of the outliers on the track.
	   for ( size_t point_iter = 0; point_iter < idx_removed_points.size(); point_iter++ ) {

             if ( point_iter == idx_removed_points.at( point_iter ) ) point_is_an_outlier = true;

           }

           if ( point_is_an_outlier == true ) continue;

           // Skip itself and the other point under consideration.
	   if ( dist_to_first_point_iter == idx_of_outer_point || dist_to_first_point_iter == idx_of_inner_point ) continue;
	 
	   double distance = TMath::Sqrt( ( x_points_above_threshold[ dist_to_first_point_iter ] - x_points_above_threshold[ idx_of_inner_point ] ) * ( x_points_above_threshold[ dist_to_first_point_iter ] - x_points_above_threshold[ idx_of_inner_point ] ) + ( y_points_above_threshold[ dist_to_first_point_iter ] - y_points_above_threshold[ idx_of_inner_point ] ) * ( y_points_above_threshold[ dist_to_first_point_iter ] - y_points_above_threshold[ idx_of_inner_point ] ) + ( z_points_above_threshold[ dist_to_first_point_iter ] - z_points_above_threshold[ idx_of_inner_point ] )	* ( z_points_above_threshold[ dist_to_first_point_iter ] - z_points_above_threshold[ idx_of_inner_point ] ) );

	   if ( distance < closest_distance_to_inner_point ) {
	     closest_distance_to_inner_point = distance;
	   }
	   
	 }

	 double distance_to_start_of_track = TMath::Sqrt( ( x_points_above_threshold[ idx_of_outer_point] - vertex_location_x ) * ( x_points_above_threshold[ idx_of_outer_point] - vertex_location_x ) + ( y_points_above_threshold[ idx_of_outer_point] - vertex_location_y ) * ( y_points_above_threshold[ idx_of_outer_point] - vertex_location_y ) + ( z_points_above_threshold[ idx_of_outer_point] - vertex_location_z ) * ( z_points_above_threshold[ idx_of_outer_point] - vertex_location_z ) );
	 double distance_to_end_of_track = TMath::Sqrt( ( x_points_above_threshold[ idx_of_outer_point] - other_end_location_x ) * ( x_points_above_threshold[ idx_of_outer_point] - other_end_location_x ) + (y_points_above_threshold[ idx_of_outer_point] - other_end_location_y ) * ( y_points_above_threshold[ idx_of_outer_point] - other_end_location_y ) + ( z_points_above_threshold[ idx_of_outer_point] - other_end_location_z ) * ( z_points_above_threshold[ idx_of_outer_point] - other_end_location_z ) );

	 std::cout << "Closest Distance To Outer Point = " << closest_distance_to_outer_point << " cm." << std::endl;
	 std::cout << "Closest Distance To Inner Point = " << closest_distance_to_inner_point << " cm." << std::endl;
	 std::cout << "Distance To Start of Track = " << distance_to_start_of_track << " cm." << std::endl;
	 std::cout << "Distance To End of Track = " << distance_to_end_of_track << " cm." << std::endl;

	 // Remove either or both of these points if they are greater than 0.6 cm from the closest point.
	 if ( closest_distance_to_outer_point > 0.6 && closest_distance_to_outer_point < 9999998. && ( distance_to_start_of_track < distance_to_end_of_track || distance_to_start_of_track > 50.0 ) ) {
	   // Append 'idx_of_outer_point' to 'idx_removed_points' and set 'outliers_found_last_iteration' to 'true'.
	   idx_removed_points.push_back( idx_of_outer_point );
	   outliers_found_last_iteration = true;
	   num_outliers_removed++;
	 }

	 if ( closest_distance_to_inner_point > 0.6 && closest_distance_to_inner_point < 9999998. && ( distance_to_start_of_track < distance_to_end_of_track || distance_to_start_of_track > 50.0 ) ) {
	   // Append 'idx_of_inner_point' to 'idx_removed_points' and set 'outliers_found_last_iteration' to 'true'.                                                                                    
           idx_removed_points.push_back( idx_of_inner_point );
           outliers_found_last_iteration = true;
	   num_outliers_removed++;
         }

       } // End of the requirement that the length of the track be greater than 3.0 cm to look for outliers.

       else {
	 outliers_found_last_iteration = false;
       }
       
     } while( outliers_found_last_iteration == true && ( num_calo_points_above_threshold - num_outliers_removed ) > 2 ); // End of the do....while statement.

     // Declare the TGraph object.
     TGraph2D * gr = new TGraph2D();

     std::cout << "The initial number of points above threshold = " << num_calo_points_above_threshold << "." << std::endl;

     for ( int initial_point_iter = 0; initial_point_iter < num_calo_points_above_threshold; initial_point_iter++ ) {

       std::cout << "Initial Point #" << initial_point_iter << ": x = " << x_points_above_threshold[ initial_point_iter ] << " cm y = " << y_points_above_threshold[ initial_point_iter ] << " cm z = " << z_points_above_threshold[ initial_point_iter ] << " cm." << std::endl;

     }
     
     int num_points_above_threshold_not_outliers  = 0;

     // Declare variables for the coordinates of the points on the track that are not outliers.
     double x_points_above_threshold_without_outliers[10000];
     double y_points_above_threshold_without_outliers[10000];
     double z_points_above_threshold_without_outliers[10000];

     // Fill these arrays with 0s.
     for ( size_t array_iter = 0; array_iter < 10000; array_iter++ ) {

       x_points_above_threshold_without_outliers[ array_iter ] = 0;
       y_points_above_threshold_without_outliers[ array_iter ] = 0;
       z_points_above_threshold_without_outliers[ array_iter ] = 0;
       
     }

     // Print out the event information.
     std::cout << "Run = " << run << " Subrun = " << subrun << " Event = " << event << "." << std::endl;
     
     // Loop through and remove all of the outliers.
     for ( size_t point_iter = 0; point_iter < num_calo_points_above_threshold; point_iter++ ) {

       bool point_is_an_outlier = false;
       
       for ( size_t outlier_iter = 0; outlier_iter < idx_removed_points.size(); outlier_iter++ ) {

	 if ( point_iter == idx_removed_points.at( outlier_iter ) ) point_is_an_outlier = true;

       }

       if ( point_is_an_outlier == true ) continue;

       std::cout << "Point #" << point_iter << ": X = " << x_points_above_threshold[point_iter] << " cm Y = " << y_points_above_threshold[point_iter] << " cm Z = " << z_points_above_threshold[point_iter] << " cm." << std::endl;

       gr->SetPoint( point_iter, x_points_above_threshold[point_iter], y_points_above_threshold[point_iter], z_points_above_threshold[point_iter] );

       // Set the point in the 'points_above_threshold_without_outliers' array.
       x_points_above_threshold_without_outliers[ num_points_above_threshold_not_outliers ] = x_points_above_threshold[point_iter];
       y_points_above_threshold_without_outliers[ num_points_above_threshold_not_outliers ] = y_points_above_threshold[point_iter];
       z_points_above_threshold_without_outliers[ num_points_above_threshold_not_outliers ] = z_points_above_threshold[point_iter];

       num_points_above_threshold_not_outliers++;

     }

     std::cout << "The number of points above threshold that are not outliers = " << num_points_above_threshold_not_outliers << "." << std::endl;

     double outer_x_coordinate          = 0.;
     double outer_y_coordinate          = 0.;
     double outer_z_coordinate          = 0.;
     double inner_x_coordinate          = 0.;
     double inner_y_coordinate          = 0.;
     double inner_z_coordinate          = 0.;

     double outer_furthest_x_coordinate = 0.;
     double outer_furthest_y_coordinate = 0.;
     double outer_furthest_z_coordinate = 0.;
     double inner_furthest_x_coordinate = 0.;
     double inner_furthest_y_coordinate = 0.;
     double inner_furthest_z_coordinate = 0.;

     double distance_between_points     =  0.;
     double furthest_distance           = -1.;

     size_t idx_of_outer_furthest_point = -1;
     size_t idx_of_inner_furthest_point = -1;
     
     // Find the pair of points with the furthest distance from one another.                                                                                                                            
     for ( size_t outer_iter = 0; outer_iter < size_t( num_points_above_threshold_not_outliers ); outer_iter++ ) {

       outer_x_coordinate = x_points_above_threshold_without_outliers[ outer_iter ];
       outer_y_coordinate = y_points_above_threshold_without_outliers[ outer_iter ];
       outer_z_coordinate = z_points_above_threshold_without_outliers[ outer_iter ];

       bool outer_point_is_an_outlier = false;

       for ( size_t outer_outlier_iter = 0; outer_outlier_iter < idx_removed_points.size(); outer_outlier_iter++ ) {
	 
         if ( outer_iter == idx_removed_points.at( outer_outlier_iter ) ) outer_point_is_an_outlier = true;
	 
       }

       if ( outer_point_is_an_outlier == true ) continue;

       for ( size_t inner_iter = outer_iter + 1; inner_iter < size_t( num_points_above_threshold_not_outliers ); inner_iter++ ) {

	 bool inner_point_is_an_outlier = false;

	 for ( size_t inner_outlier_iter = 0; inner_outlier_iter < idx_removed_points.size(); inner_outlier_iter++ ) {
	
	   if ( inner_iter == idx_removed_points.at( inner_outlier_iter ) ) inner_point_is_an_outlier = true;
	
	 }

	 if ( inner_point_is_an_outlier == true ) continue;

	 inner_x_coordinate = x_points_above_threshold_without_outliers[ inner_iter ];
	 inner_y_coordinate = y_points_above_threshold_without_outliers[ inner_iter ];
	 inner_z_coordinate = z_points_above_threshold_without_outliers[ inner_iter ];

	 distance_between_points = TMath::Sqrt( ( outer_x_coordinate - inner_x_coordinate ) * ( outer_x_coordinate - inner_x_coordinate ) + ( outer_y_coordinate - inner_y_coordinate ) * ( outer_y_coordinate - inner_y_coordinate ) + ( outer_z_coordinate - inner_z_coordinate ) * ( outer_z_coordinate - inner_z_coordinate )  );

	 if ( distance_between_points > furthest_distance ) {

	   // Set the variables for the points that are the furthest distance away.
	   furthest_distance           = distance_between_points;
	   outer_furthest_x_coordinate = outer_x_coordinate;
	   outer_furthest_y_coordinate = outer_y_coordinate;
	   outer_furthest_z_coordinate = outer_z_coordinate;
	   inner_furthest_x_coordinate = inner_x_coordinate;
           inner_furthest_y_coordinate = inner_y_coordinate;
           inner_furthest_z_coordinate = inner_z_coordinate;

	   // Set the value of the iterators that determine where the points are located.
	   idx_of_outer_furthest_point = outer_iter;
	   idx_of_inner_furthest_point = inner_iter;

	 }

       }

     }

     distance_between_furthest_calo_points = furthest_distance;

     // Print out the distance between the furthest points and the coordinates that go into calculating this value.
     std::cout << "Number of calorimetry points at the start = " << num_calo_points_above_threshold << "." << std::endl;
     std::cout << "Number of calorimetry points not outliers = " << num_points_above_threshold_not_outliers << "." << std::endl;
     std::cout << "Distance between furthest points = " << distance_between_furthest_calo_points << "." << std::endl;
     std::cout << "x-coordinate #1 = " << outer_furthest_x_coordinate << " cm." << std::endl;
     std::cout << "y-coordinate #1 = " << outer_furthest_y_coordinate << " cm." << std::endl;
     std::cout << "z-coordinate #1 = " << outer_furthest_z_coordinate << " cm." << std::endl;
     std::cout << "x-coordinate #2 = " << inner_furthest_x_coordinate << " cm." << std::endl;
     std::cout << "y-coordinate #2 = " << inner_furthest_y_coordinate << " cm." << std::endl;
     std::cout << "z-coordinate #2 = " << inner_furthest_z_coordinate << " cm." << std::endl;
	   
     double muon_momentum   = TMath::Sqrt( ( reco_muon_KE * reco_muon_KE ) + ( 2.0 * 105.7 * reco_muon_KE ) );

     // Find the direction in which to orient the track by comparing the muon endpoints to the vertex location.
     bool vertex_at_start = true;

     double distance_from_muon_start_to_location_of_vertex     = 0.;
     double distance_from_muon_end_to_location_of_vertex       = 0.;

     distance_from_muon_start_to_location_of_vertex = TMath::Sqrt( ( muon_track_first_point_x - vertex_location_x ) * ( muon_track_first_point_x - vertex_location_x ) + ( muon_track_first_point_y - vertex_location_y ) * ( muon_track_first_point_y - vertex_location_y ) + ( muon_track_first_point_z - vertex_location_z ) * ( muon_track_first_point_z - vertex_location_z ) );
     distance_from_muon_end_to_location_of_vertex   = TMath::Sqrt( ( muon_track_last_point_x - vertex_location_x ) * ( muon_track_last_point_x - vertex_location_x ) + ( muon_track_last_point_y - vertex_location_y ) * ( muon_track_last_point_y - vertex_location_y )	+ ( muon_track_last_point_z - vertex_location_z ) * ( muon_track_last_point_z - vertex_location_z ) );

     if ( distance_from_muon_end_to_location_of_vertex < distance_from_muon_start_to_location_of_vertex ) vertex_at_start = false;

     double resultant_track_length_from_points = TMath::Sqrt( ( muon_track_last_point_x - muon_track_first_point_x ) * ( muon_track_last_point_x - muon_track_first_point_x ) + ( muon_track_last_point_y - muon_track_first_point_y ) * ( muon_track_last_point_y - muon_track_first_point_y ) + ( muon_track_last_point_z - muon_track_first_point_z ) * ( muon_track_last_point_z - muon_track_first_point_z ) );
     
     muon_x_momentum = ( muon_momentum ) * ( muon_track_last_point_x - muon_track_first_point_x ) / ( resultant_track_length_from_points );
     muon_y_momentum = ( muon_momentum ) * ( muon_track_last_point_y - muon_track_first_point_y ) / ( resultant_track_length_from_points );
     muon_z_momentum = ( muon_momentum ) * ( muon_track_last_point_z - muon_track_first_point_z ) / ( resultant_track_length_from_points );

     // Printing out muon coordinates.
     std::cout << "Vertex Location: x = " << vertex_location_x << " cm y = " << vertex_location_y << " cm z = " << vertex_location_z << " cm." << std::endl;
     std::cout << "Muon Length = " << muon_candidate_track_length << " cm."  << std::endl;
     std::cout << "Reconstructed Muon KE = " << reco_muon_KE << " MeV." << std::endl;
     std::cout << "Muon Track Starting Coordinates: x = " << muon_track_first_point_x << " cm y = " << muon_track_first_point_y << " cm z = " << muon_track_first_point_z << " cm." << std::endl;
     std::cout << "Muon Track Ending Coordinates: x = " << muon_track_last_point_x << " cm y = " << muon_track_last_point_y	<< " cm z = " << muon_track_last_point_z << " cm." << std::endl;
     std::cout << "Muon Track Length From Points = " << resultant_track_length_from_points << " cm." << std::endl;
     std::cout << "Muon Momentum = " << muon_momentum << " MeV." << std::endl;
     std::cout << "Muon Momentum Components: x = " << muon_x_momentum << " y = " << muon_y_momentum << " z = " << muon_z_momentum << "." << std::endl;
     
     if ( vertex_at_start == false ) {

       std::cout << "Momentum components are flipped." << std::endl;

       muon_x_momentum *= -1.;
       muon_y_momentum *= -1.;
       muon_z_momentum *= -1.;

     }

     // fit the graph now 
     TVirtualFitter *min = TVirtualFitter::Fitter(0,4);
     min->SetObjectFit(gr);
     min->SetFCN(SumDistance2);
     
     Double_t arglist[10];
     arglist[0] = 0;
     min->ExecuteCommand("SET PRINT",arglist,1);
  
     double pStart[4] = {1,1,1,1};
     min->SetParameter(0,"x0",pStart[0],0.001,0,0);
     min->SetParameter(1,"Ax",pStart[1],0.001,0,0);
     min->SetParameter(2,"y0",pStart[2],0.001,0,0);
     min->SetParameter(3,"Ay",pStart[3],0.001,0,0);
    
     arglist[0] = 1000; // number of function calls 
     arglist[1] = 0.001; // tolerance 
     min->ExecuteCommand("MIGRAD",arglist,2);

     // Get the fit parameters from the minimizer.
     Double_t Ax_value;
     Double_t x0_value;
     Double_t Ay_value;
     Double_t y0_value;

     Double_t Ax_err;
     Double_t x0_err;
     Double_t Ay_err;
     Double_t y0_err;

     x0_value        = min->GetParameter(0);
     Ax_value        = min->GetParameter(1);
     y0_value        = min->GetParameter(2); 
     Ay_value        = min->GetParameter(3);

     x0_err          = min->GetParError(0);
     Ax_err          = min->GetParError(1);
     y0_err          = min->GetParError(2);
     Ay_err          = min->GetParError(3);

     // Set the values for the tree equal to these values.
     x_slope         = Ax_value;
     x_intercept     = x0_value;
     y_slope         = Ay_value;
     y_intercept     = y0_value;

     x_slope_err     = Ax_err;
     x_intercept_err = x0_err;
     y_slope_err     = Ay_err;
     y_intercept_err = y0_err;

     if ( num_calo_points_above_threshold > 0 ) { 

       double proton_x_starting_point = x_points_above_threshold[ idx_of_outer_furthest_point ];
       double proton_z_starting_point = z_points_above_threshold[ idx_of_outer_furthest_point ];

       double proton_x_ending_point   = x_points_above_threshold[ idx_of_inner_furthest_point ];
       double proton_z_ending_point   = z_points_above_threshold[ idx_of_inner_furthest_point ];

       double starting_t              = 0.;
       double ending_t                = 0.;

       if ( Ax_value > 0 ) {

	 starting_t = ( proton_x_starting_point - x0_value ) / ( Ax_value );
	 ending_t   = ( proton_x_ending_point - x0_value ) / ( Ax_value );

       }

       // Calculate the values of y from this.
       double proton_y_starting_point       = ( Ay_value * starting_t ) + y0_value;
       double proton_y_ending_point         = ( Ay_value * ending_t ) + y0_value;

       // Calculate the distance between the start and the end of the proton track.
       distance_using_track_projection      = TMath::Sqrt( ( proton_x_starting_point - proton_x_ending_point ) * ( proton_x_starting_point - proton_x_ending_point ) + ( proton_y_starting_point - proton_y_ending_point ) * ( proton_y_starting_point - proton_y_ending_point ) + ( proton_z_starting_point - proton_z_ending_point ) * ( proton_z_starting_point - proton_z_ending_point ) );

       // Calculae the projection of the proton along each axis.
       proton_x_projection                  = fabs( proton_x_ending_point - proton_x_starting_point );
       proton_y_projection                  = fabs( proton_y_ending_point - proton_y_starting_point );
       proton_z_projection                  = fabs( proton_z_ending_point - proton_z_starting_point );

       // Multiply the projections by each of the proton directional signs.
       proton_x_projection                 *= proton_x_direction;
       proton_y_projection                 *= proton_y_direction;
       proton_z_projection                 *= proton_z_direction;       

     } // End of the requirement that there are nonzero calorimetry points above threshold.
	 
     // Fill the tree.
     ftree->Fill();
       
   } // End of the loop over the number of events.
     
   std::cout << "The number of entries passing the additional cuts = " << num_entries_passing << "." << std::endl;
   
   file->Write();
   file->Close();

} // End of the function inside the script. 
   
int main() { 
  Updates_To_Proton_Fitting_Technique_With_All_MCWeights();
}

