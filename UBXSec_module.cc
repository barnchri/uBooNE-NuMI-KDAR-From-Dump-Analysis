////////////////////////////////////////////////////////////////////////
// Class:       UBXSec
// Plugin Type: producer
// File:        UBXSec_module.cc
//
// Generated at Thursday December 20th  12:00:00 by Christopher Barnes using cetskelgen
// from cetlib version v1_21_00.
////////////////////////////////////////////////////////////////////////

/**
 * \class UBXSec
 *
 * \ingroup UBXSec
 *
 * \brief Art producer module
 *
 *
 * \author Christopher Barnes 
 *
 * \version producer
 *
 * \date 2020/07/28
 *
 * Contact: barnchri@umich.edu
 *
 * Created on: December 20th, 2018 at 12:00:00
 *
 */

// Art include
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "canvas/Persistency/Common/FindManyP.h"

// Data products include
#include "lardataobj/MCBase/MCTrack.h"
#include "lardataobj/MCBase/MCShower.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/MCSFitResult.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"
#include "ubobj/UBXSec/FlashMatch.h"
#include "ubobj/UBXSec/MCGhost.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "larcoreobj/SummaryData/POTSummary.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "ubobj/UBXSec/TPCObject.h"
#include "lardataobj/RawData/TriggerData.h"
#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardataalg/DetectorInfo/DetectorClocks.h"
#include "larsim/EventWeight/Base/MCEventWeight.h"

#include "ubobj/UBXSec/UBXSecEvent.h"
#include "ubobj/UBXSec/SelectionResult.h"

#include "ubana/ParticleID/Algorithms/uB_PlaneIDBitsetHelperFunctions.h"

#include "ubevt/Utilities/SignalShapingServiceMicroBooNE.h"
#include "larevt/CalibrationDBI/Interface/DetPedestalService.h"
#include "larevt/CalibrationDBI/Interface/DetPedestalProvider.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusService.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusProvider.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"


// LArSoft include
#include "ubreco/UBFlashFinder/PECalib.h"
#include "larsim/MCCheater/BackTracker.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "larcore/Geometry/Geometry.h"
#include "ubobj/Trigger/ubdaqSoftwareTriggerData.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/MCBase/MCDataHolder.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

// Algorithms include
#include "ubana/UBXSec/Algorithms/UBXSecHelper.h"
#include "ubana/UBXSec/Algorithms/VertexCheck.h"
#include "ubana/UBXSec/Algorithms/FindDeadRegions.h"
#include "ubana/UBXSec/Algorithms/MuonCandidateFinder.h"
#include "ubana/UBXSec/Algorithms/FiducialVolume.h"
#include "ubana/UBXSec/Algorithms/NuMuCCEventSelection.h"
#include "ubana/UBXSec/Algorithms/TrackQuality.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/RecoBase/Slice.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "ubobj/CRT/CRTHit.hh"
#include "ubobj/RawData/DAQHeaderTimeUBooNE.h"

// Include the modules for the space charge correction.
#include "ubevt/SpaceCharge/SpaceChargeMicroBooNE.h"
#include "ubevt/SpaceChargeServices/SpaceChargeServiceMicroBooNE.h"

// Include the file that we need for the momentum calculation.                                                                                                                                          
#include "larreco/RecoAlg/TrackMomentumCalculator.h"

// event weight include file.                                                                                                                                                                            
#include "larsim/EventWeight/Base/MCEventWeight.h"

// Root include
#include "TString.h"
#include "TTree.h"
#include "TH2F.h"
#include "TH1D.h"
#include "TH2D.h"
#include <TDatabasePDG.h>
#include <TParticlePDG.h>

#include <fstream>
#include <string>
#include <iterator>

// Geometry package.
#include "larcore/Geometry/Geometry.h"

namespace ubxsec {
  struct Hit3D_t {
    double x;
    double y;
    double z;
    double q;
  };
}


class UBXSec;


class UBXSec : public art::EDProducer {
public:
  explicit UBXSec(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  UBXSec(UBXSec const &) = delete;
  UBXSec(UBXSec &&) = delete;
  UBXSec & operator = (UBXSec const &) = delete;
  UBXSec & operator = (UBXSec &&) = delete;

  // Required functions.
  void produce(art::Event & e) override;
  void endSubRun(art::SubRun &sr) override;

private:

  /// Prints MC particles from GENIE on the screen
  void PrintMC(std::vector<art::Ptr<simb::MCTruth>> mclist);

  /// Calculates flash position
  void GetFlashLocation(std::vector<double>, double&, double&, double&, double&);
  float GetTrackShowerScore(art::Ptr<recob::PFParticle> pfParticle, const lar_pandora::PFParticlesToMetadata pfParticleToMetadata);
  FindDeadRegions deadRegionsFinder;
void GetTaggedPFP(art::Event const & e, std::string cosmictag_producer, double score_cut, lar_pandora::PFParticleVector & pfpTaggedOut,std::vector<int> & tagid_v);
int PFPInCommon(lar_pandora::PFParticleVector first, lar_pandora::PFParticleVector second);
  //ubxsec::McPfpMatch mcpfpMatcher;
  ::ubana::FiducialVolume _fiducial_volume;
  ::ubana::MuonCandidateFinder _muon_finder;
  ::ubana::NuMuCCEventSelection _event_selection;
  ::pmtana::PECalib _pecalib;
  ::trkf::TrackMomentumCalculator  _trk_mom_calculator{0.1};

  // Database to understand particle pdg
  const TDatabasePDG* _database_pdg = TDatabasePDG::Instance();

  // Services
  ::detinfo::DetectorProperties const* _detector_properties;
  ::detinfo::DetectorClocks const* _detector_clocks;
  spacecharge::SpaceCharge const* _SCE;

  // To be set via fcl parameters
  std::string _hitfinderLabel;
  std::string _pfp_producer;
  std::string _geantModuleLabel;
  std::string _spacepointLabel;
  std::string _cosmic_tag_producer;
  std::string _neutrino_flash_match_producer;
  std::string _cosmic_flash_match_producer;
  std::string _opflash_producer_beam;

 std::string _cosmic_flash_tag_producer;
  std::string _cosmic_geo_tag_producer;
  std::string _cosmic_acpt_tag_producer;
  std::string _cosmic_stopmu_tag_producer;

  std::string _tpcobject_producer;
  std::string _potsum_producer;
  std::string _particle_id_producer;
  std::string _mc_ghost_producer;
  std::string _geocosmictag_producer;
  std::string _candidateconsistency_producer;
  std::string _mcsfitresult_mu_producer;
  std::string _calorimetry_producer;
  std::string _eventweight_producer;
  std::string _genie_eventweight_pm1_producer;
  std::string _genie_eventweight_multisim_producer;
  std::string _flux_eventweight_multisim_producer;
  bool _debug = false;                   ///< Debug mode
  bool _debug_cr = false;                   ///< Debug mode

  int _minimumHitRequirement;           ///< Minimum number of hits in at least a plane for a track
  double _minimumDistDeadReg;           ///< Minimum distance the track end points can have to a dead region
  bool _use_genie_info;                 ///< Turn this off if looking at cosmic only files
  double _beam_spill_start;             ///< Start time of the beam spill (us)
  double _beam_spill_end;               ///< Start time of the beam spill (us)
  double _total_pe_cut;                 ///< PE cut to be applied to beam flash
  double _geo_cosmic_score_cut;         ///< Cut on the score of the pandoraNu geo cosmic tagger
  double _tolerance_track_multiplicity; ///< Tolerance to consider a track coming from the nu reco vertex
  double _min_track_len;                ///< Min track length for momentum calculation
  bool _make_ophit_csv;                 ///< If true makea a csv file with ophit info
  bool _make_pida_csv;                  ///< If true makea a csv file with pida/tracklength info

  bool _do_opdet_swap;                  ///< If true swaps reconstructed OpDets according to _opdet_swap_map
  std::vector<int> _opdet_swap_map;     ///< The OpDet swap map for reco flashes
  double _cosmic_flash_tag_score_cut; ///< Score cut used in the analysis to consider the PFP as cosmic (applied to flash tagger)
  double _cosmic_geo_tag_score_cut;   ///< Score cut used in the analysis to consider the PFP as cosmic (applied to geo tagger)
  double _cosmic_acpt_tag_score_cut;  ///< Score cut used in the analysis to consider the PFP as cosmic (applied to acpt tagger)
  double _cosmic_stopmu_tag_score_cut;  ///< Score cut used in the analysis to consider the PFP as cosmic (applied to stopmu tagger)
  // Constants
  const simb::Origin_t NEUTRINO_ORIGIN = simb::kBeamNeutrino;

  // To be filled within module
  bool _is_data, _is_mc;
  //double _candidate_flash_time;
  double _drift_velocity;

  // Outputs trees
  TTree* _tree1;
  UBXSecEvent *ubxsec_event = new UBXSecEvent();

  TH2F * _deadRegion2P;
  TH2F * _deadRegion3P;

  // PIDA related variables
  TH1D * _h_pida_proton,     * _h_pida_muon,     * _h_pida_pion,     * _h_pida_kaon;
  TH2D * _h_pida_len_proton, * _h_pida_len_muon, * _h_pida_len_pion, * _h_pida_len_kaon;

  // Momentum Related Variables
  TH2D* _h_mom_true_mcs; ///< 2D histogram of true muon momentum VS reconstructed (using MCS)
  TH2D* _h_mom_true_mcs_contained; ///< 2D histogram of true muon momentum VS reconstructed (using MCS) (contained tracks)
  TH2D* _h_mom_true_mcs_uncontained; ///< 2D histogram of true muon momentum VS reconstructed (using MCS) (uncontained tracks)
  TH2D* _h_mom_true_range_contained; ///< 2D histogram of true muon momentum VS reconstructed (using Length) (contained tracks)
  TH2D* _h_mom_range_mcs_contained;  ///< 2D histogram of reconstructed (using MCS) muon momentum VS reconstructed (using Length) (contained tracks)
  TH1D* _h_mcs_cosmic_track_direction; ///< Track direction from cosmic origin TPCObjects as given by mcs (0: downward, 1: upward)
  TH2D* _h_mcs_cosmic_track_direction_deltall; /// Track direction from cosmic origin TPCObjects as given by mcs (0: downward, 1: upward) VS delta LL from MCS fit
  TH2D* _h_mcs_cosmic_track_direction_ratioll; /// Track direction from cosmic origin TPCObjects as given by mcs (0: downward, 1: upward) VS ratio LL from MCS fit
  TTree *_mom_tree_contained, *_mom_tree_uncontained;
  int _run, _subrun, _event;
  double _mom_true_contained;
  double _mom_mcs_contained;
  double _mom_range_contained;
  double _mom_true_uncontained;
  double _mom_mcs_uncontained;
  TTree *_mcs_cosmic_track_direction_tree;
  double _mcs_cosmic_track_direction;
  double _mcs_cosmic_track_downll;
  double _mcs_cosmic_track_upll;
  TTree *_mom_cosmic_tree;
  double _mom_cosmic_true;
  double _mom_cosmic_mcs;
  double _mom_cosmic_mcs_downforced;
  double _mom_cosmic_range;
  bool _mom_cosmic_down;

  // Subrun Tree.
  TTree* _sr_tree;
  int _sr_run, _sr_subrun;
  double _sr_begintime, _sr_endtime;
  double _sr_pot;

  // Passing Events Tree.
  TTree* _passing_events_tree;
  int    kdar_from_dump_event;
  int    num_of_neutrinos;
  double neutrino_energy;
  double truth_theta_angle_for_weighting;
  double truth_reco_vtx_distance;
  double parent_fvz;
  double nu_vtx_x_truth;
  double nu_vtx_y_truth;
  double nu_vtx_z_truth;
  double u_plane_vertex_location_x;
  double u_plane_vertex_location_y;
  double u_plane_vertex_location_z;
  double v_plane_vertex_location_x;
  double v_plane_vertex_location_y;
  double v_plane_vertex_location_z;
  double y_plane_vertex_location_x;
  double y_plane_vertex_location_y;
  double y_plane_vertex_location_z;
  double ubxsec_muon_phi;
  double ubxsec_muon_cos_theta;
  int    _nflashes_in_beamgate;
  int    _nflashes_in_beamspill;
  int    _nflashes_in_beamspill_window_passing_filter_PE_cut;
  double flash_time;
  double crt_flash_time; // Taking advantage of the different flash PE cut used by the CRT.
  double flash_PEs;
  double flash_z;
  double flash_y;
  double flash_z_width;
  double flash_y_width;
  int    nslices_in_event;
  int    number_of_tracks_in_TPCObject;
  int    num_of_tracks_originating_from_vertex;
  double muon_candidate_length;
  double muon_candidate_kinetic_energy;
  int    num_mctrack_points;
  double truth_muon_length;
  double truth_muon_kinetic_energy;
  double sum_of_TPCObject_track_lengths;
  double total_length_of_tracks_originating_from_vertex;
  int    num_calo_points;
  double TPCObject_kinetic_energy_using_muon_and_pion_assumption;
  double truncated_dQdx_of_TPCObject_muon_candidate_u_plane;
  double truncated_dQdx_of_TPCObject_muon_candidate_v_plane;
  double truncated_dQdx_of_TPCObject_muon_candidate_y_plane;
  double median_dQdx_of_TPCObject_muon_candidate_u_plane;
  double median_dQdx_of_TPCObject_muon_candidate_v_plane;
  double median_dQdx_of_TPCObject_muon_candidate_y_plane;
  double total_sum_of_slice_associated_hits_ADCs;
  double u_plane_sum_of_slice_associated_hits_ADCs;
  double v_plane_sum_of_slice_associated_hits_ADCs;
  double y_plane_sum_of_slice_associated_hits_ADCs;
  int    total_number_of_slice_associated_hits;
  int    u_plane_number_of_slice_associated_hits;
  int    v_plane_number_of_slice_associated_hits;
  int    y_plane_number_of_slice_associated_hits;
  double total_sum_of_slice_hits_in_vicinity_ADCs;
  double u_plane_sum_of_slice_hits_in_vicinity_ADCs;
  double v_plane_sum_of_slice_hits_in_vicinity_ADCs;
  double y_plane_sum_of_slice_hits_in_vicinity_ADCs;
  int    total_number_of_slice_hits_in_vicinity;
  int    u_plane_number_of_slice_hits_in_vicinity;
  int    v_plane_number_of_slice_hits_in_vicinity;
  int    y_plane_number_of_slice_hits_in_vicinity;
  int    num_of_top_pandora_crossings;
  int    num_of_bottom_pandora_crossings;
  int    num_of_front_pandora_crossings;
  int    num_of_back_pandora_crossings;
  int    NC_channel;
  double spline_fix_mcweight;
  double central_value_mcweight;
  double rootino_fix_mcweight;
  double other_universe_mcweights[100];
  double axialff_mcweights[2];
  double rpaccqe_mcweights[2];
  double xsecshape_mcweights[2];
  int    num_muminus_tracks;
  int    num_muplus_tracks;
  int    num_piplus_tracks;
  int    num_piminus_tracks;
  int    num_pi0_tracks;
  int    num_proton_tracks;
  int    num_electron_showers;
  int    num_positron_showers;
  int    num_photon_showers;

  // Truth variables used for weighting the Monte Carlo events.
  double truth_px_for_theta_for_weighting;
  double truth_py_for_theta_for_weighting;
  double truth_pz_for_theta_for_weighting;

  std::ofstream _csvfile, _csvfile2;

  // Tree For Quantities Necessary for Proton Reconstruction.
  TTree* proton_hit_identification_tree;
  int    num_points_above_threshold;
  int    num_u_plane_points_above_threshold_in_vertexing_loop;
  int    num_v_plane_points_above_threshold_in_vertexing_loop;
  int    num_y_plane_points_above_threshold_in_vertexing_loop;
  double x_points_above_threshold[10000];
  double y_points_above_threshold[10000];
  double z_points_above_threshold[10000];
  double reconstructed_muon_KE;
  double total_proton_KE;
  double leading_proton_KE;
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
  double proton_kinetic_energy_from_additional_tracks;
  double x_momentum_component_sum_additional_tracks;
  double y_momentum_component_sum_additional_tracks;
  double z_momentum_component_sum_additional_tracks;
  int    num_additional_tracks_originating_from_vertex;
  int    proton_x_direction;
  int    proton_y_direction;
  int    proton_z_direction;

  // Tree To Use For Comparing Events Across Systematic Samples.
  TTree* num_events_looping_over_tree;
  int    num_event_we_are_on;

  // Declare the variables necessary for the CRT Veto cut.                                                                                                                                              
  int    within_resolution;
  int    _nflashes_in_beamgate_passing_beamspill_and_PE_cuts;

  // Vertex info with SCE.
  double vertex_location_x_with_SCE;
  double vertex_location_y_with_SCE;
  double vertex_location_z_with_SCE;
  
  // Plane-level vertex info.
  double u_plane_truth_reco_vertex_distance;
  double v_plane_truth_reco_vertex_distance;
  double y_plane_truth_reco_vertex_distance;

  int    has_u_plane_info;
  int    has_v_plane_info;
  int    has_y_plane_info;
  int    has_u_plane_points_to_use_for_vertex;
  int    has_v_plane_points_to_use_for_vertex;
  int    has_y_plane_points_to_use_for_vertex;

  double u_plane_vertex_location_x_with_SCE;
  double u_plane_vertex_location_y_with_SCE;
  double u_plane_vertex_location_z_with_SCE;
  double v_plane_vertex_location_x_with_SCE;
  double v_plane_vertex_location_y_with_SCE;
  double v_plane_vertex_location_z_with_SCE;
  double y_plane_vertex_location_x_with_SCE;
  double y_plane_vertex_location_y_with_SCE;
  double y_plane_vertex_location_z_with_SCE;

  // Thresholds for finding calorimetry hits above threshold.
  double fUPlaneThreshold;
  double fVPlaneThreshold;
  double fYPlaneThreshold;

  // Declare variables for finding the direction of the muon candidate.
  bool   u_plane_has_more_than_20_points;
  bool   v_plane_has_more_than_20_points;
  bool   y_plane_has_more_than_20_points;
  bool   no_calo_points_on_u_plane;
  bool   no_calo_points_on_v_plane;
  bool   no_calo_points_on_y_plane;
  bool   u_plane_has_track_forward;
  bool   v_plane_has_track_forward;
  bool   y_plane_has_track_forward;

  // Calorimetry vectors.
  std::vector< float > residualrange_vector;
  std::vector< float > dedx_vector;

  // All of the variables for the hits that are off the track but above threshold.
  double wire_space_max_distance;
  double tick_space_max_distance;
  double maximum_proton_length;
  double maximum_proton_kinetic_energy;

  int    num_hits_above_threshold_u_plane;
  int    num_hits_above_threshold_v_plane;
  int    num_hits_above_threshold_y_plane;

  std::vector< double > u_plane_ticks_above_threshold;
  std::vector< double > u_plane_wires_above_threshold;
  std::vector< double > u_plane_charge_integrals_above_threshold;

  std::vector< double > v_plane_ticks_above_threshold;
  std::vector< double > v_plane_wires_above_threshold;
  std::vector< double > v_plane_charge_integrals_above_threshold;

  std::vector< double > y_plane_ticks_above_threshold;
  std::vector< double > y_plane_wires_above_threshold;
  std::vector< double > y_plane_charge_integrals_above_threshold;

  // Variables necessary for calculating start/end points of the track on each plane.
  double u_plane_muon_track_start_x;
  double u_plane_muon_track_start_y;
  double u_plane_muon_track_start_z;
  double u_plane_muon_track_end_x;
  double u_plane_muon_track_end_y;
  double u_plane_muon_track_end_z;

  double v_plane_muon_track_start_x;
  double v_plane_muon_track_start_y;
  double v_plane_muon_track_start_z;
  double v_plane_muon_track_end_x;
  double v_plane_muon_track_end_y;
  double v_plane_muon_track_end_z;

  double y_plane_muon_track_start_x;
  double y_plane_muon_track_start_y;
  double y_plane_muon_track_start_z;
  double y_plane_muon_track_end_x;
  double y_plane_muon_track_end_y;
  double y_plane_muon_track_end_z;
  
  double u_plane_wire_at_other_end;
  double u_plane_other_end_tick;
  double v_plane_wire_at_other_end;
  double v_plane_other_end_tick;
  double y_plane_wire_at_other_end;
  double y_plane_other_end_tick;

  // Variable used in more extensive two-track length cut.
  std::vector< double > lengths_of_tracks_in_TPC_object;

  // Variable used to find if other tracks fail the fiducial volume requirement.
  std::vector< size_t > idx_of_additional_tracks_originating_from_vertex;

  // Variable used in CRT distance cut.
  double closest_distance_to_CRT_tagged_track;

  // Additional variables used in proton orientation algorithm.
  bool   u_plane_has_nonzero_nearby_points_and_in_the_top_two_planes;
  bool   v_plane_has_nonzero_nearby_points_and_in_the_top_two_planes;
  bool   y_plane_has_nonzero_nearby_points_and_in_the_top_two_planes;

  bool   fails_CRT_distance_cut;
  bool   fails_new_two_track_selection;
  bool   fails_3cm_containment_cut;
  bool   fails_fiducial_volume_cut;

  // Count the total number of events, the number of events passing, and the number of events failing.                                                                                                 
  int    total_num_events;
  int    total_num_events_failing;
  int    total_num_events_passing;
  int    total_num_duplicate_events;

  // The number of events failing each of the different cuts.
  int    total_num_events_failed_crt_veto;
  int    total_num_events_failed_beam_disc_flashes;
  int    total_num_events_failed_beam_spill_flash;
  int    total_num_events_failed_has_slices;
  int    total_num_events_failed_track_length;
  int    total_num_events_failed_has_slice_tagged_as_neutrino;
  int    total_num_events_failed_fiducial_volume;
  int    total_num_events_failed_ntrack;
  int    total_num_events_failed_residuals_std_up;
  int    total_num_events_failed_perc_used_hits_in_cluster;
  int    total_num_events_failed_CRT_distance_cut;
  int    total_num_events_failing_new_two_track_selection;
  int    total_num_events_failing_3cm_containment_cut;
  int    total_num_events_failing_5cm_containment_cut;
  int    total_num_events_failing_length_of_tracks_from_vertex_requirement;
  int    total_num_events_failing_total_length_of_tracks_in_TPC_Object;
  int    total_num_events_failing_less_than_4_tracks_in_tpcobject_cut;
  int    total_num_events_failing_more_than_2_tracks_originating_from_vertex_cut;
  int    total_num_events_failing_CRT_distance_cut;

  // Count the number of events passing with a specific topology.
  int    total_num_events_passing_with_not_all_tracks_contained_in_fiducial_volume;
  
  // Count the number of signal events passing.
  int num_kdar_events_in_this_sample;
  int num_kdar_events_passing_in_this_sample;

  // Make sure there is no duplicity among events.
  bool   is_duplicate_event;
  std::vector< int > duplicate_runs;
  std::vector< int > duplicate_subruns;
  std::vector< int > duplicate_events;

  // Variables for track direction reconstruction/selection.
  int  muon_track_is_correctly_oriented;
  int  event_fails_new_two_track_selection;

  // Count the number of events you've looped over.
  int  event_counter;

};


UBXSec::UBXSec(fhicl::ParameterSet const & p) {

  ::art::ServiceHandle<geo::Geometry> geo;

  _pfp_producer                        = p.get<std::string>("PFParticleProducer");
  _hitfinderLabel                      = p.get<std::string>("HitProducer");
  _geantModuleLabel                    = p.get<std::string>("GeantModule");
  _spacepointLabel                     = p.get<std::string>("SpacePointProducer");
  _neutrino_flash_match_producer       = p.get<std::string>("NeutrinoFlashMatchProducer");
  _cosmic_flash_match_producer         = p.get<std::string>("CosmicFlashMatchProducer");
  _opflash_producer_beam               = p.get<std::string>("OpFlashBeamProducer");
  _tpcobject_producer                  = p.get<std::string>("TPCObjectProducer");
  _potsum_producer                     = p.get<std::string>("POTSummaryProducer");
  _particle_id_producer                = p.get<std::string>("ParticleIDProducer");
  _mc_ghost_producer                   = p.get<std::string>("MCGhostProducer");
  _geocosmictag_producer               = p.get<std::string>("GeoCosmicTaggerProducer");
  _candidateconsistency_producer       = p.get<std::string>("CandidateConsistencyProducer");
  _mcsfitresult_mu_producer            = p.get<std::string>("MCSFitResultMuProducer");
  _calorimetry_producer                = p.get<std::string>("CalorimetryProducer");
  _eventweight_producer                = p.get<std::string>("EventWeightProducer");
  _genie_eventweight_pm1_producer      = p.get<std::string>("GenieEventWeightPMOneProducer");
  _genie_eventweight_multisim_producer = p.get<std::string>("GenieEventWeightMultisimProducer");
  _flux_eventweight_multisim_producer  = p.get<std::string>("FluxEventWeightMultisimProducer");
  _cosmic_flash_tag_producer           = p.get<std::string>("CosmicFlashTagProducer");
  _cosmic_geo_tag_producer             = p.get<std::string>("CosmicGeoTagProducer");
  _cosmic_acpt_tag_producer            = p.get<std::string>("CosmicACPTTagProducer");
  _cosmic_stopmu_tag_producer          = p.get<std::string>("CosmicStopMuTagProducer");
  _mc_ghost_producer                   = p.get<std::string>("MCGhostProducer");
  _cosmic_flash_tag_score_cut          = p.get<double>("CosmicFlashTagScoreCut",0.99);
  _cosmic_geo_tag_score_cut            = p.get<double>("CosmicGeoTagScoreCut",0.6);
  _cosmic_acpt_tag_score_cut           = p.get<double>("CosmicACPTTagScoreCut",0.99);
  _cosmic_stopmu_tag_score_cut         = p.get<double>("CosmicStopMuTagScoreCut",0.99);

  _use_genie_info                      = p.get<bool>("UseGENIEInfo", false);
  _minimumHitRequirement               = p.get<int>("MinimumHitRequirement", 3);
  _minimumDistDeadReg                  = p.get<double>("MinimumDistanceToDeadRegion", 5.);

  _beam_spill_start                    = p.get<double>("BeamSpillStart", 3.2);
  _beam_spill_end                      = p.get<double>("BeamSpillEnd",   4.8);
  _total_pe_cut                        = p.get<double>("TotalPECut",     50);

  _do_opdet_swap                       = p.get<bool>("DoOpDetSwap", false);
  _opdet_swap_map                      = p.get<std::vector<int> >("OpDetSwapMap");

  _geo_cosmic_score_cut                = p.get<double>("GeoCosmicScoreCut", 0.6);
  _tolerance_track_multiplicity        = p.get<double>("ToleranceTrackMultiplicity", 5.);

  _make_ophit_csv                      = p.get<bool>("MakeOpHitCSV", false);
  _make_pida_csv                       = p.get<bool>("MakePIDACSV", false);

  _pecalib.Configure(p.get<fhicl::ParameterSet>("PECalib"));

  _fiducial_volume.Configure(p.get<fhicl::ParameterSet>("FiducialVolumeSettings"),
                             geo->DetHalfHeight(),
                             2.*geo->DetHalfWidth(),
                             geo->DetLength());

  _fiducial_volume.PrintConfig();

  _muon_finder.Configure(p.get<fhicl::ParameterSet>("MuonCandidateFinderSettings"));

  _muon_finder.PrintConfig();

  _event_selection.Configure(p.get<fhicl::ParameterSet>("NuMuCCSelectionSettings"));

  _event_selection.PrintConfig();

  _detector_properties = lar::providerFrom<detinfo::DetectorPropertiesService>();
  _detector_clocks = lar::providerFrom<detinfo::DetectorClocksService>();
  _SCE = lar::providerFrom<spacecharge::SpaceChargeService>();

  std::cout << "E Field: "            << _detector_properties->Efield() << std::endl;
  std::cout << "Temperature: "        << _detector_properties->Temperature() << std::endl;
  std::cout << "Drift Velocity: "     << _detector_properties->DriftVelocity(_detector_properties->Efield(), _detector_properties->Temperature())<< std::endl;
  std::cout << "Sampling Rate: "      << _detector_properties->SamplingRate() << std::endl;
  std::cout << "Beam Spill Start = "  << _beam_spill_start << " us." << std::endl;
  std::cout << "Beam Spill Ends = "   << _beam_spill_end << " us." << std::endl;

  art::ServiceHandle<art::TFileService> fs;
  _tree1 = fs->make<TTree>("tree","");

  int bufsize    = 16000;
  int splitlevel = 99;
  _tree1->Branch("ubxsec_event_split", &ubxsec_event, bufsize, splitlevel);

  _deadRegion2P                         = fs->make<TH2F>("deadRegion2P","deadRegion2P", 10350,0.0,1035.0,2300,-115.0,115.0);
  _deadRegion3P                         = fs->make<TH2F>("deadRegion3P","deadRegion3P", 10350,0.0,1035.0,2300,-115.0,115.0);

  _h_pida_muon                          = fs->make<TH1D>("h_pida_muon", "Muon tracks;PIDa [MeV/cm^{1.42}];", 50, 0, 20);
  _h_pida_proton                        = fs->make<TH1D>("h_pida_proton", "Proton tracks;PIDa [MeV/cm^{1.42}];", 50, 0, 20);
  _h_pida_pion                          = fs->make<TH1D>("h_pida_pion", "Pion tracks;PIDa [MeV/cm^{1.42}];", 50, 0, 20);
  _h_pida_kaon                          = fs->make<TH1D>("h_pida_kaon", "Kaon tracks;PIDa [MeV/cm^{1.42}];", 50, 0, 20);

  _h_pida_len_muon                      = fs->make<TH2D>("h_pida_len_muon", "Muon tracks;PIDa [MeV/cm^{1.42}];Track length [cm];", 50, 0, 20, 100, 0, 700);
  _h_pida_len_proton                    = fs->make<TH2D>("h_pida_len_proton", "Proton tracks;PIDa [MeV/cm^{1.42}];Track length [cm];", 50, 0, 20, 100, 0, 700);
  _h_pida_len_pion                      = fs->make<TH2D>("h_pida_len_pion", "Pion tracks;PIDa [MeV/cm^{1.42}];Track length [cm];", 50, 0, 20, 100, 0, 700);
  _h_pida_len_kaon                      = fs->make<TH2D>("h_pida_len_kaon", "Kaon tracks;PIDa [MeV/cm^{1.42}];Track length [cm];", 50, 0, 20, 100, 0, 700);

  _h_mom_true_mcs                       = fs->make<TH2D>("h_mom_true_mcs", ";True Muon Momentum [GeV];Reconstructed (via MCS) Muon Momentum [GeV];", 80, 0, 2, 80, 0, 2);
  _h_mom_true_mcs_contained             = fs->make<TH2D>("h_mom_true_mcs_contained", "Contained;True Muon Momentum [GeV];Reconstructed (via MCS) Muon Momentum [GeV];", 80, 0, 2, 80, 0, 2);
  _h_mom_true_mcs_uncontained           = fs->make<TH2D>("h_mom_true_mcs_uncontained", "Uncontained;True Muon Momentum [GeV];Reconstructed (via MCS) Muon Momentum [GeV];", 80, 0, 2, 80, 0, 2);
  _h_mom_true_range_contained           = fs->make<TH2D>("h_mom_true_range_contained", "Contained;True Muon Momentum [GeV];Reconstructed (via Length) Muon Momentum [GeV];", 80, 0, 2, 80, 0, 2);
  _h_mom_range_mcs_contained            = fs->make<TH2D>("h_mom_range_mcs_contained", "Contained;Reconstructed (via Length) Muon Momentum [GeV];Reconstructed (via MCS) Muon Momentum [GeV];", 80, 0, 2, 80, 0, 2);

  _h_mcs_cosmic_track_direction         = fs->make<TH1D>("h_mcs_cosmic_track_direction", "0: down, 1: up;;", 2, 0, 2);
  _h_mcs_cosmic_track_direction_deltall = fs->make<TH2D>("h_mcs_cosmic_track_direction_deltall", ";0: down, 1: up;(FWD - BWD) LL;", 2, 0, 2, 500, -1e-6, 1e-6);
  _h_mcs_cosmic_track_direction_ratioll = fs->make<TH2D>("h_mcs_cosmic_track_direction_ratioll", ";0: down, 1: up;(FWD / BWD) LL;", 2, 0, 2, 500, 0, 2);

  _mom_tree_contained              = fs->make<TTree>("mom_tree_contained","");
  _mom_tree_contained->Branch("run",                                                             &_run,                                               "run/I");
  _mom_tree_contained->Branch("subrun",                                                          &_subrun,                                            "subrun/I");
  _mom_tree_contained->Branch("event",                                                           &_event,                                             "event/I");
  _mom_tree_contained->Branch("mom_true_contained",                                              &_mom_true_contained,                                "mom_true_contained/D");
  _mom_tree_contained->Branch("mom_mcs_contained",                                               &_mom_mcs_contained,                                 "mom_mcs_contained/D");
  _mom_tree_contained->Branch("mom_range_contained",                                             &_mom_range_contained,                               "mom_range_contained/D");

  _mom_tree_uncontained            = fs->make<TTree>("mom_tree_uncontained","");
  _mom_tree_uncontained->Branch("run",                                                           &_run,                                                "run/I");
  _mom_tree_uncontained->Branch("subrun",                                                        &_subrun,                                             "subrun/I");
  _mom_tree_uncontained->Branch("event",                                                         &_event,                                              "event/I");
  _mom_tree_uncontained->Branch("mom_true_uncontained",                                          &_mom_true_uncontained,                               "mom_true_uncontained/D");
  _mom_tree_uncontained->Branch("mom_mcs_uncontained",                                           &_mom_mcs_uncontained,                                "mom_mcs_uncontained/D");

  _mcs_cosmic_track_direction_tree = fs->make<TTree>("mcs_cosmic_track_direction_tree","");
  _mcs_cosmic_track_direction_tree->Branch("run",                                                &_run,                                                "run/I");
  _mcs_cosmic_track_direction_tree->Branch("subrun",                                             &_subrun,                                             "subrun/I");
  _mcs_cosmic_track_direction_tree->Branch("event",                                              &_event,                                              "event/I");
  _mcs_cosmic_track_direction_tree->Branch("mcs_cosmic_track_direction",                         &_mcs_cosmic_track_direction,                         "mcs_cosmic_track_direction/D");
  _mcs_cosmic_track_direction_tree->Branch("mcs_cosmic_track_downll",                            &_mcs_cosmic_track_downll,                            "mcs_cosmic_track_downll/D");
  _mcs_cosmic_track_direction_tree->Branch("mcs_cosmic_track_upll",                              &_mcs_cosmic_track_upll,                              "mcs_cosmic_track_upll/D");

  _mom_cosmic_tree                 = fs->make<TTree>("mom_cosmic_tree","");
  _mom_cosmic_tree->Branch("run",                                                                &_run,                                                "run/I");
  _mom_cosmic_tree->Branch("subrun",                                                             &_subrun,                                             "subrun/I");
  _mom_cosmic_tree->Branch("event",                                                              &_event,                                              "event/I");
  _mom_cosmic_tree->Branch("mom_cosmic_true",                                                    &_mom_cosmic_true,                                    "mom_cosmic_true/D");
  _mom_cosmic_tree->Branch("mom_cosmic_mcs",                                                     &_mom_cosmic_mcs,                                     "mom_cosmic_mcs/D");
  _mom_cosmic_tree->Branch("mom_cosmic_mcs_downforced",                                          &_mom_cosmic_mcs_downforced,                          "mom_cosmic_mcs_downforced/D");
  _mom_cosmic_tree->Branch("mom_cosmic_range",                                                   &_mom_cosmic_range,                                   "mom_cosmic_range/D");
  _mom_cosmic_tree->Branch("mom_cosmic_down",                                                    &_mom_cosmic_down,                                    "mom_cosmic_down/O");

  _sr_tree = fs->make<TTree>("pottree","");
  _sr_tree->Branch("run",                                                                        &_sr_run,                                             "run/I");
  _sr_tree->Branch("subrun",                                                                     &_sr_subrun,                                          "subrun/I");
  _sr_tree->Branch("begintime",                                                                  &_sr_begintime,                                       "begintime/D");
  _sr_tree->Branch("endtime",                                                                    &_sr_endtime,                                         "endtime/D");
  _sr_tree->Branch("pot",                                                                        &_sr_pot,                                             "pot/D");

  // Declare the branches for the passing events tree.
  _passing_events_tree = fs->make<TTree>("_passing_events_tree", "A tree that contains muon quantities for the events that pass the NuMuCCInclusive Filter");
  _passing_events_tree->Branch("run",                                                            &_run,                                                "run/I");
  _passing_events_tree->Branch("subrun",                                                         &_subrun,                                             "subrun/I");
  _passing_events_tree->Branch("event",                                                          &_event,                                              "event/I");
  _passing_events_tree->Branch("kdar_from_dump_event",                                           &kdar_from_dump_event,                                "kdar_from_dump_event/I");
  _passing_events_tree->Branch("num_of_neutrinos",                                               &num_of_neutrinos,                                    "num_of_neutrinos/I");
  _passing_events_tree->Branch("neutrino_energy",                                                &neutrino_energy,                                     "neutrino_energy/D");
  _passing_events_tree->Branch("truth_theta_angle_for_weighting",                                &truth_theta_angle_for_weighting,                     "truth_theta_angle_for_weighting/D");
  _passing_events_tree->Branch("truth_reco_vtx_distance",                                        &truth_reco_vtx_distance,                             "truth_reco_vtx_distance/D");
  _passing_events_tree->Branch("parent_fvz",                                                     &parent_fvz,                                          "parent_fvz/D");
  _passing_events_tree->Branch("nu_vtx_x_truth",                                                 &nu_vtx_x_truth,                                      "nu_vtx_x_truth/D");
  _passing_events_tree->Branch("nu_vtx_y_truth",                                                 &nu_vtx_y_truth,                                      "nu_vtx_y_truth/D");
  _passing_events_tree->Branch("nu_vtx_z_truth",                                                 &nu_vtx_z_truth,                                      "nu_vtx_z_truth/D");
  _passing_events_tree->Branch("vertex_location_x",                                              &vertex_location_x,                                   "vertex_location_x/D");
  _passing_events_tree->Branch("vertex_location_y",                                              &vertex_location_y,                                   "vertex_location_y/D");
  _passing_events_tree->Branch("vertex_location_z",                                              &vertex_location_z,                                   "vertex_location_z/D");
  _passing_events_tree->Branch("u_plane_vertex_location_x",                                      &u_plane_vertex_location_x,                           "u_plane_vertex_location_x/D");
  _passing_events_tree->Branch("u_plane_vertex_location_y",                                      &u_plane_vertex_location_y,                           "u_plane_vertex_location_y/D");
  _passing_events_tree->Branch("u_plane_vertex_location_z",                                      &u_plane_vertex_location_z,                           "u_plane_vertex_location_z/D");
  _passing_events_tree->Branch("v_plane_vertex_location_x",                                      &v_plane_vertex_location_x,                           "v_plane_vertex_location_x/D");
  _passing_events_tree->Branch("v_plane_vertex_location_y",                                      &v_plane_vertex_location_y,                           "v_plane_vertex_location_y/D");
  _passing_events_tree->Branch("v_plane_vertex_location_z",                                      &v_plane_vertex_location_z,                           "v_plane_vertex_location_z/D");
  _passing_events_tree->Branch("y_plane_vertex_location_x",                                      &y_plane_vertex_location_x,                           "y_plane_vertex_location_x/D");
  _passing_events_tree->Branch("y_plane_vertex_location_y",                                      &y_plane_vertex_location_y,                           "y_plane_vertex_location_y/D");
  _passing_events_tree->Branch("y_plane_vertex_location_z",                                      &y_plane_vertex_location_z,                           "y_plane_vertex_location_z/D");
  _passing_events_tree->Branch("ubxsec_muon_phi",                                                &ubxsec_muon_phi,                                     "ubxsec_muon_phi/D");
  _passing_events_tree->Branch("ubxsec_muon_cos_theta",                                          &ubxsec_muon_cos_theta,                               "ubxsec_muon_cos_theta/D");
  _passing_events_tree->Branch("num_beam_flashes",                                               &_nflashes_in_beamgate,                               "num_beam_flashes/I"); 
  _passing_events_tree->Branch("_nflashes_in_beamspill",                                         &_nflashes_in_beamspill,                              "_nflashes_in_beamspill/I");
  _passing_events_tree->Branch("_nflashes_in_beamspill_window_passing_filter_PE_cut",            &_nflashes_in_beamspill_window_passing_filter_PE_cut, "_nflashes_in_beamspill_window_passing_filter_PE_cut/I");
  _passing_events_tree->Branch("flash_time",                                                     &flash_time,                                           "flash_time/D");
  _passing_events_tree->Branch("flash_PEs",                                                      &flash_PEs,                                            "flash_PEs/D");
  _passing_events_tree->Branch("flash_z",                                                        &flash_z,                                              "flash_z/D");
  _passing_events_tree->Branch("flash_y",                                                        &flash_y,                                              "flash_y/D");
  _passing_events_tree->Branch("flash_z_width",                                                  &flash_z_width,                                        "flash_z_width/D");
  _passing_events_tree->Branch("flash_y_width",                                                  &flash_y_width,                                        "flash_y_width/D");
  _passing_events_tree->Branch("nslices_in_event",                                               &nslices_in_event,                                     "nslices_in_event/I");
  _passing_events_tree->Branch("number_of_tracks_in_TPCObject",                                  &number_of_tracks_in_TPCObject,                        "number_of_tracks_in_TPCObject/I");
  _passing_events_tree->Branch("num_of_tracks_originating_from_vertex",                          &num_of_tracks_originating_from_vertex,                "num_of_tracks_originating_from_vertex/I");
  _passing_events_tree->Branch("muon_candidate_length",                                          &muon_candidate_length,                                "muon_candidate_length/D");
  _passing_events_tree->Branch("muon_candidate_kinetic_energy",                                  &muon_candidate_kinetic_energy,                        "muon_candidate_kinetic_energy/D");
  _passing_events_tree->Branch("num_mctrack_points",                                             &num_mctrack_points,                                   "num_mctrack_points/I");
  _passing_events_tree->Branch("truth_muon_length",                                              &truth_muon_length,                                    "truth_muon_length/D");
  _passing_events_tree->Branch("truth_muon_kinetic_energy",                                      &truth_muon_kinetic_energy,                            "truth_muon_kinetic_energy/D");
  _passing_events_tree->Branch("sum_of_TPCObject_track_lengths",                                 &sum_of_TPCObject_track_lengths,                       "sum_of_TPCObject_track_lengths/D");
  _passing_events_tree->Branch("total_length_of_tracks_originating_from_vertex",                 &total_length_of_tracks_originating_from_vertex,       "total_length_of_tracks_originating_from_vertex/D");
  _passing_events_tree->Branch("num_calo_points",                                                &num_calo_points,                                      "num_calo_points/I");
  _passing_events_tree->Branch("truncated_dQdx_of_TPCObject_muon_candidate_u_plane",             &truncated_dQdx_of_TPCObject_muon_candidate_u_plane,   "truncated_dQdx_of_TPCObject_muon_candidate_u_plane/D");
  _passing_events_tree->Branch("truncated_dQdx_of_TPCObject_muon_candidate_v_plane",             &truncated_dQdx_of_TPCObject_muon_candidate_v_plane,   "truncated_dQdx_of_TPCObject_muon_candidate_v_plane/D");
  _passing_events_tree->Branch("truncated_dQdx_of_TPCObject_muon_candidate_y_plane",             &truncated_dQdx_of_TPCObject_muon_candidate_y_plane,   "truncated_dQdx_of_TPCObject_muon_candidate_y_plane/D");
  _passing_events_tree->Branch("median_dQdx_of_TPCObject_muon_candidate_u_plane",                &median_dQdx_of_TPCObject_muon_candidate_u_plane,      "median_dQdx_of_TPCObject_muon_candidate_u_plane/D");
  _passing_events_tree->Branch("median_dQdx_of_TPCObject_muon_candidate_v_plane",                &median_dQdx_of_TPCObject_muon_candidate_v_plane,      "median_dQdx_of_TPCObject_muon_candidate_v_plane/D");
  _passing_events_tree->Branch("median_dQdx_of_TPCObject_muon_candidate_y_plane",                &median_dQdx_of_TPCObject_muon_candidate_y_plane,      "median_dQdx_of_TPCObject_muon_candidate_y_plane/D");
  _passing_events_tree->Branch("total_sum_of_slice_associated_hits_ADCs",                        &total_sum_of_slice_associated_hits_ADCs,              "total_sum_of_slice_associated_hits_ADCs/D");
  _passing_events_tree->Branch("u_plane_sum_of_slice_associated_hits_ADCs",                      &u_plane_sum_of_slice_associated_hits_ADCs,            "u_plane_sum_of_slice_associated_hits_ADCs/D");
  _passing_events_tree->Branch("v_plane_sum_of_slice_associated_hits_ADCs",                      &v_plane_sum_of_slice_associated_hits_ADCs,            "v_plane_sum_of_slice_associated_hits_ADCs/D");
  _passing_events_tree->Branch("y_plane_sum_of_slice_associated_hits_ADCs",                      &y_plane_sum_of_slice_associated_hits_ADCs,            "y_plane_sum_of_slice_associated_hits_ADCs/D");
  _passing_events_tree->Branch("total_number_of_slice_associated_hits",                          &total_number_of_slice_associated_hits,                "total_number_of_slice_associated_hits/I");
  _passing_events_tree->Branch("u_plane_number_of_slice_associated_hits",                        &u_plane_number_of_slice_associated_hits,              "u_plane_number_of_slice_associated_hits/I");
  _passing_events_tree->Branch("v_plane_number_of_slice_associated_hits",                        &v_plane_number_of_slice_associated_hits,              "v_plane_number_of_slice_associated_hits/I");
  _passing_events_tree->Branch("y_plane_number_of_slice_associated_hits",                        &y_plane_number_of_slice_associated_hits,              "y_plane_number_of_slice_associated_hits/I");
  _passing_events_tree->Branch("total_sum_of_slice_hits_in_vicinity_ADCs",                       &total_sum_of_slice_hits_in_vicinity_ADCs,             "total_sum_of_slice_hits_in_vicinity_ADCs/D");
  _passing_events_tree->Branch("u_plane_sum_of_slice_hits_in_vicinity_ADCs",                     &u_plane_sum_of_slice_hits_in_vicinity_ADCs,           "u_plane_sum_of_slice_hits_in_vicinity_ADCs/D");
  _passing_events_tree->Branch("v_plane_sum_of_slice_hits_in_vicinity_ADCs",                     &v_plane_sum_of_slice_hits_in_vicinity_ADCs,           "v_plane_sum_of_slice_hits_in_vicinity_ADCs/D");
  _passing_events_tree->Branch("y_plane_sum_of_slice_hits_in_vicinity_ADCs",                     &y_plane_sum_of_slice_hits_in_vicinity_ADCs,           "y_plane_sum_of_slice_hits_in_vicinity_ADCs/D");
  _passing_events_tree->Branch("total_number_of_slice_hits_in_vicinity",                         &total_number_of_slice_hits_in_vicinity,               "total_number_of_slice_hits_in_vicinity/I");
  _passing_events_tree->Branch("u_plane_number_of_slice_hits_in_vicinity",                       &u_plane_number_of_slice_hits_in_vicinity,             "u_plane_number_of_slice_hits_in_vicinity/I");
  _passing_events_tree->Branch("v_plane_number_of_slice_hits_in_vicinity",                       &v_plane_number_of_slice_hits_in_vicinity,             "v_plane_number_of_slice_hits_in_vicinity/I");
  _passing_events_tree->Branch("y_plane_number_of_slice_hits_in_vicinity",                       &y_plane_number_of_slice_hits_in_vicinity,             "y_plane_number_of_slice_hits_in_vicinity/I");
  _passing_events_tree->Branch("num_of_top_pandora_crossings",                                   &num_of_top_pandora_crossings,                         "num_of_top_pandora_crossings/I");  
  _passing_events_tree->Branch("num_of_bottom_pandora_crossings",                                &num_of_bottom_pandora_crossings,                      "num_of_bottom_pandora_crossings/I");
  _passing_events_tree->Branch("num_of_front_pandora_crossings",                                 &num_of_front_pandora_crossings,                       "num_of_front_pandora_crossings/I");
  _passing_events_tree->Branch("num_of_back_pandora_crossings",                                  &num_of_back_pandora_crossings,                        "num_of_back_pandora_crossings/I");
  _passing_events_tree->Branch("NC_channel",                                                     &NC_channel,                                           "NC_channel/I");
  _passing_events_tree->Branch("spline_fix_mcweight",                                            &spline_fix_mcweight,                                  "spline_fix_mcweight/D");
  _passing_events_tree->Branch("central_value_mcweight",                                         &central_value_mcweight,                               "central_value_mcweight/D");
  _passing_events_tree->Branch("rootino_fix_mcweight",                                           &rootino_fix_mcweight,                                 "rootino_fix_mcweight/D");
  _passing_events_tree->Branch("other_universe_mcweights",                                       &other_universe_mcweights,                             "other_universe_mcweights[100]/D");
  _passing_events_tree->Branch("axialff_mcweights",                                              &axialff_mcweights,                                    "axialff_mcweights[2]/D");
  _passing_events_tree->Branch("rpaccqe_mcweights",                                              &rpaccqe_mcweights,                                    "rpaccqe_mcweights[2]/D");
  _passing_events_tree->Branch("xsecshape_mcweights",                                            &xsecshape_mcweights,                                  "xsecshape_mcweights[2]/D");
  _passing_events_tree->Branch("num_muminus_tracks",                                             &num_muminus_tracks,                                   "num_muminus_tracks/I");
  _passing_events_tree->Branch("num_muplus_tracks",                                              &num_muplus_tracks,                                    "num_muplus_tracks/I");
  _passing_events_tree->Branch("num_piplus_tracks",                                              &num_piplus_tracks,                                    "num_piplus_tracks/I");
  _passing_events_tree->Branch("num_piminus_tracks",                                             &num_piminus_tracks,                                   "num_piminus_tracks/I");
  _passing_events_tree->Branch("num_pi0_tracks",                                                 &num_pi0_tracks,                                       "num_pi0_tracks/I");
  _passing_events_tree->Branch("num_proton_tracks",                                              &num_proton_tracks,                                    "num_proton_tracks/I");
  _passing_events_tree->Branch("num_electron_showers",                                           &num_electron_showers,                                 "num_electron_showers/I");
  _passing_events_tree->Branch("num_positron_showers",                                           &num_positron_showers,                                 "num_positron_showers/I");
  _passing_events_tree->Branch("num_photon_showers",                                             &num_photon_showers,                                   "num_photon_showers/I");
  _passing_events_tree->Branch("event_fails_new_two_track_selection",                            &event_fails_new_two_track_selection,                  "event_fails_new_two_track_selection/I");
  
  // Declare the branches for the proton hit identification tree.
  proton_hit_identification_tree = fs->make<TTree>("proton_hit_identification_tree", "A tree that contains the hits above threshold on the proton track");
  proton_hit_identification_tree->Branch("run",                                                  &_run,                                                 "run/I");
  proton_hit_identification_tree->Branch("subrun",                                               &_subrun,                                              "subrun/I");
  proton_hit_identification_tree->Branch("event",                                                &_event,                                               "event/I");
  proton_hit_identification_tree->Branch("neutrino_energy",                                      &neutrino_energy,                                      "neutrino_energy/D");
  proton_hit_identification_tree->Branch("num_points_above_threshold",                           &num_points_above_threshold,                           "num_points_above_threshold/I");
  proton_hit_identification_tree->Branch("num_u_plane_points_above_threshold_in_vertexing_loop", &num_u_plane_points_above_threshold_in_vertexing_loop, "num_u_plane_points_above_threshold_in_vertexing_loop/I");
  proton_hit_identification_tree->Branch("num_v_plane_points_above_threshold_in_vertexing_loop", &num_v_plane_points_above_threshold_in_vertexing_loop, "num_v_plane_points_above_threshold_in_vertexing_loop/I");
  proton_hit_identification_tree->Branch("num_y_plane_points_above_threshold_in_vertexing_loop", &num_y_plane_points_above_threshold_in_vertexing_loop, "num_y_plane_points_above_threshold_in_vertexing_loop/I");
  proton_hit_identification_tree->Branch("x_points_above_threshold",                             &x_points_above_threshold,                             "x_points_above_threshold[10000]/D");
  proton_hit_identification_tree->Branch("y_points_above_threshold",                             &y_points_above_threshold,                             "y_points_above_threshold[10000]/D");
  proton_hit_identification_tree->Branch("z_points_above_threshold",                             &z_points_above_threshold,                             "z_points_above_threshold[10000]/D");
  proton_hit_identification_tree->Branch("reconstructed_muon_KE",                                &reconstructed_muon_KE,                                "reconstructed_muon_KE/D");
  proton_hit_identification_tree->Branch("total_proton_KE",                                      &total_proton_KE,                                      "total_proton_KE/D");
  proton_hit_identification_tree->Branch("leading_proton_KE",                                    &leading_proton_KE,                                    "leading_proton_KE/D");
  proton_hit_identification_tree->Branch("muon_track_first_point_x",                             &muon_track_first_point_x,                             "muon_track_first_point_x/D");
  proton_hit_identification_tree->Branch("muon_track_first_point_y",                             &muon_track_first_point_y,                             "muon_track_first_point_y/D");
  proton_hit_identification_tree->Branch("muon_track_first_point_z",                             &muon_track_first_point_z,                             "muon_track_first_point_z/D");
  proton_hit_identification_tree->Branch("muon_track_last_point_x",                              &muon_track_last_point_x,                              "muon_track_last_point_x/D");
  proton_hit_identification_tree->Branch("muon_track_last_point_y",                              &muon_track_last_point_y,                              "muon_track_last_point_y/D");
  proton_hit_identification_tree->Branch("muon_track_last_point_z",                              &muon_track_last_point_z,                              "muon_track_last_point_z/D");
  proton_hit_identification_tree->Branch("vertex_location_x",                                    &vertex_location_x,                                    "vertex_location_x/D");
  proton_hit_identification_tree->Branch("vertex_location_y",                                    &vertex_location_y,                                    "vertex_location_y/D");
  proton_hit_identification_tree->Branch("vertex_location_z",                                    &vertex_location_z,                                    "vertex_location_z/D");
  proton_hit_identification_tree->Branch("other_end_location_x",                                 &other_end_location_x,                                 "other_end_location_x/D");
  proton_hit_identification_tree->Branch("other_end_location_y",                                 &other_end_location_y,                                 "other_end_location_y/D");
  proton_hit_identification_tree->Branch("other_end_location_z",                                 &other_end_location_z,                                 "other_end_location_z/D");
  proton_hit_identification_tree->Branch("proton_kinetic_energy_from_additional_tracks",         &proton_kinetic_energy_from_additional_tracks,         "proton_kinetic_energy_from_additional_tracks/D");
  proton_hit_identification_tree->Branch("x_momentum_component_sum_additional_tracks",           &x_momentum_component_sum_additional_tracks,           "x_momentum_component_sum_additional_tracks/D");
  proton_hit_identification_tree->Branch("y_momentum_component_sum_additional_tracks",           &y_momentum_component_sum_additional_tracks,           "y_momentum_component_sum_additional_tracks/D");
  proton_hit_identification_tree->Branch("z_momentum_component_sum_additional_tracks",           &z_momentum_component_sum_additional_tracks,           "z_momentum_component_sum_additional_tracks/D");
  proton_hit_identification_tree->Branch("proton_x_direction",                                   &proton_x_direction,                                   "proton_x_direction/I");
  proton_hit_identification_tree->Branch("proton_y_direction",                                   &proton_y_direction,                                   "proton_y_direction/I");
  proton_hit_identification_tree->Branch("proton_z_direction",                                   &proton_z_direction,                                   "proton_z_direction/I");
  proton_hit_identification_tree->Branch("muon_track_is_correctly_oriented",                     &muon_track_is_correctly_oriented,                     "muon_track_is_correctly_oriented/I");
  proton_hit_identification_tree->Branch("spline_fix_mcweight",                                  &spline_fix_mcweight,                                  "spline_fix_mcweight/D");

  // Declare the branches for the number of hits looping over tree. 
  num_events_looping_over_tree = fs->make<TTree>("num_events_looping_over_tree", "A tree that contains the number of the event that we are looking at");
  num_events_looping_over_tree->Branch("run",                                                  &_run,                                                 "run/I");
  num_events_looping_over_tree->Branch("subrun",                                               &_subrun,                                              "subrun/I");
  num_events_looping_over_tree->Branch("event",                                                &_event,                                               "event/I");

  if(_make_pida_csv) _csvfile.open ("pida_trklen.csv", std::ofstream::out | std::ofstream::trunc);
  if(_make_pida_csv) _csvfile << "pida,trklen,y" << std::endl;

  if(_make_ophit_csv) _csvfile2.open("ophit.csv", std::ofstream::out | std::ofstream::trunc);
  if(_make_ophit_csv) _csvfile2 << "ophit,opdet,time,pe" << std::endl;

  // Set the variables for the total numbers of events to 0.                                                                                                                                             
  total_num_events                                                        = 0;
  total_num_events_failing                                                = 0;
  total_num_events_passing                                                = 0;
  total_num_duplicate_events                                              = 0;

  // Set all of these variables equal to 0.                                                                                                                                                             
  total_num_events_failed_crt_veto                                        = 0;
  total_num_events_failed_beam_disc_flashes                               = 0;
  total_num_events_failed_beam_spill_flash                                = 0;
  total_num_events_failed_has_slices                                      = 0;
  total_num_events_failed_track_length                                    = 0;
  total_num_events_failed_has_slice_tagged_as_neutrino                    = 0;
  total_num_events_failed_fiducial_volume                                 = 0;
  total_num_events_failed_ntrack                                          = 0;
  total_num_events_failed_residuals_std_up                                = 0;
  total_num_events_failed_perc_used_hits_in_cluster                       = 0;
  total_num_events_failed_CRT_distance_cut                                = 0;
  total_num_events_failing_new_two_track_selection                        = 0;
  total_num_events_failing_3cm_containment_cut                            = 0;
  total_num_events_failing_length_of_tracks_from_vertex_requirement       = 0;
  total_num_events_failing_total_length_of_tracks_in_TPC_Object           = 0;
  total_num_events_failing_less_than_4_tracks_in_tpcobject_cut            = 0;
  total_num_events_failing_more_than_2_tracks_originating_from_vertex_cut = 0;
  total_num_events_failing_CRT_distance_cut                               = 0;

  fUPlaneThreshold = 5.32907;
  fVPlaneThreshold = 5.61058;
  fYPlaneThreshold = 3.34595;

  // Make sure you are only considering background.
  num_kdar_events_in_this_sample          = 0;
  num_kdar_events_passing_in_this_sample  = 0;

  num_event_we_are_on                     = 0;

  produces<std::vector<ubana::SelectionResult>>();
  produces<art::Assns<ubana::SelectionResult, ubana::TPCObject>>();

  // For the neutrino id filter
  produces<art::Assns<recob::Vertex, recob::Track>>();
  produces< art::Assns<recob::Vertex, recob::PFParticle>>();

  // Clear out the vector of duplicate event info.
  duplicate_runs.clear();
  duplicate_subruns.clear();
  duplicate_events.clear();

  event_counter  = 0;

}



void UBXSec::produce(art::Event & e) {

  std::cout << "Currently looping over event #" << event_counter << "." << std::endl;

  // Reset the event counter;
  event_counter++;

  // Reset if the event is a KDAR event.
  kdar_from_dump_event                                       = 0;

  // Reset the duplicate event info.
  is_duplicate_event                                         = false;

  // Reset the channel info.
  NC_channel                                                 = 0;

  // Reset all of the cut info.
  fails_new_two_track_selection                              = false;
  fails_CRT_distance_cut                                     = false;
  fails_3cm_containment_cut                                  = false;
  fails_fiducial_volume_cut                                  = false;

  // Reset the calorimetry/track orientation variables.
  u_plane_has_more_than_20_points                            = false;
  v_plane_has_more_than_20_points                            = false;
  y_plane_has_more_than_20_points                            = false;
  u_plane_has_track_forward                                  = false;
  v_plane_has_track_forward                                  = false;
  y_plane_has_track_forward                                  = false;
  no_calo_points_on_u_plane                                  = false;
  no_calo_points_on_v_plane                                  = false;
  no_calo_points_on_y_plane                                  = false;
  num_hits_above_threshold_u_plane                           = 0;
  num_hits_above_threshold_v_plane                           = 0;
  num_hits_above_threshold_y_plane                           = 0;

  u_plane_ticks_above_threshold.clear();
  u_plane_wires_above_threshold.clear();
  u_plane_charge_integrals_above_threshold.clear();

  v_plane_ticks_above_threshold.clear();
  v_plane_wires_above_threshold.clear();
  v_plane_charge_integrals_above_threshold.clear();

  y_plane_ticks_above_threshold.clear();
  y_plane_wires_above_threshold.clear();
  y_plane_charge_integrals_above_threshold.clear();

  // Set all of the proton momentum components to 0.
  x_momentum_component_sum_additional_tracks                  = 0;
  y_momentum_component_sum_additional_tracks                  = 0;
  z_momentum_component_sum_additional_tracks                  = 0;

  num_additional_tracks_originating_from_vertex               = 0;

  // Empty all of the variables for the hits on each of the planes.
  total_sum_of_slice_associated_hits_ADCs                     = 0.;
  u_plane_sum_of_slice_associated_hits_ADCs                   = 0.;
  v_plane_sum_of_slice_associated_hits_ADCs                   = 0.;
  y_plane_sum_of_slice_associated_hits_ADCs                   = 0.;

  total_number_of_slice_associated_hits                       = 0;
  u_plane_number_of_slice_associated_hits                     = 0;
  v_plane_number_of_slice_associated_hits                     = 0;
  y_plane_number_of_slice_associated_hits                     = 0;

  total_sum_of_slice_hits_in_vicinity_ADCs                    = 0.;
  u_plane_sum_of_slice_hits_in_vicinity_ADCs                  = 0.;
  v_plane_sum_of_slice_hits_in_vicinity_ADCs                  = 0.;
  y_plane_sum_of_slice_hits_in_vicinity_ADCs                  = 0.;

  total_number_of_slice_hits_in_vicinity                      = 0;
  u_plane_number_of_slice_hits_in_vicinity                    = 0;
  v_plane_number_of_slice_hits_in_vicinity                    = 0;
  y_plane_number_of_slice_hits_in_vicinity                    = 0;

  // Reset the variables for the number of crossings from pandora at each face of the TPC.
  num_of_top_pandora_crossings                                = 0;
  num_of_bottom_pandora_crossings                             = 0;
  num_of_front_pandora_crossings                              = 0;
  num_of_back_pandora_crossings                               = 0;
  
  // Reset the variables for orienting the proton.
  u_plane_wire_at_other_end                                   = 0;
  u_plane_other_end_tick                                      = 0;
  v_plane_wire_at_other_end                                   = 0;
  v_plane_other_end_tick                                      = 0;
  y_plane_wire_at_other_end                                   = 0; 
  y_plane_other_end_tick                                      = 0;

  u_plane_has_nonzero_nearby_points_and_in_the_top_two_planes = false;
  v_plane_has_nonzero_nearby_points_and_in_the_top_two_planes = false;
  y_plane_has_nonzero_nearby_points_and_in_the_top_two_planes = false;

  // Reset the variable for if the muon track is oriented correctly.
  muon_track_is_correctly_oriented           = 1;

  // Reset the variable for the two-track length cut.                                                                                                                                                    
  event_fails_new_two_track_selection                         = 0;

  std::vector< double > x_coordinates_above_threshold;
  std::vector< double > y_coordinates_above_threshold;
  std::vector< double > z_coordinates_above_threshold;

  // Include the object for the Space Charge effect service.
  spacecharge::SpaceChargeMicroBooNE const* sce =
    reinterpret_cast<spacecharge::SpaceChargeMicroBooNE const*>(lar::providerFrom<spacecharge::SpaceChargeService>());

  // Declare the object for calculating momentum.
  trkf::TrackMomentumCalculator p_calculator_from_length;

  // Declare the geometry object for use with the wire information.
  auto const* geom = ::lar::providerFrom<geo::Geometry>();

  // Load the wire objects that I will need to use to find charge in the vicinity of the vertex.                                                                                                       
  art::Handle<std::vector<recob::Hit> > hit_h;
  e.getByLabel("gaushit",hit_h);

  // make sure the hit objects look good                                                                                                                                                                  
  if(!hit_h.isValid()) {
    std::cerr<<"\033[93m[ERROR]\033[00m ... could not locate Hit Objects!"<<std::endl;
    throw std::exception();
  }

  // Load the wire objects that I will need to use to find charge in the vicinity of the vertex.
  art::Handle<std::vector<recob::Wire> > wire_h;
  e.getByLabel("butcher", wire_h);

  // make sure wire objects look good                                                                                                                                                                     
  if(!wire_h.isValid()) {
    std::cerr<<"\033[93m[ERROR]\033[00m ... could not locate Wire Objects!"<<std::endl;
    throw std::exception();
  }

  // Load the tracks from pandora.
  art::Handle<std::vector<recob::Track> > pandora_track_h;
  e.getByLabel("pandora",pandora_track_h);

  // make sure pandora tracks look good                                                                                                                                                                   
  if(!pandora_track_h.isValid()) {
    std::cerr<<"\033[93m[ERROR]\033[00m ... could not locate pandora track objects!"<<std::endl;
    throw std::exception();
  }

  // Load the hit information for these tracks.
  // Declare the association between tracks and hits as well.                                                                                                                                             
  art::FindMany<recob::Hit> trk_hit_assn_v(pandora_track_h, e, "pandora");

  // Load the calorimetry information for these tracks.
  art::FindMany<anab::Calorimetry> trk_calo_assn_v( pandora_track_h, e, "pandoracaliSCE" );

  // Use this to find the calorimetry objects that belong to the muon track in the event but without the space charge correction.                                                                          
  art::FindMany<anab::Calorimetry> trk_calo_assn_v_no_SCE_corrections( pandora_track_h, e, "pandoracali" );

  // Add in the information for removing the bad hits from the calorimetry information.
  auto const & track_list_ptr = e.getValidHandle<std::vector <recob::Track> >("pandora");
  art::FindManyP<recob::Hit, recob::TrackHitMeta> fmthm(track_list_ptr, e, "pandora");

  // Putting in the MCWeight info.                                                                                                                                                                     
  art::Handle<std::vector<evwgh::MCEventWeight> > eventweight_h;
  e.getByLabel( "eventweightSplines" , eventweight_h );

  if(!eventweight_h.isValid()) {
    std::cerr<<"\033[93m[ERROR]\033[00m ... could not locate 'eventweightSplines' EventWeight!"<<std::endl;
    throw std::exception();
  }

  for ( size_t i = 0; i < eventweight_h->size(); i++ ) {

    auto const& mc_weight = eventweight_h->at(i);

    std::map<std::string, std::vector<double> > weight_map = mc_weight.fWeight;

    const std::vector<double>& spline_fix_weights = weight_map.at("splines_general_Spline");

    spline_fix_mcweight = spline_fix_weights.front();

  }

  // Repeat this same procedure with the other systematic weights.
  art::Handle<std::vector<evwgh::MCEventWeight> > genie_eventweight_h;
  e.getByLabel( "eventweight" , genie_eventweight_h );

  if(!genie_eventweight_h.isValid()) {
    std::cerr<<"\033[93m[ERROR]\033[00m ... could not locate 'eventweight' EventWeight!"<<std::endl;
    throw std::exception();
  }

  
  for ( size_t i = 0; i < genie_eventweight_h->size(); i++ ) {

    auto const& mc_weight = genie_eventweight_h->at(i);

    std::map<std::string, std::vector<double> > weight_map = mc_weight.fWeight;

    for ( std::map<std::string, std::vector<double>>::iterator it = weight_map.begin(); it != weight_map.end(); it++ ) 
      {
	std::cout << it->first << std::endl;

      }

    const std::vector<double>& central_value_weights = weight_map.at("TunedCentralValue_UBGenie");

    central_value_mcweight = central_value_weights.front();

    const std::vector<double>& rootino_fix_weights = weight_map.at("RootinoFix_UBGenie");

    rootino_fix_mcweight = rootino_fix_weights.front();

    const std::vector<double>& universe_weights = weight_map.at("All_UBGenie");

    for ( size_t i = 0; i < universe_weights.size(); i++ ) {

      other_universe_mcweights[i] = universe_weights.at( i );

    }
    
    const std::vector<double>& axialff_weights = weight_map.at("AxFFCCQEshape_UBGenie");

    for( size_t i = 0;i < axialff_weights.size(); i++ ) {

      axialff_mcweights[i] = axialff_weights.at( i );

    }   

    const std::vector<double>& rpaccqe_weights = weight_map.at("RPA_CCQE_UBGenie");

    std::cout << "Size of 'RPA_CCQE' = " << rpaccqe_weights.size() << "." << std::endl;

    for( size_t i = 0;i < rpaccqe_weights.size(); i++ ) {

      rpaccqe_mcweights[i] = rpaccqe_weights.at( i );

    }   

    const std::vector<double>& xsecshape_weights = weight_map.at("XSecShape_CCMEC_UBGenie");

    for( size_t i = 0; i < xsecshape_weights.size(); i++ ) {

      xsecshape_mcweights[i] = xsecshape_weights.at( i );

    }   

  }

  // Set the four booleans to false;                                                                                                                                                                    
  bool kaon_parent_from_dump = false;
  bool kdar_energy           = false;
  bool CCNC_interaction      = false;
  bool neutrino_contained    = false;
       parent_fvz            = 0.;

  // load neutrino information.                                                                                                                                                                        
  if (_debug) { std::cout << "loading neutrino from producer generator." << std::endl; }
  art::Handle<std::vector<simb::MCTruth> > neutrino_h;
  e.getByLabel("generator", neutrino_h);

  // make sure MCTruth info looks good                                                                                                                                                                   
  if(!neutrino_h.isValid()) {
    std::cerr<<"\033[93m[ERROR]\033[00m ... could not locate Neutrino!"<<std::endl;
    throw std::exception();
  }

  // Declare a vector of pointers for the neutrino object.                   
  auto neutrino = neutrino_h->at(0).GetNeutrino();
  auto nu = neutrino.Nu();                      

  // Find the direction of the truth track to be used in the weighting. 
  truth_px_for_theta_for_weighting = nu.Px();
  truth_py_for_theta_for_weighting = nu.Py();
  truth_pz_for_theta_for_weighting = nu.Pz();

  // Variables
  TRotation RotDet2Beam;             // Rotations
  TVector3  detxyz, BeamCoords;      // Translations
  std::vector<double> rotmatrix;     // Inputs

  // input detector coordinates to translate
  detxyz = {truth_px_for_theta_for_weighting, truth_py_for_theta_for_weighting, truth_pz_for_theta_for_weighting};     

  // From beam to detector rotation matrix
  rotmatrix = {
    0.92103853804025681562, 0.022713504803924120662, 0.38880857519374290021,
    4.6254001262154668408e-05, 0.99829162468141474651, -0.058427989452906302359,
    -0.38947144863934973769, 0.053832413938664107345, 0.91946400794392302291 };

  // Return the TRotation
  TVector3 newX, newY, newZ;
  newX = TVector3(rotmatrix[0], rotmatrix[1], rotmatrix[2]);
  newY = TVector3(rotmatrix[3], rotmatrix[4], rotmatrix[5]);
  newZ = TVector3(rotmatrix[6], rotmatrix[7], rotmatrix[8]);

  RotDet2Beam.RotateAxes(newX, newY, newZ); // Return the TRotation now det to beam
  // RotDet2Beam.Invert(); // Invert back to the beam to det

  // Rotate to beam coords
  BeamCoords = RotDet2Beam * detxyz;

  TVector3 beam_dir = {0 , 0 , 1};
  truth_theta_angle_for_weighting = BeamCoords.Angle(beam_dir) * 180 / 3.1415926;

  std::vector<art::Ptr<simb::MCTruth> >NeutrinoVec;
  art::fill_ptr_vector(NeutrinoVec, neutrino_h);

  // Find out if there is a truth KDAR neutrino from the dump in the event.
  // load the truth MCFlux information.                                                                                                                                                                
  if (_debug) { std::cout << "loading MCFlux object from producer generator." << std::endl; }
  art::Handle<std::vector<simb::MCFlux> > parent_h;
  e.getByLabel("generator", parent_h);

  // make sure MCFlux info looks good                                                                                                                                                                  
  if(!parent_h.isValid()) {
    std::cerr<<"\033[93m[ERROR]\033[00m ... could not locate Parent!" << std::endl;
    throw std::exception();
  }

  // Declare a vector of pointers with the 'MCFlux' object.                                                                                                                                            
  std::vector<art::Ptr<simb::MCFlux> > ParentVec;
  art::fill_ptr_vector(ParentVec, parent_h);

  // Loop through the neutrino parents.                                                                                                                                                                   
  kaon_parent_from_dump = false;

  for (auto& parent : ParentVec){

    parent_fvz   = parent->fvz;
 
    if ( parent->fptype == 321 && parent->fvz > 72300 && parent->fvz < 72800 ) {
      kaon_parent_from_dump = true;
      break;
    }

  }

  if ( kaon_parent_from_dump == true ) {

    // loop through neutrinos themselves.                                                                                                                                                                 
    for (auto& neutrino : NeutrinoVec ) {

      // Reset 'kdar_energy', 'CCNC_interaction', and 'neutrino_contained' to true.                                                                                                                       
      kdar_energy                            = true;
      CCNC_interaction                       = true;
      neutrino_contained                     = true;

      // Unpack the neutrino object to find an MCParticle.                                                                                                                                                
      const simb::MCNeutrino& truth_neutrino = neutrino->GetNeutrino();
      const simb::MCParticle& truth_particle = truth_neutrino.Nu();

      // Unpack the coordinates for the vertex as well.                                                                                                                                                  
      nu_vtx_x_truth                         = truth_particle.Vx(0);
      nu_vtx_y_truth                         = truth_particle.Vy(0);
      nu_vtx_z_truth                         = truth_particle.Vz(0);
      
      // Only look at those events which decay via the charged current channel.                                                                                                                           
      if ( truth_neutrino.CCNC() != 0 )
	CCNC_interaction = false;

      // Continue if the particle's energy is outside of the correct peak.                                                                                                                                 
      if ( truth_particle.E(0) < 0.2355 || truth_particle.E(0) > 0.2356 ) {
	kdar_energy = false;
      }
      
      // Check to see if the neutrino vertex is contained.                                                                                                                                               
      if ( nu_vtx_x_truth < 0.0 || nu_vtx_x_truth > 256.4 || nu_vtx_y_truth < -116.5 || nu_vtx_y_truth > 116.5 || nu_vtx_z_truth < 0.0 || nu_vtx_z_truth > 1036.8 ) {
	neutrino_contained = false;
      }
      
      // Return 'true' if all of the conditions are met.                                                                                                                                                  
      if ( kdar_energy && CCNC_interaction && neutrino_contained ) {
	kdar_from_dump_event    = 1;
      }
      
    } // End of the loop over the neutrino candidates in the event.                         

  } // End of the conditional that the parent is a kaon at the dump.

  // Load in the MCTrack information.
  art::Handle<std::vector<sim::MCTrack> > mctrack_h;
  e.getByLabel("mcreco", mctrack_h);

  // make sure MCTrack info looks good                                                                                                                                                                 
  if(!mctrack_h.isValid()) {
    std::cerr<<"\033[93m[ERROR]\033[00m ... could not locate MCTrack!"<<std::endl;
    throw std::exception();
  }

  // Fill all of the truth info in a loop over the mctracks in the event.
  num_muminus_tracks   = 0;
  num_muplus_tracks    = 0;
  num_piplus_tracks    = 0;
  num_piminus_tracks   = 0;
  num_pi0_tracks       = 0;
  num_proton_tracks    = 0;

  total_proton_KE   = -1.;
  leading_proton_KE = 0.;

  truth_muon_length         = 0.;
  truth_muon_kinetic_energy = 0.;

  for ( size_t mctrack_iter = 0; mctrack_iter < mctrack_h->size(); mctrack_iter++ ) {

    // At the top of the loop, find out which type of track this is.
    if ( mctrack_h->at( mctrack_iter ).PdgCode() == 13 ) {

      num_muminus_tracks++;

      num_mctrack_points = mctrack_h->at( mctrack_iter ).size();

      if ( num_mctrack_points > 0 ) {

	std::cout << "Starting MCTrack Coordinates: x = " << mctrack_h->at( mctrack_iter ).at(0).X() << " cm y = " << mctrack_h->at( mctrack_iter ).at(0).Y() << " cm z = " << mctrack_h->at( mctrack_iter ).at(0).Z() << "." << std::endl;
	std::cout << "Ending MCTrack Coordinates: x = "<< mctrack_h->at( mctrack_iter ).at(mctrack_h->at( mctrack_iter ).size() - 1).X() << " cm y = " << mctrack_h->at( mctrack_iter ).at(mctrack_h->at( mctrack_iter ).size() - 1).Y() << " cm z = " << mctrack_h->at( mctrack_iter ).at(mctrack_h->at( mctrack_iter ).size() - 1).Z() << "." << std::endl;

	truth_muon_length = TMath::Sqrt( ( mctrack_h->at( mctrack_iter ).at(mctrack_h->at( mctrack_iter ).size() - 1).X() - mctrack_h->at( mctrack_iter ).at(0).X() ) * ( mctrack_h->at( mctrack_iter ).at(mctrack_h->at( mctrack_iter ).size() - 1).X() - mctrack_h->at( mctrack_iter ).at(0).X() ) + ( mctrack_h->at( mctrack_iter ).at(mctrack_h->at( mctrack_iter ).size() - 1).Y() - mctrack_h->at( mctrack_iter ).at(0).Y() ) * ( mctrack_h->at( mctrack_iter ).at(mctrack_h->at( mctrack_iter ).size() - 1).Y() - mctrack_h->at( mctrack_iter ).at(0).Y() ) + ( mctrack_h->at( mctrack_iter ).at(mctrack_h->at( mctrack_iter ).size() - 1).Z() - mctrack_h->at( mctrack_iter ).at(0).Z() ) * ( mctrack_h->at( mctrack_iter ).at(mctrack_h->at( mctrack_iter ).size() - 1).Z() - mctrack_h->at( mctrack_iter ).at(0).Z() ) );

	std::cout << "Truth muon length = " << truth_muon_length << "." << std::endl;

	truth_muon_kinetic_energy = ( TMath::Sqrt( ( p_calculator_from_length.GetTrackMomentum( truth_muon_length, 13 ) * 1000. ) * ( p_calculator_from_length.GetTrackMomentum( truth_muon_length, 13 ) * 1000. ) + 105.7 * 105.7 ) - 105.7 );

	std::cout << "Truth muon kinetic energy = " << truth_muon_kinetic_energy << " MeV." << std::endl;

	std::cout << "Kinetic energy of a 38 cm muon = " << ( TMath::Sqrt( ( p_calculator_from_length.GetTrackMomentum( 38.0, 13 ) * 1000. ) * ( p_calculator_from_length.GetTrackMomentum( 38.0, 13 ) * 1000. ) + 105.7 * 105.7 ) - 105.7 );

      }

    }

    if ( mctrack_h->at( mctrack_iter ).PdgCode() == -13 ) {

      num_muplus_tracks++;

    }

    if ( mctrack_h->at( mctrack_iter ).PdgCode() == 111 ) {

      num_pi0_tracks++;

    }

    if ( mctrack_h->at( mctrack_iter ).PdgCode() == 211 ) {

      num_piplus_tracks++;

    }

    if ( mctrack_h->at( mctrack_iter ).PdgCode() == -211 ) {

      num_piminus_tracks++;

    }

    if ( mctrack_h->at( mctrack_iter ).PdgCode() == 2212 ) {

      num_proton_tracks++;

    }

    if ( mctrack_h->at( mctrack_iter ).PdgCode() != 2212 || mctrack_h->at(mctrack_iter ).size() == 0 ) continue;

    total_proton_KE += ( mctrack_h->at(mctrack_iter ).at( 0 ).E() - 938.3 );

    if ( ( mctrack_h->at( mctrack_iter ).at( 0 ).E() - 938.3 ) > leading_proton_KE ) {

      leading_proton_KE = ( mctrack_h->at(mctrack_iter ).at( 0 ).E() - 938.3 );

    }

  }

  num_electron_showers = 0;
  num_positron_showers = 0;
  num_photon_showers   = 0;

  // Load in the MCShower information.
  art::Handle<std::vector<sim::MCShower> > mcshower_h;
  e.getByLabel("mcreco", mcshower_h);

  // make sure MCShower info looks good                                                                                                                                                                 
  if(!mcshower_h.isValid()) {
    std::cerr<<"\033[93m[ERROR]\033[00m ... could not locate MCShower!"<<std::endl;
    throw std::exception();
  }

  for ( size_t mcshower_iter = 0; mcshower_iter < mcshower_h->size(); mcshower_iter++ ) {

    if ( mcshower_h->at( mcshower_iter ).PdgCode() == 11 ) {

      num_electron_showers++;

    }

    if ( mcshower_h->at( mcshower_iter ).PdgCode() == -11 ) {

      num_positron_showers++;

    }

    if ( mcshower_h->at( mcshower_iter ).PdgCode() == 22  ) {

      num_photon_showers++;

    }

  }

  if(_debug) std::cout << "********** UBXSec starts" << std::endl;
  if(_debug) std::cout << "[UBXSec] Run: "           << e.id().run()    <<
                          ", subRun: "               << e.id().subRun() <<
                          ", event: "                << e.id().event()  << std::endl;

  if (_do_opdet_swap && e.isRealData()) {
    std::cout << "[UBXSec] WARNING!!! Swapping OpDets. I hope you know what you are doing." << std::endl;
  }

  // Instantiate the output
  std::unique_ptr< std::vector<ubana::SelectionResult>>                   selectionResultVector           (new std::vector<ubana::SelectionResult>);
  std::unique_ptr< art::Assns<ubana::SelectionResult, ubana::TPCObject>>  assnOutSelectionResultTPCObject (new art::Assns<ubana::SelectionResult, ubana::TPCObject>);

  std::unique_ptr<art::Assns<recob::Vertex, recob::Track>>      vertexTrackAssociations(new art::Assns<recob::Vertex, recob::Track>);
  std::unique_ptr<art::Assns<recob::Vertex, recob::PFParticle>> vertexPFParticleAssociations(new art::Assns<recob::Vertex, recob::PFParticle>);

  // Initialize the UBXSecEvent
  ubxsec_event->Init();

  // Assign values to the event information.
  _run    = ubxsec_event->run    = e.id().run();
  _subrun = ubxsec_event->subrun = e.id().subRun();
  _event  = ubxsec_event->event  = e.id().event();

  // See if this is a duplicate event.
  for ( size_t event_iter = 0; event_iter < duplicate_events.size(); event_iter++ ) {

    if ( fabs( duplicate_runs.at( event_iter ) - _run ) < 0.001 && fabs( duplicate_subruns.at( event_iter ) - _subrun ) < 0.001 && fabs( duplicate_events.at( event_iter ) - _event ) < 0.001 )
      is_duplicate_event = true;

  }

  _is_data = e.isRealData();
  _is_mc   = !_is_data;

  //::art::ServiceHandle<cheat::BackTracker> bt;
  ::art::ServiceHandle<geo::Geometry> geo;

  // Prepare the dead region finder
  //std::cout << "[UBXSec] Recreate channel status map" << std::endl;
  const lariov::ChannelStatusProvider& chanFilt = art::ServiceHandle<lariov::ChannelStatusService>()->GetProvider();
  for (unsigned int ch = 0; ch < 8256; ch++) {
    deadRegionsFinder.SetChannelStatus(ch, chanFilt.Status(ch));
  }
  //std::cout << "[UBXSec] Now force reload BWires" << std::endl;
  deadRegionsFinder.CreateBWires();

  // Use '_detp' to find 'efield' and 'temp'
  auto const* _detp = lar::providerFrom<detinfo::DetectorPropertiesService>();
  double efield = _detp -> Efield();
  double temp   = _detp -> Temperature();
  // Determine the drift velocity from 'efield' and 'temp'
  _drift_velocity = _detp -> DriftVelocity(efield,temp);
  if (_debug) std::cout << "[UBXSec] Using drift velocity = " << _drift_velocity << " cm/us, with E = " << efield << ", and T = " << temp << std::endl;

 lar_pandora::PFParticleVector pfParticleList;              //vector of PFParticles
  lar_pandora::LArPandoraHelper::CollectPFParticles(e, _pfp_producer, pfParticleList);

  // Collect tracks
  lar_pandora::TrackVector            allPfParticleTracks;
  lar_pandora::PFParticlesToTracks    pfParticleToTrackMap;
  lar_pandora::TracksToHits           trackToHitsMap;
  lar_pandora::LArPandoraHelper::CollectTracks(e, _pfp_producer, allPfParticleTracks, pfParticleToTrackMap);
  lar_pandora::LArPandoraHelper::CollectTracks(e, _pfp_producer, allPfParticleTracks, trackToHitsMap);

  // Collect showers
  lar_pandora::ShowerVector           _shower_v;
  lar_pandora::PFParticlesToShowers   _pfp_to_shower_map;
  lar_pandora::PFParticleVector _pfp_v;
  lar_pandora::LArPandoraHelper::CollectShowers(e, _pfp_producer, _shower_v, _pfp_to_shower_map);
  lar_pandora::PFParticlesToMetadata pfParticlesToMetadata;
  lar_pandora::LArPandoraHelper::CollectPFParticleMetadata(e, _pfp_producer, _pfp_v, pfParticlesToMetadata);
  lar_pandora::PFParticlesToHits recoParticlesToHits;
  lar_pandora::HitsToPFParticles recoHitsToParticles;

  lar_pandora::LArPandoraHelper::BuildPFParticleHitMaps(e, _pfp_producer, _spacepointLabel, recoParticlesToHits, recoHitsToParticles, lar_pandora::LArPandoraHelper::kUseDaughters, true);

  // Get TPCObjects from the Event
  art::Handle<std::vector<ubana::TPCObject>> tpcobj_h;
  e.getByLabel(_tpcobject_producer, tpcobj_h);
  if (!tpcobj_h.isValid()) {
    std::cout << "[UBXSec] Cannote locate ubana::TPCObject." << std::endl;
  }
  art::FindManyP<ubana::FlashMatch> tpcobjToFlashMatchAssns(tpcobj_h, e, _neutrino_flash_match_producer);
  art::FindManyP<recob::Track>      tpcobjToTrackAssns(tpcobj_h, e, _tpcobject_producer);
  art::FindManyP<recob::Shower>     tpcobjToShowerAssns(tpcobj_h, e, _tpcobject_producer);
  art::FindManyP<recob::PFParticle> tpcobjToPFPAssns(tpcobj_h, e, _tpcobject_producer);
  art::FindManyP<recob::Vertex>     tpcobjToVertexAssns(tpcobj_h, e, _tpcobject_producer);
  art::FindManyP<anab::CosmicTag>   tpcobjToCosmicTagAssns(tpcobj_h, e, _geocosmictag_producer);
  art::FindManyP<anab::CosmicTag>   tpcobjToConsistency(tpcobj_h, e, _candidateconsistency_producer);
  art::FindManyP<anab::CosmicTag>   tpcobjToStopMu(tpcobj_h, e, _cosmic_stopmu_tag_producer);

  // Get Tracks
  art::Handle<std::vector<recob::Track>> track_h;
  e.getByLabel(_pfp_producer,track_h);
  if (!track_h.isValid() || track_h->empty()) {
    std::cout << "[UBXSec] Track handle is not valid or empty." << std::endl;
    //throw std::exception();
  }

  std::vector<art::Ptr<recob::Track>> track_p_v;

  art::fill_ptr_vector(track_p_v, track_h);

  art::FindManyP<recob::PFParticle> pfp_from_track(track_h, e, _pfp_producer);

  art::FindManyP<anab::Calorimetry> calos_from_track(track_h, e, _calorimetry_producer);

  // Get PFP
  art::Handle<std::vector<recob::PFParticle> > pfp_h;
  e.getByLabel(_pfp_producer,pfp_h);
  art::FindManyP< larpandoraobj::PFParticleMetadata > pfPartToMetadataAssoc(pfp_h,e, _pfp_producer);

  if(!pfp_h.isValid()){
    std::cout << "[UBXSec] PFP product " << _pfp_producer << " not found..." << std::endl;
  }
  if(pfp_h->empty()) {
    std::cout << "[UBXSec] PFP "         << _pfp_producer << " is empty."    << std::endl;
  }
  std::vector<art::Ptr<recob::PFParticle>> pfp_v;
  art::fill_ptr_vector(pfp_v, pfp_h);
  art::FindManyP<recob::Track> tracks_from_pfp(pfp_h, e, _pfp_producer);

  ubxsec_event->n_pfp = ubxsec_event->n_pfp_primary = 0;
  for (size_t i = 0; i < pfp_h->size(); i++) {
    ubxsec_event->n_pfp++;
    if ((*pfp_h)[i].IsPrimary()&&abs((*pfp_h)[i].PdgCode()==14))
      ubxsec_event->n_pfp_primary++;

  }
 // Get Ghosts
  art::Handle<std::vector<ubana::MCGhost> > ghost_h;
  e.getByLabel(_mc_ghost_producer,ghost_h);
  if(!ghost_h.isValid()){
    std::cout << "[UBXSec] MCGhost product " << _mc_ghost_producer << " not found..." << std::endl;
  }
  art::FindManyP<ubana::MCGhost>   mcghost_from_pfp   (pfp_h,   e, _mc_ghost_producer);
  art::FindManyP<simb::MCParticle> mcpar_from_mcghost (ghost_h, e, _mc_ghost_producer);

    for (unsigned int i = 0; i < pfp_h->size(); ++i)
    {


      const std::vector< art::Ptr<larpandoraobj::PFParticleMetadata> > &pfParticleMetadataList(pfPartToMetadataAssoc.at(i));
        if (!pfParticleMetadataList.empty())
        {
            const art::Ptr<recob::PFParticle> pParticle(pfp_h, i);
            for (unsigned int j=0; j<pfParticleMetadataList.size(); ++j)
            {
                const art::Ptr<larpandoraobj::PFParticleMetadata> &pfParticleMetadata(pfParticleMetadataList.at(j));
                const larpandoraobj::PFParticleMetadata::PropertiesMap &pfParticlePropertiesMap(pfParticleMetadata->GetPropertiesMap());
                if (!pfParticlePropertiesMap.empty())

                for (larpandoraobj::PFParticleMetadata::PropertiesMap::const_iterator it = pfParticlePropertiesMap.begin(); it != pfParticlePropertiesMap.end(); ++it)
		  if(it->first=="TrackScore")
		    {
		      std::cout << " Found PFParticle " << pParticle->Self() << " with: " << std::endl;
                    std::cout << "  - " << it->first << " = " << it->second << std::endl;
		    ubxsec_event->pfp_trackscore.push_back(it->second);
		    }
            }
        }
    }

  ////
 // Get PID information
  art::FindManyP<anab::ParticleID> particleids_from_track (track_h, e, _particle_id_producer);
  if (!particleids_from_track.isValid()) {
    std::cout << "[UBXSec] anab::ParticleID is not valid." << std::endl;
  }
  // Fill a std::map Track->ParticleID
  std::map<art::Ptr<recob::Track>, art::Ptr<anab::ParticleID>> track_to_pid_map;
  for (auto track : track_p_v) {
    std::vector<art::Ptr<anab::ParticleID>> pids = particleids_from_track.at(track.key());
    if(pids.size() == 0)
      continue;
    for (auto pid : pids) {
      // Check that pid object contains plane-2 PIDA (because that's all we care about)
      std::vector<anab::sParticleIDAlgScores> AlgScoresVec = pid->ParticleIDAlgScores();
      for (size_t i_algscore=0; i_algscore<AlgScoresVec.size(); i_algscore++){
	anab::sParticleIDAlgScores AlgScore = AlgScoresVec.at(i_algscore);
	int planenum = UBPID::uB_getSinglePlane(AlgScore.fPlaneMask);
	if (AlgScore.fVariableType==anab::kPIDA && planenum==2){
	  track_to_pid_map[track] = pid;
	  continue;
	}
      }
    }
  }
   // Get MCSFitResult - Muon
   art::Handle<std::vector<recob::MCSFitResult> > mcsfitresult_mu_h;
   e.getByLabel(_mcsfitresult_mu_producer,mcsfitresult_mu_h);
   if(!mcsfitresult_mu_h.isValid()){
     std::cout << "[UBXSec] MCSFitResult product " << _mcsfitresult_mu_producer << " not found..." << std::endl;
   }
   std::vector<art::Ptr<recob::MCSFitResult>> mcsfitresult_mu_v;
   art::fill_ptr_vector(mcsfitresult_mu_v, mcsfitresult_mu_h);

  // pandoraCosmic PFPs (for cosmic removal studies)
  art::Handle<std::vector<recob::PFParticle>> pfp_cosmic_h;
  e.getByLabel("pandoraCosmic",pfp_cosmic_h);
  if(pfp_cosmic_h.isValid()){
    std::vector<art::Ptr<recob::PFParticle>> pfp_cosmic_v;
    art::fill_ptr_vector(pfp_cosmic_v, pfp_cosmic_h);
    ubxsec_event->n_primary_cosmic_pfp = 0;
    for (auto p : pfp_cosmic_v) {
      if (!p->IsPrimary()) continue;
      ubxsec_event->n_primary_cosmic_pfp++;
    }
  } else {
    std::cout << "[UBXSec] pandoraCosmic PFP product not found..." << std::endl;
  }

  // Load the DAQ time header.
  double evt_timeGPS_nsec = 0.;

  art::Handle<raw::DAQHeaderTimeUBooNE> rawHandle_DAQHeader;
  e.getByLabel("daq", rawHandle_DAQHeader);

  if(!rawHandle_DAQHeader.isValid()) {
      std::cerr << "\033[93m[ERROR]\033[00m ... could not locate DAQ header." << std::endl;
      throw std::exception();
  }
  
  raw::DAQHeaderTimeUBooNE const& my_DAQHeader(*rawHandle_DAQHeader);
  art::Timestamp evtTimeGPS = my_DAQHeader.gps_time();
  evt_timeGPS_nsec = evtTimeGPS.timeLow();
  
  // load CRT hits.                                                                                                                                                                                       
  art::Handle<std::vector<crt::CRTHit>> crthit_h;
  e.getByLabel("mixer", crthit_h);

  // make sure CRT hits look good.                                                                                                                                                                      
  if (!crthit_h.isValid()) {
    std::cerr << "\033[93m[ERROR]\033[00m ... could not locate CRT hits." << std::endl;
    throw std::exception();
  }

  // Flashes                                                                                                                                                                                             
  ::art::Handle<std::vector<recob::OpFlash>> beamflash_h;
  e.getByLabel(_opflash_producer_beam,beamflash_h);

  // Fill each of these variables in a loop over the flashes.
  flash_time                                                 = 0.;
  flash_PEs                                                  = 0.;
  flash_z                                                    = 0.;
  flash_y                                                    = 0.;
  flash_z_width                                              = 0.;
  flash_y_width                                              = 0.;

  _nflashes_in_beamgate                                      = beamflash_h->size();
  _nflashes_in_beamspill                                     = 0;
  _nflashes_in_beamspill_window_passing_filter_PE_cut        = 0;

  if ( _nflashes_in_beamgate > 0 ) { 
    
    flash_PEs = -1.;

    for ( size_t flash_iter = 0; flash_iter < beamflash_h->size(); flash_iter++ ) {

      if ( beamflash_h->at( flash_iter ).Time() > _beam_spill_start && beamflash_h->at( flash_iter ).Time() < _beam_spill_end ) { 

	_nflashes_in_beamspill++;

	if ( beamflash_h->at( flash_iter ).Time() > flash_PEs ) {

	  flash_time         = beamflash_h->at( flash_iter ).Time();
	  flash_PEs          = beamflash_h->at( flash_iter ).TotalPE();
	  flash_z            = beamflash_h->at( flash_iter ).ZCenter();
	  flash_y            = beamflash_h->at( flash_iter ).YCenter();
	  flash_z_width      = beamflash_h->at( flash_iter ).ZWidth();
	  flash_y_width      = beamflash_h->at( flash_iter ).YWidth();

	}

	// Cut for the entire KDAR filter.
	if ( beamflash_h->at( flash_iter ).TotalPE() > 50.0 ) 
	  _nflashes_in_beamspill_window_passing_filter_PE_cut++;

	// Cut for the CRT Veto.
	if ( beamflash_h->at( flash_iter ).TotalPE() > 10.0 )
	  _nflashes_in_beamgate_passing_beamspill_and_PE_cuts++;
    
      }
  
    }

  }

  // These two variables will determine if we reject the event based on CRT information.
  within_resolution                                          = 0;

  // Set the variables for the closest CRT hit time.                                                                                                                                                    
  double _dt_abs       = 100000.0;

  // Only enter this loop if there is at least one flash above threshold in the beamspill window.
  if ( _nflashes_in_beamgate_passing_beamspill_and_PE_cuts > 0 ) {
  
    // Loop over the CRT hits.                                                                                                                                                                             
    for (int j = 0; j < int( crthit_h->size() ); j++) {

      double _crt_time_temp = ((crthit_h->at(j).ts0_ns - evt_timeGPS_nsec + 68600) / 1000.);

      // Put in the 50 PEs CRT hit cut here.
      if ( crthit_h->at( j ).peshit < 50.0 ) 
	continue;
      
      if (fabs(flash_time - _crt_time_temp) < _dt_abs) {
	_dt_abs      = fabs(flash_time - _crt_time_temp);

	// set 'within_resolution' to 'true' and break the loop if 'closest_crt_diff' is less than fResolution.                                                            
	if (_dt_abs < 1.0) {
	  within_resolution = 1;
	}

      } // End of the comparison operator of the beam info to the CRT info.
    
    } // End of the loop over the CRT hits.

  } // End of conditional that there be at least one eligible flash.

  // Set the value of 'fails_crt_veto_cut'.
  if ( within_resolution == 1 && _nflashes_in_beamspill < 2 ) 
    ubxsec_event->fails_crt_veto_cut = 1;

  ubxsec_event->nbeamfls = beamflash_h->size();
  ubxsec_event->beamfls_pe.resize(ubxsec_event->nbeamfls);
  ubxsec_event->beamfls_time.resize(ubxsec_event->nbeamfls);
  ubxsec_event->beamfls_z.resize(ubxsec_event->nbeamfls);
  ubxsec_event->beamfls_spec.resize(ubxsec_event->nbeamfls);

  for (size_t n = 0; n < beamflash_h->size(); n++) {
    auto const& flash                          = (*beamflash_h)[n];
    ubxsec_event->beamfls_pe[n]                = flash.TotalPE();
    ubxsec_event->beamfls_time[n]              = flash.Time();
    ubxsec_event->beamfls_spec[n].resize(32);
    ubxsec_event->candidate_flash_time         = 0.;
    ubxsec_event->candidate_flash_z            = 0.;
    double min_pe = -1;
    for (unsigned int i = 0; i < 32; i++) {
      unsigned int opdet = geo->OpDetFromOpChannel(i);
      if (_do_opdet_swap && e.isRealData()) {
        opdet = _opdet_swap_map.at(opdet);
      }

      ubxsec_event->beamfls_spec[n][opdet] = flash.PE(i);

      if (ubxsec_event->beamfls_time[n] > _beam_spill_start && ubxsec_event->beamfls_time[n] < _beam_spill_end) {
        // Find largest flash above threshold
        if (flash.TotalPE() > _total_pe_cut && flash.TotalPE() > min_pe) {
          ubxsec_event->candidate_flash_time = flash.Time();
          ubxsec_event->candidate_flash_z    = flash.ZCenter();
          min_pe                             = flash.TotalPE();
        }
      }
    } // OpDet loop

    double Ycenter, Zcenter, Ywidth, Zwidth;
    GetFlashLocation(ubxsec_event->beamfls_spec[n], Ycenter, Zcenter, Ywidth, Zwidth);
    ubxsec_event->beamfls_z[n] = Zcenter;

  } // flash loop

  // Collecting GENIE particles
  if(_use_genie_info) {
    art::Handle< std::vector<simb::MCTruth> > mctruthListHandle;
    std::vector<art::Ptr<simb::MCTruth> > mclist;
    if (e.getByLabel("generator",mctruthListHandle))
      art::fill_ptr_vector(mclist, mctruthListHandle);

    this->PrintMC(mclist);
    ubxsec_event->fv            = 0;
    ubxsec_event->fv_sce        = 0;

    ubxsec_event->ccnc          = -1;
    ubxsec_event->nupdg         = -1;
    ubxsec_event->nu_e          = -1;
    ubxsec_event->lep_costheta  = -9999.;
    ubxsec_event->true_muon_mom = -9999.;

    ubxsec_event->ResizeGenieTruthVectors(mclist.size());
    for (size_t iList = 0; iList < mclist.size(); iList++) {

      // Check if the true neutrino vertex is in the FV
      double truth_nu_vtx[3] = {mclist[iList]->GetNeutrino().Nu().Vx(),
                                mclist[iList]->GetNeutrino().Nu().Vy(),
                                mclist[iList]->GetNeutrino().Nu().Vz()};

      if (_fiducial_volume.InFV(truth_nu_vtx)) {
        ubxsec_event->fv = 1;
      }

      // Save the vertex for all neutrinos
      ubxsec_event->tvtx_x.at(iList) = mclist[iList]->GetNeutrino().Nu().Vx();
      ubxsec_event->tvtx_y.at(iList) = mclist[iList]->GetNeutrino().Nu().Vy();
      ubxsec_event->tvtx_z.at(iList) = mclist[iList]->GetNeutrino().Nu().Vz();
 
      // Look at the space charge correction
      geo::Vector_t sce_corr    = _SCE->GetPosOffsets(geo::Point_t(mclist[iList]->GetNeutrino().Nu().Vx(),
								mclist[iList]->GetNeutrino().Nu().Vy(),
								mclist[iList]->GetNeutrino().Nu().Vz()));

      double g4Ticks            = _detector_clocks->TPCG4Time2Tick(mclist[iList]->GetNeutrino().Nu().T())
                                + _detector_properties->GetXTicksOffset(0,0,0)
                                - _detector_properties->TriggerOffset();


      double xOffset_mcc8       = _detector_properties->ConvertTicksToX(g4Ticks, 0, 0, 0) - sce_corr.X(); //(this was for mcc8)
      double xOffset            = sce_corr.X();
      double xOffset_mcc9       = _detector_properties->ConvertTicksToX(g4Ticks, 0, 0, 0);
      double yOffset            = sce_corr.Y();
      double zOffset            = sce_corr.Z();

      ubxsec_event->sce_corr_x  = xOffset;
      ubxsec_event->sce_corr_y  = yOffset;
      ubxsec_event->sce_corr_z  = zOffset;
      ubxsec_event->time_mcc9_x = xOffset_mcc9;
      ubxsec_event->timeminussce_mcc8_x = xOffset_mcc8;

      if (_fiducial_volume.InFV(mclist[iList]->GetNeutrino().Nu().Vx() + xOffset,
                                mclist[iList]->GetNeutrino().Nu().Vy() + yOffset,
                                mclist[iList]->GetNeutrino().Nu().Vz() + zOffset)) {
        ubxsec_event->fv_sce = 1;
      }

      int n_genie_particles         = 0;
      int n_genie_particles_charged = 0;
      for (int p = 0; p < mclist[iList]->NParticles(); p++) {
        const simb::MCParticle mc_par = mclist[iList]->GetParticle(p);
        if (mc_par.StatusCode() != 1) continue;
        n_genie_particles ++;
        const TParticlePDG* par_pdg   = _database_pdg->GetParticle(mc_par.PdgCode());
        if (!par_pdg) continue;
        if (par_pdg->Charge() == 0) continue;
        n_genie_particles_charged ++;
      }

      // Only save this if the true neutrino vertex is in the FV
      if (ubxsec_event->fv == 1) {
        ubxsec_event->ccnc            = mclist[iList]->GetNeutrino().CCNC();
        ubxsec_event->mode            = mclist[iList]->GetNeutrino().Mode();
        ubxsec_event->nupdg           = mclist[iList]->GetNeutrino().Nu().PdgCode();
        ubxsec_event->nu_e            = mclist[iList]->GetNeutrino().Nu().E();
        ubxsec_event->lep_costheta    = mclist[iList]->GetNeutrino().Lepton().Pz() / mclist[iList]->GetNeutrino().Lepton().P();
        ubxsec_event->lep_phi         = UBXSecHelper::GetPhi(mclist[iList]->GetNeutrino().Lepton().Px(),
                                                           mclist[iList]->GetNeutrino().Lepton().Py(),
                                                           mclist[iList]->GetNeutrino().Lepton().Pz());
        ubxsec_event->genie_mult      = n_genie_particles;
        ubxsec_event->genie_mult_ch   = n_genie_particles_charged;
      }

      ubxsec_event->nsignal = 0;
      if(ubxsec_event->nupdg==14 && ubxsec_event->ccnc==0 && ubxsec_event->fv==1) ubxsec_event->nsignal=1;

      // Also save muon momentum if is CC interaction
      ubxsec_event->true_muon_mom = -9999.;
      if (mclist[iList]->GetNeutrino().CCNC() == 0) {
        for (int p = 0; p < mclist[iList]->NParticles(); p++) {
          auto const & mcp = mclist[iList]->GetParticle(p);
          if (mcp.Mother() != 0) continue;
          if (mcp.PdgCode() != 13) continue;
          ubxsec_event->true_muon_mom = mcp.P();
        }
      }

    }
  }

  ubxsec_event->is_signal = false;
  if (ubxsec_event->ccnc == 0 && ubxsec_event->nupdg == 14 && ubxsec_event->fv == 1) {
    ubxsec_event->is_signal = true;
  }

  lar_pandora::PFParticlesToSpacePoints pfp_to_spacept;
  lar_pandora::SpacePointsToHits spacept_to_hits;

  lar_pandora::PFParticleVector temp2;
  lar_pandora::LArPandoraHelper::CollectPFParticles(e, _pfp_producer, temp2, pfp_to_spacept);

  lar_pandora::SpacePointVector temp3;
  lar_pandora::LArPandoraHelper::CollectSpacePoints(e, _pfp_producer, temp3, spacept_to_hits);

  art::Handle<std::vector<recob::OpHit>> ophit_h;
  e.getByLabel("ophitBeam", ophit_h);

  art::Handle<std::vector<recob::OpHit>> ophit_cosmic_h;
  e.getByLabel("ophitCosmic", ophit_cosmic_h);

  // Check if the muon is reconstructed
  for (auto p : pfp_v) {
    auto mcghosts = mcghost_from_pfp.at(p.key());
    if (mcghosts.size() > 0) {
      art::Ptr<simb::MCParticle> mcpar = mcpar_from_mcghost.at(mcghosts.at(0).key()).at(0);
      const auto mc_truth = UBXSecHelper::TrackIDToMCTruth(e, "largeant", mcpar->TrackId());
      if (mc_truth) {
        if (mc_truth->Origin() == simb::kBeamNeutrino
            && mcpar->PdgCode() == 13 && mcpar->Mother() == 0) {
          ubxsec_event->muon_is_reco = true;
        }
      }
    }
  }

  std::vector<lar_pandora::TrackVector     > track_v_v;
  std::vector<lar_pandora::ShowerVector    > shower_v_v;
  std::vector<lar_pandora::PFParticleVector> pfp_v_v;
  for (size_t slice = 0; slice < tpcobj_h->size(); slice++) {
    track_v_v.push_back(tpcobjToTrackAssns.at(slice));
    shower_v_v.push_back(tpcobjToShowerAssns.at(slice));
    pfp_v_v.push_back(tpcobjToPFPAssns.at(slice));
  }

  ubxsec_event->nslices = tpcobj_h->size();
  ubxsec_event->ResizeVectors(tpcobj_h->size());

  ubxsec_event->n_tpcobj_nu_origin     = 0;
  ubxsec_event->n_tpcobj_cosmic_origin = 0;

  std::vector<art::Ptr<recob::Track>> muon_candidate_track_per_slice_v;
  std::vector<art::Ptr<recob::PFParticle>> muon_candidate_pfparticle_per_slice_v;
  std::vector<art::Ptr<recob::Vertex>> neutrino_candidate_vertex_per_slice_v;
  muon_candidate_track_per_slice_v.resize(ubxsec_event->nslices);
  muon_candidate_pfparticle_per_slice_v.resize(ubxsec_event->nslices);
  neutrino_candidate_vertex_per_slice_v.resize(ubxsec_event->nslices);

  //
  // THIS IS THE MAIN LOOP OVER THE
  // TPCOBJECTS CANDIDATES IN THIS EVENT
  //

  for (unsigned int slice = 0; slice < tpcobj_h->size(); slice++){

    ubana::TPCObject tpcobj = (*tpcobj_h)[slice];

    ubxsec_event->slc_npfp[slice]    = tpcobj.GetNPFP();
    ubxsec_event->slc_ntrack[slice]  = tpcobj.GetNTracks();
    ubxsec_event->slc_nshower[slice] = tpcobj.GetNShowers();

    // Slice origin
    ubxsec_event->slc_origin[slice]  = tpcobj.GetOrigin();

    if (tpcobj.GetOrigin() == ubana::kBeamNeutrino || tpcobj.GetOrigin() == ubana::kMixed)
      ubxsec_event->n_tpcobj_nu_origin ++;
    else
      ubxsec_event->n_tpcobj_cosmic_origin ++;

    // Slice origin extra
    ubxsec_event->slc_origin_extra[slice] = tpcobj.GetOriginExtra();

    // Containment
    ubxsec_event->slc_iscontained[slice] = UBXSecHelper::TracksAreContained(tpcobj.GetTracks());

    // Cosmic tagging: stopping mu
    ubxsec_event->slc_stopmu_tagged[slice]=false;
    std::vector<art::Ptr<anab::CosmicTag>> stopmu_tags = tpcobjToStopMu.at(slice);
    if ( stopmu_tags.size() == 1 ) {
      auto smt = stopmu_tags.at(0);
      if (smt->CosmicType() == anab::CosmicTagID_t::kGeometry_Y){
	ubxsec_event->slc_stopmu_tagged[slice] = true;
      }
    }

    // Reco vertex
    double reco_nu_vtx_raw[3];
    recob::Vertex tpcobj_nu_vtx = tpcobj.GetVertex();
    //tpcobj_nu_vtx.XYZ(reco_nu_vtx_raw);
    std::vector<art::Ptr<recob::Vertex>> recob_vtx_v = tpcobjToVertexAssns.at(slice);
    if (recob_vtx_v.size() > 0) {
      recob_vtx_v.at(0)->XYZ(reco_nu_vtx_raw);
      neutrino_candidate_vertex_per_slice_v.at(slice) = recob_vtx_v.at(0);
    } else {
      reco_nu_vtx_raw[0] = reco_nu_vtx_raw[1] = reco_nu_vtx_raw[2] = -9999;
    }

    // X position correction (time offset)
    double reco_nu_vtx[3];
    UBXSecHelper::GetTimeCorrectedPoint(reco_nu_vtx_raw, reco_nu_vtx, ubxsec_event->candidate_flash_time, _drift_velocity);

    // Set these equal to the position of the vertex as I define it.
    ubxsec_event->slc_nuvtx_fv[slice] = (_fiducial_volume.InFV(reco_nu_vtx) ? 1 : 0);

    // Vertex resolution
    if (ubxsec_event->slc_origin[slice] == ubana::kBeamNeutrino) {
      ubxsec_event->vtx_resolution = sqrt( pow(ubxsec_event->slc_nuvtx_y[slice]-ubxsec_event->tvtx_y[0], 2) + pow(ubxsec_event->slc_nuvtx_z[slice]-ubxsec_event->tvtx_z[0], 2) );
    }

    // Multiplicity
    int p, t, s;
    tpcobj.GetMultiplicity(p, t, s);
    ubxsec_event->slc_mult_pfp[slice]             = p;
    ubxsec_event->slc_mult_track[slice]           = t;
    ubxsec_event->slc_mult_shower[slice]          = s;
    ubxsec_event->slc_mult_track_tolerance[slice] = tpcobj.GetNTracksCloseToVertex(_tolerance_track_multiplicity);

    // Candidate Consistency
    ubxsec_event->slc_consistency[slice]                    = true;
    std::vector<art::Ptr<anab::CosmicTag>> consistency_tags = tpcobjToConsistency.at(slice);
    auto ct = consistency_tags.at(0);
    if (ct->CosmicType() != anab::CosmicTagID_t::kNotTagged) {
      ubxsec_event->slc_consistency[slice] = false;
    }
    ubxsec_event->slc_consistency_score[slice] = ct->CosmicScore();

    // Neutrino Flash match: is this slice selected by the external flash matching as the neutrino slice?
    // Get PFPs from TPCObject
    auto pfps_from_tpcobj = tpcobjToPFPAssns.at(slice);
    // Get primary PFP (as in: the one that has no parent). Check PDG value of primary PFP. If this is the neutrino slice, it will have a neutrino PDG code (12, 14, or 16). If not, it will have a muon PDG code
    bool isnuslc = false;
    for (auto pfp : pfps_from_tpcobj){
      if (pfp->IsPrimary()){
        if (TMath::Abs(pfp->PdgCode())==12 || TMath::Abs(pfp->PdgCode()==14) || TMath::Abs(pfp->PdgCode()==16)){
          isnuslc = true;
          break;
        }
      }
    }
    ubxsec_event->slc_is_nu[slice]   = isnuslc;
  
    // Hits
    int nhits_u, nhits_v, nhits_w;
    UBXSecHelper::GetNumberOfHitsPerPlane(e, _pfp_producer, track_v_v[slice], nhits_u, nhits_v, nhits_w);
    ubxsec_event->slc_nhits_u[slice] = nhits_u;
    ubxsec_event->slc_nhits_v[slice] = nhits_v;
    ubxsec_event->slc_nhits_w[slice] = nhits_w;

    // Longest track and check boundary
    recob::Track lt;
    if (UBXSecHelper::GetLongestTrackFromTPCObj(track_v_v[slice], lt)){
      ubxsec_event->slc_longesttrack_length[slice]      = lt.Length();
      ubxsec_event->slc_longesttrack_phi[slice]         = UBXSecHelper::GetCorrectedPhi(lt, tpcobj_nu_vtx);
      ubxsec_event->slc_longesttrack_theta[slice]       = UBXSecHelper::GetCorrectedCosTheta(lt, tpcobj_nu_vtx);
      ubxsec_event->slc_longesttrack_iscontained[slice] = UBXSecHelper::TrackIsContained(lt);
      int vtx_ok;
      ubxsec_event->slc_crosses_top_boundary[slice]     = (UBXSecHelper::IsCrossingTopBoundary(lt, vtx_ok) ? 1 : 0);
    } else {
      ubxsec_event->slc_longesttrack_length[slice]      = -9999;
    }

    // Longest shower
    recob::Shower ls;
    if (UBXSecHelper::GetLongestShowerFromTPCObj(shower_v_v[slice], ls)) {
      ubxsec_event->slc_longestshower_length[slice]    = ls.Length();
      ubxsec_event->slc_longestshower_openangle[slice] = ls.OpenAngle();
      ubxsec_event->slc_longestshower_startx[slice]    = ls.ShowerStart().X();
      ubxsec_event->slc_longestshower_starty[slice]    = ls.ShowerStart().Y();
      ubxsec_event->slc_longestshower_startz[slice]    = ls.ShowerStart().Z();
      ubxsec_event->slc_longestshower_phi[slice]       = UBXSecHelper::GetPhi(ls.Direction());
      ubxsec_event->slc_longestshower_theta[slice]     = UBXSecHelper::GetCosTheta(ls.Direction());
    }

    // Track quality
    ubxsec_event->slc_kalman_chi2[slice] = -9999;

  bool goodTrack = true;
    for (auto trk : track_v_v[slice]) {
      if (deadRegionsFinder.NearDeadReg2P( (trk->Vertex()).Y(), (trk->Vertex()).Z(), _minimumDistDeadReg )  ||
          deadRegionsFinder.NearDeadReg2P( (trk->End()).Y(),    (trk->End()).Z(),    _minimumDistDeadReg )  ||
          deadRegionsFinder.NearDeadRegCollection(trk->Vertex().Z(), _minimumDistDeadReg) ||
          deadRegionsFinder.NearDeadRegCollection(trk->End().Z(),    _minimumDistDeadReg) ||
          !UBXSecHelper::TrackPassesHitRequirment(e, _pfp_producer, trk, _minimumHitRequirement) ) {
        goodTrack = false;
        break;
      }
    }

    if (goodTrack) ubxsec_event->slc_passed_min_track_quality[slice] = true;
    else ubxsec_event->slc_passed_min_track_quality[slice] = false;

    // Vertex quality
    recob::Vertex slice_vtx = tpcobj.GetVertex();
    double slice_vtx_xyz[3];
    slice_vtx.XYZ(slice_vtx_xyz);
    ubxsec_event->slc_passed_min_vertex_quality[slice] = true;
    if (deadRegionsFinder.NearDeadReg2P(slice_vtx_xyz[1], slice_vtx_xyz[2], _minimumDistDeadReg))
      ubxsec_event->slc_passed_min_vertex_quality[slice] = false;

    // Channel status
    ubxsec_event->slc_nuvtx_closetodeadregion_u[slice] = (UBXSecHelper::PointIsCloseToDeadRegion(reco_nu_vtx, 0) ? 1 : 0);
    ubxsec_event->slc_nuvtx_closetodeadregion_v[slice] = (UBXSecHelper::PointIsCloseToDeadRegion(reco_nu_vtx, 1) ? 1 : 0);
    ubxsec_event->slc_nuvtx_closetodeadregion_w[slice] = (UBXSecHelper::PointIsCloseToDeadRegion(reco_nu_vtx, 2) ? 1 : 0);

    // Vertex check
    ubxsec::VertexCheck vtxCheck(track_v_v[slice], slice_vtx);
    ubxsec_event->slc_vtxcheck_angle[slice] = vtxCheck.AngleBetweenLongestTracks();

    // OpHits                                                                                                                                                                                           
    std::vector<ubxsec::Hit3D_t> hit3d_v;
    hit3d_v.clear();
    for (auto pfp : pfp_v_v[slice]) {
      auto iter = pfp_to_spacept.find(pfp);
      if ( iter == pfp_to_spacept.end() )
	continue;
      for (auto sp_pt : (iter->second)) {
        auto iter2 = spacept_to_hits.find(sp_pt);
        if (iter2 == spacept_to_hits.end()) {
          continue;
        }
        // Save sp_pt position and hit charge for all the sp_pt you have                                                                                                                              
	auto hit = iter2->second;
	ubxsec::Hit3D_t thishit;
        thishit.x = sp_pt->XYZ()[0];
        thishit.y = sp_pt->XYZ()[1];
        thishit.z = sp_pt->XYZ()[2];
        thishit.q = hit->Integral();
        hit3d_v.emplace_back(thishit);
      }
    }

    // Now construct average position
    double sumx = 0, sumy = 0, sumz = 0;
    double totq = 0;
    for (auto hit3d : hit3d_v) {
      sumx += hit3d.q * hit3d.x;
      sumy += hit3d.q * hit3d.y;
      sumz += hit3d.q * hit3d.z;

      totq += hit3d.q;
    }
    double charge_center[3] = {sumx / totq, sumy / totq, sumz / totq};

    int this_opch = UBXSecHelper::GetClosestPMT(charge_center);

    // Look at the opHits from this pmt
    int n_intime_ophits = 0;
    double n_intime_pe  = 0;
    for (size_t oh = 0; oh < ophit_h->size(); oh++) {
      auto const & ophit = (*ophit_h)[oh];
      size_t opdet       = geo->OpDetFromOpChannel(ophit.OpChannel());
      if(_make_ophit_csv) _csvfile2 << oh << "," << opdet << "," << ophit.PeakTime() << "," << _pecalib.BeamPE(opdet,ophit.Area(),ophit.Amplitude()) << std::endl;
      if (ophit.OpChannel() != this_opch) continue;
      if (ophit.PeakTime() > _beam_spill_start && ophit.PeakTime() < _beam_spill_end) {
        n_intime_ophits ++;
        n_intime_pe     += _pecalib.BeamPE(opdet,ophit.Area(),ophit.Amplitude());
      }
    } // end loop ophit

    ubxsec_event->slc_n_intime_pe_closestpmt[slice] = n_intime_pe;

    // Muon Candidate
    _muon_finder.Reset();
    _muon_finder.SetTracks(track_v_v[slice]);
    _muon_finder.SetTrackToPIDMap(track_to_pid_map);
    art::Ptr<recob::Track> candidate_track;
    // cout<<candidate
    bool muon_cand_exists = _muon_finder.GetCandidateTrack(candidate_track);
    if (muon_cand_exists) {

      bool fully_contained = _fiducial_volume.InFV(candidate_track->Vertex<TVector3>(), candidate_track->End<TVector3>());
      
      ubxsec_event->slc_muoncandidate_exists[slice]    = true;
      ubxsec_event->slc_muoncandidate_contained[slice] = fully_contained;
      ubxsec_event->slc_muoncandidate_length[slice]    = candidate_track->Length();
      ubxsec_event->slc_muoncandidate_phi[slice]       = UBXSecHelper::GetCorrectedPhi((*candidate_track), tpcobj_nu_vtx);
      ubxsec_event->slc_muoncandidate_theta[slice]     = UBXSecHelper::GetCorrectedCosTheta((*candidate_track), tpcobj_nu_vtx);
      ubxsec_event->slc_muoncandidate_theta_xz[slice]  = UBXSecHelper::GetCorrectedCosThetaXZ((*candidate_track), tpcobj_nu_vtx);
      ubxsec_event->slc_muoncandidate_theta_yz[slice]  = UBXSecHelper::GetCorrectedCosThetaYZ((*candidate_track), tpcobj_nu_vtx);

      ubxsec_event->slc_muoncandidate_mom_range[slice] = _trk_mom_calculator.GetTrackMomentum(candidate_track->Length(), 13);
      ubxsec_event->slc_muoncandidate_mom_mcs[slice]   = _trk_mom_calculator.GetMomentumMultiScatterLLHD(candidate_track);
     
      // For MCS first check the track direction is rigth

      TVector3 temp(reco_nu_vtx[0], reco_nu_vtx[1], reco_nu_vtx[2]);
      bool track_direction_correct = (candidate_track->Vertex<TVector3>() - temp).Mag() < (candidate_track->End<TVector3>() - temp).Mag();

      if (track_direction_correct) {
	ubxsec_event->slc_muoncandidate_mom_mcs[slice] = mcsfitresult_mu_v.at(candidate_track.key())->fwdMomentum();
      } else {
	ubxsec_event->slc_muoncandidate_mom_mcs[slice] = mcsfitresult_mu_v.at(candidate_track.key())->bwdMomentum();
      }

      // Also see if the track is recon going downwards (for cosmic studies)
      bool track_going_down = candidate_track->Vertex().Y() > candidate_track->End().Y();
      
      // Look at calorimetry for the muon candidate
      std::vector<art::Ptr<anab::Calorimetry>> calos = calos_from_track.at(candidate_track.key());
      
      ubxsec_event->slc_muoncandidate_dqdx_vector[slice]      = UBXSecHelper::GetDqDxVector(calos);
      ubxsec_event->slc_muoncandidate_dqdx_trunc[slice]       = UBXSecHelper::GetDqDxTruncatedMean(calos);
      ubxsec_event->slc_muoncandidate_dqdx_u_trunc[slice]     = UBXSecHelper::GetDqDxTruncatedMean(calos, 0);
      ubxsec_event->slc_muoncandidate_dqdx_v_trunc[slice]     = UBXSecHelper::GetDqDxTruncatedMean(calos, 1);

      ubxsec_event->slc_muoncandidate_res_range_y[slice]      = UBXSecHelper::GetResRange(calos,2);
      ubxsec_event->slc_muoncandidate_res_range_u[slice]      = UBXSecHelper::GetResRange(calos, 0);
      ubxsec_event->slc_muoncandidate_res_range_v[slice]      = UBXSecHelper::GetResRange(calos, 1);

      ubxsec_event->slc_muoncandidate_dEdx_y[slice]           = UBXSecHelper::GetdEdx(calos,2);
      ubxsec_event->slc_muoncandidate_dEdx_u[slice]           = UBXSecHelper::GetdEdx(calos, 0);
      ubxsec_event->slc_muoncandidate_dEdx_v[slice]           = UBXSecHelper::GetdEdx(calos, 1);

      ubxsec_event->slc_muoncandidate_dQdx_y[slice]           = UBXSecHelper::GetdQdx(calos,2);
      ubxsec_event->slc_muoncandidate_dQdx_u[slice]           = UBXSecHelper::GetdQdx(calos, 0);
      ubxsec_event->slc_muoncandidate_dQdx_v[slice]           = UBXSecHelper::GetdQdx(calos, 1);

      ubxsec_event->slc_muoncandidate_mip_consistency[slice]  = _muon_finder.MIPConsistency(ubxsec_event->slc_muoncandidate_dqdx_trunc[slice],
                                                                                           ubxsec_event->slc_muoncandidate_length[slice]);
      ubxsec_event->slc_muoncandidate_mip_consistency2[slice] = _muon_finder.SVMPredict(ubxsec_event->slc_muoncandidate_dqdx_trunc[slice],
                                                                                        ubxsec_event->slc_muoncandidate_length[slice]);

      // Get the related PFP
      art::Ptr<recob::PFParticle> candidate_pfp = pfp_from_track.at(candidate_track.key()).at(0);
      const auto mcghosts = mcghost_from_pfp.at(candidate_pfp.key());
      if (mcghosts.size() > 0) {
        art::Ptr<simb::MCParticle> mcpar = mcpar_from_mcghost.at(mcghost_from_pfp.at(candidate_pfp.key()).at(0).key()).at(0);
        const auto mc_truth = UBXSecHelper::TrackIDToMCTruth(e, "largeant", mcpar->TrackId());
        ubxsec_event->slc_muoncandidate_truepdg[slice] = mcpar->PdgCode();
        if (mc_truth) {

          // Check the true origin of the candidate PFP
          if (mc_truth->Origin() == simb::kBeamNeutrino) {
            ubxsec_event->slc_muoncandidate_trueorigin[slice] = ubana::kBeamNeutrino;
          } else if (mc_truth->Origin() == simb::kCosmicRay) {
            ubxsec_event->slc_muoncandidate_trueorigin[slice] = ubana::kCosmicRay;
          }

          // Now make momentum distributions
          if (mc_truth->Origin() == simb::kBeamNeutrino && mcpar->PdgCode() == 13 && mcpar->Mother() == 0) {
            _h_mom_true_mcs->Fill(mcpar->P(), ubxsec_event->slc_muoncandidate_mom_mcs[slice]);
            if (fully_contained) {
              _mom_true_contained          = mcpar->P();
              _mom_mcs_contained           = ubxsec_event->slc_muoncandidate_mom_mcs[slice];
              _mom_range_contained         = ubxsec_event->slc_muoncandidate_mom_range[slice];
              _mom_tree_contained->Fill();
              _h_mom_true_mcs_contained->Fill(mcpar->P(), ubxsec_event->slc_muoncandidate_mom_mcs[slice]);
              _h_mom_true_range_contained->Fill(mcpar->P(), ubxsec_event->slc_muoncandidate_mom_range[slice]);
              _h_mom_range_mcs_contained->Fill(ubxsec_event->slc_muoncandidate_mom_range[slice],
                                               ubxsec_event->slc_muoncandidate_mom_mcs[slice]);
            } else {
              _mom_true_uncontained = mcpar->P();
              _mom_mcs_uncontained = ubxsec_event->slc_muoncandidate_mom_mcs[slice];
              _mom_tree_uncontained->Fill();
              _h_mom_true_mcs_uncontained->Fill(mcpar->P(), ubxsec_event->slc_muoncandidate_mom_mcs[slice]);
            }
          }
          if (mc_truth->Origin() == simb::kCosmicRay && (mcpar->PdgCode() == 13 || mcpar->PdgCode() == -13)) {
            _mom_cosmic_true = mcpar->P();
            _mom_cosmic_mcs = ubxsec_event->slc_muoncandidate_mom_mcs[slice];
            _mom_cosmic_mcs_downforced = track_going_down ?   mcsfitresult_mu_v.at(candidate_track.key())->fwdMomentum()
                                                            : mcsfitresult_mu_v.at(candidate_track.key())->bwdMomentum();
            _mom_cosmic_range = ubxsec_event->slc_muoncandidate_mom_range[slice];
            _mom_cosmic_down = track_going_down;
            _mom_cosmic_tree->Fill();
          }
        }
      }

      // 
      // Look at residuals
      //
      std::vector<TVector3> hit_v; // a vec of hits from coll plane
      std::vector<TVector3> track_v; // a vec of hits from coll plane

      // Collect hits
      auto iter = trackToHitsMap.find(candidate_track);
      if (iter != trackToHitsMap.end()) {
        std::vector<art::Ptr<recob::Hit>> hits = iter->second;
        for (auto hit : hits) {
          if (hit->View() == 2) {
            TVector3 h (hit->WireID().Wire, hit->PeakTime(), 0);
            hit_v.emplace_back(h);
          }
        }
      }

      // Collect track points
      for (size_t i = 0; i < candidate_track->NumberTrajectoryPoints(); i++) {
        try {
          if (candidate_track->HasValidPoint(i)) {
            TVector3 trk_pt = candidate_track->LocationAtPoint<TVector3>(i);
            double wire     = geo->NearestWire(trk_pt, 2);
            double time     = _detector_properties->ConvertXToTicks(trk_pt.X(), geo::PlaneID(0,0,2));
            TVector3 p (wire, time, 0.);
            track_v.emplace_back(p);
          }
        } catch (...) {
          continue;
        }
      }

      //if (_debug) std::cout << "[UBXSec] \t \t Hit points: " << hit_v.size() << ", track points: " << track_v.size() << std::endl;
      ubana::TrackQuality _track_quality;
      _track_quality.SetTrackPoints(track_v);
      _track_quality.SetHitCollection(hit_v);
      std::pair<double, double> residual_mean_std           = _track_quality.GetResiduals();
      std::pair<double, double> residual_truncated_mean_std = _track_quality.GetTruncatedResiduals();
      std::pair<double,int> dist_wire_pair                  = _track_quality.GetTrackGap();

      int start_wire = dist_wire_pair.second;
      int end_wire  = dist_wire_pair.second + dist_wire_pair.first;

      if (start_wire > end_wire)
        std::swap(start_wire, end_wire);

      // Create a channel to wire map
      std::map<int, int> wire_to_channel;
      for (unsigned int ch = 0; ch < 8256; ch++) {
        std::vector< geo::WireID > wire_v = geo->ChannelToWire(ch);
        wire_to_channel[wire_v[0].Wire] = ch;
      }

      int n_dead_wires = 0;

      for (int wire = start_wire; wire < end_wire; wire++) {

        int channel = wire_to_channel[wire];

        // Channel statuses: 1=dead, 3=noisy, 4=good
        if (chanFilt.Status(channel) < 4) {
          n_dead_wires++;
        }
      }

      double r = _track_quality.GetR();

      double n_hits_in_cluster = 0;
      auto it = recoParticlesToHits.find(candidate_pfp);
      if (it != recoParticlesToHits.end()) {
        for (auto h : it->second) {
          if (h->View() == 2) {
            n_hits_in_cluster++;
          }
        }
      }

      double ratio = (double)hit_v.size()/n_hits_in_cluster;

      // Also look at the scattering angle
      std::vector<TVector3> dir_v;
      for (size_t p = 0; p < candidate_track->NumberTrajectoryPoints(); p++) {
        if (!candidate_track->HasValidPoint(p)) continue;
        dir_v.push_back(candidate_track->DirectionAtPoint<TVector3>(p));
      }
      std::vector<double> angle_v;
      for (size_t p = 0; p < dir_v.size()-1; p++) {
        double angle = dir_v.at(p).Angle(dir_v.at(p+1));
        angle_v.push_back(angle);
      }
      std::sort(angle_v.begin(), angle_v.end());

      double max_angle = -1;
      if (angle_v.size() != 0)
        max_angle = angle_v.at(angle_v.size()-1) / TMath::Pi() * 180.;

      ubxsec_event->slc_muoncandidate_residuals_mean[slice]                 = residual_mean_std.first;
      ubxsec_event->slc_muoncandidate_residuals_std[slice]                  = residual_mean_std.second;
      ubxsec_event->slc_muoncandidate_residuals_truncatedmean[slice]        = residual_truncated_mean_std.first;
      ubxsec_event->slc_muoncandidate_residuals_truncatedstd[slice]         = residual_truncated_mean_std.second;
      ubxsec_event->slc_muoncandidate_wiregap[slice]                        = end_wire-start_wire;
      ubxsec_event->slc_muoncandidate_wiregap_dead[slice]                   = n_dead_wires;
      ubxsec_event->slc_muoncandidate_linearity[slice]                      = r;
      ubxsec_event->slc_muoncandidate_perc_used_hits_in_cluster[slice]      = ratio;
      ubxsec_event->slc_muoncandidate_maxscatteringangle[slice]             = max_angle;

      muon_candidate_track_per_slice_v.at(slice)                            = candidate_track;
      muon_candidate_pfparticle_per_slice_v.at(slice)                       = candidate_pfp;


    } else {

      ubxsec_event->slc_muoncandidate_exists[slice]    = false;
      ubxsec_event->slc_muoncandidate_length[slice]    = -9999;
      ubxsec_event->slc_muoncandidate_phi[slice]       = -9999;
      ubxsec_event->slc_muoncandidate_theta[slice]     = -9999;
      ubxsec_event->slc_muoncandidate_theta_xz[slice]  = -9999;
      ubxsec_event->slc_muoncandidate_theta_yz[slice]  = -9999;
      ubxsec_event->slc_muoncandidate_mom_range[slice] = -9999;
      ubxsec_event->slc_muoncandidate_mom_mcs[slice]   = -9999;

    }

    // Particle ID
    for (auto pfp : pfps_from_tpcobj){

      std::vector<art::Ptr<ubana::MCGhost>> mcghosts = mcghost_from_pfp.at(pfp.key());
      std::vector<art::Ptr<simb::MCParticle>> mcpars;
      int pdg = -1;
      if (mcghosts.size() == 0 || mcghosts.size() > 1 ) {
        continue;
      }

      mcpars = mcpar_from_mcghost.at(mcghosts[0].key());
      pdg = mcpars[0]->PdgCode();
      const auto mc_truth = UBXSecHelper::TrackIDToMCTruth(e, "largeant", mcpars[0]->TrackId());
      if (!mc_truth) {
        std::cerr << "[UBXSec] Problem with MCTruth pointer." << std::endl;
        continue;
      }
      if (mc_truth->Origin() == simb::kBeamNeutrino &&
          mcpars[0]->PdgCode() == 13 && mcpars[0]->Mother() == 0) {
        ubxsec_event->true_muon_mom_matched = mcpars[0]->P();

        if (_fiducial_volume.InFV(mcpars[0]->Vx(), mcpars[0]->Vy(), mcpars[0]->Vz()) &&
          _fiducial_volume.InFV(mcpars[0]->EndX(), mcpars[0]->EndY(), mcpars[0]->EndZ())) {
          ubxsec_event->mc_muon_contained = true;
        }
      }

      std::vector<art::Ptr<recob::Track>> tracks = tracks_from_pfp.at(pfp.key());
      for (auto track : tracks) {

        std::vector<art::Ptr<anab::ParticleID>> pids = particleids_from_track.at(track.key());
        if(pids.size() == 0){
	  continue;
	}

	std::vector<anab::sParticleIDAlgScores> AlgScoresVec = pids[0]->ParticleIDAlgScores();

	// Loop though AlgScoresVec and find the variables we want
	double tmppida=-9999;
	for (size_t i_algscore=0; i_algscore<AlgScoresVec.size(); i_algscore++){
	  anab::sParticleIDAlgScores AlgScore = AlgScoresVec.at(i_algscore);
	  int planenum = UBPID::uB_getSinglePlane(AlgScore.fPlaneMask);
	  if (planenum<0 || planenum>2) continue;
	  if (AlgScore.fVariableType==anab::kPIDA){
	    if ( planenum == 2 ) {
	      tmppida=AlgScore.fValue;
	    }
	  }
	}

	// Don't fill things if PIDA was not found
	if (tmppida==-9999) continue;

	if (pdg == 13) {
	  _h_pida_muon->Fill(tmppida);
	  _h_pida_len_muon->Fill(tmppida, track->Length());
	  if( tmppida > 0 && tmppida < 50. && _make_pida_csv) _csvfile << tmppida << "," << track->Length() << "," << "1" << std::endl;
	} else if (pdg == 2212) {
	  _h_pida_proton->Fill(tmppida);
	  _h_pida_len_proton->Fill(tmppida, track->Length());
	  if( tmppida > 0 && tmppida < 50. && _make_pida_csv) _csvfile << tmppida << "," << track->Length() << "," << "0" << std::endl;
	} else if (pdg == 211) {
	  _h_pida_pion->Fill(tmppida);
	  _h_pida_len_pion->Fill(tmppida, track->Length());
	} else if (pdg == 321) {
	  _h_pida_kaon->Fill(tmppida);
	  _h_pida_len_kaon->Fill(tmppida, track->Length());
	}
      }
    }

    // MCS - Track direction, study cosmic direction
    // Take the muon candidate track for a TPCObject
    // with cosmic origin, and check the track direction
    // (up/down) from MCS result
    if (muon_cand_exists && tpcobj.GetOrigin() == ubana::kCosmicRay) {
      bool best_fwd = mcsfitresult_mu_v.at(candidate_track.key())->isBestFwd();
      bool down_track = candidate_track->Vertex().Y() > candidate_track->End().Y();
      if (down_track && best_fwd) {               // Track is reco going down and mcs agrees (true for cosmic)
        _h_mcs_cosmic_track_direction->Fill(0);
        _mcs_cosmic_track_direction = 0;
        _mcs_cosmic_track_downll    = mcsfitresult_mu_v.at(candidate_track.key())->fwdLogLikelihood();
        _mcs_cosmic_track_upll      = mcsfitresult_mu_v.at(candidate_track.key())->bwdLogLikelihood();
      } else if (!down_track && best_fwd) {       // Track is reco going up and mcs agrees
        _h_mcs_cosmic_track_direction->Fill(1);
        _mcs_cosmic_track_direction = 1;
        _mcs_cosmic_track_downll    = mcsfitresult_mu_v.at(candidate_track.key())->bwdLogLikelihood();
        _mcs_cosmic_track_upll      = mcsfitresult_mu_v.at(candidate_track.key())->fwdLogLikelihood();
      } else if (down_track && !best_fwd) {       // Track is reco going down and mcs disagrees
        _h_mcs_cosmic_track_direction->Fill(1);
        _mcs_cosmic_track_direction = 2;
        _mcs_cosmic_track_downll    = mcsfitresult_mu_v.at(candidate_track.key())->fwdLogLikelihood();
        _mcs_cosmic_track_upll      = mcsfitresult_mu_v.at(candidate_track.key())->bwdLogLikelihood();
      } else if (!down_track && !best_fwd) {      // Track is reco going up and mcs disagrees (true for cosmic)
        _h_mcs_cosmic_track_direction->Fill(0);
        _mcs_cosmic_track_direction = 3;
        _mcs_cosmic_track_downll    = mcsfitresult_mu_v.at(candidate_track.key())->bwdLogLikelihood();
        _mcs_cosmic_track_upll      = mcsfitresult_mu_v.at(candidate_track.key())->fwdLogLikelihood();
      }
      _mcs_cosmic_track_direction_tree->Fill();
    }

    std::cout << "[UBXSec] --- SLICE INFORMATION SAVED" << std::endl;
  } // slice loop

  if (_is_mc) {

    // SW Trigger
    art::Handle<raw::ubdaqSoftwareTriggerData> softwareTriggerHandle;
    e.getByLabel("swtrigger", softwareTriggerHandle);

    if (softwareTriggerHandle.isValid()){
      if (softwareTriggerHandle->getNumberOfAlgorithms() == 1) {
        std::vector<std::string> algoNames = softwareTriggerHandle->getListOfAlgorithms();
        ubxsec_event->is_swtriggered = (softwareTriggerHandle->passedAlgo(algoNames[0]) ? 1 : 0);
      }
    }

    // MC Flash
    ::art::Handle<std::vector<recob::OpFlash> > nuMcflash_h;
    e.getByLabel("NeutrinoMCFlash",nuMcflash_h);
    if( !nuMcflash_h.isValid() || nuMcflash_h->empty() ) {
      std::cerr << "Don't have neutrino MC flashes." << std::endl;
    } else {

      auto const& flash = (*nuMcflash_h)[0];
      ubxsec_event->numc_flash_spec.resize(geo->NOpDets());
      for (unsigned int i = 0; i < geo->NOpDets(); i++) {
        unsigned int opdet = geo->OpDetFromOpChannel(i);
        ubxsec_event->numc_flash_spec[opdet] = flash.PE(i);
      }
    }

    // MCFlash vs op activity
    bool opActivityInBeamSpill = false;
    // Check if there are recon beam flashed in the beam spill window
    for (auto reco_fls_time : ubxsec_event->beamfls_time) {
      if (reco_fls_time > _beam_spill_start && reco_fls_time < _beam_spill_end) {
         opActivityInBeamSpill = true;
       }
    }
    if (nuMcflash_h->size() == 0) {
      if(opActivityInBeamSpill) {
        //std::cout << "No MCFlash but optical activity in the beam spill." << std::endl;
        ubxsec_event->no_mcflash_but_op_activity = true;
      }
    }

  }


  // POT
  art::Handle< sumdata::POTSummary > potsum_h;
  if(e.getByLabel(_potsum_producer, potsum_h))
    ubxsec_event->pot = potsum_h->totpot;
  else
    ubxsec_event->pot = 0.;

  // *********************
  // Event Selection
  // *********************

  // Use the top event loop to look at cuts that do not require more work (basic CRT, flash, track cuts).
  _event_selection.SetEvent(ubxsec_event);

  size_t slice_index=-999;
  int slice_index_1=slice_index;
  std::string reason = "no_failure";
  std::map<std::string,bool> failure_map;

  failure_map.clear();

  failure_map.emplace(std::make_pair("a_crt_veto_cut", false));
  failure_map.emplace(std::make_pair("b_beam_disc_flashes", false));
  failure_map.emplace(std::make_pair("c_beam_spill_flash", false));
  failure_map.emplace(std::make_pair("d_has_slices", false));
  failure_map.emplace(std::make_pair("e_has_slice_tagged_as_neutrino", false));
  failure_map.emplace(std::make_pair("f_ntrack", false));
  failure_map.emplace(std::make_pair("g_track_length", false));
  failure_map.emplace(std::make_pair("h_residuals_std_up", false));
  failure_map.emplace(std::make_pair("i_perc_used_hits_in_cluster", false));

  bool is_selected = _event_selection.IsSelected(slice_index_1, failure_map);
  if (_debug) std::cout << "[UBXSec] >>>>>>>>>>>>>>>>>>>>>> Is Selected? " << (is_selected ? "YES" : "NO") << std::endl;

  for (auto iter : failure_map) {

    std::cout << "[UBXSec] Cut: " << iter.first << "  >>>  " << (iter.second ? "PASSED" : "NOT PASSED") << std::endl;
    if ( !iter.second ) {
      reason = "fail_" + iter.first;
      break;
    }

  }

  std::cout << "[UBXSec] Selection Failed at Cut: " << reason << std::endl;

  ::ubana::SelectionResult selection_result;
  selection_result.SetSelectionType("numu_cc_inclusive");
  selection_result.SetFailureReason(reason);
  selection_result.SetCutFlowStatus(failure_map);

  // See if one of the last three cases fails.
  failure_map.emplace(std::make_pair("j_containment_cut", true));
  failure_map.emplace(std::make_pair("k_fiducial_volume_cut", true));
  failure_map.emplace(std::make_pair("l_length_of_tracks_from_vertex_cut", true));
  failure_map.emplace(std::make_pair("m_new_two_track_length_requirement", true));
  failure_map.emplace(std::make_pair("n_number_of_tracks_in_TPCObject_requirement", true));
  failure_map.emplace(std::make_pair("o_number_of_tracks_originating_from_vertex_requirement", true));
  failure_map.emplace(std::make_pair("p_CRT_distance_cut", true));

  // Increment the number of KDAR events if it is one.
  if ( kdar_from_dump_event == 1 ) 
    num_kdar_events_in_this_sample++;

  // *********************
  // Save Event Selection Output in the Event
  // *********************
  if ( !is_selected && is_duplicate_event == false ) {

    num_event_we_are_on++;
    
    total_num_events++;

    total_num_events_failing++;

    selection_result.SetSelectionStatus(false);

    ubxsec_event->is_selected = false;

    selectionResultVector->emplace_back(std::move(selection_result));

    // Go through the reasons that do not require additional work to find out if the event fails them.
    if ( reason == "fail_a_crt_veto_cut" ) {
      total_num_events_failed_crt_veto++;
      std::cout << "This event failed the CRT veto." << std::endl;
    }
                                                                                                                              
    if ( reason == "fail_b_beam_disc_flashes" ) {
      total_num_events_failed_beam_disc_flashes++;
      std::cout << "No beam discriminated flashes." << std::endl;
    }

    if ( reason == "fail_c_beam_spill_flash"  ) {
      total_num_events_failed_beam_spill_flash++;
      std::cout << "No beam discriminated flash > 50 PEs in the beamspill window." << std::endl;
    }

    if ( reason == "fail_d_has_slices"  ) {
      total_num_events_failed_has_slices++;
      std::cout << "There is no slice in this event." << std::endl;
    }

    if ( reason == "fail_f_ntrack"  ) {
      total_num_events_failed_ntrack++;
      std::cout << "There are less than the required number of tracks in the TPC Object." << std::endl;
    }

    if ( reason == "fail_g_track_length"  ) {
      total_num_events_failed_track_length++;
      std::cout << "The longest muon candidate track in the event fails the track length cut." << std::endl;
    }

    if ( reason == "fail_h_residuals_std_up"  ) {
      total_num_events_failed_residuals_std_up++;
      std::cout << "This event has failed residuals std up cut." << std::endl;
    }

    if ( reason == "fail_i_perc_used_hits_in_cluster"  ) {
      total_num_events_failed_perc_used_hits_in_cluster++;
      std::cout << "This event has failed the percentage of hits used in the cluster." << std::endl;
    }

  } else if ( is_selected && is_duplicate_event == false ) {

    total_num_events++;

    num_event_we_are_on++;

    // loop through the truth neutrinos themselves.                                                                                                                                                       
    num_of_neutrinos = 0;
    neutrino_energy  = -1.;

    for (auto& neutrino : NeutrinoVec ) {

      // Unpack the neutrino object to find an MCParticle.                                                                                                                                                
      const simb::MCNeutrino& truth_neutrino = neutrino->GetNeutrino();
      const simb::MCParticle& truth_particle = truth_neutrino.Nu();

      if ( ( truth_particle.E(0) * 1000.) > neutrino_energy ) 
	neutrino_energy = ( truth_particle.E(0) * 1000. );

      // Reset the coordinates of the neutrino as well from above.
      nu_vtx_x_truth               = truth_particle.Vx(0);
      nu_vtx_y_truth               = truth_particle.Vy(0);
      nu_vtx_z_truth               = truth_particle.Vz(0);

      if ( truth_neutrino.CCNC() == 1 ) 
	NC_channel = 1;

      num_of_neutrinos++;

    }

    // Save all of the information for the passing event.
    nslices_in_event = ubxsec_event->nslices;
    
    int nu_slc_idx = -1;

    // Find the slice that was tagged as a neutrino.
    for ( int slice_iter = 0; slice_iter < nslices_in_event; slice_iter++ ) {

      if ( ubxsec_event->slc_is_nu.at( slice_iter ) == true ) {
	nu_slc_idx = slice_iter;
	break;
      }

    }
 
    ubana::TPCObject tpcobj = (*tpcobj_h)[nu_slc_idx];

    auto tracks = tpcobjToTrackAssns.at( nu_slc_idx );

    // 'nu_slc_idx' necessarily has to be 0 or greater.
    std::vector<lar_pandora::TrackVector     > track_v_v_for_event_tree;
    track_v_v_for_event_tree.push_back(tpcobjToTrackAssns.at(nu_slc_idx));

    // Find the longest_track.
    recob::Track lt_track;
    UBXSecHelper::GetLongestTrackFromTPCObj(track_v_v_for_event_tree[nu_slc_idx], lt_track );

    size_t muon_candidate_idx = -1;
    
    for ( size_t pandora_track_iter = 0; pandora_track_iter < pandora_track_h->size(); pandora_track_iter++ ) {

      if ( fabs( pandora_track_h->at( pandora_track_iter ).Length() - ubxsec_event->slc_longesttrack_length[nu_slc_idx] ) < 0.001 )
	muon_candidate_idx = pandora_track_iter;

    } // End of the loop over the pandora tracks to find the one that corresponds to the muon candidate.

    // Start Neutrino Energy Calculation Code.                                                                                                                                                          
    bool first_points_of_track_are_start = true;

    // Set a variable for the calorimety object of the muon candidate on the track.                                                                                                                      
    auto Calo_v                                        = trk_calo_assn_v.at( muon_candidate_idx );
    auto Calo_v_no_SCE_corrections                     = trk_calo_assn_v_no_SCE_corrections.at( muon_candidate_idx );

    auto vhit                                          = fmthm.at( muon_candidate_idx );
    auto vmeta                                         = fmthm.data( muon_candidate_idx );

    int  count                                         = 0;

    // For each of the planes, first calculate the median values.
    median_dQdx_of_TPCObject_muon_candidate_u_plane    = 0.;
    median_dQdx_of_TPCObject_muon_candidate_v_plane    = 0.;
    median_dQdx_of_TPCObject_muon_candidate_y_plane    = 0.;

    truncated_dQdx_of_TPCObject_muon_candidate_u_plane = 0.;
    truncated_dQdx_of_TPCObject_muon_candidate_v_plane = 0.;
    truncated_dQdx_of_TPCObject_muon_candidate_y_plane = 0.;

    std::vector< float > xyz_v_x_coords;
    xyz_v_x_coords.clear();

    std::vector< float > xyz_v_y_coords;
    xyz_v_y_coords.clear();

    std::vector< float > xyz_v_z_coords;
    xyz_v_z_coords.clear();

    std::vector< float > xyz_v_no_SCE_corrections_x_coords;
    xyz_v_no_SCE_corrections_x_coords.clear();

    std::vector< float > xyz_v_no_SCE_corrections_y_coords;
    xyz_v_no_SCE_corrections_y_coords.clear();

    std::vector< float > xyz_v_no_SCE_corrections_z_coords;
    xyz_v_no_SCE_corrections_z_coords.clear();

    // These 'double's were previously 'float's.
    std::vector < float > dqdx_v;
    dqdx_v.clear();

    std::vector < float > dedx_v;
    dedx_v.clear();
    
    std::vector < float > dqdx_points_ordered;
    dqdx_points_ordered.clear();

    std::vector < float > dedx_points_ordered;
    dedx_points_ordered.clear();

    std::vector< float > dQdx_outliers_and_endpoints_removed_one_standard_deviation;
    dQdx_outliers_and_endpoints_removed_one_standard_deviation.clear();

    std::vector< float > dEdx_outliers_and_endpoints_removed;
    dEdx_outliers_and_endpoints_removed.clear();

    truncated_dQdx_of_TPCObject_muon_candidate_u_plane = 0;
    truncated_dQdx_of_TPCObject_muon_candidate_v_plane = 0;
    truncated_dQdx_of_TPCObject_muon_candidate_y_plane = 0;
    median_dQdx_of_TPCObject_muon_candidate_u_plane    = 0;
    median_dQdx_of_TPCObject_muon_candidate_v_plane    = 0;
    median_dQdx_of_TPCObject_muon_candidate_y_plane    = 0;
    
    bool fill_vector                                   = true;

    // Loop through the planes to orient the muon track.
    for ( size_t pl = 0; pl < 3; pl++ ) {

      std::cout << "Plane #" << pl << "." << std::endl;

      // Do not loop over planes with a meaningless index.
      if ( Calo_v[pl]->PlaneID().Plane != 0 && Calo_v[pl]->PlaneID().Plane != 1 && Calo_v[pl]->PlaneID().Plane != 2 ) 
	continue;

      dqdx_v.clear();
      dedx_v.clear();
      dqdx_points_ordered.clear();
      dedx_points_ordered.clear();
      dQdx_outliers_and_endpoints_removed_one_standard_deviation.clear();
      dEdx_outliers_and_endpoints_removed.clear();
      xyz_v_x_coords.clear();
      xyz_v_no_SCE_corrections_x_coords.clear();
      xyz_v_y_coords.clear();
      xyz_v_no_SCE_corrections_y_coords.clear();
      xyz_v_z_coords.clear();
      xyz_v_no_SCE_corrections_z_coords.clear();

      if ( Calo_v[pl]->dQdx().size() > 20 ) {

        if ( Calo_v[pl]->PlaneID().Plane == 0 )
	  u_plane_has_more_than_20_points = true;
      
	if ( Calo_v[pl]->PlaneID().Plane == 1 )
	  v_plane_has_more_than_20_points = true;

	if ( Calo_v[pl]->PlaneID().Plane == 2 )
	  y_plane_has_more_than_20_points = true;

      }
      
      // Loop through the vector and fill the relevant vectors.                                                                                                                                       
      for ( size_t hit_iter = 0; hit_iter < vmeta.size(); hit_iter++ ) {

	fill_vector = true;

	int ind = vmeta[hit_iter]->Index();

	// check that the traj point is in the calorimetry point                                                                     
	// and belongs to the plane we are interested in                                                                             
	if( pandora_track_h->at( muon_candidate_idx ).HasValidPoint(ind) && vhit[hit_iter]->WireID().Plane == pl && count < int( Calo_v[pl]->dQdx().size() ) ){

	  auto const& hit = vhit.at( hit_iter );

	  // Check to make sure the hit is not in one of the spikes.
	  float rms = hit->RMS();

	  // Filter on the peaks (we don't want to use these in those quantities).                                                  
	  if ( fabs( rms - 0.333 ) < 0.0001 || fabs( rms - 0.5 ) < 0.0001 || fabs( rms - 0.667 ) < 0.0001 || fabs( rms - 1.0 ) < 0.0001 || fabs( rms - 1.333 ) < 0.0001 || fabs( rms - 1.667 ) < 0.0001 || fabs( rms - 2.0 ) < 0.0001 || fabs( rms - 2.333 ) < 0.0001 || fabs( rms - 2.667 ) < 0.0001 || fabs( rms - 2.8572 ) < 0.0001 || fabs( rms - 3.0 ) < 0.0001 || fabs( rms - 3.333 ) < 0.0001 || fabs( rms - 3.667 ) < 0.0001 || fabs( rms - 4.0 ) < 0.0001 || fabs( rms - 4.333 ) < 0.0001 || fabs( rms - 4.667 ) < 0.0001 || fabs( rms - 5.0 ) < 0.0001 || fabs( rms - 5.333 ) < 0.0001 || fabs( rms - 5.667 ) < 0.0001 || fabs( rms - 6.0 ) < 0.0001 || fabs( rms - 6.333 ) < 0.0001 || fabs( rms - 6.5 ) < 0.0001 || fabs( rms - 7.0 ) < 0.0001 || fabs( rms - 7.333 ) < 0.0001 || fabs( rms - 7.5 ) < 0.0001 || fabs( rms - 7.667 ) < 0.0001 || fabs( rms - 8.0 ) < 0.0001 || fabs( rms - 8.5 ) < 0.0001 || fabs( rms - 9.0 ) < 0.0001 || fabs( rms - 9.5 ) < 0.0001 || fabs( rms - 10.0 ) < 0.0001 || fabs( rms - 10.5 ) < 0.0001 || fabs( rms - 11.0 ) < 0.0001 || fabs( rms - 11.5 ) < 0.0001 || fabs( rms - 12.0 ) < 0.0001 || fabs( rms - 12.5 ) < 0.0001 || fabs( rms - 13.0 ) < 0.0001 || fabs( rms - 13.5 ) < 0.0001 || fabs( rms - 14.0 ) < 0.0001 || fabs( rms - 14.5 ) < 0.0001 || fabs( rms - 15.0 ) < 0.0001 || fabs( rms - 15.5 ) < 0.0001 || fabs( rms - 16.0 ) < 0.0001 || fabs( rms - 16.5 ) < 0.0001 || fabs( rms - 17.0 ) < 0.0001 || fabs( rms - 17.5 ) < 0.0001 || fabs( rms - 18.0 ) < 0.0001 || fabs( rms - 18.5 ) < 0.0001 || fabs( rms - 19.0 ) < 0.0001 || fabs( rms - 19.5 ) < 0.0001 || fabs( rms - 20.0 ) < 0.0001 )
	    fill_vector = false;
	
	  if ( fill_vector == true ) {
	    dqdx_v.push_back( Calo_v[pl]->dQdx()[count] );
	    dedx_v.push_back( Calo_v[pl]->dEdx()[count] );
	    dqdx_points_ordered.push_back( Calo_v[pl]->dQdx()[count] );
	    dedx_points_ordered.push_back( Calo_v[pl]->dEdx()[count] );
	    xyz_v_x_coords.push_back( Calo_v[pl]->XYZ()[count].X() );
	    xyz_v_no_SCE_corrections_x_coords.push_back( Calo_v_no_SCE_corrections[pl]->XYZ()[count].X() );
	    xyz_v_y_coords.push_back( Calo_v[pl]->XYZ()[count].Y() );
	    xyz_v_no_SCE_corrections_y_coords.push_back( Calo_v_no_SCE_corrections[pl]->XYZ()[count].Y() );
	    xyz_v_z_coords.push_back( Calo_v[pl]->XYZ()[count].Z() );
	    xyz_v_no_SCE_corrections_z_coords.push_back( Calo_v_no_SCE_corrections[pl]->XYZ()[count].Z() );
	  }

	  count++;

	} // End of the conditional that the point actually exists along the track.

      } // End of the loop over the hits.

      if ( dedx_points_ordered.size() > 1 ) {

	// Now sort the vector (put it in order).
	std::sort( dqdx_points_ordered.begin(), dqdx_points_ordered.end() ); 
	std::sort( dedx_points_ordered.begin(), dedx_points_ordered.end() );
      
	if ( Calo_v[pl]->PlaneID().Plane == 0 )                                                                                                                                                           
	  median_dQdx_of_TPCObject_muon_candidate_u_plane = dqdx_points_ordered.at( int( dqdx_points_ordered.size() / 2 ) );                                                                              
	
	if ( Calo_v[pl]->PlaneID().Plane == 1 )                                                                                                                                                           
	  median_dQdx_of_TPCObject_muon_candidate_v_plane = dqdx_points_ordered.at( int( dqdx_points_ordered.size() / 2 ) );                                                                              
	
	if ( Calo_v[pl]->PlaneID().Plane == 2 )                                                                                                                                                           
	  median_dQdx_of_TPCObject_muon_candidate_y_plane = dqdx_points_ordered.at( int( dqdx_points_ordered.size() / 2 ) );                                                                              
	
	double dedx_median_value = dedx_points_ordered.at( int( dedx_points_ordered.size() / 2 ) );                                                                                                      
	
	// Find the standard deviation on this median.                                                                                                                                                    
	double rms_sum            = 0.;                                                                                                                                                                
	double difference_squared = 0.;                                                                                                                                                                  
	
	for ( size_t point_iter = 0; point_iter < dedx_v.size(); point_iter++ ) {                                                                                                                         
          
	  difference_squared = TMath::Power( ( dedx_median_value - dedx_v.at( point_iter ) ), 2);                                                                                                          
	  
	  rms_sum += difference_squared;                                                                                                                                                                   
                                                                                                                                                                                                          
	}                                                                                                                                                                                                  
                                                                                                                                                                                                         
	// Calculate the rms value.                                                                                                                                                                        
	double rms_value = TMath::Sqrt( ( rms_sum ) / double( dedx_v.size() ) );                                                                                                                           
      
	// Remove the points that are greater than two standard deviations from the median.                                                                                                                
	// New dE/dx vector.                                                                                                                                                                            
	std::vector< double > dEdx_outliers_and_endpoints_removed;                                                                                                                                         
	dEdx_outliers_and_endpoints_removed.clear();           
	
	for ( size_t energy_point_iter = 0; energy_point_iter < dedx_v.size(); energy_point_iter++ ) {        
	  
	  if ( fabs( dedx_median_value - dedx_v.at( energy_point_iter ) ) < ( 2.0 * rms_value ) )                                                                                                          
	    dEdx_outliers_and_endpoints_removed.push_back( dedx_v.at( energy_point_iter ) );                                                                                                              
	
	}
      
	// Do the same for the charge points.
	double dqdx_median_value = dqdx_points_ordered.at( int( dqdx_points_ordered.size() / 2 ) );
	
	// Find the standard deviation on this median.                                                                                                                                                     
	double charge_rms_sum            = 0.;
	double charge_difference_squared = 0.;
	
	for ( size_t point_iter = 0; point_iter < dqdx_v.size(); point_iter++ ) {
	
	  charge_difference_squared = TMath::Power( ( dqdx_median_value - dqdx_v.at( point_iter ) ), 2);
	  
	  charge_rms_sum            += charge_difference_squared;

	}
	
	// Calculate the rms value.                                                                                                                                                                      
	double charge_rms_value = TMath::Sqrt( ( charge_rms_sum ) / double( dqdx_v.size() ) );
	
	for ( size_t charge_point_iter = 0; charge_point_iter < dqdx_v.size(); charge_point_iter++ ) {
	  
	  if ( fabs( dqdx_median_value - dqdx_v.at( charge_point_iter ) ) < ( charge_rms_value ) )
	    dQdx_outliers_and_endpoints_removed_one_standard_deviation.push_back( dqdx_v.at( charge_point_iter ) );
	  
	}
      
	// Calculate the value of the truncated mean standard deviation on these points.
	float truncated_mean_standard_deviation = 0.;
	
	int num_of_elements                     = dQdx_outliers_and_endpoints_removed_one_standard_deviation.size();
	
	for ( int mean_iter = 0; mean_iter < num_of_elements; mean_iter++ ) {
	
	  truncated_mean_standard_deviation += dQdx_outliers_and_endpoints_removed_one_standard_deviation.at( mean_iter );
	
	}
      
	if ( Calo_v[pl]->PlaneID().Plane == 0 ) 
	  truncated_dQdx_of_TPCObject_muon_candidate_u_plane = ( truncated_mean_standard_deviation / float( num_of_elements ) );
	
	if ( Calo_v[pl]->PlaneID().Plane == 1 )
	  truncated_dQdx_of_TPCObject_muon_candidate_v_plane = ( truncated_mean_standard_deviation / float( num_of_elements ) );
	
	if ( Calo_v[pl]->PlaneID().Plane == 2 ) 
	  truncated_dQdx_of_TPCObject_muon_candidate_y_plane = ( truncated_mean_standard_deviation / float( num_of_elements ) );

	// Find the direction of the track on this plane.
	// Fit this vector to a histogram in both directions.                                                                                                                                           
	TH1D* track_forward_hist    = new TH1D("track_forward_hist", "dE/dx vs. Point on Track", dEdx_outliers_and_endpoints_removed.size(), 0, dEdx_outliers_and_endpoints_removed.size() );          
	
	TH1D* track_backward_hist   = new TH1D("track_backward_hist", "dE/dx vs. Point on Track", dEdx_outliers_and_endpoints_removed.size(), 0, dEdx_outliers_and_endpoints_removed.size() );           
	
	// Loop through the points and fill both of the histograms.                                                                                                                                     
	for ( size_t histogram_iter = 0; histogram_iter < dEdx_outliers_and_endpoints_removed.size(); histogram_iter++ ) {                                                                               
	  
	  track_forward_hist->SetBinContent( histogram_iter, dEdx_outliers_and_endpoints_removed.at( histogram_iter ) );                                                                              
	  track_backward_hist->SetBinContent( ( dEdx_outliers_and_endpoints_removed.size() - 1 - histogram_iter ), dEdx_outliers_and_endpoints_removed.at( histogram_iter ) );                         
        
	}                                                                                                                                                                                               
                                                                       	
	double forward_slope  = 0.;                                                                                                                                                                    
	double backward_slope = 0.;                 
	
	TF1 *st_y_function  = new TF1("st_y_function","[0]*x + [1]", 0, dEdx_outliers_and_endpoints_removed.size() );
	
	track_forward_hist->Fit(st_y_function,"R+");                                                                                                                                                     
	forward_slope = st_y_function->GetParameter(0);                                                                                                                                                 
	  
	track_backward_hist->Fit(st_y_function,"R+");                                                                                                                                                 
	backward_slope = st_y_function->GetParameter(0);                                                                                                                                             

	if ( backward_slope < forward_slope ) {
	
	  if ( Calo_v[pl]->PlaneID().Plane == 0 )
	    u_plane_has_track_forward = true;
	  
	  if ( Calo_v[pl]->PlaneID().Plane == 1 )
	    v_plane_has_track_forward = true;
	
	  if ( Calo_v[pl]->PlaneID().Plane == 2 )
	    y_plane_has_track_forward = true;
	
	}

      }
	
      count = 0;

    } // End of the loop over the planes.

    // Go through the nine cases of track orientation.                                                                                                                                                 
    if      ( u_plane_has_track_forward && u_plane_has_more_than_20_points && v_plane_has_track_forward && v_plane_has_more_than_20_points )
      first_points_of_track_are_start = true;

    else if (  u_plane_has_track_forward && u_plane_has_more_than_20_points && y_plane_has_track_forward && y_plane_has_more_than_20_points )
      first_points_of_track_are_start = true;

    else if ( v_plane_has_track_forward && v_plane_has_more_than_20_points && y_plane_has_track_forward && y_plane_has_more_than_20_points )
      first_points_of_track_are_start = true;

    else if ( !u_plane_has_more_than_20_points && !v_plane_has_more_than_20_points && y_plane_has_more_than_20_points && y_plane_has_track_forward )
      first_points_of_track_are_start = true;

    else if ( !u_plane_has_more_than_20_points && !y_plane_has_more_than_20_points && v_plane_has_more_than_20_points && v_plane_has_track_forward )
      first_points_of_track_are_start = true;

    else if ( !v_plane_has_more_than_20_points && !y_plane_has_more_than_20_points && u_plane_has_more_than_20_points && u_plane_has_track_forward )
      first_points_of_track_are_start = true;

    else if ( !u_plane_has_more_than_20_points && !v_plane_has_more_than_20_points && !y_plane_has_more_than_20_points )
      first_points_of_track_are_start = true;

    else
      first_points_of_track_are_start = false;

    // Orient the track based on the verdict of track orientation.
    if ( xyz_v_x_coords.size() > 0 ) {

      if ( first_points_of_track_are_start == true ) {

	vertex_location_x          = xyz_v_x_coords.at( 0 ) - ( flash_time * 0.111436 );                                                                                                                  
	vertex_location_y          = xyz_v_y_coords.at( 0 );                                                                                                                                             
	vertex_location_z          = xyz_v_z_coords.at( 0 );                                                                                                                                             

	other_end_location_x       = xyz_v_x_coords.at( xyz_v_x_coords.size() - 1 ) - ( flash_time * 0.111436 );                                                                                         
	other_end_location_y       = xyz_v_y_coords.at( xyz_v_y_coords.size() - 1 );                                                                                                                      
	other_end_location_z       = xyz_v_z_coords.at( xyz_v_z_coords.size() - 1 );                                                                                                                       
                                                                                                                                                                                                           
	vertex_location_x_with_SCE = xyz_v_no_SCE_corrections_x_coords.at( 0 ) - ( flash_time * 0.111436 );                                                                                              
	vertex_location_y_with_SCE = xyz_v_no_SCE_corrections_y_coords.at( 0 );                                                                                                                            
	vertex_location_z_with_SCE = xyz_v_no_SCE_corrections_z_coords.at( 0 ); 

      }

      else {

	vertex_location_x          = xyz_v_x_coords.at( xyz_v_x_coords.size() - 1 ) - ( flash_time * 0.111436 );                                                                                        
	vertex_location_y          = xyz_v_y_coords.at( xyz_v_y_coords.size() - 1 );                                                                                                                    
	vertex_location_z          = xyz_v_z_coords.at( xyz_v_z_coords.size() - 1 );                                                                                                                            
      
	other_end_location_x       = xyz_v_x_coords.at( 0 ) - ( flash_time * 0.111436 );                                                                                                                 
	other_end_location_y       = xyz_v_y_coords.at( 0 );                                                                                                                                              
	other_end_location_z       = xyz_v_z_coords.at( 0 );                                                                                                                                                     
                                                                                                                                                                                                          
	vertex_location_x_with_SCE = xyz_v_no_SCE_corrections_x_coords.at( xyz_v_no_SCE_corrections_x_coords.size() - 1 ) - ( flash_time * 0.111436 );                                                    
	vertex_location_y_with_SCE = xyz_v_no_SCE_corrections_y_coords.at( xyz_v_no_SCE_corrections_y_coords.size() - 1 );                                                                                
	vertex_location_z_with_SCE = xyz_v_no_SCE_corrections_z_coords.at( xyz_v_no_SCE_corrections_z_coords.size() - 1 );
	
      }

    }

    else {

      auto const& trktraj = pandora_track_h->at( muon_candidate_idx ).Trajectory();

      auto firstValid     = trktraj.FirstValidPoint();
      auto lastValid      = trktraj.LastValidPoint();

      if ( first_points_of_track_are_start == true ) {

	vertex_location_x          = pandora_track_h->at( muon_candidate_idx ).LocationAtPoint( firstValid ).X() - ( 0.111436 * flash_time );
        vertex_location_y          = pandora_track_h->at( muon_candidate_idx ).LocationAtPoint( firstValid ).Y();
	vertex_location_z          = pandora_track_h->at( muon_candidate_idx ).LocationAtPoint( firstValid ).Z();

        other_end_location_x       = pandora_track_h->at( muon_candidate_idx ).LocationAtPoint( lastValid ).X() - ( 0.111436 * flash_time );
        other_end_location_y       = pandora_track_h->at( muon_candidate_idx ).LocationAtPoint( lastValid ).Y();
        other_end_location_z       = pandora_track_h->at( muon_candidate_idx ).LocationAtPoint( lastValid ).Z();

	// Include the SCE in this case.
        vertex_location_x_with_SCE = pandora_track_h->at( muon_candidate_idx ).LocationAtPoint( firstValid ).X() - ( 0.111436 * flash_time );
	vertex_location_y_with_SCE = pandora_track_h->at( muon_candidate_idx ).LocationAtPoint( firstValid ).Y();
        vertex_location_z_with_SCE = pandora_track_h->at( muon_candidate_idx ).LocationAtPoint( firstValid ).Z();

      }

      else {

	vertex_location_x          = pandora_track_h->at( muon_candidate_idx ).LocationAtPoint( lastValid ).X() - ( 0.111436 * flash_time );
	vertex_location_y          = pandora_track_h->at( muon_candidate_idx ).LocationAtPoint( lastValid ).Y();
        vertex_location_z          = pandora_track_h->at( muon_candidate_idx ).LocationAtPoint( lastValid ).Z();

        other_end_location_x       = pandora_track_h->at( muon_candidate_idx ).LocationAtPoint( firstValid ).X() - ( 0.111436 * flash_time );
        other_end_location_y       = pandora_track_h->at( muon_candidate_idx ).LocationAtPoint( firstValid ).Y();
	other_end_location_z       = pandora_track_h->at( muon_candidate_idx ).LocationAtPoint( firstValid ).Z();

	// Include the SCE in this case.
        vertex_location_x_with_SCE = pandora_track_h->at( muon_candidate_idx ).LocationAtPoint( lastValid ).X() - ( 0.111436 * flash_time );
        vertex_location_y_with_SCE = pandora_track_h->at( muon_candidate_idx ).LocationAtPoint( lastValid ).Y();
        vertex_location_z_with_SCE = pandora_track_h->at( muon_candidate_idx ).LocationAtPoint( lastValid ).Z();

      }

    }

      muon_candidate_length         = ubxsec_event->slc_longesttrack_length[nu_slc_idx];
      reconstructed_muon_KE         = ( TMath::Sqrt( ( p_calculator_from_length.GetTrackMomentum( muon_candidate_length, 13 ) * 1000. ) * ( p_calculator_from_length.GetTrackMomentum( muon_candidate_length, 13 ) * 1000. ) + 105.7 * 105.7 ) - 105.7 );
      muon_candidate_kinetic_energy = ( TMath::Sqrt( ( p_calculator_from_length.GetTrackMomentum( muon_candidate_length, 13 ) * 1000. ) * ( p_calculator_from_length.GetTrackMomentum( muon_candidate_length, 13 ) * 1000. ) + 105.7 * 105.7 ) - 105.7 );

      // Place a requirement on the maximum length of the proton based on the muon length.                                                                                                              
      // Use 250 MeV to give a little leeway.                                                                                                                                                          
      // 105.7 MeV = Muon Mass, 40 MeV = binding energy.
      maximum_proton_kinetic_energy = ( 250 - reconstructed_muon_KE - 105.7 - 40. );
      maximum_proton_length         = 0.0;

      // Convert this number into a length based on piecewise information.                                                                                                                             
      // The maximum length is the length at the higher value of the kinetic energy.                                                                                                                   
      // The extra allowance given in 'maximum_proton_kinetic_energy' makes this ok.                                                                                                                      
      if ( maximum_proton_kinetic_energy < 30.0 )
	maximum_proton_length = 1.0;

      else if ( maximum_proton_kinetic_energy > 30.0 && maximum_proton_kinetic_energy < 45.0 )
	maximum_proton_length = 2.0;

      else if ( maximum_proton_kinetic_energy > 45.0 && maximum_proton_kinetic_energy < 57.0 )
	maximum_proton_length = 3.0;

      else if ( maximum_proton_kinetic_energy > 57.0 && maximum_proton_kinetic_energy < 67.0 )
	maximum_proton_length = 4.0;

      else if ( maximum_proton_kinetic_energy > 67.0 && maximum_proton_kinetic_energy < 77.0 )
	maximum_proton_length = 5.0;

      else if ( maximum_proton_kinetic_energy > 77.0 && maximum_proton_kinetic_energy < 85.0 )
	maximum_proton_length = 6.0;

      else if ( maximum_proton_kinetic_energy > 85.0 && maximum_proton_kinetic_energy < 93.0 )
	maximum_proton_length = 7.0;

      else if ( maximum_proton_kinetic_energy > 93.0 && maximum_proton_kinetic_energy < 101.0 )
	maximum_proton_length = 8.0;

      else if ( maximum_proton_kinetic_energy > 101.0 && maximum_proton_kinetic_energy < 108.0 )
	maximum_proton_length = 9.0;

      else if ( maximum_proton_kinetic_energy > 108.0 && maximum_proton_kinetic_energy < 115.0 )
	maximum_proton_length = 10.0;

      else if ( maximum_proton_kinetic_energy > 115.0 && maximum_proton_kinetic_energy < 122.0 )
	maximum_proton_length = 11.0;

      wire_space_max_distance    = ( maximum_proton_length / 0.3 );
      tick_space_max_distance    = ( ( maximum_proton_length / 0.055718 ) + 300.0 ); // I add an extra space onto here.

      // This part is necessary for finding the high-charge pixels at the end of the muon track.

      // Identify the points on the muon track for orientating the track in the 3D fitting algorithm.
      auto const& trktraj = pandora_track_h->at( muon_candidate_idx ).Trajectory();

      auto firstValid     = trktraj.FirstValidPoint();
      auto lastValid      = trktraj.LastValidPoint();

      // Loop through the longest track and calculate the piecewise track distance.                                                                                                                     
      for ( size_t point_iter = firstValid; point_iter < ( pandora_track_h->at( muon_candidate_idx ).NumberTrajectoryPoints() - 1 ); point_iter++ )  {

	double point0_X = track_h->at( muon_candidate_idx ).LocationAtPoint( point_iter ).X();
	double point0_Y = track_h->at( muon_candidate_idx ).LocationAtPoint( point_iter ).Y();
	double point0_Z = track_h->at( muon_candidate_idx ).LocationAtPoint( point_iter ).Z();

	double point1_X = track_h->at( muon_candidate_idx ).LocationAtPoint( point_iter + 1 ).X();
	double point1_Y = track_h->at( muon_candidate_idx ).LocationAtPoint( point_iter + 1 ).Y();
	double point1_Z = track_h->at( muon_candidate_idx ).LocationAtPoint( point_iter + 1 ).Z();

	if ( point_iter == firstValid ) {

	  muon_track_first_point_x = point0_X - ( 0.111436 * flash_time );
	  muon_track_first_point_y = point0_Y;
	  muon_track_first_point_z = point0_Z;

	  // Find the offsets on each of the points.                                                                                                                                                   
	  geo::Vector_t muon_vertex_offsets = sce->GetCalPosOffsets(geo::Point_t{muon_track_first_point_x,muon_track_first_point_y,muon_track_first_point_z});

	  muon_track_first_point_x  =  muon_track_first_point_x - muon_vertex_offsets.X();
	  muon_track_first_point_y  =  muon_track_first_point_y + muon_vertex_offsets.Y();
	  muon_track_first_point_z   = muon_track_first_point_z + muon_vertex_offsets.Z();

	}

	if ( point_iter == ( lastValid - 1 ) ) {

	  muon_track_last_point_x = point1_X - ( 0.111436 * flash_time );
	  muon_track_last_point_y = point1_Y;
	  muon_track_last_point_z = point1_Z;

	  // Find the offsets on each of the points.                                                                                                                                                  
	  geo::Vector_t muon_vertex_offsets = sce->GetCalPosOffsets(geo::Point_t{muon_track_last_point_x,muon_track_last_point_y,muon_track_last_point_z});

	  muon_track_last_point_x  =  muon_track_last_point_x - muon_vertex_offsets.X();
	  muon_track_last_point_y  =  muon_track_last_point_y + muon_vertex_offsets.Y();
	  muon_track_last_point_z  =  muon_track_last_point_z + muon_vertex_offsets.Z();

	}

      } // End of the loop over the trajectory points on the muon candidate track.
      
      // First, add in the energy from the tracks that are located at the start of the muon track.
      proton_kinetic_energy_from_additional_tracks = 0;

      idx_of_additional_tracks_originating_from_vertex.clear();

      for ( int tpcobject_track_iter = 0; tpcobject_track_iter < int(pandora_track_h->size()); tpcobject_track_iter++ ) {

	// Skip the longest track (the vertex is found with respect to that one).                                                                                                                      
	if ( tpcobject_track_iter == int( muon_candidate_idx ) )
	  continue;

	// Find the starting and ending x, y, and z coordinates of the muon track and of the track that you are currently looking at.                                                                  
	auto track = pandora_track_h->at( tpcobject_track_iter );

	auto first_valid_point_idx = track.FirstValidPoint();
	auto last_valid_point_idx  = track.LastValidPoint();
	
	auto first_valid_point     = track.LocationAtPoint( first_valid_point_idx );
	auto last_valid_point      = track.LocationAtPoint( last_valid_point_idx );

	double other_track_x0      = first_valid_point.X() - ( 0.111436 * flash_time );
	double other_track_y0      = first_valid_point.Y();
	double other_track_z0      = first_valid_point.Z();
	double other_track_x1      = last_valid_point.X() - ( 0.111436 * flash_time );
	double other_track_y1      = last_valid_point.Y();
	double other_track_z1      = last_valid_point.Z();

	// Offset these track points for the space charge effect.                                                                                                                                    
	geo::Vector_t other_track_first_point_offsets = sce->GetCalPosOffsets(geo::Point_t{other_track_x0,other_track_y0,other_track_z0});
	
	double other_track_x0_with_offset = other_track_x0 - other_track_first_point_offsets.X();
	double other_track_y0_with_offset = other_track_y0 + other_track_first_point_offsets.Y();
	double other_track_z0_with_offset = other_track_z0 + other_track_first_point_offsets.Z();
	
	geo::Vector_t other_track_second_point_offsets = sce->GetCalPosOffsets(geo::Point_t{other_track_x1,other_track_y1,other_track_z1});
	
	double other_track_x1_with_offset = other_track_x1 - other_track_second_point_offsets.X();
	double other_track_y1_with_offset = other_track_y1 + other_track_second_point_offsets.Y();
	double other_track_z1_with_offset = other_track_z1 + other_track_second_point_offsets.Z();

	// Find which point is closer to both the vertex and the end of the muon candidate.                                                                                                           
	double closest_distance_vertex = 0.;
	double closest_distance_end    = 0.;
	
	double distance_vertex_point0  = TMath::Sqrt( ( vertex_location_x - other_track_x0_with_offset ) * ( vertex_location_x - other_track_x0_with_offset ) + ( vertex_location_y - other_track_y0_with_offset ) * ( vertex_location_y - other_track_y0_with_offset ) + ( vertex_location_z - other_track_z0_with_offset ) * ( vertex_location_z - other_track_z0_with_offset ) );
	
	double distance_vertex_point1  = TMath::Sqrt( ( vertex_location_x - other_track_x1_with_offset ) * ( vertex_location_x - other_track_x1_with_offset ) + ( vertex_location_y - other_track_y1_with_offset ) * ( vertex_location_y - other_track_y1_with_offset ) + ( vertex_location_z - other_track_z1_with_offset ) * ( vertex_location_z - other_track_z1_with_offset ) );
	
	double distance_end_point0     = TMath::Sqrt( ( other_end_location_x - other_track_x0_with_offset ) * ( other_end_location_x - other_track_x0_with_offset ) + ( other_end_location_y - other_track_y0_with_offset ) * ( other_end_location_y - other_track_y0_with_offset ) + ( other_end_location_z - other_track_z0_with_offset ) * ( other_end_location_z - other_track_z0_with_offset ) );
																	
	double distance_end_point1     = TMath::Sqrt( ( other_end_location_x - other_track_x1_with_offset ) * ( other_end_location_x - other_track_x1_with_offset ) + ( other_end_location_y - other_track_y1_with_offset ) * ( other_end_location_y - other_track_y1_with_offset ) + ( other_end_location_z - other_track_z1_with_offset ) * ( other_end_location_z - other_track_z1_with_offset ) );
	
	if ( distance_vertex_point0 < distance_vertex_point1 ) {
	  closest_distance_vertex = distance_vertex_point0;
	}

	else {
	  closest_distance_vertex = distance_vertex_point1;
	}
	
	if ( distance_end_point0 < distance_end_point1 ) {
	  closest_distance_end = distance_end_point0;
	}
	
	else {
	  closest_distance_end = distance_end_point1;
	}

	// Now, determine if the track originates at the vertex or at the other end.                                                                                                                  
	if ( closest_distance_vertex < closest_distance_end && closest_distance_vertex < 3.0 ) { // The second cut is to ensure that this is a track that actually originates at the vertex.          

	  // Use the length of this track to find what its kinetic energy would be if it was a proton.                                                                                               
	  double proton_candidate_length   = track.Length();

	  double proton_x_projection         = 0.;
	  double proton_y_projection         = 0.;
	  double proton_z_projection         = 0.;

	  int    num_track_trajectory_points = track.NumberTrajectoryPoints();

	  if ( num_track_trajectory_points > 0 ) {

	    proton_x_projection       = ( track.LocationAtPoint( last_valid_point_idx ).X() - track.LocationAtPoint( first_valid_point_idx ).X() );
	    proton_y_projection       = ( track.LocationAtPoint( last_valid_point_idx ).Y() - track.LocationAtPoint( first_valid_point_idx ).Y() );
	    proton_z_projection       = ( track.LocationAtPoint( last_valid_point_idx ).Z() - track.LocationAtPoint( first_valid_point_idx ).Z() );

	    // Flip the direction of the 
	    if ( fabs( closest_distance_vertex - distance_vertex_point1 ) < 0.001 ) {
	      
	      proton_x_projection    *= -1.0;
	      proton_y_projection    *= -1.0;
	      proton_z_projection    *= -1.0;

	    } 

	  }

	  // Calculate the momentum and its components.
	  double proton_candidate_momentum = ( p_calculator_from_length.GetTrackMomentum( proton_candidate_length, 2212 ) * 1000. ); 

	  x_momentum_component_sum_additional_tracks += ( proton_candidate_momentum ) * ( proton_x_projection / proton_candidate_length );
	  y_momentum_component_sum_additional_tracks += ( proton_candidate_momentum ) * ( proton_y_projection / proton_candidate_length );
	  z_momentum_component_sum_additional_tracks += ( proton_candidate_momentum ) * ( proton_z_projection / proton_candidate_length );
	  
	  double proton_candidate_KE    = ( TMath::Sqrt( ( p_calculator_from_length.GetTrackMomentum( proton_candidate_length, 2212 ) * 1000. ) * ( p_calculator_from_length.GetTrackMomentum( proton_candidate_length, 2212 ) * 1000. ) + 938.3 * 938.3 ) - 938.3 );

	  proton_kinetic_energy_from_additional_tracks += proton_candidate_KE;

	  idx_of_additional_tracks_originating_from_vertex.push_back( tpcobject_track_iter );

	  num_additional_tracks_originating_from_vertex++;

	} // End of the requirement that the track be closer to the vertex than the other end.                                                                                                        

      } // End of the loop over the tracks in the event.

    num_of_tracks_originating_from_vertex = ( num_additional_tracks_originating_from_vertex + 1 );

    // Set the value for the number of points above threshold to 0.
    num_points_above_threshold = 0;

    num_u_plane_points_above_threshold_in_vertexing_loop = 0;
    num_v_plane_points_above_threshold_in_vertexing_loop = 0;
    num_y_plane_points_above_threshold_in_vertexing_loop = 0;

    u_plane_muon_track_start_x                           = 0.;
    u_plane_muon_track_start_y                           = 0.;
    u_plane_muon_track_start_z                           = 0.;
    u_plane_muon_track_end_x                             = 0.;
    u_plane_muon_track_end_y                             = 0.;
    u_plane_muon_track_end_z                             = 0.;

    v_plane_muon_track_start_x                           = 0.;
    v_plane_muon_track_start_y                           = 0.;
    v_plane_muon_track_start_z                           = 0.;
    v_plane_muon_track_end_x                             = 0.;
    v_plane_muon_track_end_y                             = 0.;
    v_plane_muon_track_end_z                             = 0.;

    y_plane_muon_track_start_x                           = 0.;
    y_plane_muon_track_start_y                           = 0.;
    y_plane_muon_track_start_z                           = 0.;
    y_plane_muon_track_end_x                             = 0.;
    y_plane_muon_track_end_y                             = 0.;
    y_plane_muon_track_end_z                             = 0.;
    
    // Set the vertex location variables equal to the location of the end of the track.
    // That way, when you look at the hits below, it will not matter if you have zero calorimetry points above threshold on the track.
    // U Plane.
    u_plane_vertex_location_x                            = vertex_location_x;
    u_plane_vertex_location_y                            = vertex_location_y;
    u_plane_vertex_location_z                            = vertex_location_z;
    u_plane_vertex_location_x_with_SCE                   = vertex_location_x_with_SCE;
    u_plane_vertex_location_y_with_SCE                   = vertex_location_y_with_SCE;
    u_plane_vertex_location_z_with_SCE                   = vertex_location_z_with_SCE;

    // V Plane. 
    v_plane_vertex_location_x                            = vertex_location_x;
    v_plane_vertex_location_y                            = vertex_location_y;
    v_plane_vertex_location_z                            = vertex_location_z;
    v_plane_vertex_location_x_with_SCE                   = vertex_location_x_with_SCE;
    v_plane_vertex_location_y_with_SCE                   = vertex_location_y_with_SCE;
    v_plane_vertex_location_z_with_SCE                   = vertex_location_z_with_SCE;

    // Y Plane.
    y_plane_vertex_location_x                            = vertex_location_x;
    y_plane_vertex_location_y                            = vertex_location_y;
    y_plane_vertex_location_z                            = vertex_location_z;
    y_plane_vertex_location_x_with_SCE                   = vertex_location_x_with_SCE;
    y_plane_vertex_location_y_with_SCE                   = vertex_location_y_with_SCE;
    y_plane_vertex_location_z_with_SCE                   = vertex_location_z_with_SCE;

    // Introduce the vertex coordinates and ticks here.
    double u_plane_wire_at_vertex                        = -9999.;
    double u_plane_vertex_tick                           = -9999.;
    double u_plane_wire_at_vertex_integral               = -9999.;

    double v_plane_wire_at_vertex                        = -9999.;
    double v_plane_vertex_tick                           = -9999.;
    double v_plane_wire_at_vertex_integral               = -9999.;

    double y_plane_wire_at_vertex                        = -9999.;
    double y_plane_vertex_tick                           = -9999.;
    double y_plane_wire_at_vertex_integral               = -9999.;

    // This is the conditional for whether the hit is bad (in the spikes) or not.
    fill_vector                   = true;
    count                         = 0;

    // This is going to be a loop over the 'meta' information on each plane again. 
    for ( size_t pl = 0; pl < 3; pl++ ) {

      if ( Calo_v[pl]->PlaneID().Plane != 0 && Calo_v[pl]->PlaneID().Plane != 1 && Calo_v[pl]->PlaneID().Plane != 2 )
        continue;

      // Clear out the residual range vector and the dEdx vector.                                                                                                                                
      residualrange_vector.clear();
      dedx_vector.clear();

      auto calo_orig                     = Calo_v.at( pl );
      auto xyz_v_orig                    = calo_orig->XYZ();
      auto dedx_v_orig                   = calo_orig->dEdx();
      auto rr_v_orig                     = calo_orig->ResidualRange();
      auto calo_no_SCE_corrections_orig  = Calo_v_no_SCE_corrections.at( pl );
      auto xyz_v_no_SCE_corrections_orig = calo_no_SCE_corrections_orig->XYZ();

      // Declare the new vectors that you will fill later on in the loop.
      std::vector< float > xyz_v_x_coords;
      xyz_v_x_coords.clear();

      std::vector< float > xyz_v_y_coords;
      xyz_v_y_coords.clear();
	  
      std::vector< float > xyz_v_z_coords;
      xyz_v_z_coords.clear();

      std::vector< float > dedx_v;
      dedx_v.clear();

      std::vector< float > rr_v;
      rr_v.clear();

      std::vector< float > xyz_v_no_SCE_corrections_x_coords;
      xyz_v_no_SCE_corrections_x_coords.clear();

      std::vector< float > xyz_v_no_SCE_corrections_y_coords;
      xyz_v_no_SCE_corrections_y_coords.clear();

      std::vector< float > xyz_v_no_SCE_corrections_z_coords;
      xyz_v_no_SCE_corrections_z_coords.clear();

      std::vector< float > hit_tick_values;
      hit_tick_values.clear();

      std::vector< float > hit_wire_values;
      hit_wire_values.clear();

      std::vector< float > hit_integral_values;
      hit_integral_values.clear();

      // Set the value of 'num_calorimetry_points'.                                                                                                                                                   
      num_calo_points = xyz_v_orig.size();

      // Loop through the vector and fill the relevant vectors.                                                                                                                                         
      for ( size_t hit_iter = 0; hit_iter < vmeta.size(); hit_iter++ ) {

        fill_vector = true;

	int ind = vmeta[hit_iter]->Index();

        // check that the traj point is in the calorimetry point                                                                                                                                        
	// and belongs to the plane we are interested in                                                                                                                                                
        if( pandora_track_h->at( muon_candidate_idx ).HasValidPoint(ind) && vhit[hit_iter]->WireID().Plane == pl && count < int( dedx_v_orig.size() ) ){

          auto const& hit = vhit.at( hit_iter );

          // Check to make sure the hit is not in one of the spikes.                                                                                                                                    
          float rms = hit->RMS();

          // Filter on the peaks (we don't want to use these in those quantities).                                                                                                                       
          if ( fabs( rms - 0.333 ) < 0.0001 || fabs( rms - 0.5 ) < 0.0001 || fabs( rms - 0.667 ) < 0.0001 || fabs( rms - 1.0 ) < 0.0001 || fabs( rms - 1.333 ) < 0.0001 || fabs( rms - 1.667 ) < 0.0001 || fabs( rms - 2.0 ) < 0.0001 || fabs( rms - 2.333 ) < 0.0001 || fabs( rms - 2.667 ) < 0.0001 || fabs( rms - 2.8572 ) < 0.0001 || fabs( rms - 3.0 ) < 0.0001 || fabs( rms - 3.333 ) < 0.0001 || fabs( rms - 3.667 ) < 0.0001 || fabs( rms - 4.0 ) < 0.0001 || fabs( rms - 4.333 ) < 0.0001 || fabs( rms - 4.667 ) < 0.0001 || fabs( rms - 5.0 ) < 0.0001 || fabs( rms - 5.333 ) < 0.0001 || fabs( rms - 5.667 ) < 0.0001 || fabs( rms - 6.0 ) < 0.0001 || fabs( rms - 6.333 ) < 0.0001 || fabs( rms - 6.5 ) < 0.0001 || fabs( rms - 7.0 ) < 0.0001 || fabs( rms - 7.333 ) < 0.0001 || fabs( rms - 7.5 ) < 0.0001 || fabs( rms - 7.667 ) < 0.0001 || fabs( rms - 8.0 ) < 0.0001 || fabs( rms - 8.5 ) < 0.0001 || fabs( rms - 9.0 ) < 0.0001 || fabs( rms - 9.5 ) < 0.0001 || fabs( rms - 10.0 ) < 0.0001 || fabs( rms - 10.5 ) < 0.0001 || fabs( rms - 11.0 ) < 0.0001 || fabs( rms - 11.5 ) < 0.0001 || fabs( rms - 12.0 ) < 0.0001 || fabs( rms - 12.5 ) < 0.0001 || fabs( rms - 13.0 ) < 0.0001 || fabs( rms - 13.5 ) < 0.0001 || fabs( rms - 14.0 ) < 0.0001 || fabs( rms - 14.5 ) < 0.0001 || fabs( rms - 15.0 ) < 0.0001 || fabs( rms - 15.5 ) < 0.0001 || fabs( rms - 16.0 ) < 0.0001 || fabs( rms - 16.5 ) < 0.0001 || fabs( rms - 17.0 ) < 0.0001 || fabs( rms - 17.5 ) < 0.0001 || fabs( rms - 18.0 ) < 0.0001 || fabs( rms - 18.5 ) < 0.0001 || fabs( rms - 19.0 ) < 0.0001 || fabs( rms - 19.5 ) < 0.0001 || fabs( rms - 20.0 ) < 0.0001 )
            fill_vector = false;

          if ( fill_vector == true ) {

            dedx_v.push_back( dedx_v_orig.at(count) );
            xyz_v_x_coords.push_back( xyz_v_orig.at(count).X() );
            xyz_v_no_SCE_corrections_x_coords.push_back( xyz_v_no_SCE_corrections_orig.at(count).X() );
            xyz_v_y_coords.push_back( xyz_v_orig.at(count).Y() );
            xyz_v_no_SCE_corrections_y_coords.push_back( xyz_v_no_SCE_corrections_orig.at(count).Y() );
            xyz_v_z_coords.push_back( xyz_v_orig.at(count).Z() );
            xyz_v_no_SCE_corrections_z_coords.push_back( xyz_v_no_SCE_corrections_orig.at(count).Z() );
	    hit_tick_values.push_back( ( ( hit->StartTick() + hit->EndTick() ) / 2.0 ) + 2400. );
	    hit_wire_values.push_back( hit->WireID().Wire );
	    hit_integral_values.push_back( hit->Integral() );
	    rr_v.push_back( rr_v_orig.at( count ) );

          }

          count++;

        } // End of the conditional that the point actually exists along the track.                                                                                                                       

      } // End of the loop over the hits.                                                

      int u_plane_vertex_idx = 0;
      int v_plane_vertex_idx = 0;
      int y_plane_vertex_idx = 0;

      // Continue on to find the vertex coordinates.
      if ( xyz_v_x_coords.size() > 0 ) {
		  
	// These loops are going to loop over the points on the track from the end furthest from the vertex.
	// Find the vertex using U plane info.                                                                                                                                                       
	if ( Calo_v.at( pl )->PlaneID().Plane == 0 ) {
	  
	  has_u_plane_info = 1;
	  
	  if ( first_points_of_track_are_start ) {

	    if ( xyz_v_x_coords.size() > 0 ) {

	      u_plane_vertex_idx         = 0;

	      u_plane_muon_track_start_x = ( xyz_v_x_coords.at( 0 ) - ( 0.111436 * flash_time ) );
	      u_plane_muon_track_start_y = xyz_v_y_coords.at( 0 );
	      u_plane_muon_track_start_z = xyz_v_z_coords.at( 0 );
	      u_plane_muon_track_end_x   = ( xyz_v_x_coords.at( xyz_v_x_coords.size() - 1 ) - ( 0.111436 * flash_time ) );
	      u_plane_muon_track_end_y   = xyz_v_y_coords.at( xyz_v_y_coords.size() - 1 );
	      u_plane_muon_track_end_z   = xyz_v_z_coords.at( xyz_v_z_coords.size() - 1 );

	      u_plane_wire_at_vertex     = hit_wire_values.at( 0 );
	      u_plane_vertex_tick        = hit_tick_values.at( 0 );
	      u_plane_wire_at_vertex_integral = hit_integral_values.at( 0 );

	      u_plane_wire_at_other_end  = hit_wire_values.at( xyz_v_x_coords.size() - 1 );
	      u_plane_other_end_tick     = hit_wire_values.at( xyz_v_x_coords.size() - 1 );

	    }
	    
	    double furthest_distance_from_vertex = -1.;
	    
	    for ( int calo_point_iter = 0; calo_point_iter < int(xyz_v_x_coords.size()); calo_point_iter++ ) {

	      double distance = TMath::Sqrt( ( xyz_v_x_coords.at( xyz_v_x_coords.size() - calo_point_iter - 1 ) - xyz_v_x_coords.at( 0 ) ) * ( xyz_v_x_coords.at( xyz_v_x_coords.size() - calo_point_iter - 1 ) - xyz_v_x_coords.at( 0 ) ) + ( xyz_v_y_coords.at( xyz_v_y_coords.size() - calo_point_iter - 1 ) - xyz_v_y_coords.at( 0 ) ) * ( xyz_v_y_coords.at( xyz_v_y_coords.size() - calo_point_iter - 1 ) - xyz_v_y_coords.at( 0 ) ) + ( xyz_v_z_coords.at( xyz_v_z_coords.size() - calo_point_iter - 1 ) - xyz_v_z_coords.at( 0 ) ) * ( xyz_v_z_coords.at( xyz_v_z_coords.size() - calo_point_iter - 1 ) - xyz_v_z_coords.at( 0 ) ) );

	      if ( distance < maximum_proton_length && dedx_v.at( xyz_v_x_coords.size() - calo_point_iter - 1 ) > fUPlaneThreshold ) {

                num_u_plane_points_above_threshold_in_vertexing_loop++;

              }
	      
	      if ( distance > furthest_distance_from_vertex && distance < maximum_proton_length && dedx_v.at( xyz_v_x_coords.size() - calo_point_iter - 1 ) > fUPlaneThreshold ) {
		  
		has_u_plane_points_to_use_for_vertex     = 1;
		
		u_plane_vertex_idx                       = int( xyz_v_x_coords.size() - calo_point_iter - 1 );

		furthest_distance_from_vertex            = distance;
		u_plane_vertex_location_x                = ( xyz_v_x_coords.at( xyz_v_x_coords.size() - calo_point_iter - 1 ) - ( 0.111436 * flash_time ) );
		u_plane_vertex_location_y                = xyz_v_y_coords.at( xyz_v_y_coords.size() - calo_point_iter - 1 );
		u_plane_vertex_location_z                = xyz_v_z_coords.at( xyz_v_z_coords.size() - calo_point_iter - 1 );
		
		u_plane_truth_reco_vertex_distance       = TMath::Sqrt( ( nu_vtx_x_truth - u_plane_vertex_location_x ) * ( nu_vtx_x_truth - u_plane_vertex_location_x ) + ( nu_vtx_y_truth - u_plane_vertex_location_y ) * ( nu_vtx_y_truth - u_plane_vertex_location_y ) + ( nu_vtx_z_truth - u_plane_vertex_location_z ) * ( nu_vtx_z_truth - u_plane_vertex_location_z ) );
		
		u_plane_vertex_location_x_with_SCE       = ( xyz_v_no_SCE_corrections_x_coords.at( xyz_v_no_SCE_corrections_x_coords.size() - calo_point_iter - 1 ) - ( 0.111436 * flash_time ) );
		u_plane_vertex_location_y_with_SCE       = xyz_v_no_SCE_corrections_y_coords.at( xyz_v_no_SCE_corrections_y_coords.size() - calo_point_iter - 1 );
		u_plane_vertex_location_z_with_SCE       = xyz_v_no_SCE_corrections_z_coords.at( xyz_v_no_SCE_corrections_z_coords.size() - calo_point_iter - 1 );

		u_plane_wire_at_vertex                   = hit_wire_values.at( xyz_v_no_SCE_corrections_x_coords.size() - calo_point_iter - 1 );
		u_plane_vertex_tick                      = hit_tick_values.at( xyz_v_no_SCE_corrections_x_coords.size() - calo_point_iter - 1 );
		u_plane_wire_at_vertex_integral          = hit_integral_values.at( xyz_v_no_SCE_corrections_x_coords.size() - calo_point_iter - 1 );

	      }
	      
	    }
	    
	  } // End of the forward orientation case.

	  else {
	    

	    if ( xyz_v_x_coords.size() > 0 ) {

	      u_plane_vertex_idx                        = int( xyz_v_x_coords.size() - 1 );

              u_plane_muon_track_start_x                = ( xyz_v_x_coords.at( xyz_v_x_coords.size() - 1 ) - ( 0.111436 * flash_time ) );
              u_plane_muon_track_start_y                = xyz_v_y_coords.at( xyz_v_y_coords.size() - 1 );
	      u_plane_muon_track_start_z                = xyz_v_z_coords.at( xyz_v_z_coords.size() - 1 );
              u_plane_muon_track_end_x                  = ( xyz_v_x_coords.at( 0 ) - ( 0.111436 * flash_time ) );
              u_plane_muon_track_end_y                  = xyz_v_y_coords.at( 0 );
              u_plane_muon_track_end_z                  = xyz_v_z_coords.at( 0 );

	      u_plane_wire_at_vertex                    = hit_wire_values.at( xyz_v_x_coords.size() - 1 );
              u_plane_vertex_tick                       = hit_tick_values.at( xyz_v_x_coords.size() - 1 );
	      u_plane_wire_at_vertex_integral           = hit_integral_values.at( xyz_v_x_coords.size() - 1 );

	      u_plane_wire_at_other_end                 = hit_wire_values.at( 0 );
              u_plane_other_end_tick                    = hit_wire_values.at( 0 );

	    }   

	    double furthest_distance_from_vertex = -1.;
	    
	    for ( int calo_point_iter = 0; calo_point_iter < int(xyz_v_x_coords.size()); calo_point_iter++ ) {
	      
	      double distance = TMath::Sqrt( ( xyz_v_x_coords.at( calo_point_iter ) - xyz_v_x_coords.at( xyz_v_x_coords.size() - 1 ) ) * ( xyz_v_x_coords.at( calo_point_iter ) - xyz_v_x_coords.at( xyz_v_x_coords.size() - 1 ) ) + ( xyz_v_y_coords.at( calo_point_iter ) - xyz_v_y_coords.at( xyz_v_y_coords.size() - 1 ) ) * ( xyz_v_y_coords.at( calo_point_iter ) - xyz_v_y_coords.at( xyz_v_y_coords.size() - 1 ) ) + ( xyz_v_z_coords.at( calo_point_iter ) - xyz_v_z_coords.at( xyz_v_z_coords.size() - 1 ) ) * ( xyz_v_z_coords.at( calo_point_iter ) - xyz_v_z_coords.at( xyz_v_z_coords.size() - 1 ) ) );

	      if ( distance < maximum_proton_length && dedx_v.at( calo_point_iter ) > fUPlaneThreshold ) {

                num_u_plane_points_above_threshold_in_vertexing_loop++;

              }
									
	      if ( distance > furthest_distance_from_vertex && distance < maximum_proton_length && dedx_v.at( calo_point_iter ) > fUPlaneThreshold ) {

		has_u_plane_points_to_use_for_vertex     = 1;

		u_plane_vertex_idx                       = calo_point_iter;
		
		furthest_distance_from_vertex            = distance;
		u_plane_vertex_location_x                = ( xyz_v_x_coords.at( calo_point_iter ) - ( 0.111436 * flash_time ) );
		u_plane_vertex_location_y                = xyz_v_y_coords.at( calo_point_iter );
		u_plane_vertex_location_z                = xyz_v_z_coords.at( calo_point_iter );

		u_plane_truth_reco_vertex_distance       = TMath::Sqrt(( nu_vtx_x_truth - u_plane_vertex_location_x ) * ( nu_vtx_x_truth - u_plane_vertex_location_x ) + ( nu_vtx_y_truth - u_plane_vertex_location_y ) * ( nu_vtx_y_truth - u_plane_vertex_location_y ) + ( nu_vtx_z_truth - u_plane_vertex_location_z ) * ( nu_vtx_z_truth - u_plane_vertex_location_z ) );
		
		u_plane_vertex_location_x_with_SCE       = ( xyz_v_no_SCE_corrections_x_coords.at( calo_point_iter ) - ( 0.111436 * flash_time ) );
		u_plane_vertex_location_y_with_SCE       = xyz_v_no_SCE_corrections_y_coords.at( calo_point_iter );
		u_plane_vertex_location_z_with_SCE       = xyz_v_no_SCE_corrections_z_coords.at( calo_point_iter );

		u_plane_wire_at_vertex                   = hit_wire_values.at( calo_point_iter );
		u_plane_vertex_tick                      = hit_tick_values.at( calo_point_iter );
		u_plane_wire_at_vertex_integral          = hit_integral_values.at( calo_point_iter );

	      }
	      
	    }

	  } // End of the backwards orientation case. 

	} // End of the case of the U Plane.
	
	if ( Calo_v.at( pl )->PlaneID().Plane == 1 ) {
	  
	  has_v_plane_info = 1;
	  
	  // Find the vertex using collection plane info.                                                                                                                                                
	  if ( first_points_of_track_are_start ) {

	    if ( xyz_v_x_coords.size() > 0 ) {

	      v_plane_vertex_idx                         = 0;

              v_plane_muon_track_start_x                 = ( xyz_v_x_coords.at( 0 ) - ( 0.111436 * flash_time ) );
              v_plane_muon_track_start_y                 = xyz_v_y_coords.at( 0 );
              v_plane_muon_track_start_z                 = xyz_v_z_coords.at( 0 );
              v_plane_muon_track_end_x                   = ( xyz_v_x_coords.at( xyz_v_x_coords.size() - 1 ) - ( 0.111436 * flash_time ) );
              v_plane_muon_track_end_y                   = xyz_v_y_coords.at( xyz_v_y_coords.size() - 1 );
              v_plane_muon_track_end_z                   = xyz_v_z_coords.at( xyz_v_z_coords.size() - 1 );

	      v_plane_wire_at_vertex                     = hit_wire_values.at( 0 );
              v_plane_vertex_tick                        = hit_tick_values.at( 0 );
	      v_plane_wire_at_vertex_integral            = hit_integral_values.at( 0 );

	      v_plane_wire_at_other_end                  = hit_wire_values.at( xyz_v_x_coords.size() - 1 );
              v_plane_other_end_tick                     = hit_wire_values.at( xyz_v_x_coords.size() - 1 );

            }
	      
	    double furthest_distance_from_vertex = -1.;

	    for ( int calo_point_iter = 0; calo_point_iter < int(xyz_v_x_coords.size()); calo_point_iter++ ) {

	      double distance = TMath::Sqrt( ( xyz_v_x_coords.at( xyz_v_x_coords.size() - calo_point_iter - 1 ) - xyz_v_x_coords.at( 0 ) ) * ( xyz_v_x_coords.at( xyz_v_x_coords.size() - calo_point_iter - 1 ) - xyz_v_x_coords.at( 0 ) ) + ( xyz_v_y_coords.at( xyz_v_y_coords.size() - calo_point_iter - 1 ) - xyz_v_y_coords.at( 0 ) ) * ( xyz_v_y_coords.at( xyz_v_y_coords.size() - calo_point_iter - 1 ) - xyz_v_y_coords.at( 0 ) ) + ( xyz_v_z_coords.at( xyz_v_z_coords.size() - calo_point_iter - 1 ) - xyz_v_z_coords.at( 0 ) ) * ( xyz_v_z_coords.at( xyz_v_z_coords.size() - calo_point_iter - 1 ) - xyz_v_z_coords.at( 0 ) ) );

	      if ( distance < maximum_proton_length && dedx_v.at( xyz_v_x_coords.size() - calo_point_iter - 1 ) > fVPlaneThreshold ) {
		
		num_v_plane_points_above_threshold_in_vertexing_loop++;

	      }
	      
	      if ( distance > furthest_distance_from_vertex && distance < maximum_proton_length && dedx_v.at( xyz_v_x_coords.size() - calo_point_iter - 1 ) > fVPlaneThreshold ) {
		
		has_v_plane_points_to_use_for_vertex     = 1;

		v_plane_vertex_idx                       = int( xyz_v_x_coords.size() - calo_point_iter - 1 );

		furthest_distance_from_vertex            = distance;
		v_plane_vertex_location_x                = ( xyz_v_x_coords.at( xyz_v_x_coords.size() - calo_point_iter - 1 ) - ( 0.111436 * flash_time ) );
		v_plane_vertex_location_y                = xyz_v_y_coords.at( xyz_v_y_coords.size() - calo_point_iter - 1 );
		v_plane_vertex_location_z                = xyz_v_z_coords.at( xyz_v_z_coords.size() - calo_point_iter - 1 );

		v_plane_truth_reco_vertex_distance       = TMath::Sqrt( ( nu_vtx_x_truth - v_plane_vertex_location_x ) * ( nu_vtx_x_truth - v_plane_vertex_location_x ) + ( nu_vtx_y_truth - v_plane_vertex_location_y ) *  ( nu_vtx_y_truth - v_plane_vertex_location_y ) + ( nu_vtx_z_truth - v_plane_vertex_location_z ) * ( nu_vtx_z_truth - v_plane_vertex_location_z ) );
		
		v_plane_vertex_location_x_with_SCE       = ( xyz_v_no_SCE_corrections_x_coords.at( xyz_v_no_SCE_corrections_x_coords.size() - calo_point_iter - 1 ) - ( 0.111436 * flash_time ) );
		v_plane_vertex_location_y_with_SCE       = xyz_v_no_SCE_corrections_y_coords.at( xyz_v_no_SCE_corrections_y_coords.size() - calo_point_iter - 1);
		v_plane_vertex_location_z_with_SCE       = xyz_v_no_SCE_corrections_z_coords.at( xyz_v_no_SCE_corrections_z_coords.size() - calo_point_iter - 1);
	       
		v_plane_wire_at_vertex                   = hit_wire_values.at( xyz_v_no_SCE_corrections_x_coords.size() - calo_point_iter - 1 );
                v_plane_vertex_tick                      = hit_tick_values.at( xyz_v_no_SCE_corrections_x_coords.size() - calo_point_iter - 1 );
		v_plane_wire_at_vertex_integral          = hit_integral_values.at( xyz_v_no_SCE_corrections_x_coords.size() - calo_point_iter - 1 );

	      }

	    } 

	  } // End of the forward orientation case. 

	  else {

	    if ( xyz_v_x_coords.size() > 0 ) {

	      v_plane_vertex_idx                         = int( xyz_v_x_coords.size() - 1 );

              v_plane_muon_track_start_x                 = ( xyz_v_x_coords.at( xyz_v_x_coords.size() - 1 ) - ( 0.111436 * flash_time ) );
              v_plane_muon_track_start_y                 = xyz_v_y_coords.at( xyz_v_y_coords.size() - 1 );
              v_plane_muon_track_start_z                 = xyz_v_z_coords.at( xyz_v_z_coords.size() - 1 );
              v_plane_muon_track_end_x                   = ( xyz_v_x_coords.at( 0 ) - ( 0.111436 * flash_time ) );
              v_plane_muon_track_end_y                   = xyz_v_y_coords.at( 0 );
              v_plane_muon_track_end_z                   = xyz_v_z_coords.at( 0 );

	      v_plane_wire_at_vertex                     = hit_wire_values.at( xyz_v_x_coords.size() - 1 );
              v_plane_vertex_tick                        = hit_tick_values.at( xyz_v_x_coords.size() - 1 );
	      v_plane_wire_at_vertex_integral            = hit_integral_values.at( xyz_v_x_coords.size() - 1 );

	      v_plane_wire_at_other_end                  = hit_wire_values.at( 0 );
              v_plane_other_end_tick                     = hit_wire_values.at( 0 );

            }
	    
	    double furthest_distance_from_vertex = -1.;
	    
	    for ( int calo_point_iter = 0; calo_point_iter < int(xyz_v_x_coords.size()); calo_point_iter++ ) {
	      
	      double distance = TMath::Sqrt( ( xyz_v_x_coords.at( calo_point_iter ) - xyz_v_x_coords.at( xyz_v_x_coords.size() - 1 ) ) * ( xyz_v_x_coords.at( calo_point_iter ) - xyz_v_x_coords.at( xyz_v_x_coords.size() - 1 ) ) + ( xyz_v_y_coords.at( calo_point_iter ) - xyz_v_y_coords.at( xyz_v_y_coords.size() - 1 ) ) * ( xyz_v_y_coords.at( calo_point_iter ) - xyz_v_y_coords.at( xyz_v_y_coords.size() - 1 ) ) + ( xyz_v_z_coords.at( calo_point_iter ) - xyz_v_z_coords.at( xyz_v_z_coords.size() - 1 ) ) * ( xyz_v_z_coords.at( calo_point_iter ) - xyz_v_z_coords.at( xyz_v_z_coords.size() - 1 ) ) );

	      if ( distance < maximum_proton_length && dedx_v.at( calo_point_iter ) > fVPlaneThreshold ) {
		
		num_v_plane_points_above_threshold_in_vertexing_loop++;

	      } 

	      if ( distance > furthest_distance_from_vertex && distance < maximum_proton_length && dedx_v.at( calo_point_iter ) > fVPlaneThreshold ) {

		has_v_plane_points_to_use_for_vertex      = 1;

		v_plane_vertex_idx                        = calo_point_iter;

		furthest_distance_from_vertex             = distance;
		v_plane_vertex_location_x                 = ( xyz_v_x_coords.at( calo_point_iter ) - ( 0.111436 * flash_time ) );
		v_plane_vertex_location_y                 = xyz_v_y_coords.at( calo_point_iter );
		v_plane_vertex_location_z                 = xyz_v_z_coords.at( calo_point_iter );
		
		v_plane_truth_reco_vertex_distance        = TMath::Sqrt(( nu_vtx_x_truth - v_plane_vertex_location_x ) * ( nu_vtx_x_truth - v_plane_vertex_location_x ) + ( nu_vtx_y_truth - v_plane_vertex_location_y ) *( nu_vtx_y_truth - v_plane_vertex_location_y ) + ( nu_vtx_z_truth - v_plane_vertex_location_z ) * ( nu_vtx_z_truth - v_plane_vertex_location_z ) );
		
		v_plane_vertex_location_x_with_SCE        = ( xyz_v_no_SCE_corrections_x_coords.at( calo_point_iter ) - ( 0.111436 * flash_time ) );
		v_plane_vertex_location_y_with_SCE        = xyz_v_no_SCE_corrections_y_coords.at( calo_point_iter );
		v_plane_vertex_location_z_with_SCE        = xyz_v_no_SCE_corrections_z_coords.at( calo_point_iter );

		v_plane_wire_at_vertex                    = hit_wire_values.at( calo_point_iter );
                v_plane_vertex_tick                       = hit_tick_values.at( calo_point_iter );
		v_plane_wire_at_vertex_integral           = hit_integral_values.at( calo_point_iter );

	      }

	    }

	  } // End of the backwards orientation case. 

	} // End of the case of the V Plane.

	if ( Calo_v.at( pl )->PlaneID().Plane == 2 ) {

	  has_y_plane_info = 1;

	  // Find the vertex using collection plane info.                                                                                                                                               
	  if ( first_points_of_track_are_start ) {

	    if ( xyz_v_x_coords.size() > 0 ) {

	      y_plane_vertex_idx                          = 0;

              y_plane_muon_track_start_x                  = ( xyz_v_x_coords.at( 0 ) - ( 0.111436 * flash_time ) );
              y_plane_muon_track_start_y                  = xyz_v_y_coords.at( 0 );
              y_plane_muon_track_start_z                  = xyz_v_z_coords.at( 0 );
              y_plane_muon_track_end_x                    = ( xyz_v_x_coords.at( xyz_v_x_coords.size() - 1 ) - ( 0.111436 * flash_time ) );
              y_plane_muon_track_end_y                    = xyz_v_y_coords.at( xyz_v_y_coords.size() - 1 );
              y_plane_muon_track_end_z                    = xyz_v_z_coords.at( xyz_v_z_coords.size() - 1 );

	      y_plane_wire_at_vertex                      = hit_wire_values.at( 0 );
              y_plane_vertex_tick                         = hit_tick_values.at( 0 );
	      y_plane_wire_at_vertex_integral             = hit_integral_values.at( 0 );

	      y_plane_wire_at_other_end                   = hit_wire_values.at( xyz_v_x_coords.size() - 1 );
              y_plane_other_end_tick                      = hit_wire_values.at( xyz_v_x_coords.size() - 1 );

            }

	    double furthest_distance_from_vertex = -1.;

	    for ( int calo_point_iter = 0; calo_point_iter < int(xyz_v_x_coords.size()); calo_point_iter++ ) {

	      double distance = TMath::Sqrt( ( xyz_v_x_coords.at( xyz_v_x_coords.size() - calo_point_iter - 1 ) - xyz_v_x_coords.at( 0 ) ) * ( xyz_v_x_coords.at( xyz_v_x_coords.size() - calo_point_iter - 1 ) - xyz_v_x_coords.at( 0 ) ) + ( xyz_v_y_coords.at( xyz_v_y_coords.size() - calo_point_iter - 1 ) - xyz_v_y_coords.at( 0 ) ) * ( xyz_v_y_coords.at( xyz_v_y_coords.size() - calo_point_iter - 1 ) - xyz_v_y_coords.at( 0 ) ) + ( xyz_v_z_coords.at( xyz_v_z_coords.size() - calo_point_iter - 1 ) - xyz_v_z_coords.at( 0 ) ) * ( xyz_v_z_coords.at( xyz_v_z_coords.size() - calo_point_iter - 1 ) - xyz_v_z_coords.at( 0 ) ) );

	      if ( distance < maximum_proton_length && dedx_v.at( xyz_v_x_coords.size() - calo_point_iter - 1 ) > fYPlaneThreshold ) {

                num_y_plane_points_above_threshold_in_vertexing_loop++;

              }

	      if ( distance > furthest_distance_from_vertex && distance < maximum_proton_length && dedx_v.at( xyz_v_x_coords.size() - calo_point_iter - 1 ) > fYPlaneThreshold ) {

		has_y_plane_points_to_use_for_vertex      = 1;

		y_plane_vertex_idx                        = int( xyz_v_x_coords.size() - calo_point_iter - 1 );

		furthest_distance_from_vertex             = distance;
		y_plane_vertex_location_x                 = ( xyz_v_x_coords.at( xyz_v_x_coords.size() - calo_point_iter - 1 ) - ( 0.111436 * flash_time ) );
		y_plane_vertex_location_y                 = xyz_v_y_coords.at( xyz_v_y_coords.size() - calo_point_iter - 1 );
		y_plane_vertex_location_z                 = xyz_v_z_coords.at( xyz_v_z_coords.size() - calo_point_iter - 1 );

		y_plane_truth_reco_vertex_distance        = TMath::Sqrt(( nu_vtx_x_truth - y_plane_vertex_location_x ) * ( nu_vtx_x_truth - y_plane_vertex_location_x ) + ( nu_vtx_y_truth - y_plane_vertex_location_y ) * ( nu_vtx_y_truth - y_plane_vertex_location_y ) + ( nu_vtx_z_truth - y_plane_vertex_location_z ) * ( nu_vtx_z_truth - y_plane_vertex_location_z ) );
		
		y_plane_vertex_location_x_with_SCE        = ( xyz_v_no_SCE_corrections_x_coords.at( xyz_v_no_SCE_corrections_x_coords.size() - calo_point_iter - 1 ) - ( 0.111436 * flash_time ) );
		y_plane_vertex_location_y_with_SCE        = xyz_v_no_SCE_corrections_y_coords.at( xyz_v_no_SCE_corrections_y_coords.size() - calo_point_iter - 1);
		y_plane_vertex_location_z_with_SCE        = xyz_v_no_SCE_corrections_z_coords.at( xyz_v_no_SCE_corrections_z_coords.size() - calo_point_iter - 1);

		y_plane_wire_at_vertex                    = hit_wire_values.at( xyz_v_no_SCE_corrections_x_coords.size() - calo_point_iter - 1 );
                y_plane_vertex_tick                       = hit_tick_values.at( xyz_v_no_SCE_corrections_x_coords.size() - calo_point_iter - 1 );
		y_plane_wire_at_vertex_integral           = hit_integral_values.at( xyz_v_no_SCE_corrections_x_coords.size() - calo_point_iter - 1 );

	      }
	      
	    }
	    
	  } // End of the forward orientation case. 

	  else {

	    if ( xyz_v_x_coords.size() > 0 ) {

	      y_plane_vertex_idx                          = int( xyz_v_x_coords.size() - 1 );

              y_plane_muon_track_start_x                  = ( xyz_v_x_coords.at( xyz_v_x_coords.size() - 1 ) - ( 0.111436 * flash_time ) );
              y_plane_muon_track_start_y                  = xyz_v_y_coords.at( xyz_v_y_coords.size() - 1 );
              y_plane_muon_track_start_z                  = xyz_v_z_coords.at( xyz_v_z_coords.size() - 1 );
              y_plane_muon_track_end_x                    = ( xyz_v_x_coords.at( 0 ) - ( 0.111436 * flash_time ) );
              y_plane_muon_track_end_y                    = xyz_v_y_coords.at( 0 );
              y_plane_muon_track_end_z                    = xyz_v_z_coords.at( 0 );

	      y_plane_wire_at_vertex                      = hit_wire_values.at( xyz_v_x_coords.size() - 1 );
              y_plane_vertex_tick                         = hit_tick_values.at( xyz_v_x_coords.size() - 1 );
	      y_plane_wire_at_vertex_integral             = hit_integral_values.at( xyz_v_x_coords.size() - 1 );

	      y_plane_wire_at_other_end                   = hit_wire_values.at( 0 );
              y_plane_other_end_tick                      = hit_wire_values.at( 0 );

            }

	    double furthest_distance_from_vertex = -1.;
	    
	    for ( int calo_point_iter = 0; calo_point_iter < int(xyz_v_x_coords.size()); calo_point_iter++ ) {
	      
	      double distance = TMath::Sqrt( ( xyz_v_x_coords.at( calo_point_iter ) - xyz_v_x_coords.at( xyz_v_x_coords.size() - 1 ) ) * ( xyz_v_x_coords.at( calo_point_iter ) - xyz_v_x_coords.at( xyz_v_x_coords.size() - 1 ) ) + ( xyz_v_y_coords.at( calo_point_iter ) - xyz_v_y_coords.at( xyz_v_y_coords.size() - 1 ) ) * ( xyz_v_y_coords.at( calo_point_iter ) - xyz_v_y_coords.at( xyz_v_y_coords.size() - 1 ) ) + ( xyz_v_z_coords.at( calo_point_iter ) - xyz_v_z_coords.at( xyz_v_z_coords.size() - 1 ) ) * ( xyz_v_z_coords.at( calo_point_iter ) - xyz_v_z_coords.at( xyz_v_z_coords.size() - 1 ) ) );

	      if ( distance < maximum_proton_length && dedx_v.at( calo_point_iter ) > fYPlaneThreshold ) {

                num_y_plane_points_above_threshold_in_vertexing_loop++;

              }
	      
	      if ( distance > furthest_distance_from_vertex && distance < maximum_proton_length && dedx_v.at( calo_point_iter ) > fYPlaneThreshold ) {

		has_y_plane_points_to_use_for_vertex = 1;

		y_plane_vertex_idx = int( calo_point_iter );
		
		furthest_distance_from_vertex             = distance;
		y_plane_vertex_location_x                 = ( xyz_v_x_coords.at( calo_point_iter ) - ( 0.111436 * flash_time ) );
		y_plane_vertex_location_y                 = xyz_v_y_coords.at( calo_point_iter );
		y_plane_vertex_location_z                 = xyz_v_z_coords.at( calo_point_iter );

		y_plane_truth_reco_vertex_distance        = TMath::Sqrt(( nu_vtx_x_truth - y_plane_vertex_location_x ) * ( nu_vtx_x_truth - y_plane_vertex_location_x ) + ( nu_vtx_y_truth - y_plane_vertex_location_y ) *( nu_vtx_y_truth - y_plane_vertex_location_y ) + ( nu_vtx_z_truth - y_plane_vertex_location_z ) * ( nu_vtx_z_truth - y_plane_vertex_location_z ) );
		
		y_plane_vertex_location_x_with_SCE        = ( xyz_v_no_SCE_corrections_x_coords.at( calo_point_iter ) - ( 0.111436 * flash_time ) );
		y_plane_vertex_location_y_with_SCE        = xyz_v_no_SCE_corrections_y_coords.at( calo_point_iter );
		y_plane_vertex_location_z_with_SCE        = xyz_v_no_SCE_corrections_z_coords.at( calo_point_iter );
		
		y_plane_wire_at_vertex                    = hit_wire_values.at( calo_point_iter );
                y_plane_vertex_tick                       = hit_tick_values.at( calo_point_iter );
		y_plane_wire_at_vertex_integral           = hit_integral_values.at( calo_point_iter );

	      }

	    } // End of the loop of the calorimetry points on this plane.

	  } // End of the backwards orientation case. 

	} // End of the look at the second plane.
      
	if ( first_points_of_track_are_start ) {

	  for ( int calo_point_iter = 0; calo_point_iter < int(xyz_v_x_coords.size()); calo_point_iter++ ) {
	      
	    if ( calo_orig->PlaneID().Plane == 0 ) {
	    
	      vertex_location_x = u_plane_vertex_location_x;
	      vertex_location_y = u_plane_vertex_location_y;
	      vertex_location_z = u_plane_vertex_location_z;
	      
	    }
	    
	    if ( calo_orig->PlaneID().Plane == 1 ) {
	      
	      vertex_location_x = v_plane_vertex_location_x;
	      vertex_location_y = v_plane_vertex_location_y;
	      vertex_location_z = v_plane_vertex_location_z;
		
	    }

	    if ( calo_orig->PlaneID().Plane == 2 ) {

	      vertex_location_x = y_plane_vertex_location_x;
	      vertex_location_y = y_plane_vertex_location_y;
	      vertex_location_z = y_plane_vertex_location_z;
	      
	    }

	    double distance = TMath::Sqrt( ( vertex_location_x - ( xyz_v_x_coords.at( xyz_v_x_coords.size() - calo_point_iter - 1 ) - ( 0.111436 * flash_time ) ) ) * ( vertex_location_x - ( xyz_v_x_coords.at( xyz_v_x_coords.size() - calo_point_iter - 1 ) - ( 0.111436 * flash_time ) ) ) + ( vertex_location_y - xyz_v_y_coords.at( xyz_v_y_coords.size() - calo_point_iter - 1 ) ) * ( vertex_location_y - xyz_v_y_coords.at( xyz_v_y_coords.size() - calo_point_iter - 1 ) ) + ( vertex_location_z - xyz_v_z_coords.at( xyz_v_z_coords.size() - calo_point_iter - 1 ) ) * ( vertex_location_z - xyz_v_z_coords.at( xyz_v_z_coords.size() - calo_point_iter - 1 ) ) );
	    
	    if ( calo_orig->PlaneID().Plane == 0 ) {
	      
	      if ( distance < maximum_proton_length && dedx_v.at( xyz_v_x_coords.size() - calo_point_iter - 1 ) > fUPlaneThreshold && int( xyz_v_x_coords.size() -  calo_point_iter - 1 ) <= u_plane_vertex_idx ) {

		x_coordinates_above_threshold.push_back( xyz_v_x_coords.at( xyz_v_x_coords.size() - calo_point_iter - 1 ) - 0.111436 * flash_time );
		y_coordinates_above_threshold.push_back( xyz_v_y_coords.at( xyz_v_y_coords.size() - calo_point_iter - 1 ) );
		z_coordinates_above_threshold.push_back( xyz_v_z_coords.at( xyz_v_z_coords.size() - calo_point_iter - 1 ) );
		

	      }

	    }
	    
	    if ( calo_orig->PlaneID().Plane == 1 ) {
		
	      if ( distance < maximum_proton_length && dedx_v.at( xyz_v_x_coords.size() - calo_point_iter - 1 ) > fVPlaneThreshold && int( xyz_v_x_coords.size() -  calo_point_iter - 1 ) <= v_plane_vertex_idx ) {

		x_coordinates_above_threshold.push_back( xyz_v_x_coords.at( xyz_v_x_coords.size() - calo_point_iter - 1 ) - 0.111436 * flash_time );
		y_coordinates_above_threshold.push_back( xyz_v_y_coords.at( xyz_v_y_coords.size() - calo_point_iter - 1 ) );
		z_coordinates_above_threshold.push_back( xyz_v_z_coords.at( xyz_v_z_coords.size() - calo_point_iter - 1 ) );
		
	      }

	    }
					     
	    if ( calo_orig->PlaneID().Plane == 2 ) {

	      if ( distance < maximum_proton_length && dedx_v.at( xyz_v_x_coords.size() - calo_point_iter - 1 ) > fYPlaneThreshold && int( xyz_v_x_coords.size() -  calo_point_iter - 1 ) <= y_plane_vertex_idx ) {

		x_coordinates_above_threshold.push_back( xyz_v_x_coords.at( xyz_v_x_coords.size() - calo_point_iter - 1 ) - 0.111436 * flash_time );
		y_coordinates_above_threshold.push_back( xyz_v_y_coords.at( xyz_v_y_coords.size() - calo_point_iter - 1 ) );
		z_coordinates_above_threshold.push_back( xyz_v_z_coords.at( xyz_v_z_coords.size() - calo_point_iter - 1 ) );
		
	      }

	    }

	    residualrange_vector.push_back( rr_v.at( calo_point_iter ) );
	    dedx_vector.push_back( dedx_v.at( xyz_v_x_coords.size() - calo_point_iter - 1 ) );
	    
	  } // End of the loop over the calorimetry points on the track.                                                                                                                                   
	  
	} // End of the conditional that the calorimetry points on the track are placed in the correct order.                                                                                        

	else {

	  for ( int calo_point_iter = 0; calo_point_iter < int(xyz_v_x_coords.size()); calo_point_iter++ ) {
	    
	    if ( calo_orig->PlaneID().Plane == 0 ) {
	      
	      vertex_location_x = u_plane_vertex_location_x;
	      vertex_location_y = u_plane_vertex_location_y;
	      vertex_location_z = u_plane_vertex_location_z;
	      
	    }

	    if ( calo_orig->PlaneID().Plane == 1 ) {

	      vertex_location_x = v_plane_vertex_location_x;
	      vertex_location_y = v_plane_vertex_location_y;
	      vertex_location_z = v_plane_vertex_location_z;
	      
	    }

	    if ( calo_orig->PlaneID().Plane == 2 ) {

	      vertex_location_x = y_plane_vertex_location_x;
	      vertex_location_y = y_plane_vertex_location_y;
	      vertex_location_z = y_plane_vertex_location_z;
		
	    }
	      
	    double distance = TMath::Sqrt( ( vertex_location_x - ( xyz_v_x_coords.at( calo_point_iter ) - 0.111436 * flash_time ) ) * ( vertex_location_x - ( xyz_v_x_coords.at( calo_point_iter ) - 0.111436 * flash_time ) ) + ( vertex_location_y - xyz_v_y_coords.at( calo_point_iter ) ) * ( vertex_location_y - xyz_v_y_coords.at( calo_point_iter ) ) + ( vertex_location_z - xyz_v_z_coords.at( calo_point_iter ) ) * ( vertex_location_z - xyz_v_z_coords.at( calo_point_iter ) ) );

	    if ( calo_orig->PlaneID().Plane == 0 ) {

	      if ( distance < maximum_proton_length && dedx_v.at( calo_point_iter ) > fUPlaneThreshold && int( calo_point_iter ) >= u_plane_vertex_idx ) {

		x_coordinates_above_threshold.push_back( xyz_v_x_coords.at( calo_point_iter ) - 0.111436 * flash_time );
		y_coordinates_above_threshold.push_back( xyz_v_y_coords.at( calo_point_iter ) );
		z_coordinates_above_threshold.push_back( xyz_v_z_coords.at( calo_point_iter ) );

	      }

	    }

	    if ( calo_orig->PlaneID().Plane == 1 ) {
		
	      if ( distance < maximum_proton_length && dedx_v.at( calo_point_iter ) > fVPlaneThreshold && int( calo_point_iter ) >= v_plane_vertex_idx ) {

		x_coordinates_above_threshold.push_back( xyz_v_x_coords.at( calo_point_iter ) - 0.111436 * flash_time );
		y_coordinates_above_threshold.push_back( xyz_v_y_coords.at( calo_point_iter ) );
		z_coordinates_above_threshold.push_back( xyz_v_z_coords.at( calo_point_iter ) );

	      }

	    }
	      
	    if ( calo_orig->PlaneID().Plane == 2 ) {
		
	      if ( distance < maximum_proton_length && dedx_v.at( calo_point_iter ) > fYPlaneThreshold && int( calo_point_iter ) >= y_plane_vertex_idx ) {
		
		x_coordinates_above_threshold.push_back( xyz_v_x_coords.at( calo_point_iter ) - 0.111436 * flash_time );
		y_coordinates_above_threshold.push_back( xyz_v_y_coords.at( calo_point_iter ) );
		z_coordinates_above_threshold.push_back( xyz_v_z_coords.at( calo_point_iter ) );

	      }

	    }

	    residualrange_vector.push_back( rr_v.at( calo_point_iter ) );
	    dedx_vector.push_back( dedx_v.at( calo_point_iter ) );
	      
	  } // End of the loop over the calorimetry points on the track.                                                                                                                              

	} // End of 'else' statement over correct orientation of track points.                

      } // End of the requirement that there be a nonzero number of trajectory points associated to the track through calorimetry.

      // Consider the case when there are zero points on this plane.
      else {

	if ( pl == 0 ) no_calo_points_on_u_plane = true;
	if ( pl == 1 ) no_calo_points_on_v_plane = true;
	if ( pl == 2 ) no_calo_points_on_y_plane = true;

      } // End of the case in which there are 0 points at this end of the track.
	
      count = 0;

    } // End of the loop over the calorimetry object.

    // See if these coordinates are inside/outside the FV.                                                                                                                                              
    if ( vertex_location_x < 20 || vertex_location_x > 236.35 || vertex_location_y < -96.5 || vertex_location_y > 96.5 || vertex_location_z < 20.0 || vertex_location_z > 1016.8 )
      fails_fiducial_volume_cut = true;

    // Use the y-plane info to tell if the track is properly oriented.
    double dist_truth_vertex_to_y_plane_start = TMath::Sqrt( ( nu_vtx_x_truth - y_plane_muon_track_start_x ) * ( nu_vtx_x_truth - y_plane_muon_track_start_x ) + ( nu_vtx_y_truth - y_plane_muon_track_start_y ) * ( nu_vtx_y_truth - y_plane_muon_track_start_y ) + ( nu_vtx_z_truth - y_plane_muon_track_start_z ) * ( nu_vtx_z_truth - y_plane_muon_track_start_z ) );

    double dist_truth_vertex_to_y_plane_end   = TMath::Sqrt( ( nu_vtx_x_truth -y_plane_muon_track_end_x ) * ( nu_vtx_x_truth - y_plane_muon_track_end_x ) +( nu_vtx_y_truth - y_plane_muon_track_end_y ) * ( nu_vtx_y_truth - y_plane_muon_track_end_y )+ ( nu_vtx_z_truth - y_plane_muon_track_end_z ) * ( nu_vtx_z_truth - y_plane_muon_track_end_z ) );

    if ( dist_truth_vertex_to_y_plane_end < dist_truth_vertex_to_y_plane_start ) muon_track_is_correctly_oriented = 0;
	  
    // Set the value of 'num_points_above_threshold'.
    num_points_above_threshold = x_coordinates_above_threshold.size();

    // Print out all of the 3D points above threshold.                                                                                                                                                  
    for ( size_t x_coord_iter = 0; x_coord_iter < x_coordinates_above_threshold.size(); x_coord_iter++ ) {

      x_points_above_threshold[x_coord_iter] = x_coordinates_above_threshold.at( x_coord_iter );
      y_points_above_threshold[x_coord_iter] = y_coordinates_above_threshold.at( x_coord_iter );
      z_points_above_threshold[x_coord_iter] = z_coordinates_above_threshold.at( x_coord_iter );

    }

    // Find the high-charge hits in the event that are located within the 3D radius but are not associated with the track.                                                                              
    // Loop through the hits and make sure that it is not associated with the track.                                                                                                                     
    if ( int( muon_candidate_idx ) != -1 ) {

      // Declare the association for the hits from the longest track.                                                                                                                                   
      std::vector<const recob::Hit*> longest_track_hit_v = trk_hit_assn_v.at( muon_candidate_idx);
	
      for ( size_t hit_iter = 0; hit_iter < hit_h->size(); hit_iter++ ) {

	auto hit = hit_h->at(hit_iter);

	// Put the continue statement in for the 'bad hits'.
	float rms = hit.RMS();

	// Filter on the peaks.
	if ( fabs( rms - 0.333 ) < 0.0001 || fabs( rms - 0.5 ) < 0.0001 || fabs( rms - 0.667 ) < 0.0001 || fabs( rms - 1.0 ) < 0.0001 || fabs( rms - 1.333 ) < 0.0001 || fabs( rms - 1.667 ) < 0.0001 || fabs( rms - 2.0 ) < 0.0001 || fabs( rms - 2.333 ) < 0.0001 || fabs( rms - 2.667 ) < 0.0001 || fabs( rms - 2.8572 ) < 0.0001 || fabs( rms - 3.0 ) < 0.0001 || fabs( rms - 3.333 ) < 0.0001 || fabs( rms - 3.667 ) < 0.0001 || fabs( rms - 4.0 ) < 0.0001 || fabs( rms - 4.333 ) < 0.0001 || fabs( rms - 4.667 ) < 0.0001 || fabs( rms - 5.0 ) < 0.0001 || fabs( rms - 5.333 ) < 0.0001 || fabs( rms - 5.667 ) < 0.0001 || fabs( rms - 6.0 ) < 0.0001 || fabs( rms - 6.333 ) < 0.0001 || fabs( rms - 6.5 ) < 0.0001 || fabs( rms - 7.0 ) < 0.0001 || fabs( rms - 7.333 ) < 0.0001 || fabs( rms - 7.5 ) < 0.0001 || fabs( rms - 7.667 ) < 0.0001 || fabs( rms - 8.0 ) < 0.0001 || fabs( rms - 8.5 ) < 0.0001 || fabs( rms - 9.0 ) < 0.0001 || fabs( rms - 9.5 ) < 0.0001 || fabs( rms - 10.0 ) < 0.0001 || fabs( rms - 10.5 ) < 0.0001 || fabs( rms - 11.0 ) < 0.0001 || fabs( rms - 11.5 ) < 0.0001 || fabs( rms - 12.0 ) < 0.0001 || fabs( rms - 12.5 ) < 0.0001 || fabs( rms - 13.0 ) < 0.0001 || fabs( rms - 13.5 ) < 0.0001 || fabs( rms - 14.0 ) < 0.0001 || fabs( rms - 14.5 ) < 0.0001 || fabs( rms - 15.0 ) < 0.0001 || fabs( rms - 15.5 ) < 0.0001 || fabs( rms - 16.0 ) < 0.0001 || fabs( rms - 16.5 ) < 0.0001 || fabs( rms - 17.0 ) < 0.0001 || fabs( rms - 17.5 ) < 0.0001 || fabs( rms - 18.0 ) < 0.0001 || fabs( rms - 18.5 ) < 0.0001 || fabs( rms - 19.0 ) < 0.0001 || fabs( rms - 19.5 ) < 0.0001 || fabs( rms - 20.0 ) < 0.0001 )
	  continue;

	// See if hit is associated to the longest track.                                                                                                                                                  
	bool is_associated_to_longest_track = false;

	for ( size_t longest_track_hit_iter = 0; longest_track_hit_iter < longest_track_hit_v.size(); longest_track_hit_iter++ ) {

	  const recob::Hit* longest_track_hit = longest_track_hit_v.at( longest_track_hit_iter );

	  // Compare three pieces of information to see if the hit is associated to the longest track in the event.                                                                                       
	  if ( fabs( hit.StartTick() - longest_track_hit->StartTick() ) < 0.001 && fabs( hit.EndTick() - longest_track_hit->EndTick() ) < 0.001 && fabs( hit.Integral() - longest_track_hit->Integral() ) <0.001 ) {
	    is_associated_to_longest_track = true;
	    break;
	  }
	  
	}

	if ( is_associated_to_longest_track == true ) continue;

	bool is_associated_to_another_track_near_vertex = 0.;

	// Do the same for the other tracks that are found to be coming from the vertex.                                                                                                               
	// Loop through the indices of the other tracks found to be coming from the vertex.                                                                                                             
	for ( size_t nearby_track_iter  = 0; nearby_track_iter < idx_of_additional_tracks_originating_from_vertex.size(); nearby_track_iter++ ) {

	  std::vector<const recob::Hit*> track_from_vertex_hit_v = trk_hit_assn_v.at( idx_of_additional_tracks_originating_from_vertex.at( nearby_track_iter ) );

	  for ( size_t nearby_track_hit_iter = 0; nearby_track_hit_iter < track_from_vertex_hit_v.size(); nearby_track_hit_iter++ ) {

	    const recob::Hit* nearby_track_hit = track_from_vertex_hit_v.at( nearby_track_hit_iter );

	    // Compare three pieces of information to see if the hit is associated to the nearby track in the event.                                                                                   
	    if ( fabs( hit.StartTick() - nearby_track_hit->StartTick() ) < 0.001 && fabs( hit.EndTick() - nearby_track_hit->EndTick() ) < 0.001 && fabs( hit.Integral() - nearby_track_hit->Integral() ) < 0.001 ) {
	      is_associated_to_another_track_near_vertex = true;
	      break;
	    }

	  }

	}

	if ( is_associated_to_another_track_near_vertex == true ) continue;

	// Set the value of the 'vertex_location' variables equal to 0. before setting them equal to a nonzero value.
	vertex_location_x_with_SCE = 0.;
	vertex_location_y_with_SCE = 0.;
	vertex_location_z_with_SCE = 0.;

	if ( hit.View() == 0 ) {

	  vertex_location_x_with_SCE = u_plane_vertex_location_x_with_SCE;
	  vertex_location_y_with_SCE = u_plane_vertex_location_y_with_SCE;
	  vertex_location_z_with_SCE = u_plane_vertex_location_z_with_SCE;

	}

	if ( hit.View() == 1 ) {

	  vertex_location_x_with_SCE = v_plane_vertex_location_x_with_SCE;
	  vertex_location_y_with_SCE = v_plane_vertex_location_y_with_SCE;
	  vertex_location_z_with_SCE = v_plane_vertex_location_z_with_SCE;

	}

	if ( hit.View() == 2 ) {

	  vertex_location_x_with_SCE = y_plane_vertex_location_x_with_SCE;
	  vertex_location_y_with_SCE = y_plane_vertex_location_y_with_SCE;
	  vertex_location_z_with_SCE = y_plane_vertex_location_z_with_SCE;

	}
	
	double vertex_wire                  = geom->WireCoordinate( vertex_location_y_with_SCE, vertex_location_z_with_SCE, hit.View(), 0, 0);
	double vertex_collection_plane_tick = ( ( ( vertex_location_x_with_SCE + ( 0.111436 * flash_time ) ) / ( 0.5 * 0.111436 ) ) + 3200.  );
	double hit_tick                     = ( ( hit.StartTick() + hit.EndTick() ) / 2.0 ) + 2400.;
	double hit_wire                     = hit.WireID().Wire;
	double threshold                    = 200.0; // 296.923;   

	if ( hit.View() == 1 )  {
	  threshold = 200.0; //336.971                                                                                                                                                                        
	}

	if ( hit.View() ==  2 )  {
	  threshold = 200.0; // 294.793;                                                                                                                                                                      
	}

	if ( fabs( vertex_wire - hit_wire ) < wire_space_max_distance && fabs( vertex_collection_plane_tick - hit_tick ) < tick_space_max_distance && hit.Integral() > threshold ) {

	  if ( hit.View() == 0 )  {

	    u_plane_ticks_above_threshold.push_back( hit_tick );
	    u_plane_wires_above_threshold.push_back( hit_wire );
	    u_plane_charge_integrals_above_threshold.push_back( hit.Integral() );

	    num_hits_above_threshold_u_plane++;

	  }

	  if ( hit.View() == 1 ) {

	    v_plane_ticks_above_threshold.push_back( hit_tick );
	    v_plane_wires_above_threshold.push_back( hit_wire );
	    v_plane_charge_integrals_above_threshold.push_back( hit.Integral() );

	    num_hits_above_threshold_v_plane++;

	  }

	  if ( hit.View() == 2 ) {

	    y_plane_ticks_above_threshold.push_back( hit_tick );
	    y_plane_wires_above_threshold.push_back( hit_wire );
	    y_plane_charge_integrals_above_threshold.push_back( hit.Integral() );

	    num_hits_above_threshold_y_plane++;

	  }

	}

      }

    }

    // In order to get 3D points across the planes, you have to match the points.
    // This will be used in the proton reconstruction for hits located OFF the tracks.
    // U/V Plane Matching.                                                                                                                                                                              
    std::vector< int > u_plane_matched_tick_values_with_v_plane;
    std::vector< int > u_plane_matched_wire_values_with_v_plane;
    std::vector< int > v_plane_matched_tick_values_with_u_plane;
    std::vector< int > v_plane_matched_wire_values_with_u_plane;

    // Clear out the vectors.                                                                                                                                                                            
    u_plane_matched_tick_values_with_v_plane.clear();
    u_plane_matched_wire_values_with_v_plane.clear();
    v_plane_matched_tick_values_with_u_plane.clear();
    v_plane_matched_wire_values_with_u_plane.clear();

    // U/Y Plane Matching.                                                                                                                                                                               
    std::vector< int > u_plane_matched_tick_values_with_y_plane;
    std::vector< int > u_plane_matched_wire_values_with_y_plane;
    std::vector< int > y_plane_matched_tick_values_with_u_plane;
    std::vector< int > y_plane_matched_wire_values_with_u_plane;

    // Clear out the vectors.                                                                                                                                                                           
    u_plane_matched_tick_values_with_y_plane.clear();
    u_plane_matched_wire_values_with_y_plane.clear();
    y_plane_matched_tick_values_with_u_plane.clear();
    y_plane_matched_wire_values_with_u_plane.clear();

    // V/Y Plane Matching.                                                                                                                                                                               
    std::vector< int > v_plane_matched_tick_values_with_y_plane;
    std::vector< int > v_plane_matched_wire_values_with_y_plane;
    std::vector< int > y_plane_matched_tick_values_with_v_plane;
    std::vector< int > y_plane_matched_wire_values_with_v_plane;

    // Clear out the vectors.                                                                                                                                                                            
    v_plane_matched_tick_values_with_y_plane.clear();
    v_plane_matched_wire_values_with_y_plane.clear();
    y_plane_matched_tick_values_with_v_plane.clear();
    y_plane_matched_wire_values_with_v_plane.clear();

    
    double closest_tick_difference = 10000.;
    int    idx_of_closest_tick     = -1.; // This is so we only allow one-to-one matching across the planes.                                                                                            

    // These will be filled to only take the best two-plane match for a given point.                                                                                                                    
    std::vector< int > u_plane_points_matched_in_uv_loop;
    std::vector< int > u_plane_wires_matched_in_uv_loop;
    std::vector< int > u_plane_ticks_in_uv_loop;
    std::vector< int > u_plane_tick_difference_in_uv_loop;

    std::vector< int > u_plane_points_matched_in_uy_loop;
    std::vector< int > u_plane_wires_matched_in_uy_loop;
    std::vector< int > u_plane_ticks_in_uy_loop;
    std::vector< int > u_plane_tick_difference_in_uy_loop;

    std::vector< int > v_plane_points_matched_in_uv_loop;
    std::vector< int > v_plane_wires_matched_in_uv_loop;
    std::vector< int > v_plane_ticks_in_uv_loop;
    std::vector< int > v_plane_tick_difference_in_uv_loop;

    std::vector< int > v_plane_points_matched_in_vy_loop;
    std::vector< int > v_plane_wires_matched_in_vy_loop;
    std::vector< int > v_plane_ticks_in_vy_loop;
    std::vector< int > v_plane_tick_difference_in_vy_loop;

    std::vector< int > y_plane_points_matched_in_uy_loop;
    std::vector< int > y_plane_wires_matched_in_uy_loop;
    std::vector< int > y_plane_ticks_in_uy_loop;
    std::vector< int > y_plane_tick_difference_in_uy_loop;

    std::vector< int > y_plane_points_matched_in_vy_loop;
    std::vector< int > y_plane_wires_matched_in_vy_loop;
    std::vector< int > y_plane_ticks_in_vy_loop;
    std::vector< int > y_plane_tick_difference_in_vy_loop;

    u_plane_points_matched_in_uv_loop.clear();
    u_plane_wires_matched_in_uv_loop.clear();
    u_plane_ticks_in_uv_loop.clear();
    u_plane_tick_difference_in_uv_loop.clear();
    u_plane_points_matched_in_uy_loop.clear();
    u_plane_wires_matched_in_uy_loop.clear();
    u_plane_ticks_in_uy_loop.clear();
    u_plane_tick_difference_in_uy_loop.clear();
    v_plane_points_matched_in_uv_loop.clear();
    v_plane_wires_matched_in_uv_loop.clear();
    v_plane_ticks_in_uv_loop.clear();
    v_plane_tick_difference_in_uv_loop.clear();
    v_plane_points_matched_in_vy_loop.clear();
    v_plane_wires_matched_in_vy_loop.clear();
    v_plane_ticks_in_vy_loop.clear();
    v_plane_tick_difference_in_vy_loop.clear();
    y_plane_points_matched_in_uy_loop.clear();
    y_plane_wires_matched_in_uy_loop.clear();
    y_plane_ticks_in_uy_loop.clear();
    y_plane_tick_difference_in_uy_loop.clear();
    y_plane_points_matched_in_vy_loop.clear();
    y_plane_wires_matched_in_vy_loop.clear();
    y_plane_ticks_in_vy_loop.clear();
    y_plane_tick_difference_in_vy_loop.clear();

    if ( num_hits_above_threshold_u_plane > 0 && num_hits_above_threshold_v_plane > 0 ) {

      // Loop through the hits above threshold on the u-plane and find the one with the closest tick value on the v-plane.                                                                              
      for ( int u_plane_iter = 0; u_plane_iter < num_hits_above_threshold_u_plane; u_plane_iter++ ) {

	closest_tick_difference = 10000.;
	idx_of_closest_tick     = -1;

	// Start a loop over the v-plane points.                                                                                                                                                         
	for ( int v_plane_iter = 0; v_plane_iter < num_hits_above_threshold_v_plane; v_plane_iter++ ) {

	  bool already_matched_to_u_plane_hit = false;

	  for  ( int already_matched_iter = 0; already_matched_iter < int( v_plane_points_matched_in_uv_loop.size() ); already_matched_iter++ ) {

	    if ( v_plane_iter == v_plane_points_matched_in_uv_loop.at( already_matched_iter ) ) {
	      already_matched_to_u_plane_hit = true;
	      break;
	    }

	  }

	  if ( already_matched_to_u_plane_hit == true ) continue;

	  if ( fabs( u_plane_ticks_above_threshold.at( u_plane_iter ) - v_plane_ticks_above_threshold.at( v_plane_iter ) ) < closest_tick_difference ) {

	    closest_tick_difference = fabs( u_plane_ticks_above_threshold.at( u_plane_iter ) - v_plane_ticks_above_threshold.at( v_plane_iter ) );
	    idx_of_closest_tick     = v_plane_iter;

	  }

	}

	// Now see to make sure that an actual match was made.                                                                                                                                          
	if ( closest_tick_difference < 9998 ) {

	  u_plane_points_matched_in_uv_loop.push_back( u_plane_iter );
	  v_plane_points_matched_in_uv_loop.push_back( idx_of_closest_tick );
	  u_plane_wires_matched_in_uv_loop.push_back( u_plane_wires_above_threshold.at( u_plane_iter ) );
	  v_plane_wires_matched_in_uv_loop.push_back( v_plane_wires_above_threshold.at( idx_of_closest_tick ) );
	  u_plane_ticks_in_uv_loop.push_back( u_plane_ticks_above_threshold.at( u_plane_iter ) );
	  v_plane_ticks_in_uv_loop.push_back( v_plane_ticks_above_threshold.at( idx_of_closest_tick ) );
	  u_plane_tick_difference_in_uv_loop.push_back( closest_tick_difference );
	  v_plane_tick_difference_in_uv_loop.push_back( closest_tick_difference );

	}

      }

    }

    if ( num_hits_above_threshold_u_plane > 0 && num_hits_above_threshold_y_plane > 0 ) {

      // Loop through the hits above threshold on the u-plane and find the one with the closest tick value on the y-plane.                                                                              
      for ( int u_plane_iter = 0; u_plane_iter < num_hits_above_threshold_u_plane; u_plane_iter++ ) {

	closest_tick_difference = 10000.;
	idx_of_closest_tick     = -1;

	// Start a loop over the y-plane points.                                                                                                                                                        
	for ( int y_plane_iter = 0; y_plane_iter < num_hits_above_threshold_y_plane; y_plane_iter++ ) {

	  bool already_matched_to_u_plane_hit = false;

	  for  ( int already_matched_iter = 0; already_matched_iter < int( y_plane_points_matched_in_uy_loop.size() ); already_matched_iter++ ) {

	    if ( y_plane_iter == y_plane_points_matched_in_uy_loop.at( already_matched_iter ) ) {
	      already_matched_to_u_plane_hit = true;
	      break;
	    }

	  }

	  if ( already_matched_to_u_plane_hit == true ) continue;

	  if ( fabs( u_plane_ticks_above_threshold.at( u_plane_iter ) - y_plane_ticks_above_threshold.at( y_plane_iter ) ) < closest_tick_difference ) {

	    closest_tick_difference = fabs( u_plane_ticks_above_threshold.at( u_plane_iter ) - y_plane_ticks_above_threshold.at( y_plane_iter ) );
	    idx_of_closest_tick     = y_plane_iter;

	  }

	}

	// Now see to make sure that an actual match was made.                                                                                                                                         
	if ( closest_tick_difference < 9998 ) {

	  u_plane_points_matched_in_uy_loop.push_back( u_plane_iter );
	  y_plane_points_matched_in_uy_loop.push_back( idx_of_closest_tick );
	  u_plane_wires_matched_in_uy_loop.push_back( u_plane_wires_above_threshold.at( u_plane_iter ) );
	  y_plane_wires_matched_in_uy_loop.push_back( y_plane_wires_above_threshold.at( idx_of_closest_tick ) );
	  u_plane_ticks_in_uy_loop.push_back( u_plane_ticks_above_threshold.at( u_plane_iter ) );
	  y_plane_ticks_in_uy_loop.push_back( y_plane_ticks_above_threshold.at( idx_of_closest_tick ) );
	  u_plane_tick_difference_in_uy_loop.push_back( closest_tick_difference );
	  y_plane_tick_difference_in_uy_loop.push_back( closest_tick_difference );

	}
	
      }

    }

    if ( num_hits_above_threshold_v_plane > 0 && num_hits_above_threshold_y_plane > 0 ) {

      // Loop through the hits above threshold on the v-plane and find the one with the closest tick value on the y-plane.                                                                             
      for ( int v_plane_iter = 0; v_plane_iter < num_hits_above_threshold_v_plane; v_plane_iter++ ) {

	closest_tick_difference = 10000.;
	idx_of_closest_tick     = -1;

	// Start a loop over the y-plane points.                                                                                                                                                        
	for ( int y_plane_iter = 0; y_plane_iter < num_hits_above_threshold_y_plane; y_plane_iter++ ) {

	  bool already_matched_to_v_plane_hit = false;

	  for  ( int already_matched_iter = 0; already_matched_iter < int( y_plane_points_matched_in_vy_loop.size() ); already_matched_iter++ ) {

	    if ( y_plane_iter == y_plane_points_matched_in_vy_loop.at( already_matched_iter ) ) {
	      already_matched_to_v_plane_hit = true;
	      break;
	    }

	  }

	  if ( already_matched_to_v_plane_hit == true ) continue;

	  if ( fabs( v_plane_ticks_above_threshold.at( v_plane_iter ) - y_plane_ticks_above_threshold.at( y_plane_iter ) ) < closest_tick_difference ) {

	    closest_tick_difference = fabs( v_plane_ticks_above_threshold.at( v_plane_iter ) - y_plane_ticks_above_threshold.at( y_plane_iter ) );
	    idx_of_closest_tick     = y_plane_iter;

	  }

	}

	// Now see to make sure that an actual match was made.                                                                                                                                           
	if ( closest_tick_difference < 9998 ) {

	  v_plane_points_matched_in_vy_loop.push_back( v_plane_iter );
	  y_plane_points_matched_in_vy_loop.push_back( idx_of_closest_tick );
	  v_plane_wires_matched_in_vy_loop.push_back( v_plane_wires_above_threshold.at( v_plane_iter ) );
	  y_plane_wires_matched_in_vy_loop.push_back( y_plane_wires_above_threshold.at( idx_of_closest_tick ) );
	  v_plane_ticks_in_vy_loop.push_back( v_plane_ticks_above_threshold.at( v_plane_iter ) );
	  y_plane_ticks_in_vy_loop.push_back( y_plane_ticks_above_threshold.at( idx_of_closest_tick ) );
	  v_plane_tick_difference_in_vy_loop.push_back( closest_tick_difference );
	  y_plane_tick_difference_in_vy_loop.push_back( closest_tick_difference );

	}

      }

    }

    double x_coordinate = 0.;
    double y_coordinate = 0.;
    double z_coordinate = 0.;

    // Save the indices of the v-plane and y-plane points already used to make a match.                                                                                                                 
    std::vector< int > u_plane_indices_to_skip_in_matching_with_v_plane;
    std::vector< int > u_plane_indices_to_skip_in_matching_with_y_plane;
    std::vector< int > v_plane_indices_to_skip_in_matching_with_u_plane;
    std::vector< int > v_plane_indices_to_skip_in_matching_with_y_plane;
    std::vector< int > y_plane_indices_to_skip_in_matching_with_u_plane;
    std::vector< int > y_plane_indices_to_skip_in_matching_with_v_plane;

    u_plane_indices_to_skip_in_matching_with_v_plane.clear();
    u_plane_indices_to_skip_in_matching_with_y_plane.clear();
    v_plane_indices_to_skip_in_matching_with_u_plane.clear();
    v_plane_indices_to_skip_in_matching_with_y_plane.clear();
    y_plane_indices_to_skip_in_matching_with_u_plane.clear();
    y_plane_indices_to_skip_in_matching_with_v_plane.clear();

    // Now see if there are any duplicate matched hits and only save the one with the closer match in time.                                                                                             
    // This relies on the correspondence between the entry in the 'u_plane_points_matched_in_uv_loop' and the 'u_plane_tick_difference_in_uv_loop'.                                                     
    // This will ensure that each hit has no more than one match, although some hits will get neglected.                                                                                                 
    
    // U-Plane.                                                                                                                                                                                         
    for ( size_t v_plane_matched_iter = 0; v_plane_matched_iter < u_plane_points_matched_in_uv_loop.size(); v_plane_matched_iter++ ) {

      for ( size_t y_plane_matched_iter = 0; y_plane_matched_iter < u_plane_points_matched_in_uy_loop.size(); y_plane_matched_iter++ ) {

	if ( u_plane_points_matched_in_uv_loop.at( v_plane_matched_iter ) == u_plane_points_matched_in_uy_loop.at( y_plane_matched_iter ) ) {

	  if ( u_plane_tick_difference_in_uv_loop.at( v_plane_matched_iter ) < u_plane_tick_difference_in_uy_loop.at( y_plane_matched_iter ) ) {

	    u_plane_indices_to_skip_in_matching_with_y_plane.push_back( u_plane_points_matched_in_uy_loop.at( y_plane_matched_iter ) );

	  }

	  else {

	    u_plane_indices_to_skip_in_matching_with_v_plane.push_back( v_plane_points_matched_in_uv_loop.at( v_plane_matched_iter ) );

	  }

	}

      }

    } // End Case of the U-Plane.

    // V-Plane.                                                                                                                                                                                         
    for ( size_t u_plane_matched_iter = 0; u_plane_matched_iter < v_plane_points_matched_in_uv_loop.size(); u_plane_matched_iter++ ) {

      for ( size_t y_plane_matched_iter = 0; y_plane_matched_iter < v_plane_points_matched_in_vy_loop.size(); y_plane_matched_iter++ ) {

	if ( v_plane_points_matched_in_uv_loop.at( u_plane_matched_iter ) == v_plane_points_matched_in_vy_loop.at( y_plane_matched_iter ) ) {

	  if ( v_plane_tick_difference_in_uv_loop.at( u_plane_matched_iter ) < v_plane_tick_difference_in_vy_loop.at( y_plane_matched_iter ) ) {

	    v_plane_indices_to_skip_in_matching_with_y_plane.push_back( v_plane_points_matched_in_vy_loop.at( y_plane_matched_iter ) );

	  }

	  else {

	    v_plane_indices_to_skip_in_matching_with_u_plane.push_back( v_plane_points_matched_in_uv_loop.at( u_plane_matched_iter ) );

	  }

	}

      }

    } // End Case of the V-Plane. 

    // Y-Plane.                                                                                                                                                                                         
    for ( size_t u_plane_matched_iter = 0; u_plane_matched_iter < y_plane_points_matched_in_uy_loop.size(); u_plane_matched_iter++ ) {

      for ( size_t v_plane_matched_iter = 0; v_plane_matched_iter < y_plane_points_matched_in_vy_loop.size(); v_plane_matched_iter++ ) {

	if ( y_plane_points_matched_in_uy_loop.at( u_plane_matched_iter ) == y_plane_points_matched_in_vy_loop.at( v_plane_matched_iter ) ) {

	  if ( y_plane_tick_difference_in_uy_loop.at( u_plane_matched_iter ) < y_plane_tick_difference_in_vy_loop.at( v_plane_matched_iter ) ) {

	    y_plane_indices_to_skip_in_matching_with_v_plane.push_back( y_plane_points_matched_in_vy_loop.at( v_plane_matched_iter ) );

	  }

	  else {

	    y_plane_indices_to_skip_in_matching_with_u_plane.push_back( y_plane_points_matched_in_uy_loop.at( u_plane_matched_iter ) );

	  }

	}

      }

    } // End Case of the Y-Plane. 

    // Now, do the matching and make the 3D points.                                                                                                                                                     
    // U-V Plane matching.                                                                                                                                                                              
    for ( size_t u_plane_iter = 0; u_plane_iter < u_plane_points_matched_in_uv_loop.size(); u_plane_iter++ ) {

      // Make sure that this point isn't one of the ones that should be excluded.                                                                                                                        
      bool u_plane_point_should_be_excluded = false;

      for ( size_t u_plane_index_to_skip = 0; u_plane_index_to_skip < u_plane_indices_to_skip_in_matching_with_v_plane.size(); u_plane_index_to_skip++ ) {

	if ( u_plane_points_matched_in_uv_loop.at( u_plane_iter ) == u_plane_indices_to_skip_in_matching_with_v_plane.at( u_plane_index_to_skip ) ) {

	  u_plane_point_should_be_excluded = true;
	  break;
	}

      }

      if ( u_plane_point_should_be_excluded == true ) continue;

      int v_plane_idx = v_plane_points_matched_in_uv_loop.at( u_plane_iter );

      bool v_plane_point_should_be_excluded = false;

      for ( size_t v_plane_index_to_skip = 0; v_plane_index_to_skip < v_plane_indices_to_skip_in_matching_with_u_plane.size(); v_plane_index_to_skip++ ) {

	if ( v_plane_idx == v_plane_indices_to_skip_in_matching_with_u_plane.at( v_plane_index_to_skip ) ) {
	  v_plane_point_should_be_excluded = true;
	  break;
	}

      }

      if ( v_plane_point_should_be_excluded == true ) continue;

      // Now, calculate the coordinates.                                                                                                                                                                 
      x_coordinate = ( ( ( ( u_plane_ticks_in_uv_loop.at( u_plane_iter ) + v_plane_ticks_in_uv_loop.at( u_plane_iter ) ) / 2.0 ) - 3200. ) * ( 0.5 ) * ( 0.111436 ) );

      // Call the 'intersection_point' function to get the y and z coordinates of the crossing point.                                                                                                    
      geo::WireID wid_u_plane(0, 0, 0, u_plane_wires_matched_in_uv_loop.at( u_plane_iter ));
      geo::WireID wid_v_plane(0, 0, 1, v_plane_wires_matched_in_uv_loop.at( u_plane_iter ));;

      geo::WireIDIntersection wire_intersection;

      bool do_wires_intersect = geom->WireIDsIntersect( wid_u_plane, wid_v_plane, wire_intersection );

      // Keep this print statement to avoid an error on unused variables.
      std::cout << "The verdict on if the wires intersect = " << do_wires_intersect << "." << std::endl;

      y_coordinate                                           = wire_intersection.y;
      z_coordinate                                           = wire_intersection.z;

      // Append all of the coordinates to the vector of the 3D points corresponding to the proton track.                                                                                                
      x_points_above_threshold[ num_points_above_threshold ] = x_coordinate;
      y_points_above_threshold[ num_points_above_threshold ] = y_coordinate;
      z_points_above_threshold[ num_points_above_threshold ] = z_coordinate;

      num_points_above_threshold++;

    }

    // U-Y Plane Matching.                                                                                                                                                                               
    for ( size_t u_plane_iter = 0; u_plane_iter < u_plane_points_matched_in_uy_loop.size(); u_plane_iter++ ) {

      // Make sure that this point isn't one of the ones that should be excluded.                                                                                                                        
      bool u_plane_point_should_be_excluded = false;

      for ( size_t u_plane_index_to_skip = 0; u_plane_index_to_skip < u_plane_indices_to_skip_in_matching_with_y_plane.size(); u_plane_index_to_skip++ ) {

	if ( u_plane_points_matched_in_uy_loop.at( u_plane_iter ) == u_plane_indices_to_skip_in_matching_with_y_plane.at( u_plane_index_to_skip ) ) {
	  u_plane_point_should_be_excluded = true;
	  break;
	}

      }

      if ( u_plane_point_should_be_excluded == true ) continue;

      int y_plane_idx = y_plane_points_matched_in_uy_loop.at( u_plane_iter );

      // Make sure that this point isn't one of the ones that should be excluded.                                                                                                                       
      bool y_plane_point_should_be_excluded = false;

      for ( size_t y_plane_index_to_skip = 0; y_plane_index_to_skip < y_plane_indices_to_skip_in_matching_with_u_plane.size(); y_plane_index_to_skip++ ) {

	if ( y_plane_idx == y_plane_indices_to_skip_in_matching_with_u_plane.at( y_plane_index_to_skip ) ) {
	  y_plane_point_should_be_excluded = true;
	  break;
	}

      }

      if ( y_plane_point_should_be_excluded == true ) continue;
      
      // Now, calculate the coordinates.                                                                                                                                                                
      x_coordinate = ( ( ( ( u_plane_ticks_in_uy_loop.at( u_plane_iter ) + y_plane_ticks_in_uy_loop.at( u_plane_iter ) ) / 2.0 ) - 3200. ) * ( 0.5 ) * ( 0.111436 ) );

      // Call the 'intersection_point' function to get the y and z coordinates of the crossing point.                                                                                                  
      geo::WireID wid_u_plane(0, 0, 0, u_plane_wires_matched_in_uy_loop.at( u_plane_iter ));
      geo::WireID wid_y_plane(0, 0, 2, y_plane_wires_matched_in_uy_loop.at( u_plane_iter ));;

      geo::WireIDIntersection wire_intersection;

      bool do_wires_intersect = geom->WireIDsIntersect( wid_u_plane, wid_y_plane, wire_intersection );

      // Keep this print statement to avoid an error on unused variables. 
      std::cout << "The verdict on if the wires intersect = " << do_wires_intersect << "." << std::endl;

      y_coordinate                                           = wire_intersection.y;
      z_coordinate                                           = wire_intersection.z;

      // Append all of the coordinates to the vector of the 3D points corresponding to the proton track.                                                                                                
      x_points_above_threshold[ num_points_above_threshold ] = x_coordinate;
      y_points_above_threshold[ num_points_above_threshold ] = y_coordinate;
      z_points_above_threshold[ num_points_above_threshold ] = z_coordinate;

      num_points_above_threshold++;

    }

    // V-Y Plane Matching.                                                                                                                                                                              
    for ( size_t v_plane_iter = 0; v_plane_iter < v_plane_points_matched_in_vy_loop.size(); v_plane_iter++ ) {

      // Make sure that this point isn't one of the ones that should be excluded.                                                                                                                        
      bool v_plane_point_should_be_excluded = false;

      for ( size_t v_plane_index_to_skip = 0; v_plane_index_to_skip < v_plane_indices_to_skip_in_matching_with_y_plane.size(); v_plane_index_to_skip++ ) {

	if ( v_plane_points_matched_in_vy_loop.at( v_plane_iter ) == v_plane_indices_to_skip_in_matching_with_y_plane.at( v_plane_index_to_skip ) ) {
	  v_plane_point_should_be_excluded = true;
	  break;
	}

      }

      if ( v_plane_point_should_be_excluded == true ) continue;

      int y_plane_idx = y_plane_points_matched_in_vy_loop.at( v_plane_iter );

      // Make sure that this point isn't one of the ones that should be excluded.                                                                                                                        
      bool y_plane_point_should_be_excluded = false;

      for ( size_t y_plane_index_to_skip = 0; y_plane_index_to_skip < y_plane_indices_to_skip_in_matching_with_v_plane.size(); y_plane_index_to_skip++ ) {

	if ( y_plane_idx == y_plane_indices_to_skip_in_matching_with_v_plane.at( y_plane_index_to_skip ) ) {
	  y_plane_point_should_be_excluded = true;
	  break;
	}

      }

      if ( y_plane_point_should_be_excluded == true ) continue;

      // Now, calculate the coordinates.                                                                                                                                                              
      x_coordinate = ( ( ( ( v_plane_ticks_in_vy_loop.at( v_plane_iter ) + y_plane_ticks_in_vy_loop.at( v_plane_iter ) ) / 2.0 ) - 3200. ) * ( 0.5 ) * ( 0.111436 ) );

      // Call the 'intersection_point' function to get the y and z coordinates of the crossing point.                                                                                                    
      geo::WireID wid_v_plane(0, 0, 1, v_plane_wires_matched_in_vy_loop.at( v_plane_iter ));
      geo::WireID wid_y_plane(0, 0, 2, y_plane_wires_matched_in_vy_loop.at( v_plane_iter ));;

      geo::WireIDIntersection wire_intersection;

      bool do_wires_intersect = geom->WireIDsIntersect( wid_v_plane, wid_y_plane, wire_intersection );

      // Keep this print statement to avoid an error on unused variables. 
      std::cout << "The verdict on if the wires intersect = " << do_wires_intersect << "." << std::endl;

      y_coordinate                                           = wire_intersection.y;
      z_coordinate                                           = wire_intersection.z;

      // Append all of the coordinates to the vector of the 3D points corresponding to the proton track.                                                                                                
      x_points_above_threshold[ num_points_above_threshold ] = x_coordinate;
      y_points_above_threshold[ num_points_above_threshold ] = y_coordinate;
      z_points_above_threshold[ num_points_above_threshold ] = z_coordinate;

      num_points_above_threshold++;

    }

    double u_plane_vertex_wire                               = geom->WireCoordinate( u_plane_vertex_location_y, u_plane_vertex_location_z, 0, 0, 0 );
    double u_plane_vertex_time_tick                          = ( ( ( u_plane_vertex_location_x + ( 0.111436 * flash_time ) ) / ( 0.5 * 0.111436 ) ) + 3200.  );
    
    double v_plane_vertex_wire                               = geom->WireCoordinate( v_plane_vertex_location_y, v_plane_vertex_location_z, 1, 0, 0 );
    double v_plane_vertex_time_tick                          = ( ( ( v_plane_vertex_location_x + ( 0.111436 * flash_time ) ) / ( 0.5 * 0.111436 ) ) + 3200.  );

    double y_plane_vertex_wire                               = geom->WireCoordinate( y_plane_vertex_location_y, y_plane_vertex_location_z, 2, 0, 0 );
    double y_plane_vertex_time_tick                          = ( ( ( y_plane_vertex_location_x + ( 0.111436 * flash_time ) ) / ( 0.5 * 0.111436 ) ) + 3200.  );

    // Loop through the hits to find out how many are within 2 maximum track lengths (133 wires, 718 ticks) of the TPC object.
    for ( size_t event_hit_iter = 0; event_hit_iter < hit_h->size(); event_hit_iter++ ) {

      auto hit = hit_h->at( event_hit_iter );

      // Put the continue statement in for the 'bad hits'.                                                                                                                                                
      float rms = hit.RMS();

      // Filter on the peaks.                                                                                                                                                                             
      if ( fabs( rms - 0.333 ) < 0.0001 || fabs( rms - 0.5 ) < 0.0001 || fabs( rms - 0.667 ) < 0.0001 || fabs( rms - 1.0 ) < 0.0001 || fabs( rms - 1.333 ) < 0.0001 || fabs( rms - 1.667 ) < 0.0001 || fabs( rms - 2.0 ) < 0.0001 || fabs( rms - 2.333 ) < 0.0001 || fabs( rms - 2.667 ) < 0.0001 || fabs( rms - 2.8572 ) < 0.0001 || fabs( rms - 3.0 ) < 0.0001 || fabs( rms - 3.333 ) < 0.0001 || fabs( rms - 3.667 ) < 0.0001 || fabs( rms - 4.0 ) < 0.0001 || fabs( rms - 4.333 ) < 0.0001 || fabs( rms - 4.667 ) < 0.0001 || fabs( rms - 5.0 ) < 0.0001 || fabs( rms - 5.333 ) < 0.0001 || fabs( rms - 5.667 ) < 0.0001 || fabs( rms - 6.0 ) < 0.0001 || fabs( rms - 6.333 ) < 0.0001 || fabs( rms - 6.5 ) < 0.0001 || fabs( rms - 7.0 ) < 0.0001 || fabs( rms - 7.333 ) < 0.0001 || fabs( rms - 7.5 ) < 0.0001 || fabs( rms - 7.667 ) < 0.0001 || fabs( rms - 8.0 ) < 0.0001 || fabs( rms - 8.5 ) < 0.0001 || fabs( rms - 9.0 ) < 0.0001 || fabs( rms - 9.5 ) < 0.0001 || fabs( rms - 10.0 ) < 0.0001 || fabs( rms - 10.5 ) < 0.0001 || fabs( rms - 11.0 ) < 0.0001 || fabs( rms - 11.5 ) < 0.0001 || fabs( rms - 12.0 ) < 0.0001 || fabs( rms - 12.5 ) < 0.0001 || fabs( rms - 13.0 ) < 0.0001 || fabs( rms - 13.5 ) < 0.0001 || fabs( rms - 14.0 ) < 0.0001 || fabs( rms - 14.5 ) < 0.0001 || fabs( rms - 15.0 ) < 0.0001 || fabs( rms - 15.5 ) < 0.0001 || fabs( rms - 16.0 ) < 0.0001 || fabs( rms - 16.5 ) < 0.0001 || fabs( rms - 17.0 ) < 0.0001 || fabs( rms - 17.5 ) < 0.0001 || fabs( rms - 18.0 ) < 0.0001 || fabs( rms - 18.5 ) < 0.0001 || fabs( rms - 19.0 ) < 0.0001 || fabs( rms - 19.5 ) < 0.0001|| fabs( rms - 20.0 ) < 0.0001 )
	continue;

      double vertex_wire = 0.;
      double vertex_tick = 0.;

      if ( hit.View() == 0 ) {

	vertex_wire = u_plane_vertex_wire;
	vertex_tick = u_plane_vertex_time_tick;

      }

      if ( hit.View() == 1 ) {

	vertex_wire = v_plane_vertex_wire;
	vertex_tick = v_plane_vertex_time_tick;

      } 

      if ( hit.View() == 2 ) {

	vertex_wire = y_plane_vertex_wire;
	vertex_tick = y_plane_vertex_time_tick;

      } 

      double hit_tick                     = ( ( hit.StartTick() + hit.EndTick() ) / 2.0 ) + 2400.;
      double hit_wire                     = hit.WireID().Wire;

      if ( fabs( hit_wire - vertex_wire ) < 133 && fabs( hit_tick - vertex_tick ) < 718 ) {

	// Increment the total number of hits in the slice.
	total_number_of_slice_hits_in_vicinity++;

	if ( hit.View() == 0 ) {

	  u_plane_number_of_slice_hits_in_vicinity++;

	  u_plane_sum_of_slice_hits_in_vicinity_ADCs += hit.Integral();
	  total_sum_of_slice_hits_in_vicinity_ADCs   += hit.Integral();

	}

	if ( hit.View() == 1 ) {

	  v_plane_number_of_slice_hits_in_vicinity++;

	  v_plane_sum_of_slice_hits_in_vicinity_ADCs += hit.Integral();
          total_sum_of_slice_hits_in_vicinity_ADCs   += hit.Integral();
	  
	}

	if ( hit.View() == 2 ) {

	  y_plane_number_of_slice_hits_in_vicinity++;

	  y_plane_sum_of_slice_hits_in_vicinity_ADCs += hit.Integral();
          total_sum_of_slice_hits_in_vicinity_ADCs   += hit.Integral();

        }

      } // End of the requirement that the hits are in the vicinity of the vertex. 

    } // End of the loop over the hits.

    // Work on orienting the proton track using the information that we have available from the hits.
    proton_x_direction = 1;
    proton_y_direction = 1;
    proton_z_direction = 1;

    // Find the closest hit above threshold within 3 wires from the first wire.
    // Only look at hits that are further down the muon track from the vertex hit (if they are on the track at all ).
    double u_plane_closest_hit_wire_above_threshold_to_vertex_hit = 9999.;
    double u_plane_closest_hit_tick_above_threshold_to_vertex_hit = 9999.;
    int    num_of_u_plane_hits_considered                         = 0;

    double v_plane_closest_hit_wire_above_threshold_to_vertex_hit = 9999.;
    double v_plane_closest_hit_tick_above_threshold_to_vertex_hit = 9999.;
    int   num_of_v_plane_hits_considered                          = 0;
    
    double y_plane_closest_hit_wire_above_threshold_to_vertex_hit = 9999.;
    double y_plane_closest_hit_tick_above_threshold_to_vertex_hit = 9999.;
    int    num_of_y_plane_hits_considered                         = 0;
	 
    for ( size_t vertex_hit_iter = 0; vertex_hit_iter < hit_h->size(); vertex_hit_iter++ ) {

      auto hit = hit_h->at( vertex_hit_iter );

      // Put the continue statement in for the 'bad hits'.                                                                                                                                                
      float rms = hit.RMS();

      // Filter on the peaks.                                                                                                                                                                             
      if ( fabs( rms - 0.333 ) < 0.0001 || fabs( rms - 0.5 ) < 0.0001 || fabs( rms - 0.667 ) < 0.0001 || fabs( rms - 1.0 ) < 0.0001 || fabs( rms - 1.333 ) < 0.0001 || fabs( rms - 1.667 ) < 0.0001 || fabs( rms - 2.0 ) < 0.0001 || fabs( rms - 2.333 ) < 0.0001 || fabs( rms - 2.667 ) < 0.0001 || fabs( rms - 2.8572 ) < 0.0001 || fabs( rms - 3.0 ) < 0.0001 || fabs( rms - 3.333 ) < 0.0001 || fabs( rms - 3.667 ) < 0.0001 || fabs( rms - 4.0 ) < 0.0001 || fabs( rms - 4.333 ) < 0.0001 || fabs( rms - 4.667 ) < 0.0001 || fabs( rms - 5.0 ) < 0.0001 || fabs( rms - 5.333 ) < 0.0001 || fabs( rms - 5.667 ) < 0.0001 || fabs( rms - 6.0 ) < 0.0001 || fabs( rms - 6.333 ) < 0.0001 || fabs( rms - 6.5 ) < 0.0001 || fabs( rms - 7.0 ) < 0.0001 || fabs( rms - 7.333 ) < 0.0001 || fabs( rms - 7.5 ) < 0.0001 || fabs( rms - 7.667 ) < 0.0001 || fabs( rms - 8.0 ) < 0.0001 || fabs( rms - 8.5 ) < 0.0001 || fabs( rms - 9.0 ) < 0.0001 || fabs( rms - 9.5 ) < 0.0001 || fabs( rms - 10.0 ) < 0.0001 || fabs( rms - 10.5 ) < 0.0001 || fabs( rms - 11.0 ) < 0.0001 || fabs( rms - 11.5 ) < 0.0001 || fabs( rms - 12.0 ) < 0.0001 || fabs( rms - 12.5 ) < 0.0001 || fabs( rms - 13.0 ) < 0.0001 || fabs( rms - 13.5 ) < 0.0001 || fabs( rms - 14.0 ) < 0.0001 || fabs( rms - 14.5 ) < 0.0001 || fabs( rms - 15.0 ) < 0.0001 || fabs( rms - 15.5 ) < 0.0001 || fabs( rms - 16.0 ) < 0.0001 || fabs( rms - 16.5 ) < 0.0001 || fabs( rms - 17.0 ) < 0.0001 || fabs( rms - 17.5 ) < 0.0001 || fabs( rms - 18.0 ) < 0.0001 || fabs( rms - 18.5 ) < 0.0001 || fabs( rms - 19.0 ) < 0.0001 || fabs( rms - 19.5 ) < 0.0001|| fabs( rms - 20.0 ) < 0.0001 )
	continue;

      if ( hit.View() == 0 && no_calo_points_on_u_plane == false ) {

	// Do not consider the vertex hit.
	if ( fabs( hit.WireID().Wire - u_plane_wire_at_vertex ) < 0.001 && fabs( hit.Integral() - u_plane_wire_at_vertex_integral ) < 0.001 )
	  continue;

	if ( fabs( hit.WireID().Wire - u_plane_wire_at_vertex ) < 3.0 && hit.Integral() > 200. ) {

	  // Make sure that this is not closer to the end of the track in the cheap way.
	  if ( ( vertex_location_x > other_end_location_x && ( ( ( hit.StartTick() + hit.EndTick() ) / 2.0 ) + 2400. ) > u_plane_vertex_tick ) || ( vertex_location_x < other_end_location_x && ( ( ( hit.StartTick() + hit.EndTick() ) / 2.0 ) + 2400. ) < u_plane_vertex_tick ) ) {

	    num_of_u_plane_hits_considered++;

	    if ( fabs( hit.WireID().Wire - u_plane_wire_at_vertex ) < fabs( u_plane_closest_hit_wire_above_threshold_to_vertex_hit - u_plane_wire_at_vertex ) ) {

	      u_plane_closest_hit_wire_above_threshold_to_vertex_hit = hit.WireID().Wire;
	      u_plane_closest_hit_tick_above_threshold_to_vertex_hit = ( ( ( hit.StartTick() + hit.EndTick() ) / 2.0 ) + 2400. );

	    }

	  }

	} 

      } // End of finding the U-Plane hits above threshold.
      
      if ( hit.View() == 1 && no_calo_points_on_v_plane == false ) {

	// Do not consider the vertex hit.                                                                                                                                                             
        if ( fabs( hit.WireID().Wire - v_plane_wire_at_vertex ) < 0.001 && fabs( hit.Integral() - v_plane_wire_at_vertex_integral ) < 0.001 )
	  continue;

	if ( fabs( hit.WireID().Wire - v_plane_wire_at_vertex ) < 3.0 && hit.Integral() > 200. ) {

	  // Make sure that this is not closer to the end of the track in the cheap way.                                     
	  if ( ( vertex_location_x > other_end_location_x && ( ( ( hit.StartTick() + hit.EndTick() ) / 2.0 ) + 2400. ) > v_plane_vertex_tick ) || ( vertex_location_x < other_end_location_x && ( ( ( hit.StartTick() + hit.EndTick() ) / 2.0 ) + 2400. ) < v_plane_vertex_tick ) ) {

	    num_of_v_plane_hits_considered++;
	    
	    if ( fabs( hit.WireID().Wire - v_plane_wire_at_vertex ) < fabs( v_plane_closest_hit_wire_above_threshold_to_vertex_hit - v_plane_wire_at_vertex ) ) {
	      
	      v_plane_closest_hit_wire_above_threshold_to_vertex_hit = hit.WireID().Wire;
	      v_plane_closest_hit_tick_above_threshold_to_vertex_hit = ( ( ( hit.StartTick() + hit.EndTick() ) / 2.0 ) + 2400. );
		
	    }
	    
	  }
	  
	}
	
      } // End of finding the V-Plane hits above threshold. 
      
      if ( hit.View() == 2 && no_calo_points_on_y_plane == false ) {

	// Do not consider the vertex hit.                                                                                                                                                             
        if ( fabs( hit.WireID().Wire - y_plane_wire_at_vertex ) < 0.001 && fabs( hit.Integral() - y_plane_wire_at_vertex_integral ) < 0.001 )
	  continue;
	
	if ( fabs( hit.WireID().Wire - y_plane_wire_at_vertex ) < 3.0 && hit.Integral() > 200. ) {
	   
	  // Make sure that this is not closer to the end of the track in the cheap way.                                                                                                       
	  if ( ( vertex_location_x > other_end_location_x && ( ( ( hit.StartTick() + hit.EndTick() ) / 2.0 ) + 2400. ) > y_plane_vertex_tick ) || ( vertex_location_x < other_end_location_x && ( ( ( hit.StartTick() + hit.EndTick() ) / 2.0 ) + 2400. ) < y_plane_vertex_tick ) ) {

	    num_of_y_plane_hits_considered++;

	    if ( fabs( hit.WireID().Wire - y_plane_wire_at_vertex ) < fabs( y_plane_closest_hit_wire_above_threshold_to_vertex_hit - y_plane_wire_at_vertex ) ) {
	      
	      y_plane_closest_hit_wire_above_threshold_to_vertex_hit = hit.WireID().Wire;
	      y_plane_closest_hit_tick_above_threshold_to_vertex_hit = ( ( ( hit.StartTick() + hit.EndTick() ) / 2.0 ) + 2400. );
	      
	    }
	      
	  }
	  
	}  
	
      } // End of finding the Y-Plane hits above threshold. 

    }
    
    // Find the direction in wire space in which the high-charge hits emanate from the vertex.
    double u_plane_wire_direction = ( u_plane_closest_hit_wire_above_threshold_to_vertex_hit - u_plane_wire_at_vertex );
    double v_plane_wire_direction = ( v_plane_closest_hit_wire_above_threshold_to_vertex_hit - v_plane_wire_at_vertex );
    double y_plane_wire_direction = ( y_plane_closest_hit_wire_above_threshold_to_vertex_hit - y_plane_wire_at_vertex );

    // Look for other hits in the vicinity of the vertex hit, but ONLY in the direction of the closest hit to the vertex.
    std::vector< double > u_plane_nearby_hits_tick;
    std::vector< double > u_plane_nearby_hits_wire;
    u_plane_nearby_hits_tick.clear();
    u_plane_nearby_hits_wire.clear();

    std::vector< double > v_plane_nearby_hits_tick;
    std::vector< double > v_plane_nearby_hits_wire;
    v_plane_nearby_hits_tick.clear();
    v_plane_nearby_hits_wire.clear();

    std::vector< double > y_plane_nearby_hits_tick;
    std::vector< double > y_plane_nearby_hits_wire;
    y_plane_nearby_hits_tick.clear();
    y_plane_nearby_hits_wire.clear();
    
    if ( num_of_u_plane_hits_considered > 0 ) {

      u_plane_nearby_hits_tick.push_back( u_plane_closest_hit_tick_above_threshold_to_vertex_hit );
      u_plane_nearby_hits_wire.push_back( u_plane_closest_hit_wire_above_threshold_to_vertex_hit );

    }

    if ( num_of_v_plane_hits_considered > 0 ) {

      v_plane_nearby_hits_tick.push_back( v_plane_closest_hit_tick_above_threshold_to_vertex_hit );
      v_plane_nearby_hits_wire.push_back( v_plane_closest_hit_wire_above_threshold_to_vertex_hit );

    } 

    if ( num_of_y_plane_hits_considered > 0 ) {

      y_plane_nearby_hits_tick.push_back( y_plane_closest_hit_tick_above_threshold_to_vertex_hit );
      y_plane_nearby_hits_wire.push_back( y_plane_closest_hit_wire_above_threshold_to_vertex_hit );

    }

    //  Includes the vertex and the closest hit.
    int total_num_u_plane_points = 0;
    int total_num_v_plane_points = 0;
    int total_num_y_plane_points = 0;

    if ( no_calo_points_on_u_plane == false ) total_num_u_plane_points = 1;
    if ( no_calo_points_on_v_plane == false ) total_num_v_plane_points = 1;
    if ( no_calo_points_on_y_plane == false ) total_num_y_plane_points = 1;

    if ( num_of_u_plane_hits_considered > 0 ) 
      total_num_u_plane_points++;

    if ( num_of_v_plane_hits_considered > 0 ) 
      total_num_v_plane_points++;

    if ( num_of_y_plane_hits_considered > 0 )
      total_num_y_plane_points++;

    // Now, look for nearby hits in the direction of the closest hit.
    // Loop through hits and find the ones within five wires from the vertex hit.
    for ( size_t nearby_hit_iter = 0; nearby_hit_iter < hit_h->size(); nearby_hit_iter++ ) {

      auto hit = hit_h->at( nearby_hit_iter );

      // Put the continue statement in for the 'bad hits'.                                                                                                                                                
      float rms = hit.RMS();

      // Filter on the peaks.                                                                                                                                                                             
      if ( fabs( rms - 0.333 ) < 0.0001 || fabs( rms - 0.5 ) < 0.0001 || fabs( rms - 0.667 ) < 0.0001 || fabs( rms - 1.0 ) < 0.0001 || fabs( rms - 1.333 ) < 0.0001 || fabs( rms - 1.667 ) < 0.0001 || fabs( rms - 2.0 ) < 0.0001 || fabs( rms - 2.333 ) < 0.0001 || fabs( rms - 2.667 ) < 0.0001 || fabs( rms - 2.8572 ) < 0.0001 || fabs( rms - 3.0 ) < 0.0001 || fabs( rms - 3.333 ) < 0.0001 || fabs( rms - 3.667 ) < 0.0001 || fabs( rms - 4.0 ) < 0.0001 || fabs( rms - 4.333 ) < 0.0001 || fabs( rms - 4.667 ) < 0.0001 || fabs( rms - 5.0 ) < 0.0001 || fabs( rms - 5.333 ) < 0.0001 || fabs( rms - 5.667 ) < 0.0001 || fabs( rms - 6.0 ) < 0.0001 || fabs( rms - 6.333 ) < 0.0001 || fabs( rms - 6.5 ) < 0.0001 || fabs( rms - 7.0 ) < 0.0001 || fabs( rms - 7.333 ) < 0.0001 || fabs( rms - 7.5 ) < 0.0001 || fabs( rms - 7.667 ) < 0.0001 || fabs( rms - 8.0 ) < 0.0001 || fabs( rms - 8.5 ) < 0.0001 || fabs( rms - 9.0 ) < 0.0001 || fabs( rms - 9.5 ) < 0.0001 || fabs( rms - 10.0 ) < 0.0001 || fabs( rms - 10.5 ) < 0.0001 || fabs( rms - 11.0 ) < 0.0001 || fabs( rms - 11.5 ) < 0.0001 || fabs( rms - 12.0 ) < 0.0001 || fabs( rms - 12.5 ) < 0.0001 || fabs( rms - 13.0 ) < 0.0001 || fabs( rms - 13.5 ) < 0.0001 || fabs( rms - 14.0 ) < 0.0001 || fabs( rms - 14.5 ) < 0.0001 || fabs( rms - 15.0 ) < 0.0001 || fabs( rms - 15.5 ) < 0.0001 || fabs( rms - 16.0 ) < 0.0001 || fabs( rms - 16.5 ) < 0.0001 || fabs( rms - 17.0 ) < 0.0001 || fabs( rms - 17.5 ) < 0.0001 || fabs( rms - 18.0 ) < 0.0001 || fabs( rms - 18.5 ) < 0.0001 || fabs( rms - 19.0 ) < 0.0001 || fabs( rms - 19.5 ) < 0.0001|| fabs( rms - 20.0 ) < 0.0001 )
	continue;

      if ( u_plane_nearby_hits_wire.size() > 0 && hit.View() == 0 ) {

	if ( fabs( hit.WireID().Wire - u_plane_wire_at_vertex ) < 5.0 ) {

	  if ( ( u_plane_wire_direction > 0 && ( hit.WireID().Wire - u_plane_wire_at_vertex ) > 0 ) || ( u_plane_wire_direction < 0 && ( hit.WireID().Wire - u_plane_wire_at_vertex ) < 0 ) ) {
	  
	    u_plane_nearby_hits_tick.push_back( ( ( hit.StartTick() + hit.EndTick() ) / 2.0 ) + 2400. );
	    u_plane_nearby_hits_wire.push_back( hit.WireID().Wire );
	    
	    total_num_u_plane_points++;
	    
	  }
	    
	}

      } // End of finding U-Plane nearby hits above threshold.

      if ( v_plane_nearby_hits_wire.size() > 0 && hit.View() == 1 ) {

        if ( fabs( hit.WireID().Wire - v_plane_wire_at_vertex ) < 5.0 ) {

	  if ( ( v_plane_wire_direction > 0 && ( hit.WireID().Wire - v_plane_wire_at_vertex ) > 0 ) || ( v_plane_wire_direction < 0 && ( hit.WireID().Wire - v_plane_wire_at_vertex ) < 0 ) ) {

	    v_plane_nearby_hits_tick.push_back( ( ( hit.StartTick() + hit.EndTick() ) / 2.0 ) + 2400. );
	    v_plane_nearby_hits_wire.push_back( hit.WireID().Wire );
	    
	    total_num_v_plane_points++;

	  }

        }

      } // End of finding V-Plane nearby hits above threshold. 

      if ( y_plane_nearby_hits_wire.size() > 0 && hit.View() == 2 ) {

        if ( fabs( hit.WireID().Wire - y_plane_wire_at_vertex ) < 5.0 ) {

	  if ( ( y_plane_wire_direction > 0 && ( hit.WireID().Wire - y_plane_wire_at_vertex ) > 0 ) || ( y_plane_wire_direction < 0 && ( hit.WireID().Wire - y_plane_wire_at_vertex ) < 0 ) ) {

	    y_plane_nearby_hits_tick.push_back( (  ( hit.StartTick() + hit.EndTick() ) / 2.0) + 2400. );
	    y_plane_nearby_hits_wire.push_back( hit.WireID().Wire );
	    
	    total_num_y_plane_points++;

	  }

        }

      } // End of finding Y-Plane nearby hits above threshold.  

    }

    // For each plane, go through and find the wire furthest away from the vertex wire.
    size_t furthest_away_u_plane_wire_idx = -1;
    size_t furthest_away_v_plane_wire_idx = -1;
    size_t furthest_away_y_plane_wire_idx = -1;

    double furthest_u_plane_wire_dist     = -1;
    double furthest_v_plane_wire_dist     = -1;
    double furthest_y_plane_wire_dist     = -1;

    if (  u_plane_nearby_hits_wire.size() > 0 ) {

      for ( size_t u_plane_point_iter = 0; u_plane_point_iter < u_plane_nearby_hits_wire.size(); u_plane_point_iter++ ) {

	if ( fabs( u_plane_nearby_hits_wire.at( u_plane_point_iter ) - u_plane_wire_at_vertex ) > furthest_u_plane_wire_dist )  {

	  furthest_away_u_plane_wire_idx = u_plane_point_iter;
	  furthest_u_plane_wire_dist     = fabs( u_plane_nearby_hits_wire.at( u_plane_point_iter ) - u_plane_wire_at_vertex );

	}

      }

    } // End of requirement that there is at least one nearby U-Plane hit above threshold.

    if (  v_plane_nearby_hits_wire.size() > 0 ) {

      for ( size_t v_plane_point_iter = 0; v_plane_point_iter < v_plane_nearby_hits_wire.size(); v_plane_point_iter++ ) {

        if ( fabs( v_plane_nearby_hits_wire.at( v_plane_point_iter ) - v_plane_wire_at_vertex ) > furthest_v_plane_wire_dist )  {

          furthest_away_v_plane_wire_idx = v_plane_point_iter;
          furthest_v_plane_wire_dist     = fabs( v_plane_nearby_hits_wire.at( v_plane_point_iter ) - v_plane_wire_at_vertex );

        }

      } 

    }  // End of requirement that there is at least one nearby V-Plane hit above threshold. 

    if (  y_plane_nearby_hits_wire.size() > 0 ) {

      for ( size_t y_plane_point_iter = 0; y_plane_point_iter < y_plane_nearby_hits_wire.size(); y_plane_point_iter++ ) {

        if ( fabs( y_plane_nearby_hits_wire.at( y_plane_point_iter ) - y_plane_wire_at_vertex ) > furthest_y_plane_wire_dist )  {

          furthest_away_y_plane_wire_idx = y_plane_point_iter;
          furthest_y_plane_wire_dist     = fabs( y_plane_nearby_hits_wire.at( y_plane_point_iter ) - y_plane_wire_at_vertex );

        }

      } 

    } // End of requirement that there is at least one nearby Y-Plane hit above threshold. 

    // Find which planes have the greatest number of nearby points.
    if ( u_plane_nearby_hits_wire.size() > v_plane_nearby_hits_wire.size() || u_plane_nearby_hits_wire.size() > y_plane_nearby_hits_wire.size() ) u_plane_has_nonzero_nearby_points_and_in_the_top_two_planes = true;
 
    if ( v_plane_nearby_hits_wire.size() > u_plane_nearby_hits_wire.size() || v_plane_nearby_hits_wire.size() > y_plane_nearby_hits_wire.size() ) v_plane_has_nonzero_nearby_points_and_in_the_top_two_planes = true;

    if ( y_plane_nearby_hits_wire.size() > u_plane_nearby_hits_wire.size() || y_plane_nearby_hits_wire.size() > v_plane_nearby_hits_wire.size() ) y_plane_has_nonzero_nearby_points_and_in_the_top_two_planes = true;

    // Handle the weird of case of if they are all greater than zero and have the same size.
    if ( u_plane_nearby_hits_wire.size() > 0 && u_plane_nearby_hits_wire.size() == v_plane_nearby_hits_wire.size() && u_plane_nearby_hits_wire.size() == y_plane_nearby_hits_wire.size() ) u_plane_has_nonzero_nearby_points_and_in_the_top_two_planes = true;

    if ( y_plane_nearby_hits_wire.size() > 0 && y_plane_nearby_hits_wire.size() == u_plane_nearby_hits_wire.size() && y_plane_nearby_hits_wire.size() == v_plane_nearby_hits_wire.size() ) y_plane_has_nonzero_nearby_points_and_in_the_top_two_planes = true;

    // Find the case of >= two hits on both the u and the y planes.
    if ( u_plane_has_nonzero_nearby_points_and_in_the_top_two_planes == true && y_plane_has_nonzero_nearby_points_and_in_the_top_two_planes == true ) {

      double x_projection           = ( y_plane_nearby_hits_tick.at( furthest_away_y_plane_wire_idx ) - y_plane_vertex_tick );

      if ( x_projection < 0 ) proton_x_direction = -1;

      geo::WireID wid_vertex_u_plane(0, 0, 0, u_plane_wire_at_vertex);
      geo::WireID wid_vertex_y_plane(0, 0, 2, y_plane_wire_at_vertex);

      // Make 3D points of the vertex points and the furthest points and find the signs of the direction.
      geo::WireIDIntersection vertex_wire_intersection;

      bool do_vertex_wires_intersect = geom->WireIDsIntersect( wid_vertex_u_plane, wid_vertex_y_plane, vertex_wire_intersection );

      // Keep this print statement to avoid an error on unused variables. 
      std::cout << "The verdict on if the wires intersect = " << do_vertex_wires_intersect << "." << std::endl;
      
      double vertex_y_coordinate = vertex_wire_intersection.y;
      double vertex_z_coordinate = vertex_wire_intersection.z;

      geo::WireID wid_furthest_u_plane(0, 0, 0, u_plane_nearby_hits_wire.at( furthest_away_u_plane_wire_idx ) );
      geo::WireID wid_furthest_y_plane(0, 0, 2, y_plane_nearby_hits_wire.at( furthest_away_y_plane_wire_idx ) );

      geo::WireIDIntersection furthest_wire_intersection;

      bool do_furthest_wires_intersect = geom->WireIDsIntersect( wid_furthest_u_plane, wid_furthest_y_plane, furthest_wire_intersection );

      // Keep this print statement to avoid an error on unused variables. 
      std::cout << "Verdict on if the furthest 'uy' wires intersect = " << do_furthest_wires_intersect << "." << std::endl;

      double furthest_y_coordinate = furthest_wire_intersection.y;
      double furthest_z_coordinate = furthest_wire_intersection.z;      

      double y_projection = ( furthest_y_coordinate - vertex_y_coordinate );
      double z_projection = ( furthest_z_coordinate - vertex_z_coordinate );

      if ( y_projection < 0 ) proton_y_direction = -1;
      if ( z_projection < 0 ) proton_z_direction = -1;

    } // End of the case where the U and the Y have the greatest number of points used to match.

    // Find the case of >= two hits on both the v and the y planes.                                                                                                                                      
    else if ( v_plane_has_nonzero_nearby_points_and_in_the_top_two_planes == true && y_plane_has_nonzero_nearby_points_and_in_the_top_two_planes == true ) {

      double x_projection         = ( y_plane_nearby_hits_tick.at( furthest_away_y_plane_wire_idx ) - y_plane_vertex_tick );

      if ( x_projection < 0 ) proton_x_direction = -1;

      geo::WireID wid_vertex_v_plane(0, 0, 1, v_plane_wire_at_vertex);
      geo::WireID wid_vertex_y_plane(0, 0, 2, y_plane_wire_at_vertex);

      // Make 3D points of the vertex points and the furthest points and find the signs of the direction.                                                                                               
      geo::WireIDIntersection vertex_wire_intersection;

      bool do_vertex_wires_intersect = geom->WireIDsIntersect( wid_vertex_v_plane, wid_vertex_y_plane, vertex_wire_intersection );

      std::cout << "Verdict on if the vertex 'vy' wires intersect = " << do_vertex_wires_intersect << "." << std::endl;

      double vertex_y_coordinate = vertex_wire_intersection.y;
      double vertex_z_coordinate = vertex_wire_intersection.z;

      geo::WireID wid_furthest_v_plane(0, 0, 1, v_plane_nearby_hits_wire.at( furthest_away_v_plane_wire_idx ) );
      geo::WireID wid_furthest_y_plane(0, 0, 2, y_plane_nearby_hits_wire.at( furthest_away_y_plane_wire_idx ) );

      geo::WireIDIntersection furthest_wire_intersection;

      bool do_furthest_wires_intersect = geom->WireIDsIntersect( wid_furthest_v_plane, wid_furthest_y_plane, furthest_wire_intersection );

      // Keep this print statement to avoid an error on unused variables. 
      std::cout << "Verdict on if the furthest 'vy' wires intersect = " << do_furthest_wires_intersect << "." << std::endl;

      double furthest_y_coordinate = furthest_wire_intersection.y;
      double furthest_z_coordinate = furthest_wire_intersection.z;

      double y_projection = ( furthest_y_coordinate - vertex_y_coordinate );
      double z_projection = ( furthest_z_coordinate - vertex_z_coordinate );

      if ( y_projection< 0 ) proton_y_direction = -1;
      if ( z_projection< 0 ) proton_z_direction = -1;

    }   // End of the case where the V and the Y have the greatest number of points used to match.
      
    // Find the case of >= two hits on both the u and the v planes.                                                                                                                                     
    else if ( u_plane_has_nonzero_nearby_points_and_in_the_top_two_planes == true && v_plane_has_nonzero_nearby_points_and_in_the_top_two_planes == true ) {
      
      double x_projection         = ( u_plane_nearby_hits_tick.at( furthest_away_u_plane_wire_idx ) - u_plane_vertex_tick );

      if ( x_projection < 0 ) proton_x_direction = -1;

      geo::WireID wid_vertex_u_plane(0, 0, 0, u_plane_wire_at_vertex);
      geo::WireID wid_vertex_v_plane(0, 0, 1, v_plane_wire_at_vertex);

      // Make 3D points of the vertex points and the furthest points and find the signs of the direction.                                                                                               
      geo::WireIDIntersection vertex_wire_intersection;

      bool do_vertex_wires_intersect = geom->WireIDsIntersect( wid_vertex_u_plane, wid_vertex_v_plane, vertex_wire_intersection );

      // Keep this print statement to avoid an error on unused variables. 
      std::cout << "Verdict on if the vertex 'uv' wires intersect = " << do_vertex_wires_intersect << "." << std::endl;

      double vertex_y_coordinate = vertex_wire_intersection.y;
      double vertex_z_coordinate = vertex_wire_intersection.z;

      geo::WireID wid_furthest_u_plane(0, 0, 0, u_plane_nearby_hits_wire.at( furthest_away_u_plane_wire_idx ) );
      geo::WireID wid_furthest_v_plane(0, 0, 1, v_plane_nearby_hits_wire.at( furthest_away_v_plane_wire_idx ) );

      geo::WireIDIntersection furthest_wire_intersection;

      bool do_furthest_wires_intersect = geom->WireIDsIntersect( wid_furthest_u_plane, wid_furthest_v_plane, furthest_wire_intersection );

      // Keep this print statement to avoid an error on unused variables. 
      std::cout << "Verdict on if the furthest 'uv' wires intersect = " << do_furthest_wires_intersect << "." << std::endl;

      double furthest_y_coordinate = furthest_wire_intersection.y;
      double furthest_z_coordinate = furthest_wire_intersection.z;

      double y_projection = ( furthest_y_coordinate - vertex_y_coordinate );
      double z_projection = ( furthest_z_coordinate - vertex_z_coordinate );

      if ( y_projection< 0 ) proton_y_direction = -1;
      if ( z_projection< 0 ) proton_z_direction = -1;

    }   // End of the case where the U and the V have the greatest number of points used to match.

    else {

      // See if there is more than one point on the y-plane.
      if ( y_plane_nearby_hits_wire.size() > 0 )  {

	// Find the direction that x should point in according to the relative spacing of the ticks on the y-plane.                                                                                      
	double x_projection = ( y_plane_nearby_hits_tick.at( furthest_away_y_plane_wire_idx ) - y_plane_vertex_tick );

	if ( x_projection < 0 ) proton_x_direction = -1;

	if ( ( y_plane_nearby_hits_wire.at( furthest_away_y_plane_wire_idx ) - y_plane_wire_at_vertex ) < 0 ) proton_z_direction = -1;
	
      }

      else {

	// Flip the x-direction and the z-direction based on the orientation of the muon track.
	// Keep the x-direction and y-direction positive.
	if ( ( other_end_location_z - vertex_location_z ) > 0 ) proton_z_direction = -1;

      }

    } // End of the default cases.

    double muon_vertex_x     = 0.;
    double muon_vertex_y     = 0.;
    double muon_vertex_z     = 0.;
    double muon_ending_x     = 0.;
    double muon_ending_y     = 0.;
    double muon_ending_z     = 0.;

    // Set the muon vertex and ending point from the coordinates of the track.
    auto track_for_coords    = pandora_track_h->at( muon_candidate_idx );
    size_t starting_idx      = pandora_track_h->at( muon_candidate_idx ).FirstValidPoint();
    size_t ending_idx        = pandora_track_h->at( muon_candidate_idx ).LastValidPoint();
    
    // Set the muon coordinates here. 
    muon_vertex_x            = track_for_coords.LocationAtPoint( starting_idx ).X() - ( flash_time * 0.111436 );
    muon_vertex_y            = track_for_coords.LocationAtPoint( starting_idx ).Y();
    muon_vertex_z            = track_for_coords.LocationAtPoint( starting_idx ).Z();

    muon_ending_x            = track_for_coords.LocationAtPoint( ending_idx ).X() - ( flash_time * 0.111436 );
    muon_ending_y            = track_for_coords.LocationAtPoint( ending_idx ).Y();
    muon_ending_z            = track_for_coords.LocationAtPoint( ending_idx ).Z();      

    // Flip them if necessary.
    double dist_start_to_vtx = TMath::Sqrt( ( vertex_location_x - muon_vertex_x ) * ( vertex_location_x - muon_vertex_x ) + ( vertex_location_y - muon_vertex_y ) * ( vertex_location_y - muon_vertex_y ) + ( vertex_location_z - muon_vertex_z ) * ( vertex_location_z - muon_vertex_z ) );
    double dist_end_to_vtx   = TMath::Sqrt( ( vertex_location_x - muon_ending_x ) * ( vertex_location_x - muon_ending_x ) + ( vertex_location_y - muon_ending_y ) * ( vertex_location_y - muon_ending_y ) + ( vertex_location_z - muon_ending_z ) * ( vertex_location_z - muon_ending_z ) );

    if ( dist_end_to_vtx < dist_start_to_vtx ) {

      muon_vertex_x         = track_for_coords.LocationAtPoint( ending_idx ).X() - ( flash_time * 0.111436 );
      muon_vertex_y         = track_for_coords.LocationAtPoint( ending_idx ).Y();
      muon_vertex_z         = track_for_coords.LocationAtPoint( ending_idx ).Z();
      muon_ending_x         = track_for_coords.LocationAtPoint( starting_idx ).X() - ( flash_time * 0.111436 );
      muon_ending_y         = track_for_coords.LocationAtPoint( starting_idx ).Y();
      muon_ending_z         = track_for_coords.LocationAtPoint( starting_idx ).Z();

    }

    // Offset these track points for the space charge effect.
    geo::Vector_t muon_vertex_offsets = sce->GetCalPosOffsets(geo::Point_t{muon_vertex_x,muon_vertex_y,muon_vertex_z});
    
    double muon_vertex_x_with_offset  = muon_vertex_x - muon_vertex_offsets.X();
    double muon_vertex_y_with_offset  = muon_vertex_y + muon_vertex_offsets.Y();
    double muon_vertex_z_with_offset  = muon_vertex_z + muon_vertex_offsets.Z();
    
    geo::Vector_t muon_ending_offsets = sce->GetCalPosOffsets(geo::Point_t{muon_ending_x, muon_ending_y, muon_ending_z});
    
    double muon_ending_x_with_offset  = muon_ending_x - muon_ending_offsets.X();
    double muon_ending_y_with_offset  = muon_ending_y + muon_ending_offsets.Y();
    double muon_ending_z_with_offset  = muon_ending_z + muon_ending_offsets.Z();
    
    // Test all of the points for their position in the TPC.    
    // 3 cm containment requirement.
    if ( muon_vertex_x_with_offset < 3.0 || muon_vertex_x_with_offset > 253.35 )
      fails_3cm_containment_cut = true;
    
    if ( muon_vertex_y_with_offset < -113.5 || muon_vertex_y_with_offset > 113.5 )
      fails_3cm_containment_cut = true;
    
    if ( muon_vertex_z_with_offset < 3.0 || muon_vertex_z_with_offset > 1033.8 )
      fails_3cm_containment_cut = true;
    
    if ( muon_ending_x_with_offset < 3.0 || muon_ending_x_with_offset > 253.35 )
      fails_3cm_containment_cut = true;
    
    if ( muon_ending_y_with_offset < -113.5 || muon_ending_y_with_offset > 113.5 )
      fails_3cm_containment_cut = true;
    
    if ( muon_ending_z_with_offset < 3.0 || muon_ending_z_with_offset > 1033.8 )
      fails_3cm_containment_cut = true;
    
    // Reset the variable for the total track length to '0'.
    total_length_of_tracks_originating_from_vertex = 0.;
    
    // Find the total length of the tracks originating at the vertex.
    total_length_of_tracks_originating_from_vertex += ubxsec_event->slc_longesttrack_length[nu_slc_idx];
    
    // Loop through the track in the TPC Object and find if they originate at the vertex.  Skip the longest track.
    for ( size_t tpcobject_track_iter = 0; tpcobject_track_iter < tracks.size(); tpcobject_track_iter++ ) {
      
      // Skip the longest track (the vertex is found with respect to that one).
      if ( fabs( tracks.at( tpcobject_track_iter )->Length() - ubxsec_event->slc_longesttrack_length[nu_slc_idx] ) < 0.001 )
	continue;
      
      // Find the starting and ending x, y, and z coordinates of the muon track and of the track that you are currently looking at.
      auto track                                     = *(tracks.at( tpcobject_track_iter ));
      
      auto first_valid_point_idx                     = track.FirstValidPoint();
      auto last_valid_point_idx                      = track.LastValidPoint();
      
      auto first_valid_point                         = track.LocationAtPoint( first_valid_point_idx );
      auto last_valid_point                          = track.LocationAtPoint( last_valid_point_idx );
      
      double other_track_x0                          = first_valid_point.X() - ( 0.111436 * flash_time );
      double other_track_y0                          = first_valid_point.Y();
      double other_track_z0                          = first_valid_point.Z();
      double other_track_x1                          = last_valid_point.X() - ( 0.111436 * flash_time );
      double other_track_y1                          = last_valid_point.Y();
      double other_track_z1                          = last_valid_point.Z();
      
      // Offset these track points for the space charge effect. 
      geo::Vector_t other_track_first_point_offsets  = sce->GetCalPosOffsets(geo::Point_t{other_track_x0,other_track_y0,other_track_z0});
      
      double other_track_x0_with_offset              = other_track_x0 - other_track_first_point_offsets.X();
      double other_track_y0_with_offset              = other_track_y0 + other_track_first_point_offsets.Y();
      double other_track_z0_with_offset              = other_track_z0 + other_track_first_point_offsets.Z();
      
      geo::Vector_t other_track_second_point_offsets = sce->GetCalPosOffsets(geo::Point_t{other_track_x1,other_track_y1,other_track_z1});
      
      double other_track_x1_with_offset              = other_track_x1 - other_track_second_point_offsets.X();
      double other_track_y1_with_offset              = other_track_y1 + other_track_second_point_offsets.Y();
      double other_track_z1_with_offset              = other_track_z1 + other_track_second_point_offsets.Z();
      
      // Test all of the points for their position in the TPC.                                                                                                                                            
      // 3 cm containment requirement.                                                                                                                                                                 
      if ( other_track_x0_with_offset < 3.0 || other_track_x0_with_offset > 253.35 )
	fails_3cm_containment_cut = true;
      
      if ( other_track_y0_with_offset < -113.5 || other_track_y0_with_offset > 113.5 )
	fails_3cm_containment_cut = true;
      
      if ( other_track_z0_with_offset < 3.0 || other_track_z0_with_offset > 1033.8 )
	fails_3cm_containment_cut = true;
      
      if ( other_track_x1_with_offset < 3.0 || other_track_x1_with_offset > 253.35 )
	fails_3cm_containment_cut = true;
      
      if ( other_track_y1_with_offset < -113.5 || other_track_y1_with_offset > 113.5 )
	fails_3cm_containment_cut = true;
      
      if ( other_track_z1_with_offset < 3.0 || other_track_z1_with_offset > 1033.8 )
	fails_3cm_containment_cut = true;
                  
      // Find which point is closer to both the vertex and the end of the muon candidate.
      double closest_distance_vertex = 0.;
      double closest_distance_end    = 0.;
      
      double distance_vertex_point0  = TMath::Sqrt( ( muon_vertex_x - other_track_x0 ) * ( muon_vertex_x - other_track_x0 ) + ( muon_vertex_y - other_track_y0 ) * ( muon_vertex_y - other_track_y0 ) + ( muon_vertex_z - other_track_z0 ) * ( muon_vertex_z - other_track_z0 ) );

      double distance_vertex_point1  = TMath::Sqrt( ( muon_vertex_x - other_track_x1 ) * ( muon_vertex_x - other_track_x1 ) + ( muon_vertex_y - other_track_y1 ) * ( muon_vertex_y - other_track_y1 ) + ( muon_vertex_z - other_track_z1 ) * ( muon_vertex_z - other_track_z1 ) );

      double distance_end_point0     = TMath::Sqrt( ( muon_ending_x - other_track_x0 ) * ( muon_ending_x - other_track_x0 ) + ( muon_ending_y - other_track_y0 ) * ( muon_ending_y - other_track_y0 ) + ( muon_ending_z - other_track_z0 ) * ( muon_ending_z - other_track_z0 ) );
      
      double distance_end_point1     = TMath::Sqrt( ( muon_ending_x - other_track_x1 ) * ( muon_ending_x - other_track_x1 ) + ( muon_ending_y - other_track_y1 ) * ( muon_ending_y - other_track_y1 ) + ( muon_ending_z - other_track_z1 ) * ( muon_ending_z - other_track_z1 ) );
            
      if ( distance_vertex_point0 < distance_vertex_point1 ) {
	closest_distance_vertex = distance_vertex_point0;
      }
      
      else {
	closest_distance_vertex = distance_vertex_point1;
      }
      
      if ( distance_end_point0 < distance_end_point1 ) {
	closest_distance_end = distance_end_point0;
      }
      
      else {
	closest_distance_end = distance_end_point1;
      }
      
      // Now, determine if the track originates at the vertex or at the other end.
      if ( closest_distance_vertex < closest_distance_end ) {
	
	total_length_of_tracks_originating_from_vertex += track.Length();
	
	// See if this track can be considered part of the phase space that we are looking at.
	if ( tracks.size() == 2 ) {
	  
	  // Go through all of the possibilities of impossible combinations of track lengths.
	  if ( track.Length() > 1.0 && ubxsec_event->slc_longesttrack_length[nu_slc_idx] > 36.0 ) 
	    fails_new_two_track_selection = true;
	  
	  if ( track.Length() > 2.0 && ubxsec_event->slc_longesttrack_length[nu_slc_idx] > 28.0 )
	    fails_new_two_track_selection = true;
	  
	  if ( track.Length() > 3.0 && ubxsec_event->slc_longesttrack_length[nu_slc_idx] > 24.0 )
	    fails_new_two_track_selection = true;
	  
	  if ( track.Length() > 4.0 && ubxsec_event->slc_longesttrack_length[nu_slc_idx] > 20.0 )
	    fails_new_two_track_selection = true;

	  if ( track.Length() > 5.0 && ubxsec_event->slc_longesttrack_length[nu_slc_idx] > 11.0 )
	    fails_new_two_track_selection = true;
	  
	  if ( track.Length() > 6.0 && ubxsec_event->slc_longesttrack_length[nu_slc_idx] > 8.0 )
	    fails_new_two_track_selection = true;
	  
	  if ( track.Length() > 7.0 && ubxsec_event->slc_longesttrack_length[nu_slc_idx] > 7.0 )
	    fails_new_two_track_selection = true;

	  if ( track.Length() > 11.0)
	    fails_new_two_track_selection = true;
	  
	} // End of requirement that there be two tracks in the TPC Object.
	
      } // End of the requirement that the track be closer to the vertex than the other end.
      
    } // End of the loop over the tracks in the TPC Object.

    // Set variables for the candidate muon within the TPC object.
    muon_candidate_length          = ubxsec_event->slc_longesttrack_length[nu_slc_idx];
    reconstructed_muon_KE          = ( TMath::Sqrt( ( p_calculator_from_length.GetTrackMomentum( muon_candidate_length, 13 ) * 1000. ) * ( p_calculator_from_length.GetTrackMomentum( muon_candidate_length, 13 ) * 1000. ) + 105.7 * 105.7 ) - 105.7 );
  
    number_of_tracks_in_TPCObject  = tracks.size();
    
    sum_of_TPCObject_track_lengths = 0.;
  
    // See how many of the tracks in the TPC object have a length greater than 5 cm.
    for ( size_t tpcobj_track_iter = 0; tpcobj_track_iter < tracks.size(); tpcobj_track_iter++ ) {
    
      sum_of_TPCObject_track_lengths += tracks.at( tpcobj_track_iter )->Length();
    
      // Calculate the energy of the track using the muon and pion assumption.
      if ( fabs( tracks.at( tpcobj_track_iter )->Length() - muon_candidate_length ) > 0.001 ) 
	TPCObject_kinetic_energy_using_muon_and_pion_assumption += ( TMath::Sqrt( ( p_calculator_from_length.GetTrackMomentum( tracks.at( tpcobj_track_iter )->Length(), 111 ) * 1000. ) * ( p_calculator_from_length.GetTrackMomentum( tracks.at( tpcobj_track_iter )->Length(), 111 ) * 1000. ) + 135.0 * 135.0 ) - 135.0 );											   
          
    }

    // Include a variable for the vertex position before it is t0-corrected.
    ubxsec_muon_phi                = ubxsec_event->slc_longesttrack_phi[nu_slc_idx];
    ubxsec_muon_cos_theta           = ubxsec_event->slc_longesttrack_theta[nu_slc_idx];

    // Calculate the distance between the truth vertex and the reconstructed vertex.
    truth_reco_vtx_distance          = TMath::Sqrt( ( nu_vtx_x_truth - vertex_location_x ) * ( nu_vtx_x_truth - vertex_location_x )  + ( nu_vtx_y_truth - vertex_location_y ) * ( nu_vtx_y_truth - vertex_location_y ) + ( nu_vtx_z_truth - vertex_location_z ) * ( nu_vtx_z_truth - vertex_location_z ) );

    // Find out which tracks are associated to the TPCObject and make a vector of their hits.
    art::Handle<std::vector<ubana::TPCObject> > tpcobject_h;                                                                                                                                         
    e.getByLabel("TPCObjectMaker", tpcobject_h);                                                                                                                                                       
    
    // make sure TPC Object info looks good                                                                                                                                                             
    if(!tpcobject_h.isValid()) {                                                                                                                                                                        
      std::cerr<<"\033[93m[ERROR]\033[00m ... could not locate TPC Object!"<<std::endl;                                                                                                                    
      throw std::exception();                                                                                                                                                                           
    }  

    art::FindMany<recob::Track> tpcobject_track_assn_v(tpcobject_h, e, "TPCObjectMaker");
    const std::vector< const recob::Track* >& tracks_associated_to_TPC_object = tpcobject_track_assn_v.at( nu_slc_idx );
    
    std::vector< int > TrackID_set;
    TrackID_set.clear();
    
    // Loop through the tracks and print out their ID.
    for ( size_t charge_track_iter = 0; charge_track_iter < tracks_associated_to_TPC_object.size(); charge_track_iter++ ) {

      auto const& track = tracks_associated_to_TPC_object.at( charge_track_iter );

      TrackID_set.push_back( track->ID() );

    }

    art::FindMany< recob::Hit > track_hit_assn_v(pandora_track_h, e, "pandora");

    // Find the number of boundary crossings at this point in the code.

    // Loop through the tracks in the event and print out the information about the tracks that are associated to the TPC Object.
    for ( size_t pandora_iter = 0; pandora_iter < pandora_track_h->size(); pandora_iter++ ) {

      // See how many pandora tracks are near the Y and Z boundaries of the detector.
      // You can only do this within 20 cm because you don't have X-information; you have to x-correct the track endpoints.
      auto pandora_track = pandora_track_h->at( pandora_iter );
      
      auto first_valid_point_idx = pandora_track.FirstValidPoint();
      auto last_valid_point_idx  = pandora_track.LastValidPoint();

      auto first_valid_point     = pandora_track.LocationAtPoint( first_valid_point_idx );
      auto last_valid_point      = pandora_track.LocationAtPoint( last_valid_point_idx );

      double pandora_track_y0    = first_valid_point.Y();
      double pandora_track_z0    = first_valid_point.Z();
      double pandora_track_y1    = last_valid_point.Y();
      double pandora_track_z1    = last_valid_point.Z();

      if ( pandora_track_y0 > 96.5 || pandora_track_y1 > 96.5 ) 
	num_of_top_pandora_crossings++;

      if ( pandora_track_y0 < -96.5 || pandora_track_y1 < -96.5 )
	num_of_bottom_pandora_crossings++;

      if ( pandora_track_z0 < 20.0 || pandora_track_z1 < 20.0 )
	num_of_front_pandora_crossings++;

      if ( pandora_track_z0 > 1016.8 || pandora_track_z1 > 1016.8 )
	num_of_back_pandora_crossings++;

      int track_id = pandora_track_h->at( pandora_iter ).ID();

      for ( size_t TrackID_iter = 0; TrackID_iter < TrackID_set.size(); TrackID_iter++ ) {

	if ( track_id == TrackID_set.at( TrackID_iter ) ) {

	  const std::vector<const recob::Hit*>& hits = track_hit_assn_v.at( pandora_iter ); 
	  
	  // Loop through the vector and fill the relevant vectors.                                                                                                                                  
	  for ( size_t hit_iter = 0; hit_iter < hits.size(); hit_iter++ ) {                                                                                                                           
	    
	    auto const& hit = hits.at( hit_iter );                                                                                                                                                     
                                                                                                                                                                                                
	    total_sum_of_slice_associated_hits_ADCs      += hit->Integral();                                                                                                                           
	    total_number_of_slice_associated_hits++;
                                                                                                                                                                                            
	    if ( hit->View() == 0 ) {                                                                                                                                                                   
                                                                                                                                                                                                   
	      u_plane_sum_of_slice_associated_hits_ADCs  += hit->Integral();                                                                                                                         
	      u_plane_number_of_slice_associated_hits++;
                                                                                                                                                 
	    }                                            
                                                                                                                                             
                                                                                                                                                                                                      
	    if ( hit->View() == 1 ) {                                                                                                                                                                 
                                                                                                                                                                                                       
	      v_plane_sum_of_slice_associated_hits_ADCs  += hit->Integral();                                                                                                                            
	      v_plane_number_of_slice_associated_hits++;
                                                                                                                                                            
	    }                                                                                                                                                                                          
                                                                                                                                                                                                          
	    if ( hit->View() == 2 ) {                                                                                                                                                                   
                                                                                                                                                                                                       
	      y_plane_sum_of_slice_associated_hits_ADCs  += hit->Integral();               
	      y_plane_number_of_slice_associated_hits++;
                                                                                                     
	    }                                                                                                                                                                      

	  }
	  
	  // Break the loop after you have found the track ID.
	  break;

	}

      } // End of the loop over the Track IDs.

    }

    // Put in the CRT distance cut right here - this rejects an event based on the minimum distance between a neutrino vertex & a CRT-tagged track.
    // if so, let's fetch which reco tracks have a CRT Hit association                                                                                                                                   
    std::vector<int> matchedRecoTrkKey;
    matchedRecoTrkKey.clear();
    // Get CRT Hits containers                                                                                                                                                                       
    art::Handle< std::vector<crt::CRTHit> > crtHitHandle;  // Container of crt hits                                                                                                                   
    std::vector<art::Ptr<crt::CRTHit> >     crtHits;       // Vector of crt hits                                                                                                                        
    e.getByLabel("crthitcorr", crtHitHandle);
    art::fill_ptr_vector(crtHits, crtHitHandle);
    // Find associations                                                                                                                                                                               
    art::FindManyP<recob::Track> fCRT2TPC(crtHitHandle, e , "crttrackmatch");

    int num_tracks_associated_to_CRT = 0;

    if (fCRT2TPC.isValid())
      {
	for (unsigned int indexAssn = 0; indexAssn < fCRT2TPC.size(); ++indexAssn )
	  {
	    auto assTracksPtr = fCRT2TPC.at(indexAssn);
	    if (assTracksPtr.size()==0) continue;
	    
	    num_tracks_associated_to_CRT++;

	    recob::Track const& aTrack = *assTracksPtr.front();
	    matchedRecoTrkKey.push_back(aTrack.ID());
	  }// End Loop on Association                                                                                                                                                                       
      }// if there's an associated track

    //Scope of this bit of code is to create an association between bad vtx and tracks 
    //Loop on verteces                                                                                                                                                                                  
    int _nTracks       = pandora_track_h->size();
    int _nTrackCRTAssn = matchedRecoTrkKey.size();
    int trackUsed      =  0;

    // Reset the variable for the closest approach of a cosmic tagged by the CRT.
    closest_distance_to_CRT_tagged_track = 10000.;

    // Print out the number of tracks + the number of tracks associated to the CRT.
    std::cout << "N tracks, N tracks assn2CRT: " << _nTracks << " " << _nTrackCRTAssn << " \n";

    for ( size_t track_iter = 0; track_iter < pandora_track_h->size(); track_iter++ ) {
      auto thisTrack = pandora_track_h->at( track_iter );
      int trackID    = thisTrack.ID();

      // I need to check if this track ID is in the keys                                                                                                                                                
      // if it is, I'm interested in this track. Otherwise, I'm not and I shall continue.                                                                                                               
      if(! (std::find(matchedRecoTrkKey.begin(), matchedRecoTrkKey.end(), trackID) != matchedRecoTrkKey.end())) continue;
      trackUsed++;

      // Grab this track TrajectorPoints & calculate the trajectory                                                                                                                                    
      for(size_t p = 1; p < thisTrack.NumberTrajectoryPoints(); ++p){
	const auto& pos_cur       = thisTrack.LocationAtPoint(p);
	double dx                 =  vertex_location_x_with_SCE - ( pos_cur.x() - 0.111436 * flash_time );
	double dy                 =  vertex_location_y_with_SCE - pos_cur.y();
	double dz                 =  vertex_location_z_with_SCE - pos_cur.z();
	double thisTrajPtDistance = TMath::Sqrt(dx*dx+dy*dy+dz*dz);
	
	if ( thisTrajPtDistance < closest_distance_to_CRT_tagged_track ) {
	  closest_distance_to_CRT_tagged_track  = thisTrajPtDistance;
	}

	if ( closest_distance_to_CRT_tagged_track < 15.0 ) {
	  fails_CRT_distance_cut = true;
	  break;
	}
	
      } // End of the loop over the track trajectory points.

      if ( fails_CRT_distance_cut == true ) 
	break;
	      
    } // End of the loop over the tracks.

    // Go through the cuts that required more work than what was shown originally.
    if ( fails_3cm_containment_cut == true )  { 
      failure_map["j_containment_cut"]                      = false;
      reason                                                = "fail_containment_cut";
      total_num_events_failing_3cm_containment_cut++;
      total_num_events_failing++;
    }

    else if ( fails_fiducial_volume_cut == true) {
      failure_map["k_fiducial_volume_cut"]               = false;
      reason                                                 = "fail_fiducial_volume";
      total_num_events_failed_fiducial_volume++;
      total_num_events_failing++;
    } 
    
    else if ( number_of_tracks_in_TPCObject > 3 )  {
     failure_map["n_number_of_tracks_in_TPCObject_requirement"]       = false;
     reason = "fail_less_than_4_tracks_in_tpcobject_cut";
     total_num_events_failing_less_than_4_tracks_in_tpcobject_cut++;
     total_num_events_failing++;
    }
   
    else if ( num_of_tracks_originating_from_vertex > 2 ) {
     failure_map["o_number_of_tracks_originating_from_vertex_requirement"]       = false;
     reason = "fail_more_than_2_tracks_originating_from_vertex_cut";
     total_num_events_failing_more_than_2_tracks_originating_from_vertex_cut++;
     total_num_events_failing++;
    }

    else if ( fails_CRT_distance_cut == true ) {
     failure_map["p_CRT_distance_cut"]       = false;
     reason = "fails_CRT_distance_cut";
     total_num_events_failing_CRT_distance_cut++;
     total_num_events_failing++;
    }

    // Look at the cuts that will be used in our sidebands.
    else  {

      if ( total_length_of_tracks_originating_from_vertex > 40.0 ) {
	failure_map["l_length_of_tracks_from_vertex_cut"] = false;
	reason                                                = "fail_sum_of_track_lengths_from_vertex";
	total_num_events_failing_length_of_tracks_from_vertex_requirement++;
	total_num_events_failing++;
      }

      else if ( sum_of_TPCObject_track_lengths > 65.0 ) {
        failure_map["q_length_of_tracks_in_TPCObject_cut"] = false;
        reason                                             = "fail_sum_of_TPCObject_track_lengths";
        total_num_events_failing_total_length_of_tracks_in_TPC_Object++;
        total_num_events_failing++;

      }

      else if ( fails_new_two_track_selection == true ) {
	event_fails_new_two_track_selection = 1;
	failure_map["m_new_two_track_length_requirement"]       = false;
	reason                                                = "fail_new_two_track_requirement";
	total_num_events_failing_new_two_track_selection++;
	total_num_events_failing++;
      }
      
      else 
	total_num_events_passing++;

      if ( kdar_from_dump_event == true )
	num_kdar_events_passing_in_this_sample++;
     
      // Fill both trees.                                                                                                                                                                   
      _passing_events_tree->Fill();
      proton_hit_identification_tree->Fill();
      
     
      selection_result.SetSelectionStatus(true);
     
      ubxsec_event->is_selected = true;
      
      // Grab the selected TPCObject
      std::vector<art::Ptr<ubana::TPCObject>> tpcobj_v;
      art::fill_ptr_vector(tpcobj_v, tpcobj_h);
      art::Ptr<ubana::TPCObject> new_tpcobj = tpcobj_v.at(slice_index_1);
      
      if (_debug) std::cout << "[UBXSec] >>>>>>>>>>>>>>>>>>>>>> Selected TPCObject with index " << slice_index_1 << std::endl;
     
      // Prepare the tpcobj output
      std::vector<art::Ptr<ubana::TPCObject>>  out_tpcobj_v;
      out_tpcobj_v.resize(1);
      out_tpcobj_v.at(0) = new_tpcobj;
     
      selectionResultVector->emplace_back(std::move(selection_result));
      util::CreateAssn(*this, e, *selectionResultVector, out_tpcobj_v, *assnOutSelectionResultTPCObject);
      
      // For the TPCNeutrinoID Filter
      util::CreateAssn(*this, e, muon_candidate_track_per_slice_v.at(slice_index_1), neutrino_candidate_vertex_per_slice_v.at(slice_index_1), *vertexTrackAssociations);
      util::CreateAssn(*this, e, muon_candidate_pfparticle_per_slice_v.at(slice_index_1), neutrino_candidate_vertex_per_slice_v.at(slice_index_1), *vertexPFParticleAssociations);
      
    } // End of the 'else' statement in which the event is selected.
    
  } // End of the case in which the event is selected.

  else {
   
    std::cout << "This event failed because of " << reason << "."<< std::endl;
   
  }

  if ( is_duplicate_event )
    total_num_duplicate_events++;

  // Print out the information to see where you're at.
  // Total # of Events                                                                                                                                                                                 
  std::cout << "The total number of events = "                                                                         << total_num_events                                                        << "." << std::endl;
  std::cout << "The total number of events passing = "                                                                 << total_num_events_passing                                                << "." << std::endl;
  std::cout << "The total number of events failing = "                                                                 << total_num_events_failing                                                << "." << std::endl;
  std::cout << "The total number of duplicate events = "                                                               << total_num_duplicate_events                                              << "." << std::endl;
  
  // Cut Breakdown                                                              
  std::cout << "The total number of events failing the CRT veto = "                                                    << total_num_events_failed_crt_veto                                        << "." << std::endl;
  std::cout << "The total number of events with no 'simpleFlashBeam' flashes = "                                       << total_num_events_failed_beam_disc_flashes                               << "." << std::endl;
  std::cout << "The total number of events without a 'simpleFlashBeam' flash of > 50 PEs in the beamspill window = "   << total_num_events_failed_beam_spill_flash                                << "." << std::endl;
  std::cout << "The total number of events with no TPC Objects = "                                                     << total_num_events_failed_has_slices                                      << "." << std::endl;
  std::cout << "The total number of events with no TPC Object that is tagged as the neutrino = "                       << total_num_events_failed_has_slice_tagged_as_neutrino                    << "." << std::endl;
  std::cout << "The total number of events with a TPC Object failing the > 40 cm track length cut = "                  << total_num_events_failed_track_length                                    << "." << std::endl;
  std::cout << "The total number of events failing the fiducial volume cut = "                                         << total_num_events_failed_fiducial_volume                                 << "." << std::endl;
  std::cout << "The total number of events failing the number of tracks cut = "                                        << total_num_events_failed_ntrack                                          << "." << std::endl;
  std::cout << "The total number of events failing the residuals std up cut = "                                        << total_num_events_failed_residuals_std_up                                << "." << std::endl;
  std::cout << "The total number of events failing the perc hits used in cluster = "                                   << total_num_events_failed_perc_used_hits_in_cluster                       << "." << std::endl;
  std::cout << "The total number of events failing the new two track selection = "                                     << total_num_events_failing_new_two_track_selection                        << "." << std::endl;
  std::cout << "The total number of events failing the 3 cm containment cut = "                                        << total_num_events_failing_3cm_containment_cut                            << "." << std::endl;
  std::cout << "The total number of events failing the sum of the track lengths coming from the vertex requirement = " << total_num_events_failing_length_of_tracks_from_vertex_requirement       << "." << std::endl;
  std::cout << "The total number of events failing the sum of the track lengths in the TPC Object = "                  << total_num_events_failing_total_length_of_tracks_in_TPC_Object           << "." << std::endl;
  std::cout << "The total number of events failing the less than 4 tracks in the TPC Object cut = "                    << total_num_events_failing_less_than_4_tracks_in_tpcobject_cut            << "." << std::endl;
  std::cout << "The total number of events failing the more than 2 tracks originating from the vertex cut = "          << total_num_events_failing_more_than_2_tracks_originating_from_vertex_cut << "." << std::endl;
  std::cout << "The total number of events failing the CRT distance cut = "                                            << total_num_events_failing_CRT_distance_cut                               << "." << std::endl;

  std::cout << "\n"                                                                                                    << std::endl;
  std::cout << "The number of KDAR events in this sample = "                                                           << num_kdar_events_in_this_sample                                          << "." << std::endl;
  std::cout << "The number of KDAR events passing in this sample = "                                                   << num_kdar_events_passing_in_this_sample                                  << "." << std::endl;
  
  if(_debug) std::cout << "[UBXSec] Filling tree now." << std::endl;
  _tree1->Fill();

  // Fill the duplicate event info.
  duplicate_runs.push_back( _run );
  duplicate_subruns.push_back( _subrun );
  duplicate_events.push_back( _event );
  
  e.put(std::move(selectionResultVector));
  e.put(std::move(assnOutSelectionResultTPCObject));
  
  e.put(std::move(vertexTrackAssociations));
  e.put(std::move(vertexPFParticleAssociations));

  num_events_looping_over_tree->Fill();

  if(_debug) std::cout << "********** UBXSec ends" << std::endl;

  return;

} // End of the UBXSec function.



void UBXSec::endSubRun(art::SubRun& sr) {

  if (_debug) std::cout << "[UBXSec::endSubRun] Starts" << std::endl;

  _sr_run       = sr.run();
  _sr_subrun    = sr.subRun();
  _sr_begintime = sr.beginTime().value();
  _sr_endtime   = sr.endTime().value();

  art::Handle<sumdata::POTSummary> potsum_h;

  // MC
  if (_is_mc) {
    if (_debug) std::cout << "[UBXSec::endSubRun] Getting POT for MC" << std::endl;
    if(sr.getByLabel(_potsum_producer, potsum_h)) {
      if (_debug) std::cout << "[UBXSec::endSubRun] POT are valid" << std::endl;
      _sr_pot = potsum_h->totpot;
    }
    else
      _sr_pot = 0.;
  }

  // Data - Use Zarko's script instead
  if (_is_data) {
 //   if (_debug) std::cout << "[UBXSec::endSubRun] Getting POT for DATA, producer " << _potsum_producer << ", instance " << _potsum_instance << std::endl;
    if (sr.getByLabel(_potsum_producer,  potsum_h)){
      if (_debug) std::cout << "[UBXSec::endSubRun] POT are valid" << std::endl;
      _sr_pot = potsum_h->totpot;
    }
    else
      _sr_pot = 0;
  }

  _sr_tree->Fill();

  if (_debug) std::cout << "[UBXSec::endSubRun] Ends" << std::endl;
}




void UBXSec::PrintMC(std::vector<art::Ptr<simb::MCTruth>> mclist) {

  std::cout << "[UBXSec] ================= MC Information ================= [UBXSec]" << std::endl;

  int iList = 0;
  std::cout << " NEUTRINO:" << std::endl;
  if (mclist[iList]->NeutrinoSet()) {
    std::cout << "\tPDG      " << mclist[iList]->GetNeutrino().Nu().PdgCode() << std::endl;
    std::cout << "\tCC/NC?   " << (mclist[iList]->GetNeutrino().CCNC() == 0 ? "CC" : "NC") << std::endl;
    std::cout << "\tMode     " << mclist[iList]->GetNeutrino().Mode() << std::endl;
    std::cout << "\tQSqr     " << mclist[iList]->GetNeutrino().QSqr() << std::endl;
    std::cout << "\tW        " << mclist[iList]->GetNeutrino().W() << std::endl;
    std::cout << "\tX        " << mclist[iList]->GetNeutrino().X() << std::endl;
    std::cout << "\tY        " << mclist[iList]->GetNeutrino().Y() << std::endl;
    std::cout << "\tHitNuc   " << mclist[iList]->GetNeutrino().HitNuc() << std::endl;
    std::cout << "\tE        " << mclist[iList]->GetNeutrino().Nu().E() << std::endl;
    std::cout << "\tVx       " << mclist[iList]->GetNeutrino().Nu().Vx() << std::endl;
    std::cout << "\tVy       " << mclist[iList]->GetNeutrino().Nu().Vy() << std::endl;
    std::cout << "\tVz       " << mclist[iList]->GetNeutrino().Nu().Vz() << std::endl;

  } else
    std::cout << "\t---No Neutrino information---" << std::endl;

  std::cout << std::endl;
  std::cout << " PRIMARIES (only with status code==1):" << std::endl;
  for (int p = 0; p < mclist[0]->NParticles(); p++) {
    const simb::MCParticle mc_par = mclist[0]->GetParticle(p);
    if (mc_par.StatusCode() != 1) continue;
    std::cout << "\tPDG           " << mc_par.PdgCode() << std::endl;
    std::cout << "\tStart process " << mc_par.Process() << std::endl;
    std::cout << "\tEnd process   " << mc_par.EndProcess() << std::endl;
    std::cout << "\tEnergy        " << mc_par.E() << std::endl;
    std::cout << "\tMomentum      " << mc_par.P() << std::endl;
    std::cout << "\tVertex        " << mc_par.Vx() << ", " << mc_par.Vy() << ", " << mc_par.Vz() << std::endl;
    std::cout << "\tStatus Code   " << mc_par.StatusCode() << std::endl << std::endl;
  }

  std::cout << "[UBXSec] ================= MC Information ================= [UBXSec]" << std::endl;
}


//_______________________________________________________________________________________
float UBXSec::GetTrackShowerScore(art::Ptr<recob::PFParticle> pfParticle, const lar_pandora::PFParticlesToMetadata pfParticleToMetadata)
{
  // Find the PFParticle in the metadata map
   auto itParticle = pfParticleToMetadata.find(pfParticle);
  // Check the PFParticle was found
  if (itParticle == pfParticleToMetadata.end())
    throw cet::exception("WorkshopTrackShowerHelperComplete") << "PFParticle has no metadata" << std::endl;
  // There should only be one metadata for each PFParticle
  if (itParticle->second.size() != 1)
    throw cet::exception("WorkshopTrackShowerHelperComplete") << "PFParticle has mutiple metadata" << std::endl;
  // The metadata vector has size one as required, so just take the first element
  const auto metadata = itParticle->second.front();
  // Now get the properties map - this is a map from a string (the property name) to a float (the property value)
  const auto propertiesMap = metadata->GetPropertiesMap();
  // Look for the track shower ID score
  const auto itScore = propertiesMap.find("TrackScore");
  // Check the track score was available
  if (itScore == propertiesMap.end())
    throw cet::exception("WorkshopTrackShowerHelperComplete") << "PFParticle has no track score" << std::endl;
  return itScore->second;
}
void UBXSec::GetFlashLocation(std::vector<double> pePerOpDet,
                              double& Ycenter,
                              double& Zcenter,
                              double& Ywidth,
                              double& Zwidth)
{

  // Reset variables
  Ycenter = Zcenter = 0.;
  Ywidth  = Zwidth  = -999.;
  double totalPE = 0.;
  double sumy = 0., sumz = 0., sumy2 = 0., sumz2 = 0.;

  for (unsigned int opdet = 0; opdet < pePerOpDet.size(); opdet++) {

    // Get physical detector location for this opChannel
    double PMTxyz[3];
    ::art::ServiceHandle<geo::Geometry> geo;
    geo->OpDetGeoFromOpDet(opdet).GetCenter(PMTxyz);

    // Add up the position, weighting with PEs
    sumy    += pePerOpDet[opdet]*PMTxyz[1];
    sumy2   += pePerOpDet[opdet]*PMTxyz[1]*PMTxyz[1];
    sumz    += pePerOpDet[opdet]*PMTxyz[2];
    sumz2   += pePerOpDet[opdet]*PMTxyz[2]*PMTxyz[2];

    totalPE += pePerOpDet[opdet];
  }

  Ycenter = sumy/totalPE;
  Zcenter = sumz/totalPE;

  // This is just sqrt(<x^2> - <x>^2)
  if ( (sumy2*totalPE - sumy*sumy) > 0. )
    Ywidth = std::sqrt(sumy2*totalPE - sumy*sumy)/totalPE;

  if ( (sumz2*totalPE - sumz*sumz) > 0. )
    Zwidth = std::sqrt(sumz2*totalPE - sumz*sumz)/totalPE;
}
//____________________________________________________________________________________________
void UBXSec::GetTaggedPFP(art::Event const & e, std::string cosmictag_producer, double score_cut, lar_pandora::PFParticleVector & pfpTaggedOut,std::vector<int> & tagid_v){

  pfpTaggedOut.clear();
  tagid_v.clear();

  if (_debug_cr) std::cout << "Getting cosmic tags from " << cosmictag_producer << std::endl;

  // Get the CosmicTag from the ART event
  art::Handle<std::vector<anab::CosmicTag>> cosmicTagHandle;
  e.getByLabel(cosmictag_producer, cosmicTagHandle);

  if (!cosmicTagHandle.isValid() || cosmicTagHandle->empty()){
    std::cerr << "Cosmic tag " << cosmictag_producer << " is not valid or empty." << std::endl;
    std::cout <<"valid? "<<cosmicTagHandle.isValid()<<"  empty?? "<<cosmicTagHandle->empty()<<std::endl;
    return;
  }

  // Look up the associations to PFPs
  art::FindManyP<recob::PFParticle> cosmicPFPAssns(cosmicTagHandle, e, cosmictag_producer);

  if (_debug_cr) std::cout << " cosmicPFPAssns.size(): " << cosmicPFPAssns.size() << std::endl;

  // Loop over the cosmic tags
  for (unsigned int ct = 0; ct < cosmicPFPAssns.size(); ct++) {

    // Get the cosmic tag
    art::Ptr<anab::CosmicTag> cosmicTag(cosmicTagHandle, ct);
    //  if(_debug_cr) std::cout << "This cosmic tag (" << ct << ") has type: " << cosmicTag->CosmicType() << " and score: " << cosmicTag->CosmicScore() << std::endl;

    // Get the PFP associated with this CT
    std::vector<art::Ptr<recob::PFParticle>> cosmicTagToPFP_v = cosmicPFPAssns.at(cosmicTag.key());
    //if(_debug) std::cout << "Number of PFP associated with this Cosmic Tag: " << cosmicTagToPFP_v.size() << std::endl;

    if (score_cut < 0) {
      pfpTaggedOut.emplace_back(cosmicTagToPFP_v.at(0));
      tagid_v.emplace_back(cosmicTag->CosmicType());
    } else {
      if (cosmicTag->CosmicScore() > score_cut) {
        pfpTaggedOut.emplace_back(cosmicTagToPFP_v.at(0));
        tagid_v.emplace_back(cosmicTag->CosmicType());
      }
    }
  }

}
//____________________________________________________________________________________________
int UBXSec::PFPInCommon(lar_pandora::PFParticleVector first, lar_pandora::PFParticleVector second){

  int nInCommon = 0;


  for (unsigned int f = 0; f < first.size(); f++){

    for (unsigned int s = 0; s < second.size(); s++){

      if(first.at(f) == second.at(s)) {

        nInCommon++;

        //std::cout << "In common found, flash is  " << flash << " and geo is " << geo << std::endl;

      }
    }

  }
  return nInCommon;

}



DEFINE_ART_MODULE(UBXSec)
