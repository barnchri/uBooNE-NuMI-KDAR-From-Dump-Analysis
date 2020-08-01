#include <cstdlib>
#include <iostream>
#include <map>
#include <string>
#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"

#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"

int TMVAClassification( TString myMethodList = "" ) {
  // This loads the library
  TMVA::Tools::Instance();

  // Default MVA methods to be trained + tested
  std::map<std::string,int> Use;

  // --- Boosted Decision Trees
  Use["BDT"]             = 0; // uses Adaptive Boost
  Use["BDTG"]            = 1; // uses Gradient Boost
  Use["BDTB"]            = 0; // uses Bagging
  Use["BDTD"]            = 0; // decorrelation + Adaptive Boost
  Use["BDTF"]            = 0; // allow usage of fisher discriminant for node splitting
  std::cout << std::endl;
  std::cout << "==> Start TMVAClassification" << std::endl;

  // Select methods (don't look at this code - not of interest)
  if (myMethodList != "") {

    for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) it->second = 0;
    std::vector<TString> mlist = TMVA::gTools().SplitString( myMethodList, ',' );
    
    for (UInt_t i=0; i<mlist.size(); i++) {
      std::string regMethod(mlist[i]);
      if (Use.find(regMethod) == Use.end()) {
	std::cout << "Method \"" << regMethod << "\" not known in TMVA under this name. Choose among the following:" << std::endl;
	for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) std::cout << it->first << " ";
	std::cout << std::endl;
	return 1;
      }
      Use[regMethod] = 1;
    }
  }

  // --- Here the preparation phase begins
  
  // Create a ROOT output file where TMVA will store ntuples, histograms, etc.
  TString outfileName( "TMVA.root" );
  TFile* outputFile = TFile::Open( outfileName, "RECREATE" );

  // The second argument is the output file for the training results
  //  TMVA::Factory *factory = new TMVA::Factory( "TMVAClassification_minus", outputFile,"!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P,D:AnalysisType=Classification" );
  TMVA::Factory *factory = new TMVA::Factory( "TMVAClassification_minus", outputFile,"!V:!Silent:Color:DrawProgressBar:Transformations=I:AnalysisType=Classification" );

  TMVA::DataLoader *dataloader = new TMVA::DataLoader("dataset_calibrated_all_stats_no_neutrino_direction");
  // factory options pag. 15 TMVA userg guide
  // branch, name, unit, type

  dataloader->AddVariable("reco_vtx_x", "Reco Vertex x", "", 'F');
  dataloader->AddVariable("reco_vtx_y", "Reco Vertex y", "", 'F');
  dataloader->AddVariable("reco_vtx_z", "Reco Vertex z", "", 'F');
  dataloader->AddVariable("x_component_of_muon_momentum", "Absolute Value of X-Component of Momentum", "", 'F');
  dataloader->AddVariable("y_component_of_muon_momentum", "Absolute Value of Y-Component of Momentum", "", 'F');
  dataloader->AddVariable("z_component_of_muon_momentum", "Absolute Value of Z-Component of Momentum", "", 'F');
  dataloader->AddVariable("flash_z_difference", "Z Flash Diff.", "", 'F');
  dataloader->AddVariable("flash_y_difference", "Y Flash Diff.", "", 'F');
  dataloader->AddVariable("truncated_dQdx_of_TPCObject_muon_candidate_u_plane", "U Trunc. Mean", "", 'F');
  dataloader->AddVariable("truncated_dQdx_of_TPCObject_muon_candidate_v_plane", "V Trunc. Mean", "", 'F');
  dataloader->AddVariable("truncated_dQdx_of_TPCObject_muon_candidate_y_plane", "Y Trunc. Mean", "", 'F');
  dataloader->AddVariable("truncated_dQdx_of_TPCObject_muon_candidate_sum", "Sum Trunc. Mean", "", 'F');
  dataloader->AddVariable("hit_sum", "Total Hit Sum", "", 'F');
  dataloader->AddVariable("u_plane_hit_sum",  "U Hit Sum", "", 'F');
  dataloader->AddVariable("v_plane_hit_sum",  "V Hit Sum", "", 'F');
  dataloader->AddVariable("y_plane_hit_sum",  "Y Hit Sum", "", 'F');
  dataloader->AddVariable("vicinity_hit_sum_difference", "Difference Between Hits in Vicinity and Hits on Tracks: All Planes", "", 'F');
  dataloader->AddVariable("u_plane_vicinity_hit_sum_difference", "Difference Between Hits in Vicinity and Hits on Tracks: U Plane", "", 'F');
  dataloader->AddVariable("v_plane_vicinity_hit_sum_difference", "Difference Between Hits in Vicinity and Hits on Tracks: V Plane", "", 'F');
  dataloader->AddVariable("y_plane_vicinity_hit_sum_difference", "Difference Between Hits in Vicinity and Hits on Tracks: Y Plane", "", 'F');
  dataloader->AddVariable("total_number_of_slice_hits_in_vicinity_ADCs", "Total Number of Hits in Muon Candidate Vicinity", "", 'F');
  dataloader->AddVariable("u_plane_number_of_slice_hits_in_vicinity_ADCs", "U-Plane Number of Hits in Muon Candidate Vicinity", "", 'F');
  dataloader->AddVariable("v_plane_number_of_slice_hits_in_vicinity_ADCs", "V-Plane Number of Hits in Muon Candidate Vicinity", "", 'F');
  dataloader->AddVariable("y_plane_number_of_slice_hits_in_vicinity_ADCs", "Y-Plane Number of Hits in Muon Candidate Vicinity", "", 'F');
  
  // Read training and test data
  TString fname = "./input_file_cv_training_sample.root";  // <------------------------------------------ FILE .root WITH TTREE IN FOLDER
   
  if (gSystem->AccessPathName( fname ))  // file does not exist in local directory
    gSystem->Exec("curl -O http://root.cern.ch/files/tmva_class_example.root");
   
  TFile *input = TFile::Open( fname );
   
  std::cout << "--- TMVAClassification       : Using input file: " << input->GetName() << std::endl;
   
  // --- Register the training and test trees
  TTree *signal     = (TTree*)input->Get("signal_tree");         // <-----------  TTREE IN FILE
  TTree *background = (TTree*)input->Get("background_tree");     // <-----------  TTREE IN FILE
  
  // Global event weights per tree
  Double_t signalWeight     = 1.0;
  Double_t backgroundWeight = 1.0; //0.1;
   
  // You can add an arbitrary number of signal or background trees
  dataloader->AddSignalTree(signal, signalWeight);

  if(background->GetEntries()!=0) dataloader->AddBackgroundTree(background, backgroundWeight);

  // Set the weights for the signal and the background.
  dataloader->SetSignalWeightExpression("spline_fix_mcweight * rootino_fix_mcweight * central_value_mcweight * special_kdar_weight");
  dataloader->SetBackgroundWeightExpression("spline_fix_mcweight * rootino_fix_mcweight * central_value_mcweight * special_kdar_weight * sampleweight * treeweight");

  // Apply additional cuts on the signal and background samples (can be different)
  TCut mycuts = ""; // for example: TCut mycuts = "abs(var1)<0.5 && abs(var2-0.5)<1";
  TCut mycutb = ""; // for example: TCut mycutb = "abs(var1)<0.5";

  dataloader->PrepareTrainingAndTestTree(mycuts, mycutb, "SplitMode=random:!V" );

  // Boosted Decision Trees
  if (Use["BDTG"]) // Gradient Boost
    factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDTG",
			 "!H:!V:NTrees=1000:MinNodeSize=2.5%:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=2" );

  // Train MVAs using the set of training events
  factory->TrainAllMethods();

  // ---- Evaluate all MVAs using the set of test events
  factory->TestAllMethods();

  // ----- Evaluate and compare performance of all configured MVAs
  factory->EvaluateAllMethods();

  // --------------------------------------------------------------

  // Save the output
  outputFile->Close();

  std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
  std::cout << "==> TMVAClassification is done!" << std::endl;

  delete factory;

  // Launch the GUI for the root macros
  if (!gROOT->IsBatch()) TMVA::TMVAGui( outfileName );
  return 0;
}

int main( int argc, char** argv ) {

  // Select methods (don't look at this code - not of interest)
  TString methodList; 

  for (int i=1; i<argc; i++) {
    TString regMethod(argv[i]);
    if(regMethod=="-b" || regMethod=="--batch") continue;
    if (!methodList.IsNull()) methodList += TString(","); 
    methodList += regMethod;
  }

  return TMVAClassification_edited(methodList); 

}
