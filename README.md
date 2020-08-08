# uBoone-NuMI-KDAR-From-Dump-Analysis
This is the analysis code used for my now-defunct thesis measurement, the differential cross section of Kaon-Decay-At-Rest Muon Neutrinos from the NuMI Beam Dump with MicroBooNE.

This analysis involves a pre-selection, reconstruction after preselection, a BDT, and data visualization while properly normalizing the data.  Each of the files fulfills the following purposes:

UBXSec_module.cc - This file performs the pre-selection while extracting quantities necessary for the proton length & direction reconstruction.

Updates_To_Proton_Fitting_Technique_With_All_MCWeights.C - This file passes along quantities needed for the BDT & performs proton length & direction reconstruction.

Generating_Signal_And_Background_Trees*.C - These files prepare the trees for the BDT for each of the samples.

TMVAClassification.C - This file train the BDT.

TMVAClassificationApplication*.C - These files tests the BDT on each of the samples.

Plotting_Quantities_As_Function_Of_Cut_Muon_Candidate_Only.C - This file properly normalizes & visualizes the output for the testing samples as a function of BDT score cut.
