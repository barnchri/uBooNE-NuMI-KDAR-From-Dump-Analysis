# uBoone-NuMI-KDAR-From-Dump-Analysis
This is the analysis code used for my original Ph.D. thesis measurement, the differential cross section of Kaon-Decay-At-Rest Muon Neutrinos from the NuMI Beam Dump with MicroBooNE.

This analysis involves a pre-selection, reconstruction after preselection, a BDT, and data visualization with data normalization.  Each of the files fulfills the following purpose(s):

UBXSec_module.cc - This file performs the pre-selection, removing events with binary selection requirements, while extracting quantities necessary for the proton length & direction reconstruction.

Updates_To_Proton_Fitting_Technique_With_All_MCWeights.C - This file passes along quantities needed for the BDT & performs proton length & direction reconstruction.

Generating_Signal_And_Background_Trees*.C - These files prepare the ROOT TTrees for the BDT for each of the samples.

TMVAClassification.C - This file trains the BDT.

TMVAClassificationApplication*.C - These files test the BDT on each of the samples, providing a score in the range [-1.0, 1.0] based on how confident the BDT is that the given event is a KDAR-from-the-dump muon neutrino event.

Plotting_Quantities_As_Function_Of_Cut_Muon_Candidate_Only.C - This file normalizes & visualizes the output for the testing samples as a function of BDT score cut.

counting_weighted_number_of_events.py - This script uses the weights in the tree to count the weighted number of events passing the pre-selection.
