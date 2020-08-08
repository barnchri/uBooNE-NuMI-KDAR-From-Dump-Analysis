import ROOT
import root_numpy
import pandas as pd

file    = ROOT.TFile('KDAR_Training_Sample_For_BDT.root')
tree    = file.Get('t0ana/tree')

arr     = root_numpy.tree2array(tree, branches = ['run', 'subrun', 'event', 'truth_energy_in_tree', 'spline_fix_mcweight', 'rootino_fix_mcweight', 'central_value_mcweight'], selection = 'flash_PEs < 2000 && sum_of_TPCObject_track_lengths < 65')
df      = pd.DataFrame(arr)

num_of_unweighted_events                  = 0.
total_weighted_events_without_KDAR_weight = 0.
total_weighted_events_with_KDAR_weight    = 0.

already_used_event = False

already_used_runs    = []
already_used_subruns = []
already_used_events  = []

for i in range(df.shape[0]):

    print "Looping over event #%d." %( i )

    already_used_event = False

    # Make sure that the event we're using is unique.                                                                                                                                                    
    for j in range( len(already_used_runs) ):

        if already_used_runs[j] == df['run'][i] and already_used_subruns[j] == df['subrun'][i] and already_used_events[j] == df['event'][i]:
            already_used_event = True
            break

    if already_used_event == True:
        continue

    total_weighted_events_without_KDAR_weight += ( df['spline_fix_mcweight'][i] * df['rootino_fix_mcweight'][i] * df['central_value_mcweight'][i] )

    # Add an additional factor of 5 if the event in question is a KDAR event (whether it's from the NuMI dump or not).
    if abs( df['truth_energy_in_tree'][i] - 235.532 ) < 0.001:
        total_weighted_events_with_KDAR_weight    += ( df['spline_fix_mcweight'][i] * df['rootino_fix_mcweight'][i] * df['central_value_mcweight'][i] * 5.0 )

    else:
        total_weighted_events_with_KDAR_weight    += ( df['spline_fix_mcweight'][i] * df['rootino_fix_mcweight'][i] * df['central_value_mcweight'][i] )
        
    num_of_unweighted_events                  += 1.
    
print "The number of unweighted events = %f." %( num_of_unweighted_events )
print "The number of weighted events without the factor of 5 = %f." %( total_weighted_events_without_KDAR_weight )
print "The number of weighted events with the factor of 5 = %f." %( total_weighted_events_with_KDAR_weight )

