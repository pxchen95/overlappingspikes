# overlappingspikes
Neural spike sorting algorithms for overlapping spikes. Algorithms are all written in MATLAB. Completed as a senior thesis project, advised by Professor Alexander H. Barnett. 

- **maindriver.m:** sorts spikes with and without penalty, with and without time shifts
- **decision_bound.m:** plots the decision boundaries in the 2-spike case for all algorithms in **notimeshifts**

## synthesis
Creation of synthetic spike signals. 

*Functions:*
- **synth.m:** synthesizes spike signals (without noise)
- **gen_firing.m:** generates firing times using a Poisson distribution
  
*Drivers:*
- **demo_synth.m:** tests synth.m, gen_firing.m
  
## detectspikes
Spike detection via maximum likelihood.

*Drivers:*
- **binarytest.m:** detection of one vs. no spike using linear classifiers
  
## notimeshifts
Sorting spikes without time shifts. These functions cannot handle sorting spikes with time shifts.

*Functions:*
- **brute_force.m:** brute force, NO time shifts, NO penalty
- **greedy_sort.m:** simple greedy, NO time shifts, NO penalty
- **FoBa.m:** forward-backward greedy algorithm, currently is no different than greedy_sort.m (in performance)
- **count.m:** counts the errors for sorting spikes WITHOUT time shifts

*Drivers:*
- **driver_count.m:** tests every possible case for count.m
- **driver_brute.m:** shows that brute-force gets 100% accuracy at eta = 0
- **driver_foba.m:** compares sorting by simple greedy and foba, shows that error fractions match exactly

## timeshifts
Soring spikes with time shifts. These functions can also handle sorting spikes without time shifts.

*Functions:*
- **brute_force_ts.m:** brute-force for time shifts with penalty, tests all possible combinations of time-shifted spikes
- **combine.m:** used in brute_force_ts.m to find all possible combinations of the elements of n vectors (adapted from http://stackoverflow.com/questions/21895335/generate-a-matrix-containing-all-combinations-of-elements-taken-from-n-vectors)
- **greedy_likely.m:** simple greedy with penalty, for time-shifted spikes
- **greedy_pairs.m:** greedy with pairs with penalty, for time-shifted spikes
- **count_ts2.m:** counts the errors for sorting time-shifted spikes

*Drivers:*
- **driver_count_ts2:** tests every possible case for count_ts2.m

## realdata
Run algorithms on data from Harris, et al. "Accuracy of tetrode spike separation as determined by simultaneous intracellular and extracellular measurements."
- **harrisdata.mat:** clips resulting from getharrisdata.m, ready to be sorted

*Drivers:*
- **getharrisdata.m:** written by A. H. Barnett, filters harrisdata, extracts clips, uses clustering algotihm to sort clips
- **driver_harris.m:** sorts the real clips from harrisdata.mat
- **driver_harris2.m:** creates and sorts synthetic clips using the spike templates (with time shifts) from harrisdata.mat

## other
- **paulabars.m:** written by A. H. Barnett, creates bar graphs for error fractions with error bars (std, /sqrt(N))
- **paulabars.m:** version 2 of paulabars.m, changes are only in formatting
- **terrorbar.m:** written by T. Agus (downloaded from: https://www.mathworks.com/matlabcentral/fileexchange/52367-terrorbar-m--error-bars-with-controlled-widths-post-r2014b), used in paulabars.m to address the issue of non-standard error bar caps in pre-2016b MATLAB versions
