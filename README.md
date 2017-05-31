# overlappingspikes
Neural spike sorting algorithms for overlapping spikes in MATLAB. Completed as a senior thesis project, advised by Professor Alexander H. Barnett. 

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
Sorting spikes without time shifts.

*Functions:*

*Drivers:*

## timeshifts
Soring spikes with time shifts. These functions can also handle sorting spikes without time shifts.

*Functions:*

*Drivers:*

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
- **terrorbar.m:** written by T. Agus (downloaded from: https://www.mathworks.com/matlabcentral/fileexchange/52367-terrorbar-m--error-bars-with-controlled-widths-post-r2014b) used in paulabars.m to address the issue of non-standard error bar caps pre-2016b MATLAB versions
