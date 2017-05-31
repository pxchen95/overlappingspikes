# overlappingspikes
Neural spike sorting algorithms for overlapping spikes in MATLAB. Completed as a senior thesis project, advised by Professor Alexander H. Barnett. 

## synthesis
Creation of synthetic spike signals. 
*Functions:*
  **synth.m:** synthesizes spike signals (without noise)
  **gen_firing.m:** generates firing times using a Poisson distribution
*Drivers:*
  **demo_synth.m:** tests synth.m, gen_firing.m
  
## detectspikes
Spike detection via maximum likelihood.
*Drivers:*
  **binarytest.m:** detection of one vs. no spike using linear classifiers
  
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
*Functions:*
*Drivers:*
