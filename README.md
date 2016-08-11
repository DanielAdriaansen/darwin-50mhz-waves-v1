# 50mhzwaves

Step 1:
run mask_50MHz_profiler.ncl
This creates a precipitation filter using the 920 data.

Step 2:
run mask_min_profiler.ncl
This filters based on bad 50 data and also the 920 precip filter from step 1 as
well as a user defined minimum period length for spectral analysis.

DATA:
data/masked = raw + masked 50 MHz data using precip info from 920. Also contains mask_x arrays for mask info.
data/maskedmin = raw + masked 50 MHz data using precip info + a minimum valid period of data.
data/maskedminbad = raw + masked 50 MHz data using pecip info + bad 50 MHz data + a minimum valid period of data.

IMAGES:
** images/data = 920/50 four panel with precip mask applied
** images/mask = 920/50 four panel of PRECIP/BAD/GOOD three value 2D time/height plots
images/periods = special plots visualizing data quality
images/prime = 920/50 four panel perturbation with precip mask applied --> THESE NEED WORK
** images/raw = raw data from Christopher
** images/stdata = 920/50 four panel with precip, and minimum masks applied (input to S-transform) but not perturbations
** images/stdata2 = 920/50 four panel with precip, bad, and minimum masks applied (input to S-transform) but not perturbations
