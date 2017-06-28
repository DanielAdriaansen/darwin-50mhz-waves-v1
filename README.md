# 50mhzwavesV1

Step 1:
run mask_50MHz_profiler.ncl
This creates a precipitation mask using the 920 data (only uses test 1) and also maska out "bad" 50 MHz data
Mask key: 1 = good data, 2 = missing/bad 50 MHz data, 3 = precipitation
Reads from: /d1/dadriaan/paper/data/c2/raw
Writes to: /d1/dadriaan/paper/data/c2/masked
Processes: Individual (daily) files

Step 2:
run mask_min_profiler.ncl
This creates an additional mask for periods that are less than a user defined threshold (minutes)
Mask key: 1 = good data, 2 = missing/bad 50 MHz data, 3 = precipitation, 4 = too short
Reads from: /d1/dadriaan/paper/data/c2/masked
Writes to: /d1/dadriaan/paper/data/c2/maskedminbad
Processes: Individual (daily) files
****NOTE: This should probably be re-written in MATLAB since we'll need to do this after finding outliers****

Step 3:
run profiler_driver.m
1. Concatenate all daily files into a single array
2. Save variables to a mat file
3. Mask out the regime we're not processing
4. Find outliers in the data and mask them out
5. Re-run the chunk finding code from step 2 above


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
** images/st_monsoon = S-Transform of all the chunks in the monsoon regime. This is just ST-w, integrated, absolute, and
   square of the output
** images/st_break = S-Transform of all the chunks in the break regime. This is just ST-w, integrated, absolute, and
** images/summary = summary images showing the number of points at each time and height in the Darwin day for monsoon/break
   regime available for analysis.
** images/anomalies = images used for showing Christopher anomalous data
** images/gretchenmullendore = files shared with Gretchen for meetingsi
