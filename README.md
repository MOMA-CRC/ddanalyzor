# ddanalyzor: Simple analysis of raw droplet digital PCR fluorescence data using R

## What
To calculate the concentration of target DNA in a droplet digital PCR (ddPCR)
experiment one need to set a threshold between positive and negative droplets.
This is often done subjectively by "eyeballing", or using an automatic threshold
setting option in the software.
**ddanalyzor.R** is a simple script that uses base R to turn raw ddPCR
data into concentration when the above approaches are not acceptable.

**ddanalyzor.R** uses a smoothing kernel on a positive control sample to define
positive and negative droplets and a fluorescence threshold to separate these.
This threshold is then used on unknown samples - optionally after calibrating samples
to similar median negative fluorescence - in order to calculate copies per well.


## How
1. **Export raw fluorescence data**
2. **Read raw data into named R lists**
3. **Run ddanalyzor for each ddPCR plate**


## Example
The starting point for **ddanalyzor** is a (named) list of numeric vectors of raw
fluorescence data. Run **ddanalyzor** for each plate and each fluorescence channel/fluorophore.
For instance, say you have exported raw data from the QuantaSoft software
(`options > export amplitude and cluster data`) from a single plate to csv files
into the directory _raw_data_. If each file is name _A01_Amplitude.csv_, _A02_Amplitude.csv_ etc you can do

```r

#The raw data directory
raw_data_dir <- "raw_data"

#Find all Amplitude.csv files in raw data directory
raw_data_files <- list.files(raw_data_dir, pattern = "Amplitude.csv", full.names = T)

#Then read into list of data frames (note QuantaSoft format)
raw_data <- lapply(raw_data_files,
                   read.table,
                   header = T,
                   stringsAsFactors = F)

# Name raw data "A01", "A02" etc.
names(raw_data) <- sub("_Amplitude.csv", "", basename(raw_data_files))

# Since ddanalyzor analyzes in a plate- and channel/target-wise manner extract first column (channel 1)
ch1 <- lapply(raw_data, "[[", 1)

# This is a list of numeric vectors
str(ch1)

# Now, run ddanalyzor with the positive control in "A01" expecting one negative and one positive populations
source("ddanalyzor.R")
results <-
  ddanalyzor(cases = raw_data[!names(raw_data) %in% "A01"],
             positive_control = raw_data["A01"],
             extremes = c(1,1))
head(results)             

```

A minimal example is also described in [example.R][./example.R].


## ddanalyzor parameters
### Input
* cases. A list of numeric vectors of raw data for each well for a single channel. If named, these will be used in the results.  
* positive_control. A list of raw data for the control well for a single channel. If multiple positive controls they should be combined into one vector. If a vector of length 1 is given this is used as a threshold in which case extremes, calibrate, beta, bandwidths and n are not used.
* extremes. Vector of length 2 given the number of dominant negative and positive fluorescence populations expected. I.e. for an assay with a fluorescence profile of two negative populations and a single positive population one should set this to c(2,1).   
* calibrate [TRUE]. If TRUE, data is linearly calibrated according to the positive control so that the median of the negative droplets is similar to the negative doplets of the positive control.
* output [NULL]. If NULL, results are returned as a data frame. Otherwise results are written to the file specified in output.
* beta [0]. Stringency parameter where 0 includes all data after threshold point.
* alpha [0.05]. Confidence interval alpha.
* bandwidths [seq(1, 5000, 50)]. Vector of bandwidths (fed to R's base::density) to search through from left to right to find extremes. The smallest bandwidth that results in the number of maxima satisfying the sum of the extremes argument is used.
* n [10000]. The number of equally spaced points to used for kernel estimator.
* seed [NULL]. For reproducibility.

### Output
A data frame with results from each well (when output is NULL). A table is written
to a file if specified by output. The columns are: "sample_id", "positive_count",
"negative_count", "copies_per_well", "copies_per_well_ci_low", "copies_per_well_ci_high",
"threshold", "beta", "alpha", "extreme_profile", "calibrate", "calibrate_factor",
"kernel_bandwidth", "control_maxima", "control_minima" and "time_stamp".


## Cite
**ddanalyzor** was used in the paper [Novel DNA methylation biomarkers show high
sensitivity and specificity for blood-based detection of colorectal cancer — a
clinical biomarker discovery and validation study](https://doi.org/10.1186/s13148-019-0757-3).
Please cite this.
