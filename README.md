# ddanalyzor: Simple analysis of raw droplet digital PCR fluorescence data using R

## What
To calculate the concentration of target DNA in a droplet digital PCR (ddPCR)
experiments one need to set a threshold between positive and negative droplets.
This is often done subjectively by "eyeballing", or using an automatic threshold
setting options in the software.
**ddanalyzor.R** is a simple script that uses base R to turn raw ddPCR
data into concentration when these approaches are not acceptable.

**ddanalyzor.R** uses a smoothing kernel on a positive control sample to define
positive and negative droplets and a fluorescence threshold to separate these.
This threshold is used on unknown samples - optionally after calibrating samples
to similar median negative fluorescence - and hence to calculate copies per well.


## How
1. **Export raw fluorescence data**
2. **Read raw data into named R lists**
3. **Run ddanalyzor for each ddPCR plate**

## ddanalyzor parameters
### Input
_cases_
## Cite
**ddanalyzor** was used in the paper [Novel DNA methylation biomarkers show high
sensitivity and specificity for blood-based detection of colorectal cancer — a
clinical biomarker discovery and validation study](https://doi.org/10.1186/s13148-019-0757-3).
Please cite this.
