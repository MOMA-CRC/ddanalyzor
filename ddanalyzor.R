###-------------------###
### JENSEN ET AL 2019 ###
###-------------------###


ddanalyzor <- function(cases, positive_control, extremes, calibrate = T, output = NULL,
                       beta = 0, alpha = 0.05,
                       bandwidths = seq(1, 5000, 50), n = 10000, seed = NULL) {
  
  #INPUT:
  #cases: List of fluorescence data for each wells in a single plate for a single channel
  #where each element is a numeric vector named.
  #positive_control: Named list of fluorescence data for positive control well(s). Can take values:
  #1) A list of one numeric vector (if more than one control was used combine these into one vector)
  #2) A single numeric value which is then used as cut-point defining positive and negative droplets.
  #When a manually set threshold is used (not recommended) extremes, calibrate, beta, bandwidths and n is not used.
  #extremes: Numeric of length 2 designating the number of negative and positive extremes.
  #to expect in control data (e.g. 1 negative and 1 positive peak is c(1,1) ).
  #calibrate: (default True) Data is linearly calibrated according to the positive control
  #so that the median of negative droplets is similar to positive control.
  #output: Save a file with plate results for all wells. If NULL (default) return data frame.
  #beta:
  #alpha:
  #bandwidths: Vector of bandwidths to search through from left to right to find extremes.
  #n: The number of qually spaced points to used for kernel estimator.
  #seed:
  #
  #
  #OUTPUT:
  #A table with the following columns:
  #1) Calibrated data for the channel and ctrl well and for the given 2) band width, 
  #3) mode for negative and 4) positive droplets, and 5) the minimum. 6) channel,
  #7) ctrl well, 8) analysis date and 9) raw input plate data
  
  
  
  ### FUNCTIONS ###
  
  droplets <- function(N = 10000, n = c(0.9, 0.1, 0.1), maxima = c(2000, 8000), sd = c(500, 500), seed=1) {
    #Generate vector of totally N droplets distributed accoring to n fractions where 
    #the last index is rain - around the maxima with sd sd. length(n) - length(maxima) = 1.
    # #Example: 
    # N <- 10000 #droplets
    # n <- c(0.9, 0.1, 0.1) #relative distributions in two populations plus rain
    # maxima <- c(2000, 8000) #Mean for the two populations (rain is placed from 0-8000)
    # sd <- c(500, 500)
    seed = seed
    drops <- runif(N * n[length(n)], 0, max(maxima)) #Initiate with random uniform rain from 0 to max
    for(p in 1:(length(n)-1) ) {
      drops <- c(drops, rnorm(N*n[p], maxima[p], sd[p]))
    } 
    sample(drops)
  }
  
  time_stamp <- function(){
    #Return Sys.time() in the form yyyymmddhhmmss
    gsub(" ", "-", gsub(":|-| ", "", paste(Sys.time(), collapse = " ")) ) 
  }
  
  extremator <- function(data, bw, n) {
    ### Find extremes using density kernel for a given bw ###
    # INPUT:
    # data: A vector of data to make kernel density estimate.
    # bw: Smoothing bandwitdh.
    # OUTPUT: 
    # A list of two matrices. The first gives (x,y) of the estimated maxima, the
    # second (x,y) the estimated minima. 
    
    d <- density(data, bw = bw, n = n)
    
    i_max <- which(diff(sign(diff(d$y))) == -2)+1
    i_min <- which(diff(sign(diff(d$y))) == 2)+1
    
    list(t(mapply(c, d$x[i_max], d$y[i_max], SIMPLIFY=TRUE)),
         t(mapply(c, d$x[i_min], d$y[i_min], SIMPLIFY=TRUE))
    )
  }
  
  pois_conc <- function(pos, neg, a = 0.05) {
    p <- pos/(pos+neg) #Fraction of positives, chance of 'hit'
    p_ci <- p + c(qnorm(a/2), qnorm(1-a/2)) * sqrt(p*(1-p)/(pos+neg)) 
    lamda <- -log(1- c(p,p_ci))  #cps per droplet
    return(lamda*(pos+neg) ) #cps in well (+ ci_low and ci_high)
  }
  
  threshold_splitter <- function(v, threshold) {
    list('pos' = v[v >= threshold],
         'neg' = v[v < threshold])
  }
  
  
  
  ### READ DATA ###
  options(warn = -1)
  
  positive_control_name <- names(positive_control)
  positive_control <- unlist(positive_control, use.names = F)
  
  if(length(positive_control) > 1) {
    if(extremes[1] < 1 | extremes[2] < 1) {
      cat("[ddanalyzor] In extremes(a, b), both a and b should be an integer > 0. Aborting.\n")
      return(NULL)
    }
    
    
    ### DEFINE CUT POINT IN POSITIVE CONTROL ###
    extreme_world <- lapply(bandwidths, function(x) {
      extremator(data = positive_control, bw = x, n = n)
    })
    
    #The smallest (first) bandwidth to estimate sum(maxima), number of maxima
    i <- which(sapply(sapply(extreme_world, '[[', 1), nrow) == sum(extremes))[1]
    
    if(is.na(i)) {
      cat("[ddanalyzor] Failed to find expected extreme profile. Aborting.\n")
      return(NULL)
    }
    
    #Then use this bw to get a temporary cut point between pos and neg droplets 
    #in the control - this is actually a repeat for that bw of the above.
    e0 <- extremator(data = positive_control, bw = bandwidths[i], n = n)
    
    #The X coordinate of the first minimum after the negative maxima in extremes[1]
    cut0 <- e0[[2]][ extremes[1], 1 ]
    
    ### CALIBRATE CASES ### 
    
    if(calibrate) {
      #Median of negative droplets of positive control sample
      neg0_median <- median(positive_control[ positive_control <  cut0 ])
      
      #For cases substract this to calibrate to similar median of negatives as positive control
      cal_factor <- sapply(cases, function(x){median(x[x < cut0]) - neg0_median})
    } else {
      cal_factor <- rep(0, length(cases))
    }
    
    cases_calibrated <- list(NULL)
    for(j in 1:length(cases)) {  
      cases_calibrated[[j]] <- cases[[i]] - cal_factor[[j]]
    }
    names(cases_calibrated) <- names(cases)
    
    
    ### BOOTSTRAP THRESHOLD IN POSITIVE CONTROL ###
    
    #Extract the first positive population p_bar after last minimum (cut0)
    if (extremes[2] == 1) { #If only one positive maximum use every droplet above cut
      cut1 <- Inf
    } else {
      cut1 <- e0[[2]][extremes[1]+1, 1]
    }
    p_bar <- positive_control[ positive_control >  cut0 & positive_control <= cut1 ]
    
    #Find the beta/2-th quantile from 5000 bootstrap pouplations and use mean as threshold
    if(!is.null(seed)) set.seed(seed)
    x_star <- matrix(sample(p_bar, length(p_bar)*5000, replace = T), ncol = 5000)
    threshold <- round(mean(apply(x_star, 2, function(x)quantile(x, beta/2))))
    
    
    ### CALCULATE CONCENTRATION ###
    
    #Split all wells into pos and neg given threshold, and use pois_conc
    #to add cps per well and low/high ci for this (given alpha a)
    
    cases <- c(list(positive_control), cases)
    names(cases)[1] <- positive_control_name
    
    counts <- t(sapply(lapply(lapply(cases, threshold_splitter, threshold), "[", 1:2), lengths))
    counts <- as.data.frame(cbind(counts, t(apply(counts, 1, function(x)pois_conc(x[1], x[2], a = alpha)))), stringsAsFactors = F)
    names(counts) <- c('positive_count', 'negative_count', 'copies_per_well', 
                       'copies_per_well_ci_low', 'copies_per_well_ci_high')
    counts$sample_id <- rownames(counts)
    counts$threshold <- threshold
    counts$beta <- beta
    counts$alpha <- alpha
    counts$extreme_profile <- paste(extremes, collapse = ",")
    counts$calibrate <- calibrate
    counts$calibrate_factor <- c(0, cal_factor)
    counts$kernel_bandwidth <- bandwidths[i]
    counts$control_maxima <- paste(apply(format(e0[[1]], scientific = T, digits = 4),
                                         1, paste, collapse = ","), collapse = "|")
    counts$control_minima <- paste(apply(format(e0[[2]], scientific = T, digits = 4),
                                         1, paste, collapse = ","), collapse = "|")
    counts$time_stamp <- time_stamp()
    
    
  }
  
  if(length(positive_control) == 1) {
    # A single value given, interpret as manual threshold
    threshold <- positive_control[[1]]
    counts <- t(sapply(lapply(lapply(cases, threshold_splitter, threshold), "[", 1:2), lengths))
    counts <- as.data.frame(cbind(counts, t(apply(counts, 1, function(x)pois_conc(x[1], x[2], a = alpha)))), stringsAsFactors = F)
    names(counts) <- c('positive_count', 'negative_count', 'copies_per_well', 
                       'copies_per_well_ci_low', 'copies_per_well_ci_high')
    counts$sample_id <- rownames(counts)
    counts$threshold <- threshold
    counts$beta <- NA
    counts$alpha <- alpha
    counts$extreme_profile <- NA
    counts$calibrate <- NA
    counts$calibrate_factor <- NA
    counts$kernel_bandwidth <- NA
    counts$control_maxima <- NA
    counts$control_minima <- NA
    counts$time_stamp <- time_stamp()
    
  }
  
  #Reorder columns
  counts <- counts[, c('sample_id', 'positive_count', 'negative_count', 'copies_per_well', 
                       'copies_per_well_ci_low', 'copies_per_well_ci_high',
                       'threshold', 'beta', 'alpha', 'extreme_profile', 'calibrate', 'calibrate_factor', 
                       'kernel_bandwidth', 'control_maxima', 'control_minima', 'time_stamp')]
  
  ### OUTPUT RESULT TABLE ###
  if(is.null(output)) {
    return(counts)
  } else {
    write.table(counts, output, row.names = F)
  }
}


# ### EXAMPLE ###
#Make raw data as a named list of named numeric vectors representing fluoresence
#from individual ddPCR wells for a single channel.
#Let "A01" be the the positive control sample

seed = 8
well_names <- apply(expand.grid(LETTERS[1:8], sprintf("%02d", 1:12)), 1, paste, collapse = "")
N = 10000 #Approximate no of droplets per well
maxima = c(2000, 8000) #Assume maxima around 2000 (negatives) and 8000 (positive)
sd = c(100, 150) #Standard deviation of maxima

raw_data <- list()
for(w in well_names) {
  if(w == "A01") {
    #Pos control
    raw_data[[w]] <- droplets(N, n = c(0.6, 0.4, 0.05), maxima = maxima, sd = sd, seed = seed + 1)
  } else {
    #Cases. Randomly make few, very few and none positives
    raw_data[[w]] <- droplets(N, n = c(0.99, sample(c(0.1, 0.01, 0), 1), 0.05), sd = sd, seed = seed + 1)
  }
}

ddanalyzor(cases = raw_data[!names(raw_data) %in% "A01"],
           positive_control = raw_data["A01"],
           extremes = c(1,1))
