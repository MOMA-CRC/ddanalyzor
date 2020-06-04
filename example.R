
### MINIMAL EXAMPLE ###

#Make raw data as a named list of named numeric vectors representing fluoresence
#from individual ddPCR wells for a single channel.
#Let "A01" be the the positive control sample

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

source("ddanalyzor.R")

results <- 
  ddanalyzor(cases = raw_data[!names(raw_data) %in% "A01"],
             positive_control = raw_data["A01"],
             extremes = c(1,1))