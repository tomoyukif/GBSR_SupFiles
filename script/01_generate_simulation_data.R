## This document was prepared as supplemental text file 
## accompanied with the publication of GBScleanR.
##
## Creadted on Dec. 20, 2021
## Created by Tomoyuki Furuta
## Institute of Plant Science and Resources, 
## Okayama University

## We developed a package "SimPop" that simulates a population 
## with given pedigree and rondomly assigned genotypes and 
## read counts.
## The "SimPop" package is available at the following:
## "https://github.com/tomoyukif/SimPop"
## This package is developed just for lab-use and 
## not supposed to be reused by other users.
## Therefore, we do not guarantee that the tool
## works finely in other computational environments 
## with other settings that we did not use in our paper.


########################################################
# 2 way-F2 from inbred Parents with no allele read bias#
########################################################
dir <- "../noError_homoP2_F2"
if(!dir.exists(dir)){ dir.create(dir) }
#######################
# Simulate population
library(SimPop)
xo_freq <- 2
pos <- 50
n_progeny <- c(10, 100, 1000)
n_mar <- c(50, 100, 1000)
total_read <- c(0.1, 0.25, 0.5, 0.75, 1, 2, 3)
parental_read <- c(1, 3)
condition <- expand.grid(xo_freq, pos, n_progeny, n_mar)
names(condition) <- c("xo_freq", "pos", "n_progeny", "n_mar")
missing_rate <- NULL
for(i in 1:nrow(condition)){
  pop <- makeSimPop(n_mar = condition$n_mar[i],
                    xo_freq = condition$xo_freq[i],
                    pos = condition$pos[i],
                    ploidy = 2)
  pop <- addFounder(pop, allow_het = FALSE)
  pop <- addFounder(pop, allow_het = FALSE)
  pop <- getProgenies(object = pop,
                      crosstype = "pair",
                      n_progeny = 1,
                      n_comb = c(1, 1),
                      randomComb = FALSE)
  pop <- getProgenies(object = pop,
                      crosstype = "self",
                      n_progeny = condition$n_progeny[i])
  for(o_read in total_read){
    for(f_read in parental_read){
      pop <- simRead(pop,
                     gen = 1,
                     total_read = f_read * condition$n_mar[i],
                     dp_dist = NULL,
                     ad_dist = NULL,
                     seq_error = 0.0025,
                     mismap = 0)
      pop <- simRead(pop,
                     gen = "last",
                     total_read = o_read * condition$n_mar[i],
                     dp_dist = NULL,
                     ad_dist = NULL,
                     seq_error = 0.0025,
                     mismap = 0)
      dataset <- paste0("ind", condition$n_progeny[i],
                        "_mar", condition$n_mar[i],
                        "_fread", f_read,
                        "_oread", o_read)
      vcf_fn <- paste0(dir, "/", dataset, ".vcf")
      writeVCF(pop, gen = c(1, "last"), vcf_fn = vcf_fn)
      save(pop, file = paste0(dir, "/", dataset, "_SimPop.Rdata"))
    }
  }
}

#########################################################
# 2 way-F2 from outbred Parents with no allele read bias#
#########################################################
dir <- "../noError_hetP2_F2"
if(!dir.exists(dir)){ dir.create(dir) }
#######################
# Simulate population
for(i in 1:nrow(condition)){
  pop <- makeSimPop(n_mar = condition$n_mar[i],
                    xo_freq = condition$xo_freq[i],
                    pos = condition$pos[i],
                    ploidy = 2)
  pop <- addFounder(pop, allow_het = TRUE)
  pop <- addFounder(pop, allow_het = TRUE)
  pop <- getProgenies(object = pop,
                      crosstype = "pair",
                      n_progeny = 1000,
                      n_comb = c(1, 1),
                      randomComb = FALSE)
  pop <- getProgenies(object = pop,
                      crosstype = "sibling",
                      n_comb = 100,
                      n_progeny = 10)
  pop <- selectSamples(object = pop, n_samples = condition$n_progeny[i])
  for(o_read in total_read){
    for(f_read in parental_read){
      pop <- simRead(pop,
                     gen = 1,
                     total_read = f_read * condition$n_mar[i],
                     dp_dist = NULL,
                     ad_dist = NULL,
                     seq_error = 0.0025,
                     mismap = 0)
      pop <- simRead(pop,
                     gen = "last",
                     total_read = o_read * condition$n_mar[i],
                     dp_dist = NULL,
                     ad_dist = NULL,
                     seq_error = 0.0025,
                     mismap = 0)
      dataset <- paste0("ind", condition$n_progeny[i],
                        "_mar", condition$n_mar[i],
                        "_fread", f_read,
                        "_oread", o_read)
      vcf_fn <- paste0(dir, "/", dataset, ".vcf")
      writeVCF(pop, gen = c(1, "last"), vcf_fn = vcf_fn)
      save(pop, file = paste0(dir, "/", dataset, "_SimPop.Rdata"))
    }
  }
}

#########################################################
# 8 way-RIL from inbred Parents with no allele read bias#
#########################################################
dir <- "../noError_homoP8_RIL"
if(!dir.exists(dir)){ dir.create(dir) }
#######################
# Simulate population
for(i in 1:nrow(condition)){
  pop <- makeSimPop(n_mar = condition$n_mar[i],
                    xo_freq = condition$xo_freq[i],
                    pos = condition$pos[i],
                    ploidy = 2)
  pop <- addFounder(pop, allow_het = FALSE)
  pop <- addFounder(pop, allow_het = FALSE)
  pop <- addFounder(pop, allow_het = FALSE)
  pop <- addFounder(pop, allow_het = FALSE)
  pop <- addFounder(pop, allow_het = FALSE)
  pop <- addFounder(pop, allow_het = FALSE)
  pop <- addFounder(pop, allow_het = FALSE)
  pop <- addFounder(pop, allow_het = FALSE)
  pop <- getProgenies(object = pop,
                      crosstype = "pair",
                      n_progeny = 1,
                      n_comb = c(4, 1),
                      randomComb = FALSE)
  pop <- getProgenies(object = pop,
                      crosstype = "pair",
                      n_progeny = 50,
                      n_comb = c(2, 1),
                      randomComb = FALSE)
  pop <- getProgenies(object = pop,
                      crosstype = "pair",
                      n_progeny = 50,
                      n_comb = c(1, 50),
                      randomComb = FALSE)
  pop <- getProgenies(object = pop,
                      crosstype = "self",
                      n_progeny = 1)
  pop <- getProgenies(object = pop,
                      crosstype = "self",
                      n_progeny = 1)
  pop <- getProgenies(object = pop,
                      crosstype = "self",
                      n_progeny = 1)
  pop <- getProgenies(object = pop,
                      crosstype = "self",
                      n_progeny = 1)
  pop <- getProgenies(object = pop,
                      crosstype = "self",
                      n_progeny = 1)
  pop <- selectSamples(object = pop, n_samples = condition$n_progeny[i])
  for(o_read in total_read){
    for(f_read in parental_read){
      pop <- simRead(pop,
                     gen = 1,
                     total_read = f_read * condition$n_mar[i],
                     dp_dist = NULL,
                     ad_dist = NULL,
                     seq_error = 0.0025,
                     mismap = 0)
      pop <- simRead(pop,
                     gen = "last",
                     total_read = o_read * condition$n_mar[i],
                     dp_dist = NULL,
                     ad_dist = NULL,
                     seq_error = 0.0025,
                     mismap = 0)
      dataset <- paste0("ind", condition$n_progeny[i],
                        "_mar", condition$n_mar[i],
                        "_fread", f_read,
                        "_oread", o_read)
      vcf_fn <- paste0(dir, "/", dataset, ".vcf")
      writeVCF(pop, gen = c(1, "last"), vcf_fn = vcf_fn)
      save(pop, file = paste0(dir, "/", dataset, "_SimPop.Rdata"))
    }
  }
}


########################################################
# 2 way-F2 from inbred Parents with allele read biases #
########################################################
dir <- "../realError_homoP2_F2"
if(!dir.exists(dir)){ dir.create(dir) }
#######################
# Simulate population
library(GBScleanR)
gbsrVCF2GDS(vcf_fn = "../realdata/gbs_nbolf2.vcf.gz", out_fn = "../realdata/gbs_nbolf2.gds", force = TRUE)
gds <- loadGDS("../realdata/gbs_nbolf2.gds")
parents <- grep("NB|OL|F1", getScanID(gds), value = TRUE)
gds <- setScanFilter(gds, id = parents)
gds <- countGenotype(gds)
gds <- thinMarker(gds, 75)
gds <- countRead(gds)
dp <- getCountRead(gds, "snp")
chr7 <- which(getChromosome(gds) == 7)
dp <- getCountRead(gds)[chr7]
dp_dist <- dp / sum(dp)
ad_dist <- getCountReadRef(gds, prop = TRUE)[chr7]
ad_dist[is.na(ad_dist)] <- 0.5
closeGDS(gds)

n_mar <- length(ad_dist)
condition <- expand.grid(xo_freq, pos, n_progeny, n_mar)
names(condition) <- c("xo_freq", "pos", "n_progeny", "n_mar")
for(i in 1:nrow(condition)){
  pop <- makeSimPop(n_mar = condition$n_mar[i],
                    xo_freq = condition$xo_freq[i],
                    pos = condition$pos[i],
                    ploidy = 2)
  pop <- addFounder(pop, allow_het = FALSE)
  pop <- addFounder(pop, allow_het = FALSE)
  pop <- getProgenies(object = pop,
                      crosstype = "pair",
                      n_progeny = 1,
                      n_comb = c(1, 1),
                      randomComb = FALSE)
  pop <- getProgenies(object = pop,
                      crosstype = "self",
                      n_progeny = condition$n_progeny[i])
  for(o_read in total_read){
    for(f_read in parental_read){
      pop <- simRead(pop,
                     gen = 1,
                     total_read = f_read * condition$n_mar[i],
                     dp_dist = dp_dist,
                     ad_dist = ad_dist,
                     seq_error = 0.0025,
                     mismap = 0.005)
      pop <- simRead(pop,
                     gen = "last",
                     total_read = o_read * condition$n_mar[i],
                     dp_dist = dp_dist,
                     ad_dist = ad_dist,
                     seq_error = 0.0025,
                     mismap = 0.005)
      dataset <- paste0("ind", condition$n_progeny[i],
                        "_mar", condition$n_mar[i],
                        "_fread", f_read,
                        "_oread", o_read)
      vcf_fn <- paste0(dir, "/", dataset, ".vcf")
      writeVCF(pop, gen = c(1, "last"), vcf_fn = vcf_fn)
      save(pop, file = paste0(dir, "/", dataset, "_SimPop.Rdata"))
    }
  }
}

#########################################################
# 2 way-F2 from outbred Parents with allele read biases #
#########################################################
dir <- "realError_hetP2_F2"
if(!dir.exists(dir)){ dir.create(dir) }
#######################
# Simulate population
for(i in 1:nrow(condition)){
  pop <- makeSimPop(n_mar = condition$n_mar[i],
                    xo_freq = condition$xo_freq[i],
                    pos = condition$pos[i],
                    ploidy = 2)
  pop <- addFounder(pop, allow_het = TRUE)
  pop <- addFounder(pop, allow_het = TRUE)
  pop <- getProgenies(object = pop,
                      crosstype = "pair",
                      n_progeny = 1000,
                      n_comb = c(1, 1),
                      randomComb = FALSE)
  pop <- getProgenies(object = pop,
                      crosstype = "sibling",
                      n_comb = 100,
                      n_progeny = 10)
  pop <- selectSamples(object = pop, n_samples = condition$n_progeny[i])
  for(o_read in total_read){
    for(f_read in parental_read){
      pop <- simRead(pop,
                     gen = 1,
                     total_read = f_read * condition$n_mar[i],
                     dp_dist = dp_dist,
                     ad_dist = ad_dist,
                     seq_error = 0.0025,
                     mismap = 0.005)
      pop <- simRead(pop,
                     gen = "last",
                     total_read = o_read * condition$n_mar[i],
                     dp_dist = dp_dist,
                     ad_dist = ad_dist,
                     seq_error = 0.0025,
                     mismap = 0.005)
      dataset <- paste0("ind", condition$n_progeny[i],
                        "_mar", condition$n_mar[i],
                        "_fread", f_read,
                        "_oread", o_read)
      vcf_fn <- paste0(dir, "/", dataset, ".vcf")
      writeVCF(pop, gen = c(1, "last"), vcf_fn = vcf_fn)
      save(pop, file = paste0(dir, "/", dataset, "_SimPop.Rdata"))
    }
  }
}


#########################################################
# 8 way-RIL from inbred Parents with allele read biases #
#########################################################
dir <- "realError_homoP8_RIL"
if(!dir.exists(dir)){ dir.create(dir) }
#######################
# Simulate population
for(i in 1:nrow(condition)){
  pop <- makeSimPop(n_mar = condition$n_mar[i],
                    xo_freq = condition$xo_freq[i],
                    pos = condition$pos[i],
                    ploidy = 2)
  pop <- addFounder(pop, allow_het = FALSE)
  pop <- addFounder(pop, allow_het = FALSE)
  pop <- addFounder(pop, allow_het = FALSE)
  pop <- addFounder(pop, allow_het = FALSE)
  pop <- addFounder(pop, allow_het = FALSE)
  pop <- addFounder(pop, allow_het = FALSE)
  pop <- addFounder(pop, allow_het = FALSE)
  pop <- addFounder(pop, allow_het = FALSE)
  pop <- getProgenies(object = pop,
                      crosstype = "pair",
                      n_progeny = 1,
                      n_comb = c(4, 1),
                      randomComb = FALSE)
  pop <- getProgenies(object = pop,
                      crosstype = "pair",
                      n_progeny = 100,
                      n_comb = c(2, 1),
                      randomComb = FALSE)
  pop <- getProgenies(object = pop,
                      crosstype = "pair",
                      n_progeny = 50,
                      n_comb = c(1, 50),
                      randomComb = FALSE)
  pop <- getProgenies(object = pop,
                      crosstype = "self",
                      n_progeny = 1)
  pop <- getProgenies(object = pop,
                      crosstype = "self",
                      n_progeny = 1)
  pop <- getProgenies(object = pop,
                      crosstype = "self",
                      n_progeny = 1)
  pop <- getProgenies(object = pop,
                      crosstype = "self",
                      n_progeny = 1)
  pop <- getProgenies(object = pop,
                      crosstype = "self",
                      n_progeny = 1)
  pop <- selectSamples(object = pop, n_samples = condition$n_progeny[i])
  for(o_read in total_read){
    for(f_read in parental_read){
      pop <- simRead(pop,
                     gen = 1,
                     total_read = f_read * condition$n_mar[i],
                     dp_dist = dp_dist,
                     ad_dist = ad_dist,
                     seq_error = 0.0025,
                     mismap = 0.005)
      pop <- simRead(pop,
                     gen = "last",
                     total_read = o_read * condition$n_mar[i],
                     dp_dist = dp_dist,
                     ad_dist = ad_dist,
                     seq_error = 0.0025,
                     mismap = 0.005)
      dataset <- paste0("ind", condition$n_progeny[i],
                        "_mar", condition$n_mar[i],
                        "_fread", f_read,
                        "_oread", o_read)
      vcf_fn <- paste0(dir, "/", dataset, ".vcf")
      writeVCF(pop, gen = c(1, "last"), vcf_fn = vcf_fn)
      save(pop, file = paste0(dir, "/", dataset, "_SimPop.Rdata"))
    }
  }
}
