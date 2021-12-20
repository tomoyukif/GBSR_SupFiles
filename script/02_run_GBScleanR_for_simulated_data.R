library(GBScleanR)
library(SimPop)
source("./compareGeno.R")


########################################################
# 2 way-F2 from inbred Parents with no allele read bias#
########################################################
dir <- "../noError_homoP2_F2"

# With the interative parameter optimization
files <- list.files(dir, ".vcf")
for(file_i in files){
  vcf_fn <- paste0(dir, "/", file_i)
  gds_fn <- sub(".vcf", "_gbsr.gds", vcf_fn)
  gbsrVCF2GDS(vcf_fn, gds_fn, force = T)
  gds <- loadGDS(gds_fn)
  parents <- grep("Founder", getScanID(gds), value = TRUE)
  gds <- setParents(gds, parents, flip = FALSE, mono = FALSE, bi = FALSE)
  gds <- initScheme(gds, "pair", matrix(1:2, 2))
  gds <- addScheme(gds, "self")
  gds <- estGeno(object = gds,
               recomb_rate = 0.04,
               call_threshold = 0.9,
               error_rate = 0.0025,
               iter = 4,
               het_parent = FALSE,
               optim = TRUE)
  closeGDS(gds)
}

######################
# Evaluate correction
files <- list.files(dir, pattern = "_gbsr.gds")
evaldf <- NULL
for(file_i in files){
  gds_fn <- paste0(dir, "/", file_i)
  gds <- loadGDS(gds_fn)
  parents <- grep("Founder", getScanID(gds), value = TRUE)
  gds <- setParents(gds, parents, flip = FALSE, mono = FALSE, bi = FALSE)
  rawgeno <- abs(GBScleanR::getGenotype(gds, node = "raw") - 2)
  geno <- abs(GBScleanR::getGenotype(gds, node = "cor") - 2)
  pop_file <- sub("_gbsr.gds", "_SimPop.Rdata", gds_fn)
  load(pop_file)
  true <- SimPop::getGenotype(pop, gen = "last", fam = "all", sib = "all")
  true <- t(sapply(true, colSums))
  df <- compareGeno(geno, true)
  df_ind_mean <- apply(df$ind, 2, mean, na.rm = TRUE)
  df_ind_sd <- apply(df$ind, 2, sd, na.rm = TRUE)
  names(df_ind_sd) <- paste(names(df_ind_sd), "SD", sep = "_")
  evaldf <- rbind(evaldf,
                  cbind(data.frame(condition = sub(".vcf", "", file_i)),
                        data.frame(t(df_ind_mean)),
                        data.frame(t(df_ind_sd))))
  closeGDS(gds)
}
write.csv(evaldf, file = paste0(dir, "/Eval.csv"))


# Without the interative parameter optimization
files <- list.files(dir, ".vcf")
for(file_i in files){
  vcf_fn <- paste0(dir, "/", file_i)
  gds_fn <- sub(".vcf", "_gbsr_noOptim.gds", vcf_fn)
  gbsrVCF2GDS(vcf_fn, gds_fn, force = T)
  gds <- loadGDS(gds_fn)
  parents <- grep("Founder", getScanID(gds), value = TRUE)
  gds <- setParents(gds, parents, flip = FALSE, mono = FALSE, bi = FALSE)
  gds <- initScheme(gds, "pair", matrix(1:2, 2))
  gds <- addScheme(gds, "self")
  gds <- estGeno(object = gds,
               recomb_rate = 0.04,
               call_threshold = 0.9,
               error_rate = 0.0025,
               iter = 4,
               het_parent = FALSE,
               optim = FALSE)
  closeGDS(gds)
}

######################
# Evaluate correction
files <- list.files(dir, pattern = "_gbsr_noOptim.gds")
evaldf <- NULL
for(file_i in files){
  gds_fn <- paste0(dir, "/", file_i)
  gds <- loadGDS(gds_fn)
  parents <- grep("Founder", getScanID(gds), value = TRUE)
  gds <- setParents(gds, parents, flip = FALSE, mono = FALSE, bi = FALSE)
  rawgeno <- abs(GBScleanR::getGenotype(gds, node = "raw") - 2)
  geno <- abs(GBScleanR::getGenotype(gds, node = "cor") - 2)
  pop_file <- sub("_gbsr_noOptim.gds", "_SimPop.Rdata", gds_fn)
  load(pop_file)
  true <- SimPop::getGenotype(pop, gen = "last", fam = "all", sib = "all")
  true <- t(sapply(true, colSums))
  df <- compareGeno(geno, true)
  df_ind_mean <- apply(df$ind, 2, mean, na.rm = TRUE)
  df_ind_sd <- apply(df$ind, 2, sd, na.rm = TRUE)
  names(df_ind_sd) <- paste(names(df_ind_sd), "SD", sep = "_")
  evaldf <- rbind(evaldf,
                  cbind(data.frame(condition = sub(".vcf", "", file_i)),
                        data.frame(t(df_ind_mean)),
                        data.frame(t(df_ind_sd))))
  closeGDS(gds)
}
write.csv(evaldf, file = paste0(dir, "/Eval_noOptim.csv"))


########################################################
# 2 way-F2 from inbred Parents with allele read biases #
########################################################
dir <- "../realError_homoP2_F2"

# With the iterative parameter optimization
files <- list.files(dir, ".vcf")
for(file_i in files){
  vcf_fn <- paste0(dir, "/", file_i)
  gds_fn <- sub(".vcf", "_gbsr.gds", vcf_fn)
  gbsrVCF2GDS(vcf_fn, gds_fn, force = T)
  gds <- loadGDS(gds_fn)
  parents <- grep("Founder", getScanID(gds), value = TRUE)
  gds <- setParents(gds, parents, flip = FALSE, mono = FALSE, bi = FALSE)
  gds <- initScheme(gds, "pair", matrix(1:2, 2))
  gds <- addScheme(gds, "self")
  gds <- estGeno(object = gds,
               recomb_rate = 0.04,
               call_threshold = 0.9,
               error_rate = 0.0025,
               iter = 4,
               het_parent = FALSE,
               optim = TRUE)
  closeGDS(gds)
}

######################
# Evaluate correction
evaldf <- NULL
files <- list.files(dir, pattern = "_gbsr.gds")
for(file_i in files){
  gds_fn <- paste0(dir, "/", file_i)
  gds <- loadGDS(gds_fn)
  parents <- grep("Founder", getScanID(gds), value = TRUE)
  gds <- setParents(gds, parents, flip = FALSE, mono = FALSE, bi = FALSE)
  rawgeno <- abs(GBScleanR::getGenotype(gds, node = "raw") - 2)
  geno <- abs(GBScleanR::getGenotype(gds, node = "cor") - 2)
  pop_file <- sub("_gbsr.gds", "_SimPop.Rdata", gds_fn)
  load(pop_file)
  true <- SimPop::getGenotype(pop, gen = "last", fam = "all", sib = "all")
  true <- t(sapply(true, colSums))
  df <- compareGeno(geno, true)
  df_ind_mean <- apply(df$ind, 2, mean, na.rm = TRUE)
  df_ind_sd <- apply(df$ind, 2, sd, na.rm = TRUE)
  names(df_ind_sd) <- paste(names(df_ind_sd), "SD", sep = "_")
  evaldf <- rbind(evaldf,
                  cbind(data.frame(condition = sub(".vcf", "", file_i)),
                        data.frame(t(df_ind_mean)),
                        data.frame(t(df_ind_sd))))
  closeGDS(gds)
}
write.csv(evaldf, file = paste0(dir, "/Eval.csv"))

# Without the iterative parameter optimization
files <- list.files(dir, ".vcf")
for(file_i in files){
  vcf_fn <- paste0(dir, "/", file_i)
  gds_fn <- sub(".vcf", "_gbsr_noOptim.gds", vcf_fn)
  gbsrVCF2GDS(vcf_fn, gds_fn, force = T)
  gds <- loadGDS(gds_fn)
  parents <- grep("Founder", getScanID(gds), value = TRUE)
  gds <- setParents(gds, parents, flip = FALSE, mono = FALSE, bi = FALSE)
  gds <- initScheme(gds, "pair", matrix(1:2, 2))
  gds <- addScheme(gds, "self")
  gds <- estGeno(object = gds,
               recomb_rate = 0.04,
               call_threshold = 0.9,
               error_rate = 0.0025,
               iter = 4,
               het_parent = FALSE,
               optim = FALSE)
  closeGDS(gds)
}

######################
# Evaluate correction
evaldf <- NULL
files <- list.files(dir, pattern = "_gbsr_noOptim.gds")
for(file_i in files){
  gds_fn <- paste0(dir, "/", file_i)
  gds <- loadGDS(gds_fn)
  parents <- grep("Founder", getScanID(gds), value = TRUE)
  gds <- setParents(gds, parents, flip = FALSE, mono = FALSE, bi = FALSE)
  rawgeno <- abs(GBScleanR::getGenotype(gds, node = "raw") - 2)
  geno <- abs(GBScleanR::getGenotype(gds, node = "cor") - 2)
  pop_file <- sub("_gbsr_noOptim.gds", "_SimPop.Rdata", gds_fn)
  load(pop_file)
  true <- SimPop::getGenotype(pop, gen = "last", fam = "all", sib = "all")
  true <- t(sapply(true, colSums))
  df <- compareGeno(geno, true)
  df_ind_mean <- apply(df$ind, 2, mean, na.rm = TRUE)
  df_ind_sd <- apply(df$ind, 2, sd, na.rm = TRUE)
  names(df_ind_sd) <- paste(names(df_ind_sd), "SD", sep = "_")
  evaldf <- rbind(evaldf,
                  cbind(data.frame(condition = sub(".vcf", "", file_i)),
                        data.frame(t(df_ind_mean)),
                        data.frame(t(df_ind_sd))))
  closeGDS(gds)
}
write.csv(evaldf, file = paste0(dir, "/Eval_noOptim.csv"))


#########################################################
# 2 way-F2 from outbred Parents with no allele read bias#
#########################################################
dir <- "../noError_hetP2_F2"

# With the iterative parameter optimization
files <- list.files(dir, ".vcf")
for(file_i in files){
  vcf_fn <- paste0(dir, "/", file_i)
  gds_fn <- sub(".vcf", "_gbsr.gds", vcf_fn)
  gbsrVCF2GDS(vcf_fn, gds_fn, force = T)
  gds <- loadGDS(gds_fn)
  parents <- grep("Founder", getScanID(gds), value = TRUE)
  gds <- setParents(gds, parents, flip = FALSE, mono = FALSE, bi = FALSE)
  gds <- initScheme(gds, "pairing", matrix(1:2, 2))
  gds <- addScheme(gds, "sibling")
  gds <- estGeno(object = gds,
               recomb_rate = 0.04,
               call_threshold = 0.9,
               error_rate = 0.0025,
               iter = 4,
               het_parent = TRUE,
               optim = TRUE)
  closeGDS(gds)
}

######################
# Evaluate correction
evaldf <- NULL
files <- list.files(dir, pattern = "_gbsr.gds")
for(file_i in files){
  gds_fn <- paste0(dir, "/", file_i)
  gds <- loadGDS(gds_fn)
  parents <- grep("Founder", getScanID(gds), value = TRUE)
  gds <- setParents(gds, parents, flip = FALSE, mono = FALSE, bi = FALSE)
  rawgeno <- abs(GBScleanR::getGenotype(gds, node = "raw") - 2)
  geno <- abs(GBScleanR::getGenotype(gds, node = "cor") - 2)
  pop_file <- sub("_gbsr.gds", "_SimPop.Rdata", gds_fn)
  load(pop_file)
  true <- SimPop::getGenotype(pop, gen = "last", fam = "all", sib = "all")
  true <- t(sapply(true, colSums))
  df <- compareGeno(geno, true)
  df_ind_mean <- apply(df$ind, 2, mean, na.rm = TRUE)
  df_ind_sd <- apply(df$ind, 2, sd, na.rm = TRUE)
  names(df_ind_sd) <- paste(names(df_ind_sd), "SD", sep = "_")
  evaldf <- rbind(evaldf,
                  cbind(data.frame(condition = sub(".vcf", "", file_i)),
                        data.frame(t(df_ind_mean)),
                        data.frame(t(df_ind_sd))))
  closeGDS(gds)
}
write.csv(evaldf, file = paste0(dir, "/Eval.csv"))

# Without the iterative parameter optimization
files <- list.files(dir, ".vcf")
for(file_i in files){
  vcf_fn <- paste0(dir, "/", file_i)
  gds_fn <- sub(".vcf", "_gbsr_noOptim.gds", vcf_fn)
  gbsrVCF2GDS(vcf_fn, gds_fn, force = T)
  gds <- loadGDS(gds_fn)
  parents <- grep("Founder", getScanID(gds), value = TRUE)
  gds <- setParents(gds, parents, flip = FALSE, mono = FALSE, bi = FALSE)
  gds <- initScheme(gds, "pairing", matrix(1:2, 2))
  gds <- addScheme(gds, "sibling")
  gds <- estGeno(object = gds,
               recomb_rate = 0.04,
               call_threshold = 0.9,
               error_rate = 0.0025,
               iter = 4,
               het_parent = TRUE,
               optim = FALSE)
  closeGDS(gds)
}

######################
# Evaluate correction
evaldf <- NULL
files <- list.files(dir, pattern = "_gbsr_noOptim.gds")
for(file_i in files){
  gds_fn <- paste0(dir, "/", file_i)
  gds <- loadGDS(gds_fn)
  parents <- grep("Founder", getScanID(gds), value = TRUE)
  gds <- setParents(gds, parents, flip = FALSE, mono = FALSE, bi = FALSE)
  rawgeno <- abs(GBScleanR::getGenotype(gds, node = "raw") - 2)
  geno <- abs(GBScleanR::getGenotype(gds, node = "cor") - 2)
  pop_file <- sub("_gbsr_noOptim.gds", "_SimPop.Rdata", gds_fn)
  load(pop_file)
  true <- SimPop::getGenotype(pop, gen = "last", fam = "all", sib = "all")
  true <- t(sapply(true, colSums))
  df <- compareGeno(geno, true)
  df_ind_mean <- apply(df$ind, 2, mean, na.rm = TRUE)
  df_ind_sd <- apply(df$ind, 2, sd, na.rm = TRUE)
  names(df_ind_sd) <- paste(names(df_ind_sd), "SD", sep = "_")
  evaldf <- rbind(evaldf,
                  cbind(data.frame(condition = sub(".vcf", "", file_i)),
                        data.frame(t(df_ind_mean)),
                        data.frame(t(df_ind_sd))))
  closeGDS(gds)
}
write.csv(evaldf, file = paste0(dir, "/Eval_noOptim.csv"))


#########################################################
# 2 way-F2 from outbred Parents with allele read biases #
#########################################################
dir <- "realError_hetP2_F2"
# With the iterative paramter optimization
files <- list.files(dir, ".vcf")
for(file_i in files){
  vcf_fn <- paste0(dir, "/", file_i)
  gds_fn <- sub(".vcf", "_gbsr.gds", vcf_fn)
  gbsrVCF2GDS(vcf_fn, gds_fn, force = T)
  gds <- loadGDS(gds_fn)
  parents <- grep("Founder", getScanID(gds), value = TRUE)
  gds <- setParents(gds, parents, flip = FALSE, mono = FALSE, bi = FALSE)
  gds <- initScheme(gds, "pair", matrix(1:2, 2))
  gds <- addScheme(gds, "sibling")
  gds <- estGeno(object = gds,
               recomb_rate = 0.04,
               call_threshold = 0.9,
               error_rate = 0.0025,
               iter = 4,
               het_parent = TRUE,
               optim = TRUE)
  closeGDS(gds)
}

######################
# Evaluate correction
evaldf <- NULL
files <- list.files(dir, pattern = "_gbsr.gds")
for(file_i in files){
  gds_fn <- paste0(dir, "/", file_i)
  gds <- loadGDS(gds_fn)
  parents <- grep("Founder", getScanID(gds), value = TRUE)
  gds <- setParents(gds, parents, flip = FALSE, mono = FALSE, bi = FALSE)
  rawgeno <- abs(GBScleanR::getGenotype(gds, node = "raw") - 2)
  geno <- abs(GBScleanR::getGenotype(gds, node = "cor") - 2)
  pop_file <- sub("_gbsr.gds", "_SimPop.Rdata", gds_fn)
  load(pop_file)
  true <- SimPop::getGenotype(pop, gen = "last", fam = "all", sib = "all")
  true <- t(sapply(true, colSums))
  df <- compareGeno(geno, true)
  df_ind_mean <- apply(df$ind, 2, mean, na.rm = TRUE)
  df_ind_sd <- apply(df$ind, 2, sd, na.rm = TRUE)
  names(df_ind_sd) <- paste(names(df_ind_sd), "SD", sep = "_")
  evaldf <- rbind(evaldf,
                  cbind(data.frame(condition = sub(".vcf", "", file_i)),
                        data.frame(t(df_ind_mean)),
                        data.frame(t(df_ind_sd))))
  closeGDS(gds)
}
write.csv(evaldf, file = paste0(dir, "/Eval.csv"))

# Without the iterative parameter optimization
library(GBScleanR)
files <- list.files(dir, ".vcf")
for(file_i in files){
  vcf_fn <- paste0(dir, "/", file_i)
  gds_fn <- sub(".vcf", "_gbsr_noOptim.gds", vcf_fn)
  gbsrVCF2GDS(vcf_fn, gds_fn, force = T)
  gds <- loadGDS(gds_fn)
  parents <- grep("Founder", getScanID(gds), value = TRUE)
  gds <- setParents(gds, parents, flip = FALSE, mono = FALSE, bi = FALSE)
  gds <- initScheme(gds, "pair", matrix(1:2, 2))
  gds <- addScheme(gds, "sibling")
  gds <- estGeno(object = gds,
               recomb_rate = 0.04,
               call_threshold = 0.9,
               error_rate = 0.0025,
               iter = 4,
               het_parent = TRUE,
               optim = FALSE)
  closeGDS(gds)
}

######################
# Evaluate correction
evaldf <- NULL
files <- list.files(dir, pattern = "_gbsr_noOptim.gds")
for(file_i in files){
  gds_fn <- paste0(dir, "/", file_i)
  gds <- loadGDS(gds_fn)
  parents <- grep("Founder", getScanID(gds), value = TRUE)
  gds <- setParents(gds, parents, flip = FALSE, mono = FALSE, bi = FALSE)
  rawgeno <- abs(GBScleanR::getGenotype(gds, node = "raw") - 2)
  geno <- abs(GBScleanR::getGenotype(gds, node = "cor") - 2)
  pop_file <- sub("_gbsr_noOptim.gds", "_SimPop.Rdata", gds_fn)
  load(pop_file)
  true <- SimPop::getGenotype(pop, gen = "last", fam = "all", sib = "all")
  true <- t(sapply(true, colSums))
  df <- compareGeno(geno, true)
  df_ind_mean <- apply(df$ind, 2, mean, na.rm = TRUE)
  df_ind_sd <- apply(df$ind, 2, sd, na.rm = TRUE)
  names(df_ind_sd) <- paste(names(df_ind_sd), "SD", sep = "_")
  evaldf <- rbind(evaldf,
                  cbind(data.frame(condition = sub(".vcf", "", file_i)),
                        data.frame(t(df_ind_mean)),
                        data.frame(t(df_ind_sd))))
  closeGDS(gds)
}
write.csv(evaldf, file = paste0(dir, "/Eval_noOptim.csv"))


#########################################################
# 8 way-RIL from inbred Parents with no allele read bias#
#########################################################
dir <- "../noError_homoP8_RIL"
# With the iterative parameter optimization
files <- list.files(dir, ".vcf")
for(file_i in files){
  vcf_fn <- paste0(dir, "/", file_i)
  gds_fn <- sub(".vcf", "_gbsr.gds", vcf_fn)
  gbsrVCF2GDS(vcf_fn, gds_fn, force = T)
  gds <- loadGDS(gds_fn)
  parents <- grep("Founder", getScanID(gds), value = TRUE)
  gds <- setParents(gds, parents, flip = FALSE, mono = FALSE, bi = FALSE)

  gds <- initScheme(gds, "pairing", matrix(1:8, 2))
  gds <- addScheme(gds, "pairing", matrix(9:12, 2))
  gds <- addScheme(gds, "pairing", matrix(13:14, 2))
  gds <- addScheme(gds, "self")
  gds <- addScheme(gds, "self")
  gds <- addScheme(gds, "self")
  gds <- addScheme(gds, "self")
  gds <- addScheme(gds, "self")
  gds <- estGeno(object = gds,
               recomb_rate = 0.04,
               call_threshold = 0.9,
               error_rate = 0.0025,
               iter = 4,
               het_parent = FALSE,
               optim = TRUE)
  closeGDS(gds)
}

######################
# Evaluate correction
evaldf <- NULL
files <- list.files(dir, pattern = "_gbsr.gds")
for(file_i in files){
  gds_fn <- paste0(dir, "/", file_i)
  gds <- loadGDS(gds_fn)
  parents <- grep("Founder", getScanID(gds), value = TRUE)
  gds <- setParents(gds, parents, flip = FALSE, mono = FALSE, bi = FALSE)
  rawgeno <- abs(GBScleanR::getGenotype(gds, node = "raw") - 2)
  geno <- abs(GBScleanR::getGenotype(gds, node = "cor") - 2)
  pop_file <- sub("_gbsr.gds", "_SimPop.Rdata", gds_fn)
  load(pop_file)
  true <- SimPop::getGenotype(pop, gen = "last", fam = "all", sib = "all")
  true <- t(sapply(true, colSums))
  df <- compareGeno(geno, true)
  df_ind_mean <- apply(df$ind, 2, mean, na.rm = TRUE)
  df_ind_sd <- apply(df$ind, 2, sd, na.rm = TRUE)
  names(df_ind_sd) <- paste(names(df_ind_sd), "SD", sep = "_")
  evaldf <- rbind(evaldf,
                  cbind(data.frame(condition = sub(".vcf", "", file_i)),
                        data.frame(t(df_ind_mean)),
                        data.frame(t(df_ind_sd))))
  closeGDS(gds)
}
write.csv(evaldf, file = paste0(dir, "/Eval.csv"))

# Without the iterative parameter optimization
files <- list.files(dir, ".vcf")
for(file_i in files){
  vcf_fn <- paste0(dir, "/", file_i)
  gds_fn <- sub(".vcf", "_gbsr_noOptim.gds", vcf_fn)
  gbsrVCF2GDS(vcf_fn, gds_fn, force = T)
  gds <- loadGDS(gds_fn)
  gds <- setParents(gds, parents, flip = FALSE, mono = FALSE, bi = FALSE)

  gds <- initScheme(gds, "pairing", matrix(1:8, 2))
  gds <- addScheme(gds, "pairing", matrix(9:12, 2))
  gds <- addScheme(gds, "pairing", matrix(13:14, 2))
  gds <- addScheme(gds, "self")
  gds <- addScheme(gds, "self")
  gds <- addScheme(gds, "self")
  gds <- addScheme(gds, "self")
  gds <- addScheme(gds, "self")
  gds <- estGeno(object = gds,
               recomb_rate = 0.04,
               call_threshold = 0.9,
               error_rate = 0.0025,
               iter = 4,
               het_parent = FALSE,
               optim = FALSE)
  closeGDS(gds)
}

######################
# Evaluate correction
evaldf <- NULL
files <- list.files(dir, pattern = "_gbsr_noOptim.gds")
for(file_i in files){
  gds_fn <- paste0(dir, "/", file_i)
  gds <- loadGDS(gds_fn)
  parents <- grep("Founder", getScanID(gds), value = TRUE)
  gds <- setParents(gds, parents, flip = FALSE, mono = FALSE, bi = FALSE)
  rawgeno <- abs(GBScleanR::getGenotype(gds, node = "raw") - 2)
  geno <- abs(GBScleanR::getGenotype(gds, node = "cor") - 2)
  pop_file <- sub("_gbsr_noOptim.gds", "_SimPop.Rdata", gds_fn)
  load(pop_file)
  true <- SimPop::getGenotype(pop, gen = "last", fam = "all", sib = "all")
  true <- t(sapply(true, colSums))
  df <- compareGeno(geno, true)
  df_ind_mean <- apply(df$ind, 2, mean, na.rm = TRUE)
  df_ind_sd <- apply(df$ind, 2, sd, na.rm = TRUE)
  names(df_ind_sd) <- paste(names(df_ind_sd), "SD", sep = "_")
  evaldf <- rbind(evaldf,
                  cbind(data.frame(condition = sub(".vcf", "", file_i)),
                        data.frame(t(df_ind_mean)),
                        data.frame(t(df_ind_sd))))
  closeGDS(gds)
}
write.csv(evaldf, file = paste0(dir, "/Eval_noOptim.csv"))


#########################################################
# 8 way-RIL from inbred Parents with allele read biases #
#########################################################
dir <- "realError_homoP8_RIL"
# With the iterative parameter optimization
files <- list.files(dir, ".vcf")
for(file_i in files){
  vcf_fn <- paste0(dir, "/", file_i)
  gds_fn <- sub(".vcf", "_gbsr.gds", vcf_fn)
  gbsrVCF2GDS(vcf_fn, gds_fn, force = T)
  gds <- loadGDS(gds_fn)
  parents <- grep("Founder", getScanID(gds), value = TRUE)
  gds <- setParents(gds, parents, flip = FALSE, mono = FALSE, bi = FALSE)

  gds <- initScheme(gds, "pairing", matrix(1:8, 2))
  gds <- addScheme(gds, "pairing", matrix(9:12, 2))
  gds <- addScheme(gds, "pairing", matrix(13:14, 2))
  gds <- addScheme(gds, "self")
  gds <- addScheme(gds, "self")
  gds <- addScheme(gds, "self")
  gds <- addScheme(gds, "self")
  gds <- addScheme(gds, "self")
  gds <- estGeno(object = gds,
               recomb_rate = 0.04,
               call_threshold = 0.9,
               error_rate = 0.0025,
               iter = 4,
               het_parent = FALSE,
               optim = TRUE)
  closeGDS(gds)
}

######################
# Evaluate correction
evaldf <- NULL
files <- list.files(dir, pattern = "_gbsr.gds")
for(file_i in files){
  gds_fn <- paste0(dir, "/", file_i)
  gds <- loadGDS(gds_fn)
  parents <- grep("Founder", getScanID(gds), value = TRUE)
  gds <- setParents(gds, parents, flip = FALSE, mono = FALSE, bi = FALSE)
  rawgeno <- abs(GBScleanR::getGenotype(gds, node = "raw") - 2)
  geno <- abs(GBScleanR::getGenotype(gds, node = "cor") - 2)
  pop_file <- sub("_gbsr.gds", "_SimPop.Rdata", gds_fn)
  load(pop_file)
  true <- SimPop::getGenotype(pop, gen = "last", fam = "all", sib = "all")
  true <- t(sapply(true, colSums))
  df <- compareGeno(geno, true)
  df_ind_mean <- apply(df$ind, 2, mean, na.rm = TRUE)
  df_ind_sd <- apply(df$ind, 2, sd, na.rm = TRUE)
  names(df_ind_sd) <- paste(names(df_ind_sd), "SD", sep = "_")
  evaldf <- rbind(evaldf,
                  cbind(data.frame(condition = sub(".vcf", "", file_i)),
                        data.frame(t(df_ind_mean)),
                        data.frame(t(df_ind_sd))))
  closeGDS(gds)
}
write.csv(evaldf, file = paste0(dir, "/Eval.csv"))

# Without the iterative parameter optimization
files <- list.files(dir, ".vcf")
for(file_i in files){
  vcf_fn <- paste0(dir, "/", file_i)
  gds_fn <- sub(".vcf", "_gbsr_noOptim.gds", vcf_fn)
  gbsrVCF2GDS(vcf_fn, gds_fn, force = T)
  gds <- loadGDS(gds_fn)
  parents <- grep("Founder", getScanID(gds), value = TRUE)
  gds <- setParents(gds, parents, flip = FALSE, mono = FALSE, bi = FALSE)

  gds <- initScheme(gds, "pairing", matrix(1:8, 2))
  gds <- addScheme(gds, "pairing", matrix(9:12, 2))
  gds <- addScheme(gds, "pairing", matrix(13:14, 2))
  gds <- addScheme(gds, "self")
  gds <- addScheme(gds, "self")
  gds <- addScheme(gds, "self")
  gds <- addScheme(gds, "self")
  gds <- addScheme(gds, "self")
  gds <- estGeno(object = gds,
               recomb_rate = 0.04,
               call_threshold = 0.9,
               error_rate = 0.0025,
               iter = 4,
               het_parent = FALSE,
               optim = FALSE)
  closeGDS(gds)
}

######################
# Evaluate correction
evaldf <- NULL
files <- list.files(dir, pattern = "_gbsr_noOptim.gds")
for(file_i in files){
  gds_fn <- paste0(dir, "/", file_i)
  gds <- loadGDS(gds_fn)
  parents <- grep("Founder", getScanID(gds), value = TRUE)
  gds <- setParents(gds, parents, flip = FALSE, mono = FALSE, bi = FALSE)
  rawgeno <- abs(GBScleanR::getGenotype(gds, node = "raw") - 2)
  geno <- abs(GBScleanR::getGenotype(gds, node = "cor") - 2)
  pop_file <- sub("_gbsr_noOptim.gds", "_SimPop.Rdata", gds_fn)
  load(pop_file)
  true <- SimPop::getGenotype(pop, gen = "last", fam = "all", sib = "all")
  true <- t(sapply(true, colSums))
  df <- compareGeno(geno, true)
  df_ind_mean <- apply(df$ind, 2, mean, na.rm = TRUE)
  df_ind_sd <- apply(df$ind, 2, sd, na.rm = TRUE)
  names(df_ind_sd) <- paste(names(df_ind_sd), "SD", sep = "_")
  evaldf <- rbind(evaldf,
                  cbind(data.frame(condition = sub(".vcf", "", file_i)),
                        data.frame(t(df_ind_mean)),
                        data.frame(t(df_ind_sd))))
  closeGDS(gds)
}
write.csv(evaldf, file = paste0(dir, "/Eval_noOptim.csv"))
