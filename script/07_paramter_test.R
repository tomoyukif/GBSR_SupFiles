dir <- "../ParamTest"

if(!dir.exists(dir)){ dir.create(dir) }
library(GBScleanR)
library(parallel)
source("./compareGeno.R")
src_dir <- "../realError_homoP2_F2"
files <- list.files(src_dir, ".vcf")
files <- grep("_LB.vcf", files, invert = T, value = T)
vcf_fn <- paste0(src_dir, "/ind1000_mar620_fread3_oread3.vcf")
new_vf_fn <- paste0(dir, "/ind1000_mar620_fread3_oread3.vcf")
file.copy(vcf_fn, new_vf_fn)
options(scipen = 10^5)

for(i in c(1:6)){
  gds_fn <- sub(".vcf", paste0("_iter_", i, ".gds"), new_vf_fn)
  gbsrVCF2GDS(new_vf_fn, gds_fn, force = T)
  gds <- loadGDS(gds_fn)
  parents <- grep("Founder", getScanID(gds), value = TRUE)
  gds <- setParents(gds, parents, flip = FALSE, mono = FALSE, bi = FALSE)
  gds <- initScheme(gds, "pair", matrix(1:2, 2))
  gds <- addScheme(gds, "self")
  gds <- estGeno(object = gds,
               recomb_rate = 0.04,
               call_threshold = 0.9,
               error_rate = 0.0025,
               iter = i,
               het_parent = FALSE,
               optim = TRUE)
  closeGDS(gds)
}

for(i in c(0.004, 0.01, 0.02, 0.04, 0.08, 0.1)){
  gds_fn <- sub(".vcf", paste0("_recomb_", i, ".gds"), new_vf_fn)
  gbsrVCF2GDS(new_vf_fn, gds_fn, force = T)
  gds <- loadGDS(gds_fn)
  parents <- grep("Founder", getScanID(gds), value = TRUE)
  gds <- setParents(gds, parents, flip = FALSE, mono = FALSE, bi = FALSE)
  gds <- initScheme(gds, "pair", matrix(1:2, 2))
  gds <- addScheme(gds, "self")
  gds <- estGeno(object = gds,
               recomb_rate = i,
               call_threshold = 0.9,
               error_rate = 0.0025,
               iter = 4,
               het_parent = FALSE,
               optim = TRUE)
  closeGDS(gds)
}


for(i in c(0.75, 0.8, 0.85, 0.9, 0.95, 1)){
  gds_fn <- sub(".vcf", paste0("_call_", i, ".gds"), new_vf_fn)
  gbsrVCF2GDS(new_vf_fn, gds_fn, force = T)
  gds <- loadGDS(gds_fn)
  parents <- grep("Founder", getScanID(gds), value = TRUE)
  gds <- setParents(gds, parents, flip = FALSE, mono = FALSE, bi = FALSE)
  gds <- initScheme(gds, "pair", matrix(1:2, 2))
  gds <- addScheme(gds, "self")
  gds <- estGeno(object = gds,
               recomb_rate = 0.04,
               call_threshold = i,
               error_rate = 0.0025,
               iter = 4,
               het_parent = FALSE,
               optim = TRUE)
  closeGDS(gds)
}

for(i in c(0.0005, 0.001, 0.0025, 0.005, 0.01, 0.05)){
  gds_fn <- sub(".vcf", paste0("_error_", i, ".gds"), new_vf_fn)
  gbsrVCF2GDS(new_vf_fn, gds_fn, force = T)
  gds <- loadGDS(gds_fn)
  parents <- grep("Founder", getScanID(gds), value = TRUE)
  gds <- setParents(gds, parents, flip = FALSE, mono = FALSE, bi = FALSE)
  gds <- initScheme(gds, "pair", matrix(1:2, 2))
  gds <- addScheme(gds, "self")
  gds <- estGeno(object = gds,
               recomb_rate = 0.04,
               call_threshold = 0.9,
               error_rate = i,
               iter = 4,
               het_parent = FALSE,
               optim = TRUE)
  closeGDS(gds)
}

pop_file <- sub(".vcf", "_SimPop.Rdata", vcf_fn)
load(pop_file)
true <- SimPop::getGenotype(pop, gen = "last", fam = "all", sib = "all")
true <- t(sapply(true, colSums))
files <- list.files(dir, ".gds", full.names = T)
evaldf <- NULL
for(i in files){
  value <- sub(".*_", "", sub("\\.gds", "", i))
  test <- sub("_.*", "", sub(".*fread3_oread3_", "", i))
  gds <- loadGDS(i)
  parents <- grep("Founder", getScanID(gds), value = TRUE)
  gds <- setParents(gds, parents, flip = FALSE, mono = FALSE, bi = FALSE)
  geno <- abs(getGenotype(gds, node = "cor") - 2)
  eval <- apply(compareGeno(geno, true)$ind, 2, mean)
  evaldf <- rbind(evaldf, data.frame(test = test, value = value, t(eval)))
  closeGDS(gds)
}

write.table(evaldf, "../ParamTest/ParamTest.csv", sep = ",", col.names = T, row.names = F)
