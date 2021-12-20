## magicImpute is a package for Mathematica and 
## cannot be executed via CUI.
## This script supposed users have already run 
## magicImpute for all simulated datasets and 
## the real data. 

## The file names of results created by magicImpute 
## are supposed to have "_ImputedGenotype.csv" as suffix
## that following each of file names for input VCF files.
## For example, if the input VCF file was 
## "../noError_homoP2_F2.vcf" and the converted magicsnp file
## was "../noError_homoP2_F2_magicSNP.csv", the result file
## should be "../noError_homoP2_F2_ImputedGenotype.csv"

 
# Evaluate correction by magicImpute
source("./compareGeno.R")
library(GBScleanR)
library(SimPop)
dirs <- grep("noError|realError|realdata", list.dirs("../", recursive = F), value = T)
for(dir_i in dirs){
  evaldf <- NULL
  files <- list.files(dir_i, "_ImputedGenotype.csv", full.names = T)
  for(file_i in files){
    geno <- read.table(file_i, header = FALSE, sep = ",", skip = 4, stringsAsFactors = FALSE)
    col1 <- geno[, 1]
    geno <- geno[, -1]
    founder <- grepl("Founder", col1)
    geno <- as.matrix(geno)
    geno[grepl("N", geno)] <- 3
    geno[geno == "11"] <- 2
    geno[geno %in% c("12", "21")] <- 1
    geno[geno == "22"] <- 0
    geno <- matrix(as.numeric(geno), nrow(geno), ncol(geno))

    if(grepl("realdata", file_i)){
      gds_fn <- sub("_ImputedGenotype.csv", ".gds", file_i)
      gds <- loadGDS(gds_fn)
      gds <- setParents(gds, c("01_NB", "02_OL"), flip = FALSE, mono = FALSE, bi = FALSE)
      p_index <- getParents(gds)$indexes
      smapleorder <- c(p_index, (1:nscan(gds, valid = FALSE))[-p_index])
      geno[smapleorder, ] <- geno
      gdsn <- gdsfmt::add.gdsn(gds@data@handler, "corrected.genotype",
                               val = geno, storage = "bit2", compress = "LZMA_RA", replace = TRUE)
      closeGDS(gds)

    } else {
      geno[geno == 3] <- NA
      geno <- abs(geno - 2)
      geno <- geno[!founder, ]
      pop_file <- sub("_ImputedGenotype.csv", "_SimPop.Rdata", file_i)
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
    }
  }
  if(!grepl("realdata", file_i)){
    write.csv(evaldf, file = paste0(dir_i, "/Eval_Magic.csv"))
  }
}
