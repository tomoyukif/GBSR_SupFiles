dirs <- c("../noError_homoP2_F2", "../realError_homoP2_F2")
          
for(dir_i in dirs){
  files <- list.files(dir_i, ".vcf", full.names = T)
  files <- grep("_LB.vcf", files, invert = TRUE, value = T)
  for(file_i in files){
    if(grepl("realdata", file_i)){
      lb_arg <-  paste0("-jar ", "./LB-Impute/LB-Impute.jar",
                        " -method impute",
                        " -f ", file_i,
                        " -o ", sub("\\.vcf.*", "_LB.vcf", file_i),
                        " -parents 01_NB,02_OL",
                        " -offspringimpute")

    } else {
      lb_arg <-  paste0("-jar ", "./LB-Impute/LB-Impute.jar",
                        " -method impute",
                        " -f ", file_i,
                        " -o ", sub("\\.vcf.*", "_LB.vcf", file_i),
                        " -parents Founder1,Founder2",
                        " -offspringimpute")
    }
    std_out <- system2(command = "java",
                       args = lb_arg,
                       wait = TRUE,
                       stdout = TRUE,
                       stderr = TRUE)
  }
}

# Evaluate correction
library(GBScleanR)
library(SimPop)
source("./compareGeno.R")

for(dir_i in dirs){
  evaldf <- NULL
  files <- list.files(dir_i, "_LB.vcf", full.names = T)
  for(file_i in files){
    gds_fn <- sub(".vcf", ".gds", file_i)
    gbsrVCF2GDS(file_i, gds_fn, force = T)
    gds <- loadGDS(gds_fn)

    if(grepl("realdata", file_i)){
      gds <- setParents(gds, c("01_NB", "02_OL"), flip = TRUE, mono = FALSE, bi = FALSE)

    } else {
      parents <- grep("Founder", getScanID(gds), value = TRUE)
      gds <- setParents(gds, parents, flip = FALSE, mono = FALSE, bi = FALSE)
      geno <- abs(GBScleanR::getGenotype(gds, node = "raw") - 2)
      pop_file <- sub("_LB.gds", "_SimPop.Rdata", gds_fn)
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
    closeGDS(gds)
  }
  if(!grepl("realdata", file_i)){
    write.csv(evaldf, file = paste0(dir_i, "/Eval_LB.csv"))
  }
}
