library(SimPop)
library(GBScleanR)
dirs <- grep("noError|realError|realdata", list.dirs("../", recursive = F), value = T)
dirs <- grep("realdata", dirs, value = T, invert = T)
for(dir_i in dirs){
parallel::mclapply(dirs, mc.cores = 6, FUN = function(dir_i){
  files <- list.files(dir_i, ".vcf", full.names = T)
  files <- grep("_LB.vcf", files, invert = TRUE, value = T)
  for(file_i in files){
    if(grepl("realdata", file_i)){
      if(grepl("adQtile", file_i)){ next }
      gds <- sub("\\.vcf.*", ".gds", file_i)
      gbsrVCF2GDS(file_i, gds, force = T)
      gds <- loadGDS(gds)
      header <- cbind("nFounder", 2)
      write.table(header, sub("\\.vcf.*", "_MagicSNP.csv",file_i),
                  quote = TRUE, sep = ",", row.names = F, col.names = F)

      write.table(cbind("marker", t(1:nsnp(gds))),
                  sub("\\.vcf.*", "_MagicSNP.csv",file_i),
                  quote = TRUE, sep = ",", row.names = F, col.names = F, append = T)

      write.table(cbind("chromosome", t(getChromosome(gds))),
                  sub("\\.vcf.*", "_MagicSNP.csv",file_i),
                  quote = TRUE, sep = ",", row.names = F, col.names = F, append = T)

      write.table(cbind("pos(cM)", t(getPosition(gds)*10^-6*4)),
                  sub("\\.vcf.*", "_MagicSNP.csv",file_i),
                  quote = TRUE, sep = ",", row.names = F, col.names = F, append = T)
      read <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(gds@data@handler,
                                                   "annotation/format/AD/data"))
      p_read <- rbind(read[getScanID(gds)  == "01_NB", ],
                      read[getScanID(gds)  == "02_OL", ])
      p_read <- apply(p_read, 1, function(x){
        x <- matrix(x, 2)
        apply(x, 2, paste, collapse = "|")
      })
      p_read <- cbind(c("Founder1", "Founder2"), t(p_read))
      s_read <- read[!getScanID(gds) %in% c("01_NB", "02_OL"), ]
      s_read <- apply(s_read, 1, function(x){
        x <- matrix(x, 2)
        apply(x, 2, paste, collapse = "|")
      })
      s_read <- cbind(paste("Offspring", 1:ncol(s_read), sep = ""),
                      t(s_read))
      write.table(p_read, sub("\\.vcf.*", "_MagicSNP.csv",file_i),
                  quote = TRUE, sep = ",", row.names = F, col.names = F, append = T)
      write.table(s_read, sub("\\.vcf.*", "_MagicSNP.csv",file_i),
                  quote = TRUE, sep = ",", row.names = F, col.names = F, append = T)

    } else {
      pop_file <- sub("\\.vcf.*", "_SimPop.Rdata", file_i)
      load(pop_file)
      header <- cbind("nFounder", getNumInd(pop, 1))
      write.table(header, sub("\\.vcf.*", "_MagicSNP.csv",file_i),
                  quote = TRUE, sep = ",", row.names = F, col.names = F)

      write.table(cbind("marker", t(1:getNumMar(pop))),
                  sub("\\.vcf.*", "_MagicSNP.csv",file_i),
                  quote = TRUE, sep = ",", row.names = F, col.names = F, append = T)

      write.table(cbind("chromosome", t(rep(1, getNumMar(pop)))),
                  sub("\\.vcf.*", "_MagicSNP.csv",file_i),
                  quote = TRUE, sep = ",", row.names = F, col.names = F, append = T)

      write.table(cbind("pos(cM)", t(c(0, cumsum(getGenPos(pop) * 100)))),
                  sub("\\.vcf.*", "_MagicSNP.csv",file_i),
                  quote = TRUE, sep = ",", row.names = F, col.names = F, append = T)

      p_read <- getReads(pop, 1)
      p_read <- sapply(p_read, function(x){
        apply(x, 2, paste, collapse = "|")
      })
      p_read <- cbind(paste("Founder", 1:ncol(p_read)), t(p_read))
      s_read <- getReads(pop, "last")
      s_read <- sapply(s_read, function(x){
        apply(x, 2, paste, collapse = "|")
      })
      s_read <- cbind(paste("Offspring", 1:ncol(s_read), sep = ""), t(s_read))
      write.table(p_read, sub("\\.vcf.*", "_MagicSNP.csv",file_i),
                  quote = TRUE, sep = ",", row.names = F, col.names = F, append = T)
      write.table(s_read, sub("\\.vcf.*", "_MagicSNP.csv",file_i),
                  quote = TRUE, sep = ",", row.names = F, col.names = F, append = T)
    }
  }
})


# Evaluate correction
source("./compareGeno.R")
library(GBScleanR)
library(SimPop)
library(ggplot2)
dirs <- grep("simpop", list.dirs("./ForManuscript", recursive = F), value = T)
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
      gdsfmt::readmode.gdsn(gdsn)
      pdf(sub("_ImputedGenotype.csv", "_Magic_genoRatio.pdf", file_i))
      gds <- countGenotype(gds, node = "cor")
      plot(gds, "geno", coord = c(6, 2))
      dev.off()
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
      pdf(sub("_ImputedGenotype.csv", "_Magic_genoRatio.pdf", file_i), width = 5, height = 2)
      plotGeno(geno, getPhysPos(pop))
      plotGeno(true, getPhysPos(pop))
      dev.off()

    }
  }
  if(!grepl("realdata", file_i)){
    save(evaldf, file = paste0(dir_i, "/Eval_Magic.Rdata"))
    write.csv(evaldf, file = paste0(dir_i, "/Eval_Magic.csv"))
  }
}
