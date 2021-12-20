library(GBScleanR)
source("./compareGeno.R")

dir <- "../realdata"
if(!dir.exists(dir)){ dir.create(dir) }

read_quantile <- 0.9
gds_fn <- paste0("../realdata/gbs_nbolf2_adQtile", read_quantile, ".gds")
gbsrVCF2GDS(vcf_fn = "../realdata/gbs_nbolf2.vcf.gz",
            out_fn = gds_fn,
            force = TRUE)
gdata <- loadGDS(gds_fn)
p1 <- grep("01_NB", getScanID(gdata), value = TRUE)
p2 <- grep("02_OL", getScanID(gdata), value = TRUE)
gdata <- setParents(object = gdata, parents = c(p1, p2), flip = T, mono = T, bi = T)
nsnp(gdata)
gdata <- countGenotype(gdata)
gdata <- thinMarker(gdata, 75)
gdata <- setCallFilter(gdata,
                       scan_ref_qtile = c(0, read_quantile),
                       scan_alt_qtile = c(0, read_quantile))
gbsrGDS2VCF(gdata, paste0("../realdata/gbs_nbolf2_adQtile", read_quantile, ".vcf"))
gdata <- initScheme(gdata, "pair", matrix(1:2, 2))
gdata <- addScheme(gdata, "self")
gdata <-  estGeno(object = gdata,
                  recomb_rate = 0.04,
                  call_threshold = 0.9,
                  error_rate = 0.0025,
                  iter = 4,
                  het_parent = FALSE,
                  optim = TRUE)
pdf(paste0(dir, "/realdata_genoRatio_adQtile", read_quantile, ".pdf"))
gdata <- countGenotype(gdata, node = "raw")
plotGBSR(gdata, "geno", coord = c(6, 2))
gdata <- countGenotype(gdata, node = "filt")
plotGBSR(gdata, "geno", coord = c(6, 2))
gdata <- countGenotype(gdata, node = "cor")
plotGBSR(gdata, "geno", coord = c(6, 2))
dev.off()
closeGDS(gdata)



