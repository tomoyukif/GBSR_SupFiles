dir.create("../Figures")

library(GBScleanR)
library(SimPop)
library(plyr)
library(ggplot2)
tag_facet2 <- function(p, open = "(", close = ")", tag_pool = letters, x = -Inf, y = Inf,
                       hjust = -0.5, vjust = 1.5, fontface = 2, family = "", ...) {
  gb <- ggplot_build(p)
  lay <- gb$layout$layout
  tags <- cbind(lay, label = paste0(open, tag_pool[lay$PANEL], close), x = x, y = y)
  p <- p + geom_text(data = tags, aes_string(x = "x", y = "y", label = "label"), ..., hjust = hjust,
                     vjust = vjust, fontface = fontface, family = family, inherit.aes = FALSE)
  invisible(p)
}

dirs <- grep("noError|realError", list.dirs("../", recursive = F), value = T)
file.create("../EvalSummary.csv")
header <- T
for(dir_i in dirs){
  dataset <- sub(".*/", "", dir_i)
  files <- list.files(dir_i, "Eval.*\\.csv")
  files <- grep("hap", files, invert = T, value = T)
  for(file_i in files){
    tool <- "gbsr"
    if(grepl("LB", file_i)){tool = "LB"}
    if(grepl("Magic", file_i)){tool = "magic"}
    if(grepl("noOptim", file_i)){tool = "gbsr_noOptim"}
    df <- read.csv(paste0(dir_i,"/", file_i), row.names = 1)
    if(tool %in% c("LB", "magic")){
      df[, 1] <- sub("_LB.*|_Imputed.*", "", sub(".*\\/", "", df[, 1]))
    } else {
      df[, 1] <- sub("\\.gds", "", sub(".*\\/", "", df[, 1]))
    }
    df <- cbind(dataset = dataset, tool = tool, df)
    write.table(df, "../EvalSummary.csv", sep = ",", quote = F,
                append = T, col.names = header, row.names = F)
    if(header){header <- F}
  }
}
df <- read.table("../EvalSummary.csv", sep = ",", quote = "", header = T, stringsAsFactors = T)

# Correction accuracy in the simulation datasets
for(ind in c(10, 100, 1000)){
  target_dataset <- "noError_homoP2_F2|noError_hetP2_F2|noError_homoP8_RIL"
  sub_df <- subset(df, subset = grepl(paste0("ind", ind, "_"), condition) & grepl(target_dataset, dataset))
  sub_df <- droplevels(sub_df)
  sub_df$dataset <- ordered(sub_df$dataset, levels = c("noError_homoP2_F2", "noError_hetP2_F2", "noError_homoP8_RIL"))
  levels(sub_df$dataset) <- c("homoP2_F2", "hetP2_F2", "homoP8_RIL")
  sub_df$mar <- sub("_fread.*", "", sub(".*_mar", "", as.character(sub_df$condition)))
  sub_df$mar <- factor(sub_df$mar, levels = c(50, 100, 1000))
  sub_df$oread <- sub("_.*", "", sub("ind.*_oread", "", as.character(sub_df$condition)))
  sub_df$fread <- sub("_.*", "", sub("ind.*_fread", "", as.character(sub_df$condition)))
  sub_df <- subset(sub_df, subset = !(tool == "LB" & grepl("het|P8", dataset)))
  sub_df <- sub_df[order(as.numeric(sub_df$oread)), ]
  sub_df$strip <- paste0(sub_df$dataset, "Ind=",sub_df$mar)
  sub_df$tool <- factor(sub_df$tool, levels = c("LB", "gbsr_noOptim", "magic", "gbsr"))
  sub_df <- sub_df[order(as.numeric(sub_df$tool)), ]
  pdf(paste0("../Figures/Figure_eval_simuldata_ind", ind, ".pdf"), width = 3.54*2, height = 5)
  p <- ggplot(sub_df) +
    geom_path(aes(x = oread, y = correct, color = tool, linetype = fread, group = paste(fread, tool), alpha = fread)) +
    geom_point(aes(x = oread, y = correct, color = tool, alpha = fread)) +
    facet_grid(facets = mar ~ dataset, scale = "free") +
    scale_alpha_manual(breaks = c("3", "1"), values = c(1, 0.6)) +
    scale_linetype_manual(name = "Founder read depth",
                          breaks = c("3", "1"),
                          values = c(1, 2)) +
    scale_color_manual(name = "Tool",
                       breaks = c("gbsr", "gbsr_noOptim", "LB", "magic"),
                       values = c("blue", "green", "gray", "magenta"),
                       labels = c("GBSR", "GBSR_noOpt", "LB", "Magic")) +
    xlab("Offspring read depth") +
    ylab("Correct call rate") +
    theme(legend.position = "top",
          legend.box = "vertical",
          legend.spacing.y = unit(0, "mm"),
          legend.box.margin = margin(r = 40),
          legend.box.just = "left",
          legend.spacing.x = unit(1, "mm"),
          axis.text = element_text(size = 8),
          strip.text = element_text(size = 8),
          axis.title = element_text(size = 10)) +
    guides(alpha = "none")
  p <- tag_facet2(p, tag_pool = LETTERS, hjust = -0.2)
  print(p)
  dev.off()

}

# Correction accuracy in the realdata-based simulation datasets
for(ind in c(10, 100, 1000)){
  target_dataset <- "realError_homoP2_F2|realError_hetP2_F2|realError_homoP8_RIL"
  sub_df <- subset(df, subset = grepl(paste0("ind", ind, "_"), condition) & grepl(target_dataset, dataset))
  sub_df <- droplevels(sub_df)
  sub_df$dataset <- ordered(sub_df$dataset, levels = c("realError_homoP2_F2", "realError_hetP2_F2", "realError_homoP8_RIL"))
  levels(sub_df$dataset) <- c("homoP2_F2", "hetP2_F2", "homoP8_RIL")
  sub_df$oread <- sub("_.*", "", sub("ind.*_oread", "", as.character(sub_df$condition)))
  sub_df$fread <- sub("_.*", "", sub("ind.*_fread", "", as.character(sub_df$condition)))
  sub_df <- subset(sub_df, subset = !(tool == "LB" & grepl("het|P8", dataset)))
  sub_df <- sub_df[order(as.numeric(sub_df$oread)), ]
  sub_df$strip <- paste0(sub_df$dataset, "Ind=",sub_df$mar)
  sub_df$tool <- factor(sub_df$tool, levels = c("LB", "gbsr_noOptim", "magic", "gbsr"))
  sub_df <- sub_df[order(as.numeric(sub_df$tool)), ]

  pdf(paste0("../Figures/Figure_eval_realbased_ind", ind, ".pdf"), width = 3.54, height = 5)
  p <- ggplot(sub_df) +
    geom_path(aes(x = oread, y = correct, color = tool, linetype = fread, group = paste(fread, tool), alpha = fread)) +
    geom_point(aes(x = oread, y = correct, color = tool, alpha = fread)) +
    facet_wrap(facets = ~ dataset, ncol = 1, scale = "free") +
    scale_alpha_manual(breaks = c("3", "1"), values = c(1, 0.6)) +
    scale_linetype_manual(name = "Founder read depth", breaks = c("3", "1"), values = c(1, 2)) +
    scale_color_manual(name = "Tool",
                       breaks = c("gbsr", "gbsr_noOptim", "LB", "magic"),
                       values = c("blue", "green", "gray", "magenta"),
                       labels = c("GBSR", "GBSR_noOpt", "LB", "Magic")) +
    xlab("Offspring read depth") +
    ylab("Correct call rate") +
    theme(legend.position = "top",
          legend.box = "vertical",
          legend.spacing.y = unit(0, "mm"),
          legend.box.margin = margin(r = 40),
          legend.box.just = "left",
          legend.spacing.x = unit(1, "mm"),
          axis.text = element_text(size = 8),
          strip.text = element_text(size = 8),
          axis.title = element_text(size = 10)) +
    guides(alpha = "none")
  p <- tag_facet2(p, tag_pool = LETTERS, hjust = -0.2)
  print(p)
  dev.off()
}

#################################################
# Allele read bias
library(GBScleanR)
gds <- loadGDS("../realdata/gbs_nbolf2_adQtile0.9.gds"))
offspring <- grep("NB|OL|F1", getScanID(gds), value = TRUE, invert = TRUE)
gds <- setScanFilter(gds, id = offspring)
gds <- countGenotype(gds)
gds <- thinMarker(gds, 75)
gds <- setFiltGenotype(gds)
gds <- countRead(gds)
chr <- getChromosome(gds)
ad_dist <- getCountReadRef(gds, prop = TRUE)
p <- ggplot(data.frame(x = ad_dist, y = chr), aes(x = x)) +
  geom_histogram() +
  facet_wrap(~ y, 3, 4) +
  xlab("Allele read ratio per marker") +
  scale_x_continuous(breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
  theme(axis.text = element_text(size = 10),
        strip.text = element_text(size = 10),
        axis.title = element_text(size = 10))

pdf("../Figures/Figure_readbias.pdf", width = 5, height = 5)
print(p)
dev.off()
closeGDS(gds)


# Genotype ratio per marker in the realdata-based simulation datasets
for(ind in c(10, 100, 1000)){
  gds_fn <- paste0("../realError_homoP2_F2/ind",ind,"_mar620_fread3_oread3_gbsr.gds")
  gds <- loadGDS(gds_fn)
  gbsr_geno <- abs(GBScleanR::getGenotype(gds, node = "cor") - 2)
  gbsr_geno <- gbsr_geno[-(1:2), ]
  closeGDS(gds)
  pop_file <- sub("_gbsr.gds", "_SimPop.Rdata", gds_fn)
  load(pop_file)
  pos <- SimPop::getPhysPos(pop)
  true <- SimPop::getGenotype(pop, gen = "last", fam = "all", sib = "all")
  true <- t(sapply(true, colSums))
  gbsr <- compareGeno(gbsr_geno, true)
  magic_fn <- paste0("../realError_homoP2_F2/ind",ind,"_mar620_fread3_oread3_ImputedGenotype.csv")
  magic_geno <- read.table(magic_fn, header = FALSE, sep = ",", skip = 4, stringsAsFactors = FALSE)
  founder <- grepl("Founder", magic_geno[, 1])
  magic_geno <- magic_geno[, -1]
  magic_geno <- as.matrix(magic_geno)
  magic_geno[grepl("N", magic_geno)] <- NA
  magic_geno[magic_geno == "11"] <- 2
  magic_geno[magic_geno %in% c("12", "21")] <- 1
  magic_geno[magic_geno == "22"] <- 0
  magic_geno <- matrix(as.numeric(magic_geno), nrow(magic_geno), ncol(magic_geno))
  magic_geno <- abs(magic_geno - 2)
  magic_nonallelic <- magic_geno[founder, ][1, ] == magic_geno[founder, ][2, ]
  magic_geno <- magic_geno[!founder, ]
  magic <- compareGeno(magic_geno, true)
  gds_fn <- paste0("../realError_homoP2_F2/ind", ind, "_mar620_fread3_oread3_LB.gds")
  lbgds <- loadGDS(gds_fn)
  lb_geno <- abs(GBScleanR::getGenotype(lbgds, node = "raw") - 2)
  lb_nonallelic <- lb_geno[1, ] == lb_geno[2, ]
  lb_nonallelic[is.na(lb_nonallelic)] <- TRUE
  lb_geno <- lb_geno[-(1:2), ]
  closeGDS(lbgds)
  lb <- compareGeno(lb_geno, true)

  ratio <- apply(gbsr_geno, 2, function(x)table(factor(x, levels = 0:2)) / sum(!is.na(x)))
  df1 <- data.frame(pos, t(ratio), nonallelic = FALSE, Tool = "GBSR")
  ratio <- apply(true, 2, function(x)table(factor(x, levels = 0:2)) / sum(!is.na(x)))
  df1 <- rbind(df1, data.frame(pos, t(ratio), nonallelic = FALSE, Tool = "True"))
  ratio <- apply(magic_geno, 2, function(x)table(factor(x, levels = 0:2)) / sum(!is.na(x)))
  df1 <- rbind(df1, data.frame(pos, t(ratio), nonallelic = magic_nonallelic, Tool = "Magic"))
  ratio <- apply(lb_geno, 2, function(x)table(factor(x, levels = 0:2)) / sum(!is.na(x)))
  df1 <- rbind(df1, data.frame(pos, t(ratio), nonallelic = lb_nonallelic, Tool = "LB"))
  names(df1) <- c("pos", "Ref", "Het", "Alt", "Nonallelic", "Tool")

  df1$Tool <- factor(df1$Tool, levels = c("True", "GBSR", "LB", "Magic"))
  df1 <- df1[order(as.numeric(df1$Tool)), ]

  df2 <- data.frame(pos, Tool = "GBSR", miscall = gbsr$mar$miscall)
  df2 <- rbind(df2, data.frame(pos, Tool = "Magic", miscall = magic$mar$miscall))
  df2 <- rbind(df2, data.frame(pos, Tool = "LB", miscall = lb$mar$miscall))
  df2$Tool <- factor(df2$Tool, levels = c("GBSR", "LB", "Magic"))
  df2 <- df2[order(as.numeric(df2$Tool)), ]

  df1 <- tidyr::pivot_longer(df1, -c(pos, Tool, Nonallelic), names_to = "Genotype", values_to = "ratio")
  df1$Tool <- factor(df1$Tool, levels = c("True", "GBSR", "LB", "Magic"))
  df1 <- df1[order(as.numeric(df1$Tool)), ]
  df1_nonallelic <- df1[df1$Nonallelic, ]
  df1 <- df1[!df1$Nonallelic, ]
  p <- ggplot(df1) +
    geom_vline(data = df1_nonallelic, mapping = aes(xintercept = pos * 10^-6), color = "darkgray") +
    geom_line(aes(y = ratio, x = pos * 10^-6, group = Genotype, color = Genotype)) +
    geom_line(data = df2, mapping = aes(x = pos * 10^-6, y = miscall),
              stat = "identity", col = "darkorange", alpha = 0.7) +
    facet_wrap(facets = ~ Tool, ncol = 1) +
    xlab("Physical position (Mb)") +
    ylab("Proportion of genotypes and miscall rate") +
    scale_x_continuous(expand = expansion(0, 0.5), limits = c(0, NA)) +
    scale_y_continuous(limits = c(0, 1)) +
    scale_color_manual(values = c("darkgreen", "magenta", "blue"), breaks = c("Ref", "Het", "Alt")) +
    theme(legend.spacing.y = unit(0, "mm"),
          axis.text = element_text(size = 8),
          strip.text = element_text(size = 10),
          axis.title = element_text(size = 10),
          legend.position = "top",
          legend.box = "vertical")
  p <- tag_facet2(p, tag_pool = LETTERS, hjust = -0.2)

  pdf(paste0("../Figures/Figure_genoratio_realbased_ind",ind,".pdf"), width = 3.54, height = 6)
  print(p)
  dev.off()
}

# Correction of the real data
library(GBScleanR)
gds_fn <- "../realdata/gbs_nbolf2_adQtile0.9.gds"
gds <- loadGDS(gds_fn)
gds <- setParents(gds, c("01_NB", "02_OL"), flip = T, mono = T, bi = T)
gds <- countGenotype(gds)
gds <- thinMarker(gds, 75)
gbsr_geno <- abs(GBScleanR::getGenotype(gds, node = "cor") - 2)
p_geno <- abs(GBScleanR::getGenotype(gds, node = "cor", parents = "only") - 2)

closeGDS(gds)
magic_fn <- "../realdata/gbs_nbolf2_adQtile0.9_ImputedGenotype.csv"
magic_geno <- read.table(magic_fn, header = FALSE, sep = ",", skip = 4, stringsAsFactors = FALSE)
founder <- grepl("Founder", magic_geno[, 1])
magic_geno <- magic_geno[, -1]
magic_geno <- as.matrix(magic_geno)
magic_geno[grepl("N", magic_geno)] <- NA
magic_geno[magic_geno == "11"] <- 2
magic_geno[magic_geno %in% c("12", "21")] <- 1
magic_geno[magic_geno == "22"] <- 0
magic_geno <- matrix(as.numeric(magic_geno), nrow(magic_geno), ncol(magic_geno))
magic_geno <- abs(magic_geno - 2)
magic_geno <- magic_geno[!founder, ]
gds_fn <- "../realdata/gbs_nbolf2_adQtile0.9_LB.gds"
lbgds <- loadGDS(gds_fn)
lbgds <- setParents(lbgds, c("01_NB", "02_OL"), flip = T, mono = F, bi = F)
lbgds <- flipData(lbgds)
lb_geno <- abs(GBScleanR::getGenotype(lbgds, node = "raw") - 2)
closeGDS(lbgds)
pos <- getPosition(lbgds)
chr <- getChromosome(lbgds)

ratio <- apply(gbsr_geno, 2, function(x)table(factor(x, levels = 0:2)) / sum(!is.na(x)))
df3 <- data.frame(chr, pos, t(ratio), Tool = "GBSR")
ratio <- apply(lb_geno, 2, function(x)table(factor(x, levels = 0:2)) / sum(!is.na(x)))
df3 <- rbind(df3, data.frame(chr, pos, t(ratio), Tool = "LB"))
ratio <- apply(magic_geno, 2, function(x)table(factor(x, levels = 0:2)) / sum(!is.na(x)))
df3 <- rbind(df3, data.frame(chr, pos, t(ratio), Tool = "Magic"))
names(df3) <- c("chr","pos", "Ref", "Het", "Alt", "Tool")

df3 <- tidyr::pivot_longer(df3, -c(chr, pos, Tool), names_to = "Genotype", values_to = "ratio")
p <- ggplot(subset(df3, subset = chr %in% c(1, 7))) +
  geom_line(aes(y = ratio, x = pos * 10^-6, group = Genotype, color = Genotype)) +
  facet_grid(facets = chr ~ Tool, scales = "free_x") +
  xlab("Physical position (Mb)") +
  ylab("Proportion of genotypes and miscall rate") +
  scale_x_continuous(expand = expansion(0, 0.5), limits = c(0, NA)) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_manual(values = c("darkgreen", "magenta", "blue"), breaks = c("Ref", "Het", "Alt")) +
  theme(legend.spacing.y = unit(0, "mm"),
        axis.text = element_text(size = 8),
        strip.text = element_text(size = 10),
        axis.title = element_text(size = 10),
        legend.position = "top",
        legend.box = "vertical")
p <- tag_facet2(p, tag_pool = LETTERS, hjust = -0.2)

pdf("../Figures/Figure_realdata_1_7.pdf", height = 3)
print(p)
dev.off()

p <- ggplot(df3) +
  geom_line(aes(y = ratio, x = pos * 10^-6, group = Genotype, color = Genotype)) +
  facet_grid(facets = chr ~ Tool, scales = "free_x") +
  xlab("Physical position (Mb)") +
  ylab("Proportion of genotypes and miscall rate") +
  scale_x_continuous(expand = expansion(0, 0.5), limits = c(0, NA)) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_manual(values = c("darkgreen", "magenta", "blue"), breaks = c("Ref", "Het", "Alt")) +
  theme(legend.spacing.y = unit(0, "mm"),
        axis.text = element_text(size = 8),
        strip.text = element_text(size = 10),
        axis.title = element_text(size = 10),
        legend.position = "top",
        legend.box = "vertical")

pdf("../Figures/Figure_realdata.pdf", height = 10)
print(p)
dev.off()


# Evaluate different parameter settings
library(cowplot)
options(scipen = 10^6)
df <- read.table("../ParamTest/ParamTest.csv", sep =",", header = T)
df$test <- factor(df$test, levels = c("iter", "call", "recomb", "error"))
p1 <- ggplot(subset(df, subset = test == "iter" & value != 1)) +
  geom_point(aes(x = correct, y = miscall, shape = factor(as.integer(value))), col = "blue") +
  geom_path(aes(x = correct, y = miscall, group = test), col = "blue") +
  facet_wrap(facets = ~ test, ncol = 1, scale = "free",
             labeller = labeller(test = c(iter = "Iteration number", recomb = "Recombination rate", call = "Call threshold", error = "Sequence error rate"))) +
  xlab("Correct call rate") +
  ylab("Miscall rate") +
  scale_shape_manual(name = "Iteration Number", values = 1:5) +
  theme(legend.position = "top",
        legend.title = element_text(size = 8),
        legend.spacing.y = unit(0, "mm"),
        legend.box.margin = margin(r = 40, t = 0, b = 0),
        legend.box.just = "left",
        legend.spacing.x = unit(1, "mm"),
        axis.text = element_text(size = 8),
        strip.text = element_text(size = 10),
        axis.title = element_text(size = 10)) +
  guides(shape = guide_legend(nrow = 2, byrow = T, title.position = "top", title.hjust = 0.5))

p2 <- ggplot(subset(df, subset = test == "recomb" & value != 1)) +
  geom_point(aes(x = correct, y = miscall, shape = factor(value)), col = "blue") +
  geom_path(aes(x = correct, y = miscall, group = test), col = "blue") +
  facet_wrap(facets = ~ test, ncol = 1, scale = "free",
             labeller = labeller(test = c(iter = "Iteration number", recomb = "Recombination rate", call = "Call threshold", error = "Sequence error rate"))) +
  xlab("Correct call rate") +
  ylab("Miscall rate") +
  scale_shape_manual(name = expression("Genetic distance per Mb "~italic(E^d)), values = 1:6) +
  theme(legend.position = "top",
        legend.box = "vertical",
        legend.title = element_text(size = 8),
        legend.spacing.y = unit(0, "mm"),
        legend.box.margin = margin(r = 40, t = 0, b = 0),
        legend.box.just = "left",
        legend.spacing.x = unit(1, "mm"),
        axis.text = element_text(size = 8),
        strip.text = element_text(size = 10),
        axis.title = element_text(size = 10)) +
  guides(shape = guide_legend(nrow = 2, byrow = T, title.position = "top", title.hjust = 0.5))

p3 <- ggplot(subset(df, subset = test == "call" & value != 1)) +
  geom_point(aes(x = correct, y = miscall, shape = factor(value)), col = "blue") +
  geom_path(aes(x = correct, y = miscall, group = test), col = "blue") +
  facet_wrap(facets = ~ test, ncol = 1, scale = "free",
             labeller = labeller(test = c(iter = "Iteration number", recomb = "Recombination rate", call = "Call threshold", error = "Sequence error rate"))) +
  xlab("Correct call rate") +
  ylab("Miscall rate") +
  scale_shape_manual(name = expression("Genotype call threshold "~italic(P^call)), values = 1:6) +
  theme(legend.position = "top",
        legend.box = "vertical",
        legend.title = element_text(size = 8),
        legend.spacing.y = unit(0, "mm"),
        legend.box.margin = margin(r = 40, t = 0, b = 0),
        legend.box.just = "left",
        legend.spacing.x = unit(1, "mm"),
        axis.text = element_text(size = 8),
        strip.text = element_text(size = 10),
        axis.title = element_text(size = 10)) +
  guides(shape = guide_legend(nrow = 2, byrow = T, title.position = "top", title.hjust = 0.5))

p4 <- ggplot(subset(df, subset = test == "error" & value != 1)) +
  geom_point(aes(x = correct, y = miscall, shape = factor(value)), col = "blue") +
  geom_path(aes(x = correct, y = miscall, group = test), col = "blue") +
  facet_wrap(facets = ~ test, ncol = 1, scale = "free",
             labeller = labeller(test = c(iter = "Iteration number", recomb = "Recombination rate", call = "Call threshold", error = "Sequence error rate"))) +
  xlab("Correct call rate") +
  ylab("Miscall rate") +
  scale_shape_manual(name = expression("Sequence error rate "~italic(e^seq)), values = 1:6) +
  theme(legend.position = "top",
        legend.box = "vertical",
        legend.title = element_text(size = 8),
        legend.spacing.y = unit(0, "mm"),
        legend.box.margin = margin(r = 40, t = 0, b = 0),
        legend.box.just = "left",
        legend.spacing.x = unit(1, "mm"),
        axis.text = element_text(size = 8),
        strip.text = element_text(size = 10),
        axis.title = element_text(size = 10)) +
  guides(shape = guide_legend(nrow = 2, byrow = T, title.position = "top", title.hjust = 0.5))
p <- plot_grid(p1, p2, p3, p4, labels = "AUTO")

pdf("../Figures/Figure_paramtest.pdf", width = 3.54*2, height = 7)
print(p)
dev.off()


# Evaluate founder genotype estimation
dirs <- c("../noError_homoP2_F2/",
          "../noError_hetP2_F2/",
          "../noError_homoP8_RIL/")
evaldf <- NULL
for(dir_i in dirs){
  scenario <- switch (dir_i,
                      "../noError_homoP2_F2/" = "homoP2_F2",
                      "../noError_hetP2_F2/" = "hetP2_F2",
                      "../noError_homoP8_RIL/" = "homoP8_RIL"
  )
  gbsr_files <- list.files(dir_i, "*_gbsr.gds", full.names = T)
  gbsr_no_files <- list.files(dir_i, "*_gbsr_noOptim.gds", full.names = T)
  magic_files <- list.files(dir_i, "_ImputedGenotype.csv", full.names = T)
  pop_files <- list.files(dir_i, "SimPop", full.names = T)
  datasets <- sub(".*/", "", sub("_gbsr.gds", "", gbsr_files))
  for(i_data in datasets){
    gds <- loadGDS(grep(paste0(i_data, "_"), gbsr_files, value = T))
    parents <- grep("Founder", getScanID(gds), value = T)
    gds <- setParents(gds, parents = parents, flip = F, mono = F, bi = F)
    p_gbsr <- GBScleanR::getGenotype(gds, node = "parents")
    p_gbsr <- p_gbsr[c(T, F), ] + p_gbsr[c(F, T), ]

    gds_no <- loadGDS(grep(paste0(i_data, "_"), gbsr_no_files, value = T))
    parents <- grep("Founder", getScanID(gds_no), value = T)
    gds_no <- setParents(gds_no, parents = parents, flip = F, mono = F, bi = F)
    p_gbsr_no <- GBScleanR::getGenotype(gds_no, node = "parents")
    p_gbsr_no <- p_gbsr_no[c(T, F), ] + p_gbsr_no[c(F, T), ]

    magic <- read.table(grep(paste0(i_data, "_"), magic_files, value = T),
                        sep = ",", skip = 4, row.names = 1)
    p_magic <- magic[grepl("Founder", rownames(magic)), ]
    if(grepl("het", dir_i)){
      p_magic <- plyr::alply(p_magic, 1, function(x){
        x <- do.call("rbind", strsplit(as.character(x), ""))
        x_dim <- dim(x)
        x <- matrix(as.numeric(x), x_dim[1], x_dim[2])
        return(t(x))
      })
      p_magic <- do.call("rbind", p_magic)
      p_magic <- p_magic - 1
    } else {
      p_magic <- plyr::alply(p_magic, 1, function(x){
        x <- as.numeric(x) - 1
        x <- rbind(x, x)
        return(x)
      })
      p_magic <- do.call("rbind", p_magic)
    }
    p_magic <- p_magic[c(T, F), ] + p_magic[c(F, T), ]


    load(grep(paste0(i_data, "_"), pop_files, value = T))
    p_true <- SimPop::getGenotype(pop, gen = 1, fam = "all", sib = "all", fmt = "matrix")
    p_true <- p_true[c(T, F), ] + p_true[c(F, T), ]

    eval_gbsr <- compareGeno(p_gbsr, p_true)$ind
    eval_gbsr <- subset(eval_gbsr, select = c("correct", "miscall"))
    eval_gbsr_no <- compareGeno(p_gbsr_no, p_true)$ind
    eval_gbsr_no <- subset(eval_gbsr_no, select = c("correct", "miscall"))
    eval_magic <- compareGeno(p_magic, p_true)$ind
    eval_magic <- subset(eval_magic, select = c("correct", "miscall"))
    gbsr_df <- cbind(scenario = scenario,
                     dataset = i_data,
                     tool = "gbsr",
                     data.frame(t(apply(eval_gbsr, 2, mean))),
                     data.frame(t(apply(eval_gbsr, 2, sd))))
    names(gbsr_df)[6:7] <- c("correct_sd", "miscall_sd")
    gbsr_no_df <- cbind(scenario = scenario,
                        dataset = i_data,
                        tool = "gbsr_noOptim",
                        data.frame(t(apply(eval_gbsr_no, 2, mean))),
                        data.frame(t(apply(eval_gbsr_no, 2, sd))))
    names(gbsr_no_df)[6:7] <- c("correct_sd", "miscall_sd")
    magic_df <- cbind(scenario = scenario,
                      dataset = i_data,
                      tool = "magic",
                      data.frame(t(apply(eval_magic, 2, mean))),
                      data.frame(t(apply(eval_magic, 2, sd))))
    names(magic_df)[6:7] <- c("correct_sd", "miscall_sd")
    evaldf <- rbind(evaldf,
                    gbsr_df, gbsr_no_df,
                    magic_df)
  }
}

write.csv(evaldf, file = "../Eval_founders_noError.csv")

evaldf <- read.csv(file = "../Eval_founders_noError.csv")
library(ggplot2)
evaldf$mar <- sub("_.*", "", sub(".*mar", "", evaldf$dataset))
evaldf$oread <- sub(".*_oread", "", evaldf$dataset)
evaldf$fread <- sub(".*fread", "",sub("_oread.*", "", evaldf$dataset))
evaldf$ind <- sub("_.*", "", sub("ind", "", evaldf$dataset))
evaldf$indmar <- paste0("Ind = ", evaldf$ind, "\nMar = ", evaldf$mar)
evaldf$scenario <- factor(evaldf$scenario, levels = c("homoP2_F2", "hetP2_F2", "homoP8_RIL"))
evaldf$oread <- factor(evaldf$oread, levels = c("0.1", "0.25", "0.5", "0.75", "1", "2", "3"))
evaldf$fread <- factor(evaldf$fread, levels = c("3", "1"))
evaldf <- evaldf[order(as.numeric(evaldf$ind), as.numeric(evaldf$mar)), ]
evaldf$indmar <- factor(evaldf$indmar, levels = unique(evaldf$indmar))
evaldf$tool <- factor(evaldf$tool, levels = c("gbsr_noOptim", "magic", "gbsr"))
evaldf <- evaldf[order(as.numeric(evaldf$tool)), ]
p <- ggplot(evaldf) +
  geom_path(aes(x = oread, y = correct, color = tool, group = paste(tool, fread), linetype = fread)) +
  geom_point(aes(x = oread, y = correct, color = tool)) +
  facet_grid(facets = indmar ~ scenario, scale = "free") +
  scale_color_manual(name = "Tool",
                     breaks = c("gbsr", "gbsr_noOptim", "magic"),
                     values = c("blue", "green", "magenta"),
                     labels = c("GBSR", "GBSR_noOpt", "Magic")) +
  scale_linetype_manual(name = "Founder read depth",
                        breaks = c("3", "1"),
                        values = c(1, 2)) +
  ylab("Correct call rate") +
  xlab("Read depth") +
  theme(legend.position = "top",
        legend.box = "vertical",
        legend.spacing.y = unit(0, "mm"),
        legend.box.margin = margin(r = 40, t = 0, b = 0),
        legend.box.just = "left",
        legend.spacing.x = unit(1, "mm"),
        axis.text = element_text(size = 8),
        strip.text = element_text(size = 10),
        axis.title = element_text(size = 10)) +
  guides(shape = guide_legend(nrow = 1, byrow = T))
pdf("../Figures/Figure_FounderGenotype_noError.pdf", width = 3.54 * 2, height = 10)
print(p)
dev.off()

dirs <- c("../realError_homoP2_F2/",
          "../realError_hetP2_F2/",
          "../realError_homoP8_RIL/")
evaldf <- NULL
for(dir_i in dirs){
  scenario <- switch (dir_i,
                      "../realError_homoP2_F2/" = "homoP2_F2",
                      "../realError_hetP2_F2/" = "hetP2_F2",
                      "../realError_homoP8_RIL/" = "homoP8_RIL"
  )
  gbsr_files <- list.files(dir_i, "*_gbsr.gds", full.names = T)
  gbsr_no_files <- list.files(dir_i, "*_gbsr_noOptim.gds", full.names = T)
  magic_files <- list.files(dir_i, "_ImputedGenotype.csv", full.names = T)
  pop_files <- list.files(dir_i, "SimPop", full.names = T)
  datasets <- sub(".*/", "", sub("_gbsr.gds", "", gbsr_files))
  for(i_data in datasets){
    gds <- loadGDS(grep(paste0(i_data, "_"), gbsr_files, value = T))
    parents <- grep("Founder", getScanID(gds), value = T)
    gds <- setParents(gds, parents = parents, flip = F, mono = F, bi = F)
    p_gbsr <- GBScleanR::getGenotype(gds, node = "parents")
    p_gbsr <- p_gbsr[c(T, F), ] + p_gbsr[c(F, T), ]

    gds_no <- loadGDS(grep(paste0(i_data, "_"), gbsr_no_files, value = T))
    parents <- grep("Founder", getScanID(gds_no), value = T)
    gds_no <- setParents(gds_no, parents = parents, flip = F, mono = F, bi = F)
    p_gbsr_no <- GBScleanR::getGenotype(gds_no, node = "parents")
    p_gbsr_no <- p_gbsr_no[c(T, F), ] + p_gbsr_no[c(F, T), ]

    magic <- read.table(grep(paste0(i_data, "_"), magic_files, value = T),
                        sep = ",", skip = 4, row.names = 1)
    p_magic <- magic[grepl("Founder", rownames(magic)), ]
    if(grepl("het", dir_i)){
      p_magic <- plyr::alply(p_magic, 1, function(x){
        x <- do.call("rbind", strsplit(as.character(x), ""))
        x_dim <- dim(x)
        x <- matrix(as.numeric(x), x_dim[1], x_dim[2])
        return(t(x))
      })
      p_magic <- do.call("rbind", p_magic)
      p_magic <- p_magic - 1
    } else {
      p_magic <- plyr::alply(p_magic, 1, function(x){
        x <- as.numeric(x) - 1
        x <- rbind(x, x)
        return(x)
      })
      p_magic <- do.call("rbind", p_magic)
    }
    p_magic <- p_magic[c(T, F), ] + p_magic[c(F, T), ]


    load(grep(paste0(i_data, "_"), pop_files, value = T))
    p_true <- SimPop::getGenotype(pop, gen = 1, fam = "all", sib = "all", fmt = "matrix")
    p_true <- p_true[c(T, F), ] + p_true[c(F, T), ]

    eval_gbsr <- compareGeno(p_gbsr, p_true)$ind
    eval_gbsr <- subset(eval_gbsr, select = c("correct", "miscall"))
    eval_gbsr_no <- compareGeno(p_gbsr_no, p_true)$ind
    eval_gbsr_no <- subset(eval_gbsr_no, select = c("correct", "miscall"))
    eval_magic <- compareGeno(p_magic, p_true)$ind
    eval_magic <- subset(eval_magic, select = c("correct", "miscall"))
    gbsr_df <- cbind(scenario = scenario,
                     dataset = i_data,
                     tool = "gbsr",
                     data.frame(t(apply(eval_gbsr, 2, mean))),
                     data.frame(t(apply(eval_gbsr, 2, sd))))
    names(gbsr_df)[6:7] <- c("correct_sd", "miscall_sd")
    gbsr_no_df <- cbind(scenario = scenario,
                     dataset = i_data,
                     tool = "gbsr_noOptim",
                     data.frame(t(apply(eval_gbsr_no, 2, mean))),
                     data.frame(t(apply(eval_gbsr_no, 2, sd))))
    names(gbsr_no_df)[6:7] <- c("correct_sd", "miscall_sd")
    magic_df <- cbind(scenario = scenario,
          dataset = i_data,
          tool = "magic",
          data.frame(t(apply(eval_magic, 2, mean))),
          data.frame(t(apply(eval_magic, 2, sd))))
    names(magic_df)[6:7] <- c("correct_sd", "miscall_sd")
    evaldf <- rbind(evaldf,
                    gbsr_df,
                    gbsr_no_df,
                    magic_df)
  }
}
write.csv(evaldf, file = "../Eval_founders_realError.csv")

evaldf <- read.csv(file = "../Eval_founders_realError.csv")
library(ggplot2)
evaldf$oread <- sub(".*_oread", "", evaldf$dataset)
evaldf$fread <- sub(".*fread", "",sub("_oread.*", "", evaldf$dataset))
evaldf$ind <- sub("_.*", "", sub("ind", "", evaldf$dataset))
evaldf$scenario <- factor(evaldf$scenario, levels = c("homoP2_F2", "hetP2_F2", "homoP8_RIL"))
evaldf <- evaldf[order(as.numeric(evaldf$oread)), ]
evaldf$oread <- factor(evaldf$oread, levels = c("0.1", "0.25", "0.5", "0.75", "1", "2", "3"))
evaldf$fread <- factor(evaldf$fread, levels = c("3", "1"))
evaldf$tool <- factor(evaldf$tool, levels = c("gbsr_noOptim", "magic", "gbsr"))
evaldf <- evaldf[order(as.numeric(evaldf$tool)), ]
p <- ggplot(evaldf) +
  geom_path(aes(x = oread, y = correct, color = tool, group = paste(tool, fread), linetype = fread)) +
  geom_point(aes(x = oread, y = correct, color = tool)) +
  facet_grid(facets = ind ~ scenario, scale = "free") +
  scale_color_manual(name = "Tool",
                     breaks = c("gbsr", "gbsr_noOptim", "magic"),
                     values = c("blue", "green", "magenta"),
                     labels = c("GBSR", "GBSR_noOpt", "Magic")) +
  scale_linetype_manual(name = "Founder read depth",
                        breaks = c("3", "1"),
                        values = c(1, 2)) +
  ylab("Correct call rate") +
  xlab("Read depth") +
  theme(legend.position = "top",
        legend.box = "vertical",
        legend.spacing.y = unit(0, "mm"),
        legend.box.margin = margin(r = 40, t = 0, b = 0),
        legend.box.just = "left",
        legend.spacing.x = unit(1, "mm"),
        axis.text = element_text(size = 8),
        strip.text = element_text(size = 10),
        axis.title = element_text(size = 10)) +
  guides(shape = guide_legend(nrow = 1, byrow = T))
p <- tag_facet2(p, tag_pool = LETTERS, hjust = -0.2)
pdf("../Figures/Figure_FounderGenotype_realError.pdf", width = 3.54 * 2, height = 5)
print(p)
dev.off()


# Evaluate phased genotype of offspring
dirs <- c("../noError_homoP2_F2/",
          "../noError_hetP2_F2/",
          "../noError_homoP8_RIL/")
evaldf <- NULL
for(dir_i in dirs){
  scenario <- switch (dir_i,
                      "../noError_homoP2_F2/" = "homoP2_F2",
                      "../noError_hetP2_F2/" = "hetP2_F2",
                      "../noError_homoP8_RIL/" = "homoP8_RIL"
  )
  gbsr_files <- list.files(dir_i, "*_gbsr.gds", full.names = T)
  gbsr_no_files <- list.files(dir_i, "*_gbsr_noOptim.gds", full.names = T)
  magic_files <- list.files(dir_i, "_ImputedGenotype.csv", full.names = T)
  pop_files <- list.files(dir_i, "SimPop", full.names = T)
  datasets <- sub(".*/", "", sub("_gbsr.gds", "", gbsr_files))
  for(i_data in datasets){
    gds <- loadGDS(grep(paste0(i_data, "_"), gbsr_files, value = T))
    parents <- grep("Founder", getScanID(gds), value = T)
    gds <- setParents(gds, parents = parents, flip = F, mono = F, bi = F)
    p_gbsr <- GBScleanR::getGenotype(gds, node = "parents")
    if(!grepl("het", scenario)){
      p_gbsr <- p_gbsr[c(T, F), ]
    }
    gbsr <- GBScleanR::getHaplotype(gds)
    for(i in 1:dim(gbsr)[3]){
      for(j in 1:dim(gbsr)[2]){
        gbsr[, j, i] <- p_gbsr[gbsr[,j, i], j]
      }
    }

    gds_no <- loadGDS(grep(paste0(i_data, "_"), gbsr_no_files, value = T))
    parents <- grep("Founder", getScanID(gds_no), value = T)
    gds_no <- setParents(gds_no, parents = parents, flip = F, mono = F, bi = F)
    p_gbsr_no <- GBScleanR::getGenotype(gds_no, node = "parents")
    gbsr_no <- GBScleanR::getHaplotype(gds_no)
    if(!grepl("het", scenario)){
      p_gbsr_no <- p_gbsr_no[c(T, F), ]
    }
    for(i in 1:dim(gbsr_no)[3]){
      for(j in 1:dim(gbsr_no)[2]){
        gbsr_no[, j, i] <- p_gbsr_no[gbsr_no[,j, i], j]
      }
    }

    magic <- read.table(grep(paste0(i_data, "_"), magic_files, value = T),
                        sep = ",", skip = 4, row.names = 1)
    magic <- magic[!grepl("Founder", rownames(magic)), ]
    magic <- plyr::alply(magic, 1, function(x){
      x <- do.call("rbind", strsplit(as.character(x), ""))
      x_dim <- dim(x)
      x <- matrix(as.numeric(x), x_dim[1], x_dim[2])
      na_x <- is.na(rowSums(x))
      x[na_x, ] <- NA
      return(t(x-1))
    })

    load(grep(paste0(i_data, "_"), pop_files, value = T))
    true <- SimPop::getGenotype(pop, gen = "last", fam = "all", sib = "all")

    for(i in 1:dim(gbsr)[3]){
      check1 <- sum(gbsr[1,,i] == true[[i]][1, ], na.rm = T) + sum(gbsr[2,,i] == true[[i]][2, ], na.rm = T)
      check2 <- sum(gbsr[1,,i] == true[[i]][2, ], na.rm = T) + sum(gbsr[2,,i] == true[[i]][1, ], na.rm = T)
      if(check2 > check1){
        tmp <- gbsr[1,,i]
        gbsr[1,,i] <- gbsr[2,,i]
        gbsr[2,,i] <- tmp
      }

      check1 <- sum(gbsr_no[1,,i] == true[[i]][1, ], na.rm = T) + sum(gbsr_no[2,,i] == true[[i]][2, ], na.rm = T)
      check2 <- sum(gbsr_no[1,,i] == true[[i]][2, ], na.rm = T) + sum(gbsr_no[2,,i] == true[[i]][1, ], na.rm = T)
      if(check2 > check1){
        tmp <- gbsr_no[1,,i]
        gbsr_no[1,,i] <- gbsr_no[2,,i]
        gbsr_no[2,,i] <- tmp
      }

      check1 <- sum(magic[[i]][1, ] == true[[i]][1, ], na.rm = T) + sum(magic[[i]][2, ] == true[[i]][2, ], na.rm = T)
      check2 <- sum(magic[[i]][1, ] == true[[i]][2, ], na.rm = T) + sum(magic[[i]][2, ] == true[[i]][1, ], na.rm = T)
      if(check2 > check1){
        tmp <- magic[[i]][1, ]
        magic[[i]][1, ] <- magic[[i]][2, ]
        magic[[i]][2, ] <- tmp
      }
    }
    gbsr <- apply(gbsr, 2, function(x)return(x))
    gbsr_no <- apply(gbsr_no, 2, function(x)return(x))
    true <- do.call("rbind", true)
    magic <- do.call("rbind", magic)

    eval_gbsr <- compareGeno(gbsr, true)$ind
    eval_gbsr <- subset(eval_gbsr, select = c("correct", "miscall"))
    eval_gbsr_no <- compareGeno(gbsr_no, true)$ind
    eval_gbsr_no <- subset(eval_gbsr_no, select = c("correct", "miscall"))
    eval_magic <- compareGeno(magic, true)$ind
    eval_magic <- subset(eval_magic, select = c("correct", "miscall"))
    gbsr_df <- cbind(scenario = scenario,
                     dataset = i_data,
                     tool = "gbsr",
                     data.frame(t(apply(eval_gbsr, 2, mean))),
                     data.frame(t(apply(eval_gbsr, 2, sd))))
    names(gbsr_df)[6:7] <- c("correct_sd", "miscall_sd")
    gbsr_no_df <- cbind(scenario = scenario,
                        dataset = i_data,
                        tool = "gbsr_noOptim",
                        data.frame(t(apply(eval_gbsr_no, 2, mean))),
                        data.frame(t(apply(eval_gbsr_no, 2, sd))))
    names(gbsr_no_df)[6:7] <- c("correct_sd", "miscall_sd")
    magic_df <- cbind(scenario = scenario,
                      dataset = i_data,
                      tool = "magic",
                      data.frame(t(apply(eval_magic, 2, mean))),
                      data.frame(t(apply(eval_magic, 2, sd))))
    names(magic_df)[6:7] <- c("correct_sd", "miscall_sd")
    evaldf <- rbind(evaldf,
                    gbsr_df,
                    gbsr_no_df,
                    magic_df)
  }
}

write.csv(evaldf, file = "../Eval_phasedgeno_noError.csv")

evaldf <- read.csv(file = "../Eval_phasedgeno_noError.csv")
library(ggplot2)
evaldf$mar <- sub("_.*", "", sub(".*mar", "", evaldf$dataset))
evaldf$oread <- sub(".*_oread", "", evaldf$dataset)
evaldf$fread <- sub(".*fread", "",sub("_oread.*", "", evaldf$dataset))
evaldf$ind <- sub("_.*", "", sub("ind", "", evaldf$dataset))
evaldf$indmar <- paste0("Ind = ", evaldf$ind, "\nMar = ", evaldf$mar)
evaldf$scenario <- factor(evaldf$scenario, levels = c("homoP2_F2", "hetP2_F2", "homoP8_RIL"))
evaldf <- evaldf[order(as.numeric(evaldf$oread)), ]
evaldf$oread <- factor(evaldf$oread, levels = c("0.1", "0.25", "0.5", "0.75", "1", "2", "3"))
evaldf$fread <- factor(evaldf$fread, levels = c("3", "1"))
evaldf <- evaldf[order(as.numeric(evaldf$ind), as.numeric(evaldf$mar)), ]
evaldf$indmar <- factor(evaldf$indmar, levels = unique(evaldf$indmar))
evaldf$tool <- factor(evaldf$tool, levels = c("gbsr_noOptim", "magic", "gbsr"))
evaldf <- evaldf[order(as.numeric(evaldf$tool)), ]
p <- ggplot(evaldf) +
  geom_path(aes(x = oread, y = correct, color = tool, group = paste(tool, fread), linetype = fread)) +
  geom_point(aes(x = oread, y = correct, color = tool)) +
  facet_grid(facets = indmar ~ scenario, scale = "free") +
  scale_color_manual(name = "Tool",
                     breaks = c("gbsr", "gbsr_noOptim", "magic"),
                     values = c("blue", "green", "magenta"),
                     labels = c("GBSR", "GBSR_noOpt", "Magic")) +
  scale_linetype_manual(name = "Founder read depth",
                        breaks = c("3", "1"),
                        values = c(1, 2)) +
  ylab("Correct call rate") +
  xlab("Read depth") +
  theme(legend.position = "top",
        legend.box = "vertical",
        legend.spacing.y = unit(0, "mm"),
        legend.box.margin = margin(r = 40, t = 0, b = 0),
        legend.box.just = "left",
        legend.spacing.x = unit(1, "mm"),
        axis.text = element_text(size = 8),
        strip.text = element_text(size = 10),
        axis.title = element_text(size = 10)) +
  guides(shape = guide_legend(nrow = 1, byrow = T))
pdf("../Figures/Figure_PhasedGenotype_noError.pdf", width = 3.54 * 2, height = 10)
print(p)
dev.off()

dirs <- c("../realError_homoP2_F2/",
          "../realError_hetP2_F2/",
          "../realError_homoP8_RIL/")
evaldf <- NULL
for(dir_i in dirs){
  scenario <- switch (dir_i,
                      "../realError_homoP2_F2/" = "homoP2_F2",
                      "../realError_hetP2_F2/" = "hetP2_F2",
                      "../realError_homoP8_RIL/" = "homoP8_RIL"
  )
  gbsr_files <- list.files(dir_i, "*_gbsr.gds", full.names = T)
  gbsr_no_files <- list.files(dir_i, "*_gbsr_noOptim.gds", full.names = T)
  magic_files <- list.files(dir_i, "_ImputedGenotype.csv", full.names = T)
  pop_files <- list.files(dir_i, "SimPop", full.names = T)
  datasets <- sub(".*/", "", sub("_gbsr.gds", "", gbsr_files))
  for(i_data in datasets){
    gds <- loadGDS(grep(paste0(i_data, "_"), gbsr_files, value = T))
    parents <- grep("Founder", getScanID(gds), value = T)
    gds <- setParents(gds, parents = parents, flip = F, mono = F, bi = F)
    p_gbsr <- GBScleanR::getGenotype(gds, node = "parents")
    if(!grepl("het", scenario)){
      p_gbsr <- p_gbsr[c(T, F), ]
    }
    gbsr <- GBScleanR::getHaplotype(gds)
    for(i in 1:dim(gbsr)[3]){
      for(j in 1:dim(gbsr)[2]){
        gbsr[, j, i] <- p_gbsr[gbsr[,j, i], j]
      }
    }

    gds_no <- loadGDS(grep(paste0(i_data, "_"), gbsr_no_files, value = T))
    parents <- grep("Founder", getScanID(gds_no), value = T)
    gds_no <- setParents(gds_no, parents = parents, flip = F, mono = F, bi = F)
    p_gbsr_no <- GBScleanR::getGenotype(gds_no, node = "parents")
    if(!grepl("het", scenario)){
      p_gbsr_no <- p_gbsr_no[c(T, F), ]
    }
    gbsr_no <- GBScleanR::getHaplotype(gds_no)
    for(i in 1:dim(gbsr_no)[3]){
      for(j in 1:dim(gbsr_no)[2]){
        gbsr_no[, j, i] <- p_gbsr_no[gbsr_no[,j, i], j]
      }
    }

    magic <- read.table(grep(paste0(i_data, "_"), magic_files, value = T),
                        sep = ",", skip = 4, row.names = 1, stringsAsFactors = F)
    magic <- magic[!grepl("Founder", rownames(magic)), ]
    magic <- plyr::alply(magic, 1, function(x){
        x <- do.call("rbind", strsplit(as.character(x), ""))
        x_dim <- dim(x)
        x <- matrix(as.numeric(x), x_dim[1], x_dim[2])
        na_x <- is.na(rowSums(x))
        x[na_x, ] <- NA
        return(t(x-1))
      })

    load(grep(paste0(i_data, "_"), pop_files, value = T))
    true <- SimPop::getGenotype(pop, gen = "last", fam = "all", sib = "all")

    for(i in 1:dim(gbsr)[3]){
      check1 <- sum(gbsr[1,,i] == true[[i]][1, ], na.rm = T) + sum(gbsr[2,,i] == true[[i]][2, ], na.rm = T)
      check2 <- sum(gbsr[1,,i] == true[[i]][2, ], na.rm = T) + sum(gbsr[2,,i] == true[[i]][1, ], na.rm = T)
      if(check2 > check1){
        tmp <- gbsr[1,,i]
        gbsr[1,,i] <- gbsr[2,,i]
        gbsr[2,,i] <- tmp
      }

      check1 <- sum(gbsr_no[1,,i] == true[[i]][1, ], na.rm = T) + sum(gbsr_no[2,,i] == true[[i]][2, ], na.rm = T)
      check2 <- sum(gbsr_no[1,,i] == true[[i]][2, ], na.rm = T) + sum(gbsr_no[2,,i] == true[[i]][1, ], na.rm = T)
      if(check2 > check1){
        tmp <- gbsr_no[1,,i]
        gbsr_no[1,,i] <- gbsr_no[2,,i]
        gbsr_no[2,,i] <- tmp
      }

      check1 <- sum(magic[[i]][1, ] == true[[i]][1, ], na.rm = T) + sum(magic[[i]][2, ] == true[[i]][2, ], na.rm = T)
      check2 <- sum(magic[[i]][1, ] == true[[i]][2, ], na.rm = T) + sum(magic[[i]][2, ] == true[[i]][1, ], na.rm = T)
      if(check2 > check1){
        tmp <- magic[[i]][1, ]
        magic[[i]][1, ] <- magic[[i]][2, ]
        magic[[i]][2, ] <- tmp
      }
    }
    gbsr <- apply(gbsr, 2, function(x)return(x))
    gbsr_no <- apply(gbsr_no, 2, function(x)return(x))
    true <- do.call("rbind", true)
    magic <- do.call("rbind", magic)

    eval_gbsr <- compareGeno(gbsr, true)$ind
    eval_gbsr <- subset(eval_gbsr, select = c("correct", "miscall"))
    eval_gbsr_no <- compareGeno(gbsr_no, true)$ind
    eval_gbsr_no <- subset(eval_gbsr_no, select = c("correct", "miscall"))
    eval_magic <- compareGeno(magic, true)$ind
    eval_magic <- subset(eval_magic, select = c("correct", "miscall"))
    gbsr_df <- cbind(scenario = scenario,
                     dataset = i_data,
                     tool = "gbsr",
                     data.frame(t(apply(eval_gbsr, 2, mean))),
                     data.frame(t(apply(eval_gbsr, 2, sd))))
    names(gbsr_df)[6:7] <- c("correct_sd", "miscall_sd")
    gbsr_no_df <- cbind(scenario = scenario,
                        dataset = i_data,
                        tool = "gbsr_noOptim",
                        data.frame(t(apply(eval_gbsr_no, 2, mean))),
                        data.frame(t(apply(eval_gbsr_no, 2, sd))))
    names(gbsr_no_df)[6:7] <- c("correct_sd", "miscall_sd")
    magic_df <- cbind(scenario = scenario,
                      dataset = i_data,
                      tool = "magic",
                      data.frame(t(apply(eval_magic, 2, mean))),
                      data.frame(t(apply(eval_magic, 2, sd))))
    names(magic_df)[6:7] <- c("correct_sd", "miscall_sd")
    evaldf <- rbind(evaldf,
                    gbsr_df,
                    gbsr_no_df,
                    magic_df)
  }
}

write.csv(evaldf, file = "../Eval_phasedgeno_realError.csv")

evaldf <- read.csv(file = "../Eval_phasedgeno_realError.csv")
library(ggplot2)
evaldf$oread <- sub(".*_oread", "", evaldf$dataset)
evaldf$fread <- sub(".*fread", "",sub("_oread.*", "", evaldf$dataset))
evaldf$ind <- sub("_.*", "", sub("ind", "", evaldf$dataset))
evaldf$scenario <- factor(evaldf$scenario, levels = c("homoP2_F2", "hetP2_F2", "homoP8_RIL"))
evaldf <- evaldf[order(as.numeric(evaldf$oread)), ]
evaldf$oread <- factor(evaldf$oread, levels = c("0.1", "0.25", "0.5", "0.75", "1", "2", "3"))
evaldf$fread <- factor(evaldf$fread, levels = c("3", "1"))
evaldf$tool <- factor(evaldf$tool, levels = c("gbsr_noOptim", "magic", "gbsr"))
evaldf <- evaldf[order(as.numeric(evaldf$tool)), ]
p <- ggplot(evaldf) +
  geom_path(aes(x = oread, y = correct, color = tool, group = paste(tool, fread), linetype = fread)) +
  geom_point(aes(x = oread, y = correct, color = tool)) +
  facet_grid(facets = ind ~ scenario, scale = "free") +
  scale_color_manual(name = "Tool",
                     breaks = c("gbsr", "gbsr_noOptim", "magic"),
                     values = c("blue", "green", "magenta"),
                     labels = c("GBSR", "GBSR_noOpt", "Magic")) +
  scale_linetype_manual(name = "Founder read depth",
                        breaks = c("3", "1"),
                        values = c(1, 2)) +
  ylab("Correct call rate") +
  xlab("Read depth") +
  theme(legend.position = "top",
        legend.box = "vertical",
        legend.spacing.y = unit(0, "mm"),
        legend.box.margin = margin(r = 40, t = 0, b = 0),
        legend.box.just = "left",
        legend.spacing.x = unit(1, "mm"),
        axis.text = element_text(size = 8),
        strip.text = element_text(size = 10),
        axis.title = element_text(size = 10)) +
  guides(shape = guide_legend(nrow = 1, byrow = T))
pdf("../Figures/Figure_PhasedGenotype_realError.pdf", width = 3.54 * 2, height = 7)
print(p)
dev.off()
