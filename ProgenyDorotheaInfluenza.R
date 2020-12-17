## ----setup, include=FALSE, warning=FALSE---------------------------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## ---- message=FALSE------------------------------------------------------------------------------------------------------------------------------------------
library(argparser)
library(readr)
# 
p <- arg_parser("Run DOROTHEA and PROGENy on differential expression data")
p <- add_argument(p, "dds_file",  help="Differential expression wald statistic  (VN1203 vs Mock holding timepoint 7h constant) or  (12hr vs 7hr holding VN1203 constant")
p <- add_argument(p, "norm_counts_file", help="Normalized counts")
p <- add_argument(p, "--variable", help="Intervention (e.g. VN1203)")
p <- add_argument(p, "--baseline", help="The thing to be compared against. E.g. Mock")
p <- add_argument(p, "--constant", help="The things that are held constant for both variable and baseline. For example timepoint (7h) or strain (VN1203")
p <- add_argument(p, "--organism", help="Which organism: Human or Mouse", default="Human")
p <- add_argument(p, "--outdir", help="output directory")
argv <- parse_args(p) #c("Dorothea_ICL103_Proteins_VN1203_vs_NS1_24hr.csv", 
                      #  "Progeny_ICL103_Proteins_VN1203_vs_Mock_24hr.csv",
                      #  "--variable", "VN1203", 
                      #  "--baseline", "NS1", 
                      #  "--constant", "24hr",
                      #  "--outdir", "IntermediateFiles"))

library(devtools)
# install_github("saezlab/progeny")
# install_github("saezlab/dorothea")
library(progeny)
library(dorothea)
library(tibble)
library(tidyr)
library(dplyr)
library(ggplot2)
library(pheatmap)


## ------------------------------------------------------------------------------------------------------------------------------------------------------------
## Differential expression table and Normalized counts
dds <- as.matrix(read.csv(argv$dds_file))
norm_counts <-  as.matrix(read.csv(argv$norm_counts_file))


## ---- message=FALSE------------------------------------------------------------------------------------------------------------------------------------------
dds_df <- as.data.frame(dds) %>% 
    rownames_to_column(var = "GeneID") %>% 
    dplyr::select(GeneID, stat) %>% 
    dplyr::filter(!is.na(stat)) %>% 
    column_to_rownames(var = "GeneID") 


## I also need to run progeny in such a way to have values between 1 and -1 to
## use as CARNIVAL input
pathways_activity_score_inputCarnival <- 
  t(progeny(as.matrix(dds_df), 
    scale=TRUE, organism=argv$organism, top = 100, perm = 10000, z_scores = FALSE))
colnames(pathways_activity_score_inputCarnival) <- "Activity"


## ------------------------------------------------------------------------------------------------------------------------------------------------------------
## We load Dorothea Regulons
data(dorothea_hs, package = "dorothea")
regulons <- dorothea_hs %>%
  dplyr::filter(confidence %in% c("A", "B","C"))


## ------------------------------------------------------------------------------------------------------------------------------------------------------------
dds_stat <-  as.data.frame(dds) %>% 
    rownames_to_column(var = "GeneID") %>% 
    dplyr::select(GeneID, stat) %>% 
    dplyr::filter(!is.na(stat)) %>% 
    column_to_rownames(var = "GeneID") %>%
    as.matrix()

tf_activities_stat <- 
    dorothea::run_viper(as.matrix(dds_stat), regulons,
    options =  list(minsize = 5, eset.filter = FALSE, 
    cores = 1, verbose = FALSE, nes = TRUE))


## ------------------------------------------------------------------------------------------------------------------------------------------------------------
tf_activities_top25 <- tf_activities_stat %>%
    as.data.frame() %>% 
    rownames_to_column(var = "GeneID") %>%
    dplyr::rename(NES = "stat") %>%
    dplyr::top_n(25, wt = abs(NES)) %>%
    dplyr::arrange(NES) %>% 
    dplyr::mutate(GeneID = factor(GeneID))


## ------------------------------------------------------------------------------------------------------------------------------------------------------------
tf_activities_counts <- 
    dorothea::run_viper(norm_counts, regulons,
    options =  list(minsize = 5, eset.filter = FALSE, 
    cores = 1, verbose = FALSE, method = c("scale")))

tf_activities_counts_filter <- tf_activities_counts %>% 
    as.data.frame() %>% 
    rownames_to_column(var = "GeneID") %>%
    dplyr::filter(GeneID %in% tf_activities_top25$GeneID) %>%
    column_to_rownames(var = "GeneID") %>%
    as.matrix()

tf_activities <- as.vector(tf_activities_counts_filter)

paletteLength <- 100
myColor <- 
    colorRampPalette(c("darkblue", "whitesmoke","indianred"))(paletteLength)

dorotheaBreaks <- c(seq(min(tf_activities), 0, 
    length.out=ceiling(paletteLength/2) + 1),
    seq(max(tf_activities)/paletteLength, 
    max(tf_activities), 
    length.out=floor(paletteLength/2)))


## ------------------------------------------------------------------------------------------------------------------------------------------------------------
write.csv(as.data.frame(pathways_activity_score_inputCarnival), 
     paste0(argv$outdir,"/", argv$variable, "vs", argv$baseline, "_for_", argv$constant, "_pathways_activity_score_inputCarnival.csv"))
write.csv(as.data.frame(regulons), 
    paste0(argv$outdir, "/dorothea_regulons.csv"))
write.csv(as.data.frame(tf_activities_stat), 
    paste0(argv$outdir,"/", argv$variable, "vs", argv$baseline, "_for_", argv$constant, "_tf_activities_stat_inputCarnival.csv"))

