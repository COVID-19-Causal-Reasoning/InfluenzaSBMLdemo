## ----setup, include=FALSE, warning=FALSE---------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## ------------------------------------------------------------------------------------------------------------------
# install_github("saezlab/progeny")
# install_github("saezlab/dorothea")
#install_github('COVID-19-Causal-Reasoning/CARNIVAL', ref="add_gurobi", force=TRUE)
library(argparser)
library(readr)
# 
p <- arg_parser("Run CARNIVAL on DOROTHEA and PROGENy output")
p <- add_argument(p, "pathway_activities_file",  help="Pathway activities  (from Progeny)")
p <- add_argument(p, "tf_activities_file", help="Transcription Factor Activities (From Dorothea)")
p <- add_argument(p, "--perturbation_file", help="Perturbation File", default="NoInput")
p <- add_argument(p, "--variable", help="Intervention (e.g. VN1203)")
p <- add_argument(p, "--baseline", help="The thing to be compared against. E.g. Mock")
p <- add_argument(p, "--constant", help="The things that are held constant for both variable and baseline. For example timepoint (7h) or strain (VN1203")
p <- add_argument(p, "--outdir", help="output directory")
p <- add_argument(p, "--timelimit", help="Number of seconds before timeout. Default 1200", default=1200)
p <- add_argument(p, "--mipgap", help="MIP Gap. default: 3%", default=0.03)
p <- add_argument(p, "--solver", "cplex gurobi lpsolve or cbc", default="cplex")
p <- add_argument(p, "--solver_path", "location of solver binary", default="/Applications/CPLEX_Studio_Community129/cplex/bin/x86-64_osx/cplex")
argv <- parse_args(p) #, c( "IntermediateFiles/VN1203vsMock_for_24h_pathways_activity_score_inputCarnival.csv",
                      #  "IntermediateFiles/VN1203vsMock_for_24h_tf_activities_stat_inputCarnival.csv",
                      #  "--variable", "VN1203", 
                      #  "--baseline", "Mock", 
                      #  "--constant", "24hr",
                      #  "--timelimit", "150",
                      #  "--mipgap", "0.40",
                      #  "--outdir", "ResultsCARNIVAL_VN1203vsMock_for_24hr"))



Sys.setenv(CPLEX_STUDIO_KEY='api_cos_0c913dc4-503c-4334-9b4e-ae4c7d13158b')
Sys.setenv(CPLEX_STUDIO_DIR1210 = '/Applications/CPLEX_Studio_Community129' )
carnival_dir <-"~/Projects/CausalInference/Pathogenesis/OmniPath/CARNIVAL/R"
for(filename in list.files(carnival_dir)){
  source(paste( carnival_dir, filename, sep="/"))}

#library(CARNIVAL)


## ------------------------------------------------------------------------------------------------------------------
library(rlist)
library(devtools)
library(tidyverse)
library(OmnipathR)
library(dplyr)
library(tibble)
library(openxlsx)

## We also define a function to format the CARNIVAL output to cytoscape
OutputCyto <- function(CarnivalResults, outputFile) {
    CarnivalNetwork <- 
        as.data.frame(CarnivalResults$weightedSIF, stringsAsFactors = FALSE) 
    
    CarnivalNetworkNodes <- 
        unique(c(CarnivalNetwork$Node1,CarnivalNetwork$Node2))
    
    CarnivalAttributes <- CarnivalResults$nodesAttributes %>% 
        as.data.frame() %>%
        dplyr::filter(Node %in% CarnivalNetworkNodes) %>%
        dplyr::mutate(NodeType = as.character(NodeType)) %>%
        dplyr::mutate(NodeType=if_else(NodeType =="", "I", NodeType))
            
    nameOutputNetwork <- paste0(outputFile, "Network.sif")
    nameOutputAttributes <-  paste0(outputFile, "Attributes.txt")    
    
    write.table(CarnivalNetwork[,c(1,2,3)], file = nameOutputNetwork,
        quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ")
    
    write.table(CarnivalAttributes, file = nameOutputAttributes,
        quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
}


## ---- message=FALSE------------------------------------------------------------------------------------------------
## the OmniPath PPI interaction network
ia_omnipath <- import_Omnipath_Interactions() %>% as_tibble()

## We also download the other datasets containing interactions
ia_ligrec <- import_LigrecExtra_Interactions() %>% as_tibble()
ia_pwextra <- import_PathwayExtra_Interactions() %>% as_tibble()
ia_kinaseextra <- import_KinaseExtra_Interactions() %>% as_tibble()

## We bind the datasets
interactions <- as_tibble(
    bind_rows(
        ia_omnipath %>% mutate(type = 'ppi'),
        ia_pwextra %>% mutate(type = 'ppi'),
        ia_kinaseextra %>% mutate(type = 'ppi'),
        ia_ligrec %>% mutate(type = 'ppi')))

## I am going to keep only directed interactions (consensus_direction) and 
## signed interactions (consensus_stimulation/consensus_inhibition)
## We transform to the format needed by CARNIVAL. We just keep signed and 
## directed interactions 
SignedDirectedInteractions <- 
    dplyr::filter(interactions, consensus_direction==1) %>%
    filter(consensus_stimulation == 1 | consensus_inhibition == 1)

NetworkCarnival_df <- bind_rows(
  (SignedDirectedInteractions %>%
  filter(consensus_stimulation == 1 & consensus_inhibition == 0) %>%
  transmute(source_genesymbol, interaction = 1, target_genesymbol)),   
  (SignedDirectedInteractions %>%
     filter(consensus_stimulation == 0 & consensus_inhibition == 1) %>%
     transmute(source_genesymbol, interaction = -1, target_genesymbol))) %>%  
  distinct() 

## We transform the network to an igraph object to simplify
NetworkCarnival_igraph <- 
    graph_from_data_frame(NetworkCarnival_df[c(1,3,2)], directed = TRUE) %>% 
    igraph::simplify(remove.multiple = TRUE, remove.loops = TRUE, 
        edge.attr.comb = "first")

## We transform back to the format required by CARNIVAL
NetworkCarnival_df <- igraph::as_data_frame(NetworkCarnival_igraph) %>%
    dplyr::select(from, interaction, to) %>%  
    distinct() 

## We have to be careful with the gene names with a "-". CPLEX gets crazy. 
NetworkCarnival_df$from <- gsub("-","_", NetworkCarnival_df$from)
NetworkCarnival_df$to <- gsub("-","_", NetworkCarnival_df$to)

AllNodesNetwork <- unique(c(NetworkCarnival_df$from, NetworkCarnival_df$to))


## ------------------------------------------------------------------------------------------------------------------
pathways_activity_score_inputCarnival <- as.matrix(read.csv(argv$pathway_activities_file, row.names=1))
tf_activities_stat <- as.matrix(read.csv(argv$tf_activities_file, row.names=1))


## ---- eval=TRUE, echo=TRUE-----------------------------------------------------------------------------------------
### VN1203 vs Mock
tf_activities_stat_top50 <- tf_activities_stat %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "GeneID") %>%
  dplyr::filter(GeneID %in% AllNodesNetwork) %>%
  dplyr::arrange(desc(abs(stat))) %>%
  dplyr::top_n(50, wt = abs(stat)) %>%
  column_to_rownames(var = "GeneID") %>%
  t()
write.csv(colnames(tf_activities_stat_top50), 
   paste0(argv$outdir, "/Top50_tf_activities_", argv$variable, "vs", argv$baseline, "_for_", argv$constant))




## ------------------------------------------------------------------------------------------------------------------
if ( argv$perturbation_file == "NoInput"){
  HostVirus_perturbation = NULL
} else {
 
  perturbation <- read.csv(argv$perturbation_file)

  NetworkCarnivalHostVirus_df <- 
    bind_rows(NetworkCarnival_df, perturbation) %>% 
    dplyr::distinct()

## We have to be careful with the gene names with a "-". CPLEX gets crazy. 
  NetworkCarnivalHostVirus_df$from <- 
    gsub("-","_", NetworkCarnivalHostVirus_df$from)
  NetworkCarnivalHostVirus_df$to <- 
   gsub("-","_", NetworkCarnivalHostVirus_df$to)

  Virusproteins <- unique(perturbation$from)
  Virusproteins <-  gsub("-","_", Virusproteins)
  HostVirus_perturbation <- 
   data.frame(Virusproteins = rep(1,length(Virusproteins)))
  rownames(HostVirus_perturbation) <- Virusproteins
  HostVirus_perturbation <- as.data.frame(t(as.matrix(HostVirus_perturbation)))
} 
  




## ------------------------------------------------------------------------------------------------------------------
if (is.null(perturbation_file)) {
  saveRDS(CarnivalResults, file = paste0(argv$outdir, "/", argv$variable, "vs", argv$baseline, "_for_", argv$constant, "_noPerturbation.rds"))
  OutputCyto(CarnivalResults, 
  outputFile= paste0(argv$outdir, "/", argv$variable, "vs", argv$baseline, "_for_", argv$constant, "_noPerturbation"))
} else {
   saveRDS(CarnivalResults, file = paste0(argv$outdir, "/", argv$variable, "vs", argv$baseline, "_for_", argv$constant, "_HostVirusPerturbation.rds"))
   OutputCyto(CarnivalResults, 
  outputFile= paste0(argv$outdir, "/", argv$variable, "vs", argv$baseline, "_for_", argv$constant, "_HostVirusPerturbation"))
}


