setwd("~/Box/RA Cliff Metaclustering")

library(flowStats)
library(Rtsne)

#### Load Data ####
#list all files in the FCS data folder
all_files <- list.files(pattern = ".fcs", full.names = T, recursive = T)

#Subset to total CD4 files
cd4_pre_files <- all_files[str_detect(all_files, regex("Pre_Tube", ignore_case = T))]

#Subset to tmr+ files
flu_tmr_files <- all_files[str_detect(all_files, "TMR-INFLUENZA")]
agg_tmr_files <- all_files[str_detect(all_files, "TMR-AGGRECAN")]
eno_tmr_files <- all_files[str_detect(all_files, "TMR-ENOLASE")]
cilp_tmr_files <- all_files[str_detect(all_files, "TMR-CILP")]
vf_tmr_files <- all_files[str_detect(all_files, "TMR-VIMENTIN")]

build_fcsList <- function(file_list) {
  tmp_list <- list()
  for (fcs in file_list) {
    fcs_name <- paste(gsub(".*Cohort (1[a-b]).*", "\\1", fcs),
                      str_extract(fcs, "[\\d]{3}_[\\d]{3}") %>% str_sub(1, 3),
                      sep = "_")
    tmp_list[[fcs_name]] <- flowCore::read.FCS(fcs, transformation = F)
  }
  return(tmp_list)
}

cd4_pre_files <- build_fcsList(cd4_pre_files)
flu_tmr_files <- build_fcsList(flu_tmr_files)
agg_tmr_files <- build_fcsList(agg_tmr_files)
eno_tmr_files <- build_fcsList(eno_tmr_files)
cilp_tmr_files <- build_fcsList(cilp_tmr_files)
vf_tmr_files <- build_fcsList(vf_tmr_files)

#Process data
key = read.csv("MarkerKey.csv", header = T)

process_data <- function(fcs){
  #counter <<- counter + 1 # This code used to find problematic FCS files where all events are boundary events
  #print(counter)
  
  # Tidy marker names 
  marker_key <- pData(parameters(fcs))
  pData(parameters(fcs))$desc <- key$desc[match(marker_key$name, key$name)]
  
  # This changes parameters(fcs)$name, featureNames(fcs), and colnames(fcs) - aka events colnames - all in one fell swoop.
  colnames(fcs) = make.names(pData(parameters(fcs))$desc)
  
  # Remove markers with multiple markers (i.e. all the NAs)
  fcs = fcs[,!(duplicated(colnames(fcs)) | duplicated(colnames(fcs), fromLast = T))]
  fcs = fcs[, order(colnames(fcs))]
  
  # Remove boundary events (those outside of FCS's range)
  selected_markers = colnames(fcs)[!str_detect(colnames(fcs),"Tmr")]
  fcs =  Subset(fcs, filter(fcs, boundaryFilter(x = selected_markers, tolerance = 0, side = "both")))
}

cd4_pre_files <- lapply(cd4_pre_files, process_data)
flu_tmr_files <- lapply(flu_tmr_files, process_data)
agg_tmr_files <- lapply(agg_tmr_files, process_data)
eno_tmr_files <- lapply(eno_tmr_files, process_data)
#counter = 0
cilp_tmr_files = cilp_tmr_files[-c(27:28)] # All events are boundary events for these individuals
cilp_tmr_files <- lapply(cilp_tmr_files, process_data) 
#counter = 0
vf_tmr_files = vf_tmr_files[-c(23, 27, 73, 80, 91)] # All events are boundary events for these individuals
vf_tmr_files <- lapply(vf_tmr_files, process_data) 

cd4_files <- cd4_pre_files

# Merge total and tmr+ expression data
all_merged = cd4_files
for(subj in names(cd4_files)){
  exprs(all_merged[[subj]]) = rbind(exprs(cd4_pre_files[[subj]]), exprs(flu_tmr_files[[subj]]), 
                                    exprs(agg_tmr_files[[subj]]), exprs(eno_tmr_files[[subj]]), 
                                    if(is.null(cilp_tmr_files[[subj]]) == F) exprs(cilp_tmr_files[[subj]]), 
                                    if(is.null(vf_tmr_files[[subj]]) == F) exprs(vf_tmr_files[[subj]]))
}

# Transform the data
asinh_tfm_data <- function(fcs){
  # Arcsinh transform remaining columns
  tl <- transformList(colnames(fcs), arcsinhTransform(a=0, b=1/150), transformationId="asinh")
  fcs = transform(fcs, tl) 
}

all_tfm <- lapply(all_merged, asinh_tfm_data)
fcsNames = names(all_tfm)

tfm_fs <- as(all_tfm, "flowSet")
assay_markers <- c("Agg.Tmr", "CD45ra", "CD38", "Flu.Tmr", "CCR4", "CCR6", "Eno.Tmr", "Cilp.Tmr", "VF.Tmr", "CXCR3", "CCR7", "CD4")
pheno_markers <- c("CD45ra", "CD38", "CCR4", "CCR6", "CXCR3", "CCR7", "CD4")
clust_markers <- c("CD45ra", "CD38", "CCR4", "CCR6", "CXCR3", "CCR7")

#### Extract expression data and label Tmr+ events ####
mergedExpr = setNames(data.frame(matrix(ncol = ncol(exprs(tfm_fs[[1]]))+2, nrow = 0)), c(colnames(tfm_fs), "samp", "tmr_pos"))
for(n in fcsNames){
  tmp.expr = as.data.frame(exprs(tfm_fs[[n]]))[,assay_markers]
  tmp.expr$samp = as.character(n)
  tmp.expr$tmr_pos = c(rep("none", nrow(exprs(cd4_files[[n]]))), 
                       rep("flu", nrow(exprs(flu_tmr_files[[n]]))),
                       rep("agg", nrow(exprs(agg_tmr_files[[n]]))),
                       rep("eno", nrow(exprs(eno_tmr_files[[n]]))),
                       if(is.null(cilp_tmr_files[[n]]) == F) rep("cilp", nrow(exprs(cilp_tmr_files[[n]]))),
                       if(is.null(vf_tmr_files[[n]]) == F) rep("vf", nrow(exprs(vf_tmr_files[[n]]))))
  mergedExpr = rbind(mergedExpr, tmp.expr)
}


#### Cluster using Rphenograph with kd tree-type adaptation ####
find_neighbors <- function(data, k){
  nearest <- RANN::nn2(data, data, k, treetype = "kd", searchtype = "standard")
  return(nearest[[1]])
}

Rpheno <- function(data, k=30){
  if(is.data.frame(data))
    data <- as.matrix(data)
  
  if(!is.matrix(data))
    stop("Wrong input data, should be a data frame or matrix!")
  
  if(k<1){
    stop("k must be a positive integer!")
  }else if (k > nrow(data)-2){
    stop("k must be smaller than the total number of points!")
  }
  
  message("Run Rphenograph starts:","\n", 
          "  -Input data of ", nrow(data)," rows and ", ncol(data), " columns","\n",
          "  -k is set to ", k)
  
  cat("  Finding nearest neighbors...")
  t1 <- system.time(neighborMatrix <- find_neighbors(data, k=k+1)[,-1])
  cat("DONE ~",t1[3],"s\n", " Compute jaccard coefficient between nearest-neighbor sets...")
  t2 <- system.time(links <- Rphenograph:::jaccard_coeff(neighborMatrix))
  
  cat("DONE ~",t2[3],"s\n", " Build undirected graph from the weighted links...")
  links <- links[links[,1]>0, ]
  relations <- as.data.frame(links)
  colnames(relations)<- c("from","to","weight")
  t3 <- system.time(g <- igraph::graph.data.frame(relations, directed=FALSE))
  
  cat("DONE ~",t3[3],"s\n", " Run louvain clustering on the graph ...")
  t4 <- system.time(community <- igraph::cluster_louvain(g))
  cat("DONE ~",t4[3],"s\n")
  
  message("Run Rphenograph DONE, took a total of ", sum(c(t1[3],t2[3],t3[3],t4[3])), "s.")
  cat("  Return a community class\n  -Modularity value:", igraph::modularity(community),"\n")
  cat("  -Number of clusters:", length(unique(igraph::membership(community))))
  
  return(community)
}


PhenographClust_kd = function(fcs, clustering_markers) {
  exprs_mat = as.matrix(as.data.frame(exprs(fcs))[,clustering_markers])
  RPvect = as.numeric(igraph::membership(Rpheno(data = exprs_mat)))
  return(RPvect)
}

mergedExpr$RPclust = unlist(lapply(all_tfm, PhenographClust_kd, clust_markers))

n_pheno_clusts <- mergedExpr %>%
  group_by(samp) %>%
  summarize(k_clusters = max(RPclust))
# This generated 12-22 clusters per individual (mostly in the 14-20 range).

save(mergedExpr, all_tfm, tfm_fs, agg_tmr_files, cilp_tmr_files, eno_tmr_files, flu_tmr_files, vf_tmr_files, 
     cd4_files, pheno_markers, file = "Phenograph_clusters.RData")


#### Count CD8+ and Tmr+ cells in each cluster, run enrichment ####
## For RP clustering
# Collect mean and median expression data for each cluster in each individual
RP_median <- aggregate(. ~ RPclust + samp, data = mergedExpr[,c(1:13, 15)], median)
RP_median = rbind(RP_median, data.frame(RPclust = "Total_CD4", aggregate(. ~ samp, data = mergedExpr[,c(1:13)], median)))
RP_mean <- aggregate(. ~ RPclust + samp, data = mergedExpr[,c(1:13, 15)], mean)
RP_mean = rbind(RP_mean, data.frame(RPclust = "Total_CD4", aggregate(. ~ samp, data = mergedExpr[,c(1:13)], mean)))

# Count Tmr+ cells in each cluster for each individual
RPtmr_counting = mergedExpr[,c(13:15)] %>% group_by(samp, RPclust) %>% 
  summarise(clust_size = n(), Aggrecan = sum(tmr_pos=="agg"), Cilp = sum(tmr_pos=="cilp"), 
            Enolase = sum(tmr_pos=="eno"), Vimentin_Fibrinogen = sum(tmr_pos=="vf"),
            Influenza = sum(tmr_pos=="flu"), None = sum(tmr_pos=="none"))

aggreg_counts = aggregate(. ~ samp, data = RPtmr_counting[,c(1,3:9)], sum)
colnames(aggreg_counts)[2:8] = paste(colnames(aggreg_counts)[2:8], "tot", sep="_")
RPtmr_counting = merge(RPtmr_counting, aggreg_counts)

calc_enrichment = function(tmr_count, tmr_column) {
  expected = paste("E", tmr_column, sep="_")
  tot_tmr_perSamp = paste(tmr_column, "tot", sep="_")
  tmr_count[,expected] = tmr_count[,tot_tmr_perSamp] * tmr_count[,"clust_size"] / tmr_count[,"clust_size_tot"] %>% round(digits = 3)
  tmr_count[,paste("hg_pval", tmr_column, sep="_")] = phyper(tmr_count[,tmr_column]-1, tmr_count[,tot_tmr_perSamp], 
                                                             (tmr_count[,"clust_size_tot"] - tmr_count[,tot_tmr_perSamp]), 
                                                             tmr_count[,"clust_size"], lower.tail = FALSE) %>% signif(digits = 3)
  return(tmr_count)
}


RPtmr_counting = calc_enrichment(RPtmr_counting, "Aggrecan")
RPtmr_counting = calc_enrichment(RPtmr_counting, "Cilp")
RPtmr_counting = calc_enrichment(RPtmr_counting, "Enolase")
RPtmr_counting = calc_enrichment(RPtmr_counting, "Vimentin_Fibrinogen")
RPtmr_counting = calc_enrichment(RPtmr_counting, "Influenza")


RPtmr_counting = mutate(RPtmr_counting, Aggrecan_pct_Aggrecan = Aggrecan/Aggrecan_tot * 100,
                        Aggrecan_pct_tot_tmr = Aggrecan/(clust_size_tot - None_tot) * 100,
                        Cilp_pct_Cilp = Cilp/Cilp_tot * 100,
                        Cilp_pct_tot_tmr = Cilp/(clust_size_tot - None_tot) * 100,
                        Enolase_pct_Enolase = Enolase/Enolase_tot * 100,
                        Enolase_pct_tot_tmr = Enolase/(clust_size_tot - None_tot) * 100,
                        Vimentin_Fibrinogen_pct_Vimentin_Fibrinogen = Vimentin_Fibrinogen/Vimentin_Fibrinogen_tot * 100,
                        Vimentin_Fibrinogen_pct_tot_tmr = Vimentin_Fibrinogen/(clust_size_tot - None_tot) * 100,
                        Influenza_pct_Influenza = Influenza/Influenza_tot * 100,
                        Influenza_pct_tot_tmr = Influenza/(clust_size_tot - None_tot) * 100)

for_mario = subset(RPtmr_counting, hg_pval_Aggrecan < 0.05 | hg_pval_Cilp < 0.05 | hg_pval_Enolase < 0.05 |
                     hg_pval_Vimentin_Fibrinogen < 0.05 | hg_pval_Influenza < 0.05)
for_mario$RPclust = as.character(for_mario$RPclust)
for_mario = left_join(for_mario, RP_mean, by = c("samp", "RPclust"))

RA_only = subset(for_mario, hg_pval_Aggrecan < 0.05 | hg_pval_Cilp < 0.05 | hg_pval_Enolase < 0.05 |
                   hg_pval_Vimentin_Fibrinogen < 0.05)


#### tSNE Clustering ####
tsne_flowset <- function(fcs, selected_markers){
  tsne_data <- as.data.frame(exprs(fcs))
  
  tsne_out <- tsne_data[,colnames(tsne_data) %in% selected_markers] %>%
    as.matrix %>%
    Rtsne(check_duplicates = F)
  
  tsne_data$tsne_1 = tsne_out$Y[,1]
  tsne_data$tsne_2 = tsne_out$Y[,2]
  
  return(tsne_data)
}

tsne_results <- lapply(all_tfm, tsne_flowset, selected_markers = clust_markers)

save(tsne_results, clust_markers, all_tfm, tfm_fs, file = "individual_tsne_results.Rdata")

tsne_results <- do.call("rbind", tsne_results)
if(setequal(tsne_results[,colnames(tsne_results) %in% clust_markers], mergedExpr[,colnames(mergedExpr) %in% clust_markers])==T) {
  tsne_results <- cbind(tsne_results, mergedExpr$samp, mergedExpr$RPclust, mergedExpr$tmr_pos)
  colnames(tsne_results) = str_replace(colnames(tsne_results), "mergedExpr\\$", "")
  tsne_results$RPclust = as.factor(tsne_results$RPclust)
} else {
  print("The mergedExpr and tsne_results expression data do NOT match!!")
}


#### tSNE plots Tmr+ overlays ####
cb_pal <- colorblind_pal()(8)
vm_pal <- cb_pal[-1] # This is Hannah's favorite color scheme

#vm_pal <- brewer.pal(8, "Set1")
#inf_pal <- inferno(30)
#vm_pal <- inf_pal[1:26]

tmrBySubj = lapply(unique(tsne_results$samp), function(x) {
  ggplot(data = tsne_results[tsne_results$samp == x,], aes(x=tsne_1, y = tsne_2, color = RPclust, alpha = tmr_pos, shape = tmr_pos))+
    geom_point(size=0.9)+
    geom_point(data = tsne_results[tsne_results$samp == x & tsne_results$tmr_pos != "none",], size=1.5, stroke = 2, color = "black")+
    scale_color_manual(values = colorRampPalette(vm_pal)(length(unique(tsne_results$RPclust[tsne_results$samp == x]))))+
    scale_alpha_manual(values = c(.7, .7, .7, .7, .6, .7))+
    scale_shape_manual(values = c(15, 17, 8, 1, 16, 25))+
    labs(title=x, x="tsne 1", y="tsne 2")+
    guides(color=guide_legend(ncol=2))+
    theme(text = element_text(size=14))
})
grid.arrange(grobs = tmrBySubj, ncol = 10)

tmrNoFluBySubj = lapply(unique(tsne_results$samp), function(x) {
  ggplot(data = tsne_results[tsne_results$samp == x,], aes(x=tsne_1, y = tsne_2, color = RPclust, alpha = tmr_pos, shape = tmr_pos))+
    geom_point(size=0.9)+
    geom_point(data = tsne_results[tsne_results$samp == x & tsne_results$tmr_pos != "none" & tsne_results$tmr_pos != "flu",], size=1.5, stroke = 2, color = "black")+
    scale_color_manual(values = colorRampPalette(vm_pal)(length(unique(tsne_results$RPclust[tsne_results$samp == x]))))+
    scale_alpha_manual(values = c(.7, .7, .7, .6, .6, .7))+
    scale_shape_manual(values = c(15, 17, 8, 16, 16, 25))+
    labs(title=x, x="tsne 1", y="tsne 2")+
    guides(color=guide_legend(ncol=2))+
    theme(text = element_text(size=14))
})
grid.arrange(grobs = tmrNoFluBySubj, ncol = 10)


save(for_mario, tsne_results, clust_markers, all_tfm, tfm_fs, mergedExpr, RP_mean, RP_median, RPtmr_counting, 
     agg_tmr_files, cilp_tmr_files, eno_tmr_files, flu_tmr_files, vf_tmr_files, cd4_files, file = "All_phenograph_tSNE_data.RData")

