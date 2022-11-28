rm(list = ls())

# Functions for generating box plot
compare_boxplot <- function(df, x = "Group",y= "Value", fill = "Group", outlier.colour= "white", width=0.5, title.name = NULL,
                            xlabel = NULL, ylabel = "Number of interactions", alpha.value = 1,
                            label.x = 1.5, remove.xtick = FALSE,
                            show.legend = TRUE, x.lab.rot = FALSE, size.text = 10) {
  library(ggplot2)
  gg <- ggpubr::ggboxplot(df, x, y, fill = fill, outlier.colour = outlier.colour, width = width) +
    theme_classic()
  gg <- gg + ylab(ylabel) + xlab(xlabel) +
    labs(title = title.name) +  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) +
    theme(text = element_text(size = size.text), axis.text = element_text(colour="black"))
  gg <- gg + theme(legend.title = element_blank())
  if (remove.xtick) {
    gg <- gg + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
  }
  if (!show.legend) {
    gg <- gg + theme(legend.position = "none")
  }
  if (x.lab.rot) {
    gg <- gg + theme(axis.text.x = element_text(angle = 90, hjust = 1, size=size.text))
  }
  gg
  return(gg)
}

library(Seurat)
# Figs. 1, 2, 3
load("/Users/jinsuoqin/OneDrive/projects/project_Chen/e135_higher_depth/metaData_publishedNC2022/rostral_allcells_Fig1.RData")
load("/Users/jinsuoqin/OneDrive/projects/project_Chen/e135_higher_depth/metaData_publishedNC2022/caudal_allcells_Fig1.RData")
DimPlot(rostral, reduction = "umap")
DimPlot(caudal, reduction = "umap")

# Fig. 4
# load data
load("/Users/jinsuoqin/OneDrive/projects/project_Chen/e135_higher_depth/metaData_publishedNC2022/e135_brachial_LMC_withSubclusters_Fig4.RData")
DimPlot(w10x_bLMC, reduction = "umap", group.by = "quadrant", cols= color.use.major)
DimPlot(w10x_bLMC, reduction = "umap", cols = color.use.subcluster)

# Fig. 5
# Fig. 5a
TF.all <- read.table("/Users/jinsuoqin/OneDrive/projects/project_Chen/e135_higher_depth/brachial_LMC/Mus_musculus_TF.txt", header = T,  sep = "\t")
TF.all <- as.character(unique(TF.all$Symbol))
neuroP.all <- read.table("/Users/jinsuoqin/OneDrive/projects/project_Chen/e135_higher_depth/brachial_LMC/Mus_musculus_Neuropeptide.txt", header = F,  sep = "\t")$V1
hvg <- VariableFeatures(w10x_bLMC)
TF <- intersect(TF.all, hvg)
neuroP <- intersect(neuroP.all, hvg)
df.geneset = read.table(file = "/Users/jinsuoqin/OneDrive/projects/project_Chen/e135_higher_depth/brachial_LMC/geneset_gini.txt", header = T, sep = '\t')
df.geneset <- rbind(df.geneset, data.frame(Group = "TFs", Gene = TF), data.frame(Group = "Neuropeptides", Gene = neuroP))
df.geneset$Group <- factor(df.geneset$Group, levels = c("Housekeeping","Ion channel","Cell adhesion","TFs","Ligands/receptors","Neuropeptides"))

geneset.levels = levels(df.geneset$Group)
df.gini =  data.frame()
for (i in 1:length(geneset.levels)) {
  gene.use <- df.geneset$Gene[df.geneset$Group==geneset.levels[i]]
  data.use <- AverageExpression(w10x_bLMC, features = gene.use)$RNA
  gini.i = apply(as.data.frame(data.use), 1, function(x){REAT::gini(as.numeric(x))})
  df.gini  = rbind(df.gini, data.frame(Value = as.numeric(gini.i), Group = geneset.levels[i], Gene = names(gini.i), Dataavg = apply(data.use, 1, mean)))
}
df.gini$Group <- factor(df.gini$Group, levels = geneset.levels)

compare_boxplot(df.gini, x.lab.rot = T, ylabel = "Gini index") + ylim(c(0,1))
wilcox.test(df.gini$Value[df.gini$Group == "Housekeeping"], df.gini$Value[df.gini$Group == "Neuropeptides"])$p.value


# Fig. 5g: correlation of subcluster-specific TFs and neuropeptides vs. random pairs of TFs and neuropeptides
markers.subcluster <- FindAllMarkers(w10x_bLMC, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
markers.subcluster <- subset(markers.subcluster, p_val_adj < 0.05)
cell.group <- levels(w10x_bLMC)
corr.all <- data.frame()
corr.m.all <- data.frame()
for (i in 1:length(cell.group)) {
  TF.i <- intersect(TF.all, markers.subcluster$gene[markers.subcluster$cluster == cell.group[i] & markers.subcluster$pct.2 < 0.5])
  np.i <- intersect(neuroP.all, markers.subcluster$gene[markers.subcluster$cluster == cell.group[i] & markers.subcluster$pct.2 < 0.5])
  if (length(TF.i) >0 & length(np.i) > 0) {
    data.use.tf <- AverageExpression(w10x_bLMC, features = TF.i)$RNA
    data.use.np <- AverageExpression(w10x_bLMC, features = np.i)$RNA
    corr = as.vector(stats::cor(t(data.use.tf), t(data.use.np), method = "spearman"))
    corr.m <- stats::cor(t(data.use.tf), t(data.use.np), method = "spearman")
    colnames(corr.m) <- np.i
    corr.m <- as.data.frame(as.table(corr.m))
    colnames(corr.m) <- c("TF", "neuropeptide", "spearman")
    corr.m$subcluster <- cell.group[i]
    corr.all <- rbind(data.frame(corr = corr, label = cell.group[i], group = "Combinatorial codes"), corr.all)
    corr.m.all <- rbind(corr.m, corr.m.all)
  }
}
corr.m.all$subcluster <- factor(corr.m.all$subcluster, cell.group)
corr.m.all$TF <- factor(corr.m.all$TF, levels = unique(corr.m.all$TF))
idx <- with(corr.m.all, order(TF, -spearman))
corr.m.all <- corr.m.all[idx, ]
corr.all.avg =corr.all


corr.all <- c()
for (i in 1:length(cell.group)) {
  TF.i <- intersect(TF.all, rownames(w10x_bLMC))
  np.i <- intersect(np.i, rownames(w10x_bLMC))
  if (length(TF.i) >0 & length(np.i) > 0) {
    data.use.tf <- AverageExpression(w10x_bLMC, features = TF.i)$RNA
    data.use.np <- AverageExpression(w10x_bLMC, features = np.i)$RNA
    corr = as.vector(stats::cor(t(data.use.tf), t(data.use.np), method = "spearman"))
    corr.all <- rbind(data.frame(corr = corr, label = "Random pairs", group = "Random pairs"), corr.all)
  }
}
corr.all.avg.random =corr.all
set.seed(1)
samples <- sample.int(nrow(corr.all.avg.random), size = nrow(corr.all.avg))
corr.all.avg.random.subset = corr.all.avg.random[samples, ]

df = rbind(corr.all.avg, corr.all.avg.random.subset)
compare_boxplot(df, x = "group", y= "corr", fill = "group", ylabel = "Correlation", title.name = "TFs and neuropeptides",  remove.xtick = T)
wilcox.test(corr ~ group, df)$p.value


# Fig. 7
library(harmony)
library(Seurat)
library(patchwork)
library(dplyr)
library(ggplot2)
load("/Users/jinsuoqin/OneDrive/projects/project_Chen/e135_higher_depth/metaData_publishedNC2022/e135_brachial_LMC_withSubclusters_Fig4.RData")
load("/Users/jinsuoqin/OneDrive/projects/project_Chen/e135_higher_depth/metaData_publishedNC2022/e135_lumbar_LMC_withSubclusters_S10e.RData")
w10x_bLMC$sample = "Brachial"
w10x_bLMC$subclusters = paste0("Brachial: ", as.character(w10x_bLMC$subclusters))
w10x_lLMC$subclusters = paste0("Lumbar: ", as.character(w10x_lLMC$subclusters))
combined = merge(w10x_bLMC, w10x_lLMC, add.cell.ids = c("Brachial","Lumbar"))

nPC = 30
combined <- combined %>% Seurat::NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(pc.genes = combined@var.genes, npcs = nPC, verbose = FALSE)
# Run Harmony
combined <- combined %>% RunHarmony("sample", plot_convergence = TRUE)
combined <- combined %>%
  RunUMAP(reduction = "harmony", dims = 1:nPC) %>%
  FindNeighbors(reduction = "harmony", dims = 1:nPC)
DimPlot(combined, reduction = "umap", label = F, split.by = "sample", group.by = c("subclusters"))
save(combined, file = "e135_brachial_lumbar_merged_LMC_withSubclusters.RData")

# Fig. 7d
### Codes for performing similarity analyses between pairwise cell groups
# Here is an example showing the similarity between our LMC subclusters from the brachial and lumbar segments

library(batchelor)
# load the integrated data including all LMC neurons from the brachial and lumbar segments
# load("e135_brachial_lumbar_merged_LMC_withSubclusters.RData")
# Get the cell labels given by seperate analysis of each segment
cell.levels.brachial <- c('cl1','cl2','cl3','cl4','rl1','rl2','rl3','rl4','rl5','cm1','cm2','cm3','rm1','rm2','rm3','rm4')
cell.levels.brachial <- paste0("Brachial: ", cell.levels.brachial)
cell.levels.lumbar <- c("c/rl1", "c/rl2", "c/rl3", "c/rl4", "c/rl5", "c/rl6","cm1",   "cm2" ,  "cm3","rm" )
cell.levels.lumbar <- paste0("Lumbar: ", cell.levels.lumbar)
# Get the joint UMAP space
dr.umap = combined@reductions$umap@cell.embeddings
dr.brachial = dr.umap[combined$sample == "Brachial", ]
dr.lumbar = dr.umap[combined$sample == "Lumbar", ]
# Identifying MNN and computing average overlap ratio using different numbers of neighbors
k.mnn.all = c(10, 15, 20)
sim.avg <- 0
for (k.mnn in k.mnn.all) {
  out <- findMutualNN(dr.brachial, dr.lumbar, k1=k.mnn)
  mn = data.frame(first = out$first, second = out$second)
  a = combined$subclusters[combined$sample == "Brachial"]
  a <- factor(a, levels = cell.levels.brachial)
  names(a) = 1:length(a)
  mn$first.name = a[mn$first]
  b = combined$subclusters[combined$sample == "Lumbar"]
  b <- factor(b, levels = cell.levels.lumbar)
  names(b) = 1:length(b)
  mn$second.name = b[mn$second]

  mn$first.name <- factor(mn$first.name, levels = cell.levels.brachial)
  mn$second.name <- factor(mn$second.name, levels = cell.levels.lumbar)
  sim <- c()
  for (i in cell.levels.brachial) {
    sim <- rbind(sim, as.numeric(table(mn$second.name[mn$first.name == i])))
  }
  rownames(sim) <- cell.levels.brachial
  colnames(sim) <- cell.levels.lumbar
  # Computing an average overlap ratio
  sim.avg <- sim.avg + sim/as.numeric(table(a))/2/k.mnn
}
sim.avg <- sim.avg/3

## Visualize the similarity between cell groups using heatmap
sim.scale = sim.avg
library(ComplexHeatmap)
library(colorspace)
source("/Users/jinsuoqin/OneDrive/projects/project_Chen/e135_higher_depth/colorRamp3.R")
color.heatmap = colorRamp3(c(0, 0.05, 0.15, 0.4), c('#4575b4', '#ffffbf','#fdae61','#d73027'))
ht = Heatmap(sim.scale,col = color.heatmap, cluster_rows = FALSE, cluster_columns = FALSE,show_heatmap_legend = T, row_names_side = "left",row_title = "Brachial",column_title = "Similarity", name = NULL,
             cell_fun = function(j, i, x, y, width, height, fill) {
               if(sim.scale[i, j] > 0.04)
                 grid.text(sprintf("%.2f", sim.scale[i, j]), x, y, gp = gpar(fontsize = 8))
             })
draw(ht)


# Fig. 8
load("/Users/jinsuoqin/OneDrive/projects/project_Chen/e135_higher_depth/metaData_publishedNC2022/rostral_MMC_Fig8.RData")
DimPlot(w10x.MMC, reduction  = "umap", label = T)




###### Codes for estimating the number of clusters in scRNA-seq data
rm(list = ls())

## R codes for computing consensus matrix
library(Seurat)
library(dplyr)

setwd("/Users/suoqinjin/Google Drive/projects/project_Chen/e135_higher_depth/brachial_LMC")
# load Seurat object
load("e135_brachial_LMC_cLMCl.RData")
w10x <- w10x_bLMC_cLMCl

N <- ncol(w10x@data)
C <- matrix(0, N, N)
resR <- seq(0.1,3,by = 0.05) # a range of resolution for clustering
for (res in resR) {
  w10x <- FindClusters(object = w10x, reduction.type = "pca", dims.use = 1:numPC, resolution = res, algorithm = 1,save.SNN = TRUE,print.output = 0,force.recalc = T)
  print(length(unique(w10x@ident)))
  clusIndex <- as.numeric(as.character(w10x@ident))
  adjMat <- as.numeric(outer(clusIndex, clusIndex, FUN = "==")) - outer(1:N, 1:N, "==")
  C <- C + adjMat
}
CM <- C/length(resR)
# save the consensus matrix for running the Matlab code 'determineNumClusters.m'
write.table(CM,file = "consensusMatrix_brachial_LMC_cLMCl.txt",sep = '\t', row.names = F, col.names = F)

## Matlab codes for determining the number of clusters
# Please do this in Matlab
fileName = 'brachial_LMC_cLMCl'; # the file name of the figure to save
CM = importdata(['consensusMatrix_',fileName,'.txt']);
[numCluster, numCluster0, eigenvalues] = determineNumClusters(CM,fileName);
# Output:
 #  numCluster: Number of inferred clusters
 #  numCluster0: the minimum number of inferred clusters based on the number of zero eigenvalues
 #  eigenvalues: eigenvalues


# below is the Matlab codes for determining the number of clusters
function [numCluster, numCluster0, eigenvalues] = determineNumClusters(CM,fileName)
% Estimation of the number of clusters from the consensus matrix
%
% Input:
%   CM: consensus matrix constructed by running Seurat using multiple resolutions
%   fileName: a char giving the file name of the figure to save
% Output:
%   numCluster: Number of inferred clusters
%   numCluster0: the minimum number of inferred clusters based on the number of zero eigenvalues
%   eigenvalues: eigenvalues of the graph Laplacian
%
n = size(CM,1);
CM(1:n+1:end) = 1;

tol = 0.01;
numEigs = min(100,size(CM,1));

D = diag(CM*ones(n,1));
Prw = eye(n) - D^(-1/2)*CM*D^(-1/2);
all_eigs = real(eigs(Prw,numEigs,'sm'));
ZZ = sort(abs(real(all_eigs)));
numCluster0 = length(find(ZZ<=tol));
tau = 0.3;
CM(CM <= tau) = 0;
CM = (1/2)*(CM + CM');

D = diag(CM*ones(n,1));
Prw = eye(n) - D^(-1/2)*CM*D^(-1/2);
all_eigs = real(eigs(Prw,numEigs,'sm'));

zz = sort(abs(real(all_eigs)));

gap = zz(2:end) - zz(1:end-1);
[~,numCluster] = max(gap);


numCluster0 = length(find(zz<=tol));
display('Number of cluster based on zero eigenvalues & Largest gap ');
display([numCluster0 numCluster]);

eigenvalues = zz;

figure('position', [600, 200, 250, 200])
scatter(1:min([30 size(eigenvalues,1)]),eigenvalues(1:min([30 size(eigenvalues,1)])),20,'k','filled');
hold on
scatter(numCluster,eigenvalues(numCluster),40,'r')
box off;
set(gca,'LineWidth',1);
set(gca,'FontSize',8,'FontName','Arial');
xlabel('Number of clusters','FontSize',10);
ylabel('Eigenvalue of graph Laplacian','FontSize',10);
title(['Inferred number of clusters: ', num2str(numCluster),'; Min number: ',num2str(numCluster0)],'FontSize',10)
saveas(gcf,['eigGap_',fileName,'.pdf'])






