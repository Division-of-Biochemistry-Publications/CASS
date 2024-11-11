#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Initially created on Wed Jun 21 13:02:18 2023

@author: David Rüscher (github.com/DavidRuescher95)

Keep in mind that results might differ slightly as the UMAP and Louvain algorithms are
stochastic in nature and no fixed seed was used.

Versions:
    - Python 3.9.13
        - umap-learn 0.5.3
        - numpy 1.21.5
        - pandas 1.5.2
        - matplotlib 3.6.2
        - networkx 3.0
        - scipy 1.7.3
        - scikit-learn 1.0.2

    - Linux binaries:
        - FastQC 0.11.9
        - multiQC 1.13
        - bbduk 38.97
        - STAR 2.7.10a
        - FeatureCounts 2.03

    -R 4.3.0:
        - DESeq2 1.40.1
        - clusterProfiler 4.8.1
        - limma 3.54.1
        - UpSetR 1.4.0
        - tidyverse 2.0.0
        
"""
# =============================================================================
#  import modules
# =============================================================================

import os
import subprocess
import re
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as matplotlib
import numpy as np
import itertools
import string
from scipy.stats import rankdata
from sklearn.preprocessing import StandardScaler
from Bio import SeqIO
import math

# =============================================================================
#  change working directory
# =============================================================================

RScriptPath = "/usr/local/bin/Rscript"
wd = "/Users/david/iCloud/PhD/PUBLICATIONS/CRYOSECTIONING_2024/XX_code_results/rüscher_et_al_cryosectioning_2024/"
os.chdir(wd)
dataPath = "../RNASEQ/"
scriptPath = "./scripts/cryosectioning/"

# =============================================================================
# load custom classes
# =============================================================================

from modules.SampleCommunitiesPCA import SampleCommunitiesPCA
from modules.SampleUMAP import SampleUMAP
from modules.CoRegulatoryNetworkAnalysis import CRN

# =============================================================================
# Run DESeq2 Rscript
# =============================================================================

subprocess.run([f"{RScriptPath}", f"{scriptPath}DESEQ2.R"])

# =============================================================================
#  load data
# ============================================================================

metaData = pd.read_csv(
    filepath_or_buffer=f"{dataPath}00_META_DATA/META_DATA.csv",
    index_col=0
    )
vst = pd.read_csv(
    filepath_or_buffer=f"{dataPath}08_DESEQ2/VST_LIMMA.csv",
    index_col=0
    )
vstRaw = pd.read_csv(
    filepath_or_buffer=f"{dataPath}08_DESEQ2/VST.csv",
    index_col=0
    )
norm=pd.read_csv(
    filepath_or_buffer=f"{dataPath}08_DESEQ2/NORMALIZED.csv",
    index_col=0
    )

# =============================================================================
# sample clustering
# =============================================================================
# create and define path
outputPath = f"{dataPath}09_SAMPLE_CLUSTERING/"
figurePath = f"{outputPath}FIGURES/"
os.makedirs(figurePath,exist_ok=True)


umapResRaw = []

for tissue in metaData.Tissue.unique():
    print(tissue)
    
    umap = SampleUMAP(vstRaw,norm,metaData,"Tissue",tissue)
    umap.scale()
    umap.depthFilter(50)
    umap.umap(45,0.01,"cosine")
    umap.extractData()
    # join with meta data and save file
    umapResRaw.append(umap.Export["UMAP"])
    
umapResRaw = metaData.join(pd.concat(umapResRaw))
umapResRaw.to_csv(f"{outputPath}UMAP_RAW.csv")

# plot


def umapPlot(plotData, colorValues, markerValues, figurePath):
    
    markers = ["o","s","^","P","p","v"]
    colors = plt.cm.plasma(np.linspace(0,1,len(plotData[colorValues].unique())))
    plotData["color"] = rankdata(plotData[colorValues],method="dense")-1
    plotData["marker"] = rankdata(plotData[markerValues],method="dense")-1
    f1 = lambda c: matplotlib.lines.Line2D([], [], linewidth=0, color = c, marker="o",
                                          markeredgecolor = "black", markeredgewidth=.1)
    handles1 = [f1(colors[i]) for i in plotData.color.unique()]
    labels1 = plotData[colorValues].unique()
    f2 = lambda m: matplotlib.lines.Line2D([], [], linewidth=0, color = "black", marker=m,
                                          markeredgecolor = "black", markeredgewidth=.1)
    handles2 = [f2(markers[i]) for i in plotData.marker.unique()]
    labels2 = plotData[markerValues].unique()
    plt.plasma()
    
    plotDataGrouped = plotData.groupby("Tissue")
    
    cm = 1/2.54
    fig, ax = plt.subplots(2,1)
    fig.set_size_inches(8*cm,10*cm)
    
    for index, tissue in enumerate(metaData.Tissue.unique()):
    
        plotDataTissue = plotDataGrouped.get_group(tissue)
        
        for markerIndex, marker in enumerate(plotDataTissue.marker):
        
            ax[index].scatter(x = plotDataTissue.UMAP1[markerIndex],y = plotDataTissue.UMAP2[markerIndex],
                          c = colors[plotDataTissue.color[markerIndex].astype(int)],s=32,
                          marker=markers[plotDataTissue.marker[markerIndex]],
                          linewidths=.25,edgecolor="black")
            ax[index].set_title(tissue, size = 9, weight="bold")
            ax[index].tick_params(labelsize=6,size=2)
            ax[index].spines[["top","right"]].set_visible(False)
    
    fig.supxlabel("UMAP1", fontsize = 8, weight = "bold")
    fig.supylabel("UMAP2", fontsize = 8, weight = "bold")
    legend1 = plt.legend(handles1,labels1,bbox_to_anchor = [1,1.60,0,0],
                         title=colorValues, loc = "center left",markerscale = 0.75,
                         fontsize=6,title_fontsize=7,alignment="left",handletextpad=0.1)
    legend1.get_frame().set(linewidth=.75,edgecolor="black",facecolor="None")
    legend2 = plt.legend(handles2,labels2,bbox_to_anchor = [1,.525,0,0],
                         title=markerValues, loc = "center left",markerscale = 0.75,
                         fontsize=6,title_fontsize=7,alignment="left",handletextpad=0.1)
    legend2.get_frame().set(linewidth=.75,edgecolor="black",facecolor="None")
    fig.add_artist(legend1)
    fig.add_artist(legend2)
    legend1.set_in_layout(False)
    legend2.set_in_layout(False)
    plt.tight_layout()
    plt.subplots_adjust(top = .925, bottom = .1,right = 0.82, left = 0.15,hspace=0.35)
    
    fig.savefig(figurePath,dpi=600,transparent=True)
    plt.show()

plotData = pd.read_csv(f"{outputPath}UMAP_RAW.csv",index_col=0)
umapPlot(plotData, "Distance", "Plant", f"{figurePath}UMAP_vstRaw.png")

# =============================================================================
# Stem analysis
# =============================================================================
stem = SampleCommunitiesPCA(vst,norm,metaData,"Tissue","Stem")
stem.scale()
stem.depthFilter(50)
stem.pca(nFeatures=5000)
stem.louvain(10, 90, 7, "cosine", 1.5)
stem.umap(45,0.01,"cosine")
stem.extractData()

clusterMedians = stem.metaData.groupby("Cluster").median("Distance").sort_values("Distance")
for clusterIndex, cluster in enumerate(clusterMedians.index):
    stem.metaData.loc[stem.metaData.Cluster == cluster,"Cluster"] = string.ascii_uppercase[clusterIndex]

# =============================================================================
# SR analysis
# =============================================================================
SR = SampleCommunitiesPCA(vst,norm,metaData,"Tissue","SR")
SR.scale()
SR.depthFilter(50)
SR.pca(nFeatures=5000)
SR.louvain(10, 90, 7, "cosine", 1.5)
SR.umap(45,0.01,"cosine")
SR.extractData()

clusterMedians = SR.metaData.groupby("Cluster").median("Distance").sort_values("Distance")
for clusterIndex, cluster in enumerate(clusterMedians.index):
    SR.metaData.loc[SR.metaData.Cluster == cluster,"Cluster"] = string.ascii_uppercase[clusterIndex]

# =============================================================================
# save data 
# =============================================================================

metaDataClusters = pd.concat(
    [stem.metaData,SR.metaData]
    ).loc[metaData.index,:]
metaDataClusters.to_csv(f"{outputPath}META_DATA_CLUSTERS.csv")

def concatSaveData(name,path):
    file = metaDataClusters.join(pd.concat([stem.Export[name],SR.Export[name]]
                                        )).sort_values(["Tissue","Cluster","Plant"])
    file.to_csv(path + name + ".csv")
    return file

pcaClustering = concatSaveData("PCA",outputPath)
umapClustering = concatSaveData("UMAP",outputPath)

stem.Export["VarPCA"]["Tissue"] = "Stem"
SR.Export["VarPCA"]["Tissue"] = "SR"
pcaVarClustering = pd.concat([stem.Export["VarPCA"],SR.Export["VarPCA"]])
pcaVarClustering.to_csv(f"{outputPath}VarPCA.csv",index=False)

# =============================================================================
# cluster umap plots
# =============================================================================

plotData = pd.read_csv(f"{outputPath}UMAP.csv",index_col=0).sort_values(["Distance","Cluster"])

umapPlot(plotData, "Distance", "Plant", f"{figurePath}UMAP_PLANT.png")
umapPlot(plotData, "Distance", "Cluster", f"{figurePath}UMAP_CLUSTER.png")

# =============================================================================
# Cluster distance Boxplot
# =============================================================================
metaDataClusters = pd.read_csv(f"{outputPath}META_DATA_CLUSTERS.csv", index_col=0)
plotData = metaDataClusters.copy()

markers = np.array(["o","s","^"])
colors = plt.cm.plasma(np.linspace(0,1,len(plotData.Tissue.unique())))
plotData["color"] = rankdata(plotData.Tissue,method="dense")-1
plotData["marker"] = rankdata(plotData.Plant,method="dense")-1
f1 = lambda c: matplotlib.lines.Line2D([], [], linewidth=0, color = c, marker="o",
                                      markeredgecolor = "black", markeredgewidth=.1)
handles1 = [f1(colors[i]) for i in plotData.color.unique()]
labels1 = plotData.Tissue.unique()
f2 = lambda m: matplotlib.lines.Line2D([], [], linewidth=0, color = "black", marker=m,
                                      markeredgecolor = "black", markeredgewidth=.1)
handles2 = [f2(markers[i]) for i in plotData.marker.unique()]
labels2 = plotData.Plant.unique()
plotData.Cluster=plotData.Cluster.astype(str)
plotData["color"] = rankdata(plotData.Tissue,method="dense")-1
plotDataGrouped = plotData.groupby("Tissue")

cm = 1/2.54
fig, ax = plt.subplots(2,1, sharey = True)

for i in ax: i.set_axis_off()

fig.set_size_inches(8*cm,10*cm)

for index, tissue in enumerate(metaData.Tissue.unique()): 
    plotDataTissue = plotDataGrouped.get_group(tissue)
    clusters = plotDataTissue.Cluster.unique()
    f = lambda cluster: plotDataTissue.loc[plotDataTissue.Cluster == cluster,"Distance"].to_numpy()
    dataset = [f(cluster) for cluster in plotDataTissue.Cluster.unique()]
    ax[index].set_axis_on()
    box = ax[index].boxplot(dataset,showfliers=False,widths=0.55,patch_artist=True,
                                    positions=range(len(plotDataTissue.Cluster.unique())), zorder=0)
    for i,j in zip(box["boxes"],box["medians"]):
        i.set(facecolor=colors[plotDataTissue.color][0],edgecolor="black",linewidth=.5,alpha=.75)
        j.set(c="black",linewidth=0.5)
    for i,j in zip(box["whiskers"],box["caps"]):
        i.set(c="black",linewidth=.5)
        j.set(c="black",linewidth=.5)
    ax[index].tick_params(labelsize=6,size=2)
    ax[index].set_title(tissue, size = 9, weight="bold")
    ax[index].set_xticklabels(plotDataTissue.Cluster.unique())
    ax[index].spines[["top","right"]].set_visible(False)

    plotDataTissueGrouped = plotDataTissue.groupby("Plant")
    
    for markerIndex, plant in enumerate(plotDataTissue.Plant.unique()):
        
        plotDataTissuePlant = plotDataTissueGrouped.get_group(plant)
        
        f = lambda cluster: plotDataTissuePlant.loc[plotDataTissuePlant.Cluster == cluster,"Distance"].to_numpy()
        dataset = [f(cluster) for cluster in plotDataTissuePlant.Cluster.unique()]
        
        for cluster, value in enumerate(dataset):
            group = np.random.normal(cluster, 0.125, size=len(value))
            scatter = ax[index].scatter(group,value,s=12,linewidths=.35,marker=markers[plotDataTissuePlant.marker.values[0]],
                                        c="white",edgecolor="black",zorder=1,alpha=0.9)
            
fig.supxlabel("Sample cluster", fontsize = 8, weight="bold")
fig.supylabel("Distance [µm]", fontsize = 8, weight="bold")

legend1 = plt.legend(handles1,labels1,bbox_to_anchor = [1,2.2,0,0],
                     title="Tissue", loc = "center left",markerscale = 0.75,
                     fontsize=6,title_fontsize=7,alignment="left",handletextpad=0.1)
legend1.get_frame().set(linewidth=.75,edgecolor="black",facecolor="None")

legend2 = plt.legend(handles2,labels2,bbox_to_anchor = [1,1.75,0,0],
                     title="Plant", loc = "center left",markerscale = 0.75,
                     fontsize=6,title_fontsize=7,alignment="left",handletextpad=0.1)
legend2.get_frame().set(linewidth=.75,edgecolor="black",facecolor="None")

fig.add_artist(legend1)
fig.add_artist(legend2)
legend1.set_in_layout(False)
legend2.set_in_layout(False)


plt.tight_layout()
plt.subplots_adjust(top = .925, bottom = .115,right = 0.82, left = 0.15,hspace=0.35)
fig.savefig(f"{figurePath}CLUSTER_DISTANCE.png",dpi=600,transparent=True)
plt.show()


# rename and replot 

metaDataClusters = pd.read_csv(f"{outputPath}META_DATA_CLUSTERS.csv", index_col=0)
metaDataClusters["clusterName"] = np.nan
metaDataClusters.loc[metaDataClusters.Cluster == "A","clusterName"] = "Phloem"
metaDataClusters.loc[metaDataClusters.Cluster == "B","clusterName"] = "Phloem/Cambium"
metaDataClusters.loc[metaDataClusters.Cluster == "C","clusterName"] = "Cambium/Xylem"
metaDataClusters.loc[metaDataClusters.Cluster == "D","clusterName"] = "Xylem"

metaDataClusters["clusterNameShort"] = np.nan
metaDataClusters.loc[metaDataClusters.Cluster == "A","clusterNameShort"] = "P"
metaDataClusters.loc[metaDataClusters.Cluster == "B","clusterNameShort"] = "P/C"
metaDataClusters.loc[metaDataClusters.Cluster == "C","clusterNameShort"] = "C/X"
metaDataClusters.loc[metaDataClusters.Cluster == "D","clusterNameShort"] = "X"

metaDataClusters.to_csv(f"{outputPath}META_DATA_CLUSTERS.csv")


metaDataClusters = pd.read_csv(f"{outputPath}META_DATA_CLUSTERS.csv", index_col=0)
plotData = metaDataClusters.copy()
plotData.Cluster = metaDataClusters.Cluster + "\n(" + metaDataClusters.clusterNameShort +")"

markers = np.array(["o","s","^"])
colors = plt.cm.plasma(np.linspace(0,1,len(plotData.Tissue.unique())))
plotData["color"] = rankdata(plotData.Tissue,method="dense")-1
plotData["marker"] = rankdata(plotData.Plant,method="dense")-1
f1 = lambda c: matplotlib.lines.Line2D([], [], linewidth=0, color = c, marker="o",
                                      markeredgecolor = "black", markeredgewidth=.1)
handles1 = [f1(colors[i]) for i in plotData.color.unique()]
labels1 = plotData.Tissue.unique()
f2 = lambda m: matplotlib.lines.Line2D([], [], linewidth=0, color = "black", marker=m,
                                      markeredgecolor = "black", markeredgewidth=.1)
handles2 = [f2(markers[i]) for i in plotData.marker.unique()]
labels2 = plotData.Plant.unique()
plotData.Cluster=plotData.Cluster.astype(str)
plotData["color"] = rankdata(plotData.Tissue,method="dense")-1
plotDataGrouped = plotData.groupby("Tissue")

cm = 1/2.54
fig, ax = plt.subplots(2,1, sharey = True)

for i in ax: i.set_axis_off()

fig.set_size_inches(8*cm,10*cm)

for index, tissue in enumerate(metaData.Tissue.unique()): 
    plotDataTissue = plotDataGrouped.get_group(tissue)
    clusters = plotDataTissue.Cluster.unique()
    f = lambda cluster: plotDataTissue.loc[plotDataTissue.Cluster == cluster,"Distance"].to_numpy()
    dataset = [f(cluster) for cluster in plotDataTissue.Cluster.unique()]
    ax[index].set_axis_on()
    box = ax[index].boxplot(dataset,showfliers=False,widths=0.55,patch_artist=True,
                                    positions=range(len(plotDataTissue.Cluster.unique())), zorder=0)
    for i,j in zip(box["boxes"],box["medians"]):
        i.set(facecolor=colors[plotDataTissue.color][0],edgecolor="black",linewidth=.5,alpha=.75)
        j.set(c="black",linewidth=0.5)
    for i,j in zip(box["whiskers"],box["caps"]):
        i.set(c="black",linewidth=.5)
        j.set(c="black",linewidth=.5)
    ax[index].tick_params(labelsize=6,size=2)
    ax[index].set_title(tissue, size = 9, weight="bold")
    ax[index].set_xticklabels(plotDataTissue.Cluster.unique())
    ax[index].spines[["top","right"]].set_visible(False)

    plotDataTissueGrouped = plotDataTissue.groupby("Plant")
    
    for markerIndex, plant in enumerate(plotDataTissue.Plant.unique()):
        
        plotDataTissuePlant = plotDataTissueGrouped.get_group(plant)
        
        f = lambda cluster: plotDataTissuePlant.loc[plotDataTissuePlant.Cluster == cluster,"Distance"].to_numpy()
        dataset = [f(cluster) for cluster in plotDataTissuePlant.Cluster.unique()]
        
        for cluster, value in enumerate(dataset):
            group = np.random.normal(cluster, 0.125, size=len(value))
            scatter = ax[index].scatter(group,value,s=12,linewidths=.35,marker=markers[plotDataTissuePlant.marker.values[0]],
                                        c="white",edgecolor="black",zorder=1,alpha=0.9)
            
fig.supxlabel("Sample cluster", fontsize = 8, weight="bold")
fig.supylabel("Distance [µm]", fontsize = 8, weight="bold")

legend1 = plt.legend(handles1,labels1,bbox_to_anchor = [1,2.2,0,0],
                     title="Tissue", loc = "center left",markerscale = 0.75,
                     fontsize=6,title_fontsize=7,alignment="left",handletextpad=0.1)
legend1.get_frame().set(linewidth=.75,edgecolor="black",facecolor="None")

legend2 = plt.legend(handles2,labels2,bbox_to_anchor = [1,1.75,0,0],
                     title="Plant", loc = "center left",markerscale = 0.75,
                     fontsize=6,title_fontsize=7,alignment="left",handletextpad=0.1)
legend2.get_frame().set(linewidth=.75,edgecolor="black",facecolor="None")

fig.add_artist(legend1)
fig.add_artist(legend2)
legend1.set_in_layout(False)
legend2.set_in_layout(False)


plt.tight_layout()
plt.subplots_adjust(top = .925, bottom = .115,right = 0.82, left = 0.15,hspace=0.35)
fig.savefig(f"{figurePath}CLUSTER_DISTANCE_RENAMED.png",dpi=600,transparent=True)
plt.show()

# =============================================================================
# Run DESeq2 LRT Rscript
# =============================================================================

subprocess.run([f"{RScriptPath}", f"{scriptPath}LRT.R"])


# =============================================================================
# CRN analysis
# =============================================================================

# create and define path
outputPath = f"{dataPath}10_CRN/"
figurePath = f"{outputPath}FIGURES/"
os.makedirs(figurePath,exist_ok=True)

lrt = pd.read_csv(f"{dataPath}08_DESEQ2/LRT.csv",index_col=0)
metaDataClusters = pd.read_csv(f"{dataPath}09_SAMPLE_CLUSTERING/META_DATA_CLUSTERS.csv",index_col=0)

stem_CRN = CRN(vst,norm,metaDataClusters,lrt,"Tissue","Stem")
stem_CRN.scale()
stem_CRN.pValFilter(0.001)
stem_CRN.corrNetwork(filterType="pValue",corrThrsh=0.8, runPCA = True, minPCs=10, minPctVar=90)
stem_CRN.networkGeneration(remove = True)
stem_CRN.louvain(minEdges=25, minClusterSize=50, minClusterDegreeRatio = 0, resolution=1.0)
stem_CRN.umap(25,0,"correlation")
stem_CRN.extractData()


SR_CRN = CRN(vst,norm,metaDataClusters,lrt,"Tissue","SR")
SR_CRN.scale()
SR_CRN.pValFilter(0.001)
SR_CRN.corrNetwork(filterType="pValue",corrThrsh=0.8, runPCA = True, minPCs=10, minPctVar=90)
SR_CRN.networkGeneration()
SR_CRN.louvain(minEdges=25, minClusterSize=50, minClusterDegreeRatio = 0, resolution=1.1)
SR_CRN.umap(25,0,"correlation")
SR_CRN.extractData()


# =============================================================================
# combine data 
# =============================================================================

SR_clusters = SR_CRN.Clusters
SR_clusters["Tissue"] = "SR"
stem_clusters = stem_CRN.Clusters
stem_clusters["Tissue"] = "Stem"
clusters = pd.concat([SR_clusters,stem_clusters]).sort_index().sort_values(["Tissue","Cluster"])
clusters.to_csv(f"{outputPath}CLUSTERS.csv")
SR_umap = SR_CRN.Export["UMAP"]
SR_umap["Tissue"] = "SR"
stem_umap = stem_CRN.Export["UMAP"]
stem_umap["Tissue"] = "Stem"
umapGene = pd.concat([SR_umap,stem_umap])
umapGene = clusters.reset_index().merge(umapGene.reset_index()).set_index("index")
umapGene.to_csv(f"{outputPath}umapGene.csv")
scaled = pd.concat([SR_CRN.Export["Scaled"],stem_CRN.Export["Scaled"]],axis=1)
scaled.to_csv(f"{outputPath}SCALED.csv")


# =============================================================================
# gene umapGene plot
# =============================================================================

umapGene = pd.read_csv(f"{outputPath}umapGene.csv",index_col=0)
umapGene.Cluster=umapGene.Cluster.astype(str)

plotData = umapGene.copy()
colors = plt.cm.plasma(np.linspace(0,1,len(plotData.Cluster.unique())))
colors[0] = np.array([0, 0, 0, 1])
plotData["color"] = rankdata(plotData.Cluster,method="dense")-1
f = lambda c: matplotlib.lines.Line2D([], [], linewidth=0, color = c, marker="o",
                                      markeredgecolor = "black", markeredgewidth=.1)
handles = [f(colors[i]) for i in plotData.color.unique()]
labels = plotData.Cluster.unique()

plotData = plotData.groupby(["Tissue"])
cm = 1/2.54
fig, ax = plt.subplots(2,1)
fig.set_size_inches(8*cm,10*cm)
for index, tissue in enumerate(metaData.Tissue.unique()):
    tmp = plotData.get_group(tissue)
    ax[index].scatter(x = tmp.UMAP1,y = tmp.UMAP2,
                  c = colors[tmp.color.astype(int)],marker="o",s=4,
                  linewidths=0,edgecolor="black",alpha=.1)
    ax[index].tick_params(labelsize=6,size=2)
    ax[index].set_title(tissue, size = 9,weight="bold")
    ax[index].spines[["top","right"]].set_visible(False)
    
fig.supxlabel("UMAP1", fontsize = 8,weight="bold")
fig.supylabel("UMAP2", fontsize = 8,weight="bold")

legend = plt.legend(handles,labels,bbox_to_anchor = [1,2.1,0,0],
                     title="Cluster", loc = "center left",markerscale = 0.75,
                     fontsize=6,title_fontsize=7,alignment="left",handletextpad=0.1)
legend.get_frame().set(linewidth = 0.75, edgecolor="black",facecolor="None")

fig.add_artist(legend)
legend.set_in_layout(False)

plt.tight_layout()
plt.subplots_adjust(top = .9, bottom = .10,right = 0.835, left = 0.135,hspace=.35)
fig.savefig(f"{figurePath}UMAP_CLUSTERS.png",dpi=600,transparent=True)
plt.show()

# =============================================================================
# cluster expression plot 
# =============================================================================

scaled = pd.read_csv(f"{outputPath}SCALED.csv",index_col=0)


metaDataClusters = pd.read_csv(f"{dataPath}09_SAMPLE_CLUSTERING/META_DATA_CLUSTERS.csv",index_col=0)
sampleClusters = metaDataClusters.copy().reset_index()
sampleClusters["sampleCluster"] = sampleClusters.Cluster + "\n(" + sampleClusters.clusterNameShort +")"
sampleClusters = sampleClusters.drop(columns = ["Cluster","clusterName","clusterNameShort"])

umapGene = pd.read_csv(f"{outputPath}umapGene.csv",index_col=0)
umapGene.Cluster=umapGene.Cluster.astype(str)

plotData = umapGene.copy().loc[umapGene.Cluster != "-1",:]
colors = plt.cm.plasma(np.linspace(0,1,len(plotData.Cluster.unique())))
plotData["color"] = rankdata(plotData.Cluster,method="dense")-1
f = lambda c: matplotlib.lines.Line2D([], [], linewidth=0, color = c, marker="o",
                                      markeredgecolor = "black", markeredgewidth=.1)
handles = [f(colors[i]) for i in plotData.color.unique()]
labels = plotData.Cluster.unique()

plotDataGrouped = plotData.groupby("Tissue")

cm = 1/2.54
fig, ax = plt.subplots(2,len(plotData.Cluster.unique()), sharey = True)
for i in ax: 
    for j in i: j.set_axis_off()
    
fig.set_size_inches(10*cm,10*cm)
for index1, tissue in enumerate(metaData.Tissue.unique()): 
    samples = sampleClusters.loc[sampleClusters.Tissue == tissue,"Sample"]
    plotDataTissue = plotDataGrouped.get_group(tissue).join(scaled.loc[:,samples]
        ).melt(id_vars = plotData.columns,var_name="Sample").merge(sampleClusters,how="left")
    clusters = plotDataTissue.Cluster.unique()
    plotDataTissue = plotDataTissue.groupby("Cluster")
    for index2, geneCluster in enumerate(clusters):
        plotDataTissueCluster = plotDataTissue.get_group(geneCluster)
        f = lambda cluster: plotDataTissueCluster.loc[plotDataTissueCluster.sampleCluster == cluster,"value"].to_numpy()
        dataset = [f(i) for i in plotDataTissueCluster.sampleCluster.unique()]
        ax[index1,index2].set_axis_on()
        violin = ax[index1,index2].violinplot(dataset,showextrema=False,
                                                   widths=0.85,positions=range(len(plotDataTissueCluster.sampleCluster.unique())))
        box = ax[index1,index2].boxplot(dataset,showfliers=False,widths=0.2,patch_artist=True,
                                        positions=range(len(plotDataTissueCluster.sampleCluster.unique())))
        ax[index1,index2].tick_params(labelsize=6,size=2)
        ax[index1,index2].set_title( tissue + " " +geneCluster, size = 6,weight="bold")
        ax[index1,index2].set_xticklabels(plotDataTissueCluster.sampleCluster.unique())
        ax[index1,index2]
        for i in violin["bodies"]:
            i.set(color=colors[plotDataTissueCluster.color],alpha=.75,linewidth=.35)
        for i,j in zip(box["boxes"],box["medians"]):
            i.set(facecolor="white",edgecolor="black",linewidth=.35)
            j.set(c="black",linewidth=0.35)
        for i,j in zip(box["whiskers"],box["caps"]):
            i.set(c="black",linewidth=.35)
            j.set(c="black",linewidth=0)
        plt.ylim([-6.5,6.5])
        ax[index1,index2].spines[["top","right"]].set_visible(False)
fig.supxlabel("Sample cluster", fontsize = 8,weight="bold")
fig.supylabel("Z-Score", fontsize = 8,weight="bold")
legend = plt.legend(handles,labels,bbox_to_anchor = [1,2.1,0,0],
                     title="Cluster", loc = "center left",markerscale = 0.75,
                     fontsize=6,title_fontsize=7,alignment="left",handletextpad=0.1)
legend.get_frame().set(linewidth = 0.75, edgecolor="black",facecolor="None")
fig.add_artist(legend)
legend.set_in_layout(False)
plt.tight_layout()
plt.subplots_adjust(top = .9, bottom = .10,right = 0.87, left = 0.11,hspace=.35,wspace=.1)
fig.savefig(f"{figurePath}CLUSTER_EXPRESSION.png",dpi=600,transparent=True)
plt.show()


# rename and re-plot

clusters = pd.read_csv(f"{outputPath}CLUSTERS.csv",index_col=0)
clusters["clusterName"] = np.nan
clusters.loc[clusters.Cluster == -1,"clusterName"] = "Removed"
clusters.loc[(clusters.Tissue == "Stem") & (clusters.Cluster == 0),"clusterName"] = "Phloem/Xylem"
clusters.loc[(clusters.Tissue == "Stem") & (clusters.Cluster == 1),"clusterName"] = "Phloem/Cambium"
clusters.loc[(clusters.Tissue == "Stem") & (clusters.Cluster == 2),"clusterName"] = "Cambium/Xylem"
clusters.loc[(clusters.Tissue == "Stem") & (clusters.Cluster == 3),"clusterName"] = "Xylem"
clusters.loc[(clusters.Tissue == "Stem") & (clusters.Cluster == 4),"clusterName"] = "Phloem"

clusters.loc[(clusters.Tissue == "SR") & (clusters.Cluster == 2),"clusterName"] = "Phloem/Xylem"
clusters.loc[(clusters.Tissue == "SR") & (clusters.Cluster == 0),"clusterName"] = "Phloem/Cambium"
clusters.loc[(clusters.Tissue == "SR") & (clusters.Cluster == 4),"clusterName"] = "Cambium/Xylem"
clusters.loc[(clusters.Tissue == "SR") & (clusters.Cluster == 1),"clusterName"] = "Xylem"
clusters.loc[(clusters.Tissue == "SR") & (clusters.Cluster == 3),"clusterName"] = "Phloem"

clusters["clusterNameShort"] = np.nan
clusters.loc[clusters.clusterName == "Removed","clusterNameShort"] = "Removed"
clusters.loc[clusters.clusterName == "Phloem","clusterNameShort"] = "P"
clusters.loc[clusters.clusterName == "Phloem/Cambium","clusterNameShort"] = "P/C"
clusters.loc[clusters.clusterName == "Phloem/Xylem","clusterNameShort"] = "P/X"
clusters.loc[clusters.clusterName == "Cambium/Xylem","clusterNameShort"] = "C/X"
clusters.loc[clusters.clusterName == "Xylem","clusterNameShort"] = "X"

clusters.to_csv(f"{outputPath}CLUSTERS.csv")

clusters = pd.read_csv(f"{outputPath}CLUSTERS.csv",index_col=0)

plotData = clusters.copy().loc[clusters.clusterName != "Removed",:]
clusters = ["Phloem", "Phloem/Cambium", "Phloem/Xylem", "Cambium/Xylem", "Xylem"]
#clusters = ["P", "P/C", "P/X", "C/X", "X"]
colors = plt.cm.plasma(np.linspace(0,1,len(plotData.Cluster.unique())))
plotData["color"] = rankdata(plotData.Cluster,method="dense")-1
f = lambda c: matplotlib.lines.Line2D([], [], linewidth=0, color = c, marker="o",
                                      markeredgecolor = "black", markeredgewidth=.1)
handles = [f(colors[i]) for i in plotData.color.unique()]
labels = plotData.Cluster.unique()

plotDataGrouped = plotData.groupby("Tissue")

cm = 1/2.54
fig, ax = plt.subplots(2,len(plotData.clusterName.unique()), sharey = True)
for i in ax: 
    for j in i: j.set_axis_off()
    
fig.set_size_inches(10*cm,10*cm)
for index1, tissue in enumerate(metaData.Tissue.unique()): 
    samples = sampleClusters.loc[sampleClusters.Tissue == tissue,"Sample"]
    plotDataTissue = plotDataGrouped.get_group(tissue).join(scaled.loc[:,samples]
        ).melt(id_vars = plotData.columns,var_name="Sample").merge(sampleClusters,how="left")
    plotDataTissue = plotDataTissue.groupby("clusterName")
    for index2, geneCluster in enumerate(clusters):
        plotDataTissueCluster = plotDataTissue.get_group(geneCluster)
        f = lambda cluster: plotDataTissueCluster.loc[plotDataTissueCluster.sampleCluster == cluster,"value"].to_numpy()
        dataset = [f(i) for i in plotDataTissueCluster.sampleCluster.unique()]
        ax[index1,index2].set_axis_on()
        violin = ax[index1,index2].violinplot(dataset,showextrema=False,
                                                   widths=0.85,positions=range(len(plotDataTissueCluster.sampleCluster.unique())))
        box = ax[index1,index2].boxplot(dataset,showfliers=False,widths=0.2,patch_artist=True,
                                        positions=range(len(plotDataTissueCluster.sampleCluster.unique())))
        ax[index1,index2].tick_params(labelsize=5,size=2)
        ax[index1,index2].set_title( tissue + "\n" +geneCluster, size = 5,weight="bold")
        ax[index1,index2].set_xticklabels(plotDataTissueCluster.sampleCluster.unique())
        ax[index1,index2]
        for i in violin["bodies"]:
            i.set(color=colors[plotDataTissueCluster.color],alpha=.75,linewidth=.35)
        for i,j in zip(box["boxes"],box["medians"]):
            i.set(facecolor="white",edgecolor="black",linewidth=.35)
            j.set(c="black",linewidth=0.35)
        for i,j in zip(box["whiskers"],box["caps"]):
            i.set(c="black",linewidth=.35)
            j.set(c="black",linewidth=0)
        plt.ylim([-6.5,6.5])
        ax[index1,index2].spines[["top","right"]].set_visible(False)
fig.supxlabel("Sample cluster", fontsize = 8,weight="bold")
fig.supylabel("Z-Score", fontsize = 8,weight="bold")
plt.tight_layout()
plt.subplots_adjust(top = .9, bottom = .115,right = 0.985, left = 0.09,hspace=.45,wspace=.1)
fig.savefig(f"{figurePath}CLUSTER_EXPRESSION_RENAMED.png",dpi=600,transparent=True)
plt.show()


# =============================================================================
#  generate intersection file
# =============================================================================

clusters = pd.read_csv(f"{outputPath}CLUSTERS.csv", index_col=0)

clustersClean = clusters.loc[clusters.Cluster != -1,:]
clustersClean.clusterNameShort = clustersClean.Tissue.values + " " + clustersClean.clusterNameShort.values.astype(str)
clusteredGenes = clustersClean.index.unique().to_numpy()
intersections = pd.DataFrame(data=np.zeros((len(clusteredGenes),1)),index=clusteredGenes,columns=["clusterIntersection"])


for gene in clusteredGenes:
    if clustersClean.loc[gene,:].shape == (5,):
        intersections.loc[gene,"clusterIntersection"] = clustersClean.loc[gene,"clusterNameShort"]
    else:
        intersections.loc[gene,"clusterIntersection"] = ", ".join(list(clustersClean.loc[gene,"clusterNameShort"]))

intersections.to_csv(f"{outputPath}INTERSECTIONS.csv")  

# include tissu differences


up_arrow = "↑"
down_arrow = "↓"

wald = pd.read_csv(f"{dataPath}08_DESEQ2/WALD_TISSUE.csv", index_col=0)
waldSig = wald.loc[(wald.padj < 0.001) & (np.abs(wald.log2FoldChange) > np.log2(1.5)),]
waldSig["waldGroup"] = np.nan


waldSigGenes = waldSig.index.unique().to_numpy()
waldSigGroups = pd.DataFrame(data=np.zeros((len(waldSigGenes),1)),index=waldSigGenes,columns=["waldGroup"])
waldSigGroups.loc[waldSig.log2FoldChange > 0,"waldGroup"] = f"SR{up_arrow}"
waldSigGroups.loc[waldSig.log2FoldChange < 0,"waldGroup"] = f"SR{down_arrow}"

waldClusterIntersections = intersections.join(waldSigGroups, how="outer")

waldClusteredGenes = waldClusterIntersections.index.unique().to_numpy()

waldClusterIntersections["waldClusterIntersection"] = [
    (waldClusterIntersections.loc[gene,:].waldGroup + ", " + waldClusterIntersections.loc[gene,:].clusterIntersection) 
    if (waldClusterIntersections.loc[gene,"waldGroup"] is not np.nan) & 
        (waldClusterIntersections.loc[gene,"clusterIntersection"] is not np.nan)
    else (waldClusterIntersections.loc[gene,:].clusterIntersection) 
    if waldClusterIntersections.loc[gene,"clusterIntersection"] is not np.nan
    else (waldClusterIntersections.loc[gene,:].waldGroup)
    for gene in waldClusteredGenes]

waldClusterIntersections.to_csv(f"{outputPath}WALD_CLUSTER_INTERSECTIONS.csv") 

# =============================================================================
# Marker gene expression
# =============================================================================


scaled = pd.read_csv(f"{outputPath}SCALED.csv",index_col=0)

sampleClusters = metaDataClusters.copy().reset_index()
sampleClusters["sampleCluster"] = sampleClusters.Cluster + "\n(" + sampleClusters.clusterNameShort +")"
sampleClusters = sampleClusters.drop(columns = ["Cluster","clusterName","clusterNameShort"])

waldClusterIntersections = pd.read_csv(f"{outputPath}WALD_CLUSTER_INTERSECTIONS.csv", index_col=0)

AtMarkers = pd.read_csv(f"{dataPath}00_META_DATA/AT_MARKER.csv")
blast = pd.read_csv(f"{dataPath}00_META_DATA/Mesculenta_RECIPROCAL_HIT_ARAPORT11_BLASTP.csv")

MeMarkers = AtMarkers.merge(blast).merge(waldClusterIntersections.reset_index(names="Geneid"))


MeMarkers["Symbol"] = MeMarkers.Symbol + " | " +  MeMarkers.waldClusterIntersection

scaler = StandardScaler()
scaled = pd.DataFrame(data=scaler.fit_transform(vst.T).T,
                      columns=vst.columns,index=vst.index)



plotData = scaled.copy().reset_index(names="Geneid").merge(MeMarkers)
plotData=plotData.melt(id_vars = MeMarkers.columns,var_name="Sample").merge(sampleClusters)

markers = np.array(["o","s","^"])

colors = plt.cm.plasma(np.linspace(0, 1,len(plotData.Tissue.unique())))
plotData["color"] = rankdata(plotData.Tissue,method="dense")-1
plotData["marker"] = rankdata(plotData.Plant,method="dense")-1
f1 = lambda c: matplotlib.lines.Line2D([], [], linewidth=0, color = c, marker="o",
                                      markeredgecolor = "black", markeredgewidth=.1)
handles1 = [f1(colors[i]) for i in plotData.color.unique()]
labels1 = plotData.Tissue.unique()
f2 = lambda m: matplotlib.lines.Line2D([], [], linewidth=0, color = "black", marker=m,
                                      markeredgecolor = "black", markeredgewidth=.1)
handles2 = [f2(markers[i]) for i in plotData.marker.unique()]
labels2 = plotData.Plant.unique()

cm = 1/2.54
fig, ax = plt.subplots(3,3,sharey=True)
fig.set_size_inches(18*cm,18*cm)
for roleIndex, role in enumerate(["Phloem","Cambium","Xylem"]):
    plotDataRole = plotData.loc[plotData.Role == role,:]
    
    for geneIndex, gene in enumerate(plotDataRole.Symbol.unique()):
        plotDataRoleGene = plotDataRole.loc[plotDataRole.Symbol == gene,:]
        xPos = np.array([*range(len(plotDataRoleGene.sampleCluster.unique()))])
        
        for tissue in metaData.Tissue.unique():
            plotDataRoleGeneTissue = plotDataRoleGene.loc[plotDataRoleGene.Tissue == tissue,:]
            f = lambda cluster: plotDataRoleGeneTissue.loc[plotDataRoleGeneTissue.sampleCluster == cluster,"value"].to_numpy()
            dataset = [f(i) for i in plotDataRoleGeneTissue.sampleCluster.unique()]
            
            if tissue == "SR":
                xPosNew = xPos+0.225
            else:
                xPosNew = xPos-0.225
                
            box = ax[geneIndex,roleIndex].boxplot(dataset,showfliers=False,widths=0.375,patch_artist=True,
                                            positions=xPosNew, zorder=0)
            
            for i,j in zip(box["boxes"],box["medians"]):
                i.set(facecolor=colors[plotDataRoleGeneTissue.color][0],edgecolor="black",linewidth=.5,alpha=.75)
                j.set(c="black",linewidth=0.5)
            for i,j in zip(box["whiskers"],box["caps"]):
                i.set(c="black",linewidth=.5)
                j.set(c="black",linewidth=.5)
                
            
            for markerIndex, plant in enumerate(plotDataRoleGeneTissue.Plant.unique()):
                
                plotDataRoleGeneTissuePlant = plotDataRoleGeneTissue.loc[plotDataRoleGeneTissue.Plant == plant,:]
                
                f = lambda cluster: plotDataRoleGeneTissuePlant.loc[plotDataRoleGeneTissuePlant.sampleCluster == cluster,"value"].to_numpy()
                dataset = [f(cluster) for cluster in plotDataRoleGeneTissuePlant.sampleCluster.unique()]
                
                for cluster, value in zip(xPosNew,dataset):
                    group = np.random.normal(cluster, 0.05, size=len(value))
                    scatter = ax[geneIndex,roleIndex].scatter(group,value,s=12,linewidths=.35,marker=markers[plotDataRoleGeneTissuePlant.marker.values[0]],
                                                c="white",edgecolor="black",zorder=1,alpha=0.9)

        ax[geneIndex,roleIndex].set_title(gene,fontsize=8, weight = "bold")
        ax[geneIndex,roleIndex].tick_params(labelsize=8,size=2)
        ax[geneIndex,roleIndex].set_xticks(xPos,plotData.sampleCluster.unique())
        ax[geneIndex,roleIndex].spines[["right","top"]].set_visible(False)
        ax[geneIndex,roleIndex].set_ylim(-4,4)
        ax[geneIndex,roleIndex].set_yticks(np.linspace(-3,3,5),np.linspace(-3,3,5))
legend1 = plt.legend(handles1,labels1,bbox_to_anchor = [1,3.7,0,0],
                     title="Tissue", loc = "center left",markerscale = 0.75,
                     fontsize=8,title_fontsize=9,alignment="left",handletextpad=0.1)
legend1.get_frame().set(linewidth=.75,edgecolor="black",facecolor="None")
legend2 = plt.legend(handles2,labels2,bbox_to_anchor = [1,3.25,0,0],
                     title="Plant", loc = "center left",markerscale = 0.75,
                     fontsize=8,title_fontsize=9,alignment="left",handletextpad=0.1)
legend2.get_frame().set(linewidth=.75,edgecolor="black",facecolor="None")
fig.supxlabel("Sample cluster", fontsize = 10, weight="bold")
fig.supylabel("Z-Score", fontsize = 10, weight="bold")
fig.add_artist(legend1)
fig.add_artist(legend2)
legend1.set_in_layout(False)
legend2.set_in_layout(False)
plt.tight_layout()
plt.subplots_adjust(top = .95, bottom = .1,right = 0.9, left = 0.1,hspace=.45, wspace = .1)
fig.savefig(f"{figurePath}MARKER_EXPRESSION.png",dpi=600,transparent=True)
plt.show()



# =============================================================================
# define boxplot function
# =============================================================================

           
def goiBoxPlot(figurePath,name,gene):
    
    os.makedirs(f"{figurePath}", exist_ok = True)
    cm = 1/2.54
    
    scaler = StandardScaler()
    scaled = pd.DataFrame(data=scaler.fit_transform(vst.T).T,
                          columns=vst.columns,index=vst.index)
    
    plotData = scaled.copy().loc[genes,:].reset_index(names="Geneid").merge(waldClusterIntersections)
    plotData=plotData.melt(id_vars = waldClusterIntersections.columns,var_name="Sample").merge(sampleClusters)
    
    markers = np.array(["o","s","^"])
    colors = plt.cm.plasma(np.linspace(0, 1,len(plotData.Tissue.unique())))
    
    plotData["color"] = rankdata(plotData.Tissue,method="dense")-1
    plotData["marker"] = rankdata(plotData.Plant,method="dense")-1
    
    f1 = lambda c: matplotlib.lines.Line2D([], [], linewidth=0, color = c, marker="o",
                                          markeredgecolor = "black", markeredgewidth=.1)
    handles1 = [f1(colors[i]) for i in plotData.color.unique()]
    labels1 = plotData.Tissue.unique()
    f2 = lambda m: matplotlib.lines.Line2D([], [], linewidth=0, color = "black", marker=m,
                                          markeredgecolor = "black", markeredgewidth=.1)
    handles2 = [f2(markers[i]) for i in plotData.marker.unique()]
    labels2 = plotData.Plant.unique()
    
    xPos = np.array([*range(len(plotData.sampleCluster.unique()))])
        
    geneExpression = plotData.loc[plotData.Geneid == gene,:]
    
    title = geneExpression.waldClusterIntersection.unique()[0]
    
    y_lim_low = geneExpression.value.min()*1.1
    y_lim_high = geneExpression.value.max()*1.1
    
    yBreakStart = math.floor(y_lim_low)
    yBreakStop = math.floor(y_lim_high)
    
    nPoints = int(np.ceil((yBreakStop - yBreakStart))) + 1
    
    yTickBreaks = np.linspace(yBreakStart, yBreakStop, nPoints)
    
    fig, ax = plt.subplots()
    fig.set_size_inches(8*cm,6*cm)
    
    for tissue in ["Stem","SR"]:
        geneExpressionTissue = geneExpression.loc[geneExpression.Tissue == tissue,:]
        f = lambda cluster: geneExpressionTissue.loc[geneExpressionTissue.sampleCluster == cluster,"value"].to_numpy()
        dataset = [f(i) for i in geneExpressionTissue.sampleCluster.unique()]
        if tissue == "SR":
            xPosNew = xPos-0.225
        else:
            xPosNew = xPos+0.225
        box = ax.boxplot(dataset,showfliers=False,widths=0.375,patch_artist=True,
                                        positions=xPosNew, zorder=0)
        for i,j in zip(box["boxes"],box["medians"]):
            i.set(facecolor=colors[geneExpressionTissue.color][0],edgecolor="black",linewidth=.5,alpha=.75)
            j.set(c="black",linewidth=0.5)
        for i,j in zip(box["whiskers"],box["caps"]):
            i.set(c="black",linewidth=.5)
            j.set(c="black",linewidth=.5)
            
        for plant, marker in zip(geneExpressionTissue.Plant.unique(), itertools.cycle(markers)):
            geneExpressionTissuePlant = geneExpressionTissue.loc[geneExpressionTissue. Plant == plant,:]
            f = lambda cluster: geneExpressionTissuePlant.loc[geneExpressionTissuePlant.sampleCluster == cluster,"value"].to_numpy()
            dataset = [f(i) for i in geneExpressionTissuePlant.sampleCluster.unique()]
            for i,j in zip(xPosNew,dataset):
                y = j
                x = np.random.normal(i, 0.075, size=len(y))
                ax.scatter(x,y,s=12,linewidths=.35,marker=marker,
                                            c="white",edgecolor="black",zorder=1,alpha=0.9)
            
    ax.set_title(str(title),fontsize=10, weight = "bold")
    ax.tick_params(labelsize=8,size=2)
    ax.set_xticks(xPos,plotData.sampleCluster.unique())
    ax.spines[["right","top"]].set_visible(False)
    ax.set_ylim(y_lim_low,y_lim_high)
    ax.set_yticks(yTickBreaks,yTickBreaks)
    
    legend1 = plt.legend(handles1,labels1,bbox_to_anchor = [1,.875,0,0],
                         title="Tissue", loc = "center left",markerscale = 0.75,
                         fontsize=8,title_fontsize=9,alignment="left",handletextpad=0.1)
    legend1.get_frame().set(linewidth=.75,edgecolor="black",facecolor="None")
    legend2 = plt.legend(handles2,labels2,bbox_to_anchor = [1,.475,0,0],
                         title="Plant", loc = "center left",markerscale = 0.75,
                         fontsize=8,title_fontsize=9,alignment="left",handletextpad=0.1)
    legend2.get_frame().set(linewidth=.75,edgecolor="black",facecolor="None")
    ax.set_xlabel("Sample cluster", fontsize = 10, weight="bold")
    ax.set_ylabel("Z-Score", fontsize = 10, weight="bold")
    fig.add_artist(legend1)
    fig.add_artist(legend2)
    legend1.set_in_layout(False)
    legend2.set_in_layout(False)
    plt.tight_layout()
    plt.subplots_adjust(top = .9, bottom = .225,right = 0.775, left = 0.175)
    fig.savefig(f"{figurePath}/{name}_{gene}.png",dpi=600,transparent=True)
    plt.show()


# =============================================================================
# Genes of interest expression
# =============================================================================

subprocess.run([f"{RScriptPath}", f"{scriptPath}HEATMAPS.R"])

metaDataClusters = pd.read_csv(f"{dataPath}09_SAMPLE_CLUSTERING/META_DATA_CLUSTERS.csv",index_col=0)
sampleClusters = metaDataClusters.copy().reset_index()
sampleClusters["sampleCluster"] = sampleClusters.Cluster + "\n(" + sampleClusters.clusterNameShort +")"
sampleClusters = sampleClusters.drop(columns = ["Cluster","clusterName","clusterNameShort"])

outputPath = f"{dataPath}11_GOI/"
figurePath = f"{outputPath}FIGURES/AT_GOI/BOXPLOTS/"
os.makedirs(figurePath,exist_ok=True)

waldClusterIntersections = pd.read_csv(f"{dataPath}10_CRN/WALD_CLUSTER_INTERSECTIONS.csv", index_col=0).reset_index(names="Geneid")

goi = pd.read_csv(f"{outputPath}/Me_GOI.csv").merge(waldClusterIntersections)

for ivnolvement in goi.Involvement.unique():
    
    print(ivnolvement)
    
    genes = goi.loc[goi.Involvement == ivnolvement,"Geneid"]
    names = goi.loc[goi.Involvement == ivnolvement,"Name"]
    
    for gene, name in zip(genes,names):
        
        print(f"{name}_{gene}")
    
        goiBoxPlot(f"{figurePath}{ivnolvement}", name, gene)

# add expression profiles of interest to intersection file

patternsOfInterest = waldClusterIntersections.copy() 
patternsOfInterest["Pattern"] = np.nan
patternsOfInterest.loc[patternsOfInterest.clusterIntersection == "SR P/X, Stem P", "Pattern"] = "KNOX1"
patternsOfInterest.loc[patternsOfInterest.clusterIntersection == "SR P/X", "Pattern"] = "KNOX1"

patternsOfInterest.loc[patternsOfInterest.waldClusterIntersection == "SR↓, SR X, Stem X", "Pattern"] = "MYB46"
patternsOfInterest.loc[patternsOfInterest.waldClusterIntersection == "SR↓, Stem X", "Pattern"] = "MYB46"
patternsOfInterest.loc[patternsOfInterest.waldClusterIntersection == "Stem X", "Pattern"] = "MYB46"

patternsOfInterest.loc[patternsOfInterest.clusterIntersection == "SR C/X, Stem X", "Pattern"] = "WOX14"
patternsOfInterest.loc[patternsOfInterest.clusterIntersection == "SR P/C, Stem X", "Pattern"] = "WOX14"
patternsOfInterest.loc[patternsOfInterest.clusterIntersection == "SR C/X", "Pattern"] = "WOX14"
patternsOfInterest.loc[patternsOfInterest.clusterIntersection == "SR P/C", "Pattern"] = "WOX14"

patternsOfInterest.to_csv(f"{outputPath}patternsOfInterest.csv",index=False)


# =============================================================================
# Cluster analysis and enrichment
# =============================================================================

subprocess.run([f"{RScriptPath}", f"{scriptPath}CHORDDIAGRAMS.R"])
subprocess.run([f"{RScriptPath}", f"{scriptPath}ENRICHMENT_ANALYSIS.R"])


# =============================================================================
# phylogenetic expression
# =============================================================================

for experiment in ["CRYOSECTIONING","TRANSPORT"]:
    
    print(experiment)

    outputPath = f"{dataPath}11_GOI/"
    figurePath = f"{outputPath}FIGURES/PHYLOGENY/{experiment}/"
    os.makedirs(figurePath,exist_ok=True)
    
    metaDataClusters = pd.read_csv(f"{dataPath}09_SAMPLE_CLUSTERING/META_DATA_CLUSTERS.csv",index_col=0)
    sampleClusters = metaDataClusters.copy().reset_index()
    sampleClusters["sampleCluster"] = sampleClusters.Cluster + "\n(" + sampleClusters.clusterNameShort +")"
    sampleClusters = sampleClusters.drop(columns = ["Cluster","clusterName","clusterNameShort"])
    
    
    waldClusterIntersections = pd.read_csv(f"{dataPath}10_CRN/WALD_CLUSTER_INTERSECTIONS.csv", index_col=0).reset_index(names="Geneid")
    
    phylogenies = pd.read_csv(f"{dataPath}00_META_DATA/{experiment}_NOMENCLATURE.csv").merge(waldClusterIntersections)
    
    
    for family in phylogenies.Family.unique():
        
        print(family)
        
        genes = phylogenies.loc[phylogenies.Family == family,"Geneid"]
        names = phylogenies.loc[phylogenies.Family == family,"Clade_name"]
        
        for gene, name in zip(genes,names):
            
            print(f"{name}_{gene}")
        
            goiBoxPlot(f"{figurePath}{family}", name, gene)



# =============================================================================
# phylogenetic expression
# =============================================================================

outputPath = f"{dataPath}11_GOI/"
figurePath = f"{outputPath}FIGURES/PatternsOfInterest/"
os.makedirs(figurePath,exist_ok=True)


patterns = pd.read_csv(f"{outputPath}patternsOfInterest.csv")
phylogenies = pd.read_csv(f"{dataPath}00_META_DATA/CRYOSECTIONING_NOMENCLATURE.csv").merge(patterns)
phylogenies = phylogenies.loc[(phylogenies.Pattern.notnull()),:]

metaDataClusters = pd.read_csv(f"{dataPath}09_SAMPLE_CLUSTERING/META_DATA_CLUSTERS.csv",index_col=0)
sampleClusters = metaDataClusters.copy().reset_index()
sampleClusters["sampleCluster"] = sampleClusters.Cluster + "\n(" + sampleClusters.clusterNameShort +")"
sampleClusters = sampleClusters.drop(columns = ["Cluster","clusterName","clusterNameShort"])


scaler = StandardScaler()
scaled = pd.DataFrame(data=scaler.fit_transform(vst.T).T,
                      columns=vst.columns,index=vst.index)

plotData = scaled.copy().reset_index(names="Geneid").merge(phylogenies)
plotData=plotData.melt(id_vars = phylogenies.columns,var_name="Sample").merge(sampleClusters)


markers = np.array(["o","s","^"])
colors = plt.cm.plasma(np.linspace(0, 1,len(plotData.Tissue.unique())))

plotData["color"] = rankdata(plotData.Tissue,method="dense")-1
plotData["marker"] = rankdata(plotData.Plant,method="dense")-1

f1 = lambda c: matplotlib.lines.Line2D([], [], linewidth=0, color = c, marker="o",
                                      markeredgecolor = "black", markeredgewidth=.1)
handles1 = [f1(colors[i]) for i in plotData.color.unique()]
labels1 = plotData.Tissue.unique()
f2 = lambda m: matplotlib.lines.Line2D([], [], linewidth=0, color = "black", marker=m,
                                      markeredgecolor = "black", markeredgewidth=.1)
handles2 = [f2(markers[i]) for i in plotData.marker.unique()]
labels2 = plotData.Plant.unique()

xPos = np.array([*range(len(plotData.sampleCluster.unique()))])

cm = 1/2.54

for Family in plotData.Family.unique():
    print(Family)
    FamilyPlotData = plotData.loc[plotData.Family == Family,:]
    for group in FamilyPlotData.Pattern.unique():
        groupName = group.replace(", ", "_").replace(" ","-").replace("/","")
        os.makedirs(f"{figurePath}{groupName}",exist_ok=True)
        groupPlotData = FamilyPlotData.loc[FamilyPlotData.Pattern == group,:]
        for gene in groupPlotData.Geneid.unique():
            print(gene)
            name = groupPlotData.loc[groupPlotData.Geneid == gene,"Clade_name"].unique()
            geneExpression = groupPlotData.loc[groupPlotData.Geneid == gene,:]
            title = geneExpression.Family.unique()[0] + "\n" + geneExpression.locusName.unique()[0]
            fig, ax = plt.subplots()
            fig.set_size_inches(8*cm,6*cm)
            
            for tissue in ["Stem","SR"]:
                geneExpressionTissue = geneExpression.loc[geneExpression.Tissue == tissue,:]
                f = lambda cluster: geneExpressionTissue.loc[geneExpressionTissue.sampleCluster == cluster,"value"].to_numpy()
                dataset = [f(i) for i in geneExpressionTissue.sampleCluster.unique()]
                if tissue == "SR":
                    xPosNew = xPos-0.225
                else:
                    xPosNew = xPos+0.225
                box = ax.boxplot(dataset,showfliers=False,widths=0.375,patch_artist=True,
                                                positions=xPosNew, zorder=0)
                for i,j in zip(box["boxes"],box["medians"]):
                    i.set(facecolor=colors[geneExpressionTissue.color][0],edgecolor="black",linewidth=.5,alpha=.75)
                    j.set(c="black",linewidth=0.5)
                for i,j in zip(box["whiskers"],box["caps"]):
                    i.set(c="black",linewidth=.5)
                    j.set(c="black",linewidth=.5)
                    
                for plant, marker in zip(geneExpressionTissue.Plant.unique(), itertools.cycle(markers)):
                    geneExpressionTissuePlant = geneExpressionTissue.loc[geneExpressionTissue. Plant == plant,:]
                    f = lambda cluster: geneExpressionTissuePlant.loc[geneExpressionTissuePlant.sampleCluster == cluster,"value"].to_numpy()
                    dataset = [f(i) for i in geneExpressionTissuePlant.sampleCluster.unique()]
                    for i,j in zip(xPosNew,dataset):
                        y = j
                        x = np.random.normal(i, 0.075, size=len(y))
                        scatter = ax.scatter(x,y,s=12,linewidths=.35,marker=marker,
                                                    c="white",edgecolor="black",zorder=1,alpha=0.9)
                    
            ax.set_title(str(title),fontsize=8, weight = "bold")
            ax.tick_params(labelsize=8,size=2)
            ax.set_xticks(xPos,plotData.sampleCluster.unique())
            ax.spines[["right","top"]].set_visible(False)
            ax.set_ylim(-4,4)
            ax.set_yticks(np.linspace(-3,3,5),np.linspace(-3,3,5))
            
            legend1 = plt.legend(handles1,labels1,bbox_to_anchor = [1,.875,0,0],
                                 title="Tissue", loc = "center left",markerscale = 0.75,
                                 fontsize=8,title_fontsize=9,alignment="left",handletextpad=0.1)
            legend1.get_frame().set(linewidth=.75,edgecolor="black",facecolor="None")
            legend2 = plt.legend(handles2,labels2,bbox_to_anchor = [1,.475,0,0],
                                 title="Plant", loc = "center left",markerscale = 0.75,
                                 fontsize=8,title_fontsize=9,alignment="left",handletextpad=0.1)
            legend2.get_frame().set(linewidth=.75,edgecolor="black",facecolor="None")
            ax.set_xlabel("Sample cluster", fontsize = 10, weight="bold")
            ax.set_ylabel("Z-Score", fontsize = 10, weight="bold")
            fig.add_artist(legend1)
            fig.add_artist(legend2)
            legend1.set_in_layout(False)
            legend2.set_in_layout(False)
            plt.tight_layout()
            plt.subplots_adjust(top = .875, bottom = .225,right = 0.775, left = 0.175)
            fig.savefig(f"{figurePath}{groupName}/{Family}_{gene}_{name}.png",dpi=600,transparent=True)
            plt.show()


# =============================================================================
# other goi plots 
# =============================================================================


outputPath = f"{dataPath}11_GOI/"
figurePath = f"{outputPath}FIGURES/RANDOM_GENES/"
os.makedirs(figurePath,exist_ok=True)

metaDataClusters = pd.read_csv(f"{dataPath}09_SAMPLE_CLUSTERING/META_DATA_CLUSTERS.csv",index_col=0)
sampleClusters = metaDataClusters.copy().reset_index()
sampleClusters["sampleCluster"] = sampleClusters.Cluster + "\n(" + sampleClusters.clusterNameShort +")"
sampleClusters = sampleClusters.drop(columns = ["Cluster","clusterName","clusterNameShort"])


waldClusterIntersections = pd.read_csv(f"{dataPath}10_CRN/WALD_CLUSTER_INTERSECTIONS.csv", index_col=0).reset_index(names="Geneid")



# go enriched in SR up/ SR P/X
name = "SR_Auxin"

genes = ["Manes.04G080700.v8.1","Manes.04G080800.v8.1",
         "Manes.04G081350.v8.1","Manes.04G081628.v8.1","Manes.04G081740.v8.1","Manes.04G081800.v8.1",
         "Manes.11G067901.v8.1","Manes.11G068000.v8.1","Manes.11G068100.v8.1","Manes.15G036200.v8.1"]


for gene in genes:
    
    print(f"{name}_{gene}")

    goiBoxPlot(f"{figurePath}{name}", name, gene)


# lsh enriched in bp/knat1 pattern
name = "SR_LSH"

genes = [
    "Manes.02G216200.v8.1","Manes.06G054300.v8.1","Manes.12G033200.v8.1","Manes.15G094500.v8.1",
    "Manes.17G055900.v8.1","Manes.18G142500.v8.1"]


for gene in genes:
    
    print(f"{name}_{gene}")

    goiBoxPlot(f"{figurePath}{name}", name, gene)

# =============================================================================
# # pattern expression profile
# =============================================================================


outputPath = f"{dataPath}11_GOI/"
figurePath = f"{outputPath}FIGURES/PatternsOfInterest/"
os.makedirs(figurePath,exist_ok=True)

metaDataClusters = pd.read_csv(f"{dataPath}09_SAMPLE_CLUSTERING/META_DATA_CLUSTERS.csv",index_col=0)
sampleClusters = metaDataClusters.copy().reset_index()
sampleClusters["sampleCluster"] = sampleClusters.Cluster + "\n(" + sampleClusters.clusterNameShort +")"
sampleClusters = sampleClusters.drop(columns = ["Cluster","clusterName","clusterNameShort"])

waldClusterIntersections = pd.read_csv(f"{dataPath}10_CRN/WALD_CLUSTER_INTERSECTIONS.csv", index_col=0)


MeMarkers = pd.read_csv(f"{outputPath}ME_CAMBIUM_DEVIANTS.csv", index_col=0).join(waldClusterIntersections)


MeMarkers["Symbol"] = MeMarkers.Pattern + " \n " +  MeMarkers.waldClusterIntersection

scaler = StandardScaler()
scaled = pd.DataFrame(data=scaler.fit_transform(vst.T).T,
                      columns=vst.columns,index=vst.index)



plotData = MeMarkers.join(scaled.copy()).reset_index(names="Geneid")
plotData=plotData.melt(id_vars = MeMarkers.columns,var_name="Sample").merge(sampleClusters)

markers = np.array(["o","s","^"])

colors = plt.cm.plasma(np.linspace(0, 1,len(plotData.Tissue.unique())))
plotData["color"] = rankdata(plotData.Tissue,method="dense")-1
plotData["marker"] = rankdata(plotData.Plant,method="dense")-1
f1 = lambda c: matplotlib.lines.Line2D([], [], linewidth=0, color = c, marker="o",
                                      markeredgecolor = "black", markeredgewidth=.1)
handles1 = [f1(colors[i]) for i in plotData.color.unique()]
labels1 = plotData.Tissue.unique()
f2 = lambda m: matplotlib.lines.Line2D([], [], linewidth=0, color = "black", marker=m,
                                      markeredgecolor = "black", markeredgewidth=.1)
handles2 = [f2(markers[i]) for i in plotData.marker.unique()]
labels2 = plotData.Plant.unique()

cm = 1/2.54
fig, ax = plt.subplots(2,3,sharey=True)
for i in ax: 
    for j in i: j.set_axis_off()

fig.set_size_inches(18*cm,12*cm)
for roleIndex, role in enumerate(["KNOX1","WOX14","PXL"]):
    plotDataRole = plotData.loc[plotData.Pattern == role,:]
    
    for geneIndex, gene in enumerate(plotDataRole.Symbol.unique()):
        
        ax[geneIndex,roleIndex].set_axis_on()
        
        plotDataRoleGene = plotDataRole.loc[plotDataRole.Symbol == gene,:]
        xPos = np.array([*range(len(plotDataRoleGene.sampleCluster.unique()))])
        
        for tissue in metaData.Tissue.unique():
            plotDataRoleGeneTissue = plotDataRoleGene.loc[plotDataRoleGene.Tissue == tissue,:]
            f = lambda cluster: plotDataRoleGeneTissue.loc[plotDataRoleGeneTissue.sampleCluster == cluster,"value"].to_numpy()
            dataset = [f(i) for i in plotDataRoleGeneTissue.sampleCluster.unique()]
            
            if tissue == "SR":
                xPosNew = xPos+0.225
            else:
                xPosNew = xPos-0.225
                
            box = ax[geneIndex,roleIndex].boxplot(dataset,showfliers=False,widths=0.375,patch_artist=True,
                                            positions=xPosNew, zorder=0)
            
            for i,j in zip(box["boxes"],box["medians"]):
                i.set(facecolor=colors[plotDataRoleGeneTissue.color][0],edgecolor="black",linewidth=.5,alpha=.75)
                j.set(c="black",linewidth=0.5)
            for i,j in zip(box["whiskers"],box["caps"]):
                i.set(c="black",linewidth=.5)
                j.set(c="black",linewidth=.5)
                
            
            for markerIndex, plant in enumerate(plotDataRoleGeneTissue.Plant.unique()):
                
                plotDataRoleGeneTissuePlant = plotDataRoleGeneTissue.loc[plotDataRoleGeneTissue.Plant == plant,:]
                
                f = lambda cluster: plotDataRoleGeneTissuePlant.loc[plotDataRoleGeneTissuePlant.sampleCluster == cluster,"value"].to_numpy()
                dataset = [f(cluster) for cluster in plotDataRoleGeneTissuePlant.sampleCluster.unique()]
                
                for cluster, value in zip(xPosNew,dataset):
                    group = np.random.normal(cluster, 0.05, size=len(value))
                    scatter = ax[geneIndex,roleIndex].scatter(group,value,s=12,linewidths=.35,marker=markers[plotDataRoleGeneTissuePlant.marker.values[0]],
                                                c="white",edgecolor="black",zorder=1,alpha=0.9)

        ax[geneIndex,roleIndex].set_title(gene,fontsize=8, weight = "bold")
        ax[geneIndex,roleIndex].tick_params(labelsize=8,size=2)
        ax[geneIndex,roleIndex].set_xticks(xPos,plotData.sampleCluster.unique())
        ax[geneIndex,roleIndex].spines[["right","top"]].set_visible(False)
        ax[geneIndex,roleIndex].set_ylim(-4,4)
        ax[geneIndex,roleIndex].set_yticks(np.linspace(-3,3,5),np.linspace(-3,3,5))
legend1 = plt.legend(handles1,labels1,bbox_to_anchor = [1,2.3,0,0],
                     title="Tissue", loc = "center left",markerscale = 0.75,
                     fontsize=8,title_fontsize=9,alignment="left",handletextpad=0.1)
legend1.get_frame().set(linewidth=.75,edgecolor="black",facecolor="None")
legend2 = plt.legend(handles2,labels2,bbox_to_anchor = [1,1.85,0,0],
                     title="Plant", loc = "center left",markerscale = 0.75,
                     fontsize=8,title_fontsize=9,alignment="left",handletextpad=0.1)
legend2.get_frame().set(linewidth=.75,edgecolor="black",facecolor="None")
fig.supxlabel("Sample cluster", fontsize = 10, weight="bold")
fig.supylabel("Z-Score", fontsize = 10, weight="bold")
fig.add_artist(legend1)
fig.add_artist(legend2)
legend1.set_in_layout(False)
legend2.set_in_layout(False)
plt.tight_layout()
plt.subplots_adjust(top = .925, bottom = .115,right = 0.9, left = 0.1,hspace=.45, wspace = .125)
fig.savefig(f"{figurePath}PATTERN_EXPRESSION.png",dpi=600,transparent=True)
plt.show()

# =============================================================================
# extract promoter regions
# =============================================================================

outputPath = f"{dataPath}12_PROMOTERS/"
figurePath = f"{outputPath}FIGURES"
os.makedirs(figurePath,exist_ok=True)

genomeDir = "/Users/david/Documents/GENOMES"

# first get information of total chromosome length
genome = f"{genomeDir}/Mesculenta/v8.1_w_Mt/assembly/Mesculenta_671_v8.0_MT_PT.fa"

# Initialize a dictionary to store chromosome lengths
chromosomeLengths = {}
# Iterate through the genome assembly file
with open(genome, "r") as handle:
    # iterate through fasta file (genome) as before for translation
    for record in SeqIO.parse(handle, "fasta"):
        seqname = record.id # seqname is the default heading of a gff3 file
        length = len(record.seq)
        # Store the chromosome length in the dictionary
        chromosomeLengths[seqname] = length
        

gffFile = f"{genomeDir}/Mesculenta/v8.1_w_Mt/annotation/Mesculenta_671_v8.1.gene_exons_MT_PT.gff3"
promoterRegions = []
with open(gffFile,"r") as gff:
    for line in gff:
        # Skip comment lines
        if line.startswith("#"):
            continue
        # the whole line is stored as str; split into individual fields; gff is tab (\t) separated
        fields = line.strip().split("\t") # .split() removes potential whitespaces
        # select only complete features that are gene
        if fields[2] != "gene" or len(fields) < 9 or "Manes." not in fields[8]:
            continue
        fields[3] = int(fields[3])
        fields[4] = int(fields[4])
        fields[8] = re.sub("ID=","",fields[8].strip().split(";")[0])
        if fields[6] == "+":
            fields[4] = fields[3]
            fields[3] = max(0,fields[3]-1000)
        else:
            fields[3] = fields[4]
            fields[4] = min(chromosomeLengths[fields[0]],fields[4]+2000)
        promoterRegions.append(fields)
promoterRegions = np.array(promoterRegions)[:,[0,3,4,8,5,6]]

promoters = f"{outputPath}PROMOTER_REGIONS.bed"
np.savetxt(f"{promoters}", promoterRegions, delimiter='\t', fmt='%s')


command = f"/opt/homebrew/bin/bedtools getfasta -nameOnly -s -fi {genome} -bed {promoters} -fo {outputPath}PROMOTERS.fa"
subprocess.call(command,shell=True)

# https://meme-suite.org/meme/doc/fimo-tutorial.html?man_type=web

# fimo [options] <motif file> <sequence file>

command = f"/usr/local/bin/fimo --o {outputPath}FIMO --no-pgc --thresh 2.5e-05 {outputPath}Mes_TF_binding_motifs.meme {outputPath}PROMOTERS.fa"
subprocess.call(command,shell=True)

# run enrichment

subprocess.run([f"{RScriptPath}", f"{scriptPath}TF_ENRICHMENT.R"])

# =============================================================================
# plot allf TFs
# =============================================================================

outputPath = f"{dataPath}12_PROMOTERS/"
figurePath = f"{outputPath}FIGURES/TF_EXPRESSION/"
os.makedirs(figurePath,exist_ok=True)

patterns = pd.read_csv(f"{dataPath}11_GOI/patternsOfInterest.csv")
TFs = pd.read_csv(f"{outputPath}TF_LIST_CLEAN.csv")
TFs["Geneid"] = TFs["TF"]
TFs = TFs.merge(patterns,how="left")
TFs = TFs.loc[(TFs.Pattern.notnull()),:]

metaDataClusters = pd.read_csv(f"{dataPath}09_SAMPLE_CLUSTERING/META_DATA_CLUSTERS.csv",index_col=0)
sampleClusters = metaDataClusters.copy().reset_index()
sampleClusters["sampleCluster"] = sampleClusters.Cluster + "\n(" + sampleClusters.clusterNameShort +")"
sampleClusters = sampleClusters.drop(columns = ["Cluster","clusterName","clusterNameShort"])


for pattern in TFs.Pattern.unique():
    
    print(pattern)
    
    PatternTFs = TFs.loc[TFs.Pattern == pattern,:]
    
    for family in PatternTFs.Family.unique():
        
        print(family)
        
        genes = PatternTFs.loc[PatternTFs.Family == family,"Geneid"]
        
        for gene in genes:
            
            print(f"{family}_{gene}")
        
            goiBoxPlot(f"{figurePath}{pattern}", family, gene)
