#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: David Rüscher (github.com/DavidRuescher95)

A Python script that runs the entire analyses used for the paper.

For your own use of the RNA-seq data: 
    Change the working directory (wd) to the path to this script and put the raw 
    RNA-seq data into "../data/RNASEQ/01_RAW_DATA" and the genome and meta data
    into "../data/RNASEQ/00_META_DATA". 
    Use the default folder structure from phytozome (Mesculenta/v8.1/ ...). Also specify
    the full path to your RScript binary.
    
    The META_DATA.csv file will be produced from the biosample attribute in the DESEQ2.R script.
    Save these as BIOSAMPLE_ATTRIBUTES.csv into the meta data folder. It should contain columns
    named Sample, Tissue and Stage.
    SoL, Bark (Stem_LP), Wood (Stem_LC), SR and FR were analyzed. The other available tissues
    were processed but did not add any useful information and were removed for conciseness' reasons
    from the manuscript only.
    
    For the information in Figure 9, data from a second experiment was added. The META_DATA_BULKING.csv
    file contained the corresponding Sample, Experiment, Tissue and Stage meta data.
    It uses the SR and FR samples from the main experiment plus samples starting with "A" instead of "T".
    
This script will only produce rudimentary graphs. The ones in the paper were created
    in ggplot2 (R) and aggregated in Sketch (macOS) for purely asthetic reasons.
    
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
import pandas as pd
import numpy as np
from scipy.stats import rankdata
import matplotlib
import matplotlib.pyplot as plt

# =============================================================================
# variables 
# =============================================================================
# define working directory
wd = "/Users/david/Documents/GitHub/Ruescher_et_al_2024/"
# define data path
dataPath = "../data/RNASEQ/"
# define path to scripts
scriptPath = f"{wd}scripts/"
# define RScript location
RScriptPath = "/usr/local/bin/Rscript"

# =============================================================================
#  change working directory
# =============================================================================
os.chdir(path = wd)

# =============================================================================
# load custom classes
# =============================================================================
from modules.SampleUMAP import SampleUMAP
from modules.CoRegulatoryNetworkAnalysis import CRN

# =============================================================================
# run pre-processing pipeline
# =============================================================================
# QC, trimming, QC, mapping, quantification
subprocess.call(f"{scriptPath}RNASEQ/DATA_PREPROCESSING.sh", shell = True)

# =============================================================================
# Run DESeq2 Rscript
# =============================================================================
# run DESeq2 for normalization, vst and DE Analysis
subprocess.run([f"{RScriptPath}", f"{scriptPath}RNASEQ/DESEQ2.R"])

# =============================================================================
#  load data
# ============================================================================
metaData = pd.read_csv(
    filepath_or_buffer=f"{dataPath}00_META_DATA/META_DATA.csv",
    index_col=0
    )
vst = pd.read_csv(
    filepath_or_buffer=f"{dataPath}08_DESEQ2/VST.csv",
    index_col=0
    )
norm=pd.read_csv(
    filepath_or_buffer=f"{dataPath}08_DESEQ2/NORMALIZED.csv",
    index_col=0
    )

# =============================================================================
# UMAP
# =============================================================================

figurePath = f"{dataPath}09_SAMPLE_UMAP/FIGURES/"
outputPath = f"{dataPath}09_SAMPLE_UMAP/"
os.makedirs(figurePath,exist_ok=True)

# keep in mind that UMAP results will differ slightly between runs
# class is described in ./modules/SampleUMAP
umap = SampleUMAP(vst,norm,metaData)
umap.scale()
umap.depthFilter(50)
umap.umap(30,0.01)
umap.extractData()
# join with meta data and save file
umapRes = umap.Export["UMAP"]
umapRes = metaData.join(umapRes)
umapRes.to_csv(f"{outputPath}UMAP.csv")

# =============================================================================
# sample umap plot
# =============================================================================

plotData = umapRes.copy()
markers = ["o","s","^"]
colors = plt.cm.plasma(np.linspace(0,1,len(umapRes.Tissue.unique())))
plotData["color"] = [(np.where(plotData.Tissue.unique() == i)[0][0]) for i in plotData.Tissue]
plotData["marker"] = [(np.where(plotData.Stage.unique() == i)[0][0]) for i in plotData.Stage]
f1 = lambda c: matplotlib.lines.Line2D([], [], linewidth=0, color = c, marker="o",
                                      markeredgecolor = "black", markeredgewidth=.1)
handles1 = [f1(colors[i]) for i in plotData.color.unique()]
labels1 = plotData.Tissue.unique()
f2 = lambda m: matplotlib.lines.Line2D([], [], linewidth=0, color = "black", marker=m,
                                      markeredgecolor = "black", markeredgewidth=.1)
handles2 = [f2(markers[i]) for i in plotData.marker.unique()]
labels2 = plotData.Stage.unique()
plt.plasma()

cm = 1/2.54
fig, ax = plt.subplots()
fig.set_size_inches(8*cm,6*cm)

for markerIndex, marker in enumerate(plotData.marker):

    ax.scatter(x = plotData.UMAP1[markerIndex],y = plotData.UMAP2[markerIndex],
                  c = colors[plotData.color[markerIndex].astype(int)],s=12,
                  marker=markers[plotData.marker[markerIndex]],
                  linewidths=.1,edgecolor="black",alpha=0.75)

    ax.tick_params(labelsize=6,size=2)
    ax.spines[["top","right"]].set_visible(False)

fig.supxlabel("UMAP1", fontsize = 8, weight = "bold")
fig.supylabel("UMAP2", fontsize = 8, weight = "bold")
legend1 = plt.legend(handles1,labels1,bbox_to_anchor = [1,.8,0,0],
                     title="Tissue", loc = "center left",markerscale = 0.75,
                     fontsize=6,title_fontsize=7,alignment="left",handletextpad=0.1)
legend1.get_frame().set(linewidth=.75,edgecolor="black",facecolor="None")
legend2 = plt.legend(handles2,labels2,bbox_to_anchor = [1,.425,0,0],
                     title="Stage", loc = "center left",markerscale = 0.75,
                     fontsize=6,title_fontsize=7,alignment="left",handletextpad=0.1)
legend2.get_frame().set(linewidth=.75,edgecolor="black",facecolor="None")
fig.add_artist(legend1)
fig.add_artist(legend2)
legend1.set_in_layout(False)
legend2.set_in_layout(False)

plt.tight_layout()
plt.subplots_adjust(top = .925, bottom = .16,right = 0.8, left = 0.15,hspace=0.35)
fig.savefig(f"{figurePath}UMAP.png",dpi=600,transparent=True)
plt.show()

# =============================================================================
# cluster analysis
# =============================================================================

figurePath = f"{dataPath}10_CRN/FIGURES/"
outputPath = f"{dataPath}10_CRN/"
os.makedirs(figurePath,exist_ok=True)


#  load lrt data
lrt = pd.read_csv(filepath_or_buffer=f"{dataPath}08_DESEQ2/LRT.csv",index_col=0)

clusterRes = []
umapRes = []
scaledTissue = []
corrThrsh = {}

for tissue in metaData.Tissue.unique():

    print(tissue)    
    # execute analysuis
    crn = CRN(vst,norm,metaData,lrt,"Tissue",tissue)
    crn.scale()
    crn.pValFilter(0.001)
    crn.corrNetwork(nPerm=1000,filterType="pValue",corrThrsh="permutation",quantileCutoff=0.99995)
    crn.networkGeneration()
    crn.louvain(50,0.75)
    crn.umap(30,0.1)
    crn.extractData(Network=False)
    # save permutation threshold (saves time for reruns, if needed)
    corrThrsh[tissue] = crn.CorrThrsh
    # append cluster to list
    clusters = crn.Clusters
    clusters["Tissue"] = tissue
    clusterRes.append(clusters)
    # append umap to list
    umap = crn.Export["UMAP"]
    umap["Tissue"] = tissue
    umapRes.append(umap)
    # append individually scaled expression to list
    scaledTissue.append(crn.Export["Scaled"])
    # remove class to save memory
    del crn
    
# combine and save data
clusterRes = pd.concat(clusterRes)
clusterRes.to_csv(f"{outputPath}CLUSTERS.csv")

umapRes = pd.concat(umapRes)
umapRes = clusterRes.reset_index().merge(umapRes.reset_index()).set_index("index")
umapRes.to_csv(f"{outputPath}UMAP.csv")

scaledTissue = pd.concat(scaledTissue, axis = 1)
scaledTissue.to_csv(f"{outputPath}SCALED.csv")

# cluster intersections 
# was used for comprehensive lists
clusters = clusterRes.copy()
clusters.Cluster = clusters.Tissue.values + " " + clusters.Cluster.values.astype(str)
clusteredGenes = clusters.index.unique().to_numpy()
intersections = pd.DataFrame(data=np.zeros((len(clusteredGenes),1)),index=clusteredGenes,columns=["clusterIntersection"])
i = 0
for gene in clusteredGenes:
    print(i)
    if clusters.loc[gene,:].shape == (3,):
        intersections.loc[gene,"clusterIntersection"] = clusters.loc[gene,"Cluster"]
    else:
        intersections.loc[gene,"clusterIntersection"] = ", ".join(list(clusters.loc[gene,"Cluster"]))
    i = i+1
intersections.to_csv(f"{outputPath}INTERSECTIONS.csv")  


# gene umap plot

umapRes.Cluster=umapRes.Cluster.astype(str)
plotData = umapRes.copy()
colors = plt.cm.plasma(np.linspace(0,1,len(umapRes.Cluster.unique())))
plotData["color"] = rankdata(plotData.Cluster,method="dense")-1
f = lambda c: matplotlib.lines.Line2D([], [], linewidth=0, color = c, marker="o",
                                      markeredgecolor = "black", markeredgewidth=.1)
handles = [f(colors[i]) for i in plotData.color.unique()]
labels = plotData.Cluster.unique()

plotDataGrouped = plotData.groupby(["Tissue"])

cm = 1/2.54
fig, ax = plt.subplots(2,3)

# turn axis off, until used; removes empty subplots
for i in ax: 
    for j in i: j.set_axis_off()

fig.set_size_inches(10*cm,10*cm)

for (index, tissue) in enumerate(plotData.Tissue.unique()):
    print(tissue)
    plotDataTissue = plotDataGrouped.get_group(tissue)
    
    if index <= 2:
        index1 = 0
        index2 = index
    else:
        index1 = 1
        index2 = index-3
        

    ax[index1,index2].set_axis_on()
    ax[index1,index2].scatter(x = plotDataTissue.UMAP1,y = plotDataTissue.UMAP2,
                  c = colors[plotDataTissue.color.astype(int)],marker="o",s=2,
                  linewidths=0,edgecolor="black",alpha=.1)

    ax[index1,index2].tick_params(labelsize=6,size=2)
    ax[index1,index2].spines[["top","right"]].set_visible(False)
    ax[index1,index2].set_title(tissue, size = 6,weight="bold")


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
fig.savefig(f"{figurePath}CLUSTER_UMAP.png",dpi=600,transparent=True)
plt.show()


# cluster expression plot 

metaDataPlot = metaData.reset_index()

plotData = umapRes.copy()
plotData.Cluster=plotData.Cluster.astype(str)
plotData["color"] = rankdata(plotData.Cluster,method="dense")-1
plotDataGrouped = plotData.groupby("Tissue")

cm = 1/2.54
fig, ax = plt.subplots(len(plotData.Tissue.unique()),len(plotData.Cluster.unique()), sharey = True)#
# turn off axis until used
for i in ax: 
    for j in i: j.set_axis_off()
    
fig.set_size_inches(10*cm,10*cm)

for index1, tissue in enumerate(plotData.Tissue.unique()): 
    
    samples = metaDataPlot.loc[metaDataPlot.Tissue == tissue,"Sample"]
    plotDataTissue = plotDataGrouped.get_group(tissue).join(scaledTissue.loc[:,samples]
        ).melt(id_vars = plotData.columns,var_name="Sample").merge(metaDataPlot,how="left")
    clusters = plotDataTissue.Cluster.unique()
    plotDataTissue = plotDataTissue.groupby("Cluster")
    
    for index2, geneCluster in enumerate(clusters):
        plotDataTissueCluster = plotDataTissue.get_group(geneCluster)
        f = lambda cluster: plotDataTissueCluster.loc[plotDataTissueCluster.Stage == cluster,"value"].to_numpy()
        dataset = [f(i) for i in plotDataTissueCluster.Stage.unique()]
        
        ax[index1,index2].set_axis_on()
        
        violinPlot = ax[index1,index2].violinplot(dataset,showextrema=False,
                                                   widths=0.85,positions=range(len(plotDataTissueCluster.Stage.unique())))
        
        boxPlot = ax[index1,index2].boxplot(dataset,showfliers=False,widths=0.2,patch_artist=True,
                                        positions=range(len(plotDataTissueCluster.Stage.unique())))
        
        ax[index1,index2].tick_params(labelsize=6,size=2)
        ax[index1,index2].set_title( tissue + " | "+geneCluster, size = 6,weight="bold")
        ax[index1,index2].set_xticklabels(plotDataTissueCluster.Stage.unique())
        ax[index1,index2].spines[["top","right"]].set_visible(False)
        
        for body in violinPlot["bodies"]:
            body.set(color=colors[plotDataTissueCluster.color],alpha=.75,linewidth=.35)
            
        for box,boxMedian in zip(boxPlot["boxes"],boxPlot["medians"]):
            box.set(facecolor="white",edgecolor="black",linewidth=.35)
            boxMedian.set(c="black",linewidth=0.35)
            
        for whisker,cap in zip(boxPlot["whiskers"],boxPlot["caps"]):
            whisker.set(c="black",linewidth=.35)
            cap.set(c="black",linewidth=0)
            
        plt.ylim([-3.5,3.5])
        
fig.supxlabel("Stage", fontsize = 8,weight="bold")
fig.supylabel("Z-Score", fontsize = 8,weight="bold")
legend = plt.legend(handles,labels,bbox_to_anchor = [1,8.25,0,0],
                     title="Cluster", loc = "center left",markerscale = 0.75,
                     fontsize=6,title_fontsize=7,alignment="left",handletextpad=0.1)
legend.get_frame().set(linewidth = 0.75, edgecolor="black",facecolor="None")
fig.add_artist(legend)
legend.set_in_layout(False)
plt.tight_layout()
plt.subplots_adjust(top = .95, bottom = .1,right = 0.85, left = 0.15,hspace=.95,wspace=.15)
fig.savefig(f"{figurePath}CLUSTER_EXPRESSION.png",dpi=600,transparent=True)
plt.show()

# =============================================================================
# Run cluster functional analysis Rscript
# =============================================================================
# run DESeq2 for normalization, vst and DE Analysis
subprocess.run([f"{RScriptPath}", f"{scriptPath}RNASEQ/CLUSTER_ANALYSIS.R"])


# =============================================================================
# Second bulking experiment
# =============================================================================

# define data path
dataPath = "../data/RNASEQ_BULKING/"

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
# # sample umap 
# =============================================================================


outputPath = f"{dataPath}09_SAMPLE_UMAP/"
figurePath = f"{outputPath}FIGURES/"
os.makedirs(figurePath,exist_ok=True)


def umapPlot(plotData, colorValues, markerValues, figurePath):
    
    markers = ["o","s","^","P","p","v"]
    colors = plt.cm.plasma(np.linspace(0,1,len(plotData[colorValues].unique())))
    plotData["color"] = rankdata(plotData[colorValues],method="dense")-1
    plotData["marker"] = rankdata(plotData[markerValues],method="dense")-1
    
    plotData["color"] = [(np.where(plotData[colorValues].unique() == i)[0][0]) for i in plotData[colorValues]]
    plotData["marker"] = [(np.where(plotData[markerValues].unique() == i)[0][0]) for i in plotData[markerValues]]
    
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
                          c = colors[plotDataTissue.color[markerIndex].astype(int)],s=16,
                          marker=markers[plotDataTissue.marker[markerIndex]],
                          linewidths=.1,edgecolor="black")
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
    
     
# =============================================================================
# Raw data UMAP
# =============================================================================

umapResRaw = []

for tissue in metaData.Tissue.unique():
    print(tissue)
    
    umap = SampleUMAP(vst,norm,metaData,"Tissue",tissue)
    umap.scale()
    umap.depthFilter(50)
    umap.umap(15,0.01)
    umap.extractData()
    # join with meta data and save file
    umapResRaw.append(umap.Export["UMAP"])
    
umapResRaw = metaData.join(pd.concat(umapResRaw))
umapResRaw.to_csv(f"{dataPath}UMAP.csv")


plotData = umapResRaw.copy()    
    
umapPlot(plotData, "Stage", "Experiment", f"{figurePath}UMAP_vstRaw.png")


# =============================================================================
# Corrected data UMAP
# =============================================================================


umapRes = []

for tissue in metaData.Tissue.unique():
    print(tissue)
    
    umap = SampleUMAP(vst,norm,metaData,"Tissue",tissue)
    umap.scale()
    umap.depthFilter(50)
    umap.umap(15,0.01)
    umap.extractData()
    # join with meta data and save file
    umapRes.append(umap.Export["UMAP"])
    
umapRes = metaData.join(pd.concat(umapRes))
umapRes.to_csv(f"{dataPath}UMAP.csv")

plotData = umapRes.copy()
umapPlot(plotData, "Stage", "Experiment", f"{figurePath}UMAP.png")

# =====================================
# CRN analysis
# =============================================================================


figurePath = f"{dataPath}10_CRN/FIGURES/"
outputPath = f"{dataPath}10_CRN/"
os.makedirs(figurePath,exist_ok=True)


#  load lrt data
lrt = pd.read_csv(filepath_or_buffer=f"{dataPath}08_DESEQ2/LRT.csv",index_col=0)

clusterRes = []
umapRes = []
scaledTissue = []
corrThrsh = {}

for tissue in metaData.Tissue.unique():

    print(tissue)    
    # execute analysuis
    crn = CRN(vst,norm,metaData,lrt,"Tissue",tissue)
    crn.scale()
    crn.pValFilter(0.001)
    crn.corrNetwork(nPerm=1,filterType="pValue",corrThrsh="permutation",quantileCutoff=0.99995)
    crn.networkGeneration()
    crn.louvain(50,0.75)
    crn.umap(30,0.1)
    crn.extractData(Network=False)
    # save permutation threshold (saves time for reruns, if needed)
    corrThrsh[tissue] = crn.CorrThrsh
    # append cluster to list
    clusters = crn.Clusters
    clusters["Tissue"] = tissue
    clusterRes.append(clusters)
    # append umap to list
    umap = crn.Export["UMAP"]
    umap["Tissue"] = tissue
    umapRes.append(umap)
    # append individually scaled expression to list
    scaledTissue.append(crn.Export["Scaled"])
    # remove class to save memory
    del crn
    
# combine and save data
clusterRes = pd.concat(clusterRes)
clusterRes.to_csv(f"{outputPath}CLUSTERS.csv")

umapRes = pd.concat(umapRes)
umapRes = clusterRes.reset_index().merge(umapRes.reset_index()).set_index("index")
umapRes.to_csv(f"{outputPath}UMAP.csv")

scaledTissue = pd.concat(scaledTissue, axis = 1)
scaledTissue.to_csv(f"{outputPath}SCALED.csv")

# umap plot

umapRes.Cluster=umapRes.Cluster.astype(str)
plotData = umapRes.copy()
colors = plt.cm.plasma(np.linspace(0,1,len(umapRes.Cluster.unique())))
plotData["color"] = rankdata(plotData.Cluster,method="dense")-1
f = lambda c: matplotlib.lines.Line2D([], [], linewidth=0, color = c, marker="o",
                                      markeredgecolor = "black", markeredgewidth=.1)
handles = [f(colors[i]) for i in plotData.color.unique()]
labels = plotData.Cluster.unique()

plotDataGrouped = plotData.groupby(["Tissue"])

cm = 1/2.54
fig, ax = plt.subplots(2,3)

# turn axis off, until used; removes empty subplots
for i in ax: 
    for j in i: j.set_axis_off()

fig.set_size_inches(10*cm,10*cm)

for (index, tissue) in enumerate(plotData.Tissue.unique()):
    print(tissue)
    plotDataTissue = plotDataGrouped.get_group(tissue)
    
    if index <= 2:
        index1 = 0
        index2 = index
    else:
        index1 = 1
        index2 = index-3
        

    ax[index1,index2].set_axis_on()
    ax[index1,index2].scatter(x = plotDataTissue.UMAP1,y = plotDataTissue.UMAP2,
                  c = colors[plotDataTissue.color.astype(int)],marker="o",s=2,
                  linewidths=0,edgecolor="black",alpha=.1)

    ax[index1,index2].tick_params(labelsize=6,size=2)
    ax[index1,index2].spines[["top","right"]].set_visible(False)
    ax[index1,index2].set_title(tissue, size = 6,weight="bold")


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
fig.savefig(f"{figurePath}CLUSTER_UMAP.png",dpi=600,transparent=True)
plt.show()


# cluster expression plot 

metaDataPlot = metaData.reset_index()

plotData = umapRes.copy()
plotData.Cluster=plotData.Cluster.astype(str)
plotData["color"] = rankdata(plotData.Cluster,method="dense")-1
plotDataGrouped = plotData.groupby("Tissue")

cm = 1/2.54
fig, ax = plt.subplots(len(plotData.Tissue.unique()),len(plotData.Cluster.unique()), sharey = True)#
# turn off axis until used
for i in ax: 
    for j in i: j.set_axis_off()
    
fig.set_size_inches(10*cm,10*cm)

for index1, tissue in enumerate(plotData.Tissue.unique()): 
    
    samples = metaDataPlot.loc[metaDataPlot.Tissue == tissue,"Sample"]
    plotDataTissue = plotDataGrouped.get_group(tissue).join(scaledTissue.loc[:,samples]
        ).melt(id_vars = plotData.columns,var_name="Sample").merge(metaDataPlot,how="left")
    clusters = plotDataTissue.Cluster.unique()
    plotDataTissue = plotDataTissue.groupby("Cluster")
    
    for index2, geneCluster in enumerate(clusters):
        plotDataTissueCluster = plotDataTissue.get_group(geneCluster)
        f = lambda cluster: plotDataTissueCluster.loc[plotDataTissueCluster.Stage == cluster,"value"].to_numpy()
        dataset = [f(i) for i in plotDataTissueCluster.Stage.unique()]
        
        ax[index1,index2].set_axis_on()
        
        violinPlot = ax[index1,index2].violinplot(dataset,showextrema=False,
                                                   widths=0.85,positions=range(len(plotDataTissueCluster.Stage.unique())))
        
        boxPlot = ax[index1,index2].boxplot(dataset,showfliers=False,widths=0.2,patch_artist=True,
                                        positions=range(len(plotDataTissueCluster.Stage.unique())))
        
        ax[index1,index2].tick_params(labelsize=6,size=2)
        ax[index1,index2].set_title( tissue + " | "+geneCluster, size = 6,weight="bold")
        ax[index1,index2].set_xticklabels(plotDataTissueCluster.Stage.unique())
        ax[index1,index2].spines[["top","right"]].set_visible(False)
        
        for body in violinPlot["bodies"]:
            body.set(color=colors[plotDataTissueCluster.color],alpha=.75,linewidth=.35)
            
        for box,boxMedian in zip(boxPlot["boxes"],boxPlot["medians"]):
            box.set(facecolor="white",edgecolor="black",linewidth=.35)
            boxMedian.set(c="black",linewidth=0.35)
            
        for whisker,cap in zip(boxPlot["whiskers"],boxPlot["caps"]):
            whisker.set(c="black",linewidth=.35)
            cap.set(c="black",linewidth=0)
            
        plt.ylim([-3.5,3.5])
        
fig.supxlabel("Stage", fontsize = 8,weight="bold")
fig.supylabel("Z-Score", fontsize = 8,weight="bold")
legend = plt.legend(handles,labels,bbox_to_anchor = [1,8.25,0,0],
                     title="Cluster", loc = "center left",markerscale = 0.75,
                     fontsize=6,title_fontsize=7,alignment="left",handletextpad=0.1)
legend.get_frame().set(linewidth = 0.75, edgecolor="black",facecolor="None")
fig.add_artist(legend)
legend.set_in_layout(False)
plt.tight_layout()
plt.subplots_adjust(top = .95, bottom = .1,right = 0.85, left = 0.15,hspace=.95,wspace=.15)
fig.savefig(f"{figurePath}CLUSTER_EXPRESSION.png",dpi=600,transparent=True)
plt.show()



# =============================================================================
# Run analyses of metabolite data
# =============================================================================
# download data from supplementary material
subprocess.run([f"{RScriptPath}", f"{scriptPath}METABOLITES/NSC.R"])
subprocess.run([f"{RScriptPath}", f"{scriptPath}METABOLITES/SOL_SUGAR_CONCENTRATION.R"])
subprocess.run([f"{RScriptPath}", f"{scriptPath}METABOLITES/NAF.R"])
subprocess.run([f"{RScriptPath}", f"{scriptPath}METABOLITES/PHLOEM_EXUDATES.R"])
subprocess.run([f"{RScriptPath}", f"{scriptPath}METABOLITES/METABOLITE_PROFILES.R"])


# =============================================================================
# Run phylogenetic analyses
# =============================================================================
# Download CDS and peptide data of the from phytozome for P. trichocarpa, A thaliana and
#   M esculenta. Do a PFAM scan and BLASTP analysis against araport11 proteome
#   and save all to ../data/PHYLOGENIES/00_META_DATA
#   Call the files Athaliana_*, Ptrichocarpa_* and Mesculenta_* for each.
subprocess.run([f"{RScriptPath}", f"{scriptPath}PHYLOGENIES/01_GENERATE_INPUT.R"])
subprocess.call(f"{scriptPath}PHYLOGENIES/02_CDS_MAFFT.sh", shell=True)
subprocess.run([f"{RScriptPath}", f"{scriptPath}PHYLOGENIES/03_FILTER_CDS_MSA.R"])
subprocess.call(f"{scriptPath}PHYLOGENIES/04_CDS_FILTERED_IQTREE.sh", shell=True)
