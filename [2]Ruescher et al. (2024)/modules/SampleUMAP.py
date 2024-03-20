import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from umap.umap_ import UMAP
class SampleUMAP:
    '''
    A class for performing UMAP projection on RNA-seq data
    '''
    def __init__(self,logCounts,normCounts,metaData,groupVar=None,group=None):
        '''
        Initialize the SampleCommunitiesPCA class.

        Parameters:
            - logCounts: DataFrame, log-transformed gene expression data.
            - normCounts: DataFrame, normalized gene expression data.
            - metaData: DataFrame, metadata for samples.
            - groupVar: str, variable for grouping samples (optional).
            - group: str, specific group of samples to analyze (required, if groupVar != None).

        Raises:
            - ValueError: metaData, normCount and logCount Sample IDs (index for metaData, columns for rest) have to match,
                otherwise the results will not make sense for visualization, etc; a common problem for beginners.
        '''
        if not logCounts.columns.equals(metaData.index) or not normCounts.columns.equals(metaData.index):
            raise ValueError("CountMatrix and MetaData do not match. They should have the same sample IDs.")
        self.logCounts = logCounts
        self.normCounts = normCounts
        self.group = group
        self.groupVar = groupVar
        if groupVar != None:
            if group == None:
                raise RuntimeError("No group is defined. Please set groupVar = None, or define group")
            self.Samples = metaData.loc[metaData[groupVar] == group,:].index.array
        else:self.Samples = metaData.index.array
        self.metaData = metaData.loc[self.Samples,:]
    def scale(self,with_std=True):
        '''
        Required scaling of the data using sklearns StandardScaler

        Parameters:
            - with_sd: True/False, gives user choice to center or scale the data

        Note: scaling is highly encouraged as it usually improves results.

        '''
        scaler = StandardScaler(with_mean=True,with_std=with_std)
        self.Scaled = scaler.fit_transform(self.logCounts.T.loc[self.Samples,:])
        return self
    def depthFilter(self,thrsh=50):
        '''
        A required function that selects features (genes) based on normalized read depth.

        Filtering for depth of RNA-seq data is very important, expression of 0 is most common and
            low count genes tend to have extremely high variance and worsen analysis.

        Parameters:
        - thrsh: int, threshold for depth filtering, default 50 is a good heuristic

        Returns:
            - DepthFeatures: array of int, Indices of genes with mean expression > thrsh
        '''
        tmp=self.normCounts.loc[:,self.Samples].mean(axis=1)>thrsh
        self.DepthFeatures = self.normCounts.reset_index().index.array[tmp]
        self.DepthGenes = self.normCounts.index.array[tmp]
        return self
    
    def umap(self,nNeighbors=15,minDist=0.01):
        '''
        Additional function that executes UMAP projection for visual inspection.
        Does not require PCA or louvain to be executed.

        Parameters:
            - nNeighbors: int, see umap-learn n_neighbors
            - minDist: int, see umap-learn min_dist

        Raises:
            - RuntimeError: Check, if scale and depthFilter were executed beforehand

        Returns:
            - umap-learn object
        '''
        if self.Scaled is None:
            raise RuntimeError("Scaled data not found. Did you forget to run .scale()?")
        if self.DepthFeatures is None:
            raise RuntimeError("DepthFeatures not found. Did you forget to run .depthFilter()?")
        reducer = UMAP(n_neighbors=nNeighbors,min_dist=minDist)
        self.Embedding = reducer.fit(self.Scaled[:,self.DepthFeatures])
        return self
    def extractData(self,UMAP=True):
        '''
        Extracts additionally generated data as data frame for saving, visualizations, etc.

        Parameters:
            - PCA: True/False, define, if PCA data shall be extracted
            - VarPCA: True/False, define, if explained PCA variance, etc. shall be extracted
            - UMAP: True/False, define, if UMAP data shall be extracted
        '''
        self.Export = {}
        if UMAP:
            self.Export["UMAP"] = pd.DataFrame(
                data = self.Embedding.embedding_,
                index = self.Samples,
                columns = ("UMAP" + pd.Series(np.arange(self.Embedding.n_components)+1).astype(str)))