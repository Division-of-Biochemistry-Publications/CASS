import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
import networkx as nx
import networkx.algorithms.community as nx_comm
from umap.umap_ import UMAP
class CRN:
    '''
    A class for performing Community detection on RNA-seq data based on a graph adjacency matrix created
    through pearson correlation
    '''
    def __init__(self,logCounts,normCounts,metaData,LRT=None,groupVar=None,group=None):
        '''
        Initialize the CRN class.

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
        else:
            self.logCounts = logCounts
            self.normCounts = normCounts
        self.group = group
        self.groupVar = groupVar
        if groupVar != None:
            if group == None:
                raise RuntimeError("No group is defined. Please set groupVar = None, or define group")
            self.Samples = metaData.loc[metaData[groupVar] == group,:].index.array
            if LRT is not None:
                self.lrt = LRT.loc[LRT[groupVar]==group,:]
        else:
            self.Samples = metaData.index.array
            if LRT is not  None:
                self.lrt = LRT
        self.metaData = metaData.loc[self.Samples,:]
    def scale(self,with_mean=True,with_std=True):        
        scaler = StandardScaler(with_std=with_std)        
        self.Scaled = scaler.fit_transform(self.logCounts.T.loc[self.Samples,:]).T       
        return self
        '''
        Required scaling of the data using sklearns StandardScaler

        Parameters:
            - with_sd: True/False, gives user choice to center or scale the data

        Note: scaling is highly encouraged as it usually improves results.

        '''
    def depthFilter(self,thrsh=50):
        '''
        Required function, if corrNetwork shall be depth based.
        Selects features (genes) based on normalized read depth.

        Filtering for depth of RNA-seq data is very important, expression of 0 is most common and
            low count genes tend to have extremely high variance and worsen analysis.

        Parameters:
        - thrsh: int, threshold for depth filtering, default 50 is a good heuristic

        Returns:
            - DepthFeatures: array of int, Indices of genes with mean expression > thrsh
        '''
        tmp=self.normCounts.loc[:,self.Samples].mean(axis=1)>thrsh
        self.DepthFeatures = self.normCounts.reset_index().index.array[tmp]
        return self
    def pValFilter(self,thrsh=0.001):
        '''
        Required function, if corrNetwork shall be pValue based.
        Selects features (genes) based on adjusted pValues reported by DESeq2 LRT.

        Filtering for depth of RNA-seq data is done by DESeq2. No need here.

        Parameters:
        - thrsh: int, threshold for pValue filter. DESeq2 FDR is very lenient, choose stringent cutoff 

        Returns:
            - pValFeatures: array of int, Indices of genes with pAdj < thrsh
        '''
        tmp=self.lrt.padj < thrsh
        self.pValFeatures = self.normCounts.reset_index().index.array[tmp]
        return self
    def corrNetwork(self,nPerm=1000,filterType=["pValue","depth"],corrThrsh="permutation",quantileCutoff=0.995):
        if not hasattr(self, "Scaled"):
            raise RuntimeError("Scaled not found. Did you forget to run .scale() ?")
        if filterType == "pValue":
            if not hasattr(self, "pValFeatures"):
                raise RuntimeError("pValFeatures not found. Did you forget to run .pValFilter() ?")
            features = self.pValFeatures
        elif filterType == "depth":
            if not hasattr(self, "DepthFeatures"):
                raise RuntimeError("DepthFeatures not found. Did you forget to run .depthFilter() ?")
            features = self.DepthFeatures
        if corrThrsh=="permutation":
            data=self.Scaled[features,:]
            CorrThrsh = np.zeros(nPerm)
            k=0
            while k < nPerm:
                print(k)
                tmp = data[np.arange(len(data))[:,None], np.random.randn(*data.shape).argsort(axis=1)]
                tmp = np.corrcoef(tmp)
                CorrThrsh[k] = np.quantile(tmp[np.triu_indices(tmp.shape[1],k=1)],quantileCutoff)
                k=k+1
            self.CorrThrsh = np.median(CorrThrsh)
        else: self.CorrThrsh=corrThrsh
        corr = np.corrcoef(self.Scaled[features,:])
        self.CorMatrix = corr
        G = corr-self.CorrThrsh
        G[G>0] = 1
        G[G<=0] = 0
        self.CorrNetwork = G
        self.NetworkFeatures = features
        return self
    def networkGeneration(self, remove = False):
        if not hasattr(self, "CorrNetwork"):
            raise RuntimeError("CorrNetwork not found. Did you forget to run .corrNetwork() ?")
        self.G = nx.from_numpy_array(self.CorrNetwork)
        if remove is True:
            delattr(self, "CorrNetwork")
    def louvain(self,minEdges=50,resolution=1):
        if not hasattr(self, "G"):
            raise RuntimeError("Graph not found. Did you forget to run .networkGeneration() ?")
        remove = [node for node,degree in dict(self.G.degree()).items() if degree < minEdges]
        self.G.remove_nodes_from(remove)
        clusters = nx_comm.louvain_communities(self.G, resolution=resolution)
        clusters.sort(key=len)
        tmp = []
        for index, cluster in enumerate(clusters): 
            print(index)
            g = self.G.subgraph(list(cluster))
            tmp.append(
                pd.DataFrame(data={"CommunityDegree":[val for (node, val) in g.degree()],
                                   "Cluster":index}, index = self.logCounts.index.array[self.NetworkFeatures[list(cluster)]])
                )
        self.Clusters = pd.concat(tmp)
        return self
    def umap(self,nNeighbors=15,minDist=0.01):
        reducer = UMAP(n_neighbors=nNeighbors,min_dist=minDist)
        self.Embedding = reducer.fit(self.Scaled[self.NetworkFeatures,:])
        return self
    def extractData(self,Network=False,UMAP=True,Scaled=True):
        self.Export = {}
        if Network:
            self.Export["Network"] = pd.DataFrame(
                data = self.CorrNetwork,
                index = self.logCounts.index.array[self.NetworkFeatures],
                columns = self.logCounts.index.array[self.NetworkFeatures]
                )
        if UMAP:
            self.Export["UMAP"] = pd.DataFrame(
                data = self.Embedding.embedding_,
                index = self.logCounts.index.array[self.NetworkFeatures],
                columns = ("UMAP" + pd.Series(np.arange(self.Embedding.n_components)+1).astype(str)))
        if Scaled:
            self.Export["Scaled"] = pd.DataFrame(
                data=self.Scaled,index=self.normCounts.index,columns=self.Samples)
