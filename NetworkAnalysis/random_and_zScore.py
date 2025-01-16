
import networkx as nx
import numpy as np
import pandas as pd
import concurrent.futures
from tqdm.auto import tqdm
import sys
import configparser
import random

numRandoms=5000
currencyMetaPer=2.0


if __name__ == '__main__':
    config = configparser.ConfigParser()
    print(sys.argv[1])#network

network = str(sys.argv[1])


mainPath = "../../networks/"
if network == "string850":
  networkPath = mainPath+"stringFull850_v115.gexf"
elif network == "string400":
  networkPath = mainPath+"stringFull400_v115.gexf"
elif network == "biogrid":
  networkPath = mainPath+"biogrid_human_44218.gexf"
elif network == "recon22":
  networkPath = mainPath+"networks_1Feb2023/recon2_v04"+"_CurrMet"+str(currencyMetaPer).replace(".", "")+".gexf"
elif network == "recon3D2":
  networkPath = mainPath+"networks_1Feb2023/recon3D_v301"+"_CurrMet"+str(currencyMetaPer).replace(".", "")+".gexf"


result1 = pd.read_csv("result/Experiment1_resultDF_"+network+"_08082023.csv", index_col=0)
result2 = pd.read_csv("result/Experiment2_resultDF_"+network+"_08082023.csv", index_col=0)
result3 = pd.read_csv("result/Experiment3_resultDF_"+network+"_08082023.csv", index_col=0)

numberNodesList = set(result1.numberOfNodes).union(set(result2.numberOfNodes).union(set(result3.numberOfNodes)))
numberNodesList = [i for i in numberNodesList if i >=5]


G = nx.read_gexf(networkPath)

randomsDFAll = pd.DataFrame(columns={"numberSubgraphNodes", "graphNo", "connectivity"})

for nn in numberNodesList:
    
    connectivity = []
    
    for r in range(numRandoms):
        
        randomGenes = random.sample(list(G.nodes()), nn)
        
        SR = G.subgraph(randomGenes)
        conn = 1-(nx.number_of_isolates(SR) / SR.number_of_nodes())

        connectivity.append(conn)

    randomsDF = pd.DataFrame({"numberSubgraphNodes": nn,
                                    "graphNo": list(range(1,numRandoms+1)),
                                    "connectivity": connectivity})
    randomsDFAll = randomsDFAll.append(randomsDF)
    
randomsDFAll.to_csv("result/randomsAll_"+network+".csv")



randoms_summary = randomsDFAll.groupby("numberSubgraphNodes").agg(conn_mean=('connectivity', 'mean'),
                                                                           conn_std=('connectivity', 'std'))
randoms_summary = randoms_summary.reset_index()

reconResult1 = pd.merge(result1, randoms_summary, left_on="numberOfNodes", right_on="numberSubgraphNodes", how="left")
reconResult2 = pd.merge(result2, randoms_summary, left_on="numberOfNodes", right_on="numberSubgraphNodes", how="left")
reconResult3 = pd.merge(result3, randoms_summary, left_on="numberOfNodes", right_on="numberSubgraphNodes", how="left")

reconResult = pd.concat([reconResult1, reconResult2, reconResult3], axis=0).reset_index(drop=True)

reconResult.conn_std = reconResult.conn_std.replace({0:1})
reconResult["zScore"] = (reconResult.connectivity-reconResult.conn_mean)/reconResult.conn_std

reconResult.to_csv("result/randomsSummary_"+network+".csv")









