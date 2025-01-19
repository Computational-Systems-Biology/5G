# Import necessary libraries
import networkx as nx  # For graph manipulation and analysis
import numpy as np  # For numerical operations
import pandas as pd  # For data manipulation and analysis
import concurrent.futures  # For parallel processing
from tqdm.auto import tqdm  # For progress bars
import sys  # For handling command-line arguments
import configparser  # For reading configuration files
import random  # For generating random samples

# Parameters for random sampling and file naming
numRandoms = 5000  # Number of random subgraphs to generate for each node count
currencyMetaPer = 2.0  # Parameter for network-specific file paths

if __name__ == '__main__':
    # Read the network name from the first command-line argument
    config = configparser.ConfigParser()
    print(sys.argv[1])  # Debugging print for the network name

# Store the network name for use throughout the script
network = str(sys.argv[1])

# Define the base path to network files
mainPath = "../../networks/"

# Map the network name to its corresponding file path
if network == "string850":
    networkPath = mainPath + "stringFull850_v115.gexf"
elif network == "string400":
    networkPath = mainPath + "stringFull400_v115.gexf"
elif network == "biogrid":
    networkPath = mainPath + "biogrid_human_44218.gexf"
elif network == "recon22":
    networkPath = mainPath + "networks_1Feb2023/recon2_v04" + "_CurrMet" + str(currencyMetaPer).replace(".", "") + ".gexf"
elif network == "recon3D2":
    networkPath = mainPath + "networks_1Feb2023/recon3D_v301" + "_CurrMet" + str(currencyMetaPer).replace(".", "") + ".gexf"

# Load experimental results from CSV files
result1 = pd.read_csv("result/Experiment1_resultDF_" + network + "_08082023.csv", index_col=0)
result2 = pd.read_csv("result/Experiment2_resultDF_" + network + "_08082023.csv", index_col=0)
result3 = pd.read_csv("result/Experiment3_resultDF_" + network + "_08082023.csv", index_col=0)

# Extract the unique node counts from all experiments, filtering out small values (< 5)
numberNodesList = set(result1.numberOfNodes).union(set(result2.numberOfNodes).union(set(result3.numberOfNodes)))
numberNodesList = [i for i in numberNodesList if i >= 5]

# Load the specified network graph
G = nx.read_gexf(networkPath)

# DataFrame to store results of random subgraph analyses
randomsDFAll = pd.DataFrame(columns={"numberSubgraphNodes", "graphNo", "connectivity"})

# Loop through each unique node count
for nn in numberNodesList:
    connectivity = []  # Store connectivity values for current node count

    # Generate random subgraphs and calculate connectivity
    for r in range(numRandoms):
        randomGenes = random.sample(list(G.nodes()), nn)  # Randomly select nodes
        SR = G.subgraph(randomGenes)  # Create subgraph
        conn = 1 - (nx.number_of_isolates(SR) / SR.number_of_nodes())  # Calculate connectivity
        connectivity.append(conn)

    # Create a DataFrame for the current node count and append it to the main DataFrame
    randomsDF = pd.DataFrame({
        "numberSubgraphNodes": nn,
        "graphNo": list(range(1, numRandoms + 1)),
        "connectivity": connectivity
    })
    randomsDFAll = randomsDFAll.append(randomsDF)

# Save the complete random subgraph results to a CSV file
randomsDFAll.to_csv("result/randomsAll_" + network + ".csv")

# Compute mean and standard deviation of connectivity for each node count
randoms_summary = randomsDFAll.groupby("numberSubgraphNodes").agg(
    conn_mean=('connectivity', 'mean'),
    conn_std=('connectivity', 'std')
).reset_index()

# Merge the summary statistics with experimental results
combResult1 = pd.merge(result1, randoms_summary, left_on="numberOfNodes", right_on="numberSubgraphNodes", how="left")
combResult2 = pd.merge(result2, randoms_summary, left_on="numberOfNodes", right_on="numberSubgraphNodes", how="left")
combResult3 = pd.merge(result3, randoms_summary, left_on="numberOfNodes", right_on="numberSubgraphNodes", how="left")

# Combine results from all experiments into a single DataFrame
combResult = pd.concat([combResult1, combResult2, combResult3], axis=0).reset_index(drop=True)

# Replace standard deviation of zero with one to avoid division by zero in z-score calculation
combResult.conn_std = combResult.conn_std.replace({0: 1})

# Calculate z-scores for connectivity
combResult["zScore"] = (combResult.connectivity - combResult.conn_mean) / combResult.conn_std

# Save the final results with z-scores to a CSV file
combResult.to_csv("result/randomsSummary_" + network + ".csv")
