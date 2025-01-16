import networkx as nx
import numpy as np
import pandas as pd
import concurrent.futures
from tqdm.auto import tqdm
import sys

# Main entry point for the script
if __name__ == '__main__':
    # Print the experiment name passed as a command-line argument
    print(sys.argv[1]) # experiment

# Path to the network file
networkPath = "networks/stringFull850_v115.gexf"
# Load the network graph from a GEXF file
GR3D = nx.read_gexf(networkPath)

# Load the experimental data from a CSV file
combinedDataEnsF = pd.read_csv("combinedDataEnsF.csv", index_col=0, low_memory=False)

# Get the experiment name from the command-line argument
exp = sys.argv[1]

# Extract unique values from various columns in the data
cellTypes = list(combinedDataEnsF.cellType.unique())
dataTypes = list(combinedDataEnsF.dataType.unique())
designTypes = list(combinedDataEnsF.designType.unique())

# Remove empty strings from design types if they exist
if "" in designTypes:
    designTypes.remove("")

# Extract unique values for method types and remove empty strings
methodTypes = list(combinedDataEnsF.methodType.unique())
if "" in methodTypes:
    methodTypes.remove("")

# Extract unique values for folders and combinations
folds = list(combinedDataEnsF.folderName.unique())
combinations = list(combinedDataEnsF.combination.unique())

# Initialize a result dictionary to store the processed data
result = {}

# Loop over each cell type
for cell in tqdm(cellTypes):
    for dataT in dataTypes:
        # Check if the data type is RNA-Seq
        if dataT == "RNASeq":
            for design in designTypes:
                for fold in folds:
                    for comb in combinations:
                        # Filter the data based on the current parameters
                        tempDF = combinedDataEnsF[
                            (combinedDataEnsF.expType == exp) &
                            (combinedDataEnsF.cellType == cell) &
                            (combinedDataEnsF.designType == design) &
                            (combinedDataEnsF.dataType == "RNASeq") &
                            (combinedDataEnsF.folderName == fold) &
                            (combinedDataEnsF.combination == comb)]

                        # Proceed if the filtered data is not empty
                        if len(tempDF) > 0:
                            # Combine up-regulated and down-regulated gene IDs into a set
                            genes = set(tempDF.ensemblIDfinal_up).union(
                                set(tempDF.ensemblIDfinal_down))
                            
                            # Remove NaN values from the set
                            if np.nan in genes:
                                genes.remove(np.nan)

                            # Create a subgraph of the network with the selected genes
                            SP = GR3D.subgraph(list(genes))

                            # Calculate connectivity if there are at least 5 nodes in the subgraph
                            if SP.number_of_nodes() >= 5:
                                connectivity = 1 - (
                                    nx.number_of_isolates(SP) /
                                    SP.number_of_nodes())
                            else:
                                connectivity = np.nan

                            # Add the computed metrics to the result dictionary
                            result[len(result)] = {
                                "network": "String850",
                                "experiment": exp,
                                "cellType": cell,
                                "data": "RNASeq",
                                "design": design,
                                "method": "",
                                "designANDmethod": design,
                                "expDesign": fold,
                                "combination": comb,
                                "numberOfDiffExp": len(list(genes)),
                                "numberOfNodes": SP.number_of_nodes(),
                                "connectivity": connectivity
                            }

        # Check if the data type is Methylation
        elif dataT == "Methylation":
            for method in methodTypes:
                for fold in folds:
                    for comb in combinations:
                        # Filter the data based on the current parameters
                        tempDF = combinedDataEnsF[
                            (combinedDataEnsF.expType == exp) &
                            (combinedDataEnsF.cellType == cell) &
                            (combinedDataEnsF.methodType == method) &
                            (combinedDataEnsF.dataType == "Methylation") &
                            (combinedDataEnsF.folderName == fold) &
                            (combinedDataEnsF.combination == comb)]

                        # Proceed if the filtered data is not empty
                        if len(tempDF) > 0:
                            # Combine up-regulated and down-regulated gene IDs into a set
                            genes = set(tempDF.ensemblIDfinal_up).union(
                                set(tempDF.ensemblIDfinal_down))
                            
                            # Remove NaN values from the set
                            if np.nan in genes:
                                genes.remove(np.nan)

                            # Create a subgraph of the network with the selected genes
                            SP = GR3D.subgraph(list(genes))

                            # Calculate connectivity if there are at least 5 nodes in the subgraph
                            if SP.number_of_nodes() >= 5:
                                connectivity = 1 - (
                                    nx.number_of_isolates(SP) /
                                    SP.number_of_nodes())
                            else:
                                connectivity = np.nan

                            # Add the computed metrics to the result dictionary
                            result[len(result)] = {
                                "network": "String850",
                                "experiment": exp,
                                "cellType": cell,
                                "data": "Methylation",
                                "design": "",
                                "method": method,
                                "designANDmethod": method,
                                "expDesign": fold,
                                "combination": comb,
                                "numberOfDiffExp": len(list(genes)),
                                "numberOfNodes": SP.number_of_nodes(),
                                "connectivity": connectivity
                            }

# Convert the result dictionary to a DataFrame
resultDF = pd.DataFrame.from_dict(result, orient="index")

# Save the DataFrame to a CSV file
resultDF.to_csv(exp+"_resultDF_string850.csv")
