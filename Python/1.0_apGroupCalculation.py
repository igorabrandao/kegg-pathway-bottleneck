# Data manipulation import
import pandas as pd
import glob, os

# Class import
from graph import Graph

# Select the .csv files
folder = r'../output/allGraphs'
all_files = glob.glob(os.path.join(folder, "*.csv"))

# Define the graph array
graphArray = []

# Read all csv files
for filename in all_files:
    df = pd.read_csv(filename, index_col=None, header=0)
    graphArray.append(df)

pathway = 160

# Create a graph object from the pathway list
g1 = Graph(len(graphArray[pathway]))

# Add the nodes to the graph
for index in range(len(graphArray[pathway])):
    g1.addEdge(graphArray[pathway]['node1'][index], graphArray[pathway]['node2'][index])

# Create a list of unique nodes from the current graph
nodeList = graphArray[pathway][['node1', 'node2']].stack().reset_index(level=[0,1], drop=True)
nodeList = nodeList.unique()

# Define the result array
result = []

# Remove the unique nodes from the graph and test the AP groups
for i in range(0, len(nodeList)):
    print('Remove root ' + nodeList[i])
    worker(g1, nodeList[i], result, nodeList[i])
    
print(result)