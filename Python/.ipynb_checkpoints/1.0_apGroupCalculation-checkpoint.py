# Python program to find articulation points in an undirected graph
import copy
import itertools
from collections import defaultdict
import sys
sys.setrecursionlimit(10**6)

#This code is contributed by Neelam Yadav
#https://www.geeksforgeeks.org/articulation-points-or-cut-vertices-in-a-graph/
#This class represents an undirected graph 
#using adjacency list representation
class Graph:
  
    def __init__(self,vertices):
        self.V= vertices #No. of vertices
        self.graph = defaultdict(list) # default dictionary to store graph
        self.Time = 0
  
    # function to add an edge to graph
    def addEdge(self,u,v):
        self.graph[u].append(v)
        self.graph[v].append(u)
    
    def removeEdge(self, u, v):
        self.graph[u].remove(v)
        self.graph[v].remove(u)
    
    def removeNode(self, u):
        copy_u = copy.copy(self.graph[u])
        for i in copy_u:
            self.removeEdge(u, i)
        del self.graph[u]

    def combineNode(self, u, v):
        for i in self.graph:
            if u in self.graph[i]:
                self.graph[i].remove(u)
                if v not in self.graph[i]:
                    self.graph[i].append(v)
            
  
    '''A recursive function that find articulation points 
    using DFS traversal
    u --> The vertex to be visited next
    visited[] --> keeps tract of visited vertices
    disc[] --> Stores discovery times of visited vertices
    parent[] --> Stores parent vertices in DFS tree
    ap[] --> Store articulation points'''
    def APUtil(self,u, visited, ap, parent, low, disc):
 
        #Count of children in current node 
        children =0
 
        # Mark the current node as visited and print it
        visited[u]= True
 
        # Initialize discovery time and low value
        disc[u] = self.Time
        low[u] = self.Time
        self.Time += 1
 
        #Recur for all the vertices adjacent to this vertex
        for v in self.graph[u]:
            # If v is not visited yet, then make it a child of u
            # in DFS tree and recur for it
            if visited[v] == False :
                parent[v] = u
                children += 1
                self.APUtil(v, visited, ap, parent, low, disc)
 
                # Check if the subtree rooted with v has a connection to
                # one of the ancestors of u
                low[u] = min(low[u], low[v])
 
                # u is an articulation point in following cases
                # (1) u is root of DFS tree and has two or more chilren.
                if parent[u] == -1 and children > 1:
                    ap[u] = True
 
                #(2) If u is not root and low value of one of its child is more
                # than discovery value of u.
                if parent[u] != -1 and low[v] >= disc[u]:
                    ap[u] = True   
                     
                # Update low value of u for parent function calls    
            elif v != parent[u]: 
                low[u] = min(low[u], disc[v])
 
 
    #The function to do DFS traversal. It uses recursive APUtil()
    def AP(self):
  
        # Mark all the vertices as not visited 
        # and Initialize parent and visited, 
        # and ap(articulation point) arrays
        visited = [False] * (self.V)
        disc = [float("Inf")] * (self.V)
        low = [float("Inf")] * (self.V)
        parent = [-1] * (self.V)
        ap = [False] * (self.V) #To store articulation points
 
        # Call the recursive helper function
        # to find articulation points
        # in DFS tree rooted with vertex 'i'
        for i in range(self.V):
            if visited[i] == False:
                self.APUtil(i, visited, ap, parent, low, disc)
 
        result = []
        for index, value in enumerate (ap):
            if value == True: 
                result.append(index)

        return result

      
#This is by Michael KH Tai, 31/07/2018
def worker(g, prefix, arr, u):
    
    container = g.graph[u]
    if u in g.AP():
        print prefix
        arr.append([prefix])
        del g
    else:
        for i in container:
            if str(i) not in prefix.split(' '):
                new_prefix = prefix + ' ' + str(i)
                new_g = copy.deepcopy(g)
                new_g.combineNode(u, i)
                if len(new_g.graph) > 1:
                    worker(new_g, new_prefix, arr, i)
    

 # Create a graph given in the above diagram
total = 16
'''
struct = {
    1:[2,3,4],
    2:[4,5],
    3:[5,6],
    4:[6],
    5:[6]
}
'''
struct = {
    0:[1,12,13,14,15],
    1:[2,12,14,15],
    2:[3,4,5,14],
    3:[4],
    4:[5,6],
    5:[6,9],
    6:[7,8,9],
    7:[8,10,11],
    8:[9,10,11],
    9:[10],
    10:[11],
    12:[15],
    13:[14,15]
} 
g1 = Graph (total)

for key,value in struct.iteritems():
    for child in value:
        g1.addEdge(key, child)

result = []
for i in range(0, total):
    print 'Remove root ' + str(i)
    worker(g1, str(i), result, i)
print result

'''
g1 = Graph (16)
g1.addEdge(0, 1)
g1.addEdge(0, 12)
g1.addEdge(0, 13)
g1.addEdge(0, 14)
g1.addEdge(0, 15)
g1.addEdge(1, 2)
g1.addEdge(1, 12)
g1.addEdge(1, 13)
g1.addEdge(1, 14)
g1.addEdge(1, 15)
g1.addEdge(2, 3)
g1.addEdge(2, 4)
g1.addEdge(2, 5)
g1.addEdge(2, 14)
g1.addEdge(3, 4)
g1.addEdge(4, 5)
g1.addEdge(4, 6)
g1.addEdge(5, 6)
g1.addEdge(5, 9)
g1.addEdge(6, 7)
g1.addEdge(6, 8)
g1.addEdge(6, 9)
g1.addEdge(7, 8)
g1.addEdge(7, 10)
g1.addEdge(7, 11)
g1.addEdge(8, 9)
g1.addEdge(8, 10)
g1.addEdge(8, 11)
g1.addEdge(9, 10)
g1.addEdge(10, 11)
g1.addEdge(12, 15)
g1.addEdge(13, 14)
g1.addEdge(13, 15)
print "\nArticulation points in graph "
g1.AP()
g2 = copy.deepcopy(g1)
print "\nArticulation points in 2 graph "
g2.AP()