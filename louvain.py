import numpy as np
import pandas as pd

def make_edgelist():
   nn = np.loadtxt('nn.csv', delimiter=',')

   edgelist = list()
   for v in range(nn.shape[0]):
      edgelist.append( list() )
      for i in range(nn.shape[1]):
         if i not in edgelist[v]:
            edgelist[v].append( int(nn[v,i]) )
   return edgelist 

def actual_fraction_of_internal_edges( edgelist, group ):
   n_total = 0
   n_internal = 0
   for v in range(len(edgelist)):
      neighbors = edgelist[v]
      for w in neighbors:
         n_total += 1
         if group[v] == group[w]:
            n_internal += 1
   return n_internal / n_total

def expected_fraction_of_internal_edges( edgelist, group ):
   n_groups = group.max() + 1
   edges_per_group = np.zeros( n_groups )
   for v in range(len(edgelist)):
      g = group[v]
      edges_per_group[g] += len(edgelist[v])
   return ( ( edges_per_group / edges_per_group.sum() )**2 ).sum()

def modularity_score( edgelist, group ):
   return actual_fraction_of_internal_edges( edgelist, group ) - \
      expected_fraction_of_internal_edges( edgelist, group ) 

edgelist = make_edgelist()
#print( "neighbours of vertex 17:", edgelist[17] )

tbl = pd.read_csv( "ump_cl.csv" )
group = tbl['group'].values - 1
print( modularity_score( edgelist, group ) )

group = np.arange( len(edgelist) )

n_in_group = np.ones( len(edgelist), dtype=np.integer )

for v in range(len(edgelist)):
   print( "optimizing vertex", v )
   best_score = -1
   best_group = -1
   for g in range(len(n_in_group)):
      group[v] = g
      s = modularity_score( edgelist, group )
      if s > best_score:
         best_score = s
         best_group = g
   group[v] = best_group
   print( best_score )	


print( modularity_score( edgelist, group ) )
