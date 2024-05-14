import sys
import numpy as np
import pandas as pd

# We represent our undirected graph as a list of lists.
# The outer list has one element for each vertex; this
# element is a list of that vertex's neighbors.
# We call this a "list of neighbour list" (LNL)

def make_lnl():
   """This function loads the nearest-neighbor CSV
   file that we made in R and constructs the LNL.
   """

   nn = np.loadtxt('nn.csv', delimiter=',')

   # initialize with empty list
   nln = list()
   for v in range(nn.shape[0]):
      nln.append( list() )
   # Go through cells
   for v1 in range(nn.shape[0]):
      # Go through v1's nearest neighbors
      for i in range(nn.shape[1]):
         v2 = int(nn[v1,i])
         # Add edge to both neighborhoods (if necessary)
         if v2 not in nln[v1]:
            nln[v1].append( v2 )
         if v1 not in nln[v2]:
            nln[v2].append( v1 )
   return nln 

def actual_fraction_of_internal_edges( lnl, membv ):
   """This function goes through all edges and
   counts how many of them are internal, i.e.,
   connecting vertices in the same group."""
   n_total = 0     # this counts stubs (i.e., half-edges)
   n_internal = 0  # this, too
   for v in range(len(lnl)):
      neighbors = lnl[v]
      for w in neighbors:
         n_total += 1
         if membv[v] == membv[w]:
            n_internal += 1
   return n_internal / n_total

def edge_stubs_per_group( lnl, membv ):
   n_groups = membv.max() + 1
   answer = np.zeros( n_groups )
   for v in range(len(lnl)):
      g = membv[v]
      answer[g] += len(lnl[v])
   return answer

def expected_fraction_of_internal_edges( lnl, membv ):
   espg = edge_stubs_per_group( lnl, membv )
   return ( ( espg / espg.sum() )**2 ).sum()

def modularity_score( edgelist, group ):
   return actual_fraction_of_internal_edges( lnl, membv ) - \
      expected_fraction_of_internal_edges( lnl, membv ) 

lnl = make_lnl()

tbl = pd.read_csv( "ump_cl.csv" )
membv = tbl['group'].values - 1

print( modularity_score(lnl, membv))

# How many edges are there?
n_edges = 0
for nl in lnl:
   n_edges += len(nl)
n_edges /= 2
print( n_edges )   # 21176

### Try: Move a vertex from one group to another
v = 17
g0 = membv[v]
g1 = 7

# How many actual edges will we gain or lose?
delta_n_internal_edges = 0
for nv in lnl[v]:
   if membv[nv] == g0:
      delta_n_internal_edges -= 1
   elif membv[nv] == g1:
      delta_n_internal_edges += 1
print( delta_n_internal_edges / n_edges )

af0 = actual_fraction_of_internal_edges(lnl, membv)
membv[v] = g1
print( actual_fraction_of_internal_edges(lnl, membv) - af0 )

# How does the expected fraction change?

# First count the edges per group:
n_groups = membv.max() + 1
edges_per_group = np.zeros( n_groups )
for v in range(len(lnl)):
   g = membv[v]
   edges_per_group[g] += len(lnl[v])
