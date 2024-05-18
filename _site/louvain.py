import sys
import numpy as np
import pandas as pd
import numba

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
   lnl = list()
   for v in range(nn.shape[0]):
      lnl.append( list() )
   # Go through cells
   for v1 in range(nn.shape[0]):
      # Go through v1's nearest neighbors
      for i in range(nn.shape[1]):
         v2 = int(nn[v1,i])
         # Add edge to both neighborhoods (if necessary)
         if v2 not in lnl[v1]:
            lnl[v1].append( v2 )
         if v1 not in lnl[v2]:
            lnl[v2].append( v1 )
   lnl2 = numba.typed.List()
   for l in lnl:
      lnl2.append( np.array(l) )
   return lnl2

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

@numba.njit
def edge_stubs_per_group( lnl, membv ):
   n_groups = len(lnl) #membv.max() + 1
   answer = np.zeros( n_groups )
   for v in range(len(lnl)):
      g = membv[v]
      answer[g] += len(lnl[v])
   return answer

def expected_fraction_of_internal_edges( lnl, membv ):
   espg = edge_stubs_per_group( lnl, membv )
   return ( espg**2 ).sum() / espg.sum()**2

def modularity_score( lnl, membv ):
   return actual_fraction_of_internal_edges( lnl, membv ) - \
      expected_fraction_of_internal_edges( lnl, membv ) 

@numba.njit
def delta_modularity_score( lnl, membv, vertex, new_group, espg ):
   old_group = membv[vertex]
   n_edges = espg.sum()/2
   # change in actual egde fraction
   delta_n_internal_edges = 0
   for nv in lnl[vertex]:
      if membv[nv] == old_group:
         delta_n_internal_edges -= 1
      elif membv[nv] == new_group:
         delta_n_internal_edges += 1
   delta_actual_fraction = delta_n_internal_edges / n_edges
   # change expected edge fraction
   numerator_change = ( espg[old_group] - len(lnl[vertex]) )**2 - espg[old_group]**2 + \
      ( espg[new_group] + len(lnl[vertex]) )**2 - espg[new_group]**2
   delta_expected_fraction = numerator_change / (2*n_edges)**2
   return delta_actual_fraction - delta_expected_fraction

@numba.njit
def optimize( lnl, membv ):
   
   for i in range(300):      
      v = np.random.randint(0, len(lnl))
      espg = edge_stubs_per_group( lnl, membv )
      best_delta = -1
      best_group = -1
      for g in range(len(lnl)):
         if espg[g] == 0:
            continue
         delta = delta_modularity_score( lnl, membv, v, g, espg )
         if delta > best_delta:
            best_delta = delta
            best_group = g
      if best_delta > 0:
         espg[membv[v]] -= len(lnl[v])
         membv[v] = best_group
         espg[best_group] += len(lnl[v])
   return membv


lnl = make_lnl()
membv = np.arange( len(lnl) )

tbl = pd.read_csv( "ump_cl.csv" )
#membv = tbl['group'].values - 1

ms0 = modularity_score(lnl, membv)
print( ms0 )

espg = edge_stubs_per_group( lnl, membv )
print( delta_modularity_score( lnl, membv, 17, 7, espg ) )
membv[17] = 7
print( modularity_score(lnl, membv) - ms0 )


espg = edge_stubs_per_group( lnl, membv )
n_groups = len(espg)

for i in range(2000):
   membv = optimize( lnl, membv )
   print( modularity_score(lnl, membv), np.unique(membv).size )
