# Implementing automatic state decomposition
# (from Automatic discovery of metastable states for the construction of Markov
# models of macromolecular conformational dynamics
# Chodera, Singhal, Pande, Dill, Swope, 2007, JCP)
# in Python

import numpy as np
import numpy.random as npr
import matplotlib.pyplot as plt
from optimizers import MCSA

# define some objective functions suggested in section III.C.2.

def metastability_objective(T,**kwargs):
   '''
   metastability of each state, summed over all states

   \sum_i T_{i,i}

   '''

   return np.trace(T)/len(T)

def weighted_metastability_objective(T,w,**kwargs):
   '''
   metastability of each state * weight of each state, summed over all states

   \sum_i w_i T_{i,i}

   where w is an arbitrary weight vector

   '''

   return np.dot(np.diag(T),w)

def pi_weighted_metastability_objective(T,**kwargs):
   '''
   metastability of each state * its stationary probability, summed over all states

   \sum_i \pi_i T_{i,i}

   where \pi is the stationary distribution implied by T '''

   evals,evecs = np.linalg.eigh(M)
   pi = evecs[:,np.argmax(evals)]
   return weighted_metastability_objective(T,pi)

def state_lifetimes(T,**kwargs):
   ''' \tau_i = (1-T_{i,i})^{-1} '''

   T_ii = np.diag(T)
   return (1-T_ii)**-1

def state_lifetimes_objective(T,**kwargs):
   ''' \sum_i \tau_i '''

   return np.sum(state_lifetimes(T))


def count_weighted_metastability_objective(T,counts,cg_map,**kwargs):
   '''
   weight each macrostate's metastability by the number of observed simulation
   frames corresponding to that macrostate

   '''
   macrostate_counts = np.zeros(len(T))
   for i in range(len(cg_map)):
      macrostate_counts[cg_map[i]] += counts[i]
   return weighted_metastability_objective(T,macrostate_counts)

objectives = {'metastability':metastability_objective,
              'pi_weighted_metastability':pi_weighted_metastability_objective,
              'state_lifetimes':state_lifetimes_objective,
              'count_weighted_metastability':count_weighted_metastability_objective}

class CoarseGrain():
   def __init__(self,
                n_macrostates=10,
                objective='metastability',
                max_iter=10000):

      self.n_macrostates=n_macrostates
      self.objective_name = objective
      self.objective_function = objectives[objective]
      self.max_iter=max_iter

   def proposal(self,cg_map_old):
      '''
      Creates a new array, with 2 elements of cg_map_old swapped at random.

      Parameters
      ----------
      cg_map_old : array-like
         current coarse-graining map

      Returns
      -------
      cg_map : numpy.ndarray
         proposed coarse-graining map

      '''

      cg_map = np.array(cg_map_old)
      ind1,ind2 = npr.randint(0,len(cg_map),2)
      tmp = cg_map[ind1]
      cg_map[ind1] = cg_map[ind2]
      cg_map[ind2] = tmp
      return cg_map

   def cg_T(self,microstate_T,microstate_pi, cg_map):
      ''' Coarse-grain a microstate transition matrix by applying cg_map
      Parameters
      ----------
      microstate_T : (N,N), array-like, square
         microstate transition matrix
      microstate_pi : (N,), array-like
         microstate stationary distribution
      cg_map : (N,), array-like
         assigns each microstate i to a macrostate cg_map[i]

      Returns
      -------
      T : numpy.ndarray, square
         macrostate transition matrix
      '''

      n_macrostates = np.max(cg_map)+1
      n_microstates = len(microstate_T)

      # compute macrostate stationary distribution
      macrostate_pi = np.zeros(n_macrostates)
      for i in range(n_microstates):
         macrostate_pi[cg_map[i]] += microstate_pi[i]
      macrostate_pi /= np.sum(macrostate_pi)

      # accumulate macrostate transition matrix
      T = np.zeros((n_macrostates,n_macrostates))
      for i in range(n_microstates):
         for j in range(n_microstates):
            T[cg_map[i],cg_map[j]] += microstate_pi[i] * microstate_T[i,j]

      # normalize
      for a in range(n_macrostates):
         T[a] /= macrostate_pi[a]

      return T

   def fit(self,microstate_T,microstate_pi,
           init_cg_map=None,
           microstate_counts=None):
      # if numba is installed, use JIT compilation for ~100x speed-ups
      # of cg_T (which contains a nested for loop)
      try:
         from numba import jit
         self.cg_T = jit(self.cg_T)
         print('Successfully JIT-compiled inner loop! :)')
      except:
         print('Did not JIT-compile inner loop :(')

      def objective(cg_map):
         return self.objective_function(self.cg_T(microstate_T,microstate_pi,cg_map),
                    counts=microstate_counts,cg_map=cg_map)

      mcsa = MCSA(proposal_function=self.proposal,
                  objective_function=objective,
                  annealing_schedule=np.logspace(0,4,self.max_iter))
      if init_cg_map==None:
         init_cg_map = npr.randint(0,self.n_macrostates,len(microstate_T))
      solns = mcsa.maximize(init_cg_map)
      self.solns = solns
      print('Optimization complete: coarse-grained {0} = {1:.3f}'.format(
                self.objective_name,solns[-1][1]
                ))
      self.cg_map = solns[-1][0]

      self.optimization_trace = np.array([s[1] for s in solns])
