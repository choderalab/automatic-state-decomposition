# Implementing automatic state decomposition
# (from Automatic discovery of metastable states for the construction of Markov
# models of macromolecular conformational dynamics
# Chodera, Singhal, Pande, Dill, Swope, 2007, JCP)
# in Python

import numpy as np
import numpy.random as npr
import matplotlib.pyplot as plt

# define some objective functions suggested in section III.C.2.

def metastability_objective(T):
   '''
   metastability of each state, summed over all states

   \sum_i T_{i,i}

   '''

   return np.trace(T)

def weighted_metastability_objective(T,w):
   '''
   metastability of each state * weight of each state, summed over all states

   \sum_i w_i T_{i,i}

   where \pi is the stationary distribution implied by T '''

   return np.dot(np.diag(T),w)

def pi_weighted_metastability_objective(T):
   '''
   metastability of each state * its stationary probability, summed over all states

   \sum_i \pi_i T_{i,i}

   where \pi is the stationary distribution implied by T '''

   evals,evecs = np.linalg.eigh(M)
   pi = evecs[:,np.argmax(evals)]
   return weighted_metastability_objective(T,pi)

def state_lifetimes(T):
   ''' \tau_i = (1-T_{i,i})^{-1} '''

   T_ii = np.diag(T)
   return (1-T_ii)**-1

def state_lifetimes_objective(T):
   ''' \sum_i \tau_i '''

   return np.sum(state_lifetimes(T))


def count_weighted_metastability_objective(T,counts):
   '''
   weight each macrostate's metastability by the number of observed simulation frames corresponding to that macrostate

   '''
   raise NotImplementedError
   #return weighted_metastability_objective(T,microstate_counts)

class MCSA():
   def __init__(self,
             proposal_function,
             objective_function,
             annealing_schedule):
      '''

      Monte Carlo Simulated Annealing

      Parameters
      ----------

      proposal_function : function
         accepts a solution object and returns a different solution object

      objective_function : function
         accepts a solution object and returns a real number

      annealing_schedule : array-like
         inverse temperature at each step of optimization

      '''
      self.proposal_function = proposal_function
      self.objective_function = objective_function
      self.annealing_schedule = annealing_schedule

   def maximize(self,init_solution):
      '''

      Parameters
      ----------

      init_solution : object
         object that can be passed to proposal_function or objective_function

      Returns
      -------
      solutions : list of tuples
         each tuple contains:
            (solution object, objective function value)

      '''

      solutions = [(init_solution,self.objective_function(init_solution))]


      for beta in self.annealing_schedule:
         old_soln,old_f=solutions[-1]
         proposal = self.proposal_function(old_soln)
         f_proposal = self.objective_function(proposal)
         delta_Q = f_proposal - old_f
         if npr.rand() < np.exp(delta_G * beta):
            solutions.append((proposal,f_proposal))
         else:
            solutions.append((old_soln,old_f))

      return solutions

