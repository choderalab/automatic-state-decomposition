import numpy as np
import numpy.random as npr

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
      Maximize the objective function, given an initial solution

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
         if npr.rand() < np.exp(delta_Q * beta):
            solutions.append((proposal,f_proposal))
         else:
            solutions.append((old_soln,old_f))

      return solutions
