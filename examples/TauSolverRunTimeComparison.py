
# coding: utf-8

# In[1]:


# %matplotlib
# %matplotlib inline
import numpy
import matplotlib.pyplot as plt
import time


# In[2]:


import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.getcwd(), os.pardir)))
print(sys.path)


# In[3]:


import math
import gillespy2
from gillespy2.basic_tau_hybrid_solver import BasicTauHybridSolver
from gillespy2.basic_ssa_solver import BasicSSASolver
from gillespy2.basic_tau_leaping_solver import BasicTauLeapingSolver


# In[4]:


class SimpleHybridModel(gillespy2.Model):
     def __init__(self, parameter_values=None, init_v=1):
            #initialize Model
            gillespy2.Model.__init__(self, name="Simple_Hybrid_Model")

            
            #Species
            A = gillespy2.Species(name='A', initial_value=0)
            V = gillespy2.Species(name='V', initial_value=init_v)

            self.add_species([A, V])
            
            #parameters
            rate1 = gillespy2.Parameter(name='rate1', expression= 20.0)
            rate2 = gillespy2.Parameter(name='rate2', expression= 10.0)
#             rate_rule1 = gillespy2.RateRule(V, "cos(t)")
            self.add_parameter([rate1, rate2])
#             self.add_rate_rule(rate_rule1)
            
            #reactions
            r1 = gillespy2.Reaction(name="r1",reactants={}, products={A:1},
                   propensity_function="rate1 * V")
            
            r2 = gillespy2.Reaction(name="r2",reactants={A:1}, products={},
                    rate=rate2)
            
            self.add_reaction([r1,r2])
            self.timespan(numpy.linspace(0,100, 101))


# In[5]:


model = SimpleHybridModel()


# In[6]:


get_ipython().run_line_magic('time', 'results1 = model.run(solver=BasicSSASolver(), show_labels=True)')
get_ipython().run_line_magic('time', 'results2 = model.run(solver=BasicTauLeapingSolver(), show_labels=True)')
get_ipython().run_line_magic('time', 'results3 = model.run(solver=BasicTauHybridSolver(), show_labels=True)')


# In[7]:


v_range = range(1, 15)
def run_test(solver, v_range):
    run_data = []
    for n in v_range:
        model = SimpleHybridModel(init_v=n)
        time_start = time.perf_counter()
        model.run(solver=solver, show_labels=True)
        time_end = time.perf_counter()
        run_data.append(time_end-time_start)
    return run_data


# In[8]:


timing_data = {'basic':[], 'tau':[], 'hybrid_tau':[]}
get_ipython().run_line_magic('time', "timing_data['basic'] = run_test(BasicSSASolver(), v_range)")
get_ipython().run_line_magic('time', "timing_data['tau'] = run_test(BasicTauLeapingSolver(), v_range)")
get_ipython().run_line_magic('time', "timing_data['hybrid_tau'] = run_test(BasicTauHybridSolver(), v_range)")
print(timing_data)


# In[9]:


plt.figure(figsize=(20,10))
plt.title("Time Comparison of Solvers")
plt.xlabel("initial Value of V")
plt.ylabel("Simulation Run Time")
plt.plot(v_range, timing_data['basic'], label='basic')
plt.plot(v_range, timing_data['tau'], label='tau')
plt.plot(v_range, timing_data['hybrid_tau'], label='hybrid_tau')
plt.legend(loc='best')
plt.savefig("TimeComparisonOfSolvers.pdf")

