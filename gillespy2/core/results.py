import matplotlib.pyplot as plt

import warnings

from collections import UserDict

class results(UserDict):
    """ Result object for a gillespy2 solver, extends the dict object. """

    def __init__(self,data,model_name = "Undefined model name",solver_name = "Undefined solver name"):
        self.data = data
        self.model_name = model_name
        self.solver_name = solver_name

    def plot(self,xlabel = "Time (s)",ylabel = "Species Population",style="default"):

        plt.figure(figsize=(18, 10))
        plt.suptitle(self.model_name,fontsize=18)
        plt.title(self.solver_name,fontsize=12,y=1.05)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        try:
            plt.style.use(style)
        except:
            warnings.warn("Invalid matplotlib style... Using default style.")
            plt.style.use("default")

        for key in self.data:
            if key is not 'time':
                plt.plot(self.data['time'],self.data[key],label=key)

        plt.plot([0], [11])
        plt.legend(loc='best')


    ##https: // github.com / pyurdme / pyurdme / blob / master / pyurdme / pyurdme.py