import warnings

from collections import UserDict,UserList

class results(UserDict):
    """ Result object for a gillespy2 solver, extends the UserDict object. """

    def __init__(self,data,model = None,solver_name = "Undefined solver name"):
        self.data = data
        self.model = model
        self.solver_name = solver_name

    def plot(self, xaxis_label ="Time (s)", yaxis_label ="Species Population", style="default"):

        import matplotlib.pyplot as plt

        try:
            plt.style.use(style)
        except:
            warnings.warn("Invalid matplotlib style. Try using one of the following {}".format(plt.style.available))
            plt.style.use("default")

        plt.figure(figsize=(18, 10))
        plt.title((self.model.name+ " - " +self.solver_name),fontsize=18)
        plt.xlabel(xaxis_label)
        plt.ylabel(yaxis_label)

        for key in self.data:
            if key is not 'time':
                plt.plot(self.data['time'],self.data[key],label=key)

        plt.plot([0], [11])
        plt.legend(loc='best')

    def plotplotly(self,xaxis_label = "Time (s)", yaxis_label="Species Population"):

        from plotly.offline import init_notebook_mode, iplot
        import plotly.graph_objs as go

        init_notebook_mode(connected=True)

        outputdata = []
        for key in self.data:
            if key is not 'time':
                outputdata.append(
                    go.Scatter(
                        x=self.data['time'],
                        y=self.data[key],
                        mode = 'lines',
                        name = key
                    )
                )
        layout = go.Layout(
            title= (self.model.name + " - " + self.solver_name),
            xaxis=dict(
                title=xaxis_label),
            yaxis=dict(
                title=yaxis_label)
        )
        fig = dict(data = outputdata,layout=layout)
        iplot(fig)

class ensemble_results(UserList):
    """ Result object for a gillespy2 solver, extends the UserList object. """

    def __init__(self,data):
        self.data = data

    def plot(self,multiplegraphs = True):

        if multiplegraphs is True:
            for result in self.data:
                result.plot()

        if multiplegraphs is False:
            print("Graph merging is not complete yet!")

        #TODO Add more funcionality to ensemble_results for better graphing, statistical analysis, etc.