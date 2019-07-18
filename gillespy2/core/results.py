import warnings

from collections import UserDict,UserList

def _plot_iterate(self,num = ""):
    import matplotlib.pyplot as plt

    for species in self.data:
        if species is not 'time':
            plt.plot(self.data['time'], self.data[species], label=(species + num))

def _plotplotyl_iterate(self, num = "", trace_list = None):

    if trace_list is None:
        trace_list = []

    import plotly.graph_objs as go

    for species in self.data:
        if species is not 'time':
            trace_list.append(
                go.Scatter(
                    x=self.data['time'],
                    y=self.data[species],
                    mode='lines',
                    name=species + " " + num
                )
            )
    return trace_list

class Results(UserDict):
    """ Results Dict created by a gillespy2 solver with single trajectory, extends the UserDict object.

        Attributes
        ----------
        data : UserList
            A list of Results that are created by solvers with multiple trajectories
        """

    def __init__(self,data,model = None,solver_name = "Undefined solver name"):

        self.data = data
        self.model = model
        self.solver_name = solver_name

    def plot(self, xaxis_label ="Time (s)", yaxis_label ="Species Population", title = "default", style="default",**kwargs):
        """ Plots the Results using matplotlib.

         Attributes
        ----------
        xaxis_label : str
            the label for the x-axis
        yaxis_label : str
            the label for the y-axis
        title : str
            the title of the graph
        """
        import matplotlib.pyplot as plt

        try:
            plt.style.use(style)
        except:
            warnings.warn("Invalid matplotlib style. Try using one of the following {}".format(plt.style.available))
            plt.style.use("default")

        if title is "default":
            title = (self.model.name + " - " + self.solver_name)

        plt.figure(figsize=(18, 10))
        plt.title(title,fontsize=18)
        plt.xlabel(xaxis_label)
        plt.ylabel(yaxis_label)

        _plot_iterate(self)

        plt.plot([0], [11])
        plt.legend(loc='best')

    def plotplotly(self,xaxis_label = "Time (s)", yaxis_label="Species Population",title = "default",**kwargs):
        """ Plots the Results using plotly. Can only be viewed in a Jupyter Notebook.

         Attributes
        ----------
        xaxis_label : str
            the label for the x-axis
        yaxis_label : str
            the label for the y-axis
        title : str
            the title of the graph
        """

        from plotly.offline import init_notebook_mode, iplot
        import plotly.graph_objs as go

        init_notebook_mode(connected=True)

        if title is "default":
            title = (self.model.name + " - " + self.solver_name)

        trace_list = _plotplotyl_iterate(self)

        layout = go.Layout(
            title= title,
            xaxis=dict(
                title=xaxis_label),
            yaxis=dict(
                title=yaxis_label)
        )
        fig = dict(data = trace_list,layout=layout)
        iplot(fig)

class EnsembleResults(UserList):
    """ List of Results Dicts created by a gillespy2 solver with multiple trajectories, extends the UserList object.

        Attributes
        ----------
        data : UserList
            A list of Results
        """

    def __init__(self,data):
        self.data = data

    def plot(self, xaxis_label ="Time (s)", yaxis_label ="Species Population", style="default", title = "default", multiple_graphs = False,**kwargs):
        """ Plots the Results using matplotlib.

        Attributes
        ----------
        xaxis_label : str
            the label for the x-axis
        yaxis_label : str
            the label for the y-axis
            style : str
            the matplotlib style to be used for the graph or graphs
        title : str
            the title of the graph
        multiple_graphs : bool
            if each trajectory should have its own graph or if they should overlap
            """
        import matplotlib.pyplot as plt
        results_list = self.data

        if title is "default":
            title = (self[0].model.name + " - " + self[0].solver_name)

        if len(results_list) < 2:
                multiple_graphs = False

        if multiple_graphs is True:
            for i,result in enumerate(results_list):
                result.plot(xaxis_label=xaxis_label,yaxis_label=yaxis_label,title=title + " " + str(i + 1),style=style)

        if multiple_graphs is False:
            try:
                plt.style.use(style)
            except:
                warnings.warn("Invalid matplotlib style. Try using one of the following {}".format(plt.style.available))
                plt.style.use("default")

            plt.figure(figsize=(18, 10))
            plt.title(title, fontsize=18)
            plt.xlabel(xaxis_label)
            plt.ylabel(yaxis_label)

            for i,result in enumerate(results_list):
                _plot_iterate(result,num = (" " + str(i + 1)))

            plt.plot([0], [11])
            plt.legend(loc='best')


    def plotplotly(self,xaxis_label = "Time (s)", yaxis_label="Species Population",title = "default",multiple_graphs = False,**kwargs):
        """ Plots the Results using plotly. Can only be viewed in a Jupyter Notebook.

        Attributes
        ----------
        xaxis_label : str
            the label for the x-axis
        yaxis_label : str
            the label for the y-axis
        title : str
            the title of the graph
        multiple_graphs : bool
            if each trajectory should have its own graph or if they should overlap

        """

        from plotly.offline import init_notebook_mode, iplot
        import plotly.graph_objs as go

        init_notebook_mode(connected=True)

        results_list = self.data
        number_of_trajectories =len(results_list)

        if title is "default":
            title = (self[0].model.name + " - " + self[0].solver_name)

        fig = dict(data=[], layout=[])

        if len(results_list) < 2:
            multiple_graphs = False

        if multiple_graphs is True:

            from plotly import tools

            fig = tools.make_subplots(print_grid=False,rows=int(number_of_trajectories/2) + int(number_of_trajectories%2),cols = 2)

            for i, result in enumerate(results_list):
                trace_list = _plotplotyl_iterate(result, num = str(i + 1), trace_list=[])
                for k in range(0,len(trace_list)):
                    if i%2 == 0:
                        fig.append_trace(trace_list[k], int(i/2) + 1, 1)
                    else:
                        fig.append_trace(trace_list[k], int(i/2) + 1, 2)

                fig['layout'].update(autosize=True,
                                     height=400*len(results_list),
                                     showlegend=True,title =title)

            iplot(fig)

        else:
            trace_list = []
            for i,result in enumerate(results_list):
                trace_list = _plotplotyl_iterate(result, num = str(i + 1), trace_list=trace_list)

            layout = go.Layout(
                title=title,
                xaxis=dict(
                    title=xaxis_label),
                yaxis=dict(
                    title=yaxis_label)
            )

            fig['data'] = trace_list
            fig['layout'] = layout
            iplot(fig)

    def average_ensemble(self):
        """
                Generate a single Results dictionary that is made of the mean of all trajectories' outputs
                :return: the Results dictionary
                """

        results_list = self.data
        number_of_trajectories = len(results_list)

        output = Results(data={},model=results_list[0].model,solver_name=results_list[0].solver_name)

        for species in results_list[0]: #Initialize the output to be the same size as the inputs
            output[species] = [0]*len(results_list[0][species])

        output['time'] = results_list[0]['time']

        for i in range(0,number_of_trajectories): #Add every value of every Results Dict into one output Results
            results_dict = results_list[i]
            for species in results_dict:
                if species is 'time':
                    continue
                for k in range(0,len(output[species])):
                    output[species][k] += results_dict[species][k]

        for species in output:   #Divide for mean of every value in output Results
            if species is 'time':
                continue
            for i in range(0,len(output[species])):
                output[species][i] /= number_of_trajectories

        return output

    def stddev_ensemble(self):
        """
                Generate a single Results dictionary that is made of the samplestandard deviation of all trajectories'
                outputs.
                :return: the Results dictionary
                """

        import math

        results_list = self.data
        average_list = self.average_ensemble()
        number_of_trajectories = len(results_list)
        output = Results(data={}, model=results_list[0].model, solver_name=results_list[0].solver_name)

        for species in results_list[0]: #Initialize the output to be the same size as the inputs
            output[species] = [0]*len(results_list[0][species])

        output['time'] = results_list[0]['time']

        for i in range(0,number_of_trajectories):
            results_dict = results_list[i]
            for species in results_dict:
                if species is 'time':
                    continue
                for k in range(0,len(output[species])):
                    output[species][k] += (results_dict[species][k] - average_list[species][k])\
                                          *(results_dict[species][k] - average_list[species][k])

        for species in output:   #Divide for mean of every value in output Results
            if species is 'time':
                continue
            for i in range(0,len(output[species])):
                output[species][i] /= number_of_trajectories
                output[species][i] = math.sqrt(output[species][i])

        return output



        #TODO Add more funcionality to ensemble_results for better graphing, statistical analysis, etc.