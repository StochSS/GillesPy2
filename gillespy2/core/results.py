import warnings

from collections import UserDict,UserList

def _plot_iterate(self, show_labels = True, included_species_list = []):
    import matplotlib.pyplot as plt

    for i,species in enumerate(self.data):
        if species is not 'time':

            if species not in included_species_list and included_species_list:
                continue

            line_color = 'C' + str(i%10)

            if show_labels:
                label = species
            else:
                label = ""

            plt.plot(self.data['time'], self.data[species], label=label,color = line_color)

def _plotplotyl_iterate(self, show_labels = True, trace_list = None, line_dict= None, included_species_list= []):

    # List of 50 hex color values used for plotly graphs
    common_rgb_values = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f',
                         '#bcbd22', '#17becf','#ff0000', '#00ff00', '#0000ff', '#ffff00', '#00ffff', '#ff00ff',
                         '#800000', '#808000','#008000', '#800080', '#008080', '#000080', '#ff9999', '#ffcc99',
                         '#ccff99', '#cc99ff','#ffccff', '#62666a', '#8896bb', '#77a096', '#9d5a6c', '#9d5a6c',
                         '#eabc75', '#ff9600','#885300', '#9172ad', '#a1b9c4', '#18749b', '#dadecf', '#c5b8a8',
                         '#000117', '#13a8fe','#cf0060', '#04354b', '#0297a0', '#037665', '#eed284', '#442244',
                         '#ffddee', '#702afb']

    if trace_list is None:
        trace_list = []

    import plotly.graph_objs as go

    for i,species in enumerate(self.data):
        if species is not 'time':

            if species not in included_species_list and included_species_list:
                continue

            if line_dict is None:
                line_dict = {}

            #If number of species exceeds number of available colors, loop back through colors
            line_dict['color'] = common_rgb_values[(i-1)%len(common_rgb_values)]

            if show_labels:
                trace_list.append(
                    go.Scatter(
                        x=self.data['time'],
                        y=self.data[species],
                        mode='lines',
                        name=species,
                        line = line_dict
                    )
                )
            else:
                trace_list.append(
                    go.Scatter(
                        x=self.data['time'],
                        y=self.data[species],
                        mode='lines',
                        name=species,
                        line=line_dict,
                        showlegend=False
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

    def __getitem__(self, key):
        if type(key) is type(1):
            warnings.warn("Results is of type dictionary. Use results['species'] instead of results[0]['species'] ")
            return self
        if key in self.data:
            return self.data[key]
        if hasattr(self.__class__, "__missing__"):
            return self.__class__.__missing__(self, key)
        raise KeyError(key)

    def plot(self, xaxis_label ="Time (s)", yaxis_label ="Species Population", title = None, style="default",
             show_legend=True, included_species_list=[],save_png=False,figsize = (18,10)):
        """ Plots the Results using matplotlib.

         Attributes
        ----------
        xaxis_label : str
            the label for the x-axis
        yaxis_label : str
            the label for the y-axis
        title : str
            the title of the graph
        show_legend : bool
            whether or not to display a legend which lists species
        included_species_list : list
            A list of strings describing which species to include. By default displays all species.
        save_png : bool or str
            Should the graph be saved as a png file. If True, File name is title of graph. If a string is given, file
            is named after that string.
        figsize : tuple
            the size of the graph. A tuple of the form (width,height). Is (18,10) by default.

        """
        import matplotlib.pyplot as plt

        try:
            plt.style.use(style)
        except:
            warnings.warn("Invalid matplotlib style. Try using one of the following {}".format(plt.style.available))
            plt.style.use("default")

        if title is None:
            title = (self.model.name + " - " + self.solver_name)

        plt.figure(figsize=figsize)
        plt.title(title,fontsize=18)
        plt.xlabel(xaxis_label)
        plt.ylabel(yaxis_label)

        _plot_iterate(self, included_species_list=included_species_list)

        plt.plot([0], [11])

        if show_legend:
            plt.legend(loc='best')

        if isinstance(save_png, str):
            plt.savefig(save_png)

        elif save_png:
            plt.savefig(title)


    def plotplotly(self, xaxis_label = "Time (s)", yaxis_label="Species Population", title = None, show_legend=True,
                   included_species_list=[], return_plotly_figure = False):
        """ Plots the Results using plotly. Can only be viewed in a Jupyter Notebook.

         Attributes
        ----------
        xaxis_label : str
            the label for the x-axis
        yaxis_label : str
            the label for the y-axis
        title : str
            the title of the graph
        show_legend : bool
            whether or not to display a legend which lists species
        included_species_list : list
            A list of strings describing which species to include. By default displays all species.
        return_plotly_figure : bool
            whether or not to return a figure dictionary of data(graph object traces) and layout options
            which may be edited by the user.

        """

        from plotly.offline import init_notebook_mode, iplot
        import plotly.graph_objs as go

        init_notebook_mode(connected=True)

        if title is None:
            title = (self.model.name + " - " + self.solver_name)

        trace_list = _plotplotyl_iterate(self, included_species_list = included_species_list,show_labels=True)

        layout = go.Layout(
            showlegend=show_legend,
            title= title,
            xaxis=dict(
                title=xaxis_label),
            yaxis=dict(
                title=yaxis_label)
        )
        fig = dict(data = trace_list,layout=layout)
        iplot(fig)

        if return_plotly_figure:
            return fig

class EnsembleResults(UserList):
    """ List of Results Dicts created by a gillespy2 solver with multiple trajectories, extends the UserList object.

        Attributes
        ----------
        data : UserList
            A list of Results
        """

    def __init__(self,data):
        self.data = data

    def plot(self, xaxis_label ="Time (s)", yaxis_label ="Species Population", style="default", title = None,
             show_legend=True, multiple_graphs = False, included_species_list=[],save_png=False,figsize = (18,10)):
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
        included_species_list : list
             A list of strings describing which species to include. By default displays all species.
        save_png : bool or str
            Should the graph be saved as a png file. If True, File name is title of graph. If a string is given, file
            is named after that string.
        figsize : tuple
            the size of the graph. A tuple of the form (width,height). Is (18,10) by default.


            """
        import matplotlib.pyplot as plt
        results_list = self.data

        if title is None:
            if isinstance(self[0].model.name, str):
                title = (self[0].model.name + " - " + self[0].solver_name)
            else: title=''

        if len(results_list) < 2:
                multiple_graphs = False

        if multiple_graphs:

            for i,result in enumerate(results_list):

                if isinstance(save_png, str):
                    result.plot(xaxis_label=xaxis_label, yaxis_label=yaxis_label, title=title + " " + str(i + 1), style=style,
                                                 included_species_list=included_species_list,save_png=save_png + str(i + 1),figsize=figsize)
                else:
                    result.plot(xaxis_label=xaxis_label, yaxis_label=yaxis_label, title=title + " " + str(i + 1),style=style,
                                included_species_list=included_species_list, save_png=save_png, figsize=figsize)

        else:
            try:
                plt.style.use(style)
            except:
                warnings.warn("Invalid matplotlib style. Try using one of the following {}".format(plt.style.available))
                plt.style.use("default")

            plt.figure(figsize=figsize)
            plt.title(title, fontsize=18)
            plt.xlabel(xaxis_label)
            plt.ylabel(yaxis_label)

            for i,result in enumerate(results_list):

                if i > 0:
                    _plot_iterate(result, included_species_list=included_species_list,show_labels=False)
                else:
                    _plot_iterate(result, included_species_list=included_species_list)

            if show_legend:
                plt.legend(loc='best')
            plt.plot([0], [11])

            if isinstance(save_png, str):
                plt.savefig(save_png)

            elif save_png:
                plt.savefig(title)

    def plotplotly(self, xaxis_label = "Time (s)", yaxis_label="Species Population", title = None, show_legend=True,
                   multiple_graphs = False, included_species_list=[],return_plotly_figure=False):
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
        included_species_list : list
             A list of strings describing which species to include. By default displays all species.
        return_plotly_figure : bool
            whether or not to return a figure dictionary of data(graph object traces) and layout options
            which may be edited by the user.
        """

        from plotly.offline import init_notebook_mode, iplot
        import plotly.graph_objs as go

        init_notebook_mode(connected=True)

        results_list = self.data
        number_of_trajectories =len(results_list)

        if title is None:
            title = (self[0].model.name + " - " + self[0].solver_name)

        fig = dict(data=[], layout=[])

        if len(results_list) < 2:
            multiple_graphs = False

        if multiple_graphs:

            from plotly import tools

            fig = tools.make_subplots(print_grid=False,rows=int(number_of_trajectories/2) + int(number_of_trajectories%2),
                                      cols = 2)

            for i, result in enumerate(results_list):
                if i > 0:
                    trace_list = _plotplotyl_iterate(result, trace_list=[], included_species_list= included_species_list,
                                                     show_labels=False)
                else:
                    trace_list = _plotplotyl_iterate(result, trace_list=[], included_species_list=included_species_list)

                for k in range(0,len(trace_list)):
                    if i%2 == 0:
                        fig.append_trace(trace_list[k], int(i/2) + 1, 1)
                    else:
                        fig.append_trace(trace_list[k], int(i/2) + 1, 2)

                fig['layout'].update(autosize=True,
                                     height=400*len(results_list),
                                     showlegend=show_legend,title =title)

            iplot(fig)

        else:
            trace_list = []
            for i,result in enumerate(results_list):
                if i > 0:
                    trace_list = _plotplotyl_iterate(result, trace_list=trace_list,included_species_list= included_species_list,
                                                     show_labels = False)
                else:
                    trace_list = _plotplotyl_iterate(result, trace_list=trace_list,included_species_list= included_species_list)

            layout = go.Layout(
                showlegend=show_legend,
                title=title,
                xaxis=dict(
                    title=xaxis_label),
                yaxis=dict(
                    title=yaxis_label)
            )

            fig['data'] = trace_list
            fig['layout'] = layout
            iplot(fig)

        if return_plotly_figure:
            return fig


    def average_ensemble(self):
        """
                Generate a single Results dictionary that is made of the means of all trajectories' outputs
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

    def stddev_ensemble(self,ddof = 0):
        """
                Generate a single Results dictionary that is made of the sample standard deviations of all trajectories'
                outputs.

                  Attributes
                ----------
                ddof : int
                    Delta Degrees of Freedom. The divisor used in calculations is N - ddof, where N represents
                    the number of trajectories. Sample standard deviation uses ddof of 1. Defaults to population
                    standard deviation where ddof is 0.

                :return: the Results dictionary
                """

        from math import sqrt

        results_list = self.data
        number_of_trajectories = len(results_list)

        if ddof == number_of_trajectories:
            warnings.warn("ddof must be less than the number of trajectories. Using ddof of 0")
            ddof = 0

        average_list = self.average_ensemble()

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
                output[species][i] /= (number_of_trajectories - ddof)
                output[species][i] = sqrt(output[species][i])

        return output

    def plotplotly_std_dev_range(self, xaxis_label = "Time (s)", yaxis_label="Species Population", title = None,
                                 show_legend=True, included_species_list = [],return_plotly_figure=False,ddof = 0):
        """
           Plot a plotly graph depicting standard deviation and the mean graph of an ensemble_results object

         Attributes
        ----------
        xaxis_label : str
            the label for the x-axis
        yaxis_label : str
            the label for the y-axis
        title : str
            the title of the graph
        show_legend : bool
            whether or not to display a legend which lists species
        included_species_list : list
            A list of strings describing which species to include. By default displays all species.
        return_plotly_figure : bool
            whether or not to return a figure dictionary of data(graph object traces) and layout options
            which may be edited by the user.
        ddof : int
            Delta Degrees of Freedom. The divisor used in calculations is N - ddof, where N represents
            the number of trajectories. Sample standard deviation uses ddof of 1. Defaults to population
            standard deviation where ddof is 0.

        """

        average_result = self.average_ensemble()
        stddev_result = self.stddev_ensemble(ddof= ddof)

        from plotly.offline import init_notebook_mode, iplot
        import plotly.graph_objs as go

        init_notebook_mode(connected=True)

        if title is None:
            title = (average_result.model.name + " - " + average_result.solver_name + " - Standard Deviation Range")

        trace_list=[]
        for species in average_result:
            if species is not 'time':

                if species not in included_species_list and included_species_list:
                    continue

                upper_bound = []
                for i in range(0, len(average_result[species])):
                    upper_bound.append(average_result[species][i] + stddev_result[species][i])

                trace_list.append(
                    go.Scatter(
                        name=species+ ' Upper Bound',
                        x=average_result['time'],
                        y = upper_bound,
                        mode='lines',
                        marker=dict(color="#444"),
                        line=dict(width=1,dash='dot'),
                        legendgroup="Standard Deviation",
                        showlegend=False
                    )
                )
                trace_list.append(
                    go.Scatter(
                        x=average_result['time'],
                        y=average_result[species],
                        name=species,
                        fillcolor='rgba(68, 68, 68, 0.2)',
                        fill='tonexty'
                    )
                )

                lower_bound = []
                for i in range(0, len(average_result[species])):
                    lower_bound.append(average_result[species][i] - stddev_result[species][i])

                trace_list.append(
                    go.Scatter(
                        name=species + ' Lower Bound',
                        x=average_result['time'],
                        y= lower_bound,
                        mode='lines',
                        marker=dict(color="#444"),
                        line=dict(width=1,dash='dot'),
                        fillcolor='rgba(68, 68, 68, 0.2)',
                        fill='tonexty',
                        legendgroup="Standard Deviation",
                        showlegend=False
                    )
                )
        layout = go.Layout(
            showlegend=show_legend,
            title=title,
            xaxis=dict(
                title=xaxis_label),
            yaxis=dict(
                title=yaxis_label)
        )
        fig = dict(data=trace_list, layout=layout)
        iplot(fig)

        if return_plotly_figure:
            return fig

    def plot_std_dev_range(self, xaxis_label ="Time (s)", yaxis_label ="Species Population", title = None,
                           style="default", show_legend=True, included_species_list=[],ddof=0,save_png = False,figsize = (18,10)):
        """
            Plot a matplotlib graph depicting standard deviation and the mean graph of an ensemble_results object

         Attributes
        ----------
        xaxis_label : str
            the label for the x-axis
        yaxis_label : str
            the label for the y-axis
        title : str
            the title of the graph
        show_legend : bool
            whether or not to display a legend which lists species
        included_species_list : list
            A list of strings describing which species to include. By default displays all species.
        ddof : int
            Delta Degrees of Freedom. The divisor used in calculations is N - ddof, where N represents
            the number of trajectories. Sample standard deviation uses ddof of 1. Defaults to population
            standard deviation where ddof is 0.
        save_png : bool or str
            Should the graph be saved as a png file. If True, File name is title of graph. If a string is given, file
            is named after that string.
        figsize : tuple
            the size of the graph. A tuple of the form (width,height). Is (18,10) by default.

        """

        average_result = self.average_ensemble()
        stddev_result = self.stddev_ensemble(ddof=ddof)

        import matplotlib.pyplot as plt

        try:
            plt.style.use(style)
        except:
            warnings.warn("Invalid matplotlib style. Try using one of the following {}".format(plt.style.available))
            plt.style.use("default")

        plt.figure(figsize=figsize)

        for species in average_result:
            if species is 'time':
                continue

            if species not in included_species_list and included_species_list:
                continue

            lowerBound = [a-b for a,b in zip(average_result[species], stddev_result[species])]
            upperBound = [a+b for a,b in zip(average_result[species], stddev_result[species])]

            plt.fill_between(average_result['time'], lowerBound, upperBound,color='whitesmoke')
            plt.plot(average_result['time'],lowerBound,upperBound,color='grey',linestyle='dashed')
            plt.plot(average_result['time'],average_result[species],label=species)

        if title is None:
            title = (average_result.model.name + " - " + average_result.solver_name + " - Standard Deviation Range")

        plt.title(title, fontsize=18)
        plt.xlabel(xaxis_label)
        plt.ylabel(yaxis_label)
        plt.plot([0], [11])
        if show_legend:
            plt.legend(loc='best')

        if isinstance(save_png, str):
            plt.savefig(save_png)

        elif save_png:
            plt.savefig(title)

        
