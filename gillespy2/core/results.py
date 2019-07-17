import warnings

from collections import UserDict,UserList

def _plot_iterate(self,num = ""):
    import matplotlib.pyplot as plt

    for key in self.data:
        if key is not 'time':
            plt.plot(self.data['time'], self.data[key], label=(key + num))

def _plotplotyl_iterate(self,num = "",outputdata = []):

    import plotly.graph_objs as go

    for key in self.data:
        if key is not 'time':
            outputdata.append(
                go.Scatter(
                    x=self.data['time'],
                    y=self.data[key],
                    mode='lines',
                    name=key + " " + num
                )
            )
    return outputdata


class Results(UserDict):
    """ Result object for a gillespy2 solver, extends the UserDict object. """

    def __init__(self,data,model = None,solver_name = "Undefined solver name"):
        self.data = data
        self.model = model
        self.solver_name = solver_name

    def plot(self, xaxis_label ="Time (s)", yaxis_label ="Species Population", title = "default", style="default"):

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

    def plotplotly(self,xaxis_label = "Time (s)", yaxis_label="Species Population",title = "default"):

        from plotly.offline import init_notebook_mode, iplot
        import plotly.graph_objs as go

        init_notebook_mode(connected=True)

        if title is "default":
            title = (self.model.name + " - " + self.solver_name)

        outputdata = _plotplotyl_iterate(self)

        layout = go.Layout(
            title= title,
            xaxis=dict(
                title=xaxis_label),
            yaxis=dict(
                title=yaxis_label)
        )
        fig = dict(data = outputdata,layout=layout)
        iplot(fig)

class EnsembleResults(UserList):
    """ Result object for a gillespy2 solver, extends the UserList object. """

    def __init__(self,data):
        self.data = data

    def plot(self, xaxis_label ="Time (s)", yaxis_label ="Species Population", style="default", title = "default", multiple_graphs = False):

        import matplotlib.pyplot as plt

        if title is "default":
            title = (self[0].model.name + " - " + self[0].solver_name)

        if multiple_graphs is True:
            for i,result in enumerate(self.data):
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

            for i,result in enumerate(self.data):
                _plot_iterate(result,num = (" " + str(i + 1)))

            plt.plot([0], [11])
            plt.legend(loc='best')


    def plotplotly(self,xaxis_label = "Time (s)", yaxis_label="Species Population",title = "default",multiple_graphs = False):
        from plotly.offline import init_notebook_mode, iplot
        import plotly.graph_objs as go

        init_notebook_mode(connected=True)

        if title is "default":
            title = (self[0].model.name + " - " + self[0].solver_name)

        outputdata = []
        fig = dict(data=[], layout=[])

        if multiple_graphs is True:

            from plotly import tools

            fig = tools.make_subplots(rows=len(self.data))

            for i, result in enumerate(self.data):
                out = _plotplotyl_iterate(result,num = str(i + 1))
                for b in range(0,len(out)):
                    fig.append_trace(out[b], i + 1, 1)


            iplot(fig)

        else:
            for i,result in enumerate(self.data):
                outputdata = _plotplotyl_iterate(result,num = str(i + 1),outputdata=outputdata)


            layout = go.Layout(
                title=title,
                xaxis=dict(
                    title=xaxis_label),
                yaxis=dict(
                    title=yaxis_label)
            )

            fig['data'] = outputdata
            fig['layout'] = layout
            iplot(fig)


    def average_ensemble(self):

        output = Results(data={},model=self.data[0].model,solver_name=self.data[0].solver_name)

        for key in self.data[0]:
            output[key] = [0]*len(self.data[0][key])

        output['time'] = self.data[0]['time']

        for a in range(0,len(self.data)):
            results_dict = self.data[a]
            for key in results_dict:
                if key is 'time':
                    continue
                for i in range(0,len(output[key])):
                    output[key][i] += results_dict[key][i]

        for key in output:
            if key is 'time':
                continue
            for i in range(0,len(output[key])):
                output[key][i] /= len(self.data)

        return output

        #TODO Add more funcionality to ensemble_results for better graphing, statistical analysis, etc.