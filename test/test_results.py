# GillesPy2 is a modeling toolkit for biochemical simulation.
# Copyright (C) 2019-2022 GillesPy2 developers.

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import unittest
import os
import tempfile
from gillespy2.core import Model
from gillespy2.core.results import Results, Trajectory

class TestResults(unittest.TestCase):

    def test_xscale_plot(self):
        from unittest import mock
        trajectory = Trajectory(data={'time':[0.],'foo':[1.]}, model=Model('test_model'))
        results = Results(data=[trajectory])
        with mock.patch('matplotlib.pyplot.xscale') as mock_xscale:
            results.plot(xscale='log')
        mock_xscale.assert_called_with('log')

    def test_yscale_plot(self):
        from unittest import mock
        trajectory = Trajectory(data={'time':[0.],'foo':[1.]}, model=Model('test_model'))
        results = Results(data=[trajectory])
        with mock.patch('matplotlib.pyplot.yscale') as mock_yscale:
            results.plot(yscale='log')
        mock_yscale.assert_called_with('log')

    def test_xscale_plot_std_dev_range(self):
        from unittest import mock
        trajectory = Trajectory(data={'time':[0.],'foo':[1.]}, model=Model('test_model'))
        results = Results(data=[trajectory])
        with mock.patch('matplotlib.pyplot.xscale') as mock_xscale:
            results.plot_mean_stdev(xscale='log')
        mock_xscale.assert_called_with('log')

    def test_yscale_plot_std_dev_range(self):
        from unittest import mock
        trajectory = Trajectory(data={'time':[0.],'foo':[1.]}, model=Model('test_model'))
        results = Results(data=[trajectory])
        with mock.patch('matplotlib.pyplot.yscale') as mock_yscale:
            results.plot(yscale='log')
        mock_yscale.assert_called_with('log')

    def test_plotly_layout_args(self):
        from unittest import mock
        trajectory = Trajectory(data={'time':[0.],'foo':[1.]}, model=Model('test_model'))
        results = Results(data=[trajectory])
        with mock.patch('plotly.offline.init_notebook_mode') as mock_notebook:
            with mock.patch('plotly.graph_objs.Scatter') as mock_scatter:
                with mock.patch('plotly.graph_objs.Layout') as mock_layout:
                    results.plotplotly(return_plotly_figure=True,xscale='log')
        mock_layout.assert_called_with(showlegend=True, title='', xaxis_title='Time', xscale='log', yaxis_title='Value')

    def test_plotly_std_dev_range_layout_args(self):
        from unittest import mock
        from unittest.mock import MagicMock
        trajectory = Trajectory(data={'time':[0.],'foo':[1.]}, model=Model('test_model'))
        results = Results(data=[trajectory])
        _plotly_iterate=MagicMock()
        with mock.patch('plotly.offline.init_notebook_mode') as mock_notebook:
            with mock.patch('plotly.graph_objs.Scatter') as mock_scatter:
                with mock.patch('plotly.graph_objs.Layout') as mock_layout:
                    results.plotplotly_mean_stdev(return_plotly_figure=True, xscale='log')
        mock_layout.assert_called_with(legend={'traceorder':'normal'},showlegend=True, title='Mean and Standard Deviation',
                                       xaxis_title='Time', xscale='log', yaxis_title='Value')

    def test_pickle_stable_plot_iterate(self):
        from unittest import mock
        from gillespy2.core.results import _plot_iterate
        trajectory = Trajectory(data={'time':[0.],'foo':[1.]}, model=Model('test_model'))
        import pickle
        trajectory_unpickled = pickle.loads(pickle.dumps(trajectory))
        import matplotlib
        with mock.patch('matplotlib.pyplot.plot') as mock_method_before_pickle:
            _plot_iterate(trajectory)
        with mock.patch('matplotlib.pyplot.plot') as mock_method_after_pickle:
            _plot_iterate(trajectory_unpickled)
        assert mock_method_before_pickle.call_args_list == mock_method_after_pickle.call_args_list

    def test_pickle_stable_plotplotly_iterate(self):
        from unittest import mock
        from gillespy2.core.results import _plotplotly_iterate
        trajectory = Trajectory(data={'time':[0.],'foo':[1.]}, model=Model('test_model'))
        import pickle
        trajectory_unpickled = pickle.loads(pickle.dumps(trajectory))
        with mock.patch('plotly.graph_objs.Scatter') as mock_method_before_pickle:
            _plotplotly_iterate(trajectory)
        with mock.patch('plotly.graph_objs.Scatter') as mock_method_after_pickle:
            _plotplotly_iterate(trajectory_unpickled)
        assert mock_method_before_pickle.call_args_list == mock_method_after_pickle.call_args_list

    def test_to_csv_single_result_no_data(self):
        result = Results(data=None)
        test_nametag = "test_nametag"
        test_stamp = "test_stamp"

        with tempfile.TemporaryDirectory() as tempdir:
            result.to_csv(stamp = test_stamp, nametag = test_nametag, path=tempdir)
            test_path = tempdir+"/"+test_nametag+test_stamp
            assert not os.path.isdir(test_path)
            
    def test_to_csv_single_result_directory_exists(self):
        test_data = {'time':[0]}
        result = Results(data=[test_data])
        test_nametag = "test_nametag"
        test_stamp = "test_stamp"

        with tempfile.TemporaryDirectory() as tempdir:
            result.to_csv(stamp = test_stamp, nametag = test_nametag, path=tempdir)
            test_path = tempdir+"/"+test_nametag+test_stamp
            assert os.path.isdir(test_path)

    def test_to_csv_single_result_file_exists(self):
        test_data = {'time':[0]}
        result = Results(data=[test_data])
        test_nametag = "test_nametag"
        test_stamp = "test_stamp"

        with tempfile.TemporaryDirectory() as tempdir:
            result.to_csv(stamp = test_stamp, nametag = test_nametag, path=tempdir)
            test_path = tempdir+"/"+test_nametag+test_stamp+"/"+test_nametag+"0.csv"
            assert os.path.isfile(test_path)

    def test_to_csv_single_result_no_stamp(self):
        test_data = {'time':[0]}
        result = Results(data=[test_data])
        test_nametag = "test_nametag"

        with tempfile.TemporaryDirectory() as tempdir:
            result.to_csv(nametag=test_nametag, path=tempdir)
            assert len(os.listdir(tempdir)) != 0

    def test_to_csv_single_result_no_nametag(self):
        test_model = Model('test_model')
        test_data = Trajectory(data={'time':[0]},model=test_model)

        with tempfile.TemporaryDirectory() as tempdir:
            result.to_csv(nametag=test_nametag, path=tempdir)
            assert len(os.listdir(tempdir)) != 0

    def test_to_csv_single_result_no_nametag(self):
        test_model = Model('test_model')
        test_data = Trajectory(data={'time':[0]},model=test_model)
        result = Results(data=[test_data])
        result.solver_name = 'test_solver'
        test_stamp = "test_stamp"

        with tempfile.TemporaryDirectory() as tempdir:
            result.to_csv(stamp=test_stamp, path=tempdir)
            assert len(os.listdir(tempdir)) != 0

    def test_to_csv_single_result_no_path(self):
        test_data = Trajectory({'time':[0]},model=Model('test_model'),solver_name='test_solver_name')
        result = Results(data=[test_data])
        test_nametag = "test_nametag"
        test_stamp = "test_stamp"

