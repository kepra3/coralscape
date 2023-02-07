#!/usr/anaconda3/bin/python
# -*- coding: utf-8 -*-
# remember to edit python environment if required.

"""
@author: kprata
@date created: 26/5/22
@description: Rotate each plot parallel
"""
import numpy as np
import pandas as pd

annotations = pd.read_csv("results/annotations.csv")
SB20 = annotations.iloc[453:561, :]
rest_of_plots = annotations.iloc[0:453, :]
SB05 = annotations.iloc[562:570, :]

SB20_dict = SB20.set_index('X').T.to_dict('list')
rest_of_plots_dict = rest_of_plots.set_index('X').T.to_dict('list')
SB05_dict = SB05.set_index('X').T.to_dict('list')

R_minus90 = np.matrix([[-0.000000043711, 1.000000000000, 0.000000000000, 0.000000000000],
                       [-1.000000000000, -0.000000043711, 0.000000000000, 0.000000000000],
                       [0.000000000000, 0.000000000000, 1.000000000000, 0.000000000000],
                       [0.000000000000, 0.000000000000, 0.000000000000, 1.000000000000]])

R_plus90 = np.matrix([[-0.000000043711, -1.000000000000, 0.000000000000, 0.000000000000],
                      [1.000000000000, -0.000000043711, 0.000000000000, 0.000000000000],
                      [0.000000000000, 0.000000000000, 1.000000000000, 0.000000000000],
                      [0.000000000000, 0.000000000000, 0.000000000000, 1.000000000000]])

parallel_rest_of_plots = {}
for name in rest_of_plots_dict:
    rest_of_plots_dict[name].append(1)
    parallel_rest_of_plots[name] = np.matmul(R_minus90, rest_of_plots_dict[name])
    parallel_rest_of_plots[name] = parallel_rest_of_plots[name].tolist()
    parallel_rest_of_plots[name] = parallel_rest_of_plots[name][0]

parallel_SB05 = {}
for name in SB05_dict:
    SB05_dict[name].append(1)
    parallel_SB05[name] = np.matmul(R_minus90, SB05_dict[name])
    parallel_SB05[name] = parallel_SB05[name].tolist()
    parallel_SB05[name] = parallel_SB05[name][0]

parallel_SB20 = {}
for name in SB20_dict:
    SB20_dict[name].append(1)
    parallel_SB20[name] = np.matmul(R_plus90, SB20_dict[name])
    parallel_SB20[name] = parallel_SB20[name].tolist()
    parallel_SB20[name] = parallel_SB20[name][0]

rest_of_plots_df = pd.DataFrame(parallel_rest_of_plots).T
rest_of_plots_df.columns = ['x', 'y', 'z', 'p']
rest_of_plots_df.to_csv(
    '/Users/kprata/Dropbox/agaricia_project_2019/shalo_ag/gen_project/3 - Spatial/rest_of_plots_parallel.txt')

SB05_df = pd.DataFrame(parallel_SB05).T
SB05_df.columns = ['x', 'y', 'z', 'p']
SB05_df.to_csv('/Users/kprata/Dropbox/agaricia_project_2019/shalo_ag/gen_project/3 - Spatial/SB05_parallel.txt')

SB20_df = pd.DataFrame(parallel_SB20).T
SB20_df.columns = ['x', 'y', 'z', 'p']
SB20_df.to_csv('/Users/kprata/Dropbox/agaricia_project_2019/shalo_ag/gen_project/3 - Spatial/SB20_parallel.txt')
