#!/usr/anaconda3/bin/python
# -*- coding: utf-8 -*-
# remember to edit python environment if required.

"""
@author: kprata
@date created: 3/2/22
@description: Script to test overhang directly
"""
from scripts.struc_complex import *

ply_filename = "cur_kal_20m_20200214_decvis_02.ply"
annotations_filename = "cur_kal_20m_20200214_decvis_02_KP_16-12-21_completed.txt"
subsets_filename = "subsets.json"

short_name = "_".join(ply_filename.split('_')[0:4])

print('Reading PLY file ...')
pcd = o3d.io.read_point_cloud('{}/{}'.format(PATH, ply_filename))

print('Read viscore metadata file ...')
viscore_md = Viscore_metadata(subsets_filename, short_name)

print('Rotating matrix ...')
up_vector = viscore_md.dd[0:3]
R = rotate_matrix(pcd, up_vector)
pcd_r = copy.deepcopy(pcd)  # delete pcd?
pcd_r.rotate(R, center=(0, 0, 0))

print('Read assignment file ...')
annotations = get_annotations('{}/{}'.format(PATH, annotations_filename))

print('Rotate annotations ...')
rotated_annotations = {}
for name in annotations:
    rotated_annotations[name] = np.matmul(R, annotations[name])

print('Get ranges for each sample ...')
ranges = get_ranges('{}/{}'.format(PATH, annotations_filename), annotations)

print("Searching for colony around annotations ...")
rotated_colonies = get_neighbourhood(rotated_annotations, pcd_r, pcd_tree_r, ranges)
print("Found colony!")

# Define neighbourhood range
double_range = {}
for name in ranges:
    double_range[name] = ranges[name] * 2

print('Getting environment')
rotated_environment = get_neighbourhood(rotated_annotations, pcd_r, pcd_tree_r, double_range)
print('Found environment!')

overhang_value = 0.01 / viscore_md.scale_factor
overhang = calc_overhang(rotated_annotations, overhang_value, rotated_environment, ranges)