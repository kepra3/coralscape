#!/usr/anaconda3/bin/python
# -*- coding: utf-8 -*-
# remember to edit python environment if required.

"""
@author: kprata
@date created: 27/1/22
@description: Script for everyone for visualising plots and points
"""
import open3d as o3d
import numpy as np  # very useful package that can basically do everything!
# import json
# import copy

PATH = "/Users/kprata/Dropbox/agaricia_project_2019/shalo_ag/Photogrammetry/CloudCompare/"
# V_DISTANCE = -10
ply_filename = "cur_kal_40m_20200214_dec5M.ply"
#annotations_filename = "cur_kal_20m_20200214_decvis_02_KP_16-12-21_completed.txt"
radius = 2  # choose a radius around your colony
# subsets_filename = "subsets.json" if you have an upvector!!
# short_name = "_".join(ply_filename.split('_')[0:4])

print('Reading PLY file ...')
pcd = o3d.io.read_point_cloud('{}/{}'.format(PATH, ply_filename))

print('Build KDTree from point cloud ...')
pcd_tree = o3d.geometry.KDTreeFlann(pcd)

o3d.visualization.draw_geometries([pcd])

print('Read assignment file ...')
annotations_path = '{}/{}'.format(PATH, annotations_filename)
annotations = {}
annotations_file = open(annotations_path, 'r')
# Create a dictionary where each 'key' is the sample and the values are the x, y, z of the annotated point
for line in annotations_file:
    if any(flag in line for flag in ['left', 'right', 'X']):  # don't have to include this in yours if you haven't annotated like this
        continue
    cols = line.rstrip().replace(',', ' ').split()
    annotations[cols[0]] = [float(i) for i in cols[1:4]]
annotations_file.close()


# print('Read viscore file ...')
# subsets = json.load(open('{}/{}'.format(PATH, subsets_filename)))
# Determine json key of primary model
#if '{0}/{0}'.format(short_name) in subsets['d']:
#    subsets_ortho = subsets['d']['{0}/{0}'.format(short_name)]['c']['ortho']
#elif self.short_name in subsets['d']:
#    subsets_ortho = subsets['d'][short_name]['c']['ortho']
#else:
#    print('Model not found in subsets.json!')
#dd = subsets_ortho['dd']
#scale_factor = subsets_ortho['scale_factor']
#r = subsets_ortho['vecs']['r']
#u = subsets_ortho['vecs']['u']
#n = subsets_ortho['vecs']['n']
#c = subsets_ortho['vecs']['c']
#cc = subsets_ortho['vecs']['cc']
#cam_up = subsets_ortho['vecs']['cam']['up']
#cam_eye = subsets_ortho['vecs']['cam']['eye']
#cam_target = subsets_ortho['vecs']['cam']['target']

# Create an empty dictionaries and lists
colonies = {}
# connect_points = []
# connect_colors = []
print('Searching radius around annotations ...')
for name in annotations:
    [k, idx, _] = pcd_tree.search_radius_vector_3d(annotations[name], radius=radius)
    # Store colony as separate point cloud
    colonies[name] = o3d.geometry.PointCloud()
    colonies[name].points = o3d.utility.Vector3dVector(np.asarray(pcd.points)[idx[1:], :])
    random_color = list(np.random.choice(np.linspace(0, 1, 255), size=3))
    colonies[name].paint_uniform_color(random_color)
    # if you have an upvector!
    # connect_colors.append(random_color)
    # Vertical offset along up vector
    # v_translation = np.array(dd[0:3]) * V_DISTANCE
    # colonies[name] = copy.deepcopy(colonies[name])
    # colonies[name].translate(v_translation)
    # Store points for the connecting lines
    # connect_points.append(annotations[name])
    # connect_points.append(annotations[name] + v_translation)
    # Create connected lineset
# connecting_lineset = generate_connecting_lineset(connect_points, connect_colors)
# Join all geometries
all_geoms = list(colonies.values())
all_geoms.append(pcd)
# all_geoms.append(connecting_lineset)
# Visualise
o3d.visualization.draw_geometries(all_geoms)#,
                                  #zoom=0.4,
                                  #front=cam_eye,
                                  #lookat=cam_target,
                                  #up=cam_up)