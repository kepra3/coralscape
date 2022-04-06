#!/usr/anaconda3/bin/python
# -*- coding: utf-8 -*-
# remember to edit python environment if required.

"""
@author: kprata
@date created: 19/3/22
@description: Take plot point and annotation points and export environment points for each sample.
"""

# Modules
import argparse
import numpy as np
import open3d as o3d
import copy
import json

IGNORE_ANNOTATIONS = ['left', 'right', 'X']


def get_annotations(annotations_path):
    """ Read annotations from csv file """
    annotations = {}
    annotations_file = open(annotations_path, 'r')
    for line in annotations_file:
        if any(flag in line for flag in IGNORE_ANNOTATIONS):
            continue
        cols = line.rstrip().replace(',', ' ').split()
        if cols[0] == 'x':
            continue
        else:
            annotations[cols[0]] = [float(i) for i in cols[1:5]]
    annotations_file.close()
    return annotations


def get_neighbourhood(annotations, pcd, pcd_tree, environment_distance):
    """ Find neighbouring points """
    colonies = {}
    print('Searching ...')
    for name in annotations:
        [k, idx, _] = pcd_tree.search_radius_vector_3d(annotations[name][0:3],
                                                       annotations[name][3] + environment_distance)
        # Store colony as separate point cloud with points, normals and colours!!!
        colonies[name] = o3d.geometry.PointCloud()
        colonies[name].points = o3d.utility.Vector3dVector(np.asarray(pcd.points)[idx[1:], :])
        colonies[name].normals = o3d.utility.Vector3dVector(np.asarray(pcd.normals)[idx[1:], :])
        colonies[name].colors = o3d.utility.Vector3dVector(np.asarray(pcd.colors)[idx[1:], :])
    return colonies


def main(ply_filename, annotations_filename, environment_distance, path):
    """ Exporting colony point clouds
    """
    short_name = "_".join(ply_filename.split('_')[0:4])

    # 1. PREPARATION SUBSET COLONY POINTS AND SCALE & ROTATE ALL POINTS ####
    print('Reading PLY file {} ...'.format(ply_filename))
    pcd = o3d.io.read_point_cloud('{}/{}'.format(path, ply_filename))
    print('Read assignment file ...')
    annotations = get_annotations('{}/{}'.format(path, annotations_filename))
    print('Building KDTree ...')
    pcd_tree = o3d.geometry.KDTreeFlann(pcd)

    # 2. FIND ENV ####
    print('Getting neighbourhood range')
    environment = get_neighbourhood(annotations, pcd, pcd_tree, environment_distance)

    # 3. EXPORT POINTS
    for name in environment:
        colony_env = environment[name]
        print(colony_env)
        print(type(colony_env))
        o3d.io.write_point_cloud('/Volumes/KP3/coralscape/environment_points/{}_env.ply'.format(name), colony_env)  # should save normals & rgb


if __name__ == '__main__':
    # Arguments
    parser = argparse.ArgumentParser(prog="Colony clean and measure")
    parser.add_argument('ply_filename')
    parser.add_argument('annotations_filename')
    parser.add_argument('environment_distance', type=float)
    args = parser.parse_args()

    ply_filename = args.ply_filename
    annotations_filename = args.annotations_filename
    environment_distance = args.environment_distance
    path = "/Volumes/KP3/coralscape"

    # e.g.,
    # args.ply_filename = "cur_cas_05m_20201212_rotated_scaled.ply"
    # args.annotations_filename = "scaled_annotations_cur_cas_05m_20201212.csv"
    # args.environment_distance = 0.1

    # PLOTS & ROTATIONS:
    # CA05 cur_cas_05m_20201212_rotated_scaled.ply
    # scaled_annotations_cur_cas_05m_20201212.csv
    # WP05 cur_kal_05m_20200214_rotated_scaled.ply
    # scaled_annotations_cur_kal_05m_20200214.csv
    # WP10 cur_kal_10m_20200214_rotated_scaled.ply
    # scaled_annotations_cur_kal_10m_20200214.csv
    # WP20 cur_kal_20m_20200214_rotated_scaled.ply
    # scaled_annotations_cur_kal_20m_20200214.csv
    # SB05 cur_sna_05m_20200303_rotated_scaled.ply
    # scaled_annotations_cur_sna_05m_20200303.csv
    # SB10 cur_sna_10m_20200303_rotated_scaled.ply
    # scaled_annotations_cur_sna_10m_20200303.csv
    # SB20 cur_sna_20m_20200303_rotated_scaled.ply
    # scaled_annotations_cur_sna_20m_20200303.csv
    main(ply_filename, annotations_filename, environment_distance, path)
