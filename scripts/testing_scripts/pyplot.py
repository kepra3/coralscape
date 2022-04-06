#!/usr/anaconda3/bin/python
# -*- coding: utf-8 -*-
# remember to edit python environment if required.

import sys
import argparse
import numpy as np
import open3d as o3d
import open3d.visualization.gui as gui
import copy
import json

__author__ = 'Pim Bongaerts'
__copyright__ = 'Copyright (C) 2021 Pim Bongaerts'
__license__ = 'GPL'
IGNORE_ANNOTATIONS = ['left', 'right', 'X']
SEARCH_RADIUS = 0.2
V_DISTANCE = -10


class Viscore_metadata(object):
    """  """

    def __init__(self, subsets_filename, short_name):
        subsets = json.load(open(subsets_filename))
        # Determine json key of primary model
        if '{0}/{0}'.format(short_name) in subsets['d']:
            subsets_ortho = subsets['d']['{0}/{0}'.format(short_name)]['c']['ortho']
        elif self.short_name in subsets['d']:
            subsets_ortho = subsets['d'][short_name]['c']['ortho']
        else:
            print('Model not found in subsets.json!')
        self.dd = subsets_ortho['dd']
        self.scale_factor = subsets_ortho['scale_factor']
        self.r = subsets_ortho['vecs']['r']
        self.u = subsets_ortho['vecs']['u']
        self.n = subsets_ortho['vecs']['n']
        self.c = subsets_ortho['vecs']['c']
        self.cc = subsets_ortho['vecs']['cc']
        self.cam_up = subsets_ortho['vecs']['cam']['up']
        self.cam_eye = subsets_ortho['vecs']['cam']['eye']
        self.cam_target = subsets_ortho['vecs']['cam']['target']


def get_annotations(annotations_filename):
    """ Read anotations from txt file """
    annotations = {}
    annotations_file = open(annotations_filename, 'r')
    for line in annotations_file:
        if any(flag in line for flag in IGNORE_ANNOTATIONS):
            next
        cols = line.rstrip().replace(',', ' ').split()
        annotations[cols[0]] = [float(i) for i in cols[1:4]]
    annotations_file.close()
    return annotations


def generate_connecting_lineset(connect_points, connect_colors):
    # Create lines connecting original colony and offset pointcloud
    connect_lines = []
    for i in range(0, len(connect_points), 2):
        connect_lines.append([i, i + 1])
    connecting_lineset = o3d.geometry.LineSet()
    connecting_lineset.points = o3d.utility.Vector3dVector(connect_points)
    connecting_lineset.lines = o3d.utility.Vector2iVector(connect_lines)
    connecting_lineset.colors = o3d.utility.Vector3dVector(connect_colors)
    return connecting_lineset


def main(ply_filename, annotations_filename, subsets_filename):
    short_name = '_'.join(ply_filename.split('_')[0:4])
    print('Reading PLY file ...')
    pcd = o3d.io.read_point_cloud(ply_filename)
    print('Build KDTree from point cloud ...')
    pcd_tree = o3d.geometry.KDTreeFlann(pcd)
    print('Read assignment file ...')
    annotations = get_annotations(annotations_filename)
    print('Read viscore metadata file ...')
    viscore_md = Viscore_metadata(subsets_filename, short_name)
    print('Searching radius around annotations ...')
    colonies = {}
    connect_points = []
    connect_colors = []
    for name in annotations:
        [k, idx, _] = pcd_tree.search_radius_vector_3d(annotations[name], SEARCH_RADIUS / viscore_md.scale_factor)
        # np.asarray(pcd.colors)[idx[1:], :] = np.asarray(original_colors)[idx[1:], :]

        # Store colony as separate point cloud
        colonies[name] = o3d.geometry.PointCloud()
        colonies[name].points = o3d.utility.Vector3dVector(np.asarray(pcd.points)[idx[1:], :])
        random_color = list(np.random.choice(np.linspace(0, 1, 255), size=3))
        colonies[name].paint_uniform_color(random_color)
        connect_colors.append(random_color)
        # Vertical offset along up vector
        v_translation = np.array(viscore_md.dd[0:3]) * V_DISTANCE
        colonies[name].translate(v_translation)
        # Store points for the connecting lines
        connect_points.append(annotations[name])
        connect_points.append(annotations[name] + v_translation)
    # Create connected lineset
    connecting_lineset = generate_connecting_lineset(connect_points, connect_colors)
    # Join all geometries and visualize
    all_geoms = list(colonies.values())
    all_geoms.append(pcd)
    all_geoms.append(connecting_lineset)
    o3d.visualization.draw_geometries(all_geoms,
                                      zoom=0.4,
                                      front=viscore_md.cam_eye,
                                      lookat=viscore_md.cam_target,
                                      up=viscore_md.cam_up)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('ply_filename', metavar='ply_filename',
                        help='PLY filename')
    parser.add_argument('annotations_filename', metavar='annotations_filename',
                        help='txt file with annotations (label, x, y, z)')
    parser.add_argument('subsets_filename', metavar='subsets_filename',
                        help='subset.json file from viscore')
    args = parser.parse_args()
    main(args.ply_filename, args.annotations_filename, args.subsets_filename)