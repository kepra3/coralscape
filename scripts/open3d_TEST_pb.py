#!/usr/anaconda3/envs/open3d/bin/python
# -*- coding: utf-8 -*-

"""
@author: kprata
@date created: 16/12/21
@description: TODO
"""
# TODO: statistical outlier removal

import numpy as np
import open3d as o3d

IGNORE_ANNOTATIONS = ['left', 'right', 'X']
SEARCH_RADIUS = 0.2


def get_annotations(annotations_filename):
    """ Read annotations from txt file """
    annotations = {}
    annotations_file = open(annotations_filename, 'r')
    for line in annotations_file:
        if any(flag in line for flag in IGNORE_ANNOTATIONS):
            continue
        cols = line.rstrip().replace(',', ' ').split()
        annotations[cols[0]] = [float(i) for i in cols[1:4]]
    annotations_file.close()
    return annotations


def main(ply_filename, annotations_filename, SCALE):
    print('Reading PLY file ...')
    pcd = o3d.io.read_point_cloud(ply_filename)
    print('Build KDTree from point cloud ...')
    pcd_tree = o3d.geometry.KDTreeFlann(pcd)
    print('Read assignment file ...')
    annotations = get_annotations(annotations_filename)
    print('Searching radius around annotations ...')
    colonies = {}
    connect_points = []
    connect_colors = []
    for name in annotations:
        [k, idx, _] = pcd_tree.search_radius_vector_3d(annotations[name], SEARCH_RADIUS / SCALE)
        # np.asarray(pcd.colors)[idx[1:], :] = np.asarray(original_colors)[idx[1:], :]

        # Store colony as separate point cloud
        colonies[name] = o3d.geometry.PointCloud()
        colonies[name].points = o3d.utility.Vector3dVector(np.asarray(pcd.points)[idx[1:], :])
        random_color = list(np.random.choice(np.linspace(0, 1, 255), size=3))
        colonies[name].paint_uniform_color(random_color)
        connect_colors.append(random_color)

        # Store points for the connecting lines
        connect_points.append(annotations[name])
    # Join all geometries and visualize
    all_geoms = list(colonies.values())
    all_geoms.append(pcd)
    # might have to use another function here
    o3d.visualization.draw_geometries(all_geoms,
                                      zoom=0.4)


if __name__ == '__main__':
    # parser = argparse.ArgumentParser(description=__doc__)
    # parser.add_argument('ply_filename', metavar='ply_filename',
    #                    help='PLY filename')
    # parser.add_argument('annotations_filename', metavar='annotations_filename',
    #                    help='txt file with annotations (label, x, y, z)')
    # parser.add_argument('subsets_filename', metavar='subsets_filename',
    #                    help='subset.json file from viscore')
    # args = parser.parse_args()
    # main(args.ply_filename, args.annotations_filename, args.subsets_filename)
    PATH = "/Users/kprata/Dropbox/agaricia_project_2019/shalo_ag/Photogrammetry/CloudCompare/WP20"
    ply_filename = "{}/cur_kal_20m_20200214_decvis_02_subsample01.ply".format(PATH)
    annotations_filename = "{}/cur_kal_20m_20200214_decvis_02_KP_16-12-21_completed.txt".format(PATH)

    # Set scaling factor
    if ply_filename == "{}/cur_kal_20m_20200214_decvis_02_subsample01.ply".format(PATH):
        SCALE: np.ndarray = 0.18456881059341612
    else:
        SCALE = np.around(1, 17)

    main(ply_filename, annotations_filename, SCALE)
