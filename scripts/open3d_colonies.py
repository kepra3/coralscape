#!/usr/anaconda3/envs/open3d/bin/python
# -*- coding: utf-8 -*-

"""
@author: kprata
@date created: 30/12/21
@description: NOTE: there is a PATH variable
TODO: dealing with overlapping colonies
"""

# Modules
import argparse
import numpy as np
import open3d as o3d
import copy
import json
from sklearn.linear_model import LinearRegression

__author__ = 'Pim Bongaerts and Katharine Prata'
__copyright__ = 'Copyright (C) 2021 Pim Bongaerts and Katharine Prata'
__license__ = 'GPL'

IGNORE_ANNOTATIONS = ['left', 'right', 'X']
V_DISTANCE = -10
PATH = "/Users/kprata/Dropbox/agaricia_project_2019/shalo_ag/Photogrammetry/CloudCompare/WP20"


class Viscore_metadata(object):
    """  """

    def __init__(self, subsets_filename, short_name):
        subsets = json.load(open('{}/{}'.format(PATH, subsets_filename)))
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


def get_annotations(annotations_path):
    """ Read annotations from txt file """
    annotations = {}
    annotations_file = open(annotations_path, 'r')
    for line in annotations_file:
        if any(flag in line for flag in IGNORE_ANNOTATIONS):
            continue
        cols = line.rstrip().replace(',', ' ').split()
        annotations[cols[0]] = [float(i) for i in cols[1:4]]
    annotations_file.close()
    return annotations


def calculate_euler_angle(opposite, adjacent):
    # Formula for:
    # atan2(opposite, adjacent)
    if adjacent > 0:
        theta = np.arctan(opposite / adjacent) * 180 / np.pi
    elif adjacent < 0 & opposite >= 0:
        theta = (np.arctan(opposite / adjacent) + np.pi) * 180 / np.pi
    elif adjacent < 0 & opposite < 0:
        theta = (np.arctan(opposite / adjacent) - np.pi) * 180 / np.pi
    elif adjacent == 0 & opposite > 0:
        theta = (np.arctan(opposite / adjacent) + np.pi / 2) * 180 / np.pi
    elif adjacent == 0 & opposite < 0:
        theta = (np.arctan(opposite / adjacent) - np.pi / 2) * 180 / np.pi
    else:
        print("theta is undefined")

    return theta.__float__()


def rotate_matrix(pcd, up_vector, short_name):
    # TODO: use euler angles to rotate pcd and annotations
    # up_vector = viscore_md.dd[0:3]
    origin = [0, 0, 0]

    # angle xz, phi
    x_diff = origin[0] - up_vector[0]  # opposite - reef perpendicular
    z_diff = origin[2] - up_vector[2]  # adjacent - depth
    theta_xz = calculate_euler_angle(x_diff, z_diff)
    if short_name == "cur_sna_20m_20201202":
        theta_xz = -theta_xz

    # angle yz, psi
    y_diff = origin[1] - up_vector[1]  # opposite - reef parallel, y-coord stays the same
    z_diff = origin[2] - up_vector[2]  # adjacent - depth, z-coord is changed
    psi_yz = calculate_euler_angle(y_diff, z_diff)

    R = pcd.get_rotation_matrix_from_xyz((0, theta_xz, psi_yz))
    pcd_r = copy.deepcopy(pcd)
    pcd_r.rotate(R, center=(0, 0, 0))

    return pcd_r


def get_ranges(annotations_path, annotations):
    complete_set = {}
    annotations_file = open(annotations_path, 'r')
    for line in annotations_file:
        cols = line.rstrip().replace(',', ' ').split()
        complete_set[cols[0]] = [float(i) for i in cols[1:4]]
    annotations_file.close()
    ranges = {}
    for name in annotations:
        # euclidean distance
        ranges[name] = ((complete_set['{}_left'.format(name)][0] - complete_set['{}_right'.format(name)][0]) ** 2
                        + (complete_set['{}_left'.format(name)][1] - complete_set['{}_right'.format(name)][1]) ** 2
                        + (complete_set['{}_left'.format(name)][2] - complete_set['{}_right'.format(name)][2]) ** 2) \
                       ** 0.5
    return ranges


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


def calc_attachment_angles(colonies, axes_order):
    # TODO: use the PCA method of calculating angles
    """
    :param colonies: colony point clouds
    :param axes_order: the order of axes for angle calculation, e.g., [0, 1, 2]
    would be for xy, [0, 2, 1] would be for xz and [1, 2, 0] would be for yz
    :return all_theta: all theta values
    """
    all_theta = {}
    for name in colonies:
        sample = np.asarray(colonies[name].points)
        if len(sample) <= 3:
            continue

        # Fit a line
        axes_one = np.array(sample[:, axes_order[0]]).reshape(-1, 1)
        axes_two = np.array(sample[:, axes_order[1]])
        axes_three = np.mean(sample[:, axes_order[2]])
        model = LinearRegression().fit(axes_one, axes_two)
        # r_sq = model.score(axes_one, axes_two)
        axes_two_prediction = (model.intercept_ + model.coef_ * axes_one).reshape(-1, 1)

        # Calculate theta
        axes_one_diff = abs(axes_one[1] - axes_one[0])
        axes_two_diff = abs(axes_two_prediction[1] - axes_two_prediction[0])
        theta = np.arctan(axes_two_diff / axes_one_diff) * 180 / np.pi
        all_theta[name] = theta

    return all_theta


def calc_relative_height_and_overhang(colonies, annotations, ranges, pcd, pcd_tree, overhang_value):
    heights = {}
    colonies_environment = {}
    all_relative_heights = {}
    all_overhang = {}

    for name in colonies:
        sample = np.asarray(colonies[name].points)
        if len(sample) <= 3:
            continue
        colony_mean_height = np.mean(sample[:, 2])
        heights[name] = colony_mean_height

        # Use a larger radius, i.e., environment
        [k, idx, _] = pcd_tree.search_radius_vector_3d(annotations[name], ranges[name] * 2)
        # Store colony environment as separate point cloud
        colonies_environment[name] = o3d.geometry.PointCloud()
        colonies_environment[name].points = o3d.utility.Vector3dVector(np.asarray(pcd.points)[idx[1:], :])
        sample_environment = np.asarray(colonies_environment[name].points)
        env_min_height = min(sample_environment[:, 2])
        env_max_height = max(sample_environment[:, 2])
        env_height_range = env_max_height - env_min_height

        # Proportion
        all_relative_heights[name] = (colony_mean_height - env_min_height) / env_height_range

        x_high = (annotations[name][0] + overhang_value)
        x_low = (annotations[name][0] - overhang_value)
        y_high = (annotations[name][1] + overhang_value)
        y_low = (annotations[name][1] - overhang_value)

        overhang = sample_environment[(sample_environment[:, 0] < x_high.__float__()) &
                                      (sample_environment[:, 0] > x_low.__float__()) &
                                      (sample_environment[:, 1] < y_high.__float__()) &
                                      (sample_environment[:, 1] > y_low.__float__())]
        overhang_list = []
        for i in range(0, len(overhang[:, 2])):
            if overhang[:, 2][i] > (annotations[name][2] + (overhang_value * 2)):
                overhang_list.append('Yes')
            else:
                overhang_list.append('No')

        if overhang_list.count('Yes') > 1:
            overhang_presence = 'Yes'
        else:
            overhang_presence = 'No'

        all_overhang[name] = overhang_presence
    return all_relative_heights, all_overhang


def calc_overhang(colonies, range, pcd):
    # TODO: subset a cylinder and calculate proportion of z values shaded,
    #  use an expanding cylinder from centre?
    all_overhang_percentage = {}
    return all_overhang_percentage


def calc_rugosity():
    # TODO: calculate rugosity!
    all_rugosity_colony = {}
    all_rugosity_environment = {}
    return all_rugosity_colony, all_rugosity_environment


def main(ply_filename, annotations_filename, subsets_filename):
    # Make a short name from the file name
    short_name = "_".join(ply_filename.split('_')[0:4])
    print('Reading PLY file ...')
    pcd = o3d.io.read_point_cloud('{}/{}'.format(PATH, ply_filename))
    print('Build KDTree from point cloud ...')
    pcd_tree = o3d.geometry.KDTreeFlann(pcd)
    print('Read assignment file ...')
    annotations = get_annotations('{}/{}'.format(PATH, annotations_filename))
    print('Get ranges for each sample ...')
    ranges = get_ranges('{}/{}'.format(PATH, annotations_filename), annotations)
    print('Read viscore metadata file ...')
    viscore_md = Viscore_metadata(subsets_filename, short_name)
    print('Searching radius around annotations ...')
    colonies = {}
    translated_colonies = {}
    connect_points = []
    connect_colors = []
    for name in annotations:
        [k, idx, _] = pcd_tree.search_radius_vector_3d(annotations[name], ranges[name])
        # Store colony as separate point cloud
        colonies[name] = o3d.geometry.PointCloud()
        colonies[name].points = o3d.utility.Vector3dVector(np.asarray(pcd.points)[idx[1:], :])
        random_color = list(np.random.choice(np.linspace(0, 1, 255), size=3))
        colonies[name].paint_uniform_color(random_color)
        connect_colors.append(random_color)
        # Vertical offset along up vector
        v_translation = np.array(viscore_md.dd[0:3]) * V_DISTANCE
        translated_colonies[name] = copy.deepcopy(colonies[name])
        translated_colonies[name].translate(v_translation)
        # Store points for the connecting lines
        connect_points.append(annotations[name])
        connect_points.append(annotations[name] + v_translation)
    # Create connected lineset
    connecting_lineset = generate_connecting_lineset(connect_points, connect_colors)
    # Join all geometries
    all_geoms = list(translated_colonies.values())
    all_geoms.append(pcd)
    all_geoms.append(connecting_lineset)
    # Visualise
    # o3d.visualization.draw_geometries(all_geoms,
    #                                  zoom=0.4,
    #                                  front=viscore_md.cam_eye,
    #                                  lookat=viscore_md.cam_target,
    #                                  up=viscore_md.cam_up)

    # Calculate colony angles
    theta_xy = calc_attachment_angles(colonies, axes_order=[0, 1, 2])
    print("Values for theta x-y: {}".format(theta_xy))
    theta_xz = calc_attachment_angles(colonies, axes_order=[0, 2, 1])
    print("Values for theta x-z: {}".format(theta_xz))
    theta_yz = calc_attachment_angles(colonies, axes_order=[1, 2, 0])
    print("Values for theta y-z: {}".format(theta_yz))

    # Calculate relative height and overhang
    overhang_value = 0.1 / viscore_md.scale_factor
    relative_height, overhang = calc_relative_height_and_overhang(colonies, annotations, ranges,
                                                                  pcd, pcd_tree, overhang_value)
    print(relative_height)
    print(overhang)


# TODO: COMBINE ALL THE DICTIONARIES FOR OUTPUT

if __name__ == '__main__':
    ply_filename = "cur_kal_20m_20200214_decvis_02_subsampler100000.ply"
    annotations_filename = "cur_kal_20m_20200214_decvis_02_KP_16-12-21_completed.txt"
    subsets_filename = "subsets.json"

    main(ply_filename, annotations_filename, subsets_filename)
