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
import pandas as pd
from sklearn.linear_model import LinearRegression

__author__ = 'Pim Bongaerts and Katharine Prata'
__copyright__ = 'Copyright (C) 2021 Pim Bongaerts and Katharine Prata'
__license__ = 'GPL'

IGNORE_ANNOTATIONS = ['left', 'right', 'X']
V_DISTANCE = -10
PATH = "/Users/kprata/Dropbox/agaricia_project_2019/shalo_ag/Photogrammetry/CloudCompare/WP20"


# check sample 897 in plot

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


def rotate_matrix(pcd, up_vector):
    # TODO: use euler angles to rotate pcd and annotations
    # up_vector = viscore_md.dd[0:3]
    origin = [0, 0, 0]
    # angle xz, phi
    x_diff = origin[0] - up_vector[0]  # opposite - reef perpendicular
    z_diff = origin[2] - up_vector[2]  # adjacent - depth
    theta_xz = calculate_euler_angle(x_diff, z_diff)
    print('Theta is ...', theta_xz)
    # angle yz, psi
    y_diff = origin[1] - up_vector[1]  # opposite - reef parallel, y-coord stays the same
    z_diff = origin[2] - up_vector[2]  # adjacent - depth, z-coord is changed
    psi_yz = calculate_euler_angle(y_diff, z_diff)
    print('Psi is ...', psi_yz)
    # needs radians input
    theta_xz_radians = theta_xz / 180 * np.pi
    psi_yz_radians = psi_yz / 180 * np.pi
    R = pcd.get_rotation_matrix_from_xyz((psi_yz_radians, theta_xz_radians, 0))
    return R


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
        ranges[name] = (((complete_set['{}_left'.format(name)][0] - complete_set['{}_right'.format(name)][0]) ** 2
                         + (complete_set['{}_left'.format(name)][1] - complete_set['{}_right'.format(name)][1]) ** 2
                         + (complete_set['{}_left'.format(name)][2] - complete_set['{}_right'.format(name)][2]) ** 2) \
                        ** 0.5) / 2
    return ranges


def get_neighbourhood(annotations, pcd, ranges, viscore_md=None, lines=False):

    print('Building KDTree ...')
    pcd_tree = o3d.geometry.KDTreeFlann(pcd)

    colonies = {}

    if lines:
        translated_colonies = {}
        connect_points = []
        connect_colors = []
    else:
        print("Connecting lines not used")

    print('Searching ...')
    for name in annotations:
        [k, idx, _] = pcd_tree.search_radius_vector_3d(annotations[name], ranges[name])
        # Store colony as separate point cloud
        colonies[name] = o3d.geometry.PointCloud()
        colonies[name].points = o3d.utility.Vector3dVector(np.asarray(pcd.points)[idx[1:], :])
        if lines:
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
        else:
            continue
    if lines:
        return colonies, translated_colonies, connect_points
    else:
        return colonies


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


def calc_attachment_angles(colonies, axes_order, interval_num):
    # TODO: select a point every three
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

        # Select interval
        interval_axes_one = (sample[:, 0].max() - sample[:, 0].min()) / interval_num
        intervals = [0]
        subset = {}
        representative_points = []
        for i in range(1, interval_num):
            intervals.append(sample[:, 0].min() + interval_axes_one * i)
            q = np.where(sample[:, 0] < intervals[i])
            subset[i] = q
            representative_points.append(np.random.choice(subset[i][0], size=3))
            # TODO: test different sizes and draw
        # Get rid of hierarchical structure in List
        rep_points_list = [element for sublist in representative_points for element in sublist]

        # Choose points
        axes_one = np.array(sample[rep_points_list, axes_order[0]]).reshape(-1, 1)
        axes_two = np.array(sample[rep_points_list, axes_order[1]])

        # Fit a line
        model = LinearRegression().fit(axes_one, axes_two)
        # r_sq = model.score(axes_one, axes_two)
        # calculate residuals of the model
        axes_two_prediction = (model.intercept_ + model.coef_ * axes_one).reshape(-1, 1)

        # Calculate theta
        axes_one_diff = abs(axes_one[1] - axes_one[0])
        axes_two_diff = abs(axes_two_prediction[1] - axes_two_prediction[0])
        theta = np.arctan(axes_two_diff / axes_one_diff) * 180 / np.pi
        all_theta[name] = theta.__float__()

    return all_theta


def calc_relative_height_and_overhang(rotated_colonies, rotated_annotations, ranges,
                                      pcd_r, pcd_tree_r, overhang_value):
    heights = {}
    colonies_environment_r = {}
    all_relative_heights = {}
    all_overhang = {}

    for name in rotated_colonies:
        sample = np.asarray(rotated_colonies[name].points)
        if len(sample) <= 3:
            continue
        colony_mean_height = np.mean(sample[:, 2])
        heights[name] = colony_mean_height

        # Extract points
        sample_environment_r = np.asarray(colonies_environment_r[name].points)
        env_min_height = min(sample_environment_r[:, 2])
        env_max_height = max(sample_environment_r[:, 2])
        env_height_range = env_max_height - env_min_height

        # Proportion
        all_relative_heights[name] = (colony_mean_height - env_min_height) / env_height_range

        overhang_range = []
        for i in range(1, 20):
            # Trying to make stop if exceed range for that sample.
            print('Sample is:', name)
            x_high = (rotated_annotations[name][0] + overhang_value * i)
            print("x_high:", x_high)
            x_low = (rotated_annotations[name][0] - overhang_value * i)
            print("x_low:", x_low)
            y_high = (rotated_annotations[name][1] + overhang_value * i)
            print("y_high:", y_high)
            y_low = (rotated_annotations[name][0] - overhang_value * i)
            print("y_low:", y_low)
            print("x_range:", (x_high - x_low))
            print("y_range:", (y_high - y_low))
            print("range:", ranges[name])
            print("overhang:", overhang_value * 1)
            if ranges[name] < (x_high - x_low) and ranges[name] < (y_high - y_low):
                break
            else:
                # Overhang presence TODO: see notes figure out expanding cylinder problem
                overhang = sample_environment_r[(sample_environment_r[:, 0] < x_high.__float__()) &
                                                (sample_environment_r[:, 0] > x_high.__float__() - overhang_value * i) &

                                                (sample_environment_r[:, 0] > x_low.__float__()) &
                                                (sample_environment_r[:, 1] < y_high.__float__()) &
                                                (sample_environment_r[:, 1] > y_low.__float__())]
                overhang_list = []
                for j in range(0, len(overhang[:, 2])):
                    if overhang[:, 2][j] > (rotated_annotations[name][2] + (overhang_value * 2)):
                        overhang_list.append(1)  # Yes
                    else:
                        overhang_list.append(0)  # No

                if overhang_list.count(1) > 1:
                    overhang_presence = 1
                else:
                    overhang_presence = 0

                overhang_range.append(overhang_presence)
                print('overhang range is', overhang_range)
        print('overhang range is', overhang_range)
        overhang_prop = sum(overhang_range) / len(overhang_range)
        all_overhang[name] = overhang_prop
    return all_relative_heights, all_overhang


# def calc_overhang(colonies, range, pcd):
# TODO: subset a cylinder and calculate proportion of z values shaded,
#  use an expanding cylinder from centre?
#    all_overhang_percentage = {}
#    return all_overhang_percentage


# def calc_rugosity():
#    # TODO: calculate rugosity! Look at xz and yz profiles and then calculate the distances between each point over a specific x or y distance.
#    all_rugosity_colony = {}
#    all_rugosity_environment = {}
#    return all_rugosity_colony, all_rugosity_environment


def main(ply_filename, annotations_filename, subsets_filename):
    # Make a short name from the file name
    short_name = "_".join(ply_filename.split('_')[0:4])

    print('Reading PLY file ...')
    pcd = o3d.io.read_point_cloud('{}/{}'.format(PATH, ply_filename))

    print('Read viscore metadata file ...')
    viscore_md = Viscore_metadata(subsets_filename, short_name)

    print('Rotating matrix ...')
    up_vector = viscore_md.dd[0:3]
    R = rotate_matrix(pcd, up_vector, short_name)
    pcd_r = copy.deepcopy(pcd)
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
    rotated_colonies = get_neighbourhood(rotated_annotations, pcd_r, ranges)

    # Calculate rotated colony angles
    # TODO: test optimisation of slope angle
    theta_xy = calc_attachment_angles(rotated_colonies, axes_order=[0, 1, 2], interval_num=10)
    theta_xz = calc_attachment_angles(rotated_colonies, axes_order=[0, 2, 1], interval_num=10)
    theta_yz = calc_attachment_angles(rotated_colonies, axes_order=[1, 2, 0], interval_num=10)

    # Define neighbourhood range
    double_range = {}
    for name in ranges:
        double_range[name] = ranges[name] * 2
    print('Getting neighbourhood range')
    colonies_environment_r = get_neighbourhood(rotated_annotations, pcd_r, double_range)

    # Calculate relative height and overhang
    overhang_value = 0.01 / viscore_md.scale_factor
    relative_height, overhang = calc_relative_height_and_overhang(rotated_colonies, rotated_annotations, ranges,
                                                                  pcd_r, pcd_tree_r, overhang_value)
    # Scale ranges
    scaled_ranges = copy.deepcopy(ranges)
    for name in ranges:
        scaled_ranges[name] = ranges[name] * viscore_md.scale_factor

    # Join dictionaries
    dicts = [theta_xy, theta_xz, theta_yz, relative_height, overhang, scaled_ranges]
    sample_metadata = {}
    for key in dicts[0]:
        sample_metadata[key] = [d[key] for d in dicts]

    # Scale annotations
    scaled_annotations = copy.deepcopy(rotated_annotations)
    for name in rotated_annotations:
        for i in range(3):
            scaled_annotations[name][i] = rotated_annotations[name][i] * viscore_md.scale_factor

    scaled_annotations_df = pd.DataFrame(scaled_annotations).T
    scaled_annotations_df.columns = ['x', 'y', 'z']
    scaled_annotations_df.to_csv('~/git/coralscape/results/scaled_annotations_{}.csv'.format(ply_filename))
    print('Saved scaled annotations to file ...')

    # Convert to .csv for input to R
    df = pd.DataFrame(sample_metadata).T
    df.columns = ['xy', 'xz', 'yz', 'prop', 'overhang', 'range']
    df['range'] = df['range'] * viscore_md.scale_factor  # scaled range
    df.to_csv('~/git/coralscape/results/sample_metadata_{}.csv'.format(ply_filename))
    print('Saved metadata to file ...')

    coordinates = pd.DataFrame(rotated_annotations).T
    coordinates.columns = ['x', 'y', 'z']
    coordinates.to_csv('~/git/coralscape/results/rotated_coordinates_{}.csv'.format(ply_filename))
    print('Saved coordinates to file ...')


if __name__ == '__main__':
    ply_filename = "cur_kal_20m_20200214_decvis_02.ply"
    annotations_filename = "cur_kal_20m_20200214_decvis_02_KP_16-12-21_completed.txt"
    subsets_filename = "subsets.json"

    # WP05 cur_kal_05m_20200214_decvis_02_KP
    # theta -6.91, psi 21.08
    # WP10 cur_kal_10m_20200214_decvis_02_KP_905
    # theta = -25.11, psi = 11.65
    # WP20 cur_kal_20m_20200214_decvis_02_KP_16-12-21_completed
    # theta = 9.02, psi = 19.25
    # SB05 cur_sna_05m_20200303_decvis_02_SF_HI_19-1-22
    # theta = -3.90, psi = 12.30
    # SB10 cur_sna_10m_20201202_decvisann_HI_14-12-21
    # theta = -0.72, psi = -2.71
    # SB20 cur_sna_20m_20190410_decvisann_HI_12_12
    # theta = -16.82, psi = 18.23
    main(ply_filename, annotations_filename, subsets_filename)
