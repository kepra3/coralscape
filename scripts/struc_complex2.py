#!/usr/anaconda3/envs/open3d/bin/python
# -*- coding: utf-8 -*-

# Modules
import argparse
import numpy as np
import open3d as o3d
import copy
import json
import pandas as pd
from sklearn.linear_model import LinearRegression
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from typing import Any
import alphashape
from descartes import PolygonPatch

__author__ = 'Pim Bongaerts and Katharine Prata'
__copyright__ = 'Copyright (C) 2021 Pim Bongaerts and Katharine Prata'
__license__ = 'GPL'

IGNORE_ANNOTATIONS = ['left', 'right', 'X']
V_DISTANCE = -10
PATH = "/Users/kprata/Dropbox/agaricia_project_2019/shalo_ag/Photogrammetry/CloudCompare/SB20"


class Viscore_metadata(object):
    """ Defining viscore metadata as a class object"""

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
    """ Formula for calculating euler angles """
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
        theta = None
        print("theta is undefined")
    return theta.__float__()


def rotate_matrix(pcd, up_vector):
    """ Create rotation matrix from euler angles and up vector """
    origin = [0, 0, 0]
    # angle xz, theta
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


def get_ranges(annotations_path, annotations, scale):
    """ Get longest length of colony divided by 2"""
    complete_set = {}
    annotations_file = open(annotations_path, 'r')
    for line in annotations_file:
        cols = line.rstrip().replace(',', ' ').split()
        complete_set[cols[0]] = [float(i) for i in cols[1:4]]
    annotations_file.close()
    ranges = {}
    for name in annotations:
        # euclidean distance
        x1 = complete_set['{}_left'.format(name)][0]
        x2 = complete_set['{}_right'.format(name)][0]
        y1 = complete_set['{}_left'.format(name)][1]
        y2 = complete_set['{}_right'.format(name)][1]
        z1 = complete_set['{}_left'.format(name)][2]
        z2 = complete_set['{}_right'.format(name)][2]
        ranges[name] = (((x2 - x1) ** 2 + (y2 - y1) ** 2 + (z2 - z1) ** 2) ** 0.5 / 2) * scale
    return ranges


def get_neighbourhood(annotations, pcd, pcd_tree, ranges, viscore_md=None, lines=False):
    """ Find neighbouring points """
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
        # Store colony as separate point cloud with points, normals and colours!!!
        colonies[name] = o3d.geometry.PointCloud()
        colonies[name].points = o3d.utility.Vector3dVector(np.asarray(pcd.points)[idx[1:], :])
        colonies[name].normals = o3d.utility.Vector3dVector(np.asarray(pcd.normals)[idx[1:], :])
        colonies[name].colors = o3d.utility.Vector3dVector(np.asarray(pcd.colors)[idx[1:], :])
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
    """ Create lines connecting original colony and offset point cloud """
    connect_lines = []
    for i in range(0, len(connect_points), 2):
        connect_lines.append([i, i + 1])
    connecting_lineset = o3d.geometry.LineSet()
    connecting_lineset.points = o3d.utility.Vector3dVector(connect_points)
    connecting_lineset.lines = o3d.utility.Vector2iVector(connect_lines)
    connecting_lineset.colors = o3d.utility.Vector3dVector(connect_colors)
    return connecting_lineset


def create_mesh_ball_pivot(pcd):
    # estimate radius for rolling ball
    distances = pcd.compute_nearest_neighbor_distance()
    avg_dist = np.mean(distances)
    radius = 1.5 * avg_dist
    radii = []
    for i in range(1, 10):
        radii.append(radius * i)
    mesh = o3d.geometry.TriangleMesh.create_from_point_cloud_ball_pivoting(pcd, o3d.utility.DoubleVector(radii))
    return mesh


def create_mesh_poisson(pcd):
    pcd.estimate_normals()
    poisson_mesh = o3d.geometry.TriangleMesh.create_from_point_cloud_poisson(pcd,
                                                                             depth=8,
                                                                             width=0,
                                                                             scale=1.1,
                                                                             linear_fit=False)[0]
    # bbox = pcd.get_axis_aligned_bounding_box()
    # p_mesh_crop = poisson_mesh.crop(bbox)
    return poisson_mesh


def remove_outlier_points(pcd):
    print("Statistical outlier removal")
    cl, ind = pcd.remove_statistical_outlier(nb_neighbors=20, std_ratio=2.0)
    return ind


def display_inlier_outlier(cloud, ind):
    inlier_cloud = cloud.select_by_index(ind)
    outlier_cloud = cloud.select_by_index(ind, invert=True)
    print("Showing outliers (red) and inliers (gray): ")
    outlier_cloud.paint_uniform_color([1, 0, 0])
    inlier_cloud.paint_uniform_color([0.8, 0.8, 0.8])
    o3d.visualization.draw_geometries([inlier_cloud, outlier_cloud])


def get_cluster_triangles(mesh):
    print("Clustering connected triangles ...")
    with o3d.utility.VerbosityContextManager(o3d.utility.VerbosityLevel.Debug) as cm:
        triangle_clusters, cluster_n_triangles, cluster_area = (mesh.cluster_connected_triangles())
    triangle_clusters = np.asarray(triangle_clusters)
    cluster_n_triangles = np.asarray(cluster_n_triangles)
    cluster_area = np.asarray(cluster_area)
    return triangle_clusters, cluster_n_triangles, cluster_area


def largest_cluster(mesh, cluster_n_triangles, triangle_clusters):
    large_mesh: o3d.cpu.pybind.geometry.TriangleMesh = copy.deepcopy(mesh)
    largest_cluster_idx = cluster_n_triangles.argmax()
    triangles_to_remove = triangle_clusters != largest_cluster_idx
    large_mesh.remove_triangles_by_mask(triangles_to_remove)
    return large_mesh


def fit_a_plane_ransac(pcd):
    plane_model, inliers = pcd.segment_plane(distance_threshold=0.001,
                                             ransac_n=3,
                                             num_iterations=1000)
    [a, b, c, d] = plane_model
    print(f"Plane equation: {a:.2f}x + {b:.2f}y + {c:.2f}z + {d:.2f} = 0")
    return plane_model, inliers


def fit_a_lm(pcd, axes_order):
    # Choose points
    axes_one = np.asarray(pcd.points)[:, axes_order[0]].reshape(-1, 1)
    axes_two = np.asarray(pcd.points)[:, axes_order[1]]
    # Fit a line
    model = LinearRegression().fit(axes_one, axes_two)
    r_sq = model.score(axes_one, axes_two)
    print('R squared is ...', r_sq)
    # Calculate residuals of the model
    axes_two_prediction = (model.intercept_ + model.coef_ * axes_one).reshape(-1, 1)
    # Calculate theta
    axes_one_diff = abs(axes_one[1] - axes_one[0])
    axes_two_diff = abs(axes_two_prediction[1] - axes_two_prediction[0])
    theta = np.arctan(axes_two_diff / axes_one_diff) * 180 / np.pi
    print('Lm angle is ...', theta)
    return theta.__float__()


def calc_plane_angles(plane_model):
    # using the formula
    # plane_normal dot axis of interest normal / (magnitude of plane normal * magnitude of axis of interest normal)
    # which simplifies to
    #xz_normal = [0, 1, 0] have to figure out why these don't make sense
    #yz_normal = [1, 0, 0]
    plane_normal = [plane_model[0], plane_model[1], plane_model[2]]
    #mag_xz = np.linalg.norm(xz_normal)
    #mag_yz = np.linalg.norm(yz_normal)
    #mag_plane = np.linalg.norm(plane_normal)
    #cos_theta = np.dot(xz_normal, plane_normal)/(mag_xz*mag_plane)
    #cos_psi = np.dot(yz_normal, plane_normal)/(mag_yz*mag_plane)
    #theta = np.arccos(cos_theta) * 180./np.pi
    slope_xz = plane_normal[0] / plane_normal[2]
    slope_yz = plane_normal[1] / plane_normal[2]
    theta = np.arctan(slope_xz) * 180 / np.pi
    if theta > 90:
        theta = 180 - theta
    else:
        theta = theta
    print('The angle between x and z is ...', theta)  # about y-axis (xz)
    #psi = np.arccos(cos_psi) * 180./np.pi
    psi = np.arctan(slope_yz) * 180 / np.pi
    if psi > 90:
        psi = 180 - 90
    else:
        psi = psi
    print('The angle between y and z is ...', psi)  # about x-axis (yz)
    return theta.__float__(), psi.__float__()


def calc_rugosity(pcd_r, threeD_area):
    print('Projecting points to xy plane, z=0')
    x = np.asarray(pcd_r.points)[:, 0]
    y = np.asarray(pcd_r.points)[:, 1]
    z = np.repeat(0, len(np.asarray(pcd_r.points)[:, 1]))
    print('Creating polygon for 2D points')
    points_2d = np.asarray((x, y)).transpose()
    alpha_shape = alphashape.alphashape(points_2d, 2.0)
    rotated_twoD_area = alpha_shape.area
    print('2D area is ...', rotated_twoD_area)
    print('Rugosity is', threeD_area / rotated_twoD_area)
    rugosity = threeD_area / rotated_twoD_area
    if rugosity < 1:
        print("shape not complex rugosity likely essentially 1")
        rugosity = 1
    else:
        pass
    return rugosity, rotated_twoD_area


def calc_overhang(colony_pcd, pcd_env):
    dists = pcd_env.compute_point_cloud_distance(colony_pcd)
    dists = np.asarray(dists)
    ind = np.where(dists < 0.01)[0]
    np.asarray(pcd_env.colors)[ind] = [0, 0, 1]
    colony = np.asarray(pcd_env.points)[ind]
    env = np.asarray(pcd_env.points)
    q = []
    for i in range(0, len(colony)):
        p = np.where((env[:, 0] < colony[i, 0] + 0.01) & (env[:, 0] > colony[i, 0] - 0.01) &
                     (env[:, 1] < colony[i, 1] + 0.01) & (env[:, 1] > colony[i, 1] - 0.01) &
                     (env[:, 2] > max(colony[:, 2])))
        if p[0].size:
            np.asarray(pcd_env.colors)[p[0]] = [1, 0, 0]
            if len(p[0]) > 1:
                for j in range(0, len(p)):
                    p_int = int(p[0][j])
                    q.append(p_int)
            elif len(p[0]) == 1:
                p_int = int(p[0])
                q.append(p_int)
    unique_q = list()
    unique_points = 0
    for item in q:
        if item not in unique_q:
            unique_q.append(item)
            unique_points += 1
    if colony.size:
        overhang_prop = unique_points / len(colony)  # proportion works for now but area would be better
    else:
        print('issue with colony points ...')
        overhang_prop = None
    return overhang_prop


def calc_metrics(colonies, environment):
    colonies_theta_xz = {}
    colonies_psi_yz = {}
    colonies_rugosity = {}
    overhang_prop = {}
    outcrop_prop = {}
    colonies_height = {}
    environment_rugosity = {}

    for name in colonies:
        if len(np.asarray(colonies[name].points)) <= 3:
            print('Not including {} because < 3 points')
            continue
        else:
            print('\n\n\n ***** Starting new colony ...', name, '******')
            print('Getting colony angles ...')
            mesh = create_mesh_ball_pivot(colonies[name])
            triangle_clusters, cluster_n_triangles, cluster_area = get_cluster_triangles(mesh)
            large_mesh = largest_cluster(mesh, cluster_n_triangles, triangle_clusters)
            triangle_clusters, cluster_n_triangles, cluster_area = get_cluster_triangles(large_mesh)
            # TODO: NEED TO FILL HOLES IN MESH!
            threeD_area = cluster_area
            #vis = o3d.visualization.Visualizer()
            #vis.create_window(visible=True)  # works for me with False, on some systems needs to be true
            #vis.add_geometry(large_mesh)
            #vis.update_geometry(large_mesh)
            #vis.poll_events()
            #vis.update_renderer()
            #vis.capture_screen_image('{}.png'.format(name))
            #vis.destroy_window()
            print('Cluster area is ... {} m^2 for {}'.format(threeD_area, name))  # TODO: should output this
            print('Sampling points from mesh first uniformly then with poisson ...')
            colony_pcd = large_mesh.sample_points_uniformly(number_of_points=len(np.asarray(colonies[name].points)))
            colony_pcd = large_mesh.sample_points_poisson_disk(number_of_points=len(np.asarray(colonies[name].points)),
                                                               pcl=colony_pcd)
            # NEED TO FILL MESH HOLES
            print('Fitting a plane with RAMSAC to get colony angles')
            plane_model, inliers = fit_a_plane_ransac(colony_pcd)
            colony_theta_xz, colony_psi_yz = calc_plane_angles(plane_model)
            colonies_theta_xz[name] = abs(colony_theta_xz)
            colonies_psi_yz[name] = abs(colony_psi_yz)
            # RUGOSITY
            print('Getting colony rugosity ...')
            print('Rotate colony points according to colony slope')
            colony_psi_radians = colony_psi_yz / 180. * np.pi
            colony_theta_radians = colony_theta_xz / 180. * np.pi
            colony_R = colony_pcd.get_rotation_matrix_from_xyz((colony_psi_radians, colony_theta_radians, 0))
            colony_pcd_R = copy.deepcopy(colony_pcd)
            colony_pcd_R.rotate(colony_R, center=(0, 0, 0))
            rugosity = calc_rugosity(colony_pcd_R, threeD_area)
            colonies_rugosity[name] = rugosity
            # OVERHANG_PROP
            print('Getting overhang proportion ...')
            overhang_prop[name] = calc_overhang(colony_pcd, environment[name])
            print('Overhang proportion is ...', overhang_prop[name])
            # OUTCROP_PROP
            print('Getting outcrop proportion ...')
            sample = np.asarray(colonies[name].points)
            if len(sample) <= 3:
                continue
            colony_mean_height = np.mean(sample[:, 2])
            colonies_height[name] = colony_mean_height
            environment_array = np.asarray(environment[name].points)
            env_min_height = min(environment_array[:, 2])
            env_max_height = max(environment_array[:, 2])
            env_height_range = env_max_height - env_min_height
            outcrop_prop[name] = (colony_mean_height - env_min_height) / env_height_range
            print('Outcrop proportion is ...', outcrop_prop[name])
            # ENV RUGOSITY
            print('Getting environment rugosity ...')
            mesh = create_mesh_ball_pivot(environment[name])
            triangle_clusters, cluster_n_triangles, cluster_area = get_cluster_triangles(mesh)
            threeD_area = np.sum(cluster_area)
            print('Cluster area is ... {} m^2 for {}'.format(threeD_area, name))
            print('Sampling points from mesh first uniformaly then with poisson ...')
            environment_pcd = mesh.sample_points_uniformly(number_of_points=len(np.asarray(environment[name].points)))
            environment_pcd = mesh.sample_points_poisson_disk(number_of_points=len(np.asarray(environment[name].points)),
                                                              pcl=environment_pcd)
            print('Fit a plane to correct for environment slope')
            plane_model_env, inliers_env = fit_a_plane_ransac(environment_pcd)
            env_theta_xz, env_psi_yz = calc_plane_angles(plane_model_env)
            env_theta_radians = env_theta_xz / (180 * np.pi)
            env_psi_radians = env_psi_yz / (180 * np.pi)
            env_R = environment_pcd.get_rotation_matrix_from_xyz((env_psi_radians, env_theta_radians, 0))
            environment_pcd.rotate(env_R, center=(0, 0, 0))
            rugosity = calc_rugosity(environment_pcd, threeD_area)
            environment_rugosity[name] = rugosity
    return colonies_theta_xz, colonies_psi_yz, colonies_rugosity, overhang_prop, outcrop_prop, environment_rugosity


def main(ply_filename, annotations_filename, subsets_filename):
    # Make a result files as you go through samples
    out_one = "struc_complex_results.txt"
    with open(out_one, 'a') as results_out:
        if results_out.tell() == 0:
            print('Creating a new file\n')
            results_out.write(
                "\tplot_name\tsample_name\tcloud_points\tcolony_elevation\tcolony_rugosity\toverhang_prop\toutcrop_prop"
                "\tenvironment_rugosity\trange\n")
        else:
            print('File exists, appending\n')
    out_two = "colony_angle_rotation_results.txt"
    with open(out_two, 'a') as results_out:
        if results_out.tell() == 0:
            print('Creating a new file\n')
            results_out.write(
                "sample_name\tcloud_points\tplane_i,\tplane_j\tplane_k\taxis_i1\taxis_j1\taxis_k1\t"
                "axis_i2\taxis_j2\taxis_k2\televation\n")
        else:
            print('File exists, appending\n')
    out_three = "colony_area_results.txt"
    with open(out_three, 'a') as results_out:
        if results_out.tell() == 0:
            print('Creating a new file\n')
            results_out.write(
                "sample_name\tcloud_points\tthreeD_area\ttwoD_area\n")
        else:
            print('File exists, appending\n')
    # 1. PREPARATION SUBSET COLONY POINTS AND SCALE & ROTATE ALL POINTS ####
    short_name = "_".join(ply_filename.split('_')[0:4])
    print('Reading PLY file {} ...'.format(ply_filename))
    pcd = o3d.io.read_point_cloud('{}/{}'.format(PATH, ply_filename))
    print('Read viscore metadata file ...')
    viscore_md = Viscore_metadata(subsets_filename, short_name)
    print('Rotating matrix ...')
    up_vector = viscore_md.dd[0:3]
    R = rotate_matrix(pcd, up_vector)
    pcd_r = copy.deepcopy(pcd)
    pcd_r.rotate(R, center=(0, 0, 0))
    print('Scaling point cloud ...')
    pcd_r.scale(viscore_md.scale_factor, center=(0, 0, 0))
    print('Read assignment file ...')
    annotations = get_annotations('{}/{}'.format(PATH, annotations_filename))
    print('Rotate and scale annotations ...')
    rotated_annotations = {}
    for name in annotations:
        rotated_annotations[name] = np.matmul(R, annotations[name])
        rotated_scaled_annotations = rotated_annotations
        for i in range(3):
            rotated_scaled_annotations[name][i] = rotated_annotations[name][i] * viscore_md.scale_factor
    print('Get scaled ranges for each sample...')
    ranges = get_ranges('{}/{}'.format(PATH, annotations_filename), annotations, viscore_md.scale_factor)
    print('Building KDTree ...')
    pcd_tree_r = o3d.geometry.KDTreeFlann(pcd_r)

    #print('Get angles of construct, i.e., ground slope')
    #plane_model, inliers = fit_a_plane_ransac(pcd_r)
    #theta_xz, psi_yz = calc_plane_angles(plane_model)
    #o3d.visualization.draw_geometries([pcd_r])

    # 2. FIND COLONIES AND ENV ####
    print("\nSearching for colony around annotations ...")
    colonies = get_neighbourhood(rotated_scaled_annotations, pcd_r, pcd_tree_r, ranges)
    env_range = {}
    for name in ranges:
        env_range[name] = ranges[name] + 0.4
        #print('environment range is:', env_range[name])
    print('Getting neighbourhood range')
    environment = get_neighbourhood(rotated_scaled_annotations, pcd_r, pcd_tree_r, env_range)

    # 3. STATISTICS ####
    colony_theta_xz, colony_psi_yz, colony_rugosity, \
    overhang_prop, outcrop_prop, environment_rugosity = calc_metrics(colonies, environment)
    # TODO: rugosity sometimes high (5) check for specific colonies
    # TODO: should I instead calculate rugosity within a square around colony that's not scaled by the colony size?

    # 4. SAVE METADATA ####
    dicts = [colony_theta_xz, colony_psi_yz, colony_rugosity, overhang_prop, outcrop_prop, environment_rugosity, ranges]
    sample_metadata = {}
    for key in dicts[0]:
        sample_metadata[key] = [d[key] for d in dicts]
    # Convert to .csv for input to R
    df = pd.DataFrame(sample_metadata).T
    df.columns = ['colony_theta_xz', 'colony_psi_yz', 'colony_rugosity', 'overhang_prop', 'outcrop_prop', 'environment_rugosity', 'range']
    # Save to file
    df.to_csv('~/git/coralscape/results/sample_metadata_{}.csv'.format(ply_filename))
    print('Saved metadata to file ...')

    anno_df = pd.DataFrame(rotated_scaled_annotations).T
    anno_df.columns = ['x', 'y', 'z']
    anno_df.to_csv('~/git/coralscape/results/scaled_annotations_{}.csv'.format(ply_filename))

    #for name in ['KP0294_AC_WP20', 'KP0302_AC_WP20', 'KP0306_LM_WP20', 'KP0477_AC_WP20', 'KP0571_AC_WP20',
    #             'KP0583_LM_WP20', 'KP0588_LM_WP20', 'KP0573_LM_WP20', 'KP0387_AC_WP20', 'KP0518_LM_WP20']:
    #    colony_env = environment[name]
    #    print(colony_env)
    #    print(type(colony_env))
    #    o3d.io.write_point_cloud('{}_env.ply'.format(name), colony_env)  # should save normals & rgb
    #    colony_ind = colonies[name]
    #    print(colony_ind)
    #    print(type(colony_ind))
    #    o3d.io.write_point_cloud('{}.ply'.format(name), colony_ind)  # should save normals & rgb


if __name__ == '__main__':
    # Arguments
    parser = argparse.ArgumentParser(prog="Colony clean and measure")
    parser.add_argument('ply_filename')
    parser.add_argument('annotations_filename')
    parser.add_argument('subsets_filename')
    args = parser.parse_args()
    #args.ply_filename = "cur_kal_20m_20200214_decvis_02.ply"
    #args.annotations_filename = "cur_kal_20m_20200214_decvis_02_KP_16-12-21_completed.txt"
    # args.subsets_filename = "subsets.json"
    ply_filename = args.ply_filename
    annotations_filename = args.annotations_filename
    subsets_filename = args.subsets_filename
    # WP05 cur_kal_05m_20200214_decvis_02_KP
    # theta -6.91, psi 21.08
    # WP10 cur_kal_10m_20200214_decvis_02_KP_905
    # theta = -25.11, psi = 11.65
    # WP20 cur_kal_20m_20200214_decvis_02_KP_16-12-21_completed
    # theta = 9.02, psi = 19.25
    # SB05 cur_sna_05m_20200303_decvis_02_SF_HI_19-1-22
    # theta = -3.90, psi = 12.30
    # SB10 cur_sna_10m_20200303_decvis_02
    # SB10 cur_sna_10m_20201202_decvisann_HI_14-12-21
    # theta = -0.72, psi = -2.71
    # SB20 cur_sna_20m_20200303_decvis_02
    # SB20 cur_sna_20m_20190410_decvisann_HI_12_12
    # theta = -16.82, psi = 18.23
    main(ply_filename, annotations_filename, subsets_filename)