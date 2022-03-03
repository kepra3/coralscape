#!/usr/anaconda3/bin/python
# -*- coding: utf-8 -*-
# remember to edit python environment if required.

"""
@author: kprata
@date created: 21/2/22
@description: TODO
"""

import argparse
import open3d as o3d
import numpy as np
import copy
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from typing import Any
import alphashape
from descartes import PolygonPatch
from sklearn.linear_model import LinearRegression


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
    print("Statistical oulier removal")
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
    print("Cluster connected triangles")
    with o3d.utility.VerbosityContextManager(o3d.utility.VerbosityLevel.Debug) as cm:
        triangle_clusters, cluster_n_triangles, cluster_area = (mesh.cluster_connected_triangles())
    triangle_clusters = np.asarray(triangle_clusters)
    print(triangle_clusters)
    cluster_n_triangles = np.asarray(cluster_n_triangles)
    cluster_area = np.asarray(cluster_area)
    return triangle_clusters, cluster_n_triangles, cluster_area


def remove_small_clusters(mesh, cluster_n_triangles, triangle_clusters, threshold=5):
    mesh_removed = copy.deepcopy(mesh)
    triangles_to_remove = cluster_n_triangles[triangle_clusters] < threshold
    mesh_removed.remove_triangles_by_mask(triangles_to_remove)
    print("Show mesh with small clusters removed")
    o3d.visualization.draw_geometries([mesh_removed])
    return mesh_removed


def largest_cluster(mesh, cluster_n_triangles, triangle_clusters):
    large_mesh: o3d.cpu.pybind.geometry.TriangleMesh = copy.deepcopy(mesh)
    largest_cluster_idx = cluster_n_triangles.argmax()
    triangles_to_remove = triangle_clusters != largest_cluster_idx
    large_mesh.remove_triangles_by_mask(triangles_to_remove)
    print("Show largest cluster")
    o3d.visualization.draw_geometries([large_mesh])
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
    print('Theta is ...', theta)
    return theta.__float__()


def calc_plane_angles(plane_model):
    # using the formula
    # plane_normal dot axis of interest normal / (magnitude of plane normal * magnitude of axis of interest normal)
    # which simplifies to
    # xz_normal = [0, 1, 0] have to figure out why these don't make sense
    # yz_normal = [1, 0, 0]
    plane_normal = [plane_model[0], plane_model[1], plane_model[2]]
    # mag_xz = np.linalg.norm(xz_normal)
    # mag_yz = np.linalg.norm(yz_normal)
    # mag_plane = np.linalg.norm(plane_normal)
    # cos_theta = np.dot(xz_normal, plane_normal)/(mag_xz*mag_plane)
    # cos_psi = np.dot(yz_normal, plane_normal)/(mag_yz*mag_plane)
    # theta = np.arccos(cos_theta) * 180./np.pi
    slope_xz = plane_normal[0] / plane_normal[2]
    slope_yz = plane_normal[1] / plane_normal[2]
    theta = np.arctan(slope_xz) * 180 / np.pi
    # if theta > 90:
    #    theta = 180 - theta
    # else:
    #    theta = theta
    print('The angle between x and z is ...', theta)  # about y-axis (xz)
    # psi = np.arccos(cos_psi) * 180./np.pi
    psi = np.arctan(slope_yz) * 180 / np.pi
    # if psi > 90:
    #    psi = 180 - 90
    # else:
    #    psi = psi
    print('The angle between y and z is ...', psi)  # about x-axis (yz)
    # the angle between the x-y plane... i.e., the elevation
    xy_normal = [0, 0, 1]
    mag_xy = np.linalg.norm(xy_normal)
    mag_plane = np.linalg.norm(plane_normal)
    cos_elevation = np.dot(xy_normal, plane_normal) / (mag_xy * mag_plane)
    elevation = np.arccos(cos_elevation) * 180 / np.pi
    print('the angle between plane and x-y plane is ...', elevation)
    return theta.__float__(), psi.__float__(), elevation


def calc_rugosity(pcd_r, threeD_area):
    # project points onto 2D xy axes
    print('Projecting points to xy plane, z=0')
    x = np.asarray(pcd_r.points)[:, 0]
    y = np.asarray(pcd_r.points)[:, 1]
    z = np.repeat(0, len(np.asarray(pcd_r.points)[:, 1]))
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    points = ax.scatter(xs=x,
                        ys=y,
                        zs=z,
                        c='red')
    ax.set_xlabel("Reef parallel")
    ax.set_ylabel("Reef perpendicular")
    ax.set_zlabel("Depth")
    x1 = np.linspace(min(np.asarray(pcd_r.points)[:, 0]), max(np.asarray(pcd_r.points)[:, 0]), 10)
    y1 = np.linspace(min(np.asarray(pcd_r.points)[:, 1]), max(np.asarray(pcd_r.points)[:, 1]), 10)
    X, Y = np.meshgrid(x1, y1)
    Z = np.zeros((10, 10))
    surf = ax.plot_surface(X, Y, Z, alpha=0.5)
    plt.show()

    print('Creating polygon for 2D points')
    points_2d = np.asarray((x, y)).transpose()
    fig, ax = plt.subplots()
    ax.scatter(x=points_2d[:, 0], y=points_2d[:, 1])
    alpha_shape = alphashape.alphashape(points_2d, 2.0)
    ax.add_patch(PolygonPatch(alpha_shape, alpha=0.2))
    plt.show()
    twoD_area = alpha_shape.area * (scale ** 2)
    print('2D area is ...', twoD_area)

    # print('Creating polygon for 3D points')
    # points_3d = np.asarray(pcd_r.points)
    # fig = plt.figure()
    # ax = plt.axes(projection='3d')
    # ax.scatter(xs=points_3d[:, 0],
    #           ys=points_3d[:, 1],
    #           zs=points_3d[:, 2],
    #           c='red')
    # ax.set_xlabel("Reef parallel")
    # ax.set_ylabel("Reef perpendicular")
    # ax.set_zlabel("Depth")
    # alpha_shape = alphashape.alphashape(points_3d, 1.1)
    # ax.plot_trisurf(*zip(*alpha_shape.vertices), triangles=alpha_shape.faces)
    # plt.show()
    # threeD_area_2 = alpha_shape.area * (scale ** 2)
    # print('3D area is ...', threeD_area_2)
    print('Rugosity is', threeD_area / twoD_area)
    # print('rugosity 2 is', threeD_area_2 / twoD_area)
    rugosity = threeD_area / twoD_area
    if rugosity < 1:
        print("shape complex rugosity likely essentially 1")
        rugosity = 1
    else:
        pass
    return rugosity


def calc_overhang(pcd, pcd_env):
    dists = pcd_env.compute_point_cloud_distance(pcd)
    dists = np.asarray(dists)
    ind = np.where(dists < 0.01)[0]
    not_ind = np.where(dists > 0.01)[0]
    np.asarray(pcd_env.colors)[ind] = [0, 0, 1]
    o3d.visualization.draw_geometries([pcd_env])
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

    o3d.visualization.draw_geometries([pcd_env])
    unique_q = list()
    unique_points = 0
    for item in q:
        if item not in unique_q:
            unique_q.append(item)
            unique_points += 1
    print(unique_q)

    overhang_prop = unique_points / len(colony)  # proportion works for now but area would be better
    return overhang_prop


def calc_overhang_mainpy(colony_pcd, pcd_env):
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


def calc_outcrop_proportion(colonies, environment):
    """ Calculate the proportion the colony sits at in an environment """
    heights = {}
    all_relative_heights = {}

    for name in colonies:
        sample = np.asarray(colonies[name].points)
        if len(sample) <= 3:
            continue
        colony_mean_height = np.mean(sample[:, 2])
        heights[name] = colony_mean_height

        # Extract points
        environment_array = np.asarray(environment[name].points)
        env_min_height = min(environment_array[:, 2])
        env_max_height = max(environment_array[:, 2])
        env_height_range = env_max_height - env_min_height

        # Proportion
        all_relative_heights[name] = (colony_mean_height - env_min_height) / env_height_range

    return all_relative_heights


def calc_environment_rugosity(environment):
    environment_rugosity = {}
    for name in environment:
        mesh = create_mesh_ball_pivot(environment[name])
        triangle_clusters, cluster_n_triangles, cluster_area = get_cluster_triangles(mesh)
        threeD_area = np.sum(cluster_area)
        print('Cluster area is ... {} m^2 for {}'.format(threeD_area, name))
        print('Sampling points from mesh first uniformaly then with poisson ...')
        environment_pcd = mesh.sample_points_uniformly(number_of_points=len(np.asarray(environment[name].points)))
        environment_pcd = mesh.sample_points_poisson_disk(number_of_points=len(np.asarray(environment[name].points)),
                                                          pcl=environment_pcd)
        print('Fitting a plane to correct for environment slope')
        plane_model_env, inliers_env = fit_a_plane_ransac(environment_pcd)
        env_theta_xz, env_psi_yz = calc_plane_angles(plane_model_env)
        env_theta_radians = env_theta_xz / (180 * np.pi)
        env_psi_radians = env_psi_yz / (180 * np.pi)
        env_R = environment_pcd.get_rotation_matrix_from_xyz((env_psi_radians, env_theta_radians, 0))
        center = environment_pcd.get_center()
        environment_pcd.rotate(env_R, center=center)
        rugosity = calc_rugosity(environment_pcd, threeD_area)
        environment_rugosity[name] = rugosity
    return environment_rugosity


def plot_plane(pcd, plane_model, inliers, z_adjust):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(xs=np.asarray(pcd.points)[:, 0],
               ys=np.asarray(pcd.points)[:, 1],
               zs=np.asarray(pcd.points)[:, 2],
               c='red')
    ax.set_xlabel("Reef parallel")
    ax.set_ylabel("Reef perpendicular")
    ax.set_zlabel("Depth")
    # ax.scatter(xs=np.asarray(pcd_r.points)[:, 0],
    #           ys=np.asarray(pcd_r.points)[:, 1],
    #           zs=np.asarray(pcd_r.points)[:, 2],
    #           c='green')
    # ax.scatter(xs=np.asarray(pcd_c.points)[:, 0],
    #           ys=np.asarray(pcd_c.points)[:, 1],
    #           zs=np.asarray(pcd_c.points)[:, 2],
    #           c='blue')
    plane_points = np.asarray(pcd.points)[inliers]
    x = plane_points[:, 0]
    y = plane_points[:, 1]
    X, Y = np.meshgrid(x, y)
    Z = plane_model[3] - (plane_model[0] * X + plane_model[1] * Y) / plane_model[2] + z_adjust
    surf = ax.plot_surface(X, Y, Z, alpha=0.5)
    print('Showing points and plane ...')
    plt.show()


def plot_points_together(pcd, pcd_r, sample_name, number):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(xs=np.asarray(pcd.points)[:, 0],
               ys=np.asarray(pcd.points)[:, 1],
               zs=np.asarray(pcd.points)[:, 2],
               c='red')
    ax.set_xlabel("Reef parallel")
    ax.set_ylabel("Reef perpendicular")
    ax.set_zlabel("Depth")
    ax.set_title(sample_name)
    ax.scatter(xs=np.asarray(pcd_r.points)[:, 0],
               ys=np.asarray(pcd_r.points)[:, 1],
               zs=np.asarray(pcd_r.points)[:, 2],
               c='green')
    # ax.scatter(xs=np.asarray(pcd_c.points)[:, 0],
    #           ys=np.asarray(pcd_c.points)[:, 1],
    #           zs=np.asarray(pcd_c.points)[:, 2],
    #           c='blue')
    plt.savefig('../colony_point_clouds/colony_pics/{}_{}_rotated.png'.format(sample_name, number))


def plot_points_together_colours(pcd, pcd_r, sample_name, number):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(xs=np.asarray(pcd.points)[:, 0],
               ys=np.asarray(pcd.points)[:, 1],
               zs=np.asarray(pcd.points)[:, 2],
               c=np.asarray(pcd.colors))
    ax.set_xlabel("Reef parallel")
    ax.set_ylabel("Reef perpendicular")
    ax.set_zlabel("Depth")
    ax.set_title(sample_name)
    ax.scatter(xs=np.asarray(pcd_r.points)[:, 0],
               ys=np.asarray(pcd_r.points)[:, 1],
               zs=np.asarray(pcd_r.points)[:, 2],
               c=np.asarray(pcd_r.colors))
    # ax.scatter(xs=np.asarray(pcd_c.points)[:, 0],
    #           ys=np.asarray(pcd_c.points)[:, 1],
    #           zs=np.asarray(pcd_c.points)[:, 2],
    #           c='blue')
    plt.savefig('../colony_point_clouds/colony_pics/{}_{}_colour.png'.format(sample_name, number))


def main(sample_name, largest_cluster_mesh):
    # Make results file or append ot it
    out_name = "../colony_point_clouds/results.txt"
    with open(out_name, 'a') as results_out:
        if results_out.tell() == 0:
            print('Creating a new file\n')
            results_out.write(
                "sample_name\tcloud_points\tplane_i,\tplane_j\tplane_k\taxis_i1\taxis_j1\taxis_k1\t"
                "axis_i2\taxis_j2\taxis_k2\televation\n")
        else:
            print('File exists, appending\n')
    # Import point cloud
    pcd = o3d.io.read_point_cloud("../colony_point_clouds/{}.ply".format(sample_name))
    #pcd_env = o3d.io.read_point_cloud("colony_point_clouds/{}_env.ply".format(sample_name))

    o3d.visualization.draw_geometries([pcd])
    cloud_points = len(np.asarray(pcd.points))
    #o3d.visualization.draw_geometries([pcd_env])
    # Remove outlier points
    # ind = remove_outlier_points(pcd_env)
    # display_inlier_outlier(pcd_env, ind)  # may need to do this

    # Create mesh
    mesh = create_mesh_ball_pivot(pcd)
    # o3d.visualization.draw_ geometries([mesh])  # need to clean points before meshing in some cases maybe
    triangle_clusters, cluster_n_triangles, cluster_area = get_cluster_triangles(mesh)
    mesh_removed = remove_small_clusters(mesh, cluster_n_triangles, triangle_clusters, threshold=5)
    if largest_cluster_mesh == 'Yes':
        large_mesh = largest_cluster(mesh, cluster_n_triangles, triangle_clusters)
        triangle_clusters, cluster_n_triangles, cluster_area = get_cluster_triangles(large_mesh)
        threeD_area = cluster_area
    else:
        large_mesh = mesh_removed
        threeD_area = np.sum(cluster_area)
        'Print not subsampling to largest mesh'

    print('Cluster area is ... {} m^2'.format(threeD_area))

    # Convert mesh back to point cloud within non-clustered sampling, i.e., using poisson
    print('Sampling points from mesh first uniformaly then with poisson ...')
    pcd = large_mesh.sample_points_uniformly(number_of_points=2500)
    pcd = large_mesh.sample_points_poisson_disk(number_of_points=500, pcl=pcd)
    o3d.visualization.draw_geometries([pcd])

    # overhang
    # overhang = calc_overhang(pcd, pcd_env)
    # print(overhang)

    #  Fit plane ransac
    plane_model, inliers = fit_a_plane_ransac(pcd)
    theta, psi, elevation = calc_plane_angles(plane_model)
    theta_radians = theta / 180 * np.pi
    psi_radians = psi / 180 * np.pi
    plot_plane(pcd, plane_model, inliers, -20.75)
    center = pcd.get_center()
    pcd_r1 = copy.deepcopy(pcd)
    pcd_r2 = copy.deepcopy(pcd)
    pcd_r3 = copy.deepcopy(pcd)
    n_z = [0, 0, 1]  # xy plane normal
    n_p = [plane_model[0], plane_model[1], plane_model[2]]  # plane normal
    # get the axis of rotation vector through the cross product
    axis_R1 = np.cross(n_z, n_p)  # n_z pointing finger, n_z middle
    axis_R2 = np.cross(n_p, n_z)  # n_p pointing finger, n_z middle finger
    R1 = pcd.get_rotation_matrix_from_axis_angle(axis_R1)
    R2 = pcd.get_rotation_matrix_from_axis_angle(axis_R2)
    R3 = pcd.get_rotation_matrix_from_axis_angle(n_z)
    pcd_r1.rotate(R1, center=center)
    pcd_r2.rotate(R2, center=center)
    pcd_r3.rotate(R3, center=center)
    plot_points_together(pcd, pcd_r1, sample_name, 1)
    plot_points_together(pcd, pcd_r2, sample_name, 2)
    plot_points_together(pcd, pcd_r3, sample_name, 3)
    print('Sample:', sample_name, 'plane model', plane_model, '\naxis of rotation', axis_R1, '\nangle of elevation',
          elevation)
    print('axis of rotation two:', axis_R2, '\n')
    plot_points_together_colours(pcd, pcd_r1, sample_name, 1)
    plot_points_together_colours(pcd, pcd_r2, sample_name, 2)
    plot_points_together_colours(pcd, pcd_r3, sample_name, 3)
    print('Outputting results ... :)')
    with open(out_name, "a") as results_out:
        results_out.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}"
                          "\t{8}\t{9}\t{10}\t{11}\n".format(sample_name, cloud_points, n_p[0], n_p[1], n_p[2],
                                                            axis_R1[0], axis_R1[1], axis_R1[2], axis_R2[0], axis_R2[1],
                                                            axis_R2[2], elevation))


if __name__ == '__main__':
    # 287 overhang included but removed with mesh
    # 350 is a large/difficult one, but it actually worked!
    # 479 Mesh isn't connect through whole colony thus largest cluster removes part of colony thus say No
    # 490 is small
    # 554 works well with normals
    # 558
    # Arguments
    # parser = argparse.ArgumentParser(prog="Colony clean and measure")
    # parser.add_argument('filename')  # e.g., KP0287_LM_WP20
    # parser.add_argument('largest_cluster_mesh')  # e.g., Yes or No
    # args = parser.parse_args()
    # filename = args.filename
    largest_cluster_mesh = "Yes"
    for sample in ['KP0287_LM_WP20', 'KP0294_AC_WP20', 'KP0302_AC_WP20', 'KP0306_LM_WP20', 'KP0350_LM_WP20',
                   'KP0387_AC_WP20', 'KP0477_AC_WP20', 'KP0479_AC_WP20', 'KP0490_AC_WP20', 'KP0518_LM_WP20',
                   'KP0554_AC_WP20', 'KP0558_LM_WP20', 'KP0571_AC_WP20', 'KP0573_LM_WP20', 'KP0583_LM_WP20',
                   'KP0588_LM_WP20']:
        sample_name = sample

        main(sample, largest_cluster_mesh)
    # largest_cluster_mesh = args.largest_cluster_mesh

