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


def get_rugosity(pcd_r, threeD_area):
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

    overhang_prop = unique_points/len(colony)  # proportion works for now but area would be better
    return overhang_prop


def main(filename, largest_cluster_mesh):
    # import point cloud
    pcd = o3d.io.read_point_cloud("{}.ply".format(filename))
    pcd_env = o3d.io.read_point_cloud("{}_env.ply".format(filename))

    o3d.visualization.draw_geometries([pcd])
    o3d.visualization.draw_geometries([pcd_env])
    # Remove outlier points
    #ind = remove_outlier_points(pcd_env)
    #display_inlier_outlier(pcd_env, ind)  # may need to do this

    # Create mesh
    mesh = create_mesh_ball_pivot(pcd)
    #o3d.visualization.draw_geometries([mesh])  # need to clean points before meshing in some cases maybe
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
    pcd = large_mesh.sample_points_uniformly(number_of_points=len(np.asarray(pcd.points)))
    pcd = large_mesh.sample_points_poisson_disk(number_of_points=len(np.asarray(pcd.points)), pcl=pcd)
    o3d.visualization.draw_geometries([pcd])

    # overhang
    overhang = calc_overhang(pcd, pcd_env)
    print(overhang)


if __name__ == '__main__':
    # 287 overhang included but removed with mesh
    # 350 is a large/difficult one, but it actually worked!
    # 479 Mesh isn't connect through whole colony thus largest cluster removes part of colony thus say No
    # 490 is small
    # 554 works well with normals
    # 558
    # Arguments
    parser = argparse.ArgumentParser(prog="Colony clean and measure")
    parser.add_argument('filename')  # e.g., KP0287_LM_WP20
    parser.add_argument('largest_cluster_mesh')  # e.g., Yes or No
    args = parser.parse_args()
    filename = args.filename
    largest_cluster_mesh = args.largest_cluster_mesh
    main(filename, largest_cluster_mesh)
