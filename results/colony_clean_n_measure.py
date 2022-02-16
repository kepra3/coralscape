#!/usr/anaconda3/envs/open3d/bin/python
# -*- coding: utf-8 -*-
# remember to edit python environment if required.

"""
@author: kprata
@date created: 9/2/22
@description: TODO: segmenting colony, semantic segmentation?
# For keeping normals/rgb values maybe try directly changing the pcd
pcd = o3d.geometry.PointCloud()
pcd.points = o3d.utility.Vector3dVector(point_cloud[:,:3])
pcd.colors = o3d.utility.Vector3dVector(point_cloud[:,3:6]/255)
pcd.normals = o3d.utility.Vector3dVector(point_cloud[:,6:9])

"""

import argparse
import open3d as o3d
import numpy as np
import copy
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from typing import Any


def create_mesh_ball_pivot(pcd):
    pcd.estimate_normals()
    # NEED TO GET TRUE NORMALS!

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
    # First scale pcd points
    pcd.scale(0.18, center=(0, 0, 0))
    plane_model, inliers = pcd.segment_plane(distance_threshold=0.001,
                                             ransac_n=3,
                                             num_iterations=1000)
    [a, b, c, d] = plane_model
    print(f"Plane equation: {a:.2f}x + {b:.2f}y + {c:.2f}z + {d:.2f} = 0")
    return plane_model, inliers


def plot_plane(pcd, plane_model, scale, z_adjust):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    points = ax.scatter(xs=np.asarray(pcd.points)[:, 0] * scale,
                        ys=np.asarray(pcd.points)[:, 1] * scale,
                        zs=np.asarray(pcd.points)[:, 2] * scale,
                        c='red')
    ax.set_xlabel("Reef parallel")
    ax.set_ylabel("Reef perpendicular")
    ax.set_zlabel("Depth")
    #ax.set_xlim(min(np.asarray(pcd.points)[:, 0]) * scale, max(np.asarray(pcd.points)[:, 0]) * scale)
    #ax.set_ylim(min(np.asarray(pcd.points)[:, 1]) * scale, max(np.asarray(pcd.points)[:, 1]) * scale)
    #ax.set_zlim(min(np.asarray(pcd.points)[:, 2]) * scale, max(np.asarray(pcd.points)[:, 2]) * scale)
    # NEED TO CHOOSE BETTER X AND Y VALUES! Choose based on the inliers
    # what are two furthest points once rotated
    x = np.linspace(min(np.asarray(pcd.points)[:, 0]) * scale, max(np.asarray(pcd.points)[:, 0]) * scale, 10)
    y = np.linspace(min(np.asarray(pcd.points)[:, 1]) * scale, max(np.asarray(pcd.points)[:, 1]) * scale, 10)
    X, Y = np.meshgrid(x, y)
    Z = ((plane_model[0] * X) + (plane_model[1] * Y) + (plane_model[3]*scale) / plane_model[2]*scale) + z_adjust
    surf = ax.plot_surface(X, Y, Z, alpha=0.5)
    print('Showing points and plane ...')
    plt.show()


def rotate_based_on_plane(pcd, plane_model):
    psi = np.arccos(plane_model[0] / (plane_model[0] ** 2 + plane_model[1] ** 2 + plane_model[2] ** 2))
    print(psi * 180 / np.pi)  # about x-axis (yz)
    theta = np.arccos(plane_model[1] / (plane_model[0] ** 2 + plane_model[1] ** 2 + plane_model[2] ** 2))
    print(theta * 180 / np.pi)  # about y-axis (xy)
    R = pcd.get_rotation_matrix_from_xyz((psi, theta, 0))
    pcd_r = copy.deepcopy(pcd)
    pcd_r.rotate(R, center=(0, 0, 0))
    return pcd_r


def get_rugosity(pcd, threeD_area, scale):
    # need to find the largest distances between x and y once oriented (i.e., rotate points then find max xy)
    x_diff = (max(np.asarray(pcd.points)[:, 0]) - min(np.asarray(pcd.points)[:, 0])) * scale
    y_diff = (max(np.asarray(pcd.points)[:, 1]) - min(np.asarray(pcd.points)[:, 1])) * scale
    twoD_area = (x_diff * y_diff)
    print('2D area is ...', twoD_area)
    rugosity = threeD_area / twoD_area
    print("Rugosity is ...", rugosity)
    return rugosity


def main(filename, largest_cluster_mesh):
    # import point cloud
    pcd = o3d.io.read_point_cloud("{}.ply".format(filename))

    o3d.visualization.draw_geometries([pcd])

    # Remove outlier points
    # ind = remove_outlier_points(pcd)
    # display_inlier_outlier(pcd, ind)

    # Create mesh
    mesh = create_mesh_ball_pivot(pcd)
    o3d.visualization.draw_geometries([mesh])  # need to clean points before meshing in some cases maybe
    triangle_clusters, cluster_n_triangles, cluster_area = get_cluster_triangles(mesh)
    mesh_removed = remove_small_clusters(mesh, cluster_n_triangles, triangle_clusters, threshold=5)
    if largest_cluster_mesh == 'Yes':
        large_mesh = largest_cluster(mesh, cluster_n_triangles, triangle_clusters)
        triangle_clusters, cluster_n_triangles, cluster_area = get_cluster_triangles(large_mesh)
        threeD_area = cluster_area * (0.18**2)
    else:
        large_mesh = mesh_removed
        threeD_area = np.sum(cluster_area) * (0.18**2)
        'Print not subsampling to largest mesh'

    print('Cluster area is ... {} m^2'.format(threeD_area))

    # Convert mesh back to point cloud within non-clustered sampling, i.e., using poisson
    print('Sampling points from mesh first uniformaly then with poisson ...')
    pcd = large_mesh.sample_points_uniformly(number_of_points=2500)
    pcd = large_mesh.sample_points_poisson_disk(number_of_points=500, pcl=pcd)
    o3d.visualization.draw_geometries([pcd])

    # Fit a plane to point cloud
    print('Fitting a plane')
    plane_model, inliers = fit_a_plane_ransac(pcd)
    display_inlier_outlier(pcd, inliers)

    # plot plane & point
    plot_plane(pcd, plane_model, scale=0.18)

    # rotate points
    pcd_r = rotate_based_on_plane(pcd, plane_model)
    plot_plane(pcd_r, plane_model, scale=0.18, z_adjust=-1.75)  # trial and error with this!

    # Rugosity
    rugosity = get_rugosity(pcd, threeD_area, scale=0.18)

    # mesh cone to use for caclulating overhang
    # o3d.geometry.create_mesh_cone(radius=1.0, height=2.0, resolution=20, split=1)


if __name__ == '__main__':
    # 287 overhang included but removed with mesh
    # 350 is a large/difficult one, but it actually worked!
    # 479 Mesh isn't connect through whole colony thus largest cluster removes part of colony thus say No
    # 490 is weird! Not sure of real structure here.
    # 554 Mesh doesn't work so well also need point normals, so say No
    # 558
    # Arguments
    parser = argparse.ArgumentParser(prog="Colony clean and measure")
    parser.add_argument('filename')  # e.g., KP0287_LM_WP20
    parser.add_argument('largest_cluster_mesh')  # e.g., Yes or No
    args = parser.parse_args()
    filename = args.filename
    largest_cluster_mesh = args.largest_cluster_mesh
    main(filename, largest_cluster_mesh)