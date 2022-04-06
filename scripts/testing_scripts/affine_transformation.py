#!/usr/anaconda3/bin/python
# -*- coding: utf-8 -*-
# remember to edit python environment if required.

"""
@author: kprata
@date created: 19/3/22
@description: TODO
"""
import numpy as np
import open3d as o3d
import copy
import matplotlib.pyplot as plt


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


def get_cluster_triangles(mesh):
    print("Cluster connected triangles")
    with o3d.utility.VerbosityContextManager(o3d.utility.VerbosityLevel.Debug) as cm:
        triangle_clusters, cluster_n_triangles, cluster_area = (mesh.cluster_connected_triangles())
    triangle_clusters = np.asarray(triangle_clusters)
    print(triangle_clusters)
    cluster_n_triangles = np.asarray(cluster_n_triangles)
    cluster_area = np.asarray(cluster_area)
    return triangle_clusters, cluster_n_triangles, cluster_area


def largest_cluster(mesh, cluster_n_triangles, triangle_clusters):
    large_mesh: o3d.cpu.pybind.geometry.TriangleMesh = copy.deepcopy(mesh)
    largest_cluster_idx = cluster_n_triangles.argmax()
    triangles_to_remove = triangle_clusters != largest_cluster_idx
    large_mesh.remove_triangles_by_mask(triangles_to_remove)
    print("Show largest cluster")
    o3d.visualization.draw_geometries([large_mesh])
    return large_mesh


def fit_a_plane_ransac(pcd):
    #TODO: RAMSA plane regression against PCA plane regression
    plane_model, inliers = pcd.segment_plane(distance_threshold=0.001,
                                             ransac_n=3,
                                             num_iterations=1000)
    [a, b, c, d] = plane_model
    return plane_model, inliers


def plot_points_together(pcd, pcd_r):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(xs=np.asarray(pcd.points)[:, 0],
               ys=np.asarray(pcd.points)[:, 1],
               zs=np.asarray(pcd.points)[:, 2],
               c='red')
    ax.set_xlabel("Reef parallel")
    ax.set_ylabel("Reef perpendicular")
    ax.set_zlabel("Depth")
    ax.scatter(xs=np.asarray(pcd_r.points)[:, 0],
               ys=np.asarray(pcd_r.points)[:, 1],
               zs=np.asarray(pcd_r.points)[:, 2],
               c='green')


pcd = o3d.io.read_point_cloud("colony_point_clouds/KP0287_LM_WP20.ply")

o3d.visualization.draw_geometries([pcd])

# Create mesh
mesh = create_mesh_ball_pivot(pcd)
o3d.visualization.draw_geometries([mesh])  # need to clean points before meshing in some cases maybe
triangle_clusters, cluster_n_triangles, cluster_area = get_cluster_triangles(mesh)
large_mesh = largest_cluster(mesh, cluster_n_triangles, triangle_clusters)
triangle_clusters, cluster_n_triangles, cluster_area = get_cluster_triangles(large_mesh)

# Convert mesh back to point cloud within non-clustered sampling, i.e., using poisson
print('Sampling points from mesh first uniformaly then with poisson ...')
pcd = large_mesh.sample_points_uniformly(number_of_points=2500)
pcd = large_mesh.sample_points_poisson_disk(number_of_points=500, pcl=pcd)
o3d.visualization.draw_geometries([pcd])


# Get plane
plane_model, inliers = fit_a_plane_ransac(pcd)

center = pcd.get_center()
z_n = [plane_model[0], plane_model[1], plane_model[2]]
random_vector = [1, 1, 1]
# No matter what the random vector is x_n will lie in the plane normal to z_n
x_n = np.cross(random_vector, z_n) / np.linalg.norm(np.cross(random_vector, z_n))
# then in a right handed system with x_n as the x-axis and z_n as pointing up,
# y_n will be to the left of x_n and the right of z_n
y_n = np.cross(z_n, x_n)
R = np.array([x_n, y_n, z_n])
pcd_r = copy.deepcopy(pcd)
pcd_r.rotate(R, center=center)
o3d.visualization.draw_geometries([pcd, pcd_r])
plot_points_together(pcd, pcd_r)
