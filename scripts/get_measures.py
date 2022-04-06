#!/usr/anaconda3/bin/python
# -*- coding: utf-8 -*-
# remember to edit python environment if required.

"""
@author: kprata
@date created: 6/4/22
@description: TODO
"""

# Modules
import argparse
import numpy as np
import open3d as o3d
import copy

def get_annotations(annotations_path):
    """ Read annotations from csv file """
    annotations = {}
    annotations_file = open(annotations_path, 'r')
    for line in annotations_file:
        cols = line.rstrip().replace(',', ' ').split()
        if cols[0] == 'x':
            continue
        else:
            annotations[cols[0]] = [float(i) for i in cols[1:5]]
    annotations_file.close()
    return annotations


def get_neighbourhood(annotations, pcd, pcd_tree):
    """ Find neighbouring points """
    colonies = {}
    print('Searching ...')
    for name in annotations:
        [k, idx, _] = pcd_tree.search_radius_vector_3d(annotations[name][0:3],
                                                       annotations[name][3])
        # Store colony as separate point cloud with points, normals and colours!!!
        colonies[name] = o3d.geometry.PointCloud()
        colonies[name].points = o3d.utility.Vector3dVector(np.asarray(pcd.points)[idx[1:], :])
        colonies[name].normals = o3d.utility.Vector3dVector(np.asarray(pcd.normals)[idx[1:], :])
        colonies[name].colors = o3d.utility.Vector3dVector(np.asarray(pcd.colors)[idx[1:], :])
    return colonies


def create_mesh_ball_pivot(pcd):
    """ Create a mesh from a point cloud through the ball pivot method """
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
    """ Create a mesh from a point cloud through the poisson method """
    pcd.estimate_normals()
    poisson_mesh = o3d.geometry.TriangleMesh.create_from_point_cloud_poisson(pcd,
                                                                             depth=8,
                                                                             width=0,
                                                                             scale=1.1,
                                                                             linear_fit=False)[0]
    # bbox = pcd.get_axis_aligned_bounding_box()
    # p_mesh_crop = poisson_mesh.crop(bbox)
    return poisson_mesh


def get_cluster_triangles(mesh):
    """ Calculate properties of triangular mesh clusters """
    print("Clustering connected triangles ...")
    with o3d.utility.VerbosityContextManager(o3d.utility.VerbosityLevel.Debug) as cm:
        triangle_clusters, cluster_n_triangles, cluster_area = (mesh.cluster_connected_triangles())
    triangle_clusters = np.asarray(triangle_clusters)
    cluster_n_triangles = np.asarray(cluster_n_triangles)
    cluster_area = np.asarray(cluster_area)
    return triangle_clusters, cluster_n_triangles, cluster_area


def largest_cluster(mesh, cluster_n_triangles, triangle_clusters):
    """ Select the largest cluster from a mesh """
    large_mesh: o3d.cpu.pybind.geometry.TriangleMesh = copy.deepcopy(mesh)
    largest_cluster_idx = cluster_n_triangles.argmax()
    triangles_to_remove = triangle_clusters != largest_cluster_idx
    large_mesh.remove_triangles_by_mask(triangles_to_remove)
    return large_mesh


def fit_a_plane_ransac(pcd):
    """ Fit a plane to a point cloud using the RANSAC method """
    plane_model, inliers = pcd.segment_plane(distance_threshold=0.001,
                                             ransac_n=3,
                                             num_iterations=1000)
    [a, b, c, d] = plane_model
    print(f"Plane equation: {a:.2f}x + {b:.2f}y + {c:.2f}z + {d:.2f} = 0")
    return plane_model, inliers


def calc_plane_angles(plane_model):
    """ Calculate the angles of a plane through the plane equation
        Angles are calculated by finding the angle between two vectors """
    plane_normal = [plane_model[0], plane_model[1], plane_model[2]]
    slope_xz = plane_normal[0] / plane_normal[2]
    slope_yz = plane_normal[1] / plane_normal[2]
    theta = np.arctan(slope_xz) * 180 / np.pi
    print('The angle between x and z is ...', theta)  # about y-axis (xz)
    psi = np.arctan(slope_yz) * 180 / np.pi
    print('The angle between y and z is ...', psi)  # about x-axis (yz)
    # the angle between the x-y plane... i.e., the elevation
    xy_normal = [0, 0, 1]
    mag_xy = np.linalg.norm(xy_normal)
    mag_plane = np.linalg.norm(plane_normal)
    cos_elevation = np.dot(xy_normal, plane_normal) / (mag_xy * mag_plane)
    elevation = np.arccos(cos_elevation) * 180 / np.pi
    print('the angle between plane and x-y plane is ...', elevation)
    return theta.__float__(), psi.__float__(), elevation.__float__()


def calc_rugosity(pcd_r, threeD_area):
    """ Calculate the rugosity of a point cloud, given the 3D area of the mesh and 2D area created by the
    convex hull (alphashape) """
    print('Projecting points to xy plane, z=0')
    x = np.asarray(pcd_r.points)[:, 0]
    y = np.asarray(pcd_r.points)[:, 1]
    z = np.repeat(0, len(np.asarray(pcd_r.points)[:, 1]))
    print('Creating polygon for 2D points')
    points_2d = np.asarray((x, y)).transpose()
    alpha_shape = alphashape.alphashape(points_2d, 2.0)
    rotated_twoD_area = alpha_shape.area
    print('2D area is ...', rotated_twoD_area)
    rugosity = threeD_area / rotated_twoD_area
    if rugosity < 1:
        print("shape not complex rugosity likely essentially 1")
        rugosity = 1
    else:
        pass
    return rugosity, rotated_twoD_area


def calc_overhang(colony_pcd, pcd_env):
    """ Calculate the proportion of a colony is covered by an overhang given surrounding environment points """
    dists = pcd_env.compute_point_cloud_distance(colony_pcd)
    dists = np.asarray(dists)
    ind = np.where(dists < 0.01)[0]
    np.asarray(pcd_env.colors)[ind] = [0, 0, 1]  # green
    colony = np.asarray(pcd_env.points)[ind]
    env = np.asarray(pcd_env.points)
    q = []
    for i in range(0, len(colony)):
        p = np.where((env[:, 0] < colony[i, 0] + 0.01) & (env[:, 0] > colony[i, 0] - 0.01) &
                     (env[:, 1] < colony[i, 1] + 0.01) & (env[:, 1] > colony[i, 1] - 0.01) &
                     (env[:, 2] > max(colony[:, 2])))
        if p[0].size:
            np.asarray(pcd_env.colors)[p[0]] = [1, 0, 0]  # red
            if len(p[0]) > 1:
                for j in range(0, len(p)):
                    p_int = int(p[0][j])
                    q.append(p_int)
            elif len(p[0]) == 1:
                p_int = int(p[0])
                q.append(p_int)
    # Calculate area for overhang
    if colony.size:
        # Calculate area for shadowed colony
        colony_x = colony[:, 0]
        colony_y = colony[:, 1]
        colony_2d = np.asarray((colony_x, colony_y)).transpose()
        alpha_shape_colony = alphashape.alphashape(colony_2d, 18.0)
        twoD_area_colony = alpha_shape_colony.area
        if len(q) > 0:
            unique_int = np.unique(q)
            overhang_x = np.asarray(pcd_env.points)[unique_int, 0]
            overhang_y = np.asarray(pcd_env.points)[unique_int, 1]
            overhang = np.asarray((overhang_x, overhang_y)).transpose()
            alpha_shape_overhang = alphashape.alphashape(overhang, 25.0)
            twoD_area_overhang = alpha_shape_overhang.area
            overhang_prop = twoD_area_overhang / twoD_area_colony
        else:
            print('No overhang')
            overhang_prop = 0.
            twoD_area_overhang = 0.
    else:
        print('issue with colony points ...')
        overhang_prop = None
        twoD_area_overhang = None

    return overhang_prop, twoD_area_colony, twoD_area_overhang


def calc_outcrop(colony_pcd, pcd_env):
    """ Calculate the proportion of height, the mean height of a point cloud sits at in a given environment """
    sample = np.asarray(colony_pcd.points)
    colony_mean_height = np.mean(sample[:, 2])
    colony_min_height = min(sample[:, 2])
    colony_max_height = max(sample[:, 2])
    environment_array = np.asarray(pcd_env.points)
    env_min_height = min(environment_array[:, 2])
    env_max_height = max(environment_array[:, 2])
    env_mean_height = np.mean(environment_array[:, 2])
    env_height_range = env_max_height - env_min_height
    outcrop_prop = (colony_mean_height - env_min_height) / env_height_range
    return outcrop_prop, colony_mean_height, colony_min_height, colony_max_height, \
           env_mean_height, env_min_height, env_max_height


def main(environment_points, annotations):
    """ """

if __name__ == '__main__':
    # Arguments
    parser = argparse.ArgumentParser(prog="Colony clean and measure")
    parser.add_argument('ply_filename')
    parser.add_argument('annotations_filename')
    parser.add_argument('environment_distance', type=float)
    args = parser.parse_args()

    ply_filename = args.ply_filename
    annotations_filename = args.annotations_filename
    path = "/Volumes/KP3/coralscape"
    main(ply_filename, annotations_filename, path)
