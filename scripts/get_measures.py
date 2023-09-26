#!/usr/anaconda3/bin/python
# -*- coding: utf-8 -*-
# remember to edit python environment if required.

"""
@author: kprata
@date created: 6/4/22
@description: Calculates all measures for a single colony.
"""

# Modules
import argparse
import matplotlib.pyplot as plt
import numpy as np
import open3d as o3d
import copy
import alphashape
from descartes import PolygonPatch
import os


def get_annotations(annotations_path):
    """ Read annotations from csv file """
    annotations = {}
    annotations_file = open(annotations_path, 'r')
    for line in annotations_file:
        cols = line.rstrip().replace(',', ' ').split()
        if cols[1] == 'x':
            continue
        else:
            annotations[cols[0]] = [float(i) for i in cols[1:5]]
    annotations_file.close()
    return annotations


def get_neighbourhood(annotations, ply_filename, pcd, pcd_tree):
    """ Find neighbouring points """
    print('Searching ...')
    if annotations[ply_filename][3] <= 0.:
        annotations[ply_filename][3] = 0.01
    [k, idx, _] = pcd_tree.search_radius_vector_3d(annotations[ply_filename][0:3],
                                                   annotations[ply_filename][3])
    # Store colony as separate point cloud with points, normals and colours!!!
    colony = o3d.geometry.PointCloud()
    colony.points = o3d.utility.Vector3dVector(np.asarray(pcd.points)[idx[1:], :])
    colony.normals = o3d.utility.Vector3dVector(np.asarray(pcd.normals)[idx[1:], :])
    colony.colors = o3d.utility.Vector3dVector(np.asarray(pcd.colors)[idx[1:], :])
    return colony


def find_ground_points(past_colony_pcd, colony):
    dists = past_colony_pcd.compute_point_cloud_distance(colony)
    dists = np.asarray(dists)
    ind = np.where(dists != 0)[0]
    ground_pcd = o3d.geometry.PointCloud()
    ground_pcd.points = o3d.utility.Vector3dVector(np.asarray(past_colony_pcd.points)[ind[0:], :])
    ground_pcd.normals = o3d.utility.Vector3dVector(np.asarray(past_colony_pcd.normals)[ind[0:], :])
    ground_pcd.colors = o3d.utility.Vector3dVector(np.asarray(past_colony_pcd.colors)[ind[0:], :])
    return ground_pcd


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


def mesh_alpha(pcd, alpha):
    print(f"alpha={alpha:.3f}")
    mesh = o3d.geometry.TriangleMesh.create_from_point_cloud_alpha_shape(pcd, alpha)
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


def calc_rugosity(pcd_r, threeD_area, name):
    """ Calculate the rugosity of a point cloud, given the 3D area of the mesh and 2D area created by the
    convex hull (alphashape) """
    print('Projecting points to xy plane, z=0')
    x = np.asarray(pcd_r.points)[:, 0]
    y = np.asarray(pcd_r.points)[:, 1]
    z = np.repeat(0, len(np.asarray(pcd_r.points)[:, 1]))
    print('Creating polygon for 2D points')
    points_2d = np.asarray((x, y)).transpose()
    alpha_shape = alphashape.alphashape(points_2d, 2.0)
    fig, ax = plt.subplots()
    ax.scatter(x=points_2d[:, 0], y=points_2d[:, 1])
    ax.add_patch(PolygonPatch(alpha_shape, alpha=0.2))
    name_list = name.split("_")
    if len(name_list) == 5:
        try:
            os.makedirs("{}/results/{}/polygon".format(path, "_".join(name_list[3:5])), exist_ok=True)
            print("New directory created, '%s'" % "polygon")
        except OSError as error:
            print("Directory '%s' can not be created" % "polygon")
        plt.savefig("{}/results/{}/polygon/polygonPatch_{}.png".format(path, "_".join(name_list[3:5]), name))
    else:
        plt.savefig("{}/results/colony_pics/polygonPatch_{}.png".format(path, name))
    rotated_twoD_area = alpha_shape.area
    print('2D area is ...', rotated_twoD_area)
    rugosity = threeD_area / rotated_twoD_area
    if rugosity < 1:
        print("shape not complex rugosity likely essentially 1")
        rugosity = 1
    else:
        pass
    return rugosity, rotated_twoD_area


def calc_overhang(colony_pcd, pcd_env, name):
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
        unique_int = np.unique(q)
        if len(unique_int) > 1:
            overhang_x = np.asarray(pcd_env.points)[unique_int, 0]
            overhang_y = np.asarray(pcd_env.points)[unique_int, 1]
            overhang = np.asarray((overhang_x, overhang_y)).transpose()
            alpha_shape_overhang = alphashape.alphashape(overhang, 25.0)
            fig, ax = plt.subplots()
            ax.scatter(x=colony_x, y=colony_y, color='blue')
            ax.scatter(x=overhang_x, y=overhang_y, color='red')
            ax.add_patch(PolygonPatch(alpha_shape_colony, alpha=0.2, color='blue'))
            ax.add_patch(PolygonPatch(alpha_shape_overhang, alpha=0.2, color='red'))
            name_list = name.split("_")
            try:
                os.makedirs("{}/results/{}/overhang".format(path, "_".join(name_list[3:5])), exist_ok=True)
                print("New directory created, '%s'" % "overhang")
            except OSError as error:
                print("Directory '%s' can not be created" % "overhang")
            plt.savefig('{}/results/{}/overhang/overhang_cover_{}.png'.format(path, "_".join(name_list[3:5]), name))
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


def main(ply_filename, annotations_filename, path):
    """ Getting structural complexity measures for a colony """
    name_split = ply_filename.split("_")
    environment_distance = name_split[4]
    name = "_".join(name_split[0:3])
    env = "_".join(name_split[3:5])
    # Create result files these will be saved
    try:
        os.makedirs("{}/results/{}".format(path, env), exist_ok=True)
        print("New directory created, '%s'" % env)
    except OSError as error:
        print("Directory '%s' can not be created" % env)
    all_results = "{}/results/{}/struc_complex_results_{}.txt".format(path, env, ply_filename)
    with open(all_results, 'a') as results_out:
        print("Creating a new file\n")
        results_out.write(
            "plot_name\tsample_name\tcolony_points\tcolony_elevation\tground_elevation\t"
            "colony_rugosity\toverhang_prop\toutcrop_prop"
            "\tenvironment_rugosity\tcolony_range\tenvironment_range\n")
    detailed_results = "{}/results/{}/detailed_results_{}.txt".format(path, env, ply_filename)
    with open(detailed_results, 'a') as results_out:
        results_out.write(
            "plot_name\tsample_name\tcolony_points"
            "\ttheta_ground\tpsi_ground\televation_ground"
            "\ttheta\tpsi\televation"
            "\tx_i\t_x_j\tx_k\ty_i\ty_j\ty_k\tz_i\tz_j\tz_k"
            "\tcolony_threeD_area\trot_colony_twoD_area\tcolony_rugosity"
            "\tcolony_2D_area\toverhang_2D_area\toverhang_prop"
            "\tmean_colony\tlow_colony\thigh_colony\tmean_env\tlow_env\thigh_env\toutcrop_prop"
            "\tmean_ground\tlow_ground\thigh_ground\toutcrop_prop_ground"
            "\tenv_points\tenv_threeD_area\tenv_twoD_area\tenv_rugosity\n")

    # 1. Read in files
    print('Reading PLY file {} ...'.format(name))
    pcd = o3d.io.read_point_cloud('{}/data/environment_points/{}.ply'.format(path, ply_filename))
    print('Read assignment file ...')
    annotations = get_annotations('{}/data/{}'.format(path, annotations_filename))
    # Make KDTree
    print('Making KDTree ...')
    pcd_tree = o3d.geometry.KDTreeFlann(pcd)

    # 2. Find colony within environment
    print("\nSearching for colony around annotations ...")
    colony = get_neighbourhood(annotations, name, pcd, pcd_tree)
    colony_length = len(np.asarray(colony.points))
    if colony_length <= 10:
        print('Not including {} because < 10 points'.format(name))
    else:
        with open(detailed_results, 'a') as results_out:
            results_out.write("{0}\t{1}\t{2}".format(name_split[2], name, colony_length))
        # 3. Find ground past colony
        annotations_two = copy.deepcopy(annotations)
        annotations_two[name][3] = annotations_two[name][3] + 0.04  # maybe 4cm?
        past_colony_pcd = get_neighbourhood(annotations_two, name, pcd, pcd_tree)
        ground_pcd = find_ground_points(past_colony_pcd, colony)

        # 4. Mesh the colony to remove outlier points
        print('Meshing ...')
        mesh = create_mesh_ball_pivot(colony)
        # o3d.visualization.draw_geometries([mesh])
        triangle_clusters, cluster_n_triangles, cluster_area = get_cluster_triangles(mesh)
        large_mesh = largest_cluster(mesh, cluster_n_triangles, triangle_clusters)
        # o3d.visualization.draw_geometries([large_mesh])
        print('Sampling points from mesh first uniformly then with poisson ...')
        colony_pcd = large_mesh.sample_points_uniformly(number_of_points=colony_length)
        colony_pcd = large_mesh.sample_points_poisson_disk(number_of_points=colony_length, pcl=colony_pcd)
        # maybe need an estimate point normals here.
        # alpha_mesh = mesh_alpha(colony_pcd, 0.01) not working well...
        # o3d.visualization.draw_geometries([alpha_mesh])
        # vis = o3d.visualization.Visualizer()
        # vis.create_window(visible=True)  # does not work as false for me
        # vis.add_geometry(alpha_mesh)
        # vis.update_geometry(alpha_mesh)
        # vis.poll_events()
        # vis.update_renderer()
        # vis.capture_screen_image('{}/results/colony_pics/{}.png'.format(path, name))
        # vis.destroy_window()
        triangle_clusters, cluster_n_triangles, cluster_area = get_cluster_triangles(large_mesh)
        colony_threeD_area = np.sum(cluster_area).__float__()
        print('Cluster area is ... {} m^2 for {}'.format(colony_threeD_area, name))

        # 5. Calculating colony angles and ground slope angles
        print('Fitting a plane with RAMSAC to get colony angles ...')
        plane_model, inliers = fit_a_plane_ransac(colony_pcd)
        print('Getting colony anlges ...')
        colony_theta_xz, colony_psi_yz, colony_elevation = calc_plane_angles(plane_model)
        print('Colony angles: theta {}, psi {}, elevation {}'.format(colony_theta_xz, colony_psi_yz,
                                                                     colony_elevation))
        plane_model_g, inliers_g = fit_a_plane_ransac(ground_pcd)
        print('Getting ground angles ...')
        ground_theta_xz, ground_psi_yz, ground_elevation = calc_plane_angles(plane_model_g)
        print('Ground angles: theta {}, psi {}, elevation {}'.format(ground_theta_xz, ground_psi_yz, ground_elevation))
        with open(detailed_results, 'a') as results_out:
            results_out.write('\t{0}\t{1}\t{2}'.format(ground_theta_xz, ground_psi_yz, ground_elevation))

        # 6. Calculating colony rugosity
        print('Getting colony rugosity ...')
        print('Rotate colony points according to colony slope')
        center = colony_pcd.get_center()
        z_n = [plane_model[0], plane_model[1], plane_model[2]]
        random_vector = [1, 1, 1]
        # No matter what the random vector is x_n will lie in the plane normal to z_n
        x_n = np.cross(random_vector, z_n) / np.linalg.norm(np.cross(random_vector, z_n))
        # Because of the right handed system y_n will be to the left of x_n and the right of z_n
        # And orthogonal to both
        y_n = np.cross(z_n, x_n)
        # This creates a new coordinate system with the z_n normal to the x-y plane
        R = np.array([x_n, y_n, z_n])
        with open(detailed_results, 'a') as results_out:
            results_out.write("\t{0}\t{1}\t{2}"
                              "\t{3}\t{4}\t{5}"
                              "\t{6}\t{7}\t{8}"
                              "\t{9}\t{10}\t{11}".format(colony_theta_xz, colony_psi_yz, colony_elevation,
                                                         x_n[0], x_n[1], x_n[2],
                                                         y_n[0], y_n[1], y_n[2],
                                                         z_n[0], z_n[1], z_n[2]))
        colony_pcd_r = copy.deepcopy(colony_pcd)
        colony_pcd_r.rotate(R, center=center)
        colony_rugosity, rot_colony_twoD_area = calc_rugosity(colony_pcd_r, colony_threeD_area, name)
        print('Rugosity is ...', colony_rugosity)
        with open(detailed_results, 'a') as results_out:
            results_out.write("\t{0}\t{1}\t{2}".format(colony_threeD_area, rot_colony_twoD_area, colony_rugosity))

        # 7. Calculate overhang proportion
        print('Getting overhang proportion ...')
        overhang_prop, twoD_area_colony, twoD_area_overhang = calc_overhang(colony_pcd, pcd, ply_filename)
        print('Overhang proportion is ...', overhang_prop)
        with open(detailed_results, 'a') as results_out:
            results_out.write("\t{0}\t{1}\t{2}".format(twoD_area_colony, twoD_area_overhang, overhang_prop))

        # 8. Calculate outcrop proportion
        print('Getting outcrop proportion ...')
        outcrop_prop, colony_mean_height, colony_min_height, colony_max_height, \
        env_mean_height, env_min_height, env_max_height = calc_outcrop(colony_pcd, pcd)
        print('Outcrop proportion is ...', outcrop_prop)
        with open(detailed_results, 'a') as results_out:
            results_out.write("\t{0}\t{1}\t{2}\t{3}\t{4}"
                              "\t{5}\t{6}".format(colony_mean_height, colony_min_height,
                                                  colony_max_height, env_mean_height,
                                                  env_min_height, env_max_height, outcrop_prop))
        # use this for outcrop_prop...
        outcrop_prop_ground, colony_mean_height, colony_min_height, colony_max_height, \
        ground_mean_height, ground_min_height, ground_max_height = calc_outcrop(colony_pcd, ground_pcd)
        with open(detailed_results, 'a') as results_out:
            results_out.write("\t{0}\t{1}\t{2}\t{3}".format(ground_mean_height,
                                                            ground_min_height, ground_max_height, outcrop_prop_ground))

        # 9. Calculate environment rugosity
        print('Getting environment rugosity ...')
        mesh = create_mesh_ball_pivot(pcd)
        triangle_clusters, cluster_n_triangles, cluster_area = get_cluster_triangles(mesh)
        env_threeD_area = np.sum(cluster_area)
        print('Cluster area is ... {} m^2 for {}'.format(env_threeD_area, name))
        print('Sampling points from mesh first uniformly then with poisson ...')
        environment_pcd = mesh.sample_points_uniformly(number_of_points=len(np.asarray(pcd.points)))
        environment_pcd = mesh.sample_points_poisson_disk(number_of_points=len(np.asarray(pcd.points)),
                                                          pcl=environment_pcd)
        env_center = environment_pcd.get_center()

        if name_split[2] == "WP05":
            elevation = 5.35
        elif name_split[2] == "WP10":
            elevation = 14.56
        elif name_split[2] == "WP20":
            elevation = 9.88
        elif name_split[2] == "CA05":
            elevation = 9.66
        elif name_split[2] == "CA10":
            elevation = 3.53
        elif name_split[2] == "CA20":
            elevation = 35.06
        elif name_split[2] == "SB05":
            elevation = 5.70
        elif name_split[2] == "SB10":
            elevation = 8.85
        elif name_split[2] == "SB20":
            elevation = 34.35
        elif name_split[2] == "SQ12":
            elevation = 12.23
        elif name_split[2] == "SQ20":
            elevation = 14.97
        else:
            print("Elevation not defined, find elevation of plot in plot_info.txt")
        elevation_radians = elevation / 180 * np.pi
        env_R = environment_pcd.get_rotation_matrix_from_axis_angle([0, -elevation_radians, 0])
        # TODO: need find direction vector in relation to elevation!
        environment_pcd.rotate(env_R, center=env_center)
        environment_rugosity, env_twoD_area = calc_rugosity(environment_pcd, env_threeD_area, ply_filename)
        print('Environment rugosity is  ... {}'.format(environment_rugosity))
        with open(detailed_results, 'a') as results_out:
            results_out.write("{0}\t{1}\t{2}\t{3}\n".format(len(np.asarray(environment_pcd.points)),
                                                            env_threeD_area, env_twoD_area,
                                                            environment_rugosity))
        with open(all_results, 'a') as results_out:
            results_out.write("{0}\t{1}\t{2}\t{3}\t{4}"
                              "\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\n".format(name_split[2], name, colony_length,
                                                                         colony_elevation, ground_elevation,
                                                                         colony_rugosity,
                                                                         overhang_prop, outcrop_prop,
                                                                         environment_rugosity,
                                                                         annotations[name][3], environment_distance))


if __name__ == '__main__':
    # Arguments
    parser = argparse.ArgumentParser(prog="Measures per colony")
    parser.add_argument('ply_filename', help="ply_filename is sample with '_env_0.1' with no .ply extension")
    parser.add_argument('annotations_filename')
    args = parser.parse_args()

    ply_filename = args.ply_filename
    annotations_filename = args.annotations_filename
    path = ".."
    main(ply_filename, annotations_filename, path)
