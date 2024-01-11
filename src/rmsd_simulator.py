import random
import Rhino.Geometry as rg
import csv
import os
import subprocess


translation_vector = rg.Vector3d(constant, constant, constant)
main_path = "your_path"
vertices_num = 10

os.chdir(main_path)


def main():
    original_points = generate_random_points()
    original_edges = generate_edges_from_points(original_points)
    translated_points = generate_various_points(original_points)
    moved_points = [rotate_point(pt, axis, angle) for pt in translated_points]
    moved_edges = generate_edges_from_points(moved_points)
    to_csv(original_points, "original_points")
    to_csv(moved_points, "moved_points")
    call_py_file("rmsd_calculator")
    return original_points, original_edges, moved_points, moved_edges


def to_csv(points, file_name):
    with open("%s%s.csv" % (main_path, file_name), "w") as file:
        writer = csv.writer(file)
        x_list = []
        y_list = []
        z_list = []
        for coordinate in points:
            x_list.append(coordinate[0])
            y_list.append(coordinate[1])
            z_list.append(coordinate[2])
        writer.writerow(x_list)
        writer.writerow(y_list)
        writer.writerow(z_list)


def call_py_file(file_name):
    subprocess.call(["/usr/local/bin/docker", "exec", "mt", "./%s" % (file_name)])


def generate_random_points(point_num=vertices_num):
    points = []
    random.seed(0)
    for i in range(10):
        x = random.uniform(-10, 10)
        y = random.uniform(-10, 10)
        z = random.uniform(-10, 10)
        points.append(rg.Point3d(x, y, z))
    return points


def generate_edges_from_points(points):
    edges = []
    for i in range(len(points) - 1):
        edge = rg.Line(points[i], points[i + 1])
        edges.append(edge)
    return edges


def generate_various_points(points):
    translated_points = []
    for point in points:
        translated_point = point + translation_vector
        translated_points.append(translated_point)
    return translated_points


def rotate_point(point, axis, angle, origin=rg.Point3d(0, 0, 0)):
    transform = rg.Transform.Rotation(angle, axis, origin)
    rotated_point = point.Clone()
    rotated_point.Transform(transform)
    return rotated_point


if __name__ == "__main__":
    original_points, original_edges, moved_points, moved_edges = main()
