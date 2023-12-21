import rhinoscriptsyntax as rs
import Rhino.Geometry as rg

main_path = "/Users/koyanobunsho/Desktop/architecture/membrane_tensegrity/src/"  # Change the path to your directory
x_coord_list = []
y_coord_list = []
z_coord_list = []
points = []
with open("%soriginal_coordinates.csv" % (main_path), "r") as f:
    for line in f.readlines():
        line = line.split(",")
        x = float(line[0])
        y = float(line[1])
        z = float(line[2][:-1])
        x_coord_list.append(x)
        y_coord_list.append(y)
        z_coord_list.append(z)
        points.append(rg.Point3d(x, y, z))

edges = []
with open("%soriginal_edges.csv" % (main_path), "r") as f:
    for line in f.readlines():
        line = line.split(",")
        print(line)
        vertex_from = line[0]
        for i in range(len(line) - 2):
            point_from = points[int(vertex_from)]
            vertex_to = line[i + 1]
            point_to = points[int(vertex_to)]
            edge = rg.Line(point_from, point_to)
            edges.append(edge)
