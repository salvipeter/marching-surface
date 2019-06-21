from math import sin, cos
import numpy as np
import marching, superd

def gyroid(p):
    return cos(p[0]) * sin(p[1]) + cos(p[1]) * sin(p[2]) + cos(p[2]) * sin(p[0])

def gyroid_gradient(p):
    return np.array([ cos(p[2]) * cos(p[0]) - sin(p[0]) * sin(p[1]),
                      cos(p[0]) * cos(p[1]) - sin(p[1]) * sin(p[2]),
                      cos(p[1]) * cos(p[2]) - sin(p[2]) * sin(p[0]) ])

def run():
    cell = marching.Cell(np.array([-3., -3, -3]), 6.1)
    cell.evaluate(gyroid, gyroid_gradient, 2)
    model = cell.surfaces()
    superd.write_model(model, 15, '/tmp/gyroid.stl')
    superd.write_boundaries(model, 50, '/tmp/gyroid-boundaries.obj')
