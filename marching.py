import numpy as np
import superd

#####################
# Marching Surfaces #
#####################

# Example:
# >>> import marching, numpy, superd
# >>> c = marching.Cell(numpy.array([-1.6, -1.6, -1.6]), 3)
# >>> c.evaluate(lambda p: numpy.linalg.norm(p) - 1,   # Unit sphere distance function
# ...            lambda p: p / numpy.linalg.norm(p),   # Unit sphere distance gradient
# ...            2)                                    # Number of subdivisons (-> 8^2 = 64 cells)
# >>> m = c.surfaces()
# >>> superd.write_model(m, 15, '/tmp/sphere.stl')
# >>> superd.write_boundaries(m, 50, '/tmp/sphere-boundaries.obj')

class Cell:

    def __init__(self, origin, length):
        """The constructor requires the cell's origin point and edge length.
        The cell's corner far from its origin will be `origin + [length,length,length]`.
        """
        self.origin = origin
        self.length = length
        self.values = []
        self.gradients = []

    def evaluate(self, f, df, levels, initialized = False):
        """Initializes the cell with the given implicit function.
        `f` is the distance function, `df` is its gradient.
        `levels` is the number of times subdivision should be performed.
        The variable `initialized` is used for efficient evaluation of subcells.
        """
        self.values = [f(self.vertex(i)) for i in range(8)]
        self.gradients = [df(self.vertex(i)) for i in range(8)]
        self.children = []
        if levels > 0:
            new_length = self.length / 2
            new_values = []
            new_gradients = []
            def add(i, j, k):
                p = self.origin + np.array([i, j, k], 'float64') * new_length
                new_values.append(f(p))
                new_gradients.append(df(p))
            # Edge midpoints
            add(1, 0, 0)
            add(2, 1, 0) #
            add(1, 2, 0) #    +---6----+
            add(0, 1, 0) #   11       9|
            add(1, 0, 2) #  / 7      / 5
            add(2, 1, 2) # +---2----+  |
            add(1, 2, 2) # |  +---4-|--+
            add(0, 1, 2) # 3 10     1 8
            add(2, 0, 1) # |/       |/
            add(2, 2, 1) # +---0----+
            add(0, 0, 1) #
            add(0, 2, 1)
            # Face midpoints
            add(1, 1, 0) # front
            add(1, 1, 2) # back
            add(2, 1, 1) # right
            add(0, 1, 1) # left
            add(1, 0, 1) # bottom
            add(1, 2, 1) # top
            # Center point
            add(1, 1, 1)
            index = 0
            for i in range(2):
                for j in range(2):
                    for k in range(2):
                        new_origin = self.origin + np.array([i, j, k], 'float64') * new_length
                        cell = Cell(new_origin, new_length)
                        for v in range(8):
                            r = index * 8 + v
                            if refinement[r] < 0:
                                cell.values.append(self.values[v])
                                cell.gradients.append(self.gradients[v])
                            else:
                                cell.values.append(new_values[refinement[r]])
                                cell.gradients.append(new_gradients[refinement[r]])
                        cell.evaluate(f, df, levels - 1, True)
                        self.children.append(cell)
                        index += 1

    def vertex(self, i):
        """Returns the position of the `i`-th vertex.
        For the numbering scheme, see the diagram below:
              7         6
             +--------+
            /|       /|
         3 / |    2 / |
          +--------+  |          y
          |  +-----|--+          ^
          | / 4    | / 5         | z
          |/       |/            |/
        0 +--------+ 1           +---> x
        """
        return self.origin + vertex_dirs[i] * self.length

    def surfaces(self):
        """Generates surfaces approximating the implicit function given with `evaluate`."""
        if self.children:
            return sum([child.surfaces() for child in self.children], [])

        # Find zero crossings
        crosses = []
        for i in range(12):
            if self.values[edges[i][0]] * self.values[edges[i][1]] < 0:
                crosses.append(i)
        n_crosses = len(crosses)

        if n_crosses == 0:
            return []

        # Sort the crossings and compute the corresponding corner points and normals
        sorted_crosses = []
        corners = []
        normals = []
        cross = 0
        last_cross = -1
        while True:
            i1 = edges[crosses[cross]][0]
            i2 = edges[crosses[cross]][1]
            v1 = self.values[i1]
            v2 = self.values[i2]
            alpha = abs(v1) / abs(v2 - v1)
            sorted_crosses.append(crosses[cross])
            corners.append(self.vertex(i1) * (1 - alpha) + self.vertex(i2) * alpha)
            normals.append(self.gradients[i1] * (1 - alpha) + self.gradients[i2] * alpha)
            for j in range(n_crosses):
                if j != last_cross and j != cross and \
                   same_plane(edges[crosses[cross]], edges[crosses[j]]):
                    last_cross = cross
                    cross = j
                    break
            if cross == 0:
                break
        sides = len(corners)

        if sides != n_crosses:
            return []                     # TODO - what should we do here?

        # Reverse the loop if needed, such that positive is outside
        first_edge = edges[sorted_crosses[0]]
        second_edge = edges[sorted_crosses[1]]
        left_face = faces[left_faces[sorted_crosses[0]]]
        negative = all(vertex in left_face for vertex in second_edge)
        if (negative and self.values[first_edge[0]] > 0) or \
            (not negative and self.values[first_edge[0]] < 0):
             sorted_crosses.reverse()
             corners.reverse()
             normals.reverse()

        return [ self.generate_superd(sorted_crosses, corners, normals) ]

    def generate_superd(self, crosses, points, normals):
        """Generates a SuperD patch, given a sorted list of crossed edges,
        and the corresponding points and normal vectors.
        """
        sides = len(points)

        # Generate the mid-edge control points
        pe = []
        for i in range(sides):
            ip = (i + 1) % sides
            # Compute normal vectors
            vertices = list(set([edges[crosses[i]][0], edges[crosses[ip]][0],
                                 edges[crosses[i]][1], edges[crosses[ip]][1]]))
            p1 = points[i]
            p2 = points[ip]
            n1 = normalize(normals[i])
            n2 = normalize(normals[ip])
            pn = plane_normal(vertices[0], vertices[1], vertices[2])

            # Compute edge control point positions
            A = np.array([n1, n2, pn])
            b = np.array([np.dot(p1, n1), np.dot(p2, n2), np.dot((p1 + p2) / 2, pn)])
            p = np.linalg.solve(A, b)

            # Move the computed point to a preferable position
            d1 = np.array([1., 0, 0])
            d2 = np.array([0., 1, 0])
            if abs(np.dot(d1, pn)) > 0.5:
                d1 = np.array([0., 0, 1])
            elif abs(np.dot(d2, pn)) > 0.5:
                d2 = np.array([0., 0, 1])
            base = self.origin if min(pn) < 0 else self.origin + pn * self.length
            t1 = np.dot(p - base, d1)
            t2 = np.dot(p - base, d2)
            t1 = min(max(t1, 0), self.length)
            t2 = min(max(t2, 0), self.length)
            p = base + d1 * t1 + d2 * t2
            pe.append(p)

        pf = [points[(i+1)%sides] for i in range(sides)]
        pv = sum(pe) / sides
        return dict(n = sides, f = pf, e = pe, v = pv)


# Global variables

vertex_dirs = [ np.array([0., 0, 0]), np.array([1., 0, 0]),
                np.array([1., 1, 0]), np.array([0., 1, 0]),
                np.array([0., 0, 1]), np.array([1., 0, 1]),
                np.array([1., 1, 1]), np.array([0., 1, 1]) ]

edges = [ (0, 1), (1, 2), (2, 3), (3, 0),
          (4, 5), (5, 6), (6, 7), (7, 4),
          (1, 5), (6, 2), (0, 4), (7, 3) ]

faces = [ [0, 1, 2, 3], [4, 5, 6, 7], [1, 5, 6, 2],
          [0, 4, 7, 3], [0, 1, 5, 4], [3, 2, 6, 7] ]

planes = [ np.array([0.,  0, -1]), np.array([0.,  0, 1]), np.array([1., 0, 0]),
           np.array([-1., 0,  0]), np.array([0., -1, 0]), np.array([0., 1, 0]) ]

refinement = [ -1,  0, 12,  3, 10, 16, 18, 15,
               10, 16, 18, 15, -1,  4, 13,  7,
                3, 12,  2, -1, 15, 18, 17, 11,
               15, 18, 17, 11,  7, 13,  6, -1,
                0, -1,  1, 12, 16,  8, 14, 18,
               16,  8, 14, 18,  4, -1,  5, 13,
               12,  1, -1,  2, 18, 14,  9, 17,
               18, 14,  9, 17, 13,  5, -1,  6 ]

left_faces = [ 0, 0, 0, 0, 4, 2, 5, 3, 2, 2, 4, 5 ]


# Utilities

def same_plane(e1, e2):
    """Returns whether two edges belong to the same plane."""
    vertices = set([e1[0], e1[1], e2[0], e2[1]])
    return any(all(vertex in face for vertex in vertices) for face in faces)

def plane_normal(i, j, k):
    """Returns the normal vector of the plane containing vertices `i`, `j` and `k`."""
    for plane in range(6):
        found = 0
        face = faces[plane]
        for vertex in range(4):
            if face[vertex] == i or face[vertex] == j or face[vertex] == k:
                if found == 2:
                    return planes[plane]
                found += 1
    assert False, 'No plane of the cell contains these vertices'

def normalize(v):
    """A normalized version of `v`. Does not change the original."""
    norm = np.linalg.norm(v)
    if norm == 0: 
       return v
    return v / norm

def run_example():
    c = Cell(np.array([-1.6, -1.6, -1.6]), 3)
    c.evaluate(lambda p: np.linalg.norm(p) - 1, lambda p: normalize(p), 2)
    m = c.surfaces()
    superd.write_model(m, 15, '/tmp/sphere.stl')
    superd.write_boundaries(m, 50, '/tmp/sphere-boundaries.obj')
