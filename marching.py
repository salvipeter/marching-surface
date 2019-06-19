import superd

class Cell:

    def __init__(self, origin, length):
        """The constructor requires the cell's origin point and edge length.
        The cell's corner far from its origin will be `origin + [length,length,length]`.
        """
        self.origin = origin
        self.length = length

    def evaluate(self, f, df, levels, initialized = False):
        """Initializes the cell with the given implicit function.
        `f` is the distance function, `df` is its gradient.
        `levels` is the number of times subdivision should be performed.
        The variable `initialized` is used for efficient evaluation of subcells.
        """
        values = [f(self.vertex(i)) for i in range(8)]
        gradients = [df(self.vertex(i)) for i in range(8)]
        if levels > 0:
            new_length = self.length / 2
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
            self.children = []
            for i in range(2):
                for j in range(2):
                    for k in range(2):
                        new_origin = origin + np.array([i, j, k], 'float64') * new_length
                        cell = Cell(new_origin, new_length)
                        for v in range(8):
                            r = index * 8 + v
                            if refinement[r] < 0:
                                cell.values[v] = self.values[v]
                                cell.gradients[v] = self.gradients[v]
                            else:
                                cell.values[v] = new_values[refinement[r]]
                                cell.gradients[v] = new_gradients[refinement[r]]
                        cell.evaluate(f, df, levels - 1, True)
                        self.children.append(cell)
                        index += 1

    def vertex(self, i):
        return self.origin + vertex_dirs[i] * self.length

    def plane_normal(self, i, j, k):
        for plane in range(6):
            found = 0
            face = faces[plane]
            for vertex in range(4):
                if face[vertex] == i or face[vertex] == j or face[vertex] == k:
                    if found == 2:
                        return planes[plane]
                    found += 1
        assert False, 'No plane of the cell contains these vertices'


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


# Utilities

def same_plane(e1, e2):
    vertices = set([e1[0], e1[1], e2[0], e2[1]])
    return any(all(vertex in face for vertex in vertices) for face in faces)
