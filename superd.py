import math
import numpy as np


#########################
# SuperD implementation #
#########################

# SuperD patches are simple dict objects, containing the following information:
# - n: the number of sides
# - f: a list of 3D points associated with cage faces
# - e: a list of 3D points associated with cage edges
# - v: a 3D point associated with a cage vertex

# All 3D points are represented as numpy arrays.
# Also note that e[i] is between f[i-1] and f[i].

# Example:
# >>> import superd
# >>> m = superd.read_model('trebol.sdm')
# >>> superd.write_model(m, 30, '/tmp/test.stl')

fullness = 0.5                     # Fullness is currently a global variable.


# Utilities

def normalize(v):
    """A normalized version of `v`. Does not change the original"""
    norm = np.linalg.norm(v)
    if norm == 0: 
       return v
    return v / norm

def affine_combine(p, x, q):
    """Linear interpolation between `p` and `q`, giving `p` when `x=0`"""
    return p * (1 - x) + q * x


# Parameterization

def generate_domain(n):
    v = [np.array([math.cos(a), math.sin(a)]) for a in np.linspace(0, 2 * math.pi, n + 1)][:n]
    def calc_dir(i):
        w = normalize(v[i] - v[(i-1)%n])
        dn = np.array([w[1], -w[0]])
        return -dn if np.dot(v[i], dn) > 0 else dn
    dir = [calc_dir(i) for i in range(n)]
    md = abs(np.dot(v[0] - v[1], dir[0])) # n is always <= 6 now
    return dict(v = v, dir = dir, md = md)

def parameters(domain, p):
    n = len(domain['v'])
    d = [np.dot(p - domain['v'][i], domain['dir'][i]) / domain['md'] for i in range(n)]
    def calc_s(i):
        im = (i - 1) % n
        ip = (i + 1) % n
        denom = d[im] + d[ip]
        return 0 if abs(denom) < 1.0e-8 else d[im] / denom
    s = [calc_s(i) for i in range(n)]
    return (s, d)


# Ribbon generation

def generate_quartic(points):
    x1 = (2/5 * fullness + 3/5) * fullness
    x2 = (-2/7 * fullness + 9/7) * fullness
    return [points[0],
            affine_combine(points[0], x1, points[1]),
            affine_combine(affine_combine(points[0], x2, points[1]),
                           1/2,
                           affine_combine(points[2], x2, points[1])),
            affine_combine(points[2], x1, points[1]),
            points[2]]

def generate_base(patch, i):
    im = (i - 1) % patch['n']
    return generate_quartic([patch['f'][im], patch['e'][i], patch['f'][i]])

def generate_mid(patch, i):
    im = (i - 1) % patch['n']
    ip = (i + 1) % patch['n']
    return generate_quartic([patch['e'][im], patch['v'], patch['e'][ip]])

def generate_opp(patch, i):
    n = patch['n']
    if n == 3:
        return [patch['f'][(i+1)%n] for _ in range(5)]
    if n == 4:
        return generate_base(patch, (i + 2) % n)[::-1]
    assert False, "generate_opp() should only be called with 3 or 4 sides"

def generate_ribbon(patch, i):
    result = []
    base = generate_base(patch, i)
    if patch['n'] > 4:
        left = generate_base(patch, (i - 1) % patch['n'])
        right = generate_base(patch, (i + 1) % patch['n'])
        for j in range(5):
            row = []
            for k in range(5):
                row.append(base[j] + (4-j)/4 * (left[4-k] - left[4]) + j/4 * (right[k] - right[0]))
            result.append(row)
        return result
    mid = generate_mid(patch, i)
    opp = generate_opp(patch, i)
    for j in range(5):
        result.append(generate_quartic([base[j], mid[j], opp[j]]))
    return result


# Patch definition

def blend(d, i):
    return np.prod([d[j] ** 2 for j in set(range(len(d))).difference([i])])

def bernstein(n, u):
    coeff = [1]
    for j in range(n):
        saved = 0
        for k in range(j + 1):
            tmp = coeff[k]
            coeff[k] = saved + tmp * (1 - u)
            saved = tmp * u
        coeff.append(saved)
    return coeff

def eval_bezier(cp, u, v):
    order = [len(cp), len(cp[0])]
    coeff_u = bernstein(order[0] - 1, u)
    coeff_v = bernstein(order[1] - 1, v)
    result = np.array([0., 0, 0])
    for i in range(order[0]):
        for j in range(order[1]):
            result += cp[i][j] * coeff_u[i] * coeff_v[j]
    return result

def eval_patch_impl(domain, ribbons, p):
    result = np.array([0., 0, 0])
    n = len(ribbons)
    s, d = parameters(domain, p)
    tolerance = 1.0e-5
    if sum([1 for x in d if x < tolerance]) == 2:
        for i in range(n):
            if d[i] < tolerance and d[(i-1)%n] < tolerance:
                return ribbons[i][0][0]
    blends = [blend(d, i) for i in range(n)]
    blendsum = sum(blends)
    for i in range(n):
        result += eval_bezier(ribbons[i], s[i], d[i]) * blends[i] / blendsum
    return result

def eval_patch(patch, p):
    domain = generate_domain(patch['n'])
    ribbons = [generate_ribbon(patch, i) for i in range(patch['n'])]
    return eval_patch_impl(domain, ribbons, p)


# I/O

def vertices(poly, resolution):
    result = []
    n = len(poly)
    if n == 3:
        for j in range(resolution + 1):
            u = j / resolution
            p = poly[0] * u + poly[2] * (1 - u)
            q = poly[1] * u + poly[2] * (1 - u)
            for k in range(j + 1):
                v = 1 if j == 0 else k / j
                result.append(p * (1 - v) + q * v)
        return result
    if n == 4:
        for j in range(resolution + 1):
            u = j / resolution
            p = poly[0] * (1 - u) + poly[1] * u
            q = poly[3] * (1 - u) + poly[2] * u
            for k in range(resolution + 1):
                v = k / resolution
                result.append(p * (1 - v) + q * v)
        return result
    # n > 4
    lines = [(poly[(i-1)%n], poly[i]) for i in range(n)]
    center = np.array([0, 0])
    result.append(center)
    for j in range(1, resolution + 1):
        coeff = j / resolution
        for k in range(n):
            for i in range(j):
                lp = affine_combine(lines[k][0], i / j, lines[k][1])
                result.append(affine_combine(center, coeff, lp))
    return result

def triangles(n, resolution):
    result = []
    if n == 3:
        prev = 0
        current = 1
        for i in range(resolution):
            for j in range(i):
                result.append([current + j, current + j + 1, prev + j])
                result.append([current + j + 1, prev + j + 1, prev + j])
            result.append([current + i, current + i + 1, prev + i])
            prev = current
            current += i + 2
        return result
    if n == 4:
        for i in range(resolution):
            for j in range(resolution):
                index = i * (resolution + 1) + j
                result.append([index, index + resolution + 1, index + 1])
                result.append([index + 1, index + resolution + 1, index + resolution + 2])
        return result
    # n > 4
    inner_start = 0
    outer_vert = 1
    for layer in range(1, resolution + 1):
        inner_vert = inner_start
        outer_start = outer_vert
        for side in range(1, n + 1):
            vert = 1
            while True:
                next_vert = outer_start if side == n and vert == layer else outer_vert + 1
                result.append([inner_vert, outer_vert, next_vert])
                outer_vert += 1
                vert += 1
                if vert == layer + 1:
                    break
                inner_next = inner_start if side == n and vert == layer else inner_vert + 1
                result.append([inner_vert, next_vert, inner_next])
                inner_vert = inner_next
        inner_start = outer_start
    return result

def mesh_size(n, resolution):
    if n == 3:
        return (resolution + 1) * (resolution + 2) // 2
    if n == 4:
        return (resolution + 1) ** 2
    return 1 + n * resolution * (resolution + 1) // 2

def write_obj(verts, tris, filename):
    with open(filename, 'w') as f:
        for v in verts:
            f.write("v {0} {1} {2}\n".format(v[0], v[1], v[2]))
        for t in tris:
            f.write("f {0} {1} {2}\n".format(t[0] + 1, t[1] + 1, t[2] + 1))

def write_patch(patch, resolution, filename):
    """Writes a single `patch` into an OBJ file."""
    domain = generate_domain(patch['n'])
    ribbons = [generate_ribbon(patch, i) for i in range(patch['n'])]
    points = [eval_patch_impl(domain, ribbons, p) for p in vertices(domain['v'], resolution)]
    write_obj(points, triangles(patch['n'], resolution), filename)

def read_model(filename):
    """Reads a SuperD Model (SDM) file."""
    result = []
    def read_point(f):
        return np.array(list(map(float, f.readline().split())))
    with open(filename) as f:
        n_patches = int(f.readline())
        for _ in range(n_patches):
            n = int(f.readline())
            pf = [read_point(f) for _ in range(n)]
            pe = [read_point(f) for _ in range(n)]
            pv = read_point(f)
            result.append(dict(n = n, f = pf, e = pe, v = pv))
    return result

def write_ribbon(patch, i, filename):
    """Writes the `i`-th ribbon of `patch` into an OBJ file."""
    r = generate_ribbon(patch, i)
    with open(filename, "w") as f:
        for i in range(5):
            for j in range(5):
                p = r[i][j]
                f.write("v {0} {1} {2}\n".format(p[0], p[1], p[2]))
        for i in range(5):
            for j in range(4):
                index = i * 5 + j
                f.write("l {0} {1}\n".format(index + 1, index + 2))
                index = j * 5 + i
                f.write("l {0} {1}\n".format(index + 1, index + 6))

def write_model(model, resolution, filename):
    """Writes `model` (a list of patches) into an STL file."""
    with open(filename, "wb") as f:
        f.write(b'Generated by SuperD.' + b' ' * 60)
        n_faces = sum([mesh_size(patch['n'], resolution) for patch in model])
        f.write(int(n_faces).to_bytes(4, byteorder='big', signed=False))
        for patch in model:
            domain = generate_domain(patch['n'])
            ribbons = [generate_ribbon(patch, i) for i in range(patch['n'])]
            points = [eval_patch_impl(domain, ribbons, p)
                      for p in vertices(domain['v'], resolution)]
            for tri in triangles(patch['n'], resolution):
                n = normalize(np.cross(points[tri[1]] - points[tri[0]],
                                       points[tri[2]] - points[tri[0]]))
                np.array(n, 'float32').tofile(f)
                for i in range(3):
                    np.array(points[tri[i]], 'float32').tofile(f)
                f.write(int(0).to_bytes(2, byteorder='big', signed=False))
