#include <algorithm>
#include <array>
#include <cassert>
#include <fstream>
#include <iostream>
#include <set>
#include <sstream>

#include <Eigen/LU>

#include <domain.hh>
#include <ribbon.hh>
#include <surface-generalized-bezier.hh>
#include <surface-spatch.hh>
#include <surface-superd.hh>

#include "gb-io.hh"
#include "superd-io.hh"

using namespace Geometry;
using Transfinite::Surface;
using Transfinite::SurfaceGeneralizedBezier;
using Transfinite::SurfaceSPatch;
using Transfinite::SurfaceSuperD;

enum class SurfaceType { GENERALIZED_BEZIER, SUPERD, SPATCH };

// Convenience function for testing membership
template<typename Container, typename T>
bool contains(const Container &c, const T &value) {
  return std::find(c.begin(), c.end(), value) != c.end();
}

// The main class, represents one cell (but it can generate subcells)
class Cell {
public:
  // The constructor requires the cell's origin point and edge length,
  // so the cell's corner far from its origin will be `origin + Vector3D(length,length,length)`
  Cell(const Point3D &origin, double length);

  // Initializes the cell with the given implicit function.
  //  `F` should be of the type `Point3D -> double`,   [the distance]
  // `dF` should be of the type `Point3D -> Vector3D`. [the gradient]
  // `min_levels` is the minimum number of times subdivision should be performed,
  // `max_levels` is the maximum number of times subdivision should be performed.
  // The variable `initialized` is used for efficient evaluation of subcells.
  // Subdivisions between `min_levels` and `max_levels` are performed when the cell
  // is not splittable into two (i.e., it would define multiple surfaces).
  template<typename F, typename DF>
  void init(std::pair<F,DF> fdf, size_t min_levels, size_t max_levels, bool initialized = false);

  // Returns the value at the `i`-th vertex (see `vertex`)
  double value(int i) const { return values[i]; }

  // Returns the gradient at the `i`-th vertex (see `vertex`)
  const Vector3D &gradient(int i) const { return gradients[i]; }

  // Returns the position of the `i`-th vertex.
  // For the numbering scheme, see the diagram below:
  //       7         6
  //      +--------+
  //     /|       /|
  //  3 / |    2 / |
  //   +--------+  |          y
  //   |  +-----|--+          ^
  //   | / 4    | / 5         | z
  //   |/       |/            |/
  // 0 +--------+ 1           +---> x
  Point3D vertex(int i) const;

  // Returns the normal vector of the plane containing vertices `i`, `j` and `k` (see `vertex`)
  const Vector3D &planeNormal(int i, int j, int k) const;

  // Generates surfaces of the given type, approximating the implicit function given with `init`.
  std::vector<std::unique_ptr<Surface>> surfaces(SurfaceType type) const;

private:
  // Edges and faces are represented by their vertices (see `vertex`)
  using Edge = std::pair<int,int>;
  using Face = std::array<int,4>;

  // Returns whether two edges belong to the same plane
  static bool samePlane(const Edge &e1, const Edge &e2);

  // Generates a GB patch, given a sorted list of crossed edges,
  // and the corresponding points and normal vectors.
  std::unique_ptr<SurfaceGeneralizedBezier> generateGB(const std::vector<int> &crosses,
                                                       const PointVector &points,
                                                       const VectorVector &normals) const;

  // Generates a SuperD patch, given a sorted list of crossed edges,
  // and the corresponding points and normal vectors.
  std::unique_ptr<SurfaceSuperD> generateSuperD(const std::vector<int> &crosses,
                                                const PointVector &points,
                                                const VectorVector &normals) const;

  // Generates an S-patch, given a sorted list of crossed edges,
  // and the corresponding points and normal vectors.
  std::unique_ptr<SurfaceSPatch> generateSPatch(const std::vector<int> &crosses,
                                                const PointVector &points,
                                                const VectorVector &normals) const;

  Point3D origin;
  double length;
  std::array<double, 8> values;
  std::array<Vector3D, 8> gradients;
  std::array<std::unique_ptr<Cell>,8> children;

  static const std::array<Edge,12> edges;
  static const std::array<Face,6> faces;
  static const std::array<Vector3D,6> planes;
  static const std::array<int,64> refinement;
  static const std::array<unsigned char,122> good_configs;
};

Cell::Cell(const Point3D &origin, double length)
  : origin(origin), length(length)
{
}

const std::array<Cell::Edge,12> Cell::edges = {
  {
    {0, 1}, {1, 2}, {2, 3}, {3, 0},
    {4, 5}, {5, 6}, {6, 7}, {7, 4},
    {1, 5}, {6, 2}, {0, 4}, {7, 3}
  }
};

const std::array<Cell::Face,6> Cell::faces = {
  {
    {0, 1, 2, 3}, {4, 5, 6, 7}, {1, 5, 6, 2},
    {0, 4, 7, 3}, {0, 1, 5, 4}, {3, 2, 6, 7}
  }
};

const std::array<Vector3D,6> Cell::planes = {
  {
    { 0, 0, -1}, {0,  0, 1}, {1, 0, 0},
    {-1, 0,  0}, {0, -1, 0}, {0, 1, 0}
  }
};

const std::array<int,64> Cell::refinement = {
  -1, 0, 12, 3, 10, 16, 18, 15,
  10, 16, 18, 15, -1, 4, 13, 7,
  3, 12, 2, -1, 15, 18, 17, 11,
  15, 18, 17, 11, 7, 13, 6, -1,
  0, -1, 1, 12, 16, 8, 14, 18,
  16, 8, 14, 18, 4, -1, 5, 13,
  12, 1, -1, 2, 18, 14, 9, 17,
  18, 14, 9, 17, 13, 5, -1, 6
};

const std::array<unsigned char,122> Cell::good_configs = { // 122/256 configurations
  // 0 vertex (2) - no surface
  0b00000000,
  0b11111111,
  // 1 vertex (16) - 3-sided
  0b10000000, 0b01000000, 0b00100000, 0b00010000, 0b00001000, 0b00000100, 0b00000010, 0b00000001,
  0b01111111, 0b10111111, 0b11011111, 0b11101111, 0b11110111, 0b11111011, 0b11111101, 0b11111110,
  // 2 vertices (24) - 4-sided
  0b11000000, 0b10010000, 0b10001000, 0b01100000, 0b01000100, 0b00110000, 0b00100010, 0b00010001,
  0b00111111, 0b01101111, 0b01110111, 0b10011111, 0b10111011, 0b11001111, 0b11011101, 0b11101110,
  0b00001100, 0b00001001, 0b00000110, 0b00000011,
  0b11110011, 0b11110110, 0b11111001, 0b11111100,
  // 3 vertices (48) - 5-sided
  0b11100000, 0b11010000, 0b11001000, 0b11000100, 0b10110000, 0b10011000, 0b10010001, 0b10001100,
  0b00011111, 0b00101111, 0b00110111, 0b00111011, 0b01001111, 0b01100111, 0b01101110, 0b01110011,
  0b10001001, 0b01110000, 0b01100100, 0b01100010, 0b01001100, 0b01000110, 0b00110010, 0b00110001,
  0b01110110, 0b10001111, 0b10011011, 0b10011101, 0b10110011, 0b10111001, 0b11001101, 0b11001110,
  0b00100110, 0b00100011, 0b00011001, 0b00010011, 0b00001110, 0b00001101, 0b00001011, 0b00000111,
  0b11011001, 0b11011100, 0b11100110, 0b11101100, 0b11110001, 0b11110010, 0b11110100, 0b11111000,
  // 4 vertices - A (6) - 4-sided
  0b11110000, 0b11001100, 0b10011001, 0b01100110, 0b00110011, 0b00001111,
  // 4 vertices - B (2) - 6-sided
  0b10100101, 0b01011010,
  // 4 vertices - C (24) - 6-sided
  0b11101000, 0b11100010, 0b11010100, 0b11010001, 0b11001001, 0b11000110, 0b10111000, 0b10110010,
  0b10011100, 0b10010011, 0b10001110, 0b10001011, 0b01110100, 0b01110001, 0b01101100, 0b01100011,
  0b01001101, 0b01000111, 0b00111001, 0b00110110, 0b00101110, 0b00101011, 0b00011101, 0b00010111
};

Point3D
Cell::vertex(int i) const {
  switch (i) {
  case 0: return origin;
  case 1: return origin + Vector3D(length, 0,      0);
  case 2: return origin + Vector3D(length, length, 0);
  case 3: return origin + Vector3D(0     , length, 0);
  case 4: return origin + Vector3D(0     , 0     , length);
  case 5: return origin + Vector3D(length, 0,      length);
  case 6: return origin + Vector3D(length, length, length);
  case 7: return origin + Vector3D(0     , length, length);
  default:;
  }
  assert(false && "Invalid vertex");
}

const Vector3D &
Cell::planeNormal(int i, int j, int k) const {
  for (int plane = 0; plane < 6; ++plane) {
    int found = 0;
    const auto &face = faces[plane];
    for (int vertex = 0; vertex < 4; ++vertex)
      if (face[vertex] == i || face[vertex] == j || face[vertex] == k) {
        if (found == 2)
          return planes[plane];
        ++found;
      }
  }
  assert(false && "No plane of the cell contains these vertices");
}

bool
Cell::samePlane(const Edge &e1, const Edge &e2) {
  std::set<int> vertices;
  vertices.insert(e1.first);
  vertices.insert(e1.second);
  vertices.insert(e2.first);
  vertices.insert(e2.second);
  return std::any_of(faces.begin(), faces.end(), [&vertices](const Face &face) {
      return std::all_of(vertices.begin(), vertices.end(), [&face](int vertex) {
          return contains(face, vertex);
        });
    });
}

template<typename F, typename DF>
void
Cell::init(std::pair<F,DF> fdf, size_t min_levels, size_t max_levels, bool initialized) {
  unsigned char id = 0, bit = 1;
  for (int i = 0; i < 8; ++i, bit <<= 1) {
    if (!initialized) {
      values[i] = fdf.first(vertex(i));
      gradients[i] = fdf.second(vertex(i));
    }
    if (values[i] < 0)
      id += bit;
  }
  if (min_levels > 0 || (max_levels > 0 && !contains(good_configs, id))) {
    double new_length = length / 2;

    // Compute new points to be reused by the children
    std::array<double,19> new_values;
    std::array<Vector3D,19> new_gradients;
    auto add = [&](int index, int i, int j, int k) {
                 auto p = origin + Vector3D(i, j, k) * new_length;
                 new_values[index] = fdf.first(p); new_gradients[index] = fdf.second(p);
               };
    // Edge midpoints
    add( 0, 1, 0, 0);
    add( 1, 2, 1, 0); //
    add( 2, 1, 2, 0); //    +---6----+
    add( 3, 0, 1, 0); //   11       9|
    add( 4, 1, 0, 2); //  / 7      / 5
    add( 5, 2, 1, 2); // +---2----+  |
    add( 6, 1, 2, 2); // |  +---4-|--+
    add( 7, 0, 1, 2); // 3 10     1 8
    add( 8, 2, 0, 1); // |/       |/
    add( 9, 2, 2, 1); // +---0----+
    add(10, 0, 0, 1); //
    add(11, 0, 2, 1);
    // Face midpoints
    add(12, 1, 1, 0); // front
    add(13, 1, 1, 2); // back
    add(14, 2, 1, 1); // right
    add(15, 0, 1, 1); // left
    add(16, 1, 0, 1); // bottom
    add(17, 1, 2, 1); // top
    // Center point
    add(18, 1, 1, 1);

    int index = 0;
    for (int i = 0; i <= 1; ++i)
      for (int j = 0; j <= 1; ++j)
        for (int k = 0; k <= 1; ++k) {
          auto new_origin = origin + Vector3D(i, j, k) * new_length;
          auto cell = std::make_unique<Cell>(new_origin, new_length);
          for (int v = 0; v < 8; ++v) {
            int r = index * 8 + v;
            if (refinement[r] < 0) { // same as parent's
              cell->values[v] = values[v];
              cell->gradients[v] = gradients[v];
            } else {
              cell->values[v] = new_values[refinement[r]];
              cell->gradients[v] = new_gradients[refinement[r]];
            }
          }
          cell->init(fdf, min_levels ? min_levels - 1 : 0, max_levels - 1, true);
          children[index++] = std::move(cell);
        }
  }
}

std::unique_ptr<SurfaceGeneralizedBezier>
Cell::generateGB(const std::vector<int> &crosses,
                 const PointVector &points,
                 const VectorVector &normals) const {
  int sides = points.size();
  auto surf = std::make_unique<SurfaceGeneralizedBezier>();
  surf->initNetwork(sides, 3);

  // Corner control points are already computed
  for (int i = 0; i < sides; ++i)
    surf->setControlPoint(i, 0, 0, points[i]);

  // Tangent control points
  for (int i = 0; i < sides; ++i) {
    // Find the 3 or 4 vertices defining this curve,
    // and mark the central vertex in the 3-vertex case
    int ip = (i + 1) % sides;
    std::array<int,4> vertices = { 
      edges[crosses[i]].first,
      edges[crosses[i]].second,
      edges[crosses[ip]].first,
      edges[crosses[ip]].second
    };
    std::array<bool,8> seen;
    seen.fill(false);
    int central_vertex = -1;
    for (int j = 0; j < 4; ++j) {
      if (seen[vertices[j]])
        central_vertex = vertices[j];
      else
        seen[vertices[j]] = true;
    }
    int k = 0;
    for (int j = 0; j < 12; ++j)
      if (seen[j])
        vertices[k++] = j;

    // Now we can compute the tangent control points
    auto p1 = points[i], p2 = points[ip];
    auto n1 = normals[i], n2 = normals[ip];
    auto pn = planeNormal(vertices[0], vertices[1], vertices[2]);
    auto t1 = n1.normalize() ^ pn, t2 = n2.normalize() ^ pn;

    // Reverse the sign when appropriate:
    if (central_vertex != -1) {
      auto center = vertex(central_vertex);
      if (t1 * (p2 - center) < 0)
        t1 *= -1;
      if (t2 * (p1 - center) < 0)
        t2 *= -1;
    } else {
      if (t1 * (p2 - p1) < 0)
        t1 *= -1;
      if (t2 * (p1 - p2) < 0)
        t2 *= -1;
    }

    // Scale the tangents and set the control points
    double scaling = (p1 - p2).norm() / 3.0; // kutykurutty
    surf->setControlPoint(i, 1, 0, p1 + t1 * scaling);
    surf->setControlPoint(i, 2, 0, p2 + t2 * scaling);
  }

  // Twist control points by the parallelogram rule
  Point3D p;
  for (int i = 0; i < sides; ++i) {
    p = surf->controlPoint(i, 0, 1) + surf->controlPoint(i, 1, 0) - surf->controlPoint(i, 0, 0);
    surf->setControlPoint(i, 1, 1, p);
  }

  // Central control point is the mass center of the twist control points
  p = Point3D(0, 0, 0);
  for (int i = 0; i < sides; ++i)
    p += surf->controlPoint(i, 1, 1);
  surf->setCentralControlPoint(p / sides);

  surf->setupLoop();
  return surf;
}

std::unique_ptr<SurfaceSuperD>
Cell::generateSuperD(const std::vector<int> &crosses,
                     const PointVector &points,
                     const VectorVector &normals) const {
  int sides = points.size();
  auto surf = std::make_unique<SurfaceSuperD>();
  surf->initNetwork(sides);
  Point3D cp(0, 0, 0);
  for (int i = 0; i < sides; ++i) {
    int ip = (i + 1) % sides;
    surf->setFaceControlPoint(i, points[ip]);

    // Compute normal vectors
    std::array<int, 4> vertices = { 
      edges[crosses[i]].first,
      edges[crosses[i]].second,
      edges[crosses[ip]].first,
      edges[crosses[ip]].second
    };
    std::sort(vertices.begin(), vertices.end());
    std::unique(vertices.begin(), vertices.end());
    auto p1 = points[i], p2 = points[ip];
    auto n1 = normals[i]; n1.normalize();
    auto n2 = normals[ip]; n2.normalize();
    auto pn = planeNormal(vertices[0], vertices[1], vertices[2]);

    // Compute edge control point positions
    Eigen::Matrix3d A;
    Eigen::Vector3d b;
    A << n1[0], n1[1], n1[2], n2[0], n2[1], n2[2], pn[0], pn[1], pn[2];
    b << p1 * n1, p2 * n2, (p1 + p2) / 2 * pn;
    Eigen::Vector3d x = A.fullPivLu().solve(b);
    Point3D p(x(0), x(1), x(2));

    // Move the computed point to a preferable position
    Vector3D d1(1, 0, 0), d2(0, 1, 0);
    if (std::abs(d1 * pn) > 0.5)
      d1 = Vector3D(0, 0, 1);
    else if (std::abs(d2 * pn) > 0.5)
      d2 = Vector3D(0, 0, 1);
    auto base = (pn[0] < 0 || pn[1] < 0 || pn[2] < 0) ? origin : origin + pn * length;
    double t1 = (p - base) * d1, t2 = (p - base) * d2;
    t1 = std::min(std::max(t1, 0.0), length);
    t2 = std::min(std::max(t2, 0.0), length);
    p = base + d1 * t1 + d2 * t2;

    surf->setEdgeControlPoint(i, p);
    cp += surf->edgeControlPoint(i);
  }
  surf->setVertexControlPoint(cp / sides); // kutykurutty
  surf->setupLoop();
  surf->updateRibbons();
  return surf;
}

std::unique_ptr<SurfaceSPatch>
Cell::generateSPatch(const std::vector<int> &crosses,
                     const PointVector &points,
                     const VectorVector &normals) const {
  int sides = points.size(), depth = 1; // TODO: should be 3
  auto surf = std::make_unique<SurfaceSPatch>();
  surf->initNetwork(sides, depth);
  Point3D cp(0, 0, 0);
  SurfaceSPatch::Index index(sides);
  for (int i = 0; i < sides; ++i) {
    std::fill(index.begin(), index.end(), 0);
    index[i] = depth;
    surf->setControlPoint(index, points[i]);
  }
  surf->setupLoop();
  return surf;
}

std::vector<std::unique_ptr<Surface>>
Cell::surfaces(SurfaceType type) const {
  std::vector<std::unique_ptr<Surface>> result;

  if (children[0]) {
    for (int i = 0; i < 8; ++i) {
      auto surfaces = children[i]->surfaces(type);
      std::move(surfaces.begin(), surfaces.end(), std::back_inserter(result));
    }
    return result;
  }

  std::vector<int> crosses;
  for (int i = 0; i < 12; ++i)
    if (value(edges[i].first) * value(edges[i].second) < 0)
      crosses.push_back(i);
  int n_crosses = crosses.size();

  if (n_crosses == 0)           // trivial case
    return { };

  std::vector<int> sorted_crosses;
  PointVector corners;
  VectorVector normals;
  int cross = 0, last_cross = -1;
  do {
    int i1 = edges[crosses[cross]].first, i2 = edges[crosses[cross]].second;
    double v1 = value(i1), v2 = value(i2);
    double alpha = std::abs(v1) / std::abs(v2 - v1);
    sorted_crosses.push_back(crosses[cross]);
    corners.push_back(vertex(i1) * (1 - alpha) + vertex(i2) * alpha);
    normals.push_back(gradient(i1) * (1 - alpha) + gradient(i2) * alpha);
    for (int j = 0; j < n_crosses; ++j)
      if (j != last_cross && j != cross && samePlane(edges[crosses[cross]], edges[crosses[j]])) {
        last_cross = cross;
        cross = j;
        break;
      }
  } while (cross != 0);
  int sides = corners.size();

  if (sides != n_crosses)       // ambiguous case
    return { };                 // kutykurutty

  // Reverse the loop if needed, such that positive is outside
  static const std::array<int,12> left_faces = { 0, 0, 0, 0, 4, 2, 5, 3, 2, 2, 4, 5 };
  auto const &first_edge = edges[sorted_crosses[0]];
  auto const &second_edge = edges[sorted_crosses[1]];
  auto const &left_face = faces[left_faces[sorted_crosses[0]]];
  bool negative = contains(left_face, second_edge.first) && contains(left_face, second_edge.second);
  if ((negative && values[first_edge.first] > 0) ||
      (!negative && values[first_edge.first] < 0)) {
    std::reverse(sorted_crosses.begin(), sorted_crosses.end());
    std::reverse(corners.begin(), corners.end());
    std::reverse(normals.begin(), normals.end());
  }

  std::unique_ptr<Surface> surf;
  switch (type) {
  case SurfaceType::GENERALIZED_BEZIER:
    surf = generateGB(sorted_crosses, corners, normals);
    break;
  case SurfaceType::SUPERD:
    surf = generateSuperD(sorted_crosses, corners, normals);
    break;
  case SurfaceType::SPATCH:
    surf = generateSPatch(sorted_crosses, corners, normals);
    break;
  default:
    assert(false && "Invalid patch type");
  }

  result.emplace_back(surf.release());
  return result;
}

auto sphere(const Point3D &origin, double radius) {
  return std::make_pair([=](const Point3D &p) { return (p - origin).norm() - radius; },
                        [=](const Point3D &p) { return (p - origin) / (p - origin).norm(); });
}

template<typename F1, typename DF1, typename F2, typename DF2>
auto multiply(std::pair<F1,DF1> fdf1, std::pair<F2,DF2> fdf2) {
  return std::make_pair([=](const Point3D &p) { return fdf1.first(p) * fdf2.first(p); },
                        [=](const Point3D &p) {
                          return fdf1.second(p) * fdf2.first(p) + fdf2.second(p) * fdf1.first(p);
                        });
}

auto gyroid() {
  return std::make_pair([](const Point3D &p) {
      return cos(p[0]) * sin(p[1]) + cos(p[1]) * sin(p[2]) + cos(p[2]) * sin(p[0]);
    }, [](const Point3D &p) -> Vector3D {
      return {
          cos(p[2]) * cos(p[0]) - sin(p[0]) * sin(p[1]),
          cos(p[0]) * cos(p[1]) - sin(p[1]) * sin(p[2]),
          cos(p[1]) * cos(p[2]) - sin(p[2]) * sin(p[0])
            };
    });
}

void writeBoundaries(const std::vector<std::unique_ptr<Surface>> &surfaces,
                     const std::string &filename, size_t resolution) {
  std::ofstream f(filename);
  std::vector<size_t> sizes;
  for (const auto &s : surfaces) {
    size_t n = s->domain()->size();
    sizes.push_back(n * resolution);
    for (size_t i = 0; i < n; ++i) {
      auto curve = s->ribbon(i)->curve();
      for (size_t j = 0; j < resolution; ++j) {
        double u = (double)j / resolution;
        auto p = curve->eval(u);
        f << "v " << p[0] << ' ' << p[1] << ' ' << p[2] << std::endl;
      }
    }
  }
  size_t sum = 0;
  for (size_t s : sizes) {
    for (size_t i = 1; i < s; ++i)
      f << "l " << sum + i << ' ' << sum + i + 1 << std::endl;
    f << "l " << sum + s << ' ' << sum + 1 << std::endl;
    sum += s;
  }
}

int main() {
  size_t resolution = 30;
  SurfaceType type = SurfaceType::SPATCH;
  Cell cell({ -3, -3, -3 }, 6.1);
  cell.init(gyroid(), 2, 2);
  // Cell cell({ -1.6, -1.6, -1.6 }, 3);
  // cell.init(sphere({ 0, 0, 0 }, 1), 2, 2);
  // Cell cell({ 0, 0, 0 }, 1);
  // cell.init(multiply(sphere({-0.1, 0, 0}, 0.5), sphere({1.2, 0.9, 0.1}, 0.6)), 0, 0);
  auto surfaces = cell.surfaces(type);
  std::cout << "Generated " << surfaces.size() << " surfaces." << std::endl;
  for (size_t i = 0; i < surfaces.size(); ++i) {
    std::stringstream s;
    s << "/tmp/cell-" << i;
    switch (type) {
    case SurfaceType::GENERALIZED_BEZIER:
      {
        auto *gb = dynamic_cast<SurfaceGeneralizedBezier*>(surfaces[i].get());
        // saveBezier(*gb, s.str() + ".gbp");
        writeBezierControlPoints(*gb, s.str() + "-cp.obj");
      }
      break;
    case SurfaceType::SUPERD:
      {
        auto *sd = dynamic_cast<SurfaceSuperD*>(surfaces[i].get());
        writeSuperDControlPoints(*sd, s.str() + "-cp.obj");
      }
      break;
    default:
      ;
    }
    surfaces[i]->eval(resolution).writeOBJ(s.str() + ".obj");
  }
  writeBoundaries(surfaces, "/tmp/boundaries.obj", 50);
}
