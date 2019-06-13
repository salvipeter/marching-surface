#include <algorithm>
#include <array>
#include <cassert>
#include <iostream>
#include <set>
#include <sstream>

#include <surface-generalized-bezier.hh>

#include "gb-io.hh"

using namespace Geometry;
using Surface = Transfinite::SurfaceGeneralizedBezier;

class Cell {
public:
  Cell(const Point3D &origin, double length) : origin(origin), length(length) { }
  void init(double (*f)(const Point3D &), Vector3D (*df)(const Point3D &));
  double value(int i) const { return values[i]; }
  const Vector3D &gradient(int i) const { return gradients[i]; }
  Point3D vertex(int i) const;
  const Vector3D &planeNormal(int i, int j, int k) const;
  Surface surface() const;
private:
  using Edge = std::pair<int,int>;
  using Face = std::array<int,4>;
  static bool samePlane(const Edge &e1, const Edge &e2);

  Point3D origin;
  double length;
  std::array<double, 8> values;
  std::array<Vector3D, 8> gradients;
  static const std::array<Edge,12> edges;
  static const std::array<Face,6> faces;
  static const std::array<Vector3D,6> planes;
  /*
       7         6
      +--------+
     /|       /|
  3 / |    2 / |
   +--------+  |          y
   |  +-----|--+          ^
   | / 4    | / 5         | z
   |/       |/            |/
 0 +--------+ 1           +---> x
   */
};

const std::array<Cell::Edge,12> Cell::edges = {
  {
    {0, 1}, {1, 2}, {2, 3}, {3, 0},
    {4, 5}, {5, 6}, {6, 7}, {7, 4},
    {1, 5}, {4, 0}, {2, 6}, {7, 3}
  }
};

const std::array<Cell::Face,6> Cell::faces = {
  {
    {0, 1, 2, 3}, {4, 5, 6, 7}, {0, 1, 5, 4},
    {3, 2, 6, 7}, {1, 5, 6, 2}, {0, 4, 7, 3}
  }
};

const std::array<Vector3D,6> Cell::planes = {
  {
    {0, 0, -1}, {0, 0, 1}, { 0, -1, 0},
    {0, 1,  0}, {1, 0, 0}, {-1,  0, 0}
  }
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
          return std::find(face.begin(), face.end(), vertex) != face.end();
        });
    });
}

void
Cell::init(double (*f)(const Point3D &), Vector3D (*df)(const Point3D &)) {
  for (int i = 0; i < 8; ++i) {
    values[i] = f(vertex(i));
    gradients[i] = df(vertex(i));
  }
}

Surface
Cell::surface() const {
  std::vector<int> crosses;
  for (int i = 0; i < 12; ++i)
    if (value(edges[i].first) * value(edges[i].second) < 0)
      crosses.push_back(i);
  int n_crosses = crosses.size();

  std::vector<int> sorted_crosses;
  std::vector<Point3D> corners;
  std::vector<Vector3D> normals;
  int cross = 0, last_cross = -1;
  do {
    int i1 = edges[crosses[cross]].first, i2 = edges[crosses[cross]].second;
    double v1 = value(i1), v2 = value(i2);
    double alpha = std::abs(v1) / std::abs(v2 - v1);
    sorted_crosses.push_back(cross);
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

  assert(n_crosses == sides && "Ambiguous case!");

  Surface surf;
  surf.initNetwork(sides, 3);

  // Corner control points are already computed
  for (int i = 0; i < sides; ++i)
    surf.setControlPoint(i, 0, 0, corners[i]);

  // Tangent control points
  // Questions:
  // - should we normalize the gradients?
  // - how should we scale the tangents?
  for (int i = 0; i < sides; ++i) {
    // Find the 3 or 4 vertices defining this curve,
    // and mark the central vertex in the 3-vertex case
    int ip = (i + 1) % sides;
    std::array<int,4> vertices = { 
      edges[crosses[sorted_crosses[i]]].first,
      edges[crosses[sorted_crosses[i]]].second,
      edges[crosses[sorted_crosses[ip]]].first,
      edges[crosses[sorted_crosses[ip]]].second
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
    auto p1 = corners[i], p2 = corners[ip];
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
    surf.setControlPoint(i, 1, 0, p1 + t1 * scaling);
    surf.setControlPoint(i, 2, 0, p2 + t2 * scaling);
  }

  // Twist control points by the parallelogram rule
  Point3D p;
  for (int i = 0; i < sides; ++i) {
    p = surf.controlPoint(i, 0, 1) + surf.controlPoint(i, 1, 0) - surf.controlPoint(i, 0, 0);
    surf.setControlPoint(i, 1, 1, p);
  }

  // Central control point is the mass center of the twist control points
  p = Point3D(0, 0, 0);
  for (int i = 0; i < sides; ++i)
    p += surf.controlPoint(i, 1, 1);
  surf.setCentralControlPoint(p / sides);

  surf.setupLoop();
  return surf;
}

double sphere(const Point3D &p) {
  return p.norm() - 1;
}

Vector3D sphereGradient(const Point3D &p) {
  return p / p.norm();
}

int main() {
  Cell cell({ -0.3, -0.7, -1.5 }, 1);
  cell.init(sphere, sphereGradient);
  auto surface = cell.surface();
  saveBezier(surface, "/tmp/cell.gbp");
  writeBezierControlPoints(surface, "/tmp/cell-cp.obj");
  surface.eval(100).writeOBJ("/tmp/cell.obj");
  // for (int i = 0; i <= 1; ++i)
  //   for (int j = 0; j <= 1; ++j)
  //     for (int k = 0; k <= 1; ++k) {
  //       double l = 1.2;
  //       auto p = Point3D(0.1, 0.1, 0.1) + (Point3D(i, j, k) - Point3D(1, 1, 1)) * l;
  //       Cell cell(p, l);
  //       cell.init(sphere, sphereGradient);
  //       auto surface = cell.surface();
  //       std::stringstream s;
  //       s << "/tmp/cell-" << i << j << k << ".obj";
  //       surface.eval(100).writeOBJ(s.str());
  //     }
}
