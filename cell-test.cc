#include <array>
#include <cassert>
#include <iostream>
#include <sstream>

#include <surface-generalized-bezier.hh>

#include "gb-io.hh"

using namespace Geometry;
using Surface = Transfinite::SurfaceGeneralizedBezier;

double sphere(const Point3D &p) {
  return p.norm() - 1;
}

Vector3D sphereGradient(const Point3D &p) {
  return p / p.norm();
}

class Cell {
public:
  Cell(const Point3D &origin, double length) : origin(origin), length(length) { }
  void init(double (*f)(const Point3D &), Vector3D (*df)(const Point3D &));
  double value(int i) const { return values[i]; }
  const Vector3D &gradient(int i) const { return gradients[i]; }
  Point3D vertex(int i) const;
  const Vector3D &planeNormal(int i, int j, int k) const;
  Surface GB3sided() const;
private:
  Point3D origin;
  double length;
  std::array<double, 8> values;
  std::array<Vector3D, 8> gradients;
  static const std::array<std::pair<int,int>,12> edges;
  static const std::array<std::array<int,4>,6> faces;
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

const std::array<std::pair<int,int>,12> Cell::edges = {
  {
    {0, 1}, {1, 2}, {2, 3}, {3, 0},
    {4, 5}, {5, 6}, {6, 7}, {7, 4},
    {1, 5}, {4, 0}, {2, 6}, {7, 3}
  }
};

const std::array<std::array<int,4>,6> Cell::faces = {
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

void
Cell::init(double (*f)(const Point3D &), Vector3D (*df)(const Point3D &)) {
  for (int i = 0; i < 8; ++i) {
    values[i] = f(vertex(i));
    gradients[i] = df(vertex(i));
  }
}

Surface
Cell::GB3sided() const {
  std::vector<int> crosses;
  for (int i = 0; i < 12; ++i)
    if (value(edges[i].first) * value(edges[i].second) < 0)
      crosses.push_back(i);
  assert(crosses.size() == 3);

  std::array<int, 4> vertices;              // vertices[0] is the "central" vertex
  std::array<bool, 12> seen;
  seen.fill(false);
  for (int c : crosses) {
    if (seen[edges[c].first])
      vertices[0] = edges[c].first;
    else
      seen[edges[c].first] = true;
    if (seen[edges[c].second])
      vertices[0] = edges[c].second;
    else
      seen[edges[c].second] = true;
  }
  int index = 1;
  for (int i = 0; i < 12; ++i)
    if (seen[i] && i != vertices[0])
      vertices[index++] = i;

  std::array<Point3D, 3> corners;
  std::array<Vector3D, 3> normals;
  for (int i = 1; i <= 3; ++i) {
    int i1 = vertices[0], i2 = vertices[i];
    double v1 = value(i1), v2 = value(i2);
    double length = std::abs(v2 - v1), alpha = std::abs(v1) / length;
    corners[i-1] = vertex(i1) * (1 - alpha) + vertex(i2) * alpha;
    normals[i-1] = gradient(i1) * (1 - alpha) + gradient(i2) * alpha;
  }

  Surface surf;
  surf.initNetwork(3, 3);

  // Corner control points are already computed
  for (int i = 0; i < 3; ++i)
    surf.setControlPoint(i, 0, 0, corners[i]);

  // Tangent control points
  // Questions:
  // - should we normalize the gradients?
  // - how should we scale the tangents?
  auto center = vertex(vertices[0]);
  for (int i = 0; i < 3; ++i) {
    int ip = (i + 1) % 3;
    auto p1 = corners[i], p2 = corners[ip];
    auto n1 = normals[i], n2 = normals[ip];
    auto pn = planeNormal(vertices[0], vertices[i+1], vertices[ip+1]);
    auto t1 = n1.normalize() ^ pn, t2 = n2.normalize() ^ pn;
    if (t1 * (p2 - center) < 0)
      t1 *= -1;
    if (t2 * (p1 - center) < 0)
      t2 *= -1;
    double scaling = (p1 - p2).norm() / 3.0; // kutykurutty
    surf.setControlPoint(i, 1, 0, p1 + t1 * scaling);
    surf.setControlPoint(i, 2, 0, p2 + t2 * scaling);
  }

  // Twist control points by the parallelogram rule
  Point3D p;
  for (int i = 0; i < 3; ++i) {
    p = surf.controlPoint(i, 0, 1) + surf.controlPoint(i, 1, 0) - surf.controlPoint(i, 0, 0);
    surf.setControlPoint(i, 1, 1, p);
  }

  // Central control point is the mass center of the twist control points
  p = (surf.controlPoint(0, 1, 1) + surf.controlPoint(1, 1, 1) + surf.controlPoint(2, 1, 1)) / 3;
  surf.setCentralControlPoint(p);

  surf.setupLoop();
  return surf;
}

int main() {
  // Cell cell({ 0, -1, 0.5 }, 1);
  // cell.init(sphere, sphereGradient);
  // auto surface = cell.GB3sided();
  // saveBezier(surface, "/tmp/cell.gbp");
  // writeBezierControlPoints(surface, "/tmp/cell-cp.obj");
  // surface.eval(100).writeOBJ("/tmp/cell.obj");
  for (int i = 0; i <= 1; ++i)
    for (int j = 0; j <= 1; ++j)
      for (int k = 0; k <= 1; ++k) {
        double l = 1.2;
        auto p = Point3D(0.1, 0.1, 0.1) + (Point3D(i, j, k) - Point3D(1, 1, 1)) * l;
        Cell cell(p, l);
        cell.init(sphere, sphereGradient);
        auto surface = cell.GB3sided();
        std::stringstream s;
        s << "/tmp/cell-" << i << j << k << ".obj";
        surface.eval(100).writeOBJ(s.str());
      }
}
