#include <algorithm>
#include <cassert>
#include <fstream>
#include <iostream>
#include <optional>
#include <sstream>

#include <domain.hh>
#include <ribbon.hh>
#include <surface-generalized-bezier.hh>

#include "gb-io.hh"

using namespace Geometry;
using Surface = Transfinite::SurfaceGeneralizedBezier;

class Volume {
public:
  Volume(const Point3D &origin, double length, size_t size);

  template<typename F, typename DF>
  void init(std::pair<F,DF> fdf);

  std::vector<std::shared_ptr<Surface>> getSurfaces() const;

private:
  using Index = std::array<size_t, 3>;
  Point3D position(const Index &i) const;
  size_t idx(const Index &i) const;
  size_t cellIdx(const Index &i) const;
  std::array<std::optional<std::pair<int,Volume::Index>>,6> neighbors(const Index &i) const;
  double value(const Index &i) const;
  double &value(const Index &i);
  Vector3D gradient(const Index &i) const;
  Vector3D &gradient(const Index &i);
  std::pair<Index, Index> edge(const Index &i, int e) const;
  static Index vertex(const Index &index, int i);
  int findCurve(const Index &index, int face1, int face2) const;
  std::shared_ptr<Surface> generateSurfacePass1(const Index &index,
                                                const PointVector &points,
                                                const VectorVector &normals) const;
  std::shared_ptr<Surface> generateSurfacePass2(const Index &index) const;

  Point3D origin;
  double length, step;
  size_t size;
  std::vector<double> values;
  std::vector<Vector3D> gradients;
  std::vector<std::vector<int>> crosses;
  std::vector<std::shared_ptr<Surface>> surfaces;
};

Volume::Volume(const Point3D &origin, double length, size_t size)
  : origin(origin), length(length), size(size)
{
  step = length / size;
}

Point3D
Volume::position(const Index &i) const {
  return origin + Vector3D(i[0], i[1], i[2]) * step;
}

size_t
Volume::idx(const Index &i) const {
  return i[0] * std::pow(size + 1, 2) + i[1] * (size + 1) + i[2];
}

size_t
Volume::cellIdx(const Index &i) const {
  return i[0] * std::pow(size, 2) + i[1] * size + i[2];
}

std::array<std::optional<std::pair<int,Volume::Index>>,6>
Volume::neighbors(const Index &i) const {
  std::array<std::optional<std::pair<int,Index>>,6> result;
  size_t index = 0;
  for (size_t j = 0; j < 3; ++j, index += 2) {
    if (i[j] > 0) {
      result[index] = std::make_pair(index, Index());
      std::copy(i.begin(), i.end(), result[index].value().second.begin());
      result[index].value().second[j]--;
    }
    if (i[j] < size - 1) {
      result[index+1] = std::make_pair(index + 1, Index());
      std::copy(i.begin(), i.end(), result[index+1].value().second.begin());
      result[index+1].value().second[j]++;
    }
  }
  return result;
}

double
Volume::value(const Index &i) const {
  return values[idx(i)];
}

double &
Volume::value(const Index &i) {
  return values[idx(i)];
}

Vector3D
Volume::gradient(const Index &i) const {
  return gradients[idx(i)];
}

Vector3D &
Volume::gradient(const Index &i) {
  return gradients[idx(i)];
}

std::pair<Volume::Index, Volume::Index>
Volume::edge(const Index &i, int edge_index) const {
  static const std::array<size_t,12*6> indices = {
    0,0,0, 1,0,0,
    1,0,0, 1,1,0,
    1,1,0, 0,1,0,
    0,1,0, 0,0,0,
    0,0,1, 1,0,1,
    1,0,1, 1,1,1,
    1,1,1, 0,1,1,
    0,1,1, 0,0,1,
    0,0,0, 0,0,1,
    1,0,0, 1,0,1,
    1,1,0, 1,1,1,
    0,1,0, 0,1,1
  };
  size_t k = edge_index * 6;
  size_t a = i[0] + indices[k++];
  size_t b = i[1] + indices[k++];
  size_t c = i[2] + indices[k++];
  size_t d = i[0] + indices[k++];
  size_t e = i[1] + indices[k++];
  size_t f = i[2] + indices[k];
  return { { a, b, c }, { d, e, f } };
}

namespace {

    static std::array<int,6*4> faces = {
      3,  7,  8, 11,
      1,  5,  9, 10,
      0,  4,  8,  9,
      2,  6, 10, 11,
      0,  1,  2,  3,
      4,  5,  6,  7
    };

}

static bool samePlane(int a, int b) {
  size_t index = 0;
  for (size_t i = 0; i < 6; ++i) {
    size_t found = 0;
    for (size_t j = 0; j < 4; ++j, ++index)
      if (faces[i*4+j] == a || faces[i*4+j] == b)
        ++found;
    if (found == 2)
      return true;
  }
  return false;
}

static bool edgeInFace(int edge, int face) {
  auto start = faces.begin() + face * 4, end = start + 4;
  return std::find(start, end, edge) != end;
}

static const Vector3D &
planeNormal(int i, int j, int k) {
  static const std::array<std::array<int,4>,6> faces = {
    {
      {0, 1, 2, 3}, {4, 5, 6, 7}, {1, 5, 6, 2},
      {0, 4, 7, 3}, {0, 1, 5, 4}, {3, 2, 6, 7}
    }
  };
  static const std::array<Vector3D,6> planes = {
    {
      { 0, 0, -1}, {0,  0, 1}, {1, 0, 0},
      {-1, 0,  0}, {0, -1, 0}, {0, 1, 0}
    }
  };
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

Volume::Index
Volume::vertex(const Index &index, int i) {
  static std::array<size_t,8*3> vertices = {
    0,0,0, 1,0,0, 1,1,0, 0,1,0, 0,0,1, 1,0,1, 1,1,1, 0,1,1
  };
  return { index[0] + vertices[i*3+0],
           index[1] + vertices[i*3+1],
           index[2] + vertices[i*3+2] };
}

std::shared_ptr<Surface>
Volume::generateSurfacePass1(const Index &index,
                             const PointVector &points,
                             const VectorVector &normals) const {
  int sides = points.size();
  auto surf = std::make_shared<Surface>();
  surf->initNetwork(sides, 3);

  // Corner control points are already computed
  for (int i = 0; i < sides; ++i)
    surf->setControlPoint(i, 0, 0, points[i]);

  static const std::array<std::pair<int,int>,12> edges = {
    {
      {0, 1}, {1, 2}, {2, 3}, {3, 0},
      {4, 5}, {5, 6}, {6, 7}, {7, 4},
      {0, 4}, {1, 5}, {2, 6}, {3, 7}
    }
  };

  // Tangent control points
  for (int i = 0; i < sides; ++i) {
    // Find the 3 or 4 vertices defining this curve,
    // and mark the central vertex in the 3-vertex case
    int ip = (i + 1) % sides;
    size_t ci = cellIdx(index);
    std::array<int,4> vertices = { 
      edges[crosses[ci][i]].first,
      edges[crosses[ci][i]].second,
      edges[crosses[ci][ip]].first,
      edges[crosses[ci][ip]].second
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
      auto center = position(vertex(index, central_vertex));
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

  surf->setupLoop();            // for domain()-size()

  return surf;
}

int
Volume::findCurve(const Index &index, int face1, int face2) const {
  size_t ci = cellIdx(index);
  int sides = surfaces[ci]->domain()->size();
  for (int i = 0; i < sides; ++i) {
    int ip = (i + 1) % sides;
    if (edgeInFace(crosses[ci][i], face1) &&
        edgeInFace(crosses[ci][ip], face2) &&
        !edgeInFace(crosses[ci][ip], face1) &&
        edgeInFace(crosses[ci][i], face2))
      return i + 1;
    if (edgeInFace(crosses[ci][i], face2) &&
        edgeInFace(crosses[ci][ip], face1) &&
        !edgeInFace(crosses[ci][i], face1) &&
        edgeInFace(crosses[ci][ip], face2))
      return -(i + 1);
  }
  return 0;
}

static Point3D extractCornerCP(const std::shared_ptr<Surface> &s, int c) {
  if (c < 0)
    return s->controlPoint(-c-1, 3, 0);
  return s->controlPoint(c-1, 0, 0);
}

static Point3D extractTangentCP(const std::shared_ptr<Surface> &s, int c) {
  if (c < 0)
    return s->controlPoint(-c-1, 2, 0);
  return s->controlPoint(c-1, 1, 0);
}

static void setTangentCP(std::shared_ptr<Surface> &s, int c, const Point3D &p) {
  if (c < 0)
    s->setControlPoint(-c-1, 2, 0, p);
  else
    s->setControlPoint(c-1, 1, 0, p);
}

std::shared_ptr<Surface>
Volume::generateSurfacePass2(const Index &index) const {
  auto surf = surfaces[cellIdx(index)];
  if (!surf)
    return surf;
  int sides = surf->domain()->size();

  // Modify tangent control points based on adjacent patches
  for (const auto &neighbor : neighbors(index)) {
    if (!neighbor.has_value())
      continue;
    auto [f, index2] = neighbor.value();
    if (f % 2 == 1)
      continue;                 // for symmetry reasons
    auto surf2 = surfaces[cellIdx(index2)];
    if (!surf2)
      continue;

    int f2 = f + 1;
    for (int other_face = 0; other_face < 6; ++other_face) {
      if (other_face == f || other_face == f2)
        continue;
      int c = findCurve(index, f, other_face);
      int c2 = findCurve(index2, f2, other_face);
      if (c * c2 == 0)
        continue;
      auto ccp = extractCornerCP(surf, c);
      auto tcp = extractTangentCP(surf, c);
      auto tcp2 = extractTangentCP(surf2, c2);
      auto d = (tcp2 - tcp) / 2;
      setTangentCP(surf, c, ccp - d);
      setTangentCP(surf2, c2, ccp + d);
    }
  }

  // Twist control points by the parallelogram rule
  Point3D p;
  for (int i = 0; i < sides; ++i) {
    p = surf->controlPoint(i, 0, 1) + surf->controlPoint(i, 1, 0) - surf->controlPoint(i, 0, 0);
    surf->setControlPoint(i, 1, 1, p);
  }

  // Central control point is computed from the mass center of the corner & twist control points
  p = Point3D(0, 0, 0);
  auto q = Point3D(0, 0, 0);
  for (int i = 0; i < sides; ++i) {
    p += surf->controlPoint(i, 1, 1);
    q += surf->controlPoint(i, 0, 0);
  }
  surf->setCentralControlPoint((p * 2 - q) / sides);

  surf->setupLoop();

  return surf;
}

template<typename F, typename DF>
void
Volume::init(std::pair<F,DF> fdf) {
  values.clear(); values.resize(std::pow(size + 1, 3));
  gradients.clear(); gradients.resize(std::pow(size + 1, 3));
  for (size_t i = 0; i <= size; ++i)
    for (size_t j = 0; j <= size; ++j)
      for (size_t k = 0; k <= size; ++k) {
        Index index = { i, j, k };
        auto p = position(index);
        value(index) = fdf.first(p);
        gradient(index) = fdf.second(p);
      }

  crosses.clear(); crosses.resize(std::pow(size, 3));
  surfaces.clear(); surfaces.resize(std::pow(size, 3));
  for (size_t i = 0; i < size; ++i)
    for (size_t j = 0; j < size; ++j)
      for (size_t k = 0; k < size; ++k) {
        Index index = { i, j, k };
        std::vector<int> tmp_crosses;
        for (int l = 0; l < 12; ++l) {
          auto e = edge(index, l);
          if (value(e.first) * value(e.second) < 0)
            tmp_crosses.push_back(l);
        }
        int n_crosses = tmp_crosses.size();
        if (n_crosses == 0)
          continue;

        std::vector<int> sorted_crosses;
        PointVector corners;
        VectorVector normals;
        int cross = 0, last_cross = -1;
        do {
          auto e = edge(index, tmp_crosses[cross]);
          Index i1 = e.first, i2 = e.second;
          double v1 = value(i1), v2 = value(i2);
          double alpha = std::abs(v1) / std::abs(v2 - v1);
          sorted_crosses.push_back(tmp_crosses[cross]);
          corners.push_back(position(i1) * (1 - alpha) + position(i2) * alpha);
          normals.push_back(gradient(i1) * (1 - alpha) + gradient(i2) * alpha);
          for (int j = 0; j < n_crosses; ++j)
            if (j != last_cross && j != cross &&
                samePlane(tmp_crosses[cross], tmp_crosses[j])) {
              last_cross = cross;
              cross = j;
              break;
            }
        } while (cross != 0);
        int sides = corners.size();

        if (sides != n_crosses) // ambiguous case
          continue;

        // Reverse the loop if needed, such that positive is outside
        static const std::array<int,12> left_faces = { 4, 4, 4, 4, 2, 1, 3, 0, 1, 1, 2, 3 };
        auto first_edge = edge(index, sorted_crosses[0]);
        bool negative = edgeInFace(sorted_crosses[1], left_faces[sorted_crosses[0]]);
        if ((negative && value(first_edge.first) > 0) ||
            (!negative && value(first_edge.first) < 0)) {
          std::reverse(sorted_crosses.begin(), sorted_crosses.end());
          std::reverse(corners.begin(), corners.end());
          std::reverse(normals.begin(), normals.end());
        }

        crosses[cellIdx(index)] = sorted_crosses;
        surfaces[cellIdx(index)] = generateSurfacePass1(index, corners, normals);
      }

  std::vector<std::shared_ptr<Surface>> tmp(std::pow(size, 3));
  for (size_t i = 0; i < size; ++i)
    for (size_t j = 0; j < size; ++j)
      for (size_t k = 0; k < size; ++k) {
        Index index = { i, j, k };
        tmp[cellIdx(index)] = generateSurfacePass2(index);
      }
  surfaces = tmp;
}

std::vector<std::shared_ptr<Surface>>
Volume::getSurfaces() const {
  std::vector<std::shared_ptr<Surface>> result;
  std::copy_if(surfaces.begin(), surfaces.end(), std::back_inserter(result),
               [](const std::shared_ptr<Surface> &sp) { return (bool)sp; });
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

template<typename F, typename DF>
auto normalized(std::pair<F,DF> fdf) {
  return std::make_pair([=](const Point3D &p) { return fdf.first(p) / fdf.second(p).norm(); },
                        fdf.second);
}

void writeBoundaries(const std::vector<std::shared_ptr<Surface>> &surfaces,
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
  bool generate_controlnet = false;
  size_t resolution = 30;
  // Volume volume({ -3, -3, -3 }, 6.1, 8);
  // volume.init(normalized(gyroid()));
  Volume volume({ -1.6, -1.6, -1.6 }, 3.0, 4);
  volume.init(sphere({ 0, 0, 0 }, 1));
  // Volume volume({ 0, 0, 0 }, 1.0, 2);
  // volume.init(multiply(sphere({-0.1, 0, 0}, 0.5), sphere({1.2, 0.9, 0.1}, 0.6)));
  auto surfaces = volume.getSurfaces();
  std::cout << "Generated " << surfaces.size() << " surfaces." << std::endl;
  for (size_t i = 0; i < surfaces.size(); ++i) {
    std::stringstream s;
    s << "/tmp/cell-" << i;
    surfaces[i]->eval(resolution).writeOBJ(s.str() + ".obj");
    if (!generate_controlnet)
      continue;
    // saveBezier(*surfaces[i], s.str() + ".gbp");
    writeBezierControlPoints(*surfaces[i], s.str() + "-cp.obj");
  }
  writeBoundaries(surfaces, "/tmp/boundaries.obj", 50);
}
