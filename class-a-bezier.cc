#include <cmath>
#include <exception>

#include "class-a-bezier.hh"

using namespace Geometry;

class ClassABezierException : public std::exception { };

static Vector2D
rotateVector2D(const Vector2D &v, double phi) {
  double c = cos(phi), s = sin(phi);
  return { c * v[0] - s * v[1], s * v[0] + c * v[1] };
}

static std::pair<double, double>
solveQuadraticEquation(double a, double b, double c) {
  double eps = 1.0e-15;
  if (std::abs(a) < eps) {
    if (std::abs(b) < eps) {
      if (std::abs(c) < eps)
        return { 0, 0 };
      throw ClassABezierException(); // no solution
    }
    return { -a / b, -a / b };
  }
  double D = b * b - 4 * a * c;
  if (D < -eps)
    throw ClassABezierException(); // no real solution
  if (D < 0)
    D = 0;
  D = std::sqrt(D);
  return { (-b + D) / (2 * a), (-b - D) / (2 * a) };
}

static std::pair<double, double>
solveSystem(const Vector2D &t0, const Vector2D &t1, const Vector2D &r, const Vector2D &d) {
  auto v0 = t0 * -1, v1 = r * -1, v2 = t1, w = d * -1;
  // Normalized to: v2 s^2 + v1 s + v0 = w / s0
  int x = std::abs(w[0]) > std::abs(w[1]) ? 0 : 1, y = 1 - x;
  double ww = w[y] / w[x];
  double a = v2[y] - ww * v2[x], b = v1[y] - ww * v1[x], c = v0[y] - ww * v0[x];
  auto [s1, s2] = solveQuadraticEquation(a, b, c);
  double s = std::max(s1, s2);
  double s0 = w[x] / (v2[x] * s * s + v1[x] * s + v0[x]);
  return { s0, s };
}

Point2DVector
fitCubicClassABezier(const Point2D &p0, const Vector2D &t0,
                     const Point2D &p1, const Vector2D &t1) {
  try {
    double phi = acos(std::min(std::max(-t0 * t1, -1.0), 1.0)) / 2;
    if (Vector2D(t0[1], -t0[0]) * t1 < 0)
      phi *= -1;
    auto r = rotateVector2D(t0, phi);
    auto [s0, s] = solveSystem(t0, t1, r, p1 - p0);
    double s1 = s0 * s * s;
    return { p0, p0 + t0 * s0, p1 + t1 * s1, p1 };
  } catch(ClassABezierException &) {
    return { };
  }
}

PointVector
fitCubicClassABezier(const Point3D &p0, const Vector3D &t0,
                     const Point3D &p1, const Vector3D &t1) {
  auto n = (t0 ^ t1).normalize();
  auto u = (p1 - p0).normalize();
  auto v = n ^ u;
  auto to2D = [&](const Vector3D &w) -> Vector2D { return { w * u, w * v }; };
  auto pv = fitCubicClassABezier({0, 0}, to2D(t0), to2D(p1 - p0), to2D(t1));
  PointVector result;
  for (const auto &p : pv)
    result.push_back(p0 + u * p[0] + v * p[1]);
  return result;
}
