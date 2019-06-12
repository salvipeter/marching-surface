#include "gb-io.hh"

#include <fstream>
#include <iostream>

#include <domain.hh>

using namespace Geometry;

void saveBezier(const Transfinite::SurfaceGeneralizedBezier &surf, const std::string &filename) {
  std::ofstream f(filename);
  if (!f.is_open()) {
    std::cerr << "Unable to open file: " << filename << std::endl;
    return;
  }

  size_t n = surf.domain()->vertices().size();
  size_t d = surf.degree();
  f << n << ' ' << d << std::endl;
  size_t l = (d + 1) / 2;
  size_t cp = 1 + d / 2;
  cp = n * cp * l + 1;          // # of control points

  Point3D p = surf.centralControlPoint();
  f << p[0] << ' ' << p[1] << ' ' << p[2] << std::endl;

  for (size_t i = 1, side = 0, col = 0, row = 0; i < cp; ++i, ++col) {
    if (col >= d - row) {
      if (++side >= n) {
        side = 0;
        ++row;
      }
      col = row;
    }
    p = surf.controlPoint(side, col, row);
    f << p[0] << ' ' << p[1] << ' ' << p[2] << std::endl;
  }
}

void writeBezierControlPoints(const Transfinite::SurfaceGeneralizedBezier &surf,
                              const std::string &filename) {
  // Slow but simple implementation creating a nice mesh
  size_t n = surf.domain()->vertices().size();
  size_t d = surf.degree();
  size_t l = surf.layers();
  size_t cp = 1 + d / 2;
  cp = n * cp * l + 1;

  auto findControlPoint = [n,d,cp](size_t i, size_t j, size_t k) -> size_t {
    for (size_t c = 1, side = 0, col = 0, row = 0; c < cp; ++c, ++col) {
      if (col >= d - row) {
        if (++side >= n) {
          side = 0;
          ++row;
        }
        col = row;
      }
      size_t side_m = (side + n - 1) % n, side_p = (side + 1) % n;
      if ((i == side && j == col && k == row) ||
          (i == side_m && j == d - row && k == col) ||
          (i == side_p && j == row && k == d - col))
        return c + 1;
    }
    return 0;
  };

  std::ofstream f(filename);
  Point3D p = surf.centralControlPoint();
  f << "v " << p[0] << " " << p[1] << " " << p[2] << std::endl;
  for (size_t i = 1, side = 0, col = 0, row = 0; i < cp; ++i, ++col) {
    if (col >= d - row) {
      if (++side >= n) {
        side = 0;
        ++row;
      }
      col = row;
    }
    p = surf.controlPoint(side, col, row);
    f << "v " << p[0] << " " << p[1] << " " << p[2] << std::endl;
  }

  for (size_t i = 0; i < n; ++i)
    for (size_t j = 0; j <= d / 2; ++j)
      for (size_t k = 0; k < l - 1; ++k) {
        size_t a = findControlPoint(i, j, k);
        size_t b = findControlPoint(i, j + 1, k);
        size_t c = findControlPoint(i, j + 1, k + 1);
        size_t d = findControlPoint(i, j, k + 1);
        f << "f " << a << " " << b << " " << c << " " << d << std::endl;
      }
  if (d % 2 == 0)
    for (size_t i = 0; i < n; ++i) {
      size_t im = (i + n - 1) % n;
      size_t a = findControlPoint(i, l - 1, l - 1);
      size_t b = findControlPoint(i, l, l - 1);
      size_t c = 1;
      size_t d = findControlPoint(im, l, l - 1);
      f << "f " << a << " " << b << " " << c << " " << d << std::endl;
    }
  else
    for (size_t i = 0; i < n; ++i) {
      size_t a = findControlPoint(i, l - 1, l - 1);
      size_t b = findControlPoint(i, l, l - 1);
      size_t c = 1;
      f << "f " << a << " " << b << " " << c << std::endl;
    }

  f.close();
}
