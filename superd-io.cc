#include <fstream>

#include <domain.hh>

#include "superd-io.hh"

static void writePoint(std::ostream &os, const Geometry::Point3D &p) {
  os << "v " << p[0] << ' ' << p[1] << ' ' << p[2] << std::endl;
}

void writeSuperDControlPoints(const Transfinite::SurfaceSuperD &surf, const std::string &filename) {
  std::ofstream f(filename);
  size_t n = surf.domain()->size();
  writePoint(f, surf.vertexControlPoint());
  for (size_t i = 0; i < n; ++i) {
    writePoint(f, surf.edgeControlPoint(i));
    writePoint(f, surf.faceControlPoint(i));
  }
  for (size_t i = 0; i < n; ++i) {
    size_t ip = (i + 1) % n;
    f << "l " << 2 * i + 2 << ' ' << 2 * i + 3 << std::endl;
    f << "l " << 2 * i + 3 << ' ' << 2 * ip + 2 << std::endl;
    f << "l 1 " << 2 * i + 2 << std::endl;
  }
}
