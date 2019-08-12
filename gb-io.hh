#pragma once

#include <surface-generalized-bezier.hh>

void saveBezier(const Transfinite::SurfaceGeneralizedBezier &surf, const std::string &filename);

void writeBezierControlPoints(const Transfinite::SurfaceGeneralizedBezier &surf,
                              const std::string &filename);

void writeIncompatibleBezierControlPoints(const Transfinite::SurfaceGeneralizedBezier &surf,
                                          const std::string &filename);
