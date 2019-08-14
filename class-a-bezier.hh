#pragma once

#include "geometry.hh"

// p0 and p1 are the endpoints
// t0 and t1 are the end tangents (should be unit vectors)
// Both tangents should be directed towards the curve
Geometry::Point2DVector
fitCubicClassABezier(const Geometry::Point2D &p0, const Geometry::Vector2D &t0,
                     const Geometry::Point2D &p1, const Geometry::Vector2D &t1);

// Same, but with 3D points/vectors
// Assumes that everything is in a common plane
Geometry::PointVector
fitCubicClassABezier(const Geometry::Point3D &p0, const Geometry::Vector3D &t0,
                     const Geometry::Point3D &p1, const Geometry::Vector3D &t1);
