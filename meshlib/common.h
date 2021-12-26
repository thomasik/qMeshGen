#pragma once

#ifndef COMMON_H__INCLUDED
#define COMMON_H__INCLUDED

#define __STDC_WANT_LIB_EXT1__ 1

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <cmath>
#ifndef _DEBUG
 #ifndef NDEBUG
  #define NDEBUG
 #endif
#endif
#include <cassert>
#include <chrono>
#include <ctime>
using namespace std;

#ifdef _OPENMP
#include <omp.h>
#else
#pragma warning( disable : 161 )
#endif

#ifdef NDEBUG
#define LOG_ASSERT(cond) LOG4CPLUS_ASSERT(MeshLog::logger_mesh, (cond));
#else
#define LOG_ASSERT(cond) assert(cond);
#endif

#include "MeshLog.h"

/// Globally used values
#define PI		3.14159265358979323846
#define SQRT2	1.41421356237309504880
#define SQRT3	1.73205080756887729353

typedef pair<bool, bool> bool_pair;

enum class Axis { X = 0, Y = 1, Z = 2 };

enum OctVertexWhich {
	VX3D_LSW = 0, VX3D_LSE = 1, VX3D_LNW = 2, VX3D_LNE = 3,
	VX3D_HSW = 4, VX3D_HSE = 5, VX3D_HNW = 6, VX3D_HNE = 7,
	VX3D_FIRST = 0, VX3D_LAST = 7, VX3D_COUNT = 8
};

enum OctFaceWhich {
	FC3D_SOUTH = 0, FC3D_NORTH = 1, FC3D_WEST = 2, FC3D_EAST = 3,
	FC3D_LOW = 4,   FC3D_HIGH = 5,  
	FC3D_FIRST = 0, FC3D_LAST = 5, FC3D_COUNT = 6
};

enum OctEdgeWhich {
	ED3D_LS = 0, ED3D_LN = 1, ED3D_LW = 2,  ED3D_LE = 3,
	ED3D_SW = 4, ED3D_SE = 5, ED3D_NW = 6,  ED3D_NE = 7,
	ED3D_HS = 8, ED3D_HN = 9, ED3D_HW = 10, ED3D_HE = 11,
	ED3D_FIRST = 0, ED3D_LAST = 11, ED3D_COUNT = 12
};

/// Types of used objects
enum ElementType { ELEMENT_UNKNOWN,
				   CURVE_UNKNOWN, CURVE_CIRCLE, CURVE_SEGMENT, CURVE_BSPLINE, CURVE_ANALYTIC, CURVE_LQUADRIC,
				   CURVE_3D_UNKNOWN, CURVE_3D_SEGMENT, CURVE_3D_SURFACE_PARAMETRIC, CURVE_3D_BSPLINE, CURVE_3D_LQUADRIC,
				   EDGE_SIMPLE, EDGE_CURVE, EDGE_DOMAIN, EDGE_SIMPLE_3D, EDGE_DOMAIN_3D, 
				   ELEMENT_MESH_AREA, ELEMENT_MESH_TRIANGLE, ELEMENT_MESH_QUAD,
				   FACE_UNKNOWN, FACE_DOMAIN, FACE_TRIANGLE, FACE_QUAD, FACE_POLY,
				   SURFACE_UNKNOWN, SURFACE_PLANE, SURFACE_ANALYTIC, SURFACE_PLANAR_QUADRIC, SURFACE_QUADRIC_QTREE,
				   SURFACE_MONOANALYTIC, SURFACE_BSPLINE, SURFACE_CYLINDER, SURFACE_SPHERE,
				   SURFACE_TRANSLATED, SURFACE_MULTI, SURFACE_CORRECTED,
				   BLOCK_UNKNOWN, BLOCK_DOMAIN, BLOCK_TETRA, BLOCK_HEX,
				   ELEMENT_VERTEX, ELEMENT_EDGE, ELEMENT_BORDER};

/// Square of x
template <class T> T sqr(const T& x){ return x*x; }
/// Sign of x
template <class T> int sgn(const T& x){ return x > 0 ? 1 : (x < 0 ? -1 : 0); }


template<typename T>
double diffKdValue(const T& v1, const T& v2) {
	//return std::abs(v2 - v1);
	assert(v1 > 0.0 && v2 > 0.0);
	double d = std::abs(v1 / v2 + v2 / v1 - 2.0);
	return SQRT3 * d;
}

//#define SHOW_MESH(caption, set)
//#define SHOW_MESH_NORESET(caption, set)
#define SHOW_MESH(caption, set) MeshViewSet::showViewSet(caption, set)
#define SHOW_MESH_NORESET(caption, set) MeshViewSet::showViewSetNoReset(caption, set)

#define VERY_SMALL_NUMBER 1e-40
#define SMALL_NUMBER 1e-10
#define LARGE_NUMBER 1e20
#define METRIC_SMALL_NUMBER 1e-5
#define MIN_PARAM_GRATIO 1e-10

#define MEASURE_PRECISION 0.0002

#define METRIC_LENGTH_RATIO 2.0
#define METRIC_LENGTH_RATIO2 4.0

#define AQ_UNKNOWN -1.0
#define AQ_INVALID -0.5
#define AQ_VALID_MIN 0.0
#define AQ_VALID_MAX 1.0
#define AQ_ASSERT( q ) assert( ((q) == AQ_UNKNOWN) || ((q) == AQ_INVALID) || ((q) >=AQ_VALID_MIN && (q) <= AQ_VALID_MAX) )

#endif // COMMON_H__INCLUDED
