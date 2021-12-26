/////////////////////////////////////////////////////////////////////////////
// DSurfaceFitting.cpp
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2010-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#include "common.h"

#include "DSurfaceFitting.h"
#include "Metric3dContext.h"
#include "MeshPoint3d.h"
#include "SurfaceParametric.h"

#include "MeshViewSet.h"

/// Fit local parameterical surface to surface mesh sub-region (with tolerance in local metric)
SurfaceParametric* DSurfaceFitting::fitSurface(Metric3dContext& mc, MeshPoint3d* point, double tolerance)
{
	return NULL;
}