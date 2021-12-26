/////////////////////////////////////////////////////////////////////////////
// DSurfaceFitting.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2010-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#ifndef DSURFACEFITTING_H__INCLUDED
#define DSURFACEFITTING_H__INCLUDED

class Metric3dContext;
class MeshPoint3d;
class SurfaceParametric;

/**
 * This class implements methods for fitting parametrical surfaces to discrete data.
 */
class DSurfaceFitting
{
public:
	/// Fit local parameterical surface to surface mesh sub-region (with tolerance in local metric)
	static SurfaceParametric* fitSurface(Metric3dContext& mc, MeshPoint3d* point, double tolerance);
};

#endif // !defined(DSURFACEFITTING_H__INCLUDED)
