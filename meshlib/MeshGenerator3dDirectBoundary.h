/////////////////////////////////////////////////////////////////////////////
// MeshGenerator3dDirectBoundary.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2008-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#if !defined(MESHGENERATOR3DDIRECTBOUNDARY_H__INCLUDED)
#define MESHGENERATOR3DDIRECTBOUNDARY_H__INCLUDED

#include "DPoint.h"
#include "DataVector.h"

class Metric3dContext;
class MeshContainer3d;
class MeshPoint3d;
class MeshFace;
class MeshDomainVolume;
class ControlSpace3d;
class MeshContainer3dSurface;

/**
 * This class gathers several procedures responsible for discretizing volumes 
 * into 3D tetrahedral meshes (direct 3d triangulation of surface mesh).
 */
class MeshGenerator3dDirectBoundary  
{
public:
	/// Generate boundary constrained tetrahedral mesh
	//static MeshContainer3d* createBoundaryConstrainedMesh(Metric3dContext& mc, 
	//	const DataVector<MeshFace*> & bfaces, const DataVector<MeshPoint3d*> & bpoints, 
	//	MeshDomainVolume* mdv);
	/// Check "star-shape" property of boundary mesh - by calculating minimum volume
	static double minStarTetrahedraVolume(Metric3dContext& mc, 
		MeshContainer3dSurface* surface_mesh, 
		MeshDomainVolume* mdv, DPoint3d& center);
	/// Create simple boundary constrained mesh for convex (star shaped) surface mesh
	static MeshContainer3d* createSimpleConvexBoundaryConstrainedMesh(Metric3dContext& mc,
		MeshContainer3dSurface* surface_mesh, MeshDomainVolume* mdv, const DPoint3d& center);
	/// Create boundary constrained mesh using frontal approach
	static MeshContainer3d* createFrontalBoundaryConstrainedMesh(Metric3dContext& mc, 
		MeshContainer3dSurface* surface_mesh, const MeshDomainVolume* mdv);
protected:
	/// calculate inertial center from the surface mesh (weighted)
	static DPoint3d calculateInertialCenterFromSurface(MeshContainer3dSurface* surface_mesh);
};

#endif // !defined(MESHGENERATOR3DDIRECTBOUNDARY_H__INCLUDED)
