/////////////////////////////////////////////////////////////////////////////
// SurfaceDomainFaces.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2014-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#if !defined(SURFACEDOMAINFACES_H__INCLUDED)
#define SURFACEDOMAINFACES_H__INCLUDED

#include "SurfaceDomain.h"
#include "DataList.h"
#include "DPoint.h"
#include "DTriangle.h"
#include "DataVector.h"

#include <vector>

class SurfaceDomainFaces : public SurfaceDomain
{
public:
	//SurfaceDomainFaces() {}
	SurfaceDomainFaces( const DataVector<DPoint3d> & local_points, 	const DataVector<ITriangle>& local_faces, 
		const DataVector<double>& f_dist2_tol, const DataVector<double>& p_approx_err);
public:
	virtual double getInsideQuality(const DPoint3d& point, const DPoint2d& param) const override;
public:
	void draw(MeshViewSet* set = nullptr, const SurfaceParametric* surface = nullptr, int id = -2) const override;
	/// Store XML description to stream
	virtual ostream& storeXML(ostream& os, const string& prefix = "") const override;
private:
	std::vector<DPoint3d> m_points;
	std::vector<double> m_papprox_err;
	std::vector<double> m_fdist2_tol; // permitted distance (squared) from a face for a point to be valid
	std::vector<ITriangle> m_tri_indices;
};

#endif // SURFACEDOMAINFACES_H__INCLUDED