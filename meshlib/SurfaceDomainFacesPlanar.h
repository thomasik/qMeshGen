/////////////////////////////////////////////////////////////////////////////
// SurfaceDomainFacesPlanar.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2014-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#if !defined(SURFACEDOMAINFACESPLANAR_H__INCLUDED)
#define SURFACEDOMAINFACESPLANAR_H__INCLUDED

#include "SurfaceDomain.h"
#include "DataList.h"
#include "DPoint.h"
#include "DTriangle.h"
#include "DataVector.h"

#include <vector>

class SurfaceDomainFacesPlanar : public SurfaceDomain
{
public:
	//SurfaceDomainFacesPlanar() {}
	SurfaceDomainFacesPlanar( 
		const DataVector<DPoint2d> & local_points, 
		const DataVector<ITriangle>& local_faces,
		const DataVector<double>& p_approx_quality);
public:
	virtual double getInsideQuality(const DPoint3d& point, const DPoint2d& param ) const override;
public:
	int invCount() const { return inv_count; }
	void draw(MeshViewSet* set = nullptr, const SurfaceParametric* surface = nullptr, int id = -2) const override;
	/// Store XML description to stream
	virtual ostream& storeXML(ostream& os, const string& prefix = "") const override;
	virtual DPoint2d getMiddleParam() const override;
private:
	int segmentKey( int a, int b, int n ) { return std::min( a, b) * n + std::max( a, b ); }
private:
	std::vector<DPoint2d> m_points;
	std::vector<double> m_papprox_quality;
	std::vector<ITriangle> m_tri_indices;
	double m_eps;
	int inv_count;
};

#endif // SURFACEDOMAINFACESPLANAR_H__INCLUDED