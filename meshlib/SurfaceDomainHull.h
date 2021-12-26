/////////////////////////////////////////////////////////////////////////////
// SurfaceDomainHull.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2013-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#if !defined(SURFACEDOMAINHULL_H__INCLUDED)
#define SURFACEDOMAINHULL_H__INCLUDED

#include <deque>
#include <vector>

#include "SurfaceDomain.h"
#include "DataList.h"
#include "DPoint.h"

class SurfaceDomainHull : public SurfaceDomain
{
public:
	SurfaceDomainHull() {}
	SurfaceDomainHull(const DataVector<DPoint2d> & hull);
public:
	virtual double getInsideQuality(const DPoint3d& point, const DPoint2d& param) const override;
public:
	void addPoint(const DPoint2d& pt);
	bool createHull();
	DPoint2d calculateHullMiddle( Metric3dContext & mc, const SurfaceParametric * surface ) const;
	bool createHullDist( Metric3dContext & mc, const SurfaceParametric * surface );
	void draw(MeshViewSet* set = nullptr, const SurfaceParametric* surface = nullptr, int id = -2) const override;
	/// Store XML description to stream
	virtual ostream& storeXML(ostream& os, const string& prefix = "") const override;
	virtual DPoint2d getMiddleParam() const override;
//private:
//	struct HullDist {
//		HullDist( double a = 0.0, double dl = 0.0, double dr = 0.0 ) : angle(a), dlen(dl), dratio(dr) {}
//		double angle;
//		double dlen;
//		double dratio;
//	};
private:
	std::deque<DPoint2d> m_points;
	std::vector<DPoint2d> m_hull;
	DPoint2d m_middle;
	DataVector<double> m_hdist;
	int m_hoffset;
	double m_slen2, m_tlen2;
};

#endif // SURFACEDOMAINHULL_H__INCLUDED