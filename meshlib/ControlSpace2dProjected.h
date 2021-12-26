// ControlSpace2dProjected.h: interface for the ControlSpace2dProjected class.
//
//////////////////////////////////////////////////////////////////////

#pragma once

#if !defined(CONTROLSPACEPROJECTED_H__INCLUDED)
#define CONTROLSPACEPROJECTED_H__INCLUDED

#include "ControlSpace2d.h"
#include "ControlSpace3d.h"
#include "MeshData.h"

class ControlSpace2dProjected : public ControlSpace2d  
{
public:
	/// Standard constructor
	ControlSpace2dProjected(CS3dPtr cs3d, SurfaceConstPtr surface = nullptr)
		: ControlSpace2d(surface), m_cs3d(cs3d) {}
public:
	/// Returns the type of element (should be reimplemented in derived classes)
	virtual int getType() const { return MeshData::CONTROL_PROJECTED; }
	/// Get sizing info (matrix mode) at the given point
	virtual ControlDataMatrix2d getMetricAtPoint(const DPoint2d& pt2d) const { 
		return base_surface->projectTransformationTensor(
			pt2d, m_cs3d->getMetricAtPoint( base_surface->getPoint( pt2d ) ) );
		}
private:
	CS3dPtr m_cs3d;
};

#endif // !defined(CONTROLSPACEPROJECTED_H__INCLUDED)
