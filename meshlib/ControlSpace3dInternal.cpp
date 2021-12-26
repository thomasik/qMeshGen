// ControlSpace3dInternal.cpp: implementation of the ControlSpace3dInternal class.
//
//////////////////////////////////////////////////////////////////////

#include "ControlSpace3dInternal.h"
#include "DTriangle.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

ControlSpace3dInternal::ControlSpace3dInternal()
{ }

ControlSpace3dInternal::~ControlSpace3dInternal()
{ }

/// Get sizing info (matrix mode) at the given point
ControlDataMatrix3d ControlSpace3dInternal::getMetricAtPoint(const DPoint3d& pt) const
{
	static const double X0 = -0.4;
	static const double dmax = 0.2; 
	static const double dmin = dmax * 0.01;
	static const double ro = 2.0;
	static const double log_2 = log(2.0);
	//static const double log_ro = log(ro);

	double d = 0.2;
	double dist = (pt.x < X0) ? abs(pt.y) : sqrt((pt.x - X0)*(pt.x - X0) + pt.y * pt.y);
	if( dist < 2*dmin ) d = dmin;
	else{
		double a = dist / dmin - 1.0;
		double x = pow(ro, 1.0 + log(a) / log_2);
		d = x * dmin;
		if(d > dmax) d = dmax;
	}

	return ControlDataMatrix3d(d, d, d, 0.0, 0.0, 0.0);
}
