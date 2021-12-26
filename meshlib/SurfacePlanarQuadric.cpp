// SurfacePlanarQuadric.cpp: implementation of the SurfacePlanarQuadric class.
// Tomasz Jurczyk, 2013-
// Generation of unstructured meshes
//
//////////////////////////////////////////////////////////////////////

#include "SurfacePlanarQuadric.h"
#include "DEquation.h"
//#include "Curve2dParametric.h"
//#include "DMatrix.h"
//#include "DPlane.h"

const DPoint3d SurfacePlanarQuadric::getPoint(const DPoint2d& param) const
{
	return m_pquadric.getPoint(param);
}

const DPoint2d SurfacePlanarQuadric::getParameters(const DPoint3d& point) const
{
	return m_pquadric.plane.projectToPlane(point);
	//const DPoint2d start_point = m_pquadric.plane.projectToPlane(point);
	//return SurfaceParametric::getParametersNear(point, start_point);
}

const DVector3d SurfacePlanarQuadric::getDerivative(int deriv, const DPoint2d& param) const
{
	return m_pquadric.getDerivative( deriv, param );
}

ostream& operator<<(ostream& os, const SurfacePlanarQuadric *surf)
{
	if(surf->m_valid){
		os << surf->m_pquadric;
	}else{
		os << "invalid";
	}
	return os;
}

istream& operator>>(istream& is, SurfacePlanarQuadric *surf)
{
	is >> surf->m_pquadric;
	surf->m_valid = is.good();
	return is;
}

/// Store XML description to stream
ostream& SurfacePlanarQuadric::storeXML(ostream& os, const string& prefix) const
{
	os << prefix << "<planarquadric>\n";
	os << prefix << "\t<pt0> " << m_pquadric.plane.p0 << " </pt0>\n";
	os << prefix << "\t<e0> "  << m_pquadric.plane.e0 << " </e0>\n";
	os << prefix << "\t<e1> "  << m_pquadric.plane.e1 << " </e1>\n";
	os << prefix << "\t<vq> "  << m_pquadric.vq << " </vq>\n";
	return os << prefix << "</planarquadric>\n";
}

/// invert the orientation of the surface (change diretion of normal vector) if possible
bool SurfacePlanarQuadric::invertOrientation()
{
	return m_pquadric.invertOrientation();
}

DOrientedBox SurfacePlanarQuadric::getOrientedBox() const
{
	return m_pquadric.getOrientedBox();
}

DOrientedBox SurfacePlanarQuadric::getOrientedBox( const DataVector<DPoint3d>& points ) const
{
	return m_pquadric.getOrientedBox( points );
}

DOrientedBox SurfacePlanarQuadric::getOrientedBoxOpt( const DataVector<DPoint3d>& points ) const
{
	return m_pquadric.getOrientedBoxOpt( points );
}
