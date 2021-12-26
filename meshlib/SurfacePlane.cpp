// SurfacePlane.cpp: implementation of the SurfacePlane class.
// Tomasz Jurczyk, 2006-
// Generation of unstructured meshes
//
//////////////////////////////////////////////////////////////////////

#include "SurfacePlane.h"
#include "DEquation.h"
#include "Curve2dParametric.h"
#include "DMatrix.h"
#include "DPlane.h"
#include "DataVector.h"
#include "MeshFace.h"
#include "MeshPoint3d.h"

/// Standard constructor

SurfacePlane::SurfacePlane(const DPoint3d& pnt0, const DPoint3d& pnt1, const DPoint3d& pnt2) 
	: SurfaceParametric(), m_p0(pnt0), m_v01(pnt1-pnt0), m_v02(pnt2-pnt0), 
		m_vn(m_v01.crossProduct(m_v02).normalized())
{
	double l1 = m_v01.length2();
	double l2 = m_v02.length2();
	m_dratio = std::min(l1, l2) / std::max(l1, l2);

	m_valid = true;
}

SurfacePlane::SurfacePlane(const DPlane& plane)
	: SurfaceParametric(), m_p0(plane.p0), m_v01(plane.e0), m_v02(plane.e1), m_vn(plane.vn)
{
	double l1 = m_v01.length2();
	double l2 = m_v02.length2();
	m_dratio = std::min(l1, l2) / std::max(l1, l2);

	m_valid = true;
}

SurfacePlane::SurfacePlane(const DPoint3d& pnt0, const DVector3d& v01, const DVector3d& v02) 
	: SurfaceParametric(), m_p0(pnt0), m_v01(v01), m_v02(v02),
		m_vn(m_v01.crossProduct(m_v02).normalized())
{
	double l1 = m_v01.length2();
	double l2 = m_v02.length2();
	m_dratio = std::min(l1, l2) / std::max(l1, l2);

	m_valid = true;
}

SurfacePlane::SurfacePlane(const DPoint3d& pnt0, const DVector3d& vn) 
	: SurfaceParametric(), m_p0(pnt0), m_vn(vn)
{
	vn.orthogonalVectors( m_v01, m_v02 );
	double l1 = m_v01.length2();
	double l2 = m_v02.length2();
	m_dratio = std::min(l1, l2) / std::max(l1, l2);

	m_valid = true;
}

/// Standard constructor
SurfacePlane::SurfacePlane(const DVector3d& v01, const DVector3d& v02) 
	: SurfaceParametric(), m_p0(DPoint3d::zero), m_v01(v01), m_v02(v02),
		m_vn(m_v01.crossProduct(m_v02).normalized())
{
	double l1 = m_v01.length2();
	double l2 = m_v02.length2();
	m_dratio = std::min(l1, l2) / std::max(l1, l2);

	m_valid = true;
}
/// Standard constructor

const DPoint3d SurfacePlane::getPoint(const DPoint2d& param) const
{
	return m_p0 + (m_v01 * param.x) + (m_v02 * param.y);
}

double SurfacePlane::getShapeParameters(const DPoint3d& point, Curve2dConstPtr shape, double near_t, double min_t, double max_t) const
{
	double t = shape->getParameterInRange(getParameters(point), near_t, min_t, max_t);
	assert((t >= (min_t-SMALL_NUMBER)) && (t <= (max_t+SMALL_NUMBER)));
	return std::max(min_t, std::min(max_t, t));
}

double SurfacePlane::getSegmentParameters(const DPoint3d& point, const DPoint2d& pt0, const DPoint2d& pt1, double /* near_t */, double min_t, double max_t) const
{
	const DPoint3d pt0_3d = getPoint(pt0);
	const DPoint3d pt1_3d = getPoint(pt1);
	double dist1 = pt0_3d.distance(pt1_3d);
	double dist0 = pt0_3d.distance(point);
	double t = dist0 / dist1;
	return std::max(min_t, std::min(max_t, t));
}

/// Returns the parameters of the surface for the given point (numerical approximation) and signed distance (according to normal vecor)
bool SurfacePlane::getParametersAndSignedDistance(const DPoint3d& point, DPoint2d& param, double & z ) const 
{
	DMatrix3d A;
	A.setColumn(0, m_v01);
	A.setColumn(1, m_v02);
	A.setColumn(2, m_vn);
	DVector3d res;
	if(!A.solve(point - m_p0, res)) {
		assert(false);
		return false;
	}
	
	param.x = res.x;
	param.y = res.y;
	z = res.z;

	return true;
}

const DPoint2d SurfacePlane::getParameters(const DPoint3d& point) const
{
	DMatrix3d A;
	A.setColumn(0, m_v01);
	A.setColumn(1, m_v02);
	A.setColumn(2, m_vn);
	DVector3d res;
	if(A.solve(point - m_p0, res)){
		return DPoint2d(res.x, res.y);
	}else{
		// old method
		double Wxy = m_v01.x * m_v02.y - m_v01.y * m_v02.x;
		double Wxz = m_v01.x * m_v02.z - m_v01.z * m_v02.x;
		double Wyz = m_v01.y * m_v02.z - m_v01.z * m_v02.y;
		double Ws, Wt;
		DPoint2d result;

		DVector3d point_diff = point - m_p0;

		if(abs(Wxy) > abs(Wxz)){
			if(abs(Wxy) > abs(Wyz)){
				//xy
				Ws = point_diff.x * m_v02.y - point_diff.y * m_v02.x;
				Wt = m_v01.x * point_diff.y - m_v01.y * point_diff.x;
				result.x = Ws / Wxy;
				result.y = Wt / Wxy;
			}else{
				//yz
				Ws = point_diff.y * m_v02.z - point_diff.z * m_v02.y;
				Wt = m_v01.y * point_diff.z - m_v01.z * point_diff.y;
				result.x = Ws / Wyz;
				result.y = Wt / Wyz;
			}
		}else{
			if(abs(Wxz) > abs(Wyz)){
				//xz
				Ws = point_diff.x * m_v02.z- point_diff.z* m_v02.x;
				Wt = m_v01.x * point_diff.z- m_v01.z* point_diff.x;
				result.x = Ws / Wxz;
				result.y = Wt / Wxz;
			}else{
				//yz
				Ws = point_diff.y * m_v02.z - point_diff.z * m_v02.y;
				Wt = m_v01.y * point_diff.z - m_v01.z * point_diff.y;
				result.x = Ws / Wyz;
				result.y = Wt / Wyz;
			}
		}
		return result;
	}
}

const DVector3d SurfacePlane::getDerivative(int deriv, const DPoint2d& /* param */) const
{
	switch(deriv){
	case DEquation::deriv_ds:
		return m_v01;
	case DEquation::deriv_dt:
		return m_v02;
	case DEquation::deriv_dss:
	case DEquation::deriv_dst:
	case DEquation::deriv_dtt:
		return DVector3d(0.0, 0.0, 0.0);
	case DEquation::deriv_dsss:
	case DEquation::deriv_dsst:
	case DEquation::deriv_dstt:
	case DEquation::deriv_dttt:
		return DVector3d(0.0, 0.0, 0.0);
	default:
		assert(false);
	}
	return DVector3d(0.0, 0.0, 0.0);
}

double SurfacePlane::segmentLength(const DPoint2d & param_0, const DPoint2d & param_1) const
{
	return getPoint(param_0).distance(getPoint(param_1));
}

/*
bool SurfacePlane::setData(const DPoint3d& pnt0, const DPoint3d& pnt1, const DPoint3d& pnt2)
{
	m_p0 = pnt0;
	m_v01 = pnt1 - pnt0;
	m_v02 = pnt2 - pnt0;
	// TODO : check for linearity of these points

	double l1 = m_v01.length2();
	double l2 = m_v02.length2();
	m_dratio = std::min(l1, l2) / std::max(l1, l2);

	return m_valid = true;
}
*/
ostream& operator<<(ostream& os, const SurfacePlane *surf)
{
	if(surf->m_valid){
		os << surf->m_p0 << " " << surf->m_v01 << " " << surf->m_v02;
	}else{
		os << "invalid";
	}
	return os;
}

istream& operator>>(istream& is, SurfacePlane *surf)
{
	is >> ws >> surf->m_p0 >> ws >> surf->m_v01 >> ws >> surf->m_v02;

	double l1 = surf->m_v01.length2();
	double l2 = surf->m_v02.length2();
	surf->m_dratio = std::min(l1, l2) / std::max(l1, l2);

	surf->m_valid = true;
	return is;
}
/*
bool SurfacePlane::getSegmentCrossingPointParam(const DPoint3d& pnt0, const DPoint3d& pnt1, 
												 DPoint2d& pnt) const
{
	const DVector3d pt_r = pnt1-pnt0;
	const DVector3d pt_s = m_v01;
	const DVector3d pt_t = m_v02;
	const DVector3d pt_x = m_p0-pnt0;

	double w = pt_r.det(pt_s, pt_t);
	if( abs(w) < mesh_data.relative_small_number) return false;

	double wr = pt_x.det(pt_s, pt_t);
	double r = wr / w;
	//if(r < -mesh_data.relative_small_number || r > 1+mesh_data.relative_small_number) return false;

	//double ws = pt_r.det(pt_x, pt_t);
	//double s = ws / w;
	//if(s < -mesh_data.relative_small_number || s > 1+mesh_data.relative_small_number) return false;

	//double wt = pt_r.det(pt_s, pt_x);
	//double t = wt / w;
	//if(t < -mesh_data.relative_small_number || t > 1+mesh_data.relative_small_number) return false;

	//if((s+t) < -mesh_data.relative_small_number || (s+t) > 1+mesh_data.relative_small_number) return false;

	pnt = getParameters(pnt0*(1-r) + pnt1*r);
	return true;
}
*/

/// Returns the approximation of a segment on surfaces via polyline (array of points)
void SurfacePlane::getPolyLine(DataVector<DPoint3d> & polyline, const DPoint2d& a, const DPoint2d& b) const
{
	polyline.add(getPoint(a));
	polyline.add(getPoint(b));
}

/// Returns the approximation of a segment on surfaces via polyline (array of points)
void SurfacePlane::getPolyLine(DataVector<DPoint3d> & polyline, Curve2dConstPtr shape, double t0, double t1) const
{
	DataVector<double> poly_double;
	shape->getPolyLineInRange(t0, t1, poly_double);
	for(size_t i = 0; i < poly_double.countInt(); i++){
		polyline.add(getPoint(shape->getPoint(poly_double[i])));
	}
}

/// Store XML description to stream
ostream& SurfacePlane::storeXML(ostream& os, const string& prefix) const
{
	os << prefix << "<plane>\n";
	os << prefix << "\t<pt0> " << m_p0 << " </pt0>\n";
	os << prefix << "\t<e0> "  << m_v01 << " </e0>\n";
	os << prefix << "\t<e1> "  << m_v02 << " </e1>\n";
	return os << prefix << "</plane>\n";
}

/// invert the orientation of the surface (change diretion of normal vector) if possible
bool SurfacePlane::invertOrientation()
{
	m_v01 = -m_v01;
	m_vn = -m_vn;
	return true;
}

DOrientedBox SurfacePlane::getOrientedBox() const
{
	return DOrientedBox( m_p0, m_v01, m_v02 );
}

DOrientedBox SurfacePlane::getOrientedBox( const DataVector<DPoint3d>& points ) const
{
	DOrientedBox obox( m_p0, m_v01, m_v02 );
	for(size_t i = 0; i < points.countInt(); i++) {
		obox.addPoint( points[i] );
	}
	return obox;
}

DOrientedBox SurfacePlane::getOrientedBoxOpt( const DataVector<DPoint3d>& points ) const
{
	const DVector3d v1 = m_v01.normalized();
	const DVector3d v2 = m_v02.normalized();
	DOrientedBox obox( m_p0, (v2+v1)*0.5, (v2-v1)*0.5 );
	for(size_t i = 0; i < points.countInt(); i++) {
		obox.addPoint( points[i] );
	}
	return obox;
}

std::shared_ptr<SurfaceParametric> SurfacePlane::adjustedForFaces( const DataVector< MeshFace* > & mfaces,
		const DataVector< MeshPoint3d* > & mpoints, MeshFace* central_face ) const
{
	assert( mfaces.notEmpty() );
	assert( mpoints.notEmpty() );
	if( mfaces.empty() || mpoints.empty() ) return nullptr;

	DVector3d ave_normal = mfaces[0]->getBaseNormal();
	// ... or least-squares ?
	for(size_t i = 1; i < mfaces.countInt(); i++ )
		ave_normal += mfaces[i]->getBaseNormal();
	ave_normal.normalize();

	DPoint3d ave_middle;
	size_t mpct = mpoints.countInt();
	double pf = 1.0 / mpct;
	for(size_t i = 0; i < mpct; i++) {
		ave_middle.add( mpoints[i]->getCoordinates(), pf );
	}

	auto surf = std::make_shared<SurfacePlane>( ave_middle, ave_normal );
	if( central_face != nullptr ) surf->setOrientationLikeFace( central_face );

	return surf;
}
