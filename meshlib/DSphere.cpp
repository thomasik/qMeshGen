#include "DSphere.h"
#include "DMatrix.h"
#include "SurfaceCurvature.h"
#include "DEquation.h"

#include "MeshViewSet.h"

double DCircle::implicitValue( const DPoint2d& pt ) const
{
	return sqr( pt.x - m_center.x ) + sqr( pt.y - m_center.y ) - m_radius * m_radius;
}

DPoint2d DCircle::projectToCurve(const DPoint2d& pt) const
{
	DVector2d dv = pt - m_center;
	if( dv.isZero() ) return m_center + DVector2d(m_radius, 0.0);

	return m_center + dv.normalized( m_radius );
}

DPoint2d DCircle::getCirclePoint( const DPoint2d& mid, double r, double t )
{
	double angle = 2 * PI * t;
	return DPoint2d( mid.x + r * cos(angle), mid.y + r * sin(angle) );
}

DVector2d DCircle::getCircleDerivative( double r, double t )
{
	double angle = 2 * PI * t;
	double a = 2 * PI * r;
	return DVector2d( -a * sin(angle), a * cos(angle) );
}

DVector2d DCircle::getCircleDerivative2( double r, double t )
{
	double angle = 2 * PI * t;
	double a = 4 * PI * PI * r;
	return DVector2d( -a * cos(angle), -a * sin(angle) );
}

DVector2d DCircle::getCircleDerivative3( double r, double t )
{
	double angle = 2 * PI * t;
	double a = 8 * PI * PI * PI * r;
	return DVector2d( a * sin(angle), -a * cos(angle) );
}

double DCircle::getCircleParam( const DPoint2d& cmid, double cr, const DPoint2d& pt, double ts )
{
	if(cr <= 0.0) return 0.0;
	double r = cmid.distance(pt);
	double cos_dx = (pt.x - cmid.x) / r;
	if(cos_dx > 1.0) {
		cos_dx = 1.0;
	}else if(cos_dx < -1.0){
		cos_dx = -1.0;
	}
	double sin_dy = (pt.y - cmid.y) / r;
	if(sin_dy > 1.0) {
		sin_dy = 1.0;
	}else if(sin_dy < -1.0){
		sin_dy = -1.0;
	}
	double t = acos(cos_dx) / (2*PI);	// angle -> [0,PI], t -> [0-0.5]
	if(asin(sin_dy) < 0) t = 1 - t;		// angle -> [0,2PI], t -> [0-1]

	// ts - choose t nearest
	if( t > ts )
		while( (t-ts) > 0.5 ) t -= 1.0;
	else
		while( (ts-t) > 0.5 ) t += 1.0;

	return t;
}

double DCircle::getParam( const DPoint2d& pt, double ts ) const {
	return getCircleParam( m_center, m_radius, pt, ts );
}

double DSphere::implicitValue( const DPoint3d& pt ) const
{
	return sqr( pt.x - m_center.x ) + sqr( pt.y - m_center.y ) + sqr( pt.z - m_center.z ) - 
		m_radius * m_radius;
}

double DEllipsoid::implicitValue(const DPoint3d& pt) const
{
	return sqr(pt.x - m_center.x)/sqr(m_rx) + sqr(pt.y - m_center.y)/sqr(m_ry)
		+ sqr(pt.z - m_center.z)/sqr(m_rz) - 1.0;
}

double DEllipsoid::getRadius(int i) const {
	switch (i) {
	case 0: return m_rx;
	case 1: return m_ry;
	case 2: return m_rz;
	default:
		assert(false);
	}
	return 0;
}

double & DEllipsoid::getRadius(int i) {
	switch (i) {
	case 0: return m_rx;
	case 1: return m_ry;
	case 2: return m_rz;
	default:
		assert(false);
	}
	return m_rx;
}

DPoint3d DSphere::projectToSurface(const DPoint3d& pt) const
{
	DVector3d dv = pt - m_center;
	double r = getRadius();
	if( dv.isZero() ) return m_center + DVector3d(r, 0.0, 0.0);

	return m_center + dv.normalized( r );
}

const DPoint3d DSphere::getPoint(const DPoint2d& param) const
{
	return m_center + getNormalVector(param) * m_radius;
}

const DPoint3d DEllipsoid::getPoint(const DPoint2d& param) const
{
	double a1 = 2 * PI * param.x;
	double a2 = PI * param.y;
	double sin_a2 = sin(a2);

	DVector3d dv = DVector3d::zero;

	switch (m_pole) {
	case POLE_OX:
		dv = DVector3d(
			-cos(a2),
			cos(a1) * sin_a2,
			sin(a1) * sin_a2);
		break;
	case POLE_OY:
		dv = DVector3d(
			cos(a1) * sin_a2,
			cos(a2),   // ??? TO_CHECK
			sin(a1) * sin_a2);
		break;
	case POLE_OZ:
		dv =  DVector3d(
			cos(a1) * sin_a2,
			sin(a1) * sin_a2,
			-cos(a2));
		break;
	}
	return DPoint3d(
		m_center.x + dv.x * m_rx, 
		m_center.y + dv.y * m_ry,
		m_center.z + dv.z * m_rz);
}

/// Returns the principle curvature for this sphere and parameters [s,t]
const SurfaceCurvature DSphere::getCurvature(const DPoint2d& /* pt */, double & g_ratio) const
{ 
	g_ratio = 1.0; 
	double c = 1.0 / m_radius;
	return SurfaceCurvature( c, c, 0.0);
}

const DPoint2d DSphere::getParam(const DPoint3d& point) const
{
	DVector3d dvn = (point - m_center).normalized();
	if( m_radius < 0.0 ) dvn = -dvn;

	double a2 = 0.0;
	double param_x = 0.0;

	switch( m_pole ) {
	case POLE_OX:
		{
			a2 = acos( -dvn.x );
			double r_yz = sqrt( dvn.y * dvn.y + dvn.z * dvn.z );
			param_x = DCircle::getCircleParam( DPoint2d::zero, r_yz, DPoint2d( dvn.y, dvn.z ), 0.5 );
			break;
		}
	case POLE_OY:
		{
			a2 = acos( dvn.y );
			double r_xz = sqrt( dvn.x * dvn.x + dvn.z * dvn.z );
			param_x = DCircle::getCircleParam( DPoint2d::zero, r_xz, DPoint2d( dvn.x, dvn.z ), 0.5 );
			break;
		}
	case POLE_OZ:
		{
			a2 = acos( -dvn.z );
			double r_xy = sqrt( dvn.x * dvn.x + dvn.y * dvn.y );
			param_x = DCircle::getCircleParam( DPoint2d::zero, r_xy, DPoint2d( dvn.x, dvn.y ), 0.5 );
			break;
		}
	default:
		assert(false);
	}

	return DPoint2d( param_x, a2 / PI );
}

bool DSphere::getParamAndSignedDistance(const DPoint3d& point, DPoint2d& param, double & z) const
{
	param = getParam(point);

	double dist = point.distance(m_center);
	z = ( m_radius > 0.0 ) ? ( dist - m_radius ) : (-dist - m_radius);

	return true;
}

const DVector3d DSphere::getNormalVector(const DPoint2d& param) const
{
	double a1 = 2*PI * param.x;
	double a2 = PI * param.y;
	double sin_a2 = sin(a2);
	
	switch( m_pole ) {
	case POLE_OX:
		return DVector3d(
			-cos(a2),
			cos(a1) * sin_a2, 
			sin(a1) * sin_a2);
	case POLE_OY:
		return DVector3d(
			cos(a1) * sin_a2, 
			cos(a2),   // ??? TO_CHECK
			sin(a1) * sin_a2);
	case POLE_OZ:
		return DVector3d(
			cos(a1) * sin_a2, 
			sin(a1) * sin_a2, 
			-cos(a2));
	}

	assert(false);
	return DVector3d::zero;
}

const DVector3d DEllipsoid::getNormalVector(const DPoint2d& param) const
{
	return (getPoint(param) - m_center).normalized();
}

const DVector3d DSphere::getNormalVectorDerivative(int deriv, const DPoint2d& param) const
{
	double ax = 2*PI * param.x;
	double ay = PI * param.y;

	switch (deriv) {
	case DEquation::deriv_ds:
	{
		double sin_ay = sin(ay);
		switch (m_pole) {
		case POLE_OX:
			return DVector3d( 0.0, -2*PI*sin(ax)*sin_ay, 2*PI*cos(ax)*sin_ay );
		case POLE_OY:
			return DVector3d( -2*PI*sin(ax)*sin_ay, 0.0, 2*PI*cos(ax)*sin_ay );
		case POLE_OZ:
			return DVector3d( -2*PI*sin(ax)*sin_ay, 2*PI*cos(ax)*sin_ay, 0.0 );
		}
	}
	case DEquation::deriv_dt:
	{
		double sin_ay = sin(ay);
		double cos_ay = cos(ay);
		switch (m_pole) {
		case POLE_OX:
			return DVector3d( PI*sin_ay, PI*cos(ax)*cos_ay, PI*sin(ax)*cos_ay );
		case POLE_OY:
			return DVector3d( PI*cos(ax)*cos_ay, -PI*sin_ay, PI*sin(ax)*cos_ay );
		case POLE_OZ:
			return DVector3d( PI*cos(ax)*cos_ay, PI*sin(ax)*cos_ay, PI*sin_ay );
		}
	}
	case DEquation::deriv_dss:
	{
		double sin_ay = sin(ay);
		switch (m_pole) {
		case POLE_OX:
			return DVector3d( 0.0, -4*PI*PI*cos(ax)*sin_ay, -4*PI*PI*sin(ax)*sin_ay );
		case POLE_OY:
			return DVector3d( -4*PI*PI*cos(ax)*sin_ay, 0.0, -4*PI*PI*sin(ax)*sin_ay );
		case POLE_OZ:
			return DVector3d( -4*PI*PI*cos(ax)*sin_ay, -4*PI*PI*sin(ax)*sin_ay, 0.0 );
		}
	}
	case DEquation::deriv_dst:
	{
		double sin_ay = sin(ay);
		double cos_ay = cos(ay);
		switch (m_pole) {
		case POLE_OX:
			return DVector3d( 0.0, -2*PI*PI*sin(ax)*cos_ay, 2*PI*PI*cos(ax)*cos_ay );
		case POLE_OY:
			return DVector3d( -2*PI*PI*sin(ax)*cos_ay, 0.0, 2*PI*PI*cos(ax)*cos_ay );
		case POLE_OZ:
			return DVector3d( -2*PI*PI*sin(ax)*cos_ay, 2*PI*PI*cos(ax)*cos_ay, 0.0 );
		}
	}
	case DEquation::deriv_dtt:
	{
		double sin_ay = sin(ay);
		double cos_ay = cos(ay);
		switch (m_pole) {
		case POLE_OX:
			return DVector3d( PI*PI*cos_ay, -PI*PI*cos(ax)*sin_ay, -PI*PI*sin(ax)*sin_ay );
		case POLE_OY:
			return DVector3d( -PI*PI*cos(ax)*sin_ay, -PI*PI*cos_ay, -PI*PI*sin(ax)*sin_ay );
		case POLE_OZ:
			return DVector3d( -PI*PI*cos(ax)*sin_ay, -PI*PI*sin(ax)*sin_ay, PI*PI*cos_ay );
		}
	}
	case DEquation::deriv_dsss:
	case DEquation::deriv_dsst:
	case DEquation::deriv_dstt:
	case DEquation::deriv_dttt:
	{
		assert(false);
		return DVector3d::zero;
	}
	default:
		assert(false);
	}
	return DVector3d::zero;
}

const DVector3d DSphere::getNormalVectorForPoint3d(const DPoint3d& pt) const
{
	return (pt - m_center).normalized();
}

const DVector3d DSphere::getDerivative(int deriv, const DPoint2d& param) const
{
	double a1 = 2*PI * param.x;
	double a2 = PI * param.y;

	switch(deriv){
	case DEquation::deriv_ds:
		{
			double f = m_radius * sin( a2 ) * 2 * PI;
			double d1 = -f * sin( a1 );
			double d2 =  f * cos( a1 );
			switch( m_pole ) {
			case POLE_OX:
				return DVector3d( 0.0, d1, d2 );
			case POLE_OY:
				return DVector3d( d1, 0.0, d2 );
			case POLE_OZ:
				return DVector3d( d1, d2, 0.0 );
			}
		}
	case DEquation::deriv_dt:
		{
			double f = m_radius * cos( a2 ) * PI;
			double d1 = f * cos( a1 );
			double d2 = f * sin( a1 );
			double dp = -m_radius * sin(a2) * PI;
			switch( m_pole ) {
			case POLE_OX:
				return DVector3d( -dp, d1, d2 );
			case POLE_OY:
				return DVector3d( d1, dp, d2 );
			case POLE_OZ:
				return DVector3d( d1, d2, -dp );
			}
		}
	case DEquation::deriv_dss:
	case DEquation::deriv_dst:
	case DEquation::deriv_dtt:
	case DEquation::deriv_dsss:
	case DEquation::deriv_dsst:
	case DEquation::deriv_dstt:
	case DEquation::deriv_dttt:
		{
			assert( false ); // needs to be written if required
			return DVector3d::zero;

		}
	default:
		assert(false);
	}
	return DVector3d::zero;
}

bool DSphere::invertOrientation()
{
	m_radius = -m_radius;
	return true;
}

void DSphere::test() 
{
	const double R1 = 5.5;
	DSphere sphere1( DPoint3d( 1.0, 1.0, 2.0 ), R1 );
	sphere1.setPoleAxis( DSphere::POLE_OZ );

	double s_min = 0.0;
	double s_max = 0.5;
	double t_min =  0.25;
	double t_max =  0.75;
	int s_n = 6;
	int t_n = 2;

	for(int ii = 0; ii < 2; ii++){
		MeshViewSet* set = new MeshViewSet;

		double ds = (s_max - s_min) / s_n;
		double dt = (t_max - t_min) / t_n;
		double dparam_max = 0.0;

		for(int it = 0; it < t_n; it++) {
			double t0 = t_min + it * dt;
			double t1 = t0 + dt;

			for(int is = 0; is < s_n; is++) {
				double s0 = s_min + is * ds;
				double s1 = s0 + ds;

				DPoint2d params[3] = {	DPoint2d( s0, t0 ),	DPoint2d( s1, t0),
					DPoint2d( 0.5*(s0 + s1), t1) };

				DPoint3d pts[3], pt_middle = DPoint3d::zero;
				DPoint2d check_params[3], param_middle = DPoint2d::zero;
				double dparams[3];

				for(int j = 0; j < 3; j++ ){
					pts[j] = sphere1.getPoint( params[j] );
					check_params[j] = sphere1.getParam( pts[j] );
					dparams[j] = params[j].distance( check_params[j] );
					dparam_max = std::max( dparam_max, dparams[j] );
					param_middle.add( params[j], 1.0/3.0 );
					pt_middle.add( pts[j], 1.0/3.0 );
				}

				DVector3d fnormal = DVector3d::crossProduct( pts[0], pts[1], pts[2] );
				DVector3d snormal = sphere1.getNormalVector( param_middle );
				DPoint3d spt_middle = sphere1.getPoint( param_middle );

				set->addFace( pts[0], pts[1], pts[2] );
				set->addEdge( pt_middle, pt_middle + fnormal.normalized( R1/2 ), 0 );
				set->addEdge( spt_middle, spt_middle + snormal.normalized( R1/2 ), 1 );
			}
		}

		set->addInfo("dparam_max", dparam_max);
		SHOW_MESH( "DSphere test 1", set );

		sphere1.invertOrientation();
	}
}

DCylinder::DCylinder( const DPoint3d& center, const DVector3d& axis_vt, const double& radius, bool outside )
	: m_center(center), m_axis_vt(axis_vt), m_radius(radius), m_outside(outside)
{	
	assert( m_axis_vt.isNormal() );
	m_axis_vt.orthonormalVectors(m_e0, m_e1); 
	if( !m_outside ) m_e0 = -m_e0;
}

const DPoint3d DCylinder::getPoint(const DPoint2d& param) const 
{
	DPoint2d cpt = DCircle::getCirclePoint( DPoint2d::zero, m_radius, param.x);
	return m_center + (m_axis_vt * param.y) + m_e0 * cpt.x + m_e1 * cpt.y;
}

const DPoint2d DCylinder::getParam(const DPoint3d& pt) const 
{
	DMatrix3d A;
	A.setColumn(0, m_e0);
	A.setColumn(1, m_e1);
	A.setColumn(2, m_axis_vt);
	DVector3d res;
	bool result = A.solve(pt - m_center, res);
	assert( result );
	// res.x -> param.x
	// res.y and res.z -> change to alpha (param.y),
	return DPoint2d( DCircle::getCircleParam( DPoint2d::zero, m_radius, DPoint2d(res.x, res.y), 0.5), res.z );
}

bool DCylinder::getParamAndSignedDistance(const DPoint3d& point, DPoint2d& param, double & z) const
{
	DMatrix3d A;
	A.setColumn(0, m_e0);
	A.setColumn(1, m_e1);
	A.setColumn(2, m_axis_vt);
	DVector3d res;
	bool result = A.solve(point - m_center, res);
	assert( result );
	param.x = DCircle::getCircleParam( DPoint2d::zero, m_radius, DPoint2d(res.x, res.y), 0.5);
	param.y = res.z;
	z = sqrt(res.x*res.x + res.y*res.y) - m_radius;
	if(!m_outside) z = -z;

	return true;
}

const DPoint3d DCylinder::projectToSurface(const DPoint3d& pt) const
{
	assert( m_axis_vt.isNormal() );
	//return getPoint( getParam(pt) );
	double tq = m_axis_vt.scalarProduct(pt - m_center);
	DPoint3d axis_pt = m_center + m_axis_vt * tq;
	return axis_pt + (pt - axis_pt).normalized( m_radius );
}

double DCylinder::distanceFromSurface(const DPoint3d& pt) const
{
	assert( m_axis_vt.isNormal() );
	//return pt.distance( projectToSurface(pt) );
	double tq = m_axis_vt.scalarProduct(pt - m_center);
	DPoint3d axis_pt = m_center + m_axis_vt * tq;
	return abs( (pt - axis_pt).length() - m_radius );
}

/// Returns the principle curvature for this sphere and parameters [s,t]
const SurfaceCurvature DCylinder::getCurvature(const DPoint2d& pt , double & g_ratio) const
{ 
	g_ratio = 1.0; 
	double c = m_outside ? (1.0/m_radius) : (-1.0/m_radius);
	return SurfaceCurvature(0.0, c, 0.0);
}

const DVector3d DCylinder::getNormalVector(const DPoint2d& param) const
{
	DVector2d dt = DCircle::getCircleDerivative( m_radius, param.x );
	DVector3d dv = m_e0 * dt.x + m_e1 * dt.y;
	return dv.crossProduct( m_axis_vt ).normalized();
}

const DVector3d DCylinder::getNormalVectorDerivative(int deriv, const DPoint2d& param) const
{
	switch (deriv) {
	case DEquation::deriv_dt:
		return DVector3d::zero;
	case DEquation::deriv_ds:
	{
		DVector2d dt = DCircle::getCircleDerivative(m_radius, param.x);
		DVector3d dv = m_e0 * dt.x + m_e1 * dt.y;
		double vn_len = dv.crossProduct(m_axis_vt).length();

		DVector2d dt2 = DCircle::getCircleDerivative2(m_radius, param.x);
		DVector3d dv2 = (m_e0 * dt2.x + m_e1 * dt2.y);
		DVector3d result = dv2.crossProduct(m_axis_vt) / vn_len;

		//double h = 0.01;
		//DVector3d res[3];
		//for (int i = 0; i < 3; i++) {
		//	DVector3d n1 = getNormalVector(DPoint2d(param.x-h, param.y));
		//	DVector3d n2 = getNormalVector(DPoint2d(param.x+h, param.y));
		//	res[i] = (n2 - n1) / (2 * h);
		//	LOG4CPLUS_INFO(MeshLog::logger_mesh, "res[" << i << "]\t" << res[i]);
		//	h /= 2;
		//}
		//LOG4CPLUS_INFO(MeshLog::logger_mesh, "result-n.s\t" << result);

		return result;
	}
	case DEquation::deriv_dtt:
	case DEquation::deriv_dst:
		return DVector3d::zero;
	case DEquation::deriv_dss:
	{
		DVector2d dt = DCircle::getCircleDerivative(m_radius, param.x);
		DVector3d dv = m_e0 * dt.x + m_e1 * dt.y;
		double vn_len = dv.crossProduct(m_axis_vt).length();

		DVector2d dt3 = DCircle::getCircleDerivative3(m_radius, param.x);
		DVector3d dv3 = (m_e0 * dt3.x + m_e1 * dt3.y);
		DVector3d result = dv3.crossProduct(m_axis_vt) / vn_len;

		double h = 0.01;
		DVector3d res[3];
		for (int i = 0; i < 3; i++) {
			DVector3d n1 = getNormalVectorDerivative(DEquation::deriv_ds, DPoint2d(param.x-h, param.y));
			DVector3d n2 = getNormalVectorDerivative(DEquation::deriv_ds, DPoint2d(param.x+h, param.y));
			res[i] = (n2 - n1) / (2 * h);
			LOG4CPLUS_INFO(MeshLog::logger_mesh, "res[" << i << "]\t" << res[i]);
			h /= 2;
		}
		LOG4CPLUS_INFO(MeshLog::logger_mesh, "result-n.ss\t" << result);

		return result;
	}
	case DEquation::deriv_dsst:
	case DEquation::deriv_dstt:
	case DEquation::deriv_dttt:
		return DVector3d::zero;
	case DEquation::deriv_dsss:
	{
		assert(false);
		//DVector2d dt = DCircle::getCircleDerivative4(m_radius, param.x);
		//return m_e0 * dt.x + m_e1 * dt.y;
		return DVector3d::zero;
	}
	default:
		assert(false);
	}
	return DVector3d::zero;
}

const DVector3d DCylinder::getNormalVectorForPoint3d(const DPoint3d& pt) const
{
	assert( m_axis_vt.isNormal() );
	double tq = m_axis_vt.scalarProduct(pt - m_center);
	DPoint3d axis_pt = m_center + m_axis_vt * tq;
	DVector3d vn = (pt - axis_pt).normalized();
	return m_outside ? vn : -vn;
}

const DVector3d DCylinder::getDerivative(int deriv, const DPoint2d& param) const
{
	switch(deriv){
	case DEquation::deriv_dt:
	{
		//double h = 0.01;
		//DVector3d res[3];
		//for (int i = 0; i < 3; i++) {
		//	DPoint3d pt1 = getPoint(DPoint2d(param.x, param.y-h));
		//	DPoint3d pt2 = getPoint(DPoint2d(param.x, param.y+h));
		//	res[i] = (pt2 - pt1) / (2 * h);
		//	LOG4CPLUS_INFO(MeshLog::logger_mesh, "res[" << i << "]\t" << res[i]);
		//	h /= 2;
		//}
		//LOG4CPLUS_INFO(MeshLog::logger_mesh, "result-dt\t" << m_axis_vt);

		return m_axis_vt;
	}
	case DEquation::deriv_ds:
		{
			DVector2d dt = DCircle::getCircleDerivative( m_radius, param.x );
			DVector3d result = m_e0 * dt.x + m_e1 * dt.y;

			//double h = 0.01;
			//DVector3d res[3];
			//for (int i = 0; i < 3; i++) {
			//	DPoint3d pt1 = getPoint(DPoint2d(param.x - h, param.y));
			//	DPoint3d pt2 = getPoint(DPoint2d(param.x + h, param.y));
			//	res[i] = (pt2 - pt1) / (2 * h);
			//	LOG4CPLUS_INFO(MeshLog::logger_mesh, "res["<<i<<"]\t" << res[i]);
			//	h /= 2;
			//}
			//LOG4CPLUS_INFO(MeshLog::logger_mesh, "result-ds\t" << result);

			return result;
		}
	case DEquation::deriv_dtt:
	case DEquation::deriv_dst:
		return DVector3d::zero;
	case DEquation::deriv_dss:
		{
			DVector2d dt = DCircle::getCircleDerivative2( m_radius, param.x );
			DVector3d result = m_e0 * dt.x + m_e1 * dt.y;

			double h = 0.01;
			DVector3d res[3];
			for (int i = 0; i < 3; i++) {
				DVector3d v1 = getDerivative(DEquation::deriv_ds, DPoint2d(param.x-h, param.y));
				DVector3d v2 = getDerivative(DEquation::deriv_ds, DPoint2d(param.x+h, param.y));
				res[i] = (v2 - v1) / (2 * h);
				LOG4CPLUS_INFO(MeshLog::logger_mesh, "res[" << i << "]\t" << res[i]);
				h /= 2;
			}
			LOG4CPLUS_INFO(MeshLog::logger_mesh, "result-dss\t" << result);

			return result;
		}
	case DEquation::deriv_dttt:
	case DEquation::deriv_dsst:
	case DEquation::deriv_dstt:
		return DVector3d::zero;
	case DEquation::deriv_dsss:
		{
			DVector2d dt = DCircle::getCircleDerivative3( m_radius, param.x );
			return m_e0 * dt.x + m_e1 * dt.y;
		}
	default:
		assert(false);
	}
	return DVector3d::zero;
}

bool DCylinder::invertOrientation()
{
	m_outside = !m_outside;
	m_e0 = -m_e0;
	return true;
}

void DCylinder::test() 
{
	const double R1 = 1.0; // 5.5;

	DCylinder surf1( 
		DPoint3d( 1.0, 1.0, 2.0 ), 
		DVector3d::v_oz,
		//DVector3d( 1.0, 1.0, 1.0 ).normalized(), 
		R1 );

	double s_min = 0.0;
	double s_max = 0.5;
	double t_min =  0.25;
	double t_max =  0.75;
	int s_n = 6;
	int t_n = 2;

	for(int ii = 0; ii < 2; ii++){

		MeshViewSet* set = new MeshViewSet;

		double ds = (s_max - s_min) / s_n;
		double dt = (t_max - t_min) / t_n;
		double dparam_max = 0.0;

		for(int it = 0; it < t_n; it++) {
			double t0 = t_min + it * dt;
			double t1 = t0 + dt;

			for(int is = 0; is < s_n; is++) {
				double s0 = s_min + is * ds;
				double s1 = s0 + ds;

				DPoint2d params[3] = {	DPoint2d( s0, t0 ),	DPoint2d( s1, t0),
					DPoint2d( 0.5*(s0 + s1), t1) };

				DPoint3d pts[3], pt_middle = DPoint3d::zero;
				DPoint2d check_params[3], param_middle = DPoint2d::zero;
				double dparams[3];

				for(int j = 0; j < 3; j++ ){
					pts[j] = surf1.getPoint( params[j] );
					check_params[j] = surf1.getParam( pts[j] );
					dparams[j] = params[j].distance( check_params[j] );
					dparam_max = std::max( dparam_max, dparams[j] );
					param_middle.add( params[j], 1.0/3.0 );
					pt_middle.add( pts[j], 1.0/3.0 );
				}

				DVector3d fnormal = DVector3d::crossProduct( pts[0], pts[1], pts[2] );
				DVector3d snormal = surf1.getNormalVector( param_middle );
				DPoint3d spt_middle = surf1.getPoint( param_middle );

				set->addFace( pts[0], pts[1], pts[2] );
				set->addEdge( pt_middle, pt_middle + fnormal.normalized( R1/2 ), 0 );
				set->addEdge( spt_middle, spt_middle + snormal.normalized( R1/2 ), 1 );
			}
		}

		set->addInfo("dparam_max", dparam_max);
		SHOW_MESH( "DCylinder test 1", set );

		surf1.invertOrientation();
	}
}

DBox DEllipsoid::getBoundingBox() const {
	return DBox(m_center.x - m_rx, m_center.x + m_rx,
		m_center.y - m_ry, m_center.y + m_ry,
		m_center.z - m_rz, m_center.z + m_rz);
}
