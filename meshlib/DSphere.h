#pragma once

#ifndef DSPHERE_H_INCLUDED
#define DSPHERE_H_INCLUDED

#include "DPoint.h"
#include "DRect.h"
#include "SurfaceCurvature.h"

class DCircle
{
public:
	DCircle(const DPoint2d& center = DPoint2d::zero, const double& radius = 0.0)
		: m_center(center), m_radius(radius) {}
	DRect getBoundingBox() const { return DRect( m_center.x - m_radius, m_center.x + m_radius, 
		m_center.y - m_radius, m_center.y + m_radius); }
	double implicitValue( const DPoint2d& pt ) const ;
	DPoint2d projectToCurve(const DPoint2d& pt) const;
	double getParam( const DPoint2d& pt, double ts ) const;
	static double getCircleParam( const DPoint2d& mid, double r, const DPoint2d& pt, double ts );
	static DPoint2d getCirclePoint( const DPoint2d& mid, double r, double t );
	static DVector2d getCircleDerivative( double r, double t );
	static DVector2d getCircleDerivative2( double r, double t );
	static DVector2d getCircleDerivative3( double r, double t );
public:
	DPoint2d m_center;
	double m_radius;
};

class DSphere
{
public:
	DSphere(const DPoint3d& center = DPoint3d::zero, const double& radius = 0.0)
		: m_center(center), m_radius(radius), m_pole(POLE_OZ) {}
	DBox getBoundingBox() const { 
		double r = getRadius();
		return DBox( m_center.x - r, m_center.x + r, 
					 m_center.y - r, m_center.y + r, 
					 m_center.z - r, m_center.z + r); 
	}
	double getRadius() const { return (m_radius >= 0.0) ? m_radius : -m_radius; }
	const DPoint3d& getCenter() const { return m_center; }
	double implicitValue( const DPoint3d& pt ) const ;
	DPoint3d projectToSurface(const DPoint3d& pt) const;

	const SurfaceCurvature getCurvature(const DPoint2d& pt, double & g_ratio) const;
	const DPoint3d getPoint(const DPoint2d& param) const;
	const DPoint2d getParam(const DPoint3d& point) const;
	bool getParamAndSignedDistance(const DPoint3d& point, DPoint2d& param, double & z) const;
	const DPoint2d getParamNear(const DPoint3d& point, const DPoint2d& /* near_point */) const
		{ return getParam(point); }
	const DVector3d getNormalVector(const DPoint2d& param) const;
	const DVector3d getNormalVectorDerivative(int deriv, const DPoint2d& param) const;
	const DVector3d getNormalVectorForPoint3d(const DPoint3d& pt) const;
	const DVector3d getDerivative(int deriv, const DPoint2d& param) const;
	bool invertOrientation();

	static void test();
public:
	enum PoleAxis { POLE_OX = 0, POLE_OY = 1, POLE_OZ = 2 };
	void setPoleAxis( PoleAxis pa) { m_pole = pa; }
public:
	DPoint3d m_center;
	double m_radius;
	PoleAxis m_pole;
};

class DEllipsoid
{
public:
	DEllipsoid(const DPoint3d& center = DPoint3d::zero, const double& rx = 0.0,
		const double& ry = 0.0, const double& rz = 0.0)
		: m_center(center), m_rx(rx), m_ry(ry), m_rz(rz), m_pole(POLE_OZ) {}
	DBox getBoundingBox() const;
	const DPoint3d& getCenter() const { return m_center; }
	double implicitValue(const DPoint3d& pt) const;
	//DPoint3d projectToSurface(const DPoint3d& pt) const;

	double getRadius(int i) const;
	double& getRadius(int i);

	//const SurfaceCurvature getCurvature(const DPoint2d& pt, double & g_ratio) const;
	const DPoint3d getPoint(const DPoint2d& param) const;
	//const DPoint2d getParam(const DPoint3d& point) const;
	//bool getParamAndSignedDistance(const DPoint3d& point, DPoint2d& param, double & z) const;
	//const DPoint2d getParamNear(const DPoint3d& point, const DPoint2d& /* near_point */) const
	//	{	return getParam(point);	}
	const DVector3d getNormalVector(const DPoint2d& param) const;
	//const DVector3d getNormalVectorDerivative(int deriv, const DPoint2d& param) const;
	//const DVector3d getNormalVectorForPoint3d(const DPoint3d& pt) const;
	//const DVector3d getDerivative(int deriv, const DPoint2d& param) const;
	//bool invertOrientation();

	//static void test();
public:
	enum PoleAxis { POLE_OX = 0, POLE_OY = 1, POLE_OZ = 2 };
	void setPoleAxis(PoleAxis pa) { m_pole = pa; }
public:
	DPoint3d m_center;
	double m_rx, m_ry, m_rz;
	PoleAxis m_pole;
};

class DCylinder
{
public:
	DCylinder(const DPoint3d& center = DPoint3d::zero, const DVector3d& axis_vt = DVector3d::v_ox, 
		const double& radius = 0.0, bool right_oriented = true);
public:
	void set(const DPoint3d& center, const DVector3d& axis_vt, const double& radius, bool outside = true) {
		m_center = center;
		m_axis_vt = axis_vt;
		m_radius = radius;
		m_axis_vt.orthonormalVectors(m_e0, m_e1);
		m_outside = outside;
		if(!m_outside) m_e0 = -m_e0;
	}
	const DPoint3d getPoint(const DPoint2d& param) const;
	const DPoint2d getParam(const DPoint3d& pt) const;
	bool getParamAndSignedDistance(const DPoint3d& point, DPoint2d& param, double & z) const;
	const SurfaceCurvature getCurvature(const DPoint2d& pt, double & g_ratio) const;
	const DPoint2d getParamNear(const DPoint3d& point, const DPoint2d& /* near_point */) const
		{ return getParam(point); }
	const DVector3d getNormalVector(const DPoint2d& param) const;
	const DVector3d getNormalVectorDerivative(int deriv, const DPoint2d& param) const;
	const DVector3d getNormalVectorForPoint3d(const DPoint3d& pt) const;
	const DVector3d getDerivative(int deriv, const DPoint2d& param) const;
	bool invertOrientation();

	const DPoint3d projectToSurface(const DPoint3d& pt) const;
	double distanceFromSurface(const DPoint3d& pt) const;

	static void test();
public:
	DPoint3d m_center;
	DVector3d m_axis_vt;
	double m_radius;
	bool m_outside;
public:
	DVector3d m_e0, m_e1;
};

#endif

