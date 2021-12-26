// ControlSpace2dIdentity.h: interface for the ControlSpace2dIdentity class.
//
//////////////////////////////////////////////////////////////////////

#pragma once

#if !defined(CONTROLSPACEIDENTITY_H__INCLUDED)
#define CONTROLSPACEIDENTITY_H__INCLUDED

#include "ControlSpace2d.h"
#include "ControlSpace3d.h"
#include "MeshData.h"

class ControlSpace2dIdentity : public ControlSpace2d  
{
public:
	/// Standard constructor
	ControlSpace2dIdentity(double length = 1.0, SurfaceConstPtr surface = nullptr) 
		: ControlSpace2d(surface), m_cdm(length, 0, length) {}
	/// Standard constructor
	ControlSpace2dIdentity(const ControlDataMatrix2d& cdm, SurfaceConstPtr surface = nullptr) 
		: ControlSpace2d(surface), m_cdm(cdm) {}
public:
	/// Returns the type of element (should be reimplemented in derived classes)
	virtual int getType() const { return MeshData::CONTROL_IDENTITY; }
	/// Get sizing info (matrix mode) at the given point
	virtual ControlDataMatrix2d getMetricAtPoint(const DPoint2d& /* pt */) const { 
		return m_cdm; }
public:
	void setMinLength(double ml) { m_cdm = ControlDataMatrix2d(ml, 0.0, ml); }
	double getMinLength() const { return m_cdm.m11; }
private:
	ControlDataMatrix2d m_cdm;
};

class ControlSpace3dIdentity : public ControlSpace3d
{
public:
	/// Standard constructor
	ControlSpace3dIdentity(double length = 1.0) : m_length(length) {}
public:
	/// Returns the type of element (should be reimplemented in derived classes)
	virtual int getType() const { return MeshData::CONTROL_IDENTITY_3D; }
	/// Get sizing info (matrix mode) at the given point
	virtual ControlDataMatrix3d getMetricAtPoint(const DPoint3d& /* pt */) const { 
		return ControlDataMatrix3d(m_length, m_length, m_length, 0.0, 0.0, 0.0); }
private:
	double m_length;
};

#endif // !defined(CONTROLSPACEIDENTITY_H__INCLUDED)
