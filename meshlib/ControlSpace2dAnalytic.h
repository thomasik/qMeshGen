/////////////////////////////////////////////////////////////////////////////
// ControlSpace2dAnalytic.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2003-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#if !defined(CONTROLSPACEANALYTIC_H__INCLUDED)
#define CONTROLSPACEANALYTIC_H__INCLUDED

#include "ControlSpace2d.h"
class DEquation;

/**
 * This class implements a control space given by a set of equations.
 */
class ControlSpace2dAnalytic : public ControlSpace2d  
{
public:
	/// Standard contructor
	ControlSpace2dAnalytic(SurfaceConstPtr surface, 
		const char * default_lx_str, const char * default_ly_str, const char * default_angle_str);
	/// Standard contructor
	ControlSpace2dAnalytic(SurfaceConstPtr surface, 
		DEquation* lx_eq, DEquation* ly_eq, DEquation* angle_eq);
	/// Destructor
	virtual ~ControlSpace2dAnalytic();
public:
	/// Returns the type of element (should be reimplemented in derived classes)
	virtual int getType() const { return MeshData::CONTROL_ANALYTICAL; }
	/// Get sizing info (matrix mode) at the given point
	virtual ControlDataMatrix2d getMetricAtPoint(const DPoint2d& pt) const;
public:
	/// Insert additional function-case
	void addCaseAsFirst(const char * test_str, 
		const char * lx_str, const char * ly_str, const char * angle_str);
	/// Insert additional function-case (two checks)
	void addCaseAsFirst(const char * test1_str, const char * test2_str, 
		const char * lx_str, const char * ly_str, const char * angle_str);
protected:
	/// struct for additionale function-cases
	struct FunctionNode{
		DEquation* condition;
		DEquation* condition2;
		DEquation* lx;
		DEquation* ly;
		DEquation* angle;
		FunctionNode* next;
	} *m_list;
	/// Default x-length of elements
	DEquation *m_default_lx;
	/// Default y-length of elements
	DEquation *m_default_ly;
	/// Default angle of stretching of elements
	DEquation *m_default_angle;
};

#endif // !defined(CONTROLSPACEANALYTIC_H__INCLUDED)
