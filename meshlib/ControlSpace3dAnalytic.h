/////////////////////////////////////////////////////////////////////////////
// ControlSpace3dAnalytic.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2006-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#if !defined(CONTROLSPACE3DANALYTIC_H__INCLUDED)
#define CONTROLSPACE3DANALYTIC_H__INCLUDED

#include "MeshData.h"
#include "ControlSpace3d.h"
#include "DEquation.h"
#include <functional>

/**
 * This class implements a control space given by a set of equations.
 */
class ControlSpace3dAnalytic : public ControlSpace3d
{
public:
	/// Standard contructor
	ControlSpace3dAnalytic(
		const char * default_lx_str, const char * default_ly_str, 
		const char * default_lz_str, const char * default_angle1_str, 
		const char * default_angle2_str, const char * default_angle3_str,
		DEquationConstTable* ctable);
	/// Standard contructor
	ControlSpace3dAnalytic(DEquation* lx_eq, DEquation* ly_eq, DEquation* lz_eq, 
		DEquation* angle1_eq, DEquation* angle2_eq, DEquation* angle3_eq);
	/// Destructor
	virtual ~ControlSpace3dAnalytic();
public:
	/// Returns the type of element (should be reimplemented in derived classes)
	virtual int getType() const { return MeshData::CONTROL_ANALYTICAL_3D; }
	/// Get sizing info (matrix mode) at the given point
	virtual ControlDataMatrix3d getMetricAtPoint(const DPoint3d& pt) const;
public:
	/// Insert additional function-case
	void addCaseAsFirst(
		const char * test_str1, const char * test_str2, 
		const char * lx_str, const char * ly_str, const char * lz_str, 
		const char * angle1_str, const char * angle2_str, const char * angle3_str,
		DEquationConstTable* ctable);
protected:
	/// struct for additionale function-cases
	struct FunctionNode{
		DEquation* condition[2];
		DEquation* len[3];
		DEquation* angle[3];
		FunctionNode* next;
	} *m_list;
	/// Default length of elements
	DEquation* m_default_len[3];
	/// Default angles of stretching of elements
	DEquation* m_default_angle[3];
};

class ControlSpace3dFunctional: public ControlSpace3d
{
public:
	/// Standard contructor
	ControlSpace3dFunctional(const std::function<ControlDataMatrix3d(const DPoint3d& pt)> &_fm)
		: fm(_fm) {}
public:
	/// Returns the type of element (should be reimplemented in derived classes)
	virtual int getType() const { return MeshData::CONTROL_FUNCTIONAL_3D; }
	/// Get sizing info (matrix mode) at the given point
	virtual ControlDataMatrix3d getMetricAtPoint(const DPoint3d& pt) const { return fm(pt); }
protected:
	std::function<ControlDataMatrix3d(const DPoint3d& pt)> fm;
};

#endif // !defined(CONTROLSPACE3DANALYTIC_H__INCLUDED)
