/////////////////////////////////////////////////////////////////////////////
// MeshBoundaryCondition.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2004-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#if !defined(MESHBOUNDARYCONDITION_H__INCLUDED)
#define MESHBOUNDARYCONDITION_H__INCLUDED

/**
 * This class contains boundary-condition information.
 */
class MeshBoundaryCondition
{
public:
	/// Standard constructor
	MeshBoundaryCondition(const string & cond = "") : m_cond(cond) {}
public:
	/// Returns condition text
	string getCondition() const { return m_cond; }
	/// Returns text description of this curve
	void setCondition(const string& cond) { m_cond = cond; }
protected:
	/// Condition text
	string m_cond;	
};

#endif // !defined(MESHBOUNDARYCONDITION_H__INCLUDED)
