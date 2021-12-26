/////////////////////////////////////////////////////////////////////////////
// MeshDomainEdge3d.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2008-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#if !defined(MESHDOMAINEDGE3D_H__INCLUDED)
#define MESHDOMAINEDGE3D_H__INCLUDED

#include <memory>

#include "MeshPoint3d.h"
#include "MeshEdge3d.h"
#include "DataVector.h"
#include "common.h"

class ControlSpace3dAdaptive;

/**
 * This class implements a mesh edge in 3D and several operations.
 */
class MeshDomainEdge3d : public MeshEdge3d
{
public:
	/// Standard constructor
	MeshDomainEdge3d(MeshPoint3d* point1, MeshPoint3d* point2, MeshFace* face = nullptr);
	/// Standard destructor
	virtual ~MeshDomainEdge3d();
public:
	bool updateACS(std::shared_ptr<ControlSpace3d> acs) const;
	/// Add free-point for this domain-edge
	void addFreePoint(const std::shared_ptr<MeshPoint3d> & fpoint);
	/// Returns the set of free-points
	DataVector<std::shared_ptr<MeshPoint3d>>* getFreePoints() { return &m_freepoints; }
	/// Returns the type of this object
	virtual ElementType getType() const { return EDGE_DOMAIN_3D; }
	/// Returns reference to discretization
	DataVector<std::shared_ptr<MeshPoint3d>> & getDiscretization() { return m_discretization_points; }
	/// Clear discretization
	void clearDiscretization() { m_discretization_points.clear(); m_discretization_valid = false; }
	/// Set valid flag
	void setValidDiscretization(bool valid = true) { m_discretization_valid = valid; }
	/// Check valid flag
	bool isValidDiscretization() const { return m_discretization_valid; }
protected:
	/// Vertices of this edge (array of references to mesh points)
	DataVector<std::shared_ptr<MeshPoint3d>>	m_discretization_points;
	/// Valid flag
	bool m_discretization_valid;
	/// FreePoints
	DataVector<std::shared_ptr<MeshPoint3d>> m_freepoints;
};

#endif // !defined(MESHDOMAINEDGE3D_H__INCLUDED)
