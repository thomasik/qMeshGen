/////////////////////////////////////////////////////////////////////////////
// FrontFace.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2008-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#if !defined(FRONTFACE_H__INCLUDED)
#define FRONTFACE_H__INCLUDED

#include "MeshData.h"
#include "DataVector.h"

class MeshFace;
class MeshPoint3d;
class Metric3dContext;

/**
 * This class represents a face of a "running front" used during frontal boundary-constrained meshing.
 * Each face has pointers to the neighbouring front-faces (thus forming a closed cycle) and is
 * clasified according to face-to-face angles.
 */
class FrontFace : public IndexedObject
{
public:
	/// Standard constructor
	FrontFace(MeshFace* face = nullptr, bool oriented = true, int lev = 0);
	/// Standard destructor
	~FrontFace();
public:
	/// clear some data for faster all-delete process
	void preDeleteAll() {} // nothing here acutally
	/// Gathers incidency data for this face (neighbours and angles)
	FrontFace* init(Metric3dContext& mc);
	/// Update connectivity for neighbors
	void gatherNeighbors(DataVector<FrontFace*> &ffaces);
	/// Check if the given front face is among the neighbors
	bool adjacentTo(FrontFace* fface) const;
	/// Returns the pointer to the mesh-face which forms this front-face
	MeshFace* getFace() const { return m_face; }
	/// Returns the points to the i-th neighbor front face
	FrontFace* getNeighbor(int i) const { return m_neighbors[i]; }
	/// Returns the angle to the i-th neighbor 
	double getAngle(int i) const { return m_angles[i]; }
	/// Returns the angle to the given neighbor 
	double getAngle(const FrontFace* ff) const;
	/// Whether oriented correctly
	bool isOriented() const { return m_oriented; }
	/// Clear front facee (temporarily)
	void clear();
	/// Initialize the front face by setting the assigned mesh-face
	void setFace(MeshFace* face, bool oriented);
	/// Compares the rank of this front-face with the given one (required by the class DataContainer)
	short compareTo(const FrontFace* item) const;
	/// Calculates dihedral angle between front faces
	double calculateDihedralAngle(FrontFace* ff, Metric3dContext& mc) const;
public:
	/// Returns i-th point of the face (with respect to orientation)
	MeshPoint3d* getPoint(int i) const;
	/// Returns other point from the (min-angle) incident front face
	MeshPoint3d* getIncidentPoint(int i) const;
	/// Returns face of the i-th incident front-face
	MeshFace* getIncidentFace(int i) const;
	/// Returns index of the edge with smallest face-face angle
	int getMinAngleIndex() const { return m_min_angle_index; }
protected:
	/// Pointer to the assigned mesh-face
	MeshFace *m_face;
	/// Valid orientation
	bool m_oriented;
	/// Pointers to the adjacent front-faces
	DataVector<FrontFace*> m_neighbors;
	/// Angle between this face and the adjacent front-faces
	DataVector<double> m_angles;
	/// The number of the "front layer" (starting from 0 and increasing)
	int		m_level;
	/// Minimum face-to-face angle for neighbors (index)
	int m_min_angle_index;
};

#endif // !defined(FRONTFACE_H__INCLUDED)
