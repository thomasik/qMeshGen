/////////////////////////////////////////////////////////////////////////////
// MeshEdge3d.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2002-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#if !defined(MESHEDGE3D_H__INCLUDED)
#define MESHEDGE3D_H__INCLUDED

#include <memory>

#include "DPoint.h"
#include "MeshData.h"
#include "TagBorder.h"
#include "TagExtended.h"
#include "common.h"

class MeshPoint3d;
class MeshFace;
class MeshBlock;
class MeshBoundaryCondition;
struct MeshViewEdgeData;
class Metric3dContext;
class Curve3dParametric;
class SurfaceParametric;

/**
 * This class implements a mesh edge in 3D and several operations.
 */
class MeshEdge3d : public TagBorder, public TagExtended
{
public:
	struct ActiveEdge{
		ActiveEdge(MeshEdge3d* _edge = 0, double _len = 0.0) : edge(_edge), len(_len) {}
		MeshEdge3d* edge;
		double len;
		//
		short compareTo(ActiveEdge* ae) const {
			if(len == ae->len) return 0;
			else if(len > ae->len) return 1;
			else return -1;
		}
		void preDeleteAll() {}
		int getIndex() const { return index; }
		void setIndex(int i) { index = i; }
		int index;
	};
	struct ActiveEdgeInfo {
		ActiveEdge* act;
		MeshPoint3d *p0, *p1;
		ActiveEdgeInfo(ActiveEdge* _act = nullptr, MeshPoint3d *_p0 = nullptr, MeshPoint3d *_p1 = nullptr)
			: act(_act), p0(_p0), p1(_p1) {}
	};
public:
	/// Standard constructor
	MeshEdge3d(MeshPoint3d* point1, MeshPoint3d* point2, MeshFace* face = nullptr);
	/// Standard destructor
	virtual ~MeshEdge3d();
public:
	/// Enumerates all adjacent blocks
	bool adjacentBlocks(DataVector<MeshBlock*> & blocks, bool boundary_check = true) const;
	/// Enumerates all adjacent blocks
	DataVector<const MeshBlock*> adjacentBlocks( bool boundary_check = true) const;
	DataVector<MeshBlock*> adjacentBlocks(bool boundary_check = true);
	/// Switches the given vertex in the edge with the other one
	void switchPoints(const MeshPoint3d* point1, MeshPoint3d* point2);
	/// fills the provided data vector with polyline approximation
	virtual int approximatePolyline(DataVector<DPoint3d> & polyline) const; 
	/// Removes link to this point (fo delete-all phase) - returns true if last one
	bool removePointLink(MeshPoint3d* point);
	/// Returns the type of this object
	virtual ElementType getType() const { return EDGE_SIMPLE_3D; }
	/// Returns the length of this edge
	double getLength() const;
	/// Returns the length (squared) of this edge
	double getLength2() const;
	/// Returns the length (squared) of this edge
	double getLength2(Metric3dContext& mc) const;
	/// Returns the length of this edge (for variable metric)
	double getLengthMetricAdapted(Metric3dContext& mc, double t0, double t1, int level = 0) const;
	/// Returns the length of this edge (for variable metric)
	double getLengthMetricAdapted(Metric3dContext& mc) const;
	/// Returns the metric length quality of this edge
	double getLengthQuality(Metric3dContext& mc, bool ext_metric) const;
	/// Sets the given array of points as "inner points" for this edge (previous array is additionally returned and should be deleted)
	int addInnerPoints(int ct, MeshPoint3d **points, MeshPoint3d*** old_points);
	/// Returns the coordinates of point, which divides the edge into segements with lengths in the specified ratio (
	DPoint3d getAspectPoint(double ksi) const;
	/// Returns true, if both edges have a common vertex 
	bool incidentTo(const MeshPoint3d* point) const {
		return points[0] == point || points[1] == point;
	}
	/// Returns true, if both edges have a common vertex 
	bool incidentTo(const MeshEdge3d* edge) const {
		return (points[0] != edge->points[0] && points[0] != edge->points[1] &&
			points[1] != edge->points[0] && points[1] != edge->points[1]);
	}
	/// Returns common vertex for two edges, or nullptr
	MeshPoint3d* commonVertex(const MeshEdge3d* edge) const;
	/// Returns the description of this point for visualization
	std::shared_ptr<MeshViewEdgeData> getViewData() const;
	/// Returns one of the vertices of this edge
	MeshPoint3d* getMeshPoint(int i) const { return points[i]; }
	/// Returns the index of the given face in the array of adjacent faces
	int getFaceIndex(MeshFace* const face) const;
	/// Returns the coordinates of the point at the given ratio [0,1] between vertices
	DPoint3d getPoint(double t) const;
	/// Returns the coordinates of the point at the given ratio [0,1] between vertices (for edge lying on surface)
	//DPoint3d getPointOnSurface( Metric3dContext& mc, double t ) const;
	/// Returns the length of the edge
	double getLength(Metric3dContext& mc, bool local_metric = false) const;
	/// Returns the length of the edge
	double getLengthNoMetric() const;
	/// Returns the index (0 or 1) of the point as the vertice of this edge
	int getPointIndex(const MeshPoint3d* point) const;
	/// Removes the reference to the given face from the array of adjacent faces
	bool removeFaceLink(MeshFace* const face);
	/// Removes the reference to the given face from the array of adjacent faces
	bool removeFaceLinkDisregardBoundary(MeshFace* const face);
	/// Adds the reference to the given face to the array of adjacent faces
	void addFaceLink(MeshFace* const face);
	/// Returns the number of faces adjacent to this edge
	int getFaceCount() const { return (int)faces.countInt(); }
	/// Returns the face at the given index in the adjacency array
	MeshFace* getFaceAt(int i) const { return faces[i]; }
	/// Return face incident to this edge and the given point 
	MeshFace* getFaceToPoint(const MeshPoint3d* mpt) const;
	/// Checks whether there is any face adjacent to this edge
	bool isBounded() const { return faces.countInt() > 0; }
	/// Returns the other vertex of the edge
	MeshPoint3d* getOtherPoint(const MeshPoint3d* pt) const { return (points[0]==pt) ? points[1] : points[0]; }
	/// Returns the other face of the edge (if there are only two!)
	MeshFace* getOtherFace(const MeshFace* face) const { 
		assert(faces.countInt() == 2); 
		return (faces[0] == face) ? faces[1] : faces[0]; 
	}
	/// Returns the other boundary face of the edge (supposing there are only two!)
	MeshFace* getOtherBorderFace(const MeshFace* face) const;
	/// Returns the local curve
	std::shared_ptr<const Curve3dParametric>  getLocalCurve() const { return local_curve; }
	/// Changes local curve
	void setLocalCurve(std::shared_ptr<const Curve3dParametric> curve) { local_curve = curve; }
	/// Clears local curve
	void clearLocalCurve() { local_curve = nullptr; }
	/// Checks local curve
	bool hasLocalCurve() const { return (local_curve != nullptr); }
protected:
	/// Vertices of this edge (array of references to mesh points)
	MeshPoint3d*	points[2];
	/// Adjacent faces 
	DataVector<MeshFace*>	faces;
	/// Local curve
	std::shared_ptr<const Curve3dParametric> local_curve;
};

#endif // !defined(MESHEDGE3D_H__INCLUDED)
