/////////////////////////////////////////////////////////////////////////////
// MeshFace.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2002-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#if !defined(MESHFACE_H__INCLUDED)
#define MESHFACE_H__INCLUDED

#include <memory>

#include "DRect.h"
#include "DPoint.h"
#include "DVector.h"
#include "MeshData.h"
#include "TagBorder.h"
#include "TagExtended.h"
#include "DTriangle.h"
//#include "DPlanarQuadric.h"

class MeshPoint3d;
class MeshEdge3d;
class MeshBlock;
struct MeshViewFaceData;
class Metric3dContext;
class SurfaceParametric;
class DPlane;

/**
 * This class implements a mesh face in 3D.
 */
class MeshFace : public IndexedObject, public TagBorder, public TagExtended
{
public:
	/// Orientation of adjacent area (for directed edges)
	enum BlockOrientation {BLOCK_DOWN = 0, BLOCK_UP = 1};
	/// Class for storing 3D surface auxiliary data
	class SurfaceData {
	public:
		SurfaceData() : local_surface_orientation(0),
			cv0(-1.0), cv1(-1.0) {}
	public:
		DVector3d base_normal;
		std::shared_ptr<const SurfaceParametric> local_surface;
		int local_surface_orientation;
		DVector3d ce0, ce1; // curvature main directions 
		double cv0, cv1, cmax; // curvature values
		//DPlanarQuadric pq; // ???
	};
public:
	/// Standard constructor
	MeshFace(int _count, MeshEdge3d** _edges = nullptr, MeshPoint3d**	_points = nullptr);
	/// Standard destructor
	virtual ~MeshFace() {}
public:
	virtual DVector3d getSurfaceNormalForMidEdge(const MeshPoint3d* /* mp0 */, const MeshPoint3d* /* mp1 */) const {
		return getNormalVector();
	}
	/// Return normal for i-th edge, planar with face, outer
	DVector3d getOuterNormalForEdge(int i, double nvalue = 1.0);
	/// Swaps the i-th edge of this face (actually works only for triangles) with the adjacent one (can be conditionally only)
	virtual MeshEdge3d* swapWithNeighbour(Metric3dContext& /* mc */, 
		std::shared_ptr<const SurfaceParametric> /* surface */, 
		int /* i */, bool /* improve_only */, bool /* local_metric */ = false, 
		TagExtended::TagType /* side_tag */ = TagExtended::TAG_NONE) { assert(false); return nullptr; }
	/// count inverted-rank for this face
	int getInvertedRank(int * inv_index = nullptr) const;
	/// get all points for the given set of faces
	static void getMeshPointsForFaces( const DataVector<MeshFace*> & mfaces, 
		DataVector< MeshPoint3d* > & mpoints );
	/// count base plane and 2d bounding box
	bool countPlaneBBox( DPlane& plane, DRect& box2d ) const;
	/// index of neighbouring face
	int getNeighbourIndex(const MeshFace* other_face) const;
	/// switch orientation of the face (direction of vertices)
	void switchOrientation();
	/// copy additional data for cloning
	void copyAllExtraData(const MeshFace* source);
	/// create new edge adjacent to this face
	MeshEdge3d* createEdge(MeshPoint3d* p1, MeshPoint3d* p2);
	/// returns shape quality for a face
	virtual double getShapeQuality(Metric3dContext& /* mc */) const { return 1.0; }
	/// returns shape quality for a face
	virtual double getShapeQuality() const { return 1.0; }
	/// remove point, return resulting face (the same, or a new one)
	virtual MeshFace* removePoint(const MeshPoint3d* /* point */) { assert(false); return nullptr; }
	/// Split face into (mesh) triangles
	virtual bool splitToTriangles( DataVector< MeshFace* > & split_faces ) const;
	/// Split face into triangles
	virtual bool splitToTriangles( DataVector< DTriangle3d > & triangles ) const;
	/// Creates and returns a copy of this face (+ whole connectivity)
	virtual MeshFace* clone() const;
	/// Returns the edge between the two faces, or nullptr if not adjacent
	MeshEdge3d* getCommonEdge(MeshFace* face) const;
	/// Set tag for all adjacent edges
	void setIntTagForPoints(TagExtended::TagType type, int t = 1) const;
	/// Set tag for all adjacent edges
	void setTagForEdges(TagExtended::TagType type, int t = 1) const;
	// check, whether this face is valid (second version)
	bool validDirect(Metric3dContext& mc, double last_q = 0.01) const;
	// check, whether this face is valid with respect to the given surface
	virtual bool valid(std::shared_ptr<const SurfaceParametric> /* surface */) const { assert(false); return false; }
	/// Create face-edges connections
	void attachToEdges(); 
	/// Removes face-edges connections
	void detachFromEdges();
	/// Removes face-edges connections
	void detachFromEdgesDisregardBoundary();
	/// Switches one of vertices to other mesh point (with adequate update of edges)
	void switchPointsWithEdges(const MeshPoint3d* point1, MeshPoint3d* point2);
	/// Switches one of vertices to other mesh point (with adequate update of edges)
	void switchPointsWithEdgesBoundary(const MeshPoint3d* point1, MeshPoint3d* point2);
	/// Returns edge from this face adjacent to the given edge through point
	MeshEdge3d* getOtherEdge(const MeshEdge3d* edge, const MeshPoint3d* point);
	/// Removes link to this edge (for delete-all phase) - returns true if last one
	virtual bool removeEdgeLink(MeshEdge3d* edge);
	/// Removed link to (previously) adjacent block
	bool removeBlockLink(const MeshBlock *block);
	/// Removed link to (previously) adjacent block
	bool removeBlockLinkDisregardBoundary(const MeshBlock *block);
	/// Returns the first point different than the given two
	MeshPoint3d* getOtherPoint(const MeshPoint3d* pt0, const MeshPoint3d* pt1);
	/// Returns true if the given block is connected with this face
	bool incidentToBlock(const MeshBlock* block) const { return blocks[0] == block || blocks[1] == block; }
	/// Returns true if the given edge is one of the edges of this face
	bool incidentToEdge(const MeshEdge3d* edge) const;
	/// Returns the vertex before the given one (anticlockwise)
	MeshPoint3d* getPrevPoint(const MeshPoint3d* point) const;
	/// Returns the vertex after the given one (anticlockwise)
	MeshPoint3d* getNextPoint(const MeshPoint3d* point) const;
	/// Returns true if both faces are oriented the same way (checking vertices)
	bool sameOrientation(const MeshFace* face) const;
	/// Returns true if the sequence of the two given points complies to the orientation of face
	bool properOrientation(const MeshPoint3d* mp1, const MeshPoint3d* mp2) const;
	/// Returns the index of the given edge within this face
	int getEdgeIndex(const MeshEdge3d* edge) const;
	/// Returns the ith block adjacent to this face
	MeshBlock* getBlock(int i) const  { return blocks[i]; }
	/// Returns the other adjacent block of this face
	MeshBlock* getOtherBlock(const MeshBlock* block) const  { return (blocks[0]==block)?blocks[1]:blocks[0]; }
	/// Returns the vector normal to this face if possible
	virtual bool checkAndGetNormalVector(DVector3d& vn) const;
	/// Returns the vector normal to this face
	virtual DVector3d getNormalVector() const;
	/// Returns the orientation-index of this face for the given point
	int getOrientation(const DPoint3d& pt) const;
	/// Returns the index of block in block-array (i.e. its orientation)
	int getBlockIndex(const MeshBlock* block) const;
	/// Returns true if the given point is one of the vertices of this face
	bool incidentToPoint(const MeshPoint3d* point) const;
	/// Returns true if the given points are vertices of this face
	bool incidentToPoints(const MeshPoint3d* point1, const MeshPoint3d* point2) const;
	/// Returns the middle point for this element
	DPoint3d getMiddlePoint() const;
	/// Returns the bounding cubicoid for this face
	virtual DBox getBoundingBox() const;
	double getBoundingBoxDiameter() const { return getBoundingBox().getDiameter(); }
	/// Adds a reference to the given mesh block adjacent to this face (via given order of face_points)
	void setBlockLink(MeshBlock* block, const MeshPoint3d* mpt1, const MeshPoint3d* mpt2);
	/// Adds a reference to the given mesh block adjacent to this face (from the side of the given point)
	void setBlockLink(MeshBlock* block, const DPoint3d& top_point);
	/// Adds a reference to the given mesh block adjacent to this face
	void setBlockLink(MeshBlock* block, BlockOrientation side);
	/// Adds a reference to the given mesh block adjacent to this face
	void setBlockLink(MeshBlock* block, int i) { blocks[i] = block; }
	/// Removes a reference to the given mesh block adjacent to this face
	bool clearBlockLink(MeshBlock* block);
	/// Returns the number of points for this face
	virtual int getPointCount() const { return count; }
	/// Returns i-th vertex (reference to 3D mesh point) of this face
	MeshPoint3d* getPoint(int i) const { return points[i]; }
	/// Returns the number of edges for this face
	virtual int getEdgeCount() const { return count; }
	/// Returns i-th edge (reference to 3D mesh edge) of this face
	MeshEdge3d* getEdge(int i) const { return edges[i]; }
	/// Checks whether the face is adjacent to any mesh blocks
	bool isBounded() const { return (blocks[0] || blocks[1]); }
	/// Checks whether the face is adjacent to two mesh blocks
	bool isBoundedBothSides() const { return (blocks[0] && blocks[1]); }
	/// Sets the boundary type for this face (+ adjacent edges)
	void setBorderFlagsWithEdges(char border_flags = TagBorder::OUTER);
	/// Sets the boundary type for this face (+ adjacent edges)
	void setBorderWithEdges(char border_flags = TagBorder::OUTER);
	/// Sets the boundary type for this face (+ adjacent edges and pointss)
	void setBorderWithEdgesAndPoints(char border_flags = TagBorder::OUTER);
	/// Returns the type of this object
	virtual ElementType getType() const { return FACE_UNKNOWN; }
	/// Returns the description of this face for visualization
	virtual std::shared_ptr<MeshViewFaceData> getViewData(double shrink_ratio = 1.0, bool proper_orientation = true) const;
	/// clears surface data
	void clearSurfaceData();
	/// checks surface data
	bool hasLocalSurface() const {
		return (sdata != nullptr) && (sdata->local_surface != nullptr); 
	}
	/// Returns local surface
	std::shared_ptr<const SurfaceParametric> getCurrentLocalSurface() const {
		return ( sdata!= nullptr ) ? sdata->local_surface : nullptr;
	}
	/// stored max curvature
	double getMaxCurvature() const {
		return sdata ? sdata->cmax : 0.0;
	}
	/// Returns local surface
	std::shared_ptr<const SurfaceParametric> getOptLocalSurface();
	/// Sets local surface
	void setLocalSurface(std::shared_ptr<const SurfaceParametric> surface){
		if(!sdata) sdata = new SurfaceData();
		sdata->local_surface = surface;
	}
	/// checks surface orientation data
	bool hasLocalSurfaceOrientation() const {
		return (sdata != nullptr) && (sdata->local_surface_orientation != 0); 
	}
	/// Sets local surface orientation -1/1
	void setLocalSurfaceOrientation(int orient) {
		if(!sdata) sdata = new SurfaceData();
		sdata->local_surface_orientation = orient;
	}
	/// Returns local surface orientation
	int getLocalSurfaceOrientation() const {
		return (sdata == nullptr) ? 0 : sdata->local_surface_orientation;
	}
	/// Sets base normal vector
	void setBaseNormal(const DVector3d& vn) {
		if(!sdata) sdata = new SurfaceData();
		sdata->base_normal = vn;
	}
	/// Returns base normal
	const DVector3d& getBaseNormal() const {
		return (sdata == nullptr) ? DVector3d::zero : sdata->base_normal;
	}
	/// Returns normal calculated from surface for the given surface, at the middle of the face
	DVector3d getLocalSurfaceNormal(std::shared_ptr<const SurfaceParametric> surface = nullptr ) const;
	/// Checks for base normal
	bool hasBaseNormal() const {
		return (sdata != nullptr) && !sdata->base_normal.isZero();
	}
	/// Sets surface curvature data
	void setCurvatureData(const DVector3d& e0, const DVector3d& e1, double c0, double c1) {
		if(!sdata) sdata = new SurfaceData();
		sdata->ce0 = e0;
		sdata->ce1 = e1;
		sdata->cv0 = c0;
		sdata->cv1 = c1;
		sdata->cmax = std::max( c0, c1 );
	}
	/// Sets approximation planar quqdric
	//void setApproxPQuadric( const DPlanarQuadric& pq ) { 
	//	if(!sdata) sdata = new SurfaceData();
	//	sdata->pq = pq;
	//}
	const DVector3d& getCurvatureDirection0() const { return sdata ? sdata->ce0 : DVector3d::zero; }
	const DVector3d& getCurvatureDirection1() const { return sdata ? sdata->ce1 : DVector3d::zero; }
	double getCurvature0() const { return sdata ? sdata->cv0 : 0.0; }
	double getCurvature1() const { return sdata ? sdata->cv1 : 0.0; }
	/// Copy all surface data
	void copySurfaceData(const MeshFace* source);
public:
	/// Checks, whether both points are "visible" from the same side of this face
	bool sameSide(const DPoint3d& ptA, const DPoint3d& ptB) const;
	/// Check whether the normal vector of this face points outside the given block (containing this face)
	bool outsideOrientation(MeshBlock* block) const;
	/// clear some data for faster all-delete process
	virtual void preDeleteAll() {}
	/// enums for inverted-rank
	enum { SP_GOOD = 1, SP_OK = 10, SP_BAD = 100 };
protected:
	/// Number of points/edges for this face
	int				count;
	/// Array of references to edges
	MeshEdge3d**	edges;
	/// Array of references to vertices
	MeshPoint3d**	points;
	/// Array of references to adjacent blocks
	MeshBlock*		blocks[2];
	/// Surface auxiliary data
	SurfaceData*	sdata;
};

#endif // !defined(MESHFACE_H__INCLUDED)
