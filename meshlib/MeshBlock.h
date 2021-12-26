/////////////////////////////////////////////////////////////////////////////
// MeshBlock.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2002-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#if !defined(MESHBLOCK_H__INCLUDED)
#define MESHBLOCK_H__INCLUDED

#include "DRect.h"
#include "DPoint.h"
#include "MeshData.h"
#include "TagExtended.h"

class MeshFace;
class MeshPoint3d;
class Metric3dContext;
class MeshEdge3d;
struct MeshViewBlockData;
class DMatrix3d;

/**
 * This class describes the base, abstract class for 3D mesh elements.
 */
class MeshBlock : public IndexedObject, public TagExtended
{
public:
	/// Standard constructor
	MeshBlock(int fct, int pct, MeshFace** new_faces = nullptr, MeshPoint3d** new_points = nullptr);
	/// Standard destructor
	virtual ~MeshBlock() {}
public:
	virtual bool checkIntTagForAnyPoint(TagExtended::TagType tag_type, int tag_value) const;
	/// clear some data for faster all-delete process
	virtual void preDeleteAll() {}
	/// add this element as a lump of mass (respecting control space)
	virtual void addForInertialCenter(Metric3dContext& mc, DPoint3d& inertial_center, double& total_mass) const;
	/// add this element as a lump of mass (respecting control space)
	virtual void addForInertialMoments(Metric3dContext& mc, const DPoint3d& inertial_center, DMatrix3d& inertial_moments) const;
	/// Sets and returns the metric difference of this element
	virtual double countMetricDiffQuality(Metric3dContext& /*mc*/) { return quality = 1.0; }
	/// Returns "mean ratio" coefficient quality
	virtual double getMeanRatio(Metric3dContext& /*mc*/, bool /*ext_metric*/ = false) const { assert(false); return 1.0; }
	/// Returns the metric difference of this element
	virtual double countMetricDiff(Metric3dContext& /*mc*/) const { return LARGE_NUMBER; }
	/// Return representative points (using param_point_distance)
	void getRepresentativePoints(DataVector<DPoint3d> & points) const;
	/// Return distance (squared) point-triangle (using param_point_distance)
	double distance2(const DPoint3d& pt) const;
	/// Removes tag for all adjacent edges and faces
	void removeTagForEdgesAndFaces(TagExtended::TagType tag_type);
	/// Marks (or clears) tag for all faces of this block
	void setIntTagForFaces(TagExtended::TagType tag_type, int t = 1);
	/// Checks validity of switching point of this block to the other one together with adjacent edges and faces
	virtual bool checkSwitchPointsWithFaces(const MeshPoint3d *point1, MeshPoint3d *point2) const;
	/// Switches point of this block to the other one together with adjacent edges and faces
	virtual void switchPointsWithFaces(const MeshPoint3d *point1, MeshPoint3d *point2);
	/// Switches point of this block to the other one together with adjacent edges and faces
	virtual void switchPointsWithFacesBoundary(const MeshPoint3d *point1, MeshPoint3d *point2);
	/// Returns the size of this element (can be calculated in modified metric)
	virtual double getVolume(Metric3dContext& mc, bool local_metric = true) const = 0;
	/// Returns the size of this element (with non-altered coordinates)
	virtual double getVolumeNoMetric() const = 0;
	/// Create block-faces connections
	virtual void attachToFaces() = 0; 
	/// Removes block-faces connections
	virtual void detachFromFaces() = 0;
	/// Returns the reference to the face (of this block) incident to the given one via the given edge
	MeshFace* getIncidentFace(MeshFace* face, MeshEdge3d* edge);
	/// Returns the reference to the block incident through the i-th face(can be nullptr)
	MeshBlock* getNeighbour(int i) const;
	/// Returns the middle point for this element
	DPoint3d getMiddlePoint() const;
	/// Compares this block with the given one (required by DataContainer if using heap-ordering)
	short compareTo(MeshBlock* item);
	/// Returns the pointer to one of vertices of this block
	MeshPoint3d* getPoint(int i) const { return points[i]; }
	/// Returns the pointer to ith face of this block
	MeshFace*	 getFace(int i) const { return faces[i]; }
	/// Returns the number of faces in this block
	int getFaceCount() const { return face_count; }
	/// Returns the number of faces in this block (of given type)
	int getFaceCount(int ect) const;
	/// Returns the number of edges in this block
	virtual int getEdgeCount() const { assert(false); return 0; }
	/// Returns true, if this block contains face
	bool isAdjacentTo(const MeshFace* face) const;
	/// Return i-th edge
	virtual MeshEdge3d* getEdge(int /* i */) const { assert(false); return nullptr; }
	/// Return indices of points within this block for the i-th edge
	virtual void getEdgeIndices(int /* i */, int& /* i1 */, int& /* i2 */) const { assert(false); }
	/// Returns the number of vertices in this block
	int getPointCount() const { return point_count; }
	/// Returns the description of this block for visualization
	virtual std::shared_ptr<MeshViewBlockData> getViewData(double /* shrink_ratio */ = 1.0) const 
		{ assert(false); return nullptr; }
	/// Returns the type of this block (should be implemented in derived classes)
	virtual ElementType getType() const { return BLOCK_UNKNOWN; }
	/// Checks whether the block is inverted (i.e. bad orientation or invalid geometry)
	virtual bool isInverted() const { assert(false); return false; }
	/// Returns the bounding box for this block
	virtual DBox getBoundingBox() const;
	/// Returns the minimum volume for this element (in unitary space)
	double getMinVolume() const { return 0.01; }
	/// Returns the area-ID of element
	int getAreaID() const { return area_id; }
	/// Sets the area-ID of element
	void setAreaID(int id);
	/// Checks wheterh the i-th face of the block belongs to the boundary
	bool isBorder(int i) const;
	/// Returns the current quality of element (bufferred or recalculated)
	double getQuality() const { return quality; }
	/// Sets the quality value for this element (buffered for efficiency)
	void setQuality(double q){ quality = q; }
	/// Sets and returns the quality of this element (with globally or explicitely selected criterion)
	virtual double countQuality(Metric3dContext& /* mc */, int /* criterion */ = -1) { 
		return quality = 0.0;
	}
	/// Returns true if any adjacent face is marked as border
	bool hasBoundaryFace() const;
	/// Returns true if any adjacent edge is marked as border
	bool hasBoundaryEdge() const;
	/// Returns true if any adjacent vertex is marked as border
	bool hasBoundaryVertex() const;
	/// Returns point opposite the given face
	virtual MeshPoint3d* getOppositePoint(const MeshFace* mf) const;
	/// Creates and returns a copy of this block (+ whole connectivity)
	virtual MeshBlock* clone() const { return nullptr; }
protected:
	/// Number of faces
	int face_count;
	/// Array of pointers to faces 
	MeshFace	**faces;
	/// Number of vertices
	int point_count;
	/// Array of pointers to vertices
	MeshPoint3d **points;
	/// Area identifier for this element
	int		area_id;
	/// Buffered quality coefficient (invalid if < 0.0)
	double	quality;
public:
	static int param_point_distance;
};

#endif // !defined(MESHBLOCK_H__INCLUDED)
