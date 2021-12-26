/////////////////////////////////////////////////////////////////////////////
// MeshTetrahedron.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2002-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#if !defined(MESHTETRAHEDRON_H__INCLUDED)
#define MESHTETRAHEDRON_H__INCLUDED

#include "MeshBlock.h"
#include "MeshPoint3d.h"
#include "MeshFace.h"

class MeshEdge3d;
struct MeshViewBlockData;

/**
 * This class implements a three-dimensional, tetrahedral mesh element
 *	with several required methods.
 */
class MeshTetrahedron : public MeshBlock  
{
public:
	/// Standard constructor
	MeshTetrahedron(MeshPoint3d *p1, MeshPoint3d *p2, MeshPoint3d *p3, MeshPoint3d *p4);
	/// Standard desctructor
	virtual ~MeshTetrahedron();
public:
	/// Creates and returns a copy of this block (+ whole connectivity)
	virtual MeshBlock* clone() const override;
	/// Create block-faces connections
	virtual void attachToFaces(); 
	/// Removes block-faces connections
	virtual void detachFromFaces();
	/// clear some data for faster all-delete process
	virtual void preDeleteAll();
	/// add this element as a lump of mass (respecting control space)
	static void addForInertialCenterAdaptiveSplitLongestEdge(Metric3dContext& mc, 
		const DPoint3d& pt0, const DPoint3d& pt1, const DPoint3d& pt2, const DPoint3d& pt3, 
		DPoint3d& inertial_center, double& total_mass, int level = 0);
	/// add this element as a lump of mass (respecting control space)
	virtual void addForInertialCenter(Metric3dContext& mc, DPoint3d& inertial_center, double& total_mass) const;
	/// add this element as a lump of mass (respecting control space)
	static void addForInertialMomentsAdaptiveSplitLongestEdge(Metric3dContext& mc, 
		const DPoint3d& pt0, const DPoint3d& pt1, const DPoint3d& pt2, const DPoint3d& pt3, 
		const DPoint3d& inertial_center, DMatrix3d& inertial_moments, int level = 0);
	/// add this element as a lump of mass (respecting control space)
	virtual void addForInertialMoments(Metric3dContext& mc, const DPoint3d& inertial_center, DMatrix3d& inertial_moments) const;
	/// Check whether face (three points) are consistent with this tetrahedron
	bool properOrientation(const MeshPoint3d* mpt0, const MeshPoint3d* mpt1, const MeshPoint3d* mpt2) const;
	/// Sets and returns the metric difference of this element
	virtual double countMetricDiffQuality(Metric3dContext& mc); 
	/// Returns the metric difference of this element
	virtual double countMetricDiff(Metric3dContext& mc) const;
	/// Returns "mean ratio" coefficient quality
	virtual double getMeanRatio(Metric3dContext& /*mc*/, bool /*ext_metric*/ = false) const;
	/// Returns the size of this element (can be calculated in modified metric)
	virtual double getVolume(Metric3dContext& mc, bool local_metric = true) const;
	/// Returns the size of this element (with non-altered coordinates)
	virtual double getVolumeNoMetric() const;
	/// Return indices of points within this block for the i-th face
	void getFaceIndices(int i, int& i1, int& i2, int& i3) const;
	/// Return indices of points within this block for the i-th edge
	void getEdgeIndices(int i, int& i1, int& i2) const;
	/// Returns the number of edges in this block
	int getEdgeCount() const { return 6; }
	/// Return i-th edge
	MeshEdge3d* getEdge(int i) const;
	/// Returns face opposite the given point
	MeshFace* getOppositeFace(const MeshPoint3d* mp) const;
	/// Returns point opposite the given face
	MeshPoint3d* getOppositePoint(const MeshFace* mf) const;
	/// Returns fourth point
	MeshPoint3d* getOtherPoint(const MeshPoint3d* pt0, const MeshPoint3d* pt1, const MeshPoint3d* pt2) const;
	/// Returns third and fourth points
	bool getOtherPoints(const MeshPoint3d* pt0, const MeshPoint3d* pt1, MeshPoint3d* pts[]) const;
	/// Checks validity of switching point of this block to the other one together with adjacent edges and faces
	virtual bool checkSwitchPointsWithFaces(const MeshPoint3d *point1, MeshPoint3d *point2) const;
	/// Switches point of this tetrahedron to the other one together with adjacent edges and faces
	virtual void switchPointsWithFaces(const MeshPoint3d *point1, MeshPoint3d *point2);
	/// Switches point of this tetrahedron to the other one together with adjacent edges and faces
	virtual void switchPointsWithFacesBoundary(const MeshPoint3d *point1, MeshPoint3d *point2);
	/// Counts (and stores) the quality of the tetrahedron (for global or explicite criterion)
	virtual double countQuality(Metric3dContext& mc, int criterion = -1);
	/// Checks whether the given mesh point is withing the outer sphere of this tetrahedron (different metric may be used)
	bool isPointInOuterSphere(Metric3dContext& mc, MeshPoint3d *point, bool local_metric = true) const 
		{ return isPointInOuterSphere(mc, point->getMetricCoordinates(mc), local_metric); }
	/// Checks whether the given point is withing the outer sphere of this tetrahedron (different metric may be used)
	bool isPointInOuterSphere(Metric3dContext& mc, const DMPoint3d& point, bool local_metric = true) const;
	/// Checks whether this tetrahedron contains the given point
	bool isPointInside(const DPoint3d & pt0) const;
	/// Checks whether this tetrahedron is inverted (wrongly orientated)
	virtual bool isInverted() const;
	/// Returns the description of this block for visualization
	virtual std::shared_ptr<MeshViewBlockData> getViewData(double shrink_ratio = 1.0) const;
	/// Returns the reference to neighbouring tetrahedron in the direction given by the point (from the middle of the face)
	MeshTetrahedron* getNeighbourInFaceDirection(const DPoint3d &mpt, bool mind_border = false, bool second_best = false);
	/// Returns the reference to neighbouring tetrahedron in the direction given by the point (from the middle of the block)
	MeshTetrahedron* getNeighbourInBlockDirection(const DPoint3d &point, bool mind_border = false, bool second_best = false);
	/// Returns the tetrahedron containing the given point, starting from this one
	MeshTetrahedron* findTetrahedronByNeighbours(const DPoint3d &point, bool mind_border = false);
	/// Returns the length of the radius of the outer sphere for this tetrahedron(different metric may be used)
	double getOuterSphereRadius(Metric3dContext& mc, bool local_metric = true) const;
	/// Returns the center of the outer sphere for this tetrahedron(different metric may be used)
	DMPoint3d getOuterSphereCenter(Metric3dContext& mc, bool local_metric = true) const;
	/// Returns the object-specific type
	virtual ElementType getType() const { return BLOCK_TETRA; }
	/// Returns the optimal volume of tetrahedra ( |edge| = 1.0 )
	static double getOptimumVolume() { return 0.11785; } // ~ sqrt(2.0)/12.0
	/// Returns the optimal radius length (squared) of the circumscribed sphere  ( |edge| = 1.0 )
	static double getOptimumCircumradius2() { return 0.375; }
	/// Calculates the minimum dihedral angle (sin, for 0-90o)
	double getMinDihedralAngleSin() const;
public:
	/// Method of assessing the quality of tetrahedron
	static int param_quality_criterion;
	/// Method of metric selection for assessing the quality of triangle
	static int param_quality_metric_selection;
	/// Height of perfect unitary tetrahedron;
	static const double opt_metric_height;
	/// Last tetrahedra-traverse counter
	static unsigned int traverse_counter;
public:
	/// Available sequences of vertices
	static const unsigned char PCONF[4][4];
	/// Edge to vertices relation
	static const unsigned char ECONF[6][2];
protected:
	/// Four-element array of mesh faces
	MeshFace*		face_array[4];
	/// Four-element array of vertices (mesh points)
	MeshPoint3d*	point_array[4];
};

#endif // !defined(MESHTETRAHEDRON_H__INCLUDED)
