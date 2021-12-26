/////////////////////////////////////////////////////////////////////////////
// ControlSpace2dQuadTree.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2005-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#if !defined(CONTROLSPACEQUADTREE_H__INCLUDED)
#define CONTROLSPACEQUADTREE_H__INCLUDED

#include "ControlSpace2dAdaptive.h"
#include "DataMatrix.h"
#include "DataVector.h"
#include "DataHashTable.h"
#include "DataContainer.h"

class MeshContainer2d;
class MeshEdge2d;
class EPSFile;
class ReparameterizationData;
class ControlSpace3d;

/**
 * This class implements a quadtree-based control space.
 */
class ControlSpace2dQuadTree : public ControlSpace2dAdaptive
{
public:
	/// Standard contructor
	ControlSpace2dQuadTree(SurfaceConstPtr surface, const DRect& box, 
		int nxy = 100);
	/// Destructor
	virtual ~ControlSpace2dQuadTree();
public:
	/// Returns the "screenshot" of this control space for visualization
	virtual MeshViewSet* getViewSet(MeshViewSet* set = nullptr, bool with_surface = false, 
		TagExtended::TagType tag = TagExtended::TAG_NONE) const;
	/// Returns the type of element
	virtual int getType() const { return MeshData::CONTROL_QUADTREE; }
	/// Adds new information (size and shape) for some point within the domain
	virtual void addControlNode(const ControlNode2d& node);
	/// Whether reparameterization is used in some subarea (point) of control space
	virtual ReparameterizationData* getReparameterization(const DPoint2d& pt) const;
	/// Return metric gradation ratio at point
	virtual double getMetricGradationRatio(const DPoint2d& pt) const;
	/// Get sizing info (matrix mode) at the given point
	virtual ControlDataMatrix2d getMetricAtPoint(const DPoint2d& pt) const;
	/// Get transformation and parameterization matrices at the given point
	virtual ControlDataMatrix2d getMetricAndParameterizationAtPoint(const DPoint2d& pt, ControlDataMatrix2d& p_cdm) const;
	/// Interpolates (initializes) the control space basing on the previously given sizing information in discrete points
	virtual bool interpolate();
	/// Draws the structure of the control space
	virtual void storeEPS(const char* name = "control", int id = 0);
	/// Initializes the control space basing on the curvature of surface and boundary
	virtual void setSurfaceCurvatureControlData();
	/// Refines control space at the given point
	virtual bool setMinControl(const DPoint2d& pt, const ControlDataMatrix2d& cdm, bool min_value_set = true);
	/// Returns number of control nodes in adaptive control structure
	virtual int getControlNodesCount() const override;
	/// Invoke function for all control nodes of this space (read-only)	
	virtual void forEachControlNode(const std::function<void(const ControlNode2d& node)>& fg) const override;
	/// Invoke function for all control nodes of this space
	virtual void forEachControlNode(const std::function<void(ControlNode2d& node)>& fg) override;
	/// Refines control space to parameterization variance
	virtual void adaptToParameterization();
	/// Refines control space with respect to reference 3d control space
	virtual bool adaptToControlSpace3d(CS3dConstPtr space) override;
	/// Log basic information about this control space
	virtual void logDescription() const;
	/// Smoothen variance of metric within the control space
	virtual bool smoothen();
	/// Return interpolated value of extended tag data from control nodes at some given point;
	virtual double interpolateDoubleTag(const DPoint2d& pt, TagExtended::TagType type) const;
	/// Set given tag value as minimum at the closest control element
	virtual void setMinDoubleTag(const DPoint2d& pt, TagExtended::TagType type, double q);
	/// Set given tag value as minimum at the closest control element
	virtual void setMaxDoubleTag(const DPoint2d& pt, TagExtended::TagType type, double q);
	/// Calculate inside information for nodes/elements
	virtual void markInsideNodes(const MeshContainer2d* boundary_mesh);
	/// Returns local resolution
	virtual double getLocalResolution(const DPoint2d& pt) const;
public:
	void testQuadTree();
public:
	/// maximum depth for quadtree
	static int param_max_depth;
	/// maximum depth for quadtree
	static int param_max_depth_surface_adaptation;
public:
	typedef DataVector<ControlNode2d*> QuadVertexList;
	typedef DataVector<ControlNode2d> QuadVertexArray;
	enum InsideMarkEdgeForbidded { 
		FORBIDDEN_EDGE_XL = 1, FORBIDDEN_EDGE_XH = 2, 
		FORBIDDEN_EDGE_YL = 4, FORBIDDEN_EDGE_YH = 8};
	struct ControlPrediction{
		ControlPrediction() {
			for(int i = 0; i < 3; i++)
				value[i] = 0.0;
			for(int i = 0; i < 5; i++)
				counter[i] = 0;
		}
		double value[3]; // (min/max/ave) of ave-arytm
		int counter[5]; // how many leaves of each type (i.e. number of in-material nodes)
	};
	/// Show CS-grid with inside info for debug
	void showMarkInside(const MeshContainer2d* boundary_mesh = nullptr) const;
	/// Calculate mesh prediction using CS
	ControlPrediction calculateControlPrediction() const;
	/// Class QuadLeaf for QuadTree structure
	class QuadLeaf{
		friend class ControlSpace2dQuadTree;
	public:
		QuadLeaf();
		QuadLeaf(double mid_x, double mid_y, double dx, double dy, int level = 0);
		~QuadLeaf();
		/// Whether reparameterization is used in some subarea (point) of control space
		ReparameterizationData* getReparameterization(const DPoint2d& pt) const;
		/// Return metric gradation ratio at point
		double getMetricGradationRatio(const DPoint2d& pt) const;
		/// Get sizing info (matrix mode) at the given point
		ControlDataMatrix2d getMetricAtPoint(const DPoint2d& pt) const;
		/// Get sizing info (matrix mode) at the given point
		ControlDataMatrix2d getMetricAndParameterizationAtPoint(const DPoint2d& pt, 
			ControlDataMatrix2d& p_cdm, SurfaceConstPtr surface) const;
		/// Set values of vertices indices
		void setVertices(ControlNode2d* v_sw, ControlNode2d* v_se, ControlNode2d* v_nw, ControlNode2d* v_ne){
			m_vertices[VX2D_SW] = v_sw;
			m_vertices[VX2D_SE] = v_se;
			m_vertices[VX2D_NW] = v_nw;
			m_vertices[VX2D_NE] = v_ne;
		}
		/// Set pointers to neighbouring leaves
		void setNeighbours(QuadLeaf* nd, QuadLeaf* nr, QuadLeaf* nu, QuadLeaf* nl){
			m_neighbours[NB_DOWN] = nd;
			m_neighbours[NB_RIGHT] = nr;
			m_neighbours[NB_UP] = nu;
			m_neighbours[NB_LEFT] = nl;
		}
		/// Set pointer to neighbouring leaf for this leaf and sub-leaves
		void setNeighbour(int side, QuadLeaf* nb, bool skip_first = false);
		/// Set coordinates (middle point and size)
		void setCoordinates(const DPoint2d& middle, const DVector2d& dl){
			m_middle = middle;
			m_dl = dl;
		}
		/// Adapt quad leaf to surface curvature (split if required)
		bool checkAndSetReparameterization(SurfaceConstPtr surface);
		/// Adapt quad leaf to surface curvature (split if required)
		void adaptToSurfaceCurvature(SurfaceConstPtr surface, 
			QuadVertexArray& grid_vertices,
			double c_ratio, double stretch_max_ratio,
			double model_diameter, double min_len, double max_len);
		/// Adapt quad leaf to parameterization (split if required)
		void adaptToParameterization(SurfaceConstPtr surface, 
			QuadVertexArray& grid_vertices);
		/// Adapt quad leaf with respect to reference 3d control space (split if required)
		bool adaptToControlSpace3d(SurfaceConstPtr surface, 
			QuadVertexArray& grid_vertices, CS3dConstPtr space);
		/// Adapt quad leaf to new source point and insert this point into proper collection
		void adaptAndInsertControlPoint(QuadVertexArray& grid_vertices, 
			const ControlNode2d& qv, SurfaceConstPtr surface);
		/// Adapt quad leaf to new source point and set minimum value
		bool adaptAndSetControlPoint(QuadVertexArray& grid_vertices, 
			const ControlNode2d& qv, SurfaceConstPtr surface, bool min_value_set = true);
		/// Split this leaf (+ build relations with neighbours)
		bool split(QuadVertexArray& grid_vertices, SurfaceConstPtr surface, 
			QuadVertexList* split_vertices = nullptr, bool init_values = false);
		/// Whether leaf is already split
		bool isSplit() const { return m_leaves[0] != nullptr; }
		/// Adjust levels of adjacent leaves
		void balance(QuadVertexArray& grid_vertices, SurfaceConstPtr surface, 
			QuadVertexList* split_vertices, bool init_values = false);
		/// Store leaf structure into EPS file
		void storeEPS(EPSFile* eps_file) const;
		/// Insert view data about leaf structure
		void addToViewSet(MeshViewSet* set, SurfaceConstPtr surface, 
			TagExtended::TagType tag = TagExtended::TAG_NONE) const;
		/// Find last leaf (not split) for givent point
		QuadLeaf* findLastLeaf(const DPoint2d& pt);
		/// Gather metric value from source nodes within leaf	
		void gatherDataFromSourcePoints(SurfaceConstPtr surface);
		/// Propagate control values from within leaves
		int propagateUp();
		/// Propagate control values into leaves
		int propagateDown();
		/// Returns i-th vertex of the leaf
		ControlNode2d* getVertex(int i) const { return m_vertices[i]; }
		/// Smoothen variance of metric for leaf nodes
		int smoothenMetricAtNodes(bool uptree = true);
		/// count statistics
		void statQuadTree(int & max_level, int & leaf_counter, int level = 1) const;
		/// Return interpolated value of extended tag data from control nodes at some given point;
		double interpolateDoubleTag(const DPoint2d& pt, TagExtended::TagType type) const;
		/// Propagate control values from within / into leaves (returns true if whole leaf is ready)
		bool propagateInsideMark();
		/// Propagate control values from within / into leaves (returns true if whole leaf is ready)
		bool finishInsideMark();
		/// mark inside infomation for nodes
		int markInsideNodes(const DPoint2d& pt0, const DPoint2d& pt1, MeshEdge2d* edge);
		/// mark inside infomation for nodes
		void markInsideNodesY(const DPoint2d& pt0, const DPoint2d& pt1, MeshEdge2d* edge, bool right_line);
		/// mark inside infomation for nodes
		void markInsideNodesX(const DPoint2d& pt0, const DPoint2d& pt1, MeshEdge2d* edge, bool top_line);
		/// Calculate mesh prediction using CS
		void calculateControlPrediction(ControlPrediction & cp) const;
		/// Calculate mesh prediction using CS
		void calculateControlPrediction(ControlPrediction & cp, int ind, int inside_count, int inside[]) const;
	public:
		/// whether to use midnodes for value interpolation
		static int param_use_midnodes;
		/// balance level
		static int param_balance_level;
	private:
		/// middle point of leaf (parametric)
		DPoint2d m_middle;
		enum QuadVertexWhich {VX2D_SW = 0, VX2D_SE = 1, VX2D_NW = 2, VX2D_NE = 3,
				VX2D_FIRST = 0, VX2D_LAST = 3};
		ControlNode2d* m_vertices[4]; // only indices, 2-3 / 0-1
		enum QuadLeafWhich {LF_SW = 0, LF_SE = 1, LF_NW = 2, LF_NE = 3,
				LF_FIRST = 0, LF_LAST = 3};
		QuadLeaf* m_leaves[4];
		/// Parameterization matrix for leaf middle point
		ControlDataMatrix2d m_middle_cdm;
		/// Reparameterization data (if required)
		ReparameterizationData* m_reparameterization;
		//-- extra data
		/// depth level
		int m_level;
		/// dimensions of leaf (parametric)
		DVector2d m_dl;
		/// dimensions of leaf (real)
		DVector2d m_real_dl;
		enum QuadNeighbourWhich {NB_DOWN = 0, NB_RIGHT = 1, NB_UP = 2, NB_LEFT = 3, 
			NB_FIRST = 0, NB_LAST = 3};
		QuadLeaf* m_neighbours[4]; 
		enum QuadMidVertexWhich {MVX2D_DOWN = 0, MVX2D_RIGHT = 1, MVX2D_UP = 2, MVX2D_LEFT = 3,
			MVX2D_FIRST = 0, MVX2D_LAST = 3};
		static QuadVertexWhich mid_to_nodes[4][2];
		static QuadMidVertexWhich nodes_to_mid[4][2];
		ControlNode2d* m_mid_vertices[4]; 
		//--- point sources before interpolation
		DataVector<ControlNode2d> *m_source_points;
	};
protected:
	/// Translates the point-coordinates into index of control matrix
	void countLocalCoordinates(const DPoint2d& pt, int& ix, int& iy) const;
	/// Translates the point-coordinates into (consolidated) index of control matrix
	int countLocalCoordinates(const DPoint2d& pt) const;
	/// Returns smalles QuadTree leaf for this point
	QuadLeaf* findLastLeaf(const DPoint2d& pt) const;
	/// Extrapolate control data from nodes
	int extrapolateMainControlNodes();
	/// Propagate inside mark for the main grid
	void propagateMainGridInsideMark();
protected:
	/// Number of columns in the matrix
	int m_nx;
	/// Number of rows in the matrix
	int m_ny;
	/// Values stored in matrix
	QuadLeaf* m_grid;
	/// Grid vertices with defined control data
	QuadVertexArray m_grid_vertices;
};

#endif // !defined(CONTROLSPACEQUADTREE_H__INCLUDED)
