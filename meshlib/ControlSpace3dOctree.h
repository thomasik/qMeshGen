/////////////////////////////////////////////////////////////////////////////
// ControlSpace3dOctree.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2006-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#if !defined(CONTROLSPACE3DOCTREE_H__INCLUDED)
#define CONTROLSPACE3DOCTREE_H__INCLUDED

#include "ControlSpace3dAdaptive.h"
#include "DataPtrMatrix.h"
#include "MeshData.h"
#include "DataList.h"
#include "DataHashTable.h"
#include "ControlSpace3dKdTree.h"

class MeshContainer3d;
class MeshFace;
//class MeshEdge2d;
//class EPSFile;

/**
 * This class implements an octree-based control space.
 */
class ControlSpace3dOctree : public ControlSpace3dKdTree
{
public:
	/// Standard contructor
	ControlSpace3dOctree(const DBox& box, int nxyz = 100);
	/// Destructor
	virtual ~ControlSpace3dOctree();
public:
	/// Insert additional info -> about local surface available for the given point
	virtual bool addLocalSurfaceAtPoint(const DPoint3d& pt, SurfaceConstPtr local_surface, 
		DataHashTableKeyValue<SurfaceSetPtr, SurfaceSetPtr> & local_surface_hash) override;
	/// Insert additional info -> about local surface available for the given bounding box
	virtual bool addLocalSurfaceAtBBox(const DBox& box, SurfaceConstPtr local_surface, 
		DataHashTableKeyValue<SurfaceSetPtr, SurfaceSetPtr> & local_surface_hash) override;
	/// Returns additional info -> about local surfaces (set) available for the given point
	virtual SurfaceSetConstPtr getLocalSurfaceSetAtPoint(const DPoint3d& pt) override;
	/// Returns the "screenshot" of this control space for visualization
	virtual MeshViewSet* getViewSet(MeshViewSet* set = nullptr, bool with_points = false) const override;
	/// Returns the type of element
	virtual int getType() const  override { return MeshData::CONTROL_OCTREE_3D; }
	/// Adds new information (size and shape) for some point within the domain
	virtual void addControlNode(const ControlNode3d& node) override;
	/// Get sizing info (matrix mode) at the given point
	virtual ControlDataMatrix3d getMetricAtPoint(const DPoint3d& pt) const override;
	/// Interpolates (initializes) the control space basing on the previously given sizing information in discrete points
	virtual bool interpolate() override;
	/// Refines control space at the given point
	virtual bool setMinControl(const DPoint3d& pt, const ControlDataMatrix3d& cdm, bool min_value_set = true) override;
	using ControlSpace3dAdaptive::setMinControl;
	/// Returns number of control nodes in adaptive control structure
	virtual int getControlNodesCount() const override;
	/// Invoke function for all control nodes of this space (read-only)
	virtual void forEachControlNode(const std::function<void(const ControlNode3d& node)>& fg) const override;
	/// Invoke function for all control nodes of this space
	virtual void forEachControlNode(const std::function<void(ControlNode3d& node)>& fg) override;
	/// Log basic information about this control space
	virtual void logDescription() const override;
	/// Smoothen variance of metric within the control space
	virtual bool smoothen() override;
	/// Return metric gradation ratio at point
	virtual double getMetricGradationRatio(const DPoint3d& pt) const override;
	/// optimize for size (lower memory requirements) - but also from now readonly !
	virtual void compact() override;
	/// Calculate inside information for nodes/elements
	virtual void markInsideNodes(const MeshContainer3dSurface* surface_mesh) override;
	/// Mark non-tagged boundary nodes as empty
	void markInsideNodesAtBoundary();
	/// Returns local resolution
	virtual double getLocalResolution(const DPoint3d& pt) const override;
	/// Return interpolated value of extended tag data from control nodes at some given point;
	virtual double interpolateDoubleTag(const DPoint3d& pt, TagExtended::TagType type) const override;
	/// Set given tag value as minimum at the closest control element
	virtual void setMinDoubleTag(const DPoint3d& /* pt */, TagExtended::TagType /* type */, double /* q */) override;
	/// Set given tag value as minimum at the closest control element
	virtual void setMaxDoubleTag(const DPoint3d& /* pt */, TagExtended::TagType /* type */, double /* q */) override;
public: // for kdtree-interface
	virtual int adaptToField(const std::function<CDM3d(const DPoint3d & pt)> & f, double max_error) override {
		assert(false); return 0; }
	// stat	
	virtual int getElementCount() const override;
	virtual int getValueCount() const override { return getControlNodesCount(); }
	virtual int getMaxDepth() const override;
	virtual int getMaxBalance() const override;
	virtual int getTotalBytes() const override;
	/// Returns descriptive label
	virtual string getTreeType() const override { return "octree-v-oryg"; }
public:
	/// Store octree grid to text file
	bool storeTXT(const char* fname);
	/// Load octree grid from text file
	bool loadTXT(const char* fname);
	/// Test integrity of this structure;
	void testOctree();
	/// If is valid
	bool valid() const;
public:
	/// maximum depth for octree
	static int param_max_depth;
public:
	class OctLeaf;
	typedef std::shared_ptr<ControlNode3d> CN3dPtr;
	// Extra data
	struct OctLeafExtra{
	public:
		/// constructor
		OctLeafExtra(int level = 0);
	public:
		/// leaf level
		int m_level;
		/// neighbours
		OctLeaf* m_neighbours[6]; 
		/// mid vertices for edges
		CN3dPtr m_midedge_vertices[12];
		/// mid vertices for faces
		CN3dPtr m_midface_vertices[6];
		//--- point sources before interpolation
		DataVector<std::shared_ptr<ControlNode3d>> *m_source_points;
	};
public:
	typedef DataVector<CN3dPtr> OctVertexList;
	typedef DataCompoundList<OctLeafExtra> OctExtraList;
	enum InsideMarkEdgeForbidded { 
		FORBIDDEN_EDGE_XL = 1,  FORBIDDEN_EDGE_XH = 2, 
		FORBIDDEN_EDGE_YL = 4,  FORBIDDEN_EDGE_YH = 8,
		FORBIDDEN_EDGE_ZL = 16, FORBIDDEN_EDGE_ZH = 32
	};
	enum InsideMarkEdges {
		EDGE_LSNW = 1,    EDGE_LSNE = 2, 
		EDGE_LSWE = 4,    EDGE_LNWE = 8,
		EDGE_HSNW = 16,   EDGE_HSNE = 32, 
		EDGE_HSWE = 64,   EDGE_HNWE = 128,
		EDGE_LHSW = 256,  EDGE_LHSE = 512,
		EDGE_LHNW = 1024, EDGE_LHNE = 2048
	}; 
	struct ControlPrediction{
		ControlPrediction() {
			for(int i = 0; i < 3; i++)
				value[i] = 0.0;
			for(int i = 0; i < 9; i++)
				counter[i] = 0;
		}
		double value[3]; // (min/max/ave) of ave-arytm
		int counter[9]; // how many leaves of each type (i.e. number of in-material nodes)
	};
	/// Show CS-grid with inside info for debug
	void showMarkInside(const DataVector<MeshFace*>& bfaces) const;
	/// Calculate mesh prediction using CS
	ControlPrediction calculateControlPrediction() const;
	// Class OctLeaf for Octree structure
	class OctLeaf{
	public:
		OctLeaf();
		OctLeaf(const DPoint3d& middle, const DVector3d& dl, int level = 0);
		~OctLeaf();
		/// Insert view data about leaf structure
		void addToViewSet(MeshViewSet* set) const;
		/// Get sizing info (matrix mode) at the given point
		ControlDataMatrix3d getMetricAtPoint(const DPoint3d& pt) const;
		/// Set values of vertices indices
		void setSingleVertex(int side, CN3dPtr vert);
		/// Set pointers to neighbouring leaves
		void setSingleNeighbour(int side, OctLeaf* nb){ 
			assert(m_extra); 
			m_extra->m_neighbours[side] = nb; 
		}
		/// Set pointer to neighbouring leaf for this leaf and sub-leaves
		void setNeighbour(int side, OctLeaf* nb, bool skip_first = false);
		/// Set coordinates (middle point and size)
		void init(const DPoint3d& middle, const DVector3d& dl, int level = 0);
		/// Return metric gradation ratio at point
		double getMetricGradationRatio(const DPoint3d& pt) const;
		/// Adapt quad leaf to new source point and insert this point into proper collection
		void adaptAndInsertControlPoint(OctVertexList& grid_vertices, std::shared_ptr<ControlNode3d> qv);
		/// Adapt quad leaf to new source point and set minimum value
		bool adaptAndSetControlPoint(OctVertexList& grid_vertices, 
			std::shared_ptr<ControlNode3d> qv, bool min_value_set = true);
		/// Split this leaf (+ build relations with neighbours)
		bool split(OctVertexList& grid_vertices, OctVertexList& split_vertices, bool init_values = false, bool with_balancing = true);
		/// Whether leaf is already split
		bool isSplit() const { return m_leaves[0] != nullptr; }
		/// Adjust levels of adjacent leaves
		void balance(OctVertexList& grid_vertices, OctVertexList& split_vertices, bool init_values = false);
		/// Find last leaf (not split) for givent point
		OctLeaf* findLastLeaf(const DPoint3d& pt);
		/// Gather metric value from source nodes within leaf	
		void gatherDataFromSourcePoints();
		/// Propagate control values from within leaves
		int propagateUp();
		/// Propagate control values into leaves
		int propagateDown();
		/// Returns i-th vertex of the leaf
		CN3dPtr getVertex(int i) const { return m_vertices[i]; }
		/// Smoothen variance of metric for leaf nodes
		int smoothenMetricAtNodes(bool uptree = true);
		/// Set mid-edge vertex for selected edge
		bool setMidEdgeVertex(CN3dPtr q0, CN3dPtr q1, CN3dPtr node);
		/// Check validity
		bool valid(int depth = 0) const;
		/// Stores topology (1-split, 0-not split) info about this leaf
		void storeTxtSplit(ostream& os) const;
		/// Loads topology (split/not-split)
		void loadTxtSplit(istream& is, OctVertexList& grid_vertices);
		/// count statistics
		void statOctree(int & max_level, int & leaf_counter, int level = 1) const;
		/// optimize for size (lower memory requirements)
		bool optimize(int & leaves_total, int & leaves_obsolete);
		/// Propagate control values from within / into leaves (returns true if whole leaf is ready)
		bool propagateInsideMark();
		/// Set control values for any remaining unset node
		bool setInsideMark(int value);
		/// mark inside infomation for nodes
		int markInsideNodes(const DPoint3d& pt0, const DPoint3d& pt1, const DPoint3d& pt2, 
			MeshFace* face, const double& MIN_DL);
		/// mark inside infomation for nodes
		void markInsideNodesXY(const DPoint3d& pt0, const DPoint3d& pt1, const DPoint3d& pt2, 
			MeshFace* face, bool plane_z1, const double& MIN_DL);
		/// mark inside infomation for nodes
		void markInsideNodesXZ(const DPoint3d& pt0, const DPoint3d& pt1, const DPoint3d& pt2, 
			MeshFace* face, bool plane_Y1, const double& MIN_DL);
		/// mark inside infomation for nodes
		void markInsideNodesYZ(const DPoint3d& pt0, const DPoint3d& pt1, const DPoint3d& pt2, 
			MeshFace* face, bool plane_X1, const double& MIN_DL);
		/// Calculate mesh prediction using CS
		void calculateControlPrediction(ControlPrediction & cp) const;
		/// Calculate mesh prediction using CS
		void calculateControlPrediction(ControlPrediction & cp, int ind, int inside_count, int inside[]) const;
		/// mark inside-information for boundary CS-leaves
		void markInsideNodesAtBoundary(OctFaceWhich side);
		/// draw octree element
		void drawToViewSet(MeshViewSet* set, int part = -2, double ratio = 1.0) const;
		/// Return interpolated value of extended tag data from control nodes at some given point;
		double interpolateDoubleTag(const DPoint3d& pt, TagExtended::TagType type) const;
		/// Return local surface set
		SurfaceSetPtr getLocalSurfaceSet() const { return m_local_surface_set ; }
		/// Change local sruface set
		void setLocalSurfaceSet( SurfaceSetPtr surf_set ) { m_local_surface_set = surf_set; }
		/// Insert additional info -> about local surface available for the given bounding box
		void addLocalSurfaceAtBBox(const DBox& box, SurfaceConstPtr local_surface, 
			DataHashTableKeyValue<SurfaceSetPtr, SurfaceSetPtr> & local_surface_hash);
		/// Insert the given local surface into surface set (false -> already there)
		bool addToLocalSurfaceSet(SurfaceConstPtr local_surface, 
			DataHashTableKeyValue<SurfaceSetPtr, SurfaceSetPtr> & local_surface_hash);
		int getElementCount() const;	
		int getMaxDepth() const;
		int getMaxDepthAndBalance(int & max_balance) const;
	protected:
		/// forbidden edges data
		static const int FEDATA[12][4];
	public:
		/// for optimization of octree size
		static double param_octree_reduce_threshold;
		/// balance level
		static int param_balance_level;
	private:
		/// middle point of leaf
		DPoint3d m_middle;
		/// dimensions of leaf
		DVector3d m_dl;
		/// (pointers to) vertices
		CN3dPtr m_vertices[8];
		/// (pointers to) subleaves
		OctLeaf* m_leaves[8];
		//-- extra data
		OctLeafExtra* m_extra;
		// local surfaces (for surface remeshing)
		SurfaceSetPtr m_local_surface_set;
		//-- other
		friend class ControlSpace3dOctree;
	};
protected:
	/// Translates the point-coordinates into index of control matrix
	void countLocalCoordinates(const DPoint3d& pt, int in[]) const;
	/// Translates the point-coordinates into (consolidated) index of control matrix
	int countLocalCoordinates(const DPoint3d& pt) const;
	/// Returns smallest Octree leaf for this point
	OctLeaf* findLastLeaf(const DPoint3d& pt) const;
	/// Extrapolate control data from nodes
	int extrapolateMainControlNodes();
	/// Propagate inside mark for the main grid
	void propagateMainGridInsideMark();
protected:
	/// Number of columns, rows and layers in the matrix
	int m_n[3];
	/// Values stored in matrix
	OctLeaf* m_grid;
	/// Grid vertices with defined control data
	OctVertexList m_grid_vertices;
};

#endif // !defined(CONTROLSPACE3DOCTREE_H__INCLUDED)
