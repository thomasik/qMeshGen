/////////////////////////////////////////////////////////////////////////////
// SurfaceQuadricQTree.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2015-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#if !defined(SURFACEQUADRICQTREE_H__INCLUDED)
#define SURFACEQUADRICQTREE_H__INCLUDED

#include "SurfaceParametric.h"
#include "DPlane.h"
#include "DTriangle.h"
#include "DataVector.h"
#include "DataMatrix.h"
#include "DRect.h"
#include "DPlanarQuadric.h"

class MeshContainer3dSurface;
class Metric2dContext;

class QuadricNode
{
	friend class SurfaceQuadricQTree;
	friend class QuadricLeaf;
//	typedef enum { NODE_INNER, NODE_BORDER, NODE_OUTER, NODE_UNKNOWN } NodeType;
public:
	QuadricNode(double x = 0.0, double y = 0.0) 
		: vq(0.0), approx_quality(AQ_UNKNOWN), coord(x,y), 
//			ntype(NODE_UNKNOWN), 
			fnearest(nullptr), min_dist2(LARGE_NUMBER)	{ }
	QuadricNode(const DPoint2d& pt) 
		: vq(0.0), approx_quality(AQ_UNKNOWN), coord(pt),
//			ntype(NODE_UNKNOWN), 
			fnearest(nullptr), min_dist2(LARGE_NUMBER)	{ }
public:
	const DPoint2d& getCooord() const { return coord; }
	bool valid(double min_q = AQ_VALID_MIN) const { return approx_quality >= min_q; }
	void setQuality(double q) { approx_quality = q; }
	bool checkFaceDist(const DMPoint2d &qpt, 
		const DMPoint2d &fa, const DMPoint2d &fb, const DMPoint2d &fc, MeshFace* face);
private:
	DVectorN<6> vq;
	double approx_quality;
	DPoint2d coord;
private:
//	NodeType	ntype;
	MeshFace*	fnearest;
	double		min_dist2;
};

typedef DataVector<QuadricNode*> QuadricVerticesPtr;
typedef DataVector<QuadricNode> QuadricVertices;

class QuadricLeaf{
	friend class SurfaceQuadricQTree;
public:
	QuadricLeaf();
	QuadricLeaf(double mid_x, double mid_y, double dx, double dy, int level = 0);
	~QuadricLeaf();
public:
	/// Get data at the given point
	DPlanarQuadric getDataAtPoint(const DPoint2d& pt) const;
	/// Set values of vertices indices
	void setVertices(QuadricNode* v_sw, QuadricNode* v_se, QuadricNode* v_nw, QuadricNode* v_ne){
		m_vertices[VX2D_SW] = v_sw;
		m_vertices[VX2D_SE] = v_se;
		m_vertices[VX2D_NW] = v_nw;
		m_vertices[VX2D_NE] = v_ne;
	}
	/// Set pointers to neighbouring leaves
	void setNeighbours(QuadricLeaf* nd, QuadricLeaf* nr, QuadricLeaf* nu, QuadricLeaf* nl){
		m_neighbours[NB_DOWN] = nd;
		m_neighbours[NB_RIGHT] = nr;
		m_neighbours[NB_UP] = nu;
		m_neighbours[NB_LEFT] = nl;
	}
	/// Set pointer to neighbouring leaf for this leaf and sub-leaves
	void setNeighbour(int side, QuadricLeaf* nb, bool skip_first = false);
	/// Set coordinates (middle point and size)
	void setCoordinates(const DPoint2d& middle, const DVector2d& dl){
		m_middle = middle;
		m_dl = dl;
	}
	// returns width of this cell
	double getWidth() const { return 4 * m_dl.x;  } 
	// returns height of this cell
	double getHeight() const { return 4 * m_dl.y; } 
	const DPoint3d getPoint( const DPoint2d& param, SurfaceConstPtr surface ) const;
	const DVector3d getDerivative(int deriv, const DPoint2d& param, SurfaceConstPtr surface) const;
	bool withinParamRange( const DPoint2d& param ) const;
	bool adaptForFaces( QuadricVertices & qv, SurfaceConstPtr surface,
		Metric2dContext& mc2d, const DataVector<MeshFace*> & sfaces );
	/// Adapt quad leaf to the givemnfunction (split if required)
	//void adaptToFunction( std::function <QuadDataItem(const DPoint2d& param)> f);
	/// Adapt quad leaf to new source point and insert this point into proper collection
	//void adaptAndInsertControlPoint(DataQuadVertices& grid_vertices, 
	//	const DataQuadNode& qv, SurfaceConstPtr surface);
	/// Adapt quad leaf to new source point and set minimum value
	//bool adaptAndSetControlPoint(DataQuadVertices& grid_vertices, 
	//	const DataQuadNode& qv, SurfaceConstPtr surface, bool min_value_set = true);
	/// Split this leaf (+ build relations with neighbours)
	bool split(QuadricVertices& grid_vertices, SurfaceConstPtr surface, 
		QuadricVerticesPtr* split_vertices = nullptr);
	/// Whether leaf is already split
	bool isSplit() const { return m_leaves[0] != nullptr; }
	/// Adjust levels of adjacent leaves
	void balance(QuadricVertices& grid_vertices, SurfaceConstPtr surface, 
		QuadricVerticesPtr* split_vertices);
	/// Store leaf structure into EPS file
	//void storeEPS(EPSFile* eps_file) const;
	/// Insert view data about leaf structure
	void addToViewSet(MeshViewSet* set, SurfaceConstPtr surface, 
		TagExtended::TagType tag = TagExtended::TAG_NONE) const;
	/// Find last leaf (not split) for givent point
	const QuadricLeaf* findLastLeaf(const DPoint2d& pt) const;
	/// Gather data from source nodes within leaf	
	//void gatherDataFromSourcePoints(SurfaceConstPtr surface);
	/// Propagate data from within leaves
	//int propagateUp();
	/// Propagate data into leaves
	//int propagateDown();
	/// Returns i-th vertex of the leaf
	QuadricNode* getVertex(int i) const { return m_vertices[i]; }
	/// Smoothen variance of metric for leaf nodes
	//int smoothenMetricAtNodes(bool uptree = true);
	/// count statistics
	//void statQuadTree(int & max_level, int & leaf_counter, int level = 1) const;
	/// Return interpolated value of some function of data at some given point;
	//double interpolateDoubleFunc(const DPoint2d& pt, std::function <double(const QuadDataItem& data)> f) const;
	/// Propagate control values from within / into leaves (returns true if whole leaf is ready)
	//bool propagateInsideMark();
	/// Propagate control values from within / into leaves (returns true if whole leaf is ready)
	//bool finishInsideMark();
	/// mark inside infomation for nodes
	//int markInsideNodes(const DPoint2d& pt0, const DPoint2d& pt1, MeshEdge2d* edge);
	/// mark inside infomation for nodes
	//void markInsideNodesY(const DPoint2d& pt0, const DPoint2d& pt1, MeshEdge2d* edge, bool right_line);
	/// mark inside infomation for nodes
	//void markInsideNodesX(const DPoint2d& pt0, const DPoint2d& pt1, MeshEdge2d* edge, bool top_line);
	/// Calculate mesh prediction using CS
	//void calculateControlPrediction(ControlPrediction & cp) const;
	/// Calculate mesh prediction using CS
	//void calculateControlPrediction(ControlPrediction & cp, int ind, int inside_count, int inside[]) const;
private:
	/// whether to use midnodes for value interpolation
	//bool m_use_midnodes;
	/// balance level
	//int m_balance_level;
private:
	enum QuadVertexWhich {VX2D_SW = 0, VX2D_SE = 1, VX2D_NW = 2, VX2D_NE = 3,
			VX2D_FIRST = 0, VX2D_LAST = 3};
	enum QuadLeafWhich {LF_SW = 0, LF_SE = 1, LF_NW = 2, LF_NE = 3,
			LF_FIRST = 0, LF_LAST = 3};
	enum QuadNeighbourWhich {NB_DOWN = 0, NB_RIGHT = 1, NB_UP = 2, NB_LEFT = 3, 
		NB_FIRST = 0, NB_LAST = 3};
	enum QuadMidVertexWhich {MVX2D_DOWN = 0, MVX2D_RIGHT = 1, MVX2D_UP = 2, MVX2D_LEFT = 3,
		MVX2D_FIRST = 0, MVX2D_LAST = 3};
	static QuadVertexWhich mid_to_nodes[4][2];
	static QuadMidVertexWhich nodes_to_mid[4][2];
	static QuadricLeaf empty;
private:
	/// middle point of leaf (parametric)
	DPoint2d m_middle;
	QuadricNode* m_vertices[4]; // only indices, 2-3 / 0-1
	QuadricLeaf* m_leaves[4];
	//-- extra data
	/// depth level
	int m_level;
	/// valid sub-domain?
	bool m_valid;
	/// dimensions of leaf
	DVector2d m_dl; // 1/4 of width and height of leaf 
	QuadricLeaf* m_neighbours[4]; 
	QuadricNode* m_mid_vertices[4]; 
	//--- mesh faces per leaf ?
	//DataVector<MeshFace*> *m_faces;
};

typedef DataMatrix<QuadricLeaf> QuadricLeaves;

class SurfaceQuadricQTree :	public SurfaceParametric
{
public:
	SurfaceQuadricQTree( SurfacePtr base_surface,
		const DRect& brect,	int nx = 4, int ny = 4);
public:
	static SurfacePtr fitQuadricQTreeSurface(Metric3dContext &mc,
		SurfacePtr base_surface, const DRect& brect, double tolerance,
		MeshContainer3dSurface* mesh, const DataVector<MeshFace*> &sfaces, 
		int nxy = 100);
public:
	virtual string getSimpleDescription() const override { return "quadric-qtree"; }
	/// Returns the object-specific type
	virtual ElementType getType() const { return SURFACE_QUADRIC_QTREE; }
	virtual const DPoint3d getPoint( const DPoint2d& param ) const override;
	virtual const DVector3d getDerivative(int deriv, const DPoint2d& param) const override;
	void countLocalCoordinates(const DPoint2d &pt, int &ix, int &iy) const;
	const QuadricLeaf* findLastLeaf(const DPoint2d& pt) const;
	bool initializeWithFaces( Metric3dContext &mc, MeshContainer3dSurface* mesh, 
		const DataVector<MeshFace*> & sfaces, double tolerance );
	void extrapolateGridLayer();
	/// Create wire-frame visualization for rectangular surface patch covering the given set of points
	virtual MeshViewSet * createViewSet( MeshViewSet* set ) const;
	MeshViewSet* createViewSetForBaseGrid( MeshViewSet* set, int id = -2 ) const;
	bool showNode( MeshViewSet* set, int k, 
		DataVector< int > & shown_nodes, DataVector< DPoint3d > & pt3d_nodes ) const;
	virtual bool withinParamRange( const DPoint2d& param ) const override;
	virtual bool invertOrientation() override;
private:
	SurfacePtr m_base_surface;
	// rotation matrix (+inverse?)
	// qtree
protected:
	/// Rectangle specifying the area of quad tree
	DRect m_box;
	/// Number of columns in the matrix
	int m_nx;
	/// Number of rows in the matrix
	int m_ny;
	/// Values stored in matrix
	QuadricLeaves m_grid;
	/// Grid vertices with defined control data
	QuadricVertices m_grid_vertices;
};

#endif // SURFACEQUADRICQTREE_H__INCLUDED
