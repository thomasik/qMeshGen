/////////////////////////////////////////////////////////////////////////////
// DataQuadTree.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2015-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#if !defined(DATAQUADTREE_H__INCLUDED)
#define DATAQUADTREE_H__INCLUDED

#include <functional>

#include "SurfaceParametric.h"
#include "DRect.h"
#include "MeshViewSet.h"

#include "Curve2dSegment.h"

//#include "DataMatrix.h"
#include "DataVector.h"
//#include "DataHashTable.h"

class EPSFile;

#define QuadDataItem int // -> to be typename

class DataQuadNode
{
	friend class DataQuadLeaf;
	friend class DataQuadTree;
public:
	DataQuadNode(double x = 0.0, double y = 0.0) : coord(x,y) { }
private:
	QuadDataItem data;
	DPoint2d coord;
};

typedef DataVector<DataQuadNode*> DataQuadVerticesPtr;
typedef DataVector<DataQuadNode> DataQuadVertices;

/// Class QuadLeaf for QuadTree structure
class DataQuadLeaf{
	friend class DataQuadTree;
public:
	DataQuadLeaf();
	DataQuadLeaf(double mid_x, double mid_y, double dx, double dy, int level = 0);
	~DataQuadLeaf();
public:
	/// Get data at the given point
	QuadDataItem getDataAtPoint(const DPoint2d& pt) const;
	/// Set values of vertices indices
	void setVertices(DataQuadNode* v_sw, DataQuadNode* v_se, DataQuadNode* v_nw, DataQuadNode* v_ne){
		m_vertices[VX_SW] = v_sw;
		m_vertices[VX_SE] = v_se;
		m_vertices[VX_NW] = v_nw;
		m_vertices[VX_NE] = v_ne;
	}
	/// Set pointers to neighbouring leaves
	void setNeighbours(DataQuadLeaf* nd, DataQuadLeaf* nr, DataQuadLeaf* nu, DataQuadLeaf* nl){
		m_neighbours[NB_DOWN] = nd;
		m_neighbours[NB_RIGHT] = nr;
		m_neighbours[NB_UP] = nu;
		m_neighbours[NB_LEFT] = nl;
	}
	/// Set pointer to neighbouring leaf for this leaf and sub-leaves
	void setNeighbour(int side, DataQuadLeaf* nb, bool skip_first = false);
	/// Set coordinates (middle point and size)
	void setCoordinates(const DPoint2d& middle, const DVector2d& dl){
		m_middle = middle;
		m_dl = dl;
	}
	/// Adapt quad leaf to the givemnfunction (split if required)
	void adaptToFunction( std::function <QuadDataItem(const DPoint2d& param)> f);
	/// Adapt quad leaf to new source point and insert this point into proper collection
	void adaptAndInsertControlPoint(DataQuadVertices& grid_vertices, 
		const DataQuadNode& qv, SurfaceConstPtr surface);
	/// Adapt quad leaf to new source point and set minimum value
	bool adaptAndSetControlPoint(DataQuadVertices& grid_vertices, 
		const DataQuadNode& qv, SurfaceConstPtr surface, bool min_value_set = true);
	/// Split this leaf (+ build relations with neighbours)
	bool split(DataQuadVertices& grid_vertices, SurfaceConstPtr surface, 
		DataQuadVerticesPtr* split_vertices = nullptr, bool init_values = false);
	/// Whether leaf is already split
	bool isSplit() const { return m_leaves[0] != nullptr; }
	/// Adjust levels of adjacent leaves
	void balance(DataQuadVertices& grid_vertices, SurfaceConstPtr surface, 
		DataQuadVerticesPtr* split_vertices, bool init_values = false);
	/// Store leaf structure into EPS file
	void storeEPS(EPSFile* eps_file) const;
	/// Insert view data about leaf structure
	void addToViewSet(MeshViewSet* set, SurfaceConstPtr surface, 
		TagExtended::TagType tag = TagExtended::TAG_NONE) const;
	/// Find last leaf (not split) for givent point
	DataQuadLeaf* findLastLeaf(const DPoint2d& pt);
	/// Gather data from source nodes within leaf	
	void gatherDataFromSourcePoints(SurfaceConstPtr surface);
	/// Propagate data from within leaves
	int propagateUp();
	/// Propagate data into leaves
	int propagateDown();
	/// Returns i-th vertex of the leaf
	DataQuadNode* getVertex(int i) const { return m_vertices[i]; }
	/// Smoothen variance of metric for leaf nodes
	int smoothenMetricAtNodes(bool uptree = true);
	/// count statistics
	void statQuadTree(int & max_level, int & leaf_counter, int level = 1) const;
	/// Return interpolated value of some function of data at some given point;
	double interpolateDoubleFunc(const DPoint2d& pt, std::function <double(const QuadDataItem& data)> f) const;
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
	/// middle point of leaf (parametric)
	DPoint2d m_middle;
	enum QuadVertexWhich {VX_SW = 0, VX_SE = 1, VX_NW = 2, VX_NE = 3,
			VX_FIRST = 0, VX_LAST = 3};
	DataQuadNode* m_vertices[4]; // only indices, 2-3 / 0-1
	enum QuadLeafWhich {LF_SW = 0, LF_SE = 1, LF_NW = 2, LF_NE = 3,
			LF_FIRST = 0, LF_LAST = 3};
	DataQuadLeaf* m_leaves[4];
	/// Parameterization matrix for leaf middle point
	//ControlDataMatrix2d m_middle_cdm;
	//-- extra data
	/// depth level
	int m_level;
	/// dimensions of leaf (parametric)
	DVector2d m_dl;
	/// dimensions of leaf (real)
	DVector2d m_real_dl;
	enum QuadNeighbourWhich {NB_DOWN = 0, NB_RIGHT = 1, NB_UP = 2, NB_LEFT = 3, 
		NB_FIRST = 0, NB_LAST = 3};
	DataQuadLeaf* m_neighbours[4]; 
	enum QuadMidVertexWhich {MVX_DOWN = 0, MVX_RIGHT = 1, MVX_UP = 2, MVX_LEFT = 3,
		MVX_FIRST = 0, MVX_LAST = 3};
	static QuadVertexWhich mid_to_nodes[4][2];
	static QuadMidVertexWhich nodes_to_mid[4][2];
	DataQuadNode* m_mid_vertices[4]; 
	//--- point sources before interpolation
	DataVector<DataQuadNode> *m_source_points;
};

DataQuadLeaf::QuadVertexWhich		DataQuadLeaf::mid_to_nodes[4][2] = 
			{{VX_SW, VX_SE}, {VX_SE, VX_NE},
			{VX_NW, VX_NE}, {VX_SW, VX_NW}};
DataQuadLeaf::QuadMidVertexWhich	DataQuadLeaf::nodes_to_mid[4][2] = 
			{{MVX_DOWN, MVX_LEFT}, {MVX_DOWN, MVX_RIGHT},
			{MVX_UP, MVX_LEFT}, {MVX_UP, MVX_RIGHT}};

/**
 * This class implements a general quadtree structure.
 */
class DataQuadTree
{
public:
	/// Standard contructor
	DataQuadTree(SurfaceConstPtr param_surface, const DRect& bbox, int nxy = 100);
	/// Destructor
	~DataQuadTree() { if(m_grid) delete[] m_grid; }
public:
	/// Returns the "screenshot" of this control space for visualization
	MeshViewSet* getViewSet(MeshViewSet* set = nullptr, bool with_surface = false, 
		TagExtended::TagType tag = TagExtended::TAG_NONE) const;
	/// Adds new data for some point within the domain
	void addDataQuadNode(const DataQuadNode& node);
	/// Get sizing info (matrix mode) at the given point
	QuadDataItem getDataAtPoint(const DPoint2d& pt) const;
	/// Interpolates (initializes) the control space basing on the previously given sizing information in discrete points
	bool interpolate();
	/// Draws the structure of the control space
	void storeEPS(const char* name = "data-quadtree", int id = 0);
	/// Initializes the control space basing on the curvature of surface and boundary
	//virtual void setSurfaceCurvatureControlData();
	/// Refines control space at the given point
	bool setMinData(const DPoint2d& pt, const QuadDataItem& data, bool min_value_set = true);
	/// Returns number of control nodes in adaptive control structure
	int getDataQuadNodesCount() const {	m_grid_vertices.countInt(); }
	/// Returns i-th control node
	DataQuadNode& getDataQuadNode(int i) {
		assert(i >= 0 && i < m_grid_vertices.countInt());
		return m_grid_vertices.get(i);	
	}
	/// Returns i-th control node
	const DataQuadNode& getDataQuadNode(int i) const {
		assert(i >= 0 && i < m_grid_vertices.countInt());
		return m_grid_vertices.get(i);	
	}
	/// Refines structure according to function variance
	void adaptToFunction( std::function <QuadDataItem(const DPoint2d& param)> f );
	/// Log basic information about this control space
	void logDescription() const;
	/// Smoothen variance of metric within the control space
	bool smoothen();
	/// Return interpolated value of extended tag data from control nodes at some given point;
	double interpolateDoubleFunc(const DPoint2d& pt, std::function <double(const QuadDataItem& data)> f) const;
	/// Set given tag value as minimum at the closest control element
	void applyFuncForLeaf(const DPoint2d& pt, std::function <void(const QuadDataItem& data)> f);
	/// Calculate inside information for nodes/elements
	//virtual void markInsideNodes(const MeshContainer2d* boundary_mesh);
	/// Returns local resolution
	double getLocalResolution(const DPoint2d& pt) const;
public:
	void testQuadTree();
private:
	/// maximum depth for quadtree
	//int param_max_depth;
	/// maximum depth for quadtree
	//int param_max_depth_adaptation;
public:
	/// Show CS-grid with inside info for debug
	//void showMarkInside(const MeshContainer2d* boundary_mesh = nullptr) const;
	/// Calculate mesh prediction using CS
	//ControlPrediction calculateControlPrediction() const;
protected:
	/// Translates the point-coordinates into index of control matrix
	void countLocalCoordinates(const DPoint2d& pt, int& ix, int& iy) const;
	/// Translates the point-coordinates into (consolidated) index of control matrix
	int countLocalCoordinates(const DPoint2d& pt) const;
	/// Returns smalles QuadTree leaf for this point
	DataQuadLeaf* findLastLeaf(const DPoint2d& pt) const;
	/// Extrapolate control data from nodes
	int extrapolateMainDataQuadNodes();
	/// Propagate inside mark for the main grid
	//void propagateMainGridInsideMark();
protected:
	/// Rectangle specifying the area of quad tree
	DRect m_box;
	/// Number of columns in the matrix
	int m_nx;
	/// Number of rows in the matrix
	int m_ny;
	/// Values stored in matrix
	DataQuadLeaf* m_grid;
	/// Grid vertices with defined control data
	DataQuadVertices m_grid_vertices;
};

DataQuadTree::DataQuadTree(SurfaceConstPtr surface, const DRect& box, int nxy) 
	: m_box(box), m_grid_vertices(2*nxy)
{
	assert(nxy >= 4);
	// length x
	Curve2dSegment middle_xline(
		DPoint2d(box.x0, (box.y0+box.y1)*0.5), 
		DPoint2d(box.x1, (box.y0+box.y1)*0.5));
	double real_dx = middle_xline.getLengthOnSurface(0.0, 1.0, surface);
	// length y
	Curve2dSegment middle_yline(
		DPoint2d((box.x0+box.x1)*0.5, box.y0), 
		DPoint2d((box.x0+box.x1)*0.5, box.y1));
	double real_dy = middle_yline.getLengthOnSurface(0.0, 1.0, surface);

	double ratio = real_dx / real_dy;	// adjust m_nx:m_ny ratio to real lengths of domain
	m_nx = std::max(2,(int)(sqrt(ratio*nxy)));
	m_ny = std::max(2, (int)(m_nx/ratio));
	m_nx = std::min(m_nx, nxy/2);
	m_ny = std::min(m_ny, nxy/2);

	double dx = box.getDX() / m_nx;
	double dy = box.getDY() / m_ny;

	m_grid = new DataQuadLeaf[m_nx*m_ny];
	
	// add regular grid vertices
	DataQuadNode qv(m_box.x0, m_box.y0);
	for(int j = 0; j <= m_ny; j++, qv.coord.y += dy){
		qv.coord.x = m_box.x0;
		for(int i = 0; i <= m_nx; i++, qv.coord.x += dx){
			m_grid_vertices.add(qv);
		}
	}

	// init regular grid boxes with coordinates and vertices indices
	DVector2d dpt(dx/4, dy/4);
	DPoint2d pt(m_box.x0 + dx / 2, m_box.y0 + dy / 2);
	for(int j = 0, k = 0, ix = 0; j < m_ny; j++, ix++){
		for(int i = 0; i < m_nx; i++, ix++){
			m_grid[k].setCoordinates(pt, dpt);
			pt.x += dx;
			m_grid[k].setNeighbours(
				((j>0)?(&m_grid[k-m_nx]):nullptr), 
				((i<m_nx-1)?(&m_grid[k+1]):nullptr),
				((j<m_ny-1)?(&m_grid[k+m_nx]):nullptr), 
				((i>0)?(&m_grid[k-1]):nullptr));
			m_grid[k++].setVertices(
				&(m_grid_vertices[ix]), 
				&(m_grid_vertices[ix+1]), 
				&(m_grid_vertices[ix+m_nx+1]), 
				&(m_grid_vertices[ix+m_nx+2]));
		}
		pt.y += dy;
		pt.x = m_box.x0 + dx / 2;
	}
}


#endif // !defined(DATAQUADTREE_H__INCLUDED)
