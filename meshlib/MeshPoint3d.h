/////////////////////////////////////////////////////////////////////////////
// MeshPoint3D.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2002-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#if !defined(MESHPOINT3D_H__INCLUDED)
#define MESHPOINT3D_H__INCLUDED

#include "DPoint.h"
#include "MeshData.h"
#include "DataVector.h"
#include "DataHashTable.h"
#include "DataList.h"
#include "TagBorder.h"
#include "TagExtended.h"
#include "SurfaceParametric.h"
#include "Curve3dParametric.h"

class MeshEdge3d;
class MeshFace;
class MeshBlock;
struct MeshViewPointData;
class Metric3dContext;

class DataStatistics;

struct ParamAndQuality {
	ParamAndQuality( const DPoint2d& _param = DPoint2d::zero, double aq = AQ_UNKNOWN ) : param(_param), quality(aq) {}
	//ParamAndQuality( const DPoint2d& param, double aq ) : param(param), quality(aq) {}
	ParamAndQuality( SurfaceConstPtr surface, const DPoint3d& pt );
	bool isValid() const { return quality >= AQ_VALID_MIN; }
	DPoint2d param;
	double quality;
};

//#define LOG_GETLOCALSURFACE(str) LOG4CPLUS_INFO(MeshLog::logger_mesh, "GETLOCALSURFACE=" << str)
#define LOG_GETLOCALSURFACE(str)

/**
 * This class implements a mesh point in 3D with its 
 *	coordinates, incident mesh edges and additional data.
 */
class MeshPoint3d : public IndexedObject, public TagBorder, public TagExtended
{
public:
	/// Standard constructor
	MeshPoint3d(double x, double y, double z = 0.0);
	/// Standard constructor
	MeshPoint3d(const DPoint3d& pt);
	/// Surface-point constructor
	MeshPoint3d(SurfaceConstPtr surface, const DPoint2d& param, double aq = AQ_UNKNOWN, 
		const DPoint3d * point3d = nullptr);
	/// Curve-point constturctor
	MeshPoint3d(Curve3dConstPtr curve, const double& t, 
		const DPoint3d * point3d = nullptr);
	/// Copying constructor
	MeshPoint3d(const MeshPoint3d& point);
	/// Standard destructor
	~MeshPoint3d();
public:
	/// clear some data for faster all-delete process
	void preDeleteAll();
	/// Copy geometric data from other point
	void copyDataFrom(const MeshPoint3d* point);
	/// Enumerates all adjacent blocks (with hash table)
	bool adjacentBlocksWithHash(DataVector<MeshBlock*> & blocks) const;
	/// Enumerates all adjacent blocks
	bool adjacentBlocks(DataVector<MeshBlock*> & blocks, bool boundary_check = true) const;
	/// Enumerates all adjacent blocks
	bool adjacentBlocks(DataVector<MeshBlock*> & blocks, DataHashTable<MeshBlock*> & visited_blocks,
		bool boundary_check = true) const;
	/// Enumerates all adjacent faces
	bool adjacentFaces(DataVector<MeshFace*> & faces) const;
	/// Returns the adjacent mesh face, joining this point with the given two (nullptr if doesn't exist)
	MeshFace* getFaceToPoints(const MeshPoint3d* point1, const MeshPoint3d* point2) const;
	/// Returns the description of this point for visualization
	std::shared_ptr<MeshViewPointData> getViewData() const;
	/// Adds the reference to the give mesh edge into the adjacency array
	void addEdgeLink(MeshEdge3d* edge);
	/// Removes the reference to the given mesh edge from the adjacency array
	void removeEdgeLink(MeshEdge3d* edge) { edges.remove(edge); }
	/// Returns the rank of the point (count of adjacent edges)
	int getRank() const { return (int)edges.countInt(); }
	/// Returns the number of adjacent boundary edges
	int getBorderEdgesCount() const;
	/// Returns the adjacent mesh edge at the given index
	MeshEdge3d* getEdge(int i) const { return edges[i]; }
	/// Returns the adjacent mesh edge, joining this point with the given one (nullptr if doesn't exist)
	MeshEdge3d* getEdgeToPoint(const MeshPoint3d* point) const;
	/// Returns the coordinates of this point
	const DPoint3d& getCoordinates() const { return coord; }
	/// Returns the metric coordinates of this point
	const DMPoint3d getMetricCoordinates(const Metric3dContext& mc);
	/// Sets the coordinates for this point
	void setCoordinates(const DPoint3d& pt){ coord = pt; assert( local_surface == nullptr ); }
	/// Sets the coordinates for this point
	void setCoordinates( SurfaceConstPtr surface, const DPoint2d& param, 
		double aq = AQ_UNKNOWN, const DPoint3d * point3d = nullptr);
	/// Sets the coordinates for this point
	void setCoordinates( Curve3dConstPtr curve, const double& t, 
		const DPoint3d * point3d = nullptr);
	/// Set active mode
	void setActive(int mode = 1) { active = mode; }
	/// Clear active mode
	void clearActive() { active = 0; }
	/// Check active mode
	bool isActive() const { return active != 0; }
	/// Get active mode
	int getActive() const { return active; }
	/// Check if point is planar with respect to adjacent faces (according to normals)
	bool planar(const Metric3dContext& mc) const;
	/// Get arbitrary face adjacent to this point
	MeshFace* anyAdjacentFace() const;
	/// Get local vicinity set of mesh points (with this point) - border and non-border
	DataVector<MeshPoint3d*> localVicinityPoints();
	/// Get local vicinity set of mesh points (with this point) - only border or only non-border
	DataVector<MeshPoint3d*> localVicinityPointsBorderEdges( bool border = true );
	/// Get neighboring point
	MeshPoint3d* incidentPoint( int i ) const;
public:
	/// Checks validity of local surface params
	bool checkLocalSurfaceParams() const;
	/// optimize selection of local surface 
	//bool selectBestLocalSurface( Metric3dContext& mc, int local_surface_tag = -1 );
	/// Tries to move point to new coordinates (gradually, checking for inverted faces)
	bool tryMovingPoint(Metric3dContext& mc, const DPoint3d& new_pt);
	/// Tries to move point to new coordinates (gradually, checking for inverted faces)
	bool tryMovingPoint(/* Metric3dContext& mc, */ SurfaceConstPtr surface, 
		const DPoint2d& new_pt, const DataVector< MeshFace* > & faces);
	/// Tries to move point to new coordinates (gradually, checking for inverted faces)
	bool tryMovingPoint(/* Metric3dContext& mc, */ Curve3dConstPtr curve, const double& t);
	/// Tries to move point to fit its ascribed local surface/curve
	bool moveToLocalShape(Metric3dContext& mc);
	/// Sets coordinates of the point after fitting to its ascribed local surface/curve
	//bool setToLocalShape(Metric3dContext & mc);
	/// Returns new coordinates of the point after fitting to its ascribed local surface/curve
	//DPoint3d getCoordinatesWithLocalShape(Metric3dContext & mc, const DPoint3d& new_pt);
	/// Clear all local shapes
	void clearLocalShapes();
	/// Returns the local curve
	Curve3dConstPtr getLocalCurve() const;
	/// Changes local curve
	void setLocalCurve(Curve3dConstPtr curve, const double & t);
	/// Checks local curve
	bool hasLocalCurve() const { return local_curve != nullptr; }
	/// Returns local param for surface
	double getLocalCurveParam( Curve3dConstPtr curve) const;
	/// Returns local param for surface
	double getLocalCurveParam( ) const;
	/// Returns the local surface
	SurfaceConstPtr getLocalSurface() const;
	/// Returns the local surface
	SurfaceConstPtr getLocalValidSurface() const;
	/// Returns the local surface, valid and common to all points
	SurfaceConstPtr getOptLocalCommonSurface( int ls_tag, const DataVector< MeshPoint3d* > & mpoints,
		DataHashTable< SurfaceConstPtr > & hchecked_surfaces, 
		SurfaceConstPtr & opt_surface, double & opt_quality ) const;
	/// Whether this surface is valid for this point
	double getLocalSurfaceQuality( SurfaceConstPtr surface ) const;
	/// Changes local surface
	void setLocalSurface(SurfaceConstPtr surface, const DPoint2d& param, double aq );
	/// Checks local surface
	bool hasLocalSurface() const { return local_surface != nullptr; }
	/// Checks the specific local surface
	bool hasLocalSurface(SurfaceConstPtr surface ) const;
	/// Returns local param for surface
	ParamAndQuality getLocalSurfaceParamQuality( SurfaceConstPtr surface) const;
	DPoint2d getLocalSurfaceParam( SurfaceConstPtr surface) const;
	/// Returns local param for surface
	//ParamAndQuality getLocalSurfaceParam( ) const;
	DPoint2d getLocalSurfaceParam( ) const;
	/// Get parameters for this surface
	//ParamAndQuality checkAndGetSurfaceParam( SurfaceConstPtr surface );
	bool checkAndGetSurfaceParam( SurfaceConstPtr surface, DPoint2d& param );
	/// Get parameters for this curve
	bool checkAndGetCurveParam( Curve3dConstPtr curve, double& t );
	/// Get (base) normal vector
	const DVector3d& getBaseNormal() const { return base_normal; }
	/// Set (base) normal vector
	void setBaseNormal(const DVector3d& vn) { base_normal = vn; }
	/// Returns the base (initial) coordinates of this point
	//const DPoint3d& getBaseCoordinates() const { return base_coord; }
	/// Move point to new coordinates, with constraint of the maximum metric distance from base coordinates
	//bool moveLocal(Metric3dContext& mc, const DPoint3d& new_coord, double max_dist = 1.0, bool max_cut = true);
	/// Gather local surface candidates for this point
	//bool gatherLocalSurfaceCandidates( Metric3dContext & mc, 
	//	int domain_tag, DataVector< SurfaceConstPtr > & scandidates, DPoint3d * dpt = nullptr ) const;
	/// Increment value for each (valid, no domain-check) 
	//void incForValidSurfaces( DataHashTableKeyValue< SurfaceConstPtr, int > & hsurface_counter, int local_surface_tag ) const;
	///// Get all local surfaces for this points
	//void getAllLocalSurfaces(DataVector<SurfaceConstPtr> & surfaces, int local_surface_tag ) const;
	/// Increment value for each (valid, no domain-check) 
	//void incForValidCurves( DataHashTableKeyValue< Curve3dConstPtr, int > & hcurve_counter ) const;
	/// Get all local curves for this points
	//void getAllLocalCurves(DataVector<Curve3dConstPtr> & curves ) const;
	/// check max dist2 for local surface params and count it as well
	//double checkMaxLocalSurfaceParamDist2( DataStatistics& stats ) const;
private:
	/// Recalculates coordinates according to local shape(s) - with possible change of active shape!
	//bool recalculateForLocalShape( Metric3dContext & mc );
	/// Recalculates coordinates according to local shape(s) - with possible change of active shape!
	//bool recalculateForLocalShape( Metric3dContext & mc, DPoint3d & pt );
	class SurfaceData
	{
		friend MeshPoint3d;
		struct ParamData {
			ParamData(SurfaceConstPtr _surface, const DPoint3d& _ref_pt = DPoint3d::zero,
				const DPoint2d& _param = DPoint2d::zero, double _aq = AQ_UNKNOWN ); 
			void set( const DPoint2d& _param, double aq, const DPoint3d& _pt ) { 
				param = _param; ref_pt = _pt; quality = aq; }
			double isValidGetQuality(const DPoint3d& pt) const { 
				return (ref_pt == pt) ? quality : (assert(false),AQ_UNKNOWN); 
			}
			bool isValid(const DPoint3d& pt) const { return (ref_pt == pt) ? isValid() : false; }
			bool isValid() const { return quality >= AQ_VALID_MIN; }
			double updateGetQuality( const DPoint3d& pt );
			bool update( const DPoint3d& pt ) { return updateGetQuality( pt ) >= AQ_VALID_MIN; 	}
			double updateIfNeededGetQuality( const DPoint3d& pt );
			bool updateIfNeeded( const DPoint3d& pt ) {
				return updateIfNeededGetQuality( pt ) >= AQ_VALID_MIN; }
			ParamAndQuality getParamAndQuality() const { return ParamAndQuality( param, quality ); }
		public:
			SurfaceConstPtr surface;
			DPoint3d ref_pt;
			DPoint2d param;
			double quality; // 1 - perfect, 0 - barely valid, -0.5 - invalid, -1.0 - unknown
		};
	public:
		/// Standard constructor
		SurfaceData( SurfaceConstPtr local_surface, const DPoint2d& param, double aq, const DPoint3d& pt );
		SurfaceData( ) { }
	public:
		/// Search for stored params for this suface (returns nullptr if missing)
		ParamData* findSurfaceParamData( SurfaceConstPtr surface ) const;
		/// Checks validity of local surface params
		bool checkSurfaceParams() const;
		/// Returns the local surface, valid and common to all points
		SurfaceConstPtr getOptLocalCommonSurface( int ls_tag, const MeshPoint3d* mpoint,
			const DataVector< MeshPoint3d* > & mpoints,
			DataHashTable< SurfaceConstPtr > & hchecked_surfaces,
			SurfaceConstPtr & opt_surface, double & opt_quality);
		/// Returns (any) local surface
		SurfaceConstPtr getLocalSurface( );
		/// Returns valid surface
		SurfaceConstPtr getValidSurface(const DPoint3d& pt);
		/// Whether this surface is valid for this point
		double getLocalSurfaceQuality(SurfaceConstPtr surface, const DPoint3d& pt);
		/// Checks the specific local surface (not checking if valid parameters, ie. if still in domain)
		bool hasLocalSurface(SurfaceConstPtr surface ) const;
		/// Add/change local surface
		void setSurfaceParam(SurfaceConstPtr surface, const DPoint2d& param, double aq, const DPoint3d& pt );
		/// Get parameters for this surface
		DPoint2d getSurfaceParam( SurfaceConstPtr surface, const DPoint3d& pt );
		ParamAndQuality getSurfaceParamQuality( SurfaceConstPtr surface, const DPoint3d& pt );
		/// Get parameters for this surface
		bool checkAndGetSurfaceParam( SurfaceConstPtr surface, DPoint2d& param, const DPoint3d& pt );
		ParamAndQuality checkAndGetSurfaceParamQuality( SurfaceConstPtr surface, const DPoint3d& pt );
		/// Recheck best surface
		//bool recheckBestSurface( const DPoint3d& pt );
		/// Increment value for each (valid, no domain-check) 
		//void incForValidSurfaces( DataHashTableKeyValue< SurfaceConstPtr, int > & hsurface_counter, 
		//	int local_surface_tag, const DPoint3d& pt ) const;
		///// Get all local surfaces for this points
		//void getAllLocalSurfaces(DataVector<SurfaceConstPtr> & surfaces, int local_surface_tag, 
		//	const DPoint3d& pt  );
		/// check max dist2 for local surface params and count it as well
		//double checkMaxLocalSurfaceParamDist2( DataStatistics& stats, const DPoint3d& pt ) const;
	protected: 
		/// Parameters for local surfaces (for multi-partial-surface)
		DataSimpleList<ParamData> surface_params;
		// auxiliary hash table 
		//DataHashTableKeyValue<SurfaceConstPtr, ParamData*> hsparams;
	};
	class CurveData
	{
		friend MeshPoint3d;
		struct ParamData {
			ParamData(Curve3dConstPtr _curve, const DPoint3d& _ref_pt, double _param );
			ParamData(Curve3dConstPtr _curve, const DPoint3d& _ref_pt);
			void set( const double& _param, const DPoint3d& pt);
			bool isValid(const DPoint3d& pt) const { return (ref_pt == pt) && valid; }
			bool update( const DPoint3d& pt );
			bool updateIfNeeded( const DPoint3d& pt );
		public:
			Curve3dConstPtr curve;
			DPoint3d ref_pt;
			double param;
			bool valid;
		};
	public:
		/// Standard constructor
		CurveData( Curve3dConstPtr local_curve, const double& t, const DPoint3d& pt);
	public:
		/// Returns (any) local curve
		Curve3dConstPtr getLocalCurve( );
		/// Returns valid curve
		Curve3dConstPtr getValidCurve(const DPoint3d& pt);
		/// Search for stored params for this curve (returns nullptr if missing)
		ParamData* findCurveParamData(Curve3dConstPtr curve ) const;
		/// Get parameters for this curve
		bool checkAndGetCurveParam(Curve3dConstPtr curve, const DPoint3d& pt, double& t );
		/// Set/change local curve param
		void setCurveParam(Curve3dConstPtr curve, const double& t, const DPoint3d& pt );
		/// Get parameters for this curve
		double getCurveParam(Curve3dConstPtr curve, const DPoint3d& pt );
		/// Increment value for each (valid, no domain-check) 
		//void incForValidCurves( DataHashTableKeyValue< Curve3dConstPtr, int > & hcurve_counter, 
		//	const DPoint3d& pt ) const;
		/// Get all local curves for this points
		//void getAllLocalCurves(DataVector<Curve3dConstPtr> & curves, const DPoint3d& pt  );
	protected: 
		/// Parameters for local curves (for corner-points)
		//DataHashTableKeyValue<Curve3dConstPtr, ParamData> * curve_params;
		DataSimpleList<ParamData> curve_params;
	};
protected: // general
	/// 3D coordinates
	DPoint3d	coord;
	/// Adjacent mesh edges
	DataVector<MeshEdge3d*> edges;
	/// Auxiliary flag for concurrency 
	int active;
protected: // for 3dsurface points
	/// Local surface
	SurfaceData* local_surface;
	/// Local curve
	CurveData* local_curve;
protected:
	/// Base vector (average)
	DVector3d base_normal;
	/// Base coordinates (initial - for translation control)
	//DPoint3d base_coord;
};

#endif // !defined(MESHPOINT3D_H__INCLUDED)
