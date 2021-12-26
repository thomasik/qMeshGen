/////////////////////////////////////////////////////////////////////////////
// SurfaceParametric.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2002-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#if !defined(SURFACEPARAMETRIC_H__INCLUDED)
#define SURFACEPARAMETRIC_H__INCLUDED

#include <memory>

#include "DPoint.h"
#include "DMetric2d.h"
#include "DataVector.h"
#include "TagExtended.h"
#include "MeshData.h"
#include "DOrientedBox.h"
#include "SurfaceCurvature.h"
#include "Metric3dContext.h"

class Curve2dParametric;
class ControlDataMatrix3d;
class MeshViewSet;
class SurfaceDomain;
class MeshPoint3d;
class MeshFace;

/**
 * This class defines a base,abstract parametrically-defined surface in 3D
 *	with several required methods.
 */
class SurfaceParametric : public TagExtended
{
public:
	/// Standard constructor
	SurfaceParametric() : m_valid(false), m_domain(nullptr), m_fixed(false) {}
	/// Standard destructor
	virtual ~SurfaceParametric();
public:
	static std::shared_ptr<SurfaceParametric> fitSurface(Metric3dContext& mc, 
		const DataVector<DPoint3d> & points, 
		const DataVector<DVector3d> & normals, 
		DataVector<DPoint2d> & params, 
		double tolerance, MeshFace* base_face = nullptr, 
		double plane_fit_ratio = 0.1, 
		DataVector<double> * approx_quality = nullptr);
	static bool checkSurfaceFit( Metric3dContext& mc, 
		std::shared_ptr<SurfaceParametric> &surface,
		const DataVector<DPoint3d> & points, 
		const DataVector<DVector3d> & normals, 
		DataVector<DPoint2d> & params, double tolerance, 
		DataVector<double> * approx_quality = nullptr,
		bool require_within_tolerance = true);
	static std::shared_ptr<SurfaceParametric> calculateBestBaseSurface( 
		const DataVector<MeshFace*> & faces, double cmin_eps );
	static std::shared_ptr<SurfaceParametric> calculateBestBaseSurface( 
		MeshFace* face, double cmin_eps, MeshFace* central_face = nullptr );
public:
	virtual string getSimpleDescription() const { return "unknown"; }
	/// return oriented bbox for this surface (invalid)
	virtual DOrientedBox getOrientedBox() const { assert(false); return DOrientedBox(); }
	/// return oriented bbox for this surface initialized with the given set of points
	virtual DOrientedBox getOrientedBox( const DataVector<DPoint3d>& /* points */) const { assert(false); return DOrientedBox(); }
	/// return oriented bbox for this surface initialized with the given set of points
	virtual DOrientedBox getOrientedBoxOpt( const DataVector<DPoint3d>& /* points */) const { assert(false); return DOrientedBox(); }
	// fix center and range (for parameters) to fix periodic surfaces
	virtual bool fixCenterAndRangeForPeriodic(const DPoint3d& /* center_point */) { return false; }
	/// invert the orientation of the surface (change diretion of normal vector) if possible
	virtual bool invertOrientation() { assert(false); return false; }
	/// store domain to XML file
	virtual ostream& storeDomainXML(ostream& os, const string& prefix = "") const;
	/// draw/add viewset data for the domain
	void drawDomain(MeshViewSet* set = nullptr, int id = -2) const;
	/// set new domain info
	void setDomain(SurfaceDomain* domain);
	/// domain check
	//virtual bool withinDomain(Metric3dContext& mc, MeshPoint3d * point, DPoint2d * param = nullptr, double * quality = nullptr) const;
	/// domain check (-1.0 - invalid, 0-1 - quality)
	virtual double withinDomainQuality( const DPoint3d& pt, DPoint2d & param ) const;
	virtual bool withinDomain( const DPoint3d& pt, DPoint2d & param ) const {
		return withinDomainQuality( pt, param ) >= AQ_VALID_MIN; }
	/// return domain middle param
	virtual DPoint2d getDomainMiddleParam() const;
	/// domain check
	//virtual bool withinDomain(Metric3dContext& mc, const DPoint3d& pt, const DVector3d& pt_normal, DPoint2d& param, double * quality = nullptr) const;
	/// Create wire-frame visualization for rectangular surface patch covering the given set of points
	virtual MeshViewSet * createViewSetForPoints(MeshViewSet* set, const DataVector<DPoint2d> & points) const;
	/// Store XML description to stream
	virtual ostream& storeXML(ostream& os, const string& prefix = "") const;
	/// Returns the parameterization matrix for this surface and parameters [u,v]
	virtual ControlDataMatrix2d projectTransformationTensor(const DPoint2d& pt, const ControlDataMatrix3d& cdm3d) const;
	/// Returns the parameterization matrix for this surface and parameters [u,v]
	virtual ControlDataMatrix2d countParameterizationMatrix(const DPoint2d& pt, double & g_ratio) const;
	/// Returns the principle curvature direction for this surface and parameters [u,v]
	virtual const DVector3d getCurvatureDirection(const DPoint2d& pt, double & g_ratio) const;
	/// Returns the principle curvature for this surface and parameters [u,v] (angle and principle values)
	virtual const SurfaceCurvature getCurvature(const DPoint2d& pt, double & g_ratio) const;
	/// Returns the point of the surface for the given parameters
	virtual const DPoint3d getPoint(const DPoint2d& param) const  = 0;
	/// Returns the parameters of the surface for the given point (numerical approximation)
	virtual const DPoint2d getParameters(const DPoint3d& point) const { 
		return getParametersNear(point, DPoint2d::zero); }
	/// Returns the parameters of the surface for the given point (numerical approximation) and signed distance (according to normal vecor)
	virtual bool getParametersAndSignedDistance(const DPoint3d& /* point */, DPoint2d& /*param*/, double & /* z */) const {
		assert(false); return false; }
	/// Returns the parameters of the surface for the given point with starting point (numerical approximation)
	virtual const DPoint2d getParametersNear(const DPoint3d& point, const DPoint2d& near_point) const;
	/// Returns the shape (curve) parameter of the surface for the given point with starting point (numerical approximation)
	virtual double getShapeParameters(const DPoint3d& point, 
		std::shared_ptr<const Curve2dParametric> shape, double near_t, double min_t, double max_t) const;
	/// Returns the segment parameter of the surface for the given point with starting point (numerical approximation)
	virtual double getSegmentParameters(const DPoint3d& point, const DPoint2d& pt0, const DPoint2d& pt1, double near_t, double min_t, double max_t) const;
	/// Returns the normal vector to surface for the given parameters
	virtual const DVector3d getNormalVector(const DPoint2d& param) const;
	/// Returns the derivative of normal vector for the given parameters
	virtual const DVector3d getNormalVectorDerivative(int /* deriv */, const DPoint2d& /* param */) const {
		assert(false);  return DVector3d::zero;
	}
	/// Returns the normal vector to surface for the given parameters
	virtual const DVector3d getNormalVectorForPoint3d(const DPoint3d& pt) const;
	/// Returns the derivative (specified) vector for the given parameters
	virtual const DVector3d getDerivative(int deriv, const DPoint2d& param) const = 0;
	/// Returns the object-specific type 
	virtual ElementType getType() const { return SURFACE_UNKNOWN; }
	/// Returns the length of a segment on the surface (for parametrically specified vertices)
	virtual double segmentLength(const DPoint2d& param_0, const DPoint2d& param_1) const;
	/// Returns the approximation of a segment on surfaces via polyline (array of points)
	virtual void getPolyLine(DataVector<DPoint3d> & polyline, const DPoint2d& a, const DPoint2d& b) const;
	/// Returns the approximation of a segment on surfaces via polyline (array of points)
	virtual void getPolyLine(DataVector<DPoint3d> & polyline, 
		std::shared_ptr<const Curve2dParametric> shape, double t0, double t1) const;
	virtual std::shared_ptr<SurfaceParametric> adjustedForFaces( const DataVector< MeshFace* > & /* mfaces */,
		const DataVector< MeshPoint3d* > & /* mpoints */, MeshFace* /* central_face */ = nullptr  ) const {
			assert( false ); return nullptr; }
	virtual bool setOrientationLikeFace( MeshFace* face );
	virtual bool withinParamRange( const DPoint2d& /* param */) const { return true; }
public:
	/// Checks whether the surfaces is properly definied
	bool isValid() const { return m_valid; }
	/// set fixed - not mutable
	void setFixed( bool fixed = true ) { m_fixed = fixed; }
	/// whether is fixed
	bool isFixed() const { return m_fixed; }
protected:
	/// Validity flag
	bool m_valid;
	/// Optional domain reference
	SurfaceDomain* m_domain;
	/// Whether this surface is "mutable" or "final-and-permanent"
	bool m_fixed;
};

typedef std::shared_ptr<SurfaceParametric> SurfacePtr;
typedef std::shared_ptr<const SurfaceParametric> SurfaceConstPtr;

class SurfaceParametricSet
{
public:
	SurfaceParametricSet(SurfaceConstPtr surface = nullptr);
	SurfaceParametricSet(std::shared_ptr<const SurfaceParametricSet> sset, SurfaceConstPtr surface = nullptr);
public:
	bool countCurvatureMetric(const DPoint3d& pt, double model_diameter, ControlDataMatrix3d& cdm) const;
	int addSurface(SurfaceConstPtr surface);
	bool containsSurface(SurfaceConstPtr surface) const;	
	//DPoint3d fitToNearest(Metric3dContext & mc, const DPoint3d& pt, const DVector3d& vn, SurfaceParametric** nearest_surface = nullptr) const;
	//DPoint3d fitToAll(Metric3dContext & mc, const DPoint3d& pt, const DVector3d& vn) const;
	int count() const;
	inline int countInt() const { return count(); }
	const SurfaceConstPtr& getSurface(int i) const { return m_surfaces[i]; }
private:
	DataVector< SurfaceConstPtr > m_surfaces;
};

typedef std::shared_ptr<SurfaceParametricSet> SurfaceSetPtr;
typedef std::shared_ptr<const SurfaceParametricSet> SurfaceSetConstPtr;


#endif // !defined(SURFACEPARAMETRIC_H__INCLUDED)
