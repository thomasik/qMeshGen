/////////////////////////////////////////////////////////////////////////////
// MeshArea.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2002-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#if !defined(MESHAREA_H__INCLUDED)
#define MESHAREA_H__INCLUDED

#include "MeshElement.h"
#include "DRect.h"
#include "DataVector.h"

class MeshPoint2d;
class SurfaceParametric;
class MeshViewSet;

/**
 * This class describes a user-given 2d-domain
 * used for defining the area to discretize.
 */
class MeshArea : public MeshElement  
{
public:
	/// Standard constructor creating new mesh-area described by set of points
	MeshArea(int ct = 0, MeshPoint2d** new_points = nullptr, bool filled = true, double k = 1.0);
	/// Standard constructor creating new mesh-area described by set of points
	MeshArea(DataVector<MeshPoint2d*> &new_points, bool filled = true, double k = 1.0);
	/// Standard destructor
	virtual ~MeshArea();
public:
	/// Returns the description of this face for visualization
	MeshViewSet* getViewSet(MeshViewSet* set = nullptr, 
		std::shared_ptr<const SurfaceParametric> surface = nullptr) const;
	/// Returns the lenght of smallest edge in the boundary of this area (if required 3D information is used)
	double getShortestEdge(std::shared_ptr<const SurfaceParametric> surface) const;
	/// Returns the bounding rectangle which contains the whole area
	DRect getBoundingRect() const;
	/// clear some data for faster all-delete process
	virtual void preDeleteAll();
	/// Not used in this derived class
	virtual double countQuality(Metric2dContext& /* mc */, bool /* mind_area */ = false, int /* criterion  */ = -1)
		{ assert(false); return 0.0; }
	/// Not used in this derived class
	virtual double getMinArea() const { assert(false); return 0.0; }
	/// Not used in this derived class
	virtual double getArea(Metric2dContext& /* mc */, bool /* local_metric */ = true) const { assert(false); return 0.0; }
	/// Not used in this derived class
	virtual double getAreaNoMetric() const { assert(false); return 0.0; }
	/// Not used in this derived class
	virtual double getAngle(Metric2dContext& /* mc */, int /* i */, bool /* local_metric */ = true) const { assert(false); return 0.0; }
	/// Not used in this derived class
	virtual double getAlphaQuality(Metric2dContext& /* mc */, bool /* local_metric */ = true) const { assert(false); return 0.0; }
	/// Not used in this derived class
	virtual bool isPointInside(const DPoint2d& /* mpt */) const { assert(false); return false;}
	/// Not used in this derived class
	virtual void exportToGRD(FILE* /* file */) const{ assert(false); }
	/// Not used in this derived class
	virtual bool isInverted() const { assert(false); return false; }
	/// Returns the type of this object
	virtual ElementType getType() const { return ELEMENT_MESH_AREA; }
	/// Draws this area using OpenGL routines
	virtual void drawGL() const;
	/// Checks whether the area is filled or defines a hole
	virtual bool isFilled() const { return filled; }
	/// Sets "filled" info for this area
	void setFilled(bool f = true) { filled = f; }
	/// Returns the material-coefficient for this area
	double getMaterial() const { return k_material; }
	/// Sets the material-coefficient for this area
	void setMaterial(double k){ k_material = k; }
public:
	/// Checks the orientation of the given polygon (returns true if anticlockwise)
	static bool getOrientation(const DataVector<MeshPoint2d*> &points);
protected:
	/// Boolean variable defining whether the area is treated as "filled"
	bool		filled;
	/// Material-coefficient for this area
	double		k_material;
};

#endif // !defined(MESHAREA_H__INCLUDED)
