/////////////////////////////////////////////////////////////////////////////
// SurfaceMulti.h 
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2009-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#if !defined(SURFACEMULTI_H__INCLUDED)
#define SURFACEMULTI_H__INCLUDED

#include <memory>

#include "DataVector.h"
#include "DPoint.h"
#include "SurfaceParametric.h"
#include "SurfacePlane.h"
#include "DEquation.h"

/**
 *  Surface patch combined of several surfaces with 
 *   base surface and consistent parameterization
 */
class SurfaceMulti : public SurfaceParametric
{
public:
	/// Standard constructor
	SurfaceMulti(std::shared_ptr<const SurfaceParametric> base_surface);
private:
	SurfaceMulti(const SurfaceMulti& ) = delete; // not available
	SurfaceMulti& operator=(const SurfaceMulti& ) = delete; // not available
public:
	/// Abstract domain class for multi-surface domain-check
	class Domain
	{
	public:
		virtual ~Domain() {}
	public:
		virtual double testWithin(const DPoint2d& param) const = 0;
	};
	/// Specialized domain class for a rectangular domain
	class DomainRect : public Domain
	{
	public:
		DomainRect(const DPoint2d& _middle, const DVector2d& _width, double margin = 0.1)
			: middle(_middle), width(_width), mwidth(_width * (1.0+margin)) {}
	public:
		double testWithin(const DPoint2d& param) const;
	private:
		DPoint2d middle;
		DVector2d width;
		DVector2d mwidth; // with margin
	};
	/// Abstract domain class for multi-surface domain-check
	class Reparameterization
	{
	public:
		virtual ~Reparameterization() {}
	public:
		virtual DPoint2d transform(const DPoint2d& param) const = 0;
	};
	/// Specialized domain class for planar reparameterization
	class ReparameterizationPlanar : public Reparameterization
	{
	public:
		ReparameterizationPlanar(const SurfacePlane& source, const SurfacePlane& destination);
	public:
		DPoint2d transform(const DPoint2d& param) const;
	private:
		DPoint2d transf_p;
		DVector2d transf_e0, transf_e1;
	};
public:
	virtual string getSimpleDescription() const override { return "multi-surface"; }
	/// Store XML description to stream
	virtual ostream& storeXML(ostream& os, const string& prefix = "") const;
	/// Returns the point of the surface for the given parameters
	virtual const DPoint3d getPoint(const DPoint2d& param) const;
	/// Returns the normal vector to surface for the given parameters
	virtual const DVector3d getNormalVector(const DPoint2d& param) const;
	/// Returns the derivative (specified) vector for the given parameters
	virtual const DVector3d getDerivative(int deriv, const DPoint2d& param) const; 
	/// Returns the object-specific type 
	virtual ElementType getType() const { return SURFACE_MULTI; }
public:
	/// Insert an additional sub-surface
	void insertSubSurface(std::shared_ptr<const SurfaceParametric> surface, Domain* domain, Reparameterization* reparam);
protected:
	/// Sub-surface for a sub-domain
	class SubSurface{
	public:
		SubSurface(std::shared_ptr<const SurfaceParametric> _surface = nullptr, 
			Domain* _domain = nullptr, Reparameterization* _rep = nullptr)
			: surface(_surface), domain(_domain), reparam(_rep) {}
		~SubSurface();
		/// Store XML description to stream
		ostream& storeXML(ostream& os, const string& prefix = "") const;
		/// Domain test
		double withinDomain(const DPoint2d& param) const { 
			assert(domain);
			return domain->testWithin(param); 
		}
		/// Reparameterization
		DPoint2d localParam(const DPoint2d& param) const {
			assert(reparam);
			return reparam->transform(param);
		}
		/// Get surface
		std::shared_ptr<const SurfaceParametric> getSurface() const { return surface; }
	private:
		std::shared_ptr<const SurfaceParametric> surface;
		Domain* domain;
		Reparameterization* reparam;
	};
protected:
	/// Base surface
	std::shared_ptr<const SurfaceParametric> m_base_surface;
	/// Sub-surfaces
	DataVector<std::shared_ptr<SubSurface>> m_surfaces;
};

#endif // !defined(SURFACEMULTI_H__INCLUDED)
