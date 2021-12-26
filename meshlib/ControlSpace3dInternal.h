/////////////////////////////////////////////////////////////////////////////
// ControlSpace3dInternal.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2010-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#if !defined(CONTROLSPACE3DINTERNAL_H__INCLUDED)
#define CONTROLSPACE3DINTERNAL_H__INCLUDED

#include "MeshData.h"
#include "ControlSpace3d.h"
#include "DTriangle.h"

/**
 * This class implements an internal, predefined control space
 */
class ControlSpace3dInternal : public ControlSpace3d
{
public:
	/// Standard contructor
	ControlSpace3dInternal();
	/// Destructor
	virtual ~ControlSpace3dInternal();
public:
	/// Returns the type of element (should be reimplemented in derived classes)
	virtual int getType() const { return MeshData::CONTROL_INTERNAL_3D; }
	/// Get sizing info (matrix mode) at the given point
	virtual ControlDataMatrix3d getMetricAtPoint(const DPoint3d& pt) const;
};

#endif // !defined(CONTROLSPACE3DINTERNAL_H__INCLUDED)
