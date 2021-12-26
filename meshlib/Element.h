/////////////////////////////////////////////////////////////////////////////
// Element.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2002-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#if !defined(ELEMENT_H__INCLUDED)
#define ELEMENT_H__INCLUDED

#include "common.h"

/**
 * This is an abstract class for all geometrical constructional objects.
 */
class CElement
{
public:
	/// Standard destructor (empty)
	virtual ~CElement(){}
	/// Returns, whether this element is properly defined (must be implemented in derived classes)
	virtual bool isValid() const = 0;
	/// Returns the type of element (should be reimplemented in derived classes)
	virtual int getType() const { return ELEMENT_UNKNOWN; }
};

#endif // !defined(ELEMENT_H__INCLUDED)
