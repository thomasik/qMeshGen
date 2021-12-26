/////////////////////////////////////////////////////////////////////////////
// TagBorder.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2010-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#if !defined(TAGBORDER_H__INCLUDED)
#define TAGBORDER_H__INCLUDED

/**
 * Class defines an extra data "border" for all mesh entities which requires this kind of tag
 */
class TagBorder
{
public:
	/// Constructor
	TagBorder() : border(NONE) {}
	/// Copying constructor
	TagBorder(const TagBorder& tb) : border(tb.border) { }
public:
	/// Available types of border flags
	enum BorderType { NONE = 0, OUTER = 1, INNER = 2, ANY = 3,
		RIDGE = 4, CORNER = 8, FIXED = 16 };
public:
	/// Marks this object as a boundary one
	void copyBorderFlagsFrom(const TagBorder* other){ assert(other); border = other->border; }
	/// Set the given border flag for this edge
	void setBorder(char bflag = TagBorder::OUTER){ border = bflag; }
	/// Clear all border flag from this edge
	void clearBorder(){ border = TagBorder::NONE; }
	/// Checks whether this edge belongs to the boundary
	bool isBorder(char bflag = TagBorder::ANY) const{ return (border & bflag) != 0; }
	/// Returns border flags
	char getBorderFlags() const { return border; }
	/// Add the given border flag for this edge
	void setBorderFlags(char bflag){ border |= bflag; }
	/// Clear the given border flag from this edge
	void clearBorderFlags(char bflag){ border = (char)(border & ~bflag); }
protected:
	/// border data
	char border;
};

#endif // !defined(TAGBORDER_H__INCLUDED)
