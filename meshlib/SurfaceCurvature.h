/////////////////////////////////////////////////////////////////////////////
// SurfaceCurvature.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2015-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#if !defined(SURFACECURVATURE_H__INCLUDED)
#define SURFACECURVATURE_H__INCLUDED

class SurfaceCurvature {
public:
	SurfaceCurvature(bool _valid = false) : c1(0.0), c2(0.0), s1(1), s2(1), angle(0.0), valid(_valid) {}
	SurfaceCurvature( const double _c1, const double _c2, const double _angle = 0.0)
		: c1(_c1), c2(_c2), s1(1), s2(1), angle(_angle), valid(true) 
	{
		if( c1 < 0.0 ) { c1 = -c1; s1 = -1; }
		if( c2 < 0.0 ) { c2 = -c2; s2 = -1; }
	}
public:
	double c1, c2; // curvature value
	signed char s1, s2; // sing of curvature
	double angle; // rotation angle
	bool valid;
};

#endif // !defined(SURFACECURVATURE_H__INCLUDED)
