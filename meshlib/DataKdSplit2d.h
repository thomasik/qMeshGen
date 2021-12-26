/////////////////////////////////////////////////////////////////////////////
// DataKdSplit2d.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2017-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#ifndef DATAKD2DSPLIT_H__INCLUDED
#define DATAKD2DSPLIT_H__INCLUDED

#include "common.h"

#include "DPoint.h"
#include "DRect.h"
#include "DataVector.h"
#include "DMetric2d.h"

class Kd2dSplitter
{
public:
	Kd2dSplitter(int _csize = 0) : csize(_csize), cache(NULL) {
		if (_csize > 0) createCache();
	}
	virtual ~Kd2dSplitter() {
		if (cache != NULL) delete[] cache;
	}
protected:
	void createCache() {
		assert(csize > 1);
		cache = new CDM2d[csize * csize];
	}
	CDM2d& cacheData(int i, int j) {
		assert(cache != NULL);
		return cache[i*csize + j];
	}
	void prepareCache(const DRect & box, CDM2d(*f)(const DPoint2d & pt));
protected:
	const int csize;
	CDM2d* cache;
	double cdx, cdy;
};

#endif // !defined(DATAKD2DSPLIT_H__INCLUDED)
