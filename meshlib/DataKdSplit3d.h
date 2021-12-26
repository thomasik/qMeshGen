/////////////////////////////////////////////////////////////////////////////
// DataKdSplit3d.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2017-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#ifndef DATAKD3DSPLIT_H__INCLUDED
#define DATAKD3DSPLIT_H__INCLUDED

#include "common.h"

#include "DPoint.h"
#include "DRect.h"
#include "DataVector.h"
#include "DMetric3d.h"
#include "MeshData.h"

class Kd3dSplitter
{
public:
	Kd3dSplitter(int _csize = 0) : csize(_csize), cache(NULL) {
		if (_csize > 0) createCache();
	}
	virtual ~Kd3dSplitter() {
		if (cache != NULL) delete[] cache;
	}
public:
	virtual int getSamplingSize() const { return csize; }
	virtual string getLabel() const = 0;
protected:
	void createCache() {
		assert(csize > 1);
		cache = new CDM3d[csize * csize * csize];
	}
	CDM3d& cacheData(int i, int j, int k) {
		assert(cache != NULL);
		return cache[i*csize*csize + j * csize + k];
	}
	void prepareCache(const DBox & box, const std::function<CDM3d(const DPoint3d & pt)> & f);
protected:
	const int csize;
	CDM3d* cache;
	double cdx, cdy, cdz;
};

#endif // !defined(DATAKD3DSPLIT_H__INCLUDED)
