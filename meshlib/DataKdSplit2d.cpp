#include "DataKdSplit2d.h"

void Kd2dSplitter::prepareCache(const DRect & box, CDM2d(*f)(const DPoint2d &pt)) 
{
	double fx = 1.0 / (double)csize;
	cdx = box.getDX() * fx;
	cdy = box.getDY() * fx;

	DPoint2d pt(box.x0 + 0.5*cdx, box.y0);
	int ix = 0, iy;
	for (; pt.x < box.x1; pt.x += cdx, ix++) {
		for (iy = 0, pt.y = box.y0 + 0.5*cdy; pt.y < box.y1; pt.y += cdy, iy++) {
			cacheData(ix, iy) = f(pt);
		}
	}
}
