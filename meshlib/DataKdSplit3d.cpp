#include "DataKdSplit3d.h"

void Kd3dSplitter::prepareCache(const DBox & box, const std::function<CDM3d(const DPoint3d & pt)> & f)
{
	double fx = 1.0 / (double)csize;
	cdx = box.getDX() * fx;
	cdy = box.getDY() * fx;
	cdz = box.getDZ() * fx;

	DPoint3d pt(box.x0 + 0.5*cdx, box.y0, box.z0);
	int ix = 0, iy, iz;
	for (; pt.x < box.x1; pt.x += cdx, ix++) {
		for (iy = 0, pt.y = box.y0 + 0.5*cdy; pt.y < box.y1; pt.y += cdy, iy++) {
			for (iz = 0, pt.z = box.z0 + 0.5*cdz; pt.z < box.z1; pt.z += cdz, iz++) {
				cacheData(ix, iy, iz) = f(pt);
			}
		}
	}
}
