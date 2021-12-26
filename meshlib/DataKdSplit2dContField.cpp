#include "DataKdSplit2dContField.h"

void Kd2dCFMaxDiffGradientSplitter::whereToSplit(const DRect & box, CDM2d(*f)(const DPoint2d &pt), Axis & split_axis, double & split_value)
{
	prepareCache(box, f);

	double diff, max_avg_diff = 0.0;
	for (int i = 0; i < csize; i++) {
		for (int j = 0; j < csize; j++) {
			if (i > 0) {
				diff = diffKdValue(cacheData(i, j), cacheData(i - 1, j));
				if (diff > max_avg_diff) {
					max_avg_diff = diff;
					split_axis = Axis::X;
					split_value = box.x0 + (i - 0.5) * cdx;
				}
			}
			if (j > 0) {
				diff = diffKdValue(cacheData(i, j), cacheData(i, j - 1));
				if (diff > max_avg_diff) {
					max_avg_diff = diff;
					split_axis = Axis::Y;
					split_value = box.y0 + (j - 0.5) * cdy;
				}
			}
		}
	}
}

void Kd2dCFMaxGradientSplitter::whereToSplit(const DRect & box, CDM2d(*f)(const DPoint2d &pt), Axis & split_axis, double & split_value)
{
	prepareCache(box, f);

	double max_diff = 0.0;
	for (int i = 0; i < csize; i++) {
		for (int j = 0; j < csize; j++) {
			if (i > 0) {
				double diff = diffKdValue(cacheData(i, j), cacheData(i - 1, j));
				if (diff > max_diff) {
					max_diff = diff;
					split_axis = Axis::X;
					split_value = box.x0 + i * cdx;
				}
			}
			if (j > 0) {
				double diff = diffKdValue(cacheData(i, j), cacheData(i, j - 1));
				if (diff > max_diff) {
					max_diff = diff;
					split_axis = Axis::Y;
					split_value = box.y0 + j * cdy;
				}
			}
		}
	}
}

void Kd2dCFMaxGradientSumSplitter::whereToSplit(const DRect & box, CDM2d(*f)(const DPoint2d &pt), Axis & split_axis, double & split_value)
{
	prepareCache(box, f);

	double sum[] = { 0, 0 };
	for (int i = 0; i < csize; i++) {
		for (int j = 0; j < csize; j++) {
			if (i > 0) sum[0] += diffKdValue(cacheData(i, j), cacheData(i - 1, j));
			if (j > 0) sum[1] += diffKdValue(cacheData(i, j), cacheData(i, j - 1));
		}
	}

	double d0[] = { box.x0, box.y0 };
	double step[] = { cdx, cdy };
	split_axis = (sum[0] > sum[1]) ? Axis::X : Axis::Y;
	int dim = (int)split_axis;

	double sum_half = 0;
	split_value = d0[dim] + step[dim] * csize / 2;
	for (int i = 1; i < csize; i++) {
		for (int j = 0; j < csize; j++) {
			int di, dj;
			int pdi, pdj;
			switch (split_axis) {
			case Axis::X:	di = i; dj = j; pdi = i - 1; pdj = j; break;
			case Axis::Y:	di = j; dj = i; pdi = j; pdj = i - 1; break;
			}
			sum_half += diffKdValue(cacheData(di, dj), cacheData(pdi, pdj));
		}
		if (2 * sum_half > sum[dim]) {
			split_value = d0[dim] + step[dim] * i;
			break;
		}
	}
}
