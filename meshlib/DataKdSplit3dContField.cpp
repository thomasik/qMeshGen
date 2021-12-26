#include "DataKdSplit3dContField.h"

void Kd3dCFMaxGradientSumSplitter::whereToSplit(const DBox & box, 
	const std::function<CDM3d(const DPoint3d & pt)> & f, 
	Axis & split_axis, double & split_value)
{
	prepareCache(box, f);

	double sum[] = { 0, 0, 0 };
	for (int i = 0; i < csize; i++) {
		for (int j = 0; j < csize; j++) {
			for (int k = 0; k < csize; k++) {
				if (i > 0) sum[0] += diffKdValue(cacheData(i, j, k), cacheData(i - 1, j, k));
				if (j > 0) sum[1] += diffKdValue(cacheData(i, j, k), cacheData(i, j - 1, k));
				if (k > 0) sum[2] += diffKdValue(cacheData(i, j, k), cacheData(i, j, k - 1));
			}
		}
	}

	double d0[] = { box.x0, box.y0, box.z0 };
	double step[] = { cdx, cdy, cdz };
	if (sum[0] > sum[1] && sum[0] > sum[2]) split_axis = Axis::X;
	else split_axis = (sum[1] > sum[2]) ? Axis::Y : Axis::Z;
	int dim = (int)split_axis;

	double sum_half = 0;
	split_value = d0[dim] + step[dim] * csize / 2;
	for (int i = 1; i < csize; i++) {
		for (int j = 0; j < csize; j++) {
			for (int k = 0; k < csize; k++) {
				int di, dj, dk;
				int pdi, pdj, pdk;
				switch (split_axis) {
				case Axis::X:	di = i; dj = j; dk = k; pdi = i - 1; pdj = j; pdk = k; break;
				case Axis::Y:	di = k; dj = i; dk = j; pdi = k; pdj = i - 1; pdk = j; break;
				case Axis::Z:	di = j; dj = k; dk = i; pdi = j; pdj = k; pdk = i - 1; break;
				}
				sum_half += diffKdValue(cacheData(di, dj, dk), cacheData(pdi, pdj, pdk));
			}
		}
		if (2 * sum_half > sum[dim]) {
			split_value = d0[dim] + step[dim] * i;
			break;
		}
	}
}

string Kd3dCFMaxGradientSumSplitter::getLabel() const {
	string label = "MGradSum";
	label += to_string(getSamplingSize());
	return label;
}

void Kd3dCFMaxGradientSplitter::whereToSplit(const DBox & box, 
	const std::function<CDM3d(const DPoint3d & pt)> & f, 
	Axis & split_axis, double & split_value)
{
	prepareCache(box, f);

	double max_diff = 0.0;
	for (int i = 0; i < csize; i++) {
		for (int j = 0; j < csize; j++) {
			for (int k = 0; k < csize; k++) {
				if (i > 0) {
					double diff = diffKdValue(cacheData(i, j, k), cacheData(i - 1, j, k));
					if (diff > max_diff) {
						max_diff = diff;
						split_axis = Axis::X;
						split_value = box.x0 + i * cdx;
					}
				}
				if (j > 0) {
					double diff = diffKdValue(cacheData(i, j, k), cacheData(i, j - 1, k));
					if (diff > max_diff) {
						max_diff = diff;
						split_axis = Axis::Y;
						split_value = box.y0 + j * cdy;
					}
				}
				if (k > 0) {
					double diff = diffKdValue(cacheData(i, j, k), cacheData(i, j, k - 1));
					if (diff > max_diff) {
						max_diff = diff;
						split_axis = Axis::Z;
						split_value = box.z0 + k * cdz;
					}
				}
			}
		}
	}
}

string Kd3dCFMaxGradientSplitter::getLabel() const {
	string label = "MGrad";
	label += to_string(getSamplingSize());
	return label;
}

void Kd3dCFMaxDiffGradientSplitter::whereToSplit(const DBox & box, 
	const std::function<CDM3d(const DPoint3d & pt)> & f, 
	Axis & split_axis, double & split_value)
{
	prepareCache(box, f);

	double diff, max_avg_diff = 0.0;
	for (int i = 0; i < csize; i++) {
		for (int j = 0; j < csize; j++) {
			for (int k = 0; k < csize; k++) {
				if (i > 0) {
					diff = diffKdValue(cacheData(i, j, k), cacheData(i - 1, j, k));
					if (diff > max_avg_diff) {
						max_avg_diff = diff;
						split_axis = Axis::X;
						split_value = box.x0 + (i - 0.5) * cdx;
					}
				}
				if (j > 0) {
					diff = diffKdValue(cacheData(i, j, k), cacheData(i, j - 1, k));
					if (diff > max_avg_diff) {
						max_avg_diff = diff;
						split_axis = Axis::Y;
						split_value = box.y0 + (j - 0.5) * cdy;
					}
				}
				if (k > 0) {
					diff = diffKdValue(cacheData(i, j, k), cacheData(i, j, k - 1));
					if (diff > max_avg_diff) {
						max_avg_diff = diff;
						split_axis = Axis::Z;
						split_value = box.z0 + (k - 0.5) * cdz;
					}
				}
			}
		}
	}
}

string Kd3dCFMaxDiffGradientSplitter::getLabel() const {
	string label = "MDiffGrad";
	label += to_string(getSamplingSize());
	return label;
}
