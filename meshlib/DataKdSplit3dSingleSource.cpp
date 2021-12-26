#include "DataKdSplit3dSingleSource.h"
#include "ControlSpace3dKdTree.h"

void Kd3dSSGoldenRatioLongestSplitter::whereToSplit(const DBox & box, const DPoint3d & pt, const DVector3d & dlen2, Axis & split_axis, double & split_value)
{
	// longest axis, but with metric (metric lengths stored in dlen)
	if (dlen2.x > dlen2.y && dlen2.x > dlen2.z) split_axis = Axis::X;
	else if (dlen2.y > dlen2.z) split_axis = Axis::Y;
	else split_axis = Axis::Z;
	//split_axis = box.getLongestAxis();
	const double GR1 = 2.0 / (1 + sqrt(5.0));
	const double GR0 = 1.0 - GR1;
	double dp = pt[split_axis];
	double d0 = box.getD0(split_axis);
	double d10 = box.getDLen(split_axis);
	double t = (dp - d0) / d10;
	t = (abs(t - GR0) < abs(t - GR1)) ? GR0 : GR1;
	split_value = d0 + t * d10;
}

void Kd3dSSMinLongestSplitter::whereToSplit(const DBox & box, const DPoint3d & pt, const DVector3d & dlen2, Axis & split_axis, double & split_value)
{
	// longest axis, but with metric (metric lengths stored in dlen)
	if (dlen2.x > dlen2.y && dlen2.x > dlen2.z) split_axis = Axis::X;
	else if (dlen2.y > dlen2.z) split_axis = Axis::Y;
	else split_axis = Axis::Z;
	double len = sqrt(dlen2[split_axis]);
	//split_axis = box.getLongestAxis();
	double dp = pt[split_axis];
	double d0 = box.getD0(split_axis);
	double d10 = box.getDLen(split_axis);
	double MR = KdTree3dParams::paramMetricLenRatio;
	assert(d10 > MR);
	double t = (dp - d0) / d10;
	double t_len = t * len;
	double k = MR / len;
	if (len < 2 * MR)
		split_value = d0 + 0.5*d10;
	else if (t_len <= MR)
		split_value = d0 + k*d10;
	else if (len - t_len <= MR)
		split_value = d0 + (1.0 - k)*d10;
	else if (t < 0.5)
		split_value = d0 + (t + k/2)*d10; // len==1.0 (metric) after _t_
	else
		split_value = d0 + (t - k/2)*d10; // len==1.0 (metric) before _t_
}
