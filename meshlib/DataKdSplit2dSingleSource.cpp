#include "DataKdSplit2dSingleSource.h"

void Kd2dSSGoldenRatioLongestSplitter::whereToSplit(const DRect & box, const DPoint2d & pt, const DVector2d & dlen2, Axis & split_axis, double & split_value)
{
	// longest axis, but with metric (metric lengths stored in dlen)
	split_axis = (dlen2.x > dlen2.y) ? Axis::X : Axis::Y;
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

void Kd2dSSMinLongestSplitter::whereToSplit(const DRect & box, const DPoint2d & pt, const DVector2d & dlen2, Axis & split_axis, double & split_value)
{
	// longest axis, but with metric (metric lengths stored in dlen)
	split_axis = (dlen2.x > dlen2.y) ? Axis::X : Axis::Y;
	double len = sqrt(dlen2[split_axis]);
	//split_axis = box.getLongestAxis();
	double dp = pt[split_axis];
	double d0 = box.getD0(split_axis);
	double d10 = box.getDLen(split_axis);
	double t = (dp - d0) / d10;
	double t_len = t * len;
	if (t_len <= 2.0)
		split_value = d0 + t*d10;
	else if (len - t_len <= 2.0)
		split_value = d0 + (1.0 - t)*d10;
	else if (t < 0.5)
		split_value = d0 + (t*(t_len - 1.0) / t_len)*d10; // len==1.0 (metric) before _t_
	else
		split_value = d0 + (1.0 - (1.0 - t)*(t_len - 1.0) / t_len)*d10; // len==1.0 (metric) after _t_
}
