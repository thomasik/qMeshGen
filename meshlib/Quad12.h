// Quad12.h: interface for the Quad12 class.
//
//////////////////////////////////////////////////////////////////////

#pragma once

#if !defined(QUAD12_H__INCLUDED)
#define QUAD12_H__INCLUDED

#include "DPoint.h"
class MeshQuad2d;

class Quad12  
{
public:
	Quad12(MeshQuad2d* quad = nullptr) { if(quad) setData(quad); }
public:
	static void svdfit(DPoint2d* points, double* z, int count, double* a, int a_count);
	static double approximationFunction(int id, const DPoint2d& pt);
	static double approximateValue(int count, DPoint2d* points, double *z, const DPoint2d & where, int fun_ct = 10);
	void getDerivativePoint(const DPoint2d & param, DPoint2d& ptx, DPoint2d& pty);
	void test();
	void getValueHesjan(const DPoint2d & param, double& fxx, double& fxy, double& fyy);
	double getValue(const DPoint2d & param);
	DPoint2d getPoint(const DPoint2d & param);
	double interpolateQuadValue(const DPoint2d & pt);
	//--
	void setData(MeshQuad2d* quad);
protected:
	DPoint2d points[12];
	double values[12];
};

#endif // !defined(AFX_QUAD12_H__6F1F5F27_5619_44F0_ACC9_BDBA2F211201__INCLUDED_)
