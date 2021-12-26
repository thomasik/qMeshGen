// Quad12.cpp: implementation of the Quad12 class.
//
//////////////////////////////////////////////////////////////////////

#include "Quad12.h"
#include "MeshQuad2d.h"
#include "MeshEdge2d.h"
#include "MeshPoint2d.h"
#include "nrutil.h"

void Quad12::setData(MeshQuad2d* quad)
{
	MeshPoint2d* quad_point = quad->getPoint(0);
	points[0] = quad_point->getCoordinates();
	LOG4CPLUS_ERROR(MeshLog::logger_console,   "Needs rewriting for extra weight in point");
/*
	values[0] = 0.0; // quad_point->getWeight();
	MeshEdge2d* quad_edge = quad->getEdge(0);
	MeshPoint2d* pnt = quad_edge->getInnerMeshPoint(0, quad_point);
	points[4] = pnt->getCoordinates();
	values[4] = 0.0; //pnt->getWeight();
	pnt = quad_edge->getInnerMeshPoint(1, quad_point);
	points[5] = pnt->getCoordinates();
	values[5] = 0.0; //pnt->getWeight();

	quad_point = quad->getPoint(2);
	points[2] = quad_point->getCoordinates();
	values[2] = 0.0; //quad_point->getWeight();
	quad_edge = quad->getEdge(2);
	pnt = quad_edge->getInnerMeshPoint(0, quad_point);
	points[6] = pnt->getCoordinates();
	values[6] = 0.0; //pnt->getWeight();
	pnt = quad_edge->getInnerMeshPoint(1, quad_point);
	points[7] = pnt->getCoordinates();
	values[7] = 0.0; //pnt->getWeight();

	quad_point = quad->getPoint(1);
	points[1] = quad_point->getCoordinates();
	values[1] = 0.0; //quad_point->getWeight();
	quad_edge = quad->getEdge(1);
	pnt = quad_edge->getInnerMeshPoint(0, quad_point);
	points[9] = pnt->getCoordinates();
	values[9] = 0.0; //pnt->getWeight();
	pnt = quad_edge->getInnerMeshPoint(1, quad_point);
	points[10] = pnt->getCoordinates();
	values[10] = 0.0; //pnt->getWeight();

	quad_point = quad->getPoint(3);
	points[3] = quad_point->getCoordinates();
	values[3] = 0.0; //quad_point->getWeight();
	quad_edge = quad->getEdge(3);
	pnt = quad_edge->getInnerMeshPoint(0, quad_point);
	points[11] = pnt->getCoordinates();
	values[11] = 0.0; //pnt->getWeight();
	pnt = quad_edge->getInnerMeshPoint(1, quad_point);
	points[8] = pnt->getCoordinates();
	values[8] = 0.0; //pnt->getWeight();
*/
}

double Quad12::interpolateQuadValue(const DPoint2d & pt)
{
	double	F_ERR		= 1e-10;	// Maksymalny b³¹d funkcji
	double	T_ERR		= 1e-15;	// Maksymalny b³¹d parametru
	int		MAX_STEPS	= 30;	// Maksymalna iloœæ kroków

	DPoint2d x(0.0, 0.0);
	DPoint2d ptx, pty;
	for(int k = 0; k < MAX_STEPS; k++){
		DPoint2d f = getPoint(x);
		DVector2d df = pt - f;
		if(abs(df.x) + abs(df.y) < F_ERR) return getValue(x);
		getDerivativePoint(x, ptx, pty);
		double W = ptx.x * pty.y - pty.x * ptx.y;
		if(W == 0.0){
			LOG4CPLUS_ERROR(MeshLog::logger_console,   "NR - W == 0");
			return getValue(x);	// ? error
		}
		double dx = (df.x * pty.y - df.y * pty.x) / W;
		double dy = (ptx.x * df.y - df.x * ptx.y) / W;
		x.x += dx;
		x.y += dy;
		if(abs(dx) + abs(dy) < T_ERR) return getValue(x);

		if(x.x < -1.0){ x.x = -1.0; }
		if(x.x >  1.0){ x.x =  1.0; }
		if(x.y < -1.0){ x.y = -1.0; }
		if(x.y >  1.0){ x.y =  1.0; }
	}

	return getValue(x);	// ? error
}

DPoint2d Quad12::getPoint(const DPoint2d & param)
{	
	double Xp = 1.0 + param.x;
	double Xm = 1.0 - param.x;
	double Yp = 1.0 + param.y;
	double Ym = 1.0 - param.y;

	DPoint2d value = points[0] * (Xm * Ym);
	value.add(points[1], Xp * Ym);
	value.add(points[2], Xp * Yp);
	value.add(points[3], Xm * Yp);
	value *= 9.0 * (param.x*param.x + param.y*param.y) - 10.0;

	DPoint2d v = points[4] * (Ym * (1.0 - 3.0 * param.x));
	v.add(points[5], Ym * (1.0 + 3.0 * param.x));
	v.add(points[6], Yp * (1.0 + 3.0 * param.x));
	v.add(points[7], Yp * (1.0 - 3.0 * param.x));
	value.add(v, 9.0 * (1.0 - param.x*param.x));

	v  = points[8]  * (Xm * (1.0 - 3.0 * param.y));
	v.add(points[9],  Xp * (1.0 - 3.0 * param.y));
	v.add(points[10], Xp * (1.0 + 3.0 * param.y));
	v.add(points[11], Xm * (1.0 + 3.0 * param.y));
	value.add(v, 9.0 * (1.0 - param.y*param.y));

	value *= 1.0 / 32.0;

	return value;
}

double Quad12::getValue(const DPoint2d & param)
{
	double Xp = 1.0 + param.x;
	double Xm = 1.0 - param.x;
	double Yp = 1.0 + param.y;
	double Ym = 1.0 - param.y;

	double value = Xm * Ym * values[0];
	value += Xp * Ym * values[1];
	value += Xp * Yp * values[2];
	value += Xm * Yp * values[3];
	value *= 9.0 * (param.x*param.x + param.y*param.y) - 10.0;

	double v = Ym * (1.0 - 3.0 * param.x) * values[4];
	v += Ym * (1.0 + 3.0 * param.x) * values[5];
	v += Yp * (1.0 + 3.0 * param.x) * values[6];
	v += Yp * (1.0 - 3.0 * param.x) * values[7];
	value += 9.0 * (1.0 - param.x*param.x) * v;

	v  = Xm * (1.0 - 3.0 * param.y) * values[8];
	v += Xp * (1.0 - 3.0 * param.y) * values[9];
	v += Xp * (1.0 + 3.0 * param.y) * values[10];
	v += Xm * (1.0 + 3.0 * param.y) * values[11];
	value += 9.0 * (1.0 - param.y*param.y) * v;

	value *= 1.0 / 32.0;

	return value;
}

void Quad12::getValueHesjan(const DPoint2d & param, double& fxx, double& fxy, double& fyy)
{
	double X = param.x;
	double Y = param.y;
	double Xp = 1.0 + X;
	double Xm = 1.0 - X;
	double Yp = 1.0 + Y;
	double Ym = 1.0 - Y;
	double X3m = 1.0 - 3.0 * X;
	double X3p = 1.0 + 3.0 * X;
	double Y3m = 1.0 - 3.0 * Y;
	double Y3p = 1.0 + 3.0 * Y;
	double XY = 1.0 / 32.0 * (9.0 * (X*X + Y*Y) - 10.0);
	double Xm2 = 27.0 / 32.0 * (1.0 - X*X);
	double Ym2 = 27.0 / 32.0 * (1.0 - Y*Y);

	fxx  = values[0] * (9.0 / 16.0 * Xm * Ym - 9.0 / 8.0 * X * Ym);
	fxx += values[1] * (9.0 / 16.0 * Xp * Ym + 9.0 / 8.0 * X * Ym);
	fxx += values[2] * (9.0 / 16.0 * Xp * Yp + 9.0 / 8.0 * X * Yp);
	fxx += values[3] * (9.0 / 16.0 * Xm * Yp - 9.0 / 8.0 * X * Yp);
	//
	fxx -= values[4] * (9.0 / 16.0 * X3m * Ym - 27.0 / 8.0 * X * Ym);
	fxx -= values[5] * (9.0 / 16.0 * X3p * Ym + 27.0 / 8.0 * X * Ym);
	fxx -= values[6] * (9.0 / 16.0 * X3p * Yp + 27.0 / 8.0 * X * Yp);
	fxx -= values[7] * (9.0 / 16.0 * X3m * Yp - 27.0 / 8.0 * X * Yp);

	fyy  = values[0] * (9.0 / 16.0 * Xm * Ym - 9.0 / 8.0 * Y * Xm);
	fyy += values[1] * (9.0 / 16.0 * Xp * Ym - 9.0 / 8.0 * Y * Xp);
	fyy += values[2] * (9.0 / 16.0 * Xp * Yp + 9.0 / 8.0 * Y * Xp);
	fyy += values[3] * (9.0 / 16.0 * Xm * Yp + 9.0 / 8.0 * Y * Xm);

	fyy -= values[8]  * (9.0 / 16.0 * Y3m * Xm - 27.0 / 8.0 * Y * Xm);
	fyy -= values[9]  * (9.0 / 16.0 * Y3m * Xp - 27.0 / 8.0 * Y * Xp);
	fyy -= values[10] * (9.0 / 16.0 * Y3p * Xp + 27.0 / 8.0 * Y * Xp);
	fyy -= values[11] * (9.0 / 16.0 * Y3p * Xm + 27.0 / 8.0 * Y * Xm);

	fxy  = values[0] * (-9.0 / 16.0 * Xm * X - 9.0 / 16.0 * Ym * Y + XY);
	fxy += values[1] * (-9.0 / 16.0 * Xp * X + 9.0 / 16.0 * Ym * Y - XY);
	fxy += values[2] * ( 9.0 / 16.0 * Xp * X + 9.0 / 16.0 * Yp * Y + XY);
	fxy += values[3] * ( 9.0 / 16.0 * Xm * X - 9.0 / 16.0 * Yp * Y - XY);

	fxy += values[4] * ( 9.0 / 16.0 * X3m * X + Xm2);
	fxy += values[5] * ( 9.0 / 16.0 * X3p * X - Xm2);
	fxy += values[6] * (-9.0 / 16.0 * X3p * X + Xm2);
	fxy += values[7] * (-9.0 / 16.0 * X3m * X - Xm2);

	fxy += values[8]  * ( 9.0 / 16.0 * Y3m * Y + Ym2);
	fxy += values[9]  * (-9.0 / 16.0 * Y3m * Y - Ym2);
	fxy += values[10] * (-9.0 / 16.0 * Y3p * Y + Ym2);
	fxy += values[11] * ( 9.0 / 16.0 * Y3p * Y - Ym2);
	// end.
}

void Quad12::getDerivativePoint(const DPoint2d & param, DPoint2d& ptx, DPoint2d& pty)
{
	double X = param.x;
	double Y = param.y;
	double Xp = 1.0 + X;
	double Xm = 1.0 - X;
	double Yp = 1.0 + Y;
	double Ym = 1.0 - Y;
	double X3m = 1.0 - 3.0 * X;
	double X3p = 1.0 + 3.0 * X;
	double Y3m = 1.0 - 3.0 * Y;
	double Y3p = 1.0 + 3.0 * Y;
	double XY  = 1.0 / 32.0 * (9.0 * (X*X + Y*Y) - 10.0);
	double Xm2 = 9.0 / 32.0 * (1.0 - X*X);
	double Ym2 = 9.0 / 32.0 * (1.0 - Y*Y);

	ptx  = points[0] * (9.0 / 16.0 * Xm * X * Ym - Ym * XY);
	ptx.add(points[1], (9.0 / 16.0 * Xp * X * Ym + Ym * XY));
	ptx.add(points[2], (9.0 / 16.0 * Xp * X * Yp + Yp * XY));
	ptx.add(points[3], (9.0 / 16.0 * Xm * X * Yp - Yp * XY));

	ptx.add(points[4], -(9.0/ 16.0 * X3m * X * Ym + 3.0 * Xm2 * Ym));
	ptx.add(points[5], -(9.0/ 16.0 * X3p * X * Ym - 3.0 * Xm2 * Ym));
	ptx.add(points[6], -(9.0/ 16.0 * X3p * X * Yp - 3.0 * Xm2 * Yp));
	ptx.add(points[7], -(9.0/ 16.0 * X3m * X * Yp + 3.0 * Xm2 * Yp));

	ptx.add(points[8],  -(Y3m * Ym2));
	ptx.add(points[9],   (Y3m * Ym2));
	ptx.add(points[10],  (Y3p * Ym2));
	ptx.add(points[11], -(Y3p * Ym2));

	pty  = points[0] * (9.0 / 16.0 * Xm * Y * Ym - Xm * XY);
	pty.add(points[1], (9.0 / 16.0 * Xp * Y * Ym - Xp * XY));
	pty.add(points[2], (9.0 / 16.0 * Xp * Y * Yp + Xp * XY));
	pty.add(points[3], (9.0 / 16.0 * Xm * Y * Yp + Xm * XY));

	pty.add(points[4], -(X3m * Xm2));
	pty.add(points[5], -(X3p * Xm2));
	pty.add(points[6],  (X3p * Xm2));
	pty.add(points[7],  (X3m * Xm2));

	pty.add(points[8],  -(9.0/ 16.0 * Y3m * Y * Xm + 3.0 * Ym2 * Xm));
	pty.add(points[9],  -(9.0/ 16.0 * Y3m * Y * Xp + 3.0 * Ym2 * Xp));
	pty.add(points[10], -(9.0/ 16.0 * Y3p * Y * Xp - 3.0 * Ym2 * Xp));
	pty.add(points[11], -(9.0/ 16.0 * Y3p * Y * Xm - 3.0 * Ym2 * Xm));
	// end.
}

void Quad12::test()
{
	points[0].x = 0.0; points[0].y = 0.0;
	points[1].x = 3.0; points[1].y = 0.0;
	points[2].x = 6.0; points[2].y = 3.0;
	points[3].x = 3.0; points[3].y = 3.0;
	points[4].x = 1.0; points[4].y = 0.0;
	points[5].x = 2.0; points[5].y = 0.0;
	points[6].x = 5.0; points[6].y = 3.0;
	points[7].x = 4.0; points[7].y = 3.0;
	points[8].x = 1.0; points[8].y = 1.0;
	points[9].x = 4.0; points[9].y = 1.0;
	points[10].x = 5.0; points[10].y = 2.0;
	points[11].x = 2.0; points[11].y = 2.0;
	for(int i = 0; i < 12; i++) values[i] = points[i].x * points[i].y;

	DPoint2d x(0.5, 0.5);
	//double fx = getValue(x);
	//DPoint2d pt = getPoint(x);
	DPoint2d ptx, pty;
	getDerivativePoint(x, ptx, pty);
	//LOG("Test: N(0.5, 0.5) = %0.2f\n", fx);
	//LOG("\tpt = %s\n", pt.exportString());
	//LOG("\tNx(0.5, 0.5) = %s\n", ptx.exportString());
	//LOG("\tNy(0.5, 0.5) = %s\n", pty.exportString());
	//
	x.x = 2.5;
	x.y = 0.5;
	/* double result = */ interpolateQuadValue(x);
	//LOG4CPLUS_INFO(MeshLog::logger_console, "Result of quad interpolation for pt[2.5,0.5] = %0.2f\n", result);
}

double Quad12::approximateValue(int count, DPoint2d *points, double *z, const DPoint2d & where, int fun_ct)
{
	if(fun_ct < 1) return 0.0;
	double *a = new double[fun_ct];
	svdfit(points, z, count, a, fun_ct);
	double res = 0.0;
	for(int i = 0; i < fun_ct; i++){
		res += a[i] * approximationFunction(i, where);
	}
	delete[] a;
	return res;
}

double Quad12::approximationFunction(int id, const DPoint2d &pt)
{
	switch(id){
	case 0:	return 1;
	case 1:	return pt.x;
	case 2: return pt.y;
	case 3:	return sqr(pt.x);
	case 4:	return sqr(pt.y);
	case 5: return pt.x * pt.y;
	case 6: return sqr(pt.x) * pt.x;
	case 7: return sqr(pt.y) * pt.y;
	case 8: return sqr(pt.x) * pt.y;
	case 9: return sqr(pt.y) * pt.x;
	default: assert(false); return 1.0;
	}
}

#define TOL 1.0e-5f

void Quad12::svdfit(DPoint2d *points, double *z, int count, double *a, int a_count)
{
	int j,i;
	double wmax,thresh;
	double **u = matrix(1, count, 1, a_count);
	double *w = nr_vector(1, a_count);
	double **v = matrix(1, a_count, 1, a_count);

	for(i = 1; i <= count; i++){
		for(j = 1; j <= a_count; j++){
			u[i][j] = (double)approximationFunction(j-1, points[i-1]);
		}
	}
	svdcmp(u, count, a_count, w, v);
	wmax = 0.0; 
	for(j = 1; j <= a_count; j++) if(w[j] > wmax) wmax = w[j];
	thresh = TOL*wmax;
	for(j = 1; j <= a_count; j++) if(w[j] < thresh) w[j] = 0.0;
	double *fz = nr_vector(1, count);
	double *fa = nr_vector(1, a_count);
	for(i = 1; i <= count; i++) fz[i] = (double)z[i-1];
	svbksb(u, w, v, count, a_count, fz, fa);
	for(i = 1; i <= a_count; i++) a[i-1] = (double)fa[i];


	free_vector(w, 1, a_count);
	free_matrix(v, 1, a_count, 1, a_count);
	free_matrix(u, 1, count, 1, a_count);
	free_vector(fz, 1, count);
	free_vector(fa, 1, a_count);
}
