/////////////////////////////////////////////////////////////////////////////
// Curve2dBSpline.cpp
// BSpline 2D
//	[Curve2dParametric->Curve2dBSpline]
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2004-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#include "Curve2dSegment.h"
#include "Curve2dBSpline.h"
#include "MeshData.h"
#include "DMetric2d.h"

//#include "SurfaceParametric.h"
//#include "MeshViewSet.h"

double Curve2dBSpline::BSplineMatrix[4][4] = 
	{	{-1.0/6.0,	 0.5,	-0.5,		1.0/6.0},
		{ 0.5,		-1.0,	 0.5,		0.0},
		{ -0.5,		 0.0,	 0.5,		0.0},
		{1.0/6.0,	2.0/3.0, 1.0/6.0,	0.0}
	};

double Curve2dBSpline::BDerivativeSplineMatrix[3][4] = 
	{	{-0.5,	 1.5,	-1.5,	0.5},
		{  1.0,	-2.0,	 1.0,	0.0},
		{ -0.5,	 0.0,	 0.5,	0.0}
	};

double Curve2dBSpline::BDDerivativeSplineMatrix[2][4] = 
	{	{ -1.0,	 3.0,	-3.0,	1.0},
		{  1.0,	-2.0,	 1.0,	0.0}
	};

double Curve2dBSpline::BDDDerivativeSplineMatrix[4] = 
	{ -1.0,	 3.0,	-3.0,	1.0};

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

Curve2dBSpline::Curve2dBSpline(int ct, const DPoint2d points[], bool opened) : Curve2dParametric(0.0, 1.0), m_opened(opened)
{
	for(int i=0; i<ct; i++)
		m_points.add(points[i]);
	countNodes();
}

Curve2dBSpline::Curve2dBSpline(const DataVector<DPoint2d> & points, bool opened) : Curve2dParametric(0.0, 1.0), m_opened(opened)
{
	for(size_t i=0; i< points.countInt(); i++)
		m_points.add(points[i]);
	countNodes();
}

Curve2dBSpline::Curve2dBSpline(const Curve2dBSpline& spline)
	: Curve2dParametric(0.0, 1.0)
{
	size_t i, ct = spline.m_points.countInt();
	for(i=0; i<ct; i++)
		m_points.add(spline.m_points.get(i));
	ct = spline.m_nodes.countInt();
	for(i=0; i<ct; i++)
		m_nodes.add(spline.m_nodes.get(i));
	m_opened = spline.m_opened;
}

DPoint2d Curve2dBSpline::getPoint(double t) const {
	int i = countSplineIndices(t);	// t -> ksi (0..1), i -> spline segment number
	DPoint2d eff[4];
	size_t j, k, ct = m_points.countInt()+1;
	// Obliczenie wspó³rzêdnych - najPIerw wspó³czynniki odpowiedniego wielomianu (zale¿nie od 
	//numeru i), a nastêpnie jego wartoœæ dla parametru t zredukowanego do czêœci u³amkowej.
	for(j=0; j<4; j++){
		eff[j].x = eff[j].y = 0.0;
		for(k=0; k<4; k++){
			if((i+k) > ct)
				eff[j].add(m_nodes.get(2), BSplineMatrix[j][k]);
			else
				eff[j].add(m_nodes.get(i+k), BSplineMatrix[j][k]);
		}
	}
	// Wyliczenie wartoœci wielomianu 3-go stopnia.
	for(j=1; j<4; j++){
		eff[0] *= t;
		eff[0].add(eff[j]);
	}

	return eff[0];
}

/* // Simple scanning
double Curve2dBSpline::getParameter(const DPoint2d& pt, double ts) const {
	int i, ct = m_points.countInt();
	if(ct < 2) return -1.0;

	double min_len = m_points.get(0).distance2(m_points.get(1));
	for(i = 2; i < ct; i++){
		double len = m_points.get(i-1).distance2(m_points.get(i));
		if(len < min_len) min_len = len;
	}
	double t_step = 0.01;
	double d_accuracy = min_len * 0.0001;
	for(double t = ts; t < 1.0; t += t_step){
		if(pt.distance2(getPoint(t)) <= d_accuracy) return t;
	}

	// Newton method

	return -1.0;
}
*/

/* Newton root finding
double Curve2dBSpline::getParameter(const DPoint2d& pt, double ts) const {

	double	F_ERR		= mesh_data.relative_small_number;	// Maksymalny b³¹d wyniku
	double	T_ERR		= 1e-5 * mesh_data.relative_small_number;	// Maksymalny b³¹d parametru
	int		MAX_STEPS	= 30;	// Maksymalna iloœæ kroków

	double t = ts;
	for(int i=0; i < MAX_STEPS; i++){
		const DPoint2d ft = getPoint(t);
		const DPoint2d fd = ft - pt;
		if(abs(fd.x) + abs(fd.y) < F_ERR) return t;
		DPoint2d dft = getDerivative(t);
		double dtx = (abs(dft.x) < mesh_data.relative_small_number) ? 0.0 : (fd.x / dft.x);
		double dty = (abs(dft.y) < mesh_data.relative_small_number) ? 0.0 : (fd.y / dft.y);
		t -= (i%2) ? dtx : dty;
		if(abs(dtx) + abs(dty) < T_ERR) return t;
	}

	return t;
}
*/

DVector2d Curve2dBSpline::getDerivative(double t) const {
	size_t j, k, i = countSplineIndices(t);	// t -> ksi (0..1), i -> spline segment number
	DPoint2d eff[3];
	// Pochodna wzgl. [x,y]	-> Podobnie jak wyliczanie wspó³rzêdnych punktu krzywej sklejanej
	//	najPIerw wyliczane s¹ wspó³czynniki odpowiedniego wielomianu, a nastêpnie jego
	//	wartoœæ dla podanego parametru. (Ró¿ni¹ siê zastosowane macierze funkcji kszta³tu)
	size_t ct = m_points.countInt()+1;
	for(j=0; j<3; j++){
		eff[j].x = eff[j].y = 0.0;
		for(k=0; k<4; k++){
			if((i+k) > ct)
				eff[j].add(m_nodes.get(2), BDerivativeSplineMatrix[j][k]);
			else
				eff[j].add(m_nodes.get(i+k), BDerivativeSplineMatrix[j][k]);
		}
	}
	// Wyliczenie wartoœci wielomianu
	for(j=1; j<3; j++){
		eff[0] *= t;
		eff[0].add(eff[j]);
	}

	return eff[0].fixedVector();
}

DVector2d Curve2dBSpline::getSecondDerivative(double t) const {
	size_t j, k, i = countSplineIndices(t);	// t -> ksi (0..1), i -> spline segment number
	DPoint2d eff[2];
	// Druga pochodna wzgl. [x,y]	-> Podobnie jak wyliczanie wspó³rzêdnych punktu krzywej sklejanej
	//	najPIerw wyliczane s¹ wspó³czynniki odpowiedniego wielomianu, a nastêpnie jego
	//	wartoœæ dla podanego parametru. (Ró¿ni¹ siê zastosowane macierze funkcji kszta³tu)
	size_t ct = m_points.countInt()+1;
	for(j=0; j<2; j++){
		eff[j].x = eff[j].y = 0.0;
		for(k=0; k<4; k++){
			if((i+k) > ct)
				eff[j].add(m_nodes.get(2), BDDerivativeSplineMatrix[j][k]);
			else
				eff[j].add(m_nodes.get(i+k), BDDerivativeSplineMatrix[j][k]);
		}
	}

	eff[0] *= t;
	eff[0].add(eff[1]);
	return eff[0].fixedVector();
}

DVector2d Curve2dBSpline::getThirdDerivative(double t) const {
	size_t i = countSplineIndices(t);	// t -> ksi (0..1), i -> spline segment number
	DPoint2d eff(0.0, 0.0);
	// Druga pochodna wzgl. [x,y]	-> Podobnie jak wyliczanie wspó³rzêdnych punktu krzywej sklejanej
	//	najPIerw wyliczane s¹ wspó³czynniki odpowiedniego wielomianu, a nastêpnie jego
	//	wartoœæ dla podanego parametru. (Ró¿ni¹ siê zastosowane macierze funkcji kszta³tu)
	size_t ct = m_points.countInt()+1;
	for(int k=0; k<4; k++)
		if((i+k) > ct) eff.add(m_nodes.get(2), BDDDerivativeSplineMatrix[k]);
		else eff.add(m_nodes.get(i+k), BDDDerivativeSplineMatrix[k]);

	return eff.fixedVector();
}

ostream& operator<<(ostream& os, const Curve2dBSpline *spline){
	size_t ct = spline->m_points.countInt();
	os << " OP=" << (spline->m_opened?1:0);
	os << " PT=";
	for(size_t i = 0; i < ct; i++)
		os << spline->m_points.get(i) << " ";
	return os;
}

istream& operator>>(istream& is, Curve2dBSpline *spline){
	int op;
	char z;
	is >> ws >> z >> z >> z >> op;	//	"OP="
	spline->m_opened = (op == 1);
	is >> ws >> z >> z >> z;	//	"PT="
	spline->m_points.clear();
	do{
		DPoint2d pt;
		is >> ws >> pt;
		if(is)	spline->m_points.add(pt);
	}while(is);
	spline->countNodes();
	return is;
}

bool Curve2dBSpline::countNodes()
{
	size_t n = m_points.countInt();
	DPoint2d* nodes = new DPoint2d[n+2];
	size_t i;

	double *A = new double[n];
	DPoint2d *B = new DPoint2d[n];

	if(m_opened){
		// Wersja dla krzywej sklejanej otwartej
		nodes[0] = m_points.get(0);
		nodes[n+1] = m_points.get(n-1);

		// Wspó³rzêdne x,y - uk³ad równañ Ax=b, gdzie macierz A jest trójdiagonalna,
		//	poni¿ej efektywny algorytm rozwi¹zywania tego specyficznego zadania.
		A[n-1] = -0.25;
		B[n-1].set(m_points.get(n-1), 6.0);
		B[n-1].add(nodes[n+1], -1.0);
		B[n-1] *= 0.25;
		for(i=n-1; i>1; i--){
			if(A[i] == -4.0) {
				delete[] A;
				delete[] B;
				delete[] nodes;
				return false;
			}
			A[i-1] = -1/(A[i]+4);
			B[i-1].set(m_points.get(i-1), 6.0);
			B[i-1].add(B[i], -1.0);
			B[i-1] /= (A[i]+4);
		}

		nodes[1].set(m_points.get(0), 6);
		nodes[1].add(nodes[0], -1.0);
		nodes[1].add(B[1], -1.0);
		nodes[1] /= (A[1]+4);
		for(i=1; i<n; i++){
			nodes[i+1].set(nodes[i], A[i]);
			nodes[i+1].add(B[i]);
		}
	}else{
		// Wersja dla krzywej sklejanej zamkniêtej, utrudnieniem w stosunku do poprzedniej jest
		//	fakt, ¿e macierz A w uk³adzie równañ jest "prawie" diagonalna, tzn. zawiera dodatkowo
		//	po jednym elemencie w rogach macierzy. Zastosowany poni¿ej algorytm rozwi¹zuje ten 
		//	problem poprzez wyliczenie dodatkowej poprawki i sk³ada siê z dwóch kroków, przy
		//	czym w ka¿dym stosowany jest algorytm podobny do zastosowanego powy¿ej.

		// Wspó³rzêdne x
		A[n-1] = -1.0/4.25;
		B[n-1].set(m_points.get(n-1), 6.0/4.25);
		for(i=n-1; i>1; i--){
			if(A[i] == -4.0) {
				delete[] A;
				delete[] B;
				delete[] nodes;
				return false;
			}
			A[i-1] = -1.0/(A[i]+4.0);
			B[i-1].set(m_points.get(i-1), 6.0);
			B[i-1].add(B[i], -1.0);
			B[i-1] /= (A[i]+4.0);
		}

		nodes[1].set(m_points.get(0), 6.0);
		nodes[1].add(B[1], -1.0);
		nodes[1] /= (A[1]+8.0);
		for(i=1; i<n; i++){
			nodes[i+1].set(nodes[i], A[i]);
			nodes[i+1].add(B[i]);
		}

		// Czêœæ II

		double *z = new double[n+1];
		double *C = new double[n];
		C[n-1] = 1.0/4.25;
		for(i=n-1; i>1; i--)
			C[i-1] = -C[i]/(A[i]+4.0);

		z[1] = (-4.0 - C[1])/(A[1]+8.0);
		for(i=1; i<n; i++)
			z[i+1] = A[i]*z[i] + C[i];

		DVector2d fact = (nodes[1] - nodes[n]*0.25) / (1.0+z[1]-z[n]/4.0); //Równanie v . x=(1 + v . z).
		for (i = 1; i <= n; i++) nodes[i] -= fact*z[i]; //Koñcowa postaæ wspó³rzêdnych.

		delete[] z;
		delete[] C;

		nodes[0] = nodes[n];
		nodes[n+1] = nodes[1];
	}

	m_nodes.clear();
	for(i=0; i < n+2; i++)
		m_nodes.add(nodes[i]);

	delete[] A;
	delete[] B;
	delete[] nodes;

	return true;
}

int Curve2dBSpline::countSplineIndices(double &t) const
{
	// t should be [0,1]
	if(m_opened){
		if(t < 0.0) t = 0.0;
		else if(t > 1.0) t = 1.0;
	}else{
		t = fmod(t, 1.0);
		if(t < 0.0) t += 1.0;
	}
	// [0,1] -> [0,max_ct]
	t *= getMaxParameter();
	// i - numer odcinka krzywej sklejanej
	double int_part;
	t = modf(t, &int_part);
	return (int)int_part;
}

double Curve2dBSpline::getParameterInRange( const DPoint2d& pt, double ts, double t_min, double t_max ) const {
	size_t n = m_points.countInt()-1;
	double dt = 1.0 / n;
	double t0 = 0.0;
	double min_dist2 = -1.0;
	double min_t;
	for(size_t i = 0; i < n; i++, t0+=dt){
		if(t0 > t_max || (t0+dt) < t_min) continue;
		// try this subsegment
		DMatrix2d ms;
		DVector2d vt(m_points[i], m_points[i+1]);
		ms.setColumn(0, vt);
		ms.setColumn(1, DVector2d(vt.y, -vt.x));
		DVector2d res;
		bool ok = ms.solve(pt - m_points[i], res);
		assert(ok);
		if(!ok) return Curve2dParametric::getParameterInRange( pt, ts, t_min, t_max );
		double dist2;
		double t = res.x;
		if(t < 0.0) t = 0.0;
		else if(t > 1.0) t = 1.0;
		double tx = t0 + t * dt;
		tx = std::max(t_min, std::min(t_max, tx));
		dist2 = pt.distance2(getPoint(tx));
		if(min_dist2 < 0.0 || dist2 < min_dist2){
			min_dist2 = dist2;
			min_t = tx;
		}
	}
	if(min_dist2 >= 0.0) return min_t;
	else return Curve2dParametric::getParameterInRange( pt, ts, t_min, t_max );
}

/// fit to point-cloud using direct fit, return max-dist
double Curve2dBSpline::fitToPoints(const DataVector<DPoint2d> & points, 
	const DPoint2d& ptA, const DPoint2d& ptB, const DMetric2d& dm)
{
	m_opened = true;
	m_nodes.clear();
	m_points.clear();
	m_points.add(ptA);
	m_points.add(ptB);
	if(fitToPointsAdaptive(points, 0, dm) == 0)
		countNodes();
	// calculate global error (max dist)
	return maxDistanceForSpline(points, dm);
}

/// fit to point-cloud using direct fit, return max-dist
void Curve2dBSpline::fitThroughPoints(const DataVector<DPoint2d> & points, 
	const DPoint2d& ptA, const DPoint2d& ptB, const DMetric2d& dm)
{
	m_opened = true;
	m_nodes.clear();
	m_points.clear();
	m_points.add(ptA);
	m_points.add(ptB);
	fitThroughPointsAdaptive(points, 0, dm);
	countNodes();
}

/// adaptive subroutine - fit to point-cloud using direct fit, return number of inserted points
int Curve2dBSpline::fitToPointsAdaptive(const DataVector<DPoint2d> & points, int ind, const DMetric2d& dm)
{
	size_t pct = points.countInt();
	if(pct < 5) return 0;
	assert(ind+1 < m_points.countInt());
	const DPoint2d& ptA = m_points[ind];
	const DPoint2d& ptB = m_points[ind+1];
	const DVector2d vt = ptB - ptA;
	const DPoint2d middle(ptA, ptB, 0.5);
	double c = vt.x * middle.x + vt.y * middle.y;
	// count distances from the ortoghonal line, crossing middle-point:
	// ... and select three points closest to the middle
	DataVector<double> distances(pct);
	DataVector<size_t> closest(3);
	for(size_t k = 0; k < pct; k++){
		const DPoint2d& pt = points[k];
		double dist = vt.x * pt.x + vt.y * pt.y - c;
		distances.add(dist);
		if(closest.countInt() < 3) closest.add(k);
		else{
			// check and replace if necessary...
			size_t ki = k;
			for(int i = 0; i < 3; i++){
				if(abs(distances[ki]) < abs(distances[closest[i]])){
					size_t temp = closest[i];
					closest[i] = ki;
					ki = temp;
				}
			}
		}
	}
	DPoint2d npoint = DPoint2d::zero;
	for(int i = 0; i < 3; i++){ // average of "most-middle" real-points
		npoint.add(points[closest[i]], 1.0/3.0);
	}
	// insert and recalculate nodes
	m_points.insertAt(ind+1, npoint);
	countNodes();

	//if(true){
	//	MeshViewSet* set = new MeshViewSet;
	//	for(int i = 0; i < m_points.countInt(); i++){
	//		if(i == ind+1)
	//			set->addPoint(DPoint3d(m_points[i], 0.0), 0, i);
	//		else 
	//			set->addPoint(DPoint3d(m_points[i], 0.0), 0);
	//	}
	//	for(int i = 0; i <= m_points.countInt(); i++)
	//		set->addEdge(DPoint3d(m_nodes[i], 0.0), DPoint3d(m_nodes[i+1], 0.0));
	//	for(int i = 0; i < points.countInt(); i++){
	//		if(i == closest[0] || i == closest[1] || i == closest[2])
	//			set->addPoint(DPoint3d(points[i], 0.0), 2);
	//		else
	//			set->addPoint(DPoint3d(points[i], 0.0), 1);
	//	}
	//	DataVector<double> polyline;
	//	getPolyLine(0.0, 1.0, polyline);
	//	for(int k = 0; k < polyline.countInt()-1; k++){
	//		set->addEdge(DPoint3d(getPoint(polyline[k]), 0.0), 
	//			DPoint3d(getPoint(polyline[k+1]), 0.0), 2);
	//	}
	//	SHOW_MESH("BSpline fit", set);
	//}

	// check maximum distance (left and right), excluding the three points used for calculating new node
	DataVector<DPoint2d> points_left(pct);
	DataVector<DPoint2d> points_right(pct);
	for(size_t k = 0; k < pct; k++){
		if(k == closest[0] || k == closest[1] || k == closest[2])
			continue;
		if(distances[k] < 0.0){
			points_left.add(points[k]);
		}else{
			points_right.add(points[k]);
		}
	}

//	double max_dist_left =  maxDistanceForSubSegment(points_left, ind, dm);
//	double max_dist_right = maxDistanceForSubSegment(points_right, ind+1, dm);

	// recurence left ...
	int additional_left = 0;
//	if(max_dist_left > 1.0){ // tolerance is stored within the metric 
		// run
		additional_left = fitToPointsAdaptive(points_left, ind, dm);
//	}
	// ... and right
	int additional_right = 0;
//	if(max_dist_right > 1.0){
		// run
		additional_right = fitToPointsAdaptive(points_right, ind+additional_left+1, dm);
//	}
	return 1 + additional_left + additional_right;
}

/// adaptive subroutine - fit to point-cloud using direct fit, return number of inserted points
int Curve2dBSpline::fitThroughPointsAdaptive(const DataVector<DPoint2d> & points, int ind, const DMetric2d& dm)
{
	size_t pct = points.countInt();
	if(pct < 3) return 0;
	assert(ind+1 < m_points.countInt());
	const DPoint2d& ptA = m_points[ind];
	const DPoint2d& ptB = m_points[ind+1];
	const DVector2d vt = ptB - ptA;
	const DPoint2d middle(ptA, ptB, 0.5);
	double c = vt.x * middle.x + vt.y * middle.y;
	// count distances from the ortoghonal line, crossing middle-point:
	// ... and select the point closest to the middle
	DataVector<double> distances(pct);
	size_t closest = std::string::npos;
	for(size_t k = 0; k < pct; k++){
		const DPoint2d& pt = points[k];
		double dist = vt.x * pt.x + vt.y * pt.y - c;
		distances.add(dist);
		if(closest == std::string::npos || (abs(distances[closest]) > abs(dist))) 
			closest = k;
	}
	// insert and recalculate nodes
	m_points.insertAt(ind+1, points[closest]);
	countNodes();

	// check maximum distance (left and right), excluding the closest point
	DataVector<DPoint2d> points_left(pct);
	DataVector<DPoint2d> points_right(pct);
	for(size_t k = 0; k < pct; k++){
		if(k == closest) continue;
		if(distances[k] < 0.0){
			points_left.add(points[k]);
		}else{
			points_right.add(points[k]);
		}
	}

	// recurence left ...
	int additional_left = fitThroughPointsAdaptive(points_left, ind, dm);
	// ... and right
	int additional_right = fitThroughPointsAdaptive(points_right, ind+additional_left+1, dm);

	return 1 + additional_left + additional_right;
}

/// max distance
double Curve2dBSpline::maxDistanceForSubSegment(const DataVector<DPoint2d>& points, int ind, const DMetric2d& dm) const
{
	const DPoint2d& ptA = m_points[ind];
	const DVector2d vt = m_points[ind+1] - ptA; 
	const DVector2d vn(vt.y, -vt.x);
	// matrix
	DMatrix2d m;
	m.setColumn(0, vt);
	m.setColumn(1, vn);
	double max_dist2 = 0.0;
	DVector2d res;
	for(size_t k = 0; k < points.countInt(); k++){
		const DPoint2d& ptx = points[k];
		bool ok = m.solve(ptx - ptA, res);
		assert(ok);
		if(ok){
			double t = std::max(0.0, std::min(1.0, res.x));
			double dist2 = dm.transformPStoMS(ptx - getPoint((ind+t) / getMaxParameter())).length2();
			if(dist2 > max_dist2) max_dist2 = dist2;
		}
	}
	return sqrt(max_dist2);
}

/// max distance
double Curve2dBSpline::maxDistanceForSpline(const DataVector<DPoint2d>& points, const DMetric2d& dm) const
{
	size_t n = m_points.countInt()-1;
	DataVector<DMatrix2d> ms(n, DMatrix2d());
	for(size_t i = 0; i < n; i++){
		DVector2d vt(m_points[i], m_points[i+1]);
		ms[i].setColumn(0, vt);
		ms[i].setColumn(1, DVector2d(vt.y, -vt.x));
	}

	double max_dist2 = 0.0;
	DVector2d res;
	for(size_t k = 0; k < points.countInt(); k++){
		const DPoint2d& ptx = points[k];
		double min_dist2 = -1.0;
		for(size_t i = 0; i < n; i++){
			bool ok = ms[i].solve(ptx - m_points[i], res);
			assert(ok);
			if(!ok) continue;
			double dist2;
			double t = res.x;
			if(t < 0.0){
				dist2 = dm.transformPStoMS(ptx - m_points[i]).length2();
			}else if(t > 1.0){
				dist2 = dm.transformPStoMS(ptx - m_points[i+1]).length2();
			}else{
				dist2 = dm.transformPStoMS(ptx - getPoint((i+t) / n)).length2();
			}
			if(min_dist2 < 0.0 || dist2 < min_dist2) min_dist2 = dist2;
		}
		if(min_dist2 > max_dist2) max_dist2 = min_dist2;
	}
	return sqrt(max_dist2);
}

/// Store XML description to stream
ostream& Curve2dBSpline::storeXML(ostream& os, const string& prefix) const
{
	os << prefix << "<bspline-curve>\n";
	Curve2dParametric::storeXML(os, prefix + "\t");
	if(m_opened)
		os << prefix << "\t<points>\n";
	else
		os << prefix << "\t<points opened=\"false\">\n";
	for(size_t i = 0; i < m_points.countInt(); i++){
		os << prefix << "\t\t<point> " << m_points[i] << " </point>\n";
	}
	return os << prefix << "\t</points>\n" << prefix << "</bspline-curve>\n";
}
