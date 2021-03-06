/////////////////////////////////////////////////////////////////////////////
// Curve3dBSpline.cpp
// BSpline 3D
//	[Curve3dParametric->Curve3dBSpline]
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2013-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#include "Curve3dSegment.h"
#include "Curve3dBSpline.h"
#include "MeshData.h"
#include "DMetric3d.h"
#include "Metric3dContext.h"
#include "DLeastSquaresFitting.h"
#include "DLinearQuadric.h"

#include "Curve2dBSpline.h"

//#include "Surface3dParametric.h"
//#include "MeshViewSet.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

Curve3dBSpline::Curve3dBSpline(int ct, const DPoint3d points[], bool opened) 
	: Curve3dParametric(0.0, 1.0), m_opened(opened)
{
	for(int i=0; i<ct; i++)
		m_points.add(points[i]);
	countNodes();
}

Curve3dBSpline::Curve3dBSpline(const DataVector<DPoint3d> & points, bool opened) 
	: Curve3dParametric(0.0, 1.0), m_opened(opened)
{
	for(int i=0; i< points.countInt(); i++)
		m_points.add(points[i]);
	countNodes();
}

Curve3dBSpline::Curve3dBSpline(const Curve3dBSpline& spline) : m_opened( spline.m_opened )
{
	size_t i, ct = spline.m_points.countInt();
	for(i=0; i<ct; i++)
		m_points.add(spline.m_points.get(i));
	ct = spline.m_nodes.countInt();
	for(i=0; i<ct; i++)
		m_nodes.add(spline.m_nodes.get(i));
}

DPoint3d Curve3dBSpline::getPoint(double t) const {
	int i = countSplineIndices(t);	// t -> ksi (0..1), i -> spline segment number
	DPoint3d eff[4];
	size_t j, k, ct = m_points.countInt()+1;
	// Obliczenie wsp??rz?dnych - najPIerw wsp??czynniki odpowiedniego wielomianu (zale?nie od 
	//numeru i), a nast?pnie jego warto?? dla parametru t zredukowanego do cz??ci u?amkowej.
	for(j=0; j<4; j++){
		for(k=0; k<4; k++){
			if((i+k) > ct)
				eff[j].add(m_nodes.get(2), Curve2dBSpline::BSplineMatrix[j][k]);
			else
				eff[j].add(m_nodes.get(i+k), Curve2dBSpline::BSplineMatrix[j][k]);
		}
	}
	// Wyliczenie warto?ci wielomianu 3-go stopnia.
	for(j=1; j<4; j++){
		eff[0] *= t;
		eff[0].add(eff[j]);
	}

	return eff[0];
}

/* // Simple scanning
double Curve3dBSpline::getParameter(const DPoint3d& pt, double ts) const {
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
double Curve3dBSpline::getParameter(const DPoint3d& pt, double ts) const {

	double	F_ERR		= mesh_data.relative_small_number;	// Maksymalny b??d wyniku
	double	T_ERR		= 1e-5 * mesh_data.relative_small_number;	// Maksymalny b??d parametru
	int		MAX_STEPS	= 30;	// Maksymalna ilo?? krok?w

	double t = ts;
	for(int i=0; i < MAX_STEPS; i++){
		const DPoint3d ft = getPoint(t);
		const DPoint3d fd = ft - pt;
		if(abs(fd.x) + abs(fd.y) < F_ERR) return t;
		DPoint3d dft = getDerivative(t);
		double dtx = (abs(dft.x) < mesh_data.relative_small_number) ? 0.0 : (fd.x / dft.x);
		double dty = (abs(dft.y) < mesh_data.relative_small_number) ? 0.0 : (fd.y / dft.y);
		t -= (i%2) ? dtx : dty;
		if(abs(dtx) + abs(dty) < T_ERR) return t;
	}

	return t;
}
*/

DVector3d Curve3dBSpline::getDerivative(double t) const {
	int j, k, i = countSplineIndices(t);	// t -> ksi (0..1), i -> spline segment number
	DPoint3d eff[3];
	// Pochodna wzgl. [x,y]	-> Podobnie jak wyliczanie wsp??rz?dnych punktu krzywej sklejanej
	//	najPIerw wyliczane s? wsp??czynniki odpowiedniego wielomianu, a nast?pnie jego
	//	warto?? dla podanego parametru. (R??ni? si? zastosowane macierze funkcji kszta?tu)
	size_t ct = m_points.countInt()+1;
	for(j=0; j<3; j++){
		for(k=0; k<4; k++){
			if((i+k) > ct)
				eff[j].add(m_nodes.get(2), Curve2dBSpline::BDerivativeSplineMatrix[j][k]);
			else
				eff[j].add(m_nodes.get(i+k), Curve2dBSpline::BDerivativeSplineMatrix[j][k]);
		}
	}
	// Wyliczenie warto?ci wielomianu
	for(j=1; j<3; j++){
		eff[0] *= t;
		eff[0].add(eff[j]);
	}

	return eff[0].fixedVector();
}

DVector3d Curve3dBSpline::getSecondDerivative(double t) const {
	int j, k, i = countSplineIndices(t);	// t -> ksi (0..1), i -> spline segment number
	DPoint3d eff[2];
	// Druga pochodna wzgl. [x,y]	-> Podobnie jak wyliczanie wsp??rz?dnych punktu krzywej sklejanej
	//	najPIerw wyliczane s? wsp??czynniki odpowiedniego wielomianu, a nast?pnie jego
	//	warto?? dla podanego parametru. (R??ni? si? zastosowane macierze funkcji kszta?tu)
	size_t ct = m_points.countInt()+1;
	for(j=0; j<2; j++){
		for(k=0; k<4; k++){
			if((i+k) > ct)
				eff[j].add(m_nodes.get(2), Curve2dBSpline::BDDerivativeSplineMatrix[j][k]);
			else
				eff[j].add(m_nodes.get(i+k), Curve2dBSpline::BDDerivativeSplineMatrix[j][k]);
		}
	}

	eff[0] *= t;
	eff[0].add(eff[1]);
	return eff[0].fixedVector();
}

DVector3d Curve3dBSpline::getThirdDerivative(double t) const {
	int i = countSplineIndices(t);	// t -> ksi (0..1), i -> spline segment number
	DPoint3d eff;
	// Druga pochodna wzgl. [x,y]	-> Podobnie jak wyliczanie wsp??rz?dnych punktu krzywej sklejanej
	//	najPIerw wyliczane s? wsp??czynniki odpowiedniego wielomianu, a nast?pnie jego
	//	warto?? dla podanego parametru. (R??ni? si? zastosowane macierze funkcji kszta?tu)
	size_t ct = m_points.countInt()+1;
	for(int k=0; k<4; k++)
		if((i+k) > ct) eff.add(m_nodes.get(2), Curve2dBSpline::BDDDerivativeSplineMatrix[k]);
		else eff.add(m_nodes.get(i+k), Curve2dBSpline::BDDDerivativeSplineMatrix[k]);

	return eff.fixedVector();
}

ostream& operator<<(ostream& os, const Curve3dBSpline *spline){
	size_t ct = spline->m_points.countInt();
	os << " OP=" << (spline->m_opened?1:0);
	os << " PT=";
	for(int i = 0; i < ct; i++)
		os << spline->m_points.get(i) << " ";
	return os;
}

istream& operator>>(istream& is, Curve3dBSpline *spline){
	int op;
	char z;
	is >> ws >> z >> z >> z >> op;	//	"OP="
	spline->m_opened = (op == 1);
	is >> ws >> z >> z >> z;	//	"PT="
	spline->m_points.clear();
	do{
		DPoint3d pt;
		is >> ws >> pt;
		if(is)	spline->m_points.add(pt);
	}while(is);
	spline->countNodes();
	return is;
}

bool Curve3dBSpline::countNodes()
{
	size_t n = m_points.countInt();
	DPoint3d* nodes = new DPoint3d[n+2];
	size_t i;

	double *A = new double[n];
	DPoint3d *B = new DPoint3d[n];

	if(m_opened){
		// Wersja dla krzywej sklejanej otwartej
		nodes[0] = m_points.get(0);
		nodes[n+1] = m_points.get(n-1);

		// Wsp??rz?dne x,y - uk?ad r?wna? Ax=b, gdzie macierz A jest tr?jdiagonalna,
		//	poni?ej efektywny algorytm rozwi?zywania tego specyficznego zadania.
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
		// Wersja dla krzywej sklejanej zamkni?tej, utrudnieniem w stosunku do poprzedniej jest
		//	fakt, ?e macierz A w uk?adzie r?wna? jest "prawie" diagonalna, tzn. zawiera dodatkowo
		//	po jednym elemencie w rogach macierzy. Zastosowany poni?ej algorytm rozwi?zuje ten 
		//	problem poprzez wyliczenie dodatkowej poprawki i sk?ada si? z dw?ch krok?w, przy
		//	czym w ka?dym stosowany jest algorytm podobny do zastosowanego powy?ej.

		// Wsp??rz?dne x
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

		// Cz??? II

		double *z = new double[n+1];
		double *C = new double[n];
		C[n-1] = 1.0/4.25;
		for(i=n-1; i>1; i--)
			C[i-1] = -C[i]/(A[i]+4.0);

		z[1] = (-4.0 - C[1])/(A[1]+8.0);
		for(i=1; i<n; i++)
			z[i+1] = A[i]*z[i] + C[i];

		DVector3d fact = (nodes[1] - nodes[n]*0.25) / (1.0+z[1]-z[n]/4.0); //R?wnanie v . x=(1 + v . z).
		for (i = 1; i <= n; i++) nodes[i] -= fact*z[i]; //Ko?cowa posta? wsp??rz?dnych.

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

int Curve3dBSpline::countSplineIndices(double &t) const
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

double Curve3dBSpline::getParameter( const DPoint3d& pt, double ts ) const {
	return getParameterInRange( pt, ts, 0.0, 1.0 );
}

double Curve3dBSpline::getParameterInRange(const DPoint3d& pt, double ts, double t_min, double t_max) const {

	if(t_min < 0.0) t_min = 0.0;
	if(t_max > 1.0) t_max = 1.0;

	size_t n = m_points.countInt()-1;
	double dt = 1.0 / n;
	double t0 = 0.0;
	double min_dist2 = -1.0;
	double min_t;
	for(size_t i = 0; i < n; i++, t0+=dt){
		if(t0 > t_max || (t0+dt) < t_min) continue;
		// try this subsegment
		DMatrix3d ms;
		DVector3d vt(m_points[i], m_points[i+1]);
		DVector3d vt1, vt2;
		ms.setColumn(0, vt);
		vt.orthogonalVectors(vt1, vt2);
		ms.setColumn(1, vt1);
		ms.setColumn(2, vt2);
		DVector3d res;
		bool ok = ms.solve(pt - m_points[i], res);
		assert(ok);
		if(!ok) return Curve3dParametric::getParameterInRange(pt, ts, t_min, t_max);
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
	else return Curve3dParametric::getParameterInRange(pt, ts, t_min, t_max);
}


Curve3dPtr Curve3dBSpline::fitToPointSequence(const DataVector<DPoint3d> & points, 
		Metric3dContext & mc, const double& tolerance )
{
	int pct = points.countInt();
	assert( pct > 1 );

	if( pct < 5 ) 
		return std::make_shared<Curve3dBSpline>(points);

	DataVector<DPoint3d> init_points(2);
	init_points.add( points.first() );
	init_points.add( points.last() );
	auto bspline = std::make_shared<Curve3dBSpline>(init_points);
	bspline->fitToPointSequenceAdaptive(points, mc, tolerance);
	return bspline;	
}

void Curve3dBSpline::fitToPointSequenceAdaptive(const DataVector<DPoint3d> & points, 
		Metric3dContext & mc, const double& tolerance )
{
	static const double MIN_MLEN = 2.0;
	static const int MIN_RANGE = 4;

	size_t pct = points.countInt();
	assert( pct > MIN_RANGE );

	DataVector<DPoint3d> mid_points(5);
	DLinearQuadric3d lquadric;
	DataVector<double> params( pct, 0.0 );
	DataVector<bool> bpoints(pct, false);
	bpoints.first() = true;
	bpoints.last() = true;
	DataVector<DPoint3d> spoints( pct );
	for(size_t i = 0; i < pct; i++)
		spoints.add( points[i] );
	DataVector<size_t> parts[2];
	int p_from = 0;
	int p_to = 1;
	parts[p_from].add( 0 * pct + (pct-1) );
	double tol2 = tolerance * tolerance;
	double iter_max_dist2 = 0.0;

	do{
		parts[p_to].clear();
		for(size_t i = 0; i < parts[p_from].countInt(); i++ ) {
			size_t tmp = parts[p_from][i];
			size_t a = tmp / pct;
			size_t b = tmp % pct;
			// check size
			assert( b-a >= MIN_RANGE );

			// check metric length
			double mlen = 0.0;
			for(size_t j = a; (j < b) && (mlen < MIN_MLEN); j++ ) {
				mc.countMetricAtPoint( DPoint3d( points[j], points[j+1], 0.5 ), true );
				mlen += mc.transformRStoMS( points[j+1]-points[j] ).length();
			}
			if( mlen < MIN_MLEN ) continue; // too short to be further adapted

			// check max dist
			double max_mdist2 = 0.0;
			for(size_t j = a; j <= b; j++ ) {
				mc.countMetricAtPoint( points[j], true );
				params[j] = getParameter( points[j], (j==0)?0.0:params[j-1]);
				double mdist2 = mc.transformRStoMS( points[j] - getPoint( params[j] ) ).length2();
				if( mdist2 > max_mdist2 ) max_mdist2 = mdist2;
			}
			if( max_mdist2 < tol2) continue; // good enough

			// if still needs to be split:
			// - select middle 5 points
			size_t mid = (a+b)/2;
			mid_points.clear();
			size_t mid_start = mid-2; assert( mid_start >= 0 );
			for(int i = 0; i < 5; i++) mid_points.add( points[mid_start+i] );
			// - fit plane & quadric-curve-3d
			if( DLeastSquaresFitting::fitLinearQuadric( mid_points, lquadric, false ) >= 0.0 ) {
				// - get middle point from quadric-curve-3d
				spoints[mid] = lquadric.getPoint( lquadric.getParameter( spoints[mid] ) );
				// spoints[mid] = qcurve->getPoint(tmid);
			}
			// - split
			bpoints[mid] = true;
			if( mid - a >= MIN_RANGE )
				parts[p_to].add( a * pct + mid);
			if( b - mid >= MIN_RANGE )
				parts[p_to].add( mid * pct + b);
		}
		// update bspline
		m_points.clear();
		for(size_t i = 0; i < pct; i++)
			if( bpoints[i] )
				m_points.add( spoints[i] );
		countNodes();

		//LOG4CPLUS_INFO(MeshLog::logger_mesh, "========== BSPLINE-3D FIT [" << m_points.countInt() << "]============");
		//LOG4CPLUS_INFO(MeshLog::logger_mesh, " tolerance = " << tolerance);
		iter_max_dist2 = 0.0;
		for(size_t j = 0; j < pct; j++ ) {
			mc.countMetricAtPoint( points[j], true );
			params[j] = getParameter( points[j], (j==0)?0.0:params[j-1]);
			double mdist2 = mc.transformRStoMS( points[j] - getPoint( params[j] ) ).length2();
			if(mdist2 > iter_max_dist2) iter_max_dist2 = mdist2;
			//LOG4CPLUS_INFO(MeshLog::logger_mesh, j << ( bpoints[j] ? "*" : "-" ) << "\t" << sqrt(mdist2));
		}
		//LOG4CPLUS_INFO(MeshLog::logger_mesh, "=========================================");

		// prepare next loop
		p_from = 1-p_from;
		p_to = 1-p_to;
	}while( parts[p_from].notEmpty() );

	if( iter_max_dist2 > tol2 ) { 
		//LOG4CPLUS_INFO(MeshLog::logger_mesh, "*** BSPLINE-3D *2 FIT ***");
		bool any_change = false;
		do{
			any_change = false;
			for(size_t i = 0; i < pct; i++ ) {
				if(bpoints[i]) continue;
				mc.countMetricAtPoint( points[i], true );
				params[i] = getParameter( points[i], (i==0)?0.0:params[i-1]);
				double mdist2 = mc.transformRStoMS( points[i] - getPoint( params[i] ) ).length2();
				if( mdist2 > tol2){
					mid_points.clear();
					size_t j_start = std::max((size_t)0, i);
					for(size_t j = 0; (j < 5) && (j_start+j < pct); j++) mid_points.add( points[j_start+j] );
					if( DLeastSquaresFitting::fitLinearQuadric( mid_points, lquadric, false ) >= 0.0 ) {
						spoints[i] = lquadric.getPoint( lquadric.getParameter( spoints[i] ) );
					}
					// - split
					bpoints[i] = any_change = true;
				}
			}
			if(any_change){
				// update bspline
				m_points.clear();
				for(size_t i = 0; i < pct; i++)
					if( bpoints[i] )
						m_points.add( spoints[i] );
				countNodes();
			}
		}while(any_change);
		
		//LOG4CPLUS_INFO(MeshLog::logger_mesh, "========== BSPLINE-3D FIT *2 [" << m_points.countInt() << "]============");
		//LOG4CPLUS_INFO(MeshLog::logger_mesh, " tolerance = " << tolerance);
		iter_max_dist2 = 0.0;
		for(size_t j = 0; j < pct; j++ ) {
			mc.countMetricAtPoint( points[j], true );
			params[j] = getParameter( points[j], (j==0)?0.0:params[j-1]);
			double mdist2 = mc.transformRStoMS( points[j] - getPoint( params[j] ) ).length2();
			if(mdist2 > iter_max_dist2) iter_max_dist2 = mdist2;
			//LOG4CPLUS_INFO(MeshLog::logger_mesh, j << ( bpoints[j] ? "*" : "-" ) << "\t" << sqrt(mdist2));
		}
		//LOG4CPLUS_INFO(MeshLog::logger_mesh, "=========================================");
	}
}

/// fit to point-cloud using direct fit
Curve3dPtr Curve3dBSpline::fitToPoints(const DataVector<DPoint3d> & points, 
	const DPoint3d& ptA, const DPoint3d& ptB)
{
	DataVector<DPoint3d> init_points(2);
	init_points.add(ptA);
	init_points.add(ptB);
	auto bspline = std::make_shared<Curve3dBSpline>(init_points);
	if(bspline->fitToPointsAdaptive(points, 0) == 0)
		bspline->countNodes();
	return bspline;
}

/// fit to point-cloud using direct fit
Curve3dPtr Curve3dBSpline::fitThroughPoints(const DataVector<DPoint3d> & points, 
	const DPoint3d& ptA, const DPoint3d& ptB)
{
	DataVector<DPoint3d> init_points(2);
	init_points.add(ptA);
	init_points.add(ptB);
	auto bspline = std::make_shared<Curve3dBSpline>(init_points);
	bspline->fitThroughPointsAdaptive(points, 0);
	bspline->countNodes();
	return bspline;
}

/// adaptive subroutine - fit to point-cloud using direct fit, return number of inserted points
int Curve3dBSpline::fitToPointsAdaptive(const DataVector<DPoint3d> & points, int ind)
{
	size_t pct = points.countInt();
	if(pct < 5) return 0;
	assert(ind+1 < m_points.countInt());
	const DPoint3d& ptA = m_points[ind];
	const DPoint3d& ptB = m_points[ind+1];
	const DVector3d vt = ptB - ptA;
	const DPoint3d middle(ptA, ptB, 0.5);
	double c = vt.x * middle.x + vt.y * middle.y + vt.z * middle.z;
	// count distances from the ortoghonal line, crossing middle-point:
	// ... and select three points closest to the middle
	DataVector<double> distances(pct);
	DataVector<size_t> closest(3);
	for(size_t k = 0; k < pct; k++){
		const DPoint3d& pt = points[k];
		double dist = vt.x * pt.x + vt.y * pt.y + vt.z * pt.z - c;
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
	DPoint3d npoint = DPoint3d::zero;
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
	DataVector<DPoint3d> points_left(pct);
	DataVector<DPoint3d> points_right(pct);
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
		additional_left = fitToPointsAdaptive(points_left, ind);
//	}
	// ... and right
	int additional_right = 0;
//	if(max_dist_right > 1.0){
		// run
		additional_right = fitToPointsAdaptive(points_right, ind+additional_left+1);
//	}
	return 1 + additional_left + additional_right;
}

/// adaptive subroutine - fit to point-cloud using direct fit, return number of inserted points
int Curve3dBSpline::fitThroughPointsAdaptive(const DataVector<DPoint3d> & points, int ind)
{
	size_t pct = points.countInt();
	if(pct < 3) return 0;
	assert(ind+1 < m_points.countInt());
	const DPoint3d& ptA = m_points[ind];
	const DPoint3d& ptB = m_points[ind+1];
	const DVector3d vt = ptB - ptA;
	const DPoint3d middle(ptA, ptB, 0.5);
	double c = vt.x * middle.x + vt.y * middle.y + vt.z * middle.z;
	// count distances from the ortoghonal line, crossing middle-point:
	// ... and select the point closest to the middle
	DataVector<double> distances(pct);
	size_t closest = std::string::npos;
	for(size_t k = 0; k < pct; k++){
		const DPoint3d& pt = points[k];
		double dist = vt.x * pt.x + vt.y * pt.y + vt.z * pt.z - c;
		distances.add(dist);
		if(closest == std::string::npos || (abs(distances[closest]) > abs(dist))) 
			closest = k;
	}
	// insert and recalculate nodes
	m_points.insertAt(ind+1, points[closest]);
	countNodes();

	// check maximum distance (left and right), excluding the closest point
	DataVector<DPoint3d> points_left(pct);
	DataVector<DPoint3d> points_right(pct);
	for(size_t k = 0; k < pct; k++){
		if(k == closest) continue;
		if(distances[k] < 0.0){
			points_left.add(points[k]);
		}else{
			points_right.add(points[k]);
		}
	}

	// recurence left ...
	int additional_left = fitThroughPointsAdaptive(points_left, ind);
	// ... and right
	int additional_right = fitThroughPointsAdaptive(points_right, ind+additional_left+1);

	return 1 + additional_left + additional_right;
}

/// max distance
double Curve3dBSpline::maxDistanceForSubSegment(const DataVector<DPoint3d>& points, int ind, const DMetric3d& dm) const
{
	const DPoint3d& ptA = m_points[ind];
	const DVector3d vt = m_points[ind+1] - ptA;
	DVector3d vt1, vt2;
	vt.orthogonalVectors(vt1, vt2);
	// matrix
	DMatrix3d m;
	m.setColumn(0, vt);
	m.setColumn(1, vt1);
	m.setColumn(2, vt2);
	double max_dist2 = 0.0;
	DVector3d res;
	for(int k = 0; k < points.countInt(); k++){
		const DPoint3d& ptx = points[k];
		bool ok = m.solve(ptx - ptA, res);
		assert(ok);
		if(ok){
			double t = std::max(0.0, std::min(1.0, res.x));
			double dist2 = dm.transformRStoMS(ptx - getPoint((ind+t) / getMaxParameter())).length2();
			if(dist2 > max_dist2) max_dist2 = dist2;
		}
	}
	return sqrt(max_dist2);
}

/// max distance
double Curve3dBSpline::maxDistanceForSpline(const DataVector<DPoint3d>& points, const DMetric3d& dm) const
{
	size_t n = m_points.countInt()-1;
	DataVector<DMatrix3d> ms(n, DMatrix3d());
	for(size_t i = 0; i < n; i++){
		DVector3d vt(m_points[i], m_points[i+1]);
		DVector3d vt1, vt2;
		vt.orthogonalVectors(vt1, vt2);
		ms[i].setColumn(0, vt);
		ms[i].setColumn(1, vt1);
		ms[i].setColumn(2, vt2);
	}

	double max_dist2 = 0.0;
	DVector3d res;
	for(size_t k = 0; k < points.countInt(); k++){
		const DPoint3d& ptx = points[k];
		double min_dist2 = -1.0;
		for(size_t i = 0; i < n; i++){
			bool ok = ms[i].solve(ptx - m_points[i], res);
			assert(ok);
			if(!ok) continue;
			double dist2;
			double t = res.x;
			if(t < 0.0){
				dist2 = dm.transformRStoMS(ptx - m_points[i]).length2();
			}else if(t > 1.0){
				dist2 = dm.transformRStoMS(ptx - m_points[i+1]).length2();
			}else{
				dist2 = dm.transformRStoMS(ptx - getPoint((i+t) / n)).length2();
			}
			if(min_dist2 < 0.0 || dist2 < min_dist2) min_dist2 = dist2;
		}
		if(min_dist2 > max_dist2) max_dist2 = min_dist2;
	}
	return sqrt(max_dist2);
}

/// Store XML description to stream
ostream& Curve3dBSpline::storeXML(ostream& os, const string& prefix) const
{
	os << prefix << "<bspline3d-curve>\n";
	Curve3dParametric::storeXML(os, prefix + "\t");
	if(m_opened)
		os << prefix << "\t<points>\n";
	else
		os << prefix << "\t<points opened=\"false\">\n";
	for(int i = 0; i < m_points.countInt(); i++){
		os << prefix << "\t\t<point> " << m_points[i] << " </point>\n";
	}
	return os << prefix << "\t</points>\n" << prefix << "</bspline3d-curve>\n";
}
