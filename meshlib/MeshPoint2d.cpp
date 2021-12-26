/////////////////////////////////////////////////////////////////////////////
// MeshPoint2d.cpp
//  Represents two-dimensional vertex of mesh 
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	1999-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#include "MeshPoint2d.h"
#include "MeshEdge2d.h"
#include "MeshElement.h"
#include "MeshPoint3d.h"
#include "DHesjan.h"
#include "SurfaceParametric.h"
#include "MeshViewSet.h"
#include "ControlSpace2d.h"
#include "Metric2dContext.h"

int MeshPoint2d::param_max_node_count = 10000000;

//////////////////////////////////////////////////////////////////////
// Konstruktor
MeshPoint2d::MeshPoint2d(double _x, double _y) : coord(_x, _y) { }

//////////////////////////////////////////////////////////////////////
MeshPoint2d::MeshPoint2d(const DPoint2d& pt) : coord(pt) { }

//////////////////////////////////////////////////////////////////////
MeshPoint2d::MeshPoint2d(const MeshPoint2d& point) : 
	IndexedObject(), TagExtended(), coord(point.coord) { }

MeshPoint2d::MeshPoint2d(const MeshPoint2d* point) : 
	TagBorder(*point), coord(point->coord) { }

//////////////////////////////////////////////////////////////////////
// Destruktor
MeshPoint2d::~MeshPoint2d()
{
	while(edges.notEmpty())
		delete edges[0];	// Destructor MeshEdge2d will take care of removing this association
}

void MeshPoint2d::preDeleteAll()
{
	for(size_t i = 0; i < edges.countInt(); i++){
		MeshEdge2d* edge = edges[i];
		if(edge->removePointLink(this)) delete edge;
	}
	edges.clear();
}

/** 
 * Returns edge joining the given point with this one
 */
MeshEdge2d* MeshPoint2d::getEdgeToPoint(const MeshPoint2d *point) const
{
	for(size_t i = 0; i < edges.countInt(); i++) {
		if(edges[i]->getOtherPoint(this) == point) return edges[i];
	}
	// Points are not connected
	return nullptr;
}

/** 
 * Returns the number of elements incident to this point
 * with the given number of edges
 */
int MeshPoint2d::getElementsCount(int edges_ct) const
{
	int ct = 0;
	for(size_t i = 0; i < edges.countInt(); i++){
		MeshElement* element = edges[i]->getMeshElement(this);
		if(element && element->getEdgeCount() == edges_ct)
			++ct;
	}
	return ct;
}

/* 

DHesjan MeshPoint2d::countHesjan(int ct, bool relative, double ratio, int method)
{
	DHesjan hesjan;	// invalid so far

	if(isBorder()) return hesjan;	// border point

	LOG4CPLUS_ERROR(MeshLog::logger_console,   "TODO - points weight for hesjan");

	MeshPoint2d* points[9];
	double A[9][9];
	double f[9];

	int ect = edges_count;
	//if(ect > ct) ect = ct;
	int k = 0;
	bool any_found;
	int ik = 0;
	do{
		any_found = false;
		for(int j = 0; j < ect; j++){
			int ict = edges[j]->getInnerPointsCount();
			if(ict <= ik) continue;
			points[k] = edges[j]->getInnerMeshPoint(ik, this);
			hesjan.point_nrs[k] = points[k]->getID();
			hesjan.point_type[k] = -1;	// zerowa sfera ?
			any_found = true;
			if(++k >= ct) break;
		}
		++ik;
	}while(any_found && k < ct);
	if(k < ct){
		for(int j = 0; j < ect; j++){
			points[k] = edges[j]->getOtherPoint(this);
			hesjan.point_nrs[k] = points[k]->getID();
			hesjan.point_type[k] = 0;	// PIerwsza sfera
			if(++k >= ct) break;
		}
	}
	if(k < ct){
		// Dodatkowe punkty ...
		for(int j = 0; j < ect; j++){
			MeshElement* element = edges[j]->getMeshElement(this);
			if(element){	// Jeœli to jest wierzcho³ek wewnêtrzny
				MeshElement* element2 = element->getNextEdge(edges[j])->getOtherElement(element);
				if(element2){
					points[k] = element2->getNextPoint(points[j]);
					if(points[k]){
						hesjan.point_nrs[k] = points[k]->getID();
						hesjan.point_type[k] = 1;	// druga sfera
						if(++k >= ct) break;
					}
				}
			}
		}
	}

	if(k < ct){	// wci¹¿ za ma³o ???
		// Dodatkowe punkty ...
		for(int j = 0; j < ect; j++){
			int rank = points[j]->getRank();
			for(int l = 0; l < rank; l++){
				MeshPoint2d* point1 = points[j]->getEdge(l)->getOtherPoint(points[j]);
				bool found = (point1 == this);
				if(!found){
					for(int m = 0; m < k; m++){
						if(points[m] == point1){
							found = true;
							break;
						}
					}
				}
				if(!found){
					points[k] = point1;
					hesjan.point_nrs[k] = points[k]->getID();
					hesjan.point_type[k] = 2;	// druga sfera
					++k;
					break;						
				}
			}
			if(k >= ct) break;						
		}
	}

	if(k < ct){	// wci¹¿ za ma³o ???
		return hesjan;	// still invalid
	}
		
	double PA = 0.0, PB = 0.0, PC = 0.0, PD = 0.0;
	if(relative){
		// Wyznacz œredni wektor normalny (z normalnych dla trójk¹tów incydentnych do œrodkowego punktu)
		DPoint3d pt3_0(coord.x, coord.y, ratio * weight);
		DVector3d ave_normal;
		for(int j = 0; j < edges_count; j++){
			MeshElement* element = edges[j]->getMeshElement(this);	// element na lewo od krawêdzi
			assert(element);	// jeœli punkt nie jest brzegowy ...
			MeshPoint2d* point1 = element->getNextPoint(this);
			MeshEdge2d* edge = this->getEdgeToPoint(point1);
			if(edge->getInnerPointsCount() > 0){
				point1 = edge->getInnerMeshPoint(0, this);
			}
			MeshPoint2d* point2 = element->getPrevPoint(this);
			edge = this->getEdgeToPoint(point2);
			if(edge->getInnerPointsCount() > 0){
				point2 = edge->getInnerMeshPoint(0, this);
			}
			DPoint3d pt3_1(point1->getCoordinates(), 0.0);
			pt3_1.z = ratio * point1->getWeight();
			DPoint3d pt3_2(point2->getCoordinates(), 0.0);
			pt3_2.z = ratio * point2->getWeight();
			DVector3d v_normal = pt3_0.crossProduct(pt3_1, pt3_2).normalize();
			ave_normal += v_normal;
		}
		ave_normal.normalize();
		PA = ave_normal.x;
		PB = ave_normal.y;
		PC = ave_normal.z;
		//PD = -(PA*pt3_0.x + PB*pt3_0.y + PC*fvalue);	//P³aszczyzna Ax + By + Cz + D = 0
		PD = -(PA*pt3_0.x + PB*pt3_0.y + PC*pt3_0.z);	//P³aszczyzna Ax + By + Cz + D = 0
	}

	for(int j = 0; j < ct; j++){
		DPoint2d ptx = points[j]->getCoordinates();
		if(relative){
//			double z = points[j]->getFValue();
			double z = ratio * points[j]->getWeight();
			f[j] = (PA*ptx.x + PB*ptx.y + PC*z + PD);	// odleg³oœæ punktu od p³aszczyzny
		}else{
			f[j] = points[j]->getWeight();
		}
		double hx = ptx.x - coord.x;
		double hy = ptx.y - coord.y;
		A[0][j] = hx;
		A[1][j] = hy;
		A[2][j] = hx*hx;
		A[3][j] = hy*hy;
		A[4][j] = hx*hy;
		if(ct == 9){
			A[5][j] = A[2][j] * hx;
			A[6][j] = A[3][j] * hy;
			A[7][j] = A[2][j] * hy;
			A[8][j] = A[3][j] * hx;
		}
	}
	if(relative){
		hesjan.solveBG(ct, A, f, 0.0, method);
	}else{
		hesjan.solveBG(ct, A, f, weight, method);
	}
	return hesjan;
}

bool MeshPoint2d::countHesjanCurvature(int ct, DPoint3d& curvature, int method) const
{
	if(isBorder()) return false;	// border point

	LOG4CPLUS_ERROR(MeshLog::logger_console,   "TODO - points weight for hesjan curvature");

	const int MAX_PTS = 9;
	assert(ct == 5 || ct == MAX_PTS);

	MeshPoint2d* points[MAX_PTS];
	double AM[MAX_PTS][MAX_PTS];

	DHesjan hesjan;	// invalid so far
	int ect = edges_count;
	bool any_found;
	int k = 0;
	int ik = 0;

	do{
		any_found = false;
		for(int j = 0; j < ect; j++){
			int ict = edges[j]->getInnerPointsCount();
			if(ict <= ik) continue;
			points[k] = edges[j]->getInnerMeshPoint(ik, this);
			hesjan.point_nrs[k] = points[k]->getID();
			hesjan.point_type[k] = -1;	// zerowa sfera ?
			any_found = true;
			if(++k >= ct) break;
		}
		++ik;
	}while(any_found && k < ct);
	if(k < ct){
		for(int j = 0; j < ect; j++){
			points[k] = edges[j]->getOtherPoint(this);
			hesjan.point_nrs[k] = points[k]->getID();
			hesjan.point_type[k] = 0;	// PIerwsza sfera
			if(++k >= ct) break;
		}
	}
	if(k < ct){
		// Dodatkowe punkty ...
		for(int j = 0; j < ect; j++){
			MeshElement* element = edges[j]->getMeshElement(this);
			if(element){	// Jeœli to jest wierzcho³ek wewnêtrzny
				MeshElement* element2 = element->getNextEdge(edges[j])->getOtherElement(element);
				if(element2){
					points[k] = element2->getNextPoint(points[j]);
					if(points[k]){
						hesjan.point_nrs[k] = points[k]->getID();
						hesjan.point_type[k] = 1;	// druga sfera
						if(++k >= ct) break;
					}
				}
			}
		}
	}

	if(k < ct){	// wci¹¿ za ma³o ???
		// Dodatkowe punkty ...
		for(int j = 0; j < ect; j++){
			int rank = points[j]->getRank();
			for(int l = 0; l < rank; l++){
				MeshPoint2d* point1 = points[j]->getEdge(l)->getOtherPoint(points[j]);
				bool found = (point1 == this);
				if(!found){
					for(int m = 0; m < k; m++){
						if(points[m] == point1){
							found = true;
							break;
						}
					}
				}
				if(!found){
					points[k] = point1;
					hesjan.point_nrs[k] = points[k]->getID();
					hesjan.point_type[k] = 2;	// druga sfera
					++k;
					break;						
				}
			}
			if(k >= ct) break;						
		}
	}

	if(k < ct){	// wci¹¿ za ma³o ???
		return false;	// still invalid
	}

	// Count average normal vector (from normals of triangles incident to middle points)
	DPoint3d pt3_0(coord.x, coord.y, weight);
	DVector3d ave_normal;
	for(int j = 0; j < edges_count; j++){
		MeshElement* element = edges[j]->getMeshElement(this);	// element left from the edge
		assert(element);	// must be valid for each inner node
		MeshPoint2d* point1 = element->getNextPoint(this);
		MeshEdge2d* edge = this->getEdgeToPoint(point1);
		if(edge->getInnerPointsCount() > 0) point1 = edge->getInnerMeshPoint(0, this);
		MeshPoint2d* point2 = element->getPrevPoint(this);
		edge = this->getEdgeToPoint(point2);
		if(edge->getInnerPointsCount() > 0) point2 = edge->getInnerMeshPoint(0, this);
		const DPoint3d pt3_1(point1->getCoordinates(), point1->getWeight());
		const DPoint3d pt3_2(point2->getCoordinates(), point2->getWeight());
		ave_normal += pt3_0.crossProduct(pt3_1, pt3_2).normalize();
	}

	ave_normal.normalize();
	double f[MAX_PTS];

	for(int j = 0; j < ct; j++){
		const DPoint2d ptx = points[j]->getCoordinates();
		f[j] = points[j]->getWeight();

		double hx = ptx.x - coord.x;
		double hy = ptx.y - coord.y;

		AM[0][j] = hx;
		AM[1][j] = hy;
		AM[2][j] = hx*hx;
		AM[3][j] = hy*hy;
		AM[4][j] = hx*hy;
		if(ct == 9){
			AM[5][j] = AM[2][j] * hx;
			AM[6][j] = AM[3][j] * hy;
			AM[7][j] = AM[2][j] * hy;
			AM[8][j] = AM[3][j] * hx;
		}
	}
	hesjan.solveBG(ct, AM, f, weight, method, false);

	// count curvature

	const DPoint3d fs(1, 0, hesjan.Dx);   // = getDerivative(DEquation::deriv_ds, pt);
	const DPoint3d ft(0, 1, hesjan.Dy);   // = getDerivative(DEquation::deriv_dt, pt);
	//const DPoint3d fss(0, 0, hesjan.Dxx); // = getDerivative(DEquation::deriv_dss, pt);
	//const DPoint3d fst(0, 0, hesjan.Dxy); // = getDerivative(DEquation::deriv_dst, pt);
	//const DPoint3d ftt(0, 0, hesjan.Dyy); // = getDerivative(DEquation::deriv_dtt, pt);
	//const DPoint3d fn  = fs.crossProduct(ft).normalize();
	const DVector3d fn  = ave_normal;

	double g11 = 1.0 + sqr(hesjan.Dx);  // fs.scalarProduct(fs);
	double g12 = hesjan.Dx * hesjan.Dy;  // fs.scalarProduct(ft);
	double g22 = 1.0 + sqr(hesjan.Dy); // ft.scalarProduct(ft);
		
	double b11 = fn.z * hesjan.Dxx; // fn.scalarProduct(fss);
	double b12 = fn.z * hesjan.Dxy; // fn.scalarProduct(fst);
	double b22 = fn.z * hesjan.Dyy; // fn.scalarProduct(ftt);

	double A = g11*g22 - sqr(g12);
	double B = 2*b12*g12 - b22*g11 - b11*g22;
	double C = b11*b22 - sqr(b12);

	if(abs(A) > mesh_data.relative_small_number){
		double delta = B*B - 4.0*A*C;
		if(delta > mesh_data.relative_small_number){
			delta = sqrt(delta);
			double x1 = (-B - delta) / (2.0*A);
			double x2 = (-B + delta) / (2.0*A);

			double c11 = b11-x1*g11;
			double c12 = b12-x1*g12;
			double c22 = b22-x1*g22;
			// eigenvector
			if(abs(c12) > mesh_data.relative_small_number)
				curvature = DPoint3d(abs(x1), abs(x2), atan(-c11 / c12));
			else if(abs(c22) > mesh_data.relative_small_number)
				curvature = DPoint3d(abs(x1), abs(x2), atan(-c12 / c22));
			else
				curvature = DPoint3d(abs(x2), abs(x1), 0.0);
		}else if(delta > -mesh_data.relative_small_number){
			double x1 = -B / (2.0*A);
			curvature = DPoint3d(abs(x1), abs(x1), 0.0);
		}
	}else if(abs(B) > mesh_data.relative_small_number){
		double x1 = - C / B;
		curvature = DPoint3d(abs(x1), abs(x1), 0.0);
	}
	return true;
}

bool MeshPoint2d::countHesjanCurvatureWithBorder(int ct, DPoint3d& curvature, int method) const
{
//	if(isBorder() || edges_count < 2) return false;	// punkt wewnêtrzny krawêdzi

	LOG4CPLUS_ERROR(MeshLog::logger_console,   "TODO - points weight for Hessian curvature");

	const int MAX_PTS = 9;
	assert(ct == 5 || ct == MAX_PTS);

	MeshPoint2d* points[MAX_PTS];
	double AM[MAX_PTS][MAX_PTS];

	DHesjan hesjan;	// invalid so far
	int ect = edges_count;
	bool any_found;
	int k = 0;
	int ik = 0;

	do{
		any_found = false;
		for(int j = 0; j < ect; j++){
			int ict = edges[j]->getInnerPointsCount();
			if(ict <= ik) continue;
			points[k] = edges[j]->getInnerMeshPoint(ik, this);
			hesjan.point_nrs[k] = points[k]->getID();
			hesjan.point_type[k] = -1;	// zerowa sfera ?
			any_found = true;
			if(++k >= ct) break;
		}
		++ik;
	}while(any_found && k < ct);
	if(k < ct){
		for(int j = 0; j < ect; j++){
			points[k] = edges[j]->getOtherPoint(this);
			hesjan.point_nrs[k] = points[k]->getID();
			hesjan.point_type[k] = 0;	// PIerwsza sfera
			if(++k >= ct) break;
		}
	}
	if(k < ct){
		// Dodatkowe punkty ...
		for(int j = 0; j < ect; j++){
			MeshElement* element = edges[j]->getMeshElement(this);
			if(element){	// Jeœli to jest wierzcho³ek wewnêtrzny
				MeshElement* element2 = element->getNextEdge(edges[j])->getOtherElement(element);
				if(element2){
					points[k] = element2->getNextPoint(points[j]);
					if(points[k]){
						hesjan.point_nrs[k] = points[k]->getID();
						hesjan.point_type[k] = 1;	// druga sfera
						if(++k >= ct) break;
					}
				}
			}
		}
	}

	if(k < ct){	// wci¹¿ za ma³o ???
		// Dodatkowe punkty ...
		for(int j = 0; j < ect; j++){
			int rank = points[j]->getRank();
			for(int l = 0; l < rank; l++){
				MeshPoint2d* point1 = points[j]->getEdge(l)->getOtherPoint(points[j]);
				bool found = (point1 == this);
				if(!found){
					for(int m = 0; m < k; m++){
						if(points[m] == point1){
							found = true;
							break;
						}
					}
				}
				if(!found){
					points[k] = point1;
					hesjan.point_nrs[k] = points[k]->getID();
					hesjan.point_type[k] = 2;	// druga sfera
					++k;
					break;						
				}
			}
			if(k >= ct) break;						
		}
	}

	if(k < ct){	// wci¹¿ za ma³o ???
		return false;	// still invalid
	}

	// Count average normal vector (from normals of triangles incident to middle points)
	DPoint3d pt3_0(coord.x, coord.y, weight);
	DVector3d ave_normal;
	for(int j = 0; j < edges_count; j++){
		MeshElement* element = edges[j]->getMeshElement(this);	// element left from the edge
		if(!element) continue;	// must be valid for each inner node
		MeshPoint2d* point1 = element->getNextPoint(this);
		MeshEdge2d* edge = this->getEdgeToPoint(point1);
		if(edge->getInnerPointsCount() > 0) point1 = edge->getInnerMeshPoint(0, this);
		MeshPoint2d* point2 = element->getPrevPoint(this);
		edge = this->getEdgeToPoint(point2);
		if(edge->getInnerPointsCount() > 0) point2 = edge->getInnerMeshPoint(0, this);
		const DPoint3d pt3_1(point1->getCoordinates(), point1->getWeight());
		const DPoint3d pt3_2(point2->getCoordinates(), point2->getWeight());
		ave_normal += pt3_0.crossProduct(pt3_1, pt3_2).normalize();
	}

	ave_normal.normalize();
	double f[MAX_PTS];

	for(int j = 0; j < ct; j++){
		const DPoint2d ptx = points[j]->getCoordinates();
		f[j] = points[j]->getWeight();

		double hx = ptx.x - coord.x;
		double hy = ptx.y - coord.y;

		AM[0][j] = hx;
		AM[1][j] = hy;
		AM[2][j] = hx*hx;
		AM[3][j] = hy*hy;
		AM[4][j] = hx*hy;
		if(ct == 9){
			AM[5][j] = AM[2][j] * hx;
			AM[6][j] = AM[3][j] * hy;
			AM[7][j] = AM[2][j] * hy;
			AM[8][j] = AM[3][j] * hx;
		}
	}
	hesjan.solveBG(ct, AM, f, weight, method, false);

	// count curvature

	const DVector3d fs(1, 0, hesjan.Dx);   // = getDerivative(DEquation::deriv_ds, pt);
	const DVector3d ft(0, 1, hesjan.Dy);   // = getDerivative(DEquation::deriv_dt, pt);
	//const DPoint3d fss(0, 0, hesjan.Dxx); // = getDerivative(DEquation::deriv_dss, pt);
	//const DPoint3d fst(0, 0, hesjan.Dxy); // = getDerivative(DEquation::deriv_dst, pt);
	//const DPoint3d ftt(0, 0, hesjan.Dyy); // = getDerivative(DEquation::deriv_dtt, pt);
	//const DPoint3d fn  = fs.crossProduct(ft).normalize();
	const DVector3d fn  = ave_normal;

	double g11 = 1.0 + sqr(hesjan.Dx);  // fs.scalarProduct(fs);
	double g12 = hesjan.Dx * hesjan.Dy;  // fs.scalarProduct(ft);
	double g22 = 1.0 + sqr(hesjan.Dy); // ft.scalarProduct(ft);
		
	double b11 = fn.z * hesjan.Dxx; // fn.scalarProduct(fss);
	double b12 = fn.z * hesjan.Dxy; // fn.scalarProduct(fst);
	double b22 = fn.z * hesjan.Dyy; // fn.scalarProduct(ftt);

	double A = g11*g22 - sqr(g12);
	double B = 2*b12*g12 - b22*g11 - b11*g22;
	double C = b11*b22 - sqr(b12);

	if(abs(A) > mesh_data.relative_small_number){
		double delta = B*B - 4.0*A*C;
		if(delta > mesh_data.relative_small_number){
			delta = sqrt(delta);
			double x1 = (-B - delta) / (2.0*A);
			double x2 = (-B + delta) / (2.0*A);

			double c11 = b11-x1*g11;
			double c12 = b12-x1*g12;
			//double c22 = b22-x1*g22;
			DVector2d v;
			// eigenvector
			if(abs(c12) < mesh_data.relative_small_number){
				if(abs(c11) > mesh_data.relative_small_number){
					v.x = 0.0;
					v.y = 1.0;
				}else{
					v.x = 1.0;
					v.y = 0.0;
				}
			}else{
				v.x = c12;
				v.y = -c11;
				v.normalize();
			}
			DVector3d dp = fs * v.x + ft * v.y;
			double alpha = fs.getAngle(dp);
			curvature = DPoint3d(abs(x1), abs(x2), alpha);
		}else if(delta > -mesh_data.relative_small_number){
			double x1 = -B / (2.0*A);
			curvature = DPoint3d(abs(x1), abs(x1), 0.0);
		}
	}else if(abs(B) > mesh_data.relative_small_number){
		double x1 = - C / B;
		curvature = DPoint3d(abs(x1), abs(x1), 0.0);
	}
	return true;
}
*/

std::shared_ptr<MeshViewPointData> MeshPoint2d::getViewData(SurfaceConstPtr surface) const
{
	if(!surface)
		return std::make_shared<MeshViewPointData>(DPoint3d(coord, 0.0), isBorder(), index, true);
	else{
		MeshPoint3d* point3d = (MeshPoint3d*)getPtrTag(TAG_MP_2D_3D);
		if(point3d)
			return std::make_shared<MeshViewPointData>(point3d->getCoordinates(), isBorder(), index, true);
		else 
			return std::make_shared<MeshViewPointData>(surface->getPoint(coord), isBorder(), index, true);
	}
}

const DPoint2d MeshPoint2d::getRealCoordinates(Metric2dContext& mc) const { 
	return mc.transformPStoRS(coord); 
}

const DMPoint2d MeshPoint2d::getMetricCoordinates(Metric2dContext& mc)
{
	return mc.transformPStoMS(coord);
}
