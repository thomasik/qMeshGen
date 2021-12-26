// MeshSpecialRoutinesDAT.cpp: implementation of the MeshSpecialRoutinesDAT class.
//
//////////////////////////////////////////////////////////////////////

#include "MeshData.h"
#include "MeshSpecialRoutinesDAT.h"
#include "MeshContainer2d.h"
#include "MeshContainer3d.h"
#include "MeshContainer3dSurface.h"
#include "MeshPoint2d.h"
#include "MeshPoint3d.h"
#include "MeshEdge2d.h"
#include "MeshEdge2dCurve.h"
#include "MeshEdge3d.h"
#include "MeshDomainEdge3d.h"
#include "MeshArea.h"
#include "MeshFace.h"
#include "MeshDomainSurface.h"
#include "MeshDomainVolume.h"
#include "MeshQuad2d.h"
#include "Curve2dAnalytic.h"
#include "SurfacePlane.h"
#include "ControlSpace2dAnalytic.h"
#include "ControlSpace2dMatrixUniform.h"
#include "ControlSpace2dMesh.h"
#include "ControlSpace2dIdentity.h"
#include "MeshTriangle3d.h"
#include "Quad12.h"
#include "MeshGenerator1d.h"
#include "MeshGenerator2d.h"
#include "MeshGenerator2dQuad.h"
#include "ControlSpace3dOctree.h"
#include "DHesjan.h"
#include "DataHashTable.h"
#include "DSegment.h"
#include "DTriangle.h"
#include "DLine.h"

/**
 * Creates specific geometry - triangular part of a regular pentagon:
 * \param a - length of the edge of the pentagon
 * \param r - radius of the cutting circle
 * \param dr1 - radius of the area with small elements
 * \param dr2 - radius of the transition area between small and large elements
 * \param l1 - size (edge length) of "small" elements
 * \param l2 - size (edge length) of "large" elements
 * \param with_user_control - whether special, user defined control-space should be used
 * \param shape_type - 0 for triangle, 1 for rectangle
 * \return generated boundary description with control space
MeshContainer3d* MeshSpecialRoutinesDAT::createBoundaryForPhaseTransformation(double a, double r, double dr1, double dr2, double dl1, double dl2, bool with_user_control, int shape_type)
{
	coefficient = 1.0;
	while(a > 10){
		a /= 10;
		coefficient *= 10;
	}
	while(a < 0.1){
		a *= 10;
		coefficient /= 10;
	}

	r /= coefficient;
	dr1 /= coefficient;
	dr2 /= coefficient;
	dl1 /= coefficient;
	dl2 /= coefficient;

	double alpha = 0.3*PI;
	double AB = 0.5 * a;
	double yAB = AB * sin(alpha);
	double xAB = AB * cos(alpha);
	double xAC = AB / cos(alpha);
	
	const int max_pts_ct = 4;
	int pts_ct = 3;
	int i;

	// *** 2d ***
	// points
	MeshPoint2d** points2d = new MeshPoint2d*[max_pts_ct];
	if(shape_type == 0){
		if(r <= 0.0){
			points2d[0] = new MeshPoint2d(0.0, 0.0); // id=1
			points2d[1] = new MeshPoint2d(xAC, 0.0); // id=2
			points2d[2] = new MeshPoint2d(xAB, yAB); // id=3
		}else if(r < AB){
			pts_ct = 4;
			points2d[0] = new MeshPoint2d(  r, 0.0);
			double ratio = r / AB;
			points2d[1] = new MeshPoint2d(xAC, 0.0);
			points2d[2] = new MeshPoint2d(xAB, yAB);
			points2d[3] = new MeshPoint2d(ratio*xAB, ratio*yAB);
		}else if(r < xAC){
			pts_ct = 3;
			double WA = sqr(xAB-xAC)+sqr(yAB);
			double WB = 2*(xAB*(xAC-xAB)-yAB*yAB);
			double WC = xAB*xAB+yAB*yAB-r*r;
			double delta = WB*WB-4*WA*WC;
			assert(delta > 0.0);
			delta = sqrt(delta);
			double t = (delta-WB)/(2*WA);
			assert( t >= 0.0 && t <= 1.0);
			xAB = (1-t)*xAB + t*xAC;
			yAB = (1-t)*yAB;

			points2d[0] = new MeshPoint2d(  r, 0.0);
			points2d[1] = new MeshPoint2d(xAC, 0.0);
			points2d[2] = new MeshPoint2d(xAB, yAB);

			alpha = asin(yAB/r);
		}else 
			return nullptr;
	}else{
		pts_ct = 4;
		points2d[0] = new MeshPoint2d(r, a);
		points2d[1] = new MeshPoint2d(r, 0);
		points2d[2] = new MeshPoint2d(a, 0);
		points2d[3] = new MeshPoint2d(a, a);		
	}
	//boundary & edges
	Curve2dParametric* curve = nullptr;
	if(r > 0.0){
		string str_x, str_y;
		DEquation::doubleToString(r, DEquation::v_value, str_x);
		str_x += "*COS(X)";
		DEquation::doubleToString(-r, DEquation::v_value, str_y);
		str_y += "*SIN(X)";
		curve = new Curve2dAnalytic(str_x,str_y);
	}
	MeshContainer2d* boundary = new MeshContainer2d(max_pts_ct);
	MeshEdge2d* edges2d[max_pts_ct];
	for(i = 0; i < pts_ct; i++){
		points2d[i]->setBorder(TagBorder::OUTER | TagBorder::FIXED);
		if(shape_type == 0){
			if(curve && i == pts_ct-1){
				edges2d[i] = new MeshEdge2dCurve(points2d[i], points2d[(i+1)%pts_ct],
					curve, 2*PI-alpha, 2*PI);
				edges2d[i]->setBorderType(1); // TODO -> if ever needed, should be boundary
			}else{
				edges2d[i] = new MeshEdge2d(points2d[i], points2d[(i+1)%pts_ct]);
				edges2d[i]->setBorderType(2);
			}
		}else{
			edges2d[i] = new MeshEdge2d(points2d[i], points2d[(i+1)%pts_ct]);
			edges2d[i]->setBorderType(i==0?1:2);
		}
		boundary->addMeshPoint(points2d[i]);
	}
	// surface
	SurfaceParametric* surface = new SurfacePlane(DVector3d(1.0, 0.0, 0.0), DVector3d(0.0, 1.0, 0.0));
	// boundary continued
	MeshArea* area = new MeshArea(pts_ct, points2d);
	area->setFilled();
	boundary->addMeshElement(area);

	// *** 3d ***
	MeshContainer3d* mesh = new MeshContainer3d(max_pts_ct);
	// points
	DataVector<MeshPoint3d*> points3d(pts_ct);
	DataVector<MeshDomainEdge3d*>  edges3d(pts_ct);
	for(i = 0; i < pts_ct; i++){
		MeshPoint3d* point3d = new MeshPoint3d(points2d[i]->getCoordinates().x, points2d[i]->getCoordinates().y, 0.0);
		points3d.add(point3d);
		points2d[i]->setPtrTag(TagExtended::TAG_MP_2D_3D, point3d);
		mesh->addMeshPoint(point3d);
	}
	for(i = 0; i < pts_ct; i++){
		MeshDomainEdge3d* edge3d = new MeshDomainEdge3d(points3d[i], points3d[(i+1)%pts_ct]);
		edges3d.add(edge3d);
		edges2d[i]->setPtrTag(TagExtended::TAG_ME_2D_3D, edge3d);
		edge3d->setBorder(edges2d[i]->getBorderType());
	}
	// face
	DataVector<MeshFace*> faces(1);
	boundary->setSurface(surface);
	faces.add(new MeshDomainSurface(points3d, edges3d, surface, boundary));
	// block
	MeshDomainVolume* volume = new MeshDomainVolume(faces, points3d);
	mesh->addMeshBlock(volume);

	// Special control space
	if(with_user_control && r+dr1 > 0.0){
		string len_str, angle_str = "0";
		DEquation::doubleToString(dl2, DEquation::v_length, len_str);
		ControlSpace2dAnalytic* user_control = new ControlSpace2dAnalytic(surface, len_str.c_str(), len_str.c_str(), angle_str.c_str());

		string condition_str;
		if(shape_type == 0){
			DEquation::doubleToString(sqr(r+dr1+dr2), DEquation::v_area, condition_str);
			condition_str += "-X*X-Y*Y";
			string str_l1, str_ratio, str_dr;
			DEquation::doubleToString(r+dr1, DEquation::v_length, str_dr);
			DEquation::doubleToString((dl2-dl1)/dr2, DEquation::v_value, str_ratio);
			DEquation::doubleToString(dl1, DEquation::v_length, str_l1);
			str_l1 += "+((X*X+Y*Y)^0.5-"+str_dr+")*"+str_ratio;
			user_control->addCaseAsFirst(condition_str.c_str(), str_l1.c_str(), str_l1.c_str(), angle_str.c_str());

			DEquation::doubleToString(sqr(r+dr1), DEquation::v_area, condition_str);
			condition_str += "-X*X-Y*Y";
			DEquation::doubleToString(dl1, DEquation::v_length, len_str);
			user_control->addCaseAsFirst(condition_str.c_str(), len_str.c_str(), len_str.c_str(), angle_str.c_str());
		}else{
			DEquation::doubleToString(r+dr1+dr2, DEquation::v_length, condition_str);
			condition_str += "-X";
			string str_l1, str_ratio, str_dr;
			DEquation::doubleToString(r+dr1, DEquation::v_length, str_dr);
			DEquation::doubleToString((dl2-dl1)/dr2, DEquation::v_value, str_ratio);
			DEquation::doubleToString(dl1, DEquation::v_length, str_l1);
			str_l1 += "+(X-"+str_dr+")*"+str_ratio;
			user_control->addCaseAsFirst(condition_str.c_str(), str_l1.c_str(), len_str.c_str(), angle_str.c_str());

			DEquation::doubleToString(r+dr1, DEquation::v_length, condition_str);
			condition_str += "-X";
			string str_l2;
			DEquation::doubleToString(dl1, DEquation::v_length, str_l2);
			user_control->addCaseAsFirst(condition_str.c_str(), str_l2.c_str(), len_str.c_str(), angle_str.c_str());
		}

		boundary->setControlSpace(user_control);
	}

	return mesh;
}
*/
/*
bool MeshSpecialRoutinesDAT::exportToDAT(MeshContainer2d* mesh, const string& fname, bool mark_boundary_points)
{
//	if(mesh->getInnerEdgesCount() != 2){
		LOG4CPLUS_ERROR(MeshLog::logger_console,   "No inner nodes at edges!");
		return false;
//	}

	mesh->renumeratePointsByNeighbours();
	mesh->renumeratePoints();

	ofstream file((fname+"p.dat").c_str());
	if(!file) return false;

	int pcount = mesh->getPointsCount();
	int ecount = mesh->getElementsCount();
	int p4count = 0;
	int i, j;
	double coeff = coefficient * 1000;	// [mm]
	LOG4CPLUS_ERROR(MeshLog::logger_console,   "TODO - points weight for DAT storage");
	for(i = 0; i < pcount; i++){
		MeshPoint2d* point = mesh->getPointAt(i);
		DPoint2d pt = point->getCoordinates();
		if(mark_boundary_points){
			int rank = point->getRank();
			int boundary_code = 0;
			for(j = 0; j < rank; j++){
				MeshEdge2d* edge = point->getEdge(j);
				if(edge->getBorderType() == 1){
					boundary_code = 1;
					break;
				}
			}
			file << (pt.x * coeff) << '\t' << (pt.y * coeff) << '\t' // << point->getWeight()
				<< '\t' << boundary_code << endl;
		}else{
			file << (pt.x * coeff) << '\t' << (pt.y * coeff) << '\t' // << point->getWeight() 
				<< endl;
		}
		if(point->getRank() > 1) ++p4count;
	}

	ofstream i_file((fname+"i.dat").c_str());
	if(!i_file) return false;

	StatData span = mesh->getElementsIDSpan();
	i_file << "Liczba pkt(12), liczba pkt(4), liczba elementow, szerokosc pasma\n";
	i_file << pcount << '\t' << p4count << '\t' << ecount << '\t' << ((int)span.maximum+1) << endl;
	i_file << " EL.   WEZEL     1     2     3     4     5     6     7     8     9    10    11    12           KOD\n";

	for(i = 0; i < ecount; i++){
		MeshQuad2d* quad = (MeshQuad2d*)mesh->getElementAt(i);
		if(quad->getEdgeCount() != 4){
			LOG4CPLUS_ERROR(MeshLog::logger_console,   "DAT store - triangle found.");
			return false;
		}
		assert(quad->getType() == ELEMENT_MESH_QUAD);

		i_file << (i+1);
		int has_marked_edge = 0;
		int has_boundary_edge = 0;
		for(j = 0; j < 4; j++){
			int border_type = quad->getEdge(j)->getBorderType();
			if(border_type >=0)
				++has_boundary_edge;
			if(border_type == 1)
				++has_marked_edge;
		}
		// boundary type == 0 -> powietrze  [typ krawêdzi 1 lub 2]
		// boundary type == 1 -> woda		[typ krawêdzi 3 lub 4]
		// boundary type == 2 -> izolacja   [krawêdzie brzegowe jak wewnêtrzne]
		int i_start = -1;
		if(has_marked_edge > 2){
			LOG4CPLUS_WARN(MeshLog::logger_console, "More then 2 marked edges");
		}
		if(has_marked_edge){
			// zapewniæ, ¿eby ta krawêdŸ mia³a kod 4 (lub 34 dla dwóch)
			while(quad->getEdge((++i_start+3)%4)->getBorderType() != 1);
			if(quad->getEdge(i_start%4)->getBorderType() == 1) i_start+=3;
		}else if(has_boundary_edge){
			// zapewniæ, ¿eby ten element zaczyna³ siê od 1
			do{
				++i_start;
			}while(quad->getEdge(i_start)->getBorderType() < 0 || 
				quad->getEdge((i_start+3)%4)->getBorderType() >= 0);
		}else{
			i_start = 0;
		}
		for(j = 0; j < 4; j++){
			int index = quad->getPoint((i_start+j)%4)->getIndex();
			i_file << '\t' << (index+1);
		}
		MeshEdge2d* edge = quad->getEdge(i_start%4);
		MeshPoint2d* point = quad->getPoint(i_start%4);

		int code = 0;
		for(j = 3; j >= 0; j--){
			int border_type = quad->getEdge((i_start+j)%4)->getBorderType();
			if(border_type == 0 || border_type == 1){
				if(code != 0) code *= 10;
				code += (j+1);
			}
		}

		i_file << '\t' << code << endl;
	}

	return true;
}
*/

/*
MeshContainer2d* MeshSpecialRoutinesDAT::importFromDAT(const string &fname)
{
	return importFromDAT(fname+"p.dat", fname+"i.dat");
}

MeshContainer2d* MeshSpecialRoutinesDAT::importFromDAT(const string& fname_p, const string& fname_i)
{
	ifstream file_p(fname_p.c_str());
	if(!file_p) return nullptr;
	ifstream file_i(fname_i.c_str());
	if(!file_i)	return nullptr;

	MeshContainer2d* mesh = new MeshContainer2d(100);

	// Wspó³rzêdne punktów i temperatura
	int pct = 0;
	double coeff = coefficient * 1000;	// [mm]
	while(file_p){
		double x, y, v;
		int b = 0;
		string buffer;
		getline(file_p, buffer);
		istringstream is(buffer);
		is >> x >> y >> v >> b;
		x /= coeff;	// mm
		y /= coeff;	// mm
		LOG4CPLUS_ERROR(MeshLog::logger_console,   "TODO - points weight for DAT");
		MeshPoint2d* point = new MeshPoint2d(x, y); //, v);
		mesh->addMeshPoint(point);
		//if(b > 0) point->setTag(1);
	}

	//fprintf(file, "X\tY\tV\n");
	//fprintf(file, "QUAD_NR W1-W12\n");
	// Tablica incydencji
	while(file_i){
		int quad_id;
		int w1,w2,w3,w4,w5,w6,w7,w8,w9,w10,w11,w12;
		int code;
		file_i >> quad_id >> w1 >> w2 >> w3 >> w4 >> w5 >> w6 >> w7 >> w8
			>> w9 >> w10 >> w11 >> w12 >> code;
		if(file_i){
			if(w1 > pct || w2 > pct || w3 > pct || w4 > pct ||
				w5 > pct || w6 > pct || w7 > pct || w8 > pct ||
				w9 > pct || w10 > pct || w11 > pct || w12 > pct) continue;
			MeshPoint2d* point1 = mesh->getPointAt(w1-1);
			MeshPoint2d* point2 = mesh->getPointAt(w2-1);
			MeshPoint2d* point3 = mesh->getPointAt(w3-1);
			MeshPoint2d* point4 = mesh->getPointAt(w4-1);
			MeshQuad2d* quad = new MeshQuad2d(point1, point2, point3, point4);
			quad->setAreaID(0);
			mesh->addMeshElement(quad);
			// Dodanie wierzcho³ków wewnêtrznych
			MeshEdge2d* edge = point1->getEdgeToPoint(point2);
		}
	}
//	mesh->setInnerEdgesCount(2);

	// Kontrola krawêdzi brzegowych
	for(int i = 0; i < pct; i++){
		MeshPoint2d* point = mesh->getPointAt(i);
		int rank = point->getRank();
		for(int j = 0; j < rank; j++){
			MeshEdge2d* edge = point->getEdge(j);
			if(!edge->isBorder() && (
					edge->getMeshElement(0) == nullptr ||
					edge->getMeshElement(1) == nullptr)){
				edge->setBorderType(0);
				point->setBorder();
				MeshPoint2d* other_point = edge->getOtherPoint(point);
				other_point->setBorder();
//				if(point->isTagged() && other_point->isTagged())
//					edge->setBorderType(1);
			}
		}
	}

	return mesh;
}
*/

bool MeshSpecialRoutinesDAT::interpolateValues(MeshContainer2d *dest_mesh, MeshContainer2d *source_mesh)
{
	int i, j, k, pct = dest_mesh->getPointsCount();
	int v_ect = source_mesh->getElementsCount();
	int v_pct = source_mesh->getPointsCount();

	Quad12 quad12;

	for(i = 0; i < pct; i++){
		MeshPoint2d* point = dest_mesh->getPointAt(i);
		DPoint2d pt = point->getCoordinates();
//		if(!view->controlStep(2, "Point to set value", pt)) return false;
		for(j = 0; j < v_ect; j++){
			MeshQuad2d* quad = (MeshQuad2d*)source_mesh->getElementAt(j);
			if(quad->isPointInside(point->getCoordinates())){
//				if(!view->controlStep(2, "Found value quad", quad->getMiddlePoint())) return false;
				quad12.setData(quad);
				double value = quad12.interpolateQuadValue(point->getCoordinates());
				LOG4CPLUS_ERROR(MeshLog::logger_console,   "TODO - points weight for DAT");
				//point->setWeight(value);

				break;
			}
		}
		if(j >= v_ect){	// nie nale¿y do ¿adnego elementu
			double min_distance = mesh_data.relative_infinity;
			MeshEdge2d* best_edge = nullptr;
			for(j = 0; j < v_pct; j++){
				MeshPoint2d* vpoint = source_mesh->getPointAt(j);
				int rank = vpoint->getRank();
				if(rank > 1){
					for(k = 0; k < rank; k++){
						MeshEdge2d* edge = vpoint->getEdge(k);
						if(edge->isBorder() && edge->getPointIndex(vpoint) == 0){
							MeshPoint2d* vpoint2 = edge->getMeshPoint(1);
							double distance2 = pt.distance2(vpoint->getCoordinates()) + pt.distance2(vpoint2->getCoordinates());
							if(distance2 < min_distance){
								best_edge = edge;
								min_distance = distance2;
							}
						}
					}
				}
			}
			assert(best_edge != nullptr);

//			if(!view->controlStep(2, "Found best edge", best_edge->getPoint(0.5))) return false;

			MeshQuad2d* quad = (MeshQuad2d*)best_edge->getMeshElement(0);
			if(!quad) quad = (MeshQuad2d*)best_edge->getMeshElement(1);
			assert(quad);
			//
			quad12.setData(quad);
			const DPoint2d& pt1 = best_edge->getMeshPoint(0)->getCoordinates();
			const DPoint2d& pt2 = best_edge->getMeshPoint(1)->getCoordinates();
			const DPoint2d& pt3 = point->getCoordinates();
			DPoint2d pt4(pt3.x + (pt2.y - pt1.y), pt3.y + (pt1.x - pt2.x));
			DPoint2d ptx = DLine2d::crossPoint(pt1, pt2, pt3, pt4);
			DVector2d dpt = ptx - pt3;
			pt4 = ptx + dpt * 0.2;
			double value = quad12.interpolateQuadValue(ptx);
			double value4 = quad12.interpolateQuadValue(pt4);
			LOG4CPLUS_ERROR(MeshLog::logger_console,   "TODO - points weight for DAT");
			//point->setWeight(value + 5*(value - value4));
		}
	}

	return true;
}

double MeshSpecialRoutinesDAT::coefficient = 1.0;

double MeshSpecialRoutinesDAT::countNextPhaseRadius(double old_r, double a, double ds, int shape_type)
{
	if(shape_type == 1){
		return old_r + ds / a;
	}

	assert(shape_type == 0);

	double alpha = 0.3*PI;
	double a_half = 0.5*a;
	if(old_r < a_half){
		double r_next = sqrt(2*ds/alpha + old_r*old_r);
		if(r_next <= a_half) return r_next;
		ds -= 0.5*alpha*(a_half*a_half - old_r*old_r);
		old_r = a_half;
	}
	
	double last_r = old_r;
	double last_area = getPhaseArea(last_r, a);
	double required_area = last_area+ds;
	double ratio = required_area / last_area;
	int counter = 0;

	while(abs(1.0 - ratio) > 1e-6 && counter < 100){
		++counter;
		double new_r = last_r * ratio;
		if(new_r > old_r){
			last_r = new_r;
		}else{
			last_r = 0.5*(last_r + old_r);
		}
		last_area = getPhaseArea(last_r, a);
		ratio = required_area / last_area;
	}

	return last_r;
}

double MeshSpecialRoutinesDAT::getPhaseArea(double r, double a)
{
	double alpha = 0.3*PI;
	double AB = 0.5 * a;

	if(r <= AB)	return 0.5*alpha*r*r;
	
	const DPoint2d A(0.0, 0.0);
	const DPoint2d B(AB * cos(alpha), AB * sin(alpha));
	const DPoint2d C(AB / cos(alpha), 0.0);

	if(r > C.x) return DTriangle2d::area(A, C, B);

	double WA = sqr(B.x-C.x)+sqr(B.y);
	double WB = 2*(B.x*(C.x-B.x)-B.y*B.y);
	double WC = B.x*B.x+B.y*B.y-r*r;
	double delta = WB*WB-4*WA*WC;
	assert(delta > 0.0);
	delta = sqrt(delta);
	double t = (delta-WB)/(2*WA);
	assert( t >= 0.0 && t <= 1.0);
	DPoint2d D(B, C, t);
	double alpha_d = asin(D.y/r);

	double area1 = DTriangle2d::area(A, D, B);
	double area2 = 0.5*alpha_d*r*r;

	return area1+area2;
}

bool MeshSpecialRoutinesDAT::setControlSpaceFromDATHesjan(MeshContainer2d *boundary, MeshContainer2d *value_mesh, double h_min_ratio, double h_factor)
{
	if(!boundary || !value_mesh) return false;

/* TODO

	DRect rect = boundary->getBoundingRect();
	rect.inflate(ControlSpace2d::param_inflate_box_factor);

	DHesjan::m_trim_max = true;
	DHesjan::m_trim_min = true;
	DHesjan::m_with_sqrt = true;
	if(h_factor < 0.0) h_factor = DHesjan::param_hesjan_factor;
	if(h_min_ratio < 0.0) h_min_ratio = DHesjan::param_hesjan_min_ratio;
	DHesjan::m_min_len = std::min(rect.getDX(), rect.getDY()) * h_min_ratio;
	DHesjan::m_max_len = std::min(rect.getDX(), rect.getDY()) * DHesjan::param_hesjan_max_ratio;
	DHesjan::m_factor = h_factor;

	int nx = ControlSpace2dMatrixUniform::param_uniform_nx;
	int ny = nx;
	ControlSpace2dMatrixUniform* space = new ControlSpace2dMatrixUniform(boundary->getSurface(), rect, nx, ny); 

	int i, pct = value_mesh->getPointsCount();
	for(i = 0; i < pct; i++){
		MeshPoint2d* point = value_mesh->getPointAt(i);
		DHesjan hesjan = point->countHesjan(9);
		if(hesjan.valid){
			DPoint2d pt = point->getCoordinates();
			// scale and trim
			hesjan.calculateLength();
			// insert into the control space
			space->addControlPoint(pt, DMetric2d::stretchToMatrix(hesjan.stretch));
		}
	}
	// Slight addition for boundary nodes
	for(i = 0; i < pct; i++){
		MeshPoint2d* point = value_mesh->getPointAt(i);
		if(point->isBorder()){
			int ect = point->getRank();
			double min_dist = mesh_data.relative_infinity;
			MeshPoint2d* selected_point = nullptr;
			DPoint2d dpt = point->getCoordinates();
			for(int j = 0; j < ect; j++){
				MeshPoint2d* other_point = point->getEdge(j)->getOtherPoint(point);
				if(!other_point->isBorder()){
					double d = dpt.distance2(other_point->getCoordinates());
					if(d < min_dist){
						min_dist = d;
						selected_point = other_point;
					}
				}
			}
			if(selected_point){
				DHesjan hesjan = selected_point->countHesjan(9);
				if(hesjan.valid){
					// scale and trim
					hesjan.calculateLength();
					// insert into the control space
					space->addControlPoint(dpt, DMetric2d::stretchToMatrix(hesjan.stretch));
				}
			}
		}
		DHesjan hesjan = point->countHesjan(9);
		if(hesjan.valid){
			DPoint2d pt = point->getCoordinates();
			// scale and trim
			hesjan.calculateLength();
			// insert into the control space
			space->addControlPoint(pt, DMetric2d::stretchToMatrix(hesjan.stretch));
		}
	}

	space->interpolate();
	boundary->setControlSpace(space);
*/

	return true;
}

//bool MeshSpecialRoutinesDAT::checkPhaseCorrectness(MeshContainer2d *mesh)
//{
//	int ect = mesh->getElementsCount();
//	for(int i = 0; i < ect; i++){
//		MeshQuad2d* quad = (MeshQuad2d*)mesh->getElementAt(i);
//		int edge_ct = quad->getEdgeCount();
//		if(edge_ct != 4) return false;
//		int special_count = 0;
//		for(int j = 0; j < edge_ct; j++){
//			MeshEdge2d* edge = quad->getEdge(j);
//			if(edge->getBorderType() == 1) ++special_count;
//		}
//		if(special_count > 1) return false;
//	}
//	return true;
//}
//
double MeshSpecialRoutinesDAT::getMaxSecantError(MeshContainer2d *dat_mesh)
{
	if(!dat_mesh) return 1.0;

	double max_ratio = 0.0;

	int i, j, qct = dat_mesh->getElementsCount();
	for(i = 0; i < qct; i++){
		MeshQuad2d* quad = (MeshQuad2d*)dat_mesh->getElementAt(i);
		Quad12 quad12(quad);

		DPoint2d params[] = {DPoint2d(-1.0, -1.0), DPoint2d( 1.0, -1.0), DPoint2d( 1.0,  1.0), DPoint2d(-1.0,  1.0), 
			DPoint2d(-1/3.0, -1.0), DPoint2d( 1/3.0, -1.0), DPoint2d( 1/3.0,  1.0), DPoint2d(-1/3.0,  1.0),
			DPoint2d(-1.0, -1/3.0), DPoint2d( 1.0, -1/3.0), DPoint2d( 1.0,  1/3.0), DPoint2d(-1.0,  1/3.0)};

		DPoint3d points[12];

		for(j = 0; j < 12; j++){
			DPoint2d pt2d = quad12.getPoint(params[j]);
			points[j].x = pt2d.x;
			points[j].y = pt2d.y;
			points[j].z = quad12.getValue(params[j]);
		}

		int pairs[12][2] = {{0,4},{4,5},{5,1},{0,8},{8,11},{11,3},{1,9},{9,10},{10,2},{3,7},{7,6},{6,2}};

		double local_ratio = 0.0;
		for(j = 0; j < 12; j++){
			int i0 = pairs[j][0];
			int i1 = pairs[j][1];
			DPoint2d middle_param = DPoint2d::average(params[i0], params[i1]);
			DPoint2d pt2d = quad12.getPoint(middle_param);
			DPoint3d val_middle(pt2d.x, pt2d.y, quad12.getValue(middle_param));
			DPoint3d ave_middle = DPoint3d::average(points[i0], points[i1]);
			double ratio = val_middle.distance(ave_middle) / points[i0].distance(points[i1]);
			if(ratio > local_ratio) local_ratio = ratio;
		}
		quad->setQuality(local_ratio);

		if(local_ratio > max_ratio) max_ratio = local_ratio;
	}

	return max_ratio;
}

double MeshSpecialRoutinesDAT::getMaxOscillationError(MeshContainer2d *dat_mesh)
{
	if(!dat_mesh) return 0.0;

	double max_error = 0.0;

	int i, j, qct = dat_mesh->getElementsCount();
	for(i = 0; i < qct; i++){
		MeshQuad2d* quad = (MeshQuad2d*)dat_mesh->getElementAt(i);
		Quad12 quad12(quad);

		DPoint2d params[] = {DPoint2d(-1.0, -1.0), DPoint2d( 1.0, -1.0), DPoint2d( 1.0,  1.0), DPoint2d(-1.0,  1.0), 
			DPoint2d(-1/3.0, -1.0), DPoint2d( 1/3.0, -1.0), DPoint2d( 1/3.0,  1.0), DPoint2d(-1/3.0,  1.0),
			DPoint2d(-1.0, -1/3.0), DPoint2d( 1.0, -1/3.0), DPoint2d( 1.0,  1/3.0), DPoint2d(-1.0,  1/3.0)};

		double values[12];

		for(j = 0; j < 12; j++){
			values[j] = quad12.getValue(params[j]);
		}

		int edges[4][4] = {{0,4,5,1},{1,9,10,2},{2,6,7,3},{3,11,8,0}};

		double local_error = 0.0;
		for(j = 0; j < 4; j++){
			double dv1 = values[edges[j][1]] - 0.5*(values[edges[j][0]] + values[edges[j][2]]);
			double dv2 = values[edges[j][2]] - 0.5*(values[edges[j][1]] + values[edges[j][3]]);
			if(dv1*dv2 < -mesh_data.relative_small_number){
				dv1 /= abs(values[edges[j][0]] - values[edges[j][2]]);
				dv2 /= abs(values[edges[j][1]] - values[edges[j][3]]);
				double dv = abs(dv1-dv2);
///				LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, endl << " -> " << dv1 << " * " << dv2 << " ** " << dv;
				if(dv > local_error) local_error = dv;
			}
		}
		quad->setQuality(local_error);

		if(local_error > max_error) max_error = local_error;
	}

//	if(max_error > 0.0) LOG4CPLUS_INFO(MeshLog::logger_mesh, endl);

	return max_error;
}

double MeshSpecialRoutinesDAT::getMaxDiffError(MeshContainer2d *dat_mesh, double tinit)
{
	if(!dat_mesh) return 1.0;

	double max_diff = 0.0;

	int i, j, qct = dat_mesh->getElementsCount();
	for(i = 0; i < qct; i++){
		MeshQuad2d* quad = (MeshQuad2d*)dat_mesh->getElementAt(i);
		Quad12 quad12(quad);

		DPoint2d params[] = {DPoint2d(-1.0, -1.0), DPoint2d( 1.0, -1.0), DPoint2d( 1.0,  1.0), DPoint2d(-1.0,  1.0), 
			DPoint2d(-1/3.0, -1.0), DPoint2d( 1/3.0, -1.0), DPoint2d( 1/3.0,  1.0), DPoint2d(-1/3.0,  1.0),
			DPoint2d(-1.0, -1/3.0), DPoint2d( 1.0, -1/3.0), DPoint2d( 1.0,  1/3.0), DPoint2d(-1.0,  1/3.0)};

		double local_diff = 0.0;

		for(j = 0; j < 12; j++){
			double diff = tinit - quad12.getValue(params[j]);
			if(diff > local_diff) local_diff = diff;
		}

		int pairs[12][2] = {{0,4},{4,5},{5,1},{0,8},{8,11},{11,3},{1,9},{9,10},{10,2},{3,7},{7,6},{6,2}};

		for(j = 0; j < 12; j++){
			int i0 = pairs[j][0];
			int i1 = pairs[j][1];
			DPoint2d middle_param = DPoint2d::average(params[i0], params[i1]);
			double diff = tinit - quad12.getValue(middle_param);
			if(diff > local_diff) local_diff = diff;
		}
		quad->setQuality(local_diff);

		if(local_diff > max_diff) max_diff = local_diff;
	}

	return max_diff;
}

MeshContainer3d* MeshSpecialRoutinesDAT::loadSurfaceBoundaryMesh(const string& fname)
{
	ifstream file_p((fname+"-p.txt").c_str());
	if(!file_p) return nullptr;
	ifstream file_i((fname+"-i.txt").c_str());
	if(!file_i)	return nullptr;

	// points
	int p_count;
	file_p >> p_count;

	if(p_count < 3) return nullptr;

	MeshContainer3dSurface* surface_mesh = new MeshContainer3dSurface(p_count);

	DataHashTableKeyValue<int,int> PIds_ref((unsigned int)p_count, -100);
	DBox bbox;
	for(int i = 0; i < p_count; i++){
		int id;
		double x,y,z;
		file_p >> id >> x >> y >> z;
		MeshPoint3d* point = new MeshPoint3d(x, y, z);
		surface_mesh->addMeshPoint(point);
		PIds_ref.insert(id, point->getIndex());
		bbox.addPoint(point->getCoordinates());
	}

	// faces
	int f_count;
	file_i >> f_count;

	string buffer;
	getline(file_i, buffer);
	for(int i = 0; i < f_count; i++){
		getline(file_i, buffer);
		istringstream sfile_i(buffer);
		int id, ect;
		int PIds[4];
		sfile_i >> id >> ect;
		assert(ect == 3);
		for(int j = 0; j < ect; j++)
			sfile_i >> PIds[j];
		MeshFace* face = nullptr;
		if(ect == 3){
			face = new MeshTriangle3d(
				surface_mesh->getPointAt(PIds_ref.getValue(PIds[1], 0)),
				surface_mesh->getPointAt(PIds_ref.getValue(PIds[0], 0)),
				surface_mesh->getPointAt(PIds_ref.getValue(PIds[2], 0)));
			surface_mesh->addMeshFace(face);
		}
	}

	MeshDomainVolume* domain_volume = new MeshDomainVolume(surface_mesh);
	domain_volume->setAreaID(0);

	CS3dPtr space(new ControlSpace3dOctree(bbox));
	if(((ControlSpace3dOctree*)(space.get()))->loadTXT((fname+"-cs.txt").c_str())){
		domain_volume->setUserControlSpace(space);
		domain_volume->createInitialControlSpace();
	}else{
		START_CLOCK("MG3dSpecial::PrepareCS");
		domain_volume->createInitialControlSpace();
		STOP_CLOCK("MG3dSpecial::PrepareCS");
	}

	MeshContainer3d *domain = new MeshContainer3d(10);
	domain->addMeshBlock(domain_volume);

	return domain;
}

MeshContainer3d* MeshSpecialRoutinesDAT::loadSpecialGrid(const char *fname)
{
	ifstream file(fname);
	if(!file) return nullptr;

	// Wspó³rzêdne punktów
	int pct, tct;
	file >> pct >> tct;

	if(pct < 3) return nullptr;

	MeshContainer3dSurface* surface_mesh = new MeshContainer3dSurface(pct);

	for(int i = 0; i < pct; i++){
		DPoint3d pt;
		file >> pt.x >> pt.y >> pt.z;
		MeshPoint3d *mpt = new MeshPoint3d(pt);
		surface_mesh->addMeshPoint(mpt);
	}

	// Tablica incydencji
	for(int i = 0; i < tct; i++){
		int i0, i1, i2;
		file >> i0 >> i1 >> i2;
		MeshFace* face = new MeshTriangle3d(
			surface_mesh->getPointAt(i0-1),
			surface_mesh->getPointAt(i1-1),
			surface_mesh->getPointAt(i2-1));
		surface_mesh->addMeshFace( face );
	}


	MeshDomainVolume* domain_volume = new MeshDomainVolume(surface_mesh);
	domain_volume->setAreaID(0);

	MeshContainer3d *domain = new MeshContainer3d(10);
	domain->addMeshBlock(domain_volume);
	
	return domain;
}

bool MeshSpecialRoutinesDAT::runSpecialGeometryPhase(int /* hesjan_control */, int /* from_step */)
{
	/*

	CFileDialog open_file_run(true, nullptr, nullptr, OFN_FILEMUSTEXIST | OFN_HIDEREADONLY, "Program symulacji (*.exe)|*.exe|Wszystkie pliki (*.*)|*.*||");
	if(open_file_run.DoModal()!=IDOK) return false;

	CFileDialog open_file(true, nullptr, nullptr, OFN_FILEMUSTEXIST | OFN_HIDEREADONLY, "Dane do procesu symulacji (*.sin)|*.sin|Wszystkie pliki (*.*)|*.*||");
	if(open_file.DoModal()!=IDOK) return false;
	FILE* file = fopen(open_file.GetPathName(), "rt");
	if(!file) return false;

	double a, dl, wspdyf, tferryt, tinit, c1, dr1, dr2, dl1, dl2, sim_time, l_ratio, last_radius;
	double d_surf, prev_r, last_sim_time, prev_c1, adaptation_error, time_adapted_ratio;
	int qmorph, smooth_t, smooth_q, max_steps, min_steps, max_rerun, shape_type;
	int isothermal, control_time_step, izad = 10, check_adaptation_error;
	int max_quads;
	double aaa,bbb, crate, tempcells;
	char mesh_name[100];

	fscanf(file, "isothermal=%d\n", &isothermal);
	fscanf(file, "aaa=%lf\n", &aaa);
	fscanf(file, "bbb=%lf\n", &bbb);
	fscanf(file, "crate=%lf\n", &crate);
	fscanf(file, "tempcells=%lf\n", &tempcells);
	fscanf(file, "time=%lf\n", &sim_time);
	fscanf(file, "max_steps=%d\n", &max_steps);
	fscanf(file, "wspdyf=%lf\n", &wspdyf);
	fscanf(file, "tferryt=%lf\n", &tferryt);
	fscanf(file, "tinit=%lf\n", &tinit);
	fscanf(file, "c1=%lf\n", &c1);
	fscanf(file, "a=%lf\n", &a);
	fscanf(file, "dl=%lf\n", &dl);
	fscanf(file, "nazwa=%s\n", mesh_name);
	fscanf(file, "layers_ratio=%lf\n", &l_ratio);
	fscanf(file, "last_radius=%lf\n", &last_radius);
	fscanf(file, "min_steps=%d\n", &min_steps);
	min_steps--;
	fscanf(file, "check_adaptation_error=%d\n", &check_adaptation_error);
	fscanf(file, "max_adaptation_error=%lf\n", &adaptation_error);
	fscanf(file, "max_adaptation_reruns=%d\n", &max_rerun);
	fscanf(file, "dr1=%lf\n", &dr1);
	fscanf(file, "dr2=%lf\n", &dr2);
	fscanf(file, "dl1=%lf\n", &dl1);
	fscanf(file, "dl2=%lf\n", &dl2);
	fscanf(file, "shape=%d\n", &shape_type);
	fscanf(file, "qmorph=%d\n", &qmorph);
	fscanf(file, "smooth_t=%d\n", &smooth_t);
	fscanf(file, "smooth_q=%d\n", &smooth_q);
	fscanf(file, "time_adapted_ratio=%lf\n", &time_adapted_ratio);
	fscanf(file, "control_time_step=%d\n", &control_time_step);
	if(fscanf(file, "max_quads=%d\n", &max_quads) < 1) max_quads = 1000;
	fclose(file);

	if(qmorph < 0) return false;

	double r = dl*a / (0.3*PI);
	last_radius *= a;

	int run_steps = -1;
	int last_solver_step = 0;
	int prev_last_solver_step;
	double total_time = 0.0;
	int qct;
	int rerun_count = 1;
	int prev_rerun_count;
	double prev_error = 0;

	if(from_step > 0){
		run_steps = from_step;
		CString file_name_din;
		file_name_din.Format("%s-%d-1p.din", mesh_name, run_steps);
		FILE* file_din = fopen(file_name_din, "rt");
		if(!file_din){
			CString message;
			message.Format("No %s generated !", file_name_din);
			AfxMessageBox(message, MB_OK); 
			return false;
		}
		fscanf(file_din, " dpole=%lf\n", &d_surf);
		fscanf(file_din, " promien=%lf\n", &prev_r);
		fscanf(file_din, " war.brzeg=%lf\n", &c1);
		//???
		prev_c1 = c1;
		fscanf(file_din, " czas=%lf\n", &last_sim_time);
		fscanf(file_din, " krok=%d\n", &last_solver_step);
		fclose(file_din);

		r = MeshSpecialRoutinesDAT::countNextPhaseRadius(prev_r, a, d_surf, shape_type);

		prev_last_solver_step = last_solver_step;
		prev_rerun_count = 1;
		izad = 1;
		while((izad*=10) < last_solver_step);
		izad /= 10;
	}else{
		if(isothermal == 1 && c1 != tinit){
			LOG4CPLUS_INFO(MeshLog::logger_mesh, "WARNING! isothermal = 1, but c1 is different from tinit ?");
		}
	}

	LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, "ITER\t" << "MAX_ERROR\t" << "SIZ\t" << "STEPS\t" << "QUADS\t";
	LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, "LAST_RADIUS\t" << "TOT_TIME\t" << "DELTA_SURF\t" << "SIM_TIME";
	if(isothermal == 1){
		LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, "\tBOUND_COND";
	}
	LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, endl;

	while(last_radius > r){
		MeshContainer2d* value_mesh = nullptr;
		MeshContainer2d* control_mesh = nullptr;
		++run_steps;

//		if(++run_steps % 100 == 99){
//			if(AfxMessageBox("Reached 100 steps. Continue?", MB_YESNO) == IDNO) return false;
//		}

		if(m_boundary) delete m_boundary;

		bool mesh_not_ready = true;
		double temp_dl1 = dl1;
		double h_factor = DHesjan::param_hesjan_factor;
		double h_min_ratio = DHesjan::param_hesjan_min_ratio;
		double sim_adapted_time = sim_time;

		bool use_hesjan = hesjan_control > -1 && hesjan_control <= run_steps;
		char mark = use_hesjan?'h':'u';
		if(run_steps == 0) mark = 'u';

		if(run_steps > 0){
			CString fp_name, fi_name;
			fp_name.Format("%s-%d-%dp%d.dat", mesh_name, run_steps-1, rerun_count, last_solver_step);
			fi_name.Format("%s-%d-%di.dat", mesh_name, run_steps-1, rerun_count);

			string str_fp = fp_name;
			string str_fi = fi_name;
			control_mesh = value_mesh = MeshSpecialRoutinesDAT::importFromDAT(str_fp, str_fi);
			if(value_mesh == nullptr){
				delete m_boundary;
				AfxMessageBox("Failed loading previous-iteration dat files !", MB_OK); 
				return false;		
			}

			double max_error = 0;
			switch(check_adaptation_error){
			case 1:
				max_error = MeshSpecialRoutinesDAT::getMaxSecantError(value_mesh);
				break;
			case 2:
				max_error = MeshSpecialRoutinesDAT::getMaxOscillationError(value_mesh);
				break;
			case 3:
				max_error = MeshSpecialRoutinesDAT::getMaxDiffError(value_mesh, tinit);
				break;
			case 0:
			default:
				adaptation_error = max_error = 0.0;
				break;
			}

			if(max_error > adaptation_error && rerun_count < max_rerun && (run_steps > from_step+1)){
				for(int k = 0; k < rerun_count; k++){
					temp_dl1 *= 0.9;
					h_min_ratio *= 0.9;
					sim_adapted_time *= time_adapted_ratio;	// if less than 1.0, time-step will be shortened
				}

				--run_steps;

				LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, "* " << rerun_count << " * \t";
				MESHLOG.precision(3);
				LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, max_error;
				MESHLOG.precision(5);
				LOG4CPLUS_INFO(MeshLog::logger_mesh, " *** adaptation error too large -> recalculating ...");
				LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, run_steps << "\t";
				LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, prev_error << "\t";

				// one step back ...
				++rerun_count;
				r = prev_r;
				c1 = prev_c1;
				total_time -= last_sim_time;

				fp_name.Format("%s-%d-%dp%d.dat", mesh_name, run_steps-1, prev_rerun_count, prev_last_solver_step);
				fi_name.Format("%s-%d-%di.dat", mesh_name, run_steps-1, prev_rerun_count);

				string str_fp = fp_name;
				string str_fi = fi_name;
				value_mesh = MeshSpecialRoutinesDAT::importFromDAT(str_fp, str_fi);
			}else if(control_time_step && last_solver_step <= min_steps+1 && rerun_count < max_rerun){

				max_steps *= 10;
				LOG4CPLUS_INFO(MeshLog::logger_mesh, "* Time step decreased 10 times: " << sim_time << "s, " << max_steps << " steps.");

				--run_steps;

				LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, "* " << rerun_count << " * \t";
				LOG4CPLUS_INFO(MeshLog::logger_mesh, " *** too few simulation steps -> recalculating ...");
				LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, run_steps << "\t";
				LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, prev_error << "\t";

				// one step back ...
				++rerun_count;
				r = prev_r;
				c1 = prev_c1;
				total_time -= last_sim_time;

				fp_name.Format("%s-%d-%dp%d.dat", mesh_name, run_steps-1, prev_rerun_count, prev_last_solver_step);
				fi_name.Format("%s-%d-%di.dat", mesh_name, run_steps-1, prev_rerun_count);

				string str_fp = fp_name;
				string str_fi = fi_name;
				value_mesh = MeshSpecialRoutinesDAT::importFromDAT(str_fp, str_fi);
			}else{
				if(control_time_step && last_solver_step > 100){
					if(max_steps > 100) max_steps /= 10;
					else sim_time *= 10.0;
					LOG4CPLUS_INFO(MeshLog::logger_mesh, "* Time step increased 10 times: " << sim_time << "s, " << max_steps << " steps.");
				}
				LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, run_steps << "\t" << max_error << "\t";
				prev_rerun_count = rerun_count;
				prev_last_solver_step = last_solver_step;
				rerun_count = 1;
				prev_error = max_error;
			}
		}else{
			LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, "0\t0\t";
		}

		while(mesh_not_ready){
			LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, mark;
			m_boundary = MeshSpecialRoutinesDAT::createBoundaryForPhaseTransformation(a, r, a*dr1, a*dr2, a*temp_dl1, a*dl2, run_steps==0 || !use_hesjan, shape_type);
			if(!m_boundary) {
				AfxMessageBox("Failed creating phase-boundary !", MB_OK); 
				return false;	// eg. if r is too big
			}		

			if(run_steps > 0 && use_hesjan){
				MeshContainer2d* mesh_boundary = m_boundary->getFirst2dBoundary();
				MeshSpecialRoutinesDAT::setControlSpaceFromDATHesjan(mesh_boundary, control_mesh, h_min_ratio, h_factor);
			}

			if(rerun_count > 1){
				MeshContainer2d* mesh_boundary = m_boundary->getFirst2dBoundary();
				CS2dPtr control = mesh_boundary->getControlSpace();
				assert(control);
				int i, qct = control_mesh->getElementsCount();
				double factor = 0.7;
				for(i = 2; i < rerun_count; i++, factor*=0.7);
				for(i = 0; i < qct; i++){
					MeshQuad2d* quad = (MeshQuad2d*)control_mesh->getElementAt(i);
					if(quad->getQuality() > adaptation_error){
						DPoint2d middle = quad->getMiddlePoint();
						double radius = 3.0*sqrt(quad->getAreaNoMetric() / PI);
						control->addFactor(middle, rerun_count*radius, factor);
					}
				}
			}

			mesh_not_ready = false;
			mesh_data.m_mesh_view_mode = 2;	// full-mesh
#ifdef USE_EXCEPTIONS
			try{
#endif // USE_EXCEPTIONS
				MeshGenerator1d::discretizeEdges(m_boundary);
				MeshGenerator2d::prepareBoundaryMesh(m_boundary);
				MeshGenerator2d::triangulateFaces(m_boundary);
				MeshGenerator2d::smoothenFaces(m_boundary, smooth_t);
				LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, 't';
				if(qmorph){
					MeshGenerator2dQuad::convertFacesToQuads(m_boundary, MeshData::QUADS_QMORPH);
					LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, 'q';
					int tct = m_boundary->getFirst2dMesh()->getElementsCount(3);
					if(tct > 1){
						MeshGenerator2dQuad::convertFacesToQuads(m_boundary, MeshData::QUADS_LEELO);
						LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, 'l';
					}
				}else{
					MeshGenerator2dQuad::convertFacesToQuads(m_boundary, MeshData::QUADS_LEELO);
					LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, 'l';
				}
				MeshGenerator2dQuad::smoothenFaces(m_boundary, smooth_q);
				LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, 's';

				MeshContainer2d* first_mesh = m_boundary->getFirst2dMesh();

				LOG4CPLUS_ASSERT(MeshLog::logger_mesh, first_mesh != nullptr && first_mesh->isValid(), 
					MeshingException("GDoc:invalid first mesh"), false);

				qct = first_mesh->getElementsCount();
				if(qct > max_quads){
					mesh_not_ready = true;
					delete m_boundary;
					temp_dl1 *= 1.3;
					h_factor *= 1.3;
					h_min_ratio *= 1.3;
				}else if(!MeshSpecialRoutinesDAT::checkPhaseCorrectness(first_mesh)){
					for(int i = 0; i < 30; i++){
						MeshGenerator2dQuad::smoothenFaces(m_boundary);
						if(MeshSpecialRoutinesDAT::checkPhaseCorrectness(first_mesh))
							break;
					}
					if(!MeshSpecialRoutinesDAT::checkPhaseCorrectness(first_mesh)){
						mesh_not_ready = true;
						delete m_boundary;
						temp_dl1 *= 0.98;
						h_factor *= 0.98;
						h_min_ratio *= 0.98;
					}
				}

				if(!mesh_not_ready){
					createInnerNodes(2, false);
					int ect = first_mesh->getElementsCount();
					for(int i = 0; i < ect; i++){
						MeshElement* element = first_mesh->getElementAt(i);
						int edge_ct = element->getEdgeCount();
						for(int j = 0; j < edge_ct; j++){
							if(element->getEdge(j)->getInnerPointsCount() != 2){
								mesh_not_ready = true;
								break;
							}
						}
						if(mesh_not_ready){
							delete m_boundary;
							temp_dl1 *= 0.98;
							h_factor *= 0.98;
							h_min_ratio *= 0.98;
							break;
						}
					}
				}
#ifdef USE_EXCEPTIONS				
			}catch(MeshingException){
				LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, '*';
				delete m_boundary;
				mesh_not_ready = true;
				temp_dl1 *= 0.99;
				h_factor *= 0.99;
				h_min_ratio *= 0.99;
			}
#endif // USE_EXCEPTIONS
		}
		if(run_steps == 0){
			setPointsWeight(tinit);
		}else{
			MeshContainer2d* mesh = m_boundary->getFirst2dMesh();
			assert(mesh);
			assert(mesh->isValid()); /// ttttttt

			MeshSpecialRoutinesDAT::interpolateValues(mesh, value_mesh);
			assert(mesh->isValid()); /// ttttttt
			if(control_mesh != value_mesh) delete control_mesh;
			delete value_mesh;
		}
		CString file_name;
		file_name.Format("%s-%d-%d", mesh_name, run_steps, rerun_count);
		storeDAT(file_name);
		assert(m_boundary->getFirst2dMesh()->isValid()); /// ttttttt

		file = fopen("temp.dat", "wt");
		if(file){
			fprintf(file, "0 0 8 \t igrain,nequat,mater\n");
			fprintf(file, "8 6 %d \t transf.model, imesh (3-rod, 5-rail,6-nonstr), izad\n", izad);
			if(isothermal == -1){	// old and simplified
				fprintf(file, "%G %G \t wspdyf, tferryt\n", wspdyf, tferryt);
			}else{
				fprintf(file, "%G \t tferryt\n", tferryt);
			}
			fprintf(file, "1.0 1 0 \t AX, IDANE, NREAD\n");
			fprintf(file, "319 2\n");
			fprintf(file, "-1  1  1  1 \t lewy, prawy, dol, gora\n");
			fprintf(file, "10.0  10.00\n");
			fprintf(file, "%G \t tinit\n", tinit);
			fprintf(file, "%G %d %G \t lratio min_steps promien\n", l_ratio, min_steps, r);
			fprintf(file,"\n0 \t brak punktow do wykresow\n\n");
			fprintf(file,"1 \t liczba etapow chlodzenia\n");
			fprintf(file,"%G %G %d \t TOTOCZ, CZAS [S], LICZBA KROKOW w etaPIe\n", c1, sim_adapted_time, max_steps);
			if(isothermal > -1){
				fprintf(file,"%G %G \t aaa, bbb\n", aaa, bbb);
				fprintf(file,"%d \t isothermal\n", isothermal);
				if(isothermal == 0){
					fprintf(file,"%G \t tempcells\n", tempcells);
				}else{
					fprintf(file,"%G \t crate\n", crate);
				}
			}
			fprintf(file,"2\n1. -1. 1. 1. 1\n");
			fclose(file);
		}

		file = fopen("dane.in", "wt");
		if(file){
			fprintf(file, "%s-%d-%di.dat\n", mesh_name, run_steps, rerun_count);
			fprintf(file, "%s-%d-%dp.dat\n", mesh_name, run_steps, rerun_count);
			fprintf(file, "0.0\n");
			fclose(file);
		}

		assert(m_boundary->getFirst2dMesh()->isValid()); /// ttttttt

		CString run_file_name;
		run_file_name.Format("%s < dane.in > %s-%d-%d.log", open_file_run.GetPathName(), mesh_name, run_steps, rerun_count);

		if(system(run_file_name) != 0){
			if(AfxMessageBox("Called programm returned non-zero. Continue?", MB_YESNO) == IDNO) 
				return false;
		}

		run_file_name.Format("%s-%d-%d.temp.dat", mesh_name, run_steps, rerun_count);
		rename("temp.dat", run_file_name);

/////////////////////////////////////////
		prev_c1 = c1;
		CString file_name_din;
		file_name_din.Format("%s-%d-%dp.din", mesh_name, run_steps, rerun_count);
		FILE* file_din = fopen(file_name_din, "rt");
		if(!file_din){
			AfxMessageBox("No xxx-xp.din generated !", MB_OK); 
			return false;
		}
		fscanf(file_din, " dpole=%lf\n", &d_surf);
		fscanf(file_din, " promien=%lf\n", &prev_r);
		fscanf(file_din, " war.brzeg=%lf\n", &c1);
		fscanf(file_din, " czas=%lf\n", &last_sim_time);
		fscanf(file_din, " krok=%d\n", &last_solver_step);
		fclose(file_din);

		total_time += last_sim_time;

		LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, "\t" << last_solver_step << "\t" << qct << "\t";
		LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, prev_r << "\t" << total_time << "\t";
		LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, d_surf << "\t" << last_sim_time;
		if(isothermal == 1){
			LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, "\t" << c1;
		}
		LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, endl;

		izad = 1;
		while((izad*=10) < last_solver_step);
		izad /= 10;

		r = MeshSpecialRoutinesDAT::countNextPhaseRadius(prev_r, a, d_surf, shape_type);
	}

	UpdateAllViews(nullptr);	

	*/
	return false;
}

bool MeshSpecialRoutinesDAT::testHessianCurvature(const string& fname, const string& equation, const string& equation2, bool with_quads)
{
/*
	MeshContainer3d* boundary = MeshStream::readFileMsh(fname.c_str());

//	ControlSpace2d::param_max_diameter_ratio = 0.02;

	MeshGenerator1d::discretizeEdgesMin(boundary);
	MeshGenerator2d::prepareBoundaryMesh(boundary);
	MeshGenerator2d::triangulateFaces(boundary);
	MeshGenerator2d::smoothenFaces(boundary, 2);
	if(with_quads){
		MeshGenerator2dQuad::convertFacesToQuads(boundary, MeshData::QUADS_LEELO);
		MeshGenerator2dQuad::smoothenFaces(boundary, 3, MeshData::SM_TOP_SWAP | MeshData::SM_LAPLACE_MIXED);
	}

	ControlSpace2dMesh::param_interpolation_method = MeshData::CONTROL_MIXED;

	CS2dPtr last_control_space = nullptr;

	char code = 'a';
	string answer;

	DEquation feq(equation);
	DEquation *feq2 = nullptr;
	if(equation2 != "")
		feq2 = new DEquation(equation2);

	while(true){
		MeshView::setViewMode(SHOW_MESH_SURF);
		if(MeshView::isRunning())
			MeshView::setBoundaryToUpdate(boundary);
		else
#ifdef UNIX_MODE
 #ifdef NO_OPENGL
		{
			LOG4CPLUS_WARN(MeshLog::logger_console, "ComPIled without OpenGL support!");
		}
 #else
		{
			pthread_t t;
			pthread_create(&t, nullptr, MeshView::startLoopGL, boundary);
		}
 #endif
#else
			_beginthread(MeshView::startLoopGL, 0, boundary);
#endif
		IteratorMesh2d it = boundary->getFirstValidMesh2d();
		if(!it.isValid()) return false;
		MeshContainer2d* mesh = it.getMesh();
		int pct = mesh->getPointsCount();

		mesh->storeEPS("mesh", code-'a');

		LOG4CPLUS_INFO(MeshLog::logger_console, "nodes", pct);
		LOG4CPLUS_INFO(MeshLog::logger_console, "triangles", mesh->getElementsCount());

		cout << "Next? (y/n)";
		cin >> answer;
		if(answer != "y") break;

		ControlSpace2dAdaptive::param_max_diameter_ratio = 0.05;

		DPoint3d curvature;

		LOG4CPLUS_ERROR(MeshLog::logger_console,   "TODO - points weight for DAT");
		for(int i = 0; i < pct; i++){
			MeshPoint2d* mpt = mesh->getPointAt(i);
			const DPoint2d coord = mpt->getCoordinates();
			LOG4CPLUS_ERROR(MeshLog::logger_console,   "TODO - points weight for DAT");
		}

		DRect rect = mesh->getBoundingRect();
		rect.inflate(ControlSpace2d::param_inflate_box_factor);
		CS2dPtr space = new ControlSpace2dMesh(mesh->getSurface(), rect, 
								ControlSpace2dMesh::param_interpolation_method);
	//	CS2dPtr space = new ControlSpace2dMatrixUniform(mesh->getSurface(), rect, 100, 100);


		double p = mesh_data.getModelDiameter();
		double max_len = p * ControlSpace2dAdaptive::param_max_diameter_ratio;
		double min_len = std::max(p * ControlSpace2dAdaptive::param_min_diameter_ratio, ControlSpace2dAdaptive::param_min_diameter_ratio);

		ControlDataStretch2d *control_data = new ControlDataStretch2d[pct];
		short *hessian_valid = new short[pct];

		for(int i = 0; i < pct; i++){
			hessian_valid[i] = 0;
			MeshPoint2d* mpt = mesh->getPointAt(i);
			if(mpt->isBorder()) continue;
			//DHesjan dh1 = mpt->countHesjan(9, false);
			//DHesjan dh2 = mpt->countHesjan(9, true);
			if(mpt->countHesjanCurvature(9, curvature)){
				control_data[i] = ControlSpace2dAdaptive::adjustCurvatureData(curvature, 
					ControlSpace2dAdaptive::param_curvature_ratio, p, min_len, max_len, 10);
//				space->addControlPoint(mpt->getCoordinates(), control_data[i]);
				hessian_valid[i] = 2;
			}else{
				LOG4CPLUS_INFO(MeshLog::logger_console, "Hessian counting failed.");
			}
		}
		if(feq2){
			LOG4CPLUS_ERROR(MeshLog::logger_console,   "TODO - points weight for DAT");
			for(int i = 0; i < pct; i++){
				MeshPoint2d* mpt = mesh->getPointAt(i);
				const DPoint2d coord = mpt->getCoordinates();
				mpt->setWeight(feq2->getValue(coord.x, coord.y));
			}

			for(int i = 0; i < pct; i++){
				MeshPoint2d* mpt = mesh->getPointAt(i);
				if(mpt->isBorder()) continue;
				if(mpt->countHesjanCurvature(9, curvature)){
					ControlDataStretch2d data_s = ControlSpace2dAdaptive::adjustCurvatureData(curvature, 
						ControlSpace2dAdaptive::param_curvature_ratio, p, min_len, max_len, 10);
					if(hessian_valid[i] == 0){
						control_data[i] = data_s;
						hessian_valid[i] = 1;
					}else{
						ControlDataMatrix2d dm1 = DMetric2d::stretchToMatrix(control_data[i]);
						ControlDataMatrix2d dm2 = DMetric2d::stretchToMatrix(data_s);
						dm1.setMinimum(dm2);
						control_data[i] = DMetric2d::matrixToStretch(dm1);
					}
				}
			}
		}
	
		for(int k = 2; k >= 1; k--){
			for(int i = 0; i < pct; i++){
				if(hessian_valid[i] == 0){
					MeshPoint2d* mpt = mesh->getPointAt(i);
					int rank = mpt->getRank();
					ControlDataMatrix2d dm;
					double w = 0.0;
					for(int j = 0; j < rank; j++){
						MeshPoint2d* other_pt = mpt->getEdge(j)->getOtherPoint(mpt);
						int idx = other_pt->getIndex();
						if(hessian_valid[idx] >= k){
							double wx = 1.0 / mpt->getCoordinates().distance2(other_pt->getCoordinates());
							w += wx;
							dm += DMetric2d::stretchToDxx(control_data[idx]) * wx;
						}
					}
					if(w > 0.0){
						control_data[i] = DMetric2d::dxxToStretch( dm / w);
						space->addControlPoint(mpt->getCoordinates(), control_data[i]);
						hessian_valid[i] = 1;
					}
				}
			}
		}
	
		for(int i = 0; i < pct; i++){
			if(hessian_valid[i] == 0) continue;
			MeshPoint2d* mpt = mesh->getPointAt(i);
			space->addControlPoint(mpt->getCoordinates(), DMetric2d::stretchToMatrix(control_data[i]));
		}
		space->interpolate();
		delete[] hessian_valid;
		delete[] control_data;

		space->storeEPS();
//		((ControlSpace2dMesh*)space)->storeVisualisationInfo();

		if(last_control_space){
			MeshData::StatData stats = space->compareMetricWith(last_control_space);
			LOG4CPLUS_INFO(MeshLog::logger_console, "Metric difference min = ", stats.minimum);
			LOG4CPLUS_INFO(MeshLog::logger_console, "Metric difference max = ", stats.maximum);
			LOG4CPLUS_INFO(MeshLog::logger_console, "Metric difference ave = ", stats.average);
			if(last_control_space->decRef()) delete last_control_space;
		}
		space->incRef();
		last_control_space = space;

		MeshContainer3d* boundary_adapted = MeshStream::readFileMsh(fname.c_str());
		MeshDomainVolume* domain_volume = (MeshDomainVolume*)(boundary_adapted->getBlockAt(0));
		MeshDomainSurface* domain_surface = (MeshDomainSurface*)(domain_volume->getFace(0));
	//	domain_surface->setUserControlSpace(space);
		domain_surface->getBoundary()->setControlSpace(space);

		MeshGenerator1d::discretizeEdgesMin(boundary_adapted);
		MeshGenerator2d::prepareBoundaryMesh(boundary_adapted);
		MeshGenerator2d::triangulateFaces(boundary_adapted);
		MeshGenerator2d::smoothenFaces(boundary_adapted, 2);
		if(with_quads){
			MeshGenerator2dQuad::convertFacesToQuads(boundary_adapted, MeshData::QUADS_LEELO);
			MeshGenerator2dQuad::smoothenFaces(boundary_adapted, 3, 
				MeshData::SM_TOP_SWAP | MeshData::SM_LAPLACE_MIXED);
		}

		MeshView::setViewMode(SHOW_MESH_SURF);
		if(MeshView::isRunning())
			MeshView::setBoundaryToUpdate(boundary_adapted);
		else
#ifdef UNIX_MODE
 #ifdef NO_OPENGL
		{
			LOG4CPLUS_WARN(MeshLog::logger_console, "ComPIled without OpenGL support!");
		}
 #else
		{
			pthread_t t;
			pthread_create(&t, nullptr, MeshView::startLoopGL, boundary_adapted);
		}
 #endif
#else
			_beginthread(MeshView::startLoopGL, 0, boundary_adapted);
#endif

		string str_name = "fmesh-";
		str_name += code++; 
		ofstream i_file((str_name+"-i.txt").c_str());
		ofstream p_file((str_name+"-p.txt").c_str());
		it = boundary_adapted->getFirstValidMesh2d();
		if(!it.isValid()) return false;
		MeshContainer2d* mesh_adapted = it.getMesh();

		pct = mesh_adapted->getPointsCount();
		p_file << pct << endl;
		for(int i = 0; i < pct; i++){
			MeshPoint2d* mpt = mesh_adapted->getPointAt(i);
			DPoint2d pt = mpt->getCoordinates();
			if(feq2){
				double f1 = feq.getValue(pt.x, pt.y);
				double f2 = feq2->getValue(pt.x, pt.y);
				double fx = (abs(f1) > abs(f2))?f1:f2;
				p_file << i << '\t' << pt.x << '\t' << pt.y << '\t' << fx << endl;
			}else{
				p_file << i << '\t' << pt.x << '\t' << pt.y << '\t' << feq.getValue(pt.x, pt.y) << endl;
			}
		}
		int ect = mesh_adapted->getElementsCount();
		i_file << ect << endl;
		for(int i = 0; i < ect; i++){
			MeshElement* element = mesh_adapted->getElementAt(i);
			int e_ct = element->getEdgeCount();
			i_file << i << '\t' << e_ct;
			for(int j = 0; j < e_ct; j++) i_file << '\t' << element->getPoint(j)->getIndex();
			i_file << "\t1" << endl;
		}

		delete boundary;
		boundary = boundary_adapted;
	}

	if(feq2) delete feq2;
	delete boundary;
*/
	return true;
}
