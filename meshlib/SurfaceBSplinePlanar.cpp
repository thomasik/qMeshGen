// SurfaceBSplinePlanar.cpp: implementation of the SurfaceBSplinePlanar class.
// Tomasz Jurczyk, 2009-
// Generation of unstructured meshes
//  * opened bspline surface with knots based on plane !
//
//////////////////////////////////////////////////////////////////////

#include "SurfaceBSplinePlanar.h"
#include "DataVector.h"
#include "DataMatrix.h"
#include "DLeastSquaresFitting.h"
#include "DPoint.h"
#include "DVector.h"
#include "MeshViewSet.h"
#include "DMatrix.h"
#include "DQuadric.h"
#include "Curve2dBSpline.h"
#include "DEquation.h"
#include "SurfaceMulti.h"
#include "DPlane.h"

#include "MeshGrain.h"
/// Standard constructor

//#define SHOW_BSURFACE_FIT

SurfaceBSplinePlanar::SurfaceBSplinePlanar(const SurfacePlane& plane, const DRect& rect) 
		: m_plane(  plane.getPoint(rect.getX0Y0()), 
					plane.getPoint(rect.getX1Y0()),
					plane.getPoint(rect.getX0Y1())),
		  m_nodes(nullptr), m_points(nullptr)
{
	m_valid = false;
}

SurfaceBSplinePlanar::SurfaceBSplinePlanar(const SurfacePlane& plane, DataMatrix<DPoint3d> *point_grid)
		: m_plane(plane), m_nodes(nullptr), m_points(nullptr)
{
	if(point_grid){
		int pct = point_grid->countInt();
		int rows = point_grid->rows();
		int cols = pct / rows;
		m_points = new DataMatrix<DPoint3d>(rows, cols, DPoint3d::zero);
		for(int i = 0; i < rows; i++)
			for(int j = 0; j < cols; j++)
				m_points->set(i, j, point_grid->get(i,j));
		m_valid = calculateNodes();
	}else
		m_valid = false;
}

SurfaceBSplinePlanar::~SurfaceBSplinePlanar()
{
	if(m_nodes) delete m_nodes;
	if(m_points) delete m_points;
}

const DPoint3d SurfaceBSplinePlanar::getPoint(const DPoint2d& param) const
{
	if(m_nodes){ 
		// calculate from B-spline equation
		int ui, vi;
		double u = param.x;
		double v = param.y;
		if(!countSplineIndices(u,v,ui,vi))
			return m_plane.getPoint(param);	// u/v -> ksi (0..1), ui/vi -> spline segment numbers

		double Nu[4], Nv[4];
		for(int i = 0; i < 4; i++){
			Nu[i] = Curve2dBSpline::BSplineMatrix[0][i];
			Nv[i] = Curve2dBSpline::BSplineMatrix[0][i];
			for(int j = 1; j < 4; j++){
				Nu[i] *= u;
				Nu[i] += Curve2dBSpline::BSplineMatrix[j][i];
				Nv[i] *= v;
				Nv[i] += Curve2dBSpline::BSplineMatrix[j][i];
			}
		}

		DPoint3d q;
		int n = m_nodes->rows() - 2;
		// sum Nu_i Nv_j B_i_j
		for(int uj = 0; uj < 4; uj++){
			for(int vj = 0; vj < 4; vj++){
				q.add(m_nodes->get(ui+uj, vi+vj), Nu[uj] * Nv[vj]);
			}
		}
		return q;
	}else{
		// use the plane directly
		return m_plane.getPoint(param);
	}
}

/// Returns the parameters of the surface for the given point with starting point
const DPoint2d SurfaceBSplinePlanar::getParametersNear(const DPoint3d& point, const DPoint2d& /* near_point */) const
{
	return m_plane.getParameters(point);
}

const DPoint2d SurfaceBSplinePlanar::getParameters(const DPoint3d& point) const
{
	return m_plane.getParameters(point);
}

const DVector3d SurfaceBSplinePlanar::getNormalVector(const DPoint2d& param) const
{
	if(m_nodes)
		return SurfaceParametric::getNormalVector(param);
	else
		return m_plane.getNormalVector(param);
}

const DVector3d SurfaceBSplinePlanar::getDerivative(int deriv, const DPoint2d& param) const
{
	if(m_nodes){ 
		// calculate from B-spline equation
		int ui, vi;
		double u = param.x;
		double v = param.y;
		if(!countSplineIndices(u,v,ui,vi))
			return m_plane.getDerivative(deriv, param);	// u/v -> ksi (0..1), ui/vi -> spline segment numbers

		double Nu[4], Nv[4];
		for(int i = 0; i < 4; i++){
			// -> Nu
			switch(deriv){
			case DEquation::deriv_dt:
			case DEquation::deriv_dtt: // f
				{
					Nu[i] = Curve2dBSpline::BSplineMatrix[0][i];
					for(int j = 1; j < 4; j++){
						Nu[i] *= u;
						Nu[i] += Curve2dBSpline::BSplineMatrix[j][i];
					}
					break;
				}
			case DEquation::deriv_ds:
			case DEquation::deriv_dst: // f'
				{
					Nu[i] = Curve2dBSpline::BDerivativeSplineMatrix[0][i];
					for(int j = 1; j < 3; j++){
						Nu[i] *= u;
						Nu[i] += Curve2dBSpline::BDerivativeSplineMatrix[j][i];
					}
					break;
				}
			case DEquation::deriv_dss: // f''
				{
					Nu[i] = Curve2dBSpline::BDDerivativeSplineMatrix[0][i];
					for(int j = 1; j < 2; j++){
						Nu[i] *= u;
						Nu[i] += Curve2dBSpline::BDDerivativeSplineMatrix[j][i];
					}
					break;
				}
			}
			// -> Nv
			switch(deriv){
			case DEquation::deriv_ds:
			case DEquation::deriv_dss: // f
				{
					Nv[i] = Curve2dBSpline::BSplineMatrix[0][i];
					for(int j = 1; j < 4; j++){
						Nv[i] *= v;
						Nv[i] += Curve2dBSpline::BSplineMatrix[j][i];
					}
					break;
				}
			case DEquation::deriv_dt:
			case DEquation::deriv_dst: // f'
				{
					Nv[i] = Curve2dBSpline::BDerivativeSplineMatrix[0][i];
					for(int j = 1; j < 3; j++){
						Nv[i] *= v;
						Nv[i] += Curve2dBSpline::BDerivativeSplineMatrix[j][i];
					}
					break;
				}
			case DEquation::deriv_dtt: // f''
				{
					Nv[i] = Curve2dBSpline::BDDerivativeSplineMatrix[0][i];
					for(int j = 1; j < 2; j++){
						Nv[i] *= v;
						Nv[i] += Curve2dBSpline::BDDerivativeSplineMatrix[j][i];
					}
					break;
				}
			}
		}

		DPoint3d q;
		int n = m_nodes->rows() - 2;
		// sum Nu_i Nv_j B_i_j
		for(int uj = 0; uj < 4; uj++){
			for(int vj = 0; vj < 4; vj++){
				q.add(m_nodes->get(ui+uj, vi+vj), Nu[uj] * Nv[vj]);
			}
		}
		return q.fixedVector();
	}else{
		// use the plane directly
		return m_plane.getDerivative(deriv, param);
	}
}

double SurfaceBSplinePlanar::fitAdditionally(const DataVector<DPoint3d> & points)
{
	if(!m_nodes || points.countInt() == 0) return 0.0;
	int level = m_nodes->rows();

	// ---
	if(false){
		MeshViewSet* set = getViewSet();
		for(int i = 0; i < points.countInt(); i++){
			set->addPoint(points[i], 1, i);
		}
		SHOW_MESH("BSpline::fitAdditionally - before", set);
	}
	DataMatrix<DPoint3d> prev_points(level, level, DPoint3d::zero);
	for(int i = 0; i < level; i++){
		for(int j = 0; j < level; j++){
			prev_points(i,j) = m_points->get(i,j);
		}
	}
	// ---


	// init base points from the cloud of points
	double dl = 1.0 / (level-3);

	DataMatrix<double> dh(level, level, 0.0); // -> bspline fitSquared -> (u,v,h)
	DataMatrix<int> w(level, level, 0); // -> bspline fitSquared -> (u,v,h)

	const DVector3d vn = m_plane.getNormalVector(DPoint2d::zero);
	DPoint3d pt = m_plane.getPoint(DPoint2d::zero);
	double d = vn.x * pt.x + vn.y * pt.y + vn.z * pt.z;
	for(int i = 0; i < points.countInt(); i++) {
		pt = points[i];
		DPoint2d pt_param = getParameters(pt);
		DPoint3d pts = getPoint(pt_param);
		double h_node = vn.x * pt.x  + vn.y * pt.y  + vn.z * pt.z  - d;
		double h_surf = vn.x * pts.x + vn.y * pts.y + vn.z * pts.z - d;
		double dh_pt = h_node - h_surf;
		int ix = 1 + (int) (pt_param.x / dl);
		int iy = 1 + (int) (pt_param.y / dl);

		if(ix < 0 || ix > level-2) continue;
		if(iy < 0 || iy > level-2) continue;

		dh(ix, iy) += dh_pt;
		w(ix, iy) += 1;
		dh(ix+1, iy) += dh_pt;
		w(ix+1, iy) += 1;
		dh(ix, iy+1) += dh_pt;
		w(ix, iy+1) += 1;
		dh(ix+1, iy+1) += dh_pt;
		w(ix+1, iy+1) += 1;
	}
	// apply
	for(int i = 1; i < level - 1; i++){
		for(int j = 1; j < level - 1; j++){
			if(w(i,j) == 0) continue;
			//
			m_points->set(i, j, m_points->get(i,j) + vn * (dh(i,j) / w(i,j)));
		}
	}
	// recalculate nodes
	if(!calculateNodes()) return -1.0;

	// ---
	if(false){
		MeshViewSet* set = getViewSet();
		for(int i = 0; i < points.countInt(); i++){
			set->addPoint(points[i], 1, i);
		}
		for(int i = 0; i < level; i++){
			for(int j = 0; j < level; j++){
				if(w(i,j) == 0) continue;
				set->addPoint(prev_points(i,j), 3, 0);
				set->addPoint(m_points->get(i,j), 3, 1);
			}
		}
		SHOW_MESH("BSpline::fitAdditionally - after", set);
	}
	// ---

	// check max dist
	double max_dist = -1.0;
	for(int i = 0; i < points.countInt(); i++) {
		pt = points[i];
		DPoint2d pt_param = getParameters(pt);
		DPoint3d pts = getPoint(pt_param);
		double h_node = vn.x * pt.x  + vn.y * pt.y  + vn.z * pt.z  - d;
		double h_surf = vn.x * pts.x + vn.y * pts.y + vn.z * pts.z - d;
		double dh_pt = abs(h_node - h_surf);
		if(max_dist < 0.0 || dh_pt > max_dist) max_dist = dh_pt;
	}
	return max_dist;
}

/// fit to point-cloud
double SurfaceBSplinePlanar::fitToPoints(const DataVector<DPoint3d> & points, int pct, int init_res)
{
	if(pct < FIT_PROBE_PTS) return -1.0; // or maybe create planar bspline?

	int level = init_res + 3; // size of grid + side control knots

	// create grid of knots x knots
	if(m_nodes) delete m_nodes;
	if(m_points) delete m_points;
	m_nodes  = new DataMatrix<DPoint3d>(level, level, DPoint3d());
	m_points = new DataMatrix<DPoint3d>(level, level, DPoint3d());

	// init base points from the cloud of points
	double dx = 1.0 / (level-3);
	double dy = 1.0 / (level-3);

#ifdef SHOW_BSURFACE_FIT
	DataMatrix<DataVector<int>> used_points(level, level, DataVector<int>(FIT_PROBE_PTS));
	DataMatrix<double> used_nearest(level, level, 0.0);
#endif

	DataMatrix<bool> ready_points(level, level, false);
	int not_ready_count = 0;

	const DVector3d vn = m_plane.getNormalVector(DPoint2d::zero);
	double max_local_dist2 = 0.7 * m_plane.getPoint(DPoint2d::zero).distance2(
		m_plane.getPoint(DPoint2d(dx, dy)));
	// for each knot (skipPIng border) ...
	for(int i = 1; i < level-1; i++){
		for(int j = 1; j < level-1; j++){
			DPoint2d knot((i-1) * dx, (j-1) * dy);
			DPoint3d knot3d = m_plane.getPoint(knot);
			DataVector<NodeDist> nearest(FIT_PROBE_PTS);
			// ... gather near cloud points
			double min_too_distant2 = -1.0;
			for(int k = 0; k < pct; k++){
				const DPoint3d plane_pt3d = m_plane.getPoint(m_plane.getParameters(points[k]));
				double dist2 = knot3d.distance2(plane_pt3d);
				if(dist2 > max_local_dist2) continue;
				if(nearest.countInt() < FIT_PROBE_PTS){
					nearest.add(NodeDist(k, dist2));
				}else{
					if(min_too_distant2 > 0.0 && dist2 > min_too_distant2) continue;
					int mmax = 0;
					for(int m = 1; m < nearest.countInt(); m++)
						if(nearest[m].dist2 > nearest[mmax].dist2) mmax = m;
					if(nearest[mmax].dist2 > dist2){
						min_too_distant2 = nearest[mmax].dist2;
						nearest[mmax] = NodeDist(k, dist2);
					}
				}
			}
			if(nearest.countInt() < FIT_PROBE_MIN){
				not_ready_count++;
				continue;
			}
			int knearest = 0;
			for(int k = 0; k < nearest.countInt(); k++){
				if(nearest[k].dist2 < nearest[knearest].dist2)
					knearest = k;
			}
			double dnearest = nearest[knearest].dist2;
#ifdef SHOW_BSURFACE_FIT
			DataVector<int> & used_list = used_points(i,j);
			used_nearest(i,j) = dnearest;
#endif
			if(dnearest < 0.1 * max_local_dist2){
				// count average height
				double h_ave = 0.0;
				size_t h_count = nearest.countInt();
				for(size_t k = 0; k < h_count; k++){
					const DPoint3d& real_pt3d = points[nearest[k].id];
					const DPoint2d plane_pt2d = m_plane.getParameters(real_pt3d);
					const DPoint3d plane_pt3d = m_plane.getPoint(plane_pt2d);
					const DVector3d dv = real_pt3d - plane_pt3d;
					if(dv.scalarProduct(vn) >= 0)
						h_ave += dv.length();
					else
						h_ave -= dv.length();
				}
				h_ave /= h_count;
				m_points->set(i, j, knot3d + vn * h_ave);
#ifdef SHOW_BSURFACE_FIT
				for(size_t k = 0; k < nearest.countInt(); k++){
					used_list.add(nearest[k].id);
				}
#endif
			}else{
				// replace nearest with all points from adjacent leaves ...
				DataVector<DPoint3d> selected_points(pct);
				for(int k = 0; k < pct; k++){
					const DPoint2d plane_pt2d = m_plane.getParameters(points[k]);
					const DPoint3d plane_pt3d = m_plane.getPoint(plane_pt2d);
					const DVector2d dpt = plane_pt2d - knot;
					if(abs(dpt.x) < dx && abs(dpt.y) < dy){
						selected_points.add(points[k]);
#ifdef SHOW_BSURFACE_FIT
						used_list.add(k);
#endif
					}
				}
				size_t h_count = selected_points.countInt();
				if(h_count >= FIT_PROBE_PTS){
					// ... fit plane
					DPlane local_plane;
					DLeastSquaresFitting::fitHyperplaneOrthogonal(
						selected_points, local_plane);
					// ... calculate value for (real) knot
					DPoint2d param;
					double z;
					if(local_plane.projectToPlane(knot3d, param, z)){
						m_points->set(i, j, knot3d + vn * z);
					}else{
						m_points->set(i, j, local_plane.projectToSpace(local_plane.projectToPlane(knot3d)));
					}
				}else if(h_count > 0){
					// count average height
					double h_ave = 0.0;
					for(size_t k = 0; k < h_count; k++){
						const DPoint3d& real_pt3d = selected_points[k];
						const DPoint2d plane_pt2d = m_plane.getParameters(real_pt3d);
						const DPoint3d plane_pt3d = m_plane.getPoint(plane_pt2d);
						const DVector3d dv = real_pt3d - plane_pt3d;
						if(dv.scalarProduct(vn) >= 0)
							h_ave += dv.length();
						else
							h_ave -= dv.length();
					}
					h_ave /= h_count;
					m_points->set(i, j, knot3d + vn * h_ave);
				}else{
					not_ready_count++;
					continue;
				}
			}
			// ... ready
			ready_points(i, j) = true;
		}
	}

	while(not_ready_count){
		const int off[4][2] = { { -1, -1}, {-1, 1}, {1, -1}, {1, 1} };
		bool any_change = false;
		for(int i = 1; i < level-1; i++){ // planar interpolation
			for(int j = 1; j < level-1; j++){
				if(ready_points(i, j)) continue;
				for(int k = 0; k < 4; k++){
					int ik = i + off[k][0];
					int jk = j + off[k][1];
					if( ik >= 0 && ik < level && jk >= 0 && jk < level &&
						ready_points(ik, j) && ready_points(i, jk) && ready_points(ik, jk))
					{
						m_points->set(i, j, m_points->get(i, jk) + 
							(m_points->get(ik, j) - m_points->get(ik, jk)));
						ready_points(i, j) = true;
						not_ready_count--;
						any_change = true;
						break; // the "k" loop
					}
				}
			}
		}
		if(any_change) continue;
		// .. else more straightforward solution
		for(int i = 1; i < level-1; i++){ 
			for(int j = 1; j < level-1; j++){
				if(ready_points(i, j)) continue;
				int ct = 0;
				DVector3d vt;
				if(i > 1 && ready_points(i-1, j)){
					DPoint3d knot3d = m_plane.getPoint(DPoint2d((i-1-1) * dx, (j-1) * dy));
					vt += m_points->get(i-1, j) - knot3d;
					++ct;
				}
				if(i < level-2 && ready_points(i+1, j)){
					DPoint3d knot3d = m_plane.getPoint(DPoint2d((i-1+1) * dx, (j-1) * dy));
					vt += m_points->get(i+1, j) - knot3d;
					++ct;
				}
				if(j > 1 && ready_points(i, j-1)){
					DPoint3d knot3d = m_plane.getPoint(DPoint2d((i-1) * dx, (j-1-1) * dy));
					vt += m_points->get(i, j-1) - knot3d;
					++ct;
				}
				if(j < level-2 && ready_points(i, j+1)){
					DPoint3d knot3d = m_plane.getPoint(DPoint2d((i-1) * dx, (j-1+1) * dy));
					vt += m_points->get(i, j+1) - knot3d;
					++ct;
				}
				if(ct > 0){
					vt /= ct;
					ready_points(i, j) = true;
					DPoint3d knot3d = m_plane.getPoint(DPoint2d((i-1) * dx, (j-1) * dy));
					m_points->set(i, j, knot3d + vt);
					not_ready_count--;
				}
			}
		}
	}

	// update the border knots
	fitToPointsAdjustBorder(level);

#ifdef SHOW_BSURFACE_FIT
	if(false){
		for(int ki = 1; ki < level-1; ki++){
			for(int kj = 1; kj < level-1; kj++){
				// -- cloud points
				DataVector<int> & used_list = used_points(ki, kj);
				if(used_list.empty()) continue;
				MeshViewSet* set = new MeshViewSet;
				for(int j = 0; j < pct; j++)
					set->addPoint(points[j], used_list.contains(j) ? 4 : 3);
				// -- (real) knots
				for(int i = 0; i < level; i++)
					for(int j = 0; j < level; j++)
						set->addPoint(m_points->get(i,j), (i==ki && j==kj) ? 1 : 2);
				// -- (knot) edges
				// ... horizontal
				for(int i = 1; i < level; i++)
					for(int j = 0; j < level; j++)
						set->addEdge(m_points->get(i,j), m_points->get(i-1,j), 0);
				// ... vertical
				for(int i = 0; i < level; i++)
					for(int j = 1; j < level; j++)
						set->addEdge(m_points->get(i,j), m_points->get(i,j-1), 0);
				// that's all
				ostringstream oss_fit;
				oss_fit << "initial SurfaceBSplinePlanar fit, min_dist2 = " << used_nearest(ki,kj);
				SHOW_MESH(oss_fit.str(), set);
			}
		}
	}
#endif

	// calculate base nodes from base points
	if(!calculateNodes()) return -1.0;

	// optimize
	int loops = 5; // 3;
	ostringstream oss_diff;
	while(loops-- > 0){
		DataMatrix<DVector3d> vd(level, level, DVector3d(0.0, 0.0, 0.0));
		DataMatrix<double> wd(level, level, 0.0);
		// gather corrections
		double max_dist2 = 0.0;
		for(int i = 0; i < pct; i++){
			const DPoint3d& pt = points[i];
			const DPoint2d param = getParameters(pt);
			const DPoint3d ptx = getPoint(param);
			DVector3d vdiff = pt - ptx;
			double dist2 = vdiff.length2(); // for checking only
			if(dist2 > max_dist2) max_dist2 = dist2;
			double pu = param.x / dx + 1.0;
			double pv = param.y / dy + 1.0;
			int iu = (int)pu;
			int iv = (int)pv;
			assert(iu >= 0 && iu < level-1);
			assert(iv >= 0 && iv < level-1);
			double du = pu - iu;
			double dv = pv - iv;
			double w;
			wd(iu, iv)     += (w = sqr((1.0 - du) * (1.0 - dv))); // 1.0 if close, 0.0 if far
			vd(iu, iv)     += vdiff * w;
			wd(iu+1, iv)   += (w = sqr(du * (1.0 - dv)));
			vd(iu+1, iv)   += vdiff * w;
			wd(iu, iv+1)   += (w = sqr((1.0 - du) * dv));
			vd(iu, iv+1)   += vdiff * w;
			wd(iu+1, iv+1) += (w = sqr(du * dv));
			vd(iu+1, iv+1) += vdiff * w;
		}
		oss_diff << sqrt(max_dist2) << " | ";
		// modify nodes
		for(int i = 1; i < level-1; i++){
			for(int j = 1; j < level-1; j++){
				double w = wd(i, j);
				if(w == 0.0) continue;
				m_points->set(i, j, m_points->get(i, j) + vd(i, j) / w);
			}
		}
		fitToPointsAdjustBorder(level);
		calculateNodes();
	}

	// measure difference
	double max_dist2 = 0.0;
#ifdef SHOW_BSURFACE_FIT
	DPoint3d max_pt_0;
	DPoint3d max_pt_1;
	DPoint3d max_pt_2;
#endif
	for(int i = 0; i < pct; i++){
		const DPoint3d& pt = points[i];
		DPoint2d param = getParameters(pt);
		const DPoint3d ptx = getPoint(param);
		double dist2 = pt.distance2(ptx);

		if(dist2 > max_dist2){
			max_dist2 = dist2;
#ifdef SHOW_BSURFACE_FIT
			max_pt_0 = pt;
			max_pt_1 = ptx;
			max_pt_2 = m_plane.getPoint(getParameters(pt));
#endif
		}
	}

#ifdef SHOW_BSURFACE_FIT
	if(true){
		MeshViewSet* set = new MeshViewSet;
		// -- cloud points
		for(int j = 0; j < pct; j++)
			set->addPoint(points[j], 3);
		// -- max diff points
		set->addPoint(max_pt_0, 0, 0);
		set->addPoint(max_pt_1, 0, 1);
		set->addPoint(max_pt_2, 0, 2);
		// -- (real) knots
		for(int i = 0; i < level; i++)
			for(int j = 0; j < level; j++){
				set->addPoint(m_points->get(i,j), 2);
			}
		// -- (knot) edges
		// ... horizontal
		for(int i = 1; i < level; i++)
			for(int j = 0; j < level; j++)
				set->addEdge(m_points->get(i,j), m_points->get(i-1,j), 0);
		// ... vertical
		for(int i = 0; i < level; i++)
			for(int j = 1; j < level; j++)
				set->addEdge(m_points->get(i,j), m_points->get(i,j-1), 0);
		// -- (control) knots
		for(int i = 0; i < level; i++)
			for(int j = 0; j < level; j++)
				set->addPoint(m_nodes->get(i,j), 1);
		// -- bspline points
		const int nl = 10*level;
		const double dnl = 1.0 / nl;
		for(int i = 0; i < nl; i++){
			for(int j = 0; j < nl; j++){
				DPoint2d knot(i*dnl, j*dnl);
				set->addPoint(getPoint(knot), 4);
			}
		}
		// that's all
		ostringstream oss;
		oss << "Initial SurfaceBSplinePlanar fit (max diff " << oss_diff.str() << "-> " << sqrt(max_dist2) << ")";
		SHOW_MESH(oss.str(), set);
	}
#endif

	// return max dist
	return sqrt(max_dist2);
}

/// fit to point-cloud - border adjust
void SurfaceBSplinePlanar::fitToPointsAdjustBorder(int level)
{
	// update the border knots
	for(int i = 1; i < level-1; i++){
		m_points->set(i, 0, m_points->get(i, 1) + 
			(m_points->get(i, 1) - m_points->get(i, 2)));
		m_points->set(i, level-1, m_points->get(i, level-2) + 
			(m_points->get(i, level-2) - m_points->get(i, level-3)));
		m_points->set(0, i, m_points->get(1, i) + 
			(m_points->get(1, i) - m_points->get(2, i)));
		m_points->set(level-1, i, m_points->get(level-2, i) + 
			(m_points->get(level-2, i) - m_points->get(level-3, i)));
	}
	// corners
	DPoint3d pt = DPoint3d::zero;
	pt.add(m_points->get(0, 1) + (m_points->get(0, 1) - m_points->get(0, 2)), 0.5);
	pt.add(m_points->get(1, 0) + (m_points->get(1, 0) - m_points->get(2, 0)), 0.5);
	m_points->set(0, 0, pt);
	pt = DPoint3d::zero;
	pt.add(m_points->get(0, level-2) + (m_points->get(0, level-2) - m_points->get(0, level-3)), 0.5);
	pt.add(m_points->get(1, level-1) + (m_points->get(1, level-1) - m_points->get(2, level-1)), 0.5);
	m_points->set(0, level-1, pt);
	pt = DPoint3d::zero;
	pt.add(m_points->get(level-2, 0) + (m_points->get(level-2, 0) - m_points->get(level-3, 0)), 0.5);
	pt.add(m_points->get(level-1, 1) + (m_points->get(level-1, 1) - m_points->get(level-1, 2)), 0.5);
	m_points->set(level-1, 0, pt);
	pt = DPoint3d::zero;
	pt.add(m_points->get(level-1, level-2) + (m_points->get(level-1, level-2) - m_points->get(level-1, level-3)), 0.5);
	pt.add(m_points->get(level-2, level-1) + (m_points->get(level-2, level-1) - m_points->get(level-3, level-1)), 0.5);
	m_points->set(level-1, level-1, pt);
}

/// Calculate difference (max square distance) between the surface and set of points
double SurfaceBSplinePlanar::calculateDiff2(const DataVector<DPoint3d>& points) const
{
	double max_dist2 = 0.0;
	for(int i = 0; i < points.countInt(); i++){
		const DPoint3d& pt = points[i];
		const DPoint3d ptx = getPoint(getParameters(pt));
		double dist2 = pt.distance2(ptx);
		if(dist2 > max_dist2){
			max_dist2 = dist2;
		}
	}
	return max_dist2;
}

bool SurfaceBSplinePlanar::calculateNodes()
{
	if(!m_points) return false;
	int level = m_points->rows();
	int n = level-2;
	if(m_nodes && m_nodes->rows() != level){
		delete m_nodes;
		m_nodes = nullptr;
	}
	if(!m_nodes){
		m_nodes = new DataMatrix<DPoint3d>(level, level, DPoint3d::zero);
	}

	// init boundary nodes
	for(int ri = 0; ri < level; ri++){ // left + right
		m_nodes->set(ri, 0,		m_points->get(ri, 0));	
		m_nodes->set(ri, n+1,	m_points->get(ri, n+1));
	}
	for(int ri = 1; ri < level-1; ri++){ // upper + bottom
		m_nodes->set(0,		ri,	m_points->get(0,	ri));	
		m_nodes->set(n+1,	ri,	m_points->get(n+1,	ri));
	}
	DataVector<double> A(n+1, 0.0);
	DataVector<DPoint3d> B(level, DPoint3d::zero);
	// for each row, calculate as with curve-BSpline
	for(int ri = 1; ri < level-1; ri++){
		// solving tridiagonal linear equation
		A[n] = -0.25;
		B[n].set(m_points->get(ri, n), 6.0);
		B[n].add(m_nodes->get(ri, n+1), -1.0);
		B[n] *= 0.25;
		for(int i = n; i > 2; i--){
			if(A[i] == -4.0) return false;
			A[i-1] = -1.0 / (A[i] + 4.0);
			B[i-1].set(m_points->get(ri, i-1), 6.0);
			B[i-1].add(B[i], -1.0);
			B[i-1] /= (A[i] + 4.0);
		}

		DPoint3d node;
		node.set(m_points->get(ri, 1), 6.0);
		node.add(m_nodes->get(ri, 0), -1.0);
		node.add(B[2], -1.0);
		node /= (A[2] + 4.0);
		m_nodes->set(ri, 1, node);
		for(int i = 1; i < n; i++){
			node.set(m_nodes->get(ri, i), A[i+1]);
			node.add(B[i+1]);
			m_nodes->set(ri, i+1, node);
		}
	}

	return true;
}

bool SurfaceBSplinePlanar::countSplineIndices(double& u, double& v, int& ui, int& vi) const
{
	if(!m_nodes) return false;

	// u,v should be [0,1]
	if(u < 0.0) u = 0.0;
	else if(u >= 1.0) u = 1.0 - SMALL_NUMBER;
	if(v < 0.0) v = 0.0;
	else if(v >= 1.0) v = 1.0 - SMALL_NUMBER;

	// [0,1] -> [0,n]
	int n = m_nodes->rows();
	if(n < 4) return false;

	u *= (n-3);
	v *= (n-3);

	double int_part;
	u = modf(u, &int_part);
	ui = (int)int_part;
	v = modf(v, &int_part);
	vi = (int)int_part;
	return true;
}

/// multi-fit to point-cloud using direct fit, using previous fit as a basis surface
std::shared_ptr<SurfaceMulti> SurfaceBSplinePlanar::multiFitToPointsRegular(
	std::shared_ptr<SurfaceBSplinePlanar> base_surface, 
	const DataVector<DPoint3d> & points, double& max_dist)
{
	auto msurface = std::make_shared<SurfaceMulti>(base_surface);

	// -> split in four equal parts ...
	double sub_max_diff2[2][2] = { {0.0, 0.0}, {0.0, 0.0} };
	DataVector<DPoint3d> sub_points[2][2];
	int sub_count = 0;

	// measure difference and arrange sub_points
	DataVector<DPoint2d> params(points.countInt());
	for(int i = 0; i < points.countInt(); i++){
		const DPoint3d& pt = points[i];
		DPoint2d param = base_surface->getParameters(pt);
		params.add(param); // will be used again...
		const DPoint3d ptx = base_surface->getPoint(param);
		double dist2 = pt.distance2(ptx);

		int iu = param.x < 0.5 ? 0 : 1;
		int iv = param.y < 0.5 ? 0 : 1;

		if(sub_max_diff2[iu][iv] < dist2) sub_max_diff2[iu][iv] = dist2;
		sub_points[iu][iv].add(pt);
	}

	MeshViewSet* sub_set = new MeshViewSet;
	MeshViewSet* sub_set_extra = new MeshViewSet;

	const DVector3d vn = base_surface->m_plane.getNormalVector(DPoint2d::zero);
	// create subsurfaces
	for(int ki = 0; ki < 2; ki++){
		for(int kj = 0; kj < 2; kj++){
			if(sub_points[ki][kj].countInt() < 3) continue; // too few for surface approximation
			// -> fit plane
			DPlane plane;
			double plane_max_dist = 
				DLeastSquaresFitting::fitHyperplaneOrthogonal(sub_points[ki][kj], plane);
			// -> avoid orthogonal normals 
			if(abs(vn.scalarProduct(plane.vn)) < 0.1){
				// --> skew plane a bit...
				plane.e0 += plane.vn * (0.1 * plane.e0.length());
				plane.e1 += plane.vn * (0.1 * plane.e1.length());
				plane.countNormal();
				assert(abs(vn.scalarProduct(plane.vn)) > 0.1);
			}
			// -> calculate plane and boundary rect
			DRect brect;
			double max_dist2 = 0.0;
			for(int j = 0; j < sub_points[ki][kj].countInt(); j++){
				const DPoint3d& pt = sub_points[ki][kj][j];
				DPoint2d param = plane.projectToPlane(pt);
				// - boundary
				brect.addPoint(param);
				// - distance
				const DPoint3d ptx = plane.projectToSpace(param);
				double dist2 = pt.distance2(ptx);
				if(dist2 > max_dist2) max_dist2 = dist2;
			}
			// --> calculate domain equation and local reparameterization
			SurfaceMulti::DomainRect* sub_domain = 
				new SurfaceMulti::DomainRect(
					DPoint2d(0.25 + 0.5 * ki, 0.25 + 0.5 * kj), DVector2d(0.5, 0.5), 0.15);
			// -> prepare subsurface structure
			// -> check distance and count
			if(sub_points[ki][kj].countInt() < FIT_PROBE_MIN){
				if(max_dist2 < sub_max_diff2[ki][kj]){
					// --> plane is (or has to be) sufficient as a subsurface
					SurfaceMulti::ReparameterizationPlanar* sub_reparam = 
						new SurfaceMulti::ReparameterizationPlanar(
							base_surface->m_plane, plane);
					msurface->insertSubSurface(
						std::make_shared<SurfacePlane>(plane.p0, plane.e0, plane.e1),
						sub_domain, sub_reparam);
					sub_count++;
				}
			}else{
				// --> fit bspline
				auto bsurf = std::make_shared<SurfaceBSplinePlanar>(plane, brect);
				double bmax_dist = bsurf->fitToPoints(sub_points[ki][kj], 2);
				auto bsurf_3 = std::make_shared<SurfaceBSplinePlanar>(plane, brect);
				double bmax_dist_3 = bsurf_3->fitToPoints(sub_points[ki][kj], 3);
				if(bmax_dist >= bmax_dist_3){
					bsurf = bsurf_3;
					bmax_dist = bmax_dist_3;
				}
				if(bmax_dist*bmax_dist < sub_max_diff2[ki][kj]){
					SurfaceMulti::ReparameterizationPlanar* sub_reparam = 
						new SurfaceMulti::ReparameterizationPlanar(
							base_surface->m_plane, bsurf->m_plane);
					msurface->insertSubSurface(bsurf, sub_domain, sub_reparam);
					sub_set = bsurf->getViewSet(sub_set, true);
					sub_count++;

					if(true){ // TEMPORARY
						int xrows = bsurf->m_points->rows();
						for(int xi = 0; xi < xrows; xi++){
							for(int xj = 0; xj < xrows; xj++){
								sub_set_extra->addPoint(bsurf->m_plane.getPoint(
									bsurf->m_plane.getParameters(
										bsurf->m_points->get(xi,xj))), 1);
							}
						}
					}

				}
			}
		}
	}

	//if(sub_count == 0){ // no multi-sub-surface was actually introduced, so no point in multi-surface ...
	//	delete msurface;
	//}

	// measure difference - final check
	double max_dist2 = 0.0;
	for(int i = 0; i < points.countInt(); i++){
		const DPoint3d& pt = points[i];
		DPoint2d param = base_surface->m_plane.getParameters(pt);
		const DPoint3d ptx = msurface->getPoint(param);
		double dist2 = pt.distance2(ptx);

		if(dist2 > max_dist2) max_dist2 = dist2;
	}

	if(true){
		int rows = base_surface->m_points->rows();
		for(int i = 0; i < rows; i++){
			for(int j = 0; j < rows; j++){
				sub_set_extra->addPoint(
					base_surface->m_plane.getPoint(
						base_surface->m_plane.getParameters(
							base_surface->m_points->get(i,j))), 0);
			}
		}
		for(int j = 0; j < points.countInt(); j++){
			sub_set_extra->addPoint(points[j], 3);
		}
		SHOW_MESH("multi-bsplines extra", sub_set_extra);
	}

	if(true){
		SHOW_MESH("multi-bsplines", sub_set);

		MeshViewSet* set = new MeshViewSet;
		// -- cloud points
		for(int j = 0; j < points.countInt(); j++){
			set->addPoint(points[j], 3);
		}
		// -- bspline points
		const int nl = 50;
		const double dnl = 1.0 / nl;
		for(int i = 0; i < nl; i++){
			for(int j = 1; j < nl; j++){
				DPoint2d knot1(i*dnl, j*dnl);
				DPoint2d knot2(i*dnl, (j-1)*dnl);
				set->addEdge(base_surface->getPoint(knot1), base_surface->getPoint(knot2), 0);
				set->addEdge(msurface->getPoint(knot1), msurface->getPoint(knot2), 2);
				knot1 = DPoint2d(j*dnl, i*dnl);
				knot2 = DPoint2d((j-1)*dnl, i*dnl);
				set->addEdge(base_surface->getPoint(knot1), base_surface->getPoint(knot2), 0);
				set->addEdge(msurface->getPoint(knot1), msurface->getPoint(knot2), 2);
			}
		}
		// that's all
		ostringstream oss;
		oss << "MultiSurface fit (max diff " << max_dist << " -> " << sqrt(max_dist2) << ")";
		SHOW_MESH(oss.str(), set);
	}

	// return results
	max_dist = sqrt(max_dist2);
	return msurface;
}

/// multi-fit to point-cloud using direct fit, using previous fit as a basis surface
std::shared_ptr<SurfaceMulti> SurfaceBSplinePlanar::multiFitToPoints(
	std::shared_ptr<SurfaceBSplinePlanar> base_surface, 
	const DataVector<DPoint3d> & points, 
	int min_level, int max_level, double tolerance, double& max_dist)
{
	auto msurface = std::make_shared<SurfaceMulti>(base_surface);

	// measure difference
	int level = base_surface->m_nodes->rows();
	double dx = 1.0 / (level-3);
	double dy = 1.0 / (level-3);
	DataMatrix<double> max_dist2_matrix(level, level, 0.0);
	DataVector<DPoint2d> params(points.countInt());
	for(int i = 0; i < points.countInt(); i++){
		const DPoint3d& pt = points[i];
		DPoint2d param = base_surface->m_plane.getParameters(pt);
		params.add(param); // will be used again...
		const DPoint3d ptx = base_surface->getPoint(param);
		double dist2 = pt.distance2(ptx);

		int iu = (int)(param.x / dx + 1.0);
		int iv = (int)(param.y / dy + 1.0);

		if(dist2 > max_dist2_matrix(iu, iv)){
			max_dist2_matrix(iu, iv) = dist2;
		}
	}

	const DVector3d vn = base_surface->m_plane.getNormalVector(DPoint2d::zero);
	// create subsurfaces
	double tol2 = tolerance * tolerance;
	for(int ki = 1; ki < level-2; ki++){
		for(int kj = 1; kj < level-2; kj++){
			LOG4CPLUS_INFO(MeshLog::logger_mesh, "ki=" << ki << ", kj=" <<  kj << ", dist2=" << max_dist2_matrix(ki,kj));
			if(max_dist2_matrix(ki,kj) <= tol2) continue;
			// -> select points for this subdomain (rectangle + margin
			double mid_i = ki + 0.5;
			double mid_j = kj + 0.5;
			double range = 1.0;
			DataVector<DPoint3d> selected_points(points.countInt());
			for(int i = 0; i < points.countInt(); i++){
				DPoint2d& param = params[i];
				double u = (param.x / dx + 1.0);
				double v = (param.y / dy + 1.0);
				if(abs(mid_i - u) < range && abs(mid_j - v) < range)
					selected_points.add(points[i]);
			}
			if(selected_points.countInt() < 3) continue;
			DPlane plane;
			double plane_max_dist = 
				DLeastSquaresFitting::fitHyperplaneOrthogonal(
					selected_points, plane);
			// -> avoid orthogonal normals 
			if(abs(vn.scalarProduct(plane.vn)) < 0.1){
				// --> skew plane a bit...
				plane.e0 += plane.vn * (0.1 * plane.e0.length());
				plane.e1 += plane.vn * (0.1 * plane.e1.length());
				plane.countNormal();
				assert(abs(vn.scalarProduct(plane.vn)) > 0.1);
			}
			// -> calculate plane and boundary rect
			DRect brect;
			double max_dist2 = 0.0;
			for(int j = 0; j < selected_points.countInt(); j++){
				const DPoint3d& pt = selected_points[j];
				DPoint2d param = plane.projectToPlane(pt);
				// - boundary
				brect.addPoint(param);
				// - distance
				const DPoint3d ptx = plane.projectToSpace(param);
				double dist2 = pt.distance2(ptx);
				if(dist2 > max_dist2) max_dist2 = dist2;
			}
			// --> calculate domain equation and local reparameterization
			SurfaceMulti::DomainRect* sub_domain = 
				new SurfaceMulti::DomainRect(
					DPoint2d((mid_i-1) / (level-3), (mid_j-1) / (level-3)), 
					DVector2d(0.5 / (level-3), 0.5 / (level-3)), 0.15);
			// -> prepare subsurface structure
			// -> check distance and count
			if(selected_points.countInt() < FIT_PROBE_MIN){ // || max_dist2 < tol2){
				// --> plane is (or has to be) sufficient as a subsurface
				SurfaceMulti::ReparameterizationPlanar* sub_reparam = 
					new SurfaceMulti::ReparameterizationPlanar(
						base_surface->m_plane, plane);
				msurface->insertSubSurface(
					std::make_shared<SurfacePlane>(plane.p0, plane.e0, plane.e1),
					sub_domain, sub_reparam);
			}else{
				// --> fit bspline
				auto bsurf = std::make_shared<SurfaceBSplinePlanar>(plane, brect);
				for(int lev = min_level; lev <= max_level; lev++){
					double bmax_dist = bsurf->fitToPoints(selected_points, lev);
					if(bmax_dist < tolerance) break;
				}
				SurfaceMulti::ReparameterizationPlanar* sub_reparam = 
					new SurfaceMulti::ReparameterizationPlanar(
						base_surface->m_plane, bsurf->m_plane);
				msurface->insertSubSurface(bsurf, sub_domain, sub_reparam);
			}
		}
	}

	// measure difference - final check
	double max_dist2 = 0.0;
	for(int i = 0; i < points.countInt(); i++){
		const DPoint3d& pt = points[i];
		DPoint2d param = base_surface->m_plane.getParameters(pt);
		const DPoint3d ptx = msurface->getPoint(param);
		double dist2 = pt.distance2(ptx);

		if(dist2 > max_dist2) max_dist2 = dist2;
	}

	if(true){
		MeshViewSet* set = new MeshViewSet;
		// -- cloud points
		for(int j = 0; j < points.countInt(); j++)
			set->addPoint(points[j], 3);
		// -- (real) knots
		for(int i = 0; i < level; i++)
			for(int j = 0; j < level; j++){
				set->addPoint(base_surface->m_points->get(i,j), 2);
			}
		// -- (knot) edges
		// ... horizontal
		for(int i = 1; i < level; i++)
			for(int j = 0; j < level; j++)
				set->addEdge(base_surface->m_points->get(i,j), base_surface->m_points->get(i-1,j), 0);
		// ... vertical
		for(int i = 0; i < level; i++)
			for(int j = 1; j < level; j++)
				set->addEdge(base_surface->m_points->get(i,j), base_surface->m_points->get(i,j-1), 0);
		// -- (control) knots
		for(int i = 0; i < level; i++)
			for(int j = 0; j < level; j++)
				set->addPoint(base_surface->m_nodes->get(i,j), 1);
		// -- bspline points
		const int nl = 10*level;
		const double dnl = 1.0 / nl;
		for(int i = 0; i < nl; i++){
			for(int j = 0; j < nl; j++){
				DPoint2d knot(i*dnl, j*dnl);
				DPoint3d bpoint = base_surface->getPoint(knot);
				set->addPoint(bpoint, 4);
				DPoint3d mpoint = msurface->getPoint(knot);
				if(bpoint.distance(mpoint) > SMALL_NUMBER)
					set->addPoint(mpoint, 0);
			}
		}
		// that's all
		ostringstream oss;
		oss << "MultiSurface fit (max diff " << max_dist << " -> " << sqrt(max_dist2) << ")";
		SHOW_MESH(oss.str(), set);
	}

	// return results
	max_dist = sqrt(max_dist2);
	return msurface;
}

/// Store XML description to stream
ostream& SurfaceBSplinePlanar::storeXML(ostream& os, const string& prefix) const
{
	os << prefix << "<bspline-surface type=\"plane-based\">" << endl;
	m_plane.storeXML(os, prefix+"\t");
	if(m_points){
		int rows = m_points->rows();
		int cols = m_points->countInt() / m_points->rows();
		os << prefix << "\t<rows>" << rows << "</rows>\n";
		os << prefix << "\t<columns>" << cols << "</columns>\n";
		os << prefix << "\t<points>\n";
		for(int ki = 0; ki < rows; ki++){
			for(int kj = 0; kj < cols; kj++){
				os << prefix << "\t\t<point row=\"" << ki << "\" column=\"" << kj << "\"> "
					<< m_points->get(ki, kj) << " </point>\n";
			}
		}
		os << prefix << "\t</points>\n";
	}
	return os << prefix << "</bspline-surface>" << endl;
}

/// get view set for visualization
MeshViewSet* SurfaceBSplinePlanar::getViewSet(MeshViewSet* set, bool with_knots) const
{
	if(!m_nodes) return set;

	if(!set) set = new MeshViewSet;
	int level = m_nodes->rows();

	if(with_knots){
		// -- (control) knots
		for(int i = 0; i < level; i++)
			for(int j = 0; j < level; j++)
				set->addPoint(m_nodes->get(i,j), 2);
		// -- (knot) edges
		// ... horizontal
		for(int i = 1; i < level; i++)
			for(int j = 0; j < level; j++)
				set->addEdge(m_nodes->get(i,j), m_nodes->get(i-1,j), 2);
		// ... vertical
		for(int i = 0; i < level; i++)
			for(int j = 1; j < level; j++)
				set->addEdge(m_nodes->get(i,j), m_nodes->get(i,j-1), 2);
	}

	// -- bspline points
	const int GRID = 10;
	const int nl = GRID*(level-3);
	DataMatrix<DPoint3d> pnet(nl+1, nl+1, DPoint3d::zero);
	const double dnl = 1.0 / nl;
	// ... count
	for(int i = 0; i <= nl; i++){
		for(int j = 0; j <= nl; j++){
			DPoint2d pt_param(i*dnl, j*dnl);
			pnet(i,j) = getPoint(pt_param);
		}
	}
	// ... draw
	for(int i = 0; i <= nl; i++){
		for(int j = 1; j <= nl; j++){
			set->addEdge(pnet(i,j), pnet(i, j-1), i % GRID == 0 ? 1 : 0);
			set->addEdge(pnet(j,i), pnet(j-1, i), i % GRID == 0 ? 1 : 0);
		}
	}

	return set;
}