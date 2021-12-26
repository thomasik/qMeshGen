// MeshBlock.cpp: implementation of the MeshBlock class.
//
//////////////////////////////////////////////////////////////////////

#include "MeshBlock.h"
#include "MeshPoint3d.h"
#include "DRect.h"
#include "MeshFace.h"
#include "MeshEdge3d.h"
#include "MeshData.h"
#include "MeshTriangle3d.h"
#include "MeshQuad3d.h"
#include "MeshElement.h"
#include "MeshTetrahedron.h"
#include "DMetric3d.h"

int MeshBlock::param_point_distance = MeshData::RP_VERTICES;

MeshBlock::MeshBlock(int fct, int pct, MeshFace** new_faces, MeshPoint3d** new_points) 
	: face_count(fct), faces(new_faces), point_count(pct), points(new_points), 
		area_id(-1), quality(-1.0) { }

DPoint3d MeshBlock::getMiddlePoint() const
{
	assert(point_count > 0);
	DPoint3d middle = points[0]->getCoordinates();
	for(int i = 1; i < point_count; i++)
		middle.add(points[i]->getCoordinates());
	middle /= (double)point_count;
	return middle;
}

DBox MeshBlock::getBoundingBox() const
{
	DBox box;
	for(int i = 0; i < point_count; i++){
		box.addPoint(points[i]->getCoordinates());
	}
	return box;
}

MeshFace* MeshBlock::getIncidentFace(MeshFace *face, MeshEdge3d *edge)
{
	for(int i = 0; i < face_count; i++){
		if(faces[i] != face && faces[i]->getEdgeIndex(edge) >= 0)
			return faces[i];
	}

	return nullptr;
}

short MeshBlock::compareTo(MeshBlock *item)
{
	double q1 = getQuality();
	double q2 = item->getQuality();
	if(q1 > q2){
		return 1;
	}else if(q1 < q2){
		return -1;
	}else return 0;
}

bool MeshBlock::isAdjacentTo(const MeshFace* face) const
{
	for(int i = 0; i < face_count; i++)
		if(faces[i] == face) return true;
	return false;
}

bool MeshBlock::hasBoundaryFace() const
{
	for(int i = 0; i < face_count; i++)
		if(faces[i]->isBorder()) return true;
	return false;
}

bool MeshBlock::hasBoundaryEdge() const
{
	for(int i = 0; i < face_count; i++){
		int edge_count = faces[i]->getEdgeCount();
		for(int j = 0; j < edge_count; j++)
			if(faces[i]->getEdge(j)->isBorder()) return true;
	}
	return false;
}

bool MeshBlock::hasBoundaryVertex() const
{
	for(int i = 0; i < point_count; i++)
		if(points[i]->isBorder()) return true;
	return false;
}

MeshPoint3d * MeshBlock::getOppositePoint(const MeshFace * mf) const
{
	// here, any point not incident to the face, but should be overwritten in extending classes
	for (int i = 0; i < point_count; i++) {
		if (!mf->incidentToPoint(points[i])) return points[i];
	}
	return nullptr;
}

int MeshBlock::getFaceCount(int ect) const
{
	int ct = 0;
	for(int i = 0; i < face_count; i++)
		if(faces[i]->getEdgeCount() == ect) ++ct;
	return ct;
}

void MeshBlock::switchPointsWithFacesBoundary(const MeshPoint3d *point1, MeshPoint3d *point2)
{
	switchPointsWithFaces(point1, point2); // for now
}

bool MeshBlock::checkSwitchPointsWithFaces(const MeshPoint3d *point1, MeshPoint3d *point2) const
{
	assert(false);
	return true; // implemented in MeshTetrahedron for now
}

void MeshBlock::switchPointsWithFaces(const MeshPoint3d *point1, MeshPoint3d *point2)
{
	assert(false);
	LOG4CPLUS_WARN(MeshLog::logger_console, "unstable version - needs correction!");
	for(int i = 0; i < point_count; i++){
		if(points[i] == point1){
			points[i] = point2;
			for(int j = 0; j < face_count; j++){
				if(faces[j]->incidentToPoint(point1)){
					MeshPoint3d* new_points[4] = {nullptr, nullptr, nullptr, nullptr};
					int ect = faces[j]->getEdgeCount();
					assert(ect == 3 || ect == 4);
					for(int k = 0; k < ect; k++)
						new_points[k] = (faces[j]->getPoint(k) == point1) ? point2 : faces[j]->getPoint(k);
					if(faces[j]->removeBlockLink(this)) delete faces[j];
					faces[j] = new_points[0]->getFaceToPoints(new_points[1], new_points[2]);
					if(faces[j])
						faces[j]->setBlockLink(this, new_points[0], new_points[1]);
					else if(ect == 3)
						faces[j] = new MeshTriangle3d(new_points[0], new_points[1], new_points[2], this);
					else 
						faces[j] = new MeshQuad3d(new_points[0], new_points[1], new_points[2], new_points[3], this);
				}
			}
			return;
		}
	}
	assert(false);
}

void MeshBlock::removeTagForEdgesAndFaces(TagExtended::TagType tag_type)
{
	for(int i = 0; i < face_count; i++) faces[i]->removeTag(tag_type);
	int ect = getEdgeCount();
	for(int i = 0; i < ect; i++) getEdge(i)->removeTag(tag_type);
}

void MeshBlock::setIntTagForFaces(TagExtended::TagType tag_type, int t)
{ 
	for(int i = 0; i < face_count; i++) 
		faces[i]->setIntTag(tag_type, t); 
}

double MeshBlock::distance2(const DPoint3d& pt) const
{
	double min_dist2 = mesh_data.relative_infinity;
	if((param_point_distance & MeshData::RP_MIDDLE) != 0){
		double dist2 = pt.distance2(getMiddlePoint());
		if(dist2 < min_dist2) min_dist2 = dist2;
	}
	if((param_point_distance & MeshData::RP_VERTICES) != 0)
		for(int i = 0; i < point_count; i++){
			double dist2 = pt.distance2(getPoint(i)->getCoordinates());
			if(dist2 < min_dist2) min_dist2 = dist2;
		}
	if((param_point_distance & MeshData::RP_MIDEDGES) != 0){
		int ect = getEdgeCount();
		for(int i = 0; i < ect; i++){
			double dist2 = pt.distance2(getEdge(i)->getPoint(0.5));
			if(dist2 < min_dist2) min_dist2 = dist2;
		}
	}
	if((param_point_distance & MeshData::RP_MIDFACES) != 0){
		for(int i = 0; i < face_count; i++){
			double dist2 = pt.distance2(getFace(i)->getMiddlePoint());
			if(dist2 < min_dist2) min_dist2 = dist2;
		}
	}
	return min_dist2;
}

void MeshBlock::getRepresentativePoints(DataVector<DPoint3d> & pts) const
{
	if((MeshBlock::param_point_distance & MeshData::RP_MIDDLE) != 0)
		pts.add(getMiddlePoint());
	if((MeshBlock::param_point_distance & MeshData::RP_VERTICES) != 0)
		for(int i = 0; i < point_count; i++)
			pts.add(getPoint(i)->getCoordinates());
	if((MeshBlock::param_point_distance & MeshData::RP_MIDEDGES) != 0){
		int ect = getEdgeCount();
		for(int i = 0; i < ect; i++)
			pts.add(getEdge(i)->getPoint(0.5));
	}
	if((MeshBlock::param_point_distance & MeshData::RP_MIDFACES) != 0)
		for(int i = 0; i < face_count; i++)
			pts.add(getFace(i)->getMiddlePoint());
}

/// add this element as a lump of mass (respecting control space)
void MeshBlock::addForInertialCenter(Metric3dContext& mc, DPoint3d& inertial_center, double& total_mass) const
{
	double volume = getVolume(mc);
	total_mass += volume;
	inertial_center.add(getMiddlePoint(), volume);
}

/// add this element as a lump of mass (respecting control space)
void MeshBlock::addForInertialMoments(Metric3dContext& mc, const DPoint3d& inertial_center, DMatrix3d& inertial_moments) const
{
	double volume = getVolume(mc);
	const DVector3d dv = getMiddlePoint() - inertial_center;
	inertial_moments.m[0][0] += volume * (dv.y*dv.y + dv.z*dv.z); // Ixx
	inertial_moments.m[1][1] += volume * (dv.x*dv.x + dv.z*dv.z); // Iyy
	inertial_moments.m[2][2] += volume * (dv.y*dv.y + dv.x*dv.x); // Izz
	inertial_moments.m[0][1] += -volume * (dv.x*dv.y); // Ixy
	inertial_moments.m[1][0] += -volume * (dv.x*dv.y); // Iyx = Ixy
	inertial_moments.m[0][2] += -volume * (dv.x*dv.z); // Ixz
	inertial_moments.m[2][0] += -volume * (dv.x*dv.z); // Izx = Ixz
	inertial_moments.m[1][2] += -volume * (dv.y*dv.z); // Iyz
	inertial_moments.m[2][1] += -volume * (dv.y*dv.z); // Izy = Iyz
}

/// Returns the reference to the block incident through the i-th face(can be nullptr)
MeshBlock* MeshBlock::getNeighbour(int i) const 
{ 
	return faces[i]->getOtherBlock(this); 
}

void MeshBlock::setAreaID(int id)
{ 
	area_id = id; 
}

/// Checks wheterh the i-th face of the block belongs to the boundary
bool MeshBlock::isBorder(int i) const 
{ 
	return faces[i]->isBorder(); 
}

bool MeshBlock::checkIntTagForAnyPoint(TagExtended::TagType tag_type, int tag_value) const
{
	for(int i = 0; i < point_count; i++)
		if(points[i]->checkIntTag(tag_type, tag_value)) return true;

	return false;
}
