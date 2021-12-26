// MeshBRepRSurf.cpp: implementation of the MeshBRepMsh class.
//
//////////////////////////////////////////////////////////////////////

#include "MeshBRepRSurf.h"
#include "MeshData.h"
#include "DPoint.h"

#include "DataHashTable.h"
#include "MeshBoundaryCondition.h"
#include "MeshViewSet.h"

/// Parse model description data in brep (resu.surf) .txt format
bool MeshBRepRSurf::parseFile(const string& /* fname */)
{
	return false;
	/*
	ifstream ifs(fname.c_str());
	int mode = 0; // 1 - points, 2 - triangles
	DataVector<DPoint3d> real_points(1000);
	DataVector<int> real_triangles(1000);
	DataVector<int> real_quads(1000);
	double coord[3];
	int ci = 0;
	while(ifs){
		string line;
		getline(ifs, line);
		if (line.find("END") == 0) mode = 0;
		else if(mode > 0){
			istringstream iss(line);
			if(mode == 1){
				while(iss){
					iss >> coord[ci];
					if(!iss) break;
					ci = (ci+1)%3;
					if(ci == 0)
						real_points.add(DPoint3d(coord[0], coord[1], coord[2]));
				}
			}else if(mode == 2){
				while(iss){
					int vertex;
					iss >> vertex;
					if(iss) real_triangles.add(vertex-1);
				}
			}else if(mode == 3){
				while(iss){
					int vertex;
					iss >> vertex;
					if(iss) real_quads.add(vertex-1);
				}
			}
		}else if(line.find("POINTS") == 0) mode = 1;
		else if(line.find("TRIANGLES") == 0) mode = 2;
		else if(line.find("QUADS") == 0) mode = 3;
	}
	if(real_points.countInt() < 3) {
		LOG4CPLUS_ERROR(MeshLog::logger_console, "Too few points", real_points.countInt());
		return nullptr;
	}
	if(ci != 0) {
		LOG4CPLUS_ERROR(MeshLog::logger_console,   "Inconsistent number of coordinates for points");
		return nullptr;
	}
	if(real_triangles.countInt() == 0) {
		LOG4CPLUS_ERROR(MeshLog::logger_console,   "No triangles?");
		return nullptr;
	}
	if(real_triangles.countInt() %3 > 0) {
		LOG4CPLUS_ERROR(MeshLog::logger_console,   "Inconsistent number of vertices for triangles");
		return nullptr;
	}
	if(real_quads.countInt() %3 > 0) {
		LOG4CPLUS_ERROR(MeshLog::logger_console,   "Inconsistent number of vertices for quads");
		return nullptr;
	}

	const DPoint3d spt0 = real_points[real_triangles[0]];
	const DPoint3d spt1 = real_points[real_triangles[1]];
	const DPoint3d spt2 = real_points[real_triangles[2]];

	const DVector3d e0 = (spt1 - spt0).normalized();
	DVector3d nv = DVector3d::crossProduct(spt0, spt1, spt2);
	const DVector3d e1 = nv.crossProduct(e0).normalized();
	SurfacePlane *plane = new SurfacePlane(spt0, e0, e1);
	MeshContainer2d* mesh = new MeshContainer2d(real_points.countInt());
	plane->incRefCount();
	mesh->setSurface(plane);

	// add mesh points
	int pct = real_points.countInt();
	for(int i = 0; i < pct; i++)
		mesh->addMeshPoint(new MeshPoint2d(plane->getParameters(real_points[i])));

	// add mesh elements (triangles)
	int tct = real_triangles.countInt() / 3;
	for(int i = 0; i < tct; i++)
		mesh->addMeshElement(new MeshTriangle2d(
			mesh->getPointAt(real_triangles[3*i]),
			mesh->getPointAt(real_triangles[3*i+1]),
			mesh->getPointAt(real_triangles[3*i+2])));
	// add mesh elements (triangles)
	int qct = real_quads.countInt() / 4;
	for(int i = 0; i < qct; i++)
		mesh->addMeshElement(new MeshQuad2d(
			mesh->getPointAt(real_triangles[4*i]),
			mesh->getPointAt(real_triangles[4*i+1]),
			mesh->getPointAt(real_triangles[4*i+2]),
			mesh->getPointAt(real_triangles[4*i+3])));

	// set boundary data
	for(int i = 0; i < tct; i++){
		MeshElement* element = mesh->getElementAt(i);
		element->setAreaID(0);
		if(element->isInverted()) LOG4CPLUS_WARN(MeshLog::logger_console, "Inverted element", i);
		int edge_count = element->getEdgeCount();
		for(int j = 0; j < edge_count; j++){
			MeshEdge2d* edge = element->getEdge(j);
			if(!edge->getMeshElement(0) || !edge->getMeshElement(1)){
				edge->setBorderType(0);
				edge->getMeshPoint(0)->setBorder();
				edge->getMeshPoint(1)->setBorder();
			}
		}
	}

	// check points
	DataVector<MeshPoint2d*> stranded_points;
	for(int i = 0; i < pct; i++){
		MeshPoint2d* point = mesh->getPointAt(i);
		if(point->getRank() < 2){
			LOG4CPLUS_WARN(MeshLog::logger_console, "Stranded point", i+1);
			stranded_points.add(point);
		}
	}
	for(int i = 0; i < stranded_points.countInt(); i++){
		delete mesh->removeMeshPoint(stranded_points[i]);
	}

	return mesh;
	*/
}

/// Stores the description of mesh-surface to resu.surf file
bool MeshBRepRSurf::storeFileResuSurf(const string& /* fname */, const MeshContainer2d* /* mesh */)
{
	return false;
	/*
	static char* FLOW_FILE_HEADER[] =
	{
		"FORMAT  ASCII",
			"# ##### Header #####",
			"HEADER",
			"PROJECT Diamesh",
			"COMMENT Diamesh.surf",
			"DATE    Tue Jan 28 15:15:52 1997",
			"HOUR    Tue Jan 28 15:15:52 1997",
			"CODE    DIAMESH",
			"VERSION 1.0 Alpha",
			"TYPE    SURFACE_DUMP",
			"END",
			"# ##### Prameters #####",
			"PARAMETERS",
			nullptr,
			"END",
			nullptr
	};

	ofstream file(fname.c_str());

	if(!file){
		LOG4CPLUS_ERROR(MeshLog::logger_console, "Can't open file for write", fname);
		return false;
	}

// Write the header block and the parameters 
	int pct = mesh->getPointsCount();
	int ect = mesh->getElementsCount();
	int	tct = mesh->getElementsCount(3);
	int qct = mesh->getElementsCount(4);
	char **header = FLOW_FILE_HEADER;
	while(*header) file << *header++ << endl;

	file << "SHAPE" << endl;
	file << "NBR_POINTS    " << pct << endl;
	file << "NBR_TRIANGLES " << tct << endl;
	file << "NBR_QUADS     " << qct << endl;
//   fprintf(fp,"NBR_LINES     %d\n",DiaGetNbExteriorCurves());
	file << "END" << endl;

// Write the points 
	file << "#    ############################ Points ############################" << endl;
	file << "POINTS" << endl;
	SurfaceParametric* surface = mesh->getSurface();
	for(int i = 0; i < pct; i++){
		const DPoint3d pt = surface->getPoint(mesh->getPointAt(i)->getCoordinates());		
		file << fixed << pt.x << " " << pt.y << " " << pt.z << " ";
		if(i%3 == 2) file << endl;
	}
	if(pct % 3 > 0) file << endl;
	file << "END" << endl;

// Write the triangles 
   if(ect > 0){
     file << "#    ############################ Triangles ############################\n";
     file << "SURFACE 1 surface1\n";
	 if(qct > 0){
		 file << "QUADS\n";
		 for(int i=0; i<ect; i++){
			 MeshElement* element = mesh->getElementAt(i);
			 if(element->getEdgeCount() != 4) continue;
			 for(int j = 0; j < 4; j++)
				 file << (element->getPoint(j)->getIndex() + 1) << " ";
			 file << endl;
		 }
		 file << "END" << endl;
	 }
	 if(tct > 0){
		 file << "TRIANGLES\n";
		 for(int i=0; i<ect; i++){
			 MeshElement* element = mesh->getElementAt(i);
			 if(element->getEdgeCount() != 3) continue;
			 for(int j = 0; j < 3; j++)
				 file << (element->getPoint(j)->getIndex() + 1) << " ";
			 file << endl;
		 }
		 file << "END" << endl;
	 }
	 file << "END" << endl;
   }

	return true;
	*/
}
