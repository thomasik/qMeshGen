// MeshSpecialRoutinesPaszynski.cpp: implementation of the MeshSpecialRoutinesPaszynski class.
//
//////////////////////////////////////////////////////////////////////

#include "MeshData.h"
#include "MeshSpecialRoutinesPaszynski.h"
#include "MeshContainer3d.h"
#include "MeshPoint3d.h"
#include "MeshEdge3d.h"
#include "MeshTriangle3d.h"
#include "MeshTetrahedron.h"
#include "MeshDomainSurface.h"
#include "MeshDomainVolume.h"
#include "MeshLog.h"

MeshContainer3d* MeshSpecialRoutinesPaszynski::loadTetrahedralMesh(const string& fname)
{
	ifstream file(fname.c_str());
	if(!file){
		LOG4CPLUS_ERROR(MeshLog::logger_console, "Error opening the file: " << fname);
		return nullptr;
	}
	// header
	string h1, h2, h3, line;
	getline(file, h1);
	getline(file, line);
	while(file && line != ""){
		h1.append("\n");
		h1.append(line);
		getline(file, line);
	}

	getline(file, h2);
	getline(file, line);
	while(file && line != ""){
		h2.append("\n");
		h2.append(line);
		getline(file, line);
	}

	getline(file, h3);
	getline(file, line);
	while(file && line != ""){
		h3.append("\n");
		h3.append(line);
		getline(file, line);
	}

	// vertices
	int p_count;
	file >> p_count;
	if(!file){
		LOG4CPLUS_ERROR(MeshLog::logger_console,   "Syntax error reading the file");
		return nullptr;
	}
	if(p_count < 1){
		LOG4CPLUS_ERROR(MeshLog::logger_console,   "Invalid number of points");
		return nullptr;
	}

	MeshContainer3d* volume_mesh = new MeshContainer3d(p_count);
	for(int i = 0; i < p_count; i++){
		double x,y,z;
		file >> x >> y >> z;
		if(!file){
			LOG4CPLUS_ERROR(MeshLog::logger_console,   "Syntax error reading the points");
			return nullptr;
		}
		MeshPoint3d* point = new MeshPoint3d(x, y, z);
		point->setIntTag(TagExtended::TAG_ID, i+1);
		volume_mesh->addMeshPoint(point);
	}

	// tetrahedra
	int t_count;
	file >> t_count;
	if(!file){
		LOG4CPLUS_ERROR(MeshLog::logger_console,   "Syntax error reading the file");
		return nullptr;
	}
	if(t_count < 1){
		LOG4CPLUS_ERROR(MeshLog::logger_console,   "Invalid number of tetrahedra");
		return nullptr;
	}

	for(int i = 0; i < t_count; i++){
		int id, vid0, vid1, vid2, vid3;
		file >> id >> vid0 >> vid1 >> vid2 >> vid3;
		if(!file){
			LOG4CPLUS_ERROR(MeshLog::logger_console,   "Syntax error reading the tetrahedra");
			return nullptr;
		}
		MeshTetrahedron* tetra = new MeshTetrahedron(
			volume_mesh->getPointAt(vid0-1), // since indices are 1-based
			volume_mesh->getPointAt(vid1-1),
			volume_mesh->getPointAt(vid2-1),
			volume_mesh->getPointAt(vid3-1));
		tetra->setAreaID(0);
		tetra->setIntTag(TagExtended::TAG_ID, i+1);
		volume_mesh->addMeshTetrahedron(tetra);
	}

	// ** create domain volume
	MeshDomainVolume* domain_volume = new MeshDomainVolume(volume_mesh); 
	domain_volume->setAreaID(0);
	// ?? create and set boundary mesh
	// ** set volume mesh

	domain_volume->createInitialControlSpace();

	domain_volume->copySurfaceMeshFromVolumeMesh();

	MeshContainer3d *domain = new MeshContainer3d(10);
	domain->addMeshBlock(domain_volume);

	return domain;
}

bool MeshSpecialRoutinesPaszynski::storeTetrahedralMesh(const string& fname, 
				const string& header, MeshContainer3d* mesh)
{
	ofstream file(fname.c_str());
	if(!file){
		LOG4CPLUS_ERROR(MeshLog::logger_console, "Error opening the file: " << fname);
		return false;
	}

	file << header << endl;

	int pct = mesh->getPointsCount();
	file << pct << endl << fixed;
	file.precision(6);
	for(int i = 0; i < pct; i++){
		const DPoint3d& pt = mesh->getPointAt(i)->getCoordinates();
		file << pt.x << " " << pt.y << " " << pt.z << endl;
	}

	int tct = mesh->getBlocksCount();
	file << tct << endl;
//	const int rct[] = {1, 0, 2, 3}; // for switching orientation of tetrahedra (sequence of vertices)
	const int rct[] = {0, 1, 2, 3}; // oryginal orientation of tetrahedra (sequence of vertices)
	DataHashTable<MeshFace*> hash_faces(6*tct, nullptr);
	for(int i = 0; i < tct; i++){
		MeshTetrahedron* tetra = (MeshTetrahedron*) mesh->getBlockAt(i);
		assert(tetra->getType() == BLOCK_TETRA);
		file << "  1";
		for(int j = 0; j < 4; j++){
			file << "  " << (tetra->getPoint(rct[j])->getIndex() + 1);
			MeshFace* face = tetra->getFace(j);
			if(face->isBorder())
				hash_faces.insert(face);
		}
		file << endl;
	}

	int fct = hash_faces.countInt();
	DataVector<MeshFace*> faces(fct);
	hash_faces.getValues(faces);
	file << fct << endl;
	for(int i = 0; i < fct; i++){
		file << " 1";
		for(int j = 0; j < 3; j++){
			file << "  " << (faces[i]->getPoint(j)->getIndex() + 1);
		}
		file << endl;
	}

	return true;
}
