/////////////////////////////////////////////////////////////////////////////
// MeshModel.cpp
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2008-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#include<iostream>
#include<string>
#include<sstream>
#include<algorithm>
#include<iterator>

#include "MeshModel.h"
#include "common.h"

#include "IteratorMesh2d.h"
#include "IteratorMesh3d.h"
#include "MeshContainer2d.h"
#include "MeshContainer3d.h"
#include "MeshContainer3dSurface.h"
#include "MeshDomainVolume.h"
#include "MeshDomainSurface.h"
#include "DEquation.h"
#include "MeshGenerator1d.h"
#include "MeshGenerator2d.h"
#include "MeshGenerator2dQuad.h"
#include "MeshGenerator3d.h"
#include "MeshGenerator3dQuality.h"
#include "MeshSpecialRoutinesDAT.h"
#include "MeshSpecialRoutinesPaszynski.h"
#include "MeshBRepXML.h"
#include "ControlSpace2dAdaptive.h"
#include "ControlSpace2dMesh.h"
#include "ControlSpace3dAdaptive.h"
#include "SurfaceParametric.h"
#include "MeshEdge3d.h"
#include "SurfaceAnalytic.h"
#include "SurfacePlane.h"
#include "ControlSpace2dQuadTree.h"
#include "EPSFile.h"
#include "MeshDecompositionUTC.h"
#include "DataStatistics.h"
#include "MeshGrain.h"
#include "MeshViewSet.h"
#include "Metric2dContext.h"
#include "MeshElement.h"
#include "MeshTriangle2d.h"
#include "MeshSpecialRoutinesUTC.h"
#include "MeshSplit3d.h"
#include "ControlSpace2dIdentity.h"
#include "DLeastSquaresFitting.h"
#include "MeshGenerator3dAdapt.h"
#include "ModelMultiSphere.h"
#include "MarchingCubes.h"
#include "ControlSpace3dKdTree.h"

#include "DSegment.h"
#include "DTriangle.h"

// #include "wininet.h"

MeshModel::MeshModel() : m_model_mesh(0), m_model_name("generic") {}

MeshModel::~MeshModel() {
	if(m_model_mesh) delete m_model_mesh;
}

int MeshModel::execute(const string& cmd_line)
{
	istringstream line(cmd_line);
	string command, buffer, arg1;
	line >> command;
//	cerr << "Command: [" << command << "]\n";
	command = MeshData::lowercase(command);
//	cerr << "command: [" << command << "]\n";

	static auto logger_cstat = Logger::getInstance(LOG4CPLUS_TEXT("mesh.stat.cstat"));

#ifdef USE_EXCEPTIONS
	try{
#endif
		int i = 0;
		//int from;
		//double d;

		if(command.compare("help") == 0){
			LOG4CPLUS_INFO(MeshLog::logger_console, "load <fname.xml>");
			LOG4CPLUS_INFO(MeshLog::logger_console, "load_and_prepare <fname.xml>");
			LOG4CPLUS_INFO(MeshLog::logger_console, "clear_cs");
			LOG4CPLUS_INFO(MeshLog::logger_console, "convert_to_xml <fname_from> [fname_to]");
			LOG4CPLUS_INFO(MeshLog::logger_console, "triangulate");
			LOG4CPLUS_INFO(MeshLog::logger_console, "triangulate3");
			LOG4CPLUS_INFO(MeshLog::logger_console, "quality3 [draw] [force_cs]");
			LOG4CPLUS_INFO(MeshLog::logger_console, "gather3");
			LOG4CPLUS_INFO(MeshLog::logger_console, "smoothen3");
			LOG4CPLUS_INFO(MeshLog::logger_console, "smoothen3_opt <mixed|vec|brute>");
			LOG4CPLUS_INFO(MeshLog::logger_console, "store_xml <fname>");
			LOG4CPLUS_INFO(MeshLog::logger_console, "store_grd <fname> [GR_2D_Q4|GR_2D_Q8|GR_2D_Q9|GR_2D_T3|GR_2D_T6|GR_3D_T4|GR_3D_T10]");
			LOG4CPLUS_INFO(MeshLog::logger_console, "store_txt <fname>");
			LOG4CPLUS_INFO(MeshLog::logger_console, "store_off <fname>");
			LOG4CPLUS_INFO(MeshLog::logger_console, "store_abaqus <fname>");
			LOG4CPLUS_INFO(MeshLog::logger_console, "store_surface_points <fname>");
			LOG4CPLUS_INFO(MeshLog::logger_console, "load_utc_bin <fname>");
			LOG4CPLUS_INFO(MeshLog::logger_console, "store_utc_bin <fname>");
			LOG4CPLUS_INFO(MeshLog::logger_console, "test_crack <cut|transform>");
			LOG4CPLUS_INFO(MeshLog::logger_console, "test_adapt_xml <xml_fname1> <xml_fname2> <cmp_mode> <debug_mode>");
			LOG4CPLUS_INFO(MeshLog::logger_console, "test_adapt_bin <bin_fname1> <bin_fname2> <cmp_mode> <debug_mode>");
			LOG4CPLUS_INFO(MeshLog::logger_console, "kdtree-test-cf> <params>");
			LOG4CPLUS_INFO(MeshLog::logger_console, "kdtree-test-ss> <params>");
			LOG4CPLUS_INFO(MeshLog::logger_console, "parse_grain <dir> <output_file>");
			LOG4CPLUS_INFO(MeshLog::logger_console, "load_pasz <fname>");
			LOG4CPLUS_INFO(MeshLog::logger_console, "store_pasz <fname>");
			LOG4CPLUS_INFO(MeshLog::logger_console, "validate3");
			LOG4CPLUS_INFO(MeshLog::logger_console, "exit");
			LOG4CPLUS_INFO(MeshLog::logger_console, "quit");
			#ifdef _OPENMP
			#pragma omp parallel
				{
					int th_id = omp_get_thread_num();
					if(th_id == 0){
						LOG4CPLUS_INFO(MeshLog::logger_console, "** OpenMP threads: " << omp_get_num_threads());
					}
				}
			#endif // _OPENMP
		//---------------------------------------------------
		}else if(command.compare("exit") == 0 || command.compare("quit") == 0){
			return CM_QUIT;
		//---------------------------------------------------
		}else if(command.compare("load_and_prepare") == 0){
			string fname = "";
			line >> fname;
			LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, "MeshName\t" << fname);
			MeshBRep* desc = nullptr;
			int fn_length = (int)fname.length();
			if(fname.rfind(".xml") == fn_length-4)
				desc = new MeshBRepXML;
			if(desc && desc->parseFile(fname) && desc->validate()) {
				MeshContainer3d* domain = desc->createDomainMesh();
				delete desc;

				if(domain){
					m_mesh_modified = true;
					if(m_model_mesh) delete m_model_mesh;
					m_model_mesh = domain;

					setNameFromFilename(fname);

					for(IteratorBoundary2d it = m_model_mesh->getFirstValidBoundary2d(); 
							it.isValid(); it.nextValidBoundary())
					{
						MeshDomainSurface* domain_surface = it.getDomainSurface();
						CS2dPtr ucs_2d = domain_surface->getUserControlSpace();
						MeshDomainVolume* dvolume0 = (MeshDomainVolume*)domain_surface->getBlock(0);
						MeshDomainVolume* dvolume1 = (MeshDomainVolume*)domain_surface->getBlock(1);
						CS3dPtr ucs_3d_0 =
							dvolume0 ? dvolume0->getUserControlSpace() : nullptr;
						CS3dPtr ucs_3d_1 =
							(dvolume1 && (dvolume0 != dvolume1)) ? dvolume1->getUserControlSpace() : nullptr;
						it.getBoundary()->createControlSpace(ucs_2d, ucs_3d_0, ucs_3d_1);
					}

					return CM_OK;
				}else 
					return CM_ERROR_PARSE;
			}else{
				if(desc) delete desc;
				return CM_ERROR_PARSE;
			}
		//---------------------------------------------------
		}else if(command.compare("clear_cs") == 0){
			if(!m_model_mesh){	
				LOG4CPLUS_ERROR(MeshLog::logger_console,   "No model boundary available");	
				return CM_ERROR_NOMESH;	
			}else{
				m_model_mesh->clearCS();
				return CM_OK;
			}
		//---------------------------------------------------
		}else if(command.compare("load") == 0){
			string fname = "";
			line >> fname;
			LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, "MeshName\t" << fname);
			MeshBRep* desc = nullptr;
			int fn_length = (int)fname.length();
			if(fname.rfind(".xml") == fn_length-4)
				desc = new MeshBRepXML;
			if(desc && desc->parseFile(fname) && desc->validate()) {
				MeshContainer3d* domain = desc->createDomainMesh();
				delete desc;
				if(domain){
					m_mesh_modified = true;
					if(m_model_mesh) delete m_model_mesh;
					m_model_mesh = domain;

					setNameFromFilename(fname);

					if (domain->getBlocksCount() > 0) {
						MeshDomainVolume* mdv = (MeshDomainVolume*)domain->getBlockAt(0);
						MeshContainer3d* mesh3d = mdv->getMesh();

						if (mesh3d) {
							MeshViewSet::ClipPlane cp(FVector3d(-0.460159f, 0.771446f, 0.439472f), -0.251687f);
							DataVector<MeshViewSet::ClipPlane> clip_planes(1);
							clip_planes.add(cp);
							auto mesh_set = mesh3d->getViewSetWithVisibleBlocksOnly(&clip_planes);
							SHOW_MESH("xx mesh3d BlocksOnly xx", mesh_set);
						}
					}
					return CM_OK;
				}else 
					return CM_ERROR_PARSE;
			}else{
				if(desc) delete desc;
				return CM_ERROR_PARSE;
			}
		//---------------------------------------------------
		}else if(command.compare("test_multi_sphere") == 0){
			int sphere_count = 1;
			line >> sphere_count;
			ModelMultiSphere mms( sphere_count );
			MeshContainer3dSurface* surf_mesh = MarchingCubes::createSurfaceMesh( mms );
			if( surf_mesh != nullptr ) {
				SHOW_MESH(" marching cubes - surface mesh (sub_id)", surf_mesh->getViewSet(nullptr, TagExtended::TAG_LOCAL_SURFACE_ID) );
				SHOW_MESH(" marching cubes - surface mesh", surf_mesh->getViewSet() );
				delete surf_mesh;
				return CM_OK;
			}else
				return CM_ERROR_RUN;
		//---------------------------------------------------
		}
		else if (command.compare("kdtree-test-ss") == 0) {
			vector<string> params(istream_iterator<string>{line}, istream_iterator<string>{});
			params.insert(params.begin(), "");
			int argc = (int)params.size();
			const char** argv = new const char*[argc];
			for (int i = 0; i < argc; i++)
				argv[i] = params[i].c_str();
			
			int ret = ControlSpace3dKdTree::checkKdTreeSS(0, argc, argv);
			return CM_OK;
			//---------------------------------------------------
		}
		else if (command.compare("kdtree-test-cf") == 0) {
			vector<string> params(istream_iterator<string>{line}, istream_iterator<string>{});
			params.insert(params.begin(), "");
			int argc = (int)params.size();
			const char** argv = new const char*[argc];
			for (int i = 0; i < argc; i++)
				argv[i] = params[i].c_str();

			int ret = ControlSpace3dKdTree::checkKdTreeCF(0, argc, argv);
			return CM_OK;
			//---------------------------------------------------
		}else if(command.compare("convert_to_xml") == 0){
			string fname_from, fname_to, str_threshold;
			line >> fname_from >> fname_to >> str_threshold;
			if(fname_to == "") fname_to = fname_from + ".xml";
			bool result = false;
			if(str_threshold == "")
				result = MeshBRepXML::convertFile(fname_from, fname_to);
			else{
				double threshold;
				if( DEquation::stringToDouble(str_threshold, DEquation::v_auto, threshold) )
					result = MeshBRepXML::convertFile(fname_from, fname_to, threshold);
				else
					result = MeshBRepXML::convertFile(fname_from, fname_to);
			}
			return result ? CM_OK : CM_ERROR_RUN;
		//---------------------------------------------------
		}else if(command.compare("load_pasz") == 0){
			string fname = "";
			line >> fname;
			MeshContainer3d* domain = MeshSpecialRoutinesPaszynski::loadTetrahedralMesh(fname);
			if(domain){
				m_mesh_modified = true;
				if(m_model_mesh) delete m_model_mesh;
				m_model_mesh = domain;
				setNameFromFilename(fname);
				return CM_OK;
			}else{
				LOG4CPLUS_ERROR(MeshLog::logger_console,   "ERROR loading.");
				return CM_ERROR_PARSE;
			}
		//---------------------------------------------------
		}else if(command.compare("store_xml") == 0){
			if(!m_model_mesh){ LOG4CPLUS_ERROR(MeshLog::logger_console,   "Undefined boundary");	return CM_ERROR_NOMESH;	}
			string fname;
			line >> fname;
			if(!line){ LOG4CPLUS_ERROR(MeshLog::logger_console,   "File Name required!"); return CM_ERROR_PARSE; }
			if(MeshBRepXML::storeFile(fname, m_model_mesh)){
				return CM_OK;
			}else CM_ERROR_RUN;
		//---------------------------------------------------
		}else if(command.compare("gather3") == 0){
			if(!m_model_mesh){ LOG4CPLUS_ERROR(MeshLog::logger_console,   "Undefined boundary");	return CM_ERROR_NOMESH;	}
			if(m_model_mesh->gatherVolumeMeshes()){
				m_mesh_modified = true;
				return CM_OK;
			}else
				return CM_ERROR_RUN;
		//---------------------------------------------------
		}else if(command.compare("store_abaqus") == 0){
			if(!m_model_mesh){ LOG4CPLUS_ERROR(MeshLog::logger_console,   "Undefined boundary");	return CM_ERROR_NOMESH;	}
			string fname;
			line >> fname;
			if(!line){ LOG4CPLUS_ERROR(MeshLog::logger_console,   "File Name required!"); return CM_ERROR_PARSE; }
			const MeshDomainVolume* mdv = m_model_mesh->getTotalMDV();
			if(mdv->getMesh()){
				mdv->storeAbaqus(fname, m_model_name);
			}
			int counter = -1;
			for(IteratorMesh3d it = m_model_mesh->getFirstValidMesh3d(); it.isValid(); it.nextValidMesh()){
				ostringstream ostr;
				ostr << "." << counter;
				it.getDomainVolume()->storeAbaqus(fname+ostr.str(), m_model_name+ostr.str());
			}
		//---------------------------------------------------
		}else if(command.compare("store_utc_bin") == 0){
			if(!m_model_mesh){ LOG4CPLUS_ERROR(MeshLog::logger_console,   "Undefined boundary");	return CM_ERROR_NOMESH;	}
			string fname;
			line >> fname;
			if(!line){ LOG4CPLUS_ERROR(MeshLog::logger_console,   "file name required!"); return CM_ERROR_PARSE; }
			// store meshes 3D
			const MeshDomainVolume* mdv = m_model_mesh->getTotalMDV();
			MeshContainer3d* mesh3d = mdv->getMesh();
			if(mesh3d){
				MeshSpecialRoutinesUTC::storeBinFile(fname, mesh3d);
			}
			int counter = -1;
			for(IteratorMesh3d it = m_model_mesh->getFirstValidMesh3d(); it.isValid(); it.nextValidMesh()){
				ostringstream ostr;
				ostr << fname << "." << ++counter;
				MeshSpecialRoutinesUTC::storeBinFile(ostr.str(), it.getMesh());
			}
		//---------------------------------------------------
		}else if(command.compare("validate3") == 0){
			if(!m_model_mesh){ LOG4CPLUS_ERROR(MeshLog::logger_console,   "Undefined boundary");	return CM_ERROR_NOMESH;	}
			// check meshes 3D
			const MeshDomainVolume* mdv = m_model_mesh->getTotalMDV();
			MeshContainer3d* mesh3d = mdv->getMesh();
			if(mesh3d){
				ostringstream ostr;
				ostr << "TotalMesh " << (mesh3d->isValid() ? "ok" : "invalid");
				LOG4CPLUS_INFO(MeshLog::logger_console, ostr.str());
			}
			int counter = -1;
			for(IteratorMesh3d it = m_model_mesh->getFirstValidMesh3d(); it.isValid(); it.nextValidMesh()){
				ostringstream ostr;
				ostr << "Mesh-" << ++counter << " " << (it.getMesh()->isValid() ? "ok" : "invalid");
				LOG4CPLUS_INFO(MeshLog::logger_console, ostr.str());
			}
		//---------------------------------------------------
		}else if(command.compare("load_utc_bin") == 0){
			string fname;
			line >> fname;
			if(!line){ LOG4CPLUS_ERROR(MeshLog::logger_console,   "file name required!"); return CM_ERROR_PARSE; }
			// load mesh 3D
			MeshContainer3d* volume_mesh = MeshSpecialRoutinesUTC::readBinFile(fname);
			if(!volume_mesh){
				LOG4CPLUS_ERROR(MeshLog::logger_console, "Error reading file: " << fname);
				return CM_ERROR_PARSE;
			}
			LOG4CPLUS_INFO(MeshLog::logger_mesh, "Bbox of the mesh: " << volume_mesh->getBoundingBox());
			// prepare model_mesh
			m_mesh_modified = true;
			if(m_model_mesh) delete m_model_mesh;
			m_model_mesh = new MeshContainer3d(3, volume_mesh);
			setNameFromFilename(fname);
		//---------------------------------------------------
		}else if(command.compare("store_pasz") == 0){
			if(!m_model_mesh){ LOG4CPLUS_ERROR(MeshLog::logger_console,   "Undefined boundary");	return CM_ERROR_NOMESH;	}
			string fname;
			line >> fname;
			if(!line){ LOG4CPLUS_ERROR(MeshLog::logger_console,   "File Name required!"); return CM_ERROR_PARSE; }
			// check 3D
			IteratorMesh3d it = m_model_mesh->getFirstValidMesh3d(); 
			if(it.isValid()){
				MeshSpecialRoutinesPaszynski::storeTetrahedralMesh(fname, 
					"1\n\n'Sphere'\n0.d0 0.d0 0.d0\n1.d0\n\n1\n",
					it.getMesh());
				return CM_OK;
			}else{
				LOG4CPLUS_ERROR(MeshLog::logger_console,   "No valid mesh");
				return CM_ERROR_NOMESH;
			}
		//---------------------------------------------------
		}else if(command.compare("store_surface_points") == 0){
			if(!m_model_mesh){	
				LOG4CPLUS_ERROR(MeshLog::logger_console,   "No model boundary available");	
				return CM_ERROR_NOMESH;	
			}
			line >> buffer;
			string fname = line ? buffer : "mesh";
			int bct = m_model_mesh->getBlocksCount();
			for(int i = 0; i < bct; i++){
				MeshDomainVolume* domain_volume = (MeshDomainVolume*)m_model_mesh->getBlockAt(i);
				if(domain_volume->getMesh() == nullptr)
					domain_volume->prepareBoundaryMesh();
				domain_volume->storeSurfacePoints(fname, i++);
			}
		//---------------------------------------------------
		}else if(command.compare("store_off") == 0){
			if(!m_model_mesh){ LOG4CPLUS_ERROR(MeshLog::logger_console,   "Undefined boundary");	return CM_ERROR_NOMESH;	}
			string fname;
			line >> fname;
			if(!line){ LOG4CPLUS_ERROR(MeshLog::logger_console,   "OFF File Name required!"); return CM_ERROR_PARSE; }
			int mesh_nr = 0;
			// check 3D
			int bct = m_model_mesh->getBlocksCount();
			for(int i = 0; i < bct; i++){
				MeshDomainVolume* domain_volume = (MeshDomainVolume*)m_model_mesh->getBlockAt(i);
				MeshContainer3dSurface* smesh = domain_volume->getSurfaceMesh();
				if(smesh == nullptr) {
					domain_volume->prepareBoundaryMesh();
					smesh = domain_volume->getSurfaceMesh();
				}

				if(smesh) {
					if( mesh_nr == 0 )
						smesh->storeOFF( fname );
					else
						smesh->storeOFF( fname + "-" + to_string(mesh_nr) + ".off" );
					++mesh_nr;
				}
			}
		//---------------------------------------------------
		}else if(command.compare("triangulate") == 0){
			if(!m_model_mesh){	
				LOG4CPLUS_ERROR(MeshLog::logger_console,   "No model boundary available");	
				return CM_ERROR_NOMESH;	
			}

			clock_t clock_start = clock();
			MeshGenerator2d::show_prediction = true;
			MeshGenerator2d::autoTriangulate(m_model_mesh, 2);
			double sec = (clock() - clock_start)/(double)CLOCKS_PER_SEC;
			m_mesh_modified = true;

			int nct = 0, tct = 0;
			for(IteratorMesh2d it = m_model_mesh->getFirstValidMesh2d(); it.isValid(); it.nextValidMesh()){
				MeshContainer2d* mesh = it.getMesh();
				nct += mesh->getPointsCount();
				tct += mesh->getElementsCount();
			}
			LOG4CPLUS_INFO(MeshLog::logger_console, "nodes: " << nct);
			LOG4CPLUS_INFO(MeshLog::logger_console, "triangles: " << tct);
			double speed = ((sec > 0.0)?(tct / sec):99999);
			LOG4CPLUS_INFO(MeshLog::logger_console, "Triangulation speed (triangles per sec.): " << speed);
		//---------------------------------------------------
		}else if(command.compare("test") == 0){
			int mode = 0;
			bool debug_on = true;
			line >> mode;
			if(line){
				string txt;
				line >> txt;
				debug_on = (txt.compare("release") != 0);
			}
			return test(mode, debug_on);
		//---------------------------------------------------
		}else if(command.compare("test_adapt_xml") == 0){
			string fname1, fname2;
			line >> fname1 >> fname2;
			int adapt_mode, cmp_mode, debug_mode;
			line >> adapt_mode >> cmp_mode >> debug_mode;
			if(line){
				return testAdaptXML(fname1, fname2, adapt_mode, cmp_mode, debug_mode != 0);
			}
			else return CM_ERROR_PARSE;
		//---------------------------------------------------
		}else if(command.compare("test_adapt_bin") == 0){
			string fname1, fname2;
			line >> fname1 >> fname2;
			int adapt_mode, cmp_mode, debug_mode;
			line >> adapt_mode >> cmp_mode >> debug_mode;
			if(line){
				return testAdaptBin(fname1, fname2, adapt_mode, cmp_mode, debug_mode != 0);
			}
			else return CM_ERROR_PARSE;
		//---------------------------------------------------
		}else if(command.compare("test_crack") == 0){
			m_mesh_modified = true;
			string str_mode;
			line >> str_mode;
			int mode = 0;
			if(line && str_mode == "transform") mode = 1;
			return testCrack(mode);
		//---------------------------------------------------
		}else if(command.compare("parse_grain") == 0){
			string dir = "", output_file = "grain.msh";
			line >> dir >> output_file;
			setNameFromFilename(dir);
			MeshGrain::parseGrainFiles(dir, output_file);
		//---------------------------------------------------
		}else if(command.compare("triangulate3") == 0){
			if(!m_model_mesh){	LOG4CPLUS_ERROR(MeshLog::logger_console,   "Undefined boundary");	return CM_ERROR_NOMESH;	}
			clock_t clock_start = clock();
			MeshGenerator2d::show_prediction = false;
			START_CLOCK("3D Triangulation");
			MeshGenerator3d::autoTriangulate(m_model_mesh, 2);
			STOP_CLOCK("3D Triangulation");
			double sec = (clock() - clock_start)/(double)CLOCKS_PER_SEC;
			m_mesh_modified = true;

			int nct = 0, ect = 0;
			for(IteratorMesh3d it = m_model_mesh->getFirstValidMesh3d(); it.isValid(); it.nextValidMesh()){
				MeshContainer3d* mesh = it.getMesh();
				nct += mesh->getPointsCount();
				ect += mesh->getBlocksCount();
			}
			LOG4CPLUS_INFO(logger_cstat, "NP\t" << nct);
			LOG4CPLUS_INFO(logger_cstat, "NT\t" << ect);
			LOG4CPLUS_INFO(logger_cstat, "Meshing3DTime-total\t" << sec);
			LOG4CPLUS_INFO(MeshLog::logger_console, "Generation speed (tetrahedra per sec.): " << ect / sec);
		//---------------------------------------------------
		}else if(command.compare("quality3") == 0){
			if(!m_model_mesh){	
				LOG4CPLUS_ERROR(MeshLog::logger_console,   
					"Undefined model");	
				return CM_ERROR_NOMESH;	
			}
			const MeshDomainVolume* mdv = m_model_mesh->getTotalMDV();
			MeshContainer3d* mesh3d = mdv->getMesh();
			if(!mesh3d){ 
				LOG4CPLUS_ERROR(MeshLog::logger_console,   
					"No tetrahedral mesh available");	
				return CM_ERROR_NOMESH;	
			}

			string str_mode;
			int draw_mode = 0;
			int force_cs_mode = 0;
			while(true){
				line >> str_mode;
				if(!line) break;
				if(str_mode == "draw") draw_mode = 1;
				else if(str_mode == "force_cs") force_cs_mode = 1;
			}

			DataStatistics mean_ratio_stats;
			DataStatistics edge_stats;
			DataStatistics nonconf_stats;
			DataStatistics edges_all_adjacency_stats;
			DataStatistics edges_inner_adjacency_stats;
			DataStatistics vertices_all_adjacency_stats;
			DataStatistics vertices_inner_adjacency_stats;
			DataStatistics dihedral_angles_stats;
			unsigned int np = mesh3d->getPointsCount();
			unsigned int nt = mesh3d->getBlocksCount();
				
			CS3dPtr space = mesh3d->getControlSpace();

			if(!space && (force_cs_mode == 1)){
				LOG4CPLUS_WARN(MeshLog::logger_console, "Recreating ACS from mesh");
				space = MeshGenerator3dAdapt::createACSfromMeshBlocks(mesh3d);
				mesh3d->setControlSpace(space);
				LOG4CPLUS_INFO(MeshLog::logger_console, "Done...");
			}

			if(space){
				Metric3dContext mc(space);
				mesh3d->statMeanRatio(mc, mean_ratio_stats, false);
				mesh3d->statMetricEdgeLength(mc, edge_stats, false);
				mesh3d->statMetricDifference(mc, nonconf_stats);
			}
			mesh3d->statEdgeBlockAdjacency(edges_all_adjacency_stats, false);
			mesh3d->statEdgeBlockAdjacency(edges_inner_adjacency_stats, true);
			mesh3d->statVertexBlockAdjacency(vertices_all_adjacency_stats, false);
			mesh3d->statVertexBlockAdjacency(vertices_inner_adjacency_stats, true);
			mesh3d->statMinDihedralAngles(dihedral_angles_stats);

			LOG4CPLUS_INFO(logger_cstat, "NP\t" << np);
			LOG4CPLUS_INFO(logger_cstat, "NT\t" << nt);
			if(mean_ratio_stats.calculate()){
				LOG4CPLUS_INFO(logger_cstat, "MetricMeanRatio-ave\t" << mean_ratio_stats.average());
				LOG4CPLUS_INFO(logger_cstat, "MetricMeanRatio-min\t" << mean_ratio_stats.minimum());
				LOG4CPLUS_INFO(logger_cstat, "MetricMeanRatio-dev\t" << mean_ratio_stats.stdDev());
			}
			if(edge_stats.calculate()){
				LOG4CPLUS_INFO(logger_cstat, "MetricEdgeLength-ave\t" << edge_stats.average());
				LOG4CPLUS_INFO(logger_cstat, "MetricEdgeLength-dev\t" << edge_stats.stdDev());
			}				
			if(nonconf_stats.calculate()){
				LOG4CPLUS_INFO(logger_cstat, "MetricDiff-ave\t" << nonconf_stats.average());
				LOG4CPLUS_INFO(logger_cstat, "MetricDiff-dev\t" << nonconf_stats.stdDev());
			}
			if(edges_all_adjacency_stats.calculate()){
				LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, "Edges-all\t" << edges_all_adjacency_stats.countInt());
				LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, "Edges-Block-adjacency-ave\t" << edges_all_adjacency_stats.average());
				LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, "Edges-Block-adjacency-max\t" << edges_all_adjacency_stats.maximum());
				LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, "Edges-Block-adjacency-min\t" << edges_all_adjacency_stats.minimum());
			}
			if(edges_inner_adjacency_stats.calculate()){
				LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, "Edges-inner\t" << edges_inner_adjacency_stats.countInt());
				LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, "Edges-Block-adjacency-ave\t" << edges_inner_adjacency_stats.average());
				LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, "Edges-Block-adjacency-max\t" << edges_inner_adjacency_stats.maximum());
				LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, "Edges-Block-adjacency-min\t" << edges_inner_adjacency_stats.minimum());
			}
			if(vertices_all_adjacency_stats.calculate()){
				LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, "Vertices-all\t" << vertices_all_adjacency_stats.countInt());
				LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, "Vertices-Block-adjacency-ave\t" << vertices_all_adjacency_stats.average());
				LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, "Vertices-Block-adjacency-max\t" << vertices_all_adjacency_stats.maximum());
				LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, "Vertices-Block-adjacency-min\t" << vertices_all_adjacency_stats.minimum());
			}
			if(vertices_inner_adjacency_stats.calculate()){
				LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, "Vertices-inner\t" << vertices_inner_adjacency_stats.countInt());
				LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, "Vertices-Block-adjacency-ave\t" << vertices_inner_adjacency_stats.average());
				LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, "Vertices-Block-adjacency-max\t" << vertices_inner_adjacency_stats.maximum());
				LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, "Vertices-Block-adjacency-min\t" << vertices_inner_adjacency_stats.minimum());
			}
			if(dihedral_angles_stats.calculate()){
				double angles_min = dihedral_angles_stats.minimum() * 180.0 / PI;
				int angles_05 = dihedral_angles_stats.getDataCountBottom(PI/36);
				int angles_10 = dihedral_angles_stats.getDataCountBottom(PI/18);
				int ct = dihedral_angles_stats.countInt();
				LOG4CPLUS_INFO(logger_cstat, "Angles-min\t" << angles_min);
				LOG4CPLUS_INFO(logger_cstat, 
					"Angles-below-05o\t" << angles_05 << "\t(" << (100.0*angles_05 / ct) << "%)");
				LOG4CPLUS_INFO(logger_cstat,
					"Angles-below-10o\t" << angles_10 << "\t(" << (100.0*angles_10 / ct) << "%)");
			}

			if(draw_mode == 1){
				MeshViewSet* set = new MeshViewSet(0, np, np, nt/10);
				double qt_min = sin(PI/36);
				for(int i = 0; i < (int)nt; i++){
					MeshTetrahedron* tetra = (MeshTetrahedron*)mesh3d->getBlockAt(i);
					double qt = tetra->getMinDihedralAngleSin();
					if(qt < qt_min) set->addBlockWithEdges(tetra);
				}
				for(IteratorEdge3d it = mesh3d->getFirstEdge3d(); it.isValid(); it.nextEdge()){
					if(it.getEdge()->isBorder(TagBorder::RIDGE | TagBorder::FIXED))
						set->addEdge(it.getEdge());
				}
				SHOW_MESH("Tetrahedra with dihedral angle < 5o", set);
			}
		//---------------------------------------------------
		}else if(command.compare("smoothen3_opt") == 0){
			if(!m_model_mesh){	LOG4CPLUS_ERROR(MeshLog::logger_console,   "Undefined boundary");	return CM_ERROR_NOMESH;	}
			const MeshDomainVolume* mdv = m_model_mesh->getTotalMDV();
			MeshContainer3d* mesh3d = mdv->getMesh();
			if(!mesh3d){ LOG4CPLUS_ERROR(MeshLog::logger_console,   "No tetrahedral mesh available");	return CM_ERROR_NOMESH;	}

			string str_mode;
			line >> str_mode;
			int mode = 0; // mixed
			if(line){
				if(str_mode == "vec") mode = 1;
				else if(str_mode == "brute") mode = 2;
			}

			if(mode == 0 || mode == 1){
				CS3dPtr space = mesh3d->getControlSpace();
				if(!space){
					LOG4CPLUS_WARN(MeshLog::logger_console, "Recreating ACS from mesh");
					space = MeshGenerator3dAdapt::createACSfromMeshBlocks(mesh3d);
					mesh3d->setControlSpace(space);
					LOG4CPLUS_INFO(MeshLog::logger_console, "Smoothing...");
				}
				Metric3dContext mc(space);
				if(mode == 0)
					MeshGenerator3dQuality::smoothenOptMixed(mc, mesh3d);
				else if(mode == 1)
					MeshGenerator3dQuality::smoothenOpt(mc, mesh3d);
			}else if(mode == 2){
				MeshGenerator3dQuality::optimizeTetrahedraForMinAngle(mesh3d);
			}
		//---------------------------------------------------
		}else if(command.compare("smoothen3") == 0){
			if(!m_model_mesh){	LOG4CPLUS_ERROR(MeshLog::logger_console,   "Undefined boundary");	return CM_ERROR_NOMESH;	}
			const MeshDomainVolume* mdv = m_model_mesh->getTotalMDV();
			MeshContainer3d* mesh3d = mdv->getMesh();
			if(!mesh3d){ LOG4CPLUS_ERROR(MeshLog::logger_console,   "No tetrahedral mesh available");	return CM_ERROR_NOMESH;	}

			CS3dPtr space = mesh3d->getControlSpace();
			if(!space){
				LOG4CPLUS_WARN(MeshLog::logger_console, "Recreating ACS from mesh");
				space = MeshGenerator3dAdapt::createACSfromMeshBlocks(mesh3d);
				mesh3d->setControlSpace(space);
			}
			Metric3dContext mc(space);

			MeshGenerator3dQuality::smoothen(mc, mesh3d, 1, TagExtended::TAG_NONE, 1, true, 
				MeshData::SM3_OPT_MIXED | MeshData::SM3_SWAP_COMPLETE | MeshData::SM3_BORDER_PRUNE);
		//---------------------------------------------------
		}else if(command.compare("store_txt") == 0){
			if(!m_model_mesh){ LOG4CPLUS_ERROR(MeshLog::logger_console,   "Undefined boundary");	return CM_ERROR_NOMESH;	}
			string fname;
			line >> fname;
			if(!line){ LOG4CPLUS_ERROR(MeshLog::logger_console,   "Grid File Name required!"); return CM_ERROR_PARSE; }
			int mesh_nr = 0;
			// check 3D
			for(IteratorMesh3d it = m_model_mesh->getFirstValidMesh3d(); it.isValid(); it.nextValidMesh()){
				it.getMesh()->storeTxt(fname, ++mesh_nr);
			}
			if(mesh_nr == 0){ // no 3D volume mesh
				for(IteratorMesh2d it = m_model_mesh->getFirstValidMesh2d(); it.isValid(); it.nextValidMesh()){
					//SurfaceConstPtr surface = it.getDomainSurface()->getBaseSurface();
					it.getMesh()->storeTxt(fname, ++mesh_nr);
				}
			}
		//---------------------------------------------------
		}else if(command.compare("store_grd") == 0){
			if(!m_model_mesh){ LOG4CPLUS_ERROR(MeshLog::logger_console,   "Undefined boundary");	return CM_ERROR_NOMESH;	}
			string fname;
			line >> fname;
			if(!line){ LOG4CPLUS_ERROR(MeshLog::logger_console,   "Grid File Name required!"); return CM_ERROR_PARSE; }
			int grid_type = 0;
			line >> arg1;
			if(line){ 
				for(size_t k = 0; k < arg1.length(); k++)
					arg1[k] = tolower(arg1[k]);
				if(arg1.compare("gr_2d_q4")==0) grid_type = 1;
				else if(arg1.compare("gr_2d_q8")==0) grid_type = 3;
				else if(arg1.compare("gr_2d_q9")==0) grid_type = 5;
				else if(arg1.compare("gr_2d_t3")==0) grid_type = 7;
				else if(arg1.compare("gr_2d_t6")==0) grid_type = 8;
				else if(arg1.compare("gr_3d_t4")==0) grid_type = 25;
				else if(arg1.compare("gr_3d_t10")==0) grid_type = 26;
				else {LOG4CPLUS_ERROR(MeshLog::logger_console,   "Unknown Grid Type!"); return CM_ERROR_PARSE; }
			}
			if(grid_type == 0){ // auto-recognize
				// check 3D
				for(IteratorMesh3d it = m_model_mesh->getFirstValidMesh3d(); it.isValid(); it.nextValidMesh()){
					grid_type = 25;
					it.getMesh()->storePJM(fname, grid_type, i++);
				}
				if(grid_type == 0){ // still, i.e. no 3D mesh
					for(IteratorMesh2d it = m_model_mesh->getFirstValidMesh2d(); it.isValid(); it.nextValidMesh()){
						//SurfaceConstPtr surface = it.getDomainSurface()->getBaseSurface();
						grid_type = (it.getMesh()->getElementsCount(4) > 0) ? 1 : 7; // linear for now
						it.getMesh()->storePJM(fname, grid_type, i++);
					}
				}
			}else if(grid_type < 25){ // 2d
				for(IteratorMesh2d it = m_model_mesh->getFirstValidMesh2d(); it.isValid(); it.nextValidMesh()){
					//SurfaceConstPtr surface = it.getDomainSurface()->getBaseSurface();
					it.getMesh()->storePJM(fname, grid_type, i++);
				}
			}else{
				for(IteratorMesh3d it = m_model_mesh->getFirstValidMesh3d(); it.isValid(); it.nextValidMesh()){
					it.getMesh()->storePJM(fname, grid_type, i++);
				}
			}
		//---------------------------------------------------
		}else if(command.compare("set") == 0){
			string name;
			line >> name;
			auto p = mesh_data.getProperty(name);
			if(p){
				if(p->isInt()){
					line >> i;
					if(line){
						p->setIntValue(i);
						m_properties_modified = true;
					}else{
						LOG4CPLUS_ERROR(MeshLog::logger_console, 
							"ERROR reading int value for " << name << "(" << p->description << ")");
						LOG4CPLUS_INFO(MeshLog::logger_console, "current value: " << p->getIntValue());
					}
				}else{
					double d = 0.0;
					line >> d;
					if(line){
						p->setDoubleValue(d);
						m_properties_modified = true;
					}else{
						LOG4CPLUS_ERROR(MeshLog::logger_console,
							"ERROR reading double value for " << name << "(" << p->description << ")");
						LOG4CPLUS_INFO(MeshLog::logger_console, "current value: " << p->getDoubleValue() );
					}
				}
			}else{
				LOG4CPLUS_ERROR(MeshLog::logger_console,   "ERROR unknown property name");
			}
/*
			LOG4CPLUS_INFO(MeshLog::logger_console, "CUT_MESH <fname>");
//			LOG4CPLUS_INFO(MeshLog::logger_console, "EDGE_INNER_NODES_2 <edge_node_ct> [nonlinear]");
//			LOG4CPLUS_INFO(MeshLog::logger_console, "EDGE_INNER_NODES_3 <edge_node_ct> [nonlinear]");
			LOG4CPLUS_INFO(MeshLog::logger_console, "EPS");
			LOG4CPLUS_INFO(MeshLog::logger_console, "GEN1D [pts_ct=-1] [real3d=1]");
			LOG4CPLUS_INFO(MeshLog::logger_console, "GEN1.5D");
			LOG4CPLUS_INFO(MeshLog::logger_console, "GEN2D");
			LOG4CPLUS_INFO(MeshLog::logger_console, "GEN2.5D");
			LOG4CPLUS_INFO(MeshLog::logger_console, "GEN3D");
			LOG4CPLUS_INFO(MeshLog::logger_console, "GET [property]");
			LOG4CPLUS_INFO(MeshLog::logger_console, "LEELO [max_ct=-1]");
			LOG4CPLUS_INFO(MeshLog::logger_console, "LOAD_AND_PREPARE <mesh_file_name>");
			LOG4CPLUS_INFO(MeshLog::logger_console, "LOAD_SURF <file_name>");
			LOG4CPLUS_INFO(MeshLog::logger_console, "LOG_NAME <file_name>");
			LOG4CPLUS_INFO(MeshLog::logger_console, "MAKE_QUADS <qmorph/leelo/mixed>");
			LOG4CPLUS_INFO(MeshLog::logger_console, "MAKE_QUADS_FOR_PATCH <file.resu.surf> [qmorph/leelo/mixed]");
			LOG4CPLUS_INFO(MeshLog::logger_console, "PAUSE [info]");
			LOG4CPLUS_INFO(MeshLog::logger_console, "QMORPH [max_ct=-1]");
			LOG4CPLUS_INFO(MeshLog::logger_console, "QUALITY");
			LOG4CPLUS_INFO(MeshLog::logger_console, "QUALITY_QUADS");
			LOG4CPLUS_INFO(MeshLog::logger_console, "QUALITY3");
			LOG4CPLUS_INFO(MeshLog::logger_console, "QUALITY3_DETAILED");
			LOG4CPLUS_INFO(MeshLog::logger_console, "QUALITY_EXTENDED");
			LOG4CPLUS_INFO(MeshLog::logger_console, "SET <property> <value>");
			LOG4CPLUS_INFO(MeshLog::logger_console, "SMOOTHEN [ct=1]");
			LOG4CPLUS_INFO(MeshLog::logger_console, "SMOOTHEN_DEL_SWAP [ct=1]");
			LOG4CPLUS_INFO(MeshLog::logger_console, "SMOOTHEN_QUADS [ct=1]");
			LOG4CPLUS_INFO(MeshLog::logger_console, "SMOOTHEN_QUADS_3D [ct=1]");
			LOG4CPLUS_INFO(MeshLog::logger_console, "SMOOTHEN_METRIC [ct=1]");
			LOG4CPLUS_INFO(MeshLog::logger_console, "SMOOTHEN_LAPLACE_METRIC [ct=1]");
			LOG4CPLUS_INFO(MeshLog::logger_console, "SPECIAL_PHASE_RUN [hesjan=0] [from=0]");
			LOG4CPLUS_INFO(MeshLog::logger_console, "SPLIT");
			LOG4CPLUS_INFO(MeshLog::logger_console, "STATS");
			LOG4CPLUS_INFO(MeshLog::logger_console, "STORE2.5 <fname>");
			LOG4CPLUS_INFO(MeshLog::logger_console, "STORE2 <fname>");
			LOG4CPLUS_INFO(MeshLog::logger_console, "STORE_AMIRA <fname>");
			LOG4CPLUS_INFO(MeshLog::logger_console, "STORE_SURF <fname>");
			LOG4CPLUS_INFO(MeshLog::logger_console, "STORE_DAT <fname>");
			LOG4CPLUS_INFO(MeshLog::logger_console, "TEST_HESSIAN <fname> <equation> [equation2]");
			LOG4CPLUS_INFO(MeshLog::logger_console, "TEST_HESSIAN_QUAD <fname> <equation> [equation2]");
			LOG4CPLUS_INFO(MeshLog::logger_console, "TRIANGULATE3 [sm_count]");
			LOG4CPLUS_INFO(MeshLog::logger_console, "TRIANGULATE3_SURFACE <fname>");
		//---------------------------------------------------
		}else if(command.compare("gen1d") == 0){
			if(!boundary){
				LOG4CPLUS_ERROR(MeshLog::logger_console,   "Undefined boundary");
				return CM_ERROR_NOMESH;
			}
			MeshGenerator1d::discretizeEdgesMin(boundary);
		//---------------------------------------------------
		}else if(command.compare("special_phase_run") == 0){
			str_text >> i;	if(!str_text) i = 0;
			str_text >> from;	if(!str_text) from = 0;
			MeshSpecialRoutinesDAT::runSpecialGeometryPhase(i-1, from);
		//---------------------------------------------------
		}else if(command.compare("load_surf") == 0){
			string fname = "";
			str_text >> fname;
			MeshContainer3d* domain = MeshSpecialRoutinesDAT::loadSurfaceBoundaryMesh(fname.c_str());
			if(domain){
				if(boundary) delete boundary;
				boundary = domain;
			}else
				LOG4CPLUS_ERROR(MeshLog::logger_console,   "ERROR loading.");
		//---------------------------------------------------
		}else if(command.compare("pause") == 0){
			string info = "";
			str_text >> info;
			LOG4CPLUS_INFO(MeshLog::logger_console, "pause, waiting for input...",info);
			cin >> info;
		//---------------------------------------------------
		}else if(command.compare("log_name") == 0){
			string fname = "";
			str_text >> fname;
			MeshLog::setLogFile(fname.c_str());
		//---------------------------------------------------
		}else if(command.compare("edge_inner_nodes_2") == 0){
			if(!boundary){	LOG4CPLUS_ERROR(MeshLog::logger_console,   "Undefined boundary");	return CM_ERROR_NOMESH;	}
			str_text >> i;
			if(!str_text) i = 0;
			bool use_shapes = false;
			str_text >> arg1;
			if(str_text && (arg1.compare("nonlinear") == 0)) use_shapes = true;
			for(IteratorMesh2d it = boundary->getFirstValidMesh2d(); it.isValid(); it.nextValidMesh()){
				Metric2dContext mc(it.getMesh()->getControlSpace());
				it.getMesh()->setEdgeInnerPoints(mc, i, use_shapes);
			}
		//---------------------------------------------------
		}else if(command.compare("edge_inner_nodes_3") == 0){
			if(!boundary){	LOG4CPLUS_ERROR(MeshLog::logger_console,   "Undefined boundary");	return CM_ERROR_NOMESH;	}
			str_text >> i;
			if(!str_text) i = 0;
			bool use_shapes = false;
			str_text >> arg1;
			if(str_text && (arg1.compare("nonlinear") == 0)) use_shapes = true;
			for(IteratorMesh3d it = boundary->getFirstValidMesh3d(); it.isValid(); it.nextValidMesh())
				it.getMesh()->addEdgeInnerPoints(i, use_shapes);
		//---------------------------------------------------
		}else if(command.compare("test_hessian") == 0){
			str_text >> buffer;	//
			if(!str_text){
				LOG4CPLUS_ERROR(MeshLog::logger_console,   "mesh file name required.");
				return CM_ERROR_PARSE;
			}
			string eq1, eq2 = "";
			str_text >> eq1;	// "2*sin(3*x)"
			if(!str_text){
				LOG4CPLUS_ERROR(MeshLog::logger_console,   "equation string required.");
				return CM_ERROR_PARSE;
			}
			str_text >> eq2;	// "2*sin(3*x)"

			MeshSpecialRoutinesDAT::testHessianCurvature(buffer, eq1, eq2);
		//---------------------------------------------------
		}else if(command.compare("test_hessian_quad") == 0){
			str_text >> buffer;	//
			if(!str_text){
				LOG4CPLUS_ERROR(MeshLog::logger_console,   "mesh file name required.");
				return CM_ERROR_PARSE;
			}
			string eq1, eq2 = "";
			str_text >> eq1;	// "2*sin(3*x)"
			if(!str_text){
				LOG4CPLUS_ERROR(MeshLog::logger_console,   "equation string required.");
				return CM_ERROR_PARSE;
			}
			str_text >> eq2;	// "2*sin(3*x)"

			MeshSpecialRoutinesDAT::testHessianCurvature(buffer, eq1, eq2, true);
		//---------------------------------------------------
		}else if(command.compare("cut_mesh") == 0){
			str_text >> buffer;	//
			if(!str_text){
				LOG4CPLUS_ERROR(MeshLog::logger_console,   "script file name required.");
				return CM_ERROR_PARSE;
			}
			MeshDecompositionUTC cut_task(buffer);
			cut_task.run();
		//---------------------------------------------------
		}else if(command.compare("store_dat") == 0){
			if(!boundary){	LOG4CPLUS_ERROR(MeshLog::logger_console,   "Undefined boundary");	return CM_ERROR_NOMESH;	}
			str_text >> buffer;
			string fname = str_text ? buffer : "mesh";
			for(IteratorMesh2d it = boundary->getFirstValidMesh2d(); it.isValid(); it.nextValidMesh())
				MeshSpecialRoutinesDAT::exportToDAT(it.getMesh(), fname, true);
		//---------------------------------------------------
		}else if(command.compare("store2") == 0){
			if(!boundary){	LOG4CPLUS_ERROR(MeshLog::logger_console,   "Undefined boundary");	return CM_ERROR_NOMESH;	}
			str_text >> buffer;
			string fname = str_text ? buffer : "mesh";
			for(IteratorMesh2d it = boundary->getFirstValidMesh2d(); it.isValid(); it.nextValidMesh()){
				it.getMesh()->storeTxt(fname, i++);
			}
		//---------------------------------------------------
		}else if(command.compare("store_surf") == 0){
			if(!boundary){	LOG4CPLUS_ERROR(MeshLog::logger_console,   "Undefined boundary");	return CM_ERROR_NOMESH;	}
			str_text >> buffer;
			string fname = str_text ? buffer : "mesh-surf";
			for(IteratorMesh2d it = boundary->getFirstValidMesh2d(); it.isValid(); it.nextValidMesh()){
				bool proper_orientation = it.getDomainVolume()->properOrientation((MeshFace*)it.getDomainSurface());
				it.getMesh()->storeSurf(fname, proper_orientation, i++);
			}
		//---------------------------------------------------
		}else if(command.compare("store_amira") == 0){
			if(!boundary){	LOG4CPLUS_ERROR(MeshLog::logger_console,   "Undefined boundary");	return CM_ERROR_NOMESH;	}
			str_text >> buffer;
			string fname = str_text ? buffer : "mesh";
			if(!boundary->getFirstValidBoundaryMesh3d().isValid())
				MeshGenerator3d::prepareBoundaryMesh(boundary);

			for(IteratorBoundaryMesh3d it = boundary->getFirstValidBoundaryMesh3d(); it.isValid(); it.nextValidBoundaryMesh())
				it.getDomainVolume()->storeAmiraMesh(fname, i++);
		//---------------------------------------------------
		}else if(command.compare("triangulate3_surface") == 0){

			if(boundary) delete[] boundary;
			string fname;
			str_text >> fname;
			LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, "MeshName\t" << fname);
			clock_t clock_start = clock();
			MeshGenerator2d::show_prediction = false;
			//MeshGenerator3d::prepareBoundaryMesh(boundary);
			boundary = MeshSpecialRoutinesDAT::loadSurfaceBoundaryMesh(fname);
			if(!boundary) { LOG4CPLUS_ERROR(MeshLog::logger_console,   "Error loading boundary mesh"); return CM_ERROR_NOMESH; }
			MeshGenerator3d::createTetrahedra(boundary);
			MeshGenerator3dQuality::smoothenBlocks(boundary, 2);
			double sec = (clock() - clock_start)/(double)CLOCKS_PER_SEC;
			int nct = 0, ect = 0;
			for(IteratorMesh3d it = boundary->getFirstValidMesh3d(); it.isValid(); it.nextValidMesh()){
				MeshContainer3d* mesh = it.getMesh();
				nct += mesh->getPointsCount();
				ect += mesh->getBlocksCount();
			}
			LOG4CPLUS_INFO(MeshLog::logger_console, "nodes", nct);
			LOG4CPLUS_INFO(MeshLog::logger_console, "tetrahedra", ect);
			LOG4CPLUS_INFO(MeshLog::logger_console, "Generation time - 3D (sec.)", sec);
			LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, "Meshing3DTime-total\t" << sec);
			LOG4CPLUS_INFO(MeshLog::logger_console, "Generation speed (tetrahedra per sec.)", ect / sec);
		//---------------------------------------------------
		}else if(command.compare("smoothen") == 0){
			str_text >> i;	if(!str_text) i = 1;
			clock_t clock_start = clock();
			MeshGenerator2d::smoothenFaces(boundary, i,
				MeshData::SM_LAPLACE_MIXED | MeshData::SM_TOP_SWAP | MeshData::SM_DEL_SWAP_COMPLETE);
			double sec = (clock() - clock_start)/(double)CLOCKS_PER_SEC;
			LOG4CPLUS_INFO(MeshLog::logger_console, "Smoothing time (sec.)", sec);
		//---------------------------------------------------
		}else if(command.compare("smoothen_metric") == 0){
			str_text >> i;	if(!str_text) i = 1;
			clock_t clock_start = clock();
			MeshGenerator2d::smoothenFaces(boundary, i, MeshData::SM_METRIC);
			double sec = (clock() - clock_start)/(double)CLOCKS_PER_SEC;
			LOG4CPLUS_INFO(MeshLog::logger_console, "Smoothing time (sec.)", sec);
		//---------------------------------------------------
		}else if(command.compare("smoothen3") == 0){
			str_text >> i;	if(!str_text) i = 1;
			MeshGenerator3dQuality::smoothenBlocks(boundary, i);
		//---------------------------------------------------
		}else if(command.compare("smoothen_del_swap") == 0){
			str_text >> i;	if(!str_text) i = 1;
			if(i > 0) MeshGenerator2d::smoothenFaces(boundary, i, MeshData::SM_DEL_SWAP);
			else MeshGenerator2d::smoothenFaces(boundary, 1, MeshData::SM_DEL_SWAP_COMPLETE);
		//---------------------------------------------------
		}else if(command.compare("smoothen_laplace_metric") == 0){
			str_text >> i;	if(!str_text) i = 1;
			MeshGenerator2d::smoothenFaces(boundary, i, MeshData::SM_LAPLACE_METRIC);
		//---------------------------------------------------
		}else if(command.compare("smoothen_quads") == 0){
			str_text >> i;	if(!str_text) i = 1;
			MeshGenerator2dQuad::smoothenFaces(boundary, i, MeshData::SM_TOP_SWAP | MeshData::SM_LAPLACE);
		//---------------------------------------------------
		}else if(command.compare("quality") == 0){
			if(!boundary){	LOG4CPLUS_ERROR(MeshLog::logger_console,   "Undefined boundary");	return CM_ERROR_NOMESH;	}
			DataStatistics mean_ratio_stats;
			DataStatistics edge_stats;
			DataStatistics nonconf_stats;
			for(IteratorMesh2d it = boundary->getFirstValidMesh2d(); it.isValid(); it.nextValidMesh()){
				MeshContainer2d* mesh = it.getMesh();
				CS2dPtr space = mesh->getControlSpace();
				Metric2dContext mc(space);
				mesh->statMeanRatio(mc, mean_ratio_stats, false);
				mesh->statMetricEdgeLength(mc, edge_stats, false);
				mesh->statMetricDifference(mc, nonconf_stats);
			}
			if(mean_ratio_stats.calculate()){
				LOG4CPLUS_INFO(MeshLog::logger_console, "Metric mean ratio 1p (ave)", mean_ratio_stats.average());
				LOG4CPLUS_INFO(MeshLog::logger_console, "-------------------- (dev)", mean_ratio_stats.stdDev());
			}
			if(edge_stats.calculate()){
				LOG4CPLUS_INFO(MeshLog::logger_console, "Metric edge lengths 1p (ave)", edge_stats.average());
				LOG4CPLUS_INFO(MeshLog::logger_console, "---------------------- (dev)", edge_stats.stdDev());
			}				
			if(nonconf_stats.calculate()){
				LOG4CPLUS_INFO(MeshLog::logger_console, "Metric difference 1p (ave)", nonconf_stats.average());
				LOG4CPLUS_INFO(MeshLog::logger_console, "-------------------- (dev)", nonconf_stats.stdDev());
			}				
		//---------------------------------------------------
		}else if(command.compare("quality3") == 0){
			if(!boundary){	LOG4CPLUS_ERROR(MeshLog::logger_console,   "Undefined boundary");	return CM_ERROR_NOMESH;	}
			DataStatistics mean_ratio_stats;
			DataStatistics edge_stats;
			DataStatistics nonconf_stats;
			unsigned int np = 0;
			unsigned int nt = 0;
			for(IteratorMesh3d it = boundary->getFirstValidMesh3d(); it.isValid(); it.nextValidMesh()){
				MeshContainer3d* mesh = it.getMesh();
				np += mesh->getPointsCount();
				nt += mesh->getBlocksCount();
				ControlSpace3d* space = mesh->getControlSpace();
				Metric3dContext mc(space);
				mesh->statMeanRatio(mc, mean_ratio_stats, false);
				mesh->statMetricEdgeLength(mc, edge_stats, false);
				mesh->statMetricDifference(mc, nonconf_stats);
			}
			LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, "NP\t" << np);
			LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, "NT\t" << nt);
			if(mean_ratio_stats.calculate()){
				LOG4CPLUS_INFO(MeshLog::logger_console, "Metric mean ratio (ave)", mean_ratio_stats.average());
				LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, "MetricMeanRatio-ave\t" << mean_ratio_stats.average());
				LOG4CPLUS_INFO(MeshLog::logger_console, "Metric mean ratio (min)", mean_ratio_stats.minimum());
				LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, "MetricMeanRatio-min\t" << mean_ratio_stats.minimum());
				LOG4CPLUS_INFO(MeshLog::logger_console, "------------------(dev)", mean_ratio_stats.stdDev());
				LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, "MetricMeanRatio-dev\t" << mean_ratio_stats.stdDev());
			}
			if(edge_stats.calculate()){
				LOG4CPLUS_INFO(MeshLog::logger_console, "Metric edge lengths (ave)", edge_stats.average());
				LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, "MetricEdgeLength-ave\t" << edge_stats.average());
				LOG4CPLUS_INFO(MeshLog::logger_console, "--------------------(dev)", edge_stats.stdDev());
				LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, "MetricEdgeLength-dev\t" << edge_stats.stdDev());
			}				
			if(nonconf_stats.calculate()){
				LOG4CPLUS_INFO(MeshLog::logger_console, "Metric difference  (ave)", nonconf_stats.average());
				LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, "MetricDiff-ave\t" << nonconf_stats.average());
				LOG4CPLUS_INFO(MeshLog::logger_console, "------------------ (dev)", nonconf_stats.stdDev());
				LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, "MetricDiff-dev\t" << nonconf_stats.stdDev());
			}				
		//---------------------------------------------------
		}else if(command.compare("quality_quads") == 0){
			if(!boundary){	LOG4CPLUS_ERROR(MeshLog::logger_console,   "Undefined boundary");	return CM_ERROR_NOMESH;	}
			DataStatistics angle_stats;
			DataStatistics edge_stats;
			DataStatistics alpha_stats;
			for(IteratorMesh2d it = boundary->getFirstValidMesh2d(); it.isValid(); it.nextValidMesh()){
				MeshContainer2d* mesh = it.getMesh();
				ControlSpace2dMesh* space = (ControlSpace2dMesh*)mesh->getControlSpace();
				SurfaceParametric* surface = it.getDomainSurface()->getBaseSurface();
				Metric2dContext mc(space);
				mesh->statMetricAlphaQuality(mc, alpha_stats, false);
				mesh->statMetricEdgeLength(mc, edge_stats, false);
				mesh->statRealAngles(angle_stats, surface);
			}
			if(alpha_stats.calculate()){
				LOG4CPLUS_INFO(MeshLog::logger_console, "Metric alpha quality (ave-geom)", alpha_stats.averageGeometric());
				LOG4CPLUS_INFO(MeshLog::logger_mesh, "[stat] Metric-alpha-ave-geom\t" << alpha_stats.averageGeometric());
			}
			if(edge_stats.calculate()){
				LOG4CPLUS_INFO(MeshLog::logger_console, "Metric edge lengths  (ave)", edge_stats.average());
				LOG4CPLUS_INFO(MeshLog::logger_console, "-------------------- (dev)", edge_stats.stdDev());
				LOG4CPLUS_INFO(MeshLog::logger_mesh, "[stat] Metric-edge-ave\t" << edge_stats.average());
				LOG4CPLUS_INFO(MeshLog::logger_mesh, "[stat] Metric-edge-stddev\t" << edge_stats.stdDev());
			}				
			if(angle_stats.calculate()){
				double angles_min = angle_stats.minimum() * 180.0 / PI;
				double angles_30 = angle_stats.getDataCountBottom(PI/6) / (double)angle_stats.countInt();
				LOG4CPLUS_INFO(MeshLog::logger_console, "Angles (min)", angles_min);
				LOG4CPLUS_INFO(MeshLog::logger_console, "------ (<30)", angles_30);
				LOG4CPLUS_INFO(MeshLog::logger_mesh, "[stat] Angles-min\t" << angles_min);
				LOG4CPLUS_INFO(MeshLog::logger_mesh, "[stat] Angles-below30\t" << angles_30);
			}				
		//---------------------------------------------------
		}else if(command.compare("quality3_detailed") == 0){
			if(!boundary){	LOG4CPLUS_ERROR(MeshLog::logger_console,   "Undefined boundary");	return CM_ERROR_NOMESH;	}
			for(IteratorMesh3d it = boundary->getFirstValidMesh3d(); it.isValid(); it.nextValidMesh()){
				MeshContainer3d* mesh = it.getMesh();
				// show block quality
				LOG4CPLUS_INFO(MeshLog::logger_console, "Quality: tetrahedra aspect ratio (~S/V)");
				const int RANGE1_COUNT = 5;
				double ranges1[RANGE1_COUNT] = {0.8, 0.6, 0.4, 0.2, 0.0};
				int stats1[RANGE1_COUNT+1];
				MeshData::StatData stat_data[5];
				int block_types[5];
				Metric3dContext mc(it.getMesh()->getControlSpace());
				int total_block_count = mesh->getTetrahedraQualityRange(mc, RANGE1_COUNT, ranges1, stats1, stat_data, block_types);
				for(i = 0; i <= RANGE1_COUNT; i++){
					ostringstream text;
					if(i < RANGE1_COUNT){
						text << "["; text.width(3); text << ranges1[i];
					}else text << "[-oo";
					text << ">] ";
					int ct = (int)((stats1[i] / (double)total_block_count) * 40.0);
					for(int j = 0; j <= ct; j++) text << '#';
					text << " (" << stats1[i] << ")";
					LOG4CPLUS_INFO(MeshLog::logger_console, text.str());
				}
				if(total_block_count > 0){
					ostringstream text;
					text << "Min = " << stat_data[4].minimum << ", Max = " << stat_data[4].maximum;
					text << ", Ave = " << stat_data[4].average;
					LOG4CPLUS_INFO(MeshLog::logger_console, text.str());
					const char* str_types[] = {"Inner:    ", "B-Vertex: ", "B-Edge:   ", "B-Face:   "};
					for(i = 0; i < 4; i++){
						ostringstream text_str;
						text_str << str_types[i] << block_types[i] << " Min = ";
						text_str << stat_data[i].minimum << ", Max = " << stat_data[i].maximum;
						text_str << ", Ave = " << stat_data[i].average;
						LOG4CPLUS_INFO(MeshLog::logger_console, text_str.str());
					}
				}
				LOG4CPLUS_INFO(MeshLog::logger_console, "----------------------------------------");
				// show edge metric quality
				LOG4CPLUS_INFO(MeshLog::logger_console, "Quality: metric edge length");
				const int RANGE2_COUNT = 5;
				double ranges2[RANGE2_COUNT] = {1.5, 1.3, 1.1, 0.9, 0.7};
				int stats2[RANGE2_COUNT+1];
				int total_edge_count = mesh->getMetricEdgeRange(mc, RANGE2_COUNT, ranges2, stats2);
				for(i = 0; i <= RANGE2_COUNT; i++){
					ostringstream text;
					if(i < RANGE2_COUNT){
						text << "["; text.width(3); text << ranges2[i];
					}else text << "[-oo";
					text << ">] ";
					int ct = (int)((stats2[i] / (double)total_edge_count) * 40.0);
					for(int j = 0; j <= ct; j++) text << '#';
					text << " (" << stats2[i] << ")";
					LOG4CPLUS_INFO(MeshLog::logger_console, text.str());
				}
			}
		//---------------------------------------------------
		}else if(command.compare("stats") == 0){
			if(!boundary){	LOG4CPLUS_ERROR(MeshLog::logger_console,   "Undefined boundary");	return CM_ERROR_NOMESH;	}
			for(IteratorMesh2d it = boundary->getFirstValidMesh2d(); it.isValid(); it.nextValidMesh()){
				MeshContainer2d* mesh = it.getMesh();
				LOG4CPLUS_INFO(MeshLog::logger_console, "------------------------------");
				LOG4CPLUS_INFO(MeshLog::logger_console, "Number of nodes", mesh->getPointsCount());
				LOG4CPLUS_INFO(MeshLog::logger_console, "Number of triangles", mesh->getElementsCount(3));
				LOG4CPLUS_INFO(MeshLog::logger_console, "Number of quads", mesh->getElementsCount(4));
				ControlSpace2dMesh* space = (ControlSpace2dMesh*)mesh->getControlSpace();
				if(space && (space->getType() == MeshData::CONTROL_MESH)){
					LOG4CPLUS_INFO(MeshLog::logger_console, "Control mesh nodes count", space->getControlNodesCount());
					LOG4CPLUS_INFO(MeshLog::logger_console, "Control mesh triangles count", space->getTrianglesCount());
//					LOG4CPLUS_INFO(MeshLog::logger_console, "Control mesh counter (simple)", space->getCounter(0));
//					LOG4CPLUS_INFO(MeshLog::logger_console, "Control mesh counter (triangle)", space->getCounter(1));
//					LOG4CPLUS_INFO(MeshLog::logger_console, "Control mesh counter (voronoi)", space->getCounter(2));
//					space->clearCounters();
				}
			}
		//---------------------------------------------------
		}else if(command.compare("quality_extended") == 0){
			if(!boundary){	LOG4CPLUS_ERROR(MeshLog::logger_console,   "Undefined boundary");	return CM_ERROR_NOMESH;	}
			ofstream file("mesh-quality.txt");
			for(IteratorMesh2d it = boundary->getFirstValidMesh2d(); it.isValid(); it.nextValidMesh()){
				MeshDomainSurface* domain_surface = it.getDomainSurface();
				file << "---- QUALITY - SELECTED MESH ----" << endl;
				file << "---------------------------------" << endl;
				domain_surface->assessQuality(file);
				file << "---------------------------------" << endl;
			}
		//---------------------------------------------------
		}else if(command.compare("eps") == 0){
			if(!boundary){	LOG4CPLUS_ERROR(MeshLog::logger_console,   "Undefined boundary");	return CM_ERROR_NOMESH;	}
			i = 0;
			for(IteratorMesh2d it = boundary->getFirstValidMesh2d(); it.isValid(); it.nextValidMesh()){
				it.getDomainSurface()->storeEPS(i);
				MeshContainer2d* mesh = it.getMesh();
				if(mesh){
					CS2dPtr space = (CS2dPtr)mesh->getControlSpace();
					if(space && ((space->getType() == MeshData::CONTROL_UNIFORM) ||
							(space->getType() == MeshData::CONTROL_MESH) ||
							(space->getType() == MeshData::CONTROL_QUADTREE)))
					{
						space->storeEPS("control", i);
					}
				}
				++i;
			}
		//---------------------------------------------------
		}else if(command.compare("set") == 0){
			string name;
			str_text >> name;
			i = mesh_data.findProperty(name);
			if(i >= 0){
				if(mesh_data.getPropertyType(i) == MeshData::PROP_INT){
					str_text >> i;
					if(str_text) mesh_data.setPropertyInt(name, i);
					else LOG4CPLUS_ERROR(MeshLog::logger_console,   "ERROR reading int value");
				}else{
					str_text >> d;
					if(str_text) mesh_data.setPropertyDouble(name, d);
					else LOG4CPLUS_ERROR(MeshLog::logger_console,   "ERROR reading double value");
				}
			}else{
				LOG4CPLUS_ERROR(MeshLog::logger_console,   "ERROR unknown property name");
			}
		//---------------------------------------------------
		}else if(command.compare("gen1.5d") == 0){
			if(!boundary){	LOG4CPLUS_ERROR(MeshLog::logger_console,   "Undefined boundary");	return CM_ERROR_NOMESH;	}
			MeshGenerator2d::prepareBoundaryMesh(boundary);
		//---------------------------------------------------
		}else if(command.compare("gen2d") == 0){
			if(!boundary){	LOG4CPLUS_ERROR(MeshLog::logger_console,   "Undefined boundary");	return CM_ERROR_NOMESH;	}
			MeshGenerator2d::triangulateFaces(boundary);
		//---------------------------------------------------
		}else if(command.compare("gen2.5d") == 0){
			if(!boundary){	LOG4CPLUS_ERROR(MeshLog::logger_console,   "Undefined boundary");	return CM_ERROR_NOMESH;	}
			MeshGenerator3d::prepareBoundaryMesh(boundary);
		//---------------------------------------------------
		}else if(command.compare("gen3d") == 0){
			if(!boundary){	LOG4CPLUS_ERROR(MeshLog::logger_console,   "Undefined boundary");	return CM_ERROR_NOMESH;	}
			MeshGenerator3d::createTetrahedra(boundary);
		//---------------------------------------------------
		}else if(command.compare("split") == 0){
			if(!boundary){	LOG4CPLUS_ERROR(MeshLog::logger_console,   "Undefined boundary");	return CM_ERROR_NOMESH;	}
			MeshGenerator2dQuad::convertFacesToQuads(boundary, MeshData::QUADS_SPLIT);
		//---------------------------------------------------
		}else if(command.compare("leelo") == 0){
			str_text >> i;	if(!str_text) i = -1;
			if(!boundary){	LOG4CPLUS_ERROR(MeshLog::logger_console,   "Undefined boundary");	return CM_ERROR_NOMESH;	}
			MeshGenerator2dQuad::convertFacesToQuads(boundary, MeshData::QUADS_LEELO, i);
		//---------------------------------------------------
		}else if(command.compare("qmorph") == 0){
			str_text >> i;	if(!str_text) i = -1;
			if(!boundary){	LOG4CPLUS_ERROR(MeshLog::logger_console,   "Undefined boundary");	return CM_ERROR_NOMESH;	}
			clock_t clock_start = clock();
			MeshGenerator2dQuad::convertFacesToQuads(boundary, MeshData::QUADS_QMORPH, i);
			double sec = (clock() - clock_start)/(double)CLOCKS_PER_SEC;
			int nct = 0, tct = 0, qct = 0;
			for(IteratorMesh2d it = boundary->getFirstValidMesh2d(); it.isValid(); it.nextValidMesh()){
				MeshContainer2d* mesh = it.getMesh();
				nct += mesh->getPointsCount();
				tct += mesh->getElementsCount(3);
				qct += mesh->getElementsCount(4);
			}
			LOG4CPLUS_INFO(MeshLog::logger_console, "nodes", nct);
			LOG4CPLUS_INFO(MeshLog::logger_console, "triangles", tct);
			LOG4CPLUS_INFO(MeshLog::logger_console, "quads", qct);
			LOG4CPLUS_INFO(MeshLog::logger_console, "Conversion speed (quads per sec.)", qct / sec);
		//---------------------------------------------------
		}else if(command.compare("make_quads") == 0){
			if(!boundary){	LOG4CPLUS_ERROR(MeshLog::logger_console,   "Undefined boundary");	return CM_ERROR_NOMESH;	}
			str_text >> arg1;
			clock_t clock_start = clock();			
			if(arg1.compare("leelo") == 0)
				MeshGenerator2dQuad::convertFacesToQuads(boundary, MeshData::QUADS_LEELO);
			else if(arg1.compare("mixed") == 0)
				MeshGenerator2dQuad::convertFacesToQuads(boundary, MeshData::QUADS_MIXED);
			else
				MeshGenerator2dQuad::convertFacesToQuads(boundary, MeshData::QUADS_QMORPH);

			MeshGenerator2dQuad::smoothenFaces(boundary, 3, MeshData::SM_TOP_SWAP | MeshData::SM_LAPLACE);
			double sec = (clock() - clock_start)/(double)CLOCKS_PER_SEC;
			int nct = 0, tct = 0, qct = 0;
			for(IteratorMesh2d it = boundary->getFirstValidMesh2d(); it.isValid(); it.nextValidMesh()){
				MeshContainer2d* mesh = it.getMesh();
				nct += mesh->getPointsCount();
				tct += mesh->getElementsCount(3);
				qct += mesh->getElementsCount(4);
			}
			LOG4CPLUS_INFO(MeshLog::logger_console, "nodes", nct);
			LOG4CPLUS_INFO(MeshLog::logger_console, "triangles", tct);
			LOG4CPLUS_INFO(MeshLog::logger_console, "quads", qct);
			LOG4CPLUS_INFO(MeshLog::logger_console, "Conversion speed (quads per sec.)", qct / sec);
			LOG4CPLUS_INFO(MeshLog::logger_mesh, "[stat]\tmake-quads-NT:\t" << tct);
			LOG4CPLUS_INFO(MeshLog::logger_mesh, "[stat]\tmake-quads-NQ:\t" << qct);
			LOG4CPLUS_INFO(MeshLog::logger_mesh, "[stat]\tmake-quads-NQs:\t" << ((sec > 0.0)?(qct / sec):99999));
		//---------------------------------------------------
		}else if(command.compare("make_quads_for_patch") == 0){
			//LOG4CPLUS_INFO(MeshLog::logger_console, "MAKE_QUADS_FOR_PATCH <file.resu.surf> [qmorph/leelo/mixed]");
			string fname;
			str_text >> fname;
			if(!str_text){ LOG4CPLUS_ERROR(MeshLog::logger_console, "Error reading file", fname);	return CM_ERROR_NOMESH;	}
			str_text >> arg1;
			MeshContainer2d* mesh =  // TODO read from file (+ planar surface)
				MeshStream::readFileResuSurf(fname);
			mesh->createControlSpaceFromTriangles(1.075); // create control space
			mesh->countMetricDifferenceQuality();
//			MeshView::setViewMode(SHOW_MESH_SURF_QUALITY);
//			MeshView::showViewSet("Surface mesh", mesh->getViewSet(), 2);
//			MeshView::setViewMode(SHOW_MESH_SURF);
			Metric2dContext mc(mesh->getControlSpace());
			clock_t clock_start = clock();			
			if(arg1.compare("leelo") == 0)
				MeshGenerator2dQuad::convertToQuadsMixed(mc, mesh, MeshData::QUADS_LEELO);
			else if(arg1.compare("mixed") == 0)
				MeshGenerator2dQuad::convertToQuadsMixed(mc, mesh, MeshData::QUADS_MIXED);
			else
				MeshGenerator2dQuad::convertToQuadsMixed(mc, mesh, MeshData::QUADS_QMORPH);

//			MeshView::showViewSet("Surface mesh", mesh->getViewSet(), 1);
			MeshGenerator2dQuad::improveQuads(mc, mesh, 5, MeshData::SM_TOP_SWAP | MeshData::SM_LAPLACE_MIXED);
			double sec = (clock() - clock_start)/(double)CLOCKS_PER_SEC;
			LOG4CPLUS_INFO(MeshLog::logger_console, "nodes", mesh->getPointsCount());
			int tct = mesh->getElementsCount(3);
			int qct = mesh->getElementsCount(4);
			LOG4CPLUS_INFO(MeshLog::logger_console, "triangles", tct);
			LOG4CPLUS_INFO(MeshLog::logger_console, "quads", qct);
			LOG4CPLUS_INFO(MeshLog::logger_console, "Conversion speed (quads per sec.)", qct / sec);
			LOG4CPLUS_INFO(MeshLog::logger_mesh, "[stat]\tmake-quads-NT:\t" << tct);
			LOG4CPLUS_INFO(MeshLog::logger_mesh, "[stat]\tmake-quads-NQ:\t" << qct);
			LOG4CPLUS_INFO(MeshLog::logger_mesh, "[stat]\tmake-quads-NQs:\t" << ((sec > 0.0)?(qct / sec):99999));
			if(true){
				DataStatistics angle_stats;
				DataStatistics edge_stats;
				DataStatistics alpha_stats;
				mesh->statMetricAlphaQuality(mc, alpha_stats, false);
				mesh->statMetricEdgeLength(mc, edge_stats, false);
				mesh->statRealAngles(angle_stats, mesh->getSurface());
				if(alpha_stats.calculate()){
					LOG4CPLUS_INFO(MeshLog::logger_console, "Metric alpha quality (ave-geom)", alpha_stats.averageGeometric());
					LOG4CPLUS_INFO(MeshLog::logger_mesh, "[stat] Metric-alpha-ave-geom\t" << alpha_stats.averageGeometric());
				}
				if(edge_stats.calculate()){
					LOG4CPLUS_INFO(MeshLog::logger_console, "Metric edge lengths  (ave)", edge_stats.average());
					LOG4CPLUS_INFO(MeshLog::logger_console, "-------------------- (dev)", edge_stats.stdDev());
					LOG4CPLUS_INFO(MeshLog::logger_mesh, "[stat] Metric-edge-ave\t" << edge_stats.average());
					LOG4CPLUS_INFO(MeshLog::logger_mesh, "[stat] Metric-edge-stddev\t" << edge_stats.stdDev());
				}				
				if(angle_stats.calculate()){
					double angles_min = angle_stats.minimum() * 180.0 / PI;
					double angles_30 = angle_stats.getDataCountBottom(PI/6) / (double)angle_stats.countInt();
					LOG4CPLUS_INFO(MeshLog::logger_console, "Angles (min)", angles_min);
					LOG4CPLUS_INFO(MeshLog::logger_console, "------ (<30)", angles_30);
					LOG4CPLUS_INFO(MeshLog::logger_mesh, "[stat] Angles-min\t" << angles_min);
					LOG4CPLUS_INFO(MeshLog::logger_mesh, "[stat] Angles-below30\t" << angles_30);
				}				
			}
			// TODO store quads to file
//			MeshView::showViewSet("Surface mesh", mesh->getViewSet(), 5);
			MeshStream::storeFileResuSurf(fname+"-quad.resu.surf", mesh);
			delete mesh;
		//---------------------------------------------------
		}else if(command.compare("get") == 0){
			str_text >> arg1;
			if(str_text){
				// show 
				i = mesh_data.findProperty(arg1);
				if(i < 0){
					LOG4CPLUS_ERROR(MeshLog::logger_console, "Property not found", arg1);
				}else{
					if(mesh_data.getPropertyType(i) == MeshData::PROP_INT){
						LOG4CPLUS_INFO(MeshLog::logger_console, mesh_data.getPropertyName(i), mesh_data.getPropertyInt(i));
					}else{
						LOG4CPLUS_INFO(MeshLog::logger_console, mesh_data.getPropertyName(i), mesh_data.getPropertyDouble(i));
					}
					LOG4CPLUS_INFO(MeshLog::logger_console, mesh_data.getPropertyDescription(i));
				}
			}else{
				// show all
				int ct = mesh_data.getPropertyCount();
				for(i = 0; i < ct; i++){
					ostringstream out_text;
					out_text << mesh_data.getPropertyName(i) << " = ";
					if(mesh_data.getPropertyType(i) == MeshData::PROP_INT){
						out_text << mesh_data.getPropertyInt(i);
					}else{
						out_text << mesh_data.getPropertyDouble(i);
					}
					LOG4CPLUS_INFO(MeshLog::logger_console, out_text.str());
				}
			}
*/
		//---------------------------------------------------
		}else{
			LOG4CPLUS_ERROR(MeshLog::logger_console, "ERROR parsing: " << cmd_line);
			return CM_ERROR_PARSE;
		}
#ifdef USE_EXCEPTIONS
	}catch(MeshingException& e){
		LOG4CPLUS_ERROR(MeshLog::logger_console, "Unexpected meshing error", e.getDescription());
		if(m_model_mesh) delete m_model_mesh;
		m_model_mesh = nullptr;
		LOG4CPLUS_INFO(MeshLog::logger_mesh, "[stat]\t* Unexpected meshing error:");
		LOG4CPLUS_INFO(MeshLog::logger_mesh, "[stat]\t* " << e.getDescription());
		LOG4CPLUS_INFO(MeshLog::logger_mesh, "[stat]\t***");
		LOG4CPLUS_INFO(MeshLog::logger_mesh, "[stat]\t***");
		return CM_EXCEPTION;
	}
#ifndef _DEBUG
	catch(...){
		LOG4CPLUS_ERROR(MeshLog::logger_console,   "Unexpected error");
		if(m_model_mesh) delete m_model_mesh;
		m_model_mesh = nullptr;
		LOG4CPLUS_INFO(MeshLog::logger_mesh, "[stat]\t* Unexpected error:");
		LOG4CPLUS_INFO(MeshLog::logger_mesh, "[stat]\t***");
		LOG4CPLUS_INFO(MeshLog::logger_mesh, "[stat]\t***");
		LOG4CPLUS_INFO(MeshLog::logger_mesh, "[stat]\t***");
		return CM_EXCEPTION;
	}
#endif // _DEBUG
#endif // USE_EXCEPTIONS

	return CM_OK;
}

int MeshModel::testCrack(int mode)
{
	if(!m_model_mesh){
		LOG4CPLUS_ERROR(MeshLog::logger_console,   "No meshing model available!");
		return CM_ERROR_NOMESH;
	}
	MeshContainer3d* mesh3d = m_model_mesh->getTotalMDV()->getMesh();
	if(!mesh3d){
		LOG4CPLUS_ERROR(MeshLog::logger_console,   "No tetrahedral mesh available!");
		return CM_ERROR_NOMESH;
	}
//	DBox crack_box(-1.0, -0.85, -0.0015, 0.0015, -1.0, 1.0);
//	DBox crack_box(-1.0, 0.5, -0.2, 0.2, -1.0, 1.0);
/*
	DBox crack_box(-1.0, -0.95, -0.04, 0.0, -0.2, 0.2);
	int destroyed_count = MeshSpecialRoutinesUTC::markCrack(mesh3d, crack_box);
	LOG4CPLUS_INFO(MeshLog::logger_console, "Marked (destroyed) blocks", destroyed_count);
	MeshContainer3d* cut_mesh = MeshSplit3d::splitMeshByBlocks(mesh3d, TagExtended::TAG_CUT_ADAPT, 1);
	m_model_mesh->addMeshBlock(new MeshDomainVolume(cut_mesh));
*/
	DataCompoundList<DTriangle3d> crack_planes(10);
	const double X0 = 0.0, TX = 0.05; // 0.05;
	const double TY = 0.5;
	const double Z0 = 0.0, Z1 = 1.0;
	crack_planes.append(DTriangle3d(DPoint3d(TX, TY, Z1), DPoint3d(TX, TY, Z0), DPoint3d(X0, TY, Z0)));
	crack_planes.append(DTriangle3d(DPoint3d(TX, TY, Z1), DPoint3d(X0, TY, Z0), DPoint3d(X0, TY, Z1)));
	const DVector3d crack_dv(0.045, 0.0, 0.0);
	ControlSpace2dAdaptive::param_stretch_max_ratio = 2.0;
	ControlSpace2dAdaptive::param_gradation_ratio = 2.0;
	CS3dPtr cs_next(MeshSpecialRoutinesUTC::prepareNextCS(mesh3d, crack_planes, crack_dv));
	mesh3d->setControlSpace(cs_next);
	Metric3dContext mc(cs_next);
	int destroyed_count = MeshSpecialRoutinesUTC::createCrack(mc, mesh3d, crack_planes);
	MeshSpecialRoutinesUTC::remeshCrack(mc, mesh3d, mode);

/*
	// -> mark faces  TAG_ADAPT_SURF
	MeshViewSet* set = new MeshViewSet;
	DBox mark_box(-1.0, -0.955, -0.01, 0.0, -0.19, 0.19);
	int cfct = 0;
	for(IteratorFace it = mesh3d->getFirstFace(); it.isValid(); it.nextFace()){
		MeshFace* face = it.getFace();
		if(face->isBoundedBothSides()) continue;
		if(mark_box.contains(face->getMiddlePoint())){
			face->setIntTag(TagExtended::TAG_ADAPT_SURF, 1);
			cfct++;
		}
		double max_dist2 = 0.0;
		for(int i = 0; i < face->getPointCount(); i++){
			double dist2 = abs(1.0 - DPoint3d::zero.distance2(face->getPoint(i)->getCoordinates()));
			if(dist2 > max_dist2) max_dist2 = dist2;
		}
		if(max_dist2 < 0.0001)
			face->setIntTag(TagExtended::TAG_ADAPT_SURF, 2);

		set->addFace(face, face->getIntTag(TagExtended::TAG_ADAPT_SURF, 0), MeshViewSet::param_shrink);
	}
	SHOW_MESH("Surface pre-marking", set);

	// ---> check quality at the beginning
	//set = new MeshViewSet;
	//int bct = mesh3d->getBlocksCount();
	//ControlSpace3dIdentity cs;
	//Metric3dContext mc(&cs);
	//for(int i = 0; i < bct; i++){
	//	MeshBlock* block = mesh3d->getBlockAt(i);
	//	set->addEdges(block);
	//	set->setPolygonFillMode(MeshViewSet::FILL_QUALITY);
	//	double q = block->countQuality(mc, MeshData::QUALITY3D_MEAN_RATIO);
	//	if(q < 0.2) set->addBlock(block);
	//}
	//SHOW_MESH("Mesh quality", set);

	DataMatrix<DPoint3d> crack_mpoints(cfct);
	DataCompoundList<MeshFace*> crack_faces;
	DataHashTable<MeshBlock*> hblocks(2*mesh3d->getBlocksCount(), nullptr);
	DataHashTable<MeshPoint3d*> hpoints(2*mesh3d->getPointsCount(), nullptr);
	for(IteratorFace it = mesh3d->getFirstFace(); it.isValid(); it.nextFace()){
		MeshFace* face = it.getFace();
		if(face->getIntTag(TagExtended::TAG_ADAPT_SURF, 0) == 1){
			crack_faces.append(face);
			crack_mpoints.add(face->getMiddlePoint());
			int fpct = face->getPointCount();
			for(int i = 0; i < fpct; i++){
				MeshPoint3d* point = face->getPoint(i);
				if(hpoints.contains(point)) 
					continue;
				// check, if within the surface patch
				DataVector<MeshFace*> pfaces(100);
				point->adjacentFaces(pfaces);
				bool ok = true;
				for(int j = 0; ok && (j < pfaces.countInt()); j++){
					MeshFace* face = pfaces[j];
					if(!face->isBoundedBothSides() &&
						face->getIntTag(TagExtended::TAG_ADAPT_SURF, 0) != 1)
					{
						ok = false; break;
					}
				}
				// insert
				if(ok && hpoints.insert(point)){
					DataVector<MeshBlock*> pblocks(100);
					point->adjacentBlocks(pblocks);
					for(int j = 0; j < pblocks.countInt(); j++)
						hblocks.insert(pblocks[j]);
				}
			}
		}
	}
	DataVector<MeshBlock*> crack_blocks(hblocks.countInt());
	hblocks.getValues(crack_blocks);
	DataVector<MeshPoint3d*> crack_points(hpoints.countInt());
	hpoints.getValues(crack_points);

	DPoint3d plane_pt;
	DVector3d plane_e0, plane_e1;
	double plane_max_dist = DLeastSquaresFitting::fitHyperplaneOrthogonal(crack_mpoints, 
						plane_pt, plane_e0, plane_e1);
	LOG4CPLUS_INFO(MeshLog::logger_mesh, "max point-to-plane distance = " <<  plane_max_dist);

	SurfacePlane crack_plane(plane_pt, plane_e0, plane_e1);

	DRect brect;
	for(int i = 0; i < crack_mpoints.countInt(); i++)
		brect.addPoint(crack_plane.getParameters(crack_mpoints[i]));

	if(true){
		MeshViewSet* set = new MeshViewSet;
		for (auto it = crack_faces.iterator(); it.valid(); it.moveNext()) {
		MeshFace* face = it.item();
			set->addFaceWithEdges(face, 1, MeshViewSet::param_shrink);
		}
		set->addEdge(crack_plane.getPoint(brect.getX0Y0()), crack_plane.getPoint(brect.getX1Y0()), 2);
		set->addEdge(crack_plane.getPoint(brect.getX0Y1()), crack_plane.getPoint(brect.getX1Y1()), 2);
		set->addEdge(crack_plane.getPoint(brect.getX0Y0()), crack_plane.getPoint(brect.getX0Y1()), 2);
		set->addEdge(crack_plane.getPoint(brect.getX1Y0()), crack_plane.getPoint(brect.getX1Y1()), 2);

		ControlSpace3dIdentity cs;
		Metric3dContext mc(&cs);
		for(int i = 0; i < crack_blocks.countInt(); i++){
			MeshBlock* block = crack_blocks[i];
			double q = block->countQuality(mc, MeshData::QUALITY3D_MEAN_RATIO);
			if(q < 0.25){
				set->addBlock(block);
				LOG4CPLUS_INFO(MeshLog::logger_mesh, "testCrack:: block #" << i << ", q=" << q);
			}
		}

		SHOW_MESH("Approximated plane", set);
	}

	// -> smoothen
	if(mode == 2){ // move up and move down
		LOG4CPLUS_INFO(MeshLog::logger_mesh, "=================== move up and down ===================");
		// move points
		for(int i = 0; i < crack_points.countInt(); i++){
			MeshPoint3d* point = crack_points[i];
			DPoint3d plane_coord = crack_plane.getPoint(crack_plane.getParameters(point->getCoordinates()));
			//MeshGenerator3dQuality::tryMovingPoint(point, plane_coord, false);
			point->setCoordinates(plane_coord);
		}
		// show
		MeshViewSet* set = new MeshViewSet;
		for (auto it = crack_faces.iterator(); it.valid(); it.moveNext()) {
			MeshFace* face = it.item();
			set->addFaceWithEdges(face, 1, MeshViewSet::param_shrink);
		}
		set->addEdge(crack_plane.getPoint(brect.getX0Y0()), crack_plane.getPoint(brect.getX1Y0()), 2);
		set->addEdge(crack_plane.getPoint(brect.getX0Y1()), crack_plane.getPoint(brect.getX1Y1()), 2);
		set->addEdge(crack_plane.getPoint(brect.getX0Y0()), crack_plane.getPoint(brect.getX0Y1()), 2);
		set->addEdge(crack_plane.getPoint(brect.getX1Y0()), crack_plane.getPoint(brect.getX1Y1()), 2);

		ControlSpace3dIdentity cs;
		Metric3dContext mc(&cs);
		for(int i = 0; i < crack_blocks.countInt(); i++){
			MeshBlock* block = crack_blocks[i];
			double q = block->countQuality(mc, MeshData::QUALITY3D_MEAN_RATIO);
			if(q < 0.25){
				set->addBlock(block);
				LOG4CPLUS_INFO(MeshLog::logger_mesh, "testCrack:: block #" << i << ", q=" << q);
			}
		}
		SHOW_MESH("After move down and up", set);

		// ---> clear slivers at the surface
		for(int i = 0; i < crack_blocks.countInt(); i++){
			MeshBlock* block = crack_blocks[i];
			double q = block->getQuality();
			if(q > 0.1) continue;
			double crack_fct = 0;
			for(int i = 0; i < block->getFaceCount(); i++){
				// ...
			}
		}
	}
*/
	return CM_OK;
}

int MeshModel::testAdaptXML(const std::string& fname1, const std::string& fname2,
							  int adapt_type, int cmp_type, bool debug_mode)
{
	MeshContainer3d * model1;
	MeshContainer3d * model2;

	// Read two mesh geometries
	MeshBRep* desc = new MeshBRepXML;
	if(desc->parseFile(fname1) && desc->validate()){
		model1 = desc->createDomainMesh();
		delete desc;
		if(!model1) return CM_ERROR_PARSE;
	}
	desc = new MeshBRepXML;
	if(desc->parseFile(fname2) && desc->validate()){
		model2 = desc->createDomainMesh();
		delete desc;
		if(!model2) return CM_ERROR_PARSE;
	}

	assert(model1->getBlocksCount() >= 1);
	assert(model2->getBlocksCount() >= 1);

	MeshDomainVolume* dvolumes[2] = {
		(MeshDomainVolume*)model1->getBlockAt(0),
		(MeshDomainVolume*)model2->getBlockAt(0)
	};

	CS3dPtr ucs[2] = {
		dvolumes[0]->getUserControlSpace(),
		dvolumes[1]->getUserControlSpace()
	};

	if(!ucs[0] || !ucs[0]->isAdaptive() || !ucs[1] || !ucs[1]->isAdaptive()){
		LOG4CPLUS_ERROR(MeshLog::logger_console,   "Both UCS should be adaptive!");
		return CM_ERROR_RUN;
	}

	MeshGenerator2d::show_prediction = false;
	START_CLOCK("MESH-ADAPT:initial-gen");
	if(!MeshGenerator3d::autoTriangulate(model1, 2)){
		LOG4CPLUS_ERROR(MeshLog::logger_console,   "Error triangulating the first model");
		return CM_ERROR_RUN;
	}
	STOP_CLOCK("MESH-ADAPT:initial-gen");

	//if(!MeshGenerator3d::autoTriangulate(model2, 2)){
	//	LOG4CPLUS_ERROR(MeshLog::logger_console,   "Error triangulating the second model");
	//	return CM_ERROR_RUN;
	//}

	MeshContainer3d* mesh[2] = {
		dvolumes[0]->getMesh(),
		dvolumes[1]->getMesh()
	};

	ControlSpace3dAdaptive* cs[2] = {
		mesh[0]->getControlSpace()->getAsAdaptive(),
		ucs[1]->getAsAdaptive()
	};

	ControlSpace3dAdaptive* csh[2] = { cs[0], cs[1] };

	int pct = mesh[0]->getPointsCount();
	int tct = mesh[0]->getBlocksCount();

	LOG4CPLUS_INFO(MeshLog::logger_mesh, "MESH-ADAPT: mesh-0, NT = " << tct << ", PT = " << pct);

	MeshViewSet* set = nullptr;
	//// calculate CS-diff for own CS
	//if(debug_mode) set = new MeshViewSet(pct, 20+3*tct, 2*tct, tct);

	//Metric3dContext mc(cs[0]);
	//for(int i = 0; i < tct; i++){
	//	MeshTetrahedron* t = (MeshTetrahedron*)mesh[0]->getBlockAt(i);
	//	if(t->getType() != BLOCK_TETRA) continue;
	//	t->countMetricDiffQuality(mc);
	//	if(debug_mode)
	//		set->addBlockWithEdges(t);
	//}
	//if(debug_mode){
	//	set->setPolygonFillMode(MeshViewSet::FILL_QUALITY);
	//	SHOW_MESH("Metric quality q(T,CS) for own CS", set);
	//}

	enum { ADAPT_CUT = 0, ADAPT_NEW = 1, ADAPT_TRANSFORM_ALL = 2, ADAPT_TRANSFORM_CUT = 3};
	enum { CMP_T_CS = 0, CMP_CS_CS = 1};

//	adapt_typ = ADAPT_CUT;
//	cmp_type = CMP_T_CS;

	if(adapt_type == ADAPT_CUT || adapt_type == ADAPT_TRANSFORM_CUT){
		if(cmp_type == CMP_T_CS){
			LOG4CPLUS_INFO(MeshLog::logger_console, "Testing T-CS evaluation and adaptation");
			// calculate CS-diff using mesh - CS
			//if(debug_mode) 
			//	set = new MeshViewSet(pct, 20+2*tct, 2*tct, tct);

			START_CLOCK("CUT-ADAPT:T-CS");
			Metric3dContext mc(cs[1]);
			for(int i = 0; i < tct; i++){
				MeshTetrahedron* tetra = (MeshTetrahedron*)mesh[0]->getBlockAt(i);
				if(tetra->getType() != BLOCK_TETRA) continue;
				double q = tetra->countMetricDiffQuality(mc);
				for(int j = 0; j < 4; j++){
					MeshPoint3d* point = tetra->getPoint(j);
					double qp = point->getDoubleTag(TagExtended::TAG_VISUALIZATION, 1.0);
					if(q < qp)
						point->setDoubleTag(TagExtended::TAG_VISUALIZATION, q);
				}
				//if(debug_mode)
				//	set->addBlockWithEdges(tetra);
			}
			STOP_CLOCK("CUT-ADAPT:T-CS");
			//if(debug_mode){
			//	set->setPolygonFillMode(MeshViewSet::FILL_QUALITY);
			//	SHOW_MESH("Metric quality q(T,CS) for different CS", set);
			//}
		}else{ // cut using q(CS,CS)
			LOG4CPLUS_INFO(MeshLog::logger_console, "Testing CS-CS evaluation and adaptation");

			// calculate CS-CS diff using only nodes of one CS
			START_CLOCK("CUT-ADAPT:CS-CS");
			// calculate diff-quality in CS nodes
			for(int k = 0; k < 2; k++){
				cs[k]->forEachControlNode([&](ControlNode3d& cn) {
					double diff = cn.control_data.countDifferenceRR(
						cs[1 - k]->getMetricAtPoint(cn.coord));
					double q = (1.0 - 0.25 * diff);
					cn.setDoubleTag(TagExtended::TAG_METRIC_DIFF_RATIO, q);
				});
			}
			STOP_CLOCK("CUT-ADAPT:CS-CS");
			//if(debug_mode){
			//	MeshViewSet* set_cs1 = cs[0]->getViewSet(nullptr, TagExtended::TAG_METRIC_DIFF_RATIO);
			//	if(set_cs1){
			//		set_cs1->setPolygonFillMode(MeshViewSet::FILL_NODES);
			//		SHOW_MESH("Metric quality q(CS,CS)", set_cs1);
			//	}
			//}

			START_CLOCK("CUT-ADAPT:CS-CS-mark");
			for(int i = 0; i < pct; i++){
				MeshPoint3d* mesh_point = mesh[0]->getPointAt(i);
				const DPoint3d& pt = mesh_point->getCoordinates();
				double res0 = cs[0]->getLocalResolution(pt);
				double res1 = cs[1]->getLocalResolution(pt);
				double v = cs[(res0 < res1) ? 0 : 1]->interpolateDoubleTag(pt,
					TagExtended::TAG_METRIC_DIFF_RATIO);
				mesh_point->setDoubleTag(TagExtended::TAG_VISUALIZATION, v);
			}
			STOP_CLOCK("CUT-ADAPT:CS-CS-mark");

			//if(debug_mode){
			//	MeshViewSet* set1 = new MeshViewSet(pct, 20+2*tct, 2*tct, tct);
			//	for(int i = 0; i < tct; i++){
			//		MeshTetrahedron* tetra = (MeshTetrahedron*)mesh[0]->getBlockAt(i);
			//		double v = 1.0;
			//		for(int j = 0; j < tetra->getPointCount(); j++){
			//			const DPoint3d& pt = tetra->getPoint(j)->getCoordinates();
			//			double res0 = cs[0]->getLocalResolution(pt);
			//			double res1 = cs[1]->getLocalResolution(pt);
			//			double vlocal = cs[(res0 < res1) ? 0 : 1]->interpolateDoubleTag(pt,
			//				TagExtended::TAG_METRIC_DIFF_RATIO);
			//			if(vlocal < v) v = vlocal;
			//		}
			//		if(tetra->getType() != BLOCK_TETRA) continue;
			//		tetra->setQuality(v);
			//		set1->addBlockWithEdges(tetra);
			//	}
			//	set1->setPolygonFillMode(MeshViewSet::FILL_QUALITY);
			//	SHOW_MESH("Metric quality q(CS,CS)", set1);
			//}
		}

		START_CLOCK("CUT-ADAPT:mark-layers");
		MeshContainer3d* cut_mesh;
		int cut_tct = 0;
		DataVector<int> layers(tct, 0);
		DataVector<int> p_layers(pct, 0);
		// mark elements incident to at least one marked node
		for(int i = 0; i < tct; i++){
			MeshBlock* block = mesh[0]->getBlockAt(i);
			block->removeTag(TagExtended::TAG_CUT_ADAPT);
			int bpct = block->getPointCount();
			for(int j = 0; j < bpct; j++){
				if(block->getPoint(j)->getDoubleTag(TagExtended::TAG_VISUALIZATION, 0.0) < 0.5){
					block->setIntTag(TagExtended::TAG_CUT_ADAPT, 1);
					layers[i] = 1;
					for(int k = 0; k < bpct; k++)
						p_layers[block->getPoint(k)->getIndex()] = 1;
					cut_tct++;
					break;
				}
			}
		}
		// one additional layer of elements ...
		for(int i = 0; i < tct; i++){
			if(layers[i] != 0) continue;
			MeshBlock* block = mesh[0]->getBlockAt(i);
			int bpct = block->getPointCount();
			for(int j = 0; j < bpct; j++){
				if(p_layers[block->getPoint(j)->getIndex()] > 0){
					block->setIntTag(TagExtended::TAG_CUT_ADAPT, 1);
					layers[i] = 2;
//					for(int k = 0; k < bpct; k++)
//						p_layers[block->getPoint(k)->getIndex()] = 2;
					cut_tct++;
					break;
				}
			}
		}
		STOP_CLOCK("CUT-ADAPT:mark-layers");

		if(debug_mode){
			set = new MeshViewSet(pct, 3*tct, 3*tct, tct);
			for(int i = 0; i < tct; i++)
				if(layers[i] == 1) set->addBlockWithEdges(mesh[0]->getBlockAt(i));
			SHOW_MESH("First layer", set);

			set = new MeshViewSet(pct, 3*tct, 3*tct, tct);
			for(int i = 0; i < tct; i++)
				if(layers[i] > 0) set->addBlockWithEdges(mesh[0]->getBlockAt(i), layers[i]);
			SHOW_MESH("All layers", set);
		}

		LOG4CPLUS_INFO(MeshLog::logger_mesh, "CUT-ADAPT: marked-0, NT = " << cut_tct);

		START_CLOCK("CUT-ADAPT:cut-elements");
		cut_mesh = MeshSplit3d::splitMeshByBlocks(mesh[0], TagExtended::TAG_CUT_ADAPT, 1);
		STOP_CLOCK("CUT-ADAPT:cut-elements");

		LOG4CPLUS_INFO(MeshLog::logger_mesh, "CUT MESH: NT=" << cut_mesh->getBlocksCount() << " NP=" << cut_mesh->getPointsCount());

		if(debug_mode){
			//set = mesh[0]->getViewSet();
			//SHOW_MESH("Mesh left after cut", set);

			set = cut_mesh->getViewSet();
			SHOW_MESH("Cut mesh", set);
		}

		if(debug_mode){
			set = new MeshViewSet(pct, 20+2*tct, 2*tct, tct);
			Metric3dContext mc(cs[1]);
			for(int i = 0; i < cut_tct; i++){
				MeshTetrahedron* tetra = (MeshTetrahedron*)cut_mesh->getBlockAt(i);
				if(tetra->getType() != BLOCK_TETRA) continue;
				double q = tetra->countMetricDiffQuality(mc);
				set->addBlockWithEdges(tetra);
			}
			set->setPolygonFillMode(MeshViewSet::FILL_QUALITY);
			SHOW_MESH("Metric quality q(T,CS) for cut-mesh", set);
		}

		if(adapt_type == ADAPT_CUT){
			START_CLOCK("CUT-ADAPT:cut-retriangulate");

			MeshDomainVolume* cut_domain_volume = new MeshDomainVolume(cut_mesh);
			cut_domain_volume->setAreaID(cut_mesh->getBlockAt(0)->getAreaID());

			cut_domain_volume->setUserControlSpace(ucs[1]);
			cut_domain_volume->createInitialControlSpace();

			MeshContainer3dSurface* cut_surface_mesh = cut_domain_volume->cutSurfaceMeshFromVolumeMesh();

			LOG4CPLUS_INFO(MeshLog::logger_mesh, 
				"CUT FACES: NF=" << cut_surface_mesh->getFacesCount() 
				<< " NP=" << cut_surface_mesh->getPointsCount());


			Metric3dContext cut_mc(ucs[1]);
			cut_domain_volume->discretizeUsingBoundary(cut_mc);
			MeshContainer3d *new_cut_mesh = cut_domain_volume->getMesh();
			MeshGenerator3d::addInnerNodes(cut_mc, new_cut_mesh);
			MeshGenerator3dQuality::smoothen(cut_mc, new_cut_mesh, 2);

			STOP_CLOCK("CUT-ADAPT:cut-retriangulate");

			LOG4CPLUS_INFO(MeshLog::logger_mesh, 
				"CUT-ADAPT: new-mesh-cut-0, NT = " << new_cut_mesh->getBlocksCount() 
				<< ", PT = " << new_cut_mesh->getPointsCount());

			if(debug_mode){
				set = new_cut_mesh->getViewSet();
				SHOW_MESH("Cut mesh after retriangulation", set);
			}

			START_CLOCK("CUT-ADAPT:cut-merge");
			MeshSplit3d::mergeMeshes(mesh[0], new_cut_mesh);
			STOP_CLOCK("CUT-ADAPT:cut-merge");

		}else if(adapt_type == ADAPT_TRANSFORM_CUT){
			START_CLOCK("TRANSFORM-CUT-ADAPT:transform");
			Metric3dContext mc(ucs[1]);
			MeshGenerator3d::addInnerNodes(mc, cut_mesh);
			MeshGenerator3dQuality::smoothen(mc, cut_mesh, 2);
			STOP_CLOCK("TRANSFORM-CUT-ADAPT:transform");

			LOG4CPLUS_INFO(MeshLog::logger_mesh, 
				"TRANSFORM-ADAPT: new-cut-mesh, NT=" << cut_mesh->getBlocksCount() 
				<< ", PT=" << cut_mesh->getPointsCount());

			if(debug_mode){
				set = cut_mesh->getViewSet();
				SHOW_MESH("Cut mesh after transformation", set);
			}

			START_CLOCK("TRANSFORM-CUT-ADAPT:cut-merge");
			MeshSplit3d::mergeMeshes(mesh[0], cut_mesh);
			STOP_CLOCK("TRANSFORM-CUT-ADAPT:cut-merge");
		}
	}else if(adapt_type == ADAPT_TRANSFORM_ALL){

		START_CLOCK("TRANSFORM-ADAPT:transform");
		Metric3dContext mc(ucs[1]);
//		MeshGenerator3dQuality::smoothen(mc, mesh[0], 2);
		MeshGenerator3d::addInnerNodes(mc, mesh[0]);
		MeshGenerator3dQuality::smoothen(mc, mesh[0], 2);
		STOP_CLOCK("TRANSFORM-ADAPT:transform");

		LOG4CPLUS_INFO(MeshLog::logger_mesh, 
			"TRANSFORM-ADAPT: new-mesh, NT=" << mesh[0]->getBlocksCount() 
			<< ", PT=" << mesh[0]->getPointsCount());

	}

	if(mesh[0]){
		DataStatistics mean_ratio_stats;
		DataStatistics edge_stats;
		Metric3dContext mc(ucs[1]);
		mesh[0]->statMeanRatio(mc, mean_ratio_stats, false);
		mesh[0]->statMetricEdgeLength(mc, edge_stats, false);

		static auto logger_cstat = Logger::getInstance(LOG4CPLUS_TEXT("mesh.stat.cstat"));
		if(mean_ratio_stats.calculate()){
			LOG4CPLUS_INFO(logger_cstat, "MetricMeanRatio-ave\t" << mean_ratio_stats.average());
			LOG4CPLUS_INFO(logger_cstat, "MetricMeanRatio-min\t" << mean_ratio_stats.minimum());
		}
		if(edge_stats.calculate()){
			LOG4CPLUS_INFO(logger_cstat, "Metric edge lengths (ave)" << edge_stats.average());
			LOG4CPLUS_INFO(logger_cstat, "MetricEdgeLength-ave\t" << edge_stats.average());
		}				
	}
/*

		}else{ // mode 3 or 4 -> in-place mesh modification
			for(int k = 0; k < 2; k++){
				// mark points with all elements around tagged
				for(int i = 0; i < pcts[k]; i++){
					MeshPoint2d* point = mesh[k]->getPointAt(i);
					int rank = point->getRank();
					bool within_cut = true;
					int last_tag_value = -1;
					for(int j = 0; j < rank; j++){
						MeshElement* element = point->getEdge(j)->getMeshElement(point);
						if(element){
							if(!element->availableTag(TagExtended::TAG_CUT_ADAPT)){
								within_cut = false;
								break;
							}else{
								last_tag_value = element->getIntTag(TagExtended::TAG_CUT_ADAPT);
							}
						}
					}
					if(within_cut)
						point->setIntTag(TagExtended::TAG_CUT_ADAPT, last_tag_value);
				}
			}

			Metric2dContext mc0(cs[0]);
			Metric2dContext mc1(cs[1]);

			// edge-collapse
			START_CLOCK("CUT-ADAPT:trans-collapse");
			MeshGenerator2d::collapseEdges(mc1, mesh[0], TagExtended::TAG_CUT_ADAPT, 1);
			MeshGenerator2d::collapseEdges(mc0, mesh[1], TagExtended::TAG_CUT_ADAPT, 1);
			STOP_CLOCK("CUT-ADAPT:trans-collapse");

			// smoothen
			START_CLOCK("CUT-ADAPT:trans-smoothen");
			MeshGenerator2d::smoothen(mc1, mesh[0], 2, TagExtended::TAG_CUT_ADAPT, 1);
			MeshGenerator2d::smoothen(mc0, mesh[1], 2, TagExtended::TAG_CUT_ADAPT, 1);
			STOP_CLOCK("CUT-ADAPT:trans-smoothen");

			// edge-collapse
			START_CLOCK("CUT-ADAPT:trans-collapse");
			MeshGenerator2d::collapseEdges(mc1, mesh[0], TagExtended::TAG_CUT_ADAPT, 1);
			MeshGenerator2d::collapseEdges(mc0, mesh[1], TagExtended::TAG_CUT_ADAPT, 1);
			STOP_CLOCK("CUT-ADAPT:trans-collapse");

			if(debug_mode){
				set = mesh[0]->getViewSet();
				set = mesh[1]->getViewSet(set);
				SHOW_MESH("Trans mesh after cellapse + smoothen", set);
			}

			// refine (inner nodes) -> set inf.high quality for un-marked elements ?
			START_CLOCK("CUT-ADAPT:trans-refine");
			MeshGenerator2d::addInnerNodes(mc1, mesh[0], TagExtended::TAG_CUT_ADAPT, 1);
			MeshGenerator2d::addInnerNodes(mc0, mesh[1], TagExtended::TAG_CUT_ADAPT, 1);
			STOP_CLOCK("CUT-ADAPT:trans-refine");

			// smoothen
			START_CLOCK("CUT-ADAPT:trans-smoothen");
			MeshGenerator2d::smoothen(mc1, mesh[0], 2, TagExtended::TAG_CUT_ADAPT, 1);
			MeshGenerator2d::smoothen(mc0, mesh[1], 2, TagExtended::TAG_CUT_ADAPT, 1);
			STOP_CLOCK("CUT-ADAPT:trans-smoothen");

			if(debug_mode){
				set = mesh[0]->getViewSet();
				set = mesh[1]->getViewSet(set);
				SHOW_MESH("Trans mesh after refine + smoothen", set);
			}
		}

		if(debug_mode){
			set = mesh[0]->getViewSet();
			set = mesh[1]->getViewSet(set);
			SHOW_MESH("Whole mesh after transformation", set);
		}
	}else if(mode == 2){ // whole-mesh retriangulation
		for(int k = 0; k < 2; k++){
			domain_surfaces[k]->setUserControlSpace(ucs[1-k]);
			domain_surfaces[k]->clearDiscretization();
			domain_surfaces[k]->getBoundary()->setControlSpace(cs[1-k]);
			domain_surfaces[k]->createBoundaryMesh();
		}
		START_CLOCK("CUT-ADAPT:whole-retriangulate");
		for(int k = 0; k < 2; k++){
			domain_surfaces[k]->triangulate();
			domain_surfaces[k]->smoothen(2);
		}
		STOP_CLOCK("CUT-ADAPT:whole-retriangulate");

		mesh[0] = domain_surfaces[0]->getMesh();
		mesh[1] = domain_surfaces[1]->getMesh();
	}

	LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, "CUT-ADAPT: new-mesh-0, NT = " << mesh[0]->getElementsCount() 
		<< ", PT = " << mesh[0]->getPointsCount() << endl;
	LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, "CUT-ADAPT: new-mesh-1 NT = " << mesh[1]->getElementsCount() 
		<< ", PT = " << mesh[1]->getPointsCount() << endl;

	DataStatistics stat_edges;
	DataStatistics stat_metric;

	Metric2dContext mc0(cs[0]);
	Metric2dContext mc1(cs[1]);

	mesh[0]->statMetricEdgeLength(mc1, stat_edges, false);
	mesh[1]->statMetricEdgeLength(mc0, stat_edges, false);
	mesh[0]->statMetricDifference(mc1, stat_metric);
	mesh[1]->statMetricDifference(mc0, stat_metric);

	stat_edges.calculate();
	stat_metric.calculate();

	DataVector<double> edge_ranges;
	edge_ranges.add(0.5);
	edge_ranges.add(0.7);
	edge_ranges.add(0.9);
	edge_ranges.add(1.1);
	edge_ranges.add(1.3);
	edge_ranges.add(1.5);
	stat_edges.logStats("Edge metric length", "ADAPT-EDGE-LEN", edge_ranges);

	DataVector<double> metric_ranges;
	metric_ranges.add(0.5);
	metric_ranges.add(1.0);
	metric_ranges.add(2.0);
	metric_ranges.add(4.0);
	metric_ranges.add(8.0);
	stat_metric.logStats("Metric diff quality", "ADAPT-METRIC-DIFF", metric_ranges);

	if(debug_mode){
		set = nullptr;
		for(int k = 0; k < 2; k++){
			Metric2dContext mc(cs[1-k]);
			tcts[k] = mesh[k]->getElementsCount();
			for(int i = 0; i < tcts[k]; i++){
				MeshTriangle2d* triangle = (MeshTriangle2d*)mesh[k]->getElementAt(i);
				if(triangle->getType() != ELEMENT_MESH_TRIANGLE) continue;
				triangle->countMetricDiffQuality(mc);
				//set->addElementWithEdges(triangle, mesh[k]->getSurface());
			}
			set = mesh[k]->getViewSet(set);
		}
		set->setPolygonFillMode(MeshViewSet::FILL_QUALITY);
		SHOW_MESH("Whole mesh after whole-retriangulation", set);
	}

*/

	delete model1;
	delete model2;

	return CM_OK;
}

int MeshModel::testAdaptBin(const std::string& fname1, const std::string& fname2,
							  int adapt_type, int cmp_type, bool debug_mode)
{
//	MeshContainer3d * model1;
//	MeshContainer3d * model2;

	// load mesh 3D x2
	MeshContainer3d* vmesh1 = MeshSpecialRoutinesUTC::readBinFile(fname1);
	if(!vmesh1){
		LOG4CPLUS_ERROR(MeshLog::logger_console, "Error reading file: " << fname1);
		return CM_ERROR_PARSE;
	}
	MeshContainer3d* vmesh2 = MeshSpecialRoutinesUTC::readBinFile(fname2);
	if(!vmesh2){
		LOG4CPLUS_ERROR(MeshLog::logger_console, "Error reading file: " << fname2);
		return CM_ERROR_PARSE;
	}

	assert(vmesh1->getBlocksCount() >= 1);
	assert(vmesh2->getBlocksCount() >= 1);

	CS3dPtr  cs2 = MeshGenerator3dAdapt::createACSfromMeshBlocks(vmesh2);
	vmesh2->setControlSpace(cs2);

//	enum { ADAPT_CUT = 0, ADAPT_NEW = 1, ADAPT_TRANSFORM_ALL = 2, ADAPT_TRANSFORM_CUT = 3};
//	enum { CMP_T_CS = 0, CMP_CS_CS = 1};

	// for now - only ADAPT_TRANSFORM_ALL for testing, no compare and no cut

	if(debug_mode){
		SHOW_MESH("Mesh before", vmesh1->getViewSet());
	}

	Metric3dContext mc(cs2);
	MeshGenerator3dAdapt::remeshWithLocalTransformations(mc, vmesh1);

	STOP_CLOCK("TRANSFORM-ADAPT:transform");

	if(debug_mode){
		SHOW_MESH_NORESET("Volume mesh final final", vmesh1->getViewSet());
	}

	LOG4CPLUS_INFO(MeshLog::logger_mesh, 
		"TRANSFORM-ADAPT: new-mesh, NT=" << vmesh1->getBlocksCount() 
			<< ", PT=" << vmesh1->getPointsCount());

	delete vmesh1;
	delete vmesh2;

	return CM_OK;
}

int MeshModel::test(int mode, bool debug_mode)
{
	if(!m_model_mesh){	LOG4CPLUS_ERROR(MeshLog::logger_console,   "Undefined model");	return CM_ERROR_NOMESH;	}
	const MeshDomainVolume* mdv = m_model_mesh->getTotalMDV();
	MeshContainer3d* mesh3d = mdv->getMesh();
	if(!mesh3d){ LOG4CPLUS_ERROR(MeshLog::logger_console,   "No tetrahedral mesh available");	return CM_ERROR_NOMESH;	}

	MeshGenerator3dAdapt::identifyLocalSurfaces(mesh3d);

	return CM_OK;
}

/*
int MeshingCommands::convertSurfPatchtoQuads(const string& fname)
{
	MeshContainer2d* mesh = MeshStream::readFileResuSurf(fname);
	if(!mesh){
		LOG4CPLUS_ERROR(MeshLog::logger_console,   "Error reading file.");
		return -1;
	}
	mesh->createControlSpaceFromTriangles(1.075); // create control space
	mesh->countMetricDifferenceQuality();
//	MeshView::setViewMode(SHOW_MESH_SURF_QUALITY);
//	MeshView::showViewSet("Surface mesh", mesh->getViewSet(), 2);
//	MeshView::setViewMode(SHOW_MESH_SURF);
	Metric2dContext mc(mesh->getControlSpace());
//	clock_t clock_start = clock();			
	try{
		int qct = MeshGenerator2dQuad::convertToQuadsMixed(mc, mesh, MeshData::QUADS_QMORPH);
		if(qct <= 0){
			LOG4CPLUS_ERROR(MeshLog::logger_console,   "Error converting to quads.");
			return -1;
		}

//		MeshView::showViewSet("Surface mesh", mesh->getViewSet(), 1);
		MeshGenerator2dQuad::improveQuads(mc, mesh, 3, MeshData::SM_TOP_SWAP | MeshData::SM_LAPLACE_MIXED);
		MeshGenerator2dQuad::improveQuads(mc, mesh, 3, MeshData::SM_LAPLACE_MIXED);
	}catch(...){
		LOG4CPLUS_ERROR(MeshLog::logger_console,   "Error converting to quads.");
		return -1;
	}
//	double sec = (clock() - clock_start)/(double)CLOCKS_PER_SEC;
//	LOG4CPLUS_INFO(MeshLog::logger_console, "nodes", mesh->getPointsCount());
//	int tct = mesh->getElementsCount(3);
//	int qct = mesh->getElementsCount(4);
//	LOG4CPLUS_INFO(MeshLog::logger_console, "triangles", tct);
//	LOG4CPLUS_INFO(MeshLog::logger_console, "quads", qct);
//	LOG4CPLUS_INFO(MeshLog::logger_console, "Conversion speed (quads per sec.)", qct / sec);
//	LOG4CPLUS_INFO(MeshLog::logger_mesh, "[stat]\tmake-quads-NT:\t" << tct);
//	LOG4CPLUS_INFO(MeshLog::logger_mesh, "[stat]\tmake-quads-NQ:\t" << qct);
//	LOG4CPLUS_INFO(MeshLog::logger_mesh, "[stat]\tmake-quads-NQs:\t" << ((sec > 0.0)?(qct / sec):99999));
//	MeshView::showViewSet("Surface mesh", mesh->getViewSet(), 5);
	MeshStream::storeFileResuSurf(fname+"-quad.resu.surf", mesh);
	delete mesh;
	return 0;
}
*/

void MeshModel::setNameFromFilename(const string& fname)
{
	m_model_name = fname;
	size_t pos = m_model_name.rfind(".");
	if(pos >= 0) m_model_name.erase(pos);
	pos = m_model_name.find_last_of("/\\");
	if(pos >= 0) m_model_name.erase(0, pos+1);
}
