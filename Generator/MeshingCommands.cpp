/////////////////////////////////////////////////////////////////////////////
// MeshingCommands.cpp
// Class describing instructions for meshing
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2004-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#include "MeshingCommands.h"
#include "common.h"

#include "IteratorMesh2d.h"
#include "IteratorMesh3d.h"
#include "MeshContainer2d.h"
#include "MeshContainer3d.h"
#include "MeshDomainVolume.h"
#include "MeshDomainSurface.h"
#include "DEquation.h"
#include "MeshGenerator1d.h"
#include "MeshGenerator2d.h"
#include "MeshGenerator2dQuad.h"
#include "MeshGenerator3d.h"
#include "MeshGenerator3dQuality.h"
#include "MeshSpecialRoutinesDAT.h"
#include "MeshStream.h"
#include "ControlSpaceMesh.h"
#include "ParametricSurface.h"
#include "MeshEdge3d.h"
#include "SurfaceFunction.h"
#include "SurfacePlane.h"
#include "ControlSpaceQuadTree.h"
#include "EPSFile.h"
#include "MeshDecompositionUTC.h"
#include "DataStatistics.h"
#include "MeshGrain.h"

int MeshingCommands::execute(MeshContainer3d* &boundary, const char* command_line)
{
	istringstream str_text(command_line);
	string command, buffer, arg1;
	str_text >> command;
	command = MeshData::lowercase(command);

#ifdef USE_EXCEPTIONS
	try{
#endif
		int i = 0, from;
		double d;
		if(command.compare("help") == 0){
			TLOG("CUT_MESH <fname>");
//			TLOG("EDGE_INNER_NODES_2 <edge_node_ct> [nonlinear]");
//			TLOG("EDGE_INNER_NODES_3 <edge_node_ct> [nonlinear]");
			TLOG("EPS");
			TLOG("EXIT");
			TLOG("GEN1D [pts_ct=-1] [real3d=1]");
			TLOG("GEN1.5D");
			TLOG("GEN2D");
			TLOG("GEN2.5D");
			TLOG("GEN3D");
			TLOG("GET [property]");
			TLOG("LEELO [max_ct=-1]");
			TLOG("LOAD <mesh_file_name>");
			TLOG("LOAD_AND_PREPARE <mesh_file_name>");
			TLOG("LOAD_SURF <file_name>");
			TLOG("LOG_NAME <file_name>");
			TLOG("MAKE_QUADS <qmorph/leelo/mixed>");
			TLOG("MAKE_QUADS_FOR_PATCH <file.resu.surf> [qmorph/leelo/mixed]");
			TLOG("PARSE_GRAIN <dir> <output_file>");
			TLOG("PAUSE [info]");
			TLOG("QMORPH [max_ct=-1]");
			TLOG("QUALITY");
			TLOG("QUALITY_QUADS");
			TLOG("QUALITY3");
			TLOG("QUALITY3_DETAILED");
			TLOG("QUALITY_EXTENDED");
			TLOG("QUIT");
			TLOG("SET <property> <value>");
			TLOG("SMOOTHEN [ct=1]");
			TLOG("SMOOTHEN_DEL_SWAP [ct=1]");
			TLOG("SMOOTHEN_QUADS [ct=1]");
			TLOG("SMOOTHEN_QUADS_3D [ct=1]");
			TLOG("SMOOTHEN_METRIC [ct=1]");
			TLOG("SMOOTHEN_LAPLACE_METRIC [ct=1]");
			TLOG("SPECIAL_PHASE_RUN [hesjan=0] [from=0]");
			TLOG("SPLIT");
			TLOG("STATS");
			TLOG("STORE2.5 <fname>");
			TLOG("STORE2 <fname>");
			TLOG("STORE_AMIRA <fname>");
			TLOG("STORE_SURF <fname>");
			TLOG("STORE_PJM <fname> <GR_2D_Q4|GR_2D_Q8|GR_2D_Q9|GR_2D_T3|GR_2D_T6|GR_3D_T4|GR_3D_T10>");
			TLOG("STORE_DAT <fname>");
			TLOG("TEST_HESSIAN <fname> <equation> [equation2]");
			TLOG("TEST_HESSIAN_QUAD <fname> <equation> [equation2]");
			TLOG("TRIANGULATE [sm_count]");
			TLOG("TRIANGULATE3 [sm_count]");
			TLOG("TRIANGULATE3_SURFACE <fname>");
		//---------------------------------------------------
		}else if(command.compare("gen1d") == 0){
			if(!boundary){
				TERROR("Undefined boundary");
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
				TERROR("ERROR loading.");
		//---------------------------------------------------
		}else if(command.compare("pause") == 0){
			string info = "";
			str_text >> info;
			TLOG_S("pause, waiting for input...",info);
			cin >> info;
		//---------------------------------------------------
		}else if(command.compare("load") == 0){
			string fname = "";
			str_text >> fname;
			MESHLOGSTAT << "MeshName\t" << fname << endl;
			clock_t clock_start = clock();
			MeshContainer3d* domain = MeshStream::readFileMsh(fname.c_str());
			if(domain){
				if(boundary) delete boundary;
				boundary = domain;
				double sec = (clock() - clock_start)/(double)CLOCKS_PER_SEC;
				MESHLOG << "[stat]\tpreparation-time:\t" << sec << endl;
				TLOG_S("Loading time", sec);
			}else
				TERROR("ERROR loading.");
		//---------------------------------------------------
		}else if(command.compare("load_and_prepare") == 0){
			string fname = "";
			str_text >> fname;
			MESHLOGSTAT << "MeshName\t" << fname << endl;
			clock_t clock_start = clock();
			MeshContainer3d* domain = MeshStream::readFileMsh(fname.c_str());
			if(domain){
				if(boundary) delete boundary;
				boundary = domain;
				for(IteratorBoundary2d it = boundary->getFirstValidBoundary2d(); 
						it.isValid(); it.nextValidBoundary())
				{
					MeshDomainSurface* domain_surface = it.getDomainSurface();
					ControlSpace* ucs_2d = domain_surface->getUserControlSpace();
					MeshDomainVolume* volume0 = domain_surface->getDomainBlock(0);
					MeshDomainVolume* volume1 = domain_surface->getDomainBlock(1);
					ControlSpace3d* ucs_3d_0 = volume0 ? volume0->getUserControlSpace() : NULL;
					ControlSpace3d* ucs_3d_1 = volume1 ? volume1->getUserControlSpace() : NULL;
					it.getBoundary()->createControlSpace(ucs_2d, ucs_3d_0, ucs_3d_1);
				}
				double sec = (clock() - clock_start)/(double)CLOCKS_PER_SEC;
				MESHLOG << "[stat]\tpreparation-time:\t" << sec << endl;
				TLOG_S("Loading time", sec);
			}else
				TERROR("ERROR loading.");
		//---------------------------------------------------
		}else if(command.compare("log_name") == 0){
			string fname = "";
			str_text >> fname;
			MeshLog::setLogFile(fname.c_str());
		//---------------------------------------------------
		}else if(command.compare("parse_grain") == 0){
			string dir = "", output_file = "grain.msh";
			str_text >> dir >> output_file;
			MeshGrain::parseGrainFiles(dir, output_file);
		//---------------------------------------------------
/*		}else if(command.compare("edge_inner_nodes_2") == 0){
			if(!boundary){	TERROR("Undefined boundary");	return CM_ERROR_NOMESH;	}
			str_text >> i;
			if(!str_text) i = 0;
			bool use_shapes = false;
			str_text >> arg1;
			if(str_text && (arg1.compare("nonlinear") == 0)) use_shapes = true;
			for(IteratorMesh2d it = boundary->getFirstValidMesh2d(); it.isValid(); it.nextValidMesh()){
				MetricContext mc(it.getMesh()->getControlSpace());
				it.getMesh()->setEdgeInnerPoints(mc, i, use_shapes);
			}
		//---------------------------------------------------
		}else if(command.compare("edge_inner_nodes_3") == 0){
			if(!boundary){	TERROR("Undefined boundary");	return CM_ERROR_NOMESH;	}
			str_text >> i;
			if(!str_text) i = 0;
			bool use_shapes = false;
			str_text >> arg1;
			if(str_text && (arg1.compare("nonlinear") == 0)) use_shapes = true;
			for(IteratorMesh3d it = boundary->getFirstValidMesh3d(); it.isValid(); it.nextValidMesh())
				it.getMesh()->addEdgeInnerPoints(i, use_shapes);
*/		//---------------------------------------------------
		}else if(command.compare("test_hessian") == 0){
			str_text >> buffer;	//
			if(!str_text){
				TERROR("mesh file name required.");
				return CM_ERROR_PARSE;
			}
			string eq1, eq2 = "";
			str_text >> eq1;	// "2*sin(3*x)"
			if(!str_text){
				TERROR("equation string required.");
				return CM_ERROR_PARSE;
			}
			str_text >> eq2;	// "2*sin(3*x)"

			MeshSpecialRoutinesDAT::testHessianCurvature(buffer, eq1, eq2);
		//---------------------------------------------------
		}else if(command.compare("test_hessian_quad") == 0){
			str_text >> buffer;	//
			if(!str_text){
				TERROR("mesh file name required.");
				return CM_ERROR_PARSE;
			}
			string eq1, eq2 = "";
			str_text >> eq1;	// "2*sin(3*x)"
			if(!str_text){
				TERROR("equation string required.");
				return CM_ERROR_PARSE;
			}
			str_text >> eq2;	// "2*sin(3*x)"

			MeshSpecialRoutinesDAT::testHessianCurvature(buffer, eq1, eq2, true);
		//---------------------------------------------------
		}else if(command.compare("cut_mesh") == 0){
			str_text >> buffer;	//
			if(!str_text){
				TERROR("script file name required.");
				return CM_ERROR_PARSE;
			}
			MeshDecompositionUTC cut_task(buffer);
			cut_task.run();
		//---------------------------------------------------
		}else if(command.compare("store_dat") == 0){
			if(!boundary){	TERROR("Undefined boundary");	return CM_ERROR_NOMESH;	}
			str_text >> buffer;
			string fname = str_text ? buffer : "mesh";
			for(IteratorMesh2d it = boundary->getFirstValidMesh2d(); it.isValid(); it.nextValidMesh())
				MeshSpecialRoutinesDAT::exportToDAT(it.getMesh(), fname, true);
		//---------------------------------------------------
		}else if(command.compare("store2") == 0){
			if(!boundary){	TERROR("Undefined boundary");	return CM_ERROR_NOMESH;	}
			str_text >> buffer;
			string fname = str_text ? buffer : "mesh";
			for(IteratorMesh2d it = boundary->getFirstValidMesh2d(); it.isValid(); it.nextValidMesh()){
				it.getMesh()->storeTxt(fname, i++);
			}
		//---------------------------------------------------
		}else if(command.compare("store_surf") == 0){
			if(!boundary){	TERROR("Undefined boundary");	return CM_ERROR_NOMESH;	}
			str_text >> buffer;
			string fname = str_text ? buffer : "mesh-surf";
			for(IteratorMesh2d it = boundary->getFirstValidMesh2d(); it.isValid(); it.nextValidMesh()){
				bool proper_orientation = it.getDomainVolume()->properOrientation((MeshFace*)it.getDomainSurface());
				it.getMesh()->storeSurf(fname, proper_orientation, i++);
			}
		//---------------------------------------------------
		}else if(command.compare("store_amira") == 0){
			if(!boundary){	TERROR("Undefined boundary");	return CM_ERROR_NOMESH;	}
			str_text >> buffer;
			string fname = str_text ? buffer : "mesh";
			if(!boundary->getFirstValidBoundaryMesh3d().isValid())
				MeshGenerator3d::prepareBoundaryMesh(boundary);

			for(IteratorBoundaryMesh3d it = boundary->getFirstValidBoundaryMesh3d(); it.isValid(); it.nextValidBoundaryMesh())
				it.getDomainVolume()->storeAmiraMesh(fname, i++);
		//---------------------------------------------------
		}else if(command.compare("store_pjm") == 0){
			if(!boundary){ TERROR("Undefined boundary");	return CM_ERROR_NOMESH;	}
			string fname;
			str_text >> fname;
			if(!str_text){ TERROR("Grid File Name required!"); return CM_ERROR_PARSE; }
			str_text >> arg1;
			if(!str_text){ TERROR("Grid Type required!"); return CM_ERROR_PARSE; }
			int grid_type;
			for(size_t k = 0; k < arg1.length(); k++)
				arg1[k] = tolower(arg1[k]);
			if(arg1.compare("gr_2d_q4")==0) grid_type = 1;
			else if(arg1.compare("gr_2d_q8")==0) grid_type = 3;
			else if(arg1.compare("gr_2d_q9")==0) grid_type = 5;
			else if(arg1.compare("gr_2d_t3")==0) grid_type = 7;
			else if(arg1.compare("gr_2d_t6")==0) grid_type = 8;
			else if(arg1.compare("gr_3d_t4")==0) grid_type = 25;
			else if(arg1.compare("gr_3d_t10")==0) grid_type = 26;
			else {TERROR("Unknown Grid Type!"); return CM_ERROR_PARSE; }
			if(grid_type < 25){ // 2d
				for(IteratorMesh2d it = boundary->getFirstValidMesh2d(); it.isValid(); it.nextValidMesh()){
					CParametricSurface* surface = it.getDomainSurface()->getBaseSurface();
					it.getMesh()->storePJM(fname, grid_type, i++, surface);
				}
			}else{
				for(IteratorMesh3d it = boundary->getFirstValidMesh3d(); it.isValid(); it.nextValidMesh()){
					it.getMesh()->storePJM(fname, grid_type, i++);
				}
			}
		//---------------------------------------------------
		}else if(command.compare("triangulate3") == 0){
			if(!boundary){	TERROR("Undefined boundary");	return CM_ERROR_NOMESH;	}
			int sm_count = 2;
			str_text >> sm_count;
			clock_t clock_start = clock();
			MeshGenerator2d::show_prediction = false;
			START_CLOCK("3D Triangulation");
			MeshGenerator3d::autoTriangulate(boundary, sm_count);
			STOP_CLOCK("3D Triangulation");
			double sec = (clock() - clock_start)/(double)CLOCKS_PER_SEC;
			int nct = 0, ect = 0;
			for(IteratorMesh3d it = boundary->getFirstValidMesh3d(); it.isValid(); it.nextValidMesh()){
				MeshContainer3d* mesh = it.getMesh();
				nct += mesh->getPointsCount();
				ect += mesh->getBlocksCount();
			}
			TLOG_S("nodes", nct);
			TLOG_S("tetrahedra", ect);
			TLOG_S("Generation time - 3D (sec.)", sec);
			MESHLOGSTAT << "Meshing3DTime-total\t" << sec << endl;
			TLOG_S("Generation speed (tetrahedra per sec.)", ect / sec);
		//---------------------------------------------------
		}else if(command.compare("triangulate3_surface") == 0){

			if(boundary) delete[] boundary;
			string fname;
			str_text >> fname;
			MESHLOGSTAT << "MeshName\t" << fname << endl;
			clock_t clock_start = clock();
			MeshGenerator2d::show_prediction = false;
			//MeshGenerator3d::prepareBoundaryMesh(boundary);
			boundary = MeshSpecialRoutinesDAT::loadSurfaceBoundaryMesh(fname);
			if(!boundary) { TERROR("Error loading boundary mesh"); return CM_ERROR_NOMESH; }
			MeshGenerator3d::createTetrahedra(boundary);
			MeshGenerator3dQuality::smoothenBlocks(boundary, 2);
			double sec = (clock() - clock_start)/(double)CLOCKS_PER_SEC;
			int nct = 0, ect = 0;
			for(IteratorMesh3d it = boundary->getFirstValidMesh3d(); it.isValid(); it.nextValidMesh()){
				MeshContainer3d* mesh = it.getMesh();
				nct += mesh->getPointsCount();
				ect += mesh->getBlocksCount();
			}
			TLOG_S("nodes", nct);
			TLOG_S("tetrahedra", ect);
			TLOG_S("Generation time - 3D (sec.)", sec);
			MESHLOGSTAT << "Meshing3DTime-total\t" << sec << endl;
			TLOG_S("Generation speed (tetrahedra per sec.)", ect / sec);
		//---------------------------------------------------
		}else if(command.compare("triangulate") == 0){
			if(!boundary){	TERROR("Undefined boundary");	return CM_ERROR_NOMESH;	}
			int sm_count = 2;
			str_text >> sm_count;
			clock_t clock_start = clock();
			MeshGenerator2d::show_prediction = true;
			MeshGenerator2d::autoTriangulate(boundary, sm_count);
			double sec = (clock() - clock_start)/(double)CLOCKS_PER_SEC;
			int nct = 0, tct = 0;
			for(IteratorMesh2d it = boundary->getFirstValidMesh2d(); it.isValid(); it.nextValidMesh()){
				MeshContainer2d* mesh = it.getMesh();
				nct += mesh->getPointsCount();
				tct += mesh->getElementsCount();
			}
			TLOG_S("nodes", nct);
			TLOG_S("triangles", tct);
			double speed = ((sec > 0.0)?(tct / sec):99999);
			TLOG_S("Triangulation speed (triangles per sec.)", speed);
			MESHLOG << "[stat]\ttriangulate-NT:\t" << tct << endl;
			MESHLOG << "[stat]\ttriangulate-NTs:\t" << speed << endl;
		//---------------------------------------------------
		}else if(command.compare("smoothen") == 0){
			str_text >> i;	if(!str_text) i = 1;
			clock_t clock_start = clock();
			MeshGenerator2d::smoothenFaces(boundary, i,
				MeshData::SM_LAPLACE_MIXED | MeshData::SM_TOP_SWAP | MeshData::SM_DEL_SWAP_COMPLETE);
			double sec = (clock() - clock_start)/(double)CLOCKS_PER_SEC;
			TLOG_S("Smoothing time (sec.)", sec);
		//---------------------------------------------------
		}else if(command.compare("smoothen_metric") == 0){
			str_text >> i;	if(!str_text) i = 1;
			clock_t clock_start = clock();
			MeshGenerator2d::smoothenFaces(boundary, i, MeshData::SM_METRIC);
			double sec = (clock() - clock_start)/(double)CLOCKS_PER_SEC;
			TLOG_S("Smoothing time (sec.)", sec);
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
			if(!boundary){	TERROR("Undefined boundary");	return CM_ERROR_NOMESH;	}
			DataStatistics mean_ratio_stats;
			DataStatistics edge_stats;
			DataStatistics nonconf_stats;
			for(IteratorMesh2d it = boundary->getFirstValidMesh2d(); it.isValid(); it.nextValidMesh()){
				MeshContainer2d* mesh = it.getMesh();
				ControlSpace* space = mesh->getControlSpace();
				MetricContext mc(space);
				mesh->statMeanRatio(mc, mean_ratio_stats, false);
				mesh->statMetricEdgeLength(mc, edge_stats, false);
				mesh->statMetricDifference(mc, nonconf_stats);
			}
			if(mean_ratio_stats.calculate()){
				TLOG_S("Metric mean ratio 1p (ave)", mean_ratio_stats.average());
				TLOG_S("-------------------- (dev)", mean_ratio_stats.stdDev());
			}
			if(edge_stats.calculate()){
				TLOG_S("Metric edge lengths 1p (ave)", edge_stats.average());
				TLOG_S("---------------------- (dev)", edge_stats.stdDev());
			}				
			if(nonconf_stats.calculate()){
				TLOG_S("Metric difference 1p (ave)", nonconf_stats.average());
				TLOG_S("-------------------- (dev)", nonconf_stats.stdDev());
			}				
		//---------------------------------------------------
		}else if(command.compare("quality3") == 0){
			if(!boundary){	TERROR("Undefined boundary");	return CM_ERROR_NOMESH;	}
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
			MESHLOGSTAT << "NP\t" << np << endl;
			MESHLOGSTAT << "NT\t" << nt << endl;
			if(mean_ratio_stats.calculate()){
				TLOG_S("Metric mean ratio (ave)", mean_ratio_stats.average());
				MESHLOGSTAT << "MetricMeanRatio-ave\t" << mean_ratio_stats.average() << endl;
				TLOG_S("Metric mean ratio (min)", mean_ratio_stats.minimum());
				MESHLOGSTAT << "MetricMeanRatio-min\t" << mean_ratio_stats.minimum() << endl;
				TLOG_S("------------------(dev)", mean_ratio_stats.stdDev());
				MESHLOGSTAT << "MetricMeanRatio-dev\t" << mean_ratio_stats.stdDev() << endl;
			}
			if(edge_stats.calculate()){
				TLOG_S("Metric edge lengths (ave)", edge_stats.average());
				MESHLOGSTAT << "MetricEdgeLength-ave\t" << edge_stats.average() << endl;
				TLOG_S("--------------------(dev)", edge_stats.stdDev());
				MESHLOGSTAT << "MetricEdgeLength-dev\t" << edge_stats.stdDev() << endl;
			}				
			if(nonconf_stats.calculate()){
				TLOG_S("Metric difference  (ave)", nonconf_stats.average());
				MESHLOGSTAT << "MetricDiff-ave\t" << nonconf_stats.average() << endl;
				TLOG_S("------------------ (dev)", nonconf_stats.stdDev());
				MESHLOGSTAT << "MetricDiff-dev\t" << nonconf_stats.stdDev() << endl;
			}				
		//---------------------------------------------------
		}else if(command.compare("quality_quads") == 0){
			if(!boundary){	TERROR("Undefined boundary");	return CM_ERROR_NOMESH;	}
			DataStatistics angle_stats;
			DataStatistics edge_stats;
			DataStatistics alpha_stats;
			for(IteratorMesh2d it = boundary->getFirstValidMesh2d(); it.isValid(); it.nextValidMesh()){
				MeshContainer2d* mesh = it.getMesh();
				ControlSpaceMesh* space = (ControlSpaceMesh*)mesh->getControlSpace();
				CParametricSurface* surface = it.getDomainSurface()->getBaseSurface();
				MetricContext mc(space);
				mesh->statMetricAlphaQuality(mc, alpha_stats, false);
				mesh->statMetricEdgeLength(mc, edge_stats, false);
				mesh->statRealAngles(angle_stats, surface);
			}
			if(alpha_stats.calculate()){
				TLOG_S("Metric alpha quality (ave-geom)", alpha_stats.averageGeometric());
				MESHLOG << "[stat] Metric-alpha-ave-geom\t" << alpha_stats.averageGeometric() << endl;
			}
			if(edge_stats.calculate()){
				TLOG_S("Metric edge lengths  (ave)", edge_stats.average());
				TLOG_S("-------------------- (dev)", edge_stats.stdDev());
				MESHLOG << "[stat] Metric-edge-ave\t" << edge_stats.average() << endl;
				MESHLOG << "[stat] Metric-edge-stddev\t" << edge_stats.stdDev() << endl;
			}				
			if(angle_stats.calculate()){
				double angles_min = angle_stats.minimum() * 180.0 / PI;
				double angles_30 = angle_stats.getDataCountBottom(PI/6) / (double)angle_stats.getCount();
				TLOG_S("Angles (min)", angles_min);
				TLOG_S("------ (<30)", angles_30);
				MESHLOG << "[stat] Angles-min\t" << angles_min << endl;
				MESHLOG << "[stat] Angles-below30\t" << angles_30 << endl;
			}				
		//---------------------------------------------------
		}else if(command.compare("quality3_detailed") == 0){
			if(!boundary){	TERROR("Undefined boundary");	return CM_ERROR_NOMESH;	}
			for(IteratorMesh3d it = boundary->getFirstValidMesh3d(); it.isValid(); it.nextValidMesh()){
				MeshContainer3d* mesh = it.getMesh();
				// show block quality
				TLOG("Quality: tetrahedra aspect ratio (~S/V)");
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
					TLOG(text.str());
				}
				if(total_block_count > 0){
					ostringstream text;
					text << "Min = " << stat_data[4].minimum << ", Max = " << stat_data[4].maximum;
					text << ", Ave = " << stat_data[4].average;
					TLOG(text.str());
					const char* str_types[] = {"Inner:    ", "B-Vertex: ", "B-Edge:   ", "B-Face:   "};
					for(i = 0; i < 4; i++){
						ostringstream text_str;
						text_str << str_types[i] << block_types[i] << " Min = ";
						text_str << stat_data[i].minimum << ", Max = " << stat_data[i].maximum;
						text_str << ", Ave = " << stat_data[i].average;
						TLOG(text_str.str());
					}
				}
				TLOG("----------------------------------------");
				// show edge metric quality
				TLOG("Quality: metric edge length");
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
					TLOG(text.str());
				}
			}
		//---------------------------------------------------
		}else if(command.compare("stats") == 0){
			if(!boundary){	TERROR("Undefined boundary");	return CM_ERROR_NOMESH;	}
			for(IteratorMesh2d it = boundary->getFirstValidMesh2d(); it.isValid(); it.nextValidMesh()){
				MeshContainer2d* mesh = it.getMesh();
				TLOG("------------------------------");
				TLOG_S("Number of nodes", mesh->getPointsCount());
				TLOG_S("Number of triangles", mesh->getElementsCount(3));
				TLOG_S("Number of quads", mesh->getElementsCount(4));
				ControlSpaceMesh* space = (ControlSpaceMesh*)mesh->getControlSpace();
				if(space && (space->getType() == MeshData::CONTROL_MESH)){
					TLOG_S("Control mesh nodes count", space->getControlNodesCount());
					TLOG_S("Control mesh triangles count", space->getTrianglesCount());
//					TLOG_S("Control mesh counter (simple)", space->getCounter(0));
//					TLOG_S("Control mesh counter (triangle)", space->getCounter(1));
//					TLOG_S("Control mesh counter (voronoi)", space->getCounter(2));
//					space->clearCounters();
				}
			}
		//---------------------------------------------------
		}else if(command.compare("quality_extended") == 0){
			if(!boundary){	TERROR("Undefined boundary");	return CM_ERROR_NOMESH;	}
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
			if(!boundary){	TERROR("Undefined boundary");	return CM_ERROR_NOMESH;	}
			i = 0;
			for(IteratorMesh2d it = boundary->getFirstValidMesh2d(); it.isValid(); it.nextValidMesh()){
				it.getDomainSurface()->storeEPS(i);
				MeshContainer2d* mesh = it.getMesh();
				if(mesh){
					ControlSpaceAdaptive* space = (ControlSpaceAdaptive*)mesh->getControlSpace();
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
		}else if(command.compare("store2.5") == 0){
			if(!boundary){	TERROR("Undefined boundary");	return CM_ERROR_NOMESH;	}
			str_text >> buffer;
			string fname = str_text ? buffer : "mesh";
			for(IteratorBoundaryMesh3d it = boundary->getFirstValidBoundaryMesh3d(); it.isValid(); it.nextValidBoundaryMesh())
				it.getDomainVolume()->storeSurfaceMeshTxt(fname, i++);
		//---------------------------------------------------
		}else if(command.compare("set") == 0){
			string name;
			str_text >> name;
			i = mesh_data.findProperty(name);
			if(i >= 0){
				if(mesh_data.getPropertyType(i) == MeshData::PROP_INT){
					str_text >> i;
					if(str_text) mesh_data.setPropertyInt(name, i);
					else TERROR("ERROR reading int value");
				}else{
					str_text >> d;
					if(str_text) mesh_data.setPropertyDouble(name, d);
					else TERROR("ERROR reading double value");
				}
			}else{
				TERROR("ERROR unknown property name");
			}
		//---------------------------------------------------
		}else if(command.compare("gen1.5d") == 0){
			if(!boundary){	TERROR("Undefined boundary");	return CM_ERROR_NOMESH;	}
			MeshGenerator2d::prepareBoundaryMesh(boundary);
		//---------------------------------------------------
		}else if(command.compare("gen2d") == 0){
			if(!boundary){	TERROR("Undefined boundary");	return CM_ERROR_NOMESH;	}
			MeshGenerator2d::triangulateFaces(boundary);
		//---------------------------------------------------
		}else if(command.compare("gen2.5d") == 0){
			if(!boundary){	TERROR("Undefined boundary");	return CM_ERROR_NOMESH;	}
			MeshGenerator3d::prepareBoundaryMesh(boundary);
		//---------------------------------------------------
		}else if(command.compare("gen3d") == 0){
			if(!boundary){	TERROR("Undefined boundary");	return CM_ERROR_NOMESH;	}
			MeshGenerator3d::createTetrahedra(boundary);
		//---------------------------------------------------
		}else if(command.compare("split") == 0){
			if(!boundary){	TERROR("Undefined boundary");	return CM_ERROR_NOMESH;	}
			MeshGenerator2dQuad::convertFacesToQuads(boundary, MeshData::QUADS_SPLIT);
		//---------------------------------------------------
		}else if(command.compare("leelo") == 0){
			str_text >> i;	if(!str_text) i = -1;
			if(!boundary){	TERROR("Undefined boundary");	return CM_ERROR_NOMESH;	}
			MeshGenerator2dQuad::convertFacesToQuads(boundary, MeshData::QUADS_LEELO, i);
		//---------------------------------------------------
		}else if(command.compare("test") == 0){
			//if(!boundary){	TERROR("Undefined boundary");	return CM_ERROR_NOMESH;	}
			test(boundary);
		//---------------------------------------------------
		}else if(command.compare("qmorph") == 0){
			str_text >> i;	if(!str_text) i = -1;
			if(!boundary){	TERROR("Undefined boundary");	return CM_ERROR_NOMESH;	}
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
			TLOG_S("nodes", nct);
			TLOG_S("triangles", tct);
			TLOG_S("quads", qct);
			TLOG_S("Conversion speed (quads per sec.)", qct / sec);
		//---------------------------------------------------
		}else if(command.compare("make_quads") == 0){
			if(!boundary){	TERROR("Undefined boundary");	return CM_ERROR_NOMESH;	}
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
			TLOG_S("nodes", nct);
			TLOG_S("triangles", tct);
			TLOG_S("quads", qct);
			TLOG_S("Conversion speed (quads per sec.)", qct / sec);
			MESHLOG << "[stat]\tmake-quads-NT:\t" << tct << endl;
			MESHLOG << "[stat]\tmake-quads-NQ:\t" << qct << endl;
			MESHLOG << "[stat]\tmake-quads-NQs:\t" << ((sec > 0.0)?(qct / sec):99999) << endl;
		//---------------------------------------------------
		}else if(command.compare("make_quads_for_patch") == 0){
			//TLOG("MAKE_QUADS_FOR_PATCH <file.resu.surf> [qmorph/leelo/mixed]");
			string fname;
			str_text >> fname;
			if(!str_text){ TERROR_S("Error reading file", fname);	return CM_ERROR_NOMESH;	}
			str_text >> arg1;
			MeshContainer2d* mesh =  // TODO read from file (+ planar surface)
				MeshStream::readFileResuSurf(fname);
			mesh->createControlSpaceFromTriangles(1.075); // create control space
			mesh->countMetricDifferenceQuality();
//			MeshView::setViewMode(SHOW_MESH_SURF_QUALITY);
//			MeshView::showViewSet("Surface mesh", mesh->getViewSet(), 2);
//			MeshView::setViewMode(SHOW_MESH_SURF);
			MetricContext mc(mesh->getControlSpace());
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
			TLOG_S("nodes", mesh->getPointsCount());
			int tct = mesh->getElementsCount(3);
			int qct = mesh->getElementsCount(4);
			TLOG_S("triangles", tct);
			TLOG_S("quads", qct);
			TLOG_S("Conversion speed (quads per sec.)", qct / sec);
			MESHLOG << "[stat]\tmake-quads-NT:\t" << tct << endl;
			MESHLOG << "[stat]\tmake-quads-NQ:\t" << qct << endl;
			MESHLOG << "[stat]\tmake-quads-NQs:\t" << ((sec > 0.0)?(qct / sec):99999) << endl;
			if(true){
				DataStatistics angle_stats;
				DataStatistics edge_stats;
				DataStatistics alpha_stats;
				mesh->statMetricAlphaQuality(mc, alpha_stats, false);
				mesh->statMetricEdgeLength(mc, edge_stats, false);
				mesh->statRealAngles(angle_stats, mesh->getSurface());
				if(alpha_stats.calculate()){
					TLOG_S("Metric alpha quality (ave-geom)", alpha_stats.averageGeometric());
					MESHLOG << "[stat] Metric-alpha-ave-geom\t" << alpha_stats.averageGeometric() << endl;
				}
				if(edge_stats.calculate()){
					TLOG_S("Metric edge lengths  (ave)", edge_stats.average());
					TLOG_S("-------------------- (dev)", edge_stats.stdDev());
					MESHLOG << "[stat] Metric-edge-ave\t" << edge_stats.average() << endl;
					MESHLOG << "[stat] Metric-edge-stddev\t" << edge_stats.stdDev() << endl;
				}				
				if(angle_stats.calculate()){
					double angles_min = angle_stats.minimum() * 180.0 / PI;
					double angles_30 = angle_stats.getDataCountBottom(PI/6) / (double)angle_stats.getCount();
					TLOG_S("Angles (min)", angles_min);
					TLOG_S("------ (<30)", angles_30);
					MESHLOG << "[stat] Angles-min\t" << angles_min << endl;
					MESHLOG << "[stat] Angles-below30\t" << angles_30 << endl;
				}				
			}
			// TODO store quads to file
//			MeshView::showViewSet("Surface mesh", mesh->getViewSet(), 5);
			MeshStream::storeFileResuSurf(fname+"-quad.resu.surf", mesh);
			delete mesh;
		//---------------------------------------------------
		}else if(command.compare("exit") == 0){
			return CM_QUIT;
		//---------------------------------------------------
		}else if(command.compare("quit") == 0){
			return CM_QUIT;
		//---------------------------------------------------
		}else if(command.compare("get") == 0){
			str_text >> arg1;
			if(str_text){
				// show 
				i = mesh_data.findProperty(arg1);
				if(i < 0){
					TERROR_S("Property not found", arg1);
				}else{
					if(mesh_data.getPropertyType(i) == MeshData::PROP_INT){
						TLOG_S(mesh_data.getPropertyName(i), mesh_data.getPropertyInt(i));
					}else{
						TLOG_S(mesh_data.getPropertyName(i), mesh_data.getPropertyDouble(i));
					}
					TLOG(mesh_data.getPropertyDescription(i));
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
					TLOG(out_text.str());
				}
			}
		//---------------------------------------------------
		}else{
			TERROR("ERROR parsing.");
		}
#ifdef USE_EXCEPTIONS
	}catch(MeshingException& e){
		TERROR("Unexpected meshing error");
		TERROR(e.getDescription());
		if(boundary) delete boundary;
		boundary = NULL;
		MESHLOG << "[stat]\t* Unexpected meshing error:" << endl;
		MESHLOG << "[stat]\t* " << e.getDescription() << endl;
		MESHLOG << "[stat]\t***" << endl;
		MESHLOG << "[stat]\t***" << endl;
		return CM_EXCEPTION;
	}catch(...){
		TERROR("Unexpected error");
		if(boundary) delete boundary;
		boundary = NULL;
		MESHLOG << "[stat]\t* Unexpected error:" << endl;
		MESHLOG << "[stat]\t***" << endl;
		MESHLOG << "[stat]\t***" << endl;
		MESHLOG << "[stat]\t***" << endl;
		return CM_EXCEPTION;
	}
#endif // USE_EXCEPTIONS

	return CM_OK;
}

int MeshingCommands::test(MeshContainer3d* boundary)
{
	if(!boundary) return 0;

	for(IteratorMesh3d it = boundary->getFirstValidMesh3d(); it.isValid(); it.nextValidMesh()){
		MeshContainer3d* mesh = it.getMesh();
		//ControlSpace3d* space = mesh->getControlSpace();
//		MeshView::showViewSet("boundary edges", mesh->getBoundaryEdgesViewSet());
		//Metric3dContext mc(space);
	}

	return 0;
}

int MeshingCommands::convertSurfPatchtoQuads(const string& fname)
{
	MeshContainer2d* mesh = MeshStream::readFileResuSurf(fname);
	if(!mesh){
		TERROR("Error reading file.");
		return -1;
	}
	mesh->createControlSpaceFromTriangles(1.075); // create control space
	mesh->countMetricDifferenceQuality();
//	MeshView::setViewMode(SHOW_MESH_SURF_QUALITY);
//	MeshView::showViewSet("Surface mesh", mesh->getViewSet(), 2);
//	MeshView::setViewMode(SHOW_MESH_SURF);
	MetricContext mc(mesh->getControlSpace());
//	clock_t clock_start = clock();			
	try{
		int qct = MeshGenerator2dQuad::convertToQuadsMixed(mc, mesh, MeshData::QUADS_QMORPH);
		if(qct <= 0){
			TERROR("Error converting to quads.");
			return -1;
		}

//		MeshView::showViewSet("Surface mesh", mesh->getViewSet(), 1);
		MeshGenerator2dQuad::improveQuads(mc, mesh, 3, MeshData::SM_TOP_SWAP | MeshData::SM_LAPLACE_MIXED);
		MeshGenerator2dQuad::improveQuads(mc, mesh, 3, MeshData::SM_LAPLACE_MIXED);
	}catch(...){
		TERROR("Error converting to quads.");
		return -1;
	}
//	double sec = (clock() - clock_start)/(double)CLOCKS_PER_SEC;
//	TLOG_S("nodes", mesh->getPointsCount());
//	int tct = mesh->getElementsCount(3);
//	int qct = mesh->getElementsCount(4);
//	TLOG_S("triangles", tct);
//	TLOG_S("quads", qct);
//	TLOG_S("Conversion speed (quads per sec.)", qct / sec);
//	MESHLOG << "[stat]\tmake-quads-NT:\t" << tct << endl;
//	MESHLOG << "[stat]\tmake-quads-NQ:\t" << qct << endl;
//	MESHLOG << "[stat]\tmake-quads-NQs:\t" << ((sec > 0.0)?(qct / sec):99999) << endl;
//	MeshView::showViewSet("Surface mesh", mesh->getViewSet(), 5);
	MeshStream::storeFileResuSurf(fname+"-quad.resu.surf", mesh);
	delete mesh;
	return 0;
}
