// MeshData.cpp: implementation of the MeshData class.
//
//////////////////////////////////////////////////////////////////////

#include "MeshData.h"

//#include <vld.h> // for Visual Leak Detector

#include <log4cplus/logger.h>
#include <log4cplus/loggingmacros.h>
using namespace log4cplus;

// for static parameters
#include "DPoint.h"
#include "ControlSpace2d.h"
#include "ControlSpace3d.h"
#include "ControlSpace2dMatrixUniform.h"
#include "ControlSpace3dMatrixUniform.h"
#include "ControlSpace2dQuadTree.h"
#include "ControlSpace3dOctree.h"
#include "ControlSpace2dMesh.h"
#include "MeshLog.h"
#include "QuadTree.h"
#include "OctTree.h"
#include "MeshGenerator1d.h"
#include "MeshGenerator2d.h"
#include "MeshGenerator2dQuad.h"
#include "MeshGenerator2dQMorph.h"
#include "MeshGenerator3d.h"
#include "MeshGenerator3dDelaunayBoundary.h"
#include "MeshGenerator3dQuality.h"
#include "MeshGenerator3dSurface.h"
#include "DHesjan.h"
#include "Metric2dContext.h"
#include "Metric3dContext.h"
#include "FrontEdge.h"
#include <iomanip>
#include "MeshViewSet.h"
#include "MeshGrain.h"
#include "MeshPoint3d.h"
#include "MeshContainer3dSurface.h"

DataVector<MeshingClock::ClockData> MeshingClock::mesh_clocks;

MeshData mesh_data;

/// meshlib version string
string MeshData::version()
{
	stringstream oss;
	oss << "meshlib v.2.4" << setw(3) << setfill('0') << 
		105 /* INCREASED_VERSION */
		<< ", " << __DATE__;
	return oss.str(); 
}

MeshData::MeshData()
	: m_parameters(128, "")
{
	m_model_diameter = 1.0;
	relative_small_number = SMALL_NUMBER;
	relative_infinity = LARGE_NUMBER;
	//--------
	m_view_polygon_density = 1.0;
	m_qmorph_min_ratio = 0.65;
	m_qmorph_max_ratio = sqrt(3.0);
	m_mesh_view_mode = 0;
	m_view_selected_blocks = false;

	//--- count epsilon
	double half = 0.5;
	double check = 1.0, lastcheck;
	m_epsilon = 1.0;
	do{
		lastcheck = check;
		m_epsilon *= half;
		check = 1.0 + m_epsilon;
	}while((check != 1.0) && (check != lastcheck));

	setDefaultParameters();
}

void MeshData::setDefaultParameters()
{
//	addPropertyDouble("GLVIEW_SHRINK", MeshViewSet::param_shrink, 0.9);

	//addPropertyInt("GRAIN_VIEW", MeshGrain::param_grain_view, 0,
	//	"View flags for multigrain parsing:\n"
	//	"1 - view node-fit faces\n"
	//	"2 - view fit/invalid faces\n"
	//	"4 - view parsed (initial) face data");

	//addPropertyDouble("GRAIN_TOLERANCE", MeshGrain::param_grain_tolerance, 2.0,
	//	"Tolerance for geometric approximation (plane, quadric, etc.)");

	//addPropertyDouble("GRAIN_NODE_TOLERANCE", MeshGrain::param_grain_node_identity_tolerance, 1.0,
	//	"Tolerance for how close can nodes be to be treated as different vertices");

	//addPropertyInt("EXTRA_DCMP_MESH", MeshGenerator2d::param_mesh_decomposition, 0,
	//	"Extra decomposition mode after boundray nodes triangulation\n"
	//	"0 - no decomposition (default)\n"
	//	"1 - special control space as created copy (2d)\n"
	//	"2 - special control space as filtered layer (2d)\n"
	//	"11 - special control space as created copy (3d)\n"
	//	"12 - special control space as filtered layer (3d)");
	//addPropertyDouble("EXTRA_DCMP_MESH_WIDTH", MeshGenerator2d::param_mesh_decomposition_width, 1.0,
	//	"Dense layer width for extra decomposition after boundray nodes triangulation");
	//addPropertyDouble("EXTRA_DCMP_MESH_GRADATION", MeshGenerator2d::param_mesh_decomposition_gradation, 8.0,
	//	"Metric gradation for extra decomposition after boundray nodes triangulation\n"
	//	"< 0: no metric smoothing");

	//addPropertyInt("GEN_NODES_MAX", MeshPoint2d::param_max_node_count, 10000000,
	//	"Maximum number of nodes for the created mesh");
	//addPropertyInt("GEN_NODES_INNER", MeshGenerator2d::param_triangulate_with_inner_nodes, 1,
	//	"Whether to insert inner nodes during Delaunay triangulation");

	//addPropertyDouble("QUALITY3D_REFINED_RATIO", MeshGenerator3dDelaunayBoundary::param_refined_ICD_ratio, 0.15,
	//	"Quality ratio for special refined incremental CD (test version)\n"
	//	" - blocks with (mean ratio in metric) lower quality are improved\n"
	//	" - 0 or below -> no refinement");

//	addPropertyInt("TRIANG_QTREE_THRESHOLD", QuadTree::param_qtree_threshold, 300,
//		"Threshold for quad tree structure for containing triangle search\n"
//		"(0 - no quad tree used\n");

//	addPropertyInt("TRIANG_QTREE_ELEMENT_POINT_DISTANCE", MeshElement::param_point_distance, MeshData::RP_VERTICES, 
//		"Which points of element should be used for point-element distance checking:\n"
//		"1 - element middle (flag)\n"
//		"2 - vertices (flag)\n"
//		"4 - middles of edges (flag)");

//	addPropertyInt("TRIANG_TYPE", MeshGenerator2d::param_triangulation_type, MeshData::TRIANGULATE_CIRCLE,
//		"Delaunay retriangulation method:\n"
//		" 0 - containing circle\n"
//		" 1 - edge swap");

//	addPropertyInt("TRIANG_QUALITY_IMPROVE",
//		MeshGenerator2d::param_quality_improvement, IMPROVE_OUTER_CIRCLE,
//		"Method of triangle improvement (where to insert new node)\n"
//		" 0 - longest edge\n"
//		" 1 - circumcentre\n"
//		" 2 - longest_edge and circumcentre");

	//addPropertyDouble("GEN2D_QUALITY_THRESHOLD", MeshGenerator2d::param_quality_threshold, 0.58,
	//	"Quality threshold for triangle improvement");

	//addPropertyInt("TRIANG_QUALITY_CRITERION",
	//	MeshTriangle2d::param_quality_criterion, QUALITY_CIRCLE_SPACE_AND_EDGES,
	//	"Method of assessing the quality of triangle:"
	//	" 0 - QUALITY_RADIA\n"
	//	" 1 - QUALITY_ALPHA\n"
	//	" 2 - QUALITY_QUOTIENT_r_p\n"
	//	" 3 - QUALITY_QUOTIENT_AREAS\n"
	//	" 4 - QUALITY_QUOTIENT_ppp\n"
	//	" 5 - QUALITY_QUOTIENT_r_h\n"
	//	" 6 - QUALITY_SPACE\n"
	//	" 7 - QUALITY_CIRCLE_SPACE\n"
	//	" 8 - QUALITY_CIRCLE_SPACE_AND_EDGES\n"
	//	" 9 - QUALITY_SMALL_ANGLES");

	//addPropertyDouble("TRIANG_INITIAL_SIZE_RATIO", MeshGenerator2d::param_initial_mesh_size_ratio, 0.1,
	//	"Dimension ratio for the initial mesh size (2d)");

	addPropertyInt("GEN2D_MAX_AUTO_RETRIANGULATION",
		MeshGenerator2d::param_max_auto_retriangulation_count, 1,
		"Maximum number of retriangulation for automatic control space 2d adjustment");

	addPropertyInt("GEN3D_MAX_AUTO_RETRIANGULATION",
		MeshGenerator3d::param_max_auto_retriangulation_count, 1,
		"Maximum number of retriangulation for automatic control space 3d adjustment");

	//addPropertyInt("TETRA_TYPE", MeshGenerator3d::param_triangulation_type, TRIANGULATE_CIRCLE,
	//	"Delaunay retriangulation method:\n"
	//	" 0 - containing sphere\n"
	//	" 1 - edge/face swap");

	//addPropertyDouble("GEN3D_TETRA_QUALITY_IMPROVE",
	//	MeshGenerator3d::param_quality_improvement, 0.1,
	//	"Threshold for tetrahedra improvement (where to insert new node)\n"
	//	" below TQI - longest edge, above TQI - circumcenter\n"
	//	" -1.0 -> all circumcenter\n"
	//	"  1.0 -> all longest edge");

	//addPropertyDouble("GEN3D_QUALITY_THRESHOLD", 
	//	MeshGenerator3d::param_quality_threshold, 0.57,
	//	"Quality threshold for tetrahedra improvement");

	//addPropertyInt("TETRA_QUALITY_CRITERION",
	//	MeshTetrahedron::param_quality_criterion, 
	//	MeshData::QUALITY3D_CIRCLE_SPACE_AND_EDGES,
	//	"Method of assessing the quality of tetrahedron:"
	//	" 0 - QUALITY3D_SPACE\n"
	//	" 1 - QUALITY3D_CIRCLE_SPACE\n"
	//	" 2 - QUALITY3D_ASPECT_RATIO\n"
	//	" 3 - QUALITY3D_RADIUS_RATIO\n"
	//	" 4 - QUALITY3D_MEAN_RATIO\n"
	//	" 5 - QUALITY3D_CIRCLE_SPACE_AND_EDGES");

	//addPropertyInt("TETRA_OCTREE_THRESHOLD", OctTree::param_octtree_threshold, 1000,
	//	"Threshold for oct tree structure for containing tetrahedron search\n"
	//	"(0 - no oct tree used\n");

	//addPropertyInt("TETRA_OCTREE_BLOCK_POINT_DISTANCE", 
	//	MeshBlock::param_point_distance, MeshData::RP_VERTICES, 
	//	"Which points of block should be used for point-element distance checking:\n"
	//	"1 - element middle (flag)\n"
	//	"2 - vertices (flag)\n"
	//	"4 - middles of edges (flag)\n"
	//	"8 - middles of faces (flag)");

	//addPropertyDouble("TETRA_INITIAL_SIZE_RATIO",  
	//	MeshGenerator3d::param_initial_mesh_size_ratio, 0.25,
	//	"Dimension ratio for the initial mesh size (3d)");

	//addPropertyInt("TRIANG_SWAP", MeshGenerator2d::param_swap_criterion, MeshData::SWAP_ANGLES,
	//	"Criterion used for swapPIng method of quality improvement:\n"
	//	" 0 - SWAP_ANGLES\n"
	//	" 1 - SWAP_ALPHA\n"
	//	" 2 - SWAP_RR\n"
	//	" 3 - SWAP_RP\n"
	//	" 4 - SWAP_SCST\n"
	//	" 5 - SWAP_PPP\n"
	//	" 6 - SWAP_RH");

	//addPropertyInt("GEN2D_SWAP_TRIANG_MAX", 
	//	MeshGenerator2d::param_swap_maximum, 100000,
	//	"Maximum number of edge swap per single point retriangulation with inner angles criterion");

	//addPropertyInt("QUALITY3D_SWAP_CRITERION",
	//	MeshGenerator3dQuality::param_swap_criterion, MeshData::SWAP3_MIN_VOLUME,
	//	"Criterion used for swapPIng method of quality improvement:\n"
	//	" 0 - SWAP3_ALWAYS\n"
	//	" 1 - SWAP3_MIN_VOLUME\n"
	//	" 2 - SWAP3_MIN_TETRAHEDRON_QUALITY\n"
	//	" 3 - SWAP3_AVE_TETRAHEDRON_QUALITY");

	//addPropertyInt("GEN3D_BOUNDARY_NODES_FOR_RECOVERY",
	//	MeshGenerator3dDelaunayBoundary::param_boundary_nodes_for_recovery, 1,
	//	"Whether additional boundary nodes should be inserted\n"
	//	"at the boundary for recovery of boundary edges");

	//addPropertyInt("GEN3D_BOUNDARY_TRIANGULATION_METHOD",
	//	MeshGenerator3d::param_boundary_triangulation_method, 0,
	//	"Method for creating the boundary constrained mesh:\n"
	//	" 0 - mixed (convex star / Delaunay)\n"
	//	" 1 - Delaunay only");

	addPropertyDouble("MG3S_SHARP_EDGE_THRESHOLD", 
		MeshGenerator3dSurface::param_sharp_edge_threshold, 0.5,
		"Tolerance for approximation of sharp edges and corners for surface mesh");


	addPropertyDouble("MG3S_LOCAL_SHAPE_TOLERANCE", 
		MeshGenerator3dSurface::param_local_shape_tolerance, 0.15,
		"Tolerance for approximation of local surfaces and curves for surface mesh");

	addPropertyDouble("ACS_CURVATURE_RATIO", 
		ControlSpace2dAdaptive::param_curvature_ratio, 0.15,
		"Curvature ratio for automated sizing");

	addPropertyDouble("ACS3D_DIAMETER_MAX_RATIO", 
		ControlSpace3dAdaptive::param_max_diameter_ratio, 0.5,
		"Maximum element (3d) size as a model diameter ratio for automated sizing");

	addPropertyDouble("ACS2D_DIAMETER_MAX_RATIO", 
		ControlSpace2dAdaptive::param_max_diameter_ratio, 0.5,
		"Maximum element (2d) size as a model diameter ratio for automated sizing");

	addPropertyDouble("ACS2D_INNER_BOUNDARY_RATIO", 
		ControlSpace2dAdaptive::param_inner_boundary_ratio, 1.0,
		"Maximum element (2d) size as a model diameter ratio for automated sizing for inner boundary faces (relative to normal ratio)");

	addPropertyDouble("ACS_DIAMETER_MIN_RATIO", 
		ControlSpace2dAdaptive::param_min_diameter_ratio, 0.0001,
		"Minimum element size as a model diameter ratio for automated sizing");

	addPropertyDouble("ACS_MIN_LENGTH", 
		ControlSpace2dAdaptive::param_min_length, 0.0,
		"Minimum element size (length) for automated sizing");

	addPropertyDouble("ACS_STRETCH_MAX_RATIO", 
		ControlSpace2dAdaptive::param_stretch_max_ratio, 2.0,
		"Maximum stretching (element length ratio) for automated sizing");

	//addPropertyDouble("ACS_STRETCH_CONTOUR_MAX_RATIO", 
	//	ControlSpace2dAdaptive::param_contour_stretch_max_ratio, 2.0,
	//	"Maximum stretching (element length ratio) for automated sizing (contour curvature)");

	//addPropertyDouble("ACS_INFLATE_RATIO", 
	//	ControlSpace2d::param_inflate_box_factor, 0.0,
	//	"Ratio of inflating bounding box for control space");

	//addPropertyInt("ACS_MESH_WEIGHTS_SQUARED", 
	//	ControlSpace2dMesh::param_weights_squared, 1,
	//	"Whether the calculated weights for the metric value inside the triangle should be squared");

	//addPropertyInt("ACS_MESH_INNER_NODES", 
	//	ControlSpace2dMesh::param_inner_nodes, 0,
	//	"Whether inner nodes should be inserted during preparation of this control mesh\n"
	//	" -1 - unlimited\n"
	//	"  0 - none\n"
	//	"  n - maximum n nodes\n");

	//addPropertyDouble("CONTROL_MIXED_THRESHOLD", 
	//	ControlSpace2dMesh::param_mixed_threshold, 0.1,
	//	"Threshold for metric difference for assesing mesh-control interpolation type");

	//addPropertyInt("ACS2D_TYPE", 
	//	ControlSpace2d::param_control_type, CONTROL_QUADTREE,
	//	"Which control type should be used:\n"
	//	"  0 - uniform rectangular mesh,\n"
	//	"  1 - triangular background mesh,\n"
	//	"  2 - analytical,\n"
	//	"  3 - identity,\n"
	//	"  4 - quadtree\n");

	//addPropertyInt("ACS3D_TYPE", 
	//	ControlSpace3d::param_control_type, CONTROL_OCTREE_3D,
	//	"Which control type (3D) should be used:\n"
	//	"  0 - uniform rectangular mesh,\n"
	//	"  1 - triangular background mesh,\n"
	//	"  2 - analytical,\n"
	//	"  3 - identity\n"
	//	"  4 - octree\n");

	//addPropertyInt("ACS_MESH_INTERPOLATION", 
	//	ControlSpace2dMesh::param_interpolation_method, CONTROL_TRIANGLE,
	//	"Interpolation method for control mesh:\n"
	//	"  0 - no interpolation,\n"
	//	"  1 - triangle interpolation,\n"
	//	"  2 - voronoi interpolation\n"
	//	"  3 - mixed interpolation\n"
	//	"  4 - NEM layered interpolation\n");

	//addPropertyInt("UNIFORM_CONTROL_NX", 
	//	ControlSpace2dMatrixUniform::param_uniform_nx, 100,
	//	"Resolution for uniform control grid");

	//addPropertyInt("CONTROL_UNIFORM_NXY", 
	//	ControlSpace2dAdaptive::param_control_nxy, 100,
	//	"starting X*Y resolution for quadtree control/uniform probing, ect.");

	//addPropertyInt("ACS_QTREE_USE_MIDNODES", 
	//	ControlSpace2dQuadTree::QuadLeaf::param_use_midnodes, 0,
	//	"Whether control qtree mid-nodes should be used for metric interpolation");

	//addPropertyInt("UNIFORM_CONTROL3_NX", 
	//	ControlSpace3dMatrixUniform::param_uniform_nx, 30,
	//	"Resolution for uniform control grid");

	//addPropertyInt("ACS2D_USE_SURFACE_CURVATURE", 
	//	ControlSpace2dAdaptive::param_use_surface_curvature, 1,
	//	"Whether use surface curvature for automated control space creation"); 

	//addPropertyInt("ACS2D_USE_CONTOUR_CURVATURE", 
	//	ControlSpace2dAdaptive::param_use_contour_curvature, 1,
	//	"Whether use (boundary) contour curvature for automated control space creation"); 

	//addPropertyInt("GEN_EVEN_NODES", 
	//	MeshGenerator1d::param_even_node_count, 0, 
	//	"Whether to enforce even number of nodes for each disretized boundary edge");

	//addPropertyInt("GEN_SMOOTHEN_LAST_NODE", 
	//	MeshGenerator1d::param_smoothen_last_node, 5, 
	//	"Whether to smoothen discretization of edge (for last segment, which may be smaller than required)\n"
	//	" -1: move all nodes,\n"
	//	"  0: don't move anything,\n"
	//	"  n: move n nearest nodes.");

	//addPropertyDouble("MAX_BOUNDARY_METRIC_CONFORM_RATIO", 
	//	MeshGenerator1d::param_boundary_metric_conform_rato, 1.8, 
	//	"Maximum difference for metric-length at edge vertices during discretization");

	//addPropertyDouble("HESJAN_MIN_RATIO", 
	//	DHesjan::param_hesjan_min_ratio, 0.01);
	//addPropertyDouble("HESJAN_MAX_RATIO", 
	//	DHesjan::param_hesjan_max_ratio, 0.2);
	//addPropertyDouble("HESJAN_FACTOR", 
	//	DHesjan::param_hesjan_factor, 0.01);
	
	//addPropertyInt("GEN2D_QMORPH_SIZE_CHECK", 
	//	MeshGenerator2dQMorph::param_size_check, 0,
	//	"Whether to check size of quads created durign conversion with the control space");

	//addPropertyInt("QUALITY2D_QUAD_SMOOTHING", 
	//	MeshGenerator2dQuad::param_convert_smoothing_method, 1,
	//	"Which method of smoothing (tyPIcally Laplace or Laplace-metric) should be used\n"
	//	"   for relocating nodes during qmorph conversion:\n"
	//	" 0 - Laplace\n"
	//	" 1 - Laplace-metric\n");

	//addPropertyDouble("ACS_METRIC_THRESHOLD_DIFF", 
	//	ControlSpace2dAdaptive::param_threshold_diff, 0.2, 
	//	"Metric difference threshold for quadtree (and other adaptive) control space");

	//addPropertyInt("ACS_QTREE_MAX_DEPTH", 
	//	ControlSpace2dQuadTree::param_max_depth, 20, 
	//	"Maximum depth of quadtree control space");

	//addPropertyInt("ACS_QTREE_MAX_DEPTH_SURFACE_ADAPTATION", 
	//	ControlSpace2dQuadTree::param_max_depth_surface_adaptation, 10, 
	//	"Maximum depth of quadtree control space for surface curvature adaptation");

	//addPropertyInt("ACS2D_QTREE_BALANCE", 
	//	ControlSpace2dQuadTree::QuadLeaf::param_balance_level, 1, 
	//	"Balance level for quad tree structure");

	//addPropertyInt("ACS3D_OCTREE_BALANCE", 
	//	ControlSpace3dOctree::OctLeaf::param_balance_level, 1, 
	//	"Balance level for octree structure");

	//addPropertyInt("TRIANG_QUALITY_METRIC_SELECTION", 
	//	MeshTriangle2d::param_quality_metric_selection, QM_MIDDLE,
	//	"How should be metric selected for triangle quality assessment:\n"
	//	"  0 - count at middle,\n"
	//	"  1 - intersect from middle and vertices,\n"
	//	"  2 - average from middle and vertices,\n"
	//	"  3 - intersect from midedges,\n"
	//	"  4 - average from midedges,\n");

	//addPropertyInt("TETRA_QUALITY_METRIC_SELECTION", 
	//	MeshTetrahedron::param_quality_metric_selection, QM_MIDDLE,
	//	"How should be metric selected for tetrahedron quality assessment:\n"
	//	"  0 - count at middle,\n"
	//	"  1 - intersect from middle and vertices,\n"
	//	"  2 - average from middle and vertices,\n"
	//	"  3 - intersect from midedges,\n"
	//	"  4 - average from midedges,\n");

	//addPropertyInt("ACS_RADIAL_PARTS", 
	//	ControlSpace2dAdaptive::param_radial_parts, 8, 
	//	"Number of probe points for radial metric definition in adaptive control space");
	//
	//addPropertyInt("ACS_CACHED_PARAMETERIZATION_MATRIX", 
	//	ControlSpace2d::param_cached_parameterization_matrix, 2, 
	//	"Whether the approximated (cached) parameterization matrix should be used\n"
	//	"  0 - no,\n"
	//	"  1 - from leaves,\n"
	//	"  2 - from nodes,\n");

	addPropertyDouble("ACS_GRADATION_RATIO", 
		ControlSpace2dAdaptive::param_gradation_ratio, 2.0, 
		"Gradation ratio for control space");
	
	//addPropertyDouble("RECOVERY_VOLUME_THRESHOLD", 
	//	MeshGenerator3d::param_recovery_volume_threshold, 1e-5, 
	//	"Volume threshold during 3D boundary recovery for checking line-face crossing etc.");

	//addPropertyDouble("GEN3D_VOLUME_THRESHOLD", 
	//	MeshGenerator3d::param_volume_threshold, 1e-5, 
	//	"Volume threshold during 3D generation (retriangulation, swapPIng, etc.)");

	//addPropertyDouble("GEN3D_SPECIAL_METRIC_RATIO", 
	//	MeshGenerator3d::param_special_metric_ratio, 1e-2, 
	//	"Special metric ratio for boundary recovery (for avoiding edges in undesired directions)");

	//addPropertyInt("ACS3D_OCTREE_MAX_DEPTH", 
	//	ControlSpace3dOctree::param_max_depth, 20, 
	//	"Maximum depth of octree control space");

	//addPropertyDouble("ACS3D_OCTREE_REDUCE_THRESHOLD", 
	//	ControlSpace3dOctree::OctLeaf::param_octree_reduce_threshold, 0.2, 
	//	"Metric threshold ratio for reducing octree size (number of leaves)\n"
	//	" < 0.0 - no reduction.");

	//addPropertyDouble("GEN2D_QUAD_METRIC_GRADATION_THRESHOLD", 
	//	MeshGenerator2dQuad::param_metric_gradation_threshold, 2.0, 
	//	"Metric gradation threshold for mixed quad conversion");

	//addPropertyInt("GEN2D_QUAD_LEELO_MIN_LEVEL", 
	//	MeshGenerator2dQuad::param_leelo_minimum_level, 1, 
	//	"Minimum front edge level for mixed quad conversion, where LeeLo method is allowed");

	//addPropertyDouble("GEN2D_QUAD_METRIC_LENGTH_THRESHOLD", 
	//	MeshGenerator2dQuad::param_metric_length_threshold, 4.0, 
	//	"Metric edge length (squared) threshold for mixed quad conversion");

	addPropertyInt("VIEW_SHOW",
		MeshViewSet::param_show_visualization, 1, 
		"whether any visualization should be prepared/shown");

	addPropertyInt("SURFACE_DOMAIN_TYPE",
		MeshContainer3dSurface::param_surface_domain_type, MeshContainer3dSurface::SDOMAIN_FACES_PLANAR,
		"Local surface domain type: \n"
		" 0 - hull in parametric space,\n"
		" 1 - direct 3d faces,\n"
		" 2 - faces/grid in parametric space");

	addPropertyInt("LOCAL_SURFACE_IDENTIFICATION",
		MeshGenerator3dSurface::param_local_surface_identification_method, 
		MeshGenerator3dSurface::LocalSurfaceIdentification::LSI_VIA_NORMALS,
		"Local surface identification method: \n"
		" 0 - greedy, arbitrary sequence,\n"
		" 1 - high curvature first,\n"
		" 2 - via normals using quadric-qtree");

	addPropertyInt("CS_KDTREE_MAX_DEPTH",
		KdTree3dParams::paramMaxLevel, 20, "Maximum depth for ControlSpace kd-trees");

	addPropertyInt("CS_OCTREE_MAX_DEPTH",
		KdOctree3dParams::paramMaxLevel, 7, "Maximum depth for ControlSpace octrees");
}

void MeshData::addPropertyInt(const string& name, int& value, int initial_value, const string& desc)
{
	assert( !m_parameters.contains(name) );
	auto data = std::make_shared<ParameterData>();
	data->name = lowercase(name);
	data->is_int = true;
	data->value.i_ref = &value;
	value = initial_value;
	data->description = desc;
	m_parameters.insert(name, data);
}

void MeshData::addPropertyDouble(const string& name, double& value, double initial_value, const string& desc)
{
	assert(!m_parameters.contains(name));
	auto data = std::make_shared<ParameterData>();
	data->name = lowercase(name);
	data->is_int = false;
	data->value.d_ref = &value;
	value = initial_value;
	data->description = desc;
	m_parameters.insert(name, data);
}

void MeshData::setModelDiameter(double d)
{
	m_model_diameter = d;
	relative_small_number = d*SMALL_NUMBER;
	relative_infinity = d*LARGE_NUMBER;
}

string MeshData::lowercase(const string& s)
{
	basic_string <char>::size_type len = s.length();
	string t = s;
	for(basic_string <char>::size_type i = 0; i < len; i++)
		t[i] = (char)tolower(s[i]);
	return t;
}

void MeshingClock::stopClock(const string & capt) {
	clock_t now = clock();
	ClockData data;

	static auto logger_time = log4cplus::Logger::getInstance(LOG4CPLUS_TEXT("mesh.stat.time"));
	do {
		if (mesh_clocks.empty()) return;
		data = mesh_clocks.removeLast();
		if (data.caption != capt) {
			LOG4CPLUS_WARN(MeshLog::logger_mesh, 
				"Time: " << data.caption << " =>\t???");
		}
	} while (data.caption != capt);
	clock_t diff = now - data.clock_step;

	//LOG4CPLUS_INFO(MeshLog::logger_mesh_stat,
	//	"Time: " << capt << " (sec.) =>\t" << (diff / (double)CLOCKS_PER_SEC));
	LOG4CPLUS_INFO(logger_time,
		capt << " (sec.) =>\t" << (diff / (double)CLOCKS_PER_SEC));
}




