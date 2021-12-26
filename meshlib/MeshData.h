/////////////////////////////////////////////////////////////////////////////
// MeshData.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2002-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#if !defined(MESHDATA_H__INCLUDED)
#define MESHDATA_H__INCLUDED

#include "common.h"
#include <log4cplus/initializer.h>

#include <time.h>
#include "DataVector.h"
#include "DataHashTable.h"

/// Struct used for (mostly quality-based) statistical data
struct StatData{
	StatData() : minimum(0.0), maximum(0.0), average(0.0) {}
	double minimum;
	double maximum;
	double average;
};

#define START_CLOCK(text) MeshingClock::startClock(text)
#define STOP_CLOCK(text) MeshingClock::stopClock(text)

class MeshingClock {
private:
	struct ClockData{
		ClockData(const string & capt = "", clock_t step = 0) : caption(capt), clock_step(step) {}
		string caption;
		clock_t clock_step;
	};
public:
	static void startClock(const string& capt = "?"){ 
		mesh_clocks.add(ClockData(capt, clock())); 
	}
	static void stopClock(const string& capt = "?");
private:
	static DataVector<ClockData> mesh_clocks;
};

class IndexedObject
{
public:
	IndexedObject() : index(0) {}
public:
	/// Returns the container-index 
	int getIndex() const { return index; }
	/// Sets the container-index 
	void setIndex(int i){ index = i; }
	/// Compares two objects
	short compareTo(const IndexedObject* ) const { assert(false); return 0; }
protected:
	/// Container-index
	int		index;
};

class RandomGen 
{
	unsigned long long u,v,w;
public:
	/// Constructor. Call with any integer seed (except value of v above).
	RandomGen(unsigned long long j) : v(4101842887655102017LL), w(1) {
		u = j ^ v; int64();
		v = u; int64();
		w = v; int64();
	}
	/// Return 64-bit random integer. 
	inline unsigned long long int64() {
		u = u * 2862933555777941757LL + 7046029254386353087LL;
		v ^= v >> 17; v ^= v << 31; v ^= v >> 8;
		w = 4294957665U*(w & 0xffffffff) + (w >> 32);
		unsigned long long x = u ^ (u << 21); x ^= x >> 35; x ^= x << 4;
		return (x + v) ^ w;
	}
	/// Return random double-precision floating value in the range 0. to 1.
	inline double doub() { return 5.42101086242752217E-20 * (double)int64(); }
	/// Return random double-precision floating value in the range a to b.
	inline double doub( const double& a, const double & b) { return a + (b-a)*doub(); }
	/// Return 32-bit random integer.
	inline unsigned int int32() { return (unsigned int)int64(); }
	/// Return random double-precision floating value in the range a to b.
	inline unsigned int int32( unsigned int a, unsigned int b) { return a + int32()%(b-a); }
};

/**
 * This class is used for storing some parameteters used by the whole application.
 *
 */
class MeshData  
{
public:
	/// Criteria of quality used for evaluation of triangles
	enum QualityCriterion {QUALITY_RADIA = 0, QUALITY_ALPHA = 1,
		QUALITY_QUOTIENT_r_p = 2, QUALITY_QUOTIENT_AREAS = 3,
		QUALITY_QUOTIENT_ppp = 4, QUALITY_QUOTIENT_r_h = 5,
		QUALITY_SPACE = 6, QUALITY_CIRCLE_SPACE = 7,
		QUALITY_CIRCLE_SPACE_AND_EDGES = 8, QUALITY_SMALL_ANGLES = 9};
	/// Criteria of quality used for evaluation of tetrahedra
	enum Quality3dCriterion {QUALITY3D_SPACE, QUALITY3D_CIRCLE_SPACE, 
		QUALITY3D_ASPECT_RATIO, QUALITY3D_RADIUS_RATIO, QUALITY3D_MEAN_RATIO,
		QUALITY3D_CIRCLE_SPACE_AND_EDGES};
	/// Criteria of quality used for visualization of triangles
	enum QualityViewCriterion {QVIEW_NONE, QVIEW_ALPHA_METRIC, QVIEW_ALPHA,
		QVIEW_ANGLES_METRIC, QVIEW_ANGLES, QVIEW_METRIC_LIKE, QVIEW_SPACE_METRIC};
	/// Criteria used for smoothing of triangles by edge-swapPIng
	enum SwapCriterions {SWAP_ANGLES = 0, SWAP_ALPHA = 1, SWAP_RR = 2, SWAP_RP = 3, SWAP_SCST = 4,
		SWAP_PPP = 5, SWAP_RH = 6};
	/// Criteria used for smoothing of tetrahedra by swapPIng (32 or 23)
	enum Swap3Criterions {SWAP3_ALWAYS = 0, SWAP3_MIN_VOLUME, 
		SWAP3_MIN_TETRAHEDRON_QUALITY, SWAP3_AVE_TETRAHEDRON_QUALITY};
	/// Method of Delaunay-retriangulation
	enum TriangulationType {TRIANGULATE_CIRCLE = 0, TRIANGULATE_ANGLE = 1};
	/// Method of insertion of new nodes into Delaunay triangulation
	enum InsertionType {IMPROVE_LONGEST_EDGE = 0, IMPROVE_OUTER_CIRCLE = 1, 
		IMPROVE_LONGEST_EDGE_AND_OUTER_CIRCLE = 2};
	/// Stage of Delaunay-triangulation
	enum TriangulationStage {CONSTRAIN_NONE = 0, CONSTRAIN_ACTIVE = 1, CONSTRAIN_DONE = 2};
	/// Method of triangle-to-quad conversion
	enum ConversionType {QUADS_UNKNOWN, QUADS_LEELO, QUADS_QMORPH, QUADS_SPLIT, QUADS_MIXED};
	/// Method of node-smoothing
	enum SmoothingType { SM_LAPLACE = 1, SM_DEL_SWAP = 2, SM_TOP_SWAP = 4, SM_LAPLACE_MIXED = 8, 
		SM_LAPLACE_METRIC = 16, SM_DEL_SWAP_COMPLETE = 32, SM_METRIC = 64};
	/// Method of node-smoothing
 	//enum Smoothing3Type { SM3_ALL, SM3_LAPLACE, SM3_SWAP};
	enum Smoothing3Type { SM3_LAPLACE = 1, SM3_LAPLACE_MIXED = 2, SM3_SWAP = 4, 
		SM3_SWAP_COMPLETE = 8, SM3_BORDER_PRUNE = 16, SM3_OPT_MIXED = 32, SM3_OPT_SIMPLE = 64};
	/// Method of metric calculation for triangle quality assessing
	enum MetricForTriangleQuality { QM_MIDDLE, QM_MIDDLE_AND_VERTICES_MIN, QM_MIDDLE_AND_VERTICES_AVE, 
			QM_MIDEDGES_MIN, QM_MIDEDGES_AVE, QM_VERTICES_AVE};
	/// Type of control space
	enum ControlType { CONTROL_UNIFORM, CONTROL_MESH, CONTROL_ANALYTICAL, 
		CONTROL_IDENTITY, CONTROL_QUADTREE, CONTROL_PROJECTED, 
		CONTROL_KDTREE_V, CONTROL_KDTREE_L };
	enum ControlType3d {CONTROL_UNIFORM_3D, CONTROL_MESH_3D, CONTROL_ANALYTICAL_3D, 
		CONTROL_INTERNAL_3D, CONTROL_IDENTITY_3D, 
		CONTROL_OCTREE_3D, CONTROL_FUNCTIONAL_3D, 
		CONTROL_KDTREE_3D, CONTROL_KDTREE_3D_L, CONTROL_KDTREE_3D_V, CONTROL_KDTREE_3D_LI, 
		CONTROL_KDOCTREE_3D, CONTROL_OCTREE_3D_L, CONTROL_OCTREE_3D_V, CONTROL_OCTREE_3D_LB
	};
	/// Type of control space interpolation
	enum ControlInterpolation { CONTROL_SIMPLE, CONTROL_TRIANGLE, CONTROL_VORONOI, CONTROL_MIXED,
                CONTROL_NEM_LAYERED};
	/// Type of points representative for element
	enum RepresentativePoints { RP_MIDDLE = 1, RP_VERTICES = 2, RP_MIDEDGES = 4, RP_MIDFACES = 8 };
	/// Statistical information
	struct StatData{
		StatData(double min=0.0, double max=0.0, double ave=0.0) : minimum(min), maximum(max), average(ave) {}
		double minimum, maximum, average;	// minimalna, maksymalna i œrednia wartoœæ danego parametru
	};
public:
	/// Standard constructor
	MeshData();
public:
	log4cplus::Initializer initializer;
public:
	/// "Very-small" number (scaled to the geometrical dimensions of the domain)
	double relative_small_number;
	/// "Very-large" number (scaled to the geometrical dimensions of the domain)
	double relative_infinity;
	/// Working-parameter (~model size) for scaling
	double m_model_diameter;
	/// machine epsilon
	double m_epsilon;
public:
	/// Factor for the density of view-mesh for surfaces patches
	double	m_view_polygon_density;
	/// Working-parameter used during QMorph triangle-to-quad conversion
	double m_qmorph_min_ratio;
	/// Working-parameter used during QMorph triangle-to-quad conversion
	double m_qmorph_max_ratio;
	/// Working-parameter used for determining the mode of mesh visualization (geometry, wire-mesh, full-mesh)
	int m_mesh_view_mode;
	/// Working-parameter used for showing selected blocks only
	bool m_view_selected_blocks;
public:
	/// meshlib version string
	static string version();
	/// Returns lowercase string
	static string lowercase(const string& s);
	/// sets the model diameter
	void setModelDiameter(double d);
	/// return the model diameter
	double getModelDiameter() const { return m_model_diameter; }
	/// Data for parameters
	struct ParameterData{
		/// Name of parameter
		string name;
		/// Type of property (int or double)
		bool is_int;
		/// Value of property 
		union RefValue {
			int* i_ref;
			double* d_ref;
		} value;
		/// Description of property
		string description;
	public:
		bool isInt() const		{ return is_int; }
		bool isDouble() const	{ return !is_int; }
		void setIntValue(int v)			{ *value.i_ref = v; }
		void setDoubleValue(double v)	{ *value.d_ref = v; }
		void setDoubleOrIntValue(double v) { if (is_int) setIntValue((int)v); else setDoubleValue(v); }
		int		getIntValue() const		{ return *value.i_ref; }
		double	getDoubleValue() const	{ return *value.d_ref; }
	};
	/// Set of properties
	DataHashTableKeyValue<string, std::shared_ptr<ParameterData>> m_parameters;
	/// Returns a property for a given name
	std::shared_ptr<ParameterData> getProperty(const string& name) const {
		return m_parameters.getValue(name, nullptr);
	}
	/// Sets the double-value of property at the given index
	void setPropertyDouble(const string& name, double value);
	/// Sets the int-value of property at the given index
	void setPropertyInt(const string& name, int value);
	/// Sets the description of property at the given index
	void setPropertyDescription(const string& name, string & description);
	/// Returns the double-value of property with the given name
	double getPropertyDouble(const string& name) const;
	/// Returns the description of property with the given name
	string getPropertyDescription(const string& name) const;
	/// Returns the int-value of property with the given name
	int getPropertyInt(const string& name) const;
	/// Adds a new property (name, and double-value)
	void addPropertyDouble(const string& name, double& value, double initial_value, const string& desc = "");
	/// Adds a new property (name, and int-value)
	void addPropertyInt(const string& name, int& value, int initial_value, const string& desc = "");
	/// Sets / Restores the default parameters
	void setDefaultParameters();
};

/// Global set of parameters used by whole application
extern MeshData mesh_data;

#endif // !defined(MESHDATA_H__INCLUDED)
