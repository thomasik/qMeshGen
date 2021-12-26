/////////////////////////////////////////////////////////////////////////////
// MeshBRepXML.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2009-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#if !defined(MESHBREPXML_H__INCLUDED)
#define MESHBREPXML_H__INCLUDED

#include "MeshBRep.h"

class TiXmlElement;
class DPoint2d;
class DPoint3d;
class DVector3d;

/**
 * Class responsible for parsing the description of the
 * domain, mesh, etc in XML format.
 */
class MeshBRepXML : public MeshBRep
{
public:
	/// Parse description data in brep .xml format
	bool parseFile(const string& fname);
	/// Store model/mesh/cs data into .xml format
	static bool storeFile(const string& fname, MeshContainer3d* mesh_model);
	/// Convert model/mesh/cs file to xml format
	static bool convertFile(const string& fname_from, const string& fname_to, double close_threshold = 1e-3);
	/// Process point/face data, modifies points/faces, returns number of removed points and point-faces adjacency counts in pref
	static int mergeCloseNodes(DataVector<DPoint3d> & points, DataSimpleList< DataVector<size_t> > & faces,
		double min_dist, DataVector<int> & pref);
protected:
	// Parse model brep information
	bool parseModel(TiXmlElement* modelElement);
	// Parse surface mesh information
	bool parseSurfaceMesh(TiXmlElement* smeshElement);
	// Parse surface mesh points
	int parseSurfaceMeshPoints(TiXmlElement* smeshElement);
	// Parse surface mesh faces
	int parseSurfaceMeshFaces(TiXmlElement* smeshElement);
	// Parse surface mesh local surfaces
	int parseSurfaceMeshLocalSurfaces(TiXmlElement* smeshElement);
	// Parse surface mesh local curves
	int parseSurfaceMeshLocalCurves(TiXmlElement* smeshElement);
	// Parse surface mesh border (points & edges) info
	int parseSurfaceMeshBorderInfo(TiXmlElement* smeshElement);
	// Parse model vertices
	int parseModelVertices(TiXmlElement* modelElement);
	// Parse model free points
	int parseModelFreePoints(TiXmlElement* modelElement);
	// Parse model edges
	int parseModelEdges(TiXmlElement* modelElement);
	// Parse model faces
	int parseModelFaces(TiXmlElement* modelElement);
	// Parse model blocks
	int parseModelBlocks(TiXmlElement* modelElement);
	// Parse model surfaces
	int parseModelSurfaces(TiXmlElement* modelElement);
	// Parse model curves
	int parseModelCurves(TiXmlElement* modelElement);
	// Parse sizing information
	bool parseSizing(TiXmlElement* sizingElement);
	// Parse sizing sources 2d (for faces)
	int parseSizingSources2D(TiXmlElement* sizingElement);
	// Parse sizing sources 3d (for blocks)
	int parseSizingSources3D(TiXmlElement* sizingElement);
	// Parse sizing params
	int parseSizingParams(TiXmlElement* sizingElement);
	// Parse const data
	bool parseConstants(TiXmlElement* constantsElement);
protected:
	/// Parse surface
	SurfacePtr parseSurface(TiXmlElement* element);
	/// Parse curve
	Curve2dPtr parseCurve(TiXmlElement* element);
	/// Parse curve 3d
	Curve3dPtr parseCurve3d(TiXmlElement* element);
	/// Parse plane-surface
	bool parsePlane(TiXmlElement* element, DPoint3d& pt, DVector3d& e0, DVector3d& e1);
	/// Parse point2d
	bool parsePointUV(TiXmlElement* element, DPoint2d& pt);
	/// Parse point2d
	bool parsePointUV(TiXmlElement* element, double& u, double &v);
	/// Parse point2d
	bool parsePointCT(TiXmlElement* element, int& cid, double& t);
	/// Parse point3d
	bool parsePointXYZ(TiXmlElement* element, DPoint3d& pt);
	/// Parse point3d
	bool parseVectorXYZ(TiXmlElement* element, DVector3d& vt);
	/// Parse int value
	bool parseIntValue(TiXmlElement* element, const char* label, int &id);
	/// Parse double value
	bool parseDoubleValue(TiXmlElement* element, const char* label, double &val);
	/// Parse string value
	bool parseStringValue(TiXmlElement* element, const char* label, string &val);
	/// Check attribute
	bool checkStringAttribute(TiXmlElement* element, const char* name, const char* value);
	/// Parse metric description
	bool parseMetric2d(TiXmlElement* element, Control2dItem *item);
	/// Parse metric description
	bool parseMetric3d(TiXmlElement* element, Control3dItem *item);
};

#endif // !defined(MESHBREPXML_H__INCLUDED)
