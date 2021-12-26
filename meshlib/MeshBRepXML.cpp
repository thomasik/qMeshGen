// MeshBRepXML.cpp: implementation of the MeshBRepXML class.
//
//////////////////////////////////////////////////////////////////////

#include <ctype.h>

#include <log4cplus/logger.h>
#include <log4cplus/loggingmacros.h>
using namespace log4cplus;

#include "MeshBRepXML.h"
#include "MeshData.h"
#include "DPoint.h"

#include "tinyxml.h"

#include "Curve2dCircle.h"
#include "Curve2dAnalytic.h"
#include "Curve2dBSpline.h"

#include "Curve3dParametric.h"
#include "Curve3dBSpline.h"
#include "Curve3dSegment.h"
#include "Curve3dSurfaceParametric.h"

#include "SurfacePlane.h"
#include "SurfaceAnalytic.h"
#include "SurfaceBSplinePlanar.h"
#include "SurfaceCorrected.h"
#include "SurfaceCylinder.h"
#include "SurfaceDomainHull.h"
#include "DPlanarQuadric.h"
#include "SurfacePlanarQuadric.h"

#include "MeshBlock.h"
#include "MeshDomainVolume.h"
#include "MeshContainer3d.h"

#include "MeshTriangle3d.h"
#include "MeshQuad3d.h"
#include "MeshPoly3d.h"

#include "MeshViewSet.h"
#include "MeshGenerator3dSurface.h"
#include "MeshContainer3dSurface.h"

#include "Metric3dContext.h"

#include <utility>
#include <algorithm>

/// Store model/mesh/cs data into .xml format
bool MeshBRepXML::storeFile(const string& fname, MeshContainer3d* mesh_model)
{
	ofstream ofs(fname.c_str());
	if(!ofs){
		LOG4CPLUS_ERROR(MeshLog::logger_console, 
			"Error opening file [" << fname << "] for write");
		return false;
	}
	ofs << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
		<< "<meshdoc xmlns=\"http://www.icsr.agh.edu.pl\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\">\n"
		<< "\t<header>\n"
		<< "\t\t<creator>QMeshGen</creator>\n"
		<< "\t\t<version>" << mesh_data.version() << "</version>\n"
		<< "\t</header>\n"
		<< "\t<model>\n";

	// if(store_parametric_model) -> vertices/edges/faces/surfaces

	int bct = mesh_model->getBlocksCount();
	if(bct > 0){
		ofs << "\t\t<blocks>\n";
		for(int i = 0; i < bct; i++){
			MeshDomainVolume* mdv = (MeshDomainVolume*) mesh_model->getBlockAt(i);
			assert(mdv->getType() == BLOCK_DOMAIN);
			ofs << "\t\t\t<block bid=\"" << mdv->getAreaID() << "\">\n";
			// if(store_parametric_model) -> face ids

			// if(store_mesh_3d)
			MeshContainer3d* mesh3d = mdv->getMesh();
			if(mesh3d && mesh3d->getBlocksCount() > 0){
				ofs << "\t\t\t\t<mesh>\n";

				int pct = mesh3d->getPointsCount();
				ofs << "\t\t\t\t\t<pointarray count=\"" << pct << "\">\n";
				for(int j = 0; j < pct; j++){
					const DPoint3d& pt = mesh3d->getPointAt(j)->getCoordinates();
					ofs << "\t\t\t\t\t\t" << pt.x << "\t" << pt.y << "\t" << pt.z << "\n";
				}
				ofs << "\t\t\t\t\t</pointarray>\n";

				int bct = mesh3d->getBlocksCount();
				ofs << "\t\t\t\t\t<blockarray type=\"tetrehedron\" count=\"" << bct << "\">\n";
				for(int j = 0; j < bct; j++){
					MeshBlock* block = mesh3d->getBlockAt(j);
					int bpct = block->getPointCount();
					assert(bpct == 4);
					ofs << "\t\t\t\t\t";
					for(int k = 0; k < bpct; k++)
						ofs << "\t" << block->getPoint(k)->getIndex();
					ofs << "\n";
				}
				ofs << "\t\t\t\t\t</blockarray>\n";

				ofs << "\t\t\t\t</mesh>\n";
			}

			ofs << "\t\t\t</block>\n";
		}
		ofs << "\t\t</blocks>\n";
	}

	ofs << "\t</model>\n"
		<< "</meshdoc>\n";

	return true;
}

bool MeshBRepXML::parseFile(const string& fname)
{
	TiXmlDocument doc(fname);
	bool loadOkay = doc.LoadFile();
	if(!loadOkay){
		LOG4CPLUS_ERROR(MeshLog::logger_console, 
			"Error parsing xml file: " << doc.ErrorDesc());
		return false;
	}

	TiXmlElement* meshdocElement = doc.FirstChildElement( "meshdoc" );
	if(!meshdocElement) return false;

	// -> constants
	TiXmlElement* constantsElement = meshdocElement->FirstChildElement("constants");
	if(constantsElement){
		m_ctable = new DEquationConstTable(100, "");
		if(!parseConstants(constantsElement)) return false;
	}

	// -> model (boundary reprezentation)
	// -> or surface mesh (alternative)
	TiXmlElement* modelElement = meshdocElement->FirstChildElement("model");
	TiXmlElement* surfaceMeshElement = meshdocElement->FirstChildElement("surface-mesh");

	if(modelElement){
		points_2d_only = true;
		if(!parseModel(modelElement)) return false;
	}else if(surfaceMeshElement){
		if(!parseSurfaceMesh(surfaceMeshElement)) return false;
	}else 
		return false;

	// -> sizing (control space)
	TiXmlElement* sizingElement = meshdocElement->FirstChildElement("sizing");
	if(sizingElement)
		if(!parseSizing(sizingElement)) return false;

	return true;
}

bool MeshBRepXML::parseModel(TiXmlElement* modelElement)
{
	int errors = 0;
	errors += parseModelVertices(modelElement);
	errors += parseModelEdges(modelElement);
	errors += parseModelFaces(modelElement);
	errors += parseModelBlocks(modelElement);
	errors += parseModelSurfaces(modelElement);
	errors += parseModelCurves(modelElement);
	errors += parseModelFreePoints(modelElement);

	return (errors == 0);
}

bool MeshBRepXML::parseSurfaceMesh(TiXmlElement* smeshElement)
{
	int errors = 0;
	errors += parseSurfaceMeshPoints(smeshElement);
	errors += parseSurfaceMeshFaces(smeshElement);
	errors += parseSurfaceMeshLocalSurfaces(smeshElement);
	errors += parseSurfaceMeshLocalCurves(smeshElement);
	errors += parseSurfaceMeshBorderInfo(smeshElement);

	return (errors == 0);
}

int MeshBRepXML::parseModelVertices(TiXmlElement* modelElement)
{
	int errors = 0;
	TiXmlElement* verticesElement = modelElement->FirstChildElement("vertices");
	while(verticesElement){
		TiXmlElement* vertexElement = verticesElement->FirstChildElement("vertex");
		while(vertexElement){
			int vid = -1;
			if(vertexElement->QueryIntAttribute("vid", &vid) == TIXML_SUCCESS){ // proper vid is required
				auto item = std::make_shared<PointItem>(vid);
				bool valid = true;
				if(parsePointXYZ(vertexElement, item->pt)){
					// ok, direct
				}else{  
					// via surface
					if(parseIntValue(vertexElement, "sid", item->sid)){
						points_2d_only = false;
						if(!parsePointCT(vertexElement, item->cid, item->pt.x)) // surface + contour + t
							valid = parsePointUV(vertexElement, item->pt.x, item->pt.y); // surface + (u,v)
					}else{
						// 2-point (and default surface)
						if(!parsePointCT(vertexElement, item->cid, item->pt.x)) // contour + t
							valid = parsePointUV(vertexElement, item->pt.x, item->pt.y); // (u,v)
					}
				}
				if(valid) point_list.add(item);
				else{
					LOG4CPLUS_WARN(MeshLog::logger_console, 
						"invalid vertex description, line: " << vertexElement->Row());
					++errors;
				}
			}else{
				LOG4CPLUS_WARN(MeshLog::logger_console, 
					"vertex without vid, line: " << vertexElement->Row());
				++errors;
			}
			vertexElement = vertexElement->NextSiblingElement("vertex");
		}
		verticesElement = verticesElement->NextSiblingElement("vertices");
	}
	return errors;
}

int MeshBRepXML::parseSurfaceMeshPoints(TiXmlElement* smeshElement)
{
	int errors = 0;
	TiXmlElement* smpe = smeshElement->FirstChildElement("pointarray");
	if(!smpe){
		LOG4CPLUS_WARN(MeshLog::logger_console, 
			"Missing points for surface-mesh, line: " << smeshElement->Row());
		return 1;
	}
	int pct = -1;
	if(smpe->QueryIntAttribute("count", &pct) == TIXML_SUCCESS)
		smesh_point_list.prepare(pct);
	
	const char* text = smpe->GetText();
	if(text){
		istringstream istr(text);
		DEquation eqx, eqy, eqz;
		int id = -1;
		while(true){
			string str_x, str_y, str_z;
			istr >> str_x >> str_y >> str_z;
			if(!istr) break;
			++id;
			if(eqx.parse(str_x, m_ctable) && eqy.parse(str_y, m_ctable) && eqz.parse(str_z, m_ctable)){
				DPoint3d pt(eqx.getValue(0.0), eqy.getValue(0.0), eqz.getValue(0.0));
				smesh_point_list.add( new MeshPoint3d(pt) );
			}else{
				LOG4CPLUS_WARN(MeshLog::logger_console, 
					"Error parsing coordinates for point with id: " << id);
				++errors;
			}
		}
	}else{
		LOG4CPLUS_WARN(MeshLog::logger_console, 
			"Missing points for surface-mesh, line: " << smeshElement->Row());
		return 1;
	}

	return errors;
}

int MeshBRepXML::parseSurfaceMeshFaces(TiXmlElement* smeshElement)
{
	int errors = 0;
	TiXmlElement* smfe = smeshElement->FirstChildElement("facearray");
	if(!smfe){
		LOG4CPLUS_WARN(MeshLog::logger_console, 
			"Missing faces for surface-mesh, line: " << smeshElement->Row());
		return 1;
	}
	int fct = -1;
	if(smfe->QueryIntAttribute("count", &fct) == TIXML_SUCCESS)
		smesh_face_list.prepare(fct);
	
	const char* text = smfe->GetText();
	if(text){
		istringstream istr(text);
		size_t fpct, bid0 = 0, bid1 = -1, id = -1;
		while(true){
			istr >> fpct;
			if(!istr) break;
			if(fpct < 3) {
				LOG4CPLUS_WARN(MeshLog::logger_console, 
					"Face with too low number of points: " << fpct);
				return 1;
			}
			DataVector<size_t> fpoints(fpct);
			for(size_t i = 0; i < fpct; i++){
				istr >> id;
				if(!istr) {
					LOG4CPLUS_WARN(MeshLog::logger_console, 
						"Error reading point id for a face");
					return 1;
				}else if( id < 0 || id >= smesh_point_list.countInt()){
					LOG4CPLUS_WARN(MeshLog::logger_console, 
						"Invalid point id=" << id << " for a face");
					return 1;
				}else
					fpoints.add(id);
			}
			istr >> bid0 >> bid1;
			if(!istr){
				LOG4CPLUS_WARN(MeshLog::logger_console,
					"Error reading block ids for a face");
				return 1;
			}

			MeshFace* mface = nullptr;
			if(fpct == 3){
				mface = new MeshTriangle3d(
					smesh_point_list[fpoints[0]], 
					smesh_point_list[fpoints[1]], 
					smesh_point_list[fpoints[2]] );
			}else if(fpct == 4){
				mface = new MeshQuad3d(
					smesh_point_list[fpoints[0]], 
					smesh_point_list[fpoints[1]], 
					smesh_point_list[fpoints[2]], 
					smesh_point_list[fpoints[3]] );
			}else{
				DataVector<MeshPoint3d*> spoints(fpct);
				for(size_t i = 0; i < fpct; i++)
					spoints.add( smesh_point_list[fpoints[i]] );
				mface = new MeshPoly3d(spoints);
			}
			assert(mface != nullptr);
			mface->setIntTag(TagExtended::TAG_BLOCK_0, (int)bid0);
			mface->setIntTag(TagExtended::TAG_BLOCK_1, (int)bid1);
			smesh_face_list.add(mface);
		}
	}else{
		LOG4CPLUS_WARN(MeshLog::logger_console, 
			"Missing faces for surface-mesh, line: " << smeshElement->Row());
		return 1;
	}

	return errors;
}

int MeshBRepXML::parseSurfaceMeshBorderInfo(TiXmlElement* smeshElement)
{
	int errors = 0;
	TiXmlElement* smbp = smeshElement->FirstChildElement("border-points");
	const char* text = smbp ? smbp->GetText() : nullptr;
	if(text){
		istringstream istr(text);
		int pid = -1, binfo;
		while(true){
			istr >> pid >> binfo;
			if(!istr) break;
			smesh_point_list[pid]->setBorderFlags((char)binfo);
		}
	}

	TiXmlElement* smbe = smeshElement->FirstChildElement("border-edges");
	text = smbe ? smbe->GetText() : nullptr;
	if(text){
		istringstream istr(text);
		int PId0 = -1, PId1 = -1, binfo;
		while(true){
			istr >> PId0 >> PId1 >> binfo;
			if(!istr) break;
			smesh_point_list[PId0]->getEdgeToPoint(smesh_point_list[PId1])->setBorderFlags((char)binfo);
		}
	}

	return errors;
}


int MeshBRepXML::parseModelFreePoints(TiXmlElement* modelElement)
{
	int errors = 0;
	TiXmlElement* fpse = modelElement->FirstChildElement("freepoints");
	while(fpse){
		TiXmlElement* fpe = fpse->FirstChildElement("point");
		while(fpe){
			auto item = std::make_shared<PointItem>(0);
			bool valid = true;
			if(parsePointXYZ(fpe, item->pt)){
				//parseStringValue(fpe, "label", item->label);
				fpe->QueryIntAttribute("bid", &(item->bid));
				const char* label_attr = fpe->Attribute("label");
				if(label_attr) item->label = label_attr;
				freepoint_list.add(item);
			}else{
				LOG4CPLUS_WARN(MeshLog::logger_console, 
					"invalid free-point description, line: " << fpe->Row());
				++errors;
			}
			fpe = fpe->NextSiblingElement("point");
		}
		fpe = fpse->FirstChildElement("pointarray");
		while(fpe){
			int bid = -1;
			fpe->QueryIntAttribute("bid", &bid);
			const char* label_attr = fpe->Attribute("label");
			const char* text = fpe->GetText();
			if(text){
				istringstream istr(text);
				DEquation eqx, eqy, eqz;
				while(true){
					string str_x, str_y, str_z;
					istr >> str_x >> str_y >> str_z;
					if(!istr) break;
					if(eqx.parse(str_x, m_ctable) && eqy.parse(str_y, m_ctable) && eqz.parse(str_z, m_ctable)){
						auto item = std::make_shared<PointItem>(0);
						item->bid = bid;
						if(label_attr) item->label = label_attr;
						item->pt.x = eqx.getValue(0.0);
						item->pt.y = eqy.getValue(0.0);
						item->pt.z = eqz.getValue(0.0);
						freepoint_list.add(item);
					}
				}
			}
			fpe = fpe->NextSiblingElement("pointarray");
		}
		fpse = fpse->NextSiblingElement("freepoints");
	}
	return errors;
}

int MeshBRepXML::parseModelEdges(TiXmlElement* modelElement)
{
	int errors = 0;
	// model -> edges
	TiXmlElement* edgesElement = modelElement->FirstChildElement("edges");
	while(edgesElement){
		TiXmlElement* edgeElement = edgesElement->FirstChildElement("edge");
		while(edgeElement){
			int sid = -1;
			edgeElement->QueryIntAttribute("sid", &sid); // without checking -> if not-int, then for all sids
			int vid0, vid1;
			if( edgeElement->QueryIntAttribute("vid0", &vid0) == TIXML_SUCCESS &&
				edgeElement->QueryIntAttribute("vid1", &vid1) == TIXML_SUCCESS)
			{
				auto item = std::make_shared<EdgeItem>(sid, vid0, vid1);
				// fixed ?
				item->fixed = !checkStringAttribute(edgeElement, "fixed", "false");
				// additional info
				parseIntValue(edgeElement, "cid", item->cid);
				parseDoubleValue(edgeElement, "t0", item->t0);
				parseDoubleValue(edgeElement, "t1", item->t1);
				//parseStringValue(edgeElement, "tag", item->tag);
				const char* label_attr = edgeElement->Attribute("label");
				if(label_attr) item->tag = label_attr;
				// insert
				edge_list.add(item);
			}else{
				LOG4CPLUS_WARN(MeshLog::logger_console, 
					"edge without vid0 or vid1, line: " << edgeElement->Row());
				++errors;
			}
			edgeElement = edgeElement->NextSiblingElement("edge");
		}
		edgesElement = edgesElement->NextSiblingElement("edges");
	}
	return errors;
}

int MeshBRepXML::parseModelFaces(TiXmlElement* modelElement)
{
	int errors = 0;
	// model -> faces
	TiXmlElement* facesElement = modelElement->FirstChildElement("faces");
	while(facesElement){
		// simple faces
		TiXmlElement* faceElement = facesElement->FirstChildElement("face");
		while(faceElement){
			int fid = -1, sid = 0;
			if(faceElement->QueryIntAttribute("fid", &fid) == TIXML_SUCCESS){ // proper fid is required
				faceElement->QueryIntAttribute("sid", &sid);
				auto item = std::make_shared<FaceItem>(fid, sid);
				FaceAreaItem area(fid);
				bool valid = true;
				// main outer contour of vertices
				TiXmlElement* ve = faceElement->FirstChildElement("vertex");
				while(ve){
					int vid = -1;
					if(ve->QueryIntAttribute("vid", &vid) == TIXML_SUCCESS){ // proper vid is required
						FacePoint fpt(vid);
						if( parsePointUV(ve, fpt.params) || // fpoint-> [u, v]
							parsePointCT(ve, fpt.cid, fpt.params.x)) // fpoint-> [cid, t]
						{ 
							fpt.params_available = true;
						}
						area.fpoints.add(fpt);
					}
					// next
					ve = ve->NextSiblingElement("vertex");
				}
				if(area.fpoints.countInt() < 3){
					LOG4CPLUS_WARN(MeshLog::logger_console, 
						"Face has to have at least three vertices, line: " << faceElement->Row());
					valid = false;
				}
				item->areas.add(area);
				// additional contours
				TiXmlElement* contourElement = faceElement->FirstChildElement("contour");
				while(contourElement){
					// -> attribute "sub_fid" - default is the same as fid
					int sub_fid = fid;
					contourElement->QueryIntAttribute("sub_fid", &sub_fid);
					// -> attribute "closed" - default "true"
					bool closed = !checkStringAttribute(contourElement, "closed", "false");
					// -> attribute "filled" - default "true"
					bool filled = !checkStringAttribute(contourElement, "filled", "false");

					FaceAreaItem sarea(sub_fid, closed, filled);
					// outer contour of vertices
					TiXmlElement* cve = contourElement->FirstChildElement("vertex");
					while(cve){
						int vid = -1;
						if(cve->QueryIntAttribute("vid", &vid) == TIXML_SUCCESS){ // proper vid is required
							FacePoint fpt(vid);
							if( parsePointUV(cve, fpt.params) || // fpoint-> [u, v]
								parsePointCT(cve, fpt.cid, fpt.params.x)) // fpoint-> [cid, t]
							{ 
								fpt.params_available = true;
							}
							sarea.fpoints.add(fpt);
						}
						// next
						cve = cve->NextSiblingElement("vertex");
					}
					if(sarea.closed && sarea.fpoints.countInt() < 3){
						LOG4CPLUS_WARN(MeshLog::logger_console, 
							"Closed contour has to have at least three vertices, line: " << contourElement->Row());
						valid = false;
					}else if(!sarea.closed && sarea.fpoints.countInt() < 2){
						LOG4CPLUS_WARN(MeshLog::logger_console, 
							"Opened contour has to have at least two vertices, line: " << contourElement->Row());
						valid = false;
					}else{
						item->areas.add(sarea);
					}
					contourElement = contourElement->NextSiblingElement("contour");
				}
				//parseStringValue(faceElement, "tag", item->tag);
				const char* label_attr = faceElement->Attribute("label");
				if(label_attr) item->tag = label_attr;
				// finish
				if(valid) face_list.add(item);
				else{
					LOG4CPLUS_WARN(MeshLog::logger_console, 
						"invalid face description, line: " << faceElement->Row());
					++errors;
				}
			}else{
				LOG4CPLUS_WARN(MeshLog::logger_console, 
					"face without fid, line: " << faceElement->Row());
				++errors;
			}
			faceElement = faceElement->NextSiblingElement("face");
		}
		facesElement = facesElement->NextSiblingElement("faces");
	}
	return errors;
}

int MeshBRepXML::parseModelBlocks(TiXmlElement* modelElement)
{
	int errors = 0;
	// model -> blocks
	TiXmlElement* blocksElement = modelElement->FirstChildElement("blocks");
	while(blocksElement){
		TiXmlElement* blockElement = blocksElement->FirstChildElement("block");
		while(blockElement){
			int bid = -1;
			if(blockElement->QueryIntAttribute("bid", &bid) == TIXML_SUCCESS){ // proper bid is required
				auto item = std::make_shared<BlockItem>(bid);
				// faces (topological)
				TiXmlElement* fe = blockElement->FirstChildElement("face");
				while(fe){
					int fid = -1;
					if(fe->QueryIntAttribute("fid", &fid) == TIXML_SUCCESS){ // proper fid is required
						bool oriented = !checkStringAttribute(fe, "inverted", "true");
						item->faces.add(fid);
						item->oriented.add(oriented);
					}
					fe = fe->NextSiblingElement("face");
				}
				// mesh3d ?
				TiXmlElement* me = blockElement->FirstChildElement("mesh");
				if(me){ // let's assume, it's only one mesh permitted
					TiXmlElement* parray_e = me->FirstChildElement("pointarray");
					const char* text = parray_e ? parray_e->GetText() : nullptr;
					if(text){
						istringstream istr(text);
						DEquation eqx, eqy, eqz;
						while(true){
							string str_x, str_y, str_z;
							istr >> str_x >> str_y >> str_z;
							if(!istr) break;
							if(eqx.parse(str_x, m_ctable) && eqy.parse(str_y, m_ctable) && eqz.parse(str_z, m_ctable)){
								DPoint3d pt( eqx.getValue(0.0), eqy.getValue(0.0), eqz.getValue(0.0) );
								item->mesh_vertices.append(pt);
							}
						}
					}
					TiXmlElement* barray_e = me->FirstChildElement("blockarray");
					text = barray_e ? barray_e->GetText() : nullptr;
					if(text){
						istringstream istr(text);
						int tv0, tv1, tv2, tv3;
						while(true){
							istr >> tv0 >> tv1 >> tv2 >> tv3;
							if(!istr) break;
							item->mesh_elements.append(tv0);
							item->mesh_elements.append(tv1);
							item->mesh_elements.append(tv2);
							item->mesh_elements.append(tv3);
						}
					}
				}

				block_list.add(item);
			}else{
				LOG4CPLUS_WARN(MeshLog::logger_console, 
					"block without bid, line: " << blockElement->Row());
				++errors;
			}
			blockElement = blockElement->NextSiblingElement("block");
		}
		blocksElement = blocksElement->NextSiblingElement("blocks");
	}
	return errors;
}

/// Parse surface
SurfacePtr MeshBRepXML::parseSurface(TiXmlElement* element)
{
	SurfaceDomain *domain = nullptr;
	TiXmlElement* de = element->FirstChildElement("domain-hull");
	if(de){
		DataVector<DPoint2d> hull_points;
		TiXmlElement* hpe = de->FirstChildElement("pt");
		DPoint2d pt;
		while(hpe){
			if(parsePointUV(hpe, pt.x, pt.y))
				hull_points.add(pt);
			hpe = hpe->NextSiblingElement("pt");
		}
		if(hull_points.notEmpty()){
			domain = new SurfaceDomainHull(hull_points);
		}
	}

	// -> try plane
	TiXmlElement* se = element->FirstChildElement("plane");
	if(se){
		DPoint3d plane_pt;
		DVector3d plane_e0, plane_e1;
		if(parsePlane(se, plane_pt, plane_e0, plane_e1))
		{
			auto surface = std::make_shared<SurfacePlane>(plane_pt, plane_e0, plane_e1);
			if(domain) surface->setDomain(domain);
			surface->setFixed();
			return surface;
		}else{
			LOG4CPLUS_WARN(MeshLog::logger_console, 
				"surface-plane error parsing, line: " << se->Row());
			return nullptr;
		}
	}
	// -> try analytic
	se = element->FirstChildElement("analytic");
	if(se){
		string fx, fy, fz;
		if( parseStringValue(se, "fx", fx) && parseStringValue(se, "fy", fy) && parseStringValue(se, "fz", fz))
		{
			auto surface = std::make_shared<SurfaceAnalytic>(fx.c_str(), fy.c_str(), fz.c_str(), m_ctable);
			if(domain) surface->setDomain(domain);
			surface->setFixed();
			return surface;
		}else{
			LOG4CPLUS_WARN(MeshLog::logger_console, 
				"surface-analytic error parsing, line: " << se->Row());
			return nullptr;
		}
	}
	// -> try planar quadric
	se = element->FirstChildElement("planarquadric");
	if(se){
		DPlanarQuadric pq;

		// pt0 + e0 + e1
		TiXmlElement* pte = se->FirstChildElement("pt0");
		if(!pte || !parsePointXYZ(pte, pq.plane.p0)){
			LOG4CPLUS_WARN(MeshLog::logger_console, 
				"planar-quadric error parsing, line: " << se->Row());
			return nullptr;
		}
		TiXmlElement* e0e = se->FirstChildElement("e0");
		if(!e0e || !parseVectorXYZ(e0e, pq.plane.e0)){
			LOG4CPLUS_WARN(MeshLog::logger_console, "planar-quadric error parsing, line: " << se->Row());
			return nullptr;
		}
		TiXmlElement* e1e = se->FirstChildElement("e1");
		if(!e1e || !parseVectorXYZ(e1e, pq.plane.e1)){
			LOG4CPLUS_WARN(MeshLog::logger_console, "planar-quadric error parsing, line: " << se->Row());
			return nullptr;
		}
		TiXmlElement* vqe = se->FirstChildElement("vq");
		const char* text = vqe ? vqe->GetText() : nullptr;
		if(text){
			istringstream istr(text);
			for(int i = 0; i < 6; i++){
				istr >> pq.vq[i];
			}
			if(istr){
				auto surface = std::make_shared<SurfacePlanarQuadric>(pq);
				if(domain) surface->setDomain(domain);
				surface->setFixed();
				return surface;
			}else{
				LOG4CPLUS_WARN(MeshLog::logger_console, 
					"planar-quadric error parsing, line: " << se->Row());
				return nullptr;
			}
		}else{
			LOG4CPLUS_WARN(MeshLog::logger_console, 
				"planar-quadric error parsing, line: " << se->Row());
			return nullptr;
		}
	}
	// -> try cylinder
	se = element->FirstChildElement("cylinder");
	if(se){
		DPoint3d cyl_pt;
		DVector3d cyl_vt;
		double cyl_r;

		// pt0 + axis + radius
		TiXmlElement* pte = se->FirstChildElement("pt0");
		if(!pte || !parsePointXYZ(pte, cyl_pt)){
			LOG4CPLUS_WARN(MeshLog::logger_console, 
				"surface-cylinder error parsing, line: " << se->Row());
			return nullptr;
		}
		TiXmlElement* vte = se->FirstChildElement("axis");
		if(!vte || !parseVectorXYZ(vte, cyl_vt)){
			vte = se->FirstChildElement("pt1");
			DPoint3d cyl_pt1;
			if(vte && parsePointXYZ(vte, cyl_pt1)){
				cyl_vt = cyl_pt1 - cyl_pt;
			}else{
				LOG4CPLUS_WARN(MeshLog::logger_console, 
					"surface-cylinder error parsing, line: " << se->Row());
				return nullptr;
			}
		}
		if(parseDoubleValue(se, "radius", cyl_r))
		{
			auto surface = std::make_shared<SurfaceCylinder>(cyl_pt, cyl_vt, cyl_r);
			if(domain) surface->setDomain(domain);
			surface->setFixed();
			return surface;
		}else{
			LOG4CPLUS_WARN(MeshLog::logger_console, 
				"surface-cylinder error parsing, line" << se->Row());
			return nullptr;
		}
	}
	// -> try bspline
	se = element->FirstChildElement("bspline-surface");
	if(se){ 
		// base plane
		DPoint3d plane_pt;
		DVector3d plane_e0, plane_e1;
		TiXmlElement* pe = se->FirstChildElement("plane");
		if(!pe || !parsePlane(pe, plane_pt, plane_e0, plane_e1)){
			LOG4CPLUS_WARN(MeshLog::logger_console, 
				"surface-bspline error parsing, line" << (pe ? pe->Row() : se->Row()));
			return nullptr;
		}
		// knots
		int rows, cols;
		if( !parseIntValue(se, "rows", rows) || rows <= 0 ||
			!parseIntValue(se, "columns", cols) || cols <= 0)
		{
			LOG4CPLUS_WARN(MeshLog::logger_console, 
				"surface-bspline error parsing, line: " << se->Row());
			return nullptr;
		}
		// ... list of nodes
		DataMatrix<DPoint3d> nodes(rows, cols, DPoint3d::zero);
		TiXmlElement* nodesElement = se->FirstChildElement("points");
		if(!nodesElement){
			LOG4CPLUS_WARN(MeshLog::logger_console, 
				"surface-bspline error parsing, line: " << se->Row());
			return nullptr;
		}
		TiXmlElement* ne = nodesElement->FirstChildElement("point");
		while(ne){
			int row = -1;
			int col = -1;
			DPoint3d node;
			ne->QueryIntAttribute("row", &row);
			ne->QueryIntAttribute("column", &col);
			if( row >= 0 && row < rows && col >= 0 && col < cols && parsePointXYZ(ne, node))
			{
				nodes.set(row, col, node);
			}else{
				LOG4CPLUS_WARN(MeshLog::logger_console, 
					"surface-bspline error parsing, line: " << ne->Row());
				return nullptr;
			}
			ne = ne->NextSiblingElement("point");
		}
		// ... finally, create bspline
		auto surface = std::make_shared<SurfaceBSplinePlanar>(SurfacePlane(plane_pt, plane_e0, plane_e1), &nodes);
		if(domain) surface->setDomain(domain);
		surface->setFixed();
		return surface;
	}
	// -> try corrected surface
	se = element->FirstChildElement("corrected-surface");
	if(se){
		SurfaceConstPtr base_surface = parseSurface(se);
		if(!base_surface) return nullptr;
		auto surf = std::make_shared<SurfaceCorrected>(base_surface);
		TiXmlElement* cse = se->FirstChildElement("corrections");
		while(cse){
			TiXmlElement* ce = cse->FirstChildElement("correction");
			while(ce){
				TiXmlElement* cce = ce->FirstChildElement("center");
				TiXmlElement* cte = ce->FirstChildElement("translation");
				DPoint2d ccenter;
				double cradius;
				DVector3d ctranslation;
				if(cce && cte && 
					parsePointUV(cce, ccenter) &&
					parseVectorXYZ(cte, ctranslation) &&
					parseDoubleValue(ce, "radius", cradius))
				{
					surf->insertCorrectionVector(ccenter, cradius, ctranslation);
				}else{
					LOG4CPLUS_WARN(MeshLog::logger_console, "surface-correction error parsing, line: " << ce->Row());
				}
				ce = ce->NextSiblingElement("correction");
			}
			cse = cse->NextSiblingElement("corrections");
		}
		if(domain) surf->setDomain(domain);
		surf->setFixed();
		return surf;
	}

	// unknown surface description or empty...
	return nullptr;
}

int MeshBRepXML::parseModelSurfaces(TiXmlElement* modelElement)
{
	int errors = 0;
	// model -> surfaces
	TiXmlElement* surfacesElement = modelElement->FirstChildElement("surfaces");
	while(surfacesElement){
		TiXmlElement* surfaceElement = surfacesElement->FirstChildElement("surface");
		while(surfaceElement){
			int sid = -1;
			if(surfaceElement->QueryIntAttribute("sid", &sid) == TIXML_SUCCESS){ // proper sid is required
				auto surface = parseSurface(surfaceElement);
				if(surface){
					surface->setIntTag(TagExtended::TAG_ID, sid);
					surface_list.add(surface);
				}else{
					LOG4CPLUS_WARN(MeshLog::logger_console, "error parsing surface, line: " << surfaceElement->Row());
					++errors;
				}
			}else{
				LOG4CPLUS_WARN(MeshLog::logger_console, "surface without sid, line: " << surfaceElement->Row());
				++errors;
			}
			surfaceElement = surfaceElement->NextSiblingElement("surface");
		}
		surfacesElement = surfacesElement->NextSiblingElement("surfaces");
	}
	return errors;
}

int MeshBRepXML::parseSurfaceMeshLocalSurfaces(TiXmlElement* smeshElement)
{
	int errors = 0;
	// surface mesh -> local surfaces
	TiXmlElement* surfacesElement = smeshElement->FirstChildElement("local-surfaces");
	while(surfacesElement){
		TiXmlElement* surfaceElement = surfacesElement->FirstChildElement("local-surface");
		while(surfaceElement){
			auto surface = parseSurface(surfaceElement);
			if(surface){
				TiXmlElement* lpointsElement = surfaceElement->FirstChildElement("local-points");
				const char* text = lpointsElement ? lpointsElement->GetText() : nullptr;
				if(lpointsElement && text){
					istringstream istr(text);
					int pid;
					while(true){
						istr >> pid;
						if(!istr) break;
						if(pid >= 0 && (size_t)pid < smesh_point_list.countInt()){
							// TODO - parse stored params and use instead of getParameters here
							assert(false);
							smesh_point_list[pid]->setLocalSurface( surface, 
								surface->getParameters( smesh_point_list[pid]->getCoordinates() ), 1.0 );
						}else{
							LOG4CPLUS_WARN(MeshLog::logger_console, "Local point id for a local surface outside the range");
							++errors;
						}
					}
					smesh_surface_list.add(surface);
				}else{
					LOG4CPLUS_WARN(MeshLog::logger_console, "No local point for a local surface");
					++errors;
				}
			}else{
				LOG4CPLUS_WARN(MeshLog::logger_console, "error parsing surface, line: " << surfaceElement->Row());
				++errors;
			}
			surfaceElement = surfaceElement->NextSiblingElement("local-surface");
		}
		surfacesElement = surfacesElement->NextSiblingElement("local-surfaces");
	}
	return errors;
}

/// Parse curve
Curve2dPtr MeshBRepXML::parseCurve(TiXmlElement* element)
{
	// -> try plane
	TiXmlElement* ce = element->FirstChildElement("circle");
	if(ce){
		LOG4CPLUS_ERROR(MeshLog::logger_console,   "curve-circle need to be rewritten, 0-1 to 0-2PI, auto-inner-vertices, ...");
		assert(false);

		TiXmlElement* me = ce->FirstChildElement("middle");
		DPoint2d middle;
		double radius;
		if(!me || !parsePointUV(me, middle)){
			LOG4CPLUS_WARN(MeshLog::logger_console, "curve-circle error parsing, line: " << me->Row());
			return nullptr;
		}
		if(!parseDoubleValue(ce, "radius", radius)){
			LOG4CPLUS_WARN(MeshLog::logger_console, "curve-circle error parsing, line: " << ce->Row());
			return nullptr;
		}
		auto curve = std::make_shared<Curve2dCircle>(middle, radius);
		curve->setFixed();
		return curve;
	}
	// -> try analytic
	ce = element->FirstChildElement("analytic");
	if(ce){
		string fx, fy;
		if( parseStringValue(ce, "fu", fx) && parseStringValue(ce, "fv", fy))
		{
			auto curve = std::make_shared<Curve2dAnalytic>(fx.c_str(), fy.c_str(), m_ctable);
			curve->setFixed();
			return curve;
		}else{
			LOG4CPLUS_WARN(MeshLog::logger_console, "curve-analytic error parsing, line: " << ce->Row());
			return nullptr;
		}
	}
	// -> try bspline
	ce = element->FirstChildElement("bspline-curve");
	if(ce){
		bool opened;
		DataVector<DPoint2d> nodes;
		// list of nodes
		TiXmlElement* nse = ce->FirstChildElement("points");
		if(nse){
			opened = ! checkStringAttribute(nse, "opened", "false");
			TiXmlElement* ne = nse->FirstChildElement("point");
			while(ne){
				DPoint2d node;
				if(parsePointUV(ne, node)){
					nodes.add(node);
				}else{
					LOG4CPLUS_WARN(MeshLog::logger_console, "curve-bspline error parsing, line: " << ne->Row());
					return nullptr;
				}
				ne = ne->NextSiblingElement("point");
			}
		}else{
			LOG4CPLUS_WARN(MeshLog::logger_console, "curve-bspline error parsing, line: " << ce->Row());
			return nullptr;
		}
		// ... finally, create bspline
		auto curve = std::make_shared<Curve2dBSpline>(nodes, opened);
		curve->setFixed();
		return curve;
	}

	// unknown curve description or empty...
	return nullptr;
}

int MeshBRepXML::parseModelCurves(TiXmlElement* modelElement)
{
	int errors = 0;
	// model -> curves
	TiXmlElement* curvesElement = modelElement->FirstChildElement("curves");
	while(curvesElement){
		TiXmlElement* curveElement = curvesElement->FirstChildElement("curve");
		while(curveElement){
			int cid = -1;
			if(curveElement->QueryIntAttribute("cid", &cid) == TIXML_SUCCESS){ // proper sid is required
				auto curve = parseCurve(curveElement);
				if(curve){
					curve->setIntTag(TagExtended::TAG_ID, cid);
					curve_list.add(curve);
				}else{
					LOG4CPLUS_WARN(MeshLog::logger_console, "error parsing curve, line: " << curveElement->Row());
					++errors;
				}
			}else{
				LOG4CPLUS_WARN(MeshLog::logger_console, "curve without sid or type, line: " << curveElement->Row());
				++errors;
			}
			curveElement = curveElement->NextSiblingElement("curve");
		}
		curvesElement = curvesElement->NextSiblingElement("curves");
	}
	return errors;
}

/// Parse curve 3d
Curve3dPtr MeshBRepXML::parseCurve3d(TiXmlElement* element)
{
	// -> try bspline
	TiXmlElement* ce = element->FirstChildElement("bspline3d-curve");
	if(ce){
		bool opened;
		DataVector<DPoint3d> nodes;
		// list of nodes
		TiXmlElement* nse = ce->FirstChildElement("points");
		if(nse){
			opened = ! checkStringAttribute(nse, "opened", "false");
			TiXmlElement* ne = nse->FirstChildElement("point");
			while(ne){
				DPoint3d node;
				if(parsePointXYZ(ne, node)){
					nodes.add(node);
				}else{
					LOG4CPLUS_WARN(MeshLog::logger_console, "bspline3d-curve error parsing, line: " << ne->Row());
					return nullptr;
				}
				ne = ne->NextSiblingElement("point");
			}
		}else{
			LOG4CPLUS_WARN(MeshLog::logger_console, "bspline3d-curve error parsing, line: " << ce->Row());
			return nullptr;
		}
		// ... finally, create bspline
		auto curve = std::make_shared<Curve3dBSpline>(nodes, opened);
		curve->setFixed();
		return curve;
	}

	// -> try segment
	ce = element->FirstChildElement("segment3d");
	if(ce){
		DataVector<DPoint3d> nodes;
		// list of nodes
		TiXmlElement* ept0 = ce->FirstChildElement("pt0");
		TiXmlElement* ept1 = ce->FirstChildElement("pt1");
		if(ept0 && ept1){
			DPoint3d pt0, pt1;
			if(parsePointXYZ(ept0, pt0) && parsePointXYZ(ept1, pt1)){
				// ... finally, create segment3d
				auto curve = std::make_shared<Curve3dSegment>(pt0, pt1);
				curve->setFixed();
				return curve;
			}else{
				LOG4CPLUS_WARN(MeshLog::logger_console, "segment3d error parsing, line: " << ce->Row());
				return nullptr;
			}
		}else{
			LOG4CPLUS_WARN(MeshLog::logger_console, "segment3d error parsing, line: " << ce->Row());
			return nullptr;
		}
	}

	// -> try surface + curve2d
	ce = element->FirstChildElement("curve3d-surface");
	if(ce){
		TiXmlElement* esurf = ce->FirstChildElement("surface");
		TiXmlElement* ecurv = ce->FirstChildElement("curve");

		SurfaceConstPtr surface = esurf ? parseSurface(esurf) : nullptr;
		Curve2dConstPtr curve = ecurv ? parseCurve(ecurv) : nullptr;

		if(!surface || !curve){
			LOG4CPLUS_WARN(MeshLog::logger_console, "curve3d-surface error parsing, line: " << ce->Row());
			return nullptr;
		}

		// ... finally, create segment3d
		auto curve3d = std::make_shared<Curve3dSurfaceParametric>(surface, curve);
		curve3d->setFixed();
		return curve3d;
	}

	// unknown curve3d description or empty...
	return nullptr;
}

int MeshBRepXML::parseSurfaceMeshLocalCurves(TiXmlElement* smeshElement)
{
	int errors = 0;
	// model -> curves
	TiXmlElement* curvesElement = smeshElement->FirstChildElement("local-curves");
	while(curvesElement){
		TiXmlElement* curveElement = curvesElement->FirstChildElement("local-curve");
		while(curveElement){
			auto curve = parseCurve3d(curveElement);
			if(curve){
				TiXmlElement* lpointsElement = curveElement->FirstChildElement("local-points");
				const char* text = lpointsElement ? lpointsElement->GetText() : nullptr;
				if(text){
					istringstream istr(text);
					int pid;
					double last_t = 0.0;
					while(true){
						istr >> pid;
						if(!istr) break;
						if(pid >= 0 && pid < smesh_point_list.countInt()){
							assert(false); // parameter for curve should be given as well ...
							smesh_point_list[pid]->setLocalCurve( curve, 
								last_t = curve->getParameter( smesh_point_list[pid]->getCoordinates(), last_t ));
						}else{
							LOG4CPLUS_WARN(MeshLog::logger_console, "Local point id for a local curve outside the range");
							++errors;
						}
					}
					smesh_curve_list.add(curve);
				}else{
					LOG4CPLUS_WARN(MeshLog::logger_console, "No local point for a local curve");
					++errors;
				}
			}else{
				LOG4CPLUS_WARN(MeshLog::logger_console, "error parsing curve, line: " << curveElement->Row());
				++errors;
			}
			curveElement = curveElement->NextSiblingElement("local-curve");
		}
		curvesElement = curvesElement->NextSiblingElement("local-curves");
	}
	return errors;
}

int MeshBRepXML::parseSizingParams(TiXmlElement* sizingElement)
{
	int errors = 0;
	// sizing -> params
	TiXmlElement* paramsElement = sizingElement->FirstChildElement("params");
	while(paramsElement){
		// 2D -> for faces
		TiXmlElement* paramElement = paramsElement->FirstChildElement("param");
		while(paramElement){
			string pname, pstrvalue;
			double pvalue;
			if( parseStringValue(paramElement, "name", pname) &&
				parseStringValue(paramElement, "value", pstrvalue) &&
				DEquation::stringToDouble(pstrvalue, DEquation::v_auto, pvalue))
			{
				auto p = mesh_data.getProperty(pname);
				if(p){
					p->setDoubleOrIntValue(pvalue);
					LOG4CPLUS_INFO(MeshLog::logger_console, pname << " -> " << pvalue);
				}else{
					LOG4CPLUS_WARN(MeshLog::logger_console, "sizing-param unknown name, line: " << paramElement->Row());
					++errors;
				}
			}else{
				LOG4CPLUS_WARN(MeshLog::logger_console, "sizing-param error parsing, line: " << paramElement->Row());
				++errors;
			}
			paramElement = paramElement->NextSiblingElement("param");
		}
		paramsElement = paramsElement->NextSiblingElement("params");
	}
	return errors;
}

int MeshBRepXML::parseSizingSources2D(TiXmlElement* sizingElement)
{
	int errors = 0;
	// sizing -> control sources
	TiXmlElement* sourcesElement = sizingElement->FirstChildElement("sources");
	while(sourcesElement){
		// 2D -> for faces
		TiXmlElement* sourceElement = sourcesElement->FirstChildElement("source2d");
		while(sourceElement){
			int fid = -1;
			int multi_nr = 0;
			const char* stype = sourceElement->Attribute("type");
			sourceElement->QueryIntAttribute("multi", &multi_nr);
			if(sourceElement->QueryIntAttribute("fid", &fid) == TIXML_SUCCESS && stype){ // proper fid is required
				bool valid = true;
				const string strtype(stype);
				if(strtype == "analytic"){ 
					auto item = std::make_shared<Control2dItem>(fid, Control2dItem::CS_DOMAIN);
					// metric + 2 x domain
					TiXmlElement* domain1_e = sourceElement->FirstChildElement("domain");
					TiXmlElement* domain2_e = domain1_e ? domain1_e->NextSiblingElement("domain") : nullptr;

					if(parseMetric2d(sourceElement, item.get()) &&
						item->parseDomain(	domain1_e ? domain1_e->GetText() : nullptr, 
											domain2_e ? domain2_e->GetText() : nullptr,
											m_ctable))
					{
						item->multi_nr = multi_nr; // 0 -> normal, 1 -> multi-domain
						control2d_list.add(item);
					}else{
						valid = false;
						LOG4CPLUS_WARN(MeshLog::logger_console, "source-analytic error parsing, line: " << sourceElement->Row());
						++errors;
					}
				}else if(strtype == "point"){ 
					auto item = std::make_shared<Control2dItem>(fid, Control2dItem::CS_POINT);
					// point[u,v] + metric[lx,ly,angle] + radius
					TiXmlElement* pt_e = sourceElement->FirstChildElement("point");
					if(!pt_e || !parsePointUV(pt_e, item->pt[0])){
						valid = false;
						LOG4CPLUS_WARN(MeshLog::logger_console, "source-point error parsing, line: " << pt_e->Row());
						++errors;
					}
					if(!parseMetric2d(sourceElement, item.get())){
						valid = false;
						LOG4CPLUS_WARN(MeshLog::logger_console, "source-point error parsing, line: " << sourceElement->Row());
						++errors;
					}
					if(!parseDoubleValue(sourceElement, "radius", item->radius)){
						valid = false;
						LOG4CPLUS_WARN(MeshLog::logger_console, "source-point error parsing, line: " << sourceElement->Row());
						++errors;
					}

					if(valid) control2d_list.add(item);
				}else if(strtype == "segment"){ 
					auto item = std::make_shared<Control2dItem>(fid, Control2dItem::CS_SEGMENT);
					item->directional = checkStringAttribute(sourceElement, "directional", "true");
					// point[u,v] + point[u,v] + metric[lx,ly,angle] + radius
					TiXmlElement* pt0_e = sourceElement->FirstChildElement("point0");
					if(!pt0_e || !parsePointUV(pt0_e, item->pt[0])){
						valid = false;
						LOG4CPLUS_WARN(MeshLog::logger_console, "source-point error parsing, line: " << pt0_e->Row());
						++errors;
					}
					TiXmlElement* pt1_e = sourceElement->FirstChildElement("point1");
					if(!pt1_e || !parsePointUV(pt1_e, item->pt[1])){
						valid = false;
						LOG4CPLUS_WARN(MeshLog::logger_console, "source-point error parsing, line: " << pt1_e->Row());
						++errors;
					}
					if(!parseMetric2d(sourceElement, item.get())){
						valid = false;
						LOG4CPLUS_WARN(MeshLog::logger_console, "source-point error parsing, line: " << sourceElement->Row());
						++errors;
					}
					if(!parseDoubleValue(sourceElement, "radius", item->radius)){
						valid = false;
						LOG4CPLUS_WARN(MeshLog::logger_console, "source-point error parsing, line: " << sourceElement->Row());
						++errors;
					}

					if(valid) control2d_list.add(item);
				}else if(strtype == "curve"){ 
					auto item = std::make_shared<Control2dItem>(fid, Control2dItem::CS_CURVE);
					// cid + t0 + t1 + metric[lx,ly,angle] + radius
					if( parseIntValue(sourceElement, "cid", item->curve_id) &&
						parseDoubleValue(sourceElement, "t0", item->pt[0].x) &&
						parseDoubleValue(sourceElement, "t1", item->pt[0].y) &&
						parseMetric2d(sourceElement, item.get()) &&
						parseDoubleValue(sourceElement, "radius", item->radius) )
					{
						control2d_list.add(item);
					}else{
						valid = false;
						LOG4CPLUS_WARN(MeshLog::logger_console, "source-point error parsing, line: " << sourceElement->Row());
						++errors;
					}
				}else{
					LOG4CPLUS_WARN(MeshLog::logger_console, "source with invalid type, line: " << sourceElement->Row());
					++errors;
				}
			}else{
				LOG4CPLUS_WARN(MeshLog::logger_console, "source without sid or type, line: " << sourceElement->Row());
				++errors;
			}
			sourceElement = sourceElement->NextSiblingElement("source2d");
		}
		sourcesElement = sourcesElement->NextSiblingElement("sources");
	}
	return errors;
}

int MeshBRepXML::parseSizingSources3D(TiXmlElement* sizingElement)
{
	int errors = 0;
	// sizing -> control sources
	TiXmlElement* sourcesElement = sizingElement->FirstChildElement("sources");
	while(sourcesElement){
		// 3D -> for blocks
		TiXmlElement* sourceElement = sourcesElement->FirstChildElement("source3d");
		while(sourceElement){
			int bid = -1;
			int multi_nr = 0;
			const char* stype = sourceElement->Attribute("type");
			sourceElement->QueryIntAttribute("multi", &multi_nr);
			if(sourceElement->QueryIntAttribute("bid", &bid) == TIXML_SUCCESS && stype){ // proper fid is required
				bool valid = true;
				const string strtype(stype);
				if(strtype == "internal"){
					control3d_list.add(std::make_shared<Control3dItem>(bid, Control3dItem::CS_INTERNAL));
				}else if(strtype == "analytic"){ 
					auto item = std::make_shared<Control3dItem>(bid, Control3dItem::CS_DOMAIN);
					// metric + 2 x domain
					TiXmlElement* domain1_e = sourceElement->FirstChildElement("domain");
					TiXmlElement* domain2_e = domain1_e ? domain1_e->NextSiblingElement("domain") : nullptr;

					TiXmlElement* box_e = sourceElement->FirstChildElement("bbox");
					if(box_e &&
						parseDoubleValue(box_e, "x0", item->bbox.x0) &&
						parseDoubleValue(box_e, "x1", item->bbox.x1) &&
						parseDoubleValue(box_e, "y0", item->bbox.y0) &&
						parseDoubleValue(box_e, "y1", item->bbox.y1) &&
						parseDoubleValue(box_e, "z0", item->bbox.z0) &&
						parseDoubleValue(box_e, "z1", item->bbox.z1))
					{
						item->bbox.valid = true;
					}

					if(parseMetric3d(sourceElement, item.get()) &&
						item->parseDomain(	domain1_e ? domain1_e->GetText() : nullptr, 
											domain2_e ? domain2_e->GetText() : nullptr,
											m_ctable))
					{
						item->multi_nr = multi_nr; // 0 -> normal, 1 -> multi-domain
						control3d_list.add(item);
					}else{
						valid = false;
						LOG4CPLUS_WARN(MeshLog::logger_console, "source-analytic error parsing, line: " << sourceElement->Row());
						++errors;
					}
				}else if(strtype == "point"){ 
					auto item = std::make_shared<Control3dItem>(bid, Control3dItem::CS_POINT);
					// point[x,y,z] or pointarray + metric + radius
					TiXmlElement* pt_e = sourceElement->FirstChildElement("point");
					DPoint3d pt;
					while(pt_e){
						if(parsePointXYZ(pt_e, pt)){
							item->pts.add(pt);
						}else{
							LOG4CPLUS_WARN(MeshLog::logger_console, "source-point error parsing, line: " << pt_e->Row());
							++errors;
						}
						pt_e = pt_e->NextSiblingElement("point");
					}
					pt_e = sourceElement->FirstChildElement("pointarray");
					const char* text = pt_e ? pt_e->GetText() : nullptr;
					while(text){
						istringstream istr(text);
						DEquation eqx, eqy, eqz;
						while(true){
							string str_x, str_y, str_z;
							istr >> str_x >> str_y >> str_z;
							if(!istr) break;
							if(eqx.parse(str_x) && eqy.parse(str_y) && eqz.parse(str_z)){
								pt.x = eqx.getValue(0.0);
								pt.y = eqy.getValue(0.0);
								pt.z = eqz.getValue(0.0);
								item->pts.add(pt);
							}
						}
						pt_e = pt_e->NextSiblingElement("pointarray");
					}
					if(!parseMetric3d(sourceElement, item.get())){
						valid = false;
						LOG4CPLUS_WARN(MeshLog::logger_console, "source-point error parsing, line: " << sourceElement->Row());
						++errors;
					}
					if(!parseDoubleValue(sourceElement, "radius", item->radius)){
						valid = false;
						LOG4CPLUS_WARN(MeshLog::logger_console, "source-point error parsing, line: " << sourceElement->Row());
						++errors;
					}
					if(item->pts.empty()){
						valid = false;
						LOG4CPLUS_WARN(MeshLog::logger_console, "at least one point required for source, line: " << sourceElement->Row());
					}

					if(valid) control3d_list.add(item);
				}else if(strtype == "segment"){ 
					auto item = std::make_shared<Control3dItem>(bid, Control3dItem::CS_SEGMENT);
					item->directional = checkStringAttribute(sourceElement, "directional", "true");
					// point[x,y,z] + point[x,y,z] + metric + radius
					TiXmlElement* pt0_e = sourceElement->FirstChildElement("point0");
					DPoint3d pt;
					if(!pt0_e || !parsePointXYZ(pt0_e, pt)){
						valid = false;
						LOG4CPLUS_WARN(MeshLog::logger_console, "source-point error parsing, line: " << pt0_e->Row());
						++errors;
					}else item->pts.add(pt);
					TiXmlElement* pt1_e = sourceElement->FirstChildElement("point1");
					if(!pt1_e || !parsePointXYZ(pt1_e, pt)){
						valid = false;
						LOG4CPLUS_WARN(MeshLog::logger_console, "source-point error parsing, line: " << pt1_e->Row());
						++errors;
					}else item->pts.add(pt);
					if(!parseMetric3d(sourceElement, item.get())){
						valid = false;
						LOG4CPLUS_WARN(MeshLog::logger_console, "source-point error parsing, line: " << sourceElement->Row());
						++errors;
					}
					if(!parseDoubleValue(sourceElement, "radius", item->radius)){
						valid = false;
						LOG4CPLUS_WARN(MeshLog::logger_console, "source-point error parsing, line: " << sourceElement->Row());
						++errors;
					}

					if(valid) control3d_list.add(item);
				}else if(strtype == "curve"){ 
					auto item = std::make_shared<Control3dItem>(bid, Control3dItem::CS_CURVE);
					// cid + t0 + t1 + metric + radius
					DPoint3d pt;
					if( parseIntValue(sourceElement, "cid", item->curve_id) &&
						parseDoubleValue(sourceElement, "t0", pt.x) &&
						parseDoubleValue(sourceElement, "t1", pt.y) &&
						parseMetric3d(sourceElement, item.get()) &&
						parseDoubleValue(sourceElement, "radius", item->radius) )
					{
						item->pts.add(pt);
						control3d_list.add(item);
					}else{
						valid = false;
						LOG4CPLUS_WARN(MeshLog::logger_console, "source-point error parsing, line: " << sourceElement->Row());
						++errors;
					}
				}else if(strtype == "triangle"){ 
					auto item = std::make_shared<Control3dItem>(bid, Control3dItem::CS_TRIANGLE);
					item->directional = checkStringAttribute(sourceElement, "directional", "true");
					// point[x,y,z] + point[x,y,z] + point[x,y,z] + metric + radius
					DPoint3d pt;
					TiXmlElement* pt0_e = sourceElement->FirstChildElement("point0");
					if(!pt0_e || !parsePointXYZ(pt0_e, pt)){
						valid = false;
						LOG4CPLUS_WARN(MeshLog::logger_console, "source-point error parsing, line: " << pt0_e->Row());
						++errors;
					}else item->pts.add(pt);
					TiXmlElement* pt1_e = sourceElement->FirstChildElement("point1");
					if(!pt1_e || !parsePointXYZ(pt1_e, pt)){
						valid = false;
						LOG4CPLUS_WARN(MeshLog::logger_console, "source-point error parsing, line: " << pt1_e->Row());
						++errors;
					}else item->pts.add(pt);
					TiXmlElement* pt2_e = sourceElement->FirstChildElement("point2");
					if(!pt2_e || !parsePointXYZ(pt2_e, pt)){
						valid = false;
						LOG4CPLUS_WARN(MeshLog::logger_console, "source-point error parsing, line: " << pt2_e->Row());
						++errors;
					}else item->pts.add(pt);
					if(!parseMetric3d(sourceElement, item.get())){
						valid = false;
						LOG4CPLUS_WARN(MeshLog::logger_console, "source-point error parsing, line: " << sourceElement->Row());
						++errors;
					}
					if(!parseDoubleValue(sourceElement, "radius", item->radius)){
						valid = false;
						LOG4CPLUS_WARN(MeshLog::logger_console, "source-point error parsing, line: " << sourceElement->Row());
						++errors;
					}

					if(valid) control3d_list.add(item);
				}else{
					LOG4CPLUS_WARN(MeshLog::logger_console, "source with invalid type, line: " << sourceElement->Row());
					++errors;
				}
			}else{
				LOG4CPLUS_WARN(MeshLog::logger_console, "source without sid or type, line: " << sourceElement->Row());
				++errors;
			}
			sourceElement = sourceElement->NextSiblingElement("source3d");
		}
		sourcesElement = sourcesElement->NextSiblingElement("sources");
	}
	return errors;
}

bool MeshBRepXML::parseSizing(TiXmlElement* sizingElement)
{
	int errors = 0;
	errors += parseSizingParams(sizingElement);
	errors += parseSizingSources2D(sizingElement);
	errors += parseSizingSources3D(sizingElement);

	return (errors == 0);
}

bool MeshBRepXML::parseConstants(TiXmlElement* constantsElement)
{
	int errors = 0;
	// sizing -> control sources
	TiXmlElement* constElement = constantsElement->FirstChildElement("const");
	while(constElement){
		TiXmlElement* cname_e = constElement->FirstChildElement("name");
		TiXmlElement* cvalue_e = constElement->FirstChildElement("value");
		if(!cname_e || !cvalue_e){
			LOG4CPLUS_WARN(MeshLog::logger_console, "missing name or value for constant, line: " << constElement->Row());
			++errors;
		}else{
			const char* text = cname_e->GetText();
			string cname = text ? text : "";
			std::transform(cname.begin(), cname.end(), cname.begin(), 
				[](unsigned char c) { return toupper(c); });
			//size_t clen = cname.length();
			//for(size_t i = 0; i < clen; i++)
			//	cname[i] = toupper(cname[i]);
			if(m_ctable->contains(cname)){
				LOG4CPLUS_WARN(MeshLog::logger_console, "duplicate name for constant, line: " << constElement->Row());
				++errors;
			}else{
				DEquation eq;
				//-------
				//DataVector<string> keys;
				//m_ctable->getKeys(keys);
				//LOG4CPLUS_INFO(MeshLog::logger_mesh, "Parsing value for NAME=" << cname << ", string value = " << cvalue_e->GetText());
				//LOG4CPLUS_INFO(MeshLog::logger_mesh, "* Current hashtable: ");
				//for(int i = 0; i < keys.countInt(); i++){
				//	LOG4CPLUS_INFO(MeshLog::logger_mesh, " NAME=" << keys[i] << " VALUE=" << m_ctable->getValue(keys[i], -1.0));
				//}
				//-------
				if(!eq.parse(cvalue_e->GetText(), m_ctable)){
					LOG4CPLUS_WARN(MeshLog::logger_console, "error parsing value for constant, line: " << constElement->Row());
					++errors;
				}else{
					// ok, add new constant
					double val = eq.getValue(0.0);
//					LOG4CPLUS_INFO(MeshLog::logger_mesh, "* adding value: " << val);
					m_ctable->insert(cname, val);
//					LOG4CPLUS_INFO(MeshLog::logger_mesh, "* contains: " << m_ctable->contains(cname));
//					LOG4CPLUS_INFO(MeshLog::logger_mesh, "* getValue: " << m_ctable->getValue(cname, -2.0));
				}
			}
		}
		// next element
		constElement = constElement->NextSiblingElement("const");
	}
	return (errors == 0);
}

/// Parse point2d
bool MeshBRepXML::parsePointUV(TiXmlElement* element, double& u, double &v)
{
	if(!element) return false;

	TiXmlElement* ue = element->FirstChildElement("u");
	TiXmlElement* ve = element->FirstChildElement("v");
	if(!ue || !ve){ // try "x"+"y"
		ue = element->FirstChildElement("x");
		ve = element->FirstChildElement("y");
	}
	if(!ue || !ve) return false;

	DEquation eq;
	if(eq.parse(ue->GetText(), m_ctable)) 
		u = eq.getValue(0.0);
	else return false;
	if(eq.parse(ve->GetText(), m_ctable)) 
		v = eq.getValue(0.0);
	else return false;

	return true;
}

/// Parse point2d
bool MeshBRepXML::parsePointUV(TiXmlElement* element, DPoint2d& pt)
{
	return parsePointUV(element, pt.x, pt.y);
}

/// Parse point2d
bool MeshBRepXML::parsePointCT(TiXmlElement* element, int &cid, double& t)
{
	if(!element) return false;

	TiXmlElement* ce = element->FirstChildElement("cid");
	TiXmlElement* te = element->FirstChildElement("t");
	if(!ce || !te) return false;

	DEquation eq;
	if(eq.parse(ce->GetText(), m_ctable)) 
		cid = (int)eq.getValue(0.0);
	else return false;
	if(eq.parse(te->GetText(), m_ctable)) 
		t = eq.getValue(0.0);
	else return false;

	return true;
}

/// Parse point3d
bool MeshBRepXML::parsePointXYZ(TiXmlElement* element, DPoint3d& pt)
{
	if(!element) return false;

	TiXmlElement* xe = element->FirstChildElement("x");
	TiXmlElement* ye = element->FirstChildElement("y");
	TiXmlElement* ze = element->FirstChildElement("z");
	if(!xe || !ye || !ze) return false;

	DEquation eq;
	if(eq.parse(xe->GetText(), m_ctable)) 
		pt.x = eq.getValue(0.0);
	else return false;
	if(eq.parse(ye->GetText(), m_ctable)) 
		pt.y = eq.getValue(0.0);
	else return false;
	if(eq.parse(ze->GetText(), m_ctable))
		pt.z = eq.getValue(0.0);
	else return false;

	return true;
}

/// Parse point3d
bool MeshBRepXML::parseVectorXYZ(TiXmlElement* element, DVector3d& vt)
{
	if(!element) return false;

	TiXmlElement* xe = element->FirstChildElement("x");
	TiXmlElement* ye = element->FirstChildElement("y");
	TiXmlElement* ze = element->FirstChildElement("z");
	if(!xe || !ye || !ze) return false;

	DEquation eq;
	if(eq.parse(xe->GetText(), m_ctable)) 
		vt.x = eq.getValue(0.0);
	else return false;
	if(eq.parse(ye->GetText(), m_ctable)) 
		vt.y = eq.getValue(0.0);
	else return false;
	if(eq.parse(ze->GetText(), m_ctable))
		vt.z = eq.getValue(0.0);
	else return false;

	return true;
}

/// Parse int id
bool MeshBRepXML::parseIntValue(TiXmlElement* element, const char* label, int &id)
{
	if(!element) return false;

	TiXmlElement* ide = element->FirstChildElement(label);
	if(!ide) return false;

	DEquation eq;
	if(eq.parse(ide->GetText(), m_ctable)) 
		id = (int)eq.getValue(0.0);
	else return false;

	return true;
}

/// Parse double value
bool MeshBRepXML::parseDoubleValue(TiXmlElement* element, const char* label, double &val)
{
	if(!element) return false;

	TiXmlElement* ve = element->FirstChildElement(label);
	if(!ve) return false;

	DEquation eq;
	if(eq.parse(ve->GetText(), m_ctable)) 
		val = eq.getValue(0.0);
	else return false;

	return true;
}

/// Parse string value
bool MeshBRepXML::parseStringValue(TiXmlElement* element, const char* label, string &val)
{
	if(!element) return false;

	TiXmlElement* ve = element->FirstChildElement(label);
	if(!ve) return false;

	const char* val_cstr = ve->GetText();
	if(val_cstr){
		val = val_cstr;
		return true;
	}else 
		return false;
}

bool MeshBRepXML::checkStringAttribute(TiXmlElement* element, const char* name, const char* value)
{
	const char* avalue = element->Attribute(name);
	if(!avalue) return false;
	return strcmp(avalue, value) == 0;
}

bool MeshBRepXML::parseMetric2d(TiXmlElement* element, Control2dItem *item)
{
	TiXmlElement* me = element->FirstChildElement("metric");
	if(!me) return false;

	TiXmlElement* lxye = me->FirstChildElement("lxy");
	if(lxye){
		const char* len = lxye->GetText();
		return item->parseCDS(len, len, "0", m_ctable);
	}

	TiXmlElement* lxe = me->FirstChildElement("lx");
	TiXmlElement* lye = me->FirstChildElement("ly");
	TiXmlElement* ane = me->FirstChildElement("angle");

	if(!lxe || !lye || !ane) return false;

	return item->parseCDS(lxe->GetText(), lye->GetText(), ane->GetText(), m_ctable);
}

bool MeshBRepXML::parseMetric3d(TiXmlElement* element, Control3dItem *item)
{
	TiXmlElement* me = element->FirstChildElement("metric");
	if(!me) return false;

	TiXmlElement* lxyze = me->FirstChildElement("lxyz");
	if(lxyze){
		const char* len = lxyze->GetText();
		return item->parseCDS(len, len, len, "0", "0", "0", m_ctable);
	}

	TiXmlElement* lxe = me->FirstChildElement("lx");
	TiXmlElement* lye = me->FirstChildElement("ly");
	TiXmlElement* lze = me->FirstChildElement("lz");
	TiXmlElement* axe = me->FirstChildElement("ax");
	TiXmlElement* aye = me->FirstChildElement("ay");
	TiXmlElement* aze = me->FirstChildElement("az");

	if(!lxe || !lye || !lze || !axe || !aye || !aze) return false;

	return item->parseCDS(
		lxe->GetText(), lye->GetText(), lze->GetText(),
		axe->GetText(), aye->GetText(), aze->GetText(),
		m_ctable);
}

/// Parse plane-surface
bool MeshBRepXML::parsePlane(TiXmlElement* element, DPoint3d& pt, DVector3d& e0, DVector3d& e1)
{
	// pt0 + e0 + e1
	TiXmlElement* pt0e = element->FirstChildElement("pt0");
	if(!pt0e || !parsePointXYZ(pt0e, pt)) return false;

	TiXmlElement* e0e = element->FirstChildElement("e0");
	if(!e0e || !parseVectorXYZ(e0e, e0)) return false;

	TiXmlElement* e1e = element->FirstChildElement("e1");
	if(!e1e || !parseVectorXYZ(e1e, e1)) return false;

	return true;
}

/// Convert model/mesh/cs file to xml format
bool MeshBRepXML::convertFile(const string& fname_from, const string& fname_to,	double close_threshold)
{
	LOG4CPLUS_INFO(MeshLog::logger_mesh, 
		"MeshBRepXML::convertFile(from->" << fname_from << ", to-> " << fname_to 
		<< ", threshold->" << close_threshold << ")");
	size_t dot_index = fname_from.rfind('.');
	if(dot_index == string::npos){
		LOG4CPLUS_ERROR(MeshLog::logger_console, "Input file [" << fname_from << "] should have an extension");
		return false;
	}
	const string extension = fname_from.substr(dot_index);
	if(extension == ".bv"){
		ifstream fin(fname_from.c_str());
		if(!fin){
			LOG4CPLUS_ERROR(MeshLog::logger_console, "Error opening input file [" << fname_from << "]");
			return false;
		}
		ofstream fout(fname_to.c_str());
		if(!fout){
			LOG4CPLUS_ERROR(MeshLog::logger_console, "Error opening output file [" << fname_to << "]");
			return false;
		}

		fout << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << endl;
		fout << "<meshdoc xmlns=\"http://www.icsr.agh.edu.pl\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\">" << endl;
		fout << "\t" << "<header>" << endl;
		fout << "\t\t" << "<creator>Tomasz Jurczyk</creator>" << endl;
		fout << "\t\t" << "<version>1.0</version>" << endl;
		fout << "\t\t" << "<description>Converted from " << fname_from << "</description>" << endl;
		fout << "\t" << "</header>" << endl;
		fout << "\t" << "<model>" << endl;
		fout << "\t\t" << "<vertices>" << endl;

		string line;
		getline(fin, line);	// group 0 SurfLab_2009
		getline(fin, line);	// 1
		int pct, fct;
		fin >> pct >> fct; // 450 494
		for(int i = 0; i < pct; i++){
			DPoint3d pt;
			fin >> pt.x >> pt.y >> pt.z;	// -1.0e-06 0.003676 1.448846
			fout << "\t\t\t" << "<vertex vid=\"" << i << "\">	<x>" << pt.x << "</x> <y>" 
				<< pt.y << "</y> <z>" << pt.z << "</z> </vertex>" << endl;
		}
		fout << "\t\t" << "</vertices>" << endl;
		fout << "\t\t" << "<faces>" << endl;
		for(int i = 0; i < fct; i++){
			int ect, pid;
			fin >> ect;
			fout << "\t\t\t" << "<face fid=\"" << i << "\">" << endl; 
			fout << "\t\t\t\t";
			for(int j = 0; j < ect; j++){
				fin >> pid;
				if((j > 0) && (j%4 == 0)) fout << endl << "\t\t\t\t";
				fout << "<vertex vid=\"" << pid << "\"/> ";
			}
			fout << endl << "\t\t\t" << "</face>" << endl; 
		}
		fout << "\t\t" << "</faces>" << endl;
		fout << "\t\t" << "<blocks>" << endl;
		fout << "\t\t\t" << "<block bid=\"1\">" << endl;
		for(int i = 0; i < fct; i++)
			fout << "\t\t\t\t" << "<face fid=\"" << i << "\"/>" << endl;
		fout << "\t\t\t" << "</block>" << endl;
		fout << "\t\t" << "</blocks>" << endl;
		fout << "\t" << "</model>" << endl;
		fout << "</meshdoc>" << endl;
		if(!fin){
			LOG4CPLUS_ERROR(MeshLog::logger_console, "Error reading input file: " << fname_from);
			return false;
		}
		if(!fout){
			LOG4CPLUS_ERROR(MeshLog::logger_console, "Error writing output file: " << fname_to);
			return false;
		}
		return true;
	}else if(extension == ".off"){
		ifstream fin(fname_from.c_str());
		if(!fin){
			LOG4CPLUS_ERROR(MeshLog::logger_console, "Error opening input file: " << fname_from);
			return false;
		}

		string line;
		getline(fin, line);	// OFF
		if(line != "OFF"){
			LOG4CPLUS_ERROR(MeshLog::logger_console, "Wrong header for an OFF file: " << line);
			return false;
		}
		int pct, fct, bct;
		fin >> pct >> fct >> bct; // 7840 5238 0
		if(pct < 3 || fct < 1){
			LOG4CPLUS_ERROR(MeshLog::logger_console,   "Required at least three points and one face");
			return false;
		}
		if(bct > 0){
			LOG4CPLUS_ERROR(MeshLog::logger_console,   "Non-zero number of blocks? Not implemented.");
			return false;
		}

		DataVector<DPoint3d> points(pct);
		DataSimpleList< DataVector<int> > faces;

		// read data
		for(int i = 0; i < pct; i++){
			DPoint3d pt;
			fin >> pt.x >> pt.y >> pt.z;	// -1.0e-06 0.003676 1.448846
			if(!fin) { 
				LOG4CPLUS_ERROR(MeshLog::logger_console, "Error reading point with id=" << i);
				return false;
			}
			points.add(pt);
		}
		for(int i = 0; i < fct; i++){
			int ect, pid;
			fin >> ect;
			if(!fin) { 
				LOG4CPLUS_ERROR(MeshLog::logger_console, "Error reading face with id=" << i);
				return false;
			}
			DataVector<int> face(ect);
			for(int j = 0; j < ect; j++){
				fin >> pid;
				if(!fin) { 
					LOG4CPLUS_ERROR(MeshLog::logger_console, "Error reading point-ids for face with id=" << i);
					return false;
				}
				face.add(pid);
			}
			faces.append(face);
		}

		if(!fin){
			LOG4CPLUS_ERROR(MeshLog::logger_console, "Error reading input file: " << fname_from);
			return false;
		}

		//DSphere::test();
		//DCylinder::test();

		MeshContainer3dSurface* surface_mesh = MeshGenerator3dSurface::createMeshFromFaces(points, faces);
		MeshGenerator3dSurface::checkAndFixInvertedFacesTopological(surface_mesh);
		// surface_mesh -> poly-mesh

		DBox bbox;
		for(int i = 0; i < pct; i++) bbox.addPoint(points[i]);
		double diameter = bbox.getDiameter();
		double min_dist = diameter * close_threshold;

//		LOG4CPLUS_DEBUG(MeshLog::logger_console,"Setting boundary by topology ...");
		surface_mesh->clearBoundaryFlags();
		surface_mesh->setBoundaryFeatureEdges();

		if(true){
			SHOW_MESH("OFF: create + top.boundary", surface_mesh->getViewSet());
		}

		MeshGenerator3dSurface::calculateNormals(surface_mesh);

//		LOG4CPLUS_DEBUG(MeshLog::logger_console,"Mesh validation & removing tiny edges ...");
		//int vct = 
		MeshGenerator3dSurface::validateSurfaceMesh(surface_mesh, min_dist);
		//if(vct > 0){
		//	string v_fname_from = fname_from;
		//	v_fname_from.insert(v_fname_from.length()-4,"-v");
		//	surface_mesh->storeOFF(v_fname_from);
		//}

		if( !surface_mesh->isValid() ) LOG4CPLUS_WARN(MeshLog::logger_console, "Invalid mesh after validation");

//		LOG4CPLUS_DEBUG(MeshLog::logger_console,"Setting boundary by sharp edges ...");
		surface_mesh->setBoundarySharpEdges(MeshGenerator3dSurface::param_sharp_edge_threshold);

#ifdef _DEBUG
		if(true){
			int bp_count = 0, bpc_count = 0;
			int be_count = 0;
			MeshViewSet* bset = new MeshViewSet();
			pct = surface_mesh->getPointsCount();
			for(int i = 0; i < pct; i++){
				MeshPoint3d* point = surface_mesh->getPointAt(i);
				if(point->isBorder()) { bset->addPoint(point); bp_count++; }
				if(point->isBorder(TagBorder::CORNER)) bpc_count++;
				int rank = point->getRank();
				for(int j = 0; j < rank; j++){
					MeshEdge3d* edge = point->getEdge(j);
					if(edge->getPointIndex(point) != 0) continue;
					if(edge->isBorder()) { bset->addEdge(edge); be_count++; }
				}
			}
			SHOW_MESH("OFF: create + validate", surface_mesh->getViewSet());
			if(bp_count > 0) {
				bset->addInfo("boundary points", bp_count);
				bset->addInfo("boundary corners", bpc_count);
				bset->addInfo("boundary edges", be_count);
				SHOW_MESH("OFF: boundary", bset);
			}else
				delete bset;
		}
#endif

		LOG4CPLUS_DEBUG(MeshLog::logger_console,"Creating ACS from approximate curvature ...");
		CS3dPtr cs = MeshGenerator3dSurface::createACSFromApproxCurvatureAverage(surface_mesh);
		//CS3dPtr cs = MeshGenerator3dSurface::createACSFromApproxCurvatureDirect(surface_mesh);
		surface_mesh->setControlSpace(cs);
		Metric3dContext mc(cs);

		surface_mesh->clearLocalShapes();
		LOG4CPLUS_DEBUG(MeshLog::logger_console,"Identifying local surfaces & curvesz...");
		//surface_mesh->identifyLocalSurfaces(mc, MeshGenerator3dSurface::param_local_shape_tolerance);
		MeshGenerator3dSurface::identifyLocalSurfaces( mc, surface_mesh );
		//LOG4CPLUS_DEBUG(MeshLog::logger_console,"Identifying local curves ...");
		//surface_mesh->identifyLocalCurves(mc,MeshGenerator3dSurface::param_local_shape_tolerance);

		if( !surface_mesh->isValid() ) 
			LOG4CPLUS_WARN(MeshLog::logger_console, 
				"Invalid mesh before poly-to-triangle");
		surface_mesh->checkLocalSurfaceParams();

		LOG4CPLUS_DEBUG(MeshLog::logger_console,"Converting polys to triangles ...");
		surface_mesh->convertPolysToTriangles();

		if( !surface_mesh->isValid() ) 
			LOG4CPLUS_WARN(MeshLog::logger_console, 
				"Invalid mesh right after poly-to-triangle");
		surface_mesh->checkLocalSurfaceParams();

		MeshGenerator3dSurface::smoothen(mc, surface_mesh, 2, TagExtended::TAG_NONE, 0, 
			MeshData::SM_LAPLACE_MIXED | MeshData::SM_TOP_SWAP | MeshData::SM_DEL_SWAP_COMPLETE); 

		//surface_mesh->checkSurfaceSetCounters( mc, false );

		if( !surface_mesh->isValid() ) LOG4CPLUS_WARN(MeshLog::logger_console, "Invalid mesh after poly-to-triangle");
		surface_mesh->checkLocalSurfaceParams();

		LOG4CPLUS_DEBUG(MeshLog::logger_console,"Surface remeshing with local transformations ...");
		MeshGenerator3dSurface::remeshSurfaceMeshWithLocalTransformations(mc, surface_mesh);

		//surface_mesh->checkSurfaceSetCounters( mc, true );
		//int removed_points = MeshGenerator3dSurface::validateSurfaceMesh(surface_mesh, min_dist);

		SHOW_MESH("OFF: create + remesh", surface_mesh->getViewSet(nullptr, TagExtended::TAG_LOCAL_SURFACE_DOMAIN));
		// check & process

		if( !surface_mesh->isValid() ) LOG4CPLUS_WARN(MeshLog::logger_console, "Invalid mesh after remeshing");

		LOG4CPLUS_DEBUG(MeshLog::logger_console,"Storing final mesh ...");
		// store data - OFF
		surface_mesh->storeOFF( fname_from + ".remeshed.off");
		// store data - XML
		if(surface_mesh->storeXML(fname_to, string("Converted and remeshed from ") + fname_from)){
			pct = surface_mesh->getPointsCount();
			fct = surface_mesh->getFacesCount();
			LOG4CPLUS_INFO(MeshLog::logger_console, "Stored " << pct << " nodes and " << fct << " faces");
			delete surface_mesh;
			return true;
		}else{
			LOG4CPLUS_ERROR(MeshLog::logger_console, "Error writing output file: " << fname_to);
			delete surface_mesh;
			return false;
		}
	}else{
		LOG4CPLUS_ERROR(MeshLog::logger_console, "Unknown extension for input file: " << extension);
	}
	return false;
}


// Join points of mesh, which are too close, plus remove resulting obsolete elements
int MeshBRepXML::mergeCloseNodes(DataVector<DPoint3d> & points, DataSimpleList< DataVector<size_t> > & faces,
								 double min_dist, DataVector<int> & pref)
{
	// show oryginal mesh / .off
	int pct = points.countInt();
	MeshViewSet *set = new MeshViewSet(pct, 0, faces.countInt(), 0);
	for(auto it = faces.iterator(); it.valid(); it.moveNext()){
		DataVector<size_t>& f = it.item();
		int fpct = f.countInt();
		DataVector<DPoint3d> polygon(fpct);
		for(int i = 0; i < fpct; i++) polygon.add(points[f[i]]);
		set->addPolygonConvex(polygon, (fpct == 3) ? 0 : 1);
	}
	SHOW_MESH(".OFF mesh, before", set);

	double min_d2 = min_dist * min_dist;
	int merged = 0;
	int removed_faces = 0;
	DataVector<int> checked_points(pct, -1);

	//update elements & remove edges or elements if degenerated
	for(auto it = faces.iterator(); it.valid(); it.moveNext()){
		auto& f = it.item();
		int ct = f.countInt();
		int i = 0;
		while(i < ct) {
			if( checked_points[f[i]] == 1 ) {
				f.removeOrderedAt(i);
				ct--;
			}else i++;
		}

		i = 0;
		while(ct > 1 && i < ct){
			int i_n = (i+1) % ct;
			int f_n = (int)f[i_n];
			if( points[f[i]].distance2(points[f_n]) < min_d2) {
				checked_points[f_n] = 1; // (to be) removed in other faces
				f.removeOrderedAt(i_n);
				ct--;
				merged++;
			} else{
				checked_points[f_n] = -1; // checked and set as OK
				i++;
			}

		}
		if(ct < 3){
			it.mark();
			removed_faces++;
		}else{
			for(int j = 0; j < ct; j++) pref[f[j]]++;
		}
	}
	if(removed_faces > 0)
		faces.removeMarked();

	int pstat[3] = { 0, 0, 0};
	for(int i = 0; i < pct; i++) pstat[ checked_points[i] + 1]++;
	LOG4CPLUS_DEBUG(MeshLog::logger_console, "pstat[-1]: " << pstat[0]);
	LOG4CPLUS_DEBUG(MeshLog::logger_console, "pstat[ 0]: " << pstat[1]);
	LOG4CPLUS_DEBUG(MeshLog::logger_console, "pstat[ 1]: " << pstat[2]);

	set = new MeshViewSet(pct, 0, faces.countInt(), 0);
	for(auto it = faces.iterator(); it.valid(); it.moveNext()){
		auto& f = it.item();
		int fpct = f.countInt();
		DataVector<DPoint3d> polygon(fpct);
		for(int i = 0; i < fpct; i++) polygon.add(points[f[i]]);
		set->addPolygonConvex(polygon, (fpct == 3) ? 0 : 1);
	}
	SHOW_MESH(".OFF mesh, after", set);

	return merged;
}
