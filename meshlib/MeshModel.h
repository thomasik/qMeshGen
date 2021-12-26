/////////////////////////////////////////////////////////////////////////////
// MeshModel.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2008-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#if !defined(MESHMODEL_H__INCLUDED)
#define MESHMODEL_H__INCLUDED

class MeshContainer3d;

#include <string>

class MeshModel
{
public:
	enum CmdResult {CM_OK, CM_QUIT, CM_ERROR_RUN, CM_ERROR_PARSE, CM_ERROR_NOMESH, CM_EXCEPTION};
public:
	MeshModel();
	~MeshModel();
public:
	int execute(const std::string& cmd);
	MeshContainer3d* getMesh() const { return m_model_mesh; }
//	int convertSurfPatchtoQuads(const string& fname);
	bool modifiedMesh() const { return m_mesh_modified; }
	void setModifiedMesh(bool m = true) { m_mesh_modified = m; }
	bool modifiedProperties() const { return m_properties_modified; }
	void setModifiedProperties(bool m = true) { m_properties_modified = m; }
private:
	int test(int mode, bool debug_mode = true);
	int testAdaptXML(const std::string& fname1, const std::string& fname2, 
		int adapt_type, int cmp_type, bool debug_mode);
	int testAdaptBin(const std::string& fname1, const std::string& fname2, 
		int adapt_type, int cmp_type, bool debug_mode);
	int testCrack(int mode = 0);
	void setNameFromFilename(const std::string& fname);
private:
	MeshContainer3d* m_model_mesh;
	std::string m_model_name;
	bool m_mesh_modified;
	bool m_properties_modified;
};

#endif // !defined(MESHMODEL_H__INCLUDED)
