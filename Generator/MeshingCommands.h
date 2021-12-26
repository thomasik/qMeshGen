/////////////////////////////////////////////////////////////////////////////
// MeshingCommands.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2004-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#if !defined(MESHINGCOMMANDS_H__INCLUDED)
#define MESHINGCOMMANDS_H__INCLUDED

#include "MeshData.h"
class MeshContainer3d;

enum CommandResult {CM_OK, CM_QUIT, CM_ERROR_RUN, CM_ERROR_PARSE, CM_ERROR_NOMESH, CM_EXCEPTION};

class MeshingCommands
{
public:
	static int execute(MeshContainer3d* &boundary, const char* command_line);
	static int convertSurfPatchtoQuads(const string& fname);
private:
	static int test(MeshContainer3d* boundary);
};

#endif // !defined(MESHINGCOMMANDS_H__INCLUDED)
