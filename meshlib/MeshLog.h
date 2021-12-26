/////////////////////////////////////////////////////////////////////////////
// MeshLog.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2008-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#if !defined(MESHLOG_H__INCLUDED)
#define MESHLOG_H__INCLUDED

#include <log4cplus/logger.h>
#include <log4cplus/loggingmacros.h>
using namespace log4cplus;

#include "common.h"

class MeshLog{
public:
	static bool initLog4(const string& log_name = "mesh");
	static bool initLog4(const log4cplus::helpers::SharedObjectPtr<Appender>& append_console, const string& log_name = "mesh");
	static void initLoggerConsole(const log4cplus::helpers::SharedObjectPtr<Appender>& append_console);
	static void initLoggerMesh();
	static void initLoggerMeshStat();
	static void initLoggerMeshCStat(const log4cplus::helpers::SharedObjectPtr<Appender>& append_console);
public:
	static string base_log_name;
	static Logger logger_console; 
	static Logger logger_mesh;
	static Logger logger_mesh_stat;
};

#endif // !defined(MESHLOG_H__INCLUDED)
