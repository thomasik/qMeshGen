// MeshLog.cpp: implementation of the MeshLog class.
//
//////////////////////////////////////////////////////////////////////

#include <log4cplus/logger.h>
#include <log4cplus/loggingmacros.h>
#include <log4cplus/fileappender.h>
#include <log4cplus/consoleappender.h>
#include <log4cplus/layout.h>
using namespace log4cplus;
using namespace log4cplus::helpers;


#include "MeshLog.h"

string MeshLog::base_log_name = "mesh";
Logger MeshLog::logger_mesh;
Logger MeshLog::logger_mesh_stat;
Logger MeshLog::logger_console;

void MeshLog::initLoggerConsole(const SharedObjectPtr<Appender>& append_console)
{
	logger_console = Logger::getInstance(LOG4CPLUS_TEXT("mesh.console"));
	logger_console.removeAllAppenders();

	logger_console.addAppender(append_console);
	logger_console.setLogLevel(DEBUG_LOG_LEVEL);
}

void MeshLog::initLoggerMesh()
{
	logger_mesh = Logger::getInstance(LOG4CPLUS_TEXT("mesh"));
	logger_mesh.removeAllAppenders();

	SharedFileAppenderPtr append_meshlog(new FileAppender(base_log_name + ".log"));
	append_meshlog->setName(LOG4CPLUS_TEXT("MeshLog"));
	auto pattern_meshlog = LOG4CPLUS_TEXT("%7r %-5p [%7c{1}] %-40m \t/ %b:%L%n");
	append_meshlog->setLayout(std::make_unique<PatternLayout>(pattern_meshlog));

	logger_mesh.addAppender(SharedAppenderPtr(append_meshlog.get()));
	logger_mesh.setLogLevel(DEBUG_LOG_LEVEL);
}

void MeshLog::initLoggerMeshStat()
{
	logger_mesh_stat = Logger::getInstance(LOG4CPLUS_TEXT("mesh.stat"));
	logger_mesh_stat.removeAllAppenders();

	SharedFileAppenderPtr append_meshstatlog(new FileAppender(base_log_name + ".stat.log"));
	append_meshstatlog->setName(LOG4CPLUS_TEXT("MeshStatLog"));
	auto pattern_meshstatlog = LOG4CPLUS_TEXT("%5c{1}| %m%n");
	append_meshstatlog->setLayout(std::make_unique<PatternLayout>(pattern_meshstatlog));

	logger_mesh_stat.addAppender(SharedAppenderPtr(append_meshstatlog.get()));
	logger_mesh_stat.setLogLevel(INFO_LOG_LEVEL);
}

void MeshLog::initLoggerMeshCStat(const SharedObjectPtr<Appender>& append_console)
{
	auto logger_cstat = Logger::getInstance(LOG4CPLUS_TEXT("mesh.stat.cstat"));
	logger_cstat.removeAllAppenders();

	logger_cstat.addAppender(append_console);
	logger_cstat.setLogLevel(INFO_LOG_LEVEL);
}

bool MeshLog::initLog4(const SharedObjectPtr<Appender>& append_console, const string& log_name)
{
	base_log_name = log_name;

	initLoggerConsole(append_console);
	initLoggerMeshStat();
	initLoggerMeshCStat(append_console);
	initLoggerMesh();

	return true;
}

bool MeshLog::initLog4(const string& log_name)
{
	SharedObjectPtr<Appender> append_console(new ConsoleAppender());
	append_console->setName(LOG4CPLUS_TEXT("MeshConsole"));
	auto pattern_console = LOG4CPLUS_TEXT("%m%n");
	append_console->setLayout(std::make_unique<PatternLayout>(pattern_console));
	//append_console->setLayout(std::make_unique<SimpleLayout>());

	return initLog4(append_console, log_name);
}
