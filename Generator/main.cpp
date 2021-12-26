
#include <iomanip>
#include "common.h"

#include "MeshData.h"
#include "MeshModel.h"
#include "DEquation.h"
#include "MeshGenerator2d.h"

#include "ControlSpace3dKdTree.h"

#include <log4cplus/logger.h>
#include <log4cplus/loggingmacros.h>
using namespace log4cplus;
using namespace log4cplus::helpers;

//ofstream tex_file;

int main(int argc, const char* argv[])
{
	cout << "Mesh Generator (" << mesh_data.version() << "), Tomasz Jurczyk" << endl;

	//double v1 = SMALL_NUMBER;
	//ControlDataMatrix3d cdm1(v1);

	log4cplus::Initializer initializer;

	int arg_start = 1;
	if (argc > 2 && strcmp(argv[1], "-logname") == 0) {
		MeshLog::initLog4(argv[2]);
		arg_start += 2;
	}
	else {
		MeshLog::initLog4();
	}

	MeshGenerator2d::param_triangulation_type  = 0;
	MeshGenerator2d::param_quality_improvement = 1;

	//checkKdTreeDouble();
	MeshModel model;
	for(int i = arg_start; i < argc; i++){
		if(argv[i][0] == '-'){
			string param(argv[i]+1);
			if(param.compare("logname") == 0){
				LOG4CPLUS_ERROR(MeshLog::logger_console, "-logname has to be the first option.");
				if(i+1 < argc){
					++i;
					//MeshLog::setLogFile(argv[++i]);
				}else{
					LOG4CPLUS_ERROR(MeshLog::logger_console, "Missing file name for -logname option.");
				}
			}else if (param.compare("kdtree-test-cf") == 0) {
				return ControlSpace3dKdTree::checkKdTreeCF(i, argc, argv);
			}else if (param.compare("kdtree-test-ss") == 0) {
				return ControlSpace3dKdTree::checkKdTreeSS(i, argc, argv);
			}
			else {
				auto p = mesh_data.getProperty(param);
				if(!p){
					LOG4CPLUS_WARN(MeshLog::logger_console, "Unknown property: "<< param);
					return -1;
				}
				if(i+1 < argc){
					string value(argv[++i]);
					double d;
					if (DEquation::stringToDouble(value, DEquation::v_auto, d)) {
						p->setDoubleOrIntValue(d);
						LOG4CPLUS_INFO(MeshLog::logger_console, "Param " << param << " set to: " << d);
					}
					else LOG4CPLUS_ERROR(MeshLog::logger_console, "Error parsing value for " << param);
				}else{
					LOG4CPLUS_ERROR(MeshLog::logger_console, "Missing value for " << param);
				}
			}
		}else{
			LOG4CPLUS_INFO(MeshLog::logger_console, "Executing instructions from the file: " << argv[i]);
			ifstream file(argv[i]);
			string buffer;
			getline(file, buffer);
			while(file){
				if(buffer[0] != '#'){
					clock_t clock_start = clock();
					int res = model.execute(buffer);
					double sec = (clock() - clock_start)/(double)CLOCKS_PER_SEC;

					switch(res){
					case MeshModel::CM_OK:
						cout << "** ok\t(" << sec << "s) " << buffer << endl;	break;
					case MeshModel::CM_QUIT:
						cout << "** quiting." << endl;	
						return 0;
					case MeshModel::CM_ERROR_NOMESH:
						cout << "** error - no mesh." << endl;	break;
					case MeshModel::CM_ERROR_RUN:
						cout << "** error executing." << endl;	break;
					case MeshModel::CM_ERROR_PARSE:
						cout << "** error parsing." << endl;	break;
					case MeshModel::CM_EXCEPTION:
						cout << "** error - mesh exception." << endl;	break;
					}
				}
				getline(file, buffer);
			}
		}
	}

	while(true){
		string buffer;
		cout << "> ";
		getline(cin, buffer);
		if(cin){
			while(buffer != ""){
				basic_string <char>::size_type pos = buffer.find(';');
				string comm;
				if(pos != std::string::npos){
					comm = buffer.substr(0, pos);
					buffer.erase(0, pos+1);
				}else{
					comm = buffer;
					buffer = "";
				}

				clock_t clock_start = clock();
				int res = model.execute(comm);
				double sec = (clock() - clock_start)/(double)CLOCKS_PER_SEC;

				switch(res){
				case MeshModel::CM_OK:
					cout << "** ok\t(" << sec << "s) " << comm << endl;	break;
				case MeshModel::CM_QUIT:
					cout << "** quiting." << endl;	
					return 0;
				case MeshModel::CM_ERROR_NOMESH:
					cout << "** error - no mesh." << endl;	break;
				case MeshModel::CM_ERROR_RUN:
					cout << "** error executing." << endl;	break;
				case MeshModel::CM_ERROR_PARSE:
					cout << "** error parsing." << endl;	break;
				case MeshModel::CM_EXCEPTION:
					cout << "** error - mesh exception." << endl;	break;
				}
			}
		}else{
			return 0;
		}
	}
}

