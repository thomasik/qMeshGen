#include <iostream>
#include <fstream>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <cassert>
using namespace std;

int main(int argc, char* argv[])
{
	if(argc < 3){
		cerr << "usage: " << argv[0] << " <file> filter|translate|scale|rotate|cut_sphere" << endl;
		return -1;
	}
	
	ifstream f(argv[1]);
	if(!f){
		cerr << "can't open " << argv[1] << endl;
		return -2;
	}
	
	// read two first lines
	const int LINE_MAX = 1000;
	char line[LINE_MAX];
	f >> line;
	assert( strcmp(line, "OFF") == 0 );
	int pct, fct, bct;
	f >> pct >> fct >> bct;
	
	cout << "OFF" << endl;
	cout << pct << " " << fct << " " << bct << endl;
	
	if(strcmp(argv[2], "filter") == 0){
		// filtering
		if(argc < 6) {
			cerr << "too few parameters for <filter>: x|y|z <min>|inf <max>|inf" << endl;
			return -4;
		}
		int coord = -1;
		if(*argv[3] == 'x') coord = 0;
		else if(*argv[3] == 'y') coord = 1;
		else if(*argv[3] == 'z') coord = 2;
		else{
			cerr << "unknown coordinate: " << argv[3] << endl;
			return -5;
		}
		bool ismin = true, ismax = true;
		double cmin, cmax;
		if(strcmp(argv[4], "inf") == 0) ismin = false;
		else cmin = atof(argv[4]); 
		if(strcmp(argv[5], "inf") == 0) ismax = false;
		else cmax = atof(argv[5]); 
		// go
		double pt[3];
		int itotal = 0, ileft = 0;
		for(int i = 0; i < pct; i++){
			f >> pt[0] >> pt[1] >> pt[2];
			if(!f) break;
			++itotal;
			if(ismin && pt[coord] < cmin) continue;
			if(ismax && pt[coord] > cmax) continue;
			++ileft;
			cout << pt[0] << ' ' << pt[1] << ' ' << pt[2] << endl;
		}
		cerr << ileft << " points left out of total " << itotal << endl;
	}else if(strcmp(argv[2], "translate") == 0){
		// translating
		if(argc < 6) {
			cerr << "too few parameters for <translate>: dx dy dz" << endl;
			return -4;
		}
		double dx = atof(argv[3]);
		double dy = atof(argv[4]);
		double dz = atof(argv[5]);
		// go
		double pt[3];
		int itotal = 0;
		for(int i = 0; i < pct; i++){
			f >> pt[0] >> pt[1] >> pt[2];
			if(!f) break;
			++itotal;
			cout << (pt[0]+dx) << ' ' 
				 << (pt[1]+dy) << ' ' 
				 << (pt[2]+dz) << endl;
		}
		cerr << itotal << " points translated" << endl;
	}else if(strcmp(argv[2], "cut_sphere") == 0){
		// translating
		if(argc < 7) {
			cerr << "too few parameters for <cut_sphere>: x y z r" << endl;
			return -4;
		}
		double x = atof(argv[3]);
		double y = atof(argv[4]);
		double z = atof(argv[5]);
		double r = atof(argv[6]);
		// go
		double pt[3];
		int itotal = 0, ileft = 0;
		for(int i = 0; i < pct; i++){
			f >> pt[0] >> pt[1] >> pt[2];
			if(!f) break;
			++itotal;
			if( (pt[0]-x)*(pt[0]-x) + (pt[1]-y)*(pt[1]-y) + (pt[2]-z)*(pt[2]-z) < r*r)
				continue;
			++ileft;
			cout << pt[0] << ' ' 
				 << pt[1] << ' ' 
				 << pt[2] << endl;
		}
		cerr << ileft << " points left out of " << itotal << " after cut-sphere " << endl;
	}else if(strcmp(argv[2], "scale") == 0){
		// translating
		if(argc < 6) {
			cerr << "too few parameters for <scale>: sx sy sz" << endl;
			return -4;
		}
		double sx = atof(argv[3]);
		double sy = atof(argv[4]);
		double sz = atof(argv[5]);
		// go
		double pt[3];
		int itotal = 0;
		for(int i = 0; i < pct; i++){
			f >> pt[0] >> pt[1] >> pt[2];
			if(!f) break;
			++itotal;
			cout << (pt[0]*sx) << ' ' 
				 << (pt[1]*sy) << ' ' 
				 << (pt[2]*sz) << endl;
		}
		cerr << itotal << " points rescaled" << endl;
	}else if(strcmp(argv[2], "rotate") == 0){
		// rotation around axes
		if(argc < 6) {
			cerr << "too few parameters for <rotate>: ax ay az" << endl;
			return -4;
		}
		// convert to radians
		double ax = atof(argv[3]) * M_PI / 180.0;
		double ay = atof(argv[4]) * M_PI / 180.0;
		double az = atof(argv[5]) * M_PI / 180.0;
		double cx = cos(ax), sx = sin(ax);
		double cy = cos(ay), sy = sin(ay);
		double cz = cos(az), sz = sin(az);
		// go
		double pt[3];
		int itotal = 0;
		for(int i = 0; i < pct; i++){
			f >> pt[0] >> pt[1] >> pt[2];
			if(!f) break;
			++itotal;
			if(ax != 0.0){
				double y = cx * pt[1] - sx * pt[2];
				double z = sx * pt[1] + cx * pt[2];
				pt[1] = y;
				pt[2] = z;
			}
			if(ay != 0.0){
				double x = cy * pt[0] - sy * pt[2];
				double z = sy * pt[0] + cy * pt[2];
				pt[0] = x;
				pt[2] = z;
			}
			if(az != 0.0){
				double x = cz * pt[0] - sz * pt[1];
				double y = sz * pt[0] + cz * pt[1];
				pt[0] = x;
				pt[1] = y;
			}
			cout << pt[0] << ' ' << pt[1] << ' ' << pt[2] << endl;
		}
		cerr << itotal << " points rotated" << endl;
	}else{
		cerr << "unknown command: " << argv[2] << endl;
		return -3;
	}

	// copy rest of the file, unchanged
	
	while(f) {
		f.getline( line, LINE_MAX );
		if(f) cout << line;
	}
	
	return 0;
}
