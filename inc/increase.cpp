#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
using namespace std;

int main(int argc, char* argv[])
{
	if(argc != 3){
		cerr << "Use: " << argv[0] << " <cpp file> <prev-line pattern>" << endl;
		return -1;
	}
	
	ifstream f(argv[1]);
	if(!f){
		cerr << "Error opening file " << argv[1] << endl;
		return -2;
	}
	int change_count = 0;
	while(f){
		string line;
		getline(f, line);
		if(line.find(argv[2]) != -1){
			istringstream is(line);
			int number;
			is >> number;
			string ending;
			getline(is, ending);
			cout << "\t\t" << ++number << ending << endl;
			++change_count;
		} else {
			cout << line << endl;
		}
	}
	if(change_count == 0){
		cerr << "No patter [" << argv[2] << "] found" << endl;
		return -3;
	}
	return 0;
}
