#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <unistd.h>
#include <vector>
#include "parser.h"

using namespace std;

void usage(char **argv) {
	printf(" Usage: %s [options] <cktname>\n", argv[0]);
	printf(" Options:-\n");
	printf("    -h\tDisplay this message\n\n");
	exit(-1);
}

int main(int argc, char **argv) {
	char *filename ;
	char c ;
    vector<vector<int> > gates;
    vector<vector<int> > pins;

	// ------------------------------------------------------------
	// Options/command line parsing
	// ------------------------------------------------------------	
	// Default variable settings

	while(1) {
		c = getopt(argc, argv, "h");
		if(c == -1) break ;
		switch(c) {
		case 'h':
		default:
			usage(argv);
		}
	}
	if(argc < 2) 
            usage(argv) ; // correct number of arguments
	//filename = argv[1];
        filename = "bin/lecture.netlist";


	// ------------------------------------------------------------
	// END Options/command line parsing
	// ------------------------------------------------------------	
    parse_netlist_file(gates, pins, filename);

    for(unsigned int i = 0; i < gates.size(); i++){
	cout << i << ": ";
    	for(unsigned int j = 0; j < gates[i].size(); j++){
		cout << gates[i][j] << " ";
	}
	cout << endl;
    }

    return 0;
}
    

