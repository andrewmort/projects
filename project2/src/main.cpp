#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <unistd.h>
#include <vector>
#include "parser.h"
#include "placer.h"

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
    vector<vector<int> > nets;
    vector<pin_t> pins;

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
    //parse_netlist_file(gates, pins, filename);
    double unit;
    parse_file(gates, nets, pins, unit, filename);

    printf("Gates:\n");
    for(unsigned i = 1; i < gates.size(); i++){
        printf("%d: ", i);
    	for(unsigned j = 0; j < gates[i].size(); j++){
            printf("%d ", gates[i][j]);
	    }
        printf("\n");
    }

    printf("\nNets:\n");
    for(unsigned i = 1; i < nets.size(); i++){
        printf("%d: ", i);
    	for(unsigned j = 0; j < nets[i].size(); j++){
            printf("%d ", nets[i][j]);
	    }
        printf("\n");
    }


    printf("\nPins:\n");

    for(unsigned i = 1; i < pins.size(); i++){
        printf("%d: %d @ (%d, %d)\n", i, pins[i].net, pins[i].x, pins[i].y);
    }

    return 0;
}
    

