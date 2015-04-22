#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <unistd.h>
#include <vector>
#include "iobitch.h"
#include "placer.h"

using namespace std;

void usage(char **argv) {
	printf(" Usage: %s [options] <netlist> <output>\n", argv[0]);
	printf(" Options:-\n");
	printf("    -h\tDisplay this message\n\n");
	exit(-1);
}

int main(int argc, char **argv) {
	char *in_filename, *out_filename;
	char c ;
    vector<vector<int> > gates;
    vector<net_t> nets;
    vector<pin_t> pins;
    vector<point_t> gate_location;
    double chipx, chipy, unit;

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
	if(argc < 3) 
            usage(argv) ; // correct number of arguments
	in_filename = argv[1];
	out_filename = argv[2];


	// ------------------------------------------------------------
	// END Options/command line parsing
	// ------------------------------------------------------------	
    //parse_netlist_file(gates, pins, in_filename);

    printf("Parsing Netlist: %s\n", in_filename);
    if (parse_netlist(gates, nets, pins, chipx, chipy, unit, in_filename) != 0){
	    exit(-2);
    }


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
        printf("%d: Pin %d - ", i, nets[i].pin);
    	for(unsigned j = 0; j < nets[i].gates.size(); j++){
            printf("%d ", nets[i].gates[j]);
	    }
        printf("\n");
    }


    printf("\nPins:\n");

    for(unsigned i = 1; i < pins.size(); i++){
        printf("%d: %d @ (%d, %d)\n", i, pins[i].net, pins[i].x, pins[i].y);
    }

    place(gate_location, gates, nets, pins, chipx, chipy, unit);

    printf("\nLocations:\n");

    for(unsigned i = 1; i < gate_location.size(); i++){
        printf("%d: (%f, %f)\n", i, gate_location[i].x, gate_location[i].y);
    }

    printf("Writing Output: %s\n", out_filename);

    if (write_output(gate_location, out_filename) != 0) {
        exit(-3);
    }


    return 0;
}
    

