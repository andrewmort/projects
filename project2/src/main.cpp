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

    place(gate_location, gates, nets, pins, chipx, chipy, unit);

    printf("Writing Output: %s\n", out_filename);

    if (write_output(gate_location, out_filename) != 0) {
        exit(-3);
    }


    return 0;
}
    

