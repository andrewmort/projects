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
    vector<vector<int> > clauses;
    int maxVarIndex;

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
	filename = argv[1];
        
	// ------------------------------------------------------------
	// END Options/command line parsing
	// ------------------------------------------------------------	

    parse_DIMACS_CNF(clauses, maxVarIndex, filename);

    return 0;
}
    

