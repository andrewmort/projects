#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <unistd.h>
#include <vector>
#include "parser.h"
#include "sat.h"

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
    vector<vector<int> > clauses;
    int maxVarIndex;

    parse_DIMACS_CNF(clauses, maxVarIndex, filename);

    printf("maxVarIndex: %d\n", maxVarIndex);
    
    for(int i = 0; i < (int) clauses.size(); i++) {
        vector<int> t = clauses[i];
        printf("i = %d\n", i);

        for(int j = 0; j < (int) t.size(); j++) {
            printf("\tj = %d: %d\n", j, t[j]);
        }
    }

    printf("\nsolve: %d\n", (int) solve(clauses));

    return 0;
}
    

