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
    vector<clause *> clauses;
    int maxVarIndex;
    bool ret;

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

#ifdef DEBUG
    printf("maxVarIndex: %d\n", maxVarIndex);

    for(unsigned int i = 0; i < clauses.size(); i++) {
        clause *c = clauses[i];
        printf("i = %d\n", i);

        for(unsigned int j = 0; j < c->vars.size(); j++) {
            printf("\tj = %d: %d @ %d\n", j, c->vars[j].index, 
                c->vars[j].value);
        }
    }
    printf("\n");
#endif

    ret = solve(clauses, maxVarIndex);

#ifdef DEBUG
    for(unsigned int i = 0; i < clauses.size(); i++) {
        clause *c = clauses[i];
        printf("i = %d\n", i);

        for(unsigned int j = 0; j < c->vars.size(); j++) {
            printf("\tj = %d: %d @ %d\n", j, c->vars[j].index, 
                c->vars[j].value);
        }
    }
#endif

    // Print solution line
    if (ret) {
        printf("s SATISFIABLE\n");
        print_solution();
    } else {
        printf("s UNSATISFIABLE\n");
    }

    free_vars(clauses);

    return 0;
}
    

