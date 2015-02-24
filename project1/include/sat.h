#ifndef H_SAT
#define H_SAT

#include <vector>

using namespace std;

#define TRUE 1
#define FALSE 0
#define UNASSIGNED -1


// Variable has index (or name) and value
typedef struct variable {
    int index;
    int value;
} variable;

// A clause contains a list of variables
typedef struct clause {
    vector<variable> vars;
} clause;

bool solve(vector<clause *> &clauses);

void print_solution();
void free_vars(vector<clause *> &clauses);

#endif
