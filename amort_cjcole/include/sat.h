#ifndef H_SAT
#define H_SAT

#include <vector>

using namespace std;

#define TRUE 1
#define FALSE -1
#define UNASSIGNED 0


bool DPLL(vector<vector<int> > *clauses, int max_var);

void print_solution();

#endif
