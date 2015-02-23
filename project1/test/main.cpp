#include <vector>
#include <cstdlib>
#include <cstdio>

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

int main() {
    vector<clause *> clauses;

    for(int i = 0; i < 5; i++) {
        clause *new_clause = new clause();
        clauses.push_back(new_clause);

        for(int j = 0; j < 5; j++) {
            variable v;

            v.index = j;
            v.value = UNASSIGNED;
            printf("v address: %p\n", &v);

            new_clause->vars.push_back(v);
            printf("insert address: %p\n", &(new_clause->vars[j]));
        }
    }

    for(int i = 0; i < clauses.size(); i++) {
        printf("i: %d\n", i);
        for(int j = 0; j < clauses[i]->vars.size(); j++) {
            printf("\tj: %d - %d\n", j, clauses[i]->vars[j].index);
        }
    }

    for(int i = 0; i < clauses.size(); i++) {
        delete clauses[i];
    }

    return 0;
}


