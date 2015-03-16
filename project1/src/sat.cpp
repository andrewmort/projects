#include <vector>
#include <queue>
#include <cstdlib>
#include <set>
#include <cstdio>
#include "sat.h"

typedef struct watch_list_t {
    set<vector<int> *> pos;
    set<vector<int> *> neg;
} watch_list_t;

// Vector of variable assignemnts (initially set to 0 = UNASSIGNED)
vector<int> variable_assignment;

// Vector of variables that are assigned
vector<int> assigned_variables;

// Vector of pointers to the first variable assigned in that level
vector<int *> level_assignment;

// Watched variable list
vector<watch_list_t> watch_list;

// Pointer to clauses
vector<vector<int> > *clauses;

// Function Prototypes
vector<int> *assign_variable(int var);
int backtrack(vector<int> *conflicting_clause);
int next_assignment();

bool DPLL(vector<vector<int> > *clauses_local, int max_var) {
    int level = 0;
    int var = 0;
    vector<int> *conflicting_clause;

    // Assign global variable to clauses to simplify function calls
    clauses = clauses_local;

    // Size watch list and variable assignment vectors to hold all vars
    watch_list.resize(max_var + 1);
    variable_assignment.resize(max_var + 1);

    while (true) {
        // Perform unit propagation with the new assignment var = val
        conflicting_clause = assign_variable(var);

        // When we have a conflict, conflicting clause will not be NULL
        if (conflicting_clause != NULL) {
            /*
            if (level == 0) {
                return false;
            }

            level = backtrack(conflicting_clause);
            */
            printf("Conflict\n");
        } else {
            // Get next variable assignment
            var = next_assignment();

            // When 0 is returned, all variables have been assigned
            if (var == 0) {
                return true;
            }

            // Increment level
            level++;
        }
    }
}

/*
 * Assign the variable var in it's current polarity to 1 (If var is negative
 * the variable abs(var) will be set to 0). Once the variable is assigned,
 * 2-literal watching will be used to find any unit clauses. The routine
 * will continue assigning variables until no more unit clauses exist.
 * When a conflict occurs, a pointer to the conflict clause will be returned.
 * Otherwise, the function will return NULL. If the fucntion is called with
 * var = 0, the 2-literal watch list will be generated using the first two
 * variables of each clause.
 *
 * Return:  pointer to conflict clause or NULL if no conflict
 */
vector<int> *assign_variable(int var) {
    queue<int> assignment_queue;

    //printf("Assign %d\n", var);

    // If var is 0, this is first call to unit_propagation
    if (var == 0) {
        // Go through all clauses and add variables to watch list
        for (unsigned int i = 0; i < clauses->size(); i++) {
            // When clause has only one variable it is unit so assign
            if (clauses->at(i).size() == 1) {
                assignment_queue.push((*clauses)[i][0]);

                continue;
            } 

            // Add first two variables of each clause to watch list
            for (int j = 0; j < 2; j++) {
                vector<int> *clause;

                clause = &((*clauses)[i]);
                var = (*clauses)[i][j];

                // Push pointer to clause to end of watch list for that var
                if (var > 0) {
                    // Add to neg watch list when variable is not inverted
                    watch_list[var].pos.insert(clause);
                } else {
                    // Add to pos watch list when variable is inverted
                    watch_list[-var].neg.insert(clause);
                }
            }
        }
    } else {
        // Add variable to assignemnt queue
        assignment_queue.push(var);
    }

    // Go through all assignemnts on the queue
    while (!assignment_queue.empty()) {
        set<vector<int> *> *wlist;
        int val;

        // Remove next variable assignemnt from queue and set var 
        var = assignment_queue.front();
        assignment_queue.pop();

        // Get value of variable assignment and pointer to watch list
        if (var > 0) {
            val = TRUE;
            wlist = &(watch_list[var].neg);
        } else {
            val = FALSE;
            wlist = &(watch_list[-var].pos);
        }

        // Assign variable and place it in the assigned variable list
        variable_assignment[abs(var)] = val;
        assigned_variables.push_back(abs(var));

        set<vector<int> *>::iterator it, next_it;

        // Go through all clauses in current watched list
        next_it = wlist->begin();
        while(next_it != wlist->end()) {
            // Must update iterator this way to avoid issues with erase
            it = next_it++;

            vector<int> *clause;
            bool clause_conflict = true;
            bool clause_unit = false;
            int new_var = 0;

            // Get the clause which will be operated on
            clause = *it;

            // Go through all variables in clause
            for(unsigned int i = 0; i < clause->size(); i++) {
                int i_val;
                set<vector<int> *> *i_wlist;

                // Get variable assignement
                i_val = variable_assignment[abs(clause->at(i))];

                // Find watch list of this lit that could contain this clause
                if (clause->at(i) > 0) {
                    i_wlist = &(watch_list[clause->at(i)].pos);
                } else {
                    i_wlist = &(watch_list[-clause->at(i)].neg);
                    i_val = -i_val;
                }

                // If watch list contains this clause, this is watched variable
                if(i_wlist->find(clause) != i_wlist->end()) {
                    // If watched var is true, we don't need to do anything else
                    if (i_val == TRUE) {
                        clause_unit = false;
                        clause_conflict = false;
                        new_var = 0;
                        break;
                    }

                    // If watched var is unassigned clause may be unit
                    if (i_val == UNASSIGNED) {
                        clause_unit = true;
                        clause_conflict = false;
                        new_var = clause->at(i);
                    }
                } else {
                    // If non watched var is unassigned it becomes watched 
                    if (i_val == UNASSIGNED) {
                        clause_unit = false;
                        clause_conflict = false;
                        new_var = clause->at(i);
                        break;
                    } 
                    
                    // We don't have conflict if any clause is true
                    if (i_val == TRUE) {
                        clause_unit = false;
                        clause_conflict = false;
                        new_var = 0;
                    }
                }
            }

            // Stop every if there is a conflict and return clause
            if (clause_conflict) {
                return clause;
            }

            if (clause_unit) {
                // Add variable to assignemnt queue
                assignment_queue.push(new_var);

                //printf("Unit: %d\n", new_var);

                // Continue to next clause
                continue;
            } 

            
            if (new_var != 0) {
                // Remove clause from current watch list
                wlist->erase(it);

                // Add clause to new watch list
                if (new_var > 0) {
                    // Add to neg watch list when variable is not inverted
                    watch_list[new_var].pos.insert(clause);
                } else {
                    // Add to pos watch list when variable is inverted
                    watch_list[-new_var].neg.insert(clause);
                }
            }
        }
    }

    return NULL;
}
        
/*
 * Find the conflict clause using 1UIP and backtrack to the level determined
 * by the algorithm.
 *
 * Return:  level to backtrack 
 *
 * Note: Should unassign all assigned variables that disappear due to backtrack
 */
int backtrack(vector<int> *conflicting_clause) {
    /*
    while(conflicting_clause has more than one literal from current level) {
        p = variable in C most recently assigned
        C = Resolve(C, antecedent(p), p);
    }

    add C to clause list
    l = max level in C - {p} 
    Unassigned assignments made at decision levels (l + 1) but not at l
    */

    return 0;
}

/*
 * Return variable that should be assigned next or 0 if all variables are 
 * assigned.
 */
int next_assignment() {
    for (unsigned int i = 1; i < variable_assignment.size(); i++) {
        if (variable_assignment[i] == UNASSIGNED) {
            return i;
        }
    }

    return 0;
}

void print_solution() {
    printf("v");

    for (unsigned int i = 1; i < variable_assignment.size(); i++) {
        int var = i;

        if (variable_assignment[i] == FALSE) {
            var = -var;
        }

        printf(" %d", var);
    }

    printf(" 0\n");
}
