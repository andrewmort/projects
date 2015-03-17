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

// Vector of antecedents (pointer to clause that went unit to imply variable)
vector<vector<int> *> variable_antecedent;

// Vector of variables in order they are assigned
vector<int> assigned_variables;

// Vector of indexes to the first variable assigned in that level
vector<unsigned int> level_assignment;

// Watched variable list
vector<watch_list_t> watch_list;

// Pointer to clauses
vector<vector<int> > *clauses;

// Current Level
unsigned int level = 0;

// Function Prototypes
vector<int> *assign_variable(int var);
int backtrack(vector<int> *conflicting_clause);
int next_assignment();
bool check_clause(vector<int> *, unsigned int , int *, unsigned int *);
void resolve(vector<int> *clause, vector<int> *antecedent, int recent_var);

bool DPLL(vector<vector<int> > *clauses_local, int max_var) {
    int var = 0;
    vector<int> *conflicting_clause;

    // Push initial value 0 so we have something to point to initially
    assigned_variables.push_back(0);

    // Assign global variable to clauses to simplify function calls
    clauses = clauses_local;

    // Size watch list and variable assignment vectors to hold all vars
    watch_list.resize(max_var + 1);
    variable_assignment.resize(max_var + 1, UNASSIGNED);
    variable_antecedent.resize(max_var + 1, NULL);

    while (true) {
        // Perform unit propagation with the new assignment var = val
        conflicting_clause = assign_variable(var);

        // When we have a conflict, conflicting clause will not be NULL
        if (conflicting_clause != NULL) {
            if (level == 0) {
                return false;
            }

            var = backtrack(conflicting_clause);
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

            // Ensure there currently as many elements as levels
            if (level_assignment.size() != level) {
                level_assignment.resize(level);
            }

            // Add pointer to last assigned variable to level assignment vector
            level_assignment.push_back(assigned_variables.size() - 1);
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
                    // Add to pos watch list when variable is not inverted
                    watch_list[var].pos.insert(clause);
                } else {
                    // Add to neg watch list when variable is inverted
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

                // Add clause causing variable to go unit as antecedent
                variable_antecedent[abs(new_var)] = clause;

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
 * Return:  Variable to assign
 *
 * Note: Should unassign all assigned variables that disappear due to backtrack
 */
int backtrack(vector<int> *conflicting_clause) {
    vector<int> clause;
    vector<int> *ptr;
    int most_recent; 
    unsigned int next_most_recent_idx, i;
    int ret;

    // Copy conflicting clause into clause vector
    for (i = 0; i < conflicting_clause->size(); i++) {
        clause.push_back(conflicting_clause->at(i));
    }
    
    while (check_clause(&clause, level, &most_recent, &next_most_recent_idx)) {
        resolve(&clause, variable_antecedent[most_recent], most_recent);
    }


    // Find the level containing the next most recent variable assignement
    for (i = level; i > 0; i--) {
        if (next_most_recent_idx > level_assignment[i]) {
            level = i;
        }
    }

    // Get variable assignment for level
    ret = assigned_variables[level_assignment[level] + 1];

    // Set ret to opposite polarity of initial variable value
    ret = variable_assignment[ret] == FALSE ? ret : -ret;

    // Unassign all variables before current level
    for (i = assigned_variables.size() - 1; i > level_assignment[level]; i--){
        variable_assignment[i] = UNASSIGNED;
        assigned_variables.pop_back();
    }

    // Push conflict clause to back of clauses and get pointer to it
    clauses->push_back(clause);
    ptr = &(clauses->back());

    // Update watched list
    unsigned int w;
    for (i = 0, w = 2; i < ptr->size() - w && w > 0; i++) {
        // When variable is unassigned, add it to watched list
        if (variable_assignment[abs(ptr->at(i))] == UNASSIGNED) {
            if (ptr->at(i) > 0) {
                // Add to pos watch list when variable is not inverted
                watch_list[ptr->at(i)].pos.insert(ptr);
            } else {
                // Add to neg watch list when variable is inverted
                watch_list[-ptr->at(i)].neg.insert(ptr);
            }

            // We need one less watched variable now
            w--;
        }
    }

    // Set remaining values in clause as watched variables if we need more
    for(; i < ptr->size() && w > 0; i++, w--) {
        if (ptr->at(i) > 0) {
            // Add to pos watch list when variable is not inverted
            watch_list[ptr->at(i)].pos.insert(ptr);
        } else {
            // Add to neg watch list when variable is inverted
            watch_list[-ptr->at(i)].neg.insert(ptr);
        }
    }

    return ret;
}

/*
 * Return true if the clause contains more than one variable at the current
 * level and sets recent_var to the most recently assigned variable from the
 * clause.
 */
bool check_clause(vector<int> *clause, unsigned int level, int *most_recent,
        unsigned int *next_most_recent_idx) {

    unsigned int largest_idx = 0;
    unsigned int next_largest_idx = 0;

    // Look for most recently assigned variable from clause
    for(unsigned int i = 0; i < clause->size(); i++) {
        int var;
        
        // Get current variable
        var = abs(clause->at(i));

        // Take variable from clause and find location in assigned_variables
        for(unsigned int j = assigned_variables.size() - 1; j >= 0; j--) {

            // Update indexes when j is larger than saved indexes
            if (var == assigned_variables[j]) {
                if (j > largest_idx) {
                    next_largest_idx = largest_idx;
                    largest_idx = j;
                } else if (j > next_largest_idx) {
                    next_largest_idx = j;
                }
            }
        }
    }

    // Set variables with most recent var and idx of next most recent var
    *most_recent = assigned_variables[largest_idx];
    *next_most_recent_idx = next_largest_idx;

    // Return true when most recently assigned variable is in current level
    if (largest_idx > level_assignment[level]) {
        return true;
    }

    return false;
}

/*
 * Set clause to the resolvent determined by clause, antecedent, and the
 * variable recent_var.
 */
void resolve(vector<int> *clause, vector<int> *antecedent, int recent_var) {

    for (unsigned int i = 0; i < clause->size(); i++) {
        if (abs(clause->at(i)) == recent_var) {
            clause->at(i) = clause->back();
            clause->pop_back();

            break;
        }
    }

    //TODO may need to remove duplicate
    for(unsigned int i = 0; i < antecedent->size(); i++) {
        if (abs(clause->at(i)) == recent_var) {
            continue;
        }

        clause->push_back(antecedent->at(i));
    }
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
