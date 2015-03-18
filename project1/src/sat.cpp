#include <vector>
#include <queue>
#include <cstdlib>
#include <set>
#include <cstdio>
#include "sat.h"

typedef struct watch_list_t {
    set<int> pos;
    set<int> neg;
} watch_list_t;

// Vector of variable assignemnts (initially set to 0 = UNASSIGNED)
vector<int> variable_assignment;

// Vector of antecedents (index of clause that went unit to imply variable)
vector<int> variable_antecedent;

// Vector of variables in order they are assigned
vector<int> assigned_variables;

// Vector of indexes to the first variable assigned in that level
vector<unsigned int> level_assignment;

// Vector of times we've tried to backtrack for each level
vector<int> level_backtrack;

// Watched variable list
vector<watch_list_t> watch_list;

// Pointer to clauses
vector<vector<int> > *clauses;

// Current Level
unsigned int level = 0;

bool watched_vars_set;

// Function Prototypes
int assign_variable(int var);
int backtrack(int);
int next_assignment();
bool check_clause(vector<int> *, unsigned int *, unsigned int *);
void resolve(vector<int> *clause, int antecedent, int recent_var);

bool DPLL(vector<vector<int> > *clauses_local, int max_var) {
    int var = 0;
    int conflicting_clause;

    // Push initial value 0 so we have something to point to initially
    assigned_variables.push_back(0);
    level_assignment.push_back(0);
    level_backtrack.push_back(0);

    // Assign global variable to clauses to simplify function calls
    clauses = clauses_local;

    // Size watch list and variable assignment vectors to hold all vars
    watch_list.resize(max_var + 1);
    variable_assignment.resize(max_var + 1, UNASSIGNED);
    variable_antecedent.resize(max_var + 1, -1);

    // We haven't set watched variabels yet
    watched_vars_set = false;

    //int count = 0;

    while (true) {
        /*
        if(count++ > 10) {
            return false;
        }
        */


        // Perform unit propagation with the new assignment var = val
        conflicting_clause = assign_variable(var);

        // When we have a conflict, conflicting clause will not be NULL
        if (conflicting_clause != -1) {
#ifdef DEBUG
            printf("Conflict\n");
#endif
            if (level == 0) {
                return false;
            }

            var = backtrack(conflicting_clause);
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
                level_backtrack.resize(level);
            }

            // Add pointer to last assigned variable to level assignment vector
            level_assignment.push_back(assigned_variables.size() - 1);
            level_backtrack.push_back(0);

#ifdef DEBUG
            printf("level %d\n", level);
#endif
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
 * Return:  index to conflict clause or -1 if no conflict
 */
int assign_variable(int var) {
    queue<int> assignment_queue;
    queue<int> antecedent_queue;

#ifdef DEBUG
    printf("Assign %d\n", var);
#endif

    // If var is 0, this is first call to unit_propagation
    if (var == 0) {
        if (!watched_vars_set) {
            // Go through all clauses and add variables to watch list
            for (unsigned int i = 0; i < clauses->size(); i++) {
                // When clause has only one variable it is unit so assign
                if (clauses->at(i).size() == 1) {
                    assignment_queue.push((*clauses)[i][0]);
                    antecedent_queue.push(static_cast<int>(i));

                    continue;
                } 

                // Add first two variables of each clause to watch list
                for (int j = 0; j < 2; j++) {
                    var = (*clauses)[i][j];

                    // Push pointer to clause to end of watch list for that var
                    if (var > 0) {
                        // Add to pos watch list when variable is not inverted
                        watch_list[var].pos.insert(static_cast<int>(i));
                    } else {
                        // Add to neg watch list when variable is inverted
                        watch_list[-var].neg.insert(static_cast<int>(i));
                    }
                }
            } 

            watched_vars_set = true;
        }
    } else {
        // Add variable to assignemnt queue
        assignment_queue.push(var);
        antecedent_queue.push(-1);
    }

    // Go through all assignemnts on the queue
    while (!assignment_queue.empty()) {
        set<int> *wlist;
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

#ifdef DEBUG
        printf("Assignment Queue: %d\n", var);
#endif

        // When variable is already set to value, we don't need to do anything
        if(variable_assignment[abs(var)] == val) {
            continue;
        }

        // When variable is already set, we have a conflict
        if(variable_assignment[abs(var)] != UNASSIGNED) {
            // TODO should never see this.
#ifdef DEBUG
            printf("Var already set: %d\n", variable_assignment[abs(var)]);
#endif

            // Return antecedent since this is conflicting clause
            val = antecedent_queue.front();
            antecedent_queue.pop();

            return val;
        }

        // Set antecedent for current variable
        variable_antecedent[abs(var)] = antecedent_queue.front();
        antecedent_queue.pop();

        // Assign variable and place it in the assigned variable list
        variable_assignment[abs(var)] = val;
        assigned_variables.push_back(abs(var));


#ifdef DEBUG
        printf("Assign Var: %d\n", var);
#endif

        set<int>::iterator it, next_it;

        // Go through all clauses in current watched list
        next_it = wlist->begin();
        while(next_it != wlist->end()) {
            // Must update iterator this way to avoid issues with erase
            it = next_it++;

            int idx = *it;
            bool clause_conflict = true;
            bool clause_unit = false;
            int new_var = 0;

            // Go through all variables in clause
            for(unsigned int i = 0; i < clauses->at(idx).size(); i++) {
                int i_val;
                int i_var;
                set<int> *i_wlist;

                // Get variable assignement
                i_var = clauses->at(idx).at(i);
                i_val = variable_assignment[abs(i_var)];

                // Find watch list of this var that could contain this clause
                if (i_var > 0) {
                    i_wlist = &(watch_list[i_var].pos);
                } else {
                    i_wlist = &(watch_list[-i_var].neg);
                    i_val = -i_val;
                }

                // If watch list contains this clause, this is watched variable
                if(i_wlist->find(idx) != i_wlist->end()) {
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
                        new_var = i_var;
                    }
                } else {
                    // If non watched var is unassigned it becomes watched 
                    if (i_val == UNASSIGNED) {
                        clause_unit = false;
                        clause_conflict = false;
                        new_var = i_var;
                        break;
                    } 
                    
                    // We don't have conflict if any clause is true
                    if (i_val == TRUE) {
                        clause_unit = false;
                        clause_conflict = false;
                        new_var = 0;

                        break;
                    }
                }
            }

            // Stop every if there is a conflict and return clause
            if (clause_conflict) {
#ifdef DEBUG
                printf("return conflict clause %d\n", idx);
#endif
                return idx;
            }

            if (clause_unit) {
                // Add variable to assignemnt queue
                assignment_queue.push(new_var);
                antecedent_queue.push(idx);

                
#ifdef DEBUG
                printf("Unit: %d\n", new_var);
                printf("Antecedent %d: ", idx);

                for(unsigned int a = 0; a < clauses->at(idx).size(); a++) {
                    printf("%d ", clauses->at(idx).at(a));
                }
                printf("\n");
#endif
               

                // Continue to next clause
                continue;
            } 

            if (new_var != 0) {
                // Remove clause from current watch list
                wlist->erase(it);

                // Add clause to new watch list
                if (new_var > 0) {
                    // Add to neg watch list when variable is not inverted
                    watch_list[new_var].pos.insert(idx);
                } else {
                    // Add to pos watch list when variable is inverted
                    watch_list[-new_var].neg.insert(idx);
                }
            }
        }
    }

    return -1;
}
        
/*
 * Find the conflict clause using 1UIP and backtrack to the level determined
 * by the algorithm.
 *
 * Return:  Variable to assign
 *
 * Note: Should unassign all assigned variables that disappear due to backtrack
 */
int backtrack(int conflicting_clause) {
    vector<int> clause;
    int idx;
    unsigned int most_recent_idx, next_most_recent_idx, i;
    int ret;

    
#ifdef DEBUG
    printf("Assigned Vars: \n");
    for (i = 0; i < level_assignment.size(); i++) {
        printf("\tLevel %d: ", i);

        unsigned int j;
        for (j = level_assignment[i]; j < assigned_variables.size(); j++) {
            if (i < level_assignment.size() - 1) {
                if (j >= level_assignment[i + 1]) {
                    break;
                }
            }

            printf(" (%u: %d)", j, assigned_variables[j]);
        }
        printf("\n");
    }

    printf("Conflicting Clause %d\n", conflicting_clause);
#endif
    

    // Copy conflicting clause into clause vector
    for (i = 0; i < clauses->at(conflicting_clause).size(); i++) {
        clause.push_back(clauses->at(conflicting_clause).at(i));
    }

    // Loop while clause contains more than 2 vars from most recent level
    while (check_clause(&clause, &most_recent_idx, &next_most_recent_idx)) {
        // Generate the resolvant and store in clause
        resolve(&clause,
            variable_antecedent[assigned_variables[most_recent_idx]], 
            assigned_variables[most_recent_idx]);
    }

    // After running FirstUIP, as implemented by the while loop above, 
    //  - the final value of clause should be the resolvant
    //  - most_recent_idx should be index of most recently set variable
    //      in the resolvant (also should be in the current level)
    //  - next_most_recent_idx should be index of second to most recently
    //      set variable in resolvant (this is used for backtracking
    //      and we backtrack to the level of this variable)

    // If most_recent is not from current level update next_most_recent
    if (most_recent_idx <= level_assignment[level]) {
        // We use next_most_recent index to determine the level to backtrack
        next_most_recent_idx = most_recent_idx;
    }

    if (next_most_recent_idx != 0) {
        // Find the level containing the next most recent variable assignement
        for (i = level; i > 0; i--) {
            if (next_most_recent_idx > level_assignment[i]) {
                level = i;
                break;
            }
        }

        // We only want to set each var to pos and neg before we give up
        while(level_backtrack[level] > 1) {
            // Reduce level until we find a level we haven't been to twice
            level--;
        }
        
        if (level != 0) {
            // Increase count for number of times we've backtracked this level
            level_backtrack[level]++;
    
            // Get variable assignment for level
            ret = assigned_variables[level_assignment[level] + 1];

            // Set ret to opposite polarity of initial variable value
            ret = variable_assignment[ret] == FALSE ? ret : -ret;
        } else {
            ret = 0;
        }
    } else {
        // Set level to 0 since there's only one var in conflict clause
        level = 0;
        ret = clause.at(0);
    }

    // Unassign all variables after current level
    for (i = assigned_variables.size() - 1; i > level_assignment[level]; i--){
        variable_assignment[assigned_variables[i]] = UNASSIGNED;
        assigned_variables.pop_back();
    }

    // Unassign all level pointers after current level
    for(i = level_assignment.size() - 1; i > level; i--) {
        level_assignment.pop_back();
        level_backtrack.pop_back();
    }

    // Push conflict clause to back of clauses and get pointer to it
    clauses->push_back(clause);
    idx = clauses->size() - 1;
    vector<int> *ptr = &(clauses->at(idx));

    // Update watched list
    unsigned int w;
    for (i = 0, w = 2; i < ptr->size() - w && w > 0; i++) {
        // When variable is unassigned, add it to watched list
        if (variable_assignment[abs(ptr->at(i))] == UNASSIGNED) {
            if (ptr->at(i) > 0) {
                // Add to pos watch list when variable is not inverted
                watch_list[ptr->at(i)].pos.insert(idx);
#ifdef DEBUG
                printf("Watch list 1 pos %d\n", ptr->at(i));
#endif
            } else {
                // Add to neg watch list when variable is inverted
                watch_list[-ptr->at(i)].neg.insert(idx);
#ifdef DEBUG
                printf("Watch list 1 neg %d\n", -ptr->at(i));
#endif
            }

            // We need one less watched variable now
            w--;
        }
    }

    // Set remaining values in clause as watched variables if we need more
    for(; i < ptr->size() && w > 0; i++, w--) {
        if (ptr->at(i) > 0) {
            // Add to pos watch list when variable is not inverted
            watch_list[ptr->at(i)].pos.insert(idx);
#ifdef DEBUG
            printf("Watch list 2 pos %d\n", -ptr->at(i));
#endif
        } else {
            // Add to neg watch list when variable is inverted
            watch_list[-ptr->at(i)].neg.insert(idx);
#ifdef DEBUG
            printf("Watch list 2 neg %d\n", -ptr->at(i));
#endif
        }
    }

#ifdef DEBUG
    printf("Backtrack ret: %d\n", ret);
#endif

    return ret;

}

/*
 * Return true if the clause contains more than one variable at the current
 * level and sets recent_var to the most recently assigned variable from the
 * clause.
 */
bool check_clause(vector<int> *clause, unsigned int *most_recent_idx,
        unsigned int *next_most_recent_idx) {

    unsigned int largest_idx = 0;
    unsigned int next_largest_idx = 0;

    // Look for most recently assigned variable from clause
    for(unsigned int i = 0; i < clause->size(); i++) {
        int var;
        
        // Get current variable
        var = abs(clause->at(i));

        // Take variable from clause and find location in assigned_variables
        for(unsigned int j = assigned_variables.size(); j > 0; j--) {

            // Update indexes when j is larger than saved indexes
            if (var == assigned_variables[j - 1]) {
                if (j - 1 > largest_idx) {
                    next_largest_idx = largest_idx;
                    largest_idx = j - 1;
                } else if (j - 1 != largest_idx && j - 1 > next_largest_idx) {
                    next_largest_idx = j - 1;
                }
            }
        }
    }

    
#ifdef DEBUG
    printf("Check Clause: Level %u, level_idx %u, largest_idx %u, next_largest_idx %u\n",
        level, level_assignment[level], largest_idx, next_largest_idx);
#endif
    

    // Set variables with most recent var and idx of next most recent var
    *most_recent_idx = largest_idx;
    *next_most_recent_idx = next_largest_idx;

    // Return true when next most recently assigned variable is in current level
    if (next_largest_idx > level_assignment[level]) {
#ifdef DEBUG
        printf("ret true\n");
#endif
        return true;
    }

#ifdef DEBUG
    printf("ret false\n");
#endif
    return false;
}

/*
 * Set clause to the resolvent determined by clause, antecedent, and the
 * variable recent_var.
 */
void resolve(vector<int> *clause, int antecedent, int recent_var) {
    set<int> clause_vars;

#ifdef DEBUG
    printf("Resolve clause: ");
    for(unsigned int i = 0; i < clause->size(); i++) {
        printf(" %d", clause->at(i));
    }
    printf("\nResolve antecedent: ");
    for(unsigned int i = 0; i < clauses->at(antecedent).size(); i++) {
        printf(" %d", clauses->at(antecedent).at(i));
    }
    printf("\nResolve var: %d\n", recent_var);
#endif


    for (unsigned int i = 0; i < clause->size(); i++) {
        if (abs(clause->at(i)) == recent_var) {
            // Move last element in clause to take the place of recent var
            clause->at(i) = clause->back();

            // Remove last element since it is now at i
            clause->pop_back();
        }

        // Ensure we have no duplicates
        while (i < clause->size()) {
            // Ensure we have no duplicates
            if (clause_vars.find(clause->at(i)) != clause_vars.end()) {
                // Move last element in clause to take the place of duplicate
                clause->at(i) = clause->back();

                // Remove last element since it is now at i
                clause->pop_back();
            } else {
                // Add variable to set to prevent duplicates
                clause_vars.insert(clause->at(i));

                break;
            }
        }
    }

    // Add all elements from antecedent except recent_var to clause
    for(unsigned int i = 0; i < clauses->at(antecedent).size(); i++) {
        // Remove the recent variable from the antecedent
        if (abs(clauses->at(antecedent).at(i)) == recent_var) {
            continue;
        }

        // Ensure the variable is not already in the clause
        if (clause_vars.find(clauses->at(antecedent).at(i)) 
                != clause_vars.end()) {
            continue;
        }

        // Add variable to clause and to set to prevent duplicates
        clause->push_back(clauses->at(antecedent).at(i));
        clause_vars.insert(clauses->at(antecedent).at(i));
    }

#ifdef DEBUG
    printf("Resolve: ");
    for(unsigned int i = 0; i < clause->size(); i++) {
        printf(" %d", clause->at(i));
    }
    printf("\n");
#endif
}

/*
 * Return variable that should be assigned next or 0 if all variables are 
 * assigned.
 */
int next_assignment() {
    unsigned int idx, size;
    int ret;

#ifdef DEBUG
    printf("Unassigned Vars: ");
    for (unsigned int i = 1; i < variable_assignment.size(); i++) {
        if(variable_assignment[i] == UNASSIGNED) {
            printf(" %u", i);
        }
    }
    printf("\n");
#endif


    size = variable_assignment.size() - 1;
    idx = rand() % (size);

    for (unsigned int i = 0; i < size; i++) {
        ret = static_cast<int>((i + idx) % size + 1);
#ifdef DEBUG
        printf("Try %d\n", ret);
#endif

        if (variable_assignment[ret] == UNASSIGNED) {
            if (rand() % 2) {
                ret = -ret;
            }

#ifdef DEBUG
            printf("Choose %d\n", ret);
#endif

            return ret;
        }
    }

    return 0;

    /*
    for (unsigned int i = 1; i < variable_assignment.size(); i++) {
        if (variable_assignment[i] == UNASSIGNED) {
            return i;
        }
    }

    return 0;
    */
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
