#include "sat.h"
#include <cstdio>
#include <cstdlib>
#include "stdint.h"
#include <set>

#define SAT     -32767
#define UNSAT   0
#define NONE    32767

// Contains each asigned var and pointers to those assigned in clauses
typedef struct assign_var {
    variable var;
    vector <variable *> ptrs;
} assign_var;

// Contains all all assigned variables
vector <assign_var *> assigned_list;

// Set of variables unused
set<int> unused_vars;

int assign(vector<clause *> &clauses, int index, int value);
void unassign(vector<clause *> &clauses, int index, int value);
bool dpll(vector<clause *> &clauses, int index, int value);

// Solve the satisfiability problem
bool solve(vector<clause *> &clauses, int max_var){
    // Create set of unused variables
    for(int i = 1; i <= max_var; i++) {
        unused_vars.insert(i);
    }

    // Call dpll 
    if (dpll(clauses, 1, TRUE)) {
        return true;
    } else {
        unassign(clauses,  1, TRUE);
        return dpll(clauses, 1, FALSE);
    }
}

// Print the satisfiability solution
void print_solution() {
    int val;

    printf("v");

    for(unsigned int i = 0; i < assigned_list.size(); i++) {
        val = assigned_list[i]->var.index;

        if (assigned_list[i]->var.value == FALSE) {
            val = -val;
        }

        printf(" %d", val);
    }

    printf(" 0\n");
}

// Free all allocated variables
void free_vars(vector<clause *> &clauses) {
    assign_var *a;
    clause *c;

    while(assigned_list.size() > 0) {
        a = assigned_list.back();
        assigned_list.pop_back();
        delete a;
    }

    while(clauses.size() > 0) {
        c = clauses.back();
        clauses.pop_back();
        delete c;
    }
}

// Go through all clauses and assign the variable a value. Returns SAT when
// all clauses have been satisfied, UNSAT when there is a conflict, a
// variable name (either complemented or not) that should be assigned 1 to
// make a clause unit, or NONE if nothing changed. A new assigned_list 
// entry is created for the assigned variable and all variables assigned
// in each clause is added to the pointers list in the entry.
int assign(vector<clause *> &clauses, int index, int value) {
    bool unit_decided = false;
    unsigned int unit = 0;
    unsigned int sat_count = 0;
    assign_var *assigned;

#ifdef DEBUG
    printf("\nAssign %d to %d\n", index, value);
#endif

    // Try to remove the current variable from the list of unused vars
    unused_vars.erase(abs(index));

    // Create new assign var and add to assigned list
    assigned = new assign_var();
    assigned_list.push_back(assigned);

    // Set value and index of assigned
    assigned->var.index = index;
    assigned->var.value = value;

    // Go through all clauses to assign variable values
    for (unsigned int i = 0; i < clauses.size(); i++) {
#ifdef DEBUG
        printf("i = %d\n", i);
#endif
        unsigned int unsat_count = 0;

        // Go through clause and assign variable as necessary
        for (unsigned int j = 0; j < clauses[i]->vars.size(); j++) {
#ifdef DEBUG
        printf("\tj = %d: %d @ %d", j, clauses[i]->vars[j].index, 
            clauses[i]->vars[j].value);
#endif
            
            // Assign desired variable value
            if (clauses[i]->vars[j].index == index) {
                // Assign uncomplemented value
                clauses[i]->vars[j].value = value;

                // Add clause variable to assigned ptr list
                assigned->ptrs.push_back(&(clauses[i]->vars[j]));
#ifdef DEBUG
                printf(", assign %d", value);
#endif
            } else if (clauses[i]->vars[j].index == -index) {
                // Assign complemented value
                clauses[i]->vars[j].value = 1 - value;

                // Add clause variable to assigned ptr list
                assigned->ptrs.push_back(&(clauses[i]->vars[j]));
#ifdef DEBUG
                printf(", assign %d", 1 - value);
#endif
            }

            // Determine if clause is sat, unsat, or unit
            if (clauses[i]->vars[j].value == TRUE) {
                // Clause is sat if variable is true
                sat_count++;
#ifdef DEBUG
                printf(", sat\n");
#endif

                // We are done with this clause so go to next one
                break;

            } else if (clauses[i]->vars[j].value == FALSE ) {
                // Increment unsat count if value is false
                unsat_count++;
#ifdef DEBUG
                printf(", unsat count: %d", unsat_count);
#endif

            } else {
                // Otherwise, value is unassigned so update unit

                if (!unit_decided) {
                    if (!unit) {
                        // Unit is not set, so make unit current variable
                        unit = clauses[i]->vars[j].index;

                        // Don't set unit decided yet since we may find another
                        // unassigned variable in clause
#ifdef DEBUG
                        printf(", unit: %d", unit);
#endif
                    } else {
                        // We know clause is not unit, since unit was not 0
                        unit_decided = true;
                        
                        // Reset unit to 0 to indicate there is no unit found
                        unit = 0;
#ifdef DEBUG
                        printf(", not unit: %d", unit);
#endif
                    }
                }
            }
#ifdef DEBUG
            printf("\n");
#endif
        }

        // If all variables in clause are 0, then return UNSAT
        if (unsat_count == clauses[i]->vars.size()) {
#ifdef DEBUG
            printf("\t return UNSAT\n");
#endif
            return UNSAT;
        }

        // The current clause is unit so set is_unit to prevent changing unit
        if (unit != 0 && !unit_decided) {
#ifdef DEBUG
            printf("\t unit = true\n");
#endif
            unit_decided = true;
        } else if (unit_decided && unit == 0) {
            // Reset unit and unit_decided variables
            unit_decided = false;
            unit = 0;
        }
    }

    // If all clauses are satisified, then return SAT
    if (sat_count == clauses.size()) {
        return SAT;
    }

    // If we have a unit clause, return the variable that should be set
    if (unit != 0 && unit_decided) {
        return unit;
    }

    return NONE;
}

// Unassign the last assigned variable
void unassign(vector<clause *> &clauses, int index, int value) {
    assign_var *assigned;
    int cur_idx, cur_val;

#ifdef DEBUG
    printf("Call unassign: %d @ %d\n", index, value);

    print_solution();
#endif


    do {
        // Get last assignment and pop from list
        assigned = assigned_list.back();
        assigned_list.pop_back();

        // Get index and value
        cur_idx = assigned->var.index;
        cur_val = assigned->var.value;

#ifdef DEBUG
        printf("\tUnassign: %d @ %d\n", cur_idx, cur_val);
#endif

        // Reassign all values to UNASSIGNED
        for(unsigned int i = 0; i < assigned->ptrs.size(); i++) {
            assigned->ptrs[i]->value = UNASSIGNED;
        }

        // Delete allocated struct
        delete assigned;

        // Add variable back to unused vars
        unused_vars.insert(cur_idx);

    } while(cur_idx != index || cur_val != value);
}

bool dpll(vector<clause *> &clauses, int index, int value) {
    int ret;

    ret = assign(clauses, index, value); 

#ifdef DEBUG
    printf("Index: %d, Value: %d, ", index, value);
    
    switch(ret) {
        case SAT:
            printf("Return: SAT\n");
            break;
        case UNSAT:
            printf("Return: UNSAT\n");
            break;
        case NONE:
            printf("Return: NONE\n");
            break;
        default:
            printf("Return: Unit, %d\n", ret);
    }
#endif

    if (ret == SAT) {
        // Returns SAT, so return true
        return true;

    } else if (ret == UNSAT) {
        // Return UNSAT, so return false;
        return false;

    } else if (ret == NONE) {

        if(unused_vars.empty()) {
            printf("Unused Vars Empty!\n");
            return false;
        }

        index = *(unused_vars.begin());
#ifdef DEBUG
        printf("** Get index: %d\n", index);
#endif

        // Nothing new happened, so choose new variable
        if (dpll(clauses, index, TRUE)) {
            // DPLL with new variable returns SAT, so return true
            return true;
        } 
#ifdef DEBUG
        printf("** Index %d @ %d Failed\n", index, 1);
#endif
        
        // Undo dpll assignment with new index
        unassign(clauses, index, TRUE);

        // Try dpll with index at false
        if (dpll(clauses, index, FALSE)) {
            return true;
        } 

#ifdef DEBUG
        printf("** Index %d @ %d Failed\n", index, 0);
#endif

        // Undo assignment since it failed
        unassign(clauses, index, FALSE);
        return false;
    } else {
        // Returned a variable, so we have a unit clause
        return dpll(clauses, ret, TRUE);
    }
}

