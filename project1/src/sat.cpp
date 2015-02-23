#include "sat.h"
#include <cstdio>
#include "stdint.h"

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

int assign(vector<clause *> &clauses, int index, int value);
void unassign(vector<clause *> &clauses, int index, int value);
bool dpll(vector<clause *> &clauses, int index, int value);

bool solve(vector<clause *> &clauses){
    return dpll(clauses, 1, TRUE);
}

// Go through all clauses and assign the variable a value. Returns SAT when
// all clauses have been satisfied, UNSAT when there is a conflict, a
// variable name (either complemented or not) that should be assigned 1 to
// make a clause unit, or NONE if nothing changed. A new assigned_list 
// entry is created for the assigned variable and all variables assigned
// in each clause is added to the pointers list in the entry.
int assign(vector<clause *> &clauses, int index, int value) {
    bool is_unit = false;
    unsigned int unit = 0;
    unsigned int sat_count = 0;
    assign_var *assigned;

    // Create new assign var and add to assigned list
    assigned = new assign_var();
    assigned_list.push_back(assigned);

    // Set value and index of assigned
    assigned->var.index = index;
    assigned->var.value = value;

    // Go through all clauses to assign variable values
    for (unsigned int i = 0; i < clauses.size(); i++) {
        unsigned int unsat_count = 0;

        // Go through clause and assign variable as necessary
        for (unsigned int j = 0; j < clauses[i]->vars.size(); j++) {
            
            // Assign desired variable value
            if (clauses[i]->vars[j].index == index) {
                // Assign uncomplemented value
                clauses[i]->vars[j].value = value;

                // Add clause variable to assigned ptr list
                assigned->ptrs.push_back(&(clauses[i]->vars[j]));

            } else if (clauses[i]->vars[j].index == -index) {
                // Assign complemented value
                clauses[i]->vars[j].value = 1 - value;

                // Add clause variable to assigned ptr list
                assigned->ptrs.push_back(&(clauses[i]->vars[j]));
            }

            // Determine if clause is sat, unsat, or unit
            if (clauses[i]->vars[j].value == TRUE) {
                // Clause is sat if variable is true
                sat_count++;
                unit = 0;

                // We are done with this clause so go to next one
                continue;

            } else if (clauses[i]->vars[j].value == FALSE ) {
                // Increment unsat count if value is false
                unsat_count++;

            } else {
                // Otherwise, value is unassigned so update unit
                if (!unit) {
                    // Unit is not set, so make unit current variable
                    unit = clauses[i]->vars[j].index;
                } else if (!is_unit) {
                    // This clause has multiple unassigned values so reset unit
                    unit = 0;
                }
            }
        }

        // If all variables in clause are 0, then return UNSAT
        if (unsat_count == clauses[i]->vars.size()) {
            return UNSAT;
        }

        // The current clause is unit so set is_unit to prevent changing unit
        if (unit) {
            is_unit = true;
        }
    }

    // If all clauses are satisified, then return SAT
    if (sat_count == clauses.size()) {
        return SAT;
    }

    // If we have a unit clause, return the variable that should be set
    if (is_unit) {
        return unit;
    }

    return NONE;
}

// Unassign the last assigned variable
void unassign(vector<clause *> &clauses, int index, int value) {
    assign_var *assigned;

    // Get last assignment
    assigned = assigned_list.back();

    // Ensure that we are suppose to unassigning the last assignment
    if (assigned_list.back()->var.index != index ||
            assigned_list.back()->var.value != value) {
        printf("Unassign return early\n");
        return;
    }

    // Reassign all values to UNASSIGNED
    for(unsigned int i = 0; i < assigned->ptrs.size(); i++) {
        assigned->ptrs[i]->value = UNASSIGNED;
    }

    // Delete allocated struct
    delete assigned;
}

bool dpll(vector<clause *> &clauses, int index, int value) {
    int ret;

    ret = assign(clauses, index, value); 

    if (ret == SAT) {
        // Returns SAT, so return true
        return true;

    } else if (ret == UNSAT) {
        // Return UNSAT, so return false;
        return false;

    } else if (ret == NONE) {
        // Nothing new happened, so choose new variable
        if (dpll(clauses, index++, TRUE)) {
            // DPLL with new variable returns SAT, so return true
            return true;

        } else {
            // Undo dpll assignment with index++ at true and try false
            unassign(clauses, index++, TRUE);
            return dpll(clauses, index++, FALSE);
        }
    } else {
        // Returned a variable, so we have a unit clause
        return dpll(clauses, ret, TRUE);
    }
}

/*
 * Data Structures required
 *  - List of all clauses
 *  - Clause
 *      * List of literal name (index)
 *      * List of current literal values
 *      * Whether clause was learned or part of original clause set
 *      * Activitiy of clause, related to how often the clause has
 *          been in conflict
 *  - Keep track of variabels assigned at each level
 *      * Don't want to find all variables to unassign them
 *      * Idea: Keep pointers to varaiables so it's easy to unassign them.
 */

/*
DPLL (cnf f, assignment a) {
    if (a satisfies f)
        return SAT
    else if (some clause is empty)
        return UNSAT
    else if (some clause {p} is unit)
        return DPLL(f, a|p)
    else
        pick an unassigned variable x
        if(DPLL(f, a|x) == SAT)
            return SAT
        else
            return DPLL(f, a|x')
}
*/






