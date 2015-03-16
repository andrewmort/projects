int DPLL() {
    int level = 0;

    while (true) {
        if (!unit_propagation()) {
            if (level == 0) {
                return UNSAT;
            }

            level = backtrack();
        } else {
            if (vars_unassigned() == 0) {
                return SAT;
            }

            assign_var();
            
            level++;
        }
    }
}

/*
 * Find any unit clauses and set them until no more unit clauses exist.
 * Implemented using 2-literal watching. Adds all assigned variables to the
 * current level's assigned variable list.
 *
 * Return:  false - if a conflict clause is detected
 *          true - otherwise
 *
 * Note: May need to return pointer to conflict clause or null otherwise
 */
bool unit_propagation();

/*
 * Find the conflict clause using 1UIP and backtrack to the level determined
 * by the algorithm.
 *
 * Return:  level to backtrack 
 *
 * Note: Should unassign all assigned variables that disappear due to backtrack
 */
int backtrack();

/*
 * Return number of variables that still need to be assigned.
 */
int vars_unassigned();

/*
 * Choose and assign the next unassigned variable
 *
 * Note: Add to next level's assigned variable list.
 */
void assign_var();
