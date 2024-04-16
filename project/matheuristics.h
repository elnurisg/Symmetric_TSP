#include "tsp.h"

/////////////////////////////////// Hard-fixing /////////////////////////////////////////

int hard_fixing(instance *inst, double pfix);
void fixing(CPXENVptr env, CPXLPptr lp, instance *inst, double pfix);
void unfixing(CPXENVptr env, CPXLPptr lp, instance *inst);


/////////////////////////////////// Local branching /////////////////////////////////////////

int local_branching(instance *inst, int k);
void add_local_branching_constraints(CPXENVptr env, CPXLPptr lp, instance *inst, int k);
void remove_last_constraints(CPXENVptr env, CPXLPptr lp);