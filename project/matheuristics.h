#include "tsp.h"

int hard_fixing(instance *inst, double pfix);
void fixing(CPXENVptr env, CPXLPptr lp, instance *inst, double pfix);
void unfixing(CPXENVptr env, CPXLPptr lp, instance *inst);