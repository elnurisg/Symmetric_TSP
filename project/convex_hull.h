#ifndef CONVEX_HULL_H
#define CONVEX_HULL_H

// #include "tsp.h"

#include <stdio.h>
#include <stdlib.h>
#include "utilities.h"

typedef struct {
    double x, y;
    int id; // the id in order to identify the node
} Point;


void swap_Pointers(Point* a, Point* b);
int orientation(Point p, Point q, Point r);
Point* grahamScan(instance *inst, int* hullSize);


#endif // CONVEX_HULL_H