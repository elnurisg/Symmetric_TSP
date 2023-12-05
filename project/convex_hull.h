#include "tsp.h"
#include <stdio.h>
#include <stdlib.h>

typedef struct {
    double x, y;
    int id; // the id in order to identify the node
} Point;


void swap(Point* a, Point* b);
int orientation(Point p, Point q, Point r);
Point* grahamScan(instance *inst, int* hullSize);
