#include "convex_hull.h"

// function to swap_Pointers two points
void swap_Pointers(Point* a, Point* b) {
    Point temp = *a;
    *a = *b;
    *b = temp;
}

// function to find the orientation of triplet (p, q, r)
int orientation(Point p, Point q, Point r) {
    double val = (q.y - p.y) * (r.x - q.x) - (q.x - p.x) * (r.y - q.y);
    if (val == 0) return 0;  // colinear
    return (val > 0) ? 1 : 2; // clock or counterclockwise
}

// function to perform Graham's Scan and return the convex hull
Point* grahamScan(instance *inst, int* hullSize) {
    if (inst->nnodes < 3) {
        // convex hull is not possible with less than 3 points
        return NULL;
    }

    // create an array of Points to represent the input points
    Point* points = (Point*)malloc(inst->nnodes * sizeof(Point));
    for (int i = 0; i < inst->nnodes; i++) {
        points[i].x = inst->xcoord[i];
        points[i].y = inst->ycoord[i];
        points[i].id = i;
    }

    // find the point with the lowest y-coordinate (and leftmost if tied)
    int pivot = 0;
    for (int i = 1; i < inst->nnodes; i++) {
        if ((points[i].y < points[pivot].y) || (points[i].y == points[pivot].y && points[i].x < points[pivot].x)) {
            pivot = i;
        }
    }

    // place the pivot at the beginning of the array
    swap_Pointers(&points[0], &points[pivot]);

    // sort the rest of the array based on polar angle from pivot
    qsort(&points[1], inst->nnodes - 1, sizeof(Point), 
          (int (*)(const void *, const void *))orientation);

    // initialize the convex hull with the pivot and the first two sorted points
    Point* hull = (Point*)malloc((inst->nnodes) * sizeof(Point));
    int hullCount = 0;
    hull[hullCount++] = points[0];
    hull[hullCount++] = points[1];
    hull[hullCount++] = points[2];

    // process the rest of the sorted points
    for (int i = 3; i < inst->nnodes; i++) {
        // keep removing top while the angle formed by points[i], hull[hullCount-1], and hull[hullCount-2] makes a non-left turn
        while (hullCount > 1 && orientation(hull[hullCount-2], hull[hullCount-1], points[i]) != 2) {
            hullCount--;
        }
        hull[hullCount++] = points[i];
    }

    // update the size of the convex hull
    *hullSize = hullCount;

    // free memory allocated for points array
    free(points);

    return hull;
}
