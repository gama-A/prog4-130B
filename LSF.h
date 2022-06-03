#ifndef LSF_H
#define LSF_H

#include <iostream>
#include <vector>
#include <random>

using namespace std;

extern "C" {
#include "matrix.c"
}

struct Point {
    double x, y;
};

void findLS(vector<Point> points) {
    int numPoints, i, j;
    int max_iterations = 1000;

    vector<Point> selected;
    
    i = 0;
    while(i < max_iterations) {
        for(j = 0; j < numPoints; j++) {

        }
        
    }
}

#endif
