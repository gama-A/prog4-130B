#ifndef LSF_H
#define LSF_H

#include <iostream>
#include <vector>
#include <random>
#include <algorithm>
#include <cmath>
#include <tuple>

extern "C" {
#include "matrix.c"
}

using namespace std;

struct Point {
    double x, y;
};

tuple<double,double,double> calculateLine(vector<Point> s) {
    double A, B, C;
    int i;
    
    double mat[9], inv[9];
    double y[3], sol[3];
    for(i = 0; i < 9; i++) {
        if(i < 3)
            mat[i] += pow(s[0].x, 2.0-i);
        if(i >= 3 && i < 6)
            mat[i] += pow(s[1].x, 5.0-1);
        if(i >= 6)
            mat[i] += pow(s[2].x, 8.0-i);
    }
    for(i = 0; i < 3; i++) {
        y[i] = s[i].y;
    }
    inverse(3, mat, inv);
    multiplication(3, 3, 1, inv, y, sol);
    return make_tuple(sol[0], sol[1], sol[2]);
}

double findError(Point p, double A, double B, double C) {
    double d;
    d = abs(A*pow(p.x, 2.) + B*p.x + C - p.y) / sqrt(A*A + B*B);
    return d;
}

void findLS(vector<Point> points) {
    int numPoints, i, j, N;
    numPoints = points.size();
    double A, B, C, A_ans, B_ans, C_ans;
    double error, prevError;

    // N = The number of samples to get a 99% chance the line has no outliers
    //N = ceil( log10(1 - 0.99) / log10(1 - pow(1-0.3, 3.0)) );
    N = 2500;
    j = 0;

    // Declaration of vectors to store 3 random points and the remaining n-3 points
    vector<Point> randomSamples(3);
    vector<Point> remainders = points;


    // Finding best coeffs from 1000 random iterations
    while(j < N) {
        // Selection of 3 random points
        random_shuffle(points.begin(), points.end());
        copy(points.begin(), points.begin() + 3, randomSamples.begin());
        copy(points.begin() + 3, points.end(), remainders.begin());
        
        // Calculate quadratic from 3 points
        tie(A, B, C) = calculateLine(randomSamples);

        error = 0.0;
        for(auto p : remainders) {
            error += findError(p, A, B, C);
        }
        error = error / (numPoints - 3);

        if(j != 0) {
            if(error < prevError)
                A_ans = A, B_ans = B, C_ans = C;
        }
        prevError = error;

        j++;
    }
    //cout << A << " " << B << " " << C << "\n";
    double e;
    int count;
    vector<Point> retained;
    
    for(i = 0; i < numPoints; i++) {
        e = findError(points[i], A_ans, B_ans, C_ans);
        if(count == numPoints/2)
            break;
        if(e < prevError) {
            retained.push_back(points[i]);
            count++;
        }
    }
    double mat[9], inv[9];
    double vec[3], sol[3];
    for(i = 0; i < 9; i++) {
        for(auto j : retained) {
            if(i < 3)
                mat[i] += pow(j.x, i);
            if(i >= 3 && i < 6)
                mat[i] += pow(j.x, i-2.0);
            if(i >= 6 && i < 9)
                mat[i] += pow(j.x, i-4.0);
        }
    }
    for(i = 0; i < 3; i++) {
        for(auto j : retained) {
            vec[i] += j.y*pow(j.x, i);
        }
    }
    inverse(3, mat, inv);
    multiplication(3, 3, 1, inv, vec, sol);

    cout << sol[2] << " " << sol[1] << " " << sol[0] << endl;

}

#endif
