#ifndef LSF_H
#define LSF_H

#include <iostream>
#include <vector>
#include <random>
#include <algorithm>
#include <cmath>
#include <omp.h>

using namespace std;

struct Point {
    double x, y;
};

double calculateLine(Point p0, Point p1, Point p2) {
    double A, B, C, denom;

    denom = ((p0.x - p1.x)*(p0.x - p2.x))*(p1.x - p2.x);
    A = (p2.x * (p1.y - p0.y) + p1.x * (p0.y - p2.y) + p0.x * (p2.y - p1.y)) / denom;
    B = (pow(p2.x, 2.0) * (p0.y - p1.y) + pow(p1.x, 2.0) * (p2.y - p0.y) + pow(p0.x, 2.0) * (p1.y - p2.y)) / denom;
    C = (p1.x * p2.x * (p1.x - p2.x) * p0.y + p2.x * p0.x * (p2.x - p0.x) * p1.y + p0.x * p1.x * (p0.x - p1.x) * p2.y) / denom;

    return A, B, C;
}

double findDist(Point p, double A, double B, double C) {
    double d;
    return d;
}

void findLS(vector<Point> points) {
    int numPoints, i, j, N;
    numPoints = points.size();
    double A, B, C;

    // N = The number of samples to get a 99% chance the line has no outliers
    N = ceil( log10(1 - 0.99) / log10(1 - pow(1-0.3, 3.0)) );

    j = 0;
    while(j < N) {
        // Selection of 3 random points
        vector<Point> randomSamples(3);
        vector<Point> remainders = points;
        vector<std::mt19937> m_RandEngines;
        Point p0, p1, p2;

        shuffle(remainders.begin(), remainders.end(), m_RandEngines[omp_get_thread_num()]); // To avoid picking the same element more than once
		copy(remainders.begin(), remainders.begin() + 3, randomSamples.begin());
        
        // Calculate quadratic from 3 points
        for(i = 0; i < 3; i++) {
            if(i == 0) p0 = randomSamples[i];

            if(i == 1) p1 = randomSamples[i];

            if(i == 2) p2 = randomSamples[i];
        }

        A, B, C = calculateLine(p0, p1, p2);

        double error = 0.0;
        for(auto p : remainders) {
            error += findDist(p, A, B, C);
        }
        error = error / (numPoints - 3);

        j++;
    }
}

#endif
