#ifndef LSF_H
#define LSF_H

#include <iostream>
#include <vector>
#include <random>
#include <algorithm>
#include <cmath>
#include <tuple>
#include <omp.h>

using namespace std;

struct Point {
    double x, y;
};

tuple<double,double,double> calculateLine(vector<Point> s) {
    double A, B, C, denom;

    denom = (s[0].x - s[1].x)*(s[0].x - s[2].x)*(s[1].x - s[2].x);
    A = (s[2].x * (s[1].y - s[0].y) + s[1].x * (s[0].y - s[2].y) + s[0].x * (s[2].y - s[1].y)) / denom;
    B = (pow(s[2].x, 2.0) * (s[0].y - s[1].y) + pow(s[1].x, 2.0) * (s[2].y - s[0].y) + pow(s[0].x, 2.0) * (s[1].y - s[2].y)) / denom;
    C = (s[1].x * s[2].x * (s[1].x - s[2].x) * s[0].y + s[2].x * s[0].x * (s[2].x - s[0].x) * s[1].y + s[0].x * s[1].x * (s[0].x - s[1].x) * s[2].y) / denom;

    return make_tuple(A, B, C);
}

double findDist(Point p, double A, double B, double C) {
    double d;
    return d;
}

void findLS(vector<Point> points) {
    int numPoints, i, j, N;
    numPoints = points.size();
    double A, B, C, A_ans, B_ans, C_ans;

    // N = The number of samples to get a 99% chance the line has no outliers
    N = ceil( log10(1 - 0.99) / log10(1 - pow(1-0.3, 3.0)) );

    j = 0;
    while(j < N) {
        // Selection of 3 random points
        vector<Point> randomSamples(3);
        vector<Point> remainders = points;
        vector<std::mt19937> m_RandEngines;

        int nThreads = max(1, omp_get_max_threads());
		for (int i = 0; i < nThreads; ++i)
		{				
            random_device SeedDevice;
			m_RandEngines.push_back(std::mt19937(SeedDevice()));
		}

        shuffle(remainders.begin(), remainders.end(), m_RandEngines[omp_get_thread_num()]); // To avoid picking the same element more than once
		copy(remainders.begin(), remainders.begin() + 3, randomSamples.begin());
        remainders.erase(remainders.begin(), remainders.begin() + 3);
        
        // Calculate quadratic from 3 points
        tie(A, B, C) = calculateLine(randomSamples);

        double error = 0.0, prevError;
        for(auto p : remainders) {
            error += findDist(p, A, B, C);
        }
        error = error / (numPoints - 3);

        if(j != 0) {
            if(error < prevEror)
                A_ans = A, B_ans = B, C_ans = C;
        }
        prevError = error

        j++;
    }
}

#endif
