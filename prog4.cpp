// Gamaliel Aristnodo, 8404071

#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>

#include "LSF.h"

using namespace std;

int main() {
    vector<Point> points;
    double x, y;

    string line, input;
    ifstream file("test_prog4/input_1.txt");
    if(file.is_open()) {
        Point temp;
        while(getline(file, line)) {
            stringstream ss(line);
            getline(ss, input, ' ');
            istringstream(input) >> x;
            getline(ss, input);
            istringstream(input) >> y;

            temp.x = x;
            temp.y = y;

            points.push_back(temp);
        }
    }
    
    findLS(points);

    return 0;
}