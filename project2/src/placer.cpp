#include "placer.h"
#include <vector>
#include <cstdlib>
#include <cstdio>

using namespace std;

double uniform_double();

void place(vector<point_t> &gate_location, vector<vector<int> > &gates, 
    vector<vector<int> > &nets, vector<pin_t> &pins, 
    double chipx, double chipy, double unit) {

    unsigned i;

    // Set initial gate locations
    gate_location.resize(gates.size());
    for (i = 1; i < gate_location.size(); i++) {
        gate_location[i].x = uniform_double() * chipx;
        gate_location[i].y = uniform_double() * chipy;
    }


}

double uniform_double() {
    return rand()/double(RAND_MAX);
}
