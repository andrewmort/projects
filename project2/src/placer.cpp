extern "C"{
#include "cg_user.h"
}

#include "placer.h"
#include <vector>
#include <cstdlib>
#include <cstdio>
#include <cmath>

using namespace std;

double uniform_double();
double calc_length();
double calc_density();
double calc_boundary();
double calc_cost(double *x, long int n);
void calc_gradient(double *g, double *x, long int n);

// Global vectors
vector<point_t> *locations;
vector<vector<int> > *gates;
vector<net_t> *nets;
vector<pin_t> *pins;
double chipx;
double chipy;
double unit;
double grid, alpha;   // Gridlength and alpha
double w_wl, w_dp, w_bp;    // Weights (wirelength, density, boundary)
int radius;                 // Radius size

void place(vector<point_t> &loc_locations, vector<vector<int> > &loc_gates, 
    vector<net_t> &loc_nets, vector<pin_t> &loc_pins, 
    double loc_chipx, double loc_chipy, double loc_unit) {

    unsigned i;

    // Set global variables
    locations = &loc_locations;
    gates = &loc_gates;
    nets = &loc_nets;
    pins = &loc_pins;
    chipx = loc_chipx;
    chipy = loc_chipy;
    unit = loc_unit;

    // Set initial values
    grid= 10;
    radius = 2;
    alpha = grid*radius;
    w_bp = 1;
    w_dp = 1;
    w_wl = 1;

    // Set initial gate locations
    locations->resize(gates->size());
    for (i = 1; i < locations->size(); i++) {
        locations->at(i).x = uniform_double() * chipx;
        locations->at(i).y = uniform_double() * chipy;
    }


    // Call optimizer to minimize cost function
    double *x = &(locations->at(1).x);
    long int n = 2*(locations->size() - 1);

    // Optimize cost function
    cg_descent(x, n, NULL, NULL, 1, calc_cost, calc_gradient, NULL, NULL);

}

double calc_cost(double *x, long int n) {
    return w_wl * calc_length() + w_dp * calc_density() + w_bp * calc_boundary();
}

void calc_gradient(double *g, double *x, long int n) {
}

double calc_length() {
    double length = 0;
    double xmax, xmin, ymax, ymin;
    point_t *loc;
    unsigned i,j;

    // Calculate smooth half-perimeter wirelength for each net
    for (i = 1; i < nets->size(); i++) {
        xmax = 0; xmin = 0; ymax = 0; ymin = 0;

        int pin = nets->at(i).pin;

        if (pin != 0) {
            xmax += exp(pins->at(pin).x / alpha);
            xmin += exp(-pins->at(pin).x / alpha);
            ymax += exp(pins->at(pin).y / alpha);
            ymin += exp(-pins->at(pin).y / alpha);
        }

        for (j = 0; j < nets->at(i).gates.size(); j++) {
            loc = &(locations->at(nets->at(i).gates[j]));
            //loc = &(locations->at(1));

            xmax += exp(loc->x / alpha);
            xmin += exp(-loc->x / alpha);
            ymax += exp(loc->y / alpha);
            ymin += exp(-loc->y / alpha);
        }

        // Find max and min
        xmax = alpha * log(xmax);
        xmin = -alpha * log(xmin);
        ymax = alpha * log(ymax);
        ymin = -alpha * log(ymin);

        // Determine half perimeter length 
        length += xmax - xmin + ymax - ymin;
            
    }

    printf("\nLength: %f, alpha: %f\n", length, alpha);
    return length;
    
}

double calc_density() {
    return 0;

}

double calc_boundary() {
    return 0;

}

double uniform_double() {
    return rand()/double(RAND_MAX);
}
