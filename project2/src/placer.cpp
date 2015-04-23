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
double potential(double d);
void calc_gradient(double *g, double *x, long int n);
void area_grid_points();
double delta_length(unsigned idx, int dimen, double dist);
double delta_density(unsigned idx, int dimen, double dist);
double delta_boundary(unsigned idx, int dimen, double dist);

// Global vectors
vector<point_t> *locations;
vector<vector<int> > *gates;
vector<net_t> *nets;
vector<pin_t> *pins;
double chipx;
double chipy;
double unit;
double grid, alpha;   		// Gridlength and alpha
double w_wl, w_dp, w_bp;    	// Weights (wirelength, density, boundary)
double area, gridpts;		// Sum of area of gates and number of gridpoints
int radius;                 	// Radius size

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
    radius = 10;
    alpha = grid*radius;
    w_bp = 1;
    w_dp = 1;
    w_wl = 1;

    area_grid_points();

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
    double cost;

    cost = w_wl * calc_length(); 
    cost += w_dp * calc_density();
    cost += w_bp * calc_boundary();

    printf("Cost: %f\n", cost);
    return cost;
}

void calc_gradient(double *g, double *x, long int n) {
    double h = grid * H_FACTOR;
    double delta;
    unsigned i;

    for (i = 1; i < gates->size(); i++) {
        delta = delta_length(i, X_DIM, h) + 
                delta_density(i, X_DIM, h) +
                delta_boundary(i, X_DIM, h);

        g[2*i - 2] = delta/h;

        delta = delta_length(i, Y_DIM, h) + 
                delta_density(i, Y_DIM, h) +
                delta_boundary(i, Y_DIM, h);

        g[2*i - 1] = delta/h;
    }

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

double delta_length(unsigned idx, int dimen, double dist) {
    double xnew, ynew;
    double xmax, xmin, ymax, ymin;
    double xmax_new, xmin_new, ymax_new, ymin_new;
    unsigned i,j;
    double initial = 0, final = 0;
    point_t *loc;

    // Set coordinates for current gate
    xnew = locations->at(idx).x;
    ynew = locations->at(idx).y;

    // Update dimension that is changing
    if (dimen == X_DIM) xnew += dist;
    else ynew += dist;
    

    // Go through all nets connected to gate
    for (i = 0; i < gates->at(idx).size(); i++) {
        xmax = 0; xmin = 0; ymax = 0; ymin = 0;

        // Find current net and pin value
        int net = gates->at(idx).at(i);
        int pin = nets->at(net).pin;

        // Store locations of pins
        if (pin != 0) {
            xmax = exp(pins->at(pin).x / alpha);
            xmin = exp(-pins->at(pin).x / alpha);
            ymax = exp(pins->at(pin).y / alpha);
            ymin = exp(-pins->at(pin).y / alpha);
        }

        // New values are same as regular values
        xmax_new = xmax;
        xmin_new = xmin;
        ymax_new = ymax;
        ymin_new = ymin;

        // Go through all gates connected to current net
        for (j = 0; j < nets->at(net).gates.size(); j++) {
            loc = &(locations->at(nets->at(net).gates[j]));

            xmax += exp(loc->x / alpha);
            xmin += exp(-loc->x / alpha);
            ymax += exp(loc->y / alpha);
            ymin += exp(-loc->y / alpha);

            // Find new max and min values
            if (nets->at(net).gates[j] == static_cast<int>(idx)) {
                // Use new values when we are looking at the current gate
                xmax_new += exp(xnew / alpha);
                xmin_new += exp(-xnew / alpha);
                ymax_new += exp(ynew / alpha);
                ymin_new += exp(-ynew / alpha);
            } else {
                // Use current values for all other gates
                xmax_new += exp(loc->x / alpha);
                xmin_new += exp(-loc->x / alpha);
                ymax_new += exp(loc->y / alpha);
                ymin_new += exp(-loc->y / alpha);
            }
        }

        // Determine initial half perimeter length
        initial += alpha * (log(xmax) + log(xmin) + log(ymax) + log(ymin));

        // Determine final half perimeter length
        final += alpha * (log(xmax_new) + log(xmin_new) 
            + log(ymax_new) + log(ymin_new));
    }

    return final - initial;
}

double calc_density() {
	double cg = area/gridpts;
	int rowlength = static_cast<int>(chipx/grid);
	double cost = 0;
	double temp_cost = 0;
    	
	for(int i = 0; i < gridpts; i++){
		temp_cost = 0;
		for(unsigned j = 1; j < locations->size(); j++) {
			double xdist = (i % rowlength) * grid - locations->at(j).x;
			double ydist = (i / rowlength) * grid - locations->at(j).y;
            double area_factor = gates->at(j).size()*unit / pow(radius, 2);
			temp_cost += area_factor * potential(xdist) * potential(ydist);
		}
		cost += pow(temp_cost - cg, 2);
	}

    printf("Density: %f\n", cost);
	return cost;
}

double delta_density(unsigned idx, int dimen, double dist) {
	double cg = area/gridpts;
	double xnew = locations->at(idx).x;
	double ynew = locations->at(idx).y;

	if(dimen == X_DIM) xnew += dist;
	if(dimen == Y_DIM) ynew += dist;

	int rowlength = static_cast<int>(chipx/grid);
	double initial = 0.0;
	double final = 0.0;

	for(int i = 0; i < gridpts; i++){
		double xdist = (i % rowlength) * grid - locations->at(idx).x;
		double ydist = (i / rowlength) * grid - locations->at(idx).y;
		double xdist_new = (i % rowlength) * grid - xnew;
		double ydist_new = (i / rowlength) * grid - ynew;

		initial += pow(((potential(xdist) * potential(ydist)) - cg), 2);
		final += pow(((potential(xdist_new) * potential(ydist_new)) - cg), 2);
	}
	return final - initial;
}

double calc_boundary() {
	double cost = 0.0;
	for(unsigned i = 1; i < locations->size(); i++) {
		double xpos = locations->at(i).x;
		double ypos = locations->at(i).y;

		if(xpos < 0) cost += pow(xpos / alpha, 2);
		if(ypos < 0) cost += pow(ypos / alpha, 2);
		if(xpos > chipx) cost += pow((xpos - chipx) / alpha, 2);
		if(ypos > chipy) cost += pow((ypos - chipy) / alpha, 2);
	}
    
    printf("Boundary: %f\n", cost);
    return cost;
}

double delta_boundary(unsigned idx, int dimen, double dist) {
    return 0;
}

double uniform_double() {
    return rand()/double(RAND_MAX);
}

double potential(double d) {
	if(-radius/2 <= d && d <= radius/2) return (1-2*pow(d,2)/pow(radius,2));
	else if(radius <= d || d <= -radius) return 0;
    else return (2*pow(abs(d) - radius,2)/pow(radius,2));
}

void area_grid_points() {
	area = 0;
	for(unsigned i = 1; i < gates->size(); i++) {
		area += gates->at(i).size()*unit;
	}

	gridpts = (chipx / grid) * (chipy / grid);

}
