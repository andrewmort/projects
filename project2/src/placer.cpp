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
double p(double d);
void calc_gradient(double *g, double *x, long int n);
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
double area;     		// Sum of area of gates 
int grid_points;
double radius;                 	// Radius size

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
    grid= 1.8;
    radius = 1.8;
    //alpha = grid*radius;
    alpha = 0.2;
    w_bp = 1;
    w_dp = 1;
    w_wl = 1;

    // Total number of grid points
    grid_points = (static_cast<int>(chipx/grid) + 1)
                  *  (static_cast<int>(chipy/grid) + 1);

    // Total area used by cells on chip
    area = 0;
    for (i = 1; i < gates->size(); i++) {
        area += gates->at(i).size()*unit;
    }


    // Set initial gate locations
    /*
    locations->resize(gates->size());
    for (i = 1; i < locations->size(); i++) {
        locations->at(i).x = uniform_double() * chipx;
        locations->at(i).y = uniform_double() * chipy;
    }
    */

    locations->resize(gates->size());
    locations->at(1).x = 2;
    locations->at(1).y = 2;
    locations->at(2).x = 3;
    locations->at(2).y = 3;

    // Call optimizer to minimize cost function
    double *x = &(locations->at(1).x);
    long int n = 2*(locations->size() - 1);
    double g[200];

    calc_cost(x, n);
    calc_gradient(g, x, n);

    // Optimize cost function
    //cg_descent(x, n, NULL, NULL, 1, calc_cost, calc_gradient, NULL, NULL);

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
    double delta, h;
    unsigned i;

    h = grid * H_FACTOR;
    printf("\nh %f\n", h);


    for (i = 1; i < gates->size(); i++) {
        printf("\nCell %d, X Dimension\n", i);
        delta = delta_length(i, X_DIM, h) + 
                delta_density(i, X_DIM, h) +
                delta_boundary(i, X_DIM, h);

        g[2*i - 2] = delta/h;

        printf("\nCell %d, Y Dimension\n", i);
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

    printf("Length Delta: %f\n", final - initial);
    return final - initial;
}

double calc_density() {
    double cg, cost, potential;
    double x_pt, y_pt, xdist, ydist, norm_area;
    unsigned i;

    cost = 0;
	cg = area/grid_points;
    	
    printf("area %f, grid_points %d, cg %f\n", area, grid_points, cg);

    for (x_pt = 0; x_pt <= chipx; x_pt += grid) {
        for (y_pt = 0; y_pt <= chipy; y_pt += grid) {
            potential = 0;

            for (i = 1; i < locations->size(); i++) {
                xdist = abs(x_pt - locations->at(i).x);
                ydist = abs(y_pt - locations->at(i).y);
                norm_area = (gates->at(i).size() * unit) / pow(radius, 2);

                //printf("xdist %f, p(xdist) %f, ydist %f, p(ydist) %f\n", 
                //    xdist, p(xdist), ydist, p(ydist));

                potential += norm_area * p(xdist) * p(ydist);
            }

            //printf("x %f, y %f, potential %f\n", x_pt, y_pt, potential);

		    cost += pow(potential - cg, 2);
        }
	}

    printf("Density: %f\n", cost);
	return cost;
}

double delta_density(unsigned idx, int dimen, double dist) {
    double cg, initial, final, potential, potential_new;
    double x_pt, y_pt, xdist, ydist, xdist_new, ydist_new, norm_area;
    unsigned i;

    initial = 0; final = 0;
	cg = area/grid_points;

    	
    for (x_pt = 0; x_pt <= chipx; x_pt += grid) {
        for (y_pt = 0; y_pt <= chipy; y_pt += grid) {
            potential = 0; potential_new = 0;
            
            for (i = 1; i < locations->size(); i++) {
                xdist = abs(x_pt - locations->at(i).x);
                ydist = abs(y_pt - locations->at(i).y);
                xdist_new = xdist;
                ydist_new = ydist;

                norm_area = (gates->at(i).size() * unit) / pow(radius, 2);

                if (i == idx) {
                    if (dimen == X_DIM) {
                        xdist_new = abs(x_pt - locations->at(i).x - dist);
                    } else {
                        ydist_new = abs(y_pt - locations->at(i).y - dist);
                    }
                } 

                potential += norm_area * p(xdist) * p(ydist);
                potential_new += norm_area * p(xdist_new) * p(ydist_new);
            }

            initial += pow(potential - cg, 2);
            final += pow(potential_new - cg, 2);
        }
	}

    printf("Density Delta: %f\n", final - initial);
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
    double initial = 0;
	double final = 0;
    double xpos, ypos, xpos_new, ypos_new;

    xpos = locations->at(idx).x;
    ypos = locations->at(idx).y;

    if (dimen == X_DIM) {
        xpos_new = xpos + dist;
        ypos_new = ypos;
    } else {
        xpos_new = xpos;
        ypos_new = ypos + dist;
    }

    // Initial X
    if (xpos < 0) initial += pow(xpos / alpha, 2);
    else if (xpos > chipx) initial += pow((xpos - chipx) / alpha, 2);

    // Initial Y
    if (ypos < 0) initial += pow(ypos / alpha, 2);
    else if (ypos > chipy) initial += pow((ypos - chipy) / alpha, 2);

    // Final X
    if (xpos_new < 0) final += pow(xpos_new / alpha, 2);
    else if (xpos_new > chipx) final += pow((xpos_new - chipx) / alpha, 2);

    // Final Y
    if (ypos_new < 0) final += pow(ypos_new / alpha, 2);
    else if (ypos_new > chipy) final += pow((ypos_new - chipy) / alpha, 2);

    printf("Boundary Delta: %f\n", final - initial);
    return final - initial;
}

double uniform_double() {
    return rand()/double(RAND_MAX);
}

double p(double d) {
	if(0 <= d && d <= radius/2) 
        return (1 - 2 * pow(d,2) / pow(radius,2));
	else if(radius / 2 <= d && d <= radius) 
        return (2 * pow(d - radius,2) / pow(radius,2));
    else 
        return 0;
}
