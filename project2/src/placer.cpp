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
double calc_length(point_t *loc);
double calc_density(point_t *loc);
double calc_boundary(point_t *loc);
double delta_length(unsigned cell, int dimen, double dist, point_t *loc);
double delta_density(unsigned cell, int dimen, double dist, double **grid_vals,
        unsigned num_x, unsigned num_y, point_t *loc);
double delta_boundary(unsigned cell, int dimen, double dist, point_t *loc);
double calc_cost(double *x, long int n);
void calc_gradient(double *g, double *x, long int n);
double **get_grid(unsigned &num_x, unsigned &num_y, point_t *loc);
void free_grid(double ** grid_vals, unsigned num_x, unsigned num_y);
double p(double d);

// Global vectors
//vector<point_t> *locations;
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

void place(vector<point_t> &locations, vector<vector<int> > &loc_gates, 
    vector<net_t> &loc_nets, vector<pin_t> &loc_pins, 
    double loc_chipx, double loc_chipy, double loc_unit) {

    unsigned i;

    // Set global variables
    //locations = &loc_locations;
    gates = &loc_gates;
    nets = &loc_nets;
    pins = &loc_pins;
    chipx = loc_chipx;
    chipy = loc_chipy;
    unit = loc_unit;


    // Set initial gate locations
    locations.resize(gates->size());
    for (i = 1; i < locations.size(); i++) {
        locations[i].x = uniform_double() * chipx;
        locations[i].y = uniform_double() * chipy;
    }

    // Call optimizer to minimize cost function
    double *x = &(locations[1].x);
    long int n = 2*(locations.size() - 1);

    //calc_cost(x, n);
    //calc_gradient(g, x, n);
    point_t *loc = reinterpret_cast<point_t *>(x);

    // Set initial values
    grid = loc_chipx/sqrt(gates->size());
    grid = 50;
    radius = 60;
    alpha = grid*radius;

    // Total number of grid points
    grid_points = (static_cast<int>(chipx/grid) + 1)
                  *  (static_cast<int>(chipy/grid) + 1);

    // Total area used by cells on chip
    area = 0;
    for (i = 1; i < gates->size(); i++) {
        area += gates->at(i).size()*unit;
    }

    printf("Before\n");
    w_bp = 10000;
    w_dp = 1;
    w_wl = calc_density(loc)/calc_length(loc);
    printf("After\n");

    return;

    // Optimize cost function
    cg_descent(x, n, NULL, NULL, 1, calc_cost, calc_gradient, NULL, NULL);
 
}

// g - calculated gradient for each cell
// x - location of cells
// n - number of cells
// return cost
double valgrad (double *g, double *x, long int n){
    return 0;
}

double calc_cost(double *x, long int n) {
    double cost, length, density, boundary;
    point_t *loc; 
    loc = reinterpret_cast<point_t *>(x);

    /*
    for(long int i = 0; i < n/2; i++){
        printf("i %ld, x %f, y %f\n", i, loc[i].x, loc[i].y);
    }
    */

    length = calc_length(loc);
    density = calc_density(loc);
    boundary = calc_boundary(loc);

    cost = w_wl * length;
    cost += w_dp * density;
    cost += w_bp * boundary;
    printf("Length: %f, Density: %f, Boundary: %f\n", 
        length, density, boundary);

    printf("Cost: %f\n", cost);
    return cost;
}

void calc_gradient(double *g, double *x, long int n) {
    double delta, h;
    unsigned i, num_x, num_y;
    point_t *loc; 
    double **grid_vals;
    
    loc = reinterpret_cast<point_t *>(x);
    h = grid * H_FACTOR;

    grid_vals = get_grid(num_x, num_y, loc);

    //printf("h %f\n", h);


    for (i = 1; i < gates->size(); i++) {
        //printf("\nCell %d, X Dimension\n", i);
        delta = w_wl * delta_length(i, X_DIM, h, loc) + 
                w_dp * delta_density(i, X_DIM, h, grid_vals, num_x, num_y, loc)+ 
                w_bp * delta_boundary(i, X_DIM, h, loc);

        g[2*i - 2] = delta/h;

        //printf("\nCell %d, Y Dimension\n", i);
        delta = w_wl * delta_length(i, Y_DIM, h, loc) + 
                w_dp * delta_density(i, Y_DIM, h, grid_vals, num_x, num_y, loc)+ 
                w_bp * delta_boundary(i, Y_DIM, h, loc);

        g[2*i - 1] = delta/h;
    }

    free_grid(grid_vals, num_x, num_y);
}

double calc_length(point_t *loc) {
    double length = 0;
    double xmax, xmin, ymax, ymin;
    point_t *cur_loc;
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
            cur_loc = loc + (nets->at(i).gates[j]) - 1;
            //cur_loc = &(locations->at(nets->at(i).gates[j]));

            xmax += exp(cur_loc->x / alpha);
            xmin += exp(-cur_loc->x / alpha);
            ymax += exp(cur_loc->y / alpha);
            ymin += exp(-cur_loc->y / alpha);
        }

        // Find max and min
        xmax = alpha * log(xmax);
        xmin = -alpha * log(xmin);
        ymax = alpha * log(ymax);
        ymin = -alpha * log(ymin);

        // Determine half perimeter length 
        length += xmax - xmin + ymax - ymin;
            
    }

    //printf("\nLength: %f, alpha: %f\n", length, alpha);
    return length;
    
}

double delta_length(unsigned cell, int dimen, double dist, point_t *loc) {
    double xnew, ynew;
    double xmax, xmin, ymax, ymin;
    double xmax_new, xmin_new, ymax_new, ymin_new;
    unsigned i,j;
    double initial = 0, final = 0;
    point_t *cur_loc;

    // Set coordinates for current gate
    xnew = loc[cell - 1].x;
    ynew = loc[cell - 1].y;

    // Update dimension that is changing
    if (dimen == X_DIM) xnew += dist;
    else ynew += dist;
    

    // Go through all nets connected to gate
    for (i = 0; i < gates->at(cell).size(); i++) {
        xmax = 0; xmin = 0; ymax = 0; ymin = 0;

        // Find current net and pin value
        int net = gates->at(cell).at(i);
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
            cur_loc = loc + nets->at(net).gates[j] - 1;

            xmax += exp(cur_loc->x / alpha);
            xmin += exp(-cur_loc->x / alpha);
            ymax += exp(cur_loc->y / alpha);
            ymin += exp(-cur_loc->y / alpha);

            // Find new max and min values
            if (nets->at(net).gates[j] == static_cast<int>(cell)) {
                // Use new values when we are looking at the current gate
                xmax_new += exp(xnew / alpha);
                xmin_new += exp(-xnew / alpha);
                ymax_new += exp(ynew / alpha);
                ymin_new += exp(-ynew / alpha);
            } else {
                // Use current values for all other gates
                xmax_new += exp(cur_loc->x / alpha);
                xmin_new += exp(-cur_loc->x / alpha);
                ymax_new += exp(cur_loc->y / alpha);
                ymin_new += exp(-cur_loc->y / alpha);
            }
        }

        // Determine initial half perimeter length
        initial += alpha * (log(xmax) + log(xmin) + log(ymax) + log(ymin));

        // Determine final half perimeter length
        final += alpha * (log(xmax_new) + log(xmin_new) 
            + log(ymax_new) + log(ymin_new));
    }

    //printf("Length Delta: %f\n", final - initial);
    return final - initial;
}

double **get_grid(unsigned &num_x, unsigned &num_y, point_t *loc) {
    double xdist, ydist, norm_area;
    unsigned i;
    double **grid_vals;

    num_y = static_cast<unsigned>(chipy/grid) + 1;
    num_x = static_cast<unsigned>(chipx/grid) + 1;

    // Allocate grid x dimension
    grid_vals = static_cast<double **>(calloc(num_x, sizeof(double *)));

    // Allocate grid y dimension
    for (i = 0; i < num_x; i++) {
        grid_vals[i] = static_cast<double *>(calloc(num_y, sizeof(double)));
    }

    // Go through gates and set grid values
    for (i = 0; i < gates->size() - 1; i++) {
        unsigned idx_x, idx_y;

        // Get index of lower left grid point within radius
        idx_x = ceil((loc[i].x - radius)/grid);
        idx_y = ceil((loc[i].y - radius)/grid);

        // Ensure x index is in valid range
        if (idx_x >= num_x) idx_x = num_x - 1;
        else if (idx_x < 0) idx_x = 0;

        // Ensure y index is in valid range
        if (idx_y >= num_y) idx_y = num_y - 1;
        else if (idx_y < 0) idx_y = 0;

        // Go through grid cells within radius around gate and find potential
        for (; idx_x * grid <= loc[i].x + radius; idx_x++) {
            if (idx_x >= num_x) break;

            for (; idx_y * grid <= loc[i].y + radius; idx_y++) {
                if (idx_y >= num_y) break;

                // Get for potential calculation
                xdist = abs(idx_x * grid - loc[i].x);
                ydist = abs(idx_y * grid - loc[i].y);
                norm_area = (gates->at(i).size() * unit) / pow(radius, 2);

                // Add potential to grid points
                grid_vals[idx_x][idx_y] += norm_area * p(xdist) * p(ydist);
            }
        }
    }

    return grid_vals;
}

void free_grid(double ** grid_vals, unsigned num_x, unsigned num_y) {
    unsigned i;

    for(i = 0; i < num_x; i++) {
        free(grid_vals[i]);
    }

    free(grid_vals);
}

double calc_density(point_t *loc) {
    double cg, cost, potential;
    double x_pt, y_pt, xdist, ydist, norm_area;
    double **grid_vals;

    unsigned num_y;
    unsigned num_x;
    unsigned i, j;

    cost = 0;
	cg = area/grid_points;

    // Set up grid for density
    grid_vals = get_grid(num_x, num_y, loc);

    // Calculate total cost
    for (i = 0; i < num_x; i++) {
        for (j = 0; j < num_y; j++) {
            printf(" x: %f, y: %f, new potential: %f\n",i*grid, j*grid, grid_vals[i][j]);
            cost += pow(grid_vals[i][j] - cg, 2);
        }
    }
    //printf("area %f, grid_points %d, cg %f\n", area, grid_points, cg);

    double old_cost = 0;

    for (x_pt = 0; x_pt <= chipx; x_pt += grid) {
        for (y_pt = 0; y_pt <= chipy; y_pt += grid) {
            potential = 0;

            for (i = 0; i < gates->size() - 1; i++) {
                xdist = abs(x_pt - loc[i].x);
                ydist = abs(y_pt - loc[i].y);
                norm_area = (gates->at(i).size() * unit) / pow(radius, 2);

                //printf("xdist %f, p(xdist) %f, ydist %f, p(ydist) %f\n", 
                //    xdist, p(xdist), ydist, p(ydist));

                potential += norm_area * p(xdist) * p(ydist);
            }

            printf("x %f, y %f, potential %f\n", x_pt, y_pt, potential);

		    old_cost += pow(potential - cg, 2);
        }
	}

    free_grid(grid_vals, num_x, num_y);

    if (static_cast<int>(cost) != static_cast<int>(old_cost)) {
        printf("Cost (%f) does not equal old cost (%f)\n", cost, old_cost);
    }

    //printf("Density: %f\n", cost);
	return cost;
}

double delta_density(unsigned cell, int dimen, double dist, double **grid_vals,
        unsigned num_x, unsigned num_y, point_t *loc) {

    double cg, initial, final, potential, potential_new;
    double x_pt, y_pt, xdist, ydist, xdist_new, ydist_new, norm_area;
    unsigned i;

    initial = 0; final = 0;
	cg = area/grid_points;
    i = cell - 1;

    double x_min, x_max, y_min, y_max; 

    double x_new = loc[i].x;
    double y_new = loc[i].y;

    double delta = 0;

    // Get box around cell
    x_min = loc[i].x - radius;
    x_max = loc[i].x + radius;
    y_min = loc[i].y - radius;
    y_max = loc[i].y + radius;

    // Expand box for gradient
    if (dimen == X_DIM) {
        if (dist > 0) x_max += dist;
        else x_min += dist;

        x_new += dist;
    } else {
        if (dist > 0) y_max += dist;
        else y_min += dist;

        y_new += dist;
    }

    unsigned idx_min_x, idx_max_x, idx_min_y, idx_max_y, idx_x, idx_y;

    // Get indexes of x and y
    idx_min_x = static_cast<unsigned>(ceil(x_min/grid));
    idx_max_x = static_cast<unsigned>(floor(x_max/grid));
    idx_min_y = static_cast<unsigned>(ceil(y_min/grid));
    idx_max_y = static_cast<unsigned>(floor(y_max/grid));

    // Ensure index is in valid range
    if (idx_max_x >= num_x) idx_max_x = num_x - 1;
    if (idx_max_y >= num_y) idx_max_y = num_y - 1;
    if (idx_min_x < 0) idx_min_x = 0;
    if (idx_min_y < 0) idx_min_y = 0;

    // Go through grid cells within radius around gate and find potential
    for (idx_x = idx_min_x; idx_x <= idx_max_x; idx_x++) {
        for (idx_y = idx_min_y; idx_y <= idx_max_y; idx_y++) {
            double pot, pot_new;

            norm_area = (gates->at(i).size() * unit) / pow(radius, 2);

            // Find initial contribution to potential
            xdist = abs(idx_x * grid - loc[i].x);
            ydist = abs(idx_y * grid - loc[i].y);
            pot = norm_area * p(xdist) * p(ydist);

            // Find new contribution to potential
            xdist_new = abs(idx_x * grid - x_new);
            ydist_new = abs(idx_y * grid - y_new);
            pot_new = norm_area * p(xdist_new) * p(ydist_new);

            // Find delta for this change distance
            delta += pow(grid_vals[idx_x][idx_y] + (pot_new - pot) - cg, 2) - 
                     pow(grid_vals[idx_x][idx_y] - cg, 2);
        }
    }

    	
    for (x_pt = 0; x_pt <= chipx; x_pt += grid) {
        for (y_pt = 0; y_pt <= chipy; y_pt += grid) {
            potential = 0; potential_new = 0;
            
            for (i = 1; i < gates->size(); i++) {
                xdist = abs(x_pt - loc[i - 1].x);
                ydist = abs(y_pt - loc[i - 1].y);
                xdist_new = xdist;
                ydist_new = ydist;

                norm_area = (gates->at(i - 1).size() * unit) / pow(radius, 2);

                if (i == cell) {
                    if (dimen == X_DIM) {
                        xdist_new = abs(x_pt - loc[i - 1].x - dist);
                    } else {
                        ydist_new = abs(y_pt - loc[i - 1].y - dist);
                    }
                } 

                potential += norm_area * p(xdist) * p(ydist);
                potential_new += norm_area * p(xdist_new) * p(ydist_new);
            }

            initial += pow(potential - cg, 2);
            final += pow(potential_new - cg, 2);
        }
	}

    if (static_cast<int>(final - initial) != static_cast<int>(delta)) {
        printf("Delta density new (%f) != old (%f)\n", delta, final - initial);
    }
    //printf("Density Delta: %f\n", final - initial);
    return delta;
}

double calc_boundary(point_t *loc) {
	double cost = 0.0;

	for(unsigned i = 1; i < gates->size(); i++) {
		double xpos = loc[i - 1].x;
		double ypos = loc[i - 1].y;

		if(xpos < 0) cost += pow(xpos / alpha, 2);
		if(ypos < 0) cost += pow(ypos / alpha, 2);
		if(xpos > chipx) cost += pow((xpos - chipx) / alpha, 2);
		if(ypos > chipy) cost += pow((ypos - chipy) / alpha, 2);
	}
    
    //printf("Boundary: %f\n", cost);
    return cost;
}

double delta_boundary(unsigned cell, int dimen, double dist, point_t *loc) {
    double initial = 0;
	double final = 0;
    double xpos, ypos, xpos_new, ypos_new;

    xpos = loc[cell - 1].x;
    ypos = loc[cell - 1].y;

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

    //printf("Boundary Delta: %f\n", final - initial);
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
