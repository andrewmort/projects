#ifndef __PLACER_H__
#  define __PLACER_H__

#include <vector>
using std::vector;

typedef struct pin_t {
    int net;
    int x;
    int y;
} pin_t;

typedef struct point_t {
    double x;
    double y;
} point_t;

void place(vector<point_t> &gate_location, vector<vector<int> > &gates,
    vector<vector<int> > &nets, vector<pin_t> &pins, 
    double chipx, double chipy, double unit);

#endif
