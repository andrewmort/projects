#ifndef __PLACER_H__
#  define __PLACER_H__

#include <vector>
using std::vector;

#define X_DIM 0
#define Y_DIM 1
#define H_FACTOR 0.01F

typedef struct pin_t {
    int net;
    int x;
    int y;
} pin_t;

typedef struct net_t {
    vector<int> gates;
    int pin;

    //net_t() {
     ////   pin = 0;
    //}
} net_t;

typedef struct point_t {
    double x;
    double y;
} point_t;

void place(vector<point_t> &locations, vector<vector<int> > &gates,
    vector<net_t> &nets, vector<pin_t> &pins, 
    double chipx, double chipy, double unit);

#endif
