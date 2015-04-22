/*********************************************************************
			   ================
			   CNF Parsing code
				18-760
			      Fall 2006
			   ================
**********************************************************************/

#ifndef __PARSER_H__
#  define __PARSER_H__
#include <vector>
#include "placer.h"
using std::vector;

int parse_netlist(vector<vector<int> > &gates, vector<net_t> &nets,
    vector<pin_t> &pins, double &chipx, double &chipy, double &unit, 
    const char *filename);

int write_output(vector<point_t> &gate_location, const char *filename);

#endif
