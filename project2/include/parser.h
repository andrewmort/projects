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

int parse_file(vector<vector<int> > &gates, vector<net_t> &nets,
    vector<pin_t> &pins, double &chipx, double &chipy, double &unit, 
    const char *netlist_filename);

#endif
