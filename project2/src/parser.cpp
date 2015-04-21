/*********************************************************************
			   ================
			   Parsing code
				18-760
			      Fall 2006
			   ================
**********************************************************************/

#include "parser.h"
#include "placer.h"
#include <cstdlib>
#include <cstdio>
#include <string>
#include <cstring>
using namespace std;

#define STR_BUF 1024

int parse_file(vector<vector<int> > &gates, vector<net_t> &nets,
    vector<pin_t> &pins, double &chipx, double &chipy, double &unit, 
    const char *netlist_filename) {

    // Try to open tester file
    FILE *netlist_file = fopen(netlist_filename, "r");
    if(netlist_file == NULL) {
        printf("Error opening %s\n", netlist_filename);
        return 1;
    }

    char str[STR_BUF];
    string buf;
    int line_pos = 0;
    unsigned line = 0;
    int cur_num = 0;
    int net;


    while(fgets(str, STR_BUF, netlist_file) != NULL) {
        // Clear the string buffer for a new line
        buf.clear();
        line_pos = 0;
        line++;

        // Go through line
        for(int i = 0; i < STR_BUF; i++) {

            // Parse string and store result at delimiter
            if(str[i] >= '0' && str[i] <= '9') {
                // Add all numbers to current buf
                buf.push_back(str[i]);

            // Use space or end of line as delimiter
            } else if(str[i] == ' ' || str[i] == '\n') {
                // When buf is not empty, save data
                if(!buf.empty()) {

                    // Line 1 contains general info
                    if (line == 1) {
                        if (line_pos == 0 ) {
                            chipx = atof(buf.c_str());
                        } else if (line_pos == 1) {
                            chipy = atof(buf.c_str());
                        } else if (line_pos == 2){
                            unit = atof(buf.c_str());
                        }

                    // Line 2 contains number of gates and nets
                    } else if (line == 2) {
                        if (line_pos == 0) {
                            gates.resize(atoi(buf.c_str()) + 1);
                        }  else if (line_pos == 1) {
                            nets.resize(atoi(buf.c_str()) + 1);
                            
                        }

                    // The next several lines mark which gates and nets connect
                    } else if (line > 2 && line < 2 + gates.size()) {
                        if (line_pos == 0) {
                            cur_num = atoi(buf.c_str());
                        } else if (line_pos == 1) {
                            // This is just number of nets connected to gate
                        } else {
                            net = atoi(buf.c_str()); 
                            gates[cur_num].push_back(net);
                            nets[net].gates.push_back(cur_num);
                        }

                    // The next line has the number of pins
                    } else if (line == 2 + gates.size()) {
                        pins.resize(atoi(buf.c_str()) + 1);

                    // The remaining lines contain pin information
                    } else {
                        if (line_pos == 0) {
                            cur_num = atoi(buf.c_str());
                        } else if (line_pos == 1) {
                            net = atoi(buf.c_str());
                            pins[cur_num].net = net;
                            nets[net].pin = cur_num;
                        } else if (line_pos == 2) {
                            pins[cur_num].x = atoi(buf.c_str());
                        } else if (line_pos == 3) {
                            pins[cur_num].y = atoi(buf.c_str());
                        }
                    }

                    // Update line position and clear buffer
                    line_pos++;
                    buf.clear();
                }
            }

            // If we happen to see comment or new line go to next line
            if (str[i] == '#' || str[i] == '\n') {
                break;
            }
        }
    }

    fclose(netlist_file);
    return 0;
}

