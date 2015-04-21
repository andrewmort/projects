/*********************************************************************
			   ================
			   Parsing code
				18-760
			      Fall 2006
			   ================
**********************************************************************/

#include "parser.h"
#include "placer.h"
#include <iostream>
#include <zlib.h>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <cstring>
using namespace std;

#define STR_BUF 512

//=====================================================================
// DIMACS Parser:

#define CHUNK_LIMIT 1048576

class StreamBuffer {
  gzFile  in;
  char    buf[CHUNK_LIMIT];
  int     pos;
  int     size;

  void assureLookahead() {
    if (pos >= size) {
      pos  = 0;
      size = gzread(in, buf, sizeof(buf)); } }

public:
  StreamBuffer(gzFile i) : in(i), pos(0), size(0) {
    assureLookahead(); }

  int  operator *  () { return (pos >= size) ? EOF : buf[pos]; }
  void operator ++ () { pos++; assureLookahead(); }
};

//-  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void skipWhitespace(StreamBuffer &in) {
    while ((*in >= 9 && *in <= 13) || *in == 32) 
        ++in;
}

void skipLine(StreamBuffer &in) {
    while (true) {
        if(*in == EOF) break;
        if(*in == '\n') {
		++in;
		break;
	}
        ++in;
    }
}

int parseInt(StreamBuffer &in) {
  int val = 0;
  skipWhitespace(in);
  //if (*in < '0' || *in > '\n')
    //fprintf(stderr, "PARSE ERROR! Unexpected char: %c\n", *in), exit(3);
  while (*in >= '0' && *in <= '9') {
    val = val*10 + (*in - '0');
    ++in;
  }
  return val;
}


void readClause(StreamBuffer &in, vector<vector<int> > &clauses) {
  int parsed_lit;
  vector<int> newClause;
  while (true) {
    parsed_lit = parseInt(in);
    if (parsed_lit == 0) break;
    newClause.push_back(parsed_lit);
  }
  clauses.push_back(newClause);
}

void readPin(StreamBuffer &in,  vector<vector<int> > &pins) {
        vector<int> newClause;
        int parsed_lit = parseInt(in);
        if(parsed_lit == 0) return;
        newClause.push_back(parsed_lit);
        pins.push_back(newClause);
}

void parse_netlist_main(StreamBuffer &in, vector<vector<int> > &gates, vector<vector<int> > &pins) {
    int num_gates = 0;
    int num_wires = 0;
    int num_pins = 0;

    // Handle the Gates
    skipLine(in);       // Skip first line
    if(*in == EOF) fprintf(stderr, "PARSE ERROR! File format wrong\n");
    num_gates = parseInt(in);
    skipLine(in);       // Skip remainder of line
    printf("Num gate: %d\n", num_gates);

    for(int i = 0; i < num_gates; i++) {
        skipWhitespace(in);
        if (*in == EOF) break;
        else {
            ++in; // Skip gate name
            num_wires = parseInt(in);
            for(int j = 0; j < num_wires; j++)
                readClause(in, gates);
        }
    }

    // Handle the Pins
    num_pins = parseInt(in);
    skipLine(in);
    for(int i = 0; i < num_pins; i++){
        ++in;
        readPin(in, pins);
    }

}

void parse_netlist(gzFile input_stream, vector<vector<int> > &gates, vector<vector<int> > &pins)
{
	StreamBuffer in(input_stream);
	parse_netlist_main(in, gates, pins);
}

void parse_netlist_file(vector<vector<int> > &gates, 
		vector<vector<int> > &pins, 
		const char *netlist_file)
{
	unsigned int i, j;
	int candidate;
	gzFile in = gzopen(netlist_file, "rb");
	if(in == NULL) {
		fprintf(stderr, "ERROR! Could not open file: %s\n", netlist_file);
	exit(1);
	}
	parse_netlist(in, gates, pins);
	gzclose(in);
}

int parse_file(vector<vector<int> > &gates, vector<vector<int> > &nets,
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
                            int net = atoi(buf.c_str()); 
                            gates[cur_num].push_back(net);
                            nets[net].push_back(cur_num);
                        }

                    // The next line has the number of pins
                    } else if (line == 2 + gates.size()) {
                        pins.resize(atoi(buf.c_str()) + 1);

                    // The remaining lines contain pin information
                    } else {
                        if (line_pos == 0) {
                            cur_num = atoi(buf.c_str());
                        } else if (line_pos == 1) {
                            pins[cur_num].net = atoi(buf.c_str());
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

