/*********************************************************************
			   ================
			   Parsing code
				18-760
			      Fall 2006
			   ================
**********************************************************************/

#include "parser.h"
#include <iostream>
using std::ifstream;
#include <zlib.h>
#include <cstdlib>
#include <cstdio>

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

