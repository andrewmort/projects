# Get all source .cpp files
SRC = $(wildcard src/*.cpp)

# Get all .h files
INC = $(wildcard include/*.h)

# Get all object files (use SRC and replace src/%.cpp with obj/%.o)
OBJ = $(patsubst src/%.cpp, obj/%.o, $(SRC))

# Set output name
EXENAME = bin/760sat

# Set compiler
CPP = g++ -m32

# Optimizing flags
FLAGS= -Wno-deprecated -Wno-write-strings -Wall -Wold-style-cast \
-Wformat=2 -ansi -pedantic -O3 -I include

all: run
#all: $(EXENAME)

# Compile main program by linking all object files (depends on object files)
$(EXENAME): $(OBJ) 
	@$(CPP) $(FLAGS) $^ -lz -o $@ 

# Compile each source file into an object file (depends on cpp file and h files)
obj/%.o: src/%.cpp $(INC)
	@$(CPP) $(FLAGS) -c $< -o $@


run: $(EXENAME)
	@$(EXENAME) tests/smallBenchmarks/SAT/sanity/sanity2.cnf

.PHONY: clean run

clean:
	@rm -rf $(OBJ) $(EXENAME)
