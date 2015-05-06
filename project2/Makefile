# Get all source .cpp files
SRC = $(wildcard src/*.cpp)
SRC_C = $(wildcard src/*.c)

# Get all .h files
INC = $(wildcard include/*.h)

# Get all object files (use SRC and replace src/%.cpp with obj/%.o)
OBJ = $(patsubst src/%.cpp, obj/%.o, $(SRC))
OBJ_C = $(patsubst src/%.c, obj/%.o, $(SRC_C))

# Set tests
#TEST = benchmarks
TEST = benchmarks/toy/toy1

# Find tests
TEST_FILES = $(shell find $(TEST) -type f)

# Set output name
EXENAME = bin/760placer
CHECKER = bin/760proj2checker

# Set debug variable
#DEBUG = -D DEBUG

# Set compiler
CPP = g++ -m32
CC = gcc -m32 -lm -O3 -I include

# Optimizing flags
FLAGS= -Wno-deprecated -Wno-write-strings -Wall -Wold-style-cast \
-Wformat=2 -ansi -pedantic -O3 -I include $(DEBUG)

all: run
#all: $(EXENAME)

# Compile main program by linking all object files (depends on object files)
$(EXENAME): $(OBJ) $(OBJ_C)
	@$(CPP) $(FLAGS) $^ -lz -o $@ 

$(OBJ_C): $(SRC_C) $(INC) Makefile
	@$(CC) -c $< -o $@

# Compile each source file into an object file (depends on cpp file and h files)
obj/%.o: src/%.cpp $(INC) Makefile
	@$(CPP) $(FLAGS) -c $< -o $@

# Compile each source file into an object file (depends on cpp file and h files)
obj/%.oc: src/%.c $(INC) Makefile
	@$(CC) -c $< -o $@

run: $(TEST_FILES)

$(TEST_FILES): $(EXENAME)
	$(EXENAME) $@ output/$(notdir $@).out
	@echo
	$(CHECKER) $@ output/$(notdir $@).out



.PHONY: clean run $(TEST_FILES)

clean:
	@rm -rf $(OBJ) $(EXENAME)
