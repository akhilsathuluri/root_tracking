## ----------------------------------------------------------------------
## The current make file generates the executable 'TrackedRoots' for the
## root-tracking problem setup in the main.cpp file
## ----------------------------------------------------------------------
# Add color and helpfile later if needed

# Replace with clang later
CC = g++ # clang++
CFLAGS = -Ofast # -Wall
LIBS = -lgsl -lgslcblas -lm
SRC = main.cpp root_tracker.cpp utils.h inputs.h
OBJ = $(SRC:.cpp = .o)

TrackedRoots: $(OBJ)	## Produces an executable for the root tracking problem formulated in main.cpp
	$(CC) $(CFLAGS) -o TrackedRoots $(OBJ) $(LIBS)

clean:	## Removes the generated *.o files
	rm -f core *.o

docs:	## Generates documentation using doxygen and opens the generated html files
	@echo Building documentation ...
	@doxygen RTdoc.conf
	@firefox html/index.html

.PHONY: clean
help:	## Shows current help
	@sed -ne '/@sed/!s/## //p' $(MAKEFILE_LIST)
