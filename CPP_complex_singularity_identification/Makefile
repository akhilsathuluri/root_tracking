## ----------------------------------------------------------------------
## The current make file generates the executable 'TrackedRoots' for the
## root-tracking problem setup in the main.cpp file
## ----------------------------------------------------------------------
# Add color and helpfile later if needed

# Replace with clang later
CC =  clang++ # g++
CFLAGS = -Ofast # -Wall -Werror
LIBS = -lgsl -lgslcblas -lm
SRC = main.cc root_tracker.cc # utils.hh inputs.hh
TIM = timings.cc root_tracker.cc # utils.hh inputs.hh
OBJ = $(SRC:.cpp = .o)
OBJ_TIM = $(TIM:.cpp = .o)
# For CDPR
SRC_CDPR = cdpr_sei.cc root_tracker.cc
OBJ_CDPR = $(SRC_CDPR:.cpp = .o)

TrackedRoots: $(OBJ)	## Produces an executable for the root tracking problem of an SRSPM
	$(CC) $(CFLAGS) -o TrackedRoots $(OBJ) $(LIBS)

Timings: $(OBJ_TIM) ## Produces an executable measureing the time taken for the root tracking methods formulated in timings.cpp
	$(CC) $(CFLAGS) -o Timings $(OBJ_TIM) $(LIBS)

CDPR: $(OBJ_SRC)	## Produces an executable for the root tracking problem of a CDPR
		$(CC) $(CFLAGS) -o CDPR $(OBJ_CDPR) $(LIBS)

clean:	## Removes the generated *.o files
	rm -f core *.o

docs:	## Generates documentation using doxygen and opens the generated html files
	@echo Building documentation ...
	@doxygen RTdoc.conf
	@firefox html/index.html

.PHONY: clean
help:	## Shows current help
	@sed -ne '/@sed/!s/## //p' $(MAKEFILE_LIST)
