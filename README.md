# Root tracking algorithms
## Background
A system of parametric nonlinear equations generally produce more than one set of feasible roots when solved. These roots change with the input parameter values. There are methods like homotopy continuation or paramotopy for tracking the branches of the roots of polynomial systems with the change in the parameter. 

## Root-tracking methods
This repository contains three implementation of the following algorithms used for root-tracking.
1. Nearest neighbour method
2. Newton-Raphson method based
3. Davidenko's integration method based 

![Tracking of a given path by the 3-3 CDPR](/assets/graphics/33cdpr_tracking.png)

## Illustrative example
The repository contains example usage of the methods to solve a path following problem of a semi-regular Stewart platform manipulator.
Further, an implementation for a 3-3 cable driven parallel robot (CDPR) is also given. 

<p float="middle">
  <img src="/assets/animations/r1.gif" width="33%" />
  <img src="/assets/animations/r2.gif" width="33%" /> 
  <img src="/assets/animations/r3.gif" width="33%" />
  <em>Different branches of the SRSPM tracked using the NR method</em>
</p>

##Usage
The problem formulation and execution is a two step process:
# Mathematica
The problem, i.e., the set of parametric non-linear equations are initially formulated in Mathematica and are stored in the file named Mathematica. These equations are then ported as `C++` expressions using the `FileTemplateApply` function. 

# `C++`
1. To compile and launch the doxygen documentation use `make docs`
2. To know the capabilities of the make file use `make help`
3. To use the codes as is simply use the command `make`


