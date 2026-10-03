

GRAPE- TWO-DIMENSIONAL GRIDS ABOUT AIRFOILS AND             /grape/readme.txt
  OTHER SHAPES BY THE USE OF POISSON'S EQUATION

The files for this program are in the directory grape 
  readme.txt      this file of general description
  input.txt       guide to preparing an input file
  grape.f90       the complete source code
  arc11379.txt    the original program description from COSMIC
  original.src    the original copy of the source code (from COSMIC)

Sample cases for this program are:
  case1.inp       input data
  case1.out       output data for case1.inp
  case2.inp       input data
  case2.out       output data for case2.inp
  case3.inp       input data
  case3.out       output data for case3.inp
  case4.inp       input data
  case4.out       output data for case4.inp

To compile this program for your computer, use the command
   gfortran grape.f90  -o grape.exe
Linux and Macintosh users may prefer to omit the .exe in the file name.
 
This program first asks for the name of the input file. This must
be a file written to conform to input.txt.  After calculating 
the solution, the program produces a file called grape.out that contains
a wealth of information concerning the flow problem to be solved.
