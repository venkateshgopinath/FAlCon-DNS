Sept 2017
IPGP
Geomagnetism Group
Code for 2 D convection in an annulus.
Contains a suite of IMEX, IMEXRK and RK time integration schemes.
---------------------------------------
For compiling:

Check the compilers in the Makefile and
go to source directory:

$ cd src
$ make clean
$ make all

For running the code: 

$ cd runs
$ cp main ./runs/

Change the input namelist as per requirement and run:
$ ./main input.txt
