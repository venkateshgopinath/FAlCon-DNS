FAlCon (Framework of time schemes for Annular Convection) 
This code is written to solve for 2 D convection in an annulus and analyse different time integration schemes.
Contains a suite of IMEX, IMEXRK and RK time integration schemes.

It was created as part of the PhD work at the Geomagnetism Group, Institute de Physique du Globe de Paris (IPGP), Universite de Paris.

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
