FAlCon (Framework of time schemes for Annular Convection) 
This code is written to solve for 2 D convection in an annulus and analyse different time integration schemes.
Contains a suite of IMEX, IMEXRK and RK time integration schemes.
For parallelization, it has OpenMP option. 

It was created from scratch as part of the PhD work at the Geomagnetism Group, Institute de Physique du Globe de Paris (IPGP), Universite de Paris.

The research work using the code is published at: 
V. Gopinath, A. Fournier and T. Gastine, An assessment of implicit-explicit time integrators for the pseudo-spectral approximation of Boussinesq 
thermal convection in an annulus, 460, 110965, JCP, 2022. 

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

For questions on usage, please contact: gopinath.venkatesh2@in.bosch.com or venkateshgk.j@gmail.com
