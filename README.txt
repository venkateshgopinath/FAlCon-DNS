Sept 2017
IPGP
Geomagnetism Group
Code for 2 D convection in an annulus.
---------------------------------------
For compiling:

Check the compilers in the Makefile and
go to source directory:

$ cd ./annulus_2d/src
$ make clean
$ make all

For consistency check for a particular case, for example no-slip and cnab2: 

$ cd ./annulus_2d/autotest
$ ./autotest_ns_cnab2.bash 

For running the code: 

$ cd ./annulus_2d/runs
$ cp ./annulus_2d/src/main ./annulus_2d/runs/

Change the input namelist as per requirement and run:
$ ./main input.txt

(OR)

$ ./run.bash
