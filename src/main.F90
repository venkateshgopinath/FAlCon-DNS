! Code for simulation of thermal convection in a 2D annulus geometry 

! Numerical methods used is a fully spectral method with Chebyshev collocation 
! method in the radial direction and Fourier method in the azimuthal direction.
 
 program maincode
   !$ use OMP_LIB
 use init, only: startmain, finishmain
 use solver, only: solver_init, solver_run, solver_exit, solver_log 
  
 implicit none

 call solver_init   ! Initialize solver

 !call cpu_time(startmain)
 startmain = OMP_GET_WTIME () 

 call solver_run    ! Run the solver

 !call cpu_time(finishmain) 
 finishmain = OMP_GET_WTIME () 
         
 call solver_exit   ! Exit the solver and deallocate

 call solver_log    ! Stores the log file containing run information

 end program maincode
