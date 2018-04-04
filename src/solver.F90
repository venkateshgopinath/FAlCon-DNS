module solver

   use double
   use namelists
   use chebyshev, only: cheballoc, dcheb, chebdealloc
   use init, only: finishmain, startmain, finishsteptime, startsteptime, &
                   & time_matbuild, time_matsolveT, time_matsolveW, time_tran, &
                   & timeNm_maxloop, timeNr_maxloop, finishsteptime, dr, fields_alloc, &
                   & init_grid, init_fields, init_perturbation, fields_dealloc, &
#if FFTW
                   & init_all_fftw_plans, destroy_all_fftw_plans, allocate_restart, &
#elif GG
                   & allocate_restart, &
#endif
                   & deallocate_restart
   use mat_assembly, only: allocate_mat, deallocate_mat
   use nonlin, only: allocate_nonlin, deallocate_nonlin

   use fourierloop, only: allocate_fourierloop_imex, deallocate_fourierloop_imex
   use steptime, only: allocate_steptime_imex, deallocate_steptime_imex, timeloop_imex
   use rhs_create, only: allocate_rhs_imex, deallocate_rhs_imex

   use fourierloop_rk, only: allocate_fourierloop_rk, deallocate_fourierloop_rk
   use steptime_rk, only: allocate_steptime_rk, deallocate_steptime_rk, timeloop_rk
   use rhs_create_rk, only: allocate_rhs_rk, deallocate_rhs_rk

   use fourierloop_imexrk, only: allocate_fourierloop_imexrk, deallocate_fourierloop_imexrk
   use steptime_imexrk, only: allocate_steptime_imexrk, deallocate_steptime_imexrk, timeloop_imexrk
   use rhs_create_imexrk, only: allocate_rhs_imexrk, deallocate_rhs_imexrk

   use timeschemes, only:  init_time_schemes, deallocate_time_schemes

   implicit none

   private

   public :: solver_init, solver_run, solver_exit, solver_log

contains

   subroutine solver_init

      call readNamelists()  ! Read the input parameters from namelist 

      call cheballoc(Nr_max)        ! Allocate Chebyshev matrices and variables
      call dcheb(Nr_max)            ! Make Chebyshev differentiation matrices
      call init_time_schemes(Nm_max,Nr_max,time_scheme_imp,time_scheme_exp,time_scheme_type) ! Initialize time scheme (Allocate rhs and weights)
      call fields_alloc(Nm_max,Np_max,Nr_max,lm,time_scheme_type,time_scheme_imp)    ! Allocate main variables
      call allocate_restart(Nm_max,Nr_max)       ! Allocate restart variables
#if FFTW
      call init_all_fftw_plans(Np_max,Nr_max)    ! Initialize all fftw plans
#endif
      call init_grid(Np_max,Nr_max,rmin)         ! Initialize grid
      call init_fields(Nm_max,Np_max,Nr_max,rmin,rmax,l_restart,n_restart_point,l_add_pert,ampT,lm,time_scheme_type, &
           & time_scheme_imp) ! Initialize main variables/get from checkpoint
      call init_perturbation(Nm_max,Np_max,Nr_max,ampT,l_restart,rmin,rmax,n_init)     ! Initialize perturbation on temperature
      call allocate_mat(Nm_max,Nr_max)           ! Allocate main operator matrices
      call allocate_nonlin(Np_max,Nr_max)        ! Allocate variables in nonlinear products used for multiplication in physical space
      if (time_scheme_type=='IMEXRK') then 
         call allocate_fourierloop_imexrk(Nm_max,Nr_max)   ! Allocate variables in 'fourier loop' for all the modes
         call allocate_rhs_imexrk(Nm_max,Nr_max)           ! Allocate rhs variables
         call allocate_steptime_imexrk(Nm_max,Nr_max)      ! Allocate timeloop variables  
      elseif (time_scheme_type=='RK') then
         call allocate_fourierloop_rk(Nm_max,Nr_max)   ! Allocate variables in 'fourier loop' for all the modes
         call allocate_rhs_rk(Nm_max,Nr_max)           ! Allocate rhs variables
         call allocate_steptime_rk(Nm_max,Nr_max)      ! Allocate timeloop variables  
      elseif (time_scheme_type=='IMEX') then
         call allocate_fourierloop_imex(Nr_max)        ! Allocate variables in 'fourier loop' for all the modes
         call allocate_rhs_imex(Nm_max,Nr_max)         ! Allocate rhs variables
         call allocate_steptime_imex(Nm_max,Nr_max)    ! Allocate timeloop variables  
      end if 
   end subroutine solver_init

   subroutine solver_run

      if (time_scheme_type=='IMEXRK') then
         call timeloop_imexrk(Nm_max,Np_max,Nr_max,eta,CFL,n_time_steps,n_checkpoint,n_snapshot,&
                    & dt,Ra,Pr,l_restart,n_restart,n_restart_point,n_snapshot_point,&
                    & n_KE,n_KEspec,time_scheme_type,time_scheme_imp,time_scheme_exp,tag,dt_coef,dt_max,mBC) 
      elseif (time_scheme_type=='RK') then
         call timeloop_rk(Nm_max,Np_max,Nr_max,eta,CFL,n_time_steps,n_checkpoint,n_snapshot,&
                    & dt,Ra,Pr,l_restart,n_restart,n_restart_point,n_snapshot_point,&
                    & n_KE,n_KEspec,time_scheme_imp,time_scheme_exp,tag,lm,dt_coef,dt_max) 
      elseif (time_scheme_type=='IMEX') then 
         call timeloop_imex(Nm_max,Np_max,Nr_max,eta,CFL,n_time_steps,n_checkpoint,n_snapshot,&
                       & dt,Ra,Pr,mBC,l_restart,n_restart,n_restart_point,n_snapshot_point,&
                       & n_KE,n_KEspec,time_scheme_imp,time_scheme_exp,tag,dt_coef,dt_max) 
      end if

      call cpu_time(finishsteptime)

   end subroutine solver_run

   subroutine solver_exit

      call deallocate_restart()       ! Deallocate restart variables
      if (time_scheme_type=='IMEXRK') then
         call deallocate_steptime_imexrk()      ! Deallocate timestep loop 
         call deallocate_rhs_imexrk()           ! Deallocate rhs variables
         call deallocate_fourierloop_imexrk()   ! Deallocate fourierloop variables
      elseif (time_scheme_type=='RK') then
         call deallocate_steptime_rk()      ! Deallocate timestep loop 
         call deallocate_rhs_rk()           ! Deallocate rhs variables
         call deallocate_fourierloop_rk()   ! Deallocate fourierloop variables
      elseif (time_scheme_type=='IMEX') then
         call deallocate_steptime_imex()      ! Deallocate timestep loop 
         call deallocate_rhs_imex()           ! Deallocate rhs variables
         call deallocate_fourierloop_imex()   ! Deallocate fourierloop variables
      end if   
      call deallocate_nonlin()        ! Deallocate nonlinear variables
      call deallocate_mat()           ! Deallocate operator matrices
      call deallocate_time_schemes(time_scheme_type)  ! Deallocate time schemes
      call fields_dealloc(time_scheme_type)           ! Deallocate variables
#if FFTW
      call destroy_all_fftw_plans()   ! Destroy all fftw plans
#endif
      call chebdealloc()              ! Deallocate Chebyshev matrices and variables

   end subroutine

   subroutine solver_log

      integer :: logunit
      !open(newunit=logunit,file="solver_log.txt",status="unknown",form="formatted", action="write")

      write (logunit,*) "Total time taken =",finishmain-startmain,"seconds."
      write (logunit,*) "Average time taken for time loop =",(finishsteptime-startsteptime)/real(n_time_steps-1,kind=dp),"seconds."
      write (logunit,*) "Time taken for Nr_max loop =",timeNr_maxloop,"seconds."
      write (logunit,*) "Time taken for Nm_max loop =",timeNm_maxloop,"seconds."
      write (logunit,*) "Time taken for Matrix building inside Nm_max loop =",time_matbuild,"seconds."
      write (logunit,*) "Time taken for Solver AX=b for temperature =",time_matsolveT,"seconds."
      write (logunit,*) "Time taken for Solver AX=b for vorticity and streamfunction =",time_matsolveW,"seconds."
      write (logunit,*) "Time taken for 1 inverse transform =",time_tran,"seconds."
      
      !close(logunit)

   end subroutine solver_log

end module solver 






