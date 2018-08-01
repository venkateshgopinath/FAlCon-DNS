module namelists

   use double
   implicit none

   private

   integer, public :: Nr_max
   integer, public :: Nm_max
   integer, public :: Np_max 
   integer, public :: n_checkpoint, n_snapshot, n_restart, &
                      & n_restart_point,n_snapshot_point, n_KE, n_KEspec
   logical, public :: l_restart
   real(kind=dp), public :: rmin
   real(kind=dp), public :: rmax
   character(len=100), public :: time_scheme_type
   character(len=100), public :: time_scheme_imp
   character(len=100), public :: time_scheme_exp
   logical, public :: l_imexrk_started 
   real(kind=dp), public :: CFL
   integer, public :: n_time_steps
   real(kind=dp), public :: dt 
   real(kind=dp), public :: dt_max 
   real(kind=dp), public :: dt_coef 
   real(kind=dp), public :: Ra
   real(kind=dp), public :: Pr
   real(kind=dp), public :: eta
   character(len=100), public :: mBC
   real(kind=dp), public :: ampT
   logical, public :: l_add_pert
   integer, public :: lm
   integer, public :: lagpts
   integer, public :: n_init
   character(len=100), public :: buo_tscheme

   character(len=100), public :: tag
   integer :: argument_count
   integer :: res
   integer :: inputHandle
   character(len=100) :: input_filename
   logical :: nml_exist

   public :: readNamelists, defaultNamelists

contains

   subroutine readNamelists()

      !Define namelists --------------------------
      namelist/grid/Nr_max,Nm_max

      namelist/physics/eta,Ra,Pr,mBC,ampT,l_add_pert,lagpts,n_init,buo_tscheme

      namelist/timecontrol/n_time_steps,dt,time_scheme_type,time_scheme_imp,time_scheme_exp,l_imexrk_started, &
               & dt_coef,dt_max,CFL,l_restart,n_restart,n_restart_point,n_snapshot_point

      namelist/output/tag,n_checkpoint,n_snapshot,n_KE,n_KEspec

      !-- Set default values of control parameters:
      call defaultNamelists

      ! Get the filename of the input file as first argument from the command line

      argument_count = command_argument_count()
      if (argument_count == 0) then
         print *, "'The filename of the input file must be provided as first argument'"
         else
            call get_command_argument(1,input_filename)

            inquire(file = input_filename, exist = nml_exist)

            if (.not. nml_exist) then
               print *, "'Input namelist file not found'"
            end if

            open(newunit=inputHandle,file=trim(input_filename))
            read(inputHandle,nml=grid,iostat=res)
            if ( res /= 0 ) then
               print *, '! No grid namelist found!'
            end if
            close(inputHandle)

            open(newunit=inputHandle,file=trim(input_filename))
            read(inputHandle,nml=physics,iostat=res)
            if ( res /= 0 ) then
               print *, '! No physics namelist found!'
            end if
            close(inputHandle)

            open(newunit=inputHandle,file=trim(input_filename))
            read(inputHandle,nml=timecontrol,iostat=res)
            if ( res /= 0 ) then
               print *, '! No timecontrol namelist found!'
            end if
            close(inputHandle)

            open(newunit=inputHandle,file=trim(input_filename))
            read(inputHandle,nml=output,iostat=res)
            if ( res /= 0 ) then
               print *, '! No output namelist found!'
            end if
            close(inputHandle)
      end if
      Np_max = 3*(Nm_max)      ! Set number of points along the Fourier expansion direction
      rmin=eta/(1.0_dp-eta)    ! Set inner radius using radius ratio 'eta'
      rmax=1.0_dp/(1.0_dp-eta) ! Set outer radius using radius ratio 'eta'
       
      lm=lagpts ! Number of points for interpolation
   end subroutine readNamelists

   subroutine defaultNamelists()

      Nr_max = 20              ! No. of radial grid points (Chebyshev-Gauss-Lobatto points)
      Nm_max = 20              ! No. of azimuthal modes (Fourier modes)
      eta = 0.35_dp            ! Radius ratio
      Ra = 1000.0_dp           ! Rayleigh number
      Pr = 1.0_dp              ! Prandtl number
      mBC = "SF"               ! Mechanical boundary condition: Give "NS" for no-slip and "SF" for stress-free
      ampT = 0.00001_dp        ! Amplitude of the perturbation given on Temperature
      l_add_pert = .FALSE.     ! Give 'True' if a perturbation on the variables is required when restarting
      lagpts = 14              ! Number of lagrange interpolation points for vorticity boundary condition
      buo_tscheme = "IMP"
      n_init = 0               ! Specify type of Gaussian perturbation 
      n_time_steps = 100       ! No. of time steps
      dt = 0.01_dp             ! Time step size 
      time_scheme_type = "IMEXRK" ! Specify type of time integration scheme
      time_scheme_imp = "ARS222" ! State the implicit time scheme
      time_scheme_exp = "expARS222"  ! State the explicit time scheme
      l_imexrk_started = .TRUE.
      dt_coef=2.0_dp           ! Coefficient for dt while calculating min dt 
      dt_max=3.e-5             ! Threshold for max dt  
      CFL=0.9_dp               ! Courant-Friedrichs-Lewy (CFL) condition
      l_restart=.FALSE.        ! Give 'restart=1' if you want to restart from an already stored data 
      n_restart=0              ! If 'restart=1' you have to specify the latest iteration which is saved in the first column of KE.txt
      n_restart_point=0        ! If 'restart=1' you have to specify the latest checkpoint no.
      n_snapshot_point=0       ! If 'restart=1' you have to specify the latest snapshot no. 
      n_checkpoint=500000      ! Checkpoint write interval 
      n_snapshot=300000        ! Snapshot write interval
      n_KE=1                   ! Kinetic energy write interval
      n_KEspec=100000          ! Kinetic energy spectra write interval  

   end subroutine defaultNamelists

end module namelists

 
