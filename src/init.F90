module init
   use double
   use constants, only: one, two, three, half, pi, ii, onefourth
#if FFTW
   use fourier, only: init_fftwplan, destroy_fftwplan, forfft, invfft
   use chebyshev, only: init_chebfftwplan, destroy_chebfftwplan, lagrange, dlagrange, d2lagrange, chebtransform, chebinvtran 
#elif GG
   use fourier, only: forfft, invfft
   use chebyshev, only: lagrange, dlagrange, d2lagrange, chebtransform, chebinvtran 
#endif
   use timeschemes, only: rhs_update_wts_imp, wt_lhs_tscheme_imp, wt_rhs_tscheme_imp, n_order_tscheme_imp, &
                          & n_order_tscheme_exp, n_order_tscheme_max, dt_array, rhs_imp_temp, rhs_exp_temp, &
                          & rhs_imp_vort, rhs_exp_vort, rhs_imp_uphi_bar, rhs_exp_uphi_bar, n_order_tscheme_exp

   implicit none
   private

   real(kind=dp), allocatable, public :: radius(:), r_radius(:), r_radius2(:), phi(:), dr(:)
   real(kind=dp), allocatable, public :: tt(:,:), up(:,:), ur(:,:), psi(:,:), omg(:,:), tcon(:,:), t_an(:,:)
   complex(kind=dp), allocatable, public :: omg_spec(:,:), temp_spec(:,:), upFC(:,:), tFR(:,:), omgFR(:,:), uphi_bar_spec(:)
   complex(kind=dp), allocatable, public :: omgFR_check(:,:)
   complex(kind=dp), allocatable, public :: upFR(:,:), urFR(:,:), psii(:,:) ! FR -> quantities in Fourier-Real space
   complex(kind=dp), allocatable, public :: upFR1(:,:), urFR1(:,:), tFR1(:,:), omgFR1(:,:) ! FR -> quantities in Fourier-Real space
   complex(kind=dp), allocatable, public :: upFR2(:,:), urFR2(:,:), tFR2(:,:), omgFR2(:,:) ! FR -> quantities in Fourier-Real space
   complex(kind=dp), allocatable, public :: upFR3(:,:), urFR3(:,:), tFR3(:,:), omgFR3(:,:) ! FR -> quantities in Fourier-Real space
   complex(kind=dp), allocatable, public :: upFR4(:,:), urFR4(:,:), tFR4(:,:), omgFR4(:,:) ! FR -> quantities in Fourier-Real space
   complex(kind=dp), allocatable, public :: TFC(:)
   complex(kind=dp), allocatable, public :: t2FR(:,:)
   complex(kind=dp), allocatable, public :: upFR_prev(:,:), urFR_prev(:,:)
   real(kind=dp), public :: startmain, finishmain, startNr_maxloop, finishNr_maxloop, startNm_maxloop
   real(kind=dp), public :: finishNm_maxloop, startsteptime, finishsteptime
   real(kind=dp), public :: timeNr_maxloop, timeNm_maxloop, startmatbuild, finishmatbuild, time_matbuild
   real(kind=dp), public :: start_matsolveT, finish_matsolveT, time_matsolveT
   real(kind=dp), public :: start_matsolveW, finish_matsolveW, time_matsolveW
   real(kind=dp), public :: start_tran, finish_tran, time_tran
   real(kind=dp), public :: dt_old, dt_new, tot_time 

   complex(kind=dp), allocatable :: tFRc(:,:), omgFRc(:,:), urFRc(:,:), upFRc(:,:) ! Previously saved data (small 'c' denotes latest timestep from checkpoint)
   complex(kind=dp), allocatable :: tFCc(:,:), omgFCc(:,:), urFCc(:,:), upFCc(:,:) ! Previously saved data (small 'c' denotes latest timestep from checkpoint)
   complex(kind=dp), allocatable :: tFCn(:,:), omgFCn(:,:), urFCn(:,:), upFCn(:,:) ! Previously saved data (small 'n' denotes latest timestep from checkpoint adapted if Nm or Nr different to current )
   complex(kind=dp), allocatable :: rhs_imp_temp_old(:,:,:), rhs_exp_temp_old(:,:,:)
   complex(kind=dp), allocatable :: rhs_imp_vort_old(:,:,:), rhs_exp_vort_old(:,:,:) 
   complex(kind=dp), allocatable :: rhs_imp_uphi_bar_old(:,:), rhs_exp_uphi_bar_old(:,:) 
   real(kind=dp), allocatable :: dt_array_old(:)
   complex(kind=dp), allocatable, public :: tmp_rhs_imp_temp(:,:,:),tmp_rhs_exp_temp(:,:,:)  
   complex(kind=dp), allocatable, public :: tmp_rhs_imp_vort(:,:,:),tmp_rhs_exp_vort(:,:,:) 
   complex(kind=dp), allocatable, public :: tmp_rhs_imp_uphi_bar(:,:),tmp_rhs_exp_uphi_bar(:,:) 
   complex(kind=dp), allocatable, public :: tmp_rhs_buo_term(:,:,:)

   real(kind=dp), allocatable, public :: w_rmin(:)
   real(kind=dp), allocatable, public :: dw_rmin(:)
   real(kind=dp), allocatable, public :: d2w_rmin(:)
   real(kind=dp), allocatable, public :: w_rmax(:)
   real(kind=dp), allocatable, public :: dw_rmax(:)
   real(kind=dp), allocatable, public :: d2w_rmax(:)
#if FFTW
   public :: allocate_restart, deallocate_restart, get_checkpoint_data, &
             & fields_alloc, fields_dealloc, init_grid, init_fields, init_perturbation, &
             & init_all_fftw_plans, destroy_all_fftw_plans, read_checkpoint
#elif GG
   public :: allocate_restart, deallocate_restart, get_checkpoint_data, &
             & fields_alloc, fields_dealloc, init_grid, init_fields, init_perturbation, &
             read_checkpoint
#endif

contains  
   subroutine fields_alloc(Nm_max,Np_max,Nr_max,lm,time_scheme_type,time_scheme_imp) ! main variables allocation

      integer, intent(in) :: lm 
      integer, intent(in) :: Nm_max
      integer, intent(in) :: Np_max   
      integer, intent(in) :: Nr_max
      character(len=100), intent(in) :: time_scheme_type
      character(len=100), intent(in) :: time_scheme_imp

      allocate( radius(Nr_max), r_radius(Nr_max), r_radius2(Nr_max),phi(Np_max), dr(Nr_max) )
      allocate( tt(Np_max,Nr_max), up(Np_max,Nr_max), ur(Np_max,Nr_max), psi(Np_max,Nr_max), &
                & omg(Np_max,Nr_max), tcon(Np_max,Nr_max), t_an(Np_max,Nr_max) )
      allocate( omg_spec(Nm_max+1,Nr_max), temp_spec(Nm_max+1,Nr_max), upFC(Nm_max+1,Nr_max), &
                & tFR(Nm_max+1,Nr_max), omgFR(Nm_max+1,Nr_max), uphi_bar_spec(Nr_max) )
      allocate( omgFR_check(Nm_max+1,Nr_max) )
      allocate( psii(Nm_max+1,Nr_max), upFR(Nm_max+1,Nr_max), urFR(Nm_max+1,Nr_max) )
      allocate( TFC(Nr_max) )
      allocate( t2FR(Nm_max+1,Nr_max), upFR_prev(Nm_max+1,Nr_max), urFR_prev(Nm_max+1,Nr_max))
      allocate( tmp_rhs_imp_temp(n_order_tscheme_imp,Nm_max+1,Nr_max), & 
                & tmp_rhs_exp_temp(n_order_tscheme_exp,Nm_max+1,Nr_max), &
                & tmp_rhs_imp_vort(n_order_tscheme_imp,Nm_max+1,Nr_max), &
                & tmp_rhs_exp_vort(n_order_tscheme_exp,Nm_max+1,Nr_max), &   
                & tmp_rhs_imp_uphi_bar(n_order_tscheme_imp,Nr_max), &
                & tmp_rhs_exp_uphi_bar(n_order_tscheme_exp,Nr_max), &   
                & tmp_rhs_buo_term(n_order_tscheme_imp,Nm_max+1,Nr_max) )
      allocate( upFR1(Nm_max+1,Nr_max), urFR1(Nm_max+1,Nr_max), tFR1(Nm_max+1,Nr_max), omgFR1(Nm_max+1,Nr_max) )
      if (time_scheme_type=='RK' .or. time_scheme_type=='IMEXRK' ) then 
         allocate( w_rmin(lm), dw_rmin(lm), d2w_rmin(lm), w_rmax(lm), dw_rmax(lm), d2w_rmax(lm) )
      end if

   end subroutine fields_alloc

   subroutine fields_dealloc(time_scheme_type) ! main variables deallocation

      character(len=100), intent(in) :: time_scheme_type

      deallocate( tmp_rhs_imp_temp, tmp_rhs_exp_temp, tmp_rhs_imp_vort, &
                  & tmp_rhs_exp_vort, tmp_rhs_exp_uphi_bar, tmp_rhs_imp_uphi_bar, &
                  & tmp_rhs_buo_term ) 
      deallocate( t2FR, upFR_prev, urFR_prev )
      deallocate( TFC )
      deallocate( psii, upFR, urFR )
      deallocate( omg_spec, temp_spec, upFC, tFR, omgFR, uphi_bar_spec )
      deallocate( omgFR_check )
      deallocate( tt, up, ur, psi, omg, tcon,t_an )
      deallocate( radius, r_radius, r_radius2, dr, phi )
      if (time_scheme_type=='RK' .or. time_scheme_type=='IMEXRK') then 
         deallocate( w_rmin, dw_rmin, d2w_rmin, w_rmax, dw_rmax, d2w_rmax )
      end if
   end subroutine fields_dealloc

   subroutine init_grid(Np_max,Nr_max,rmin)

      integer :: i
      integer, intent(in) :: Np_max   
      integer, intent(in) :: Nr_max
      real(kind=dp), intent(in) :: rmin

      do i=1,Nr_max
         radius(i)=(rmin + (0.5_dp*(1.0_dp+cos(real(Nr_max-i,kind=dp)*pi/real(Nr_max-1,kind=dp)))))
         r_radius(i)=1.0_dp/radius(i)
      end do
      
      do i=1,Nr_max-1
         dr(i+1)=radius(i+1)-radius(i)
      end do

      dr(1)=dr(2)

      r_radius2(:)=r_radius(:)*r_radius(:) ! remember to use products instead of power

      do i=1,Np_max
         phi(i)=0.0_dp
      end do

      do i=2,Np_max
         phi(i)=phi(i-1)+2.0_dp*pi/real(Np_max,kind=dp) ! 
      end do
      

   end subroutine init_grid

   subroutine init_fields(Nm_max,Np_max,Nr_max,rmin,rmax,l_restart,Nrestart_point,l_add_pert, &
                          & ampT,lm,time_scheme_type,time_scheme_imp,l_imexrk_started)

      integer :: i,j
      integer, intent(in) :: lm 
      integer, intent(in) :: Nm_max
      integer, intent(in) :: Np_max   
      integer, intent(in) :: Nr_max
      real(kind=dp), intent(in) :: rmin
      real(kind=dp), intent(in) :: rmax
      logical, intent(in) :: l_restart  
      integer, intent(in) :: Nrestart_point
      logical, intent(in) :: l_add_pert
      logical, intent(in) :: l_imexrk_started 
      real(kind=dp), intent(in) :: ampT
      character(len=100), intent(in) :: time_scheme_type
      character(len=100), intent(in) :: time_scheme_imp
     
     !--------- If restarting -----------------------------------------------
      if (l_restart) then   
         print *, "restarting"
         call get_checkpoint_data(Nm_max,Nr_max,Nrestart_point,dt_new,tot_time,dt_array,time_scheme_type, &
                                  & time_scheme_imp,l_imexrk_started) 
         t2FR=tFR
         if (l_add_pert) then
            tFR(:,:) = ampT + tFR(:,:) 
         end if
         temp_spec=tFR
         omg_spec=omgFR
         uphi_bar_spec=upFR(1,:)

         if (time_scheme_type=='RK' .or. time_scheme_type=='IMEXRK') then
            call lagrange(lm,radius(2:lm+1),rmin,w_rmin)
            call dlagrange(lm,radius(2:lm+1),rmin,dw_rmin)
            call d2lagrange(lm,radius(2:lm+1),rmin,d2w_rmin)
            call lagrange(lm,radius(Nr_max-1:Nr_max-lm:-1),rmin,w_rmax)
            call dlagrange(lm,radius(Nr_max-1:Nr_max-lm:-1),rmax,dw_rmax)
            call d2lagrange(lm,radius(Nr_max-1:Nr_max-lm:-1),rmax,d2w_rmax)
         end if

      else
     !--------- If starting from t=0 ----------------------------------------
         ! Define the conducting state 'tcon'
         do j=1,Nr_max
            do i=1,Np_max
               tcon(i,j)=log(radius(j))/log(rmin/rmax) - log(rmax)/log(rmin/rmax)
            end do
         end do  

         ! Define psi, up, ur
         do j=1,Nr_max
            do i=1,Np_max
               psi(i,j)=0.0_dp
            end do
         end do         
         ! Initialize velocity in phi direction

         do j=1,Nr_max
            do i=1,Np_max
               up(i,j)=0.0_dp
            end do
         end do

         ! Initialize velocity in r direction
         do j=1,Nr_max
            do i=1,Np_max
               ur(i,j)=0.0_dp
            end do
         end do

         ! Initialize velocities in (Fourier-Real) space
         do i=1,Nr_max
            call forfft(Nm_max,Np_max,psi(:,i),psii(:,i))
         end do

         do i=1,Nr_max
            call forfft(Nm_max,Np_max,ur(:,i),urFR(:,i))
         end do
  
         do i=1,Nr_max
            call forfft(Nm_max,Np_max,up(:,i),upFR(:,i))
         end do

         ! Initialize vorticity in z­direction and in (Fourier-Real) space

         omg=up 
         omg_spec=upFR
         omgFR=upFR
         upFC=upFR
         upFR_prev=upFR
         urFR_prev=upFR
         uphi_bar_spec(:) = 0.0_dp


         timeNr_maxloop = 0.0_dp
         timeNm_maxloop = 0.0_dp
         time_matbuild = 0.0_dp
         time_matsolveT = 0.0_dp
         time_matsolveW = 0.0_dp
         time_tran = 0.0_dp 

         tot_time=0.0_dp
         
         if (time_scheme_type=='RK' .or. time_scheme_type=='IMEXRK') then
            call dlagrange(lm,radius(2:lm+1),rmin,dw_rmin)
            call d2lagrange(lm,radius(2:lm+1),rmin,d2w_rmin)
            call dlagrange(lm,radius(Nr_max-1:Nr_max-lm:-1),rmax,dw_rmax)
            call d2lagrange(lm,radius(Nr_max-1:Nr_max-lm:-1),rmax,d2w_rmax)
         end if

      end if
   end subroutine init_fields

   subroutine init_perturbation(Nm_max,Np_max,Nr_max,ampT,l_restart,rmin,rmax,n_init)

      !-------------------------------------------------------------
      integer :: i,j
      integer, intent(in) :: Nm_max
      integer, intent(in) :: Np_max   
      integer, intent(in) :: Nr_max
      integer, intent(in) :: n_init
      real(kind=dp), intent(in) :: ampT
      logical, intent(in) :: l_restart
      real(kind=dp), intent(in) :: rmin, rmax
      real(kind=dp) :: sigma_p, sigma_r
      real(kind=dp) :: rc, c1, c2, L, gasp_c, ph, phi0
      real(kind=dp) :: c_r 
      real(kind=dp) :: gasp(Nr_max)
      real(kind=dp) :: phi_func(Np_max)
      integer :: Np,Nr,ir

      if (.not. l_restart) then
         print *, "NO RESTART"
         if (n_init==0) then ! If n_init=0, then initialize with a Gaussian perturbation (from Gaspari 1999 paper) on top of conducting state
            sigma_r = 0.1_dp/sqrt(two) 
            sigma_p = 0.2
            !sigma_r = 0.25 
            !sigma_p = 0.25
            rc = half*(rmax+rmin) 
            !-- Find the closest point to the middle radius 
            ir=minloc(abs(radius-rc),1) 
            !-- Overwrite rc to be on the grid 
            rc = radius(ir) 
            c1 = half*(radius(Nr_max)-rc) 
            c2 = half*(rc-radius(1))
            L = half * sigma_r 

            !-- Piecewise-definition of a compact-support Gaussian-like profile 
            !-- From Eq. (4.7) from Gaspari et al. (1999) 
            !-- This ensures that the BCs will always be fulfilled 
            do Nr=1,Nr_max 
               if ( Nr == ir ) then ! Middle point 
                  gasp(Nr)=one 
               else ! Middle and Edge parts 
                  if ( radius(Nr) < rc-c2 ) then 
                     call gausslike_compact_edge(rc-radius(Nr),c2,L,gasp(Nr))
                  else if ( radius(Nr) >= rc-c2 .and. radius(Nr) < rc ) then 
                     call gausslike_compact_middle(rc-radius(Nr),c2,L,gasp(Nr))
                  else if ( radius(Nr) > rc .and. radius(Nr) <= rc+c1 ) then 
                     call gausslike_compact_middle(radius(Nr)-rc,c1,L,gasp(Nr))
                  else if ( radius(Nr) > rc+c1 ) then 
                     call gausslike_compact_edge(radius(Nr)-rc,c1,L,gasp(Nr))
                  end if 
               end if 
            end do 

            !-- Normalisation of the two branches by the middle value 
            do Nr=1,Nr_max 
               if ( Nr < ir ) then
                  call gausslike_compact_center(c1,L,gasp_c) 
                  gasp(Nr) = gasp(Nr)/gasp_c 
               else if ( Nr > ir ) then 
                  call gausslike_compact_center(c2,L,gasp_c) 
                  gasp(Nr) = gasp(Nr)/gasp_c
               end if 
            end do

            !-- Now finally define the bubble 
            phi0 = pi
            do Nr=1,Nr_max
               c_r = ampT*gasp(Nr) 
               do Np=1,Np_max 
                  ph = (Np-1)*two*pi/(Np_max) 
                  phi_func(Np)=c_r*exp(-(ph-phi0)**2/(sigma_p)**2) 
                  tt(Np,Nr)=c_r*exp(-(ph-phi0)**2/(sigma_p)**2) + tcon(Np,Nr) ! Add perturbation + conducting state 'tcon' 
               end do 
               call forfft(Nm_max,Np_max,tt(:,Nr),tFR(:,Nr)) 
            end do

         elseif (n_init<0) then ! If n_init<0, initialize with the standard Gaussian profile perturbation on top of contucting state

            sigma_r = 0.25_dp
            sigma_p = 0.25_dp

            do j=1,Nr_max
               do i=1,Np_max
                  if (phi(i)>1.0_dp .and. phi(i)<3.0_dp) then
                     tt(i,j) = -1.0_dp*(rmax-radius(j))*(rmin-radius(j))*ampT* &
                               & (1.0_dp/(2.0*pi*sigma_p*sigma_r))*exp(-1.0_dp*((radius(j)-(rmax-(rmax-rmin)/2.0_dp))* & 
                               & (radius(j)-(rmax-(rmax-rmin)/2.0_dp))/ (sigma_r*sigma_r)+(phi(i)-2.0_dp)* &
                               & (phi(i)-2.0_dp)/(sigma_p*sigma_p))) + tcon(i,j) 
                  else
                     tt(i,j)=tcon(i,j) 
                  end if
               end do
            end do

            do i=1,Nr_max
               call forfft(Nm_max,Np_max,tt(:,i),tFR(:,i))
            end do

         elseif (n_init>0) then ! If n_init>0, intialize with perturbation on the mode='n_init' on top of conducting state 

            do j=1,Nr_max
               do i=1,Np_max
                  tt(i,j) = ampT*(rmax-radius(j))*(rmin-radius(j))*cos(real(n_init,kind=dp)*phi(i)) + tcon(i,j) 
               end do
            end do

            do i=1,Nr_max
               call forfft(Nm_max,Np_max,tt(:,i),tFR(:,i))
            end do
         
         end if

         print *, "maxval of T on the outerwall is", maxval(tt) ! Print to check if BC is properly satisfied after adding perturbation 
         print *, "minval of T on the outerwall is", minval(tt) ! Print to check if BC is properly satisfied after adding perturbation 

         !--------------
         temp_spec=tFR
         omg_spec=omgFR
         upFC=upFR
         !--------------

      end if
   end subroutine init_perturbation

#if FFTW
   subroutine init_all_fftw_plans(Np_max,Nr_max)
      
      integer, intent(in) :: Np_max,Nr_max
      !--------------------------------------------------------------------------------------
      call init_fftwplan(Np_max)  ! Make the plans for USEGG Forward and Backward Transforms 
      call init_chebfftwplan(Nr_max) ! Make the plans for Fast Chebyshev Transforms 
      !--------------------------------------------------------------------------------------
 
   end subroutine init_all_fftw_plans

   subroutine destroy_all_fftw_plans()

      !------------------------------------------------------------------------------------
      call destroy_chebfftwplan() ! Destroy the plans for Fast Chebyshev transform
      call destroy_fftwplan()     ! Destroy the plans for USEGG forward and backward transforms
      !------------------------------------------------------------------------------------

   end subroutine destroy_all_fftw_plans
#endif

   subroutine allocate_restart(Nm_max,Nr_max) !--------Allocate----------

      integer, intent(in) :: Nm_max
      integer, intent(in) :: Nr_max 

      allocate( tFRc(Nm_max+1,Nr_max), omgFRc(Nm_max+1,Nr_max), urFRc(Nm_max+1,Nr_max), &
                & upFRc(Nm_max+1,Nr_max) )
      allocate( tFCc(Nm_max+1,Nr_max), omgFCc(Nm_max+1,Nr_max), urFCc(Nm_max+1,Nr_max), &
                & upFCc(Nm_max+1,Nr_max) )
      allocate( tFCn(Nm_max+1,Nr_max), omgFCn(Nm_max+1,Nr_max), urFCn(Nm_max+1,Nr_max), &
                & upFCn(Nm_max+1,Nr_max) )
      allocate( rhs_imp_temp_old(n_order_tscheme_imp,Nm_max+1,Nr_max),rhs_exp_temp_old( &
                & n_order_tscheme_exp, Nm_max+1, Nr_max),rhs_imp_vort_old(n_order_tscheme_imp, &
                & Nm_max+1, Nr_max),rhs_exp_vort_old(n_order_tscheme_exp,Nm_max+1,Nr_max), &
                & rhs_imp_uphi_bar_old(n_order_tscheme_imp,Nr_max), rhs_exp_uphi_bar_old( &
                & n_order_tscheme_exp, Nr_max) )  
      allocate( dt_array_old(n_order_tscheme_max) )

   end subroutine allocate_restart

   subroutine deallocate_restart() !----------Deallocate--------

      deallocate( dt_array_old )
      deallocate( rhs_imp_temp_old,rhs_exp_temp_old,rhs_imp_vort_old,rhs_exp_vort_old )    
      deallocate( rhs_imp_uphi_bar_old,rhs_exp_uphi_bar_old )
      deallocate( tFCn, omgFCn, urFCn, upFCn )
      deallocate( tFCc, omgFCc, urFCc, upFCc )
      deallocate( tFRc, omgFRc, urFRc, upFRc )

   end subroutine deallocate_restart

   subroutine get_checkpoint_data(Nm_max,Nr_max,Nrestart_point,dt_new,tot_time,dt_array,time_scheme_type, &
                                  & time_scheme_imp,l_imexrk_started)

      integer, intent(in) :: Nm_max
      integer, intent(in) :: Nr_max 
      integer, intent(in) :: Nrestart_point
      logical, intent(in) :: l_imexrk_started 
      real(kind=dp), intent(out) :: dt_new
      real(kind=dp), intent(out) :: tot_time
      integer :: Nm_max_old, Nr_max_old
      integer :: n_order_tscheme_imp_old, n_order_tscheme_exp_old, n_order_tscheme_max_old 
      integer :: i,j,Nm,Nr
      
      real(kind=dp), intent(out) :: dt_array(n_order_tscheme_max)
      character(len=100), intent(in) :: time_scheme_type
      character(len=100), intent(in) :: time_scheme_imp

      if ((time_scheme_imp=='CN' .or. time_scheme_imp=='BDF2') .and. l_imexrk_started) then
         call read_checkpoint(dt_new,tot_time,Nm_max_old,Nr_max_old,n_order_tscheme_imp_old, &           
                              & n_order_tscheme_exp_old,Nrestart_point-1,tFR1,omgFR1,urFR1,upFR1, &    
                              & rhs_imp_temp_old,rhs_exp_temp_old,rhs_imp_vort_old,rhs_exp_vort_old, &
                              & rhs_imp_uphi_bar_old,rhs_exp_uphi_bar_old,dt_array_old, &
                              & n_order_tscheme_max_old) 

         dt_array(2)=dt_array_old(1)
         
         call read_checkpoint(dt_new,tot_time,Nm_max_old,Nr_max_old,n_order_tscheme_imp_old, &           
                              & n_order_tscheme_exp_old,Nrestart_point,tFR2,omgFR2,urFR2,upFR2, &    
                              & rhs_imp_temp_old,rhs_exp_temp_old,rhs_imp_vort_old,rhs_exp_vort_old, &
                              & rhs_imp_uphi_bar_old,rhs_exp_uphi_bar_old,dt_array_old, &
                              & n_order_tscheme_max_old)

         dt_array(1)=dt_array_old(1)
         
         if (Nm_max_old==Nm_max .and. Nr_max_old==Nr_max) then
            tFR=tFR2
            omgFR=omgFR2
            urFR=urFR2
            upFR=upFR2
         end if  

      elseif (time_scheme_imp=='BDF3' .and. l_imexrk_started) then
         call read_checkpoint(dt_new,tot_time,Nm_max_old,Nr_max_old,n_order_tscheme_imp_old, &           
                              & n_order_tscheme_exp_old,Nrestart_point-2,tFR1,omgFR1,urFR1,upFR1, &    
                              & rhs_imp_temp_old,rhs_exp_temp_old,rhs_imp_vort_old,rhs_exp_vort_old, &
                              & rhs_imp_uphi_bar_old,rhs_exp_uphi_bar_old,dt_array_old, &
                              & n_order_tscheme_max_old) 

         dt_array(3)=dt_array_old(1)
         
         call read_checkpoint(dt_new,tot_time,Nm_max_old,Nr_max_old,n_order_tscheme_imp_old, &           
                              & n_order_tscheme_exp_old,Nrestart_point-1,tFR2,omgFR2,urFR2,upFR2, &    
                              & rhs_imp_temp_old,rhs_exp_temp_old,rhs_imp_vort_old,rhs_exp_vort_old, &
                              & rhs_imp_uphi_bar_old,rhs_exp_uphi_bar_old,dt_array_old, &
                              & n_order_tscheme_max_old) 

         dt_array(2)=dt_array_old(1)
         
         call read_checkpoint(dt_new,tot_time,Nm_max_old,Nr_max_old,n_order_tscheme_imp_old, &           
                              & n_order_tscheme_exp_old,Nrestart_point,tFR3,omgFR3,urFR3,upFR3, &    
                              & rhs_imp_temp_old,rhs_exp_temp_old,rhs_imp_vort_old,rhs_exp_vort_old, &
                              & rhs_imp_uphi_bar_old,rhs_exp_uphi_bar_old,dt_array_old, &
                              & n_order_tscheme_max_old)

         dt_array(1)=dt_array_old(1)
         
         if (Nm_max_old==Nm_max .and. Nr_max_old==Nr_max) then
            tFR=tFR3
            omgFR=omgFR3
            urFR=urFR3
            upFR=upFR3
         end if  

      elseif (time_scheme_imp=='BDF4' .and. l_imexrk_started) then
         
         call read_checkpoint(dt_new,tot_time,Nm_max_old,Nr_max_old,n_order_tscheme_imp_old, &           
                              & n_order_tscheme_exp_old,Nrestart_point-3,tFR1,omgFR1,urFR1,upFR1, &    
                              & rhs_imp_temp_old,rhs_exp_temp_old,rhs_imp_vort_old,rhs_exp_vort_old, &
                              & rhs_imp_uphi_bar_old,rhs_exp_uphi_bar_old,dt_array_old, &
                              & n_order_tscheme_max_old) 

         dt_array(4)=dt_array_old(1)
         
         call read_checkpoint(dt_new,tot_time,Nm_max_old,Nr_max_old,n_order_tscheme_imp_old, &           
                              & n_order_tscheme_exp_old,Nrestart_point-2,tFR2,omgFR2,urFR2,upFR2, &    
                              & rhs_imp_temp_old,rhs_exp_temp_old,rhs_imp_vort_old,rhs_exp_vort_old, &
                              & rhs_imp_uphi_bar_old,rhs_exp_uphi_bar_old,dt_array_old, &
                              & n_order_tscheme_max_old) 

         dt_array(3)=dt_array_old(1)
         
         call read_checkpoint(dt_new,tot_time,Nm_max_old,Nr_max_old,n_order_tscheme_imp_old, &           
                              & n_order_tscheme_exp_old,Nrestart_point-1,tFR3,omgFR3,urFR3,upFR3, &    
                              & rhs_imp_temp_old,rhs_exp_temp_old,rhs_imp_vort_old,rhs_exp_vort_old, &
                              & rhs_imp_uphi_bar_old,rhs_exp_uphi_bar_old,dt_array_old, &
                              & n_order_tscheme_max_old)

         dt_array(2)=dt_array_old(1)

         call read_checkpoint(dt_new,tot_time,Nm_max_old,Nr_max_old,n_order_tscheme_imp_old, &           
                              & n_order_tscheme_exp_old,Nrestart_point,tFR4,omgFR4,urFR4,upFR4, &    
                              & rhs_imp_temp_old,rhs_exp_temp_old,rhs_imp_vort_old,rhs_exp_vort_old, &
                              & rhs_imp_uphi_bar_old,rhs_exp_uphi_bar_old,dt_array_old, &
                              & n_order_tscheme_max_old)

         dt_array(1)=dt_array_old(1)
         
         if (Nm_max_old==Nm_max .and. Nr_max_old==Nr_max) then
            tFR=tFR4
            omgFR=omgFR4
            urFR=urFR4
            upFR=upFR4
         end if  

      elseif(time_scheme_type=='IMEXRK') then
         !---- Get previously saved data ----
         call read_checkpoint(dt_new,tot_time,Nm_max_old,Nr_max_old,n_order_tscheme_imp_old, &           
                                 & n_order_tscheme_exp_old,Nrestart_point,tFRc,omgFRc,urFRc,upFRc, &    
                                 & rhs_imp_temp_old,rhs_exp_temp_old,rhs_imp_vort_old,rhs_exp_vort_old, &
                                 & rhs_imp_uphi_bar_old,rhs_exp_uphi_bar_old,dt_array_old, &
                                 & n_order_tscheme_max_old) 
!----    If order of the schemes are different from before (dt Array) ---------------- 
         if (n_order_tscheme_max > n_order_tscheme_max_old) then
            do i=1,n_order_tscheme_max_old
               dt_array(i)=dt_array_old(i)
            end do
            dt_array(n_order_tscheme_max_old+1:n_order_tscheme_max)=dt_array(n_order_tscheme_max_old)

         else if (n_order_tscheme_max < n_order_tscheme_max_old) then
            do i=1,n_order_tscheme_max
               dt_array(i)=dt_array_old(i)
            end do

!----    If order of the schemes are same as before (dt Array) ---------------- 
         else if (n_order_tscheme_max == n_order_tscheme_max_old) then
            dt_array=dt_array_old

         end if
         !tot_time=tot_time+dt_new

         if (Nm_max_old==Nm_max .and. Nr_max_old==Nr_max) then
            tFR=tFRc
            omgFR=omgFRc
            urFR=urFRc
            upFR=upFRc
         end if  
          
      else 
         !---- Get previously saved data ----
         call read_checkpoint(dt_new,tot_time,Nm_max_old,Nr_max_old,n_order_tscheme_imp_old, &           
                                 & n_order_tscheme_exp_old,Nrestart_point,tFRc,omgFRc,urFRc,upFRc, &    
                                 & rhs_imp_temp_old,rhs_exp_temp_old,rhs_imp_vort_old,rhs_exp_vort_old, &
                                 & rhs_imp_uphi_bar_old,rhs_exp_uphi_bar_old,dt_array_old, &
                                 & n_order_tscheme_max_old) 
!----    If order of the schemes are different from before (dt Array) ---------------- 
         if (n_order_tscheme_max > n_order_tscheme_max_old) then
            do i=1,n_order_tscheme_max_old
               dt_array(i)=dt_array_old(i)
            end do
            dt_array(n_order_tscheme_max_old+1:n_order_tscheme_max)=dt_array(n_order_tscheme_max_old)

         else if (n_order_tscheme_max < n_order_tscheme_max_old) then
            do i=1,n_order_tscheme_max
               dt_array(i)=dt_array_old(i)
            end do

!----    If order of the schemes are same as before (dt Array) ---------------- 
         else if (n_order_tscheme_max == n_order_tscheme_max_old) then
            dt_array=dt_array_old

         end if
! ----   ----------------------------------------------------------------------- 

!----    If order of the schemes are different from before (Time scheme implicit) ----- 
         if (n_order_tscheme_imp > n_order_tscheme_imp_old) then
            do i=1,n_order_tscheme_imp_old
               do Nm=0,Nm_max_old
                  do Nr=1,Nr_max_old
                     rhs_imp_temp(i,Nm+1,Nr)=rhs_imp_temp_old(i,Nm+1,Nr)
                     rhs_imp_vort(i,Nm+1,Nr)=rhs_imp_vort_old(i,Nm+1,Nr)
                     rhs_imp_uphi_bar(i,Nr)=rhs_imp_uphi_bar_old(i,Nr)
                  end do
               end do 
            end do
            do j=n_order_tscheme_imp_old+1,n_order_tscheme_imp
               do Nm=0,Nm_max_old
                  do Nr=1,Nr_max_old
                     rhs_imp_temp(j,Nm+1,Nr)=rhs_imp_temp(n_order_tscheme_imp_old,Nm+1,Nr)
                     rhs_imp_vort(j,Nm+1,Nr)=rhs_imp_vort(n_order_tscheme_imp_old,Nm+1,Nr)
                     rhs_imp_uphi_bar(j,Nr)=rhs_imp_uphi_bar(n_order_tscheme_imp_old,Nr)
                  end do
               end do 
            end do

         else if (n_order_tscheme_imp < n_order_tscheme_imp_old) then
            do i=1,n_order_tscheme_imp
               do Nm=0,Nm_max_old
                  do Nr=1,Nr_max_old
                     rhs_imp_temp(i,Nm+1,Nr)=rhs_imp_temp_old(i,Nm+1,Nr)
                     rhs_imp_vort(i,Nm+1,Nr)=rhs_imp_vort_old(i,Nm+1,Nr)
                     rhs_imp_uphi_bar(i,Nr)=rhs_imp_uphi_bar_old(i,Nr)
                  end do
               end do 
            end do

!----    If order of the schemes are same as before (Time scheme implicit) ------------ 
         else if (n_order_tscheme_imp == n_order_tscheme_imp_old) then
            do i=1,n_order_tscheme_imp
               do Nm=0,Nm_max_old
                  do Nr=1,Nr_max_old
                     rhs_imp_temp(i,Nm+1,Nr)=rhs_imp_temp_old(i,Nm+1,Nr)
                     rhs_imp_vort(i,Nm+1,Nr)=rhs_imp_vort_old(i,Nm+1,Nr)
                     rhs_imp_uphi_bar(i,Nr)=rhs_imp_uphi_bar_old(i,Nr)
                  end do
               end do 
            end do
         end if
               
!----    If order of the schemes are different from before (Time scheme explicit) ----- 
         if (n_order_tscheme_exp > n_order_tscheme_exp_old) then
            do i=1,n_order_tscheme_exp_old
               do Nm=0,Nm_max_old
                  do Nr=1,Nr_max_old
                     rhs_exp_temp(i,Nm+1,Nr)=rhs_exp_temp_old(i,Nm+1,Nr)
                     rhs_exp_vort(i,Nm+1,Nr)=rhs_exp_vort_old(i,Nm+1,Nr)
                     rhs_exp_uphi_bar(i,Nr)=rhs_exp_uphi_bar_old(i,Nr)
                  end do
               end do 
            end do
            do j=n_order_tscheme_exp_old+1,n_order_tscheme_exp
               do Nm=0,Nm_max_old
                  do Nr=1,Nr_max_old
                     rhs_exp_temp(j,Nm+1,Nr)=rhs_imp_temp(n_order_tscheme_exp_old,Nm+1,Nr)
                     rhs_exp_vort(j,Nm+1,Nr)=rhs_imp_vort(n_order_tscheme_exp_old,Nm+1,Nr)
                     rhs_exp_uphi_bar(j,Nr)=rhs_imp_uphi_bar(n_order_tscheme_exp_old,Nr)
                  end do
               end do 
            end do
         
         else if (n_order_tscheme_exp < n_order_tscheme_exp_old) then
            do i=1,n_order_tscheme_exp
               do Nm=0,Nm_max_old
                  do Nr=1,Nr_max_old
                     rhs_exp_temp(i,Nm+1,Nr)=rhs_exp_temp_old(i,Nm+1,Nr)
                     rhs_exp_vort(i,Nm+1,Nr)=rhs_exp_vort_old(i,Nm+1,Nr)
                     rhs_exp_uphi_bar(i,Nr)=rhs_exp_uphi_bar_old(i,Nr)
                  end do
               end do 
            end do

!----    If order of the schemes are same as before (Time scheme explicit) ------------ 
         else if (n_order_tscheme_exp == n_order_tscheme_exp_old) then
            do i=1,n_order_tscheme_exp
               do Nm=0,Nm_max_old
                  do Nr=1,Nr_max_old
                     rhs_exp_temp(i,Nm+1,Nr)=rhs_exp_temp_old(i,Nm+1,Nr)
                     rhs_exp_vort(i,Nm+1,Nr)=rhs_exp_vort_old(i,Nm+1,Nr)
                     rhs_exp_uphi_bar(i,Nr)=rhs_exp_uphi_bar_old(i,Nr)
                  end do
               end do 
            end do
         end if
!-----   -------------------------------------------------------------------------------

         do Nm=0,Nm_max_old 
            call chebtransform(Nr_max,tFRc(Nm+1,:),tFCc(Nm+1,:))
            call chebtransform(Nr_max,omgFRc(Nm+1,:),omgFCc(Nm+1,:))
            call chebtransform(Nr_max,urFRc(Nm+1,:),urFCc(Nm+1,:))
            call chebtransform(Nr_max,upFRc(Nm+1,:),upFCc(Nm+1,:))
         end do
            
         if (Nm_max_old==Nm_max .and. Nr_max_old==Nr_max) then
            tFR=tFRc
            omgFR=omgFRc
            urFR=urFRc
            upFR=upFRc
         else if (Nm_max>Nm_max_old .and. Nr_max==Nr_max_old) then
            do Nm=0,Nm_max_old
               do Nr=1,Nr_max_old
                  tFCn(Nm+1,Nr)=tFCc(Nm+1,Nr) 
                  omgFCn(Nm+1,Nr)=omgFCc(Nm+1,Nr) 
                  urFCn(Nm+1,Nr)=urFCc(Nm+1,Nr) 
                  upFCn(Nm+1,Nr)=upFCc(Nm+1,Nr) 
               end do
            end do
            do Nm=Nm_max_old+1,Nm_max
               do Nr=1,Nr_max_old
                  tFCn(Nm+1,Nr)=0.0_dp
                  omgFCn(Nm+1,Nr)=0.0_dp
                  urFCn(Nm+1,Nr)=0.0_dp
                  upFCn(Nm+1,Nr)=0.0_dp
               end do
            end do
            do Nm=0,Nm_max
               call chebinvtran(Nr_max,tFCn(Nm+1,:),tFR(Nm+1,:))
               call chebinvtran(Nr_max,omgFCn(Nm+1,:),omgFR(Nm+1,:))
               call chebinvtran(Nr_max,urFCn(Nm+1,:),urFR(Nm+1,:))
               call chebinvtran(Nr_max,upFCn(Nm+1,:),upFR(Nm+1,:))
            end do
         else if (Nm_max==Nm_max_old .and. Nr_max>Nr_max_old) then
            do Nm=0,Nm_max_old
               do Nr=1,Nr_max_old
                  tFCn(Nm+1,Nr)=tFCc(Nm+1,Nr) 
                  omgFCn(Nm+1,Nr)=omgFCc(Nm+1,Nr) 
                  urFCn(Nm+1,Nr)=urFCc(Nm+1,Nr) 
                  upFCn(Nm+1,Nr)=upFCc(Nm+1,Nr) 
               end do
            end do
            do Nm=0,Nm_max_old
               do Nr=Nr_max_old+1,Nr_max
                  tFCn(Nm+1,Nr)=0.0_dp
                  omgFCn(Nm+1,Nr)=0.0_dp
                  urFCn(Nm+1,Nr)=0.0_dp
                  upFCn(Nm+1,Nr)=0.0_dp
               end do
            end do
            do Nm=0,Nm_max
               call chebinvtran(Nr_max,tFCn(Nm+1,:),tFR(Nm+1,:))
               call chebinvtran(Nr_max,omgFCn(Nm+1,:),omgFR(Nm+1,:))
               call chebinvtran(Nr_max,urFCn(Nm+1,:),urFR(Nm+1,:))
               call chebinvtran(Nr_max,upFCn(Nm+1,:),upFR(Nm+1,:))
            end do
         end if

          
         !  ----- Create temporary rhs from previously stored rhs -----------  
         tmp_rhs_imp_temp=rhs_imp_temp 
         tmp_rhs_exp_temp=rhs_exp_temp
         tmp_rhs_imp_vort=rhs_imp_vort
         tmp_rhs_exp_vort=rhs_exp_vort
         tmp_rhs_imp_uphi_bar=rhs_imp_uphi_bar
         tmp_rhs_exp_uphi_bar=rhs_exp_uphi_bar

         !  -----------------------------------------------------------------
      end if
   end subroutine get_checkpoint_data

   subroutine read_checkpoint(dt_new,tot_time,Nm_max,Nr_max,n_order_tscheme_imp_old, &
                                 & n_order_tscheme_exp_old,Nrestart_point,tFRn,omgFRn,urFRn,upFRn, &
                                 & rhs_imp_temp_old,rhs_exp_temp_old,rhs_imp_vort_old,rhs_exp_vort_old,rhs_imp_uphi_bar_old, &
                                 & rhs_exp_uphi_bar_old, dt_array_old, &
                                 & n_order_tscheme_max_old)  

      real(kind=dp), intent(out) :: dt_new
      real(kind=dp), intent(out) :: tot_time
      integer, intent(out) :: Nm_max   
      integer, intent(out) :: Nr_max 
      integer, intent(out) :: n_order_tscheme_imp_old, n_order_tscheme_exp_old 
      integer, intent(out) :: n_order_tscheme_max_old 
      integer, intent(in) :: Nrestart_point
      complex(kind=dp), allocatable, intent(out) :: tFRn(:,:)
      complex(kind=dp), allocatable, intent(out) :: omgFRn(:,:)
      complex(kind=dp), allocatable, intent(out) :: urFRn(:,:)
      complex(kind=dp), allocatable, intent(out) :: upFRn(:,:)
      complex(kind=dp), allocatable, intent(out) :: rhs_imp_temp_old(:,:,:)
      complex(kind=dp), allocatable, intent(out) :: rhs_exp_temp_old(:,:,:)
      complex(kind=dp), allocatable, intent(out) :: rhs_imp_vort_old(:,:,:)
      complex(kind=dp), allocatable, intent(out) :: rhs_exp_vort_old(:,:,:)
      complex(kind=dp), allocatable, intent(out) :: rhs_imp_uphi_bar_old(:,:)
      complex(kind=dp), allocatable, intent(out) :: rhs_exp_uphi_bar_old(:,:)
      real(kind=dp), allocatable, intent(out) :: dt_array_old(:)

      character(len=72) :: datafile1
      integer :: inunit
      !----------------------------- Read latest files -------------------------------
      print *, Nrestart_point
      write(datafile1,'("checkpoint_",I5.5)') Nrestart_point
      open(newunit=inunit, status='old',file=datafile1,form='unformatted')
      read(inunit) dt_new,tot_time,Nm_max,Nr_max,n_order_tscheme_imp_old, &
                   & n_order_tscheme_exp_old, n_order_tscheme_max_old
      allocate(dt_array_old(n_order_tscheme_max_old))
      allocate(tFRn(Nm_max+1,Nr_max))
      allocate(omgFRn(Nm_max+1,Nr_max))
      allocate(urFRn(Nm_max+1,Nr_max))
      allocate(upFRn(Nm_max+1,Nr_max))
      allocate(rhs_imp_temp_old(n_order_tscheme_imp_old,Nm_max+1,Nr_max)) 
      allocate(rhs_exp_temp_old(n_order_tscheme_exp_old,Nm_max+1,Nr_max)) 
      allocate(rhs_imp_vort_old(n_order_tscheme_imp_old,Nm_max+1,Nr_max)) 
      allocate(rhs_exp_vort_old(n_order_tscheme_exp_old,Nm_max+1,Nr_max)) 
      allocate(rhs_imp_uphi_bar_old(n_order_tscheme_imp_old,Nr_max)) 
      allocate(rhs_exp_uphi_bar_old(n_order_tscheme_exp_old,Nr_max)) 

      !allocate(rhs_imp_temp_old(1,Nm_max+1,Nr_max)) 
      !allocate(rhs_exp_temp_old(1,Nm_max+1,Nr_max)) 
      !allocate(rhs_imp_vort_old(1,Nm_max+1,Nr_max)) 
      !allocate(rhs_exp_vort_old(1,Nm_max+1,Nr_max)) 
      !allocate(rhs_imp_uphi_bar_old(1,Nr_max)) 
      !allocate(rhs_exp_uphi_bar_old(1,Nr_max)) 

      read(inunit) dt_array_old
      read(inunit) tFRn
      read(inunit) omgFRn
      read(inunit) urFRn
      read(inunit) upFRn
      
      read(inunit) rhs_imp_temp_old   
      read(inunit) rhs_exp_temp_old   
      read(inunit) rhs_imp_vort_old   
      read(inunit) rhs_exp_vort_old   
      read(inunit) rhs_imp_uphi_bar_old   
      read(inunit) rhs_exp_uphi_bar_old   

      close(inunit)

   end subroutine read_checkpoint

   subroutine gausslike_compact_middle(r,c,L,gasp1)
   !
   ! This defines the first part of a compact-support Gaussian-like 
   ! profile. This is adapted from a Gaspari & Cohn, QJRMS, 1999, Eq. 4.7b 
   ! 
      !-- Input variables 
      real(kind=dp), intent(in) :: r 
      real(kind=dp), intent(in) :: c 
      real(kind=dp), intent(in) :: L

      real(kind=dp), intent(out) :: gasp1 

      gasp1 = pi/three*L*r*(r+three*L)*exp(-r/L)+(two*pi*L*L*(c+L)**2)/r* & 
              & exp(-two*c/L)*(one-(one-r/(c+L))*exp(r/L))+pi*L*L*L*(exp(-r/L)- & 
              & exp((r-two*c)/L)) 

   end subroutine gausslike_compact_middle 

   subroutine gausslike_compact_edge(r,c,L,gasp2)
   !
   ! This defines the second part of a compact-support Gaussian-like 
   ! profile. This is adapted from a Gaspari & Cohn, QJRMS, 1999, Eq. 4.7a 
   ! 
      !-- Input variables 
      real(kind=dp), intent(in) :: r 
      real(kind=dp), intent(in) :: c 
      real(kind=dp), intent(in) :: L

      real(kind=dp), intent(out) :: gasp2 

      gasp2 = two*pi*L/r*exp(-r/L)*(half*(r*(r+L)*(two*c-r))+one/three*( & 
              & (r-c)**3-c**3)-L*(c+L)*(r-c+L)+L*(c+L)**2*exp((r-two*c)/L)) 

   end subroutine gausslike_compact_edge 

   subroutine gausslike_compact_center(c,L,gasp3)
   ! 
   ! This defines the central point of a compact-support Gaussian-like 
   ! profile. This is adapted from a Gaspari & Cohn, QJRMS, 1999, Eq. 4.7c 
   ! 
      !-- Input variables 
      real(kind=dp), intent(in) :: c 
      real(kind=dp), intent(in) :: L 

      real(kind=dp), intent(out) :: gasp3 

      gasp3 = pi*L*L*L*(one-exp(-two*c/L))-two*pi*c*L*(c+L)*exp(-two*c/L) 

   end subroutine gausslike_compact_center 

end module init

 
 
