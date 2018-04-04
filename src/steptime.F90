module steptime

   use double
#if FFTW
   use chebyshev, only: init_chebfftwplan, destroy_chebfftwplan, chebtransform, chebinvtran
#elif GG
   use chebyshev, only: chebtransform, chebinvtran
#endif
   use init, only: TFC, t2FR, upFR, urFR, tFR, omgFR, startNm_maxloop, finishNm_maxloop, & 
                   & startNr_maxloop, finishNr_maxloop, timeNm_maxloop, timeNr_maxloop, startsteptime, &
                   & dt_new, tot_time, tmp_rhs_imp_temp, tmp_rhs_exp_temp, tmp_rhs_imp_vort, tmp_rhs_exp_vort, &
                   & tmp_rhs_imp_uphi_bar, tmp_rhs_exp_uphi_bar, tmp_rhs_buo_term, upFR_prev, urFR_prev, up, &
                   & urFR3,upFR3, tFR3, omgFR3, urFR2,upFR2, tFR2, omgFR2, urFR1,upFR1, tFR1, omgFR1  
   use nonlin, only: Nr_maxloop, dtval_p, dtval_r
   use output, only: init_output, calculate_spectra, final_output, writeke_spectral, &
                     & store_checkpoint, store_snapshot, writeke_physical
   use fourierloop, only: Nm_maxloop
   use timeschemes, only: rhs_update_wts_imp, wt_lhs_tscheme_imp, wt_rhs_tscheme_imp, n_order_tscheme_imp, &
                          & n_order_tscheme_exp, n_order_tscheme_max, dt_array, rhs_imp_temp, rhs_exp_temp, &
                          & rhs_imp_vort, rhs_exp_vort, rhs_buo_term, rhs_imp_uphi_bar, rhs_exp_uphi_bar, n_order_tscheme_exp
   use rhs_create, only: rhs_temp_from_restart, rhs_vort_from_restart, rhs_uphibar_from_restart 

   implicit none

   private

   complex(kind=dp), allocatable, public :: uphi_temp_FR(:,:), ur_temp_FR(:,:)
   complex(kind=dp), allocatable, public :: uphi_omg_FR(:,:), ur_omg_FR(:,:)

   public :: allocate_steptime_imex, deallocate_steptime_imex, timeloop_imex

contains

   subroutine allocate_steptime_imex(Nm_max,Nr_max) !----Allocate--------
      integer, intent(in) :: Nm_max
      integer, intent(in) :: Nr_max 

      allocate( uphi_temp_FR(Nm_max+1,Nr_max), ur_temp_FR(Nm_max+1,Nr_max) )
      allocate( uphi_omg_FR(Nm_max+1,Nr_max), ur_omg_FR(Nm_max+1,Nr_max) )

   end subroutine allocate_steptime_imex

   subroutine deallocate_steptime_imex() !--------------Deallocate-------

      deallocate( uphi_omg_FR, ur_omg_FR  )
      deallocate( uphi_temp_FR, ur_temp_FR )

   end subroutine deallocate_steptime_imex

   subroutine timeloop_imex(Nm_max,Np_max,Nr_max,eta,CFL,n_time_steps,n_checkpoint,n_snapshot,dt, Ra,Pr,mBC, &
                       & l_restart,n_restart,n_restart_point,n_snapshot_point,n_KE,n_KEspec, &
                       & time_scheme_imp,time_scheme_exp,tag,dt_coef,dt_max)
       
      integer :: n_step  
      character(len=100), intent(in) :: tag
      integer, intent(in) :: Nm_max
      integer, intent(in) :: Np_max   
      integer, intent(in) :: Nr_max 
      logical, intent(in) :: l_restart
      integer, intent(in) :: n_restart
      integer, intent(in) :: n_restart_point, n_snapshot_point
      integer, intent(in) :: n_KE, n_KEspec
      real(kind=dp), intent(in) :: eta, CFL
      real(kind=dp), intent(in) :: dt_coef, dt_max
      character(len=100), intent(in) :: time_scheme_imp
      character(len=100), intent(in) :: time_scheme_exp
      integer, intent(in) :: n_time_steps
      integer, intent(in) :: n_checkpoint
      integer, intent(in) :: n_snapshot
      real(kind=dp), intent(in) :: dt
      real(kind=dp), intent(in) :: Ra
      real(kind=dp), intent(in) :: Pr
      character(len=100), intent(in) :: mBC
      integer :: count_snap
      integer :: count_chkpnt
      

      count_snap = n_snapshot_point
      count_chkpnt = n_restart_point
      
      !--------- If restarting -----------------------------------------------
      if (l_restart) then
         if (time_scheme_imp=='BDF1') then
            ! n-1 step
            call Nr_maxLOOP(Nm_max,Np_max,Nr_max,tFR1,omgFR1,upFR1,urFR1,uphi_temp_FR,ur_temp_FR,uphi_omg_FR, &
                               & ur_omg_FR)
            call rhs_temp_from_restart(Nm_max,Nr_max,uphi_temp_FR,ur_temp_FR,tFR1,time_scheme_imp,n_order_tscheme_imp, &
                                       n_order_tscheme_exp)   
            call rhs_vort_from_restart(Nm_max,Nr_max,Ra,Pr,tFR1,uphi_omg_FR,ur_omg_FR,omgFR1,time_scheme_imp, &
                                       & n_order_tscheme_imp, n_order_tscheme_exp)
            call rhs_uphibar_from_restart(Nm_max,Nr_max,upFR1,urFR1,omgFR1,time_scheme_imp,n_order_tscheme_imp, &
                                       n_order_tscheme_exp) 
           
            ! n step
            call Nr_maxLOOP(Nm_max,Np_max,Nr_max,tFR2,omgFR2,upFR2,urFR2,uphi_temp_FR,ur_temp_FR,uphi_omg_FR, &
                               & ur_omg_FR)
            call rhs_temp_from_restart(Nm_max,Nr_max,uphi_temp_FR,ur_temp_FR,tFR2,time_scheme_imp,n_order_tscheme_imp-1, &
                                       n_order_tscheme_exp-1)   
            call rhs_vort_from_restart(Nm_max,Nr_max,Ra,Pr,tFR2,uphi_omg_FR,ur_omg_FR,omgFR2,time_scheme_imp, &
                                       & n_order_tscheme_imp-1, n_order_tscheme_exp-1)
            call rhs_uphibar_from_restart(Nm_max,Nr_max,upFR2,urFR2,omgFR2,time_scheme_imp,n_order_tscheme_imp-1, &
                                       n_order_tscheme_exp-1) 
            !tmp_rhs_buo_term(1,:,:)=rhs_buo_term(n_order_tscheme_imp,:,:)

         elseif (time_scheme_imp=='BDF1') then

            ! n-2 step
            call Nr_maxLOOP(Nm_max,Np_max,Nr_max,tFR1,omgFR1,upFR1,urFR1,uphi_temp_FR,ur_temp_FR,uphi_omg_FR, &
                               & ur_omg_FR)
            call rhs_temp_from_restart(Nm_max,Nr_max,uphi_temp_FR,ur_temp_FR,tFR1,time_scheme_imp,n_order_tscheme_imp, &
                                       n_order_tscheme_exp)   
            call rhs_vort_from_restart(Nm_max,Nr_max,Ra,Pr,tFR1,uphi_omg_FR,ur_omg_FR,omgFR1,time_scheme_imp, &
                                       & n_order_tscheme_imp, n_order_tscheme_exp)
            call rhs_uphibar_from_restart(Nm_max,Nr_max,upFR1,urFR1,omgFR1,time_scheme_imp,n_order_tscheme_imp, &
                                       n_order_tscheme_exp) 

            ! n-1 step
            call Nr_maxLOOP(Nm_max,Np_max,Nr_max,tFR2,omgFR2,upFR2,urFR2,uphi_temp_FR,ur_temp_FR,uphi_omg_FR, &
                               & ur_omg_FR)
            call rhs_temp_from_restart(Nm_max,Nr_max,uphi_temp_FR,ur_temp_FR,tFR2,time_scheme_imp,n_order_tscheme_imp-1, &
                                       n_order_tscheme_exp-1)   
            call rhs_vort_from_restart(Nm_max,Nr_max,Ra,Pr,tFR2,uphi_omg_FR,ur_omg_FR,omgFR2,time_scheme_imp, &
                                       & n_order_tscheme_imp-1, n_order_tscheme_exp-1)
            call rhs_uphibar_from_restart(Nm_max,Nr_max,upFR2,urFR2,omgFR2,time_scheme_imp,n_order_tscheme_imp-1, &
                                       n_order_tscheme_exp-1) 
            tmp_rhs_buo_term(2,:,:)=rhs_buo_term(n_order_tscheme_imp-1,:,:)
           
            ! n step
            call Nr_maxLOOP(Nm_max,Np_max,Nr_max,tFR3,omgFR3,upFR3,urFR3,uphi_temp_FR,ur_temp_FR,uphi_omg_FR, &
                               & ur_omg_FR)
            call rhs_temp_from_restart(Nm_max,Nr_max,uphi_temp_FR,ur_temp_FR,tFR3,time_scheme_imp,n_order_tscheme_imp-2, &
                                       n_order_tscheme_exp-2)   
            call rhs_vort_from_restart(Nm_max,Nr_max,Ra,Pr,tFR3,uphi_omg_FR,ur_omg_FR,omgFR3,time_scheme_imp, &
                                       & n_order_tscheme_imp-2, n_order_tscheme_exp-2)
            call rhs_uphibar_from_restart(Nm_max,Nr_max,upFR3,urFR3,omgFR3,time_scheme_imp,n_order_tscheme_imp-2, &
                                       n_order_tscheme_exp-2) 
            tmp_rhs_buo_term(1,:,:)=rhs_buo_term(n_order_tscheme_imp-2,:,:)


         else 

            call Nr_maxLOOP(Nm_max,Np_max,Nr_max,tFR,omgFR,upFR,urFR,uphi_temp_FR,ur_temp_FR,uphi_omg_FR, &
                               & ur_omg_FR)
         end if
      end if 
      !----------------------------------------------------------------------- 
!------------------------- Time loop begins ----------------------------------------------------------   
      do n_step=1+n_restart,n_time_steps+n_restart 
          
         if (n_step-n_restart==1) then
            dt_new=dt
         end if
            dt_new=dt
            dt_array(:)=dt
         if (n_step-n_restart>1) then
                   
            if (n_step-n_restart==2) then
                call cpu_time(startsteptime)
            end if
            call cpu_time(startNr_maxloop)
            !-------------------- Call Nr_max loop ----------------------------------------------------
               call Nr_maxLOOP(Nm_max,Np_max,Nr_max,tFR,omgFR,upFR,urFR,uphi_temp_FR,ur_temp_FR,uphi_omg_FR, &
                               & ur_omg_FR)
            !------------------------------------------------------------------------------------------
            upFR_prev = upFR
            urFR_prev = urFR

            call cpu_time(finishNr_maxloop) 
            timeNr_maxloop = timeNr_maxloop + finishNr_maxloop-startNr_maxloop

         end if 

         !------------ Compute New dt by enforcing CFL ----------------- 
        ! call compute_new_dt(n_step,n_restart,l_restart,CFL,dt_new,dt_coef,dt_max) 
         !--------------------------------------------------------------

         if ( (l_restart .and. n_step>1+n_restart) .or. (.not. l_restart) ) then
            tot_time=tot_time+dt_new 
         end if

         if (n_step-n_restart>1) then
             tmp_rhs_imp_temp=rhs_imp_temp
             tmp_rhs_exp_temp=rhs_exp_temp
             tmp_rhs_imp_vort=rhs_imp_vort
             tmp_rhs_exp_vort=rhs_exp_vort
             tmp_rhs_buo_term=rhs_buo_term
             tmp_rhs_imp_uphi_bar=rhs_imp_uphi_bar
             tmp_rhs_exp_uphi_bar=rhs_exp_uphi_bar
         end if

         call cpu_time(startNm_maxloop) 
         !-------------------- Call Nm_max loop ------------------------------------------------------------ 
         call Nm_maxLOOP(Nm_max,Nr_max,Ra,Pr,mBC,uphi_temp_FR,ur_temp_FR,uphi_omg_FR,ur_omg_FR,&
                     & n_step,n_restart,time_scheme_imp,time_scheme_exp)
         !--------------------------------------------------------------------------------------------------
         call cpu_time(finishNm_maxloop)
         timeNm_maxloop = timeNm_maxloop + finishNm_maxloop-startNm_maxloop

            
         
         !--------------------- Store Snapshots -----------------------------
         if (mod(n_step,n_snapshot)==0) then
            count_snap = count_snap + 1 
            call store_snapshot(Nm_max,Nr_max,Ra,Pr,eta,tot_time,dt_new,count_snap,tFR,omgFR,urFR,upFR)

         end if
         !------------------------------------------------------------------- 
         
         !-------------------- Store checkpoints ----------------------------
         if (mod(n_step,n_checkpoint)==0) then
            print *, tot_time, "checkpoint save"
            count_chkpnt = count_chkpnt + 1
            call store_checkpoint(Nm_max,Nr_max,count_chkpnt,dt_new,tot_time,tFR,omgFR,urFR,upFR, &
                                  & n_order_tscheme_imp,n_order_tscheme_exp, rhs_imp_temp, &
                                  & rhs_exp_temp,rhs_imp_vort,rhs_exp_vort,rhs_imp_uphi_bar,rhs_exp_uphi_bar, &
                                  & dt_array,n_order_tscheme_max)
 
         end if
         !------------------------------------------------------------------- 

         !-------------------- Store Kinetic Energy (KE) --------------------
         if (mod(n_step,n_KE)==0) then
            call init_output(tag)  ! Open output files 
            !call writeKE_spectral(Nm_max,Nr_max,TFC,tot_time,eta,n_step,omgFR,Ra,Pr,dt_new,tFR)
            call writeKE_spectral(Nm_max,Nr_max,Np_max,TFC,tot_time,eta,n_step,omgFR,Ra,Pr,dt_new,tFR)
            !call writeKE_physical(Np_max,Nr_max,Nm_max,tot_time,Ra,Pr) ! Uncomment for KE calc in physical space 
            call final_output() ! Close output files
         end if
         !-------------------- Store Kinetic Energy (KE) spectra ------------
         if (mod(n_step,n_KEspec)==0) then
            call init_output(tag)  ! Open output files 
            call calculate_spectra(Nm_max,Nr_max,urFR,upFR,n_step)
            call final_output() ! Close output files
         end if
         !------------------------------------------------------------------- 
      
      end do 
!------------------------- Time loop ends ----------------------------------------------------------   
        
   end subroutine timeloop_imex

   subroutine compute_new_dt(n_step,n_restart,l_restart,CFL,dt_new,dt_coef,dt_max)
   
      integer, intent(in) :: n_step, n_restart
      logical, intent(in) :: l_restart
      real(kind=dp), intent(in) :: CFL
      real(kind=dp), intent(inout) :: dt_new
      real(kind=dp), intent(in) :: dt_coef, dt_max
      integer :: i_order
      real(kind=dp) :: dt_cal, dt2, dt_n

      !--------------------- Enforce CFL constraint here ----------------------------
      dt_cal = min(minval(dtval_r),minval(dtval_p))
      dt_n = CFL*dt_cal
      dt2 = min(0.5_dp*(1.0_dp/dt_coef+1.0_dp)*dt_n,dt_max)
      !------------------------------------------------------------------------------
      !------------------- Optimize dt here -----------------------------------------
      if (dt_array(2) > dt_max .and. n_step>1+n_restart) then
         dt_new = dt_max
      elseif (dt_array(2) >= dt_cal .and. n_step>1+n_restart) then
         dt_new = dt2
      elseif (dt_coef*dt_array(2) < dt_cal .and. dt_array(2) < dt_max .and. n_step>1+n_restart) then
         if (dt_array(1)<dt2) then
            dt_new = dt2
         else
            dt_new = dt_array(1)
         end if
      end if
      !------------------------------------------------------------------------------

      !--------------------- Construct dt array here---- ----------------------------
      if ( (.not. l_restart) .and. (n_step-n_restart==1) ) then
            dt_array(:) = dt_new
      else if ( (l_restart) .and. (n_step-n_restart==1) ) then
            dt_new = dt_array(1)  
      end if

      if (n_step-n_restart > 1 ) then
         if (n_order_tscheme_max==1) then
            dt_array(1) = dt_new
         else if (n_order_tscheme_max>1) then
            do i_order=n_order_tscheme_max,2,-1
               dt_array(i_order) = dt_array(i_order-1)
            end do
            dt_array(1) = dt_new
         end if
      end if
      !------------------------------------------------------------------------------
   end subroutine compute_new_dt

end module steptime
