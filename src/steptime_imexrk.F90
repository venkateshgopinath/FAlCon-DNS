module steptime_imexrk

   use double
#if FFTW
   use chebyshev, only: init_chebfftwplan, destroy_chebfftwplan, chebtransform, chebinvtran, D, D2
#elif GG
   use chebyshev, only: chebtransform, chebinvtran, D, D2
#endif
   use init, only: TFC, upFR, urFR, tFR, omgFR, startNm_maxloop, finishNm_maxloop, & 
                   & startNr_maxloop, finishNr_maxloop, timeNm_maxloop, timeNr_maxloop, startsteptime, &
                   dt_old, dt_new, tot_time, psii, ur, up, omg, tt, omgFR_check, finishsteptime, looptime 
   use nonlin, only: Nr_maxloop, dtval_r, dtval_p, dtval_rkr, dtval_rkp
   use output, only: init_output, calculate_spectra, final_output, writeke_spectral, &
                     & store_checkpoint_imexrk, store_snapshot_imexrk
   use fourierloop_imexrk, only: Get_stage_var, RHS_construct_stage, Assembly_stage, Assembly_stage_SA
   use timeschemes, only: rhs_update_wts_imp, wt_lhs_tscheme_imp, wt_rhs_tscheme_imp, n_order_tscheme_imp, &
                          & n_order_tscheme_exp, n_order_tscheme_max, dt_array, rhs_imp_temp, rhs_exp_temp, &
                          & rhs_imp_vort, rhs_exp_vort, rhs_imp_uphi_bar, rhs_exp_uphi_bar, &
                          & n_order_tscheme_exp, butcher_aA, butcher_bA, rhs_update_wts_exp, &
                          & wt_rhs_tscheme_exp, butcher_aD, ars_eqn_check_A, ars_eqn_check_D 
   use mat_assembly, only: mat_build_rk, mat_build_uphibar

   implicit none

   private

   complex(kind=dp), allocatable, public :: uphi_temp_FR(:,:), ur_temp_FR(:,:)
   complex(kind=dp), allocatable, public :: uphi_omg_FR(:,:), ur_omg_FR(:,:)
   
   public :: allocate_steptime_imexrk, deallocate_steptime_imexrk, timeloop_imexrk, solver_log2

contains

   subroutine allocate_steptime_imexrk(Nm_max,Nr_max) !----Allocate--------
      integer, intent(in) :: Nm_max
      integer, intent(in) :: Nr_max 

      allocate( uphi_temp_FR(Nm_max+1,Nr_max), ur_temp_FR(Nm_max+1,Nr_max) )
      allocate( uphi_omg_FR(Nm_max+1,Nr_max), ur_omg_FR(Nm_max+1,Nr_max) )
      
   end subroutine allocate_steptime_imexrk

   subroutine deallocate_steptime_imexrk() !--------------Deallocate-------
     
      deallocate( uphi_omg_FR, ur_omg_FR  )
      deallocate( uphi_temp_FR, ur_temp_FR )

   end subroutine deallocate_steptime_imexrk

   subroutine timeloop_imexrk(Nm_max,Np_max,Nr_max,eta,CFL,n_time_steps,n_checkpoint,n_snapshot,dt,Ra,Pr, &
                       & l_restart,l_optimizedt,n_restart,n_restart_point,n_snapshot_point,n_KE,n_KEspec, &
                       & time_scheme_type,time_scheme_imp,time_scheme_exp,tag,dt_coef,dt_max,mBC, &
                       & lm,buo_tscheme,totaltime,l_vartimestep)
       
      integer :: n_step
      integer, intent(in) :: lm 
      integer, intent(in) :: Nm_max
      integer, intent(in) :: Np_max   
      integer, intent(in) :: Nr_max 
      logical, intent(in) :: l_restart, l_optimizedt
      integer, intent(in) :: n_restart
      integer, intent(in) :: n_restart_point, n_snapshot_point
      integer, intent(in) :: n_KE, n_KEspec
      real(kind=dp), intent(in) :: eta, CFL
      real(kind=dp), intent(in) :: dt_coef, dt_max
      character(len=100), intent(in) :: time_scheme_type
      character(len=100), intent(in) :: time_scheme_imp
      character(len=100), intent(in) :: time_scheme_exp
      character(len=100), intent(in) :: mBC  
      character(len=100), intent(in) :: buo_tscheme  
      integer, intent(in) :: n_time_steps
      integer, intent(in) :: n_checkpoint
      integer, intent(in) :: n_snapshot
      real(kind=dp), intent(in) :: dt
      real(kind=dp), intent(in) :: Ra
      real(kind=dp), intent(in) :: Pr
      real(kind=dp), intent(in) :: totaltime
      logical, intent(in) :: l_vartimestep
      integer :: count_snap
      integer :: count_chkpnt
      integer :: rk_stage, Nm
      character(len=100), intent(in) :: tag
      real(kind=dp) :: C

      C=((2.0_dp/real(Nr_max-1,kind=dp))**(0.5_dp))

      count_snap = n_snapshot_point
      count_chkpnt = n_restart_point

      call rhs_update_wts_imp(time_scheme_imp,wt_lhs_tscheme_imp,wt_rhs_tscheme_imp,n_order_tscheme_imp)
      call rhs_update_wts_exp(time_scheme_imp,time_scheme_exp,wt_rhs_tscheme_exp,n_order_tscheme_exp)

      do Nm=0,Nm_max
            call mat_build_rk(Nr_max,Nm,mBC) ! Build the operator matrix solving for psi and factorize them 
      end do
      finishsteptime=0.0_dp
      startsteptime=0.0_dp 
      looptime=0.0_dp
!------------------------- Time loop begins ----------------------------------------------------------   
      do n_step=1+n_restart,n_time_steps+n_restart 
          
         if (n_step-n_restart==1) then
            dt_old=dt
            dt_new=dt
            dt_array(:)=dt 
         end if

         if (l_vartimestep .eqv. .TRUE.) then ! If variable timestep using CFL is on
            if (n_step-n_restart==1 .and. l_restart) then
               !-------------------- Call Nr_max loop ---------------------------------------------------------
               call Nr_maxLOOP(Nm_max,Np_max,Nr_max,tFR,omgFR,upFR,urFR,uphi_temp_FR,ur_temp_FR,uphi_omg_FR, &
                                     & ur_omg_FR)
               !-----------------------------------------------------------------------------------------------
               !------------ Compute New dt by enforcing CFL ----------------- 
               call compute_new_dt(n_step,n_restart,l_restart,l_optimizedt,CFL,dt_new,dt_coef,dt_max,Pr) 
               !--------------------------------------------------------------
            end if
            if (n_step-n_restart>1) then
               if ( ( .not. ars_eqn_check_A ) .or. ( .not. ars_eqn_check_D ) ) then
                  !-------------------- Call Nr_max loop ---------------------------------------------------------
                  call Nr_maxLOOP(Nm_max,Np_max,Nr_max,tFR,omgFR,upFR,urFR,uphi_temp_FR,ur_temp_FR,uphi_omg_FR, &
                                        & ur_omg_FR)
                  !-----------------------------------------------------------------------------------------------
                  !------------ Compute New dt by enforcing CFL ----------------- 
                  call compute_new_dt(n_step,n_restart,l_restart,l_optimizedt,CFL,dt_new,dt_coef,dt_max,Pr) 
                  !--------------------------------------------------------------
               else
                  !------------ Compute New dt by enforcing CFL ----------------- 
                  call compute_new_dt(n_step,n_restart,l_restart,l_optimizedt,CFL,dt_new,dt_coef,dt_max,Pr) 
                  !--------------------------------------------------------------
               end if
            end if  
         else ! If constant timestep is on 
            dt_old=dt
            dt_new=dt
            dt_array(:)=dt 
         end if

         if (n_step-n_restart==2) then
            call cpu_time(startsteptime)
         end if

         do rk_stage=1,n_order_tscheme_exp ! LOOP for RK stages (rk_stage) 
               if (rk_stage>1) then
                  
                  !------ Solve for and update stage variables temp1, omg1, psi1 at each stage---- 
                  call Get_stage_var(Nm_max,Nr_max,rk_stage,tFR,omgFR,psii, &
                                     & upFR,urFR,n_step,n_restart,Ra,Pr,mBC,buo_tscheme)
                  !-------------------------------------------------------------------------------

               end if

               if ( (n_step-n_restart>1 .or. rk_stage > 1 .or. l_restart )  .and. & 
                  &   (sum(butcher_aA(:,rk_stage)) + butcher_bA(rk_stage) /=0) ) then 
                      ! Call Nrloop only if atleast one element in the column of Butcher's table for advection is non-zero)  
                  !-------------------- Call Nr_max loop ---------------------------------------------------------
                  call Nr_maxLOOP(Nm_max,Np_max,Nr_max,tFR,omgFR,upFR,urFR,uphi_temp_FR,ur_temp_FR,uphi_omg_FR, &
                                  & ur_omg_FR)
                  !-----------------------------------------------------------------------------------------------
               end if

               !-------------------- Construct RHS for all stages --------------------------------------------- 
               call RHS_construct_stage(Nm_max,Nr_max,Ra,Pr,uphi_temp_FR,ur_temp_FR,uphi_omg_FR,ur_omg_FR,& 
                                & time_scheme_type,rk_stage,tFR,omgFR,urFR,upFR)
               !-------------------------------------------------------------------------------------------

         end do ! End LOOP for RK stages 

         ! ASSEMBLY STAGE ------------------------------------------------------
         if ( ( .not. ars_eqn_check_A ) .or. ( .not. ars_eqn_check_D ) ) then ! If not satisfying equation (2.3) in ARS_97 paper
         !   !--- Assembly for non-stiffly accurate schemes ---------------------
            call Assembly_stage(Nm_max,Nr_max,dt_new,tFR,omgFR,upFR,urFR,lm,mBC)
            !call Assembly_stage(Nm_max,Nr_max,dt_new,tFR,omgFR,upFR,urFR,lm,mBC,rmin,rmax)
         else
            !--- Assembly for stiffly accurate schemes -------------------------
            call Assembly_stage_SA(Nm_max,Nr_max,tFR,omgFR,upFR)
         end if
         !----------------------------------------------------------------------     

         tot_time=tot_time+dt_new 
 
         !looptime=0.0_dp
         !call cpu_time(looptime)
         !finishsteptime=finishsteptime+looptime
         call cpu_time(finishsteptime)

         if (tot_time<totaltime) then

            ! Calculate temperature in FC space at 0-mode --       
            call chebtransform(Nr_max,tFR(1,:),TFC)
            !-----------------------------------------------

            !--------------------- Store Snapshots -----------------------------
            if (mod(n_step,n_snapshot)==0) then
               count_snap = count_snap + 1 
               call store_snapshot_imexrk(Nm_max,Nr_max,Ra,Pr,eta,tot_time,dt_new,count_snap,tFR,omgFR,urFR,upFR,omgFR_check)

            end if
            !------------------------------------------------------------------- 
            
            !-------------------- Store checkpoints ----------------------------
            if (mod(n_step,n_checkpoint)==0) then
               print *, tot_time, "checkpoint save"
               count_chkpnt = count_chkpnt + 1
               call store_checkpoint_imexrk(Nm_max,Nr_max,count_chkpnt,dt_new,tot_time,tFR,omgFR,urFR,upFR, &
                                     & n_order_tscheme_imp,n_order_tscheme_exp, rhs_imp_temp, &
                                     & rhs_exp_temp,rhs_imp_vort,rhs_exp_vort,rhs_imp_uphi_bar, &
                                     & rhs_exp_uphi_bar,dt_array,n_order_tscheme_max)
    
            end if
            !------------------------------------------------------------------- 
            !-------------------- Store Kinetic Energy (KE) --------------------
            if (mod(n_step,n_KE)==0) then
               call init_output(tag)  ! Open output files 
               call writeKE_spectral(Nm_max,Nr_max,Np_max, TFC,tot_time,eta,omgFR,Ra,Pr,dt_new,tFR,mBC)
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
         else 
            call solver_log2(n_step)
            call exit()
         end if
      
      end do 
!------------------------- Time loop ends ----------------------------------------------------------   
   end subroutine timeloop_imexrk

   subroutine solver_log2(n_step)
    integer, intent(in) :: n_step
    integer :: logunit
      open(newunit=logunit,file="solver_log2.txt",status="unknown",form="formatted", action="write")

      call cpu_time(finishsteptime)
      write (logunit,*) "Total time taken =",finishsteptime-startsteptime,"seconds."
      write (logunit,*) "Average time taken for time loop =",(finishsteptime-startsteptime)/real(n_step-1,kind=dp),"seconds."
     
      close(logunit)

   end subroutine solver_log2

   subroutine compute_new_dt(n_step,n_restart,l_restart,l_optimizedt,CFL,dt_new,dt_coef,dt_max,Pr)
   
      integer, intent(in) :: n_step, n_restart
      logical, intent(in) :: l_restart, l_optimizedt
      real(kind=dp), intent(in) :: CFL
      real(kind=dp), intent(in) :: Pr
      real(kind=dp), intent(in) :: dt_coef, dt_max
      real(kind=dp), intent(out) :: dt_new
      integer :: i_order
      real(kind=dp) :: dt_cal, dt2, dt_n

      if (n_step>1+n_restart .or. l_restart) then 
         !--------------------- Enforce CFL constraint here ----------------------------------------
         dt_cal = min(minval(dtval_r),minval(dtval_p))
         !print *, minval(dtval_r), (2.0_dp/3.0_dp)*minval(dtval_p), minval(dtval_p),"here" 
         if (Pr==1 .or. Pr>1) then 
            dt_n = max(1.0_dp/Pr,1.0_dp)*CFL*dt_cal
         else
            dt_n = min(1.0_dp/Pr,1.0_dp)*CFL*dt_cal
         end if
         dt2 = min(0.5_dp*(1.0_dp/dt_coef+1.0_dp)*dt_n,dt_max)
         !------------------------------------------------------------------------------------------   
      end if
 
      if (l_optimizedt) then
         !------------------- Optimize dt here -------------------------------------------------------
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
         !------------------------------------------------------------------------------------------   
      else
         dt_new=dt2
      end if
      !--------------------- Construct dt array here---- ----------------------------
      
      if ( (.not. l_restart) .and. (n_step-n_restart==1) ) then
         dt_array(:) = dt_new
      else if ( (l_restart) .and. (n_step-n_restart==1) ) then
         dt_array(1) = dt_new
      end if 

      if (n_step-n_restart > 1 ) then
         if (n_order_tscheme_max==1) then
            dt_array(1) = dt_new
         else if (n_order_tscheme_max>1) then
            do i_order=2,2,-1
               dt_array(i_order) = dt_array(i_order-1)
            end do
            dt_array(1) = dt_new
         end if
      end if
      !------------------------------------------------------------------------------
   end subroutine compute_new_dt

end module steptime_imexrk
