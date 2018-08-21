module fourierloop_imexrk
   !$ use OMP_LIB

   use double
   use constants, only: ii
   use chebyshev, only: chebtransform, chebinvtran, chebinvtranD1, chebinvtranD2, chebinvtranD1D2, t, D, D2
   use init, only: TFC, temp_spec,omg_spec,startmatbuild, &
                   & finishmatbuild, finishmatbuild, time_matbuild, radius, r_radius, r_radius2, &
                   & dw_rmin, d2w_rmin, dw_rmax, d2w_rmax, uphi_bar_spec, upFC, ur, omg, omgFR_check
   use mat_assembly, only: LAPpsi_all, AT_all_IRK, AF_all_IRK, A_uphi_all_IRK, IPIV_uphi_IRK, IPIV1_IRK, IPIV2_IRK, &
                           & mat_build_IRK, mat_build_uphibar_IRK, IPIV1_lap
   use algebra, only: matsolve
   use timeschemes, only: wt_lhs_tscheme_imp, wt_rhs_tscheme_imp, n_order_tscheme_imp, wt_rhs_tscheme_exp, & 
                          & n_order_tscheme_exp, rhs_update_wts_imp, rhs_update_wts_exp, &
                          & irk_max, dt_array, butcher_aA, butcher_bA, butcher_aD, butcher_bD, &
                          & ars_eqn_check_A, ars_eqn_check_D, rhs_imp_temp, rhs_exp_temp, rhs_imp_vort, &
                          & rhs_exp_vort, rhs_imp_uphi_bar, rhs_exp_uphi_bar, diag_diff, diag_index
   use rhs_create_imexrk, only: rhs_construct_temp_a, rhs_construct_vort_a, rhs_construct_temp_d, &
                                & rhs_construct_uphibar_a, rhs_construct_uphibar_d, &
                                & rhs_construct_vort_d, rhs_construct_buo   

   implicit none

   private 

   complex(kind=dp), allocatable :: real_rhs_temp(:)
   complex(kind=dp), allocatable :: rhs_omg(:), real_rhs_omg(:) 
   complex(kind=dp), allocatable :: real_rhs_psi(:), real_d_rhs_psi(:), real_d2_rhs_psi(:) 
   complex(kind=dp), allocatable :: rhs1(:), rhs2(:), rhs(:), rhsf(:), rhs_vort(:), rhs_b(:)
   real(kind=dp), allocatable :: rhs_r(:), rhs_i(:), rhsf_r(:), rhsf_i(:)
   real(kind=dp), allocatable :: rhs_uphi(:)

   complex(kind=dp), allocatable, public :: dtempdt1_a(:,:,:), dtempdt1_d(:,:,:)
   complex(kind=dp), allocatable, public :: duphibar_dt1_a(:,:), duphibar_dt1_d(:,:)
   complex(kind=dp), allocatable, public :: domgdt1_a(:,:,:), domgdt1_d(:,:,:)
   complex(kind=dp), allocatable, public :: rhs_buo(:,:,:)

   character :: TRANS='N'

   public :: allocate_fourierloop_imexrk, deallocate_fourierloop_imexrk, Get_stage_var, &
             & RHS_construct_stage, Assembly_stage, Assembly_stage_SA

contains

   subroutine allocate_fourierloop_imexrk(Nm_max,Nr_max)

      integer, intent(in) :: Nr_max 
      integer, intent(in) :: Nm_max 

      allocate( real_rhs_temp(Nr_max) )
      allocate( rhs_omg(Nr_max),real_rhs_omg(Nr_max),real_rhs_psi(Nr_max) )
      allocate( real_d_rhs_psi(Nr_max), real_d2_rhs_psi(Nr_max) )
      allocate( rhs1(Nr_max), rhs2(Nr_max), rhs(Nr_max), rhsf(2*Nr_max), rhs_vort(Nr_max), rhs_b(Nr_max) )  
      allocate( rhs_r(Nr_max), rhs_i(Nr_max), rhsf_r(2*Nr_max), rhsf_i(2*Nr_max) )
      allocate( rhs_uphi(Nr_max) )
      allocate( dtempdt1_a(n_order_tscheme_exp,Nm_max+1,Nr_max), &
                & dtempdt1_d(n_order_tscheme_exp,Nm_max+1,Nr_max) )
      allocate( duphibar_dt1_a(n_order_tscheme_exp,Nr_max), &
                & duphibar_dt1_d(n_order_tscheme_exp,Nr_max) )
      allocate( domgdt1_a(n_order_tscheme_exp,Nm_max+1,Nr_max), &
                & domgdt1_d(n_order_tscheme_exp,Nm_max+1,Nr_max) )
      allocate( rhs_buo(n_order_tscheme_exp,Nm_max+1,Nr_max) )
   
   end subroutine allocate_fourierloop_imexrk

   subroutine deallocate_fourierloop_imexrk
   
      deallocate( dtempdt1_a, dtempdt1_d, duphibar_dt1_a, duphibar_dt1_d, domgdt1_a, domgdt1_d )
      deallocate( rhs_r, rhs_i, rhsf_r, rhsf_i )  
      deallocate( rhs1, rhs2, rhs,rhsf, rhs_uphi, rhs_vort, rhs_buo )  
      deallocate( rhs_omg,real_rhs_omg,real_rhs_psi,real_d_rhs_psi, real_d2_rhs_psi )
      deallocate( real_rhs_temp )
    
   end subroutine deallocate_fourierloop_imexrk

   subroutine Get_stage_var(Nm_max,Nr_max,rk_stage,tFR,omgFR,psii,upFR,urFR,n_step,n_restart,Ra,Pr,mBC,lm,buo_tscheme)


      integer, intent(in) :: Nm_max
      integer, intent(in) :: Nr_max 
      real(kind=dp), intent(in) :: Ra
      real(kind=dp), intent(in) :: Pr
      character(len=100), intent(in) :: mBC  
      character(len=100), intent(in) :: buo_tscheme  
      integer, intent(in) :: rk_stage, n_step, n_restart
      complex(kind=dp), intent(out) :: tFR(Nm_max+1,Nr_max),omgFR(Nm_max+1,Nr_max) 
      complex(kind=dp), intent(out) :: psii(Nm_max+1,Nr_max)
      complex(kind=dp), intent(out) :: upFR(Nm_max+1,Nr_max), urFR(Nm_max+1,Nr_max)
      integer, intent(in) :: lm 
      !Local variables ---------------------------------------------------------
      complex(kind=dp) :: rhs_psi(Nr_max), F_dtemp_d(Nr_max), F_domg_d(Nr_max)
      complex(kind=dp) :: F_dtemp_a(Nr_max), F_domg_a(Nr_max), F_buo(Nr_max)
      complex(kind=dp) :: F_duphibar_a(Nr_max), F_duphibar_d(Nr_max)
      complex(kind=dp) :: rhs_stage_temp(Nr_max), rhs_stage_uphibar(Nr_max), rhs_stage_omg(Nr_max)  
      complex(kind=dp) :: rhs_temp_FC(Nr_max), rhs_omg_FC(2*Nr_max)     
      complex(kind=dp) :: rhs_stage_omg_full(2*Nr_max)
      integer :: i,j,Nm,INFO1,INFO2,Nr_max2 ! Nm -> azimuthal 'n' loop over Fourier modes
      complex(kind=dp) :: D1upFR(Nr_max)
      real(kind=dp) :: t_ref, t_final
      complex(kind=dp) :: dpsi_rmin
      complex(kind=dp) :: d2psi_rmin
      complex(kind=dp) :: dpsi_rmax
      complex(kind=dp) :: d2psi_rmax

      rhs1(:)=0.0_dp 
      rhsf(:)=0.0_dp 
      
      Nr_max2=2*Nr_max
      
      if ( diag_diff .eqv. .FALSE. ) then ! If implicit Butcher's table has all diagonal elements equal to each other
         do Nm=0,Nm_max  
            if ( (n_step-n_restart==1) .or. (dt_array(1)/=dt_array(2)) ) then 
               call cpu_time(startmatbuild) 
               call mat_build_IRK(Nr_max,dt_array(1),Nm,mBC,butcher_aD(diag_index,diag_index),diag_index,Pr) ! Build the operator matrices for (T, omega) and factorize them 
               call mat_build_uphibar_IRK(Nr_max,dt_array(1),mBC,butcher_aD(diag_index,diag_index),diag_index,Pr) ! Build the operatir matrix for uphibar and factorize them
               call cpu_time(finishmatbuild) 
               time_matbuild=time_matbuild + finishmatbuild - startmatbuild
            end if
         end do
         AT_all_IRK(:,:,:,rk_stage)=AT_all_IRK(:,:,:,diag_index)
         AF_all_IRK(:,:,:,rk_stage)=AF_all_IRK(:,:,:,diag_index)
         A_uphi_all_IRK(:,:,:,rk_stage)=A_uphi_all_IRK(:,:,:,diag_index)
         IPIV1_IRK(:,:,rk_stage)=IPIV1_IRK(:,:,diag_index) 
         IPIV2_IRK(:,:,rk_stage)=IPIV2_IRK(:,:,diag_index) 
         IPIV_uphi_IRK(:,:,rk_stage)=IPIV_uphi_IRK(:,:,diag_index) 
      elseif ( diag_diff .eqv. .TRUE. ) then ! If implicit Butcher's table has all diagonal elements unequal to each other
         do Nm=0,Nm_max  
            if ( (n_step-n_restart==1) .or. (dt_array(1)/=dt_array(2)) ) then 
               call cpu_time(startmatbuild) 
               call mat_build_IRK(Nr_max,dt_array(1),Nm,mBC,butcher_aD(diag_index,diag_index),diag_index,Pr) ! Build the operator matrices for (T, omega) and factorize them 
               call mat_build_uphibar_IRK(Nr_max,dt_array(1),mBC,butcher_aD(diag_index,diag_index),diag_index,Pr) ! Build the operatir matrix for uphibar and factorize them
               call cpu_time(finishmatbuild) 
               time_matbuild=time_matbuild + finishmatbuild - startmatbuild
            end if
         end do
         if ( rk_stage/=diag_index ) then
            if ( (n_step-n_restart==1) .or. (dt_array(1)/=dt_array(2)) ) then 
               if ( butcher_aD(rk_stage,rk_stage)==butcher_aD(diag_index,diag_index) ) then
                  AT_all_IRK(:,:,:,rk_stage)=AT_all_IRK(:,:,:,diag_index)
                  AF_all_IRK(:,:,:,rk_stage)=AF_all_IRK(:,:,:,diag_index)
                  A_uphi_all_IRK(:,:,:,rk_stage)=A_uphi_all_IRK(:,:,:,diag_index)
                  IPIV1_IRK(:,:,rk_stage)=IPIV1_IRK(:,:,diag_index) 
                  IPIV2_IRK(:,:,rk_stage)=IPIV2_IRK(:,:,diag_index) 
                  IPIV_uphi_IRK(:,:,rk_stage)=IPIV_uphi_IRK(:,:,diag_index) 
               else
                  do Nm=0,Nm_max  
                        call mat_build_IRK(Nr_max,dt_array(1),Nm,mBC,butcher_aD(rk_stage,rk_stage),rk_stage,Pr) ! Build the operator matrices for (T, omega) and factorize them 
                        call mat_build_uphibar_IRK(Nr_max,dt_array(1),mBC,butcher_aD(rk_stage,rk_stage),rk_stage,Pr) ! Build the operatir matrix for uphibar and factorize them
                  end do
               end if
            end if
         end if
      end if
      
      !-------------- Loop over the Fourier modes ---------------------------------------------------------------------------
      !$omp parallel & 
      !$omp private(Nm,i,j,F_dtemp_d,F_dtemp_a,F_duphibar_d,F_duphibar_a,F_domg_d,F_domg_a,F_buo, &
      !$omp & rhs_stage_temp,rhs_stage_uphibar,rhs_stage_omg,rhs_stage_omg_full,rhs_r,rhs_i,rhsf_r,rhsf_i, &
      !$omp & D1upFR,rhs_temp_FC,real_rhs_temp,rhs_omg_FC,rhs_omg,real_rhs_omg,rhs_psi, &
      !$omp & real_rhs_psi,real_d_rhs_psi,real_d2_rhs_psi,dpsi_rmin,d2psi_rmin,dpsi_rmax,d2psi_rmax) default(shared)  
         t_ref= OMP_GET_WTIME ()
      !$omp do    
      do Nm=0,Nm_max  

         !if ( (n_step-n_restart==1) .or. (dt_array(1)/=dt_array(2)) ) then !.or. &
         !   !& (butcher_aD(rk_stage,rk_stage)/=butcher_aD(rk_stage-1,rk_stage-1)) ) then ! Have to generalize 19/3/18 
         !   ! Call 'mat_build' only if dt_current and dt_previous are different OR if the butcher table for implicit part has unequal diagonal values
         !   call cpu_time(startmatbuild) 
         !   call mat_build(Nr_max,dt_array(1),Nm,mBC,butcher_aD(rk_stage,rk_stage),Pr) ! Build the operator matrices and factorize them 
         !   call mat_build_uphibar(Nr_max,dt_array(1),mBC,butcher_aD(rk_stage,rk_stage),Pr)
         !   call cpu_time(finishmatbuild) 
         !   time_matbuild=time_matbuild + finishmatbuild - startmatbuild
         !end if

         F_dtemp_d(:)=0.0_dp 
         F_dtemp_a(:)=0.0_dp 
         F_duphibar_d(:)=0.0_dp 
         F_duphibar_a(:)=0.0_dp 
         F_domg_d(:)=0.0_dp  
         F_domg_a(:)=0.0_dp 
         F_buo(:)=0.0_dp 
         rhs_stage_temp(:)=0.0_dp  
         rhs_stage_uphibar(:)=0.0_dp  
         rhs_stage_omg(:)=0.0_dp  
         dpsi_rmin = 0.0_dp
         d2psi_rmin = 0.0_dp
         dpsi_rmax = 0.0_dp
         d2psi_rmax = 0.0_dp

         !-- Apply weights to RHS using implicit Butcher's table: summation of (a_i F_i) 
         do i=1,rk_stage-1 
            F_dtemp_d(:)=F_dtemp_d(:)+butcher_aD(rk_stage,i)*dtempdt1_d(i,Nm+1,:) 
         end do

         !-- Apply weights to RHS using explicit Butcher's table: summation of (\hat{a}_i F_i) 
         do i=1,rk_stage-1 
            F_dtemp_a(:)=F_dtemp_a(:)+butcher_aA(rk_stage,i)*dtempdt1_a(i,Nm+1,:) 
         end do

         !-- RHS temp = Y_0 + dt * summation of (a_i * F_i) + summation of (\hat{a}_i F_i) 
         rhs_stage_temp(:) = temp_spec(Nm+1,:) + dt_array(1)*(F_dtemp_d + F_dtemp_a)

         !---------------------------------------------------------  

         if(Nm==0) then ! Apply BC for temperature
           rhs_stage_temp(1)=1.0_dp
           rhs_stage_temp(Nr_max)=0.0_dp
         else 
           rhs_stage_temp(1)=0.0_dp
           rhs_stage_temp(Nr_max)=0.0_dp
         end if 

         rhs_r=real(rhs_stage_temp)
         rhs_i=aimag(rhs_stage_temp)
         
         !******Solve for stage temperature CALL DGETRS A*vt=rhs ****************************
         call matsolve(TRANS, Nr_max, AT_all_IRK(:,:,Nm+1,rk_stage), IPIV1_IRK(:,Nm+1,rk_stage), rhs_r, rhs_i,INFO1)
         !*********************************************************************************** 
         do i=1,Nr_max
           rhs_temp_FC(i)=cmplx(rhs_r(i),rhs_i(i),kind=dp)
         end do
         
         call chebinvtran(Nr_max,rhs_temp_FC,real_rhs_temp)
         tFR(Nm+1,:)=real_rhs_temp(:) ! Update stage variable temperature

         ! Solve for uphi_bar here -------------------------------------------------------------
         if (Nm==0) then 

            !-- Apply weights to RHS using implicit Butcher's table: summation of (a_i F_i) 
            do i=1,rk_stage-1 
               !F_duphibar_d(:)=F_duphibar_d(:)+butcher_aD(rk_stage,i)*duphibar_dt1_d(i,:) 
               F_duphibar_d(:)=F_duphibar_d(:)+butcher_aD(rk_stage,i)*domgdt1_d(i,1,:) 
            end do
            
            !-- Apply weights to RHS using explicit Butcher's table: summation of (\hat{a}_i F_i) 
            do i=1,rk_stage-1 
               !F_duphibar_a(:)=F_duphibar_a(:)+butcher_aA(rk_stage,i)*duphibar_dt1_a(i,:) 
               F_duphibar_a(:)=F_duphibar_a(:)+butcher_aA(rk_stage,i)*domgdt1_a(i,1,:) 
            end do

            !-- RHS temp = Y_0 + dt * summation of (a_i * F_i) + summation of (\hat{a}_i F_i) 
            rhs_stage_uphibar(:) = uphi_bar_spec(:) + dt_array(1)*(F_duphibar_d + F_duphibar_a)

            !---------------------------------------------------------  

            ! Apply BC for uphi bar
            rhs_stage_uphibar(1)=0.0_dp
            rhs_stage_uphibar(Nr_max)=0.0_dp

            rhs_r=real(rhs_stage_uphibar)
            rhs_i=0.0_dp

            !******Solve for stage temperature CALL DGETRS A*vt=rhs ****************************
            call matsolve(TRANS, Nr_max, A_uphi_all_IRK(:,:,1,rk_stage), IPIV_uphi_IRK(:,1,rk_stage), rhs_r, rhs_i,INFO1)
            !*********************************************************************************** 

            do i=1,Nr_max
              !upFC(1,i)=cmplx(rhs_r(i),rhs_i(i),kind=dp)
              upFC(1,i)=cmplx(rhs_r(i),0.0_dp,kind=dp)
            end do

            call chebinvtran(Nr_max,upFC(1,:),upFR(1,:)) ! Update stage variable uphi_bar
            call chebinvtranD1(Nr_max,upFC(1,:),D1upFR(:))
            urFR(1,:) = 0.0_dp 

            do i=1,Nr_max
               omgFR(1,i)=upFR(1,i)*r_radius(i) + D1upFR(i)
            end do
            
         !--------------------------------------------------------------------------------------------
         else

            ! --- Construct RHS Buoyancy term which goes to vorticity solve ----- 
            do i=1,Nr_max
               rhs_buo(rk_stage,Nm+1,i)= - (Ra/Pr)*(ii*real(Nm,kind=dp)*(r_radius(i))*tFR(Nm+1,i))
                       ! buoyancy term
            end do
            !-------------------------------------------------------------------- 

            !-- Apply weights to RHS using implicit Butcher's table: summation of (a_i F_i) 
            do i=1,rk_stage-1 
               F_domg_d(:)=F_domg_d(:)+butcher_aD(rk_stage,i)*domgdt1_d(i,Nm+1,:) 
            end do
            
            !-- Apply weights to RHS using explicit Butcher's table: summation of (\hat{a}_i F_i) 
            do i=1,rk_stage-1 
               F_domg_a(:)=F_domg_a(:)+butcher_aA(rk_stage,i)*domgdt1_a(i,Nm+1,:) 
            end do
            
            if (buo_tscheme=='IMP') then ! IMPLICIT TREATMENT
               !-- Apply weights to RHS BUO using implicit Butcher's table: summation of (a_i Fbuo_i)  
               do i=1,rk_stage ! Here the loop for full rk_stage as we use the current stage Temperature to construct buoyancy
                 F_buo(:)=F_buo(:)+butcher_aD(rk_stage,i)*rhs_buo(i,Nm+1,:) 
               end do
            elseif (buo_tscheme=='EXP') then ! EXPLICIT TREATMENT
               !-- Apply weights to RHS BUO using implicit Butcher's table: summation of (a_i Fbuo_i) 
               do i=1,rk_stage-1 ! Here the loop for full rk_stage as we use the current stage Temperature to construct buoyancy
                 F_buo(:)=F_buo(:)+butcher_aA(rk_stage,i)*rhs_buo(i,Nm+1,:) 
               end do
            end if

            !-- RHS omg = Y_0 + dt * summation of (a_i * F_i) + summation of (\hat{a}_i F_i) + summation of (a_i * Fbuo_i)
            rhs_stage_omg(:) = omg_spec(Nm+1,:) + dt_array(1)*(F_domg_d + F_domg_a + F_buo) ! bug found here  was put "F_dtemp_d + F_dtemp_a + F_buo" Dec 1 2017

            rhs_stage_omg(1)=0.0_dp ! Apply BC for stream function in the vorticity equation
            rhs_stage_omg(Nr_max)=0.0_dp ! Apply BC for stream function in the vorticity equation       

            do i=1,Nr_max
               rhs_stage_omg_full(i)=0.0_dp
               rhs_stage_omg_full(i+Nr_max)=rhs_stage_omg(i)
            end do
            !------------------------------------------------------- 

            rhsf_r=real(rhs_stage_omg_full)
            rhsf_i=aimag(rhs_stage_omg_full)

            !**************************** CALL DGETRS AF*vf=rhsf **************************
            call matsolve(TRANS, Nr_max2, AF_all_IRK(:,:,Nm+1,rk_stage), IPIV2_IRK(:,Nm+1,rk_stage), rhsf_r, rhsf_i,INFO2)
            !******************************************************************************

            do i=1,2*Nr_max
              rhs_omg_FC(i)=cmplx(rhsf_r(i),rhsf_i(i),kind=dp)
            end do
            
            do i=1,Nr_max

               rhs_omg(i)=rhs_omg_FC(Nr_max+i)
               rhs_psi(i)=rhs_omg_FC(i)

            end do

            call chebinvtran(Nr_max,rhs_omg,real_rhs_omg)
            omgFR(Nm+1,:)=real_rhs_omg(:)
            call chebinvtran(Nr_max,rhs_psi,real_rhs_psi)
            call chebinvtranD1(Nr_max,rhs_psi,real_d_rhs_psi)
            call chebinvtranD2(Nr_max,rhs_psi,real_d2_rhs_psi)
            ! OMEGA FIX and CHECK start  
            !-------------------------- Johnston strategy for BC ----------------------------------------------
            !do i=1,lm ! Evaluate at rmin
            !      d2psi_rmin=d2psi_rmin + d2w_rmin(i)*real_rhs_psi(i+1) ! Summation 2nd derivative
            !      dpsi_rmin=dpsi_rmin + dw_rmin(i)*real_rhs_psi(i+1)    ! Summation 1st derivative
            !end do
               
            !do i=1,lm ! Evaluate at rmax
            !      d2psi_rmax=d2psi_rmax + d2w_rmax(i)*real_rhs_psi(Nr_max-i)  ! Summation 2nd derivative
            !      dpsi_rmax=dpsi_rmax + dw_rmax(i)*real_rhs_psi(Nr_max-i)     ! Summation 1st derivative
            !end do 

            ! ------------ Use Lagrange polynomials for approximation of No-slip boundary conditions ------------------- 
            !omgFR(Nm+1,1)=-1.0_dp*(d2psi_rmin - 2.0_dp*dw_rmin(1)*dpsi_rmin)            ! Apply omega BC ! Johnston strategy
            !omgFR(Nm+1,Nr_max)=-1.0_dp*(d2psi_rmax - 2.0_dp*dw_rmax(1)*dpsi_rmax)       ! Apply omega BC ! Johnston strategy
            ! --------------------------------------------------------------------------------------------------
            !----------------------------------------------------------------------------------------------------    

            !do i=1,Nr_max
            !   omgFR(Nm+1,i) = -1.0_dp*(r_radius(i)*real_d_rhs_psi(i)+real_d2_rhs_psi(i)- &
            !                    & real(Nm,kind=dp)*real(Nm,kind=dp)*r_radius2(i)*real_rhs_psi(i))
            !end do

            !omgFR(Nm+1,1) = -1.0_dp*(real_d2_rhs_psi(1))
            !omgFR(Nm+1,Nr_max) = -1.0_dp*(real_d2_rhs_psi(Nr_max))

            !do i=1,Nr_max
            !   omgFR_check(Nm+1,i) = -1.0_dp*(r_radius(i)*real_d_rhs_psi(i)+real_d2_rhs_psi(i)- &
            !                    & real(Nm,kind=dp)*real(Nm,kind=dp)*r_radius2(i)*real_rhs_psi(i))
            !end do
            !print *, ((abs(omgFR_check(2)))), ((abs(omgFR(Nm+1,2)))), & 
            !         & abs((omgFR_check(2)))-abs((omgFR(Nm+1,2))),   "comparison" , Nm 
            ! OMEGA FIX and CHECK end  

            do i=1,Nr_max
               psii(Nm+1,i)=real_rhs_psi(i)
               upFR(Nm+1,i)=-1.0_dp*real_d_rhs_psi(i)
               urFR(Nm+1,i)=ii*real(Nm,kind=dp)*r_radius(i)*real_rhs_psi(i)
            end do

            call chebtransform(Nr_max,upFR(Nm+1,:),upFC(Nm+1,:))

         end if 
      end do
      !$omp end do
         t_final= OMP_GET_WTIME ()
      !$omp end parallel
      !print *, t_final - t_ref
      
      !-------------- End loop over the Fourier modes --------------------------------------------------------------

   end subroutine Get_stage_var

   subroutine RHS_construct_stage(Nm_max,Nr_max,Ra,Pr,uphi_temp_FR,ur_temp_FR,uphi_omg_FR,ur_omg_FR,& 
                                  & time_scheme_type,rk_stage,tFR,omgFR,urFR,upFR)

      integer, intent(in) :: Nm_max
      integer, intent(in) :: Nr_max 
      real(kind=dp), intent(in) :: Ra
      real(kind=dp), intent(in) :: Pr  
      character(len=100), intent(in) :: time_scheme_type
      complex(kind=dp), intent(in) :: tFR(Nm_max+1,Nr_max),omgFR(Nm_max+1,Nr_max) 
      complex(kind=dp), intent(in) :: urFR(Nm_max+1,Nr_max),upFR(Nm_max+1,Nr_max) 
      complex(kind=dp), intent(in) :: uphi_temp_FR(Nm_max+1,Nr_max),ur_temp_FR(Nm_max+1,Nr_max) 
      complex(kind=dp), intent(in) :: uphi_omg_FR(Nm_max+1,Nr_max),ur_omg_FR(Nm_max+1,Nr_max)
      integer, intent(in) :: rk_stage
      integer :: Nm, i ! Nm -> azimuthal 'n' loop over Fourier modes
      real(kind=dp) :: t_ref, t_final

      rhs(:)=0.0_dp 
      rhs1(:)=0.0_dp 
      rhsf(:)=0.0_dp 
      rhs_uphi(:)=0.0_dp 
      rhs_vort(:)=0.0_dp
      rhs_b(:)=0.0_dp
      !-------------- Loop over the Fourier modes ---------------------------------------------------------------------------
      !!$omp parallel & 
      !!$omp private(Nm,rhs,rhs_uphi,rhs_vort,rhs_b) default(shared)  
      !   t_ref= OMP_GET_WTIME ()
      !!$omp do    
      do Nm=0,Nm_max  

         !----------- Call RHS construct for temperature advection ----------
         call rhs_construct_temp_a(Nm_max,Nr_max,uphi_temp_FR,ur_temp_FR,tFR(Nm+1,:),Nm,rhs, &
                                 & time_scheme_type)
         
         dtempdt1_a(rk_stage,Nm+1,:) = rhs(:)
         !---------------------------------------------------------

         rhs(:)=0.0_dp 
         !----------- Call RHS construct for temperature diffusion ----------
         call rhs_construct_temp_d(Nm_max,Nr_max,tFR(Nm+1,:),Nm,rhs, &
                                 & time_scheme_type,Pr)
         
         dtempdt1_d(rk_stage,Nm+1,:) = rhs(:)
         !---------------------------------------------------------
         
         if (Nm==0) then

            !----------- Call RHS construct for uphi_bar advection ----------
            call rhs_construct_uphibar_a(Nm_max,Nr_max,urFR,upFR,upFC,rhs_uphi, &
                                    & time_scheme_type,ur,omg)
            
            !duphibar_dt1_a(rk_stage,:) = rhs_uphi(:)
            domgdt1_a(rk_stage,Nm+1,:) = rhs_uphi(:)
            !---------------------------------------------------------

            rhs_uphi(:)=0.0_dp 
            !----------- Call RHS construct for uphi_bar diffusion ----------
            call rhs_construct_uphibar_d(Nm_max,Nr_max,upFR,upFC,rhs_uphi, &
                                    & time_scheme_type)
            
            !duphibar_dt1_d(rk_stage,:) = rhs_uphi(:)
            domgdt1_d(rk_stage,Nm+1,:) = rhs_uphi(:)
            !---------------------------------------------------------

         else

            !----------- Call RHS construct for vorticity advection ------------
            call rhs_construct_vort_a(Nm_max,Nr_max,uphi_omg_FR,ur_omg_FR,omgFR(Nm+1,:), &
                                    & Nm,rhs1,rhsf,time_scheme_type,rhs_vort) 

            domgdt1_a(rk_stage,Nm+1,:) = rhs_vort(:)
            !-------------------------------------------------------------------
            rhs_vort(:)=0.0_dp

            !----------- Call RHS construct for vorticity diffusion ------------
            call rhs_construct_vort_d(Nm_max,Nr_max,omgFR(Nm+1,:), &
                                    & Nm,rhs1,rhsf,time_scheme_type,rhs_vort) 

            domgdt1_d(rk_stage,Nm+1,:) = rhs_vort(:)
            !-------------------------------------------------------------------
            rhs_b(:)=0.0_dp
            if (rk_stage==1) then
               do i=1,Nr_max
                  rhs_buo(rk_stage,Nm+1,i)= - (Ra/Pr)*(ii*real(Nm,kind=dp)*(r_radius(i))*tFR(Nm+1,i))
                          ! buoyancy term
               end do
            end if
         end if
      end do
      !!$omp end do
      !  t_final= OMP_GET_WTIME ()
      !!$omp end parallel
      !!print *, t_final - t_ref

      !-------------- End loop over the Fourier modes --------------------------------------------------------------

   end subroutine RHS_construct_stage

subroutine Assembly_stage(Nm_max,Nr_max,dt,tFR,omgFR,upFR,urFR,lm,mBC)

      integer, intent(in) :: Nm_max
      integer, intent(in) :: Nr_max 
      real(kind=dp), intent(in) :: dt
      integer, intent(in) :: lm  
      character(len=100), intent(in) :: mBC  
      complex(kind=dp), intent(inout) :: tFR(Nm_max+1,Nr_max),omgFR(Nm_max+1,Nr_max) 
      complex(kind=dp), intent(inout) :: upFR(Nm_max+1,Nr_max), urFR(Nm_max+1,Nr_max)

      integer :: i,j,k,Nm,INFO1 ! Nm -> azimuthal 'n' loop over Fourier modes
      complex(kind=dp) :: rhs_psi(Nr_max), F_dtemp(Nr_max), F_duphibar(Nr_max), F_domg(Nr_max)
      complex(kind=dp) :: D1upFR(Nr_max) 
      real(kind=dp) :: t_ref, t_final
      complex(kind=dp) :: dpsi_rmin
      complex(kind=dp) :: d2psi_rmin
      complex(kind=dp) :: dpsi_rmax
      complex(kind=dp) :: d2psi_rmax

         !!-------------- Loop over the Fourier modes ---------------------------------------------------------------------------
         !$omp parallel & 
         !$omp private(Nm,i,j,k,F_dtemp,F_duphibar,F_domg, &
         !$omp & rhs_r,rhs_i,rhsf_r,rhsf_i, &
         !$omp & D1upFR,rhs_psi, &
         !$omp & real_rhs_psi,real_d_rhs_psi,dpsi_rmin,d2psi_rmin,dpsi_rmax,d2psi_rmax) default(shared)  
         t_ref= OMP_GET_WTIME ()
         !$omp do     
         do Nm=0,Nm_max  
            F_dtemp(:)=0.0_dp 
            F_duphibar(:)=0.0_dp 
            F_domg(:)=0.0_dp 
            dpsi_rmin = 0.0_dp
            d2psi_rmin = 0.0_dp
            dpsi_rmax = 0.0_dp
            d2psi_rmax = 0.0_dp

            !-- Apply weights to RHS using Butcher's table: summation of (b_i F_i) + summation of (\hat{b}_i F_i) 
            do i=1,n_order_tscheme_exp
               F_dtemp(:)=F_dtemp(:) + butcher_bA(i)*dtempdt1_a(i,Nm+1,:) + butcher_bD(i)*dtempdt1_d(i,Nm+1,:)
            end do

            !-- y_n = y_(n-1) + dt * summation of (b_i * F_i) -------------------------

            temp_spec(Nm+1,:) = temp_spec(Nm+1,:) + dt*F_dtemp(:)

            if(Nm==0) then ! Apply BC for temperature
               temp_spec(Nm+1,1)=1.0_dp
               temp_spec(Nm+1,Nr_max)=0.0_dp
            else 
               temp_spec(Nm+1,1)=0.0_dp
               temp_spec(Nm+1,Nr_max)=0.0_dp
            end if 

            if (Nm==0) then

               !-- Apply weights to RHS using Butcher's table: summation of (b_i F_i) + summation of (\hat{b}_i F_i) 
               do i=1,n_order_tscheme_exp
                  F_duphibar(:)=F_duphibar(:) + butcher_bA(i)*domgdt1_a(i,Nm+1,:) + butcher_bD(i)*domgdt1_d(i,Nm+1,:)
               end do

                  !F_duphibar(:)=F_duphibar(:) + duphibar_dt1_d(n_order_tscheme_exp,:)
               !-- y_n = y_(n-1) + dt * summation of (b_i * F_i) -------------------------
               uphi_bar_spec(:) = uphi_bar_spec(:) + dt*F_duphibar
               uphi_bar_spec(1)=0.0_dp
               uphi_bar_spec(Nr_max)=0.0_dp
               upFR(1,:) = uphi_bar_spec(:) ! update upFR here

               call chebtransform(Nr_max,uphi_bar_spec,upFC(1,:))
               call chebinvtranD1(Nr_max,upFC(1,:),D1upFR(:))

               urFR(1,:) = 0.0_dp 

               do i=1,Nr_max
                  omg_spec(1,i)=uphi_bar_spec(i)*r_radius(i) + D1upFR(i)
               end do

            else

               !-- Apply weights to RHS using Butcher's table: summation of (b_i F_i) + summation of (\hat{b}_i F_i) 
               do i=1,n_order_tscheme_exp
                  do j=2,Nr_max-1                   
                     F_domg(j)=F_domg(j) + butcher_bA(i)*domgdt1_a(i,Nm+1,j) + butcher_bD(i)*domgdt1_d(i,Nm+1,j) &
                               &         + butcher_bD(i)*rhs_buo(i,Nm+1,j)      
                  end do
               end do

               !-- y_n = y_(n-1) + dt * summation of (b_i * F_i) -------------------------
               do k=2,Nr_max-1                   
                  omg_spec(Nm+1,k) = omg_spec(Nm+1,k) + dt*F_domg(k)
               end do


               rhs_r=-1.0_dp*real(omg_spec(Nm+1,:))
               rhs_i=-1.0_dp*aimag(omg_spec(Nm+1,:))

               rhs_r(1)=0.0_dp  
               rhs_r(Nr_max)=0.0_dp  
               rhs_i(1)=0.0_dp  
               rhs_i(Nr_max)=0.0_dp  

               !************************** CALL DGETRS A*vt=rhs ****************************
               call matsolve(TRANS, Nr_max, LAPpsi_all(:,:,Nm+1), IPIV1_lap(:,Nm+1), rhs_r, rhs_i,INFO1)
               !****************************************************************************  

               do i=1,Nr_max
                 rhs_psi(i)=cmplx(rhs_r(i),rhs_i(i),kind=dp)
               end do
               call chebinvtran(Nr_max,rhs_psi,real_rhs_psi)
               call chebinvtranD1(Nr_max,rhs_psi,real_d_rhs_psi)
               call chebinvtranD2(Nr_max,rhs_psi,real_d2_rhs_psi)

               !-------------------------- Johnston strategy for BC ----------------------------------------------
               do i=1,lm ! Evaluate at rmin
                     d2psi_rmin=d2psi_rmin + d2w_rmin(i)*real_rhs_psi(i+1) ! Summation 2nd derivative
                     dpsi_rmin=dpsi_rmin + dw_rmin(i)*real_rhs_psi(i+1)    ! Summation 1st derivative
                  end do
                  
                  do i=1,lm ! Evaluate at rmax
                     d2psi_rmax=d2psi_rmax + d2w_rmax(i)*real_rhs_psi(Nr_max-i)  ! Summation 2nd derivative
                     dpsi_rmax=dpsi_rmax + dw_rmax(i)*real_rhs_psi(Nr_max-i)     ! Summation 1st derivative
               end do 

               do i=1,Nr_max
                  !psii(Nm+1,i)=real_rhs_psi(i)
                  upFR(Nm+1,i)=-1.0_dp*real_d_rhs_psi(i)
                  urFR(Nm+1,i)=1.0_dp*ii*real(Nm,kind=dp)*r_radius(i)*real_rhs_psi(i)
               end do

               if (mBC=='NS') then ! Johnston's strategy
                  ! ------------ Use Lagrange polynomials for approximation of No-slip boundary conditions ------------------- 
                  omg_spec(Nm+1,1)=-1.0_dp*(d2psi_rmin - 2.0_dp*dw_rmin(1)*dpsi_rmin)            ! Apply omega rmin BC 
                  omg_spec(Nm+1,Nr_max)=-1.0_dp*(d2psi_rmax - 2.0_dp*dw_rmax(1)*dpsi_rmax)       ! Apply omega rmax BC 
                  ! --------------------------------------------------------------------------------------------------
               elseif (mBC=='SF') then
                  ! ------------ Use Lagrange polynomials for approximation of Stress-free boundary conditions ------------------- 
                  omg_spec(Nm+1,1)=-2.0_dp/radius(1)*dpsi_rmin               ! Apply omega rmin BC 
                  omg_spec(Nm+1,Nr_max)=-2.0_dp/radius(Nr_max)*dpsi_rmax     ! Apply omega rmax BC 
                  ! --------------------------------------------------------------------------------------------------
               end if 

               call chebtransform(Nr_max,upFR(Nm+1,:),upFC(Nm+1,:))

            end if
         end do
         !$omp end do
         t_final= OMP_GET_WTIME ()
         !$omp end parallel
         !print *, t_final - t_ref
          
         tFR=temp_spec ! update tFR here
         omgFR=omg_spec ! update omgFR here

   end subroutine Assembly_stage

   subroutine Assembly_stage_SA(Nm_max,Nr_max,tFR,omgFR,upFR)

   integer, intent(in) :: Nm_max
   integer, intent(in) :: Nr_max 
   complex(kind=dp), intent(in) :: tFR(Nm_max+1,Nr_max),omgFR(Nm_max+1,Nr_max) 
   complex(kind=dp), intent(in) :: upFR(Nm_max+1,Nr_max)

   temp_spec = tFR
   omg_spec = omgFR
   uphi_bar_spec = upFR(1,:)

   end subroutine Assembly_stage_SA

end module fourierloop_imexrk
