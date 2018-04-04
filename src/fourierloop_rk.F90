module fourierloop_rk

   use double
   use constants, only: ii
   use chebyshev, only: chebtransform, chebinvtran, chebinvtranD1, chebinvtranD2, chebinvtranD1D2, t, D, D2
   use init, only: TFC, temp_spec, omg_spec, uphi_bar_spec, startmatbuild, &
                   & finishmatbuild, finishmatbuild, time_matbuild, r_radius, r_radius2, &
                   & dw_rmin, d2w_rmin, dw_rmax, d2w_rmax, upFC
   use mat_assembly, only: LAPpsi_all, IPIV1_lap
   use algebra, only: matsolve
   use timeschemes, only: wt_rhs_tscheme_exp,n_order_tscheme_exp, rhs_update_wts_exp, irk_max, &
                          & dt_array, butcher_a, butcher_b, rhs_exp_temp, rhs_exp_vort, &
                          & rhs_exp_uphi_bar 
   use rhs_create_rk, only: rhs_construct_temp, rhs_construct_vort, rhs_construct_uphi_bar   

   implicit none

   private 

   complex(kind=dp), allocatable :: real_rhs_temp(:)
   complex(kind=dp), allocatable :: rhs_omg(:), rhs_psi(:), real_rhs_omg(:)
   real(kind=dp), allocatable :: rhs_uphi(:)
   complex(kind=dp), allocatable :: real_rhs_psi(:), real_d_rhs_psi(:), real_d2_rhs_psi(:) 
   complex(kind=dp), allocatable :: rhs1(:), rhs(:), rhsf(:), rhs_vort(:)
   real(kind=dp), allocatable :: rhs_r(:), rhs_i(:), rhf_r(:), rhf_i(:)

   complex(kind=dp), allocatable, public :: dtempdt1(:,:,:), dtempdt2(:,:), domgdt1(:,:,:), duphibar_dt1(:,:)
   character :: TRANS='N'

   complex(kind=dp) :: dpsi_rmin
   complex(kind=dp) :: d2psi_rmin
   complex(kind=dp) :: dpsi_rmax
   complex(kind=dp) :: d2psi_rmax

   public :: allocate_fourierloop_rk, deallocate_fourierloop_rk, Get_stage_var, &
             & RHS_construct_stage, Assembly_stage

contains

   subroutine allocate_fourierloop_rk(Nm_max,Nr_max)

      integer, intent(in) :: Nr_max 
      integer, intent(in) :: Nm_max 

      allocate( real_rhs_temp(Nr_max) )
      allocate( rhs_omg(Nr_max),rhs_psi(Nr_max),real_rhs_omg(Nr_max),real_rhs_psi(Nr_max) )
      allocate( rhs_uphi(Nr_max) )
      allocate( real_d_rhs_psi(Nr_max), real_d2_rhs_psi(Nr_max) )
      allocate( rhs1(Nr_max), rhs(Nr_max), rhsf(2*Nr_max), rhs_vort(Nr_max) )  
      allocate( rhs_r(Nr_max), rhs_i(Nr_max), rhf_r(2*Nr_max), rhf_i(2*Nr_max) )
      allocate( dtempdt1(n_order_tscheme_exp,Nm_max+1,Nr_max),dtempdt2(Nm_max+1,Nr_max), &
                & domgdt1(n_order_tscheme_exp,Nm_max+1,Nr_max), &
                & duphibar_dt1(n_order_tscheme_exp,Nr_max) )
   
   end subroutine allocate_fourierloop_rk

   subroutine deallocate_fourierloop_rk
   
      deallocate( dtempdt1, dtempdt2, domgdt1, duphibar_dt1 )
      deallocate( rhs_r, rhs_i, rhf_r, rhf_i )  
      deallocate( rhs1, rhs,rhsf, rhs_vort )  
      deallocate( rhs_uphi )
      deallocate( rhs_omg,rhs_psi,real_rhs_omg,real_rhs_psi,real_d_rhs_psi, real_d2_rhs_psi )
      deallocate( real_rhs_temp )
    
   end subroutine deallocate_fourierloop_rk

   subroutine Get_stage_var(Nm_max,Nr_max,time_scheme_imp,time_scheme_exp,rk_stage,lm,tFR,omgFR,psii,upFR,urFR)

      integer, intent(in) :: lm 
      integer, intent(in) :: Nm_max
      integer, intent(in) :: Nr_max 
      character(len=100), intent(in) :: time_scheme_imp
      character(len=100), intent(in) :: time_scheme_exp
      integer, intent(in) :: rk_stage
      complex(kind=dp), intent(out) :: tFR(Nm_max+1,Nr_max),omgFR(Nm_max+1,Nr_max) 
      complex(kind=dp), intent(out) :: psii(Nm_max+1,Nr_max)
      complex(kind=dp), intent(out) :: upFR(Nm_max+1,Nr_max), urFR(Nm_max+1,Nr_max)

      complex(kind=dp) :: D1upFR(Nr_max)
      integer :: i,Nm,INFO1 ! Nm -> azimuthal 'n' loop over Fourier modes
      complex(kind=dp) :: rhs_psi(Nr_max), F_dtemp(Nr_max), F_domg(Nr_max), F_duphi_bar(Nr_max)

         call rhs_update_wts_exp(time_scheme_imp,time_scheme_exp,wt_rhs_tscheme_exp,n_order_tscheme_exp)
     
!-------------- Loop over the Fourier modes ---------------------------------------------------------------------------
         do Nm=0,Nm_max   
            F_dtemp(:)=0.0_dp 
            F_duphi_bar(:)=0.0_dp 
            F_domg(:)=0.0_dp  
            dpsi_rmin = 0.0_dp
            d2psi_rmin = 0.0_dp
            dpsi_rmax = 0.0_dp
            d2psi_rmax = 0.0_dp

            !-- Apply weights to RHS using Butcher's table: summation of (a_i F_i) 
            do i=1,rk_stage-1 
               F_dtemp(:)=F_dtemp(:)+butcher_a(rk_stage,i)*dtempdt1(i,Nm+1,:) 
            end do
            
            !-- Y_stage = Y_0 + dt * summation of (a_i * F_i) -------------------------
            tFR(Nm+1,:) = temp_spec(Nm+1,:) + dt_array(1)*F_dtemp

            if(Nm==0) then ! Apply BC for temperature
                 tFR(Nm+1,1)=1.0_dp
                 tFR(Nm+1,Nr_max)=0.0_dp
            else 
                 tFR(Nm+1,1)=0.0_dp
                 tFR(Nm+1,Nr_max)=0.0_dp
            end if

            if (rk_stage==n_order_tscheme_exp) then
               rhs_exp_temp(1,Nm+1,:) = dt_array(1)*F_dtemp
            end if
            
            if (Nm==0) then

               !-- Apply weights to RHS using Butcher's table: summation of (a_i F_i) 
               do i=1,rk_stage-1 
                  F_duphi_bar(:)=F_duphi_bar(:)+butcher_a(rk_stage,i)*duphibar_dt1(i,:) 
               end do
               
               !-- Y_stage = Y_0 + dt * summation of (a_i * F_i) -------------------------
               upFR(Nm+1,:) = uphi_bar_spec(:) + dt_array(1)*F_duphi_bar

               ! Apply BC for uphi_bar
               upFR(Nm+1,1)=0.0_dp
               upFR(Nm+1,Nr_max)=0.0_dp

               call chebtransform(Nr_max,upFR(Nm+1,:),upFC(Nm+1,:))
               call chebinvtranD1(Nr_max,upFC(Nm+1,:),D1upFR(:))

               urFR(Nm+1,:) = 0.0_dp 

               do i=1,Nr_max
                  omgFR(Nm+1,i)=upFR(Nm+1,i)*r_radius(i) + D1upFR(i)
               end do

               if (rk_stage==n_order_tscheme_exp) then
                  rhs_exp_uphi_bar(1,:) = dt_array(1)*F_duphi_bar
               end if

            else

               !-- Apply weights to RHS using Butcher's table: summation of (a_i F_i) 
               do i=1,rk_stage-1
                  F_domg(:)=F_domg(:)+butcher_a(rk_stage,i)*domgdt1(i,Nm+1,:)
               end do

               !-- Y_stage = Y_0 + dt * summation of (a_i * F_i) -------------------------
               omgFR(Nm+1,:) = omg_spec(Nm+1,:) + dt_array(1)*F_domg
               
               rhs_r=-1.0_dp*real(omgFR(Nm+1,:))
               rhs_i=-1.0_dp*aimag(omgFR(Nm+1,:))

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
                     d2psi_rmax=d2psi_rmax + d2w_rmax(i)*real_rhs_psi(Nr_max-i) ! Summation 2nd derivative
                     dpsi_rmax=dpsi_rmax + dw_rmax(i)*real_rhs_psi(Nr_max-i)    ! Summation 1st derivative
               end do 

               ! ------------ Use Lagrange polynomials for approximation of No₋slip boundary conditions -----------
               omgFR(Nm+1,1)=-1.0_dp*(d2psi_rmin - 2.0_dp*dw_rmin(1)*dpsi_rmin)            ! Apply omega BC ! Johnston strategy
               omgFR(Nm+1,Nr_max)=-1.0_dp*(d2psi_rmax - 2.0_dp*dw_rmax(1)*dpsi_rmax)       ! Apply omega BC ! Johnston strategy
               ! --------------------------------------------------------------------------------------------------
               !---------------------------------------------------------------------------------------------------

               do i=1,Nr_max
                  psii(Nm+1,i)=real_rhs_psi(i)
                  upFR(Nm+1,i)=-1.0_dp*real_d_rhs_psi(i)
                  urFR(Nm+1,i)=ii*real(Nm,kind=dp)*r_radius(i)*real_rhs_psi(i)
               end do

               call chebtransform(Nr_max,upFR(Nm+1,:),upFC(Nm+1,:))

               if (rk_stage==n_order_tscheme_exp) then
                  rhs_exp_vort(1,Nm+1,:) = dt_array(1)*F_domg
               end if

            end if

         end do

       
!-------------- End loop over the Fourier modes --------------------------------------------------------------

   end subroutine Get_stage_var

   subroutine RHS_construct_stage(Nm_max,Nr_max,Ra,Pr,uphi_temp_FR,ur_temp_FR,uphi_omg_FR,ur_omg_FR,& 
                                  & time_scheme_exp,rk_stage,urFR,upFR,tFR,omgFR)

      integer, intent(in) :: Nm_max
      integer, intent(in) :: Nr_max 
      real(kind=dp), intent(in) :: Ra
      real(kind=dp), intent(in) :: Pr  
      character(len=100), intent(in) :: time_scheme_exp
      complex(kind=dp), intent(in) :: urFR(Nm_max+1,Nr_max),upFR(Nm_max+1,Nr_max) 
      complex(kind=dp), intent(in) :: tFR(Nm_max+1,Nr_max),omgFR(Nm_max+1,Nr_max) 
      complex(kind=dp), intent(in) :: uphi_temp_FR(Nm_max+1,Nr_max),ur_temp_FR(Nm_max+1,Nr_max) 
      complex(kind=dp), intent(in) :: uphi_omg_FR(Nm_max+1,Nr_max),ur_omg_FR(Nm_max+1,Nr_max)
      integer, intent(in) :: rk_stage
      integer :: Nm ! Nm -> azimuthal 'n' loop over Fourier modes

      rhs(:)=0.0_dp 
      rhs_uphi(:)=0.0_dp
      rhs_vort(:)=0.0_dp
!-------------- Loop over the Fourier modes ---------------------------------------------------------------------------

      do Nm=0,Nm_max  

         !----------- Call RHS construct for temperature ----------
         call rhs_construct_temp(Nm_max,Nr_max,uphi_temp_FR,ur_temp_FR,tFR,Nm,rhs, &
                                 & time_scheme_exp,Pr)
         
         dtempdt1(rk_stage,Nm+1,:) = rhs(:)
         !---------------------------------------------------------
         
         if (Nm==0) then

            !----------- Call RHS construct for uphi_bar ----------
            call rhs_construct_uphi_bar(Nm_max,Nr_max,upFR,urFR,upFC,omgFR,rhs_uphi, &
                                    & time_scheme_exp)
            
            duphibar_dt1(rk_stage,:) = rhs_uphi(:)
            !---------------------------------------------------------
         
         else

            !----------- Call RHS construct for vorticity ------------
            call rhs_construct_vort(Nm_max,Nr_max,Ra,Pr,uphi_omg_FR,ur_omg_FR,omgFR, &
                                    & Nm,rhs1,rhsf,time_scheme_exp,rhs_vort) 

            domgdt1(rk_stage,Nm+1,:) = rhs_vort(:)

         end if 

      end do

!-------------- End loop over the Fourier modes --------------------------------------------------------------

   end subroutine RHS_construct_stage 

   subroutine Assembly_stage(Nm_max,Nr_max,dt,lm,tFR,omgFR,psii,upFR,urFR)

      integer, intent(in) :: lm 
      integer, intent(in) :: Nm_max
      integer, intent(in) :: Nr_max 
      real(kind=dp), intent(in) :: dt
      complex(kind=dp), intent(inout) :: tFR(Nm_max+1,Nr_max),omgFR(Nm_max+1,Nr_max) 
      complex(kind=dp), intent(inout) :: psii(Nm_max+1,Nr_max)
      complex(kind=dp), intent(out) :: upFR(Nm_max+1,Nr_max), urFR(Nm_max+1,Nr_max)

      integer :: i,Nm,INFO1 ! Nm -> azimuthal 'n' loop over Fourier modes
      complex(kind=dp) :: D1upFR(Nr_max)
      complex(kind=dp) :: rhs_psi(Nr_max), F_dtemp(Nr_max), F_domg(Nr_max), F_duphi_bar(Nr_max)
      

!-------------- Loop over the Fourier modes ---------------------------------------------------------------------------
      do Nm=0,Nm_max  
         F_dtemp(:)=0.0_dp 
         F_duphi_bar(:)=0.0_dp 
         F_domg(:)=0.0_dp 
         dpsi_rmin = 0.0_dp
         d2psi_rmin = 0.0_dp
         dpsi_rmax = 0.0_dp
         d2psi_rmax = 0.0_dp

         !-- Apply weights to RHS using Butcher's table: summation of (b_i F_i) 
         do i=1,n_order_tscheme_exp
            F_dtemp(:)=F_dtemp(:)+butcher_b(i)*dtempdt1(i,Nm+1,:)
         end do

         !-- y_n = y_(n-1) + dt * summation of (b_i * F_i) -------------------------
         temp_spec(Nm+1,:) = temp_spec(Nm+1,:) + dt*F_dtemp

         if(Nm==0) then ! Apply BC for temperature
              temp_spec(Nm+1,1)=1.0_dp
              temp_spec(Nm+1,Nr_max)=0.0_dp
         else 
              temp_spec(Nm+1,1)=0.0_dp
              temp_spec(Nm+1,Nr_max)=0.0_dp
         end if

         if (Nm==0) then

            !-- Apply weights to RHS using Butcher's table: summation of (b_i F_i) 
            do i=1,n_order_tscheme_exp
               F_duphi_bar(:)=F_duphi_bar(:)+butcher_b(i)*duphibar_dt1(i,:)
            end do

            !-- y_n = y_(n-1) + dt * summation of (b_i * F_i) -------------------------
            uphi_bar_spec(:) = uphi_bar_spec(:) + dt*F_duphi_bar

            ! Apply BC for uphi_bar 
            uphi_bar_spec(1)=0.0_dp
            uphi_bar_spec(Nr_max)=0.0_dp

            call chebtransform(Nr_max,uphi_bar_spec,upFC(Nm+1,:))
            call chebinvtranD1(Nr_max,upFC(Nm+1,:),D1upFR(:))

            urFR(Nm+1,:) = 0.0_dp 

            do i=1,Nr_max
               omg_spec(Nm+1,i)=uphi_bar_spec(i)*r_radius(i) + D1upFR(i)
            end do

         else

            !-- Apply weights to RHS using Butcher's table: summation of (b_i F_i) 
            do i=1,n_order_tscheme_exp
               F_domg(:)=F_domg(:)+butcher_b(i)*domgdt1(i,Nm+1,:)
            end do

            !-- y_n = y_(n-1) + dt * summation of (b_i * F_i) -------------------------
            omg_spec(Nm+1,:) = omg_spec(Nm+1,:) + dt*F_domg

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

            ! ------------ Use Lagrange polynomials for approximation of No-slip boundary conditions ------------------- 
            omg_spec(Nm+1,1)=-1.0_dp*(d2psi_rmin - 2.0_dp*dw_rmin(1)*dpsi_rmin)            ! Apply omega BC ! Johnston strategy
            omg_spec(Nm+1,Nr_max)=-1.0_dp*(d2psi_rmax - 2.0_dp*dw_rmax(1)*dpsi_rmax)       ! Apply omega BC ! Johnston strategy
            ! --------------------------------------------------------------------------------------------------
            !---------------------------------------------------------------------------------------------------- 

            do i=1,Nr_max
               psii(Nm+1,i)=real_rhs_psi(i)
               upFR(Nm+1,i)=-1.0_dp*real_d_rhs_psi(i)
               urFR(Nm+1,i)=ii*real(Nm,kind=dp)*r_radius(i)*real_rhs_psi(i)
            end do

            call chebtransform(Nr_max,upFR(Nm+1,:),upFC(Nm+1,:))

            end if

      end do

      tFR=temp_spec
      omgFR=omg_spec
      upFR(1,:)=uphi_bar_spec(:) 
!-------------- End loop over the Fourier modes --------------------------------------------------------------

   end subroutine Assembly_stage

end module fourierloop_rk
