module fourierloop
   !$ use OMP_LIB
   use double
   use constants, only: ii
   use chebyshev, only: chebtransform, chebinvtran, chebinvtranD1, chebinvtranD1D2, chebinvtranD2 
   use init, only: TFC, tFR, t2FR, omgFR, temp_spec, omg_spec, upFC, psii, upFR, urFR, startmatbuild, &
                   & finishmatbuild, finishmatbuild, time_matbuild, r_radius,dt_new, tmp_rhs_exp_uphi_bar, &
                   & r_radius, r_radius2, omgFR_check
   use mat_assembly, only: AT_all, AF_all, IPIV1, IPIV2, mat_build, mat_build_uphibar, A_uphi_all, IPIV_uphi
   use algebra, only: matsolve, matsolve_real
   use timeschemes, only: wt_lhs_tscheme_imp, wt_rhs_tscheme_imp, wt_rhs_tscheme_exp, n_order_tscheme_imp, &
                          & n_order_tscheme_exp, rhs_update_wts_imp, rhs_update_wts_exp, dt_array, &
                          & rhs_exp_uphi_bar
   use rhs_create, only: rhs_construct_temp, rhs_construct_vort, rhs_construct_uphi_bar   

   implicit none

   private 

   !complex(kind=dp), allocatable :: real_rhs_temp(:)
   !complex(kind=dp), allocatable :: rhs_omg(:), real_rhs_omg(:) 
   !complex(kind=dp), allocatable :: real_rhs_psi(:), real_d_rhs_psi(:) 
   !real(kind=dp), allocatable :: rhs_r(:), rhs_i(:), rhsf_r(:), rhsf_i(:)
   character :: TRANS='N'

   public :: allocate_fourierloop_imex, deallocate_fourierloop_imex, Nm_maxLOOP

contains

   subroutine allocate_fourierloop_imex(Nr_max)

      integer, intent(in) :: Nr_max 
        
      !allocate( real_rhs_temp(Nr_max) )
      !allocate( rhs_omg(Nr_max),real_rhs_omg(Nr_max),real_rhs_psi(Nr_max),real_d_rhs_psi(Nr_max) )
      !allocate( rhs_r(Nr_max), rhs_i(Nr_max), rhsf_r(2*Nr_max), rhsf_i(2*Nr_max) )  
   
   end subroutine allocate_fourierloop_imex

   subroutine deallocate_fourierloop_imex
      
      !deallocate( rhs_r, rhs_i, rhsf_r, rhsf_i )  
      !deallocate( rhs_omg,real_rhs_omg,real_rhs_psi,real_d_rhs_psi )
      !deallocate( real_rhs_temp )
    
   end subroutine deallocate_fourierloop_imex

   subroutine Nm_maxLOOP(Nm_max,Nr_max,Ra,Pr,mBC,uphi_temp_FR,ur_temp_FR,uphi_omg_FR,ur_omg_FR,n_step,n_restart,& 
                     & time_scheme_imp, time_scheme_exp,l_restart)

      integer, intent(in) :: Nm_max
      integer, intent(in) :: Nr_max 
      integer, intent(in) :: n_step, n_restart
      real(kind=dp), intent(in) :: Ra
      real(kind=dp), intent(in) :: Pr  
      logical, intent(in) :: l_restart
      character(len=100), intent(in) :: mBC  
      character(len=100), intent(in) :: time_scheme_imp
      character(len=100), intent(in) :: time_scheme_exp
      complex(kind=dp), intent(in) :: uphi_temp_FR(Nm_max+1,Nr_max),ur_temp_FR(Nm_max+1,Nr_max) 
      complex(kind=dp), intent(in) :: uphi_omg_FR(Nm_max+1,Nr_max),ur_omg_FR(Nm_max+1,Nr_max)
      !Local variables -----------------------------  
      complex(kind=dp) :: D1upFR(Nr_max)
      complex(kind=dp) :: rhs(Nr_max),rhsf(2*Nr_max),rhs_psi(Nr_max),real_rhs_temp(Nr_max)
      real(kind=dp) :: rhs_uphi(Nr_max)
      real(kind=dp) :: rhs_r(Nr_max), rhs_i(Nr_max), rhsf_r(2*Nr_max), rhsf_i(2*Nr_max)
      complex(kind=dp) :: rhs_omg(Nr_max),real_rhs_omg(Nr_max),real_rhs_psi(Nr_max),real_d_rhs_psi(Nr_max)
      complex(kind=dp) :: real_d2_rhs_psi(Nr_max)
      complex(kind=dp) :: omgFR_check(Nr_max)
      integer :: i,Nm,INFO1,INFO2,Nr_max2 ! Nm -> azimuthal 'n' loop over Fourier modes
      real(kind=dp) :: t_ref, t_final

      ! Update RHS weights ----------
      call rhs_update_wts_imp(time_scheme_imp,wt_lhs_tscheme_imp,wt_rhs_tscheme_imp,n_order_tscheme_imp)
      call rhs_update_wts_exp(time_scheme_imp,time_scheme_exp,wt_rhs_tscheme_exp,n_order_tscheme_exp)
      !------------------------------
      !----------------------------------------------------------------------------------------------
      if (n_step-n_restart==1) then ! For working with variables in Fourier-Chebyshev space 
         do i=0,Nm_max
            call chebtransform(Nr_max,tFR(i+1,:),temp_spec(i+1,:))
            call chebtransform(Nr_max,omgFR(i+1,:),omg_spec(i+1,:))
            call chebtransform(Nr_max,upFR(i+1,:),upFC(i+1,:))
         end do
      end if !---------------------------------------------------------------------------------------
      
      Nr_max2=2*Nr_max
      rhs(:)=0.0_dp
      rhs_r(:)=0.0_dp
      rhs_i(:)=0.0_dp
      rhs_uphi(:)=0.0_dp
!-------------- Loop over the Fourier modes ---------------------------------------------------------------------------
      !$omp parallel & 
      !$omp private(Nm,rhs,rhs_r,rhs_i,rhsf_r,rhsf_i,real_rhs_temp,rhs_uphi,D1upFR,real_rhs_omg,real_rhs_psi, & 
      !$omp & real_d_rhs_psi,rhsf,rhs_omg,rhs_psi,rhs_exp_uphi_bar) default(shared)  
         t_ref= OMP_GET_WTIME ()
      !$omp do   
      do Nm=0,Nm_max   
           
         if (n_step-n_restart==1 .or. dt_array(1)/=dt_array(2)) then ! Call 'mat_build' only if dt_current and dt_previous are different 
            call cpu_time(startmatbuild)
            call mat_build(Nr_max,dt_array(1),Nm,mBC,wt_lhs_tscheme_imp,Pr) ! Build the operator matrices and factorize them 
            call mat_build_uphibar(Nr_max,dt_array(1),mBC,wt_lhs_tscheme_imp,Pr)
            call cpu_time(finishmatbuild) 
            time_matbuild=time_matbuild + finishmatbuild - startmatbuild
         end if
      
! -------------------------------------------------------------------------------------------------

         !----------- Call RHS construct for temperature ----------
         call rhs_construct_temp(Nm_max,Nr_max,dt_new,uphi_temp_FR,ur_temp_FR,n_step,temp_spec(Nm+1,:),Nm,rhs,n_restart, &
                                & wt_rhs_tscheme_imp, wt_rhs_tscheme_exp,n_order_tscheme_imp, &
                                 & n_order_tscheme_exp, time_scheme_imp,Pr,l_restart)
         !---------------------------------------------------------  
         
         rhs_r(:)=real(rhs(:))
         rhs_i(:)=aimag(rhs(:))

         !************************** CALL DGETRS A*vt=rhs ****************************
         call matsolve(TRANS, Nr_max, AT_all(:,:,Nm+1), IPIV1(:,Nm+1), rhs_r, rhs_i,INFO1)
         !**************************************************************************** 
         do i=1,Nr_max
           rhs(i)=cmplx(rhs_r(i),rhs_i(i),kind=dp)
         end do
         
         if (Nm==0) then
            TFC=rhs
         end if

         call chebinvtran(Nr_max,rhs,real_rhs_temp)
         temp_spec(Nm+1,:)=rhs(:)
         tFR(Nm+1,:)=real_rhs_temp(:)

         ! SOLVE FOR UPHI_BAR HERE -----------------------------------------------------------------------------
         if (Nm==0) then
            !----------- Call RHS construct for uphi_bar -------------
            call rhs_construct_uphi_bar(Nm_max,Nr_max,dt_new,upFR,urFR,omgFR,n_step,upFC,rhs_uphi,n_restart, &
                                       & wt_rhs_tscheme_imp,wt_rhs_tscheme_exp,n_order_tscheme_imp, &
                                       & n_order_tscheme_exp,time_scheme_imp,l_restart)
            !---------------------------------------------------------   
            rhs_r=rhs_uphi
            rhs_i=0.0_dp
            ! --------------------------------------------------
               
            !************************** CALL DGETRS A_uphi*u_phi=rhs ****************************
            call matsolve_real(TRANS, Nr_max, A_uphi_all(:,:,Nm+1), IPIV_uphi(:,Nm+1), rhs_r, INFO1)
            !**************************************************************************** 

            do i=1,Nr_max
               upFC(Nm+1,i)=cmplx(rhs_r(i),rhs_i(i),kind=dp)
            end do

            call chebinvtran(Nr_max,upFC(Nm+1,:),upFR(Nm+1,:))
            call chebinvtranD1(Nr_max,upFC(Nm+1,:),D1upFR(:))
            urFR(Nm+1,:) = 0.0_dp 

            do i=1,Nr_max
               omgFR(Nm+1,i)=upFR(Nm+1,i)*r_radius(i) + D1upFR(i)
            end do
         ! ----------------------------------------------------------------------------------------------------   
         else
            !----------- Call RHS construct for vorticity ----------
            call rhs_construct_vort(Nm_max,Nr_max,dt_new,Ra,Pr,uphi_omg_FR,ur_omg_FR,n_step,omg_spec(Nm+1,:), &
                                    & Nm,rhsf,n_restart,wt_rhs_tscheme_imp, &
                                    & wt_rhs_tscheme_exp,n_order_tscheme_imp,n_order_tscheme_exp, &
                                    & time_scheme_imp,tFR(Nm+1,:),l_restart)
       
            !------------------------------------------------------- 

            rhsf_r=real(rhsf)
            rhsf_i=aimag(rhsf)

            !**************************** CALL DGETRS AF*vf=rhsf **************************
            call matsolve(TRANS, Nr_max2, AF_all(:,:,Nm+1), IPIV2(:,Nm+1), rhsf_r, rhsf_i,INFO2)
            !******************************************************************************
            do i=1,2*Nr_max
              rhsf(i)=cmplx(rhsf_r(i),rhsf_i(i),kind=dp)
            end do

            do i=1,Nr_max

               omg_spec(Nm+1,i)=rhsf(Nr_max+i)
               rhs_omg(i)=rhsf(Nr_max+i)
               rhs_psi(i)=rhsf(i)

            end do

            call chebinvtran(Nr_max,rhs_omg,real_rhs_omg)
            omgFR(Nm+1,:)=real_rhs_omg(:)
            call chebinvtran(Nr_max,rhs_psi,real_rhs_psi)
            call chebinvtranD1(Nr_max,rhs_psi,real_d_rhs_psi)
            call chebinvtranD2(Nr_max,rhs_psi,real_d2_rhs_psi)

            do i=1,Nr_max
               psii(Nm+1,i)=real_rhs_psi(i)
               upFR(Nm+1,i)=-1.0_dp*real_d_rhs_psi(i)
               urFR(Nm+1,i)=ii*real(Nm,kind=dp)*r_radius(i)*real_rhs_psi(i)
            end do

            ! OMEGA FIX - CHECK start
            !omgFR(Nm+1,1) = r_radius(1)*real_d_rhs_psi(1)+real_d2_rhs_psi(1)- &
            !                & real(Nm,kind=dp)*real(Nm,kind=dp)*r_radius2(1)*real_rhs_psi(1)
            !omgFR(Nm+1,Nr_max) = r_radius(Nr_max)*real_d_rhs_psi(Nr_max)+real_d2_rhs_psi(Nr_max)- &
            !                & real(Nm,kind=dp)*real(Nm,kind=dp)*r_radius2(Nr_max)*real_rhs_psi(Nr_max)
            !omgFR_check(1) = -1.0_dp*(real_d2_rhs_psi(1))
            !omgFR_check(Nr_max) = -1.0_dp*(real_d2_rhs_psi(Nr_max))
            !do i=1,Nr_max
            !   omgFR_check(i) = r_radius(i)*real_d_rhs_psi(i)+real_d2_rhs_psi(i)- &
            !                    & real(Nm,kind=dp)*real(Nm,kind=dp)*r_radius2(i)*real_rhs_psi(i)
            !end do
            !if (n_step==1000) then 
            !print *, (abs(aimag(omgFR_check(1)))), (abs(aimag(omgFR(Nm+1,1)))), & 
            !         & abs(aimag(omgFR_check(1)))-abs(aimag(omgFR(Nm+1,1))),   "comparison", Nm  
            !end if
            !print *, real_d2_rhs_psi(5), omgFR(Nm+1,5), "comparison"  
            ! OMEGA FIX - CHECK end
            call chebtransform(Nr_max,upFR(Nm+1,:),upFC(Nm+1,:))

         end if
                  
      end do  
      !$omp end do
         t_final= OMP_GET_WTIME ()
      !$omp end parallel
      !print *, maxval(aimag(omgFR_check))
        print *, t_final - t_ref
!-------------- End loop over the Fourier modes -------------------------------------------------------------------------

   end subroutine Nm_maxLOOP

end module fourierloop

!SHARED(uphi_temp_FR,ur_temp_FR,uphi_omg_FR,ur_omg_FR, &
!            !$OMP & Nm_max,Nr_max,n_step,n_restart,wt_rhs_tscheme_imp,wt_rhs_tscheme_exp,n_order_tscheme_imp, &
!            !$OMP & n_order_tscheme_exp, dt_array, &
!            !$OMP & mBC,Ra,Pr,AT_all,A_uphi,AF_all,PIV_uphi,IPIV1,IPIV2,INFO1,INFO2,dt_new) 
