module rhs_create
   !$ use omp_lib
   use double
   use constants, only: ii, pi
   use fourier, only: invfft
   use chebyshev, only: chebtransform, chebinvtran, chebinvtranD1, chebinvtranD2, chebinvtranD1D2
   use init, only: r_radius, r_radius2, tmp_rhs_imp_temp, tmp_rhs_exp_temp, &
                   & tmp_rhs_imp_vort, tmp_rhs_exp_vort, tmp_rhs_imp_uphi_bar, tmp_rhs_exp_uphi_bar, &
                   & tmp_rhs_buo_term, tFR,omgFR, ur, up, omg 
   use timeschemes, only:  wt_lhs_tscheme_imp, wt_rhs_tscheme_imp, rhs_exp_vort, &
                          & rhs_imp_temp, rhs_imp_vort, rhs_imp_uphi_bar, rhs_exp_uphi_bar, &
                          & dt_array, rhs_buo_term, rhs_exp_temp

   implicit none

   private

   !complex(kind=dp), allocatable :: uphi_temp_rad(:), ur_temp_rad_spec(:), ur_temp_rad(:)
   !complex(kind=dp), allocatable :: uphi_omg_rad(:), ur_omg_rad_spec(:), ur_omg_rad(:) 
   !complex(kind=dp), allocatable :: real_d_omg_rad(:), real_d2_omg_rad(:) 
   !complex(kind=dp), allocatable :: omg_real_rad(:)  
   !complex(kind=dp), allocatable :: real_temp_rad(:), real_d_temp_rad(:), real_d2_temp_rad(:)
   !complex(kind=dp), allocatable :: real_d_ur_temp_rad(:,:), real_d_ur_omg_rad(:,:)
   !complex(kind=dp), allocatable, public :: rhs2(:)
   public :: allocate_rhs_imex, deallocate_rhs_imex, rhs_construct_temp, rhs_construct_vort, &
             & rhs_construct_uphi_bar, rhs_temp_from_restart, rhs_vort_from_restart, &
             & rhs_uphibar_from_restart

contains

   subroutine allocate_rhs_imex(Nm_max,Nr_max)

      integer, intent(in) :: Nm_max
      integer, intent(in) :: Nr_max 
           
      !allocate( uphi_temp_rad(Nr_max), ur_temp_rad_spec(Nr_max), ur_temp_rad(Nr_max), &
      !          & uphi_omg_rad(Nr_max), ur_omg_rad_spec(Nr_max) )
      !allocate( omg_real_rad(Nr_max), &
      !          & ur_omg_rad(Nr_max), real_d_omg_rad(Nr_max), real_d2_omg_rad(Nr_max) ) 
      !allocate( real_temp_rad(Nr_max), real_d_temp_rad(Nr_max), real_d2_temp_rad(Nr_max) )
      !allocate( real_d_ur_temp_rad(Nm_max+1,Nr_max), real_d_ur_omg_rad(Nm_max+1,Nr_max) )
      !allocate( rhs2(Nr_max) )

   end subroutine allocate_rhs_imex

   subroutine deallocate_rhs_imex

      !deallocate( rhs2 )
      !deallocate( real_d_ur_temp_rad, real_d_ur_omg_rad )
      !deallocate( real_temp_rad, real_d_temp_rad, real_d2_temp_rad )
      !deallocate( omg_real_rad, ur_omg_rad, real_d_omg_rad, real_d2_omg_rad )
      !deallocate( uphi_temp_rad, ur_temp_rad_spec, ur_temp_rad, uphi_omg_rad, ur_omg_rad_spec )

   end subroutine deallocate_rhs_imex

   subroutine rhs_construct_temp(Nm_max,Nr_max,dt,uphi_temp_FR,ur_temp_FR,n_step,temp_specp,Nm,rhs_temp,n_restart, &
                                 & wt_rhs_tscheme_imp, wt_rhs_tscheme_exp,n_order_tscheme_imp, &
                                 & n_order_tscheme_exp, time_scheme_imp,Pr)
             
      integer, intent(in) :: Nm_max
      integer, intent(in) :: Nr_max 
      real(kind=dp), intent(in) :: Pr  
      integer, intent(in) :: n_step, n_restart
      integer, intent(in) :: n_order_tscheme_imp
      integer, intent(in) :: n_order_tscheme_exp
      real(kind=dp), intent(in) :: wt_rhs_tscheme_imp(n_order_tscheme_imp)
      real(kind=dp), intent(in) :: wt_rhs_tscheme_exp(n_order_tscheme_exp)
      !complex(kind=dp), intent(in) :: tmp_exp_temp(n_order_tscheme_exp,Nr_max)
      character(len=100), intent(in) :: time_scheme_imp
      real(kind=dp), intent(in) :: dt
      complex(kind=dp), intent(in) :: temp_specp(Nr_max)
      integer, intent(in) :: Nm
      complex(kind=dp), intent(in) :: uphi_temp_FR(Nm_max+1,Nr_max),ur_temp_FR(Nm_max+1,Nr_max)
      complex(kind=dp), intent(out) :: rhs_temp(Nr_max)
      ! Local variables --------------------------
      integer :: i, i_order
      complex(kind=dp) :: ur_temp_rad(Nr_max)
      complex(kind=dp) :: real_d_ur_temp_rad(Nr_max)
      complex(kind=dp) :: ur_temp_rad_spec(Nr_max)
      complex(kind=dp) :: uphi_temp_rad(Nr_max)
      complex(kind=dp) :: real_temp_rad(Nr_max)
      complex(kind=dp) :: real_d_temp_rad(Nr_max)
      complex(kind=dp) :: real_d2_temp_rad(Nr_max)
     
      if (n_step-n_restart>1 .or. n_restart/=0) then
         do i_order=1,n_order_tscheme_exp-1
            rhs_exp_temp(i_order+1,Nm+1,:) = tmp_rhs_exp_temp(i_order,Nm+1,:)
         end do
         do i_order=1,n_order_tscheme_imp-1
            rhs_imp_temp(i_order+1,Nm+1,:) = tmp_rhs_imp_temp(i_order,Nm+1,:)
         end do
      end if
      ! CANT MAKE A COPY HERE IF THE VARIABLE temp_specp is private
      do i=1,Nr_max
         uphi_temp_rad(i)=uphi_temp_FR(Nm+1,i)
         ur_temp_rad(i)=ur_temp_FR(Nm+1,i)
         rhs_temp(i)=0.0_dp
         real_temp_rad(i)=0.0_dp
      end do

      call chebtransform(Nr_max,ur_temp_rad,ur_temp_rad_spec)
      
      call chebinvtran(Nr_max,temp_specp,real_temp_rad)
      
      call chebinvtranD1D2(Nr_max,temp_specp,real_d_temp_rad,real_d2_temp_rad)
      
      call chebinvtranD1(Nr_max,ur_temp_rad_spec,real_d_ur_temp_rad)

      !----- If Crank-Nicholson (CN) is used as the implicit scheme -------------------------
      if (time_scheme_imp=='CN') then 

         do i=1,Nr_max

            rhs_imp_temp(1,Nm+1,i)=real_temp_rad(i)

            rhs_temp(i) = rhs_imp_temp(1,Nm+1,i)+(1.0_dp-wt_rhs_tscheme_imp(1))*dt* &
                                   & (1.0_dp/Pr)*(r_radius(i)*real_d_temp_rad(i) + real_d2_temp_rad(i) &
                                   & -real(Nm,kind=dp)*real(Nm,kind=dp)*r_radius2(i) &
                                   & *real_temp_rad(i))
            
            rhs_exp_temp(1,Nm+1,i)=(ii*real(Nm,kind=dp)*r_radius(i)*uphi_temp_rad(i)+ &
                                   & real_d_ur_temp_rad(i)+r_radius(i)*ur_temp_rad(i))


            do i_order=1,n_order_tscheme_exp
               rhs_temp(i) = rhs_temp(i) - dt*(wt_rhs_tscheme_exp(i_order)*rhs_exp_temp(i_order,Nm+1,i)) 
            end do

         end do
      
      !----- If Backward-Difference Formula (BDF) is used as the implicit scheme --------------
      else if (time_scheme_imp=='BDF2' .or. time_scheme_imp=='BDF3' .or. time_scheme_imp=='BDF4') then 

         do i=1,Nr_max

            rhs_imp_temp(1,Nm+1,i)=real_temp_rad(i)
            !rhs_imp_temp(1,Nm+1,i)=tFR(Nm+1,i)
           
            do i_order=1,n_order_tscheme_imp
               rhs_temp(i) = rhs_temp(i) + wt_rhs_tscheme_imp(i_order)*rhs_imp_temp(i_order,Nm+1,i)
            end do

            rhs_exp_temp(1,Nm+1,i)=(ii*real(Nm,kind=dp)*r_radius(i)*uphi_temp_rad(i)+ &
                                   & real_d_ur_temp_rad(i)+r_radius(i)*ur_temp_rad(i))

            do i_order=1,n_order_tscheme_exp
               rhs_temp(i) = rhs_temp(i) - dt*(wt_rhs_tscheme_exp(i_order)*rhs_exp_temp(i_order,Nm+1,i)) 
            end do

         end do

      end if  
         
      if(Nm==0) then ! Apply BC for temperature
              rhs_temp(1)=1.0_dp
              rhs_temp(Nr_max)=0.0_dp
      else 
              rhs_temp(1)=0.0_dp
              rhs_temp(Nr_max)=0.0_dp
      end if

   end subroutine rhs_construct_temp

   subroutine rhs_construct_vort(Nm_max,Nr_max,dt,Ra,Pr,uphi_omg_FR,ur_omg_FR,n_step,omg_specp,Nm,rhs1,rhsf, &
                                 & n_restart,wt_rhs_tscheme_imp,wt_rhs_tscheme_exp,n_order_tscheme_imp, &
                                 & n_order_tscheme_exp,time_scheme_imp,tFRp) 

      integer, intent(in) :: Nm_max
      integer, intent(in) :: Nr_max 
      integer, intent(in) :: n_step, n_restart
      integer, intent(in) :: n_order_tscheme_exp
      integer, intent(in) :: n_order_tscheme_imp
      real(kind=dp), intent(in) :: wt_rhs_tscheme_imp(n_order_tscheme_imp)
      real(kind=dp), intent(in) :: wt_rhs_tscheme_exp(n_order_tscheme_exp)
      character(len=100), intent(in) :: time_scheme_imp
      real(kind=dp), intent(in) :: dt
      real(kind=dp), intent(in) :: Ra
      real(kind=dp), intent(in) :: Pr  
      complex(kind=dp), intent(in) :: omg_specp(Nm_max+1,Nr_max)
      integer, intent(in) :: Nm
      complex(kind=dp), intent(in) :: uphi_omg_FR(Nm_max+1,Nr_max),ur_omg_FR(Nm_max+1,Nr_max)
      integer :: i, i_order
      complex(kind=dp), intent(in) :: rhs1(Nr_max) 
      complex(kind=dp), intent(in) :: tFRp(Nm_max+1,Nr_max) 
      complex(kind=dp), intent(out) :: rhsf(2*Nr_max) 
      ! Local variables --------------------------
      complex(kind=dp) :: ur_omg_rad(Nr_max)
      complex(kind=dp) :: real_d_ur_omg_rad(Nr_max)
      complex(kind=dp) :: ur_omg_rad_spec(Nr_max)
      complex(kind=dp) :: uphi_omg_rad(Nr_max)
      complex(kind=dp) :: real_omg_rad(Nr_max)
      complex(kind=dp) :: real_d_omg_rad(Nr_max)
      complex(kind=dp) :: real_d2_omg_rad(Nr_max)
      complex(kind=dp) :: rhs2(Nr_max)
      
      if (n_step-n_restart>1 .or. n_restart/=0) then
      !if (n_step-n_restart>1) then
         do i_order=1,n_order_tscheme_exp-1
            rhs_exp_vort(i_order+1,Nm+1,:) = tmp_rhs_exp_vort(i_order,Nm+1,:)
         end do
         do i_order=1,n_order_tscheme_imp-1
            rhs_imp_vort(i_order+1,Nm+1,:) = tmp_rhs_imp_vort(i_order,Nm+1,:)
         end do
         !if (time_scheme_imp=='BDF2' .or. time_scheme_imp=='BDF3' .or. time_scheme_imp=='BDF4') then
         !   do i_order=1,n_order_tscheme_imp-1
         !      rhs_buo_term(i_order+1,Nm+1,:) = tmp_rhs_buo_term(i_order,Nm+1,:)
         !   end do
         !end if
      end if 
     
      !if (n_step-n_restart==1 .and. n_restart/=0) then
      !   do i_order=1,n_order_tscheme_imp-1
      !      rhs_buo_term(i_order+1,Nm+1,:) = tmp_rhs_buo_term(i_order,Nm+1,:)
      !   end do
      !end if 

      do i=1,Nr_max
         ur_omg_rad(i)=ur_omg_FR(Nm+1,i)
         uphi_omg_rad(i)=uphi_omg_FR(Nm+1,i)
         rhs2(i)=0.0_dp
         real_omg_rad(i)=0.0_dp
      end do

      call chebtransform(Nr_max,ur_omg_rad,ur_omg_rad_spec)

      call chebinvtran(Nr_max,omg_specp(Nm+1,:),real_omg_rad)

      call chebinvtranD1D2(Nr_max,omg_specp(Nm+1,:),real_d_omg_rad,real_d2_omg_rad)
      
      call chebinvtranD1(Nr_max,ur_omg_rad_spec,real_d_ur_omg_rad(:)) 
      
      !----- If Crank-Nicholson (CN) is used as the implicit scheme -------------------------
      if (time_scheme_imp=='CN') then 

         do i=1,Nr_max

            rhs_imp_vort(1,Nm+1,i) = real_omg_rad(i)
          
               rhs2(i) = rhs_imp_vort(1,Nm+1,i) + (1.0_dp-wt_rhs_tscheme_imp(1))*dt*(r_radius(i)* &
                                     & real_d_omg_rad(i)+real_d2_omg_rad(i)-real(Nm,kind=dp)* &
                                     & real(Nm,kind=dp)*r_radius2(i)*real_omg_rad(i)) - (Ra/Pr)*dt*0.5_dp* &
                                     &(ii*real(Nm,kind=dp)*(r_radius( &
                                        & i))*tFRp(Nm+1,i) + ii*real(Nm,kind=dp)*(r_radius(i))*rhs_imp_temp(1,Nm+1,i)) 
                                     ! Buoyancy is treated implicitly as Ra*Pr*dt*0.5*(dT/dphi_current + dT/dpi_old)
                                     ! (tFRp is current as it is coming after solving the temperature equation and
                                     !  rhs_imp_temp is the previous time step's temperature)  

            rhs_exp_vort(1,Nm+1,i) = (ii*real(Nm,kind=dp)*r_radius(i)*uphi_omg_rad(i)+real_d_ur_omg_rad(i) &
                                           & +r_radius(i)*ur_omg_rad(i)) 
            
                  
            do i_order=1,n_order_tscheme_exp
                rhs2(i) = rhs2(i) - dt*(wt_rhs_tscheme_exp(i_order)*rhs_exp_vort(i_order,Nm+1,i)) 
            end do

         end do


      !----- If Backward-Difference Formula (BDF) is used as the implicit scheme --------------
      else if (time_scheme_imp=='BDF2' .or. time_scheme_imp=='BDF3' .or. time_scheme_imp=='BDF4') then
            
         do i=1,Nr_max
            rhs_imp_vort(1,Nm+1,i)=omgFR(Nm+1,i) 

            do i_order=1,n_order_tscheme_imp
               rhs2(i) = rhs2(i) + wt_rhs_tscheme_imp(i_order)*(rhs_imp_vort(i_order,Nm+1,i))
            end do

            rhs_buo_term(1,Nm+1,i)=-(Ra/Pr)*(ii*real(Nm,kind=dp)*(r_radius( &
                                      & i))*tFRp(Nm+1,i))
            !rhs_buo_term(1,Nm+1,i)=(Ra/Pr)*(ii*real(Nm,kind=dp)*(r_radius( &
            !                          & i))*rhs_imp_temp(1,Nm+1,i))
            
            !do i_order=1,n_order_tscheme_imp
            !   rhs2(i) = rhs2(i) + wt_lhs_tscheme_imp*dt*wt_rhs_tscheme_imp(i_order)* &
            !           & rhs_buo_term(i_order,Nm+1,i)
            !end do
            ! Buoyancy is treated implicitly - Ra*Pr*dt*weight_of_scheme*(dT/dphi_current)


            rhs_exp_vort(1,Nm+1,i)=(ii*real(Nm,kind=dp)*r_radius(i)*uphi_omg_rad(i)+ &
                                   & real_d_ur_omg_rad(i) + r_radius(i)*ur_omg_rad(i)) !&
                                   !& + rhs_buo_term(1,Nm+1,i) 

            do i_order=1,n_order_tscheme_exp
               rhs2(i)=rhs2(i) - dt*(wt_rhs_tscheme_exp(i_order)*rhs_exp_vort(i_order,Nm+1,i)) 
            end do

            rhs2(i) = rhs2(i) + dt*wt_lhs_tscheme_imp*rhs_buo_term(1,Nm+1,i)

         end do
        
      end if

      rhs2(1)=0.0_dp ! Apply BC for stream function in the vorticity equation
      rhs2(Nr_max)=0.0_dp ! Apply BC for stream function in the vorticity equation

      do i=1,Nr_max

         rhsf(i)=rhs1(i)
         rhsf(i+Nr_max)=rhs2(i)

      end do

   end subroutine rhs_construct_vort

   subroutine rhs_construct_uphi_bar(Nm_max,Nr_max,dt,upFR_p,urFR_p,omgFR_p,n_step,upFCp,rhs_uphi,n_restart, &
                                 & wt_rhs_tscheme_imp,wt_rhs_tscheme_exp,n_order_tscheme_imp, &
                                 & n_order_tscheme_exp,time_scheme_imp)

      integer, intent(in) :: Nm_max
      integer, intent(in) :: Nr_max 
      integer, intent(in) :: n_step, n_restart
      integer, intent(in) :: n_order_tscheme_imp
      integer, intent(in) :: n_order_tscheme_exp
      real(kind=dp), intent(in) :: wt_rhs_tscheme_imp(n_order_tscheme_imp)
      real(kind=dp), intent(in) :: wt_rhs_tscheme_exp(n_order_tscheme_exp)
      character(len=100), intent(in) :: time_scheme_imp
      real(kind=dp), intent(in) :: dt
      complex(kind=dp), intent(in) :: upFCp(Nm_max+1,Nr_max)
      complex(kind=dp), intent(in) :: urFR_p(Nm_max+1,Nr_max),upFR_p(Nm_max+1,Nr_max)
      complex(kind=dp), intent(in) :: omgFR_p(Nm_max+1,Nr_max)
      real(kind=dp), intent(out) :: rhs_uphi(Nr_max)

      integer :: i, i_order, Nr, Nm, Np, Npmax
      real(kind=dp) :: ur_up_FR(Nr_max)
      real(kind=dp) :: ur_up(Nr_max)
      real(kind=dp) :: ur_omg(Nr_max)
      real(kind=dp) :: ur_omg_FR(Nr_max)
      real(kind=dp) :: ur_d1up(Nr_max)
      real(kind=dp) :: d1up(3*Nm_max,Nr_max)
      real(kind=dp) :: ur_d1up_FR(Nr_max)
      complex(kind=dp) :: u_diff(Nr_max)
      complex(kind=dp) :: d1upFR(Nm_max+1,Nr_max)
      complex(kind=dp) :: d2upFR(Nm_max+1,Nr_max)

      if (n_step-n_restart>1 .or. n_restart/=0) then
      !if (n_step-n_restart>1) then
         do i_order=1,n_order_tscheme_exp-1
            rhs_exp_uphi_bar(i_order+1,:) = tmp_rhs_exp_uphi_bar(i_order,:)
         end do
         do i_order=1,n_order_tscheme_imp-1
            rhs_imp_uphi_bar(i_order+1,:) = tmp_rhs_imp_uphi_bar(i_order,:)
         end do
      end if
      
      do i=1,Nr_max
         rhs_uphi(i)=0.0_dp
      end do
       
      do Nm=0,Nm_max 
         call chebinvtranD1(Nr_max,upFCp(Nm+1,:),d1upFR(Nm+1,:))
         call chebinvtranD2(Nr_max,upFCp(Nm+1,:),d2upFR(Nm+1,:))
      end do

      ur_up_FR(:)=0.0_dp
      ur_d1up_FR(:)=0.0_dp
      ur_up(:)=0.0_dp
      ur_omg(:)=0.0_dp
      ur_omg_FR(:)=0.0_dp
      ur_d1up(:)=0.0_dp

      ! ---- Product in Fourier space using transformed variables urFR and upFR --------------------
      do Nr=1,Nr_max
         do Nm=0,Nm_max
            if (Nm==0) then
               ur_up_FR(Nr) = ur_up_FR(Nr) + real(urFR_p(Nm+1,Nr)*(upFR_p(Nm+1,Nr)))
               ur_omg_FR(Nr) = ur_omg_FR(Nr) + real(urFR_p(Nm+1,Nr)*(omgFR_p(Nm+1,Nr)))
               ur_d1up_FR(Nr) = ur_d1up_FR(Nr) + real(urFR_p(Nm+1,Nr)*(d1upFR(Nm+1,Nr)))
            else
               ur_up_FR(Nr)=ur_up_FR(Nr) + real(urFR_p(Nm+1,Nr)*conjg(upFR_p(Nm+1,Nr)) &
                           & + conjg(urFR_p(Nm+1,Nr))*upFR_p(Nm+1,Nr))
               ur_omg_FR(Nr)=ur_omg_FR(Nr) + real(urFR_p(Nm+1,Nr)*conjg(omgFR_p(Nm+1,Nr)) &
                           & + conjg(urFR_p(Nm+1,Nr))*omgFR_p(Nm+1,Nr))
               ur_d1up_FR(Nr)=ur_d1up_FR(Nr) + real(urFR_p(Nm+1,Nr)*conjg(d1upFR(Nm+1,Nr)) &
                           & + conjg(urFR_p(Nm+1,Nr))*d1upFR(Nm+1,Nr))
            end if 
         end do
      end do
    
      ! ----- Product in physical space  -----------------
      Npmax=3*Nm_max
      do Nr=1,Nr_max
         call invfft(Nm_max,Npmax,d1upFR(:,Nr),d1up(:,Nr)) 
      end do

      do Nr=1,Nr_max
         do Np=1,Npmax
            ur_up(Nr)=ur_up(Nr) + ur(Np,Nr)*up(Np,Nr)
            ur_omg(Nr)=ur_omg(Nr) + ur(Np,Nr)*omg(Np,Nr)
            ur_d1up(Nr)=ur_d1up(Nr) + ur(Np,Nr)*d1up(Np,Nr)
         end do
      end do 

      ur_up=ur_up*(1.0_dp/Npmax)
      ur_d1up=ur_d1up*(1.0_dp/Npmax)
      ur_omg=ur_omg*(1.0_dp/Npmax)
      ! ---------------------------------------------------
      do Nr=1,Nr_max
         u_diff(Nr) = d2upFR(1,Nr) - r_radius2(Nr)*upFR_p(1,Nr) + r_radius(Nr)*d1upFR(1,Nr)
      end do

      !----- If Crank-Nicholson (CN) is used as the implicit scheme -------------------------
      if (time_scheme_imp=='CN') then 

         do i=1,Nr_max

            rhs_imp_uphi_bar(1,i)=(upFR_p(1,i))

            rhs_uphi(i) = real(rhs_imp_uphi_bar(1,i))+(1.0_dp-wt_rhs_tscheme_imp(1))*dt* &
                                   & real(u_diff(i)) ! Diffusive part
            
            rhs_exp_uphi_bar(1,i)=(ur_d1up_FR(i)+r_radius(i)*ur_up_FR(i)) ! Advective part in Fourier-Real space
            !rhs_exp_uphi_bar(1,i)=(ur_d1up(i)+r_radius(i)*ur_up(i)) ! Advective part in Physical space (uncomment for checking)
            !rhs_exp_uphi_bar(1,i)=(ur_omg_FR(i)) ! Advective part in Fourier-Real space
            !rhs_exp_uphi_bar(1,i)=(ur_omg(i)) ! Advective part in Physical space (uncomment for checking)
            
            do i_order=1,n_order_tscheme_exp
               rhs_uphi(i) = rhs_uphi(i) - dt*(wt_rhs_tscheme_exp(i_order)*real(rhs_exp_uphi_bar(i_order,i))) 
            end do

         end do
      !----- If Backward-Difference Formula (BDF) is used as the implicit scheme --------------
      else if (time_scheme_imp=='BDF2' .or. time_scheme_imp=='BDF3' .or. time_scheme_imp=='BDF4' ) then 
         
         do i=1,Nr_max

            rhs_imp_uphi_bar(1,i)=(upFR_p(1,i))
           
            do i_order=1,n_order_tscheme_imp
               rhs_uphi(i) = rhs_uphi(i) + wt_rhs_tscheme_imp(i_order)*real(rhs_imp_uphi_bar(i_order,i))
            end do

            rhs_exp_uphi_bar(1,i)=((ur_d1up_FR(i)+r_radius(i)*ur_up_FR(i)))

            do i_order=1,n_order_tscheme_exp
               rhs_uphi(i) = rhs_uphi(i) - dt*(wt_rhs_tscheme_exp(i_order)*real(rhs_exp_uphi_bar(i_order,i))) 
            end do

         end do

      end if  
         
      ! Apply BC for uphi_bar
      rhs_uphi(1)=0.0_dp
      rhs_uphi(Nr_max)=0.0_dp
         
   end subroutine rhs_construct_uphi_bar 

   subroutine rhs_temp_from_restart(Nm_max,Nr_max,uphi_temp_FR,ur_temp_FR,temp_specp,time_scheme_imp, &
                                       & n_order_tscheme_imp, n_order_tscheme_exp)
             
      integer, intent(in) :: Nm_max
      integer, intent(in) :: Nr_max 
      complex(kind=dp), intent(in) :: temp_specp(Nm_max+1,Nr_max)
      integer :: Nm
      complex(kind=dp), intent(in) :: uphi_temp_FR(Nm_max+1,Nr_max),ur_temp_FR(Nm_max+1,Nr_max)
      integer :: i
      character(len=100), intent(in) :: time_scheme_imp
      integer, intent(in) :: n_order_tscheme_imp, n_order_tscheme_exp
      ! Local variables --------------------------
      complex(kind=dp) :: ur_temp_rad(Nr_max)
      complex(kind=dp) :: real_d_ur_temp_rad(Nr_max)
      complex(kind=dp) :: ur_temp_rad_spec(Nr_max)
      complex(kind=dp) :: uphi_temp_rad(Nr_max)
      complex(kind=dp) :: real_temp_rad(Nr_max)
      complex(kind=dp) :: real_d_temp_rad(Nr_max)
      complex(kind=dp) :: real_d2_temp_rad(Nr_max)

      if (time_scheme_imp=='BDF2' .or. time_scheme_imp=='BDF3' .or. time_scheme_imp=='BDF4') then 
         do Nm=0,Nm_max
            do i=1,Nr_max
               uphi_temp_rad(i)=uphi_temp_FR(Nm+1,i)
               ur_temp_rad(i)=ur_temp_FR(Nm+1,i)
            end do

            call chebtransform(Nr_max,ur_temp_rad,ur_temp_rad_spec)

            call chebinvtranD1(Nr_max,ur_temp_rad_spec,real_d_ur_temp_rad(:))

            !----- If Backward-Difference Formula (BDF) is used as the implicit scheme --------------
            do i=1,Nr_max
               
               rhs_imp_temp(n_order_tscheme_imp,Nm+1,i)=temp_specp(Nm+1,i) 
             
               rhs_exp_temp(n_order_tscheme_exp,Nm+1,i)=(ii*real(Nm,kind=dp)*r_radius(i)*uphi_temp_rad(i)+ &
                                      & real_d_ur_temp_rad(i)+r_radius(i)*ur_temp_rad(i))

            end do
         end do
      end if  
         
     
         
   end subroutine rhs_temp_from_restart

   subroutine rhs_vort_from_restart(Nm_max,Nr_max,Ra,Pr,tFRp,uphi_omg_FR,ur_omg_FR,omg_specp,time_scheme_imp, &
                                       & n_order_tscheme_imp, n_order_tscheme_exp) 

      integer, intent(in) :: Nm_max
      integer, intent(in) :: Nr_max 
      real(kind=dp), intent(in) :: Ra
      real(kind=dp), intent(in) :: Pr
      complex(kind=dp), intent(in) :: omg_specp(Nm_max+1,Nr_max)
      complex(kind=dp), intent(in) :: tFRp(Nm_max+1,Nr_max)
      integer :: Nm
      complex(kind=dp), intent(in) :: uphi_omg_FR(Nm_max+1,Nr_max),ur_omg_FR(Nm_max+1,Nr_max)
      integer :: i
      character(len=100), intent(in) :: time_scheme_imp
      integer, intent(in) :: n_order_tscheme_imp, n_order_tscheme_exp

      ! Local variables --------------------------
      complex(kind=dp) :: ur_omg_rad(Nr_max)
      complex(kind=dp) :: real_d_ur_omg_rad(Nr_max)
      complex(kind=dp) :: ur_omg_rad_spec(Nr_max)
      complex(kind=dp) :: uphi_omg_rad(Nr_max)
      complex(kind=dp) :: real_omg_rad(Nr_max)
      complex(kind=dp) :: real_d_omg_rad(Nr_max)
      complex(kind=dp) :: real_d2_omg_rad(Nr_max)

      !----- If Backward-Difference Formula (BDF) is used as the implicit scheme --------------
      if (time_scheme_imp=='BDF2' .or. time_scheme_imp=='BDF3' .or. time_scheme_imp=='BDF4') then
         do Nm=0,Nm_max      
            do i=1,Nr_max
               ur_omg_rad(i)=ur_omg_FR(Nm+1,i)
               uphi_omg_rad(i)=uphi_omg_FR(Nm+1,i)
            end do

            call chebtransform(Nr_max,ur_omg_rad,ur_omg_rad_spec)

            call chebinvtranD1(Nr_max,ur_omg_rad_spec,real_d_ur_omg_rad(:)) 
            
            do i=1,Nr_max
               rhs_imp_vort(n_order_tscheme_imp,Nm+1,i)=omg_specp(Nm+1,i) 
                        
               rhs_exp_vort(n_order_tscheme_exp,Nm+1,i)=(ii*real(Nm,kind=dp)*r_radius(i)*uphi_omg_rad(i)+ &
                                      & real_d_ur_omg_rad(i) + r_radius(i)*ur_omg_rad(i)) &
                                      & + (Ra/Pr)*(ii*real(Nm,kind=dp)*(r_radius( &
                                      & i))*tFRp(Nm+1,i))

               !rhs_buo_term(n_order_tscheme_imp,Nm+1,i)=-(Ra/Pr)*(ii*real(Nm,kind=dp)*(r_radius( &
               !                       & i))*tFRp(Nm+1,i))
                     
            end do
         end do 
      end if

   end subroutine rhs_vort_from_restart

   subroutine rhs_uphibar_from_restart(Nm_max,Nr_max,upFR_p,urFR_p,omgFR_p,time_scheme_imp, &
                                       & n_order_tscheme_imp, n_order_tscheme_exp)

      integer, intent(in) :: Nm_max
      integer, intent(in) :: Nr_max 
      character(len=100), intent(in) :: time_scheme_imp
      complex(kind=dp), intent(in) :: urFR_p(Nm_max+1,Nr_max),upFR_p(Nm_max+1,Nr_max)
      complex(kind=dp), intent(in) :: omgFR_p(Nm_max+1,Nr_max)
      integer, intent(in) :: n_order_tscheme_imp, n_order_tscheme_exp

      integer :: i, Nr, Nm, Np, Npmax
      real(kind=dp) :: ur_up_FR(Nr_max)
      real(kind=dp) :: ur_up(Nr_max)
      real(kind=dp) :: ur_omg(Nr_max)
      real(kind=dp) :: ur_omg_FR(Nr_max)
      real(kind=dp) :: ur_d1up(Nr_max)
      real(kind=dp) :: d1up(3*Nm_max,Nr_max)
      real(kind=dp) :: ur_d1up_FR(Nr_max)
      complex(kind=dp) :: u_diff(Nr_max)
      complex(kind=dp) :: d1upFR(Nm_max+1,Nr_max)
      complex(kind=dp) :: d2upFR(Nm_max+1,Nr_max)
      complex(kind=dp) :: upFCp(Nm_max+1,Nr_max)
     
      !----- If Backward-Difference Formula (BDF) is used as the implicit scheme --------------
      if (time_scheme_imp=='BDF2' .or. time_scheme_imp=='BDF3' .or. time_scheme_imp=='BDF4' ) then 

         do Nm=0,Nm_max 
            call chebtransform(Nr_max,upFR_p(Nm+1,:),upFCp(Nm+1,:))
            call chebinvtranD1(Nr_max,upFCp(Nm+1,:),d1upFR(Nm+1,:))
            call chebinvtranD2(Nr_max,upFCp(Nm+1,:),d2upFR(Nm+1,:))
         end do

         ur_up_FR(:)=0.0_dp
         ur_d1up_FR(:)=0.0_dp
         ur_up(:)=0.0_dp
         ur_omg(:)=0.0_dp
         ur_omg_FR(:)=0.0_dp
         ur_d1up(:)=0.0_dp

         ! ---- Product in Fourier space using transformed variables urFR and upFR --------------------
         do Nr=1,Nr_max
            do Nm=0,Nm_max
               if (Nm==0) then
                  ur_up_FR(Nr) = ur_up_FR(Nr) + real(urFR_p(Nm+1,Nr)*(upFR_p(Nm+1,Nr)))
                  ur_omg_FR(Nr) = ur_omg_FR(Nr) + real(urFR_p(Nm+1,Nr)*(omgFR_p(Nm+1,Nr)))
                  ur_d1up_FR(Nr) = ur_d1up_FR(Nr) + real(urFR_p(Nm+1,Nr)*(d1upFR(Nm+1,Nr)))
               else
                  ur_up_FR(Nr)=ur_up_FR(Nr) + real(urFR_p(Nm+1,Nr)*conjg(upFR_p(Nm+1,Nr)) &
                              & + conjg(urFR_p(Nm+1,Nr))*upFR_p(Nm+1,Nr))
                  ur_omg_FR(Nr)=ur_omg_FR(Nr) + real(urFR_p(Nm+1,Nr)*conjg(omgFR_p(Nm+1,Nr)) &
                              & + conjg(urFR_p(Nm+1,Nr))*omgFR_p(Nm+1,Nr))
                  ur_d1up_FR(Nr)=ur_d1up_FR(Nr) + real(urFR_p(Nm+1,Nr)*conjg(d1upFR(Nm+1,Nr)) &
                              & + conjg(urFR_p(Nm+1,Nr))*d1upFR(Nm+1,Nr))
               end if 
            end do
         end do
    
         ! ----- Product in physical space  -----------------
         Npmax=3*Nm_max
         do Nr=1,Nr_max
            call invfft(Nm_max,Npmax,d1upFR(:,Nr),d1up(:,Nr)) 
         end do

         do Nr=1,Nr_max
            do Np=1,Npmax
               ur_up(Nr)=ur_up(Nr) + ur(Np,Nr)*up(Np,Nr)
               ur_omg(Nr)=ur_omg(Nr) + ur(Np,Nr)*omg(Np,Nr)
               ur_d1up(Nr)=ur_d1up(Nr) + ur(Np,Nr)*d1up(Np,Nr)
            end do
         end do 

         ur_up=ur_up*(1.0_dp/Npmax)
         ur_d1up=ur_d1up*(1.0_dp/Npmax)
         ur_omg=ur_omg*(1.0_dp/Npmax)
         ! ---------------------------------------------------

         do Nr=1,Nr_max
            u_diff(Nr) = d2upFR(1,Nr) - r_radius2(Nr)*upFR_p(1,Nr) + r_radius(Nr)*d1upFR(1,Nr)
         end do

            
         do i=1,Nr_max

            rhs_imp_uphi_bar(n_order_tscheme_imp,i)=(upFR_p(1,i))
           
            rhs_exp_uphi_bar(n_order_tscheme_exp,i)=((ur_d1up_FR(i)+r_radius(i)*ur_up_FR(i)))

         end do

      end if  
         
   end subroutine rhs_uphibar_from_restart 


end module rhs_create
