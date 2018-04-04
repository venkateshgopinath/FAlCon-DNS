module rhs_create_rk

   use double
   use constants, only: ii
   use fourier, only: invfft
   use chebyshev, only: chebtransform, chebinvtran, chebinvtranD1, chebinvtranD1D2, chebinvtranD2 
   use init, only: r_radius, r_radius2, tFR, t2FR, omgFR, ur, up, omg
   use timeschemes, only:  wt_lhs_tscheme_imp, wt_rhs_tscheme_imp, rhs_exp_temp, rhs_exp_vort, &
                          & rhs_imp_temp, rhs_imp_vort, dt_array

   implicit none

   private

   complex(kind=dp), allocatable :: uphi_temp_rad(:), ur_temp_rad_spec(:), ur_temp_rad(:)
   complex(kind=dp), allocatable :: uphi_omg_rad(:), ur_omg_rad_spec(:), ur_omg_rad(:) 
   complex(kind=dp), allocatable :: real_omg_spec(:), temp_real_rad(:), temp_spec_rad(:) 
   complex(kind=dp), allocatable :: real_d_omg_rad(:), real_d2_omg_rad(:)
   complex(kind=dp), allocatable :: omg_real_rad(:), omg_spec_rad(:)  
   complex(kind=dp), allocatable :: real_temp_rad(:), real_d_temp_rad(:), real_d2_temp_rad(:)
   complex(kind=dp), allocatable :: real_d_ur_temp_rad(:,:), real_d_ur_omg_rad(:,:)
   complex(kind=dp), allocatable, public :: rhs_vort_buo_term(:,:)
   complex(kind=dp), allocatable, public :: rhs2(:)

   public :: allocate_rhs_rk, deallocate_rhs_rk, rhs_construct_temp, rhs_construct_vort, rhs_construct_uphi_bar

contains

   subroutine allocate_rhs_rk(Nm_max,Nr_max)

      integer, intent(in) :: Nm_max
      integer, intent(in) :: Nr_max 
           
      allocate( uphi_temp_rad(Nr_max), ur_temp_rad_spec(Nr_max), ur_temp_rad(Nr_max), &
                & uphi_omg_rad(Nr_max), ur_omg_rad_spec(Nr_max) )
      allocate( real_omg_spec(Nr_max), temp_real_rad(Nr_max), temp_spec_rad(Nr_max), &
                & omg_real_rad(Nr_max), omg_spec_rad(Nr_max), &
                & ur_omg_rad(Nr_max), real_d_omg_rad(Nr_max), real_d2_omg_rad(Nr_max) )
      allocate( real_temp_rad(Nr_max), real_d_temp_rad(Nr_max), real_d2_temp_rad(Nr_max) )
      allocate( real_d_ur_temp_rad(Nm_max+1,Nr_max), real_d_ur_omg_rad(Nm_max+1,Nr_max) )
      allocate( rhs_vort_buo_term(Nm_max+1,Nr_max) )
      allocate( rhs2(Nr_max) )
     

   end subroutine allocate_rhs_rk

   subroutine deallocate_rhs_rk

      deallocate( rhs2 )
      deallocate( rhs_vort_buo_term )
      deallocate( real_d_ur_temp_rad, real_d_ur_omg_rad )
      deallocate( real_temp_rad, real_d_temp_rad, real_d2_temp_rad )
      deallocate( real_omg_spec, temp_real_rad, temp_spec_rad,omg_real_rad, omg_spec_rad, &
                  & ur_omg_rad, real_d_omg_rad, real_d2_omg_rad )
      deallocate( uphi_temp_rad, ur_temp_rad_spec, ur_temp_rad, uphi_omg_rad, ur_omg_rad_spec )

   end subroutine deallocate_rhs_rk

   subroutine rhs_construct_temp(Nm_max,Nr_max,uphi_temp_FR,ur_temp_FR,temp_spec,Nm,rhs, &
                                 & time_scheme_exp,Pr)
             
      integer, intent(in) :: Nm_max
      integer, intent(in) :: Nr_max 
      real(kind=dp), intent(in) :: Pr
      character(len=100), intent(in) :: time_scheme_exp
      complex(kind=dp), intent(in) :: temp_spec(Nm_max+1,Nr_max)
      integer, intent(in) :: Nm
      complex(kind=dp), intent(in) :: uphi_temp_FR(Nm_max+1,Nr_max),ur_temp_FR(Nm_max+1,Nr_max)
      integer :: i
      complex(kind=dp), intent(out) :: rhs(Nr_max)
      
      if (time_scheme_exp=='RK4' .or. time_scheme_exp=='RK2') then
         do i=1,Nr_max
            uphi_temp_rad(i)=uphi_temp_FR(Nm+1,i)
            ur_temp_rad(i)=ur_temp_FR(Nm+1,i)
            temp_real_rad(i)=temp_spec(Nm+1,i)
            rhs(i)=0.0_dp
         end do

         call chebtransform(Nr_max,ur_temp_rad,ur_temp_rad_spec)
         
         call chebtransform(Nr_max,temp_real_rad,temp_spec_rad)
                        
         call chebinvtranD1D2(Nr_max,temp_spec_rad,real_d_temp_rad,real_d2_temp_rad)
         
         call chebinvtranD1(Nr_max,ur_temp_rad_spec,real_d_ur_temp_rad(Nm+1,:))

         do i=1,Nr_max

            rhs(i)= (1.0_dp/Pr)*(r_radius(i)*real_d_temp_rad(i) + real_d2_temp_rad(i) & 
                    & -real(Nm,kind=dp)*real(Nm,kind=dp)*r_radius2(i)*temp_real_rad(i)) &
                    & - (ii*real(Nm,kind=dp)*r_radius(i)*uphi_temp_rad(i)+ &
                    & real_d_ur_temp_rad(Nm+1,i)+r_radius(i)*ur_temp_rad(i))
                    ! Sum of linear + non-linear parts

         end do

      end if  
         
   end subroutine rhs_construct_temp

   subroutine rhs_construct_vort(Nm_max,Nr_max,Ra,Pr,uphi_omg_FR,ur_omg_FR,omg2,Nm,rhs1,rhsf, &
                                    & time_scheme_exp,rhs2) 

      integer, intent(in) :: Nm_max
      integer, intent(in) :: Nr_max 
      character(len=100), intent(in) :: time_scheme_exp
      real(kind=dp), intent(in) :: Ra
      real(kind=dp), intent(in) :: Pr  
      complex(kind=dp), intent(in) :: omg2(Nm_max+1,Nr_max)
      integer, intent(in) :: Nm
      complex(kind=dp), intent(in) :: uphi_omg_FR(Nm_max+1,Nr_max),ur_omg_FR(Nm_max+1,Nr_max)
      integer :: i
      complex(kind=dp), intent(in) :: rhs1(Nr_max) 
      complex(kind=dp), intent(out) :: rhsf(2*Nr_max) 
      complex(kind=dp), intent(out) :: rhs2(Nr_max) 

      if (time_scheme_exp=='RK4' .or. time_scheme_exp=='RK2') then
         do i=1,Nr_max
            ur_omg_rad(i)=ur_omg_FR(Nm+1,i)
            uphi_omg_rad(i)=uphi_omg_FR(Nm+1,i)
            rhs2(i)=0.0_dp
            omg_spec_rad(i)=0.0_dp
            omg_real_rad(i)=omg2(Nm+1,i)
         end do

         call chebtransform(Nr_max,ur_omg_rad,ur_omg_rad_spec)
         
         call chebtransform(Nr_max,omg_real_rad,omg_spec_rad)
         
         call chebinvtranD1D2(Nr_max,omg_spec_rad,real_d_omg_rad,real_d2_omg_rad)
         
         call chebinvtranD1(Nr_max,ur_omg_rad_spec,real_d_ur_omg_rad(Nm+1,:))
     
         do i=1,Nr_max
            rhs2(i)= (r_radius(i)*real_d_omg_rad(i)+real_d2_omg_rad(i)-real(Nm,kind=dp)* &
                    & real(Nm,kind=dp)*r_radius2(i)*omg_real_rad(i)) &
                    & - (ii*real(Nm,kind=dp)*r_radius(i)*uphi_omg_rad(i)+real_d_ur_omg_rad(Nm+1,i) &
                    & + r_radius(i)*ur_omg_rad(i)) - (Ra/Pr)*(ii*real(Nm,kind=dp)*(r_radius(i))*tFR(Nm+1,i))
                    ! Sum of linear + non-linear parts + buoyancy term

         end do
         
      end if

      do i=1,Nr_max

         rhsf(i)=rhs1(i)
         rhsf(i+Nr_max)=rhs2(i)

      end do

   end subroutine rhs_construct_vort 

   subroutine rhs_construct_uphi_bar(Nm_max,Nr_max,upFR_p,urFR_p,upFC,omgFR_p,rhs_uphi, &
                                    & time_scheme_exp)
      integer, intent(in) :: Nm_max
      integer, intent(in) :: Nr_max 
      character(len=100), intent(in) :: time_scheme_exp
      complex(kind=dp), intent(in) :: upFC(Nm_max+1,Nr_max)
      complex(kind=dp), intent(in) :: urFR_p(Nm_max+1,Nr_max),upFR_p(Nm_max+1,Nr_max)
      complex(kind=dp), intent(in) :: omgFR_p(Nm_max+1,Nr_max)
      real(kind=dp), intent(out) :: rhs_uphi(Nr_max)

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
     
      do i=1,Nr_max
         rhs_uphi(i)=0.0_dp
      end do
       
      do Nm=0,Nm_max 
         call chebinvtranD1(Nr_max,upFC(Nm+1,:),d1upFR(Nm+1,:))
         call chebinvtranD2(Nr_max,upFC(Nm+1,:),d2upFR(Nm+1,:))
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
                                
      if (time_scheme_exp=='RK4' .or. time_scheme_exp=='RK2') then

         do i=1,Nr_max

            rhs_uphi(i) = real(u_diff(i)) - (ur_d1up_FR(i)+r_radius(i)*ur_up_FR(i)) 
                         ! Sum of linear + non-linear parts 
            
            !(ur_d1up(i)+r_radius(i)*ur_up(i)) ! Advective part in Physical space (uncomment for checking)
            !(ur_omg_FR(i)) ! Advective part in Fourier space (Uncomment for checking)
            !(ur_omg(i)) ! Advective part in Physical space (uncomment for checking)
            
         end do

      end if  
         
   end subroutine rhs_construct_uphi_bar 

end module rhs_create_rk
