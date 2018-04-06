module rhs_create_imexrk

   use double
   use constants, only: ii
   use chebyshev, only: chebtransform, chebinvtran, chebinvtranD1, chebinvtranD1D2, chebinvtranD2 
   use init, only: r_radius, r_radius2, ur, up
   use timeschemes, only:  wt_lhs_tscheme_imp, wt_rhs_tscheme_imp, rhs_exp_temp, rhs_exp_vort, &
                          & rhs_imp_temp, rhs_imp_vort, dt_array

   implicit none

   private

   public :: rhs_construct_temp_a, rhs_construct_vort_a, &
             & rhs_construct_uphibar_a, rhs_construct_uphibar_d, & 
             & rhs_construct_temp_d, rhs_construct_vort_d, rhs_construct_buo    

contains

   subroutine rhs_construct_temp_a(Nm_max,Nr_max,uphi_temp_FR,ur_temp_FR,temprhs,Nm,rhs, &
                                 & time_scheme_type)
             
      integer, intent(in) :: Nm_max
      integer, intent(in) :: Nr_max 
      character(len=100), intent(in) :: time_scheme_type
      complex(kind=dp), intent(in) :: temprhs(Nr_max)
      integer, intent(in) :: Nm
      complex(kind=dp), intent(in) :: uphi_temp_FR(Nm_max+1,Nr_max),ur_temp_FR(Nm_max+1,Nr_max)
      integer :: i
      complex(kind=dp), intent(out) :: rhs(Nr_max)
      ! Local variables --------------------------
      complex(kind=dp) :: ur_temp_rad(Nr_max)
      complex(kind=dp) :: real_d_ur_temp_rad(Nr_max)
      complex(kind=dp) :: ur_temp_rad_spec(Nr_max)
      complex(kind=dp) :: uphi_temp_rad(Nr_max)
      complex(kind=dp) :: real_temp_rad(Nr_max)
      complex(kind=dp) :: temp_spec_rad(Nr_max)
      complex(kind=dp) :: real_d_temp_rad(Nr_max)
      complex(kind=dp) :: real_d2_temp_rad(Nr_max) 
      
      if (time_scheme_type=='IMEXRK') then
         do i=1,Nr_max
            uphi_temp_rad(i)=uphi_temp_FR(Nm+1,i)
            ur_temp_rad(i)=ur_temp_FR(Nm+1,i)
            rhs(i)=0.0_dp
         end do

         call chebtransform(Nr_max,ur_temp_rad,ur_temp_rad_spec)
         
         call chebinvtranD1(Nr_max,ur_temp_rad_spec,real_d_ur_temp_rad(:))

         do i=1,Nr_max

            rhs(i)= - (ii*real(Nm,kind=dp)*r_radius(i)*uphi_temp_rad(i)+ &
                    & real_d_ur_temp_rad(i)+r_radius(i)*ur_temp_rad(i))
                    ! non-linear part

         end do

      end if  

   end subroutine rhs_construct_temp_a

   subroutine rhs_construct_temp_d(Nm_max,Nr_max,temprhs,Nm,rhs, &
                                 & time_scheme_type,Pr)
             
      integer, intent(in) :: Nm_max
      integer, intent(in) :: Nr_max 
      real(kind=dp), intent(in) :: Pr
      character(len=100), intent(in) :: time_scheme_type
      complex(kind=dp), intent(in) :: temprhs(Nr_max)
      integer, intent(in) :: Nm
      integer :: i
      complex(kind=dp), intent(out) :: rhs(Nr_max)
      ! Local variables --------------------------
      complex(kind=dp) :: real_temp_rad(Nr_max)
      complex(kind=dp) :: temp_spec_rad(Nr_max)
      complex(kind=dp) :: real_d_temp_rad(Nr_max)
      complex(kind=dp) :: real_d2_temp_rad(Nr_max) 

      if (time_scheme_type=='IMEXRK') then
         do i=1,Nr_max
            rhs(i)=0.0_dp
         end do

         call chebtransform(Nr_max,temprhs,temp_spec_rad)
                        
         call chebinvtranD1D2(Nr_max,temp_spec_rad,real_d_temp_rad,real_d2_temp_rad)
         
         do i=1,Nr_max

            rhs(i)= (1.0_dp/Pr)*(r_radius(i)*real_d_temp_rad(i) + real_d2_temp_rad(i) & 
                    & -real(Nm,kind=dp)*real(Nm,kind=dp)*r_radius2(i)*temprhs(i)) 
                    ! linear part

         end do
            
      end if 

   end subroutine rhs_construct_temp_d  

   subroutine rhs_construct_uphibar_a(Nm_max,Nr_max,urFR_p,upFR_p,upFC,rhs_uphi, &
                                 & time_scheme_type)
             
      integer, intent(in) :: Nm_max
      integer, intent(in) :: Nr_max 
      character(len=100), intent(in) :: time_scheme_type
      complex(kind=dp), intent(in) :: urFR_p(Nm_max+1,Nr_max)
      complex(kind=dp), intent(in) :: upFR_p(Nm_max+1,Nr_max)
      complex(kind=dp), intent(in) :: upFC(Nm_max+1,Nr_max)
      real(kind=dp), intent(out) :: rhs_uphi(Nr_max)
      
      integer :: i, Nr, Nm
      real(kind=dp) :: ur_up_FR(Nr_max)
      real(kind=dp) :: ur_d1up_FR(Nr_max)
      complex(kind=dp) :: d1upFR(Nm_max+1,Nr_max)

      if (time_scheme_type=='IMEXRK') then

         do i=1,Nr_max
            rhs_uphi(i)=0.0_dp
         end do

         do Nm=0,Nm_max
            call chebinvtranD1(Nr_max,upFC(Nm+1,:),d1upFR(Nm+1,:))
         end do

         ur_up_FR(:)=0.0_dp
         ur_d1up_FR(:)=0.0_dp

         ! ---- Product in Fourier space using transformed variables urFR and upFR --------------------
         do Nr=1,Nr_max
            do Nm=0,Nm_max
               if (Nm==0) then
                  ur_up_FR(Nr) = ur_up_FR(Nr) + real(urFR_p(Nm+1,Nr)*(upFR_p(Nm+1,Nr)))
                  ur_d1up_FR(Nr) = ur_d1up_FR(Nr) + real(urFR_p(Nm+1,Nr)*(d1upFR(Nm+1,Nr)))
               else
                  ur_up_FR(Nr)=ur_up_FR(Nr) + real(urFR_p(Nm+1,Nr)*conjg(upFR_p(Nm+1,Nr)) &
                              & + conjg(urFR_p(Nm+1,Nr))*upFR_p(Nm+1,Nr))
                  ur_d1up_FR(Nr)=ur_d1up_FR(Nr) + real(urFR_p(Nm+1,Nr)*conjg(d1upFR(Nm+1,Nr)) &
                              & + conjg(urFR_p(Nm+1,Nr))*d1upFR(Nm+1,Nr))
               end if
            end do
         end do

         do i=1,Nr_max

            rhs_uphi(i)=-(ur_d1up_FR(i)+r_radius(i)*ur_up_FR(i)) ! Advective part in Fourier-Real space

         end do

      end if
      

   end subroutine rhs_construct_uphibar_a

   subroutine rhs_construct_uphibar_d(Nm_max,Nr_max,upFR_p,upFC,rhs_uphi, &
                                 & time_scheme_type)
             
      integer, intent(in) :: Nm_max
      integer, intent(in) :: Nr_max 
      character(len=100), intent(in) :: time_scheme_type
      complex(kind=dp), intent(in) :: upFC(Nm_max+1,Nr_max)
      complex(kind=dp), intent(in) :: upFR_p(Nm_max+1,Nr_max)
      real(kind=dp), intent(out) :: rhs_uphi(Nr_max)
      
      integer :: i, Nm
      real(kind=dp) :: ur_up_FR(Nr_max)
      real(kind=dp) :: ur_omg_FR(Nr_max)
      real(kind=dp) :: ur_d1up_FR(Nr_max)
      complex(kind=dp) :: d1upFR(Nm_max+1,Nr_max)
      complex(kind=dp) :: d2upFR(Nm_max+1,Nr_max)

      if (time_scheme_type=='IMEXRK') then
     
         do i=1,Nr_max
            rhs_uphi(i)=0.0_dp
         end do
          
         do Nm=0,Nm_max 
            call chebinvtranD1(Nr_max,upFC(Nm+1,:),d1upFR(Nm+1,:))
            call chebinvtranD2(Nr_max,upFC(Nm+1,:),d2upFR(Nm+1,:))
         end do

         ur_up_FR(:)=0.0_dp
         ur_d1up_FR(:)=0.0_dp
         ur_omg_FR(:)=0.0_dp

         do i=1,Nr_max

            rhs_uphi(i)= real(d2upFR(1,i) - r_radius2(i)*upFR_p(1,i) + r_radius(i)*d1upFR(1,i),kind=dp) ! Diffusive part in Fourier-Real space

         end do

      end if  

   end subroutine rhs_construct_uphibar_d


   subroutine rhs_construct_vort_a(Nm_max,Nr_max,uphi_omg_FR,ur_omg_FR,omgrhs,Nm,rhs1,rhsf, &
                                    & time_scheme_type,rhs2) 

      integer, intent(in) :: Nm_max
      integer, intent(in) :: Nr_max 
      character(len=100), intent(in) :: time_scheme_type
      complex(kind=dp), intent(in) :: omgrhs(Nr_max)
      integer, intent(in) :: Nm
      complex(kind=dp), intent(in) :: uphi_omg_FR(Nm_max+1,Nr_max),ur_omg_FR(Nm_max+1,Nr_max)
      integer :: i
      complex(kind=dp), intent(in) :: rhs1(Nr_max) 
      complex(kind=dp), intent(out) :: rhsf(2*Nr_max) 
      complex(kind=dp), intent(out) :: rhs2(Nr_max) 
      ! Local variables --------------------------
      complex(kind=dp) :: ur_omg_rad(Nr_max)
      complex(kind=dp) :: real_d_ur_omg_rad(Nr_max)
      complex(kind=dp) :: ur_omg_rad_spec(Nr_max)
      complex(kind=dp) :: uphi_omg_rad(Nr_max)
      complex(kind=dp) :: real_omg_rad(Nr_max)

      if (time_scheme_type=='IMEXRK') then
         do i=1,Nr_max
            ur_omg_rad(i)=ur_omg_FR(Nm+1,i)
            uphi_omg_rad(i)=uphi_omg_FR(Nm+1,i)
            rhs2(i)=0.0_dp
         end do

         call chebtransform(Nr_max,ur_omg_rad,ur_omg_rad_spec)
         
         call chebinvtranD1(Nr_max,ur_omg_rad_spec,real_d_ur_omg_rad(:))
     
         do i=1,Nr_max
            rhs2(i)= - (ii*real(Nm,kind=dp)*r_radius(i)*uphi_omg_rad(i)+real_d_ur_omg_rad(i) &
                    & + r_radius(i)*ur_omg_rad(i)) 
                    ! non-linear part

         end do
         
      end if

      rhs2(1)=0.0_dp ! Apply BC for stream function in the vorticity equation
      rhs2(Nr_max)=0.0_dp ! Apply BC for stream function in the vorticity equation

      do i=1,Nr_max

         rhsf(i)=rhs1(i)
         rhsf(i+Nr_max)=rhs2(i)

      end do

   end subroutine rhs_construct_vort_a

   subroutine rhs_construct_vort_d(Nm_max,Nr_max,omgrhs,Nm,rhs1,rhsf, &
                                    & time_scheme_type,rhs2) 

      integer, intent(in) :: Nm_max
      integer, intent(in) :: Nr_max 
      character(len=100), intent(in) :: time_scheme_type
      complex(kind=dp), intent(in) :: omgrhs(Nr_max)
      integer, intent(in) :: Nm
      integer :: i
      complex(kind=dp), intent(in) :: rhs1(Nr_max) 
      complex(kind=dp), intent(out) :: rhsf(2*Nr_max) 
      complex(kind=dp), intent(out) :: rhs2(Nr_max) 
      ! Local variables --------------------------
      complex(kind=dp) :: real_d_omg_rad(Nr_max)
      complex(kind=dp) :: real_d2_omg_rad(Nr_max)
      complex(kind=dp) :: omg_spec_rad(Nr_max)

      if (time_scheme_type=='IMEXRK') then
         do i=1,Nr_max
            rhs2(i)=0.0_dp
            omg_spec_rad(i)=0.0_dp
         end do

         call chebtransform(Nr_max,omgrhs,omg_spec_rad)
         
         call chebinvtranD1D2(Nr_max,omg_spec_rad,real_d_omg_rad,real_d2_omg_rad)
         
         do i=1,Nr_max
            rhs2(i)= (r_radius(i)*real_d_omg_rad(i)+real_d2_omg_rad(i)-real(Nm,kind=dp)* &
                    & real(Nm,kind=dp)*r_radius2(i)*omgrhs(i)) 
                    ! linear part

         end do
         
      end if

      rhs2(1)=0.0_dp ! Apply BC for stream function in the vorticity equation
      rhs2(Nr_max)=0.0_dp ! Apply BC for stream function in the vorticity equation

      do i=1,Nr_max

         rhsf(i)=rhs1(i)
         rhsf(i+Nr_max)=rhs2(i)

      end do

   end subroutine rhs_construct_vort_d

   subroutine rhs_construct_buo(Nm_max,Nr_max,Ra,Pr,uphi_omg_FR,ur_omg_FR,omg2,Nm, &
                                    & time_scheme_type,rhs2,temprhs) 

      integer, intent(in) :: Nm_max
      integer, intent(in) :: Nr_max 
      character(len=100), intent(in) :: time_scheme_type
      real(kind=dp), intent(in) :: Ra
      real(kind=dp), intent(in) :: Pr  
      complex(kind=dp), intent(in) :: omg2(Nm_max+1,Nr_max)
      complex(kind=dp), intent(in) :: temprhs(Nr_max)
      integer, intent(in) :: Nm
      complex(kind=dp), intent(in) :: uphi_omg_FR(Nm_max+1,Nr_max),ur_omg_FR(Nm_max+1,Nr_max)
      integer :: i
      complex(kind=dp), intent(out) :: rhs2(Nr_max) 

      if (time_scheme_type=='IMEXRK') then
         do i=1,Nr_max
            rhs2(i)=0.0_dp
         end do

         do i=1,Nr_max
            rhs2(i)= - (Ra/Pr)*(ii*real(Nm,kind=dp)*(r_radius(i))*temprhs(i))
                    ! buoyancy term

         end do
         
      end if

      rhs2(1)=0.0_dp ! Apply BC for stream function in the vorticity equation
      rhs2(Nr_max)=0.0_dp ! Apply BC for stream function in the vorticity equation

   end subroutine rhs_construct_buo

end module rhs_create_imexrk
