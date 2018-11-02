module output

   use double
   use constants, only: pi, zero, ii
   use integ, only: pInt, radInt, radInt2 
#if FFTW
   use init, only: radius, r_radius, r_radius2, phi , ur, up, tt, omg, tFR, urFR, upFR, t_an, TFC, &
                   read_checkpoint, init_all_fftw_plans, destroy_all_fftw_plans, tFR, upFR_prev, upFC
#elif GG
   use init, only: radius, r_radius, r_radius2, phi , ur, up, tt, omg, tFR, urFR, upFR, t_an, TFC, &
                   read_checkpoint, tFR, upFR_prev, upFC
#endif
   use chebyshev, only: xp, chebinvtranD1, chebinvtranD2, chebtransform, chebinvtran
   use fourier, only: invfft, forfft
   use timeschemes, only: rhs_update_wts_imp, wt_lhs_tscheme_imp, wt_rhs_tscheme_imp, n_order_tscheme_imp, &
                          & n_order_tscheme_exp, n_order_tscheme_max, dt_array, rhs_imp_temp, rhs_exp_temp, &
                          & rhs_imp_vort, rhs_exp_vort, rhs_imp_uphi_bar, rhs_exp_uphi_bar
   implicit none

   private

   integer :: ke_unit, ke_sp_unit, t_sp_unit, nu_unit, pw_unit, up_unit
   integer, public :: count_snap
   integer, public :: count_chkpnt
   public :: init_output, final_output, writeKE_physical, writeKE_spectral, calculate_spectra, &
              store_checkpoint, store_checkpoint_imexrk, store_snapshot,store_snapshot_imexrk 

contains

   subroutine init_output(tag)
      character(len=100), intent(in) :: tag

      open (newunit=ke_unit,file="KE"//"_"//tag,status="unknown", position="append", action="write")
      open (newunit=ke_sp_unit,file="KEspec"//"_"//tag,status="unknown", position="append", action="write")
      open (newunit=t_sp_unit,file="Temp_spec"//"_"//tag,status="unknown", position="append", action="write")
      open (newunit=nu_unit,file="Nusselt"//"_"//tag,status="unknown", position="append", action="write")
      open (newunit=pw_unit,file="power"//"_"//tag,status="unknown", position="append", action="write")
      open (newunit=up_unit,file="upbalance"//"_"//tag,status="unknown", position="append", action="write")

   end subroutine init_output

   subroutine final_output
      
      close (up_unit)
      close (pw_unit)
      close (ke_sp_unit)
      close (t_sp_unit)
      close (nu_unit)
      close (ke_unit)

   end subroutine final_output
   
   !-------------------------------------------------------------------------------------------------------
   subroutine writeKE_physical(Np_max,Nr_max,Nm_max,tot_time,Ra,Pr) ! 

      integer, intent(in) :: Np_max, Nr_max, Nm_max
      real(kind=dp), intent(in) :: tot_time, Ra, Pr
      real(kind=dp) :: KE_r(Np_max), vis_term_r(Np_max), buo_term_r(Np_max)
      real(kind=dp) :: KE_tot, vis_term_tot, buo_term_tot
      real(kind=dp) :: KE(Np_max,Nr_max), vis_term(Np_max,Nr_max), buo_term(Np_max,Nr_max)
      integer :: i, j, Nr
! alex
      real(kind=dp) :: KE_radial_tot, KE_azimuthal_tot
      real(kind=dp) :: lupbar, KE_radial(Np_max,Nr_max), KE_azimuthal(Np_max,Nr_max)
      real(kind=dp) :: KE_radial_r(Np_max), KE_azimuthal_r(Np_max)
      real(kind=dp) :: l_upbar(Nr_max), e_upbar(Nr_max), KE_0
!end alex
!
      do Nr=1,Nr_max  
         call invfft(Nm_max,Np_max,urFR(:,Nr),ur(:,Nr)) 
         call invfft(Nm_max,Np_max,tFR(:,Nr),tt(:,Nr))
      end do             
      do j=1,Nr_max
         do i=1,Np_max
!alex
            KE_radial(i,j) = ur(i,j)*ur(i,j)  
            KE_azimuthal(i,j) = up(i,j)*up(i,j)
            l_upbar(j) = l_upbar(j) + up(i,j)
            KE(i,j)=KE_radial(i,j) + KE_azimuthal(i,j)
!end alex
            vis_term(i,j)=omg(i,j)*omg(i,j)
            buo_term(i,j)=ur(i,j)*tt(i,j)
         end do
         l_upbar(j) = l_upbar(j) / real(Np_max)
         e_upbar(j) = l_upbar(j) * l_upbar(j)
      end do
      do i=1,Np_max
         call radInt(Nr_max,KE(i,:),KE_r(i)) 
!alex
         call radInt(Nr_max,KE_radial(i,:),KE_radial_r(i)) 
         call radInt(Nr_max,KE_azimuthal(i,:),KE_azimuthal_r(i)) 
         call radInt(Nr_max,e_upbar(:),KE_0)
!end alex
         call radInt(Nr_max,vis_term(i,:),vis_term_r(i)) 
         call radInt(Nr_max,buo_term(i,:),buo_term_r(i)) 
      end do
      
      call pInt(Np_max,phi,KE_r,KE_tot)
!alex
      call pInt(Np_max,phi,KE_radial_r,KE_radial_tot)
      call pInt(Np_max,phi,KE_azimuthal_r,KE_azimuthal_tot)
!end alex
      call pInt(Np_max,phi,vis_term_r,vis_term_tot)
      call pInt(Np_max,phi,buo_term_r,buo_term_tot)
      buo_term_tot = Ra/Pr*buo_term_tot
      KE_tot=0.5_dp*KE_tot
!alex
      KE_radial_tot=0.5_dp*KE_radial_tot
      KE_azimuthal_tot=0.5_dp*KE_azimuthal_tot
      KE_0 = 0.5_dp *  KE_0
!end alex

!     write(ke_unit,'(i8, 4es20.12)') n_step, tot_time, KE_tot, vis_term_tot, buo_term_tot
      write(ke_unit,'(6es20.12)') tot_time, KE_tot, KE_radial_tot, KE_azimuthal_tot, KE_0 ! Write KE data to a text file
!     print *, n_step, tot_time, KE_tot, vis_term_tot, buo_term_tot ! Print KE for checking purposes only
      write(6,'(6es20.12)')  tot_time, KE_tot,  KE_radial_tot, KE_azimuthal_tot, KE_0 ! Print for checking purposes only

   end subroutine writeKE_physical 
   !-------------------------------------------------------------------------------------------------------

   !-------------------------------------------------------------------------------------------------------
   subroutine writeKE_spectral(Nm_max,Nr_max,Np_max,TFC,tot_time,eta,n_step,omgFR,Ra,Pr,dt_new,tFR,mBC) ! 

      integer, intent(in) :: Nm_max, Nr_max, Np_max
      real(kind=dp), intent(in) :: eta, tot_time
      complex(kind=dp), intent(in) :: TFC(Nr_max)
      complex(kind=dp), intent(in) :: omgFR(Nm_max+1,Nr_max), tFR(Nm_max+1,Nr_max)
      integer, intent(in) :: n_step
      real(kind=dp), intent(in) :: Ra, Pr
      real(kind=dp), intent(in) :: dt_new
      character(len=100), intent(in) :: mBC
      real(kind=dp) :: KE_tot
      real(kind=dp) :: Nus_b, Nus_t
      real(kind=dp) :: urk, upk, vis_term, buo_term
      real(kind=dp) :: ur2(Nr_max), up2(Nr_max), w2(Nr_max), URT(Nr_max) 
!alex
      real(kind=dp) :: upbar2(Nr_max), upbark
!end alex
      complex(kind=dp) :: dtFR(Nr_max) 
      integer :: i,j
      complex(kind=dp) :: d1upFR(Nm_max+1,Nr_max)
      complex(kind=dp) :: d2upFR(Nm_max+1,Nr_max)

      complex(kind=dp) :: uFa_eqn(Nr_max)
      complex(kind=dp) :: u_adv(Nr_max)
      complex(kind=dp) :: u_diff(Nr_max)

      real(kind=dp) :: omg(Np_max,Nr_max)
      real(kind=dp) :: d1up(Np_max,Nr_max)
      real(kind=dp) :: ur_t(Nr_max)
      complex(kind=dp) :: ur_up(Nr_max), upavg(Nr_max)
      complex(kind=dp) :: ur_omg(Nr_max)
      complex(kind=dp) :: ur_up_FR(Nr_max)
      complex(kind=dp) :: ur_d1up_FR(Nr_max)
      integer :: Nr,Nm, Np

      complex(kind=dp), allocatable :: tFR_ref(:,:), omgFR_ref(:,:), urFR_ref(:,:), upFR_ref(:,:) ! Previously saved data (small 'c' denotes current timestep)
      real(kind=dp), allocatable :: tt_comp(:,:)
      complex(kind=dp), allocatable :: tFC_comp(:,:)
      complex(kind=dp), allocatable :: tFR_comp(:,:)
      real(kind=dp) :: eta_ref, Ra_ref, Pr_ref 
      integer :: ref_point,ref_point2, Np_max_ref
      real(kind=dp), allocatable :: tt_ref(:,:)
      real(kind=dp), allocatable :: temp_diff_phi(:), ur_diff_phi(:)
      real(kind=dp), allocatable :: ref_diff_phi(:)
      real(kind=dp) :: L2error_temp, L2error_ur
      complex(kind=dp) :: ttFC(Nm_max+1,Nr_max)
      real(kind=dp), allocatable :: rad(:)
      real(kind=dp) :: dt_new_ref, tot_time_ref
      integer :: Nm_max_ref, Nr_max_ref

      ! Check u_phi force balance  
      do Nm=1,Nm_max+1 
         call chebinvtranD1(Nr_max,upFC(Nm,:),d1upFR(Nm,:))
         call chebinvtranD2(Nr_max,upFC(Nm,:),d2upFR(Nm,:))
      end do

      ur_up_FR(:)=0.0_dp
      ur_d1up_FR(:)=0.0_dp
      u_adv(:)=0.0_dp

      ! ---- Product using transformed variables urFR and upFR --------------------

      do Nr=1,Nr_max
         do Nm=0,Nm_max
            if (Nm==0) then
               ur_up_FR(Nr) = ur_up_FR(Nr) + (urFR(Nm+1,Nr)*(upFR(Nm+1,Nr)))
               ur_d1up_FR(Nr) = ur_d1up_FR(Nr) + (urFR(Nm+1,Nr)*(d1upFR(Nm+1,Nr)))
            else
               ur_up_FR(Nr)=ur_up_FR(Nr) + (urFR(Nm+1,Nr)*conjg(upFR(Nm+1,Nr)) &
                           & + conjg(urFR(Nm+1,Nr))*upFR(Nm+1,Nr))
               ur_d1up_FR(Nr)=ur_d1up_FR(Nr) + urFR(Nm+1,Nr)*conjg(d1upFR(Nm+1,Nr)) &
                           & + conjg(urFR(Nm+1,Nr))*d1upFR(Nm+1,Nr)
            end if 
         end do
      end do

      ! ---- Product using physical space variables ur and up
      do Nr=1,Nr_max
         call invfft(Nm_max,Np_max,upFR(:,Nr),up(:,Nr)) 
         call invfft(Nm_max,Np_max,urFR(:,Nr),ur(:,Nr)) 
         call invfft(Nm_max,Np_max,omgFR(:,Nr),omg(:,Nr)) 
         call invfft(Nm_max,Np_max,d1upFR(:,Nr),d1up(:,Nr)) 
      end do

      upavg(:)=0.0_dp
      ur_up(:)=0.0_dp
      ur_t(:)=0.0_dp
      ur_omg(:)=0.0_dp

      do Nr=1,Nr_max
         do Np=1,Np_max
            upavg(Nr)=upavg(Nr) + up(Np,Nr)
            ur_up(Nr)=ur_up(Nr) + ur(Np,Nr)*up(Np,Nr)
            ur_t(Nr)=ur_t(Nr) + ur(Np,Nr)*tt(Np,Nr)
            ur_omg(Nr)=ur_omg(Nr) + ur(Np,Nr)*omg(Np,Nr)
         end do
      end do
         upavg=upavg*(1.0_dp/Np_max)
         ur_up=ur_up*(1.0_dp/Np_max)
         ur_t=ur_t*(1.0_dp/Np_max)
         ur_omg=ur_omg*(1.0_dp/Np_max)
      ! ---------------------------------------------------------------------------

      do Nr=1,Nr_max
         uFa_eqn(Nr)=-(ur_d1up_FR(Nr)+r_radius(Nr)*ur_up_FR(Nr))
      end do

      do Nr=1,Nr_max
         u_diff(Nr) = d2upFR(1,Nr) - r_radius2(Nr)*upFR(1,Nr) + r_radius(Nr)*d1upFR(1,Nr)
      end do

      !print *, real(u_adv(20)) , real(u_diff(20)), upFR(1,20), "check ur omg"

      if (mod(n_step,1000)==0) then
         do Nr=1,Nr_max
            write(up_unit,'(4ES20.12)') radius(Nr), real((upFR(1,Nr)-upFR_prev(1,Nr))/(dt_new)), &
                  & real(uFa_eqn(Nr)), real(u_diff(Nr))
            !write(up_unit,'(4ES20.12)') radius(Nr), (real(ur_up_FR(Nr))), real(ur_up(Nr)), real(ur_up_FR(Nr))-real(ur_up(Nr)) ! Check if products are same in Physical and Fourier space
         end do
      end if

      !------ Calculate KE ------------  
      do j=1,Nr_max
          ur2(j)=zero
          up2(j)=zero
          do i=1,Nm_max+1
             ur2(j)=ur2(j)+abs(urFR(i,j)*urFR(i,j))
             up2(j)=up2(j)+abs(upFR(i,j)*upFR(i,j))
          end do
!alex
          upbar2(j) = 0.0_dp
          upbar2(j) = upavg(j)*upavg(j)
!end alex
      end do

      call radInt(Nr_max,ur2,urk)
      call radInt(Nr_max,up2,upk)
!alex
      call radInt(Nr_max,upbar2,upbark)
!end alex
       
      KE_tot=2.0_dp*pi*(urk+upk)
      
      !------------------------ Write KE and Nusselt number data in a text file -------------------
!     write(ke_unit,'(6es20.12)') dt_new, tot_time, KE_tot, maxval(real(upFR(1,:)))  ! Write KE data to a text file
!     write(ke_unit,'(6es20.12)') dt_new, tot_time, KE_tot, maxval(real(tFR(1,2:62)))  ! Write KE data to a text file (test for adv-diff)
      write(ke_unit,'(6es20.12)')  dt_new, tot_time, KE_tot, 2.0_dp * pi * urk, 2.0_dp * pi * upk, 2.0_dp * pi * upbark  ! Write KE data to a text file
      !-----------------------------------

      !------ Calculate Nusselt ---------
      call chebinvtranD1(Nr_max,TFC,dtFR)
                 
      Nus_b = real(dtFR(1),kind=dp)/((1.0_dp-eta)/(eta*log(eta)))
      Nus_t = real(dtFR(Nr_max),kind=dp)/((1.0_dp-eta)/(log(eta)))

      write(nu_unit,'(3ES20.12)') tot_time, Nus_t, Nus_b ! Write Nusselt no. data to text file 
      
      !-----------------------------------
      
      !------ Calculate Power budget ---------
      ! Viscous diffusion -------------------- 
      do j=1,Nr_max
         w2(j)=zero
         do i=0,Nm_max
            if (i==0) then
               w2(j) = w2(j)+real(omgFR(i+1,j)*omgFR(i+1,j))
            else 
               w2(j)=w2(j)+real(omgFR(i+1,j)*conjg(omgFR(i+1,j))+conjg(omgFR(i+1,j))*omgFR(i+1,j))
            end if
         end do
      end do   
      call radInt(Nr_max,w2,vis_term)
      if (mBC=='NS') then
         vis_term=2.0_dp*pi*vis_term
      elseif (mBC=='SF') then
         vis_term=2.0_dp*pi*vis_term
      end if     
      ! Buoyancy -----------------------------
      do j=1,Nr_max
         URT(j)=zero
         do i=0,Nm_max
            if (i==0) then
               URT(j) = URT(j)+real(urFR(i+1,j)*tFR(i+1,j))
            else 
               URT(j)=URT(j)+real(urFR(i+1,j)*conjg(tFR(i+1,j))+conjg(urFR(i+1,j)*tFR(i+1,j)))
            end if
         end do
      end do  
      
      !call radInt(Nr_max,URT,buo_term) ! Calculating from Spectral space
      call radInt(Nr_max,ur_t,buo_term) ! Calculating from Physical space

      buo_term=2.0_dp*pi*Ra/Pr*buo_term

      write(pw_unit,'(3ES20.12)') tot_time, vis_term, buo_term ! Write power budget data to text file
      !----------------------------------------
!-----------------------------------------------------------------------------------
     print *, dt_new, tot_time, KE_tot, maxval(aimag(upFR(:,1))), maxval(real(upFR(:,1))) ! Print for checking purposes only
!      write(6,'(6es20.12)')  tot_time, KE_tot, 2.0_dp * pi * urk, 2.0_dp * pi * upk, 2.0_dp * pi * upbark  ! Write KE data to a text file
!#####################################  

   end subroutine writeKE_spectral 
   !--------------------------------------------------------------------------------------------------------

   subroutine calculate_spectra(Nm_max,Nr_max,urFR,upFR,n_step) ! 

      integer, intent(in) :: Nm_max, Nr_max, n_step
      complex(kind=dp), intent(in) :: urFR(Nm_max+1,Nr_max),upFR(Nm_max+1,Nr_max)
      real(kind=dp) ::  KEspec(Nm_max+1), Ekin(Nr_max), Tspec(Nm_max+1), Tsp_rad(Nr_max)
      integer :: i,m,r

      !-------------------------- KE at each mode Nm_max (Spectra)  -------------------------------- 

      do m=0,Nm_max
         do r=1,Nr_max
            Ekin(r) = abs((urFR(m+1,r)*urFR(m+1,r)) + (upFR(m+1,r)*upFR(m+1,r)))
         end do
         call radInt(Nr_max,Ekin,KEspec(m+1))
      end do

      do m=0,Nm_max
         do r=1,Nr_max
            Tsp_rad(r) = abs(tFR(m+1,r))
         end do
         call radInt(Nr_max,Tsp_rad,Tspec(m+1))
      end do

      do i=0,Nm_max
         write(ke_sp_unit,'(2i8, 1ES20.12)') n_step, i, 4.0_dp*pi*KEspec(i+1) ! 
         write(t_sp_unit,'(2i8, 1ES20.12)') n_step, i, 4.0_dp*pi*Tspec(i+1) ! 
      end do
                 
   end subroutine calculate_spectra
   !--------------------------------------------------------------------------------------------------------

   !--------------------------------------------------------------------------------------------------------
   subroutine store_checkpoint(Nm_max,Nr_max,count_chkpnt,dt_new,tot_time,tFR,omgFR,urFR,upFR, &
                                  & n_order_tscheme_imp,n_order_tscheme_exp, rhs_imp_temp, &
                                  & rhs_exp_temp,rhs_imp_vort,rhs_exp_vort,rhs_imp_uphi_bar, rhs_exp_uphi_bar, &
                                  & dt_array,n_order_tscheme_max,time_scheme_type)

      integer, intent(in) :: Nm_max   
      integer, intent(in) :: Nr_max 
      integer, intent(in) :: count_chkpnt 
      integer, intent(in) :: n_order_tscheme_imp,n_order_tscheme_exp,n_order_tscheme_max 
      character(len=100), intent(in) :: time_scheme_type
      real(kind=dp), intent(in) :: tot_time, dt_new
      real(kind=dp), intent(in) :: dt_array(n_order_tscheme_max)
      complex(kind=dp), intent(in) :: tFR(Nm_max+1,Nr_max)
      complex(kind=dp), intent(in) :: omgFR(Nm_max+1,Nr_max)
      complex(kind=dp), intent(in) :: urFR(Nm_max+1,Nr_max)
      complex(kind=dp), intent(in) :: upFR(Nm_max+1,Nr_max)
      complex(kind=dp), intent(in) :: rhs_imp_temp(n_order_tscheme_imp,Nm_max+1,Nr_max)
      complex(kind=dp), intent(in) :: rhs_exp_temp(n_order_tscheme_exp,Nm_max+1,Nr_max)
      complex(kind=dp), intent(in) :: rhs_imp_vort(n_order_tscheme_imp,Nm_max+1,Nr_max)
      complex(kind=dp), intent(in) :: rhs_exp_vort(n_order_tscheme_exp,Nm_max+1,Nr_max)
      complex(kind=dp), intent(in) :: rhs_imp_uphi_bar(n_order_tscheme_imp,Nr_max)
      complex(kind=dp), intent(in) :: rhs_exp_uphi_bar(n_order_tscheme_exp,Nr_max)

      character(len=72) :: datafile1
      integer :: outunit 
      
     
      !-----------------------------------------Write current time step -----------------------------------------
      write(datafile1,fmt='(a,I5.5)') "checkpoint_", count_chkpnt  
      open(newunit=outunit, file=datafile1, action="write", status="replace", form='unformatted' )
      write(outunit) dt_new,tot_time,Nm_max,Nr_max,n_order_tscheme_imp,n_order_tscheme_exp, &
                     & n_order_tscheme_max
      write(outunit) dt_array
      write(outunit) tFR
      write(outunit) omgFR
      write(outunit) urFR
      write(outunit) upFR
      write(outunit) rhs_imp_temp
      write(outunit) rhs_exp_temp
      write(outunit) rhs_imp_vort
      write(outunit) rhs_exp_vort
      write(outunit) rhs_imp_uphi_bar
      write(outunit) rhs_exp_uphi_bar
      close(outunit)
   
   end subroutine store_checkpoint
   !-------------------------------------------------------------------------------------------------------
   !--------------------------------------------------------------------------------------------------------
   subroutine store_checkpoint_imexrk(Nm_max,Nr_max,count_chkpnt,dt_new,tot_time,tFR,omgFR,urFR,upFR, &
                                  & n_order_tscheme_imp,n_order_tscheme_exp, rhs_imp_temp, &
                                  & rhs_exp_temp,rhs_imp_vort,rhs_exp_vort,rhs_imp_uphi_bar, rhs_exp_uphi_bar, &
                                  & dt_array,n_order_tscheme_max,time_scheme_type)

      integer, intent(in) :: Nm_max   
      integer, intent(in) :: Nr_max 
      integer, intent(in) :: count_chkpnt 
      integer, intent(in) :: n_order_tscheme_imp,n_order_tscheme_exp,n_order_tscheme_max 
      character(len=100), intent(in) :: time_scheme_type
      real(kind=dp), intent(in) :: tot_time, dt_new
      real(kind=dp), intent(in) :: dt_array(n_order_tscheme_max)
      complex(kind=dp), intent(in) :: tFR(Nm_max+1,Nr_max)
      complex(kind=dp), intent(in) :: omgFR(Nm_max+1,Nr_max)
      complex(kind=dp), intent(in) :: urFR(Nm_max+1,Nr_max)
      complex(kind=dp), intent(in) :: upFR(Nm_max+1,Nr_max)
      complex(kind=dp), intent(in) :: rhs_imp_temp(1,Nm_max+1,Nr_max)
      complex(kind=dp), intent(in) :: rhs_exp_temp(1,Nm_max+1,Nr_max)
      complex(kind=dp), intent(in) :: rhs_imp_vort(1,Nm_max+1,Nr_max)
      complex(kind=dp), intent(in) :: rhs_exp_vort(1,Nm_max+1,Nr_max)
      complex(kind=dp), intent(in) :: rhs_imp_uphi_bar(1,Nr_max)
      complex(kind=dp), intent(in) :: rhs_exp_uphi_bar(1,Nr_max)

      character(len=72) :: datafile1
      integer :: outunit 
      
     
      !-----------------------------------------Write current time step -----------------------------------------
      write(datafile1,fmt='(a,I5.5)') "checkpoint_", count_chkpnt  
      open(newunit=outunit, file=datafile1, action="write", status="replace", form='unformatted' )
      write(outunit) dt_new,tot_time,Nm_max,Nr_max,n_order_tscheme_imp,n_order_tscheme_exp, &
                     & n_order_tscheme_max
      write(outunit) dt_array
      write(outunit) tFR
      write(outunit) omgFR
      write(outunit) urFR
      write(outunit) upFR
      write(outunit) rhs_imp_temp
      write(outunit) rhs_exp_temp
      write(outunit) rhs_imp_vort
      write(outunit) rhs_exp_vort
      write(outunit) rhs_imp_uphi_bar
      write(outunit) rhs_exp_uphi_bar
      close(outunit)
   
   end subroutine store_checkpoint_imexrk
   !-------------------------------------------------------------------------------------------------------

   !---------------------------------------------------------------------------------------------------------
   subroutine store_snapshot_imexrk(Nm_max,Nr_max,Ra,Pr,eta,tot_time,dt_new,count_snap,tFRs,omgFRs,urFRs,upFRs,omgFR_check_s)

      integer, intent(in) :: Nm_max
      integer, intent(in) :: Nr_max
      real(kind=dp), intent(in) :: Ra, Pr, eta
      real(kind=dp), intent(in) :: tot_time,dt_new
      integer, intent(in) :: count_snap
      complex(kind=dp), intent(in) :: tFRs(Nm_max+1,Nr_max)
      complex(kind=dp), intent(in) :: omgFRs(Nm_max+1,Nr_max)
      complex(kind=dp), intent(in) :: urFRs(Nm_max+1,Nr_max)
      complex(kind=dp), intent(in) :: upFRs(Nm_max+1,Nr_max)
      complex(kind=dp), intent(in) :: omgFR_check_s(Nm_max+1,Nr_max)

      character(len=72) :: datafile1
      integer :: outunit, i 
       
      !do i=1,1
      !   print *, radius(i), maxval(abs(omgFRs(:,1)))
      !end do
      
      !-----------------------------------------Write current time step ---------------------------------------
      write(datafile1,fmt='(a,I5.5)') "snapshot_plot_", count_snap
      open(newunit=outunit, file=datafile1, action="write", status="replace", form='unformatted' )
      write(outunit) Ra, Pr, eta
      write(outunit) dt_new,tot_time
      write(outunit) Nm_max,Nr_max
      write(outunit) radius
      write(outunit) tFRs
      write(outunit) omgFRs
      write(outunit) urFRs
      write(outunit) upFRs
      write(outunit) omgFR_check_s
      close(outunit)

   end subroutine store_snapshot_imexrk
   !-------------------------------------------------------------------------------------------------------

   !---------------------------------------------------------------------------------------------------------
   subroutine store_snapshot(Nm_max,Nr_max,Ra,Pr,eta,tot_time,dt_new,count_snap,tFRs,omgFRs,urFRs,upFRs)

      integer, intent(in) :: Nm_max
      integer, intent(in) :: Nr_max
      real(kind=dp), intent(in) :: Ra, Pr, eta
      real(kind=dp), intent(in) :: tot_time,dt_new
      integer, intent(in) :: count_snap
      complex(kind=dp), intent(in) :: tFRs(Nm_max+1,Nr_max)
      complex(kind=dp), intent(in) :: omgFRs(Nm_max+1,Nr_max)
      complex(kind=dp), intent(in) :: urFRs(Nm_max+1,Nr_max)
      complex(kind=dp), intent(in) :: upFRs(Nm_max+1,Nr_max)

      character(len=72) :: datafile1
      integer :: outunit, i 
       
      !do i=1,1
      !   print *, radius(i), maxval(abs(omgFRs(:,1)))
      !end do
      
      !-----------------------------------------Write current time step ---------------------------------------
      write(datafile1,fmt='(a,I5.5)') "snapshot_plot_", count_snap
      open(newunit=outunit, file=datafile1, action="write", status="replace", form='unformatted' )
      write(outunit) Ra, Pr, eta
      write(outunit) dt_new,tot_time
      write(outunit) Nm_max,Nr_max
      write(outunit) radius
      write(outunit) tFRs
      write(outunit) omgFRs
      write(outunit) urFRs
      write(outunit) upFRs
      close(outunit)

   end subroutine store_snapshot
   !-------------------------------------------------------------------------------------------------------

  subroutine read_snapshot(dt_new,tot_time,eta,Ra,Pr,Nm_max,Nr_max,rad,tFRn,omgFRn,urFRn,upFRn,ref_point)  

      integer, intent(in) :: ref_point
      real(kind=dp), intent(out) :: dt_new
      real(kind=dp), intent(out) :: tot_time, eta, Ra, Pr
      integer, intent(out) :: Nm_max   
      integer, intent(out) :: Nr_max 
      complex(kind=dp), allocatable, intent(out) :: tFRn(:,:)
      complex(kind=dp), allocatable, intent(out) :: omgFRn(:,:)
      complex(kind=dp), allocatable, intent(out) :: urFRn(:,:)
      complex(kind=dp), allocatable, intent(out) :: upFRn(:,:)
      real(kind=dp), allocatable, intent(out) :: rad(:)

      character(len=72) :: datafile1
      integer :: inunit
      !----------------------------- Read latest files -------------------------------
      write(datafile1,'("snapshot_plot_ref_",I5.5)') ref_point
      open(newunit=inunit, status='old',file=datafile1,form='unformatted')
      read(inunit) Ra, Pr, eta
      read(inunit) dt_new, tot_time 
      read(inunit) Nm_max,Nr_max
      allocate(rad(Nr_max))
      allocate(tFRn(Nm_max+1,Nr_max))
      allocate(omgFRn(Nm_max+1,Nr_max))
      allocate(urFRn(Nm_max+1,Nr_max))
      allocate(upFRn(Nm_max+1,Nr_max))

      read(inunit) rad
      read(inunit) tFRn
      read(inunit) omgFRn
      read(inunit) urFRn
      read(inunit) upFRn
      close(inunit)

      end subroutine read_snapshot


      subroutine read_snapshot2(tFRn,omgFRn,urFRn,upFRn,ref_point)  

      integer, intent(in) :: ref_point
      real(kind=dp):: dt_new
      real(kind=dp) :: tot_time, eta, Ra, Pr
      integer :: Nm_max   
      integer :: Nr_max 
      complex(kind=dp), allocatable, intent(out) :: tFRn(:,:)
      complex(kind=dp), allocatable, intent(out) :: omgFRn(:,:)
      complex(kind=dp), allocatable, intent(out) :: urFRn(:,:)
      complex(kind=dp), allocatable, intent(out) :: upFRn(:,:)

      character(len=72) :: datafile1
      integer :: inunit
      !----------------------------- Read latest files -------------------------------
      write(datafile1,'("snapshot_plot_",I5.5)') ref_point
      open(newunit=inunit, status='old',file=datafile1,form='unformatted')
      read(inunit) Ra, Pr, eta
      read(inunit) dt_new, tot_time 
      read(inunit) Nm_max,Nr_max
      allocate(tFRn(Nm_max+1,Nr_max))
      allocate(omgFRn(Nm_max+1,Nr_max))
      allocate(urFRn(Nm_max+1,Nr_max))
      allocate(upFRn(Nm_max+1,Nr_max))

      read(inunit) tFRn
      read(inunit) omgFRn
      read(inunit) urFRn
      read(inunit) upFRn
      close(inunit)

      end subroutine read_snapshot2
end module output

