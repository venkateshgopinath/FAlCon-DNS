module nonlin
   !$ use OMP_LIB

   use double
   use constants, only: pi
   use fourier, only: forfft, invfft
   use init, only: up, ur, omg, tt, radius, dr

   implicit none

   private

   real(kind=dp), allocatable :: uphi_temp(:), ur_temp(:)
   real(kind=dp), allocatable :: uphi_omg(:), ur_omg(:)
   real(kind=dp), allocatable, public :: dtval_p(:), dtval_r(:), dtval_rkr(:), dtval_rkp(:)

   public ::  allocate_nonlin, deallocate_nonlin, Nr_maxLOOP

contains

   subroutine allocate_nonlin(Np_max,Nr_max) !---------------Allocate-------------------------------------
      integer, intent(in) :: Np_max   
      integer, intent(in) :: Nr_max 

      allocate( dtval_r(Nr_max), dtval_p(Nr_max) )
      allocate( dtval_rkr(Nr_max), dtval_rkp(Nr_max) )

   end subroutine allocate_nonlin

   subroutine deallocate_nonlin() !------------------Deallocate-----------------------------------

      deallocate( dtval_r, dtval_p )
      deallocate( dtval_rkr, dtval_rkp )

   end subroutine deallocate_nonlin

   !----------------------------Get nonlinear products--------------------------------------------
   subroutine  Nr_maxLOOP(Nm_max,Np_max,Nr_max,tFR,omgFR,upFR,urFR,uphi_temp_FR,ur_temp_FR,uphi_omg_FR,ur_omg_FR)

      integer :: i,Nr ! Nr-> radial 'r' loop over Chebyshev points
      integer, intent(in) :: Nr_max 
      integer, intent(in) :: Np_max 
      integer, intent(in) :: Nm_max

      complex(kind=dp), intent(in) :: tFR(Nm_max+1,Nr_max),omgFR(Nm_max+1,Nr_max) 
      complex(kind=dp), intent(in) :: upFR(Nm_max+1,Nr_max), urFR(Nm_max+1,Nr_max)
      complex(kind=dp), intent(out) :: uphi_temp_FR(Nm_max+1,Nr_max),ur_temp_FR(Nm_max+1,Nr_max)
      complex(kind=dp), intent(out) :: uphi_omg_FR(Nm_max+1,Nr_max), ur_omg_FR(Nm_max+1,Nr_max)
      ! Local variables -------------------------------------- 
      real(kind=dp) :: uphi_temp(Np_max), ur_temp(Np_max)
      real(kind=dp) :: uphi_omg(Np_max), ur_omg(Np_max)
      !****************** PRODUCTS in real space (U.grad T) *************************** 
      
      !-------------Nr_max Loop begins -------------------------------------------------------------      
      !!$omp parallel & 
      !!$omp private(Nr,uphi_temp,ur_temp,uphi_omg,ur_omg) default(shared)  
      !!$omp do    
      do Nr=1,Nr_max  
                          
         !Inverse fourier to get back real space variables
         call invfft(Nm_max,Np_max,tFR(:,Nr),tt(:,Nr)) 
         call invfft(Nm_max,Np_max,omgFR(:,Nr),omg(:,Nr)) 
         call invfft(Nm_max,Np_max,upFR(:,Nr),up(:,Nr)) 
         call invfft(Nm_max,Np_max,urFR(:,Nr),ur(:,Nr)) 
          
      !-------------- Get dt values from grid and velocity --------- 

            call get_dtvals(Nr,Nm_max,dtval_r,dtval_p,dtval_rkr,dtval_rkp)

      !-------------------------------------------------------------
         do i=1,Np_max
            uphi_temp(i)=up(i,Nr)*tt(i,Nr)
            ur_temp(i)=ur(i,Nr)*tt(i,Nr)
         end do
              
         call forfft(Nm_max,Np_max,uphi_temp,uphi_temp_FR(:,Nr))
         call forfft(Nm_max,Np_max,ur_temp,ur_temp_FR(:,Nr))
                              
      !***************** PRODUCTS in real space (U.grad Omega)) ********************

         do i=1,Np_max
            uphi_omg(i)=up(i,Nr)*omg(i,Nr)
            ur_omg(i)=ur(i,Nr)*omg(i,Nr)
         end do

         call forfft(Nm_max,Np_max,uphi_omg,uphi_omg_FR(:,Nr))
         call forfft(Nm_max,Np_max,ur_omg,ur_omg_FR(:,Nr))

      end do
      !!$omp end do
      !!$omp end parallel
      !-------------Nr_max Loop ends-----------------------------------------------------------------      

   end subroutine Nr_maxLOOP

   subroutine get_dtvals(Nr,Nm_max,dtval_r,dtval_p,dtval_rkr,dtval_rkp)
   
      integer, intent(in) :: Nr 
      integer, intent(in) :: Nm_max
      real(kind=dp), intent(out) :: dtval_r(:), dtval_p(:)
      real(kind=dp), intent(out) :: dtval_rkr(:), dtval_rkp(:)
      
      dtval_r(Nr) = dr(Nr)/(maxval(abs(real(ur(:,Nr),kind=dp)))) ! (delta(length_radius)/max(velocity))_measure for each radial segment 
      dtval_p(Nr) = (radius(Nr)*pi)/(real(Nm_max,kind=dp)*maxval(abs(real(up(:,Nr),kind=dp)))) !(delta(length_phi)/max(velocity))_measure for each azimuthal segment 
      dtval_rkr(Nr) = dr(Nr)*dr(Nr) ! dr^2  
      dtval_rkp(Nr) = (radius(Nr)*pi)/(real(Nm_max,kind=dp)) * (radius(Nr)*pi)/(real(Nm_max,kind=dp)) ! dphi^2 

   end subroutine get_dtvals

end module nonlin
