module fourier

   use, intrinsic :: iso_c_binding

   use double
   use constants, only: ii, pi

   implicit none

   private
   type(C_PTR) :: plan_forward, plan_backward
  
   public :: init_fftwplan, destroy_fftwplan, forfft, invfft

contains

   subroutine init_fftwplan(Np_max)
      include "fftw3.f03"
      integer, intent(in) :: Np_max
      real(kind=dp) :: tf1(Np_max) 
      complex(kind=dp) :: tf2(Np_max/2+1)
      complex(kind=dp) :: tf3(Np_max/2+1)
      real(kind=dp) :: tf4(Np_max)   
        
      plan_forward = fftw_plan_dft_r2c_1d(Np_max,tf1,tf2,FFTW_ESTIMATE)
      plan_backward = fftw_plan_dft_c2r_1d(Np_max,tf3,tf4,FFTW_ESTIMATE)

   end subroutine init_fftwplan

   subroutine destroy_fftwplan()
      include "fftw3.f03"
        
      call fftw_destroy_plan(plan_forward)
      call fftw_destroy_plan(plan_backward)

   end subroutine destroy_fftwplan

   subroutine forfft(Nm_max,Np_max,ff,f)
      include "fftw3.f03"
      integer, intent(in) :: Nm_max
      integer, intent(in) :: Np_max   
      real(kind=dp), intent(in) :: ff(Np_max)
      complex(kind=dp), intent(out) :: f(Nm_max+1)
      integer :: i
      complex(kind=dp) :: tmp1(Np_max/2+1) 
      real(kind=dp) :: fftemp(Np_max)
         
      fftemp=ff ! Make temporary copy as the input array gets destroyed

      call fftw_execute_dft_r2c(plan_forward,fftemp,tmp1)

      do i=1,(Nm_max+1)
         f(i)=tmp1(i)/(real(Np_max,kind=dp))
      end do

   end subroutine forfft

   subroutine forfftO(Nm_max,Np_max,ff,f)
      include "fftw3.f03"
      integer, intent(in) :: Nm_max
      integer, intent(in) :: Np_max   
      real(kind=dp), intent(in) :: ff(Np_max)
      complex(kind=dp), intent(out) :: f(Nm_max+1)
      complex(kind=dp) :: tmp1(Np_max/2+1)
      real(kind=dp) :: fftemp(Np_max)
         
      fftemp=ff

      call fftw_execute_dft_r2c(plan_forward,fftemp,tmp1)
        
      f=tmp1

   end subroutine forfftO

   subroutine invfft(Nm_max,Np_max,f,ff)
      include "fftw3.f03"
      integer, intent(in) :: Nm_max
      integer, intent(in) :: Np_max   
      complex(kind=dp), intent(in) :: f(Nm_max+1)
      real(kind=dp), intent(out) :: ff(Np_max)
      complex(kind=dp) :: ftempfull(Np_max/2+1) 
      integer :: i
      real(kind=dp) :: fftemp1(Np_max)

      do i=1,Nm_max+1
         ftempfull(i)=f(i)
      end do
      do i=Nm_max+2,(Np_max/2+1)
         ftempfull(i)=0.0_dp
      end do
      
      call fftw_execute_dft_c2r(plan_backward,ftempfull,fftemp1)

      ff=fftemp1

   end subroutine invfft 

end module fourier
