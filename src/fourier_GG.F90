module fourier

   use double
   use constants, only: ii, pi

   implicit none

   private

   public :: forfft, invfft

contains

   subroutine forfft(Nm_max,Np_max,ff,f)
          
      integer, intent(in) :: Nm_max
      integer, intent(in) :: Np_max   
      
      real(kind=dp), intent(in) :: ff(Np_max)
      complex(kind=dp), intent(out) :: f(Nm_max+1)
      complex(kind=dp) :: cc 
      integer :: i,j,m
      
      do i=1,Nm_max+1
         f(i)=0.0_dp
      end do
      
      
      do m=0,Nm_max
         do j=1,Np_max
            cc=ff(j)*exp(-2.0_dp*pi*ii*real(m,kind=dp)*((real(j,kind=dp)-1.0_dp)/real(Np_max,kind=dp)))
            f(m+1)=f(m+1)+cc
         end do
      end do
       
      f=(1.0_dp/real(Np_max,kind=dp))*f 
          
   end subroutine forfft

   subroutine invfft(Nm_max,Np_max,f, ff)
          
      integer, intent(in) :: Nm_max
      integer, intent(in) :: Np_max   
      
      complex(kind=dp), intent(in) :: f(Nm_max+1)
      real(kind=dp), intent(out) :: ff(Np_max)
      real(kind=dp) :: cc
      integer :: i,j,m
          

      do i=1,Np_max
         ff(i)=0.0_dp
      end do
   

   
      do m=0,Nm_max
         do j=1,Np_max
            if (m==0) then
               cc=0.5_dp*real(((f(m+1))*exp(2.0_dp*pi*ii*(m)*(j-1)/(Np_max))),kind=dp)
               ff(j)=ff(j)+cc
            else
               cc=real(((f(m+1))*exp(2.0_dp*pi*ii*real(m)*(j-1)/(Np_max))),kind=dp)
               ff(j)=ff(j)+cc 
            end if
         end do
      end do
   
      ff=2.0_dp*ff 
           
   end subroutine invfft 

end module fourier
