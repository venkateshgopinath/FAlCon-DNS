module integ

   use double
   use constants, only: pi
   use namelists, only: rmin
   use chebyshev, only: chebtransform
   !use init, only:radius
   implicit none
   private 

   real(kind=dp), public :: rs

   public :: pInt, radInt,radInt2 

contains

   subroutine pInt(Nr_max, x, y, rs) ! Regular integration on uniform grid 
   
      integer, intent(in) :: Nr_max
      real(kind=dp), intent(in) :: x(Nr_max)        
      real(kind=dp), intent(in) :: y(Nr_max)   
      real(kind=dp), intent(out) :: rs 
      integer :: i
          

      rs=0.0_dp
      do i=1,Nr_max-1
         rs= rs + (x(i+1)-x(i))*(y(i+1)+y(i))
      end do
      rs=0.5_dp*rs
        
   end subroutine pInt



   subroutine radInt(Nr_max,f,ftot) ! Chebyshev integration along radius

      integer, intent(in) :: Nr_max
      real(kind=dp), intent(in) :: f(Nr_max) 
      complex(kind=dp) :: ft(Nr_max)
      real(kind=dp) :: ftr(Nr_max)
      complex(kind=dp) :: ff(Nr_max)
      integer :: i
      real(kind=dp) :: ftot, Qm
      real(kind=dp) :: cc
      real(kind=dp) :: rad(Nr_max)
      
      do i=1,Nr_max
         rad(i)=(rmin + (0.5_dp*(1.0_dp+cos(real(Nr_max-i,kind=dp)*pi/real(Nr_max-1,kind=dp)))))
      end do
      
      do i=1,Nr_max
         ff(i)=f(i)*rad(i)
      end do
      call chebtransform(Nr_max,ff,ft)
      ftr=real(ft,kind=dp)

      ftot=0.0_dp
      do i=0,Nr_max-1
         if (mod(i,2)==0) then
            Qm=-1.0_dp/real(i**2-1,kind=dp)
         else
            Qm=0.0_dp
         end if
         if ((i==0) .or. (i==(Nr_max-1))) then
            cc=0.5_dp*ftr(i+1)*Qm
            ftot=ftot+cc
         else
            cc=ftr(i+1)*Qm
            ftot=ftot+cc
         end if
      end do
  
      ftot=(2.0_dp/real(Nr_max-1,kind=dp))**(0.5_dp)*ftot 

   end subroutine radInt


!------------------------------------------------------------------------------
   subroutine radInt2(Nr_max,f,ftot) ! Chebyshev integration along radius

      integer, intent(in) :: Nr_max
      real(kind=dp), intent(in) :: f(Nr_max) 
      complex(kind=dp) :: ft(Nr_max)
      real(kind=dp) :: ftr(Nr_max)
      complex(kind=dp) :: ff(Nr_max)
      integer :: i
      real(kind=dp) :: ftot, Qm
      real(kind=dp) :: cc
      
      do i=1,Nr_max
         ff(i)=f(i)!*radius(i)*radius(i)
      end do
      call chebtransform(Nr_max,ff,ft)
      ftr=real(ft,kind=dp)
      ftot=0.0_dp
      do i=0,Nr_max-1
         if (mod(i,2)==0) then
            Qm=-1.0_dp/real(i**2-1,kind=dp)
         else
            Qm=0.0_dp
         end if
         if ((i==0) .or. (i==(Nr_max-1))) then
            cc=0.5_dp*ftr(i+1)*Qm
            ftot=ftot+cc
         else
            cc=ftr(i+1)*Qm
            ftot=ftot+cc
         end if
      end do
  
      ftot=(2.0_dp/real(Nr_max-1,kind=dp))**(0.5_dp)*ftot 

   end subroutine radInt2
end module integ
