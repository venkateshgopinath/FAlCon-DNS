module chebyshev

   use double
   use constants, only: pi

   implicit none

   private

   real(kind=dp), allocatable, public :: D(:,:), D2(:,:), t(:,:) 
   complex(kind=dp), allocatable, public :: xp(:)

   public :: cheballoc, chebdealloc, dcheb, chebtransform, chebinvtran, chebinvtranD1, &
             & chebinvtranD2, chebinvtranD1D2, lagrange, dlagrange, d2lagrange 

contains

   subroutine cheballoc(Nr_max) 
        
      integer, intent(in) :: Nr_max 
      
      allocate( D(Nr_max,Nr_max), D2(Nr_max,Nr_max), t(Nr_max,Nr_max), xp(Nr_max) )

   end subroutine cheballoc

   subroutine chebdealloc()

      deallocate( D, D2, t, xp )

   end subroutine chebdealloc
     
   subroutine dcheb(Nr_max)

      integer, intent(in) :: Nr_max 
      real(kind=dp) :: dt(Nr_max,Nr_max), d2t(Nr_max,Nr_max)
      integer, parameter :: dm=1.0_dp
      integer :: i,j

      do i=1,Nr_max
         xp(i) = cos((real(Nr_max-i,kind=dp)*pi)/real(Nr_max-1,kind=dp)) 
      end do
      
      do i=1,Nr_max
         t(1,i)=1.0_dp
         t(2,i)=xp(i)
         dt(1,i)=0.0_dp
         dt(2,i)=1.0_dp
         d2t(1,i)=0.0_dp
         d2t(2,i)=0.0_dp
         do j=3,Nr_max
            t(j,i) = 2.0_dp*xp(i)*t(j-1,i)-t(j-2,i)
            dt(j,i) = 2.0_dp*t(j-1,i) + 2.0_dp*xp(i)*dt(j-1,i) - dt(j-2,i)
            d2t(j,i) = 4.0_dp*dt(j-1,i) + 2.0_dp*xp(i)*d2t(j-1,i)- d2t(j-2,i)
         end do
      end do

      D=2.0_dp/dm * dt

      D=transpose(D)
      t=transpose(t)
      D2 = (2.0_dp/dm)*(2.0_dp/dm)*d2t
      D2 = transpose(D2)

   end subroutine dcheb

   subroutine chebtransform(Nr_max,f,ft)

      integer, intent(in) :: Nr_max 
      integer :: m,k,i
      complex(kind=dp), intent(in) :: f(Nr_max)
      complex(kind=dp), intent(out) :: ft(Nr_max)
      complex(kind=dp), allocatable :: d(:)
      complex(kind=dp) :: Tm, c
           
      allocate( d(Nr_max) )

      do i=1,Nr_max
         d(i)=0.0_dp
      end do

      do m=1,Nr_max
         do k=1,Nr_max
            Tm=cos((m-1)*real(Nr_max-k,kind=dp)*pi/real(Nr_max-1,kind=dp))
            if ((k==1) .or. (k==Nr_max)) then
               c=0.5_dp*f(k)*Tm
               d(m)=d(m)+c
            else
               c=f(k)*Tm
               d(m)=d(m)+c
            end if
         end do
      end do 

      ft=((2.0_dp/real(Nr_max-1.0_dp,kind=dp))**(0.5_dp))*d

   end subroutine chebtransform

   subroutine chebinvtran(Nr_max,f,ff)

      integer :: m,k,i,j
      integer, intent(in) :: Nr_max 
      complex(kind=dp), intent(in) :: f(Nr_max)
      complex(kind=dp), intent(out) :: ff(Nr_max)
      complex(kind=dp), allocatable :: d(:), kk(:), d2(:), xp(:)
      complex(kind=dp), allocatable :: t(:,:), dt(:,:), d2t(:,:)
      complex(kind=dp) :: dm, Tm, c, cc, c2

      allocate( d(Nr_max), d2(Nr_max), kk(Nr_max), xp(Nr_max) )
      allocate( t(Nr_max,Nr_max), dt(Nr_max,Nr_max), d2t(Nr_max,Nr_max) )

      dm=1.0_dp

      do i=1,Nr_max
         kk(i)=i
      end do

      do i=1,Nr_max
         xp(i)=(cos((real(Nr_max-kk(i),kind=dp)*pi)/real(Nr_max-1,kind=dp)))
      end do
      do k=1,Nr_max
         t(1,k)=1.0_dp
         t(2,k)=xp(k)
         dt(1,k)=0.0_dp
         dt(2,k)=1.0_dp
         d2t(1,k)=0.0_dp
         d2t(2,k)=0.0_dp
         do m=3,Nr_max
            t(m,k)=cos((m-1)*real(Nr_max-k,kind=dp)*pi/real(Nr_max-1,kind=dp))
            dt(m,k)=2.0_dp*t(m-1,k) + 2.0_dp*xp(k)*dt(m-1,k) - dt(m-2,k)
            d2t(m,k) = 4.0_dp*dt(m-1,k) + 2.0_dp*xp(k)*d2t(m-1,k) - d2t(m-2,k)
         end do
      end do

      dt=2.0_dp/dm*dt
      d2t=(2.0_dp/dm)**2*d2t

      do j=1,Nr_max
         d(j)=0.0_dp
         ff(j)=0.0_dp
         d2(j)=0.0_dp
      end do

      do k=1,Nr_max
         do m=1,Nr_max
            if  ((m==1) .or. (m==Nr_max)) then
               c2=0.5_dp*(f(m)*d2t(m,k))
               d2(k)=d2(k)+c2
               c=0.5_dp*(f(m)*dt(m,k))
               d(k)=d(k)+c
               cc=0.5_dp*(f(m)*t(m,k))
               ff(k)=ff(k)+cc
            else
               c2=f(m)*d2t(m,k)
               d2(k)=d2(k)+c2
               c=f(m)*dt(m,k)
               d(k)=d(k)+c
               cc=(f(m)*t(m,k))
               ff(k)=ff(k)+cc
            end if
         end do
      end do

      ff=((2.0_dp/real(Nr_max-1,kind=dp))**(0.5_dp)) * ff

   end subroutine chebinvtran

   subroutine chebinvtranD1(Nr_max,f,df)

      integer :: m,k,i,j
      integer, intent(in) :: Nr_max 
      complex(kind=dp), intent(in) :: f(Nr_max)
      complex(kind=dp), intent(out) :: df(Nr_max)
      complex(kind=dp), allocatable :: d(:), kk(:), d2(:), xp(:)
      complex(kind=dp), allocatable :: t(:,:), dt(:,:), d2t(:,:)
      complex(kind=dp) :: dm, Tm, c, cc, c2

      allocate( d(Nr_max), d2(Nr_max), kk(Nr_max), xp(Nr_max) )
      allocate( t(Nr_max,Nr_max), dt(Nr_max,Nr_max), d2t(Nr_max,Nr_max) )

      dm=1.0_dp

      do i=1,Nr_max
         kk(i)=i
      end do

      do i=1,Nr_max
         xp(i)=(cos((real(Nr_max-kk(i),kind=dp)*pi)/real(Nr_max-1,kind=dp)))
      end do
      do k=1,Nr_max
         t(1,k)=1.0_dp
         t(2,k)=xp(k)
         dt(1,k)=0.0_dp
         dt(2,k)=1.0_dp
         d2t(1,k)=0.0_dp
         d2t(2,k)=0.0_dp
         do m=3,Nr_max
            t(m,k)=cos((m-1)*real(Nr_max-k,kind=dp)*pi/real(Nr_max-1,kind=dp))
            dt(m,k)=2.0_dp*t(m-1,k) + 2.0_dp*xp(k)*dt(m-1,k) - dt(m-2,k)
            d2t(m,k) = 4.0_dp*dt(m-1,k) + 2.0_dp*xp(k)*d2t(m-1,k) - d2t(m-2,k)
         end do
      end do

      dt=2.0_dp/dm*dt
      d2t=(2.0_dp/dm)**2*d2t

      do j=1,Nr_max
         d(j)=0.0_dp
         d2(j)=0.0_dp
      end do

      do k=1,Nr_max
         do m=1,Nr_max
            if  ((m==1) .or. (m==Nr_max)) then
               c2=0.5_dp*(f(m)*d2t(m,k))
               d2(k)=d2(k)+c2
               c=0.5_dp*(f(m)*dt(m,k))
               d(k)=d(k)+c
               cc=0.5_dp*(f(m)*t(m,k))
            else
               c2=f(m)*d2t(m,k)
               d2(k)=d2(k)+c2
               c=f(m)*dt(m,k)
               d(k)=d(k)+c
               cc=(f(m)*t(m,k))
            end if
         end do
      end do

      df=((2.0_dp/real(Nr_max-1,kind=dp))**(0.5_dp)) * d

   end subroutine chebinvtranD1

   subroutine chebinvtranD2(Nr_max,f,d2f)

      integer :: m,k,i,j
      integer, intent(in) :: Nr_max 
      complex(kind=dp), intent(in) :: f(Nr_max)
      complex(kind=dp), intent(out) :: d2f(Nr_max)
      complex(kind=dp), allocatable :: d(:), kk(:), d2(:), xp(:)
      complex(kind=dp), allocatable :: t(:,:), dt(:,:), d2t(:,:)
      complex(kind=dp) :: dm, Tm, c, cc, c2

      allocate( d(Nr_max), d2(Nr_max), kk(Nr_max), xp(Nr_max) )
      allocate( t(Nr_max,Nr_max), dt(Nr_max,Nr_max), d2t(Nr_max,Nr_max) )

      dm=1.0_dp

      do i=1,Nr_max
         kk(i)=i
      end do

      do i=1,Nr_max
         xp(i)=(cos((real(Nr_max-kk(i),kind=dp)*pi)/real(Nr_max-1,kind=dp)))
      end do
      do k=1,Nr_max
         t(1,k)=1.0_dp
         t(2,k)=xp(k)
         dt(1,k)=0.0_dp
         dt(2,k)=1.0_dp
         d2t(1,k)=0.0_dp
         d2t(2,k)=0.0_dp
         do m=3,Nr_max
            t(m,k)=cos((m-1)*real(Nr_max-k,kind=dp)*pi/real(Nr_max-1,kind=dp))
            dt(m,k)=2.0_dp*t(m-1,k) + 2.0_dp*xp(k)*dt(m-1,k) - dt(m-2,k)
            d2t(m,k) = 4.0_dp*dt(m-1,k) + 2.0_dp*xp(k)*d2t(m-1,k) - d2t(m-2,k)
         end do
      end do

      dt=2.0_dp/dm*dt
      d2t=(2.0_dp/dm)**2*d2t

      do j=1,Nr_max
         d(j)=0.0_dp
         d2(j)=0.0_dp
      end do

      do k=1,Nr_max
         do m=1,Nr_max
            if  ((m==1) .or. (m==Nr_max)) then
               c2=0.5_dp*(f(m)*d2t(m,k))
               d2(k)=d2(k)+c2
               c=0.5_dp*(f(m)*dt(m,k))
               d(k)=d(k)+c
               cc=0.5_dp*(f(m)*t(m,k))
            else
               c2=f(m)*d2t(m,k)
               d2(k)=d2(k)+c2
               c=f(m)*dt(m,k)
               d(k)=d(k)+c
               cc=(f(m)*t(m,k))
            end if
         end do
      end do

      d2f=((2.0_dp/real(Nr_max-1,kind=dp))**(0.5_dp)) * d2

   end subroutine chebinvtranD2

   subroutine chebinvtranD1D2(Nr_max,f,df,d2f)

      integer :: m,k,i,j
      integer, intent(in) :: Nr_max 
      complex(kind=dp), intent(in) :: f(Nr_max)
      complex(kind=dp), intent(out) :: df(Nr_max), d2f(Nr_max)
      complex(kind=dp), allocatable :: d(:), kk(:), d2(:), xp(:)
      complex(kind=dp), allocatable :: t(:,:), dt(:,:), d2t(:,:)
      complex(kind=dp) :: dm, Tm, c, cc, c2

      allocate( d(Nr_max), d2(Nr_max), kk(Nr_max), xp(Nr_max) )
      allocate( t(Nr_max,Nr_max), dt(Nr_max,Nr_max), d2t(Nr_max,Nr_max) )

      dm=1.0_dp

      do i=1,Nr_max
         kk(i)=i
      end do

      do i=1,Nr_max
         xp(i)=(cos((real(Nr_max-kk(i),kind=dp)*pi)/real(Nr_max-1,kind=dp)))
      end do
      do k=1,Nr_max
         t(1,k)=1.0_dp
         t(2,k)=xp(k)
         dt(1,k)=0.0_dp
         dt(2,k)=1.0_dp
         d2t(1,k)=0.0_dp
         d2t(2,k)=0.0_dp
         do m=3,Nr_max
            t(m,k)=cos((m-1)*real(Nr_max-k,kind=dp)*pi/real(Nr_max-1,kind=dp))
            dt(m,k)=2.0_dp*t(m-1,k) + 2.0_dp*xp(k)*dt(m-1,k) - dt(m-2,k)
            d2t(m,k) = 4.0_dp*dt(m-1,k) + 2.0_dp*xp(k)*d2t(m-1,k) - d2t(m-2,k)
         end do
      end do

      dt=2.0_dp/dm*dt
      d2t=(2.0_dp/dm)**2*d2t

      do j=1,Nr_max
         d(j)=0.0_dp
         d2(j)=0.0_dp
      end do

      do k=1,Nr_max
         do m=1,Nr_max
            if  ((m==1) .or. (m==Nr_max)) then
               c2=0.5_dp*(f(m)*d2t(m,k))
               d2(k)=d2(k)+c2
               c=0.5_dp*(f(m)*dt(m,k))
               d(k)=d(k)+c
               cc=0.5_dp*(f(m)*t(m,k))
            else
               c2=f(m)*d2t(m,k)
               d2(k)=d2(k)+c2
               c=f(m)*dt(m,k)
               d(k)=d(k)+c
               cc=(f(m)*t(m,k))
            end if
         end do
      end do

      d2f=((2.0_dp/real(Nr_max-1,kind=dp))**(0.5_dp)) * d2
      df=((2.0_dp/real(Nr_max-1,kind=dp))**(0.5_dp)) * d

   end subroutine chebinvtranD1D2

   subroutine lagrange(n,xp,x,w)

      integer, intent(in) :: n
      real(kind=dp), intent(in) :: xp(n)
      real(kind=dp), intent(in) :: x
      real(kind=dp) :: w(n)
      integer :: i,j

      w(:) = 1.0_dp

      do i=1,n 
         do j=1,n
               if (i /= j) then
                  w(i) = w(i)*(x - xp(j))/(xp(i) - xp(j))
               end if
         end do
      end do

   end subroutine lagrange

   subroutine dlagrange(n,xp,x,dw)

      integer, intent(in) :: n
      real(kind=dp), intent(in) :: xp(n)
      real(kind=dp), intent(in) :: x
      real(kind=dp), intent(out) :: dw(n)
      real(kind=dp) :: w(n)
      real(kind=dp) :: ww(n)
      integer :: i,j,m

      dw(:) = 0.0_dp
      w(:) = 1.0_dp
      ww(:) = 0.0_dp

      do j=1,n 
         do i=1,n
            w(i) = 1.0_dp
               if (i /= j) then
                  do m=1,n
                     if (m /= i .and. m /= j) then
                        w(i) = w(i)*(x - xp(m))/(xp(j) - xp(m))
                     end if
                  end do
                  ww(i) = w(i)/(xp(j)-xp(i))
                  dw(j) = dw(j) + ww(i)
               end if
         end do
      end do

   end subroutine dlagrange

   subroutine d2lagrange(n,xp,x,d2w)

      integer, intent(in) :: n
      real(kind=dp), intent(in) :: xp(n)
      real(kind=dp), intent(in) :: x
      real(kind=dp), intent(out) :: d2w(n)
      real(kind=dp) :: w(n)
      real(kind=dp) :: ww(n)
      real(kind=dp) :: mw(n)
      integer :: i,j,m,l

      d2w(:) = 0.0_dp
      w(:) = 1.0_dp
      ww(:) = 0.0_dp
      mw(:) = 0.0_dp

      do j=1,n 
         do i=1,n
               if (i /= j) then
                  mw(i) = 0.0_dp
                  do m=1,n
                     if (m /= i .and. m /= j) then
                        w(i) = 1.0_dp
                        do l=1,n
                           if ( l /= j .and. l /= i .and. l /= m ) then
                              w(i) = w(i)*(x - xp(l))/(xp(j) - xp(l))
                           end if
                        end do
                        ww(i) = w(i)/(xp(j)-xp(m))
                        mw(i) = mw(i) + ww(i)
                     end if
                  end do
                  d2w(j) = d2w(j) + mw(i)/(xp(j)-xp(i)) 
               end if
         end do
      end do

   end subroutine d2lagrange

end module chebyshev
