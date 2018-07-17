module chebyshev

   use, intrinsic :: iso_c_binding
   use double
   use constants, only: pi

   implicit none

   private
   real(kind=dp), allocatable, public :: D(:,:), D2(:,:), t(:,:) 

   real(kind=dp), allocatable, public :: xp(:)

   type(C_PTR), public :: plan_cheb_fft, plan_cheb_ifft, plan_redft 


   public :: init_chebfftwplan, destroy_chebfftwplan, cheballoc, chebdealloc, dcheb, chebtransform, chebinvtran, &
        chebinvtranD1, chebinvtranD2, chebinvtranD1D2, lagrange, dlagrange, d2lagrange

contains

   subroutine init_chebfftwplan(Nr_max) ! Create plan for Discrete Fourier Transform
      include "fftw3.f03"
      integer, intent(in) :: Nr_max
      
      complex(kind=dp),  allocatable :: tf1(:)
      complex(kind=dp), allocatable :: tf2(:)
      real(kind=dp), allocatable :: fr(:)
      real(kind=dp), allocatable :: ffr(:)


      allocate( tf1(2*Nr_max-2), tf2(2*Nr_max-2), fr(Nr_max), ffr(Nr_max) )
         
      plan_cheb_fft = fftw_plan_dft_1d(2*Nr_max-2,tf1,tf2,FFTW_BACKWARD,FFTW_ESTIMATE)
      plan_cheb_ifft = fftw_plan_dft_1d(2*Nr_max-2,tf1,tf2,FFTW_FORWARD,FFTW_ESTIMATE)
      plan_redft = fftw_plan_r2r_1d(Nr_max,fr,ffr,FFTW_REDFT00,FFTW_ESTIMATE)

   end subroutine init_chebfftwplan 
   !-----------------------------------------------------------------------------------------

   subroutine destroy_chebfftwplan() ! Destroy plan
      include "fftw3.f03"
     
      call fftw_destroy_plan(plan_cheb_fft)
      call fftw_destroy_plan(plan_cheb_ifft)
      call fftw_destroy_plan(plan_redft)

   end subroutine destroy_chebfftwplan
   !-----------------------------------------------------------------------------------------

   subroutine cheballoc(Nr_max) ! Allocate Chebyshev differentiation matrix 
     
      integer, intent(in) :: Nr_max 
      
      allocate( D(Nr_max,Nr_max), D2(Nr_max,Nr_max), t(Nr_max,Nr_max), xp(Nr_max) )


   end subroutine cheballoc
   !-----------------------------------------------------------------------------------------

   subroutine chebdealloc() ! Deallocate Chebyshev differentiation matrix

      deallocate( D, D2, t, xp )

   end subroutine chebdealloc
   !-----------------------------------------------------------------------------------------

   subroutine dcheb(Nr_max) ! Construct Chebyshev differentiation matrix
    
      integer, intent(in) :: Nr_max 
      real(kind=dp) :: dt(Nr_max,Nr_max), d2t(Nr_max,Nr_max)
      integer :: dm
      integer :: i,j

      dm=1
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
   !-----------------------------------------------------------------------------------------

   subroutine chebtransform(Nr_max,f,ft)
      include "fftw3.f03"
      integer, intent(in) :: Nr_max 
      complex(kind=dp), intent(in) :: f(Nr_max)
      complex(kind=dp), intent(out) :: ft(Nr_max)
      real(kind=dp) :: fac
      complex(kind=dp), allocatable :: f0(:), f00(:), ff(:), f2(:)
      real(kind=dp), allocatable :: fnr(:), fnt(:)

      allocate( f0(Nr_max), f00(2*Nr_max-2), ff(2*Nr_max-2), f2(Nr_max), fnr(Nr_max), fnt(Nr_max) )
      !---------------------- FAST CHEBYSHEV TRANSFORM --------------------------------------
      f0(:)=f(Nr_max:1:-1) 
      ! PRE-PROCESSING
      f00(1:Nr_max)=f0(1:Nr_max)
      f00(Nr_max+1:2*Nr_max-2)=f(2:Nr_max-1)

      call fftw_execute_dft(plan_cheb_fft,f00,ff) ! Execute the plan

      ! POST-PROCESSING
      ff=ff/real(2*Nr_max-2,kind=dp)

      f2(1)=ff(1)
      f2(Nr_max)=ff(Nr_max)
      f2(2:Nr_max-1)=2.0_dp*ff(2:Nr_max-1) ! WHAT COMES OUT TO BE SAME AS THE ONE SENT FOR PRE-PROCESSING FOR INVERSE TRANSFORM
      
      ! POST-PROCESSING TO HAVE SAME TRANSFORM COEFFICIENTS AS A STANDARD CHEBYSHEV TRANSFORM WOULD DO
      fac=(2.0_dp/real(Nr_max-1,kind=dp))**0.5_dp
      ft=f2/fac
      ft(1)=2.0_dp*ft(1)
      ft(Nr_max)=2.0_dp*ft(Nr_max) ! IT IS THE ACTUAL TRANFORM NEEDED which is same as obtained by standard Chebyshev transform

   end subroutine chebtransform


   subroutine chebinvtran(Nr_max,ft,ff)
      include "fftw3.f03"
      integer :: i
      integer, intent(in) :: Nr_max
      complex(kind=dp), intent(in) :: ft(Nr_max)
      complex(kind=dp), intent(out) :: ff(Nr_max)
      complex(kind=dp) :: ftf(Nr_max), f(Nr_max)
      real(kind=dp) :: fact
      real(kind=dp) :: fr(Nr_max), fi(Nr_max), ffr(Nr_max), ffi(Nr_max), ff_r(Nr_max), ff_i(Nr_max)
      !--------------------------- INVERSE TRANSFORMS ---------------------------------

      f=ft 
      fact=(2.0_dp/real(Nr_max-1,kind=dp))**0.5_dp
      
      ftf=f

      fr=real(f) 
      fi=aimag(f)

      !------------------------ FAST INVERSE CHEBYSHEV TRANSFORM (For the function) ------------
      
      call fftw_execute_r2r(plan_redft, fr, ffr) !Execute the plan
      
      !ff_r=cmplx(ffr(1:Nr_max), kind=dp)
      ff_r=ffr(1:Nr_max)
      ff_r = ff_r * 0.5_dp * fact 
      
      call fftw_execute_r2r(plan_redft, fi, ffi) !Execute the plan
      
      !ff_i=cmplx(ffi(1:Nr_max), kind=dp)
      ff_i=ffi(1:Nr_max)
      ff_i = ff_i * 0.5_dp * fact

      do i=1,Nr_max
         ff(i)=cmplx(ff_r(i), ff_i(i),kind=dp) 
      end do

      ff(:)=ff(Nr_max:1:-1) 

   end subroutine chebinvtran

   subroutine chebinvtranD1(Nr_max,ft,df)
      include "fftw3.f03"
      integer :: i
      integer, intent(in) :: Nr_max
      
      complex(kind=dp), intent(in) :: ft(Nr_max)
      complex(kind=dp), intent(out) :: df(Nr_max)
      complex(kind=dp) :: f(Nr_max), f2c(2*Nr_max-2), f2r(2*Nr_max-2), beta1(Nr_max)
      real(kind=dp) :: fact 
      !--------------------------- INVERSE TRANSFORMS ------------------------------------------

      f=ft 

      fact=(2.0_dp/real(Nr_max-1,kind=dp))**0.5_dp

      ! Recurrence for the 1st derivative coefficients:
      beta1(Nr_max) = 0.0_dp
      beta1(Nr_max-1) = 2.0_dp * real(Nr_max-1,kind=dp) * f(Nr_max)
      do i = Nr_max-1, 2, -1
      beta1(i-1) = beta1(i+1) + 4.0_dp * real(i-1,kind=dp) * f(i)
      end do

      beta1=beta1*fact
      beta1(1)=beta1(1)/2.0_dp
      beta1(Nr_max)=beta1(Nr_max)/2.0_dp

      !- FAST INVERSE CHEBYSHEV TRANSFORM (For the 1st derivative of the function) -------------

      ! INVERSE FAST CHEBYSHEV TRANSFORM 
      
      ! PRE-PROCESSING 
      f2c(1)=beta1(1)
      f2c(2:Nr_max)=beta1(2:Nr_max)/2.0_dp
      f2c(Nr_max+1:2*Nr_max-2)=beta1(Nr_max-1:2:-1)/2.0_dp

      call fftw_execute_dft(plan_cheb_ifft,f2c,f2r)  !Execute the plan      
      
      ! POST-PROCESSING
      
      df(1:Nr_max)=f2r(1:Nr_max) 

      df(:)=df(Nr_max:1:-1) 
      !-----------------------------------------------------------------------------------------

   end subroutine chebinvtranD1

   subroutine chebinvtranD2(Nr_max,ft,d2f)
      include "fftw3.f03"
      integer :: i
      integer, intent(in) :: Nr_max
      
      complex(kind=dp), intent(in) :: ft(Nr_max)
      complex(kind=dp), intent(out) :: d2f(Nr_max)
      complex(kind=dp) :: f(Nr_max), f2c(2*Nr_max-2), f2r(2*Nr_max-2), beta1(Nr_max), beta2(Nr_max)
      real(kind=dp) :: fact 
      !--------------------------- INVERSE TRANSFORMS -------------------------------------------

      f=ft 

      fact=(2.0_dp/real(Nr_max-1,kind=dp))**0.5_dp
      !------------------------------------------------------------------------------------------


      ! Recurrence for the 1st derivative coefficients:
      beta1(Nr_max) = 0.0_dp
      beta1(Nr_max-1) = 2.0_dp * real(Nr_max-1,kind=dp) * f(Nr_max)
      do i = Nr_max-1, 2, -1
         beta1(i-1) = beta1(i+1) + 4.0_dp * real(i-1,kind=dp) * f(i)
      end do

      ! Recurrence for the 2nd derivative coefficients:
      beta2(Nr_max) = 0.0_dp
      beta2(Nr_max-1) = 2.0_dp * real(Nr_max-1,kind=dp) * beta1(Nr_max)
      do i = Nr_max-1, 2, -1
         beta2(i-1) = beta2(i+1) + 4.0_dp * real(i-1,kind=dp) * beta1(i)
      end do 

      beta1=beta1*fact
      beta1(1)=beta1(1)/2.0_dp
      beta1(Nr_max)=beta1(Nr_max)/2.0_dp

      beta2=beta2*fact
      beta2(1)=beta2(1)/2.0_dp
      beta2(Nr_max)=beta2(Nr_max)/2.0_dp


      !-- FAST INVERSE CHEBYSHEV TRANSFORM (For the 2nd derivative of the function) -------------

      ! INVERSE FAST CHEBYSHEV TRANSFORM 

      ! PRE-PROCESSING 
      f2c(1)=beta2(1)
      f2c(2:Nr_max)=beta2(2:Nr_max)/2.0_dp
      f2c(Nr_max+1:2*Nr_max-2)=beta2(Nr_max-1:2:-1)/2.0_dp

      call fftw_execute_dft(plan_cheb_ifft,f2c,f2r)  !Execute the plan     
      
      ! POST-PROCESSING
      
      d2f(1:Nr_max)=f2r(1:Nr_max) 

      d2f(:)=d2f(Nr_max:1:-1)  

      !-------------------------------------------------------------------------------------------

   end subroutine chebinvtranD2

   subroutine chebinvtranD1D2(Nr_max,ft,df,d2f)
      include "fftw3.f03"
      integer :: i
      integer, intent(in) :: Nr_max
      
      complex(kind=dp), intent(in) :: ft(Nr_max)
      complex(kind=dp), intent(out) :: d2f(Nr_max)
      complex(kind=dp), intent(out) :: df(Nr_max)
      complex(kind=dp) :: f(Nr_max), f2c(2*Nr_max-2), f2r(2*Nr_max-2), beta1(Nr_max), beta2(Nr_max)
      real(kind=dp) :: fact 
      !--------------------------- INVERSE TRANSFORMS -------------------------------------------

      f=ft 

      fact=(2.0_dp/real(Nr_max-1,kind=dp))**0.5_dp
      !------------------------------------------------------------------------------------------


      ! Recurrence for the 1st derivative coefficients:
      beta1(Nr_max) = 0.0_dp
      beta1(Nr_max-1) = 2.0_dp * real(Nr_max-1,kind=dp) * f(Nr_max)
      do i = Nr_max-1, 2, -1
         beta1(i-1) = beta1(i+1) + 4.0_dp * real(i-1,kind=dp) * f(i)
      end do

      ! Recurrence for the 2nd derivative coefficients:
      beta2(Nr_max) = 0.0_dp
      beta2(Nr_max-1) = 2.0_dp * real(Nr_max-1,kind=dp) * beta1(Nr_max)
      do i = Nr_max-1, 2, -1
         beta2(i-1) = beta2(i+1) + 4.0_dp * real(i-1,kind=dp) * beta1(i)
      end do 

      beta1=beta1*fact
      beta1(1)=beta1(1)/2.0_dp
      beta1(Nr_max)=beta1(Nr_max)/2.0_dp

      beta2=beta2*fact
      beta2(1)=beta2(1)/2.0_dp
      beta2(Nr_max)=beta2(Nr_max)/2.0_dp


      !--FAST INVERSE CHEBYSHEV TRANSFORM (For the 2nd derivative of the function) --------------

      ! INVERSE FAST CHEBYSHEV TRANSFORM 
      ! PRE-PROCESSING 

      f2c(1)=beta1(1)
      f2c(2:Nr_max)=beta1(2:Nr_max)/2.0_dp
      f2c(Nr_max+1:2*Nr_max-2)=beta1(Nr_max-1:2:-1)/2.0_dp

      call fftw_execute_dft(plan_cheb_ifft,f2c,f2r)  !Execute the plan     
      
      ! POST-PROCESSING
      
      df(1:Nr_max)=f2r(1:Nr_max) 

      df(:)=df(Nr_max:1:-1)   

      ! PRE-PROCESSING 
      f2c(1)=beta2(1)
      f2c(2:Nr_max)=beta2(2:Nr_max)/2.0_dp
      f2c(Nr_max+1:2*Nr_max-2)=beta2(Nr_max-1:2:-1)/2.0_dp

      call fftw_execute_dft(plan_cheb_ifft,f2c,f2r)  !Execute the plan     
      
      ! POST-PROCESSING
      
      d2f(1:Nr_max)=f2r(1:Nr_max) 

      d2f(:)=d2f(Nr_max:1:-1)  

      !------------------------------------------------------------------------------------------

   end subroutine chebinvtranD1D2

   subroutine lagrange(n,xp,x,w)

      integer, intent(in) :: n
      real(kind=dp), intent(in) :: xp(n)
      real(kind=dp), intent(in) :: x
      real(kind=dp), intent(out) :: w(n)
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
