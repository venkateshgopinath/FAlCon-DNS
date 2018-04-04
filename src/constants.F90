module constants

   use double

   implicit none

   real(kind=dp), parameter :: pi=4.0_dp*atan(1.0_dp)
   complex(kind=dp), parameter :: ii=(0.0_dp, 1.0_dp)
   real(kind=dp), parameter :: zero=0.0_dp
   real(kind=dp), parameter :: one=1.0_dp
   real(kind=dp), parameter :: two=2.0_dp
   real(kind=dp), parameter :: three=3.0_dp
   real(kind=dp), parameter :: half=0.5_dp
   real(kind=dp), parameter :: onesixth = 1.0_dp/6.0_dp
   real(kind=dp), parameter :: onethird = 1.0_dp/3.0_dp
   real(kind=dp), parameter :: onefourth = 1.0_dp/4.0_dp
   real(kind=dp), parameter :: threefourth = 3.0_dp/4.0_dp
   real(kind=dp), parameter :: sevenfourth = 7.0_dp/4.0_dp
   real(kind=dp), parameter :: fivesixth = 5.0_dp/6.0_dp
   real(kind=dp), parameter :: oneninth = 1.0_dp/9.0_dp
   real(kind=dp), parameter :: twoninth = 2.0_dp/9.0_dp
   real(kind=dp), parameter :: fourninth = 4.0_dp/9.0_dp
   real(kind=dp), parameter :: threehalf = 3.0_dp/2.0_dp
   real(kind=dp), parameter :: oneeighteenth = 1.0_dp/18.0_dp
   real(kind=dp), parameter :: eleveneighteenth = 11.0_dp/18.0_dp
   real(kind=dp), parameter :: fiveeighteenth = 5.0_dp/18.0_dp

end module constants
