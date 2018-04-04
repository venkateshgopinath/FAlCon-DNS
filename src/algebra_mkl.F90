module algebra
 
   use lapack95
   use double

   implicit none

   private

   public :: factorize, matsolve

contains

   !------------------------------ Matrix factorization -----------------
   subroutine factorize(Nr_max,AT,PIV1,INFO1)

      integer, intent(in) :: Nr_max
      real(kind=dp), intent(inout) :: AT(Nr_max,Nr_max)
      integer, intent(inout) :: PIV1(Nr_max)
      integer, intent(out) :: INFO1

      call getrf( AT, PIV1, INFO1) 

   end subroutine factorize

   !------------------------------ Matrix solve Ax=b -----------------
   subroutine matsolve(TRANS,Nr_max,AT_all,IPIV1,rhs_r,rhs_i,INFO1)
                        
      character, intent(in):: TRANS
      integer, intent(in) :: Nr_max
      real(kind=dp), intent(in) :: AT_all(Nr_max,Nr_max)
      integer, intent(in) :: IPIV1(Nr_max)
      real(kind=dp), intent(inout) :: rhs_r(Nr_max)
      real(kind=dp), intent(inout) :: rhs_i(Nr_max)
      integer, intent(in) :: INFO1

      call getrs( AT_all, IPIV1, rhs_r) 
      call getrs( AT_all, IPIV1, rhs_i) 

   end subroutine matsolve

end module algebra
