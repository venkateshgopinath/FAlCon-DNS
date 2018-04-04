module algebra

   use double

   implicit none

   private

   public :: factorize,matsolve,matsolve_real

contains

   !------------------------------ Matrix factorization -----------------
   subroutine factorize(Nr_max,AT,PIV1,INFO1)

      integer, intent(in) :: Nr_max
      real(kind=dp), intent(inout) :: AT(Nr_max,Nr_max)
      integer, intent(inout) :: PIV1(Nr_max)
      integer, intent(in) :: INFO1

      call dgetrf( Nr_max, Nr_max, AT, Nr_max, PIV1, INFO1) 

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

      call dgetrs( TRANS, Nr_max,int(1),AT_all, Nr_max, IPIV1, rhs_r, Nr_max, INFO1) 
      call dgetrs( TRANS, Nr_max,int(1),AT_all, Nr_max, IPIV1, rhs_i, Nr_max, INFO1) 

   end subroutine matsolve

   !------------------------------ Matrix solve Ax=b -----------------
   subroutine matsolve_real(TRANS,Nr_max,AT_all,IPIV1,rhs_r,INFO1)
                     
      character, intent(in):: TRANS
      integer, intent(in) :: Nr_max
      real(kind=dp), intent(in) :: AT_all(Nr_max,Nr_max)
      integer, intent(in) :: IPIV1(Nr_max)
      real(kind=dp), intent(inout) :: rhs_r(Nr_max)
      integer, intent(in) :: INFO1

      call dgetrs( TRANS, Nr_max,int(1),AT_all, Nr_max, IPIV1, rhs_r, Nr_max, INFO1) 

   end subroutine matsolve_real

end module algebra


