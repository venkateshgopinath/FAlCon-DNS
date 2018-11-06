module mat_assembly
   !$ use omp_lib

   use double
   use constants, only: ii
   use chebyshev, only: t, D, D2
   use init, only: radius, r_radius, r_radius2, startmatbuild, finishmatbuild, time_matbuild
   use algebra, only: factorize
   use timeschemes, only: n_order_tscheme_imp

   implicit none

   private

   
   real(kind=dp), allocatable, public :: AT_all_IRK(:,:,:,:), AF_all_IRK(:,:,:,:)
   real(kind=dp), allocatable, public :: A_uphi_all_IRK(:,:,:,:)
   integer, allocatable, public :: IPIV1_IRK(:,:,:), IPIV2_IRK(:,:,:)
   integer, allocatable, public :: IPIV_uphi_IRK(:,:,:)
   real(kind=dp), allocatable, public :: LAPpsi(:,:) 
   real(kind=dp), allocatable, public :: LAPpsi_all(:,:,:) 
   integer, allocatable, public :: IPIV1_lap(:,:)

   real(kind=dp), allocatable, public :: AT_all(:,:,:), AF_all(:,:,:)
   real(kind=dp), allocatable, public :: A_uphi_all(:,:,:)
   integer, allocatable, public :: IPIV1(:,:), IPIV2(:,:)
   integer, allocatable, public :: IPIV_uphi(:,:)

   public :: allocate_mat, deallocate_mat, mat_build_IRK, mat_build_uphibar_IRK, mat_build, mat_build_uphibar, mat_build_rk

contains

   subroutine allocate_mat(Nm_max,Nr_max) !--------------------ALLOCATE MATRICES----------------------

      integer, intent(in) :: Nm_max
      integer, intent(in) :: Nr_max 

      allocate( IPIV1_IRK(Nr_max,Nm_max+1,n_order_tscheme_imp), IPIV2_IRK(2*Nr_max,Nm_max+1,n_order_tscheme_imp) )
      allocate( IPIV1_lap(Nr_max,Nm_max+1) )
      allocate( IPIV_uphi_IRK(Nr_max,Nm_max+1,n_order_tscheme_imp) )   
      allocate( A_uphi_all_IRK(Nr_max,Nr_max,Nm_max+1,n_order_tscheme_imp) )
      allocate( AT_all_IRK(Nr_max,Nr_max,Nm_max+1,n_order_tscheme_imp), AF_all_IRK(2*Nr_max,2*Nr_max,Nm_max+1,n_order_tscheme_imp) )
      allocate( LAPpsi(Nr_max,Nr_max) )
      allocate( LAPpsi_all(Nr_max,Nr_max,Nm_max+1) )
      allocate( IPIV1(Nr_max,Nm_max+1), IPIV2(2*Nr_max,Nm_max+1) )
      allocate( IPIV_uphi(Nr_max,Nm_max+1) ) 
      allocate( A_uphi_all(Nr_max,Nr_max,Nm_max+1) )
      allocate( AT_all(Nr_max,Nr_max,Nm_max+1), AF_all(2*Nr_max,2*Nr_max,Nm_max+1) )

   end subroutine allocate_mat

   subroutine deallocate_mat() !----------------DEALLOCATE MATRICES----------------------

      deallocate( AT_all_IRK, AF_all_IRK )
      deallocate( A_uphi_all_IRK )
      deallocate( IPIV_uphi_IRK )   
      deallocate( IPIV1_IRK, IPIV2_IRK )
      deallocate( IPIV1_lap )
      deallocate( LAPpsi )
      deallocate( LAPpsi_all )
      deallocate( AT_all, AF_all )
      deallocate( A_uphi_all )
      deallocate( IPIV_uphi )   
      deallocate( IPIV1, IPIV2 )

   end subroutine deallocate_mat

   subroutine mat_build_IRK(Nr_max,dt,n,mBC,wt_lhs_tscheme_imp,rk_stage,Pr)

      integer, intent(in) :: Nr_max 
      integer, intent(in) :: n
      character(len=100), intent(in) :: mBC 
      real(kind=dp), intent(in) :: dt
      real(kind=dp), intent(in) :: Pr
      real(kind=dp) :: C
      integer :: Nr_max2
      integer :: i,j 
      integer :: INFO1, INFO2
      real(kind=dp), intent(in) :: wt_lhs_tscheme_imp
      integer, intent(in) :: rk_stage
       
      real(kind=dp) :: AT(Nr_max,Nr_max)
      real(kind=dp) :: AF(2*Nr_max,2*Nr_max)
      real(kind=dp) :: bctemp(2,Nr_max), bcpsi1(2,Nr_max), bcpsi2(2,Nr_max)
      integer :: PIV1(Nr_max), PIV2(2*Nr_max)
      
      AT(:,:)=0.0_dp
      AF(:,:)=0.0_dp
      PIV1(:)=0
      PIV2(:)=0
                           ! OPERATOR MATRIX ASSEMBLY and LU Factorization 
      !**************** AT is the operator matrix for temperature equation *************************************
      C=((2.0_dp/real(Nr_max-1,kind=dp))**(0.5_dp))
      Nr_max2=2*Nr_max

      !-----Temperature BC-----------------------------------------------------
      bctemp(1,:) = t(1,:)      ! bottom BC
      bctemp(2,:) = t(Nr_max,:) ! top BC
      !------------------------------------------------------------------------

      !-----Vorticity and streamfunction BC------------------------------------
      bcpsi1(1,:) = t(1,:)      ! 1st BC is psi=0 (bottom BC) 
      bcpsi1(2,:) = t(Nr_max,:) ! 1st BC is psi=0 (top BC)
 
      if (mBC=='NS') then ! Enforce no-slip boundary condition
               bcpsi2(1,:)=D(1,:)        ! 2nd BC for psi (bottom BC)
               bcpsi2(2,:)=D(Nr_max,:)   ! 2nd BC for psi (top BC)
      elseif (mBC=='SF') then ! Enforce stress-free boundary condition
               bcpsi2(1,:)=D2(1,:)-r_radius(1)*D(1,:)                ! 2nd BC for psi (bottom BC)
               bcpsi2(2,:)=D2(Nr_max,:)-r_radius(Nr_max)*D(Nr_max,:) ! 2nd BC for psi (top BC)
      end if
      !------------------------------------------------------------------------

      do j=1,Nr_max
         do i=1,Nr_max
            AT(i,j)=t(i,j)-wt_lhs_tscheme_imp*dt*(1.0_dp/Pr)*(r_radius(i)*D(i,j)+D2(i,j)-real(n,kind=dp)* & 
                    & real(n,kind=dp)*r_radius2(i)*t(i,j))
         end do
            AT(1,j)=bctemp(1,j)
            AT(Nr_max,j)=bctemp(2,j)
      end do

      do i=1,Nr_max
         AT(i,1)=0.5_dp*AT(i,1)
         AT(i,Nr_max)=0.5_dp*AT(i,Nr_max)
      end do 
                            
      AT=C*AT

      !****** AF is coupled operator matrix for psi and omega *****************************************************
      do j=1,Nr_max
         do i=1,Nr_max
            AF(i,j)=r_radius(i)*D(i,j)+D2(i,j)-real(n,kind=dp)*real(n,kind=dp)*r_radius2(i)*t(i,j)

            AF(i,Nr_max+j)=t(i,j)

            AF(Nr_max+i,j)=0.0_dp

            AF(Nr_max+i,Nr_max+j)=t(i,j)-wt_lhs_tscheme_imp*dt*(r_radius(i)*D(i,j)+D2(i,j)-real(n,kind=dp)*real(n,kind=dp)* &
                          & r_radius2(i)*t(i,j))

         end do
         AF(1,j)=bcpsi1(1,j)
         AF(Nr_max,j)=bcpsi1(2,j)

         AF(1,Nr_max+j)=0.0_dp
         AF(Nr_max,Nr_max+j)=0.0_dp

         AF(Nr_max+1,j)=bcpsi2(1,j)
         AF(Nr_max+Nr_max,j)=bcpsi2(2,j)

         AF(Nr_max+1,Nr_max+j)=0.0_dp
         AF(Nr_max+Nr_max,Nr_max+j)=0.0_dp

      end do

      do i=1,Nr_max
         AF(i,1)=0.5_dp*AF(i,1)
         AF(i,Nr_max)=0.5_dp*AF(i,Nr_max)

         AF(i,Nr_max+1)=0.5_dp*AF(i,Nr_max+1)
         AF(i,Nr_max+Nr_max)=0.5_dp*AF(i,Nr_max+Nr_max)

         AF(Nr_max+i,1)=0.5_dp*AF(Nr_max+i,1)
         AF(Nr_max+i,Nr_max)=0.5_dp*AF(Nr_max+i,Nr_max)

         AF(Nr_max+i,Nr_max+1)=0.5_dp*AF(Nr_max+i,Nr_max+1)
         AF(Nr_max+i,Nr_max+Nr_max)=0.5_dp*AF(Nr_max+i,Nr_max+Nr_max)

      end do

      AF=C*AF

      !***** CALL DGETRF factorization for AT matrix
      call factorize(Nr_max,AT,PIV1,INFO1)
                           
      AT_all_IRK(:,:,n+1,rk_stage)=AT
      IPIV1_IRK(:,n+1,rk_stage)=PIV1
           
      !***** CALL DGETRF factorization for AF matrix
      call factorize(2*Nr_max,AF,PIV2,INFO2)

      AF_all_IRK(:,:,n+1,rk_stage)=AF
      IPIV2_IRK(:,n+1,rk_stage)=PIV2
                                    
   end subroutine mat_build_IRK

   subroutine mat_build_uphibar_IRK(Nr_max,dt,mBC,wt_lhs_tscheme_imp,rk_stage,Pr)

      integer, intent(in) :: Nr_max 
      character(len=100), intent(in) :: mBC 
      real(kind=dp), intent(in) :: dt
      real(kind=dp), intent(in) :: Pr
      real(kind=dp) :: C
      integer :: i,j 
      integer :: INFO3
      real(kind=dp), intent(in) :: wt_lhs_tscheme_imp
      real(kind=dp) :: A_uphi(Nr_max,Nr_max)
      real(kind=dp) :: bcuphibar(2,Nr_max)
      integer :: PIV_uphi(Nr_max)
      integer, intent(in) :: rk_stage
      
      A_uphi(:,:)=0.0_dp
      PIV_uphi(:)=0

      !print *, wt_lhs_tscheme_imp, "wt imp"
                           ! OPERATOR MATRIX ASSEMBLY and LU Factorization 
      !**************** AT is the operator matrix for temperature equation *************************************
      C=((2.0_dp/real(Nr_max-1,kind=dp))**(0.5_dp))

      if (mBC=='NS') then ! Enforce no-slip boundary condition
               bcuphibar(1,:)=t(1,:)        ! BC for uphibar (bottom BC)
               bcuphibar(2,:)=t(Nr_max,:)   ! BC for uphibar (top BC)
      elseif (mBC=='SF') then ! Enforce stress-free boundary condition
               bcuphibar(1,:)=D(1,:)-r_radius(1)*t(1,:)     ! BC for uphibar (bottom BC)
               bcuphibar(2,:)=D(Nr_max,:)-r_radius(Nr_max)*t(Nr_max,:)   ! BC for uphibar (top BC)
      end if
      !------------------------------------------------------------------------

      do j=1,Nr_max
         do i=1,Nr_max
            A_uphi(i,j)=t(i,j)-wt_lhs_tscheme_imp*dt*(r_radius(i)*D(i,j)+D2(i,j) - r_radius2(i)*t(i,j))
         end do
            A_uphi(1,j)=bcuphibar(1,j)
            A_uphi(Nr_max,j)=bcuphibar(2,j)
      end do

      do i=1,Nr_max
         A_uphi(i,1)=0.5_dp*A_uphi(i,1)
         A_uphi(i,Nr_max)=0.5_dp*A_uphi(i,Nr_max)
      end do 
                            
      A_uphi=C*A_uphi


      !***** CALL DGETRF factorization for A_uphi matrix
      call factorize(Nr_max,A_uphi,PIV_uphi,INFO3)
      A_uphi_all_IRK(:,:,1,rk_stage)=A_uphi
      IPIV_uphi_IRK(:,1,rk_stage)=PIV_uphi
                                    
   end subroutine mat_build_uphibar_IRK

   subroutine mat_build(Nr_max,dt,n,mBC,wt_lhs_tscheme_imp,Pr)

      integer, intent(in) :: Nr_max 
      integer, intent(in) :: n
      character(len=100), intent(in) :: mBC 
      real(kind=dp), intent(in) :: dt
      real(kind=dp), intent(in) :: Pr
      real(kind=dp) :: C
      integer :: Nr_max2
      integer :: i,j 
      integer :: INFO1, INFO2
      real(kind=dp), intent(in) :: wt_lhs_tscheme_imp
       
      real(kind=dp) :: AT(Nr_max,Nr_max)
      real(kind=dp) :: AF(2*Nr_max,2*Nr_max)
      real(kind=dp) :: bctemp(2,Nr_max), bcpsi1(2,Nr_max), bcpsi2(2,Nr_max)
      integer :: PIV1(Nr_max), PIV2(2*Nr_max)
      
      AT(:,:)=0.0_dp
      AF(:,:)=0.0_dp
      PIV1(:)=0
      PIV2(:)=0
                           ! OPERATOR MATRIX ASSEMBLY and LU Factorization 
      !**************** AT is the operator matrix for temperature equation *************************************
      C=((2.0_dp/real(Nr_max-1,kind=dp))**(0.5_dp))
      Nr_max2=2*Nr_max

      !-----Temperature BC-----------------------------------------------------
      bctemp(1,:) = t(1,:)      ! bottom BC
      bctemp(2,:) = t(Nr_max,:) ! top BC
      !------------------------------------------------------------------------

      !-----Vorticity and streamfunction BC------------------------------------
      bcpsi1(1,:) = t(1,:)      ! 1st BC is psi=0 (bottom BC) 
      bcpsi1(2,:) = t(Nr_max,:) ! 1st BC is psi=0 (top BC)
 
      if (mBC=='NS') then ! Enforce no-slip boundary condition
               bcpsi2(1,:)=D(1,:)        ! 2nd BC for psi (bottom BC)
               bcpsi2(2,:)=D(Nr_max,:)   ! 2nd BC for psi (top BC)
      elseif (mBC=='SF') then ! Enforce stress-free boundary condition
               bcpsi2(1,:)=D2(1,:)-r_radius(1)*D(1,:)                ! 2nd BC for psi (bottom BC)
               bcpsi2(2,:)=D2(Nr_max,:)-r_radius(Nr_max)*D(Nr_max,:) ! 2nd BC for psi (top BC)
      end if
      !------------------------------------------------------------------------

      do j=1,Nr_max
         do i=1,Nr_max
            AT(i,j)=t(i,j)-wt_lhs_tscheme_imp*dt*(1.0_dp/Pr)*(r_radius(i)*D(i,j)+D2(i,j)-real(n,kind=dp)* & 
                    & real(n,kind=dp)*r_radius2(i)*t(i,j))
         end do
            AT(1,j)=bctemp(1,j)
            AT(Nr_max,j)=bctemp(2,j)
      end do

      do i=1,Nr_max
         AT(i,1)=0.5_dp*AT(i,1)
         AT(i,Nr_max)=0.5_dp*AT(i,Nr_max)
      end do 
                            
      AT=C*AT

      !****** AF is coupled operator matrix for psi and omega *****************************************************
      do j=1,Nr_max
         do i=1,Nr_max
            AF(i,j)=r_radius(i)*D(i,j)+D2(i,j)-real(n,kind=dp)*real(n,kind=dp)*r_radius2(i)*t(i,j)

            AF(i,Nr_max+j)=t(i,j)

            AF(Nr_max+i,j)=0.0_dp

            AF(Nr_max+i,Nr_max+j)=t(i,j)-wt_lhs_tscheme_imp*dt*(r_radius(i)*D(i,j)+D2(i,j)-real(n,kind=dp)*real(n,kind=dp)* &
                          & r_radius2(i)*t(i,j))

         end do
         AF(1,j)=bcpsi1(1,j)
         AF(Nr_max,j)=bcpsi1(2,j)

         AF(1,Nr_max+j)=0.0_dp
         AF(Nr_max,Nr_max+j)=0.0_dp

         AF(Nr_max+1,j)=bcpsi2(1,j)
         AF(Nr_max+Nr_max,j)=bcpsi2(2,j)

         AF(Nr_max+1,Nr_max+j)=0.0_dp
         AF(Nr_max+Nr_max,Nr_max+j)=0.0_dp

      end do

      do i=1,Nr_max
         AF(i,1)=0.5_dp*AF(i,1)
         AF(i,Nr_max)=0.5_dp*AF(i,Nr_max)

         AF(i,Nr_max+1)=0.5_dp*AF(i,Nr_max+1)
         AF(i,Nr_max+Nr_max)=0.5_dp*AF(i,Nr_max+Nr_max)

         AF(Nr_max+i,1)=0.5_dp*AF(Nr_max+i,1)
         AF(Nr_max+i,Nr_max)=0.5_dp*AF(Nr_max+i,Nr_max)

         AF(Nr_max+i,Nr_max+1)=0.5_dp*AF(Nr_max+i,Nr_max+1)
         AF(Nr_max+i,Nr_max+Nr_max)=0.5_dp*AF(Nr_max+i,Nr_max+Nr_max)

      end do

      AF=C*AF

      !***** CALL DGETRF factorization for AT matrix
      call factorize(Nr_max,AT,PIV1,INFO1)
                           
      AT_all(:,:,n+1)=AT
      IPIV1(:,n+1)=PIV1
           
      !***** CALL DGETRF factorization for AF matrix
      call factorize(2*Nr_max,AF,PIV2,INFO2)

      AF_all(:,:,n+1)=AF
      IPIV2(:,n+1)=PIV2
                                    
   end subroutine mat_build

   subroutine mat_build_uphibar(Nr_max,dt,mBC,wt_lhs_tscheme_imp,Pr)

      integer, intent(in) :: Nr_max 
      character(len=100), intent(in) :: mBC 
      real(kind=dp), intent(in) :: dt
      real(kind=dp), intent(in) :: Pr
      real(kind=dp) :: C
      integer :: i,j 
      integer :: INFO3
      real(kind=dp), intent(in) :: wt_lhs_tscheme_imp
      real(kind=dp) :: A_uphi(Nr_max,Nr_max)
      real(kind=dp) :: bcuphibar(2,Nr_max)
      integer :: PIV_uphi(Nr_max)
      
      A_uphi(:,:)=0.0_dp
      PIV_uphi(:)=0
      A_uphi_all(:,:,1)=0.0_dp
      IPIV_uphi(:,1)=0

      !print *, wt_lhs_tscheme_imp, "wt imp"
                           ! OPERATOR MATRIX ASSEMBLY and LU Factorization 
      !**************** AT is the operator matrix for temperature equation *************************************
      C=((2.0_dp/real(Nr_max-1,kind=dp))**(0.5_dp))

      if (mBC=='NS') then ! Enforce no-slip boundary condition
               bcuphibar(1,:)=t(1,:)        ! BC for uphibar (bottom BC)
               bcuphibar(2,:)=t(Nr_max,:)   ! BC for uphibar (top BC)
      elseif (mBC=='SF') then ! Enforce stress-free boundary condition
               bcuphibar(1,:)=D(1,:)-r_radius(1)*t(1,:)        ! BC for uphibar (bottom BC)
               bcuphibar(2,:)=D(Nr_max,:)-r_radius(Nr_max)*t(Nr_max,:)   ! BC for uphibar (top BC)
      end if
      !------------------------------------------------------------------------

      do j=1,Nr_max
         do i=1,Nr_max
            A_uphi(i,j)=t(i,j)-wt_lhs_tscheme_imp*dt*(r_radius(i)*D(i,j)+D2(i,j) - r_radius2(i)*t(i,j))
         end do
            A_uphi(1,j)=bcuphibar(1,j)
            A_uphi(Nr_max,j)=bcuphibar(2,j)
      end do

      do i=1,Nr_max
         A_uphi(i,1)=0.5_dp*A_uphi(i,1)
         A_uphi(i,Nr_max)=0.5_dp*A_uphi(i,Nr_max)
      end do 
                            
      A_uphi=C*A_uphi


      !***** CALL DGETRF factorization for A_uphi matrix
      call factorize(Nr_max,A_uphi,PIV_uphi,INFO3)
      A_uphi_all(:,:,1)=A_uphi
      IPIV_uphi(:,1)=PIV_uphi
                                    
   end subroutine mat_build_uphibar

   subroutine mat_build_rk(Nr_max,n,mBC)

      integer, intent(in) :: Nr_max 
      integer, intent(in) :: n
      character(len=100), intent(in) :: mBC 
      real(kind=dp) :: C
      integer :: i,j 
      integer :: INFO1
      integer :: PIV1(Nr_max)

                           ! OPERATOR MATRIX ASSEMBLY and LU Factorization 
      !**************** AT is the operator matrix for temperature equation *************************************
      C=((2.0_dp/real(Nr_max-1,kind=dp))**(0.5_dp))
      
      do j=1,Nr_max
         do i=1,Nr_max
            LAPpsi(i,j)=r_radius(i)*D(i,j)+D2(i,j)-real(n,kind=dp)*real(n,kind=dp)*r_radius2(i)*t(i,j) 
         end do
            !if (mBC=='NS') then
               LAPpsi(1,j)=1.0_dp*t(1,j)
               LAPpsi(Nr_max,j)=1.0_dp*t(Nr_max,j)
            !   ! Attempt at solving with 4 BCs instead of Johnstion's fix for assembly type IMEXRK cases
            !   LAPpsi(2,j)=1.0_dp*D(1,j)
            !   LAPpsi(Nr_max-1,j)=1.0_dp*D(Nr_max,j)
            !   ! ----------------------------------
            !else if (mBC=='SF')  then
            !   LAPpsi(1,j)=1.0_dp*t(1,j)
            !   LAPpsi(Nr_max,j)=1.0_dp*t(Nr_max,j)
            !   ! Attempt at solving with 4 BCs instead of Johnstion's fix for assembly type IMEXRK cases
            !   LAPpsi(2,j)=D2(1,j)-r_radius(1)*D(1,j)
            !   LAPpsi(Nr_max-1,j)=D2(Nr_max,j)-r_radius(Nr_max)*D(Nr_max,j)
            !end if
      end do

      do i=1,Nr_max
         LAPpsi(i,1)=0.5_dp*LAPpsi(i,1)
         LAPpsi(i,Nr_max)=0.5_dp*LAPpsi(i,Nr_max)
      end do
                   
      LAPpsi=1.0_dp*C*LAPpsi
     
      !***** CALL DGETRF factorization for AT matrix
      call factorize(Nr_max,LAPpsi,PIV1,INFO1)
                           
      LAPpsi_all(:,:,n+1)=LAPpsi
      IPIV1_lap(:,n+1)=PIV1
        
   end subroutine mat_build_rk

end module mat_assembly
 
