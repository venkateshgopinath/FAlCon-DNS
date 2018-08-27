module timeschemes

   use double
   use constants, only: one, onesixth, onethird, onefourth, threefourth, threehalf, &
                        & half, eleveneighteenth, oneeighteenth, sevenfourth, fivesixth, &
                        & fiveeighteenth, oneninth, twoninth, fourninth, twothird
   use namelists, only: dt
   implicit none

   private

   real(kind=dp), public :: wt_lhs_tscheme_imp
   real(kind=dp), allocatable, public :: wt_rhs_tscheme_imp(:), wt_rhs_tscheme_exp(:)
   integer, public :: n_order_tscheme_imp, n_order_tscheme_exp, n_order_tscheme_max, irk_max
   complex(kind=dp), allocatable, public :: rhs_imp_temp(:,:,:), rhs_exp_temp(:,:,:) 
   complex(kind=dp), allocatable, public :: rhs_imp_vort(:,:,:), rhs_exp_vort(:,:,:)
   complex(kind=dp), allocatable, public :: rhs_imp_uphi_bar(:,:), rhs_exp_uphi_bar(:,:)
   complex(kind=dp), allocatable, public :: rhs_buo_term(:,:,:)
   real(kind=dp), allocatable, public :: dt_array(:)
   real(kind=dp), allocatable, public :: butcher_a(:,:)
   real(kind=dp), allocatable, public :: butcher_b(:)
   real(kind=dp), allocatable, public :: butcher_aD(:,:),butcher_aA(:,:) 
   real(kind=dp), allocatable, public :: butcher_bD(:),butcher_bA(:)
   logical, public :: ars_eqn_check_A, ars_eqn_check_D, diag_diff
   integer, public :: diag_index
   real(kind=dp) :: gammaa, delta, alphaa, betaa, etaa, lambda
   real(kind=dp) :: a31, a32, a41, a42, a43, a44, a53, a54, a55, bb1, bb2, bb3, bb4, bb5, bb6, c3, c4 
   real(kind=dp) :: aa31, aa32, aa41, aa42, aa43, aa51, aa52, aa53, aa54, aa64, aa65 
   real(kind=dp) :: a21, a22, a33, c2, aa21, c1 

   public :: init_time_schemes, deallocate_time_schemes, rhs_update_wts_imp, &
             & rhs_update_wts_exp

contains

   subroutine init_time_schemes(Nm_max,Nr_max,time_scheme_imp,time_scheme_exp,time_scheme_type)

      integer, intent(in) :: Nm_max,Nr_max
      character(len=100), intent(in) :: time_scheme_imp, time_scheme_exp, time_scheme_type

      select case (time_scheme_imp)

         case ('CN')
            n_order_tscheme_imp=2  
            allocate( wt_rhs_tscheme_imp(n_order_tscheme_imp) )
            allocate( rhs_imp_temp(n_order_tscheme_imp,Nm_max+1,Nr_max) )
            allocate( rhs_imp_vort(n_order_tscheme_imp,Nm_max+1,Nr_max) )
            allocate( rhs_imp_uphi_bar(n_order_tscheme_imp,Nr_max) )
            allocate( rhs_buo_term(n_order_tscheme_imp,Nm_max+1,Nr_max) )
            rhs_imp_temp(2,:,:)=0.0_dp
            rhs_imp_vort(2,:,:)=0.0_dp
            rhs_imp_uphi_bar(2,:)=0.0_dp

            select case (time_scheme_exp)
               case ('AB2')
                  n_order_tscheme_exp=2  
                  allocate( wt_rhs_tscheme_exp(n_order_tscheme_exp) )
                  allocate( rhs_exp_temp(n_order_tscheme_exp,Nm_max+1,Nr_max) )
                  allocate( rhs_exp_vort(n_order_tscheme_exp,Nm_max+1,Nr_max) )
                  allocate( rhs_exp_uphi_bar(n_order_tscheme_exp,Nr_max) )
                  rhs_exp_temp(2,:,:)=0.0_dp 
                  rhs_exp_vort(2,:,:)=0.0_dp
                  rhs_exp_uphi_bar(2,:)=0.0_dp 
            
               case ('AB3')
                  n_order_tscheme_exp=3  
                  allocate( wt_rhs_tscheme_exp(n_order_tscheme_exp) )
                  allocate( rhs_exp_temp(n_order_tscheme_exp,Nm_max+1,Nr_max) )
                  allocate( rhs_exp_vort(n_order_tscheme_exp,Nm_max+1,Nr_max) )
                  allocate( rhs_exp_uphi_bar(n_order_tscheme_exp,Nr_max) )
                  rhs_exp_temp(2:3,:,:)=0.0_dp 
                  rhs_exp_vort(2:3,:,:)=0.0_dp
                  rhs_exp_uphi_bar(2:3,:)=0.0_dp 
               
            end select

         case ('BDF2') 
            n_order_tscheme_imp=2  
            allocate( wt_rhs_tscheme_imp(n_order_tscheme_imp) )
            allocate( rhs_imp_temp(n_order_tscheme_imp,Nm_max+1,Nr_max) )
            allocate( rhs_imp_vort(n_order_tscheme_imp,Nm_max+1,Nr_max) )
            allocate( rhs_imp_uphi_bar(n_order_tscheme_imp,Nr_max) )
            allocate( rhs_buo_term(n_order_tscheme_imp,Nm_max+1,Nr_max) )
            rhs_imp_temp(2,:,:)=0.0_dp
            rhs_imp_vort(2,:,:)=0.0_dp
            rhs_imp_uphi_bar(2,:)=0.0_dp
            rhs_buo_term(2,:,:)=0.0_dp

            select case (time_scheme_exp)
               case ('AB2')
                  n_order_tscheme_imp=2  
                  n_order_tscheme_exp=2  
                  allocate( wt_rhs_tscheme_exp(n_order_tscheme_exp) )
                  allocate( rhs_exp_temp(n_order_tscheme_exp,Nm_max+1,Nr_max) )
                  allocate( rhs_exp_vort(n_order_tscheme_exp,Nm_max+1,Nr_max) )
                  allocate( rhs_exp_uphi_bar(n_order_tscheme_exp,Nr_max) )
                  rhs_exp_temp(2,:,:)=0.0_dp 
                  rhs_exp_vort(2,:,:)=0.0_dp 
                  rhs_exp_uphi_bar(2,:)=0.0_dp 
               
            end select

         case ('BDF3') 
            n_order_tscheme_imp=3  
            allocate( wt_rhs_tscheme_imp(n_order_tscheme_imp) )
            allocate( rhs_imp_temp(n_order_tscheme_imp,Nm_max+1,Nr_max) )
            allocate( rhs_imp_vort(n_order_tscheme_imp,Nm_max+1,Nr_max) )
            allocate( rhs_imp_uphi_bar(n_order_tscheme_imp,Nr_max) )
            allocate( rhs_buo_term(n_order_tscheme_imp,Nm_max+1,Nr_max) )
            rhs_imp_temp(2:3,:,:)=0.0_dp 
            rhs_imp_vort(2:3,:,:)=0.0_dp 
            rhs_imp_uphi_bar(2:3,:)=0.0_dp
            rhs_buo_term(2:3,:,:)=0.0_dp

            select case (time_scheme_exp)
               case('AB3')
                  n_order_tscheme_imp=3  
                  n_order_tscheme_exp=3  
                  allocate( wt_rhs_tscheme_exp(n_order_tscheme_exp) )
                  allocate( rhs_exp_temp(n_order_tscheme_exp,Nm_max+1,Nr_max) )
                  allocate( rhs_exp_vort(n_order_tscheme_exp,Nm_max+1,Nr_max) )
                  allocate( rhs_exp_uphi_bar(n_order_tscheme_exp,Nr_max) )
                  rhs_exp_temp(2:3,:,:)=0.0_dp 
                  rhs_exp_vort(2:3,:,:)=0.0_dp 
                  rhs_exp_uphi_bar(2:3,:)=0.0_dp 

            end select

         case ('BDF4') 
            n_order_tscheme_imp=4  
            allocate( wt_rhs_tscheme_imp(n_order_tscheme_imp) )
            allocate( rhs_imp_temp(n_order_tscheme_imp,Nm_max+1,Nr_max) )
            allocate( rhs_imp_vort(n_order_tscheme_imp,Nm_max+1,Nr_max) )
            allocate( rhs_imp_uphi_bar(n_order_tscheme_imp,Nr_max) )
            allocate( rhs_buo_term(n_order_tscheme_imp,Nm_max+1,Nr_max) )
            rhs_imp_temp(2:4,:,:)=0.0_dp 
            rhs_imp_vort(2:4,:,:)=0.0_dp 
            rhs_imp_uphi_bar(2:4,:)=0.0_dp
            rhs_buo_term(2:4,:,:)=0.0_dp

            select case (time_scheme_exp)
               case('AB4')
                  n_order_tscheme_imp=4  
                  n_order_tscheme_exp=4  
                  allocate( wt_rhs_tscheme_exp(n_order_tscheme_exp) )
                  allocate( rhs_exp_temp(n_order_tscheme_exp,Nm_max+1,Nr_max) )
                  allocate( rhs_exp_vort(n_order_tscheme_exp,Nm_max+1,Nr_max) )
                  allocate( rhs_exp_uphi_bar(n_order_tscheme_exp,Nr_max) )
                  rhs_exp_temp(2:4,:,:)=0.0_dp 
                  rhs_exp_vort(2:4,:,:)=0.0_dp 
                  rhs_exp_uphi_bar(2:4,:)=0.0_dp 

            end select

         case ('NONE')

            select case (time_scheme_exp)

              case ('RK2')
                  n_order_tscheme_imp=2 
                  allocate( rhs_imp_temp(1,Nm_max+1,Nr_max) )
                  allocate( rhs_imp_vort(1,Nm_max+1,Nr_max) )
                  allocate( rhs_imp_uphi_bar(1,Nr_max) )
                  n_order_tscheme_exp=2
                  allocate( butcher_a(n_order_tscheme_exp,n_order_tscheme_exp) )
                  allocate( butcher_b(n_order_tscheme_exp) )
                  allocate( rhs_exp_temp(n_order_tscheme_exp,Nm_max+1,Nr_max) )
                  allocate( rhs_exp_vort(n_order_tscheme_exp,Nm_max+1,Nr_max) )
                  allocate( rhs_exp_uphi_bar(n_order_tscheme_exp,Nr_max) )
                  allocate( wt_rhs_tscheme_imp(1) ) 
                  allocate( wt_rhs_tscheme_exp(1) )  

              case ('RK4')
                  n_order_tscheme_imp=4 
                  allocate( rhs_imp_temp(1,Nm_max+1,Nr_max) )
                  allocate( rhs_imp_vort(1,Nm_max+1,Nr_max) )
                  allocate( rhs_imp_uphi_bar(1,Nr_max) )
                  n_order_tscheme_exp=4
                  allocate( butcher_a(n_order_tscheme_exp,n_order_tscheme_exp) )
                  allocate( butcher_b(n_order_tscheme_exp) )
                  allocate( rhs_exp_temp(n_order_tscheme_exp,Nm_max+1,Nr_max) )
                  allocate( rhs_exp_vort(n_order_tscheme_exp,Nm_max+1,Nr_max) )
                  allocate( rhs_exp_uphi_bar(n_order_tscheme_exp,Nr_max) )
                  allocate( wt_rhs_tscheme_imp(1) ) 
                  allocate( wt_rhs_tscheme_exp(1) )  

            end select 
         ! FOR IMEX RK SCHEMES below, n_order_tscheme_imp or exp = number of stages -------
         case ('EUL111')
            n_order_tscheme_imp=2 
            allocate( butcher_aD(n_order_tscheme_imp,n_order_tscheme_imp) )
            allocate( butcher_bD(n_order_tscheme_imp) )
            allocate( rhs_imp_temp(1,Nm_max+1,Nr_max) )
            allocate( rhs_imp_vort(1,Nm_max+1,Nr_max) )
            allocate( rhs_imp_uphi_bar(1,Nr_max) )
            allocate( rhs_buo_term(n_order_tscheme_imp,Nm_max+1,Nr_max) )
            allocate( wt_rhs_tscheme_imp(1) ) 

            select case (time_scheme_exp)

              case ('expEUL111')
                  n_order_tscheme_exp=2
                  allocate( butcher_aA(n_order_tscheme_exp,n_order_tscheme_exp) )
                  allocate( butcher_bA(n_order_tscheme_exp) )
                  allocate( rhs_exp_temp(1,Nm_max+1,Nr_max) )
                  allocate( rhs_exp_vort(1,Nm_max+1,Nr_max) )
                  allocate( rhs_exp_uphi_bar(1,Nr_max) )
                  allocate( wt_rhs_tscheme_exp(1) )  

            end select   

         case ('ARS222') ! Ascher, Ruuth and Spiteri (ARS222) 2nd order, Appl. Numer. Math. 1997 paper
            n_order_tscheme_imp=3 
            allocate( butcher_aD(n_order_tscheme_imp,n_order_tscheme_imp) )
            allocate( butcher_bD(n_order_tscheme_imp) )
            allocate( rhs_imp_temp(1,Nm_max+1,Nr_max) )
            allocate( rhs_imp_vort(1,Nm_max+1,Nr_max) )
            allocate( rhs_imp_uphi_bar(1,Nr_max) )
            allocate( wt_rhs_tscheme_imp(1) ) 
            gammaa = (2.0_dp-sqrt(2.0_dp))/2.0_dp
            delta = 1.0_dp-1.0_dp/gammaa/2.0_dp
            select case (time_scheme_exp)

              case ('expARS222')
                  n_order_tscheme_exp=3
                  allocate( butcher_aA(n_order_tscheme_exp,n_order_tscheme_exp) )
                  allocate( butcher_bA(n_order_tscheme_exp) )
                  allocate( rhs_exp_temp(1,Nm_max+1,Nr_max) )
                  allocate( rhs_exp_vort(1,Nm_max+1,Nr_max) )
                  allocate( rhs_exp_uphi_bar(1,Nr_max) )
                  allocate( wt_rhs_tscheme_exp(1) )  

            end select 

         case ('ARS232') ! Ascher, Ruuth and Spiteri (ARS232) 2nd order, Appl. Numer. Math. 1997 paper
            n_order_tscheme_imp=3 
            allocate( butcher_aD(n_order_tscheme_imp,n_order_tscheme_imp) )
            allocate( butcher_bD(n_order_tscheme_imp) )
            allocate( rhs_imp_temp(1,Nm_max+1,Nr_max) )
            allocate( rhs_imp_vort(1,Nm_max+1,Nr_max) )
            allocate( rhs_imp_uphi_bar(1,Nr_max) )
            allocate( wt_rhs_tscheme_imp(1) ) 
            gammaa = (2.0_dp-sqrt(2.0_dp))/2.0_dp
            delta = -2.0_dp*sqrt(2.0_dp)/3.0_dp

            select case (time_scheme_exp)

              case ('expARS232')
                  n_order_tscheme_exp=3
                  allocate( butcher_aA(n_order_tscheme_exp,n_order_tscheme_exp) )
                  allocate( butcher_bA(n_order_tscheme_exp) )
                  allocate( rhs_exp_temp(1,Nm_max+1,Nr_max) )
                  allocate( rhs_exp_vort(1,Nm_max+1,Nr_max) )
                  allocate( rhs_exp_uphi_bar(1,Nr_max) )
                  allocate( wt_rhs_tscheme_exp(1) )  

            end select  

         case ('ARS233') ! Ascher, Ruuth and Spiteri (ARS233) 3rd order, Appl. Numer. Math. 1997 paper
            n_order_tscheme_imp=3 
            allocate( butcher_aD(n_order_tscheme_imp,n_order_tscheme_imp) )
            allocate( butcher_bD(n_order_tscheme_imp) )
            allocate( rhs_imp_temp(1,Nm_max+1,Nr_max) )
            allocate( rhs_imp_vort(1,Nm_max+1,Nr_max) )
            allocate( rhs_imp_uphi_bar(1,Nr_max) )
            allocate( wt_rhs_tscheme_imp(1) ) 
            gammaa = (3.0_dp+sqrt(3.0_dp))/6.0_dp

            select case (time_scheme_exp)

              case ('expARS233')
                  n_order_tscheme_exp=3
                  allocate( butcher_aA(n_order_tscheme_exp,n_order_tscheme_exp) )
                  allocate( butcher_bA(n_order_tscheme_exp) )
                  allocate( rhs_exp_temp(1,Nm_max+1,Nr_max) )
                  allocate( rhs_exp_vort(1,Nm_max+1,Nr_max) )
                  allocate( rhs_exp_uphi_bar(1,Nr_max) )
                  allocate( wt_rhs_tscheme_exp(1) )  

            end select  

         case ('ARS343') ! Ascher, Ruuth and Spiteri (ARS343) 3rd order, Appl. Numer. Math. 1997 paper
            n_order_tscheme_imp=4 
            allocate( butcher_aD(n_order_tscheme_imp,n_order_tscheme_imp) )
            allocate( butcher_bD(n_order_tscheme_imp) )
            allocate( rhs_imp_temp(1,Nm_max+1,Nr_max) )
            allocate( rhs_imp_vort(1,Nm_max+1,Nr_max) )
            allocate( rhs_imp_uphi_bar(1,Nr_max) )
            allocate( wt_rhs_tscheme_imp(1) ) 

            select case (time_scheme_exp)

              case ('expARS343')
                  n_order_tscheme_exp=4
                  allocate( butcher_aA(n_order_tscheme_exp,n_order_tscheme_exp) )
                  allocate( butcher_bA(n_order_tscheme_exp) )
                  allocate( rhs_exp_temp(1,Nm_max+1,Nr_max) )
                  allocate( rhs_exp_vort(1,Nm_max+1,Nr_max) )
                  allocate( rhs_exp_uphi_bar(1,Nr_max) )
                  allocate( wt_rhs_tscheme_exp(1) )  

            end select   

         case ('MARS343') ! Modified ARS343 Boscarino et. al. (MARS343) 3rd order, Communications to SIMAI Congress. 2007 paper
            n_order_tscheme_imp=4 
            allocate( butcher_aD(n_order_tscheme_imp,n_order_tscheme_imp) )
            allocate( butcher_bD(n_order_tscheme_imp) )
            allocate( rhs_imp_temp(1,Nm_max+1,Nr_max) )
            allocate( rhs_imp_vort(1,Nm_max+1,Nr_max) )
            allocate( rhs_imp_uphi_bar(1,Nr_max) )
            allocate( wt_rhs_tscheme_imp(1) ) 

            select case (time_scheme_exp)

              case ('expMARS343')
                  n_order_tscheme_exp=4
                  allocate( butcher_aA(n_order_tscheme_exp,n_order_tscheme_exp) )
                  allocate( butcher_bA(n_order_tscheme_exp) )
                  allocate( rhs_exp_temp(1,Nm_max+1,Nr_max) )
                  allocate( rhs_exp_vort(1,Nm_max+1,Nr_max) )
                  allocate( rhs_exp_uphi_bar(1,Nr_max) )
                  allocate( wt_rhs_tscheme_exp(1) )  

            end select   

         case ('IJ3') ! Izzo, Jackiewicz (IMEX-RK33pi/2) 3rd order, Appl. Numer. Math. 2016 paper
            n_order_tscheme_imp=4 
            allocate( butcher_aD(n_order_tscheme_imp,n_order_tscheme_imp) )
            allocate( butcher_bD(n_order_tscheme_imp) )
            allocate( rhs_imp_temp(1,Nm_max+1,Nr_max) )
            allocate( rhs_imp_vort(1,Nm_max+1,Nr_max) )
            allocate( rhs_imp_uphi_bar(1,Nr_max) )
            allocate( wt_rhs_tscheme_imp(1) ) 
            lambda = 0.7886866510998523_dp

            select case (time_scheme_exp)

              case ('expIJ3')
                  n_order_tscheme_exp=4
                  allocate( butcher_aA(n_order_tscheme_exp,n_order_tscheme_exp) )
                  allocate( butcher_bA(n_order_tscheme_exp) )
                  allocate( rhs_exp_temp(1,Nm_max+1,Nr_max) )
                  allocate( rhs_exp_vort(1,Nm_max+1,Nr_max) )
                  allocate( rhs_exp_uphi_bar(1,Nr_max) )
                  allocate( wt_rhs_tscheme_exp(1) )  

            end select   

         case ('PC2') ! Predictor-corrector 2nd order, given by N. Schaeffer 
            n_order_tscheme_imp=4 
            allocate( butcher_aD(n_order_tscheme_imp,n_order_tscheme_imp) )
            allocate( butcher_bD(n_order_tscheme_imp) )
            allocate( rhs_imp_temp(1,Nm_max+1,Nr_max) )
            allocate( rhs_imp_vort(1,Nm_max+1,Nr_max) )
            allocate( rhs_imp_uphi_bar(1,Nr_max) )
            allocate( wt_rhs_tscheme_imp(1) ) 

            select case (time_scheme_exp)

              case ('expPC2')
                  n_order_tscheme_exp=4
                  allocate( butcher_aA(n_order_tscheme_exp,n_order_tscheme_exp) )
                  allocate( butcher_bA(n_order_tscheme_exp) )
                  allocate( rhs_exp_temp(1,Nm_max+1,Nr_max) )
                  allocate( rhs_exp_vort(1,Nm_max+1,Nr_max) )
                  allocate( rhs_exp_uphi_bar(1,Nr_max) )
                  allocate( wt_rhs_tscheme_exp(1) )  

            end select      

         case ('SSP332') ! Shu-Chao Duan, https://arxiv.org/pdf/1606.02053.pdf  
            n_order_tscheme_imp=4 
            allocate( butcher_aD(n_order_tscheme_imp,n_order_tscheme_imp) )
            allocate( butcher_bD(n_order_tscheme_imp) )
            allocate( rhs_imp_temp(1,Nm_max+1,Nr_max) )
            allocate( rhs_imp_vort(1,Nm_max+1,Nr_max) )
            allocate( rhs_imp_uphi_bar(1,Nr_max) )
            allocate( wt_rhs_tscheme_imp(1) ) 

            select case (time_scheme_exp)

              case ('expSSP332')
                  n_order_tscheme_exp=4
                  allocate( butcher_aA(n_order_tscheme_exp,n_order_tscheme_exp) )
                  allocate( butcher_bA(n_order_tscheme_exp) )
                  allocate( rhs_exp_temp(1,Nm_max+1,Nr_max) )
                  allocate( rhs_exp_vort(1,Nm_max+1,Nr_max) )
                  allocate( rhs_exp_uphi_bar(1,Nr_max) )
                  allocate( wt_rhs_tscheme_exp(1) )  

            end select      

         case ('SSP543') ! Shu-Chao Duan, https://arxiv.org/pdf/1606.02053.pdf
            n_order_tscheme_imp=6 
            allocate( butcher_aD(n_order_tscheme_imp,n_order_tscheme_imp) )
            allocate( butcher_bD(n_order_tscheme_imp) )
            allocate( rhs_imp_temp(1,Nm_max+1,Nr_max) )
            allocate( rhs_imp_vort(1,Nm_max+1,Nr_max) )
            allocate( rhs_imp_uphi_bar(1,Nr_max) )
            allocate( wt_rhs_tscheme_imp(1) ) 

            select case (time_scheme_exp)

              case ('expSSP543')  
                  n_order_tscheme_exp=6
                  allocate( butcher_aA(n_order_tscheme_exp,n_order_tscheme_exp) )
                  allocate( butcher_bA(n_order_tscheme_exp) )
                  allocate( rhs_exp_temp(1,Nm_max+1,Nr_max) )
                  allocate( rhs_exp_vort(1,Nm_max+1,Nr_max) )
                  allocate( rhs_exp_uphi_bar(1,Nr_max) )
                  allocate( wt_rhs_tscheme_exp(1) )  

            end select 

         case ('ARS443') ! Ascher, Ruuth and Spiteri (ARS443) 3rd order, Appl. Numer. Math. 1997 paper
            n_order_tscheme_imp=5 
            allocate( butcher_aD(n_order_tscheme_imp,n_order_tscheme_imp) )
            allocate( butcher_bD(n_order_tscheme_imp) )
            allocate( rhs_imp_temp(1,Nm_max+1,Nr_max) )
            allocate( rhs_imp_vort(1,Nm_max+1,Nr_max) )
            allocate( rhs_imp_uphi_bar(1,Nr_max) )
            allocate( wt_rhs_tscheme_imp(1) ) 

            select case (time_scheme_exp)

              case ('expARS443')
                  n_order_tscheme_exp=5
                  allocate( butcher_aA(n_order_tscheme_exp,n_order_tscheme_exp) )
                  allocate( butcher_bA(n_order_tscheme_exp) )
                  allocate( rhs_exp_temp(1,Nm_max+1,Nr_max) )
                  allocate( rhs_exp_vort(1,Nm_max+1,Nr_max) )
                  allocate( rhs_exp_uphi_bar(1,Nr_max) )
                  allocate( wt_rhs_tscheme_exp(1) )  

            end select   

         case ('PR433') ! Pareschi and Russo (PR433) 3rd order, 2004 paper
            n_order_tscheme_imp=5 
            allocate( butcher_aD(n_order_tscheme_imp,n_order_tscheme_imp) )
            allocate( butcher_bD(n_order_tscheme_imp) )
            allocate( rhs_imp_temp(1,Nm_max+1,Nr_max) )
            allocate( rhs_imp_vort(1,Nm_max+1,Nr_max) )
            allocate( rhs_imp_uphi_bar(1,Nr_max) )
            allocate( wt_rhs_tscheme_imp(1) ) 
            alphaa=0.24169426078821_dp 
            betaa=0.06042356519705_dp
            etaa=0.12915286960590_dp

            select case (time_scheme_exp)

              case ('expPR433')
                  n_order_tscheme_exp=5
                  allocate( butcher_aA(n_order_tscheme_exp,n_order_tscheme_exp) )
                  allocate( butcher_bA(n_order_tscheme_exp) )
                  allocate( rhs_exp_temp(1,Nm_max+1,Nr_max) )
                  allocate( rhs_exp_vort(1,Nm_max+1,Nr_max) )
                  allocate( rhs_exp_uphi_bar(1,Nr_max) )
                  allocate( wt_rhs_tscheme_exp(1) )  

            end select   

         case ('BPR353') ! Boscarino, Pareschi and Russo (BPR353) 3rd order, SIAM Journal on Numerical Analysis 2013 paper
            n_order_tscheme_imp=5 
            allocate( butcher_aD(n_order_tscheme_imp,n_order_tscheme_imp) )
            allocate( butcher_bD(n_order_tscheme_imp) )
            allocate( rhs_imp_temp(1,Nm_max+1,Nr_max) )
            allocate( rhs_imp_vort(1,Nm_max+1,Nr_max) )
            allocate( rhs_imp_uphi_bar(1,Nr_max) )
            allocate( wt_rhs_tscheme_imp(1) ) 

            select case (time_scheme_exp)

              case ('expBPR353')
                  n_order_tscheme_exp=5
                  allocate( butcher_aA(n_order_tscheme_exp,n_order_tscheme_exp) )
                  allocate( butcher_bA(n_order_tscheme_exp) )
                  allocate( rhs_exp_temp(1,Nm_max+1,Nr_max) )
                  allocate( rhs_exp_vort(1,Nm_max+1,Nr_max) )
                  allocate( rhs_exp_uphi_bar(1,Nr_max) )
                  allocate( wt_rhs_tscheme_exp(1) )  

            end select    

         case ('KC343') ! Additive Runge Kutta (ARK3(2)4L[2]SA) 3rd order, Carpenter and Kennedy, Appl. Numer. Math 2003 paper
            n_order_tscheme_imp=4 
            allocate( butcher_aD(n_order_tscheme_imp,n_order_tscheme_imp) )
            allocate( butcher_bD(n_order_tscheme_imp) )
            allocate( rhs_imp_temp(1,Nm_max+1,Nr_max) )
            allocate( rhs_imp_vort(1,Nm_max+1,Nr_max) )
            allocate( rhs_imp_uphi_bar(1,Nr_max) )
            allocate( wt_rhs_tscheme_imp(1) ) 

            select case (time_scheme_exp)

              case ('expKC343')
                  n_order_tscheme_exp=4
                  allocate( butcher_aA(n_order_tscheme_exp,n_order_tscheme_exp) )
                  allocate( butcher_bA(n_order_tscheme_exp) )
                  allocate( rhs_exp_temp(1,Nm_max+1,Nr_max) )
                  allocate( rhs_exp_vort(1,Nm_max+1,Nr_max) )
                  allocate( rhs_exp_uphi_bar(1,Nr_max) )
                  allocate( wt_rhs_tscheme_exp(1) )  

            end select   

         case ('KC564') ! Additive Runge Kutta (ARK4(3)6L[2]SA) 4th order, Carpenter and Kennedy, Appl. Numer. Math 2003 paper
            n_order_tscheme_imp=6 
            allocate( butcher_aD(n_order_tscheme_imp,n_order_tscheme_imp) )
            allocate( butcher_bD(n_order_tscheme_imp) )
            allocate( rhs_imp_temp(1,Nm_max+1,Nr_max) )
            allocate( rhs_imp_vort(1,Nm_max+1,Nr_max) )
            allocate( rhs_imp_uphi_bar(1,Nr_max) )
            allocate( wt_rhs_tscheme_imp(1) ) 

            select case (time_scheme_exp)

              case ('expKC564')  
                  n_order_tscheme_exp=6
                  allocate( butcher_aA(n_order_tscheme_exp,n_order_tscheme_exp) )
                  allocate( butcher_bA(n_order_tscheme_exp) )
                  allocate( rhs_exp_temp(1,Nm_max+1,Nr_max) )
                  allocate( rhs_exp_vort(1,Nm_max+1,Nr_max) )
                  allocate( rhs_exp_uphi_bar(1,Nr_max) )
                  allocate( wt_rhs_tscheme_exp(1) )  

            end select 

         case ('FW53') ! Fritzen and Wittekindt (FW53) 3rd order, Math. Meth. Appl. Sci. 1997 paper 
            n_order_tscheme_imp=5 
            allocate( butcher_aD(n_order_tscheme_imp,n_order_tscheme_imp) )
            allocate( butcher_bD(n_order_tscheme_imp) )
            allocate( rhs_imp_temp(1,Nm_max+1,Nr_max) )
            allocate( rhs_imp_vort(1,Nm_max+1,Nr_max) )
            allocate( rhs_imp_uphi_bar(1,Nr_max) )
            allocate( wt_rhs_tscheme_imp(1) ) 

            select case (time_scheme_exp)

              case ('expFW53')
                  n_order_tscheme_exp=5
                  allocate( butcher_aA(n_order_tscheme_exp,n_order_tscheme_exp) )
                  allocate( butcher_bA(n_order_tscheme_exp) )
                  allocate( rhs_exp_temp(1,Nm_max+1,Nr_max) )
                  allocate( rhs_exp_vort(1,Nm_max+1,Nr_max) )
                  allocate( rhs_exp_uphi_bar(1,Nr_max) )
                  allocate( wt_rhs_tscheme_exp(1) )  

         end select   

         case ('CK222') ! Boscarino, Pareschi, Russo  (CK222) 2nd order, SIAM Journal on Numerical Analysis 2017 paper 
            n_order_tscheme_imp=3 
            allocate( butcher_aD(n_order_tscheme_imp,n_order_tscheme_imp) )
            allocate( butcher_bD(n_order_tscheme_imp) )
            allocate( rhs_imp_temp(1,Nm_max+1,Nr_max) )
            allocate( rhs_imp_vort(1,Nm_max+1,Nr_max) )
            allocate( rhs_imp_uphi_bar(1,Nr_max) )
            allocate( wt_rhs_tscheme_imp(1) ) 

            select case (time_scheme_exp)

              case ('expCK222')
                  n_order_tscheme_exp=3
                  allocate( butcher_aA(n_order_tscheme_exp,n_order_tscheme_exp) )
                  allocate( butcher_bA(n_order_tscheme_exp) )
                  allocate( rhs_exp_temp(1,Nm_max+1,Nr_max) )
                  allocate( rhs_exp_vort(1,Nm_max+1,Nr_max) )
                  allocate( rhs_exp_uphi_bar(1,Nr_max) )
                  allocate( wt_rhs_tscheme_exp(1) )  

            end select 

         case ('BPR442') ! Boscarino, Pareschi, Russo  (BPR442) 2nd order, SIAM Journal on Numerical Analysis 2017 paper 
            n_order_tscheme_imp=5 
            allocate( butcher_aD(n_order_tscheme_imp,n_order_tscheme_imp) )
            allocate( butcher_bD(n_order_tscheme_imp) )
            allocate( rhs_imp_temp(1,Nm_max+1,Nr_max) )
            allocate( rhs_imp_vort(1,Nm_max+1,Nr_max) )
            allocate( rhs_imp_uphi_bar(1,Nr_max) )
            allocate( wt_rhs_tscheme_imp(1) ) 

            select case (time_scheme_exp)

              case ('expBPR442')
                  n_order_tscheme_exp=5
                  allocate( butcher_aA(n_order_tscheme_exp,n_order_tscheme_exp) )
                  allocate( butcher_bA(n_order_tscheme_exp) )
                  allocate( rhs_exp_temp(1,Nm_max+1,Nr_max) )
                  allocate( rhs_exp_vort(1,Nm_max+1,Nr_max) )
                  allocate( rhs_exp_uphi_bar(1,Nr_max) )
                  allocate( wt_rhs_tscheme_exp(1) )  

            end select

         case ('CFN343') ! Calvo, Frutos, Novo  (CFN343) 3rd order, Applied Numerical Analysis 2001 paper 
            n_order_tscheme_imp=4 
            allocate( butcher_aD(n_order_tscheme_imp,n_order_tscheme_imp) )
            allocate( butcher_bD(n_order_tscheme_imp) )
            allocate( rhs_imp_temp(1,Nm_max+1,Nr_max) )
            allocate( rhs_imp_vort(1,Nm_max+1,Nr_max) )
            allocate( rhs_imp_uphi_bar(1,Nr_max) )
            allocate( wt_rhs_tscheme_imp(1) ) 
            gammaa = 0.435866521508459_dp 

            select case (time_scheme_exp)

              case ('expCFN343')
                  n_order_tscheme_exp=4
                  allocate( butcher_aA(n_order_tscheme_exp,n_order_tscheme_exp) )
                  allocate( butcher_bA(n_order_tscheme_exp) )
                  allocate( rhs_exp_temp(1,Nm_max+1,Nr_max) )
                  allocate( rhs_exp_vort(1,Nm_max+1,Nr_max) )
                  allocate( rhs_exp_uphi_bar(1,Nr_max) )
                  allocate( wt_rhs_tscheme_exp(1) )  

            end select   

         case ('BHR553') ! Boscarino  (BHR553) 3rd order, Applied Numerical Mathematics 2009 paper 
            n_order_tscheme_imp=5 
            allocate( butcher_aD(n_order_tscheme_imp,n_order_tscheme_imp) )
            allocate( butcher_bD(n_order_tscheme_imp) )
            allocate( rhs_imp_temp(1,Nm_max+1,Nr_max) )
            allocate( rhs_imp_vort(1,Nm_max+1,Nr_max) )
            allocate( rhs_imp_uphi_bar(1,Nr_max) )
            allocate( wt_rhs_tscheme_imp(1) ) 
            gammaa = 0.435866521508482_dp
            aa31 = gammaa
            aa32 = aa31
            aa41 = -475883375220285986033264.0_dp/594112726933437845704163.0_dp
            aa42 = 0.0_dp
            aa43 = 1866233449822026827708736.0_dp/594112726933437845704163.0_dp
            aa51 = 62828845818073169585635881686091391737610308247.0_dp/176112910684412105319781630311686343715753056000.0_dp
            aa52 = -302987763081184622639300143137943089.0_dp/1535359944203293318639180129368156500.0_dp
            aa53 = 262315887293043739337088563996093207.0_dp/297427554730376353252081786906492000.0_dp
            aa54 = -987618231894176581438124717087.0_dp/23877337660202969319526901856000.0_dp
            a31 = aa31
            a32 = -31733082319927313.0_dp/455705377221960889379854647102.0_dp
            a41 = -3012378541084922027361996761794919360516301377809610.0_dp/ &
                & 45123394056585269977907753045030512597955897345819349.0_dp
            a42 = -62865589297807153294268.0_dp/102559673441610672305587327019095047.0_dp
            a43 = 418769796920855299603146267001414900945214277000.0_dp/212454360385257708555954598099874818603217167139.0_dp
            bb1 = 487698502336740678603511.0_dp/1181159636928185920260208.0_dp
            bb2 = 0.0_dp
            bb3 = -1.0_dp*aa52
            bb4 = -105235928335100616072938218863.0_dp/2282554452064661756575727198000.0_dp
            c3 = 902905985686.0_dp/1035759735069.0_dp
            c4 = 2684624.0_dp/1147171.0_dp
            
            select case (time_scheme_exp)

              case ('expBHR553')
                  n_order_tscheme_exp=5
                  allocate( butcher_aA(n_order_tscheme_exp,n_order_tscheme_exp) )
                  allocate( butcher_bA(n_order_tscheme_exp) )
                  allocate( rhs_exp_temp(1,Nm_max+1,Nr_max) )
                  allocate( rhs_exp_vort(1,Nm_max+1,Nr_max) )
                  allocate( rhs_exp_uphi_bar(1,Nr_max) )
                  allocate( wt_rhs_tscheme_exp(1) )  

            end select

         case ('CB3e') ! Cavaglieri and Bewley (IMEXRKCB3e) 3rd order, JCP 2015  
            n_order_tscheme_imp=4 
            allocate( butcher_aD(n_order_tscheme_imp,n_order_tscheme_imp) )
            allocate( butcher_bD(n_order_tscheme_imp) )
            allocate( rhs_imp_temp(1,Nm_max+1,Nr_max) )
            allocate( rhs_imp_vort(1,Nm_max+1,Nr_max) )
            allocate( rhs_imp_uphi_bar(1,Nr_max) )
            allocate( wt_rhs_tscheme_imp(1) ) 

            select case (time_scheme_exp)

              case ('expCB3e')
                  n_order_tscheme_exp=4
                  allocate( butcher_aA(n_order_tscheme_exp,n_order_tscheme_exp) )
                  allocate( butcher_bA(n_order_tscheme_exp) )
                  allocate( rhs_exp_temp(1,Nm_max+1,Nr_max) )
                  allocate( rhs_exp_vort(1,Nm_max+1,Nr_max) )
                  allocate( rhs_exp_uphi_bar(1,Nr_max) )
                  allocate( wt_rhs_tscheme_exp(1) )  

            end select      

         case ('CB3f') ! Cavaglieri and Bewley (IMEXRKCB3f) 3rd order, JCP 2015  
            n_order_tscheme_imp=4 
            allocate( butcher_aD(n_order_tscheme_imp,n_order_tscheme_imp) )
            allocate( butcher_bD(n_order_tscheme_imp) )
            allocate( rhs_imp_temp(1,Nm_max+1,Nr_max) )
            allocate( rhs_imp_vort(1,Nm_max+1,Nr_max) )
            allocate( rhs_imp_uphi_bar(1,Nr_max) )
            allocate( wt_rhs_tscheme_imp(1) ) 
            c2=49.0_dp/50.0_dp
            a21=c2/2.0_dp
            a22=c2/2.0_dp
            a31=-785157464198.0_dp/1093480182337.0_dp
            a32=-30736234873.0_dp/978681420651.0_dp
            a33=983779726483.0_dp/1246172347126.0_dp
            bb1=-2179897048956.0_dp/603118880443.0_dp
            bb2=99189146040.0_dp/891495457793.0_dp
            bb3=6064140186914.0_dp/1415701440113.0_dp
            bb4=146791865627.0_dp/668377518349.0_dp
            aa21=c2
            aa31=13244205847.0_dp/647648310246.0_dp
            aa32=13419997131.0_dp/686433909488.0_dp
            aa42=231677526244.0_dp/1085522130027.0_dp
            aa43=3007879347537.0_dp/683461566472.0_dp

            select case (time_scheme_exp)

              case ('expCB3f')
                  n_order_tscheme_exp=4
                  allocate( butcher_aA(n_order_tscheme_exp,n_order_tscheme_exp) )
                  allocate( butcher_bA(n_order_tscheme_exp) )
                  allocate( rhs_exp_temp(1,Nm_max+1,Nr_max) )
                  allocate( rhs_exp_vort(1,Nm_max+1,Nr_max) )
                  allocate( rhs_exp_uphi_bar(1,Nr_max) )
                  allocate( wt_rhs_tscheme_exp(1) )  

            end select      

         case ('CB4') ! Cavaglieri and Bewley (IMEXRKCB4) 3rd order, JCP 2015
            n_order_tscheme_imp=6 
            allocate( butcher_aD(n_order_tscheme_imp,n_order_tscheme_imp) )
            allocate( butcher_bD(n_order_tscheme_imp) )
            allocate( rhs_imp_temp(1,Nm_max+1,Nr_max) )
            allocate( rhs_imp_vort(1,Nm_max+1,Nr_max) )
            allocate( rhs_imp_uphi_bar(1,Nr_max) )
            allocate( wt_rhs_tscheme_imp(1) ) 
            c2=0.25_dp
            a21=c2/2.0_dp
            a22=c2/2.0_dp
            a31=216145252607.0_dp/961230882893.0_dp
            a32=257479850128.0_dp/1143310606989.0_dp
            a33=30481561667.0_dp/101628412017.0_dp
            a42=-381180097479.0_dp/1276440792700.0_dp
            a43=-54660926949.0_dp/461115766612.0_dp
            a44=344309628413.0_dp/552073727558.0_dp
            a53=-100836174740.0_dp/861952129159.0_dp
            a54=-250423827953.0_dp/1283875864443.0_dp
            a55=0.5_dp
            aa21=c2 
            aa31=153985248130.0_dp/1004999853329.0_dp
            aa32=902825336800.0_dp/1512825644809.0_dp
            aa42=99316866929.0_dp/820744730663.0_dp
            aa43=82888780751.0_dp/969573940619.0_dp
            aa53=57501241309.0_dp/765040883867.0_dp
            aa54=76345938311.0_dp/676824576433.0_dp
            aa64=-4099309936455.0_dp/6310162971841.0_dp
            aa65=1395992540491.0_dp/933264948679.0_dp
            bb1=232049084587.0_dp/1377130630063.0_dp
            bb2=322009889509.0_dp/2243393849156.0_dp
            bb3=-195109672787.0_dp/1233165545817.0_dp 
            bb4=-340582416761.0_dp/705418832319.0_dp
            bb5=463396075661.0_dp/409972144477.0_dp
            bb6=323177943294.0_dp/1626646580633.0_dp

            select case (time_scheme_exp)

              case ('expCB4')  
                  n_order_tscheme_exp=6
                  allocate( butcher_aA(n_order_tscheme_exp,n_order_tscheme_exp) )
                  allocate( butcher_bA(n_order_tscheme_exp) )
                  allocate( rhs_exp_temp(1,Nm_max+1,Nr_max) )
                  allocate( rhs_exp_vort(1,Nm_max+1,Nr_max) )
                  allocate( rhs_exp_uphi_bar(1,Nr_max) )
                  allocate( wt_rhs_tscheme_exp(1) )  

            end select 

         case ('LZ3L1') ! Liu, Zou  (LZ3L1) 3rd order, Journal of Computational and Applied Math 2005 paper 
            n_order_tscheme_imp=5 
            allocate( butcher_aD(n_order_tscheme_imp,n_order_tscheme_imp) )
            allocate( butcher_bD(n_order_tscheme_imp) )
            allocate( rhs_imp_temp(1,Nm_max+1,Nr_max) )
            allocate( rhs_imp_vort(1,Nm_max+1,Nr_max) )
            allocate( rhs_imp_uphi_bar(1,Nr_max) )
            allocate( wt_rhs_tscheme_imp(1) ) 

            select case (time_scheme_exp)

              case ('expLZ3L1')
                  n_order_tscheme_exp=5
                  allocate( butcher_aA(n_order_tscheme_exp,n_order_tscheme_exp) )
                  allocate( butcher_bA(n_order_tscheme_exp) )
                  allocate( rhs_exp_temp(1,Nm_max+1,Nr_max) )
                  allocate( rhs_exp_vort(1,Nm_max+1,Nr_max) )
                  allocate( rhs_exp_uphi_bar(1,Nr_max) )
                  allocate( wt_rhs_tscheme_exp(1) )  

            end select

         case ('LZ3A4a') ! Liu, Zou  (LZ3A4a) 3rd order, Journal of Computational and Applied Math 2005 paper 
            n_order_tscheme_imp=5 
            allocate( butcher_aD(n_order_tscheme_imp,n_order_tscheme_imp) )
            allocate( butcher_bD(n_order_tscheme_imp) )
            allocate( rhs_imp_temp(1,Nm_max+1,Nr_max) )
            allocate( rhs_imp_vort(1,Nm_max+1,Nr_max) )
            allocate( rhs_imp_uphi_bar(1,Nr_max) )
            allocate( wt_rhs_tscheme_imp(1) ) 

            select case (time_scheme_exp)

              case ('expLZ3A4a')
                  n_order_tscheme_exp=5
                  allocate( butcher_aA(n_order_tscheme_exp,n_order_tscheme_exp) )
                  allocate( butcher_bA(n_order_tscheme_exp) )
                  allocate( rhs_exp_temp(1,Nm_max+1,Nr_max) )
                  allocate( rhs_exp_vort(1,Nm_max+1,Nr_max) )
                  allocate( rhs_exp_uphi_bar(1,Nr_max) )
                  allocate( wt_rhs_tscheme_exp(1) )  

            end select


         case ('CFN564') ! Calvo, Frutos, Novo  (CFN564) 4th order, Applied Numerical Analysis 2001 paper 
            n_order_tscheme_imp=6 
            allocate( butcher_aD(n_order_tscheme_imp,n_order_tscheme_imp) )
            allocate( butcher_bD(n_order_tscheme_imp) )
            allocate( rhs_imp_temp(1,Nm_max+1,Nr_max) )
            allocate( rhs_imp_vort(1,Nm_max+1,Nr_max) )
            allocate( rhs_imp_uphi_bar(1,Nr_max) )
            allocate( wt_rhs_tscheme_imp(1) ) 

            select case (time_scheme_exp)

              case ('expCFN564')  
                  n_order_tscheme_exp=6
                  allocate( butcher_aA(n_order_tscheme_exp,n_order_tscheme_exp) )
                  allocate( butcher_bA(n_order_tscheme_exp) )
                  allocate( rhs_exp_temp(1,Nm_max+1,Nr_max) )
                  allocate( rhs_exp_vort(1,Nm_max+1,Nr_max) )
                  allocate( rhs_exp_uphi_bar(1,Nr_max) )
                  allocate( wt_rhs_tscheme_exp(1) )  

            end select 

         case ('LZ4A1')  ! Liu, Zou  (LZ4A1) 3rd order, Journal of Computational and Applied Math 2005 paper 
            n_order_tscheme_imp=7 
            allocate( butcher_aD(n_order_tscheme_imp,n_order_tscheme_imp) )
            allocate( butcher_bD(n_order_tscheme_imp) )
            allocate( rhs_imp_temp(1,Nm_max+1,Nr_max) )
            allocate( rhs_imp_vort(1,Nm_max+1,Nr_max) )
            allocate( rhs_imp_uphi_bar(1,Nr_max) )
            allocate( wt_rhs_tscheme_imp(1) ) 

            select case (time_scheme_exp)

              case ('expLZ4A1')  
                  n_order_tscheme_exp=7
                  allocate( butcher_aA(n_order_tscheme_exp,n_order_tscheme_exp) )
                  allocate( butcher_bA(n_order_tscheme_exp) )
                  allocate( rhs_exp_temp(1,Nm_max+1,Nr_max) )
                  allocate( rhs_exp_vort(1,Nm_max+1,Nr_max) )
                  allocate( rhs_exp_uphi_bar(1,Nr_max) )
                  allocate( wt_rhs_tscheme_exp(1) )  

            end select 

      end select
       
      !------- Allocate dt array for time schemes of order > 1 --------- 
      n_order_tscheme_max = max(n_order_tscheme_imp,n_order_tscheme_exp)
      !if (time_scheme_type=='IMEX') then
         allocate( dt_array(n_order_tscheme_max) )
      !else
      !   allocate( dt_array(2) )
      !end if
      !-----------------------------------------------------------------

   end subroutine init_time_schemes

   subroutine deallocate_time_schemes(time_scheme_type)
      character(len=100), intent(in) :: time_scheme_type

      deallocate( dt_array )
      if (time_scheme_type/='RK') then
         deallocate( rhs_imp_vort )
         deallocate( rhs_imp_temp )
      end if
      if (time_scheme_type=='IMEX') then
         deallocate( rhs_buo_term )
      end if
         deallocate( rhs_exp_vort )
         deallocate( rhs_exp_temp )
         deallocate( wt_rhs_tscheme_imp )
         deallocate( wt_rhs_tscheme_exp )

   end subroutine deallocate_time_schemes

   subroutine rhs_update_wts_imp(time_scheme_imp,wt_lhs_tscheme_imp,wt_rhs_tscheme_imp,n_order_tscheme_imp)

      integer, intent(in) :: n_order_tscheme_imp
      character(len=100), intent(in) :: time_scheme_imp
      real(kind=dp), intent(out) :: wt_rhs_tscheme_imp(n_order_tscheme_imp)
      real(kind=dp), intent(out) :: wt_lhs_tscheme_imp
      real(kind=dp) :: a0,a1,a2,a3,a4,r_dt1,r_dt2,r_dt3
      integer :: i,j

      select case (time_scheme_imp)

         case ('CN') 
            wt_rhs_tscheme_imp(1) = half
            wt_lhs_tscheme_imp=half

         case ('BDF2')
            r_dt1=dt_array(2)/dt_array(1)
            a0=(2.0_dp+r_dt1)/(1.0_dp+r_dt1)
            a1=1.0_dp+1.0_dp/r_dt1
            a2=-1.0_dp*(1.0_dp/(r_dt1)/(1.0_dp+r_dt1)) ! Thomas found mistake here to fix on Sept 20th 2017 (fixed)
            wt_lhs_tscheme_imp=1.0_dp/a0
            wt_rhs_tscheme_imp(1) = a1/a0
            wt_rhs_tscheme_imp(2) = a2/a0

         case ('BDF3')
            r_dt1=dt_array(2)/dt_array(1)
            r_dt2=dt_array(3)/dt_array(1)
            a0=(1.0_dp + 1.0_dp/(1.0_dp+r_dt1) + 1.0_dp/(1.0_dp+r_dt1+r_dt2))
            a1=(((1.0_dp+r_dt1)*(1.0_dp+r_dt1+r_dt2))/(r_dt1*(r_dt1+r_dt2)))
            a2=(-1.0_dp*(1.0_dp+r_dt1+r_dt2)/(r_dt1*r_dt2*(1.0_dp+r_dt1)))
            a3=((1.0_dp+r_dt1)/(r_dt2*(r_dt1+r_dt2)*(1.0_dp+r_dt1+r_dt2)))
            wt_lhs_tscheme_imp=1.0_dp/a0
            wt_rhs_tscheme_imp(1) = a1/a0
            wt_rhs_tscheme_imp(2) = a2/a0
            wt_rhs_tscheme_imp(3) = a3/a0

         case ('BDF4')
            r_dt1=dt_array(1)/dt_array(2)
            r_dt2=dt_array(2)/dt_array(3)
            r_dt3=dt_array(3)/dt_array(4)

            c1 = one+r_dt3*(one+r_dt2) 
            c2 = one+r_dt2*(one+r_dt1) 
            c3 = one+r_dt3*c2 

            a0= one + r_dt1/(one+r_dt1)+r_dt2*r_dt1/c2+r_dt3* &
            & r_dt2*r_dt1/c3
            a1 = -one-r_dt1*(one+r_dt2*(one+r_dt1)*(one+r_dt3*c2/c1)/ &
            &    (one+r_dt2))
            a2 = r_dt1*(r_dt1/(one+r_dt1)+r_dt2*r_dt1 * (c3+r_dt3)/&
            &    (one+r_dt3))
            a3 = -r_dt2**3*r_dt1**2*(one+r_dt1)/(one+r_dt2) * c3/c2
            a4 = (one+r_dt1)/(one+r_dt3) * c2/c1/c3 *(r_dt3**4* &
            &    r_dt2**3*r_dt1**2)

            wt_lhs_tscheme_imp = 1.0_dp/a0
            wt_rhs_tscheme_imp(1) = -a1/a0
            wt_rhs_tscheme_imp(2) = -a2/a0
            wt_rhs_tscheme_imp(3) = -a3/a0
            wt_rhs_tscheme_imp(4) = -a4/a0
            !print *,  wt_lhs_tscheme_imp, wt_rhs_tscheme_imp(1), wt_rhs_tscheme_imp(2), wt_rhs_tscheme_imp(3), &
            !& wt_rhs_tscheme_imp(4)
         case ('EUL111')

            butcher_aD = reshape([0.0_dp, 0.0_dp, &
                                  0.0_dp, 1.0_dp], &
                                  [2,2],order=[2,1]) 

            butcher_bD = reshape([0.0_dp, 1.0_dp],[2]) 

            diag_diff = .FALSE.
            diag_index = 2

         case ('ARS222')

            butcher_aD = reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                                  0.0_dp, gammaa, 0.0_dp, &
                                  0.0_dp, 1.0_dp-gammaa, gammaa], &
                                  [3,3],order=[2,1])  

            butcher_bD = reshape([0.0_dp, 1.0_dp-gammaa, gammaa],[3]) 

            do i=1,n_order_tscheme_imp
               if (butcher_aD(n_order_tscheme_imp,i)==butcher_bD(i)) then
                  ars_eqn_check_D=.TRUE.
               end if
            end do 

            diag_diff = .FALSE.
            diag_index = 2

         case ('ARS232')

            butcher_aD = reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                                  0.0_dp, gammaa, 0.0_dp, &
                                  0.0_dp, 1.0_dp-gammaa, gammaa], &
                                  [3,3],order=[2,1])  

            butcher_bD = reshape([0.0_dp, 1.0_dp-gammaa, gammaa],[3]) 

            do i=1,n_order_tscheme_imp
               if (butcher_aD(n_order_tscheme_imp,i)==butcher_bD(i)) then
                  ars_eqn_check_D=.TRUE.
               end if
            end do  

            diag_diff = .FALSE.
            diag_index = 2

         case ('ARS233')

            butcher_aD = reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                                  0.0_dp, gammaa, 0.0_dp, &
                                  0.0_dp, 1.0_dp-2.0_dp*gammaa, gammaa], &
                                  [3,3],order=[2,1])  

            butcher_bD = reshape([0.0_dp, 0.5_dp, 0.5_dp],[3]) 

            do i=1,n_order_tscheme_imp
               if (butcher_aD(n_order_tscheme_imp,i)==butcher_bD(i)) then
                  ars_eqn_check_D=.TRUE.
               end if
            end do  

            diag_diff = .FALSE.
            diag_index = 2

         case ('ARS343')

            butcher_aD = reshape([0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                                  0.0_dp, 0.4358665215_dp, 0.0_dp, 0.0_dp, &
                                  0.0_dp, 0.2820667392_dp, 0.4358665215_dp, 0.0_dp, &
                                  0.0_dp, 1.208496649_dp, -0.644363171_dp, 0.4358665215_dp], &
                                  [4,4],order=[2,1])  

            butcher_bD = reshape([0.0_dp, 1.208496649_dp, -0.644363171_dp, 0.4358665215_dp],[4]) 

            do i=1,n_order_tscheme_imp
               if (butcher_aD(n_order_tscheme_imp,i)==butcher_bD(i)) then
                  ars_eqn_check_D=.TRUE.
               end if
            end do   

            diag_diff = .FALSE.
            diag_index = 2

         case ('MARS343')

            butcher_aD = reshape([0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                                  0.0_dp, 0.435866521508458_dp, 0.0_dp, 0.0_dp, &
                                  0.0_dp, 0.28206673924577_dp, 0.435866521508458_dp, 0.0_dp, &
                                  0.0_dp, 1.20849664917601_dp, -0.64436317068446_dp, 0.435866521508458_dp], &
                                  [4,4],order=[2,1])  

            butcher_bD = reshape([0.0_dp, 1.20849664917601_dp, -0.64436317068446_dp, 0.435866521508458_dp],[4]) 

            ars_eqn_check_D=.TRUE.

            diag_diff = .FALSE.
            diag_index = 2

         case ('IJ3')

            butcher_aD = reshape([0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                                  0.0_dp, lambda, 0.0_dp, 0.0_dp, &
                  0.0_dp, lambda/(3.0_dp*(1.0_dp-2.0_dp*lambda)), (1.0_dp-3.0_dp*lambda)/(3.0_dp*(1.0_dp-2.0_dp*lambda)), 0.0_dp, &
                  0.0_dp, -1.0_dp*lambda/(1.0_dp-2.0_dp*lambda), (1.0_dp-lambda)/(1.0_dp-2.0_dp*lambda), 0.0_dp], &
                                  [4,4],order=[2,1])  

            butcher_bD = reshape([0.0_dp, 0.0_dp, 0.75_dp, 0.25_dp],[4]) 
            ars_eqn_check_D=.FALSE.

            diag_diff = .FALSE.
            diag_index = 2

         case ('PC2')

            butcher_aD = reshape([0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                                  0.5_dp, 0.5_dp, 0.0_dp, 0.0_dp, &
                                  0.5_dp, 0.0_dp, 0.5_dp, 0.0_dp, &
                                  0.5_dp, 0.0_dp, 0.0_dp, 0.5_dp], &
                                  [4,4],order=[2,1])  

            butcher_bD = reshape([0.5_dp, 0.0_dp, 0.0_dp, 0.5_dp],[4]) 

            ars_eqn_check_D=.TRUE.

            diag_diff = .FALSE.
            diag_index = 2

         case ('SSP332')

            butcher_aD = reshape([0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                                  0.0_dp, 0.5_dp, 0.0_dp, 0.0_dp, &
                                  0.0_dp, 0.5_dp, 0.5_dp, 0.0_dp, &
                                  0.0_dp, 1.0_dp, -0.5_dp, 0.5_dp], &
                                  [4,4],order=[2,1])  

            butcher_bD = reshape([0.0_dp, 1.0_dp, -0.5_dp, 0.5_dp],[4]) 

            ars_eqn_check_D=.TRUE.

            diag_diff = .FALSE.
            diag_index = 2

         case ('SSP543')
  
            butcher_aD = reshape([0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                                  0.0_dp, 0.5_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                                  0.0_dp, 0.5_dp, 0.5_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                                  0.0_dp, 0.5_dp, -0.5_dp, 0.5_dp, 0.0_dp, 0.0_dp, &
                                  0.0_dp, 2.0_dp, 0.5_dp, -3.0_dp, 0.5_dp, 0.0_dp, &
                                  0.0_dp, 2.0_dp/3.0_dp, -1.0_dp/3.0_dp, 0.0_dp, 1.0_dp/6.0_dp, 0.5_dp], &
                                  [6,6],order=[2,1]) 

            butcher_bD = reshape([0.0_dp, 2.0_dp/3.0_dp, -1.0_dp/3.0_dp, 0.0_dp, 1.0_dp/6.0_dp, 0.5_dp],[6]) 

            ars_eqn_check_D=.TRUE.

            diag_diff = .FALSE.
            diag_index = 2
 
         case ('ARS443')
            
            butcher_aD = reshape([ 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                                   0.0_dp, 0.5_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                                   0.0_dp,1.0_dp/6.0_dp, 0.5_dp, 0.0_dp,0.0_dp, &
                                   0.0_dp, -0.5_dp, 0.5_dp, 0.5_dp, 0.0_dp, &
                                   0.0_dp, 1.5_dp, -1.5_dp, 0.5_dp, 0.5_dp], &
                                   [5,5],order=[2,1]) 

            butcher_bD = reshape([0.0_dp, 1.5_dp, -1.5_dp, 0.5_dp, 0.5_dp],[5]) 

            do i=1,n_order_tscheme_imp
               if (butcher_aD(n_order_tscheme_imp,i)==butcher_bD(i)) then
                  ars_eqn_check_D=.TRUE.
               end if
            end do

            diag_diff = .FALSE.
            diag_index = 2

         case ('PR433')
            
            butcher_aD = reshape([ 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                                   0.0_dp, alphaa, 0.0_dp, 0.0_dp, 0.0_dp, &
                                   0.0_dp,-alphaa, alphaa, 0.0_dp,0.0_dp, &
                                   0.0_dp, 0.0_dp, 1.0_dp-alphaa, alphaa, 0.0_dp, &
                                   0.0_dp, betaa, etaa, 0.5_dp-betaa-etaa-alphaa, alphaa], &
                                   [5,5],order=[2,1]) 

            !butcher_aD = reshape([ alphaa, 0.0_dp, 0.0_dp, 0.0_dp, &
            !                       -alphaa, alphaa, 0.0_dp,0.0_dp, &
            !                       0.0_dp, 1.0_dp-alphaa, alphaa, 0.0_dp, &
            !                       betaa, etaa, 0.5_dp-betaa-etaa-alphaa, alphaa], &
            !                       [4,4],order=[2,1]) 

            !butcher_bD = reshape([0.0_dp, 1.0_dp/6.0_dp, 1.0_dp/6.0_dp, 2.0_dp/3.0_dp],[4]) 
            butcher_bD = reshape([0.0_dp,0.0_dp, 1.0_dp/6.0_dp, 1.0_dp/6.0_dp, 2.0_dp/3.0_dp],[5]) 

            ars_eqn_check_D=.FALSE.

            diag_diff = .FALSE.
            diag_index = 2

         case ('BPR353')

            butcher_aD = reshape([0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                                  0.5_dp, 0.5_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                                  5.0_dp/18.0_dp, -1.0_dp/9.0_dp, 0.5_dp, 0.0_dp, 0.0_dp, &
                                  0.5_dp, 0.0_dp, 0.0_dp, 0.5_dp,0.0_dp, &
                                  0.25_dp,0.0_dp,0.75_dp,-0.5_dp,0.5_dp], &
                                  [5,5],order=[2,1])

            butcher_bD = reshape([0.25_dp,0.0_dp,0.75_dp,-0.5_dp,0.5_dp],[5])

            do i=1,n_order_tscheme_imp
               if (butcher_aD(n_order_tscheme_imp,i)==butcher_bD(i)) then
                  ars_eqn_check_D=.TRUE.
               end if
            end do

            diag_diff = .FALSE.
            diag_index = 2

         case ('KC343')
 
            butcher_aD = reshape([0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                   1767732205903.0_dp/4055673282236.0_dp, 1767732205903.0_dp/4055673282236.0_dp, 0.0_dp, 0.0_dp, &
                   2746238789719.0_dp/10658868560708.0_dp, -640167445237.0_dp/6845629431997.0_dp, &
                  & 1767732205903.0_dp/4055673282236.0_dp, 0.0_dp, &
            1471266399579.0_dp/7840856788654.0_dp, -4482444167858.0_dp/7529755066697.0_dp, &
            & 11266239266428.0_dp/11593286722821.0_dp,1767732205903.0_dp/4055673282236.0_dp], &
                                  [4,4],order=[2,1]) 

            butcher_bD = reshape([1471266399579.0_dp/7840856788654.0_dp, -4482444167858.0_dp/7529755066697.0_dp, &
            & 11266239266428.0_dp/11593286722821.0_dp,1767732205903.0_dp/4055673282236.0_dp],[4]) 
            !butcher_bD = reshape([2756255671327.0_dp/12835298489170.0_dp, -10771552573575.0_dp/22201958757719.0_dp, &
            !& 9247589265047.0_dp/10645013368117.0_dp,2193209047091.0_dp/5459859503100.0_dp],[4]) 
            !butcher_bD = reshape([-215264564351.0_dp/13552729205753.0_dp,17870216137069.0_dp/13817060693119.0_dp, &
            !& -28141676662227.0_dp/17317692491321.0_dp,2508943948391.0_dp/7218656332882.0_dp],[4]) 

            do i=1,n_order_tscheme_imp
               if (butcher_aD(n_order_tscheme_imp,i)==butcher_bD(i)) then
                  ars_eqn_check_D=.TRUE.
               end if
            end do

            diag_diff = .FALSE.
            diag_index = 2

         case ('KC564')
  
            butcher_aD = reshape([0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                                  0.25_dp, 0.25_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                                  8611.0_dp/62500.0_dp, -1743.0_dp/31250.0_dp, 0.25_dp, 0.0_dp,0.0_dp, 0.0_dp, &
            5012029.0_dp/34652500.0_dp, -654441.0_dp/2922500.0_dp, 174375.0_dp/388108.0_dp, 0.25_dp, 0.0_dp, 0.0_dp, &
            15267082809.0_dp/155376265600.0_dp, -71443401.0_dp/120774400.0_dp, 730878875.0_dp/902184768.0_dp, & 
            & 2285395.0_dp/8070912.0_dp, 0.25_dp, 0.0_dp, &
      82889.0_dp/524892.0_dp, 0.0_dp, 15625.0_dp/83664.0_dp, 69875.0_dp/102672.0_dp, -2260.0_dp/8211.0_dp, 0.25_dp], &
                                  [6,6],order=[2,1]) 

            butcher_bD = reshape([82889.0_dp/524892.0_dp, 0.0_dp, 15625.0_dp/83664.0_dp, 69875.0_dp/102672.0_dp, & 
                                  & -2260.0_dp/8211.0_dp, 0.25_dp],[6]) 

            !do i=1,n_order_tscheme_imp
            !   if (butcher_aD(n_order_tscheme_imp,i)==butcher_bD(i)) then
                  ars_eqn_check_D=.TRUE.
            !      ars_eqn_check_D=.TRUE.
            !   end if
            !end do

            diag_diff = .FALSE.
            diag_index = 2
            
         case ('FW53')
            
            butcher_aD = reshape([0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                                        0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                                        0.5_dp, 0.0_dp, 0.5_dp, 0.0_dp, 0.0_dp, &
                                        0.25_dp, -0.25_dp, 0.0_dp, 0.5_dp, 0.0_dp, &
                                        1.0_dp/6.0_dp, 1.0_dp/6.0_dp, -1.0_dp/6.0_dp, 4.0_dp/6.0_dp, 1.0_dp/6.0_dp], &
                                        [5,5],order=[2,1])

            butcher_bD = reshape([1.0_dp/6.0_dp, 1.0_dp/6.0_dp, -1.0_dp/6.0_dp, 4.0_dp/6.0_dp, 1.0_dp/6.0_dp],[5]) 

            do i=1,n_order_tscheme_imp
               if (butcher_aD(n_order_tscheme_imp,i)==butcher_bD(i)) then
                  ars_eqn_check_D=.TRUE.
               end if
            end do

            diag_diff = .FALSE.
            diag_index = 2

         case ('CK222')

            butcher_aD = reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                                  -onethird+sqrt(2.0_dp)*0.5_dp, 1.0_dp-sqrt(2.0_dp)*0.5_dp, 0.0_dp, &
                   threefourth-sqrt(2.0_dp)*0.25_dp, -threefourth+sqrt(2.0_dp)*0.75_dp , 1.0_dp-sqrt(2.0_dp)*0.5_dp], &
                                  [3,3],order=[2,1])  

            butcher_bD = reshape([threefourth-sqrt(2.0_dp)*0.25_dp, -threefourth+sqrt(2.0_dp)*0.75_dp , &
                  1.0_dp-sqrt(2.0_dp)*0.5_dp],[3]) 

            do i=1,n_order_tscheme_imp
               if (butcher_aD(n_order_tscheme_imp,i)==butcher_bD(i)) then
                  ars_eqn_check_D=.TRUE.
               end if
            end do 

            diag_diff = .FALSE.
            diag_index = 2

         case ('BPR442')
            
            butcher_aD = reshape([ 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                                   0.0_dp, 0.25_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                                   0.0_dp, 0.0_dp, 0.25_dp, 0.0_dp, 0.0_dp, &
                                   0.0_dp, 1.0_dp/24.0_dp, 11.0_dp/24.0_dp, 0.25_dp, 0.0_dp, &
                                   0.0_dp, 11.0_dp/24.0_dp, 1.0_dp/6.0_dp, 1.0_dp/8.0_dp, 0.25_dp], &
                                   [5,5],order=[2,1]) 

            butcher_bD = reshape([0.0_dp, 11.0_dp/24.0_dp, 1.0_dp/6.0_dp, 1.0_dp/8.0_dp, 0.25_dp],[5]) 

            ars_eqn_check_D=.TRUE.

            diag_diff = .FALSE.
            diag_index = 2

         case ('CFN343') 

            butcher_aD = reshape([0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                                  0.0_dp, gammaa, 0.0_dp, 0.0_dp, &
                                  0.0_dp, (1.0_dp-gammaa)/2.0_dp, gammaa, 0.0_dp, &
                                  0.0_dp, 1.208496649176010_dp, -0.644363170684469_dp, gammaa], &
                                  [4,4],order=[2,1])  

            butcher_bD = reshape([0.0_dp, 1.208496649176010_dp, -0.644363170684469_dp, gammaa],[4]) 

            ars_eqn_check_D=.TRUE.

            diag_diff = .FALSE.
            diag_index = 2

         case ('BHR553')
            
            butcher_aD = reshape([ 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                                   gammaa, gammaa, 0.0_dp, 0.0_dp, 0.0_dp, &
                                   a31, a32, gammaa, 0.0_dp, 0.0_dp, &
                                   a41, a42, a43, gammaa, 0.0_dp, &
                                   bb1, 0.0_dp, bb3, bb4, gammaa], &
                                   [5,5],order=[2,1]) 

            butcher_bD = reshape([bb1, 0.0_dp, bb3, bb4, gammaa],[5]) 

            ars_eqn_check_D=.TRUE.

            diag_diff = .FALSE.
            diag_index = 2

         case ('CB3e')

            butcher_aD = reshape([0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                                  0.0_dp, 1.0_dp/3.0_dp, 0.0_dp, 0.0_dp, &
                                  0.0_dp, 0.5_dp, 0.5_dp, 0.0_dp, &
                                  0.0_dp, 0.75_dp, -0.25_dp, 0.5_dp], &
                                  [4,4],order=[2,1])  

            butcher_bD = reshape([0.0_dp, 0.75_dp, -0.25_dp, 0.5_dp],[4]) 

            ars_eqn_check_D=.TRUE.

            diag_diff = .TRUE.
            diag_index = 3

         case ('CB3f')

            butcher_aD = reshape([0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                                  a21, a22, 0.0_dp, 0.0_dp, &
                                  a31, a32, a33, 0.0_dp, &
                                  bb1, bb2, bb3, bb4], &
                                  [4,4],order=[2,1])  

            butcher_bD = reshape([bb1, bb2, bb3, bb4],[4]) 

            ars_eqn_check_D=.TRUE.

            diag_diff = .TRUE.
            diag_index = 2

         case ('CB4')
  
            butcher_aD = reshape([0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                                  a21, a22, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                                  a31, a32, a33, 0.0_dp, 0.0_dp, 0.0_dp, &
                                  bb1, a42, a43, a44, 0.0_dp, 0.0_dp, &
                                  bb1, bb2, a53, a54, a55, 0.0_dp, &
                                  bb1, bb2, bb3, bb4, bb5, bb6], &
                                  [6,6],order=[2,1]) 

            butcher_bD = reshape([bb1, bb2, bb3, bb4, bb5, bb6],[6]) 
          
            ars_eqn_check_D=.TRUE.

            diag_diff = .TRUE.
            diag_index = 2

         case ('LZ3L1')
            
            butcher_aD = reshape([ 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                                   3.0_dp/20.0_dp, 1.0_dp/10.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                                   9.0_dp/10.0_dp, -13.0_dp/10.0_dp, 9.0_dp/10.0_dp, 0.0_dp, 0.0_dp, &
                                   17.0_dp/10.0_dp, -11.0_dp/4.0_dp, 3.0_dp/2.0_dp, 3.0_dp/10.0_dp, 0.0_dp, &
                                   1.0_dp, -10.0_dp/3.0_dp, 17.0_dp/3.0_dp, -10.0_dp/3.0_dp, 1.0_dp], &
                                   [5,5],order=[2,1]) 

            butcher_bD = reshape([1.0_dp, -10.0_dp/3.0_dp, 17.0_dp/3.0_dp, -10.0_dp/3.0_dp, 1.0_dp],[5]) 

            ars_eqn_check_D=.TRUE.

            diag_diff = .TRUE.
            diag_index = 2

         case ('LZ3A4a')
            
            butcher_aD = reshape([ 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                                   0.5_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                                   0.25_dp, -0.75_dp, 1.0_dp, 0.0_dp, 0.0_dp, &
                                   0.0_dp, -3.0_dp, 4.0_dp, 0.0_dp, 0.0_dp, &
                                   1.0_dp/6.0_dp, 0.0_dp, 2.0_dp/3.0_dp, -0.5_dp, 2.0_dp/3.0_dp], &
                                   [5,5],order=[2,1]) 

            butcher_bD = reshape([1.0_dp/6.0_dp, 0.0_dp, 2.0_dp/3.0_dp, -0.5_dp, 2.0_dp/3.0_dp],[5]) 

            ars_eqn_check_D=.TRUE.

            diag_diff = .TRUE.
            diag_index = 3

         case ('CFN564')
  
            butcher_aD = reshape([0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                                  0.0_dp, 0.25_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                                  0.0_dp, 0.5_dp, 0.25_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                                  0.0_dp, 17.0_dp/50.0_dp, -1.0_dp/25.0_dp, 0.25_dp, 0.0_dp, 0.0_dp, &
                                  0.0_dp, 371.0_dp/1360.0_dp, -137.0_dp/2720.0_dp, 15.0_dp/544.0_dp, 0.25_dp, 0.0_dp, &
                                  0.0_dp, 25.0_dp/24.0_dp, -49.0_dp/48.0_dp, 125.0_dp/16.0_dp, -85.0_dp/12.0_dp, 0.25_dp], &
                                  [6,6],order=[2,1]) 

            butcher_bD = reshape([0.0_dp, 25.0_dp/24.0_dp, -49.0_dp/48.0_dp, 125.0_dp/16.0_dp, -85.0_dp/12.0_dp, 0.25_dp],[6]) 
          
            ars_eqn_check_D=.TRUE.

            diag_diff = .FALSE.
            diag_index = 2

         case ('LZ4A1')
  
            butcher_aD = reshape([0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                                  -1.0_dp/6.0_dp, 0.5_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                                  1.0_dp/6.0_dp, -1.0_dp/3.0_dp, 0.5_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                                  3.0_dp/8.0_dp, -3.0_dp/8.0_dp, 0.0_dp, 0.5_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                                  1.0_dp/8.0_dp, 0.0_dp, 3.0_dp/8.0_dp, -0.5_dp, 0.5_dp, 0.0_dp, 0.0_dp, &
                                  -0.5_dp, 0.0_dp, 3.0_dp, -2.0_dp, 0.0_dp, 0.5_dp, 0.0_dp, &
                                  1.0_dp/6.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 2.0_dp/3.0_dp, -0.5_dp, 2.0_dp/3.0_dp], &
                                  [7,7],order=[2,1]) 

            butcher_bD = reshape([1.0_dp/6.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 2.0_dp/3.0_dp, -0.5_dp, 2.0_dp/3.0_dp],[7]) 
          
            ars_eqn_check_D=.TRUE.

            diag_diff = .TRUE.
            diag_index = 2

      end select            
      
   end subroutine rhs_update_wts_imp

   subroutine rhs_update_wts_exp(time_scheme_imp,time_scheme_exp,wt_rhs_tscheme_exp,n_order_tscheme_exp)

      integer, intent(in) :: n_order_tscheme_exp
      character(len=100), intent(in) :: time_scheme_imp
      character(len=100), intent(in) :: time_scheme_exp
      real(kind=dp), intent(out) :: wt_rhs_tscheme_exp(n_order_tscheme_exp)
      real(kind=dp) :: a0,b0,b1,b2,b3,r_dt1,r_dt2,r_dt3
      integer :: i
      
      select case (time_scheme_imp)

         case ('CN')

            select case (time_scheme_exp)
               case('AB2')  
                  wt_rhs_tscheme_exp(1) = 1.0_dp + half*dt_array(1)/dt_array(2)
                  wt_rhs_tscheme_exp(2) = -half*dt_array(1)/dt_array(2)

               case ('AB3')
                  wt_rhs_tscheme_exp(2) = (dt_array(1)*dt_array(3)*half + onesixth*dt_array(1)*dt_array(1)) &
                                          & /(dt_array(2)*dt_array(2)*half - dt_array(2)*dt_array(3))
                  wt_rhs_tscheme_exp(3) = (-dt_array(1)*half - dt_array(2)*wt_rhs_tscheme_exp(2))*half &
                                          & /(dt_array(3))
                  wt_rhs_tscheme_exp(1) = 1.0_dp - wt_rhs_tscheme_exp(2) - wt_rhs_tscheme_exp(3)
               
            end select

         case ('BDF2')

            select case (time_scheme_exp)
               case ('AB2')
                  r_dt1=dt_array(2)/dt_array(1)
                  a0=(2.0_dp+r_dt1)/(1.0_dp+r_dt1)
                  b0=1.0_dp+1.0_dp/r_dt1
                  b1=-1.0_dp/r_dt1
                  wt_rhs_tscheme_exp(1) = b0/a0 
                  wt_rhs_tscheme_exp(2) = b1/a0         
            
            end select
            
         case ('BDF3')

            select case (time_scheme_exp)
               case ('AB3')
                  r_dt1=dt_array(2)/dt_array(1)
                  r_dt2=dt_array(3)/dt_array(1)
                  a0=(1.0_dp + 1.0_dp/(1.0_dp+r_dt1) + 1.0_dp/(1.0_dp+r_dt1+r_dt2))
                  b0=1.0_dp*((1.0_dp + r_dt1)*(1.0_dp+r_dt1+r_dt2))/(r_dt1*(r_dt1+r_dt2))
                  b1=1.0_dp*(-(1.0_dp+r_dt1+r_dt2)/(r_dt1*r_dt2))
                  b2=1.0_dp*(1.0_dp+r_dt1)/(r_dt2*(r_dt1+r_dt2))
                  wt_rhs_tscheme_exp(1) = b0/a0 
                  wt_rhs_tscheme_exp(2) = b1/a0
                  wt_rhs_tscheme_exp(3) = b2/a0

            end select

         case ('BDF4')

            select case (time_scheme_exp)
               case ('AB4')
                  r_dt1=dt_array(1)/dt_array(2)
                  r_dt2=dt_array(2)/dt_array(3)
                  r_dt3=dt_array(3)/dt_array(4)

                  c1 = one+r_dt3*(one+r_dt2) 
                  c2 = one+r_dt2*(one+r_dt1) 
                  c3 = one+r_dt3*c2 

                  a0= one + r_dt1/(one+r_dt1)+r_dt2*r_dt1/c2+r_dt3* &
                  & r_dt2*r_dt1/c3

                  b3 = -r_dt3**3*r_dt2**2*r_dt1*(one+r_dt1)/(one+r_dt3)* & 
                      & c2/c1 
                  b2 = r_dt2**2*r_dt1*(one+r_dt1)/(one+r_dt2)*c3 
                  b1 = -c2*c3 * r_dt1/(one+r_dt3) 
                  b0 = r_dt2*(one+r_dt1)/(one+r_dt2) * ((one+r_dt1) * & 
                      & (c3+r_dt3)+(one+r_dt3)/r_dt2)/c1 
                  
                  wt_rhs_tscheme_exp(1) = b0/a0 
                  wt_rhs_tscheme_exp(2) = b1/a0
                  wt_rhs_tscheme_exp(3) = b2/a0
                  wt_rhs_tscheme_exp(4) = b3/a0

            end select

         case ('NONE')
 
            select case (time_scheme_exp)

               case ('RK2')
                 
                  butcher_a = reshape([0.0_dp, 0.0_dp, &
                                       0.5_dp, 0.0_dp], &
                                       [2,2],order=[2,1]) 

                  butcher_b = reshape([0.0_dp, 1.0_dp],[2]) 

               case ('RK4')

                  butcher_a = reshape([0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                                        0.5_dp, 0.0_dp, 0.5_dp, 0.0_dp, &
                                        0.0_dp, 0.5_dp, 0.0_dp, 0.0_dp, &
                                        0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp], &
                                        [4,4],order=[2,1])

                  butcher_b = reshape([1.0_dp/6.0_dp, 1.0_dp/3.0_dp, 1.0_dp/3.0_dp, 1.0_dp/6.0_dp],[4]) 
                  
            end select

         case ('EUL111')
 
            select case (time_scheme_exp)

               case ('expEUL111')
                 
                  butcher_aA = reshape([0.0_dp, 0.0_dp, &
                                        1.0_dp, 0.0_dp], &
                                        [2,2],order=[2,1]) 

                  butcher_bA = reshape([1.0_dp, 0.0_dp],[2]) 

             end select

         case ('ARS222')
 
            select case (time_scheme_exp)

               case ('expARS222')

                  butcher_aA = reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                                        gammaa, 0.0_dp, 0.0_dp, &
                                        delta, 1.0_dp-delta, 0.0_dp], &
                                        [3,3],order=[2,1])  

                  butcher_bA = reshape([delta, 1.0_dp-delta, 0.0_dp],[3]) 

                  do i=1,n_order_tscheme_exp
                     if (butcher_aA(n_order_tscheme_exp,i)==butcher_bA(i)) then
                        ars_eqn_check_A=.TRUE.
                     end if
                  end do

            end select

         case ('ARS232')
 
            select case (time_scheme_exp)

               case ('expARS232')

                  butcher_aA = reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                                        gammaa, 0.0_dp, 0.0_dp, &
                                        delta, 1.0_dp-delta, 0.0_dp], &
                                        [3,3],order=[2,1])  

                  butcher_bA = reshape([0.0_dp, 1.0_dp-gammaa, gammaa],[3]) 

                  do i=1,n_order_tscheme_exp
                     if (butcher_aA(n_order_tscheme_exp,i)==butcher_bA(i)) then
                        ars_eqn_check_A=.TRUE.
                     end if
                  end do

            end select 

         case ('ARS233')
 
            select case (time_scheme_exp)

               case ('expARS233')

                  butcher_aA = reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                                        gammaa, 0.0_dp, 0.0_dp, &
                                        gammaa-1.0_dp, 2.0_dp*(1.0_dp-gammaa), 0.0_dp], &
                                        [3,3],order=[2,1])  

                  butcher_bA = reshape([0.0_dp, 0.5_dp, 0.5_dp],[3]) 

                  do i=1,n_order_tscheme_exp
                     if (butcher_aA(n_order_tscheme_exp,i)==butcher_bA(i)) then
                        ars_eqn_check_A=.TRUE.
                     end if
                  end do

            end select 

         case ('ARS343')
 
            select case (time_scheme_exp)

               case ('expARS343')

                  butcher_aA = reshape([0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                                        0.4358665215_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                                        0.321278886_dp, 0.3966543747_dp, 0.0_dp, 0.0_dp, &
                                        -0.105858296_dp, 0.5529291479_dp, 0.5529291479_dp, 0.0_dp], &
                                        [4,4],order=[2,1])  

                  butcher_bA = reshape([0.0_dp, 1.208496649_dp, -0.644363171_dp, 0.4358665215_dp],[4]) 

                  do i=1,n_order_tscheme_exp
                     if (butcher_aA(n_order_tscheme_exp,i)==butcher_bA(i)) then
                        ars_eqn_check_A=.TRUE.
                     end if
                  end do

            end select  

         case ('MARS343')
 
            select case (time_scheme_exp)

               case ('expMARS343')

                  butcher_aA = reshape([0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                                        0.435866521508458_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                                        0.535396540307354_dp, 0.182536720446875_dp, 0.0_dp, 0.0_dp, &
                                        0.63041255815287_dp, -0.83193390106308_dp, 1.20152134291021_dp, 0.0_dp], &
                                        [4,4],order=[2,1])  

                  butcher_bA = reshape([0.0_dp, 1.20849664917601_dp, -0.64436317068446_dp, 0.435866521508458_dp],[4]) 

                  ars_eqn_check_A=.FALSE.

            end select  

         case ('IJ3')
 
            select case (time_scheme_exp)

               case ('expIJ3')

                  butcher_aA = reshape([0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                                        0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                                        0.0_dp, onethird, 0.0_dp, 0.0_dp, &
                                        0.0_dp, -1.0_dp, 2.0_dp, 0.0_dp], &
                                        [4,4],order=[2,1])  

                  butcher_bA = reshape([0.0_dp, 0.0_dp, 0.75_dp, 0.25_dp],[4]) 

                  ars_eqn_check_A=.FALSE.

            end select  

          case ('PC2')
 
            select case (time_scheme_exp)

               case ('expPC2')

                  butcher_aA = reshape([0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                                        1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                                        0.5_dp, 0.5_dp, 0.0_dp, 0.0_dp, &
                                        0.5_dp, 0.0_dp, 0.5_dp, 0.0_dp], &
                                        [4,4],order=[2,1])  

                  butcher_bA = reshape([0.5_dp, 0.0_dp, 0.5_dp, 0.0_dp],[4]) 

                  ars_eqn_check_A=.TRUE.

            end select  

          case ('SSP332')
 
            select case (time_scheme_exp)

               case ('expSSP332')

                  butcher_aA = reshape([0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                                        0.5_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                                        0.5_dp, 0.5_dp, 0.0_dp, 0.0_dp, &
                                        1.0_dp/3.0_dp, 1.0_dp/3.0_dp, 1.0_dp/3.0_dp, 0.0_dp], &
                                        [4,4],order=[2,1])  

                  butcher_bA = reshape([1.0_dp/3.0_dp, 1.0_dp/3.0_dp, 1.0_dp/3.0_dp, 0.0_dp],[4]) 

                  ars_eqn_check_A=.TRUE.

            end select  
 
          case ('SSP543')
 
            select case (time_scheme_exp)

               case ('expSSP543')

                  butcher_aA = reshape([0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                                  0.5_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                                  0.5_dp, 0.5_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                                  1.0_dp/6.0_dp, 1.0_dp/6.0_dp, 1.0_dp/6.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                                  0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                                  1.0_dp/6.0_dp, 1.0_dp/6.0_dp, 1.0_dp/6.0_dp, 0.5_dp, 0.0_dp, 0.0_dp], &
                                  [6,6],order=[2,1]) 

                  butcher_bA = reshape([1.0_dp/6.0_dp, 1.0_dp/6.0_dp, 1.0_dp/6.0_dp, 0.5_dp, 0.0_dp, 0.0_dp],[6]) 

                  ars_eqn_check_A=.TRUE.

            end select  
 

         case ('ARS443')
 
            select case (time_scheme_exp)

               case ('expARS443')

                  butcher_aA = reshape([ 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                                         0.5_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                                         11.0_dp/18.0_dp,1.0_dp/18.0_dp, 0.0_dp, 0.0_dp,0.0_dp, &
                                         5.0_dp/6.0_dp, -5.0_dp/6.0_dp, 0.5_dp, 0.0_dp, 0.0_dp, &
                                         0.25_dp, 7.0_dp/4.0_dp, 0.75_dp, -7.0_dp/4.0_dp, 0.0_dp], &
                                         [5,5],order=[2,1]) 

                  butcher_bA = reshape([0.25_dp, 7.0_dp/4.0_dp, 0.75_dp, -7.0_dp/4.0_dp, 0.0_dp],[5]) 

                  do i=1,n_order_tscheme_exp
                     if (butcher_aA(n_order_tscheme_exp,i)==butcher_bA(i)) then
                        ars_eqn_check_A=.TRUE.
                     end if
                  end do

            end select  

         case ('PR433')
 
            select case (time_scheme_exp)

               case ('expPR433')

                  !butcher_aA = reshape([ 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                  !                       0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                  !                       0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                  !                       0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                  !                       0.0_dp, 0.25_dp, 0.25_dp, 0.0_dp, 0.0_dp], &
                  !                       [5,5],order=[2,1]) 
                  butcher_aA = reshape([ 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                                         0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                                         0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                                         0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, &
                                         0.0_dp, 0.0_dp, 0.25_dp, 0.25_dp, 0.0_dp], &
                                         [5,5],order=[2,1]) 
                  !butcher_aA = reshape([ 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                  !                       0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                  !                       0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, &
                  !                       0.0_dp, 0.25_dp, 0.25_dp, 0.0_dp], &
                  !                       [4,4],order=[2,1]) 

                  !butcher_bA = reshape([0.0_dp, 1.0_dp/6.0_dp, 1.0_dp/6.0_dp, 2.0_dp/3.0_dp],[4]) 
                  butcher_bA = reshape([0.0_dp,0.0_dp, 1.0_dp/6.0_dp, 1.0_dp/6.0_dp, 2.0_dp/3.0_dp],[5]) 
                  !butcher_bA = reshape([0.0_dp, 1.0_dp/6.0_dp, 1.0_dp/6.0_dp, 2.0_dp/3.0_dp,0.0_dp],[5]) 

                  ars_eqn_check_A=.FALSE.

            end select  

         case ('BPR353')
 
            select case (time_scheme_exp)

               case ('expBPR353')

                  butcher_aA = reshape([0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                                        1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                                        4.0_dp/9.0_dp, 2.0_dp/9.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                                        0.25_dp, 0.0_dp, 0.75_dp, 0.0_dp, 0.0_dp, &
                                        0.25_dp, 0.0_dp, 0.75_dp, 0.0_dp, 0.0_dp], &
                                        [5,5],order=[2,1])

                  butcher_bA = reshape([0.25_dp, 0.0_dp, 0.75_dp, 0.0_dp, 0.0_dp],[5]) 

                  do i=1,n_order_tscheme_exp
                     if (butcher_aA(n_order_tscheme_exp,i)==butcher_bA(i)) then
                        ars_eqn_check_A=.TRUE.
                     end if
                  end do

            end select 

         case ('KC343')
 
            select case (time_scheme_exp)

               case ('expKC343')
                  
                  butcher_aA = reshape([0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                                  1767732205903.0_dp/2027836641118.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
               5535828885825.0_dp/10492691773637.0_dp, 788022342437.0_dp/10882634858940.0_dp,0.0_dp, 0.0_dp, &
               6485989280629.0_dp/16251701735622.0_dp,-4246266847089.0_dp/9704473918619.0_dp, &
              & 10755448449292.0_dp/10357097424841.0_dp, 0.0_dp], &
                                  [4,4],order=[2,1]) 

            butcher_bA = reshape([1471266399579.0_dp/7840856788654.0_dp, -4482444167858.0_dp/7529755066697.0_dp, &
            & 11266239266428.0_dp/11593286722821.0_dp,1767732205903.0_dp/4055673282236.0_dp],[4]) 
            !butcher_bA = reshape([2756255671327.0_dp/12835298489170.0_dp, -10771552573575.0_dp/22201958757719.0_dp, &
            !& 9247589265047.0_dp/10645013368117.0_dp,2193209047091.0_dp/5459859503100.0_dp],[4]) 
            !butcher_bA = reshape([4655552711362.0_dp/22874653954995.0_dp, -18682724506714.0_dp/9892148508045.0_dp, &
            !& 34259539580243.0_dp/13192909600954.0_dp,584795268549.0_dp/6622622206610.0_dp],[4]) 
            !      butcher_bA = reshape([0.187641024346724_dp, -0.595297473576955_dp, 0.971789927721772_dp, &
            !                            & 0.435866521508459_dp],[4]) 

                  do i=1,n_order_tscheme_exp
                     if (butcher_aA(n_order_tscheme_exp,i)==butcher_bA(i)) then
                        ars_eqn_check_A=.TRUE.
                     end if
                  end do

            end select  

         case ('KC564')
 
            select case (time_scheme_exp)

               case ('expKC564')

                  butcher_aA = reshape([ 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                                         0.5_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                                         13861.0_dp/62500.0_dp,6889.0_dp/62500.0_dp, 0.0_dp, 0.0_dp,0.0_dp, 0.0_dp, &
                  -116923316275.0_dp/2393684061468.0_dp, -2731218467317.0_dp/15368042101831.0_dp,  & 
                  & 9408046702089._dp/11113171139209.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                  -451086348788.0_dp/2902428689909.0_dp, -2682348792572.0_dp/7519795681897.0_dp, & 
                  & 12662868775082.0_dp/11960479115383.0_dp, 3355817975965.0_dp/11060851509271.0_dp, 0.0_dp, 0.0_dp, &
      647845179188.0_dp/3216320057751.0_dp, 73281519250.0_dp/8382639484533.0_dp, 552539513391.0_dp/3454668386233.0_dp, &
      & 3354512671639.0_dp/8306763924573.0_dp, 4040.0_dp/17871.0_dp, 0.0_dp], &
                                         [6,6],order=[2,1]) 

                  butcher_bA = reshape([82889.0_dp/524892.0_dp, 0.0_dp, 15625.0_dp/83664.0_dp, 69875.0_dp/102672.0_dp, & 
                                  & -2260.0_dp/8211.0_dp, 0.25_dp],[6]) 

                  !do i=1,n_order_tscheme_exp
                  !   if (butcher_aA(n_order_tscheme_exp,i)==butcher_bA(i)) then
                  !      ars_eqn_check_A=.TRUE.
                        ars_eqn_check_A=.FALSE.
                  !   end if
                  !end do

            end select   
         
         case ('FW53')
 
            select case (time_scheme_exp)

               case ('expFW53')
                  butcher_aA = reshape([ 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                                   1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                                   0.5_dp, 0.5_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                                   0.25_dp, 0.25_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                                   1.0_dp/6.0_dp, 1.0_dp/6.0_dp, 0.0_dp, 4.0_dp/6.0_dp, 0.0_dp], &
                                   [5,5],order=[2,1]) 

                  butcher_bA = reshape([1.0_dp/6.0_dp, 1.0_dp/6.0_dp, 0.0_dp, 4.0_dp/6.0_dp, 0.0_dp],[5]) 

                  do i=1,n_order_tscheme_exp
                     if (butcher_aA(n_order_tscheme_exp,i)==butcher_bA(i)) then
                        ars_eqn_check_A=.TRUE.
                     end if
                  end do

            end select    

         case ('CK222')
 
            select case (time_scheme_exp)

               case ('expCK222')

                  butcher_aA = reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                                        twothird, 0.0_dp, 0.0_dp, &
                                        0.25_dp, 0.75_dp, 0.0_dp], &
                                        [3,3],order=[2,1])  

                  butcher_bA = reshape([0.25_dp, 0.75_dp, 0.0_dp],[3]) 

                  do i=1,n_order_tscheme_exp
                     if (butcher_aA(n_order_tscheme_exp,i)==butcher_bA(i)) then
                        ars_eqn_check_A=.TRUE.
                     end if
                  end do

            end select

         case ('BPR442')
 
            select case (time_scheme_exp)

               case ('expBPR442')

                  butcher_aA = reshape([0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                                        0.25_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                                        13.0_dp/4.0_dp, -3.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                                        0.25_dp, 0.0_dp, 0.5_dp, 0.0_dp, 0.0_dp, &
                                        0.0_dp, 1.0_dp/3.0_dp, 1.0_dp/6.0_dp, 0.5_dp, 0.0_dp], &
                                        [5,5],order=[2,1])

                  butcher_bA = reshape([0.0_dp, 1.0_dp/3.0_dp, 1.0_dp/6.0_dp, 0.5_dp, 0.0_dp],[5]) 

                  do i=1,n_order_tscheme_exp
                     if (butcher_aA(n_order_tscheme_exp,i)==butcher_bA(i)) then
                        ars_eqn_check_A=.TRUE.
                     end if
                  end do

            end select 

          case ('CFN343')
 
            select case (time_scheme_exp)

               case ('expCFN343')

                  butcher_aA = reshape([0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                                        gammaa, 0.0_dp, 0.0_dp, 0.0_dp, &
                                        (1.0_dp+gammaa)/2.0_dp+0.35_dp, -0.35_dp, 0.0_dp, 0.0_dp, &
                                        0.0_dp, 1.0_dp+0.989175724679846_dp, -0.989175724679846_dp, 0.0_dp], &
                                        [4,4],order=[2,1])  

                  butcher_bA = reshape([0.0_dp, 1.208496649176010_dp, -0.644363170684469_dp, gammaa],[4]) 

                  ars_eqn_check_A=.FALSE.

            end select              

         case ('BHR553')
 
            select case (time_scheme_exp)

               case ('expBHR553')

                  butcher_aA = reshape([0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                                        2.0_dp*gammaa, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                                        aa31, aa32, 0.0_dp, 0.0_dp, 0.0_dp, &
                                        aa41, aa42, aa43, 0.0_dp, 0.0_dp, &
                                        aa51, aa52, aa53, aa54, 0.0_dp], &
                                        [5,5],order=[2,1])

                  butcher_bA = reshape([bb1, 0.0_dp, bb3, bb4, gammaa],[5]) 

                  ars_eqn_check_A=.FALSE.

            end select 

         case ('CB3e')
 
            select case (time_scheme_exp)

               case ('expCB3e')

                  butcher_aA = reshape([0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                                        1.0_dp/3.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                                        0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, &
                                        0.0_dp, 0.75_dp, 0.25_dp, 0.0_dp], &
                                        [4,4],order=[2,1])  

                  butcher_bA = reshape([0.0_dp, 0.75_dp, -0.25_dp, 0.5_dp],[4]) 

                  ars_eqn_check_A=.FALSE.

            end select              

         case ('CB3f')
 
            select case (time_scheme_exp)

               case ('expCB3f')

                  butcher_aA = reshape([0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                                        aa21, 0.0_dp, 0.0_dp, 0.0_dp, &
                                        aa31, aa32, 0.0_dp, 0.0_dp, &
                                        bb1, aa42, aa43, 0.0_dp], &
                                        [4,4],order=[2,1])  

                  butcher_bA = reshape([bb1, bb2, bb3, bb4],[4]) 

                  ars_eqn_check_A=.FALSE.

            end select              

         case ('CB4')
 
            select case (time_scheme_exp)

               case ('expCB4')
                  butcher_aA = reshape([0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                                       aa21, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                                       aa31, aa32, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                                       bb1, aa42, aa43, 0.0_dp, 0.0_dp, 0.0_dp, &
                                       bb1, bb2, aa53, aa54, 0.0_dp, 0.0_dp, &
                                       bb1, bb2, bb3, aa64, aa65, 0.0_dp], &
                                       [6,6],order=[2,1]) 

                  butcher_bA = reshape([bb1, bb2, bb3, bb4, bb5, bb6],[6]) 
 
                  ars_eqn_check_A=.FALSE.

            end select              

         case ('LZ3L1')
 
            select case (time_scheme_exp)

               case ('expLZ3L1')

                  butcher_aA = reshape([0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                                        0.25_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                                        0.0_dp, 0.5_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                                        -0.5_dp, 5.0_dp/4.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                                        0.0_dp, 2.0_dp/3.0_dp, -1.0_dp/3.0_dp, 2.0_dp/3.0_dp, 0.0_dp], &
                                        [5,5],order=[2,1])

                  butcher_bA = reshape([0.0_dp, 2.0_dp/3.0_dp, -1.0_dp/3.0_dp, 2.0_dp/3.0_dp, 0.0_dp],[5]) 

                  ars_eqn_check_A=.TRUE.

            end select 

         case ('LZ3A4a')
 
            select case (time_scheme_exp)

               case ('expLZ3A4a')

                  butcher_aA = reshape([0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                                        0.5_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                                        0.25_dp, 0.25_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                                        0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                                        1.0_dp/6.0_dp, 0.0_dp, 2.0_dp/3.0_dp, 1.0_dp/6.0_dp, 0.0_dp], &
                                        [5,5],order=[2,1])

                  butcher_bA = reshape([1.0_dp/6.0_dp, 0.0_dp, 2.0_dp/3.0_dp, 1.0_dp/6.0_dp, 0.0_dp],[5]) 

                  ars_eqn_check_A=.TRUE.

            end select 



         case ('CFN564')
 
            select case (time_scheme_exp)

               case ('expCFN564')

                  butcher_aA = reshape([0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                                  0.25_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                                  -0.25_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                                  -13.0_dp/100.0_dp, 43.0_dp/75.0_dp, 8.0_dp/75.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                                  -6.0_dp/85.0_dp, 42.0_dp/85.0_dp, 179.0_dp/1360.0_dp, -15.0_dp/272.0_dp, 0.0_dp, 0.0_dp, &
                                  0.0_dp, 79.0_dp/24.0_dp, -5.0_dp/8.0_dp, 25.0_dp/2.0_dp, -85.0_dp/6.0_dp, 0.0_dp], &
                                  [6,6],order=[2,1]) 

                  butcher_bA = reshape([0.0_dp, 25.0_dp/24.0_dp, -49.0_dp/48.0_dp, 125.0_dp/16.0_dp, -85.0_dp/12.0_dp, 0.25_dp],[6]) 

                  ars_eqn_check_A=.FALSE.

            end select  

         case ('LZ4A1')
 
            select case (time_scheme_exp)

               case ('expLZ4A1')

                  butcher_aA = reshape([0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                                        1.0_dp/3.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                                        1.0_dp/6.0_dp, 1.0_dp/6.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                                        1.0_dp/8.0_dp, 0.0_dp, 3.0_dp/8.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                                        1.0_dp/8.0_dp, 0.0_dp, 3.0_dp/8.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                                        0.5_dp, 0.0_dp, -1.5_dp, 1.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, &
                                        1.0_dp/6.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 2.0_dp/3.0_dp, 1.0_dp/6.0_dp, 0.0_dp], &
                                        [7,7],order=[2,1]) 

                  butcher_bA = reshape([1.0_dp/6.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 2.0_dp/3.0_dp, 1.0_dp/6.0_dp, 0.0_dp],[7]) 

                  ars_eqn_check_A=.TRUE.

            end select  

      end select 

   end subroutine rhs_update_wts_exp

end module timeschemes
