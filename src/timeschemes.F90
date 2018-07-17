module timeschemes

   use double
   use constants, only: onesixth, onethird, onefourth, threefourth, threehalf, &
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
   logical, public :: ars_eqn_check_A, ars_eqn_check_D
   real(kind=dp) :: gammaa, delta, alphaa, betaa, etaa

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

         case ('PC2') ! Ascher, Ruuth and Spiteri (ARS343) 3rd order, Appl. Numer. Math. 1997 paper
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
            n_order_tscheme_imp=4 
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
                  n_order_tscheme_exp=4
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

         case ('CFN343') ! Calvo 3nd order, 2001 paper 
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
      real(kind=dp) :: a0,a1,a2,a3,a4,r_dt1,r_dt2
      integer :: i,j

      select case (time_scheme_imp)

         case ('CN') 
            wt_rhs_tscheme_imp(1) = half
            wt_lhs_tscheme_imp=half

         case ('BDF2')
            r_dt1=dt_array(2)/dt_array(1)
            a0=(2.0_dp+r_dt1)/(1.0_dp+r_dt1)
            a1=1.0_dp+1.0_dp/r_dt1
            a2=-1.0_dp*(1.0_dp/(r_dt1)/(1.0_dp+r_dt1)) ! Thomas found mistake here to fix on Sept 20th 2017
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
            a0=25.0_dp/12.0_dp
            a1=4.0_dp
            a2=-3.0_dp
            a3=4.0_dp/3.0_dp
            a4=-1.0_dp/4.0_dp
            wt_lhs_tscheme_imp=1.0_dp/a0
            wt_rhs_tscheme_imp(1) = a1/a0
            wt_rhs_tscheme_imp(2) = a2/a0
            wt_rhs_tscheme_imp(3) = a3/a0
            wt_rhs_tscheme_imp(4) = a4/a0

         case ('EUL111')

            butcher_aD = reshape([0.0_dp, 0.0_dp, &
                                  0.0_dp, 1.0_dp], &
                                  [2,2],order=[2,1]) 

            butcher_bD = reshape([0.0_dp, 1.0_dp],[2]) 

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

         case ('PC2')

            butcher_aD = reshape([0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                                  0.5_dp, 0.5_dp, 0.0_dp, 0.0_dp, &
                                  0.5_dp, 0.0_dp, 0.5_dp, 0.0_dp, &
                                  0.5_dp, 0.0_dp, 0.0_dp, 0.5_dp], &
                                  [4,4],order=[2,1])  

            butcher_bD = reshape([0.5_dp, 0.0_dp, 0.0_dp, 0.5_dp],[4]) 

            do i=1,n_order_tscheme_imp
               if (butcher_aD(n_order_tscheme_imp,i)==butcher_bD(i)) then
                  ars_eqn_check_D=.TRUE.
               end if
            end do   

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

         case ('PR433')
            
            !butcher_aD = reshape([ 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
            !                       0.0_dp, alphaa, 0.0_dp, 0.0_dp, 0.0_dp, &
            !                       0.0_dp,-alphaa, alphaa, 0.0_dp,0.0_dp, &
            !                       0.0_dp, 0.0_dp, 1.0_dp-alphaa, alphaa, 0.0_dp, &
            !                       0.0_dp, betaa, etaa, 0.5_dp-betaa-etaa-alphaa, alphaa], &
            !                       [5,5],order=[2,1]) ! padded version

            butcher_aD = reshape([ alphaa, 0.0_dp, 0.0_dp, 0.0_dp, &
                                   -alphaa, alphaa, 0.0_dp,0.0_dp, &
                                   0.0_dp, 1.0_dp-alphaa, alphaa, 0.0_dp, &
                                   betaa, etaa, 0.5_dp-betaa-etaa-alphaa, alphaa], &
                                   [4,4],order=[2,1]) 

            butcher_bD = reshape([0.0_dp, 1.0_dp/6.0_dp, 1.0_dp/6.0_dp, 2.0_dp/3.0_dp],[4]) 

            ars_eqn_check_D=.FALSE.

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

         case ('BPR442')
            
            butcher_aD = reshape([ 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                                   0.0_dp, 0.25_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                                   0.0_dp, 0.0_dp, 0.25_dp, 0.0_dp, 0.0_dp, &
                                   0.0_dp, 1.0_dp/24.0_dp, 11.0_dp/24.0_dp, 0.25_dp, 0.0_dp, &
                                   0.0_dp, 11.0_dp/24.0_dp, 1.0_dp/6.0_dp, 1.0_dp/8.0_dp, 0.25_dp], &
                                   [5,5],order=[2,1]) 

            butcher_bD = reshape([0.0_dp, 11.0_dp/24.0_dp, 1.0_dp/6.0_dp, 1.0_dp/8.0_dp, 0.25_dp],[5]) 

            do i=1,n_order_tscheme_imp
               if (butcher_aD(n_order_tscheme_imp,i)==butcher_bD(i)) then
                  ars_eqn_check_D=.TRUE.
               end if
            end do

         case ('CFN343')

            butcher_aD = reshape([0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                                  0.0_dp, gammaa, 0.0_dp, 0.0_dp, &
                                  0.0_dp, (1.0_dp-gammaa)/2.0_dp, gammaa, 0.0_dp, &
                                  0.0_dp, 1.208496649176010_dp, -0.644363170684469_dp, gammaa], &
                                  [4,4],order=[2,1])  

            butcher_bD = reshape([0.0_dp, 1.208496649176010_dp, -0.644363170684469_dp, gammaa],[4]) 

            do i=1,n_order_tscheme_imp
               if (butcher_aD(n_order_tscheme_imp,i)==butcher_bD(i)) then
                  ars_eqn_check_D=.TRUE.
               end if
            end do   


      end select            
      
   end subroutine rhs_update_wts_imp

   subroutine rhs_update_wts_exp(time_scheme_imp,time_scheme_exp,wt_rhs_tscheme_exp,n_order_tscheme_exp)

      integer, intent(in) :: n_order_tscheme_exp
      character(len=100), intent(in) :: time_scheme_imp
      character(len=100), intent(in) :: time_scheme_exp
      real(kind=dp), intent(out) :: wt_rhs_tscheme_exp(n_order_tscheme_exp)
      real(kind=dp) :: a0,b0,b1,b2,b3,r_dt1,r_dt2
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
                  a0=25.0_dp/12.0_dp
                  b0=4.0_dp
                  b1=-6.0_dp
                  b2=4.0_dp
                  b3=-1.0_dp
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

          case ('PC2')
 
            select case (time_scheme_exp)

               case ('expPC2')

                  butcher_aA = reshape([0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                                        1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                                        0.5_dp, 0.5_dp, 0.0_dp, 0.0_dp, &
                                        0.5_dp, 0.0_dp, 0.5_dp, 0.0_dp], &
                                        [4,4],order=[2,1])  

                  butcher_bA = reshape([0.5_dp, 0.0_dp, 0.5_dp, 0.0_dp],[4]) 

                  do i=1,n_order_tscheme_exp
                     if (butcher_aA(n_order_tscheme_exp,i)==butcher_bA(i)) then
                        ars_eqn_check_A=.TRUE.
                     end if
                  end do

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

                 ! butcher_aA = reshape([ 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                 !                        0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                 !                        0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                 !                        0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, &
                 !                        0.0_dp, 0.0_dp, 0.25_dp, 0.25_dp, 0.0_dp], &
                 !                        [5,5],order=[2,1]) ! padded version 
                  butcher_aA = reshape([ 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                                         0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                                         0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, &
                                         0.0_dp, 0.25_dp, 0.25_dp, 0.0_dp], &
                                         [4,4],order=[2,1]) 

                  butcher_bA = reshape([0.0_dp, 1.0_dp/6.0_dp, 1.0_dp/6.0_dp, 2.0_dp/3.0_dp],[4]) 

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

      end select 

   end subroutine rhs_update_wts_exp

end module timeschemes
