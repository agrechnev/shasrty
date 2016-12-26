! Global data and routines for triangle
module ttt
 implicit none

 real(8), parameter :: pi=3.14159265358979_8

 !--- Global data----
 ! Input data
 real(8) :: J1,J2,Bu,H  ! Hamiltonian parameters

 integer :: mode_alpha  !What to look for?
 ! 0 = read alpha1,2
 ! .ne. 0 = automatic alpha1,2
 ! 1 = Neel : alpha1,2=1/2
 ! 2 = Flip (FM) : alpha1,2=1/2 (alphas do not matter)
 ! 3 = Dimer : alpha1,2=0
 ! 4 = UUD  : alpha1,2= (2-rho+hh)/(3*hh)
 ! 5 = Spin Flop : alpha1,2=1/2
 ! 6 = Umbrella : alpha1,2=1/(rho+1) Is this right ?
 !
 ! mode_alpha is also used to print the exact result if possible 
 
 real(8) :: alpha1, alpha2

 integer :: maxit  ! Maximum number of iterations
 
 real(8) :: lambda ! The relaxation parameter

 real(8) :: tol_g, tol_e ! The desired accuracy of the gradients and energy/site

 integer :: nshuffle ! # of shuffles

 
 real(8) :: rho ! rho=J2/J1
 real(8) :: beta ! beta=Bu/J1
 real(8) :: hh ! hh=H/J1


 !----- The spin configuration
 
 ! Angles phi, theta (x,y,z)
 real(8), dimension(1:3) :: phi, theta

 ! The "gradient" dE/d phi , dE/d theta

 real(8), dimension(1:3) :: de_dphi,de_dtheta

 real(8) :: maxgrad_t, maxgrad_p ! The maximum component of each grad

 ! The "best" values for the shuffle algorithm
 real(8), dimension(1:3) :: best_phi, best_theta
 real(8) :: best_e, best_m ! Energy and magnetic moment

 integer :: best_shuf ! Shuffle iteration number


 ! The total energy

 real(8) :: energy

 real(8) :: energy_old, energy_diff ! Comparison to the previous energy

 ! The magnetic moment
 real(8) :: mag_mom

 !-------------------------------------------------------------
contains

 !------------------------
 ! Seed the random number generator

 
subroutine init_random_seed()
 ! Example from GNU manual, modified
 ! Does not work well in date_and_time mode somehow
            implicit none
            integer, allocatable :: seed(:)
            integer :: i, n, un, istat, dt(8), pid, t(2), s
            integer(8) :: count, tms
          
            call random_seed(size = n)
            allocate(seed(n))
            ! First try if the OS provides a random number generator
            open(21, file="/dev/urandom", access="stream", &
                 form="unformatted", action="read", status="old", iostat=istat)
            if (istat == 0) then
               write (6,*) "Seed from OS"
               read(21) seed
               close(21)         
            else
               ! Fallback to XOR:ing the current time and pid. The PID is
               ! useful in case one launches multiple instances of the same
               ! program in parallel.
               call system_clock(count)
               if (count /= 0) then
                  write (6,*) "Seed from system_clock"
                  t = transfer(count, t)
               else
                  write (6,*) "Seed from date_and_time"
                  call date_and_time(values=dt)
                  tms = (dt(1) - 1970) * 365_8 * 24 * 60 * 60 * 1000 &
                       + dt(2) * 31_8 * 24 * 60 * 60 * 1000 &
                       + dt(3) * 24 * 60 * 60 * 60 * 1000 &
                       + dt(5) * 60 * 60 * 1000 &
                       + dt(6) * 60 * 1000 + dt(7) * 1000 &
                       + dt(8)
                  t = transfer(tms, t)
               end if
               s = ieor(t(1), t(2))
  
               ! getpid is not portable
               !pid = getpid() + 1099279 ! Add a prime
               pid = 2 + 1099279 ! Add a prime
               s = ieor(s, pid)
               if (n >= 3) then
                  seed(1) = t(1) + 36269
                  seed(2) = t(2) + 72551
                  seed(3) = pid
                  if (n > 3) then
                     seed(4:) = s + 37 * (/ (i, i = 0, n - 4) /)
                  end if
               else
                  seed = s + 37 * (/ (i, i = 0, n - 1 ) /)
               end if
            end if
            call random_seed(put=seed)
          end subroutine init_random_seed

 ! --- And my one ---

 subroutine ttt_seed
  ! Seed the random number generator
  ! Funny, this works beautifully, unlike the example above ;)

  integer :: car(1:8),seed(1:8)
  integer :: i

  call date_and_time(VALUES=car)

  do i=1,8
   seed(i)=car(8)+i
  end do

  !write (6,*) "system clock: car=",car
  !write (6,*) "seed=",seed
  
  call random_seed(put=seed)
 end subroutine ttt_seed
 
 !------------------------
 ! Randomize the spin configuration
 subroutine ttt_random

  integer :: i ! Loop counter over 3 spins

  real(8) :: p0,t0 ! Random phi, theta

  real(8) :: harv ! The random number

  do i=1,3
   call random_number(harv)
   p0=2*pi*harv   
   call random_number(harv)
   t0=acos(2.0_8*harv-1.0_8)  ! Proper distribution of angles
   
   phi(i)=p0
   theta(i)=t0
  end do ! i

 end subroutine ttt_random


!--------------------------------
! Calculate the contribution of one isotropic exchange bond

 subroutine ex_is(J,n1,n2,gp,gt,e)
  real(8), intent(in) :: J  ! The exchange interaction
  integer, intent(in) :: n1,n2  ! The 2 nodes
  real(8), intent(inout) :: gp,gt,e ! The "gradient", energy

  ! Local data
  real(8) :: p1,p2,t1,t2 ! The two theta, phi
  real(8) :: gp0,gt0,e0 ! The contributions

  real(8) :: a1,a2,b1,b2,c,d

  !--- Code
  ! The 2 spins

  t1=theta(n1)
  p1=phi(n1)
  t2=theta(n2)
  p2=phi(n2)

  ! The trig functions
  a1=sin(t1)
  a2=sin(t2)
  b1=cos(t1)
  b2=cos(t2)
  c=sin(p1-p2)
  d=cos(p1-p2)

  ! The dE/dphi1, dE/dtheta1, E
  gp0=-a1*a2*c
  gt0=b1*a2*d-a1*b2
  e0=a1*a2*d+b1*b2
 
  ! Add the data
  gp=gp+J*gp0
  gt=gt+J*gt0
  e=e+J*e0/2 ! 1/2 to avoid double counting

end subroutine ex_is


!------------------
! Calculate gradient and energy for a given spin config

 subroutine ttt_grad(nstep)
  integer, intent (in) :: nstep
   
  integer :: i

  ! The contributions to gradient and energy
  real(8) :: dgp, dgt, de

  ! Misc
  real(8) ::a,b
  
  ! ---- Start the code ---
  energy =0
  mag_mom=0
  maxgrad_t=0
  maxgrad_p=0
    
  do i=1,3
   
     ! Contributions for each node
     dgp=0
     dgt=0
     de=0

     ! Add the neighbors
     ! Note J2/2 instead of J2 for the triangle
     
     select case (i) ! J2 bond is 1-3, the rest is J1
      case (1)
       call ex_is(J1,1,2,dgp,dgt,de)
       call ex_is(J2/2,1,3,dgp,dgt,de)
      
      case (2)
       call ex_is(J1,2,1,dgp,dgt,de)
       call ex_is(J1,2,3,dgp,dgt,de)
    
      case (3)
       call ex_is(J1,3,2,dgp,dgt,de)
       call ex_is(J2/2,3,1,dgp,dgt,de)

     end select

     ! Add the Zeeman term and anisotropy term
     
     a=sin(theta(i))
     b=cos(theta(i))
     if (i==2) then ! Use alpha1,alpha2

      de  = de  - H*alpha1*b + Bu*alpha2*a*a/2
      dgt = dgt + H*alpha1*a + Bu*alpha2*a*b

     else ! use (1-alpha1,2)/2

      de  = de  - H*((1-alpha1)/2)*b+ Bu*((1-alpha2)/2)*a*a/2
      dgt = dgt + H*((1-alpha1)/2)*a+ Bu*((1-alpha2)/2)*a*b

     end if ! i==1
    
     ! Store this data
     
     mag_mom=mag_mom+cos(theta(i))
     de_dphi(i)=dgp
     de_dtheta(i)=dgt
     energy=energy+de
     
     ! Max gradient
     if (abs(dgp)>maxgrad_p) maxgrad_p=abs(dgp)
     if (abs(dgt)>maxgrad_t) maxgrad_t=abs(dgt)
  end do  !i

  ! Normalize magnetic moment per one spin 
  mag_mom = mag_mom/3

  if (nstep>1) then ! Calculate the energy difference
   energy_diff=energy - energy_old
   !if (energy_diff > 0) then
    !lambda=lambda*0.9_8  ! Cut the step if the energy INCREASED in the last step
    !write (6,*) "ttt_grad: lambda decreased to ",lambda
   !end if
  end if

  energy_old=energy ! Store the energy
  
 end subroutine ttt_grad

!-----------------
! Make one step

subroutine ttt_step

  ! local data
  integer :: i
  real(8) :: delta_p, delta_t, a, t,p
 

  ! code
  
  ! Change theta, phi based on gradients, lambda

  do i=1,3
   a=sin(theta(i))

   delta_t = - lambda*de_dtheta(i)
   
   if (abs(a) < 1.e-15) then  ! Safety cutoff for theta=0,pi case to avoid division by 0
    delta_p=0
   else
    delta_p = - lambda*de_dphi(i)/(a*a)
   end if

   ! Rotate the spins
   t = theta(i) + delta_t
   p = phi(i) + delta_p

   ! Enforce the limits on theta
   if (t < 0) then
    t= - t
    p = p + pi
   end if

   if (t > pi) then
    t= 2* pi - t
    p = p + pi
   end if
  
   ! Enforce the limits on phi

   if (p < 0) p=p+2*pi
   if (p> 2*pi) p=p-2*pi

   ! Put the spins back
   theta(i) = t
   phi(i) = p

  end do ! i

end subroutine ttt_step
!---------------------------
! Write some data in the end
subroutine ttt_write
 integer :: i

 ! Code
 write (6,*) "Angles"
 do i=1,3
  write (6,"(i1,a,2f20.15)") i,":",best_theta(i),best_phi(i)
 end do
end subroutine ttt_write

 
end module ttt
