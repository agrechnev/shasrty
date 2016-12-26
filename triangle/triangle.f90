! Just like Eowyn, but relax one triangle (J1,J1,J2/2 bonds), Bu and H
! B (uniaxial anisotropy) under construction

program triangle 

 use ttt
 implicit none

 integer :: iter ! Iteration loop counter

 integer :: shuf ! Shuffle loop counter

 integer :: pstep  ! How often do we print result

 real(8) :: lll


 !----------- Start the code ----------
 ! Read data
 open(21,file='triangle.in',status='old')
 read (21,*) J1,J2,Bu,H  ! Hamiltonian parameters
 read (21,*) mode_alpha, alpha1, alpha2
 read (21,*) maxit,lambda,tol_g,tol_e  ! Convergence parameters
 read (21,*) nshuffle ! # of shuffles
 close(21)

 write (6,*) "I'm just a little triangle ..."
 write (6,*) "But am I not perfect?"

 ! Unitless params
 rho=J2/J1
 hh=H/J1
 beta=Bu/J1

 ! Set up alpha if automatic

 select case (mode_alpha) 
  ! If mode_alpha=0 then use the alpha1,2 from the input file

  case (1)  ! 1 = Neel : alpha1,2=1/2
   write (6,*) "1: Neel"
   alpha1=0.5_8
   alpha2=alpha1

  case(2) ! 2 = Flip (FM) : alpha1,2=1/2 (alphas do not matter)
   write (6,*) "2: Spin Flip (FM)"
   alpha1=0.5_8
   alpha2=alpha1

  case(3) ! 3 = Dimer : alpha1,2=0
   write (6,*) "3: Dimer"
   alpha1=0
   alpha2=0

  case(4) ! 4 = UUD : alpha1,2= (2-rho+hh)/(3*hh)
   write (6,*) "4: UUD"
   alpha1= (2-rho+hh)/(3*hh)
   ! Skip for testing purposes
   !alpha2=alpha1  

  case(5) ! 5 = Spin Flop : alpha1,2=1/2
   write (6,*) "5: Spin Flop"
   alpha1=0.5_8
   alpha2=alpha1

  case(6) ! 6 = Umbrella : alpha1,2=1/(rho+1) Is this right ?
   write (6,*) "6: Umbrella"
   alpha1=1/(rho+1)
   alpha2=alpha1
 end select

 write (6,*) "alpha1=",alpha1
 write (6,*) "(1-alpha1)=",1-alpha1


 write (6,*) "alpha2=",alpha2
 write (6,*) "1-alpha2=",1-alpha2

 pstep = max(maxit/20,1) ! Set up pstep

 ! Seed the random number generator
 !call ttt_seed
 call init_random_seed
 
 best_e=0
 ! Loop over shuffles
 do shuf=1,nshuffle
  
  write (6,*) "shuf = ", shuf  

  ! Randomize the configuration
   call ttt_random

  ! For temp debugging: read trig.in instead of randomization
  
   !open(21,file='trig.dat',status='old')
   !read (21,*) theta(2),phi(2)
   !read (21,*) theta(1),phi(1)
   !read (21,*) theta(3),phi(3)
   !close(21)
  
  ! Data template

  write (6,"(a)") "iter : energy, mag_mom, energy_diff, maxgrad_t, maxgrad_p "
  write (6,"(a)") "----------------------------------------------------------"

  ! The main loop
  do iter=1,maxit

   call ttt_grad(iter)  ! Find the gradient, energy, etc.

   ! Write the data - every pstep iterations

   if ((iter/pstep)*pstep==iter .or. iter==1 .or. iter==maxit) then  ! write, but not  everything
    write (6,"(i6,a,5f15.5)") iter,":", energy, mag_mom, energy_diff, maxgrad_t, maxgrad_p
   end if
  
   ! Check for convergence
   if ((iter > 1) .and. (abs(energy_diff) < tol_e) .and. (maxgrad_t < tol_g) .and. (maxgrad_p < tol_g)) then
     write (6,"(i6,a,5f15.5)") iter,":", energy, mag_mom, energy_diff, maxgrad_t, maxgrad_p
     write (6,*) "Convergence reached at iteration : ", iter
     exit
   end if

  
   call ttt_step  ! Make one relaxation step (rotate spins)
 
   ! Check for no convergence
   if (iter == maxit) then
    write (6,*) "No convergence after ", maxit, " iterations"
   end if

  end do ! iter

  ! Find the best energy value algorithm
  if ((shuf==1) .or. (energy < best_e-tol_e/2)) then
   !write (6,*) "Found !",energy,best_e,best_e-tol_e/2

   best_phi=phi
   best_theta=theta
   best_e=energy
   best_m=mag_mom
   best_shuf=shuf
  end if
 end do ! shuf

 ! Write the output
 write (6,*) "Found at shuffle : ",best_shuf
 write (6,*) "E=",best_e
 write (6,*) "M=",best_m
  
 ! Write theoretic energy (all cases)
 
 select case (mode_alpha)
  ! Nothing for mode_alpha=0, obviously
  
  case(1) ! Neel
   write (6,*) "Theory Neel:"
   write (6,*) "E=",-2*J1+J2/2
   write (6,*) "M=",0
  case(2) ! Flip
   write (6,*) "Theory Flip:"
   write (6,*) "E=",+2*J1+J2/2 - H
   write (6,*) "M=",1
  case(3) ! Dimer
   write (6,*) "Theory Dimer:"
   write (6,*) "E=",-J2/2
   write (6,*) "M=",0
  case(4) ! UUD
   write (6,*) "Theory UUD:"
   write (6,*) "E=",-J1*2/3 - J2/6 - H/3
   write (6,*) "M=",1.0_8/3
  case(5) ! Flop
   write (6,*) "Theory Spin Flop:"
   write (6,*) "E=",-2*J1+J2/2 + Bu/2 - H*H/(2*(8*J1-Bu))
   write (6,*) "M=", H/(8*J1-Bu)
  case(6) ! Umbrella
   write (6,*) "Theory Umbrella:"
   write (6,*) "E=",-J2/2-J1*J1/J2+ Bu/2-H*H*J2/(4*(J1+J2)*(J1+J2)-2*Bu*J2)
   write (6,*) "M=",H*J2/(4*(J1+J2)*(J1+J2)-2*Bu*J2)
 end select

 !--------------
 
 call ttt_write  !Write the 3 angles

end program triangle
