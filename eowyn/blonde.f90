module blonde
 ! All global data and subroutines for eowyn

 implicit none


 real(8), parameter :: pi=3.14159265358979_8

 ! Global data
 !----- Input data

 real(8) :: J1,J2  ! Exchange integrals of the SS lattice
 real(8) :: J3     ! And in the z-direction
 real(8) :: Bu,Be  ! Uniaxial and exchange anisotropy (>0 for easy axis)

 real(8) :: H1,H2  ! Field range
 integer :: numh   ! Number of H-points

 integer :: nx,ny,nz ! Number of sites in the 3 directions
 ! nx,ny must be EVEN for smooth periodic boundary conditions

 integer :: maxit  ! Maximum number of iterations
 
 integer :: nshuffle ! # of shuffles
 
 real(8) :: lambda,lambda0 ! The relaxation parameter

 real(8) :: tol_g, tol_e ! The desired accuracy of the gradients and energy/site

 ! Input data for blonde_trigen
  
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


 !----------- More data
 real(8) :: H      ! External field (current)
 integer :: pstep  ! How often do we print result

 logical :: e_anis  ! If .true. then exchange anistotropy is used
 integer :: nsite  ! Number of sites: nsite=nx*ny*nz

 ! The exchange interactions (with possible exchange anisotropy)
 integer :: j1x,j1y,j1z,j2x,j2y,j2z,j3x,j3y,j3z

 !----- The spin configuration
 
 ! Angles phi, theta (x,y,z)
 real(8), dimension(:,:,:), allocatable:: phi, theta

 ! The "gradient" dE/d phi , dE/d theta

 real(8), dimension(:,:,:), allocatable :: de_dphi,de_dtheta

 real(8) :: maxgrad_t, maxgrad_p ! The maximum component of each grad

 ! The total energy

 real(8) :: energy

 real(8) :: energy_old, energy_diff,energy_diff2 ! Comparison to the previous energy

 ! The magnetic moment
 real(8) :: mag_mom

!-------The "best" values-----------------
 real(8), dimension(:,:,:), allocatable :: best_phi, best_theta
 real(8) :: best_e, best_m ! Energy and magnetic moment

 integer :: best_shuf ! Shuffle iteration number

!----------------------------------------------

contains

subroutine blonde_seed
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
 end subroutine blonde_seed

 subroutine blonde_init
  ! Init everything

  ! Allocate the arrays
  allocate(phi(nx,ny,nz),theta(nx,ny,nz))
  allocate(best_phi(nx,ny,nz),best_theta(nx,ny,nz))
  allocate(de_dphi(nx,ny,nz),de_dtheta(nx,ny,nz))

  nsite=nx*ny*nz  ! Number of sites

  if (abs(Be)>1.e-10) then
   ! Exchange anisotropy used
   e_anis=.true.
   ! Set the xyz exchange interactions

   j1x=J1
   j1y=J1
   j1z=J1+Be

   j2x=J2
   j2y=J2
   j2z=J2

   j3x=J3
   j3y=J3
   j3z=J3
  else ! No exchange anisotropy
   e_anis=.false.
  end if

 
  ! Seed the random number generator
  call blonde_seed
  
  
 end subroutine blonde_init

!-------------

 subroutine blonde_finish
  deallocate(phi,theta,de_dphi,de_dtheta,best_phi,best_theta)
 end subroutine blonde_finish

!--------------

 ! Randomize the spin configuration
 subroutine blonde_random
  integer :: ix, iy, iz

  real(8) :: p0,t0 ! Random phi, theta

  real(8) :: harv ! The random number

  do ix=1,nx
  do iy=1,ny
  do iz=1,nz
   call random_number(harv)
   p0=2*pi*harv   
   call random_number(harv)
   t0=acos(2.0_8*harv-1.0_8)  ! Proper distribution of angles
   
   phi(ix,iy,iz)=p0
   theta(ix,iy,iz)=t0
  end do ! ix, iy, iz
  end do
  end do
 end subroutine blonde_random

!--------------------------------
! Calculate the contribution of one isotropic exchange bond

subroutine ex_is(J,x1,y1,z1,x2,y2,z2,gp,gt,e)
  real(8), intent(in) :: J  ! The exchange interaction
  integer, intent(in) :: x1,y1,z1,x2,y2,z2  ! The 2 nodes
  real(8), intent(inout) :: gp,gt,e ! The "gradient", energy

  ! Local data
  real(8) :: p1,p2,t1,t2 ! The two theta, phi
  real(8) :: gp0,gt0,e0 ! The contributions

  real(8) :: a1,a2,b1,b2,c,d

  !--- Code
  ! The 2 spins
  t1=theta(x1,y1,z1)
  p1=phi(x1,y1,z1)
  t2=theta(x2,y2,z2)
  p2=phi(x2,y2,z2)

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

!-------------------------------------------
! ex_is2 : simpler version which uses theta, phi and calculates energy only

real(8) function ex_is2(J,t1,p1,t2,p2)
  ! Input parameters
  real(8) :: J  ! The exchange interaction
  real(8) :: t1,p1,t2,p2 ! The 2 spins

  ! Other data
  real(8) :: a1,a2,b1,b2,c,d ! Trig functions

  ! Code
  
  ! The trig functions
  a1=sin(t1)
  a2=sin(t2)
  b1=cos(t1)
  b2=cos(t2)
  c=sin(p1-p2)
  d=cos(p1-p2)

  ! Energy
  
  ex_is2 = J*(a1*a2*d+b1*b2)

end function ex_is2
 


!------------------
! Calculate gradient and energy for a given spin config

 subroutine blonde_grad(nstep)
  integer, intent (in) :: nstep
   
  integer :: ix, iy, iz

   ! The xyz neighbors
  integer :: nex1,nex2,ney1,ney2,nez1,nez2

  ! The diagonal neighbors
  integer :: dix,diy

  ! The contributions to gradient and energy
  real(8) :: dgp, dgt, de

  ! Misc
  real(8) ::a,b
  
  ! ---- Start the code ---
  energy =0
  mag_mom=0
  maxgrad_t=0
  maxgrad_p=0
    
  do ix=1,nx
   ! Find the x neighbors
   nex1=ix-1
   nex2=ix+1
   ! Periodic boundary conditions
   if (nex1==0) nex1=nx
   if (nex2==nx+1) nex2=1

   do iy=1,ny
    ! Find the y neighbors
    ney1=iy-1
    ney2=iy+1
    ! Periodic boundary conditions
    if (ney1==0) ney1=ny
    if (ney2==ny+1) ney2=1

    do iz=1,nz
     ! Find the z neighbors
     nez1=iz-1
     nez2=iz+1
     ! Periodic boundary conditions
     if (nez1==0) nez1=nz
     if (nez2==nz+1) nez2=1
      
     ! Find the diagonal neighbor
     ! odd ix => diy=iy+1, even ix => diy=iy-1
     ! odd iy => dix=ix+1, even iy => dix=ix-1
     ! + wrap around the edges
     ! In particular:
     ! (1,1) <-> (2,2)
     ! (1,2) <-> (0,3)
     ! (2,3) <-> (3,2)
    
     diy=iy-1 + 2* mod(ix,2)
     if (diy==0) diy=ny
     if (diy==ny+1) diy=1

     dix=ix-1 + 2* mod(iy,2)
     if (dix==0) dix=nx
     if (dix==nx+1) dix=1
   
     ! Contributions for each node
     dgp=0
     dgt=0
     de=0

     ! Check for exchange anis
     if (e_anis) then
      print *,"Error: exchange anisotropy not supported yet"
      print *,"Be=",Be
      stop
     end if

     ! Add 4 in-plane neigbors
     call ex_is(J1,ix,iy,iz,nex1,iy,iz,dgp,dgt,de)

     call ex_is(J1,ix,iy,iz,nex2,iy,iz,dgp,dgt,de)

     call ex_is(J1,ix,iy,iz,ix,ney1,iz,dgp,dgt,de)

     call ex_is(J1,ix,iy,iz,ix,ney2,iz,dgp,dgt,de)

     ! Add the diagonal neigbor
     call ex_is(J2,ix,iy,iz,dix,diy,iz,dgp,dgt,de)

     ! Add z-neighbors
     if (nz>1) then
      call ex_is(J3,ix,iy,iz,ix,iy,nez1,dgp,dgt,de)

      call ex_is(J3,ix,iy,iz,ix,iy,nez2,dgp,dgt,de)
     end if ! nz > 1
 
     ! Add uniaxial anisotropy and the Zeeman term
     a=sin(theta(ix,iy,iz))
     b=cos(theta(ix,iy,iz))
     de=de + Bu*a*a/2 - H*b
     dgt=dgt + Bu*a*b + H*a
    
     ! Store this data
     
     mag_mom=mag_mom+b
     de_dphi(ix,iy,iz)=dgp
     de_dtheta(ix,iy,iz)=dgt
     energy=energy+de
     
     ! Max gradient
     if (abs(dgp)>maxgrad_p) maxgrad_p=abs(dgp)
     if (abs(dgt)>maxgrad_t) maxgrad_t=abs(dgt)

    end do ! ix, iy ,iz
   end do
  end do

  ! Normalize per one spin
  energy = energy/nsite
  mag_mom = mag_mom/nsite

  if (nstep>1) then ! Calculate the energy difference
   if (nstep>2) energy_diff2=energy_diff
   energy_diff=energy - energy_old
   if (energy_diff > 0) then
    !lambda=lambda/2  ! Cut the step if the energy INCREASED in the last step
    !write (6,*) nstep,"blonde_grad: lambda decreased to ",lambda
   !else if (abs(energy_diff2/energy_diff-1) < 0.001_8) then
    ! Disable for now
    !lambda=lambda*1.2  ! Increase the step if too slow
    !write (6,*) nstep,"blonde_grad: lambda increased to ",lambda
   end if
  end if

  energy_old=energy ! Store the energy

  
 end subroutine blonde_grad

!-----------------
! Make one step

subroutine blonde_step

  ! local data
  integer :: ix, iy, iz
  real(8) :: delta_p, delta_t, a, t,p

  !real(8) :: lambda0  ! Experimental: randomize lambda a bit
  !real(8) :: harv
 

  ! code
  
  ! Change theta, phi based on gradients, lambda

  do ix=1,nx
  do iy=1,ny
  do iz=1,nz
   a=sin(theta(ix,iy,iz))

   ! Randomize lambda -- not needed yet
   ! call random_number(harv)

   ! lambda0=lambda*(1-harv/3)

   delta_t = - lambda*de_dtheta(ix,iy,iz)
   
   if (abs(a) < 1.e-10) then  ! Safety cutoff for theta=0,pi case to avoid division by 0
    delta_p=0
   else
    delta_p = - lambda*de_dphi(ix,iy,iz)/(a*a)
   end if

   ! Rotate the spins
   t = theta(ix,iy,iz) + delta_t
   p = phi(ix,iy,iz) + delta_p

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
   p = mod(p,2*pi)
   if (p < 0) p=p+2*pi
   
   ! Put the spins back
   theta(ix,iy,iz) = t
   phi(ix,iy,iz) = p

  end do ! ix,iy,iz
  end do
  end do

end subroutine blonde_step

!--------------------------------
! Write the spin config
subroutine blonde_write
 integer :: ix,iy,iz

 open(22,file='spin.dat',status='replace')

 do ix=1,nx
 do iy=1,ny
 do iz=1,nz
  write (22,"(3i4,2f20.14)") ix,iy,iz,best_theta(ix,iy,iz),best_phi(ix,iy,iz)
 end do
 end do 
 end do

 close(22)
end subroutine blonde_write

!--------------------------
! Calculate the energies of all triangles -- to compare with the progam triangle
! Because of the strange bug (?) with the anisotropy and UUD region

subroutine blonde_trigen
  ! Set alpha1,2 if automatic

 integer :: ix,iy           ! Coordinates of the spins = right angle of the triangle
 integer :: ix1,iy1,ix2,iy2 ! The other 2 corners

 real(8) :: rho,hh ! Unitless

 real(8) :: th0,th1,th2,ph0,ph1,ph2 ! 3 spins in polar coords, 0 = tip
  
 real(8) :: energy_trig ! Energy of a triangle
 
 real(8) :: esum ! Total energy

 ! Start the code
 write (6,*) "-----"
 write (6,*) "blonde_trigen:"

 rho=J2/J1
 hh=H/J1

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
   alpha2=alpha1

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
 write (6,*) "alpha2=",alpha2

 esum=0
 open(22,file='trigen.dat',status='replace')
 ! Loop over all spins
 ! Triangles are identified by its right angle
 do ix=1,nx
 do iy=1,ny
     
  ! Find coordinates of 2 other corners

  ix1=ix
  iy1=iy - 2* mod(ix,2)+1
  
  ix2=ix - 2*mod(iy,2)+1
  iy2=iy

  ! Wrap them around
  if (ix1==0) then
   ix1=nx
  else if (ix1==nx+1)  then
   ix1=1
  end if ! ix1

  if (ix2==0) then
   ix2=nx
  else if (ix2==nx+1) then
   ix2=1
  end if ! ix2
  
  if (iy1==0) then
   iy1=ny
  else if (iy1==ny+1) then
   iy1=1
  end if ! iy1
  
  if (iy2==0) then
   iy2=ny
  else if (iy2==ny+1) then
   iy2=1
  end if ! iy2

  ! Find theta, phi of these spins
  th0=best_theta(ix,iy,1)
  ph0=best_phi(ix,iy,1)
  th1=best_theta(ix1,iy1,1)
  ph1=best_phi(ix1,iy1,1)
  th2=best_theta(ix2,iy2,1)
  ph2=best_phi(ix2,iy2,1)

  ! Find the energy : exchange 

  energy_trig = ex_is2(J2/2,th1,ph1,th2,ph2) + &  ! Diagonal bond
    ex_is2(J1,th0,ph0,th1,ph1)+ex_is2(J1,th0,ph0,th2,ph2)
  
  ! External field
  energy_trig = energy_trig - H*alpha1*cos(th0) - H*(1-alpha1)*(cos(th1)+cos(th2))/2
  
  ! Anisotropy
  energy_trig = energy_trig + Bu*alpha2*sin(th0)**2/2 + &
     Bu*(1-alpha2)*(sin(th1)**2 + sin(th2)**2)/4
  
  esum=esum+energy_trig
 
  ! Write the energy
  write (22,"(6i4,f20.14)") ix,iy,ix1,iy1,ix2,iy2,energy_trig
  
 end do !iy
 end do !ix
 close(22)

 write (6,*) "Average E=",esum/(nx*ny)

end subroutine blonde_trigen

! Make one complete run for given parameters (for each shuffle)
subroutine blonde_run(verb)
 logical :: verb ! print data ?
 

 integer :: iter ! Iteration number
 !--------------------

 
    call blonde_random ! Start with random spin config

    ! The main interation loop
    do iter=1,maxit

    call blonde_grad(iter)  ! Find the gradient, energy, etc.

    ! Write the data - every pstep iterations

    if (verb  .and. ((iter/pstep)*pstep==iter .or. iter==1 .or. iter==maxit)) then  ! write, but not  everything
      write (6,"(i6,a,5f15.5)") iter,":", energy, mag_mom, energy_diff, maxgrad_t, maxgrad_p
    end if
    
    ! Check for convergence
    if (verb .and. (iter > 1) .and. (abs(energy_diff) < tol_e) .and. (maxgrad_t < tol_g) .and. (maxgrad_p < tol_g)) then
      write (6,"(i6,a,5f15.5)") iter,":", energy, mag_mom, energy_diff, maxgrad_t, maxgrad_p
      write (6,*) "Convergence reached at iteration : ", iter
      exit
    end if

    
    call blonde_step  ! Make one relaxation step (rotate spins)
  
    ! Check for no convergence
    if (verb .and. iter == maxit) then
      write (6,*) "No convergence after ", maxit, " iterations"
    end if
    
    end do ! iter
end subroutine blonde_run

end module blonde