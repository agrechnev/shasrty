! Read spin.dat and eowyn.in from eowyn and analyze pairs and "trinagles" of spins
program spin23
 use spin23_run
 implicit none

 ! Input data
 integer :: nx,ny,nz ! Lattice size

 ! Angles phi, theta (x,y,z)
 real(8), dimension(:,:), allocatable:: phi, theta

 ! Other data
 integer :: choice ! Choose between 2 and 3 spins
 
 integer :: ix,iy ! Loop counters
 integer :: iix,iiy,iiz ! Dummy

 integer :: x1,y1,x2,y2,x3,y3 ! Coordinates of chosen nodes

 !-------- Start the code --------------

 ! Read eowyn.in
 open(21,file='eowyn.in',status='old')
 read (21,*) nx,ny,nz  ! Lattice size

 ! Line 2,3 not needed!
 ! read (21,*) J1,J2,J3,Bu,Be,H  ! Hamiltonian parameters
 ! read (21,*) maxit,lambda,tol_g,tol_e  ! Convergence parameters
 close(21)

 if (nz>1) then 
  stop "Error : nz>1 : 3D not supported yet !"
 end if

 ! Allocate the arrays
 allocate(phi(nx,ny),theta(nx,ny))

 ! Read spin.dat
 open(22,file='spin.dat',status='old')

 do ix=1,nx
 do iy=1,ny
  read (22,*) iix,iiy,iiz,theta(ix,iy),phi(ix,iy)
  if ((ix .ne. iix) .or. (iy .ne. iiy)) then
   stop "Error : index mismatch !"
  end if
 end do
 end do 

 close(22)

 ! Do the main part
 do
  print *,"What to do?"
  print *," 1 - exit"
  print *," 2 - 2 spins"
  print *," 3 - 3 spins"
 
  read (5,*) choice
  
  select case(choice)
   case(1)
    exit

   case(2) ! 2 spins
    print *, "Enter x1,y1,x2,y2"
    read (5,*) x1,y1,x2,y2
   
    call run_2spins(theta(x1,y1),phi(x1,y1),theta(x2,y2),phi(x2,y2))

   case(3) ! 3 spins
    print *, "Enter x1,y1,x2,y2,x3,y3"
    read (5,*) x1,y1,x2,y2,x3,y3
   
    call run_3spins(theta(x1,y1),phi(x1,y1),theta(x2,y2),phi(x2,y2),theta(x3,y3),phi(x3,y3))

  end select 
 end do


end program spin23