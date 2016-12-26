! Read spin.dat and eowyn.in from eowyn and count the number of types
program count
 use ccc
 implicit none

 ! Input data
 integer :: nx,ny,nz ! Lattice size

 ! Angles phi, theta (x,y,z)
 real(8), dimension(:,:), allocatable:: phi, theta

 ! Other data

 integer :: ix,iy,i ! Loop counters
 integer :: iix,iiy,iiz ! Dummy

 ! Data for counting

 integer :: ntype  ! Number of types

 real(8), dimension(:,:), allocatable :: stype  ! stype(type,i) Vector s for each type

 logical :: foundflag ! .true. = spin blongs to an existing type

 logical :: planar ! Planar on not ?
 real(8) :: phi0 ! The first phi

 integer :: type ! Type for current ix,iy
  
 real(8) :: s(1:3)  ! Spin for current ix,iy as a vector

 ! Data for xfig
 integer,parameter :: rr=250 ! radius

 integer,parameter :: max_xy=8000 ! Max coordinates in xfig plot

 integer :: stepx,stepy,offsx,offsy ! steps, offsets
 
 integer :: xx,yy, xx2,yy2 ! Coords of circle center, line etc

 integer ::xxmax,yymax ! Max coordinates (max_xy - roundabout errows)
 
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
 allocate(stype(1:nx*ny,1:3)) ! Allocate the type array
 
 open(23,file='type.dat',status='replace')

 ! Data for xfig
 open(24,file='tail',status='replace')

 stepx=max_xy/max(nx,ny)
 offsx=stepx
 stepy=stepx
 offsy=stepy

 ntype=0
 planar=.true. 
 phi0=phi(1,1)

 do ix=1,nx   ! Loop over all spins
 do iy=1,ny

  ! Check phi's for planar
  if (abs(sin(phi0-phi(ix,iy))) > 1.0e-4) planar=.false.

  foundflag=.false. ! No existing types  -- to start with
  s=tp2v(theta(ix,iy),phi(ix,iy)) ! Current spin (ix,iy) as vector
    
  do i=1,ntype ! Check for existing types

    if (match(s,stype(i,1:3))) then ! Found !
     foundflag=.true.
     type=i
    end if !

  end do

  if (.not. foundflag) then ! Not found: a new type
   ntype=ntype+1
   type=ntype
   stype(ntype,1:3)=s
  end if

  ! Print the data
  write (23,"(3i10)") ix,iy,type

  ! Data for xfig
  xx=(ix-1)*stepx+offsx
  yy=(ny-iy)*stepy+offsy

  write (24,"(a,2i2,a,8i5)") "1 3 0 1 ",type,type," 50 -1 20 0.0 1 0.0 ",xx,yy,rr,rr,xx,yy,xx+rr,yy

 end do !iy
 end do !ix

 close(23)

 
 ! Write lattice info into the xfig plot
 
 xxmax=(nx-1)*stepx+offsx
 yymax=(ny-1)*stepx+offsx
 
 ! Vertical lines
 do ix=1,nx
  xx=(ix-1)*stepx+offsx
  write (24,"(a)") "2 1 0 2 0 7 60 -1 -1 0.0 0 0 -1 0 0 2"  ! Line header
  write (24,"(a,4i7)") "     ",xx,offsy,xx,yymax
 end do

 
 ! Horizontal lines
 do iy=1,ny
  yy=(ny-iy)*stepy+offsy
  write (24,"(a)") "2 1 0 2 0 7 60 -1 -1 0.0 0 0 -1 0 0 2"  ! Line header
  write (24,"(a,4i7)") "     ",offsx,yy,xxmax,yy
 end do

 ! Diagonal lines : the difficult part

 do iy=1,ny-1 ! Loop over "rows", minus the last one
  yy=(ny-iy)*stepy+offsy
  yy2=yy-stepy
 
  do ix=1,nx
   xx=(ix-1)*stepx+offsx
   
   if (mod(ix,2)==1) then ! Odd ix only
    xx2=xx+stepx*(2*mod(iy,2)-1)

    if ((xx2 .ge. offsx) .and. (xx2 .le. xxmax)) then ! Add diag line
      write (24,"(a)") "2 1 0 2 0 7 60 -1 -1 0.0 0 0 -1 0 0 2"  ! Line header
      write (24,"(a,4i7)") "     ",xx,yy,xx2,yy2
    end if

   end if  ! mod(ix,2)=1

  end do ! ix
 end do ! iy

 close(24)

 call system("cat head tail > temp.fig") 

 print *,"Found ",ntype," types !"

 if (planar) then
  print *,"PLANAR"
 else
  print *,"NOT PLANAR"
 end if !planar

 deallocate(stype)
end program count