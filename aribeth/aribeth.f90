! Create a nice xfig file for the SS lattice
! With arrows showing spins from spin.dat as arrows
! Optional: plot the circles

program aribeth
 implicit none
 
 ! Input data
 integer :: nx, ny ! Number of spins along x and y directions
 
 real(8), dimension(:,:), allocatable:: phi, theta  ! Spins from spin.dat

 ! More data

 integer,parameter :: max_xy=8000 ! Max coordinates in xfig plot - not used

 integer,parameter :: max_rad=900 ! Max circle radius

 integer :: rr  ! Circle radius
 integer :: circ_col ! Circle color
 integer :: del_x, del_y ! Arrow direction

 integer,parameter :: ar_len=900 ! Max arrow half-length

 integer :: stepx,stepy,offsx,offsy ! steps, offsets

 integer ::xxmax,yymax ! Max coordinates (max_xy - roundabout errows)

 integer :: xx,yy, xx2,yy2 ! Coords of line ends

 integer :: ix,iy ! Loop counters

 integer :: dum1, dum2, dum3 ! Dummy vars

 

 ! Start the code
 
 ! Read configuration data - nx, ny only
 open(21,file='eowyn.in',status='old')
 read (21,*) nx, ny
 close(21)
 
 allocate(phi(nx,ny),theta(nx,ny))

 ! Read the spins from spin.dat
 open(22,file='spin.dat',status='old')

 do ix=1,nx
 do iy=1,ny
   read (22,"(3i4,2f20.14)") dum1,dum2,dum3,theta(ix,iy),phi(ix,iy)
 end do
 end do

 close(22)
 ! Read the file spin.dat

 ! Data for xfig
 open(24,file='tail',status='replace')

 stepx=2000
 offsx=stepx
 stepy=2000
 offsy=stepy

 xxmax=(nx-1)*stepx+offsx
 yymax=(ny-1)*stepx+offsx

 ! Vertical lines 
 do ix=1,nx
  xx=(ix-1)*stepx+offsx
  write (24,"(a)") "2 1 0 3 2 7 60 -1 -1 0.0 0 0 -1 0 0 2"  ! Line header
  write (24,"(a,4i7)") "     ",xx,offsy,xx,yymax
 end do

 ! Horizontal lines 
 do iy=1,ny
  yy=(ny-iy)*stepy+offsy
  write (24,"(a)") "2 1 0 3 2 7 60 -1 -1 0.0 0 0 -1 0 0 2"  ! Line header
  write (24,"(a,4i7)") "     ",offsx,yy,xxmax,yy
 end do

! Diagonals
do iy=1,ny-1 ! Loop over "rows", minus the last one
  yy=(ny-iy)*stepy+offsy
  yy2=yy-stepy
 
  do ix=1,nx
   xx=(ix-1)*stepx+offsx
   
   if (mod(ix,2)==1) then ! Odd ix only
    xx2=xx+stepx*(2*mod(iy,2)-1)

    if ((xx2 .ge. offsx) .and. (xx2 .le. xxmax)) then ! Add diag line
      write (24,"(a)") "2 1 1 5 2 7 60 -1 -1 12.0 0 0 -1 0 0 2"  ! Line header
      write (24,"(a,4i7)") "     ",xx,yy,xx2,yy2
    end if

   end if  ! mod(ix,2)=1

  end do ! ix
 end do ! iy

 ! Circles and arrows

 do ix=1,nx   ! Loop over all spins
 do iy=1,ny

  ! Centers of the circles/arrow for a given spin
  xx=(ix-1)*stepx+offsx
  yy=(ny-iy)*stepy+offsy

  rr=max_rad*abs(cos(theta(ix,iy)))
  del_x=ar_len*sin(theta(ix,iy))*cos(phi(ix,iy))
  del_y=ar_len*sin(theta(ix,iy))*sin(phi(ix,iy))

  if (cos(theta(ix,iy)) .ge. 0) then
   circ_col=4 ! Spin up color
  else
   circ_col=6 ! Spin down color
  end if !

  ! Circle
  write (24,"(a,2i2,a,8i7)") "1 3 0 5 ",circ_col,circ_col," 50 -1 20 0.0 1 0.0 ",xx,yy,rr,rr,xx,yy,xx+rr,yy

  ! Now the arrow
  write (24,"(a)") "2 1 0 10 0 7 40 -1 -1 12.0 0 0 -1 1 0 2"  ! Line header
  write (24,"(a)") " 	1 1 7.00 200.00 400.00"
  write (24,"(a,4i7)") "     ",xx-del_x,yy-del_y,xx+del_x,yy+del_y  

 end do ! iy
 end do ! ix
 
 close(24)

 call system("cat head tail > temp.fig") 

end program aribeth