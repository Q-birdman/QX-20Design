program plot
  implicit none
  integer ::n,i,j,k,l
  real(4),dimension(4,3) ::coord, coord1, coord2
  real(4),dimension(:),allocatable::phi,theta,psi,x0,y0,z0
  real(4),dimension(4,4) ::f1, f2, f3, ff1, ff2
  real(4)::c1,s1,c2,s2,c3,s3,pi
  character ::filename*128

  pi=4*atan(1.0)


  
  
  open(10,file = "log.dat", status = 'old')
  read(10,*) n
  
  allocate(x0(0:n))
  allocate(y0(0:n))
  allocate(z0(0:n))
  allocate(phi(0:n))
  allocate(theta(0:n))
  allocate(psi(0:n))
  
  ! Initial state
  read(10,*) x0(0),y0(0),z0(0),phi(0),theta(0),psi(0)
  
  ! Read datas
  do i = 1, n
     read(10,*) x0(i),y0(i),z0(i),phi(i),theta(i),psi(i)
  end do
  
  close(10)
  
   open(11,file = "plot.dat", status = 'old')
  
  
  do i = 0, n, 300
  
  coord(1,1) = 5
  coord(2:3,1) = 0
  coord(4,1:3) = 1
  coord(1,2:3) = 0
  coord(3,2:3) = 0
  coord(2,2) = 10
  coord(2,3) = -10
  
  coord1 = 0.0
  coord2 = 0.0
  
  c1 = cos(phi(i)*pi/180); c2 = cos(theta(i)*pi/180); c3 = cos(psi(i)*pi/180)
  s1 = sin(phi(i)*pi/180); s2 = sin(theta(i)*pi/180); s3 = sin(psi(i)*pi/180)

  !call rotate(phi(i),theta(i),psi(i),0.0,0.0,0.0,3,coord,coord1)
  coord1(1,1) = 5*c2*c3
  coord1(2,1) = 5*c2*s3
  coord1(3,1) = 5*s2
  coord1(1,2) = 10*(-c1*s3)
  coord1(2,2) = 10*c1*c3
  coord1(3,2) = 10*s1  
  coord1(1,3) = -10*(-c1*s3)
  coord1(2,3) = -10*c1*c3
  coord1(3,3) = -10*s1 

  coord2(1,1) = coord1(1,1) + x0(i)  
  coord2(2,1) = coord1(2,1) + y0(i)
  coord2(3,1) = coord1(3,1) + z0(i)
  coord2(4,1) = 1
  coord2(1,2) = coord1(1,2) + x0(i) 
  coord2(2,2) = coord1(2,2) + y0(i)
  coord2(3,2) = coord1(3,2) + z0(i)
  coord2(4,2) = 1
  coord2(1,3) = coord1(1,3) + x0(i) 
  coord2(2,3) = coord1(2,3) + y0(i)
  coord2(3,3) = coord1(3,3) + z0(i)
  coord2(4,3) = 1

  !write(filename, '("plot", i3.3, ".dat")')i
  !open(10,"/plots/")
  !open(11,file = filename,status = 'replace')
  !open(11,file = "plot.dat", status = 'replace')
  !write(11,*) coord2(1,1), coord2(2,1), coord2(3,1)
  !write(11,*) coord2(1,2), coord2(2,2), coord2(3,2)
  !write(11,*) coord2(1,3), coord2(2,3), coord2(3,3)
  do j = 1,100
     write(11,*) (coord2(1,1)*j/100)+(coord2(1,2)*(100-j)/100)&
	             ,(coord2(2,1)*j/100)+(coord2(2,2)*(100-j)/100)&
				 ,(coord2(3,1)*j/100)+(coord2(3,2)*(100-j)/100)
				 
	 write(11,*) (coord2(1,2)*j/100)+(coord2(1,3)*(100-j)/100)&
	             ,(coord2(2,2)*j/100)+(coord2(2,3)*(100-j)/100)&
				 ,(coord2(3,2)*j/100)+(coord2(3,3)*(100-j)/100)
				 
	 write(11,*) (coord2(1,3)*j/100)+(coord2(1,1)*(100-j)/100)&
	             ,(coord2(2,3)*j/100)+(coord2(2,1)*(100-j)/100)&
				 ,(coord2(3,3)*j/100)+(coord2(3,1)*(100-j)/100) 
  end do
  write(11,'(a)')
  
  end do
  
  close(11)
  !close(10)
  
  open(12,file = "plotlog.dat",status = 'replace')
  
  do i = 0, n
  write(12,*) x0(i), y0(i), z0(i)
  end do

  close(12)
  
  
end program
  



subroutine rotate(phi,theta,psi,x0,y0,z0,n,coord,coord1) !x0,y0,z0を回転中心として回転させる
  implicit none
  real(4),intent(IN)::phi,theta,psi,x0,y0,z0
  integer,intent(IN)::n
  real(4),dimension(4,4)::r
  real(4),dimension(4,n),intent(IN)::coord
  real(4),dimension(4,n),intent(OUT)::coord1
  real(4),dimension(4,n)::coord_otw1,coord_otw2 !on the way
  real(4)::l1,m1,n1,l2,m2,n2,l3,m3,n3
  integer::i, j, k
  real(4)::s

  !回転行列の定義
  l1 = cos(phi); l2 = cos(theta); l3 = cos(psi)
  m1 = sin(phi); m2 = sin(theta); m3 = sin(psi)
  r = 0.0
  r(1,1) = l1*l2
  r(2,1) = m1*l2
  r(3,1) = -m2
  r(1,2) = l1*m2*m3-m1*l3
  r(2,2) = m1*m2*m3+l1*l3
  r(3,2) = l2*m3
  r(1,3) = l1*m2*l3+m1*m3
  r(2,3) = m1*m2*l3-l1*m3
  r(3,3) = l2*l3
  r(4,4) = 1
  
  !原点に移動→回転→元の位置に移動
  call parallel_translation(-x0,-y0,-z0,n,coord,coord_otw1)
  coord_otw2 = 0.0
  do i = 1, 4
     do j = 1, n
        s = 0
        do k = 1, 4
           s = s+r(i,k)*coord_otw1(k,j)
        end do
        coord_otw2(i,j) = s
     end do
  end do
  call parallel_translation(x0,y0,z0,n,coord_otw2,coord1)
  
end subroutine rotate


subroutine scaling(a,b,c,n,coord,coord1) !拡大縮小
  implicit none
  real(4),intent(IN)::a,b,c
  integer,intent(IN)::n
  real(4),dimension(4,4)::s
  real(4),dimension(4,n),intent(IN)::coord
  real(4),dimension(4,n),intent(OUT)::coord1
  integer::i, j, k
  real(4)::t

  !拡大縮小行列の定義
  s = 0.0
  s(1,1) = a
  s(2,2) = b
  s(3,3) = c
  s(4,4) = 1

  !拡大縮小
  coord1 = 0
  do i = 1, 4
     do j = 1, n
        t = 0
        do k = 1, 4
           t = t+s(i,k)*coord(k,j)
        end do
        coord1(i,j) = t
     end do
  end do
     
end subroutine scaling


subroutine parallel_translation(p,q,r,n,coord,coord1) !平行移動
  implicit none
  real(4),intent(IN)::p,q,r
  integer,intent(IN)::n
  real(4),dimension(4,4)::pt
  real(4),dimension(4,n),intent(IN)::coord
  real(4),dimension(4,n),intent(OUT)::coord1
  integer::i, j, k
  real(4)::s

  !平行移動行列の定義
  pt = 0.0
  do i = 1, 4
     pt(i,i) = 1
  end do
  pt(1,4) = p
  pt(2,4) = q
  pt(3,4) = r

  !平行移動
  coord1 = 0.0
  do i = 1, 4
     do j = 1, n
        s = 0
        do k = 1, 4
           s = s+pt(i,k)*coord(k,j)
        end do
        coord1(i,j) = s
     end do
  end do
   
end subroutine parallel_translation