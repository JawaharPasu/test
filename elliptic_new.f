
	implicit none
	integer M,N
	parameter(M=21)
	parameter(N=41)
	real e,q,dtime,etime,t(2),tim(2)
	integer i,j,k,choice,t1,t2,iter
	double precision l,h,delx,dely,beta,x,y,eva,om
	double precision TEMP(M,N),TEMP1(M,N),TEMP2(M,N),Val(M,N)

	common /yval1/ i,j,beta,eva,TEMP,TEMP1,TEMP2,om

	open(10,file='elliptic_t.dat')
	l = 1.0d0		!x-length
	h = 2.0d0		!y-length
	!om = 1.2d0		!relaxation parameter
	iter=0			!initialize iteration
	delx = l/(M-1)		!step length in x-direction
	dely = h/(N-1)		!step length in y-direction
	beta = (delx**2)/(dely**2)	!beta value
	eva = 1/(2*(1+beta**2))
	print *,beta
!-----------------------------------------------------------------------
	!initializing boundary conditions
9	do i = 1,M
		do j = 1,N
		TEMP(i,j) = 0.d0
		TEMP1(i,j) = TEMP(i,j)
		TEMP2(i,j) = TEMP(i,j)
		end do
	end do
10	do i = 1,M
	TEMP(i,1) = 100.d0
	TEMP1(i,1)= 100.d0
	TEMP2(i,1)= 100.d0
	end do
	TEMP(1,1)= 50.d0
	TEMP(21,1)=50.d0
	TEMP1(1,1)= 50.d0
	TEMP1(21,1)=50.d0
	TEMP2(1,1)= 50.d0
	TEMP2(21,1)=50.d0
!------------------------------------------------------------------------
	write(*,*) "Jacobi - 1"
	write(*,*) "Point Gauss seidel - 2"
	write(*,*) "line Gauss seidel - 3"
	write(*,*) "PSOR - 4"
	write(*,*) "LSOR - 5"
	write(*,*) "ADI - 6"
	write(*,*) "ADI SOR - 7"
11	write(*,*) "Enter a value="
	read(*,*) choice
	if (choice.eq. 1) then
	call jacobi
	else if (choice.eq. 2) then
	call pgs
	else if (choice.eq. 3) then
	call lgsi
	else if (choice.eq. 4) then
	write(*,*) "enter relaxation parameter value less than 1.8"
	read(*,*) om
	call psor
	else if (choice.eq. 5) then
	write(*,*) "enter relaxation parameter value less than 1.3"
	read(*,*) om
	call lsor
	else if (choice.eq. 6) then
	call adi
	else if (choice.eq. 7) then
	write(*,*) "enter relaxation parameter value greater than 1"
	read(*,*) om
	call adisor
	else
	write(*,*) "wrong choice,pls enter correct value (1-7)"
	goto 11
	end if

12	if (choice.eq. 1) then
	call jacobi
	else if (choice.eq. 2) then
	call pgs
	else if (choice.eq. 3) then
	call lgsi
	else if (choice.eq. 4) then
	call psor
	else if (choice.eq. 5) then
	call lsor
	else if (choice.eq. 6) then
	call adi
	else if (choice.eq. 7) then
	call adisor
	else
	end if
	val = TEMP1 - TEMP
	TEMP = TEMP1
13 	do i = 2,M-1
		do j = 2,N-1
14 			if (abs(val(i,j)) .le. 0.001) then
				exit
			else
				iter=iter+1
				goto 12
			end if
		end do
	end do
	TEMP = TEMP - val
	x = 0.d0
	y = 0.d0
	write(10,*) 'TITLE="TEMP DIST"'
	write(10,*) 'VARIABLES = "X","Y" "TEMP"'
	write(10,*) 'ZONE T="ONLY ZONE", I=21,J=41, F=POINT'
	write(*,*) "x , ", "y, ", "Temperatures"
	do j = 1,N
		do i = 1,M
		write(10,*) x,y,TEMP1(i,j)
		write(*,*) x,y,TEMP1(i,j)
		x = x + delx
		end do
		x = 0.d0
		y = y + dely
	end do
	write(*,*) "No. of iterations=",iter
	e=dtime(t)
	q=etime(tim)
	write(*,*) "CPU time=",e,"total time=",q
	stop
	end
!-------------------------------------------------------------------------
	subroutine jacobi
	implicit none
	integer M,N
	parameter(M=21)
	parameter(N=41)
	integer i,j,k
	double precision l,h,delx,dely,beta,x,y,eva,tend,om
	double precision TEMP(M,N),TEMP1(M,N),TEMP2(M,N),Val(M,N)

	common /yval1/ i,j,beta,eva,TEMP,TEMP1,TEMP2,om

	do i = 2,M-1
		do j = 2,N-1
		tend = (beta**2)*((TEMP(i,j+1))+(TEMP(i,j-1)))
		TEMP1(i,j) = eva*(TEMP((i+1),j) + TEMP((i-1),j) + tend)
		end do
	end do
	return
	end subroutine
!-------------------------------------------------------------------------
	subroutine pgs
	implicit none
	integer M,N
	parameter(M=21)
	parameter(N=41)
	integer i,j,k
	double precision l,h,delx,dely,beta,x,y,eva,tend,om
	double precision TEMP(M,N),TEMP1(M,N),TEMP2(M,N),Val(M,N)

	common /yval1/ i,j,beta,eva,TEMP,TEMP1,TEMP2,om

	do i = 2,M-1
		do j = 2,N-1
		tend = (beta**2)*((TEMP(i,j+1))+(TEMP1(i,j-1)))
		TEMP1(i,j) = eva*(TEMP((i+1),j) + TEMP1((i-1),j) + tend)
		end do
	end do
	return
	end subroutine
!--------------------------------------------------------------------------
	subroutine psor
	implicit none
	integer M,N
	parameter(M=21)
	parameter(N=41)
	integer i,j,k
	double precision l,h,delx,dely,beta,x,y,eva,tend,om,tat
	double precision TEMP(M,N),TEMP1(M,N),TEMP2(M,N),Val(M,N)

	common /yval1/ i,j,beta,eva,TEMP,TEMP1,TEMP2,om
	do i = 2,M-1
		do j = 2,N-1
		tend = (beta**2)*((TEMP(i,j+1))+(TEMP1(i,j-1)))
		tat = eva*(TEMP((i+1),j) + TEMP1((i-1),j) + tend)
		TEMP1(i,j)=((1-om)*TEMP(i,j))+(om*tat)
		end do
	end do
	return
	end subroutine
!--------------------------------------------------------------------------
	subroutine lgsi
	implicit none
	integer M,N
	parameter(M=21)
	parameter(N=41)
	integer i,j,k
	double precision l,h,delx,dely,beta,x,y,eva,tend,om,tat
	double precision TEMP(M,N),TEMP1(M,N),TEMP2(M,N),Val(M,N)
	double precision a(M),b(M),c(M),d(M),nan,cas
	
	common /yval1/ i,j,beta,eva,TEMP,TEMP1,TEMP2,om

	do j=2,N-1

	
	d(2)=-(2*(1+beta**2))
	
	c(2)=-(beta**2)*(TEMP1(2,j-1)+TEMP(2,j+1))
	
	do i=3,M-1
	a(i-1)=1.d0
	b(i)=1.d0
	d(i)=-(2*(1+beta**2))
	d(i)=d(i)-(b(i)*a(i-1))/d(i-1)
	
	nan=-(beta**2)*(TEMP1(i,j-1)+TEMP(i,j+1))
	
	c(i)=nan-(b(i)*c(i-1))/d(i-1)
	
	
	end do
	
	TEMP1(20,j)=c(M-1)/d(M-1)  
	
	do i=M-2,2,-1

	TEMP1(i,j)=(c(i)-(a(i)*TEMP1(i+1,j)))/d(i)
	
	end do
	end do
	return
	end subroutine
!---------------------------------------------------------------------------

!---------------------------------------------------------------------------
	subroutine adi
	implicit none
	integer M,N
	parameter(M=21)
	parameter(N=41)
	integer i,j,k,can,cud,lo,jo
	double precision l,h,delx,dely,beta,x,y,eva,tend,om,tat
	double precision TEMP(M,N),TEMP1(M,N),TEMP2(M,N),Val(M,N)
	double precision a1(M),b1(M),c1(M),d1(M)
	double precision a2(N),b2(N),c2(N),d2(N)
	double precision nan,cas,tak

	common /yval1/ i,j,beta,eva,TEMP,TEMP1,TEMP2,om

	do j=2,N-1
	a1(2)=1.d0
	b1(2)=1.d0
	d1(2)=-2*(1+beta**2)
	c1(2)=-(beta**2)*(TEMP2(2,j-1)+TEMP(2,j+1))
	do i=3,M-1
	a1(i-1)=1.d0
	b1(i)=1.d0
	d1(i)=-(2*(1+beta**2))
	d1(i)=d1(i)-(b1(i)*a1(i-1))/d1(i-1)
	
	nan=-(beta**2)*(TEMP2(i,j-1)+TEMP(i,j+1))
	
	c1(i)=nan-(b1(i)*c1(i-1))/d1(i-1)
	!write(*,*) i
	!write(*,*) c1(i)
	end do
	
	TEMP2(20,j)=c1(M-1)/d1(M-1) 
	!nan=19
	do i=M-2,2,-1
	!write(*,*) i
	TEMP2(i,j)=(c1(i)-(a1(i)*TEMP2(i+1,j)))/d1(i)
	!nan=nan-1
	!write(*,*) i,j,TEMP2(i,j)
	end do
	end do


	do i=2,M-1
	a2(2)=beta**2
	b2(2)=beta**2
	d2(2)=-2*(1+beta**2)
	c2(2)=-(TEMP1(i-1,2)+TEMP2(i+1,2))
	do j=3,N-1
	a2(j-1)=beta**2
	b2(j)=beta**2
	d2(j)=-(2*(1+beta**2))
	d2(j)=d2(j)-(b2(j)*a2(j-1))/d2(j-1)
	
	nan=-(TEMP1(i-1,j)+TEMP2(i+1,j))
	
	c2(j)=nan-(b2(j)*c2(j-1))/d2(j-1)
	!write(*,*) j
	end do
	!write(*,*) c2(j)
	TEMP1(i,40)=c2(N-1)/d2(N-1) 
	
	do j=N-2,2,-1

	TEMP1(i,j)=(c2(j)-(a2(j)*TEMP1(i,j+1)))/d2(j)
	!write(*,*) i,j,TEMP1(i,j)
	
	end do
	end do
	
	return
	end subroutine
!------------------------------------------------------------------------
	subroutine adisor
	implicit none
	integer M,N
	parameter(M=21)
	parameter(N=41)
	integer i,j,k,can,cud,lo,jo
	double precision l,h,delx,dely,beta,x,y,eva,tend,om,tat
	double precision TEMP(M,N),TEMP1(M,N),TEMP2(M,N),Val(M,N)
	double precision a1(M),b1(M),c1(M),d1(M)
	double precision a2(N),b2(N),c2(N),d2(N)
	double precision nan,cas,tak

	common /yval1/ i,j,beta,eva,TEMP,TEMP1,TEMP2,om

	do j=2,N-1
	a1(2)=om
	b1(2)=om
	d1(2)=-2*(1+beta**2)
	eva =-om*((beta**2)*(TEMP2(2,j-1)+TEMP(2,j+1)))
	c1(2)=-((1-om)*(2*(1+beta**2))*TEMP(2,j))+eva
	do i=3,M-1
	a1(i-1)=1.d0
	b1(i)=1.d0
	d1(i)=-(2*(1+beta**2))
	d1(i)=d1(i)-(b1(i)*a1(i-1))/d1(i-1)
	eva =-om*(beta**2)*(TEMP2(i,j-1)+TEMP(i,j+1))
	nan=-((1-om)*(2*(1+beta**2))*TEMP(i,j))+eva
	
	c1(i)=nan-(b1(i)*c1(i-1))/d1(i-1)
	!write(*,*) i
	!write(*,*) c1(i)
	end do
	
	TEMP2(20,j)=c1(M-1)/d1(M-1) 
	!nan=19
	do i=M-2,2,-1
	!write(*,*) i
	TEMP2(i,j)=(c1(i)-(a1(i)*TEMP2(i+1,j)))/d1(i)
	!nan=nan-1
	!write(*,*) i,j,TEMP2(i,j)
	end do
	end do


	do i=2,M-1
	a2(2)=om*beta**2
	b2(2)=om*beta**2
	d2(2)=-2*(1+beta**2)
	eva=-om*(TEMP1(i-1,2)+TEMP2(i+1,2))
	c2(2)=-((1-om)*(2*(1+beta**2))*TEMP2(i,2))+eva
	do j=3,N-1
	a2(j-1)=beta**2
	b2(j)=beta**2
	d2(j)=-(2*(1+beta**2))
	d2(j)=d2(j)-(b2(j)*a2(j-1))/d2(j-1)
	eva=-om*(TEMP1(i-1,j)+TEMP2(i+1,j))
	nan=-((1-om)*(2*(1+beta**2))*TEMP2(i,j))+eva
	
	c2(j)=nan-(b2(j)*c2(j-1))/d2(j-1)
	!write(*,*) j
	end do
	!write(*,*) c2(j)
	TEMP1(i,40)=c2(N-1)/d2(N-1) 
	
	do j=N-2,2,-1

	TEMP1(i,j)=(c2(j)-(a2(j)*TEMP1(i,j+1)))/d2(j)
	!write(*,*) i,j,TEMP1(i,j)
	
	end do
	end do
	
	return
	end subroutine
	
