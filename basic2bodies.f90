program nbody
	implicit none

	integer :: i,j,n
	double precision :: G, dt, rad
	double precision, DIMENSION(3, 2) :: r,v,a
	double precision, DIMENSION(2) :: m
	double precision, DIMENSION(3) :: dvector
	dt = 1 ! time interval
	n = 2
	G = 6.67e-11 ! G Nm**2/kg**2
	r = 0.0 ! coordinates
	v = 0.0 ! velocities
	a = 0.0 ! acceleration
	dvector = 0.0 ! direction of the gravitational force
	m(1) = 1.99e30 ! mass sun kg
	m(2) = 5.97e24 ! mass earth kg
	r(1,2) = 1.496e11 ! earth sun distance m
	v(1,2) = 29.8e3 ! earth velocity km/s
	
	do i = 1,n
		do j = 1,n
			if (j == i) CYCLE
			dvector(1:3) = r(1:3,i)-r(1:3,j)
			rad = sqrt((r(1,i)-r(1,j))**2 + (r(2,i)-r(2,j))**2 + (r(3,i)-r(3,j))**2)
			a(1:3,j) = ((G * m(i))/rad**3)*dvector(1:3)
			v(1:3,j) = v(1:3,j) + a(1:3,j)*dt
			r(1:3,j) = r(1:3,j) + v(1:3,j)*dt
		end do
	end do 
	print *, rad
end program nbody
