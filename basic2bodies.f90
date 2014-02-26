program nbody
	implicit none

	integer :: i,j,n, counter
	double precision :: G, dt, rad, year, t
	! vectors in format (x1,y1,z1), (x2,y2,z2)
	double precision, DIMENSION(3, 2) :: pos,v,a
	double precision, DIMENSION(2) :: m
	double precision, DIMENSION(3) :: r
	counter = 0
	t = 0.0
	year = 3.15570e7 !number of seconds in a year
	dt = 1 ! time interval (s)
	n = 2
	G = 6.67e-11 ! G (Nm**2/kg**2)
	pos = 0.0 ! position coordinates (m)
	v = 0.0 ! velocities (m/s)
	a = 0.0 ! acceleration (m/s*2)
	r = 0.0 ! direction of the gravitational force
	m(1) = 1.99e30 ! mass sun (kg)
	m(2) = 5.97e24 ! mass earth (kg)
	pos(1,2) = 1.496e11 ! earth sun distance (m)
	v(2,2) = 29.8e3 ! earth velocity (m/s)

	do while (t < year)
		do i = 1,n
			do j = 1,n
				if (j == i) CYCLE
				r(1:3) = (pos(1:3,i)-pos(1:3,j))
				rad = sqrt((pos(1,i)-pos(1,j))**2 + (pos(2,i)-pos(2,j))**2 + (pos(3,i)-pos(3,j))**2)
				a(1:3,j) = -((G * m(i))/rad**3)*r(1:3)
				v(1:3,j) = v(1:3,j) - a(1:3,j)*dt
			end do
			pos(1:3,i) = pos(1:3,i) + v(1:3,i)*dt
			if (MOD(counter,1000000) == 0) print *, pos(1:2,1)
			counter = counter + 1
		end do
		t = t + dt
	end do 
end program nbody
