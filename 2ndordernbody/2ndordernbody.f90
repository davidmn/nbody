program nbody
	implicit none

	integer :: i,j,n, counter, l, k
	double precision :: G, dt, rad, year, t, e_0, e_t, ke_tot,ke_temp, u_temp, u_tot, m_tot
	double precision :: vel_temp
	! vectors in format (x1,y1,z1), (x2,y2,z2)
	double precision, DIMENSION(3, 3) :: pos,v,v2,a
	double precision, DIMENSION(3) :: m
	double precision, DIMENSION(3) :: r,com,cov
	counter = 0 ! used for printing
	t = 0.0 ! the time
	year = 3.15570e7 !number of seconds in a year
	dt = 1 ! time interval (s)
	n = 3 ! number of bodies
	G = 6.67e-11 ! G (Nm**2/kg**2)
	pos = 0.0 ! position coordinates (m)
	v = 0.0 ! velocities (m/s)
	a = 0.0 ! acceleration (m/s*2)
	r = 0.0 ! direction of the gravitational force
	m(1) = 1.99e30 ! mass sun (kg)
	m(2) = 5.97e24 ! mass earth (kg)
	m(2) = 1.8986e27 ! mass jupiter (kg)
	pos(1,2) = 1.496e11 ! earth sun distance (m)
	pos(1,3) = 778547200e3 ! juptier sun distance (m)
	v(2,2) = 29.8e3 ! earth velocity (m/s)
	v(2,3) = 13.07e3 ! jupiter velocity (m/s)
	com = 0.0 ! centre of mass array
	cov = 0.0 ! centre of velocity array
	m_tot = 0.0 ! total mass of system
	e_0 = 0.0 ! total energy at time zero
	e_t = 0.0 ! total energy at later time
	u_temp = 0.0 ! temp potential energy
	u_tot = 0.0 ! total potential energy
	ke_tot = 0.0 ! total kinetic energy
	ke_temp = 0.0 ! temp kinetic energy

	! calculate adjustments for changing frames
	do i = 1, n, 1
		com(1:3) = com(1:3) + m(i)*pos(1:3,i)
		cov(1:3) = cov(1:3) +m(i)*v(1:3,i)
		m_tot = m_tot + m(i)
	end do
	com(1:3) = com(1:3) / m_tot
	cov(1:3) = cov(1:3) / m_tot
	! apply corrections to change frame
 	do i = 1, n, 1
		pos(1:3,i) = pos(1:3,i) - com(1:3)
		v(1:3,i) = v(1:3,i) - cov(1:3)
	end do	

	! first calculate kinetic energy
	do i = 1, n, 1
		vel_temp = sqrt(v(1,i)*v(1,i)+v(2,i)*v(2,i)+v(3,i)*v(3,i))
		ke_temp = 0.5*m(i)*vel_temp*vel_temp
		ke_tot = ke_tot +  ke_temp
	end do
	print *, "total initial kinetic energy =",ke_tot

	! then calculate potential energy
	do i = 1, n-1, 1
		u_temp = 0.0
		do j = i+1, n, 1
			r(1:3) = (pos(1:3,i)-pos(1:3,j))
			rad = sqrt((pos(1,i)-pos(1,j))**2 + (pos(2,i)-pos(2,j))**2 + (pos(3,i)-pos(3,j))**2)
			u_temp = u_temp + (G*m(i)*m(j))/rad
		end do
		u_tot = u_tot + u_temp
	end do

	! total initial energy is ke_tot + u_tot
	e_0 = ke_tot + u_tot

	! main loop, loops over time
	do while (t < 12*year)
		!in the first loop (nested) calculated the acceleration on each body due to each other body
		a = 0.0
		do i = 1,n
			do j = 1,n
				if (j == i) CYCLE
				r(1:3) = (pos(1:3,i)-pos(1:3,j)) ! vector between objects i and j 
				rad = sqrt((pos(1,i)-pos(1,j))**2 + (pos(2,i)-pos(2,j))**2 + (pos(3,i)-pos(3,j))**2) ! length of r
				a(1:3,i) = a(1:3,i) - ((G * m(j))/rad**3)*r(1:3) ! acceleration ma = Gmm/r^3 * rvector
			end do
		end do
		! once the acceleration on each body is calculated the positions are updated
		do i = 1, n, 1
			pos(1:3,i) = pos(1:3,i) + v(1:3,i)*dt ! move the objects
			v(1:3,i) = v(1:3,i) + a(1:3,i)*dt ! update the objects velocities
			! counter loop to check accuracy and print data
			counter = counter + 1
			if ( counter >= 10000000 ) then
				ke_tot = 0.0
				u_tot = 0.0
				! find kinetic energy
				do k = 1, n, 1
					vel_temp = sqrt(v(1,k)*v(1,k)+v(2,k)*v(2,k)+v(3,i)*v(3,k))
					ke_temp = 0.5*m(k)*vel_temp*vel_temp
					ke_tot = ke_tot +  ke_temp
				end do
				! find potential energy
				do k = 1, n-1, 1
					u_temp = 0.0
						do l = k+1, n, 1
							r(1:3) = (pos(1:3,k)-pos(1:3,l))
							rad = sqrt((pos(1,k)-pos(1,l))**2 + (pos(2,k)-pos(2,l))**2 + (pos(3,k)-pos(3,l))**2)
							u_temp = u_temp + (G*m(k)*m(l))/rad
						end do
					u_tot = u_tot + u_temp
				end do
				! sum energies	
				e_t = u_tot + ke_tot
				! energy difference
				print *, 1-e_t/e_0
				counter = 0
			end if
		end do
		t = t + dt ! step time
	end do 
end program nbody
