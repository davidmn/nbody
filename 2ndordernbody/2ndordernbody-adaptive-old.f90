program nbody
	implicit none

	integer :: i,j,n, counter, l, k
	double precision :: G, dt, rad, year, t, e_0, e_t, ke_tot,ke_temp, u_temp, u_tot, m_tot
	double precision :: vel_temp, current_time, criteria, AU, uranus_distance
	! vectors in format (x1,y1,z1), (x2,y2,z2)
	double precision, DIMENSION(3, 7) :: pos,v,a,a2
	double precision, DIMENSION(7) :: m
	double precision, DIMENSION(3) :: r,com,cov
	integer, DIMENSION(7) :: class
	uranus_distance = 0.0
	AU = 1.5e11 ! astronomical unit
	criteria = 0.0
	class = 1 ! start all planets having acc calculated every step, higher class are calculaed less often
	counter = 0 ! used for printing
	t = 0.0 ! the time
	year = 3.15570e7 !number of seconds in a year
	dt = 30 ! time interval (s)
	n = 7 ! number of bodies
	G = 6.67e-11 ! G (Nm**2/kg**2)
	pos = 0.0 ! position coordinates (m)
	v = 0.0 ! velocities (m/s)
	a2 = 0.0 
	a = 0.0 ! acceleration (m/s*2)
	r = 0.0 ! direction of the gravitational force
	m(1) = 1.99e30 ! mass sun (kg)
	m(2) = 4.87e24 ! mass venus(kg)
	m(3) = 5.97e24 ! mass earth (kg)
	m(4) = 6.39e23 ! mass mars (kg)
	m(5) = 1.90e27 ! mass jupiter (kg)
	m(6) = 5.683e27 ! mass saturn (kg)
	m(7) = 8.68e25 ! uranus mass (kg)
	pos(1,2) = 1.089e11 ! venus sun distance (m)
	pos(1,3) = 1.496e11 ! earth sun distance (m)
	pos(1,4) = 2.387e11 ! mars sun distance (m)
	pos(1,5) = 7.835e11 ! juptier sun distance (m)
	pos(1,6) = 1.485e12 ! saturn sun distance (m)
	pos(1,7) = 3.124e12 ! uranus sun distance (m)
	v(2,2) = 3.5e4 ! venus velocity (m/s)
	v(2,3) = 2.98e4 ! earth velocity (m/s)
	v(2,4) = 2.41e4 ! mars velocity (m/s)
	v(2,5) = 1.307e4 ! jupiter velocity (m/s)
	v(2,6) = 9.64e3 ! saturn velocity (m/s)
	V(2,7) = 6.7e3 ! uranus velocity (m/s)
	com = 0.0 ! centre of mass array
	cov = 0.0 ! centre of velocity array
	m_tot = 0.0 ! total mass of system
	e_0 = 0.0 ! total energy at time
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
			u_temp = u_temp - (G*m(i)*m(j))/rad
		end do
		u_tot = u_tot + u_temp
	end do

	! total initial energy is ke_tot + u_tot
	e_0 = ke_tot + u_tot

	ke_tot = 0.0
	u_tot = 0.0

	! kick
	do i = 1,n
		do j = 1,n
			if (j == i) CYCLE
			r(1:3) = (pos(1:3,i)-pos(1:3,j)) ! vector between objects i and j 
			rad = sqrt((pos(1,i)-pos(1,j))**2 + (pos(2,i)-pos(2,j))**2 + (pos(3,i)-pos(3,j))**2) ! length of r
			a(1:3,i) = a(1:3,i) - ((G * m(j))/rad**3)*r(1:3) ! acceleration ma = Gmm/r^3 * rvector
		end do
	end do

	! works out what integer time each body should be in
	do i = 1, n, 1
		criteria = 0.01 / sqrt(a(1,i)*a(1,i) + a(2,i)*a(2,i) + a(3,i))
		print *, criteria
		if ( criteria < 2.0 ) class(i) = 1
		if ( criteria > 2.0 .and. criteria < 4.0 ) class(i) = 2
		if ( criteria > 4.0 .and. criteria < 8.0 ) class(i) = 4
		if ( criteria > 8.0 .and. criteria < 16.0 ) class(i) = 8
		if ( criteria > 16.0 .and. criteria < 32.0 ) class(i) = 16
		if ( criteria > 32.0 ) class(i) = 32		
	end do

	do i = 1, n, 1
		print *, class(i)
	end do

	! main loop, loops over time
	do while (t < 10000*year)
		!in the first loop (nested) calculated the acceleration on each body due to each other body
		a2 = 0.0
		do i = 1,n
			do j = 1,n
				if (j == i) CYCLE
				r(1:3) = (pos(1:3,i)-pos(1:3,j)) ! vector between objects i and j 
				rad = sqrt((pos(1,i)-pos(1,j))**2 + (pos(2,i)-pos(2,j))**2 + (pos(3,i)-pos(3,j))**2) ! length of r
				a2(1:3,i) = a2(1:3,i) - ((G * m(j))/rad**3)*r(1:3) ! acceleration ma = Gmm/r^3 * rvector
			end do
		end do
		! once the acceleration on each body is calculated the positions are updated
		do i = 1, n, 1
			pos(1:3,i) = pos(1:3,i) + v(1:3,i)*dt + (0.5 * a2(1:3,i) * dt*dt) ! move the objects
			v(1:3,i) = v(1:3,i) + 0.5 * (a(1:3,i)+a2(1:3,i)) * dt ! update the objects velocities
			! counter loop to check accuracy and print data
			counter = counter + 1
			if ( counter >= 10000000 ) then
				ke_tot = 0.0
				u_tot = 0.0
				e_t = 0.0
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
							u_temp = u_temp - (G*m(k)*m(l))/rad
						end do
					u_tot = u_tot + u_temp
				end do
				! sum energies	
				e_t = u_tot + ke_tot
				! energy difference
				current_time = t/year
				uranus_distance = (sqrt((pos(1,1)-pos(1,7))**2 + (pos(2,1)-pos(2,7))**2 + (pos(3,1)-pos(3,7))**2))/AU
				print *, "time =", current_time, ((e_t-e_0)/e_0)*100
				print *, "time =", current_time, uranus_distance
				counter = 0
			end if
		end do
		a = a2
		t = t + dt ! step time
	end do 
end program nbody

