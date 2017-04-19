program main

use mtmod

implicit none

real*8, dimension(:), allocatable		:: m
real*8, dimension(:,:), allocatable	    :: pos, g
real*8                                  :: L, ener, tau, ener2, max_bb, max_cm
real*8                                  :: bbrat, T0, dist1, amaitu, hasi
integer                                 :: N_part, N_links, MCS, N_meas
integer                                 :: i, acc, rej, j, clock, seed

! CONSTANS
!************************************************************
real*8, parameter				       :: hbar = 1.0      ,&
                                          pi = acos(-1.0)
!************************************************************


!Reading parameters
!************************************************************
open(unit=188,file='pathintegral_input.txt',status='old',action='read')
read(unit=188,fmt=*) T0, MCS, max_bb, max_cm, bbrat, n_meas
!************************************************************


!Reading initial configuration
!************************************************************
open(unit=123,file='initial_configuration.dat',status='old',action='read')
open(unit=124,file='interaction.dat',status='old',action='read')
read(unit=123,fmt=*) N_part, N_links, L

tau = 1.0/(T0*N_links)
allocate(m(N_part),pos(N_part,N_links),g(N_part,N_part))

read(unit=123, fmt=*) m(:)

do i = 1, N_links
    read(unit=123,fmt=*) pos(:,i)
enddo
print*, "Number of particles", N_part
print*, "Number of beads", N_links
do i = 1, N_part
    read(unit=124,fmt=*) g(i,:)
enddo

close(unit=123)
close(unit=124)
!****************************************************************



!MONTE CARLO LOOP
!****************************************************************
open(unit=155,file="energy.dat",status="replace",action="write")
open(unit=156,file="position.dat",status="replace",action="write")

call energy(pos,g,m,tau,ener)

acc = 0
rej = 0

CALL SYSTEM_CLOCK(COUNT=clock)
seed = clock + 37
call sgrnd(seed)
!call sgrnd(10)
call cpu_time(hasi)
do i = 1, MCS

    call update(pos,g,m, max_bb, max_cm, bbrat,tau,acc,rej)
    call energy(pos,g,m,tau,ener)

    if (real(i,8)/real(N_meas,8) == real(i/N_meas,8)) then
        write(unit=156,fmt='(i10,100f20.10)') i, pos(1,:), pos(2,:), pos(3,:)
        write(unit=155,fmt='(i10,f20.10)') i , ener
    endif

enddo
call cpu_time(amaitu)
print*, "time", amaitu-hasi
close(unit=155)
close(unit=156)

!****************************************************************

!SAVE LAST CONFIGURATION
!****************************************************************
open(unit=144,file="positions_final2.dat",status="replace",action="write")
write(unit=144,fmt='(2i10,f20.10)') N_part, N_links, L
write(unit=144,fmt='(100f20.10)') m(:)

do j = 1, N_links
    write(unit=144,fmt='(100f20.10)') pos(:,j)
enddo
close(unit=144)
!****************************************************************

print*, "Accepted moves",acc, "Rejected moves",rej

!****************************************************************
!****************************************************************
!****************************************************************
contains

!********************************************
!					UPDATE
!********************************************
subroutine update(pos,g,m,max_bb,max_cm,bbrat,tau,acc,rej)
    real*8, dimension(:,:), intent(inout)		:: pos
    real*8, dimension(:,:), intent(in)			:: g
    real*8, dimension(:), intent(in)			:: m
    real*8, intent(in)							:: max_bb, max_cm, bbrat, tau
    integer, intent(inout)                      :: acc, rej
    real*8										:: random, random2, V0, V
    real*8, dimension(size(pos,1),size(pos,2))	:: pos_new
    integer										:: N_part, N_links, i, j, k, k2

    N_part  = size(pos,1)
    N_links = size(pos,2)

    random = grnd()

    if ( random < bbrat ) then


    	do k = 1, N_links*N_part
!            print*, "****PARTICLE POSITION BEFORE*****"
 !           do k2 = 1, N_links
  !              print*, pos(:,k2)
   !         enddo
            random2 = grnd()
        	i = int(random2*N_part + 1.0) !Choose particle
            random2 = grnd()
            j = int(random2*N_links + 1.0) !Choose beed
!            print*, "selected bead", i, j

    		V0 = density_dif_bed(pos,g,m,j,i,tau)
 !           print*, "V0=", V0
    		random2 = grnd()
    		random2 = (random2 - 0.5) * 2.0 * max_bb

    		pos_new = pos
    		pos_new(i,j) = pos_new(i,j) + random2
  !          print*, "****PARTICLE POSITION AFTER*****"
    !        do k2 = 1, N_links
   !             print*, pos_new(:,k2)
     !       enddo
    		!PBC
    		if ( pos_new(i,j) < 0.0 ) then
    			pos_new(i,j) = pos_new(i,j) + L
    		elseif ( pos_new(i,j) > L ) then
    			pos_new(i,j) = pos_new(i,j) - L
    		endif

    		V = density_dif_bed(pos_new,g,m,j,i,tau)
      !      print*, "V", V
    		random2 = grnd()
       !     print*, "Random", random2, V/V0


    		!ACCEPTANCE
    		if ( random2 < V/V0 ) then
    			pos = pos_new
                acc = acc + 1
        !        print*, "ACEPTED"
    		else
                rej = rej + 1
         !       print*, "REJECTED"
    		endif

    	enddo

    else

    	do i = 1, N_part

    		V0 = density_dif_cm(pos,g,m,i,tau)

    		random2 = grnd()
    		random2 = (random2 - 0.5) * 2.0 * max_cm

    		pos_new = pos
    		pos_new(i,:) = pos_new(i,:) + random2

    		!PBC
    		where ( pos_new < 0.0)
    			pos_new = pos_new + L
    		elsewhere ( pos_new > L)
    			pos_new = pos_new - L
    		endwhere

    		V = density_dif_cm(pos_new,g,m,i,tau)

    		!ACCEPTANCE
    		random2 = grnd()
    		if (random2 < V/V0) then
    			pos = pos_new
                acc = acc + 1
    		else
                rej = rej + 1
    		endif

    	enddo

    endif

end subroutine


!******************************************
!	BEAD-BEAD DENSITY MATRIX DIFFERENCE
!******************************************
function density_dif_bed(pos,g,m,j,i,tau) result(ans)
real*8, dimension(:,:), intent(in)		:: pos, g
real*8, dimension(:), intent(in)		:: m
real*8, intent(in)						:: tau
integer, intent(in)						:: i, j
real*8									:: ans, single, rel, xrel1_sp, xrel2_sp, x_lehen1, x_lehen2
real*8									:: xrel1, xrel2, xrel3, xrel4, xrel5, xrel6, nu1, nu2, g1, g2
integer									:: N, N_links, k, j2, j3, i2, i3

N = size(pos,1)
N_links = size(pos,2)

!CHAIN
j2 = j + 1 - N_links*int(real(j)/real(N_links))
j3 = j - 1 + N_links*int(abs(j-1-N_links)/real(N_links))

xrel1_sp = pos(i,j) - pos(i,j2)
xrel1_sp = xrel1_sp - L*nint(xrel1_sp/L)
xrel2_sp = pos(i,j) - pos(i,j3)
xrel2_sp = xrel2_sp - L*nint(xrel2_sp/L)
!print*, "xrel",xrel1_sp,xrel2_sp
single = sp_dens(xrel1_sp,m(i))*sp_dens(xrel2_sp,m(i))
!print*, "single", single, sp_dens(pos(i,j),pos(i,j+1),m(i)/2.0)
!BEDS
i2 = i + 1 - N*int(real(i)/real(N))
i3 = i - 1 + N*int(abs(i-1-N)/real(N))

rel = 1.0d0
do k = 1, N
	if (k /= i) then

		nu1 = m(i)*m(k)/(m(i) + m(k))
		g1 = g(i,k)

		xrel1 = pos(i,j) - pos(k,j)
		xrel2 = pos(i,j2) - pos(k,j2)
		call minimum_image(xrel1,xrel2,L)
		rel = rel_dens(xrel1,xrel2,nu1,g1)*rel

		xrel1 = pos(i,j) - pos(k,j)
		xrel3 = pos(i,j3) - pos(k,j3)
		call minimum_image(xrel1,xrel3,L)
		rel = rel_dens(xrel1,xrel3,nu1,g1)*rel

	endif
enddo

ans = single*rel
end function

!*******************************************
!		CENTER OF MASS MATRIX DIFFERENCE
!*******************************************

function density_dif_cm(pos,g,m,i,tau) result(ans)
real*8, dimension(:,:), intent(in)		:: pos, g
real*8, dimension(:), intent(in)		:: m
integer, intent(in)						:: i
real*8, intent(in)						:: tau
real*8									:: ans, xrel1, xrel2, xrel4, xrel5, nu1, nu2
integer									:: N_part, N_links, j, j2, i2, i3, k

N_part = size(pos,1)
N_links = size(pos,2)

ans = 1.0

do j = 1, N_links

	j2 = j + 1 - N_links*int(real(j)/real(N_links))

	do k = 1, N_part

		if (k /= i) then

			xrel1 = pos(i,j) - pos(k,j)
			xrel2 = pos(i,j2) - pos(k,j2)
			call minimum_image(xrel1,xrel2,L)

			nu1 = m(i)*m(k)/(m(i) + m(k))

			ans = ans * rel_dens(xrel1,xrel2,nu1,g(i,k))
		endif

	enddo

enddo

end function

!*******************************************
!		SINGLE PARTICLE DENSITY MATRIX
!*******************************************
function sp_dens(xrel,m) result(ans)
real*8, intent(in)		:: xrel, m
real*8					:: ans

ans = exp(-m*xrel**2/(2.0*tau*hbar**2))*sqrt(m/(2.0d0*pi*tau*hbar**2))

!ans = (m/(2.0*pi*tau*hbar**2))**(1.0/2.0)*exp(-m*xrel**2/(2.0*tau*hbar**2))
!print*, "single", ans, xrel
end function

!*******************************************
!			REL DENSITY MATRIX
!*******************************************
function rel_dens(x1,x2,nu,g) result(ans)
real*8, intent(in)		:: x1, x2, nu, g
real*8					:: ans, u

u = nu*(abs(x2)+abs(x1)+g*tau)/sqrt(2.0d0*nu*tau*hbar**2)
ans = 1.0d0 - exp(-nu*(x1*x2 + abs(x1*x2))/(tau*hbar)) &
    *sqrt(pi*nu*tau/2.0d0)*(g/hbar)*f(u)
!print*, "rel", x1, x2, exp(-nu*(x1*x2 + abs(x1*x2))/(tau*hbar))
end function

!********************************************
!				ENERGY MEASURE
!********************************************
subroutine energy(pos,g,m,tau,ener)
    real*8, dimension(:,:), intent(in)			:: g, pos
    real*8, dimension(:), intent(in)			:: m
    real*8, intent(in)							:: tau
    real*8, intent(out)					    	:: ener
    integer										:: i,j,N_links, N_part, k, j2
    real*8										:: sum1, sum2, a1, a2, a3, a4, u, sp
    real*8										:: nu1, xrel1, xrel2, dp, xrel
    N_part  = size(pos,1)
    N_links = size(pos,2)

    ener = 0.0


    sp = 0.0d0
	dp = 0.0d0
    do i = 1, N_part



    	do j = 1, N_links
    		j2 = j + 1 - N_links*int(real(j)/real(N_links))
    		xrel = pos(i,j) - pos(i,j2)
    		xrel = xrel - L*nint(xrel/L)

    		sp = sp + m(i)*xrel**2/(2.0d0*hbar**2*tau**2)
    		!print*, "sp",j, m(i)*xrel**2/(2.0d0*hbar**2*tau**2), xrel
    		sum1 = 0.0d0
    		sum2 = 0.0d0
    		do k = i+1, N_part

    			xrel1 = pos(i,j) - pos(k,j)
    			xrel2 = pos(i,j2) - pos(k,j2)
    			!print*, i,j,k, pos(i,j), pos(k,j), pos(i,j2), pos(k,j2)
    			call minimum_image(xrel1,xrel2,L)
                !print*, i, j, k, xrel1, xrel2
    			nu1 = m(i)*m(k)/(m(i) + m(k))

    			u = sqrt(nu1/(2.0d0*tau))*(abs(xrel1) + abs(xrel2) + g(i,k)*tau)/hbar
    			!print*, "u", u
    			a1 = nu1*(xrel1*xrel2 + abs(xrel1*xrel2))/(hbar**2*tau**(3.0d0/2.0d0))*f(u)
    			a2 = f(u)/(2.0d0*sqrt(tau))
    			a3 = -nu1*(2.0d0*u*f(u)-2.0d0/pi)*(abs(xrel1)+abs(xrel2)-g(i,k)*tau)&
    				/(hbar*tau*sqrt(8.0d0*nu1))
    			a4 = (g(i,k)/hbar)*sqrt(pi*nu1*tau/2.0d0)*f(u) - &
    				exp(nu1*(xrel1*xrel2 + abs(xrel1*xrel2))/(tau*hbar**2))
    			!print*, a4, (g(i,k)/hbar)*sqrt(pi*nu1*tau/2.0d0)*erfc(u)*exp(u**2), exp(nu1*(xrel1*xrel2 + abs(xrel1*xrel2))/(tau*hbar**2))
    			sum1 = sum1 + (g(i,k)/hbar)*sqrt(pi*nu1/2.0d0)*(a1+a2+a3)/a4
    			sum2 = sum2 + log(1.0d0 - (g(i,k)/hbar)*sqrt(pi*nu1*tau/2.0d0)*f(u)&
    				*exp(-nu1*(xrel1*xrel2 + abs(xrel1*xrel2))/(tau*hbar**2)))
				!print*, i, j, k, 1.0d0 - (g(i,k)/hbar)*sqrt(pi*nu1*tau/2.0d0)*f(u)&
    			!	*exp(-nu1*(xrel1*xrel2 + abs(xrel1*xrel2))/(tau*hbar**2)), xrel1, xrel2
				!print*, i, j, k, (g(i,k)/hbar)*sqrt(pi*nu1/2.0d0)
    		enddo

    		dp = dp + exp(sum2)*sum1
			!print*, "dp", i, j, sum1, sum2, dp
    		!print*, i, j, exp(sum2)*sum1, sum1, sum2
    	enddo
        !ener = 1.0d0/(2.0d0*tau) - sp/real(N_links,8) + dp/real(N_links,8) + ener
!        print*, "part",i,1.0d0/(2.0d0*tau) - sp/(real(N_part,8)*real(N_links,8)) + dp/(real(N_part,8)*real(N_links,8))
        !print*, dp/real(N_links,8)

    enddo

    !ener = 1.0d0/(2.0d0*tau) - sp/(real(N_part,8)*real(N_links,8)) + dp/(real(N_part,8)*real(N_links,8))
    !ener = ener/real(N_part,8)
	ener = 1.0d0/(2.0d0*tau) - sp/real(N_part*N_links,8) + dp/real(N_links,8)
	!print*, "single_part = ", sp/real(N_part*N_links,8)
	!print*, "double_part = ", dp/real(N_part*N_links,8)
    !print*, "energy:",ener
end subroutine

subroutine energy2(pos,ener)
real*8, dimension(:,:), intent(in)      :: pos
real*8, intent(out)                     :: ener
integer                                 :: M, j

M = size(pos,2)
ener = 0.0d0
do j = 1, M-1
    ener = (pos(1,j) - pos(1,j+1))**2 + ener
enddo
ener = ener + (pos(1,1)-pos(1,M))**2
ener = ener/(real(M,8)*4.0*tau**2)
end subroutine

!********************************************
!				minimum_image
!********************************************
! Funtzio honen berezitasuna da partikula bereko
! bi bola ditugunean eta bakarrik aplikatzen diogunean
! pbc bateri ez du hori egiten eta zegoen lekuan uzten du

subroutine minimum_image(x1,x2,L)
real*8, intent(in)		:: L
real*8, intent(inout)	:: x1, x2
real*8                  :: x1_old, x2_old
x1_old = x1
x2_old = x2
x1 = x1 - L*nint(x1/L)
x2 = x2 - L*nint(x2/L)

if (x1*x2 < 0.0d0 .and. (x1_old > L/2.0d0 .or. x2_old > L/2.0d0)) then
    x1 = x1_old - L*nint(x1_old/L)
    x2 = x2_old - L*nint(x2_old/L)
endif
end subroutine

function f(u) result(ans)
real*8, intent(in)  :: u
real*8              :: ans

ans = 2.0d0/(sqrt(pi)*(u+sqrt(u**2+4.0d0/pi)))
end function
end program
