program energy_kaka

implicit none

real*8, dimension(3,2)		:: pos
real*8, dimension(3,3)		:: g
real*8, dimension(3)		:: m
real*8						:: tau, ener, L, hbar, pi
integer						:: i

hbar = 1.0d0
pi = acos(-1.0d0)
L = 10.0d0
g = 1.0d0
m = 1.0d0
tau = 1.0d0

pos(1,:) = (/2.0d0,3.0d0/)
pos(2,:) = (/3.0d0,4.0d0/)
pos(3,:) = (/5.0d0,6.0d0/)
do i = 1, 3
	print*, pos(i,:)
enddo

call energy(pos,g,m,tau,ener)

print*, ''
print*, "energy = ", ener

contains

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
