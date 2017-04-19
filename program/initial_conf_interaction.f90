program initial_conf_interaction

use mtmod

implicit none

integer                             ::  N_A, N_B, N_links, i, j
real*8                              ::  g_A, g_B, g_AB, L, rho_A, rho_B, M_A, M_B
real*8, dimension(:,:), allocatable :: pos, g
real*8, dimension(:), allocatable   :: m
!**************************************************

open(unit=123,file='initial_conf_interaction.txt',status='old',action='read')
read(unit=123, fmt=*) N_A, N_B, M_A, M_B, g_A, g_B, g_AB, L, N_links
close(unit=123)


allocate(pos(N_A+N_B,N_links),m(N_A+N_B),g(N_A+N_B,N_A+N_B))

rho_A = reaL(N_A,8)/L
rho_B = real(N_B,8)/L

call initialize_pos(rho_A,rho_B,L,pos)
call initialize_mix(rho_A,rho_B,M_A,M_B,g_A,g_B,g_AB,m,g)

open(unit=124,file='initial_configuration.dat',status='replace',action='write')
open(unit=125,file='interaction.dat',status='replace',action='write')

write(unit=124,fmt='(2i10,f20.10)') N_A+N_B, N_links, L
write(unit=124,fmt='(100f20.10)') m(:)
do j = 1, N_links
    write(unit=124,fmt='(100f20.10)') pos(:,j)
enddo

do i = 1, N_A+N_B
    write(unit=125,fmt='(100f20.10)') g(i,:)
enddo

close(unit=124)
close(unit=125)
contains

!********************************************
!			INITIALIZE POSITION
!********************************************
subroutine initialize_pos(rho_A,rho_B,L,pos)
real*8, intent(in)						:: rho_A, rho_B, L
real*8, dimension(:,:), intent(out)	:: pos
integer									:: N_links, N_part, i

N_part  = size(pos,1)
N_links = size(pos,2)

do i = 1, N_part
	pos(i,:) = (2.0d0*i - 1.0d0)*L/(2.0d0*real(N_part,8))
enddo

end subroutine

!******************************************
!			INITIALIZE MIX
!******************************************

subroutine initialize_mix(rho_A,rho_B,M_A,M_B,g_A,g_B,g_AB,m,g)
real*8, intent(in)						:: rho_A, rho_B, M_A, M_B, g_A, g_B, g_AB
real*8, dimension(:), allocatable		:: m
real*8, dimension(:,:), allocatable	    :: g
real*8									:: x
integer									:: N, N_B, i, k, j

N = size(m)
N_B = rho_B*N/(rho_A+rho_B)

m = M_A
g = 0.0

if ( M_A /= M_B ) then
	do i = 1, N_B

		x = grnd()
		j = int(x*N + 1.0)
		do while(m(j)==M_B)
			x = grnd()
			j = int(x*N + 1.0)
		enddo
		m(j) = M_B
	enddo
endif

do i = 1, N
	do j = i, N

		if (m(i) == M_A .and. m(j) == M_A .and. i /= j) then
			g(i,j) = g_A
			g(j,i) = g_A
		elseif (m(i) == M_B .and. m(j) == M_B .and. i /=j) then
			g(i,j) = g_B
			g(j,i) = g_B
		else
			g(i,j) = g_AB
			g(j,i) = g_AB
		endif

	enddo
enddo


end subroutine

end program
