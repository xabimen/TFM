program radial

integer                                 :: N_part, N_links, i, j, N_his
real*16, dimension(:,:), allocatable    :: pos
real*16, dimension(:), allocatable      :: gr
real*16                                 :: L
N_part = 10
N_links = 20
N_his = 20
L = 10.0
allocate(pos(N_part,N_links),gr(N_his))

open(unit=123,file='positions_pbc.dat',status='old',action='read')

do j = 1, 20
    read(unit=123,fmt=*) pos(:,j)
enddo


call radial_distribution(pos,L,N_his,gr)

open(unit=199, file="radial.dat", status="replace", action="write")


do i = 0, N_his
	write(unit=199, fmt='(3f20.10)') i*L/(2.0d0*real(N_his,16)), gr(i)
enddo

close(unit=199)

contains

subroutine radial_distribution(pos,L,N_his,gr)
real*16, dimension(:,:), intent(in)	    :: pos
real*16, intent(in)					    :: L
integer, intent(in)					    :: N_his
real*16, dimension(:), intent(out)		:: gr
real*16								    :: rel, x, dx
integer								    :: N_part, N_links, i1, i2, j, k

N_part = size(pos,1)
N_links = size(pos,2)

dx = L/(2.0d0*real(N_his,16))

gr = 0.0d0

!j = 1
do j = 1, N_links
do i1 = 1, N_part
	do i2 = 1, N_part

		rel = pos(i1,j) - pos(i2,j)
		rel = abs(rel - L*nint(rel/L))
        print*, rel

		do k = 1, N_his
			x = (k-1)*L/(2.0d0*real(N_his,16))

			if (rel > x .and. rel < x + dx) then
				gr(k) = gr(k) + 1.0d0

			endif
		enddo

	enddo
enddo
enddo

gr = gr*N_his/real(N_part**2*N_links,16)




end subroutine

end program
