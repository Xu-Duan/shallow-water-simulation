module helper
    implicit none
    integer, parameter :: dp = kind(1d0)

    public :: init, get_input, write_output, get_init_depth
    private :: inquire_dir
contains

    subroutine init(nx, ny, h, u, v, d)
        integer, intent(in) :: nx, ny
        real(dp), dimension(0:nx+1, 0:ny+1), intent(inout):: h, u, v, d

        ! local variables
        integer :: dunit, i, j
        character(7) :: dir = 'output/'
        logical :: ex
        integer :: status
        ex = inquire_dir(dir)
        if (.not. ex) status = system('mkdir '//dir)

        h = 0d0
        u = 0d0
        v = 0d0

        open(newunit=dunit, file=dir//'depth.txt')
        do i = 1, nx
            do j = 1, ny
                write(dunit, '(ES20.5, $)') d(i, j)
            end do
            write(dunit, *)
        enddo
        close(dunit)
    end subroutine init

    subroutine get_init_depth(nx, ny, d)
        integer, intent(in) :: nx, ny
        real(dp), dimension(0:nx+1, 0:ny+1), intent(inout):: d

        ! local variables
        integer :: i, j, nx_, ny_
        integer :: dunit

        open(newunit=dunit, file='depth.txt', status='old')
        read(dunit, *) nx_, ny_
        if (nx_ /= nx .or. nx_ /= nx) stop 'incompatible between depth file and input file'
        do i = 1, nx
            read(dunit, *) (d(i, j), j=1,ny)
        enddo
        close(dunit)

        ! left boundary condition, linear extrapolation
        do j = 1, ny
            d(0, j) = max(0d0, 2d0*d(1, j) - d(2, j))
        end do
        ! right boundary condition, linear extrapolation
        do j = 1, ny
            d(nx+1, j) = max(0d0, 2d0*d(nx, j) - d(nx-1, j))
        end do
        ! upper boundary condition, linear extrapolation
        do i = 1, nx
            d(i, ny+1) = max(0d0, 2d0*d(i, ny) - d(i, ny-1))
        end do
        ! down boundary condition, linear extrapolation
        do i = 1, nx
            d(i, 0) = max(0d0, 2d0*d(i, 1) - d(i, 2))
        end do
    end subroutine get_init_depth

    subroutine get_input(dx, dy, nx, ny, dt, t_total, odt, omega, A, imethod)
        real(dp), intent(inout) :: dx, dy, dt, t_total, odt, omega, A
        integer, intent(inout) :: nx, ny, imethod

        ! local variables
        integer :: inunit

        open(newunit=inunit, file='input.txt', status='old')
        read(inunit, *) dx
        read(inunit, *) dy
        read(inunit, *) nx
        read(inunit, *) ny
        read(inunit, *) dt
        read(inunit, *) odt
        read(inunit, *) t_total
        read(inunit, *) omega
        read(inunit, *) A
        read(inunit, *) imethod
        close(inunit)
    end subroutine get_input

    subroutine write_output(nx, ny, dt, t_current, odt, t_total, h, u, v, ioutput)
        integer, intent(in) :: nx, ny
        real(dp), intent(in) :: dt, t_current, odt, t_total
        real(dp), dimension(0:nx+1, 0:ny+1), intent(in) :: h, u, v
        integer, intent(inout) :: ioutput

        ! local variables
        integer :: unit_h, unit_u, unit_v
        character(7) :: dir = 'output/'
        character(13) :: file_h, file_u, file_v
        integer :: i, j

        if (ioutput * odt <= t_current) then
            ! write output
            write(file_h, '(a, I5.5, a)') 'h_', ioutput, '.out'
            write(file_u, '(a, I5.5, a)') 'u_', ioutput, '.out'
            write(file_v, '(a, I5.5, a)') 'v_', ioutput, '.out'
            open(newunit=unit_h, file=dir//trim(file_h))
            open(newunit=unit_u, file=dir//trim(file_u))
            open(newunit=unit_v, file=dir//trim(file_v))
            do i = 1, nx
                do j = 1, ny
                    write(unit_h, '(ES20.5e3, $)') h(i, j)
                    write(unit_u, '(ES20.5e3, $)') u(i, j)
                    write(unit_v, '(ES20.5e3, $)') v(i, j)
                end do
                write(unit_h, *)
                write(unit_u, *)
                write(unit_v, *)
            enddo
            close(unit_h)
            close(unit_u)
            close(unit_v)
            ioutput = ioutput + 1
        endif
    end subroutine write_output

    subroutine show_progress(t_current, t_total, iprogress)
        integer, intent(inout) :: iprogress
        real(dp), intent(in) :: t_current, t_total
        ! show progress
        if (t_total * iprogress/10d0 < t_current) then
            print '(a, I3, a)', ' calculated  ', iprogress*10, '%'
            iprogress = iprogress + 1
        end if
    end subroutine show_progress

    FUNCTION inquire_dir(dir) result(res)
        CHARACTER(LEN=*), INTENT(IN)     :: dir
        CHARACTER(LEN=*), PARAMETER      :: FN = '/calnskd.a' ! Random
        INTEGER                          :: stat
        LOGICAL                          :: res
        Inquire(file=Trim(dir)//FN, exist=res)
        IF(res) RETURN
        Open(unit=1 , file=Trim(dir)//FN, status='new', action='read', iostat=stat)
        IF (stat .EQ. 0) res = .true.
        close(1, status='delete')
    END FUNCTION inquire_dir
end module helper
