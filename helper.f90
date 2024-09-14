module helper
    implicit none
    integer, parameter :: dp = kind(1d0)

    public :: init, get_input, write_output
    private :: inquire_dir
contains

    subroutine init(nx, ny, h, u, v, d)
        integer, intent(in) :: nx, ny
        real(dp), dimension(0:nx+1, 0:ny+1), intent(inout):: h, u, v, d
        
        ! local variables
        integer :: dunit, i, j, nx_, ny_
        
        h = 0d0
        u = 0d0
        v = 0d0

        open(newunit=dunit, file='depth.txt', status='old')
        read(dunit, *) nx_, ny_
        if (nx_ /= nx .or. nx_ /= nx) stop 'incompatible between depth file and input file'
        do i = 1, nx
            read(dunit, *) (d(i, j), j=1,ny)
        enddo
        close(dunit)
    end subroutine init

    subroutine get_input(dx, dy, nx, ny, dt, t_total, omega, A)
        real(dp), intent(inout) :: dx, dy, dt, t_total, omega, A
        integer, intent(inout) :: nx, ny

        ! local variables
        integer :: inunit

        open(newunit=inunit, file='input.txt', status='old')
        read(inunit, *) dx
        read(inunit, *) dy
        read(inunit, *) nx
        read(inunit, *) ny
        read(inunit, *) dt
        read(inunit, *) t_total
        read(inunit, *) omega
        read(inunit, *) A
        close(inunit)
    end subroutine get_input

    subroutine write_output(nx, ny, nt, h, u, v)
        integer, intent(in) :: nx, ny, nt
        real(dp), dimension(0:nx+1, 0:ny+1), intent(in) :: h, u, v

        ! local variables
        integer :: unit_h, unit_u, unit_v
        character(7) :: dir = 'output/'
        character(13) :: file_h, file_u, file_v

        logical :: tmp
        integer :: status
        tmp = inquire_dir(dir)
        if (.not. tmp) status = system('mkdir '//dir)

        write(file_h, '(a, I5.5, a)') 'h_', nt, '.out'
        write(file_u, '(a, I5.5, a)') 'u_', nt, '.out'
        write(file_v, '(a, I5.5, a)') 'v_', nt, '.out'
        open(newunit=unit_h, file=dir//trim(file_h))
        open(newunit=unit_u, file=dir//trim(file_u))
        open(newunit=unit_v, file=dir//trim(file_v))

        close(unit_h)
        close(unit_u)
        close(unit_v)
    end subroutine write_output

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
