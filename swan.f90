program main
    use helper
    use timestep
    implicit none

    real(dp) :: dx, dy, omega, A
    integer :: nx, ny, nt=0
    real(dp), dimension(:, :), allocatable :: h, u, v, d

    real(dp) :: t_current = 0d0, t_total, dt

    print *, '|----------- SWAN ------------|'
    print *, '|---- by Duan Xu from SJTU----|'
    print *, '|-----------------------------|'
    read(*, *)

    call get_input(dx, dy, nx, ny, dt, t_total, omega, A) ! basic input parameter
    allocate(h(0:nx+1, 0:ny+1), u(0:nx+1, 0:ny+1), v(0:nx+1, 0:ny+1), d(0:nx+1, 0:ny+1))
    
    call init(nx, ny, h, u, v, d) ! initialization, initial condition and depth
    call bc(nx, ny, h, u, v, t_current, omega, A)   ! implement boundary condition
    call write_output(nx, ny, nt, h, u, v)

    do while (t_current < t_total)
        call diffcalc(nx, ny, dx, dy, dt, h, u, v, d)
        call bc(nx, ny, h, u, v, t_current, omega, A)
        t_current = t_current + dt
        nt = nt + 1
        call write_output(nx, ny, nt, h, u, v)
    end do

    print *, '|-------- end of SWAN --------|'
    read(*, *)
end program main