module timestep
    use helper, only : dp
    implicit none
    real(dp), parameter :: C = 25d0 ! chezy coefficient
    real(dp), parameter :: pi = acos(-1d0)
    real(dp), parameter :: phi = 2d0 * pi / 24 / 3600 ! self-rotating angle velocity rad/s
    real(dp), parameter :: latitude = 31 ! 31 deg North
    real(dp), parameter :: F = 2 * phi * sin(latitude/180d0*pi) ! coriolis coefficient
    real(dp), parameter :: g = 9.81
contains
    subroutine diffcalc(nx, ny, dx, dy, dt, h, u, v, d)
        integer, intent(in) :: nx, ny
        real(dp), intent(in) :: dx, dy, dt
        real(dp), dimension(0:nx+1, 0:ny+1), intent(in):: h, u, v, d

        ! local variables
        integer :: i, j
        real(dp) :: D1X, D2X, D1Y, D2Y, TU1, TU2, MV, THX, SX, TV1, TV2, MU, THY, SY
        ! D1X: the depth plus wave height for u(i+1, j) in x direction (1st eq)
        ! D2X: the depth plus wave height for u(i, j) in x direction (1st eq)
        ! D1Y: the depth plus wave height for v(i, j+1) in y direction (1st eq)
        ! D2Y: the depth plus wave height for v(i, j) in y direction (1st eq)
        ! TU1: the difference of velocity in x direction (2nd eq)
        ! TU2: the difference of velocity in y direction (2nd eq)
        ! MV: the mean value of velocity in y direction (2nd eq)
        ! THX: the difference of wave height in x direction (2nd eq)
        ! SX: the shear force of bottom friction in x direction (2nd eq)
        ! TV1: the difference of velocity in x direction (3rd eq)
        ! TV2: the difference of velocity in y direction (3rd eq)
        ! MU: the mean value of velocity in y direction (3rd eq)
        ! THY: the difference of wave height in y direction (3rd eq)
        ! SY: the shear force of bottom friction in y direction (3rd eq)
        do j = 1,ny
            do i = 1, nx
                ! 1st eq
                if (u(i+1, j) > 0d0) then
                    D1X = d(i, j) + h(i, j)
                else
                    D1X = d(i+1, j) + h(i+1, j)
                endif
                if (u(i, j) > 0d0) then
                    D2X = d(i-1, j) + h(i-1, j)
                else
                    D2X = d(i, j) + h(i, j)
                endif
                if (v(i, j+1) > 0d0) then
                    D1Y = d(i, j) + h(i, j)
                else
                    D1Y = d(i, j+1) + h(i, j+1)
                endif
                if (v(i, j) > 0d0) then
                    D2Y = d(i, j-1) + h(i, j-1)
                else
                    D2Y = d(i, j) + h(i, j)
                endif
                
                ! 2nd eq
                if (u(i, j) > 0d0) then
                    TU1 = u(i, j) - u(i-1, j)
                else
                    TU1 = u(i+1, j) - u(i, j)
                endif
                MV = (v(i-1, j) + v(i+1, j) + v(i, j+1) + v(i, j-1))/4d0
                if (MV > 0d0) then
                    TU2 = u(i, j) - u(i, j-1)
                else
                    TU2 = u(i, j+1) - u(i, j)
                endif
                THX = h(i, j) - h(i-1, j)
                SX = g * u(i, j) * sqrt(u(i, j)**2 + v(i, j)**2) / C**2 / (d(i, j)+h(i, j))

                ! 3rd eq
                MU = (u(i-1, j) + u(i+1, j) + u(i, j+1) + u(i, j-1))/4d0
                if (MU > 0d0) then
                    TV1 = v(i, j) - v(i-1, j)
                else
                    TV1 = v(i+1, j) - v(i, j)
                endif
                if (v(i, j) > 0d0) then
                    TV2 = v(i, j) - v(i, j-1)
                else
                    TV2 = v(i, j+1) - v(i, j)
                endif
                THY = h(i, j) - h(i, j-1)
                SY = g * v(i, j) * sqrt(u(i, j)**2 + v(i, j)**2) / C**2 / (d(i, j)+h(i, j))
            end do
        enddo
    end subroutine diffcalc

    subroutine bc(nx, ny, h, u, v, t_current, omega, A)
        integer, intent(in) :: nx, ny
        real(dp), dimension(0:nx+1, 0:ny+1), intent(inout):: h, u, v
        real(dp), intent(in) :: omega, A, t_current

        ! local varibales
        integer :: i, j

        ! TODO : left boundary condition, input airy wave condition
        do j = 1, ny
            h(0, j) = A * sin(omega * t_current)
            u(0, j) = 2 * A * sin(omega * t_current)
            v(0, j) = 0d0
        end do

        ! right boundary condition, linear extrapolation
        do j = 1, ny
            h(nx+1, j) = 2d0*h(nx, j) - h(nx-1, j)
            u(nx+1, j) = 2d0*u(nx, j) - u(nx-1, j)
            v(nx+1, j) = 2d0*v(nx, j) - v(nx-1, j)
        end do

        ! upper boundary condition, linear extrapolation
        do i = 1, nx
            h(i, ny+1) = 2d0*h(i, ny) - h(i, ny-1)
            u(i, ny+1) = 2d0*u(i, ny) - u(i, ny-1)
            v(i, ny+1) = 2d0*v(i, ny) - v(i, ny-1)
        end do

        ! down boundary condition, linear extrapolation
        do i = 1, nx
            h(i, 0) = 2d0*h(i, 1) - h(i, 2)
            u(i, 0) = 2d0*u(i, 1) - u(i, 2)
            v(i, 0) = 2d0*v(i, 1) - v(i, 2)
        end do

    end subroutine bc
end module timestep