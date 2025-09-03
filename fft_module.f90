module fft_module
  use, intrinsic :: iso_c_binding
  implicit none
    
contains

  subroutine initial_condition(n, ux, uy, uz, vortx, vorty, vortz)
    real(8), allocatable, dimension(:, :, :), intent(out) :: ux, uy, uz
    real(8), allocatable, dimension(:, :, :), intent(out) :: vortx, vorty, vortz
    integer, intent(in) :: n
    integer :: i,j,k
    real(8) :: dx, dy, dz

    real(8), parameter :: pi=3.1415927410125732
    real(8), parameter :: pi2=2*pi

    dx = pi2/n
    dy = pi2/n
    dz = pi2/n

    print *, "In initial conditons"
    print *, "Value of n=", n

    allocate(ux(n, n, n))
    allocate(uy(n, n, n))
    allocate(uz(n, n, n))

    allocate(vortx(n, n, n))
    allocate(vorty(n, n, n))
    allocate(vortz(n, n, n))

    do k=1, n
       do j=1, n
          do i=1, n
           ux(i,j,k) = sin((i-0.5)*dx)*cos((j-0.5)*dy)*cos((k-0.5)*dz)
           uy(i,j,k) = -cos((i-0.5)*dx)*sin((j-0.5)*dy)*cos((k-0.5)*dz)
           uz(i,j,k) = 0.0 !
           !vortx(i,j,k)=-cos((i-0.5)*dx)*sin((j-0.5)*dy)*cos((k-0.5)*dz)                                    
           !vorty(i,j,k)=-sin((i-0.5)*dx)*cos((j-0.5)*dy)*sin((k-0.5)*dz)                                    
           !vortz(i,j,k) = sin((i-0.5)*dx)*sin((j-0.5)*dy)*cos((k-0.5)*dz) //&                               
           !               -sin((i-0.5)*dx)*sin((j-0.5)*dy)*cos((k-0.5)*dz)                                  
        enddo
     enddo
  enddo

  end subroutine initial_condition



  subroutine hit_rhs()
    print *, "In semi discrete form"

  end subroutine hit_rhs

  subroutine hit_time_step()
    print *, "in time integration"
  end subroutine hit_time_step


end module fft_module

