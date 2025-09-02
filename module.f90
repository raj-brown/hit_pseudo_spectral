module fft_module
  use, intrinsic :: iso_c_binding
  include 'fftw3.f03'
  implicit none

contains

  subroutine intial_condition(n, ux, uy, uz, vortx, vorty, vortz)
    real(8), allocatable, dimension(:, :, :), intent(out) :: ux, uy, uz
    real(8), allocatable, dimension(:, :, :), intent(out) :: vortx, vorty, vortz
    integer, parameter, intent(in) :: n, i, j, k
    real(8) :: dx, dy, dz
    print *, "In initial conditons"

    do k=1, n
       do j=1, n
          do i=1, n
           ux(i,j,k) = sin((i-0.5)*dx)*cos((j-0.5)*dy)*cos((k-0.5)*dz)
           uy(i,j,k) = -cos((i-0.5)*dx)*sin((j-0.5)*dy)*cos((k-0.5)*dz)
           uz(i,j,k) = 0
           !vortx(i,j,k)=-cos((i-0.5)*dx)*sin((j-0.5)*dy)*cos((k-0.5)*dz)                                    
           !vorty(i,j,k)=-sin((i-0.5)*dx)*cos((j-0.5)*dy)*sin((k-0.5)*dz)                                    
           !vortz(i,j,k) = sin((i-0.5)*dx)*sin((j-0.5)*dy)*cos((k-0.5)*dz) //&                               
           !               -sin((i-0.5)*dx)*sin((j-0.5)*dy)*cos((k-0.5)*dz)                                  
        enddo
     enddo
  enddo

  end subroutine intial_condition



  subroutine hit_rhs()

  end subroutine hit_rhs

  subroutine hit_time_step()

  end subroutine time_step


end module fft_module

