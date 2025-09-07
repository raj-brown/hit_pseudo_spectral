module fft_module
  use, intrinsic :: iso_c_binding
  implicit none
  include 'fftw3.f03'
  private
  public :: fft_init, initial_condition, compute_wavenumbers
  public :: fft_forward_execute, fft_inverse_execute, fft_cleanup
  

  public :: hit_rhs
  public :: ab_time_integrator, RK45_time_integrator, SSP_time_integrator
  public :: vtk_output
  integer::n=0, nx=0, ny=0, nz=0
  integer:: ncube
  real(C_DOUBLE), pointer :: in_ptr(:, :, :)
  complex(C_DOUBLE_COMPLEX), pointer :: out_ptr(:, :, :)

contains
  
  subroutine initial_condition(n, ux, uy, uz)
    real(8), allocatable, dimension(:, :, :), intent(inout) :: ux, uy, uz
    integer, intent(in) :: n
    integer :: i,j,k
    real(8) :: dx, dy, dz

    real(8), parameter :: pi=3.1415927410125732
    real(8), parameter :: pi2=2*pi

    dx = pi2/n
    dy = pi2/n
    dz = pi2/n

    allocate(ux(n, n, n))
    allocate(uy(n, n, n))
    allocate(uz(n, n, n))

    do k=1, n
       do j=1, n
          do i=1, n
             ux(i,j,k) = sin((i-0.5)*dx)*cos((j-0.5)*dy)*cos((k-0.5)*dz)
             uy(i,j,k) = -cos((i-0.5)*dx)*sin((j-0.5)*dy)*cos((k-0.5)*dz)
             uz(i,j,k) = 0.0 !
          enddo
       enddo
    enddo
  end subroutine initial_condition

  subroutine compute_wavenumbers(n, nhp, kx, ky, kz, k2, k2inv, kabs, ind_k2)
    integer, intent(in) :: n
    integer, intent(in) :: nhp
    real(8), allocatable, dimension(:, :, :), intent(inout) :: kx, ky, kz
    real(8), allocatable, dimension(:, :, :), intent(inout) :: k2, k2inv, kabs
    integer, allocatable, dimension(:, :, :), intent(inout):: ind_k2
    integer :: i, j, k

    ! Wavenumbers 
    allocate(kx(nhp, n, n))
    allocate(ky(nhp, n, n))
    allocate(kz(nhp, n, n))

    ! Wavenumber based parametrs
    allocate(k2(nhp, n, n))
    allocate(k2inv(nhp, n, n))
    allocate(kabs(nhp, n, n))
    allocate(ind_k2(nhp, n, n))

    
    do i=1,nhp
       if (i <= n/2) then
          kx(i, :, :) = i-1
       else
          kx(i, :, :) = -n+i-1
       endif
    enddo

    do j=1,n
       if (j <= n/2) then
          ky(:, j, :) = j-1
       else
          ky(:, j, :) = -n+j-1
       endif
    enddo

    do k=1,n
       if (k <= n/2) then
          kz(:, :, k) = k-1
       else
          kz(:, :, k) = -n+k-1
       endif
    enddo


    do k=1,n
       do j=1,n
          do i=1, nhp
             k2(i,j,k) = kx(i,j,k)*kx(i, j, k) + ky(i,j,k)*ky(i, j, k) + kz(i,j,k)*kz(i, j, k)
             ind_k2(i, j, k) = int(sqrt(k2(i, j, k)) + 0.5)
             kabs(i,j,k) = sqrt(k2(i, j, k))
             if (k2(i,j,k)/=0) then
                k2inv(i, j, k) = 1.0/k2(i, j, k)
             else
                k2inv(i, j, k) = 1.0
             endif
          enddo
       enddo
    enddo

  end subroutine compute_wavenumbers



  ! Subroutine for initializing the infrastructure for FFTW
  ! Which includes the creating forward and inverse plan and relevant data structures
  subroutine fft_init(nx_in, ny_in, nz_in, in_array, out_array, plan_forward, plan_inverse)
    type(C_PTR), intent(out) :: plan_forward
    type(C_PTR), intent(out) :: plan_inverse
    integer, intent(in),target :: nx_in
    integer, intent(in) :: ny_in, nz_in
    real(C_DOUBLE), intent(inout), target :: in_array(nx_in, ny_in, nz_in)
    complex(C_DOUBLE_COMPLEX), intent(inout), target :: out_array((nx_in/2) + 1, ny_in, nz_in)
    plan_forward=C_NULL_PTR
    plan_inverse=C_NULL_PTR

    nx = nx_in
    ny = ny_in
    nz = nz_in
    ncube = nx*ny*nz

    in_ptr => in_array
    out_ptr => out_array
    
    plan_forward = fftw_plan_dft_r2c_3d(nx, ny, nz, in_ptr, out_ptr, FFTW_ESTIMATE)
    plan_inverse = fftw_plan_dft_c2r_3d(nx, ny, nz, out_ptr, in_ptr, FFTW_ESTIMATE)
  
    if (c_associated(plan_forward)) then
       print *, "Forward plan success!"
    end if

    if (c_associated(plan_inverse)) then
       print *, "Inverse plan success!"
    end if
    
     print *, "In fft subroutine "
   end subroutine fft_init


  ! Forward FFT
   subroutine fft_forward_execute(data_in_forward, data_out_forward, plan_forward)
     type(C_PTR), intent(in) :: plan_forward
     real(C_DOUBLE), intent(inout), target:: data_in_forward(nx, ny, nz)
     complex(C_DOUBLE_COMPLEX), intent(inout), target :: data_out_forward((nx/2)+1, ny, nz)
     in_ptr => data_in_forward
     out_ptr => data_out_forward
     call fftw_execute_dft_r2c(plan_forward, in_ptr, out_ptr)
 end subroutine fft_forward_execute

  
  ! Inverse FFT 
  subroutine fft_inverse_execute(data_in_inverse, data_out_inverse, plan_inverse)
    type(C_PTR), intent(in) :: plan_inverse
    real(C_DOUBLE), intent(inout), target :: data_out_inverse(nx, ny, nz)
    complex(C_DOUBLE_COMPLEX), intent(inout), target :: data_in_inverse((nx/2)+1, ny, nz)
    out_ptr => data_in_inverse
    in_ptr => data_out_inverse
    call fftw_execute_dft_c2r(plan_inverse, out_ptr, in_ptr)
    data_out_inverse(:, :, :) = data_out_inverse(:, :, :)/ncube
  end subroutine fft_inverse_execute



  ! clean up the fft analogous data structure and overheads
  subroutine fft_cleanup()
  !  if(c_associated(plan_forward)) then
  !     call fftw_destroy_plan(plan_forward)
  !     plan_forward = C_NULL_PTR
  !     print *, "Destroyed Forward Plan!"
  !  end if
  !  if (c_associated(plan_inverse)) then
  !     call fftw_destroy_plan(plan_inverse)
  !     plan_inverse=C_NULL_PTR
  !     print *, "Destroyed Inverse Plan!"
  !  end if
  end subroutine fft_cleanup


  


  
  subroutine vel_x_omega_physical()

  end subroutine vel_x_omega_physical

  subroutine vel_x_omega_forward_spectral()

  end subroutine vel_x_omega_forward_spectral
  
  
  
  subroutine hit_rhs()
    print *, "In semi discrete form"
  end subroutine hit_rhs

  subroutine ab_time_integrator()
    print *, "in time integration"
  end subroutine ab_time_integrator


  subroutine RK45_time_integrator()
    print *, "in time integration"
  end subroutine RK45_time_integrator

  subroutine SSP_time_integrator()
    print *, "in time integration"
  end subroutine SSP_time_integrator

  subroutine vtk_output(field, dx, dy, dz)
    real(8), allocatable, dimension(:, :, :), intent(in)::field
    character(len=50) :: filename    
    real(8), intent(in):: dx, dy, dz
    integer :: nx, ny, nz, i, j, k
    integer :: unit
    character(len=50) :: ext
    character(len=50) :: output_file

    ext=".vts"
    output_file = filename // ext
    
   
    nx = size(field, 1)
    ny = size(field, 2)
    nz = size(field, 3)

    write(filename, "('field',i10.10,'.plt')") 10
    open(unit=100, file=filename, status="replace")

    write(100, *) 'TITLE = "3D-Volume Data"'
    write(100, *) 'VARIABLES = "X", "Y", "Z", "U"'
    write(100, "('ZONE I=',i4,1x,', J=',i4,1x,', K=',i4,1x,'F=POINT')") nx, ny, nz
    do k=1,nx
      do j=1,ny
         do i=1,nz
            write(100, '(3(f8.5,2x), 1pe15.7, 2x)') &
                 (i-0.5)*dx, (j-0.5)*dy, (k-0.5)*dz, &
                 field(i,j,k)
         enddo
      enddo
   enddo
    
  close(unit)









    close(unit)

    print *, "finished writing file", trim(filename)
    
    
  end subroutine vtk_output
  
  


end module fft_module

