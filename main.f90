! Thi program is to compute the solution of INS Equation for HIT problem
! in a periodic box using pseudo spectral method, with Adam bashforth and Crank Nicholsan
! for time integration.
program main
  use, intrinsic :: iso_c_binding
  use fft_module
  implicit none
  include 'fftw3.f03'

  ! Mesh resolution
  integer, parameter :: n=64
  integer, parameter :: nx=n, ny=n, nz=n
  integer, parameter :: nxyz=nx*ny*nz
  integer, parameter :: nhp=(n/2)+1
  integer, parameter :: ncube=n**3
  character(len=16) :: filename
  type(C_PTR) :: plan_f
  type(C_PTR) :: plan_b
  
  ! Arrays in physical space
  real(8), allocatable, dimension(:, :, :) :: ux, uy, uz
  real(8), allocatable, dimension(:, :, :) :: vortx_ini, vorty_ini, vortz_ini

  !real(8), allocatable, dimension(:, :, :) :: vortx, vorty, vortz
  !real(8), allocatable, dimension(:, :, :) :: velvortx, velvorty, velvortz


  ! Arrays in spectral space
  complex(KIND=C_DOUBLE_COMPLEX), allocatable, dimension(:, :, :) :: ux_hat_1, uy_hat_1, uz_hat_1
  !complex(KIND=C_DOUBLE_COMPLEX), allocatable, dimension(:, :, :) :: ux_hat2, uy_hat2, uz_hat2
   complex(KIND=C_DOUBLE_COMPLEX), allocatable, dimension(:, :, :) :: omegax_hat, omegay_hat, omegaz_hat !vorticity: ω
  !complex(KIND=C_DOUBLE_COMPLEX), allocatable, dimension(:, :, :) :: vomegax_hat, vomegay_hat, vomegaz_hat ! u × ω
  !complex(KIND=C_DOUBLE_COMPLEX), allocatable, dimension(:, :, :) :: rhsx, rhsy, rhsz, rhsx0, rhsy0, rhsz0    ! RHS = F(u)
  
  ! Wavenumber related arrays
  real(8), allocatable, dimension(:, :, :) :: kx, ky, kz, k2, k2inv, kabs, ksqr
  integer, allocatable, dimension(:, :, :) :: id_absk
  real(8), parameter :: pi=3.1415927410125732
  real(8), parameter :: pi2=2*pi

  ! Kinematic viscosity
  real(8), parameter :: visc = 1/300.
  real(8) :: dx, dy, dz
  integer :: istep

  ! Number of timesteps
  integer, parameter :: nsteps=1

  ! output to be written
  integer, parameter :: ndump=50
  ! Output to ASCII files every ... steps
  integer, parameter :: noutput=50
  ! \Delta
  real(8), parameter :: dt=0.005

  !Spectral truncation  for dealiasing
  real(8) :: klimit=n/3.

  ! Benchnark timing and counters
  real :: start, finish
  integer :: i, j, k

    
  dx=pi2/nx
  dy=pi2/ny
  dz=pi2/nz

  allocate(vortx_ini(nx, ny, nz))
  allocate(vorty_ini(nx, ny, nz))
  allocate(vortz_ini(nx, ny, nz))
  allocate(ux_hat_1(nhp, ny, nz))
  allocate(uy_hat_1(nx, ny, nz))
  allocate(uz_hat_1(nx, ny, nz))

  allocate(omegax_hat(nx, ny, nz))
  allocate(omegay_hat(nx, ny, nz))
  allocate(omegaz_hat(nx, ny, nz))


  ! Initialize velocity field according to the Taylor-Green Vortex case
  call initial_condition(n, ux, uy, uz)
  print *, 'size in physical space', size(ux), size(uy), size(uz)
   ! initial vorticity 
  do k=1, n
     do j=1, n
        do i=1, n
           vortx_ini(i,j,k) = -cos((i-0.5)*dx)*sin((j-0.5)*dy)*sin((k-0.5)*dz) 
           vorty_ini(i,j,k) = -sin((i-0.5)*dx)*cos((j-0.5)*dy)*sin((k-0.5)*dz)
           vortz_ini(i,j,k) = sin((i-0.5)*dx)*sin((j-0.5)*dy)*cos((k-0.5)*dz)+sin((i-0.5)*dx)*sin((j-0.5)*dy)*cos((k-0.5)*dz) 
        enddo
     enddo
  enddo

  print *, "analytical vort:", vorty_ini(1:2, 1:2, 1:2)
  
  !call vtk_output(vortx_ini, dx, dy, dz)
 
  call compute_wavenumbers(n, nhp, kx, ky, kz, k2, k2inv, kabs, id_absk)
  call fft_init(nhp, ny, nz, uy, uy_hat_1, plan_f, plan_b)

  call fft_forward_execute(ux, ux_hat_1, plan_f)
  call fft_forward_execute(uy, uy_hat_1, plan_f)
  call fft_forward_execute(uz, uz_hat_1, plan_f)

  do k=1,nz
     do j=1,ny
        do i=1, nhp
           omegax_hat(i,j,k) = cmplx(0,1)*(ky(i,j,k)*uz_hat_1(i,j,k)-kz(i,j,k)*uy_hat_1(i,j,k))
           omegay_hat(i,j,k) = cmplx(0,1)*(kz(i,j,k)*ux_hat_1(i,j,k)-kx(i,j,k)*uz_hat_1(i,j,k))
           omegaz_hat(i,j,k) = cmplx(0,1)*(kx(i,j,k)*uy_hat_1(i,j,k)-ky(i,j,k)*ux_hat_1(i,j,k))
        end do
     end do
  end do

     call fft_inverse_execute(omegax_hat, ux, plan_b)

     print *, 'rec1', omegax_hat(1:2, 1:2, 1:2)
     !print *, 'rec2', uz(1:2, 1:2, 1:2)



  
  !------------------------------------------------------------
  ! Driver loop
  !------------------------------------------------------------
  
  call cpu_time(start)

  do istep=1, nsteps
     
     write(*,"(I10,1X,1pe14.6)") istep, istep*dt
     call fft_forward_execute(ux, ux_hat_1, plan_f)
     call fft_forward_execute(uy, uy_hat_1, plan_f)
     call fft_forward_execute(uz, uz_hat_1, plan_f)


     ! call fft_inverse_execute_half_plane(ux_hat_1, ux, plan_b_1)
     ! call fft_inverse_execute(uz_hat_1, uz, plan_b_2)
     print *,"kyh", ky(1:2, 1:2, 1:2)

     do k=1, nz
        do j=1, ny
           do i=1, nhp
              omegax_hat(i,j,k) = cmplx(0,1)*(ky(i,j,k)*uz_hat_1(i,j,k)-kz(i,j,k)*uy_hat_1(i,j,k))
              omegay_hat(i,j,k) = cmplx(0,1)*(kz(i,j,k)*ux_hat_1(i,j,k)-kx(i,j,k)*uz_hat_1(i,j,k))
              omegaz_hat(i,j,k) = cmplx(0,1)*(kx(i,j,k)*uy_hat_1(i,j,k)-ky(i,j,k)*ux_hat_1(i,j,k)) 
           end do
        end do
     end do

     call fft_inverse_execute(omegax_hat, ux, plan_b)

     !print *, 'rec1', ux(1:2, 1:2, 1:2)
     !print *, 'rec2', uz(1:2, 1:2, 1:2)

  enddo

  deallocate(ux, uy, uz)
  deallocate(vortx_ini, vorty_ini, vortz_ini)
  
  call cpu_time(finish)
  write(*, "(A, 1X, 1pe14.6)") 'elased_time:',  finish - start

endprogram main

