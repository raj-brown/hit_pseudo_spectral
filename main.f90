! Thi program is to compute the solution of INS Equation for HIT problem
! in a periodic box using pseudo spectral method, with Adam bashforth and Crank Nicholsan
! for time integration.
program main
  use, intrinsic :: iso_c_binding
  use fft_module
  implicit none

  ! Mesh resolution
  integer, parameter :: n=64, nx=n, ny=n, nz=n
  integer, parameter :: nxnynz=nx*ny*nz
  integer, parameter :: nhp1=n/2+1
  integer, parameter :: ncube=n**3
  character(len=16) :: filename
  
  ! Arrays in physical space
  real(8), allocatable, dimension(:, :, :) :: ux, uy, uz
  real(8), allocatable, dimension(:, :, :) :: vortx, vorty, vortz
  !real(8), allocatable, dimension(:, :, :) :: velvortx, velvorty, velvortz
  !real(8), allocatable, dimension(:, :, :) :: uxg, uyg, uzg

  ! Arrays in spectral space
  !complex(KIND=C_DOUBLE_COMPLEX), allocatable, dimension(:, :, :) :: ux_hat, uy_hat, uz_hat
  !complex(KIND=C_DOUBLE_COMPLEX), allocatable, dimension(:, :, :) :: ux_hat2, uy_hat2, uz_hat2
  !complex(KIND=C_DOUBLE_COMPLEX), allocatable, dimension(:, :, :) :: omegax_hat, omegay_hat, omegaz_hat !vorticity: ω
  !complex(KIND=C_DOUBLE_COMPLEX), allocatable, dimension(:, :, :) :: vomegax_hat, vomegay_hat, vomegaz_hat ! u × ω
  !complex(KIND=C_DOUBLE_COMPLEX), allocatable, dimension(:, :, :) :: rhsx, rhsy, rhsz, rhsx0, rhsy0, rhsz0    ! RHS = F(u)
  
  ! Wavenumber related arrays
  !real(8), allocatable, dimension(:, :, :) :: kx, ky, kz, k2, k2inv, kabs, ksqr
  !integer, allocatable, dimension(:, :, :) :: ik2
  real(8), parameter :: pi=3.1415927410125732
  real(8), parameter :: pi2=2*pi

  ! Kinematic viscosity
  real(8), parameter :: visc = 1/300.
  real(8) :: dx, dy, dz
  integer :: istep

  ! Number of timesteps
  integer, parameter :: nsteps=2000

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

  !allocate(ux(nx, ny, nz))
  !allocate(uy(nx, ny, nz))
  !allocate(uz(nx, ny, nz))

  ! Initialize velocity field according to the Taylor-Green Vortex case

  call initial_condition(n, ux, uy, uz, vortx, vorty, vortz)

  !do k=1, n
  !   do j=1, n
  !      do i=1, n
  !         ux(i,j,k) = sin((i-0.5)*dx)*cos((j-0.5)*dy)*cos((k-0.5)*dz)
  !         uy(i,j,k) = -cos((i-0.5)*dx)*sin((j-0.5)*dy)*cos((k-0.5)*dz)
  !         uz(i,j,k) = 0
  !      enddo
  !   enddo
  !enddo

  print *, "ux", ux(1:5, 1:5, 1:5)
  
  !------------------------------------------------------------
  ! Main time loop
  !------------------------------------------------------------
  
  !call cpu_time(start)

  !do istep=1, nsteps
  !   if (mod(istep, 10)==0) then
  !      write(*,"(I7,1X,1pe14.6)") istep, istep*dt
  !   endif
  !enddo ! timeloop

 deallocate(ux, uy, uz)

  
  call cpu_time(finish)
  write(*, "(1X, 1pe14.6)") finish - start

endprogram main

