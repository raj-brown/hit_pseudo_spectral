! Thi program is to compute the solution of INS Equation for HIT problem
! in a periodic box using pseudo spectral method, with Adam bashforth and Crank Nicholsan
! for time integration.
program main
   use, intrinsic :: iso_c_binding
   use fft_module
   implicit none
   include 'fftw3.f03'

   ! Mesh resolution
   integer, parameter :: n = 64
   integer, parameter :: nx = n, ny = n, nz = n
   integer, parameter :: nxyz = nx*ny*nz
   integer, parameter :: nhp = (n/2) + 1
   integer, parameter :: ncube = n**3
   !character(len=16) :: filename
   type(C_PTR) :: plan_f
   type(C_PTR) :: plan_b

   ! Arrays in physical space
   real(8), allocatable, dimension(:, :, :) :: ux, uy, uz
   real(8), allocatable, dimension(:, :, :) :: vortx_ini, vorty_ini, vortz_ini

   real(8), allocatable, dimension(:, :, :) :: vortx, vorty, vortz
   real(8), allocatable, dimension(:, :, :) :: velvortx, velvorty, velvortz

   ! Arrays in spectral space
   complex(KIND=C_DOUBLE_COMPLEX), allocatable, dimension(:, :, :) :: ux_hat, uy_hat, uz_hat
   complex(KIND=C_DOUBLE_COMPLEX), allocatable, dimension(:, :, :) :: ux_hat_1, uy_hat_1, uz_hat_1
   complex(KIND=C_DOUBLE_COMPLEX), allocatable, dimension(:, :, :) :: omegax_hat, omegay_hat, omegaz_hat !vorticity: ω
   complex(KIND=C_DOUBLE_COMPLEX), allocatable, dimension(:, :, :) :: vomegax_hat, vomegay_hat, vomegaz_hat ! u × ω
   complex(KIND=C_DOUBLE_COMPLEX), allocatable, dimension(:, :, :) :: rhsx_prev_hat, rhsy_prev_hat, rhsz_prev_hat  ! RHS = F(u); stage Storage
   complex(KIND=C_DOUBLE_COMPLEX), allocatable, dimension(:, :, :) :: rhsx_hat, rhsy_hat, rhsz_hat    ! RHS = F(u) Running storage

   ! Wavenumber related arrays
   real(8), allocatable, dimension(:, :, :) :: kx, ky, kz, k2, k2inv, kabs
   integer, allocatable, dimension(:, :, :) :: id_absk
   real(8), parameter :: pi = 4.0*atan(1.0)
   real(8), parameter :: pi2 = 2*pi

   ! Kinematic viscosity
   real(8), parameter :: visc = 1/300.
   real(8) :: dx, dy, dz
   integer :: istep

   ! Number of timesteps
   integer, parameter :: nsteps = 2000
   real(8) :: energy, energy_0, norm_energy

   ! output to be written
   integer, parameter :: ndump = 50
   ! Output to ASCII files every ... steps
   integer, parameter :: noutput = 50
   ! \Delta
   real(8), parameter :: dt = 0.005

   !Spectral truncation  for dealiasing
   real(8) :: klimit = nx/3.

   ! Benchnark timing and counters
   real :: start, finish
   integer :: i, j, k

   dx = pi2/nx
   dy = pi2/ny
   dz = pi2/nz

   allocate (vortx_ini(nx, ny, nz))
   allocate (vorty_ini(nx, ny, nz))
   allocate (vortz_ini(nx, ny, nz))

   allocate (vortx(nx, ny, nz))
   allocate (vorty(nx, ny, nz))
   allocate (vortz(nx, ny, nz))

   allocate (velvortx(nx, ny, nz))
   allocate (velvorty(nx, ny, nz))
   allocate (velvortz(nx, ny, nz))

   allocate (ux_hat(nhp, ny, nz))
   allocate (uy_hat(nhp, ny, nz))
   allocate (uz_hat(nhp, ny, nz))

   allocate (ux_hat_1(nhp, ny, nz))
   allocate (uy_hat_1(nhp, ny, nz))
   allocate (uz_hat_1(nhp, ny, nz))

   allocate (omegax_hat(nhp, ny, nz))
   allocate (omegay_hat(nhp, ny, nz))
   allocate (omegaz_hat(nhp, ny, nz))

   allocate (vomegax_hat(nhp, ny, nz))
   allocate (vomegay_hat(nhp, ny, nz))
   allocate (vomegaz_hat(nhp, ny, nz))

   allocate (rhsx_prev_hat(nhp, ny, nz))
   allocate (rhsy_prev_hat(nhp, ny, nz))
   allocate (rhsz_prev_hat(nhp, ny, nz))

   allocate (rhsx_hat(nhp, ny, nz))
   allocate (rhsy_hat(nhp, ny, nz))
   allocate (rhsz_hat(nhp, ny, nz))

   ! Initialize velocity field according to the Taylor-Green Vortex case
   call initial_condition(n, ux, uy, uz)

   ! initial vorticity
   do k = 1, n
      do j = 1, n
         do i = 1, n
            vortx_ini(i, j, k) = -cos((i - 0.5)*dx)*sin((j - 0.5)*dy)*sin((k - 0.5)*dz)
            vorty_ini(i, j, k) = -sin((i - 0.5)*dx)*cos((j - 0.5)*dy)*sin((k - 0.5)*dz)
            vortz_ini(i, j, k) = sin((i - 0.5)*dx)*sin((j - 0.5)*dy)*cos((k - 0.5)*dz) &
                                 + sin((i - 0.5)*dx)*sin((j - 0.5)*dy)*cos((k - 0.5)*dz)
         end do
      end do
   end do

   !call vtk_output(vortx_ini, dx, dy, dz)

   call compute_wavenumbers(n, nhp, kx, ky, kz, k2, k2inv, kabs, id_absk)
   call fft_init(nx, ny, nz, ux, ux_hat, plan_f, plan_b)

   call fft_forward_execute(ux, ux_hat, plan_f)
   call fft_forward_execute(uy, uy_hat, plan_f)
   call fft_forward_execute(uz, uz_hat, plan_f)

   energy_0 = 0.5*sum(ux**2 + uy**2 + uz**2)/nxyz
   print *, minval(id_absk)
   print *, "klimit", klimit
   ! Open text file for writing out energy
   open (unit=10, file='ke.txt', status='replace', action='write')

   !write(10, '(A)') 'i    t    ke'

   !------------------------------------------------------------
   ! Driver loop
   !------------------------------------------------------------

   call cpu_time(start)

   do istep = 1, nsteps
      print '(A, 1X, I4, 1X, A, I4)', 'Time Step', istep, 'Total steps', nsteps
      ux_hat_1 = ux_hat
      uy_hat_1 = uy_hat
      uz_hat_1 = uz_hat

      call fft_inverse_execute(ux_hat_1, ux, plan_b)
      call fft_inverse_execute(uy_hat_1, uy, plan_b)
      call fft_inverse_execute(uz_hat_1, uz, plan_b)

      energy = 0.5*sum(ux**2 + uy**2 + uz**2)/nxyz

      norm_energy = energy/energy_0
      write (10, '(I4, 1X, 1PE14.6, 1X, 1PE14.6)') istep, istep*dt, norm_energy
      write (*, '(I4, 1X, 1PE14.6, 1X, 1PE14.6)') istep, istep*dt, energy_0
    
      do k = 1, nz
         do j = 1, ny
            do i = 1, nhp
               omegax_hat(i, j, k) = cmplx(0, 1)*(ky(i, j, k)*uz_hat(i, j, k) - kz(i, j, k)*uy_hat(i, j, k))
               omegay_hat(i, j, k) = cmplx(0, 1)*(kz(i, j, k)*ux_hat(i, j, k) - kx(i, j, k)*uz_hat(i, j, k))
               omegaz_hat(i, j, k) = cmplx(0, 1)*(kx(i, j, k)*uy_hat(i, j, k) - ky(i, j, k)*ux_hat(i, j, k))
            end do
         end do
      end do

      call fft_inverse_execute(omegax_hat, vortx, plan_b)
      call fft_inverse_execute(omegay_hat, vorty, plan_b)
      call fft_inverse_execute(omegaz_hat, vortz, plan_b)

      do k = 1, nz
         do j = 1, ny
            do i = 1, nx
               velvortx(i, j, k) = uy(i, j, k)*vortz(i, j, k) - uz(i, j, k)*vorty(i, j, k)
               velvorty(i, j, k) = uz(i, j, k)*vortx(i, j, k) - ux(i, j, k)*vortz(i, j, k)
               velvortz(i, j, k) = ux(i, j, k)*vorty(i, j, k) - uy(i, j, k)*vortx(i, j, k)
            end do
         end do
      end do

      !print *, 'rec_curl', velvortx(1:2, 1:2, 1:2)

      call fft_forward_execute(velvortx, vomegax_hat, plan_f)
      call fft_forward_execute(velvorty, vomegay_hat, plan_f)
      call fft_forward_execute(velvortz, vomegaz_hat, plan_f)

      ! Dealias Mode
      do k = 1, nz
         do j = 1, ny
            do i = 1, nhp
               if (id_absk(i, j, k) > klimit) then
                  vomegax_hat(i, j, k) = 0.0
                  vomegay_hat(i, j, k) = 0.0
                  vomegaz_hat(i, j, k) = 0.0
               end if
            end do
         end do
      end do

      ! Compute RHS

      do k = 1, nz
         do j = 1, ny
            do i = 1, nhp
               rhsx_hat(i, j, k) = &
                  vomegax_hat(i, j, k) - &
                  visc*k2(i, j, k)*ux_hat(i, j, k) - &
                  kx(i, j, k)*k2inv(i, j, k)* &
                  (kx(i, j, k)*vomegax_hat(i, j, k) + ky(i, j, k)*vomegay_hat(i, j, k) &
                   + kz(i, j, k)*vomegaz_hat(i, j, k))

               rhsy_hat(i, j, k) = &
                  vomegay_hat(i, j, k) - &
                  visc*k2(i, j, k)*uy_hat(i, j, k) - &
                  ky(i, j, k)*k2inv(i, j, k)* &
                  (kx(i, j, k)*vomegax_hat(i, j, k) + ky(i, j, k)*vomegay_hat(i, j, k) &
                   + kz(i, j, k)*vomegaz_hat(i, j, k))

               rhsz_hat(i, j, k) = &
                  vomegaz_hat(i, j, k) - &
                  visc*k2(i, j, k)*uz_hat(i, j, k) - &
                  kz(i, j, k)*k2inv(i, j, k)* &
                  (kx(i, j, k)*vomegax_hat(i, j, k) + ky(i, j, k)*vomegay_hat(i, j, k) &
                   + kz(i, j, k)*vomegaz_hat(i, j, k))
            end do
         end do
      end do

      ! Stage storage for time integrator at first step

      if (istep == 1) then
         rhsx_prev_hat = vomegax_hat
         rhsy_prev_hat = vomegay_hat
         rhsz_prev_hat = vomegaz_hat
      end if

      ! Time integrate

      call ab_time_integrator(nhp, nx, ny, dt, ux_hat, uy_hat, uz_hat, &
                              rhsx_prev_hat, rhsy_prev_hat, rhsz_prev_hat, &
                              rhsx_hat, rhsy_hat, rhsz_hat)

      ! Stage Storage for AB Integrator
      rhsx_prev_hat = rhsx_hat
      rhsy_prev_hat = rhsy_hat
      rhsz_prev_hat = rhsz_hat

   end do

   call fft_cleanup(plan_f)
   call fft_cleanup(plan_b)
   deallocate (ux, uy, uz)
   deallocate (vortx_ini, vorty_ini, vortz_ini)

   call cpu_time(finish)
   write (*, "(A, 1X, 1pe14.6)") 'elased_time:', finish - start

   close (10)
end program main

