!**************************************************************
!  Copyright Euratom-CEA
!  Authors :
!     Cartier-Michaud Thomas (thomas.cartier-michaud@cea.fr)
!     Ghizzo Alain (Alain.Ghizzo@ijl.nancy-universite.fr)
!     Grandgirard Virginie (virginie.grandgirard@cea.fr)
!     Passeron (chantal.passeron@cea.fr)
!     Sarazin Yanick (yanick.sarazin@cea.fr)
!
!  This code TERESA (Trapped Element REduction in
!  Semi lagrangian Approach) solve Vlasov equation
!  for trapped ions in 4-D F(Ksi,Psi,Kappa,E,t) parametrised
!  by the Energy E and the kappa parameter
!     0 < kappa  <1   for trapped particles
!  Numerical scheme = semi-Lagrangian
!  Reference:
!     G.Depret, X.Garbet, P.Bertrand, A.Ghizzo, PPCF 42 (2000) 949
!
!  This software is governed by the CeCILL-B license
!  under French law and abiding by the rules of distribution
!  of free software.  You can  use, modify and redistribute
!  the software under the terms of the CeCILL-B license as
!  circulated by CEA, CNRS and INRIA at the following URL
!  "http://www.cecill.info".
!**************************************************************
!
! System:
!
!	d/dt \bar{F} - [\bar{Phi}, \bar{F}] + Omega_D E \partial_ksi F = 0
!
!	C2 (Phi-<Phi>) - C1 (rho^2 d^2/dpsi^2 + deltab^2 d^2/dksi^2) Phi
!        =   INTEGRAL { JO F jacob(E) jacob(K) dE dK}/n_eq - 1
!
! With:
!    F(Ksi,Psi,Ka,E)    = Distribution Function
!    Phi(Ksi,Psi)       = Electric Potential
!    Ksi                = phi - q*theta			        (~ precession angle)
!    Psi                = Poloidal Flux		                (~ radial coordinate)
!    Vad1(Ksi,Psi,Ka,E) = + d/dPsi [J0 . Phi] (+OmegaD.E)  ( Advection Field along Ksi )
!    Vad2(Ksi,Psi,Ka,E) = - d/dKsi [J0 . Phi]	                ( Advection Field along Psi )
!    E                  = Energy
!    J0                 = approximation of Bessel Function (accounts for gyro-averages)
!
! F and Phi are periodic along the first coordinate :  g(Nksi,...) = g(1,...)
!
!=======================================================================
! SCHEME : MIDPOINT (AG)
!
! * INIT F
!
! * F = krook(F,dt/2)
! * F = dissip(F,dt/2)
! * F = linear_advection(F,dt/2) ! optional
!
! MIDPOINT : advection field at t+dt/2
! * F_half_step = F
! * PhiBAR = compute_potentiel(F)
! * Vad = compute_advection_field(PhiBAR)
! * F = advection(F,PhiBAR,dt/2)
! * PhiBAR = compute_potentiel(F)
! * Vad = compute_advection_field(PhiBAR)
!
! NEWPOINT : F at t+dt
! * F = F_half_step
! * F = advection(F,PhiBAR,dt)
! * F = linear_advection(F,dt/2) ! optional
! * F = dissip(F,dt/2)
! * F = krook(F,dt/2)
!
!=======================================================================

#include "debug.h"

PROGRAM TERESA

  use prec_const
  use globals
  use mem_alloc_module
  use geometry_module
  use advection_module
  use clock_module
  use fftNRF90_module
  use HDF5
  use HDF5_io_module
  ! added from lesur's version, line 86, lili ***
  use eqn_mvt !08/07/2016
  implicit none
  include 'mpif.h'


  !*********************************************
  !	Definition of variables
  !*********************************************

  type (geometry)     :: geom
  type (advec3ts)     :: advec

  integer      :: iter
  integer      :: glob_iter
  real(RKIND)  :: time
  real(RKIND)  :: glob_time
  real(RKIND)  :: max_memory

  real(RKIND)  :: Hnorm, Hnorm_origin, Fmin, Fmin_origin
  real(RKIND)  :: dt_ksi, dt_psi, dt_t, dt_constraint, drift_ksi
  real(RKIND)  :: dt_ksi_tmp, dt_psi_tmp, dt_t_tmp, dt_tmp
  real(RKIND)  :: max_dt_t, max_vad_x1, max_vad_x2, max_vad_x1_tmp, max_vad_x2_tmp
  real(RKIND),dimension(1:3):: dt_com, dt_com_tmp
  real(RKIND)  :: delta_DIAG_f0
  real(RKIND)  :: delta_DIAG_f2_4, delta_DIAG_f2_3, delta_DIAG_f2_1, delta_DIAG_f2_dt
  real(RKIND)  :: delta_DIAG_f4_4, delta_DIAG_f4_3, delta_DIAG_f4_1, delta_DIAG_f4_dt
  real(RKIND)  :: delta_DIAG_f5
  real(RKIND)  :: delta_DIAG

  integer      :: ierr
  integer      :: i, j, k, l, n
!  integer      :: status
  logical      :: fexist

  ! get number of openmp threads
  external omp_get_num_threads
  integer omp_get_num_threads
#ifdef MPI2
  integer :: obtained
#endif

  ! variables used to test conservation during the temporal splitting
  test_density=.false.
  Hnorm = ZE
  Hnorm_origin = ZE
  Fmin = ZE
  Fmin_origin = ZE
!  particles_out_tmp = ZE
!  particles_out = ZE

  !*********************************************
  ! INITIALISATION OF THE MPI WORLD
  !*********************************************

#ifdef MPI2
  call mpi_init_thread(MPI_THREAD_FUNNELED,obtained,ierr)
!  call mpi_init_thread(MPI_THREAD_MULTIPLE,obtained,ierr)
#else
  call MPI_INIT(ierr)
#endif
  call mpi_comm_size(MPI_COMM_WORLD, nproc_glob, ierr) ! nproc_glob = total number of MPI processes
  call mpi_comm_rank(MPI_COMM_WORLD, mype_glob, ierr)  ! tag of each MPI process

  call initialization_time()
  call clck_time(bclock_total,1) ! TOTAL TIME OF EXECUTION


  call clck_time(bclock_init,1) ! TOTAL TIME OF INIT

  !*********************************************
  ! READ PARAMETERS & ALLOCATE MEMORY
  !*********************************************

  !*** write version info ***
  call write_version()

  call avoid_unwanted_restart()

  call memory_requierement()
  call read_data()
  call HDF5_data_saving()
  call create_array()

  !*********************************************
  !   INITIAL CONDITIONS
  !*********************************************

  call init_all(geom,advec)
  max_memory = max_allocate

  call print_memory_size(max_memory)

  call clck_time(eclock_init,0) ! TOTAL TIME OF INIT
  call clck_diff(bclock_init,eclock_init,global_time_init)

  !*********************************************
  !   FIRST STEP
  !*********************************************

  if (mype_glob.eq.0) then
    print*,'************************************'
    print*,'*       ....................       *'
    print*,'    ............................    '
    print*,'  ................................  '
    print*,'...          FIRST STEP          ...'
    print*,'  ................................  '
    print*,'    ............................    '
    print*,'*       ....................       *'
    print*,'************************************'
  end if

  call clck_time(bclock_first_step,1) ! TOTAL OF FIRST STEP

  call MPI_allreduce(RST_T,glob_time,1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,ierr)

  if (RST) then

     glob_iter = RST_iter

     ! in order to determine an accurate value of dt at the first iteration
     !call Compute_potential(.true.)
    !Added fom lesur's version, line 207, lili***
    call Compute_potential(.true.,time)

     call Compute_advection_field(.true.)

  else

     glob_iter = 0
     iter_f0=-1
     iter_f1=-1
     iter_f2=-1
     iter_f3=-1
     iter_f4=-1
     !added from lesur's version, line 218, lili***
     iter_f5=-1
     !Added from Medina, 259,lili
     iter_f6=-1

    ! Normalisation (taking into account the numerical error : Hnorm ~= 1)

     call filters()

     call TestDensity(F,Hnorm,Fmin,0)

     !     do l = lstart,lend
     !        do k = kstart,kend
     !           do j = 1,Npsi
     !              do i = 1,Nksi
     !                 F(i,j,k,l) = F(i,j,k,l)/Hnorm
     !              enddo
     !           enddo
     !        enddo
     !     enddo


     ! computations for diagnotics at t0 and calculate dt


     !call Compute_potential(.true.)
     ! Added from Lesur's version, line 238, lili***
     call Compute_potential(.true.,time)
     call Compute_advection_field(.true.)

     ! monitoring the time of the first saving
     call clck_time(bclock_write,1) ! TOTAL TIME OF WRITE

     call HDF5_init_saving()

     if ( compute_dt ) then
        ! save time derivative of diag
        call save_dt(ZE,glob_time, &
             ZE,ZE,ZE,ZE,&
             ZE,ZE,ZE,ZE)
     end if

     ! save diag
     call save(glob_time)

     call clck_time(eclock_write,0) ! TOTAL TIME OF WRITE
     call clck_diff(bclock_write,eclock_write,global_time_write)

  endif

  call TestDensity(F,Hnorm_origin,Fmin_origin,1) ! ******************************** OFF SET

  call clck_time(eclock_first_step,0) ! TOTAL OF FIRST STEP
  call clck_diff(bclock_first_step,eclock_first_step,global_time_first_step)

  !*********************************************
  !   START TEMPORAL LOOP
  !*********************************************

  if (mype_glob.eq.0) then
     print*,'************************************'
     print*,'*       ....................       *'
     print*,'    ............................    '
     print*,'  ................................  '
     print*,'...     START TEMPORAL LOOP      ...'
     print*,'  ................................  '
     print*,'    ............................    '
     print*,'*       ....................       *'
     print*,'************************************'
  end if

  iter = 0
  time = 0.0_RKIND

  dt_t_tmp = dt_origin
  dt_ksi_tmp = dt_origin
  dt_psi_tmp = dt_origin
  dt_constraint = dt_origin

  dt_t = dt_origin
  dt_ksi = dt_origin
  dt_psi = dt_origin

  TEMPORAL_LOOP : do while (time.lt.(Tmax-dt_eps/2.0_RKIND))

     !**************************************!
     !  Estimation of an optimal time step  !
     !**************************************!

     call clck_time(bclock_compute_dt,1) ! TOTAL TIME OF INIT

     !Added from lesur's version, line 302;303, lili***
     !call EvolvingC1(glob_time)
     !call EvolvingBoundTemperature(glob_time)

        VARIABLE_DT : if (ratio_dt.ne.1) then

        ! compute d Phi / dt
        dPhidt(:,:)=(Phi(:,:)-Phi_old(:,:))/dt

        ! compute criterion on t
        max_dt_t = maxval(abs(dPhidt(:,:))) + 1.0e-12_RKIND
        dt_t_tmp = length_t / max_dt_t

        max_vad_x1=0.0_RKIND
        max_vad_x2=0.0_RKIND


        ! depending on the way to compute the advection
        if (uncouple_linear .or. linear) then

           ! browse advection field
           do l = lstart,lend
              do k = kstart,kend
                 do j = 1,Npsi

                    ! compute criterion on psi, Vad_x2_0 is d Phi / d ksi, motion along ksi
                    max_vad_x2_tmp = maxval(abs(Vad_x2(:,j,k,l)))

                    ! compute criterion on ksi, Vad_x1_0 is d Phi / d psi, motion along psi
                    drift_ksi = OmegaD(j,k)*Energy(l) / Z_species
                    max_vad_x1_tmp = maxval(abs(Vad_x1(:,j,k,l)+drift_ksi))

                    max_vad_x1=max(max_vad_x1,max_vad_x1_tmp) + 1.0e-12_RKIND
                    max_vad_x2=max(max_vad_x2,max_vad_x2_tmp) + 1.0e-12_RKIND

                 end do
              end do
           end do

        else

           ! browse advection field
           do l = lstart,lend
              do k = kstart,kend
                 do j = 1,Npsi

                    ! compute criterion on psi, Vad_x2_0 is d Phi / d ksi, motion along ksi
                    max_vad_x2_tmp = maxval(abs(Vad_x2(:,j,k,l)))

                    ! compute criterion on ksi, Vad_x1_0 is d Phi / d psi, motion along psi
                    max_vad_x1_tmp = maxval(abs(Vad_x1(:,j,k,l)))

                    max_vad_x1=max(max_vad_x1,max_vad_x1_tmp) + 1.0e-12_RKIND
                    max_vad_x2=max(max_vad_x2,max_vad_x2_tmp) + 1.0e-12_RKIND

                 end do
              end do
           end do

        end if


        dt_ksi_tmp = length_ksi * dKsi / max_vad_x1
        dt_psi_tmp = length_psi * dPsi / max_vad_x2

        ! determine the new dt

        ! in order to perform only one communication
        dt_com_tmp(1) = dt_t_tmp
        dt_com_tmp(2) = dt_ksi_tmp
        dt_com_tmp(3) = dt_psi_tmp

        call MPI_allreduce(dt_com_tmp  ,dt_com  ,3,MPI_REAL8,MPI_MIN,MPI_COMM_WORLD,ierr)

        dt_t = dt_com(1)
        dt_ksi = dt_com(2)
        dt_psi = dt_com(3)

        dt_constraint = min(dt_t,dt_ksi,dt_psi)

        n = 1
        dt_tmp = dt_origin

        do while ( (dt_tmp.gt.dt_constraint) .and. (n.lt.ratio_dt) )
           dt_tmp = dt_origin / real(n,RKIND)
           n = n * 2
        end do

        ! logical constraint : criterion to compute the systeme at all steps requiring  diag saving

        delta_DIAG_f0 = dt_DIAG_f0 - mod(glob_time,dt_DIAG_f0)

        delta_DIAG_f5 = dt_DIAG_f5 - mod(glob_time,dt_DIAG_f5)

        delta_DIAG_f4_4 = dt_DIAG_f4 - mod(glob_time+2.0_RKIND*dt_DIAG_dt,dt_DIAG_f4)
        delta_DIAG_f4_3 = dt_DIAG_f4 - mod(glob_time+1.0_RKIND*dt_DIAG_dt,dt_DIAG_f4)
        delta_DIAG_f4_1 = dt_DIAG_f4 - mod(glob_time-1.0_RKIND*dt_DIAG_dt,dt_DIAG_f4)
        delta_DIAG_f4_dt= dt_DIAG_f4 - mod(glob_time-2.0_RKIND*dt_DIAG_dt,dt_DIAG_f4)

        delta_DIAG_f2_4 = dt_DIAG_f2 - mod(glob_time+2.0_RKIND*dt_DIAG_dt,dt_DIAG_f2)
        delta_DIAG_f2_3 = dt_DIAG_f2 - mod(glob_time+1.0_RKIND*dt_DIAG_dt,dt_DIAG_f2)
        delta_DIAG_f2_1 = dt_DIAG_f2 - mod(glob_time-1.0_RKIND*dt_DIAG_dt,dt_DIAG_f2)
        delta_DIAG_f2_dt= dt_DIAG_f2 - mod(glob_time-2.0_RKIND*dt_DIAG_dt,dt_DIAG_f2)

        if(delta_DIAG_f0.le.dt_eps/2.0_RKIND) delta_DIAG_f0 = dt_origin

        if(delta_DIAG_f5.le.dt_eps/2.0_RKIND) delta_DIAG_f5 = dt_origin

        if(delta_DIAG_f4_4.le.dt_eps/2.0_RKIND) delta_DIAG_f4_4 = dt_origin
        if(delta_DIAG_f4_3.le.dt_eps/2.0_RKIND) delta_DIAG_f4_3 = dt_origin
        if(delta_DIAG_f4_1.le.dt_eps/2.0_RKIND) delta_DIAG_f4_1 = dt_origin
        if(delta_DIAG_f4_dt.le.dt_eps/2.0_RKIND) delta_DIAG_f4_dt = dt_origin

        if(delta_DIAG_f2_4.le.dt_eps/2.0_RKIND) delta_DIAG_f2_4 = dt_origin
        if(delta_DIAG_f2_3.le.dt_eps/2.0_RKIND) delta_DIAG_f2_3 = dt_origin
        if(delta_DIAG_f2_1.le.dt_eps/2.0_RKIND) delta_DIAG_f2_1 = dt_origin
        if(delta_DIAG_f2_dt.le.dt_eps/2.0_RKIND) delta_DIAG_f2_dt = dt_origin

        ! test sur le nombre d'iterations, il faut avoir assez avance dans la simulation

        ! delta_DIAG_f0 ! verifier le type,

        delta_DIAG = min(delta_DIAG_f0,&
             delta_DIAG_f5,&
             delta_DIAG_f4_4,delta_DIAG_f4_3,delta_DIAG_f4_1,delta_DIAG_f4_dt,&
             delta_DIAG_f2_4,delta_DIAG_f2_3,delta_DIAG_f2_1,delta_DIAG_f2_dt)

!          if(mype.eq.0)then
!             open(4000+mype_glob,position="append")
!             write(4000+mype_glob,*) "---------"
!             write(4000+mype_glob,*) "delta_DIAG = ",delta_DIAG
!             write(4000+mype_glob,*) "delta_DIAG_f5_1 = ",delta_DIAG_f5_1
!             write(4000+mype_glob,*) "delta_DIAG_f5_dt = ",delta_DIAG_f5_dt
!             write(4000+mype_glob,*) "delta_DIAG_f5_3 = ",delta_DIAG_f5_3
!             write(4000+mype_glob,*) "delta_DIAG_f5_4 = ",delta_DIAG_f5_4
!             write(4000+mype_glob,*) "dt_eps = ",dt_eps
!             write(4000+mype_glob,*) "dt_tmp = ",dt_tmp
!             write(4000+mype_glob,*) "glob_time = ",glob_time
!             write(4000+mype_glob,*) "dt_DIAG_dt = ",dt_DIAG_dt
!             close(4000+mype_glob)
!          end if

        if (delta_DIAG.le.dt_tmp) then
           dt_tmp = delta_DIAG
        end if

!        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!        print*,'mype_glob =,'mype_glob,' dt_tmp = ',dt_tmp

        dt = dt_tmp

     end if VARIABLE_DT


     time = time + dt
     glob_time = glob_time + dt
     iter = iter + 1
     glob_iter = glob_iter + 1
     temps(iter)=time

     if (mype.eq.0) then
        dt_list(iter)=dt
        dt_t_list(iter)=dt_t
        dt_ksi_list(iter)=dt_ksi
        dt_psi_list(iter)=dt_psi
     end if

     if (mype_glob.eq.0) then
        write(6,'(A)') ' -----------------------------------------------------------'
        write(6,'(A24,A8,3A12)') '. iter / [ i_min i_max ]','Gstep','Time','TimeStep','GlobalTime'
        write(6,'(4I8,3(1pe12.3))') iter,iter_min,iter_max,glob_iter,time,dt,glob_time
        if (dt.gt.dt_ksi) then
           write(6,'(A,1pe12.3,A,1pe12.3)') ' WARNING : dt > dt_ksi : ',dt,'>',dt_ksi
        end if
        if (dt.gt.dt_psi) then
           write(6,'(A,1pe12.3,A,1pe12.3)') ' WARNING : dt > dt_psi : ',dt,'>',dt_psi
        end if
        if (dt.gt.dt_t) then
           write(6,'(A,1pe12.3,A,1pe12.3)') ' WARNING : dt > dt_t : ',dt,'>',dt_t
        end if
     end if


     delta_DIAG_f4_4 = dt_DIAG_f4 - mod(glob_time+2.0_RKIND*dt_DIAG_dt,dt_DIAG_f4)
     delta_DIAG_f4_3 = dt_DIAG_f4 - mod(glob_time+1.0_RKIND*dt_DIAG_dt,dt_DIAG_f4)
     delta_DIAG_f4_1 = dt_DIAG_f4 - mod(glob_time-1.0_RKIND*dt_DIAG_dt,dt_DIAG_f4)
     delta_DIAG_f4_dt= dt_DIAG_f4 - mod(glob_time-2.0_RKIND*dt_DIAG_dt,dt_DIAG_f4)

     delta_DIAG_f2_4 = dt_DIAG_f2 - mod(glob_time+2.0_RKIND*dt_DIAG_dt,dt_DIAG_f2)
     delta_DIAG_f2_3 = dt_DIAG_f2 - mod(glob_time+1.0_RKIND*dt_DIAG_dt,dt_DIAG_f2)
     delta_DIAG_f2_1 = dt_DIAG_f2 - mod(glob_time-1.0_RKIND*dt_DIAG_dt,dt_DIAG_f2)
     delta_DIAG_f2_dt= dt_DIAG_f2 - mod(glob_time-2.0_RKIND*dt_DIAG_dt,dt_DIAG_f2)

     delta_DIAG_f4_4 = min(delta_DIAG_f4_4,dt_DIAG_f4-delta_DIAG_f4_4)
     delta_DIAG_f4_3 = min(delta_DIAG_f4_3,dt_DIAG_f4-delta_DIAG_f4_3)
     delta_DIAG_f4_1 = min(delta_DIAG_f4_1,dt_DIAG_f4-delta_DIAG_f4_1)
     delta_DIAG_f4_dt= min(delta_DIAG_f4_dt,dt_DIAG_f4-delta_DIAG_f4_dt)

     delta_DIAG_f2_4 = min(delta_DIAG_f2_4,dt_DIAG_f2-delta_DIAG_f2_4)
     delta_DIAG_f2_3 = min(delta_DIAG_f2_3,dt_DIAG_f2-delta_DIAG_f2_3)
     delta_DIAG_f2_1 = min(delta_DIAG_f2_1,dt_DIAG_f2-delta_DIAG_f2_1)
     delta_DIAG_f2_dt= min(delta_DIAG_f2_dt,dt_DIAG_f2-delta_DIAG_f2_dt)

!        if(mype.eq.0)then
!            print*,delta_DIAG_f4_4
!            print*,delta_DIAG_f4_3
!            print*,delta_DIAG_f4_1
!            print*,delta_DIAG_f4_dt
!           print*,'***',delta_DIAG,dt_tmp,dt_eps
!           print*,delta_DIAG_f2_4
!           print*,delta_DIAG_f2_3
!           print*,delta_DIAG_f2_1
!           print*,delta_DIAG_f2_dt
!           print*,'**********************************************'
!        end if

     call clck_time(eclock_compute_dt,0)
     call clck_diff(bclock_compute_dt,eclock_compute_dt,global_time_compute_dt)

     !*********************************!
     !  Start of the time integration  !
     !*********************************!

     if (nu.ne.0.0) call collision(0.5_RKIND*dt)

     call Dissipation(0.5_RKIND*dt,1)

!     if (uncouple_linear) call L_advance_F(0.5_RKIND*dt)

     call clck_time(bclock_fncopy,1) ! TOTAL TIME OF COPY
     do l = lstart,lend
        do k = kstart,kend
           do j = 1,Npsi
              do i = 1,Nksi
                 F_half_step(i,j,k,l)=F(i,j,k,l)
              end do
           end do
        end do
     end do
     call clck_time(eclock_fncopy,0) ! TOTAL TIME OF COPY
     call clck_diff(bclock_fncopy,eclock_fncopy,global_time_fncopy)

     !**************************!
     !  Start of the mid point  !
     !**************************!

     if (filter_choice.ge.3) call filters()

     !call Compute_potential(.false.)
     ! Added from LEsur's version, line 542, lili ***
     call Compute_potential(.false.,time)
     if(.not. linear) call Compute_advection_field(.false.)
! Added from Medina, 617,620,lili
     if(use_particles) then !save vad_n for rk2 JM
        Vad_x1_deb=Vad_x1
        Vad_x2_deb=Vad_x2
    endif

     call advance_F(advec,geom,0.5_RKIND*dt)
     if (uncouple_linear .or. linear) call L_advance_F(0.5_RKIND*dt)


!     call MPI_allreduce(particles_out_tmp, particles_out,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
!     if (mype.eq.0) then
!        particles_out_0_list(iter)=particles_out
!     end if

     call Source(0.5_RKIND*dt,0)
     if (mype.eq.0) then
        S0_0_list(iter)=S0
     end if

     if (filter_choice.ge.2) call filters()

!     call Compute_potential(.false.)
    !Added from lesur's version, line 559, lili ***
    call Compute_potential(.false.,time-0.5_RKIND*dt)
    if(.not. linear) call Compute_advection_field(.false.)
! Added from Medina, 639,642,lili
 if(use_particles) then !save vad_n+1/2 for rk2 JM
        Vad_x1_demi=Vad_x1
        Vad_x2_demi=Vad_x2
    endif

     if (ratio_dt.ne.1) then
        ! save electric field
        Phi_old(:,:)=Phi(:,:)
     end if


     !************************!
     !  End of the mid point  !
     !************************!

     call clck_time(bclock_fncopy,1) ! TOTAL TIME OF COPY
     do l = lstart,lend
        do k = kstart,kend
           do j = 1,Npsi
              do i = 1,Nksi
                 F(i,j,k,l)=F_half_step(i,j,k,l)
              end do
           end do
        end do
     end do
     call clck_time(eclock_fncopy,0) ! TOTAL TIME OF COPY
     call clck_diff(bclock_fncopy,eclock_fncopy,global_time_fncopy)

     call Source(0.5_RKIND*dt,1)
     if (mype.eq.0) then
        S0_1_list(iter)=S0
     end if

     if (filter_choice.ge.2) call filters()
     if (uncouple_linear .or. linear) call L_advance_F(dt)
     call advance_F(advec,geom,dt)
!     call MPI_allreduce(particles_out_tmp, particles_out,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
!     if (mype.eq.0) then
!        particles_out_1_list(iter)=particles_out
!     end if

     call Source(0.5_RKIND*dt,1)
     if (mype.eq.0) then
        S0_2_list(iter)=S0
     end if

!     if (uncouple_linear) call L_advance_F(0.5_RKIND*dt)

     call Dissipation(0.5_RKIND*dt,-1)

     if (nu.ne.0.0) call collision(0.5_RKIND*dt)

     if (filter_choice.ge.1) call filters()
! Added from Lesur's version, lines 608;622, lili***




     call clck_time(bclock_test_part,1) ! TOTAL TIME OF test part.

!     if(use_particles.and.(mype_glob.eq.0)) then
     if(use_particles) then
!         call Compute_potential(.false.,time)
         !call Compute_advection_field(.false.)
     ! Added from Medina, 697,701, lili,
     if(rungekutta.eq.4) then
             call Compute_advection_field(.false.)
             Vad_x1_end=Vad_x1  !save vad_n+1 for rk4 JM
             Vad_x2_end=Vad_x2
         endif

         !call calcul_P_list(advec,geom,iter) ! JM 11/07/2016 calcule la positiondes particules test a chaque pas de temps !TestParticles
      call calcul_P_list(advec,geom,glob_time) ! added from medina, lili

   endif

     call clck_time(eclock_test_part,0) ! TOTAL TIME OF test part.
     call clck_diff(bclock_test_part,eclock_test_part,global_time_test_part)



     call clck_time(bclock_write,1) ! TOTAL TIME OF WRITE

     !**********!
     !  SAVING  !
     !**********!

     if ( compute_dt ) then
        ! save time derivative of diag
        call save_dt(time,glob_time, &
             delta_DIAG_f4_4,delta_DIAG_f4_3,delta_DIAG_f4_1,delta_DIAG_f4_dt,&
             delta_DIAG_f2_4,delta_DIAG_f2_3,delta_DIAG_f2_1,delta_DIAG_f2_dt)
     end if

     ! save diag
     call save(glob_time)

     call clck_time(eclock_write,0) ! TOTAL TIME OF WRITE
     call clck_diff(bclock_write,eclock_write,global_time_write)

     if (mype_glob.eq.0) then
        inquire(file='TERESA.stop', exist=fexist)
        if (fexist) then
           write(6,*) ' '
           write(6,'(A)') &
                ' ===> The signal file "TERESA.stop" exist '
           write(6,'(A,E12.5,A,E12.5)') &
                ' ===> So! stop the run at glob_time = ', &
                glob_time,' iter = ',iter
        end if
        call MPI_BCAST(fexist,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
     else
        call MPI_BCAST(fexist,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
     end if

     if(fexist) then
        !       status = 1
        exit
     end if

 end do TEMPORAL_LOOP

 !*********************************************
 ! EXIT PROGRAM
 !*********************************************

 call clck_time(bclock_exit,1) ! TOTAL TIME TO END

 call write_trace(iter,glob_iter,time,glob_time,Hnorm_origin,Fmin_origin) ! two warning in hdf5
 call write_restart(advec,geom,glob_time,glob_iter)

 call del_geometry(geom)
 call del_advec3ts(advec)
 call delete_array()

 call clck_time(eclock_exit,0) ! TOTAL TIME TO END
 call clck_diff(bclock_exit,eclock_exit,global_time_exit)

 call clck_time(eclock_total,0)! TOTAL TIME OF EXECUTION
 call clck_diff(bclock_total,eclock_total,global_time_total)

 call write_time() ! one warning by mype in write

 if (mype_glob.eq.0) then
    print*,'************************************'
    print*,'*       ....................       *'
    print*,'    ............................    '
    print*,'  ................................  '
    print*,'...         END OF TERESA        ...'
    print*,'  ................................  '
    print*,'    ............................    '
    print*,'*       ....................       *'
    print*,'************************************'
 end if

 call mpi_finalize(ierr)

END PROGRAM TERESA

! PRINT_DEBUG

! do i = 0,nproc_glob
!    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!    if (mype_glob.eq.i) then
!       PRINT_DEBUG , XXX
!    end if
!    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
! end do
