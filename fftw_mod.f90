
MODULE fftw

!****************************************************************
!INCLUDES:
! SUBROUTINE init_fftw-mpi
! SUBROUTINE stop_fftw_mpi
!****************************************************************

!****************************************************************
! SIZE and mapping info for parallel FFTW
 INTEGER,public :: local_nz, local_z_start
 INTEGER,public :: local_nx, local_x_start
 INTEGER,public :: local_ny_after_transpose, local_nz_after_transpose
 INTEGER,public :: local_y_start_after_transpose, local_z_start_after_transpose
 INTEGER,public :: total_local_size
! INTEGER,public :: local_size

!***************Declarations of FFTW plans**********************
	Integer*8 planf     		! forward plan for use of fft => xyz
	Integer*8 planb     		! backward plan for use of ifft => xyz
	Integer*8 plancf     		! forward plan for use of fft => xyz
	Integer*8 plancb     		! backward plan for use of ifft => xyz
!****************************************************************
!     This file contains PARAMETER statements for various constants
!     that can be passed to FFTW routines.  You should include
!     this file in any FORTRAN program that calls the fftw_f77
!     routines (either directly or with an #include statement
!     if you use the C preprocessor).

      integer FFTW_FORWARD,FFTW_BACKWARD
      parameter (FFTW_FORWARD=-1,FFTW_BACKWARD=1)

      integer FFTW_REAL_TO_COMPLEX,FFTW_COMPLEX_TO_REAL
      parameter (FFTW_REAL_TO_COMPLEX=-1,FFTW_COMPLEX_TO_REAL=1)

      integer FFTW_ESTIMATE,FFTW_MEASURE
      parameter (FFTW_ESTIMATE=0,FFTW_MEASURE=1)

      integer FFTW_OUT_OF_PLACE,FFTW_IN_PLACE,FFTW_USE_WISDOM
      parameter (FFTW_OUT_OF_PLACE=0)
      parameter (FFTW_IN_PLACE=8,FFTW_USE_WISDOM=16)

      integer FFTW_THREADSAFE
      parameter (FFTW_THREADSAFE=128)

!     Constants for the MPI wrappers:
      integer,public :: FFTW_TRANSPOSED_ORDER, FFTW_NORMAL_ORDER
      integer,public :: FFTW_SCRAMBLED_INPUT, FFTW_SCRAMBLED_OUTPUT
      parameter (FFTW_TRANSPOSED_ORDER=1, FFTW_NORMAL_ORDER=0)
      parameter(FFTW_SCRAMBLED_INPUT=8192)
      parameter(FFTW_SCRAMBLED_OUTPUT=16384)

CONTAINS

 SUBROUTINE init_fftw_mpi
!***************************************************************************
! Needs:
!   nx,ny,nz --> global system sizes in module global (parameters)
!   Code must be compiled with the fftw libs
!   MPI must be already started (ie call start_mpi)
!***************************************************************************
  USE global_para, ONLY: nx=>SIDEx, ny=>SIDEy,nz=>SIDEz
  USE control, ONLY: myid, nprocs
  USE mpi
  IMPLICIT NONE
  INTEGER :: i, j,k,ierror, problm, problm_t
!DOUBLE PRECISION, DIMENSION(:),ALLOCATABLE :: data,work
 !double complex  work,data1
! double precision ::  work,data1
!dimension data1(0:nx-1,0:ny-1,0:nz-1), work(0:nx-1,0:ny-1,0:nz-1)
! dimension data1(nx,ny,nz), work(nx,ny,nz)
!dimension data1(32,32,32), work(32,32,32)
!nx=32
!ny=32
!!nz=32
 !
 !
!  !     write(6,*)
!
!	CALL rfftw3d_f77_mpi_create_plan (planf, mpi_comm_world, nx, ny, nz,&	! Plan for ->FFT of an array of size (nx,ny,nz)
!    				            & FFTW_FORWARD, FFTW_ESTIMATE)
!	CALL rfftw3d_f77_mpi_create_plan (planb, mpi_comm_world, nx, ny, nz,&	! Plan for <-FFT of an array of size (nx,ny,nz)
!    				            & FFTW_BACKWARD,FFTW_ESTIMATE)
!	CALL rfftwnd_f77_mpi_local_sizes (planf, local_nz,local_z_start,&	! Get local partition information
!                                    & local_ny_after_transpose, local_y_start_after_transpose, total_local_size)		      

	CALL fftw3d_f77_mpi_create_plan (plancf, mpi_comm_world, nz, ny, nx,&	! Plan for ->FFT of an array of size (nz,ny,nx)
    				            & FFTW_FORWARD, FFTW_ESTIMATE)
	CALL fftw3d_f77_mpi_create_plan (plancb, mpi_comm_world, nz, ny, nx,&	! Plan for <-FFT of an array of size (nx,ny,nz)
    				            & FFTW_BACKWARD,FFTW_ESTIMATE)
	CALL fftwnd_f77_mpi_local_sizes (plancf, local_nx,local_x_start,&	! Get local partition information
                                    & local_ny_after_transpose, local_y_start_after_transpose, total_local_size)		      











  !******** Check FFTW partitioning*********
    problm = 0
    if ((local_nx == 0) .or. (local_ny_after_transpose == 0) .or. &
        (total_local_size == 0) )  then
       if (myid==0)  write(6,*) "Zero size found in FFTW partitioning !"
       problm = 1
    else if ((local_nx /= nx/nprocs) .or. (local_ny_after_transpose /= ny/nprocs)) then   
       if (myid==0)  write(6,*) "Uneven FFTW partitioning !"
       problm = 1
!    else if ((total_local_size /= 2*(nx/2+1)*ny*local_nz))  then
!       if (myid==0)  write(6,*) "Mismatched total_local_size !"
!       problm = 1
    end if


        print *,"address of fwd plan =   ",plancf
        print *,"address of bckwd plan = ",plancb
    call mpi_barrier(mpi_comm_world,ierror)
    call mpi_allreduce (problm, problm_t, 1, mpi_integer, mpi_sum, mpi_comm_world, ierror)
    
    if (problm_t > 0)  then
       do i = 0,nprocs-1
          if (i == myid) then
            write(6,*) "Node = ", myid
            write(6,*) "local_nx = ",local_nx
            write(6,*) "local_x_start = ",local_x_start
            write(6,*) "local_ny_after_transpose = ",local_ny_after_transpose
            write(6,*) "total_local_size = ",total_local_size
            write(6,*)
            write(6,*)
          endif
          call mpi_barrier(mpi_comm_world,ierror)
       end do
       call mpi_finalize(ierror)
       call exit(0)
      else  
       do i = 0,nprocs-1
          if (i == myid) then
            write(6,*) "Node = ", myid
            write(6,*) "local_nx = ",local_nx
            write(6,*) "local_x_start = ",local_x_start
            write(6,*) "ny,nz = ",ny,nz
            write(6,*) "local_ny_after_transpose = ",local_ny_after_transpose
            write(6,*) "total_local_size = ",total_local_size
            write(6,*)
            write(6,*)
          endif
        enddo
    endif


! ALLOCATE (data(1:total_local_size),work(1:total_local_size))
!data=1.0
!work=0.3
!!do i=1,32
!!do j=1,32
!!do k=1,32
!! data1(i,j,k) = 1.0
!! work(i,j,k) = 0.5
!!!work(i,j,k) = data1(i,j,k)
!!enddo
!!enddo
!!enddo
!!
!
!call mpi_barrier(mpi_comm_world,ierror)
!CALL rfftwnd_f77_mpi (planf, 1, data, work, 1, FFTW_TRANSPOSED_ORDER)
!CALL rfftwnd_f77_mpi (planb, 1, data, work, 1, FFTW_TRANSPOSED_ORDER)
!call mpi_barrier(mpi_comm_world,ierror)
!data=data/(nx*ny*nz)
!write(*,*) "data(2)",data(2),"on",myid
!!call fftwnd_f77_mpi(planf, 1, data1, work, 0, fftw_normal_order)
!call mpi_barrier(mpi_comm_world,ierror)
!write(*,*) "done FFT in fft_mod ","on",myid
!    DEALLOCATE(data)
!    DEALLOCATE(work)

!stop
END SUBROUTINE init_fftw_mpi

!subroutine rfftw3D_mpi_f(data,work,nx,ny,nz)
!use mpi
!implicit none
!!DOUBLE PRECISION, DIMENSION(:),ALLOCATABLE :: data,work
!DOUBLE PRECISION, intent(inout) :: data(total_local_size),work(total_local_size)
!integer n, nx, ny, nz, nit,ii
!integer i,j,k
!
! !ALLOCATE (data(1:total_local_size),work(1:total_local_size))
!
!
!
!call mpi_barrier(mpi_comm_world,ierror)
!CALL rfftwnd_f77_mpi (planf, 1, data, work, 1, FFTW_TRANSPOSED_ORDER)
!CALL rfftwnd_f77_mpi (planb, 1, data, work, 1, FFTW_TRANSPOSED_ORDER)
!call mpi_barrier(mpi_comm_world,ierror)
!data=data/(nx*ny*nz)
!!write(*,*) "data(2)",data(2),"on",myid
!!call fftwnd_f77_mpi(planf, 1, data1, work, 0, fftw_normal_order)
!call mpi_barrier(mpi_comm_world,ierror)
!write(*,*) "done rFFT in fft_mod ","on",myid
! !   DEALLOCATE(data)
!  !  DEALLOCATE(work)
!
!end subroutine rfftw3D_mpi_f

!************complex fft***********************

!subroutine fftw3D_mpi_f(data,work,nx,ny,nz)
!use mpi
!implicit none
!!DOUBLE PRECISION, DIMENSION(:),ALLOCATABLE :: data,work
!DOUBLE complex, intent(inout) :: data(total_local_size),work(total_local_size)
!integer n, nx, ny, nz, nit,ii
!integer i,j,k
!
! !ALLOCATE (data(1:total_local_size),work(1:total_local_size))
!
!
!
!call mpi_barrier(mpi_comm_world,ierror)
!CALL rfftwnd_f77_mpi (plancf, 1, data, work, 1, FFTW_TRANSPOSED_ORDER)
!CALL rfftwnd_f77_mpi (plancb, 1, data, work, 1, FFTW_TRANSPOSED_ORDER)
!call mpi_barrier(mpi_comm_world,ierror)
!data=data/(nx*ny*nz)
!write(*,*) "data(2).re",real(data(2)),"on",myid
!!call fftwnd_f77_mpi(planf, 1, data1, work, 0, fftw_normal_order)
!call mpi_barrier(mpi_comm_world,ierror)
!write(*,*) "done rFFT in fft_mod ","on",myid
! !   DEALLOCATE(data)
!  !  DEALLOCATE(work)
!
!end subroutine fftw3D_mpi_f
!



 subroutine stop_fftw_mpi
!***************************************************************************
! Destroy FFTW plans
!***************************************************************************
    USE mpi

    IMPLICIT NONE
 !   CALL rfftwnd_f77_mpi_destroy_plan(planf)
 !   CALL rfftwnd_f77_mpi_destroy_plan(planb)
    CALL fftwnd_f77_mpi_destroy_plan(plancf)
    CALL fftwnd_f77_mpi_destroy_plan(plancb)

  !  if (node == 0) then
  !     write(6,*)
  !     write(6,*) "FFTW plans destroyed"
  !     write(6,*)
  !  endif

  END  SUBROUTINE stop_fftw_mpi

END  MODULE fftw
