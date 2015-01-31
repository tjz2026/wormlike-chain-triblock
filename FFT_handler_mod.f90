!!!-----------------------------------------------------------------------
!!! source  : m_FFT_handler.f90
!!! type    : module
!!! author  : Jiuzhou Tang 
!!! purpose : the purpose of this module is to implement 1d,2d,3d FFT
!!!           transform for C2C and (R2C,C2R) using several public avaible 
!!!           efficient parallel fft packages. (FFTW,MKL-FFTW, and 2decomp&FFT etc..)
!!! status  : unstable
!!! comment : only support real(dp) and complex(dp) data types
!!!           FFTW3 is chosen by default,while 2decomp&FFT is compiled with -D2dpencl
!!!-----------------------------------------------------------------------
!! ============
 module fftw3
   use, intrinsic :: iso_c_binding
   include 'fftw3.f03'
 end module fftw3

 module fftw3_mpi
   use, intrinsic :: iso_c_binding
   include 'fftw3-mpi.f03'
 end module fftw3_mpi

 module FFT_handler
  use fftw3_mpi
  implicit none
  private
! fftw3 serial version, 1d fft
!-----------------------------------------------------------------------
!    Wrapper subroutines for fftw3 serial version
! COMMENTS
!    Consider use even/odd fftw, and r-to-r fftw in the future
! SOURCE
!-----------------------------------------------------------------------
   PUBLIC :: s_fft_plan
   PUBLIC :: create_s_fft_plan    ! initialize an fft_plan
   PUBLIC :: s_fftc               ! complex FFT for 1, 2, or 3D
   PUBLIC :: s_fft                ! Forward FFT for 1, 2, or 3D (Real to Complex)
   PUBLIC :: s_ifft               ! Inverse FFT for 1, 2, or 3D (Complex to Real)
! Parameters required by fftw3 supplied in fftw3.f
   integer, parameter :: FFTW_ESTIMATE=64
   integer, parameter :: FFTW_FORWARD=-1 ,FFTW_BACKWARD=1
!-------------------------------------------------------------------
! TYPE
!    s_fft_plan 
! PURPOSE
!    Contains grid dimensions for FFT grid and integer pointers to
!    the "plan" structures used by the serial FFTW package
! SOURCE
!-------------------------------------------------------------------
   type s_fft_plan
      integer    ::  n(3)   ! grid dimensions, 0 for unused dimensions
      integer*8  ::  f      ! fftw plan object for forward transform
      integer*8  ::  r      ! fftw plan object for inverse transform
   end type s_fft_plan

!  FFT wrapers for MPI FFTW version ( to be implemented )
!--------------------------------------------------------------------------------------
   !PUBLIC :: mpi_fft_plan                                                            - 
   !PUBLIC :: create_mpi_fft_plan    ! initialize an fft_plan                         -
   !PUBLIC :: mpi_fftc  ! complex FFT for 1, 2, or 3D                                 - 
   !PUBLIC :: mpi_fft   ! Forward FFT for 1, 2, or 3D (Real to Complex)               -
   !PUBLIC :: mpi_ifft  ! Inverse FFT for 1, 2, or 3D (Complex to Real)               -
!--------------------------------------------------------------------------------------
! Warning: we assume that (N1,N2,N3) correspongs to the dimension grid of X,Y,Z axis by default.
! This implies that for 3d case,we will do slab decomposition on the Z axis and on Y axis for 2d case, 
! which is different from old version code.
   integer(C_INTPTR_T), public :: N1 ! the first dimension in fortran array, varies fastest
   integer(C_INTPTR_T), public :: N2 
   integer(C_INTPTR_T), public :: N3
   type(C_PTR),public :: fplan2d,bplan2d, cdata
   type(C_PTR),public :: fplan3d,bplan3d, cdata_r
   type(C_PTR),public :: fplan_r2c_2d,fplan_c2r_2d
   type(C_PTR),public :: fplan_r2c_3d,fplan_c2r_3d

! Issue : if multigrid scheme is adopted in future, we must use a set of grids, thus 
! FFT plan should be defined as variable size arrays. like :
! type(C_PTR),allocatable, public :: fplan3d(:)

! defining C2C transform variable and work array
   complex(C_DOUBLE_COMPLEX), pointer,public :: local_data2d(:,:)
   complex(C_DOUBLE_COMPLEX), pointer,public :: local_data3d(:,:,:)
   integer(C_INTPTR_T),public ::  alloc_local, local_N2, local_N2_offset
   integer(C_INTPTR_T),public ::  local_N3_offset,local_N3
   integer,public :: local_nz,local_z_start
   integer,public :: local_ny,local_y_start
! defining R2C,C2R transform variable and work array
   real(C_DOUBLE), pointer,public  :: local_data2d_r(:,:)
   real(C_DOUBLE), pointer,public  :: local_data3d_r(:,:,:)

!---------------------------------------------------------------------------------
! FFTW public routines for 2d,3d,c2c and r2c/c2r transforms
! Issue: better overloading these routines into one public interface in future.
   public :: fftw3_mpi_2d_dft_init
   public :: fftw3_mpi_3d_dft_init
   public :: fftw3_mpi_2d_r2c_init
   public :: fftw3_mpi_3d_r2c_init
   public :: destroy_mpi_fftw3_3d
   public :: destroy_mpi_fftw3_2d
   public :: destroy_mpi_fftw3_r2c_3d
   public :: destroy_mpi_fftw3_r2c_2d
   
   public :: fftw_2d_c2c
   public :: fftw_3d_c2c
   public :: fftw_2d_r2c
   public :: fftw_3d_r2c
   public :: fftw_2d_c2r
   public :: fftw_3d_c2r
!---------------------------------------------------------------------------------
!Parallel FFT with  Pencil decomposition ( current engine : 2decomp&FFT)

contains

   !-------------------------------------------------------------------
   !****p fft_mod/create_fft_plan
   ! SUBROUTINE
   !    create_fft_plan
   ! PURPOSE
   !    Creates an fft_plan object for grids with dimensions 
   !    ngrid(1),..,ngrid(dim)
   !-------------------------------------------------------------------
   subroutine create_s_fft_plan(ngrid,plan,fft_c2c)
   integer,intent(IN)             :: ngrid(3) ! dimensions of grid
   type(s_fft_plan),intent(OUT)     :: plan
   logical, optional, intent(IN)  :: fft_c2c
   !***
   plan%n=ngrid
   end subroutine create_s_fft_plan
   !===================================================================
   !-------------------------------------------------------------------
   !****p fft_mod/fft
   ! SUBROUTINE
   !     fft(plan,in,out)
   ! PURPOSE
   !    Calculates forward fft of in, returns result in out.
   !    Wrapper for 1, 2, & 3 dimensional real -> complex transforms
   ! ARGUMENTS
   !    plan    - fft plan object
   !    in, out - real(long) 3D arrays 
   ! COMMENT
   !    in and out are dimensioned 0:ngrid(i)-1 for all i <= dim, 
   !    and 0:0 for any unused dimensions with dim < i <= 3
   !-------------------------------------------------------------------
   subroutine s_fft(plan,in,out)
   integer, parameter :: dp    = kind(1.0d0)
   type(s_fft_plan),intent(IN)   :: plan
   real(dp), intent(IN)      :: in(0:,0:,0:)
   complex(dp), intent(OUT)  :: out(0:,0:,0:)
   !***
   call dfftw_plan_dft_r2c_3d(plan%f,plan%n(1),plan%n(2),plan%n(3),&!
                              in,out,FFTW_ESTIMATE)
   call dfftw_execute(plan%f)
   call dfftw_destroy_plan(plan%f)
   end subroutine s_fft
   !===================================================================


   !-------------------------------------------------------------------
   !****p fft_mod/ifft
   ! SUBROUTINE
   !     ifft(plan,in,out)
   ! PURPOSE
   !    Calculates inverse fft of real array in, returns in out.
   !    Wrapper for 1, 2, & 3 dimensional complex -> real transforms
   ! ARGUMENTS
   !    plan - fft plan object
   !    in   - complex(long) 3D input array
   !    out  - real(long) 3D input array
   ! COMMENT
   !    in and out are dimensioned 0:ngrid(i)-1 for all i <= dim, 
   !    and 0:0 for any unused dimensions with dim < i <= 3
   !-------------------------------------------------------------------
   subroutine s_ifft(plan,in,out)
   integer, parameter :: dp    = kind(1.0d0)
   type(s_fft_plan),intent(IN)   :: plan
   complex(dp), intent(IN)   :: in(0:,0:,0:)
   real(dp), intent(OUT)     :: out(0:,0:,0:)
   !***
   call dfftw_plan_dft_c2r_3d(plan%r,plan%n(1),plan%n(2),plan%n(3),&!
                              in,out,FFTW_ESTIMATE)
   call dfftw_execute(plan%r)
   call dfftw_destroy_plan(plan%r)
   end subroutine s_ifft
   !===================================================================


   !-------------------------------------------------------------------
   !****p fft_mod/fftc
   ! SUBROUTINE
   !     fftc(direction,plan,in,out)
   ! PURPOSE
   !    Calculates forward fft of in, returns result in out.
   !    Wrapper for 1, 2, & 3 dimensional real -> complex transforms
   ! ARGUMENTS
   !    plan    - fft plan object
   !    in, out - real(long) 3D arrays 
   ! COMMENT
   !    in and out are dimensioned 0:ngrid(i)-1 for all i <= dim, 
   !    and 0:0 for any unused dimensions with dim < i <= 3
   !-------------------------------------------------------------------
   subroutine s_fftc(direction,plan,in,out)
   integer, parameter :: dp    = kind(1.0d0)
   integer,intent(IN)          :: direction
   type(s_fft_plan),intent(IN)   :: plan
   complex(dp), intent(IN)   :: in(0:,0:,0:)
   complex(dp), intent(OUT)  :: out(0:,0:,0:)
   !***
   if (direction == 1) then
      call dfftw_plan_dft_3d(plan%f,plan%n(1),plan%n(2),plan%n(3),&!
             in,out,FFTW_FORWARD,FFTW_ESTIMATE)
   else
      call dfftw_plan_dft_3d(plan%f,plan%n(1),plan%n(2),plan%n(3),&!
             in,out,FFTW_BACKWARD,FFTW_ESTIMATE)
   end if
   call dfftw_execute(plan%f)
   call dfftw_destroy_plan(plan%f)

   if (direction == +1) out = out/dcmplx( dble(plan%n(1)*plan%n(2)*plan%n(3)) , 0.0d0)
   end subroutine s_fftc


!Slab decomposition on Y axis with Ny grid number.
! Warning : tips for fortran interface of fftw3:
!Note that when we called fftw_mpi_local_size_2d and fftw_mpi_plan_dft_2d with the dimensions in reversed order, 
!since a L × M Fortran array is viewed by FFTW in C as a M × L array. This means that the array was distributed
! over the M dimension, the local portion of which is a L × local_M array in Fortran. (You must not use an allocate
! statement to allocate an L × local_M array, however; you must allocate alloc_local complex numbers, which may be 
!greater than L * local_M, in order to reserve space for intermediate steps of the transform.) Finally, we mention 
!that because C's array indices are zero-based, the local_j_offset argument can conveniently be interpreted as an 
!offset in the 1-based j index (rather than as a starting index as in C).
  subroutine fftw3_mpi_2d_dft_init(Nx,Ny)
   use mpi
   use fftw3_mpi
   use mpi_control, ONLY: myid, nprocs
   implicit none
   integer,intent(in) :: Nx,Ny

   N1=Nx
   N2=Ny

    call fftw_mpi_init()
!   get local data size and allocate (note dimension reversal)
  alloc_local = fftw_mpi_local_size_2d( N2,N1, MPI_COMM_WORLD, &
                                       local_N2, local_N2_offset)

  cdata = fftw_alloc_complex(alloc_local)

  call c_f_pointer(cdata, local_data2d, [N1,local_N2])

!   create MPI plan for in-place forward DFT (note dimension reversal)
  fplan2d = fftw_mpi_plan_dft_2d(N2,N1, local_data2d, local_data2d, MPI_COMM_WORLD, &
                              FFTW_FORWARD, FFTW_EXHAUSTIVE)

  bplan2d = fftw_mpi_plan_dft_2d(N2,N1, local_data2d, local_data2d, MPI_COMM_WORLD, &
                              FFTW_BACKWARD, FFTW_EXHAUSTIVE)

  write(*,*) "Local_N2=,local_N2_offset=",Local_N2,Local_N2_offset,"on",myid 

  local_ny=Local_N2
  Local_y_start=Local_N2_offset

! compute transform (as many times as desired)
!  call fftw_mpi_execute_dft(fplan, data, data)
!  call fftw_mpi_execute_dft(bplan, data, data)
! write(*,*) "size of local_datain fftw3",size(local_data),"on",myid

end subroutine fftw3_mpi_2d_dft_init

  subroutine fftw3_mpi_3d_dft_init(Nx,Ny,Nz)
    use mpi
    use fftw3_mpi
    USE mpi_control, ONLY: myid, nprocs
    implicit none
    integer,intent(in) :: Nx,Ny,Nz
    N1=Nx  
    N2=Ny 
    N3=Nz   

    call fftw_mpi_init()
!   get local data size and allocate (note dimension reversal)
  alloc_local = fftw_mpi_local_size_3d(N3, N2,N1, MPI_COMM_WORLD, &
                                       local_N3, local_N3_offset)
  write(*,*) "alloc_local",alloc_local,"on",myid

  cdata = fftw_alloc_complex(alloc_local)

  call c_f_pointer(cdata, local_data3d, [N1,N2,local_N3])

!   create MPI plan for in-place forward DFT (note dimension reversal)
!  fplan = fftw_mpi_plan_dft_3d(N3,N2,N1, local_data, local_data, MPI_COMM_WORLD, &
!                              FFTW_FORWARD, FFTW_EXHAUSTIVE)
  fplan3d = fftw_mpi_plan_dft_3d(N3,N2,N1, local_data3d, local_data3d, MPI_COMM_WORLD, &
                              FFTW_FORWARD, FFTW_MEASURE)

  bplan3d = fftw_mpi_plan_dft_3d(N3,N2,N1, local_data3d, local_data3d, MPI_COMM_WORLD, &
                              FFTW_BACKWARD, FFTW_MEASURE)

  write(*,*) "Local_N3=,local_N3_offset=",Local_N3,Local_N3_offset,"on",myid 

  local_nz=Local_N3
  Local_z_start=Local_N3_offset

! compute transform (as many times as desired)
!  call fftw_mpi_execute_dft(fplan, data, data)
!  call fftw_mpi_execute_dft(bplan, data, data)
!  data=data/(32*32)
! write(*,*) "data(1,1)af fftw3",data(1,1),"on",myid
! write(*,*) "size of local_datain fftw3",size(local_data),"on",myid
end subroutine fftw3_mpi_3d_dft_init


   !Real to Complex 2d slab decomposition on Y axis with Ny grid number.
  subroutine fftw3_mpi_2d_r2c_init(Nx,Ny)
   use mpi
   use fftw3_mpi
   use mpi_control, ONLY: myid, nprocs
   implicit none
   integer,intent(in) :: Nx,Ny

   N1=Nx
   N2=Ny
    call fftw_mpi_init()
!   get local data size and allocate (note dimension reversal)
  alloc_local = fftw_mpi_local_size_2d( N2,N1/2+1, MPI_COMM_WORLD, &
                                       local_N2, local_N2_offset)

  cdata_r = fftw_alloc_real(2*alloc_local)
  cdata   = fftw_alloc_complex(alloc_local)

  call c_f_pointer(cdata_r, local_data2d_r, [N1+2,local_N2])
  call c_f_pointer(cdata,   local_data2d, [N1/2+1,local_N2])

!   create MPI plan for in-place forward R2C (note dimension reversal)
  fplan_r2c_2d = fftw_mpi_plan_r2c_2d(N2,N1, local_data2d_r, local_data2d, MPI_COMM_WORLD, &
                              FFTW_FORWARD, FFTW_EXHAUSTIVE)

  bplan_c2r_2d = fftw_mpi_plan_c2r_2d(N2,N1, local_data2d, local_data2d_r, MPI_COMM_WORLD, &
                              FFTW_BACKWARD, FFTW_EXHAUSTIVE)

  write(*,*) "Local_N2=,local_N2_offset=",Local_N2,Local_N2_offset,"on",myid 

  local_ny=Local_N2
  Local_y_start=Local_N2_offset

! compute transform (as many times as desired)
!  call fftw_mpi_execute_dft(fplan, data, data)
!  call fftw_mpi_execute_dft(bplan, data, data)
! write(*,*) "size of local_datain fftw3",size(local_data),"on",myid

  end subroutine fftw3_mpi_2d_r2c_init

!out-of-place r2c transform of L × M × N real data [padded to L × M × 2(N/2+1)],
! resulting in L × M × N/2+1 complex data.
  subroutine fftw3_mpi_3d_r2c_init(Nx,Ny,Nz)
    use mpi
    use fftw3_mpi
    USE mpi_control, ONLY: myid, nprocs
    implicit none
    integer,intent(in) :: Nx,Ny,Nz
    N1=Nx  
    N2=Ny 
    N3=Nz   

    call fftw_mpi_init()
!   get local data size and allocate (note dimension reversal)
  alloc_local = fftw_mpi_local_size_3d(N3, N2,N1/2+1, MPI_COMM_WORLD, &
                                       local_N3, local_N3_offset)
  write(*,*) "alloc_local",alloc_local,"on",myid

  cdata_r = fftw_alloc_real(2*alloc_local)
  cdata = fftw_alloc_complex(alloc_local)

  call c_f_pointer(cdata_r, local_data3d_r, [N1+2,local_N2])
  call c_f_pointer(cdata, local_data3d, [N1/2+1,N2,local_N3])

  fplan_r2c_3d = fftw_mpi_plan_dft_3d(N3,N2,N1, local_data3d_r, local_data3d, MPI_COMM_WORLD, &
                              FFTW_FORWARD, FFTW_MEASURE)

  bplan_c2r_3d = fftw_mpi_plan_dft_3d(N3,N2,N1, local_data3d, local_data3d_r, MPI_COMM_WORLD, &
                              FFTW_BACKWARD, FFTW_MEASURE)

  write(*,*) "Local_N3=,local_N3_offset=",Local_N3,Local_N3_offset,"on",myid 

  local_nz=Local_N3
  Local_z_start=Local_N3_offset

! compute transform (as many times as desired)
!  call fftw_mpi_execute_dft(fplan, data, data)
!  call fftw_mpi_execute_dft(bplan, data, data)
!  data=data/(32*32)
! write(*,*) "data(1,1)af fftw3",data(1,1),"on",myid
! write(*,*) "size of local_datain fftw3",size(local_data),"on",myid
end subroutine fftw3_mpi_3d_r2c_init


! routines to performe 2d and 3d C2C/R2C/C2R fft transform 

 subroutine fftw_2d_c2c(dir)
   implicit none
   integer,intent(in) :: dir
   select case (dir)
   case (1)
   call fftw_mpi_execute_dft(fplan_2d, local_data2d, local_data2d)
   case (-1)
   call fftw_mpi_execute_dft(bplan_2d, local_data2d, local_data2d)
   case default
   write(*,*) "error in fft c2c direction, dir must be 1 or -1"
   stop
   end select
 end subroutine fftw_2d_c2c

 subroutine fftw_3d_c2c(dir)
   implicit none
   integer,intent(in) :: dir
   select case (dir)
   case (1)
   call fftw_mpi_execute_dft(fplan_3d, local_data3d, local_data3d)
   case (-1)
   call fftw_mpi_execute_dft(bplan_3d, local_data3d, local_data3d)
   case default
   write(*,*) "error in fft c2c direction, dir must be 1 or -1"
   stop
   end select
 end subroutine fftw_3d_c2c

 subroutine fftw_2d_r2c()
   implicit none
   call fftw_mpi_execute_r2c(fplan_r2c_2d, local_data2d_r, local_data2d)
 end subroutine fftw_2d_r2c

 subroutine fftw_2d_c2r()
   implicit none
   call fftw_mpi_execute_c2r(fplan_c2r_2d, local_data2d, local_data2d_r)
 end subroutine fftw_2d_c2r

 subroutine fftw_3d_r2c()
   implicit none
  call fftw_mpi_execute_r2c(fplan_r2c_3d, local_data3d_r, local_data3d)
 end subroutine fftw_3d_r2c

 subroutine fftw_3d_c2r()
   implicit none
  call fftw_mpi_execute_c2r(fplan_c2r_3d, local_data3d, local_data3d_r)
 end subroutine fftw_3d_c2r


subroutine destroy_mpi_fftw3_2d()
  use fftw3_mpi
  implicit none
  call fftw_destroy_plan(fplan2d)
  call fftw_destroy_plan(bplan2d)
  call fftw_free(cdata)
end subroutine destroy_mpi_fftw3_2d

subroutine destroy_mpi_fftw3_3d()
  use fftw3_mpi
  implicit none
  call fftw_destroy_plan(fplan3d)
  call fftw_destroy_plan(bplan3d)
  call fftw_free(cdata)
end subroutine destroy_mpi_fftw3_3d

subroutine destroy_mpi_fftw3_r2c_2d()
  use fftw3_mpi
  implicit none
  call fftw_destroy_plan(fplan_r2c_2d)
  call fftw_destroy_plan(bplan_c2r_2d)
  call fftw_free(cdata)
  call fftw_free(cdata_r)
end subroutine destroy_mpi_fftw3_r2c_2d

subroutine destroy_mpi_fftw3_r2c_3d()
  use fftw3_mpi
  implicit none
  call fftw_destroy_plan(fplan_r2c_3d)
  call fftw_destroy_plan(bplan_c2r_3d)
  call fftw_free(cdata)
  call fftw_free(cdata_r)
end subroutine destroy_mpi_fftw3_r2c_3d


end module FFT_handler

# if defined (PENCL)
!2decomp&FFT package (Parallel FFT transform using 2d pencil decomposition)
!Warning: to use decomp_fft_mod, the grid must be 3 dimension.(Can't do pencil for 2d data) 
  module decomp_fft_mod
  use decomp_2d
  use decomp_2d_fft
  implicit none
  !define public variables 
  integer, public :: Nx ! grid size on X axis
  integer, public :: Ny 
  integer, public :: Nz
! two dimensional grid (p_row, p_col ) of cpu processors
  integer, public :: p_row 
  integer, public :: p_col
! defining the start and end index of local data,more details can be found on
! http://www.2decomp.org/
  integer, dimension(3),public :: fft_start, fft_end, fft_size
  complex(mytype), allocatable, dimension(:,:,:) :: local_in, local_out
  real(mytype), allocatable, dimension(:,:,:)    :: local_in_r
  complex(mytype), allocatable, dimension(:,:,:) :: local_out_c
 contains

 
 subroutine decomp2d_mpi_3d_dft_c2c_init(N1,N2,N3,n_prow,n_pcol)
  use mpi
  USE mpi_control, ONLY: myid, nprocs
  implicit none
  integer,intent(in) :: N1,N2,N3
  integer,intent(in) :: n_prow,n_pcol
  Nx=N1
  Ny=N2
  Nz=N3

  nrank=myid
  nproc=nprocs

  p_row=n_prow
  p_col=n_pcol

  call decomp_2d_init(nx,ny,nz,p_row,p_col)
  call decomp_2d_fft_init()
  
  allocate (local_in(xstart(1)-1:xend(1)-1,xstart(2)-1:xend(2)-1,xstart(3)-1:xend(3)-1))
  allocate (local_out(zstart(1)-1:zend(1)-1,zstart(2)-1:zend(2)-1,zstart(3)-1:zend(3)-1))
  
  end subroutine decomp2d_mpi_3d_dft_c2c_init

 subroutine decomp2d_mpi_3d_dft_r2c_init(loc_hf_size)
  use mpi
  USE mpi_control, ONLY: myid, nprocs
  implicit none
  integer,intent(in) :: N1,N2,N3
  integer,intent(in) :: n_prow,n_pcol
  integer,intent(inout) :: loc_hf_size

  Nx=N1
  Ny=N2
  Nz=N3

  nrank=myid
  nproc=nprocs

  p_row=n_prow
  p_col=n_pcol

  call decomp_2d_init(nx,ny,nz,p_row,p_col)
  call decomp_2d_fft_init

  allocate (local_in_r(xstart(1)-1:xend(1)-1,xstart(2)-1:xend(2)-1,xstart(3)-1:xend(3)-1))
  allocate (local_out(zstart(1)-1:zend(1)-1,zstart(2)-1:zend(2)-1,zstart(3)-1:zend(3)-1))

  call decomp_2d_fft_get_size(fft_start,fft_end,fft_size)

  allocate (local_out_c(fft_start(1)-1:fft_end(1)-1, &
       fft_start(2)-1:fft_end(2)-1, &
       fft_start(3)-1:fft_end(3)-1))

   loc_hf_size=fft_size(1)*fft_size(2)*fft_size(3)
  
  if(myid==0) then
  write(*,*) "decomp_2d FFT for r2c transform initialized"
  write(*,*) "xsize,zsize,fftsize=",xsize(1)*xsize(2)*xsize(3),&
             zsize(1)*zsize(2)*zsize(3), loc_hf_size
  endif

end subroutine decomp2d_mpi_3d_dft_r2c_init

subroutine destroy_mpi_decomp2d()
  implicit none
! deallocate the fft input and output array here if necessary, to be done in future!

call decomp_2d_fft_finalize
call decomp_2d_finalize
end subroutine destroy_mpi_decomp2d

subroutine assemble_global(ndir,local,global,nx,ny,nz)
  
  use mpi
  USE mpi_control, ONLY: myid, nprocs
  implicit none
  
  integer, intent(IN) :: ndir  ! 1 = X-pencil; 3 = Z-pencil
  integer, intent(IN) :: nx,ny,nz
  !complex(mytype), dimension(:,:,:), intent(IN) :: local
  !complex(mytype), dimension(nx,ny,nz), intent(OUT) :: global
  
  !complex(mytype), allocatable, dimension(:,:,:) :: rbuf
  real(mytype), dimension(:,:,:), intent(IN) :: local
  real(mytype), dimension(nx,ny,nz), intent(OUT) :: global
  
  real(mytype), allocatable, dimension(:,:,:) :: rbuf
  integer, dimension(9) :: sbuf1, rbuf1
  
  integer :: ierror, i,j,k,m, i1,i2,j1,j2,k1,k2, count
  integer, dimension(MPI_STATUS_SIZE) :: status
  
  nrank=myid
  nproc=nprocs

  if (nrank==0) then
     ! master writes its own data to a global array
     if (ndir==3) then  ! Z-pencil 
        i1 = zstart(1)
        i2 = zend(1)
        j1 = zstart(2)
        j2 = zend(2)
        k1 = zstart(3)
        k2 = zend(3)
     else if (ndir==1) then  ! X-pencil
        i1 = xstart(1)
        i2 = xend(1)
        j1 = xstart(2)
        j2 = xend(2)
        k1 = xstart(3)
        k2 = xend(3)
     end if
     do k=k1,k2
        do j=j1,j2
           do i=i1,i2
              ! 'local' is assumbed shape array
              ! but it is OK as starting index for rank 0 always 1
              global(i,j,k)=local(i,j,k)
           end do
        end do
     end do
     ! then loop through all other ranks to collect data
     do m=1,nproc-1
        CALL MPI_RECV(rbuf1,9,MPI_INTEGER,m,m,MPI_COMM_WORLD, &
             status,ierror)
        allocate(rbuf(rbuf1(1):rbuf1(2),rbuf1(4):rbuf1(5), &
             rbuf1(7):rbuf1(8)))
        CALL MPI_RECV(rbuf,rbuf1(3)*rbuf1(6)*rbuf1(9),real_type,m, &
             m+nproc,MPI_COMM_WORLD,status,ierror)
        do k=rbuf1(7),rbuf1(8)
           do j=rbuf1(4),rbuf1(5)
              do i=rbuf1(1),rbuf1(2)
                 global(i,j,k)=rbuf(i,j,k)
              end do
           end do
        end do
        deallocate(rbuf)
     end do
  else
     ! slaves send data to mater
     if (ndir==3) then  ! Z-pencil
        sbuf1(1) = zstart(1)
        sbuf1(2) = zend(1)
        sbuf1(3) = zsize(1)
        sbuf1(4) = zstart(2)
        sbuf1(5) = zend(2)
        sbuf1(6) = zsize(2)
        sbuf1(7) = zstart(3)
        sbuf1(8) = zend(3)
        sbuf1(9) = zsize(3)
        count = zsize(1)*zsize(2)*zsize(3)
     else if (ndir==1) then  ! X-pencil
        sbuf1(1) = xstart(1)
        sbuf1(2) = xend(1)
        sbuf1(3) = xsize(1)
        sbuf1(4) = xstart(2)
        sbuf1(5) = xend(2)
        sbuf1(6) = xsize(2)
        sbuf1(7) = xstart(3)
        sbuf1(8) = xend(3)
        sbuf1(9) = xsize(3)
        count = xsize(1)*xsize(2)*xsize(3)
     end if
     ! send partition information
     CALL MPI_SEND(sbuf1,9,MPI_INTEGER,0,nrank,MPI_COMM_WORLD,ierror)
     ! send data array
     CALL MPI_SEND(local,count,real_type,0, &
          nrank+nproc,MPI_COMM_WORLD,ierror)
  end if
  
  return
end subroutine assemble_global


subroutine assembly_local_1d_to_global_3d(ndir,local_1d,global_3d,nx,ny,nz)
  implicit none
  integer :: k,k_k,k_j,k_i
  integer, intent(IN) :: ndir  ! 1 = X-pencil; 3 = Z-pencil
  integer, intent(IN) :: nx,ny,nz
  real(mytype), dimension(:), intent(IN) :: local_1d
  real(mytype), allocatable, dimension(:,:,:) :: local_3d
  real(mytype), dimension(nx,ny,nz), intent(OUT) :: global_3d

if(ndir==1) then
allocate (local_3d(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)))
else if(ndir==3) then
allocate (local_3d(zstart(1):zend(1),zstart(2):zend(2),zstart(3):zend(3)))
else
write(*,*) "error input of integer variable ndir, must be either 1 or 3, exit"
stop
endif

local_3d=0.0
global_3d=0.0

   if(ndir==1) then
          k=0
        do K_i=xstart(3),xend(3)
          do k_j=xstart(2),xend(2)
            do K_k=xstart(1),xend(1)
             local_3d(K_k,k_j,k_i)=local_1d(k+1)
             k=k+1
            enddo
          enddo
        enddo
    else if(ndir==3) then
          k=0
        do K_i=zstart(3),zend(3)
          do k_j=zstart(2),zend(2)
            do K_k=zstart(1),zend(1)
             local_3d(K_k,k_j,k_i)=local_1d(k+1)
             k=k+1
            enddo
          enddo
        enddo
    endif

  call assemble_global(ndir,local_3d,global_3d,nx,ny,nz) 

deallocate(local_3d)

end subroutine assembly_local_1d_to_global_3d

end module decomp_fft_mod
# endif /* PENCL */
