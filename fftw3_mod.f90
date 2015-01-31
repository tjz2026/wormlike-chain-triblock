
   module fftw3
   use, intrinsic :: iso_c_binding
   include 'fftw3.f03'
   end module fftw3

  module single_fftw3_operation
  use fftw3
  ! note that currently the single core fftw3 doesn't support SIMD optimization,soon this
  ! issue will be dealing with. and this is fftw module is for legacy fortran ,not modern fortran
  !, which could be very different from each other, be careful!
  
  type(C_PTR) :: plan_1d_f,plan_1d_b,plan_2d_f,plan_2d_b
  type(C_PTR) :: plan_3d_f,plan_3d_b
  complex(C_DOUBLE_COMPLEX), dimension(:),allocatable :: in_1d, out_1d
  complex(C_DOUBLE_COMPLEX), dimension(:,:),allocatable :: in_2d, out_2d
  complex(C_DOUBLE_COMPLEX), dimension(:,:,:),allocatable :: in_3d, out_3d

   public :: fftw3_c2c_init
     interface fftw3_c2c_init
         module procedure fftw3_1d_c2c_init
         module procedure fftw3_2d_c2c_init
         module procedure fftw3_3d_c2c_init
     end interface fftw3_c2c_init

     public :: fftw3_c2c_clean
     interface fftw3_c2c_clean
         module procedure fftw3_1d_c2c_clean
         module procedure fftw3_2d_c2c_clean
         module procedure fftw3_3d_c2c_clean
     end interface fftw3_c2c_clean



  contains

  subroutine fftw3_1d_c2c_init(N1)
  implicit none 
  integer :: N1
  allocate(in_1d(1:N1))
  allocate(out_1d(1:N1))

 call dfftw_plan_dft_1d(plan_1d_f ,N1,in_1d,out_1d,FFTW_FORWARD,FFTW_ESTIMATE)
 call dfftw_plan_dft_1d(plan_1d_b ,N1,in_1d,out_1d,FFTW_BACKWARD,FFTW_ESTIMATE)

  end subroutine fftw3_1d_c2c_init

  subroutine fftw3_2d_c2c_init(N1,N2)
  implicit none 
  integer :: N1,N2
  allocate(in_2d(1:N1,1:N2))
  allocate(out_2d(1:N1,1:N2))
 call dfftw_plan_dft_2d(plan_2d_f ,N1,N2,in_2d,out_2d,FFTW_FORWARD,FFTW_ESTIMATE)
 call dfftw_plan_dft_2d(plan_2d_b ,N1,N2,in_2d,out_2d,FFTW_BACKWARD,FFTW_ESTIMATE)

  end subroutine fftw3_2d_c2c_init

  subroutine fftw3_3d_c2c_init(N1,N2,N3)
  implicit none 
  integer :: N1,N2,N3
  allocate(in_3d(1:N1,1:N2,1:N3))
  allocate(out_3d(1:N1,1:N2,1:N3))
 call dfftw_plan_dft_2d(plan_3d_f ,N1,N2,in_3d,out_3d,FFTW_FORWARD,FFTW_ESTIMATE)
 call dfftw_plan_dft_2d(plan_3d_b ,N1,N2,in_3d,out_3d,FFTW_BACKWARD,FFTW_ESTIMATE)

  end subroutine fftw3_3d_c2c_init


  subroutine fftw3_1d_c2c_clean(N1)
  integer :: N1
  if (allocated(in_1d)) deallocate(in_1d)
  if (allocated(out_1d)) deallocate(out_1d)
  call dfftw_destroy_plan( plan_1d_f )
  call dfftw_destroy_plan( plan_1d_b )
  end subroutine fftw3_1d_c2c_clean

  subroutine fftw3_2d_c2c_clean(N1,N2)
  integer :: N1,N2
  if (allocated(in_2d)) deallocate(in_2d)
  if (allocated(out_2d)) deallocate(out_2d)
  call dfftw_destroy_plan( plan_2d_f )
  call dfftw_destroy_plan( plan_2d_b )
  end subroutine fftw3_2d_c2c_clean

  subroutine fftw3_3d_c2c_clean(N1,N2,N3)
  integer :: N1,N2,N3
  if (allocated(in_3d)) deallocate(in_3d)
  if (allocated(out_3d)) deallocate(out_3d)
  call dfftw_destroy_plan( plan_3d_f )
  call dfftw_destroy_plan( plan_3d_b )
  end subroutine fftw3_3d_c2c_clean

  end module single_fftw3_operation



# if defined (SLAB)
# if defined (Mem_OP)   
  module fftw3_mpi
   use, intrinsic :: iso_c_binding
   include 'fftw3-mpi.f03'
   end module fftw3_mpi

  module mpi_fftw3_operation
  use fftw3_mpi

  integer(C_INTPTR_T), public :: N1
  integer(C_INTPTR_T), public :: N2
  integer(C_INTPTR_T), public :: N3
  type(C_PTR),public :: fplan,bplan, cdata
  complex(C_DOUBLE_COMPLEX), pointer,public :: local_data(:)
  integer(C_INTPTR_T),public ::  alloc_local, local_N3, local_N3_offset
  integer :: local_nx,local_x_start
  type(C_PTR),public :: fplan_r2c_2d,fplan_c2r_2d
  type(C_PTR),public :: fplan_r2c_3d,fplan_c2r_3d
  contains

  subroutine fftw3_mpi_3d_dft_init()
  use mpi
  use fftw3_mpi
  use global_para,only :SIDEx,SIDEy,SIDEz
  USE control, ONLY: myid, nprocs
  implicit none
  N1=SIDEz   
  N2=SIDEy   
  N3=SIDEx   


    call fftw_mpi_init()

!   get local data size and allocate (note dimension reversal)
  alloc_local = fftw_mpi_local_size_3d(N3, N2,N1, MPI_COMM_WORLD, &
                                       local_N3, local_N3_offset)
  cdata = fftw_alloc_complex(alloc_local)

  call c_f_pointer(cdata, local_data, [N1,N2,local_N3])

!   create MPI plan for in-place forward DFT (note dimension reversal)
  fplan = fftw_mpi_plan_dft_3d(N3,N2,N1, local_data, local_data, MPI_COMM_WORLD, &
                              FFTW_FORWARD, FFTW_EXHAUSTIVE)

  bplan = fftw_mpi_plan_dft_3d(N3,N2,N1, local_data, local_data, MPI_COMM_WORLD, &
                              FFTW_BACKWARD, FFTW_EXHAUSTIVE)

  write(*,*) "Local_N3=,local_N3_offset=",Local_N3,Local_N3_offset,"on",myid 

  local_nx=Local_N3
  Local_x_start=Local_N3_offset

! compute transform (as many times as desired)
!  call fftw_mpi_execute_dft(fplan, data, data)
!  call fftw_mpi_execute_dft(bplan, data, data)
!  data=data/(32*32)
! write(*,*) "data(1,1)af fftw3",data(1,1),"on",myid

end subroutine fftw3_mpi_3d_dft_init

subroutine destroy_3d_mpi_fftw3()
  use fftw3_mpi
  implicit none
  call fftw_destroy_plan(fplan)
  call fftw_destroy_plan(bplan)
  call fftw_free(cdata)
end subroutine destroy_3d_mpi_fftw3

end module mpi_fftw3_operation

# else  /* Mem_OP */   
   module fftw3_mpi
   use, intrinsic :: iso_c_binding
   include 'fftw3-mpi.f03'
   end module fftw3_mpi

  module mpi_fftw3_operation
  use fftw3_mpi
   
  integer(C_INTPTR_T), public :: N1
  integer(C_INTPTR_T), public :: N2
  integer(C_INTPTR_T), public :: N3
  type(C_PTR),public :: fplan,bplan, cdata
# if defined (Dim2)
  complex(C_DOUBLE_COMPLEX), pointer,public :: local_data(:,:)
# else  /* Dim2 */
  complex(C_DOUBLE_COMPLEX), pointer,public :: local_data(:,:,:)
# endif /* Dim2 */
  integer(C_INTPTR_T),public ::  alloc_local, local_N2, local_N2_offset
  integer(C_INTPTR_T),public ::  local_N3_offset,local_N3
  integer,public :: local_nx,local_x_start
  contains
# if defined (Dim2)
  subroutine fftw3_mpi_2d_dft_init()
  use mpi
  use fftw3_mpi
  use global_para,only :SIDEx,SIDEy,SIDEz
  USE control, ONLY: myid, nprocs
  implicit none
  N1=SIDEy   
  N2=SIDEx   


    call fftw_mpi_init()

!   get local data size and allocate (note dimension reversal)
  alloc_local = fftw_mpi_local_size_2d( N2,N1, MPI_COMM_WORLD, &
                                       local_N2, local_N2_offset)

  cdata = fftw_alloc_complex(alloc_local)

  call c_f_pointer(cdata, local_data, [N1,local_N2])

!   create MPI plan for in-place forward DFT (note dimension reversal)
  fplan = fftw_mpi_plan_dft_2d(N2,N1, local_data, local_data, MPI_COMM_WORLD, &
                              FFTW_FORWARD, FFTW_EXHAUSTIVE)

  bplan = fftw_mpi_plan_dft_2d(N2,N1, local_data, local_data, MPI_COMM_WORLD, &
                              FFTW_BACKWARD, FFTW_EXHAUSTIVE)

  write(*,*) "Local_N2=,local_N2_offset=",Local_N2,Local_N2_offset,"on",myid 

  local_nx=Local_N2
  Local_x_start=Local_N2_offset

! compute transform (as many times as desired)
!  call fftw_mpi_execute_dft(fplan, data, data)
!  call fftw_mpi_execute_dft(bplan, data, data)
!  data=data/(32*32)
! write(*,*) "data(1,1)af fftw3",data(1,1),"on",myid
 write(*,*) "size of local_datain fftw3",size(local_data),"on",myid

end subroutine fftw3_mpi_2d_dft_init
# else  /* Dim2 */
  subroutine fftw3_mpi_3d_dft_init()
  use mpi
  use fftw3_mpi
  use global_para,only :SIDEx,SIDEy,SIDEz
  USE control, ONLY: myid, nprocs
  implicit none
  N1=SIDEz   
  N2=SIDEy   
  N3=SIDEx   


    call fftw_mpi_init()

!   get local data size and allocate (note dimension reversal)
  alloc_local = fftw_mpi_local_size_3d(N3, N2,N1, MPI_COMM_WORLD, &
                                       local_N3, local_N3_offset)
  write(*,*) "alloc_local",alloc_local,"on",myid

  cdata = fftw_alloc_complex(alloc_local)

  call c_f_pointer(cdata, local_data, [N1,N2,local_N3])

!   create MPI plan for in-place forward DFT (note dimension reversal)
!  fplan = fftw_mpi_plan_dft_3d(N3,N2,N1, local_data, local_data, MPI_COMM_WORLD, &
!                              FFTW_FORWARD, FFTW_EXHAUSTIVE)
  fplan = fftw_mpi_plan_dft_3d(N3,N2,N1, local_data, local_data, MPI_COMM_WORLD, &
                              FFTW_FORWARD, FFTW_MEASURE)

  bplan = fftw_mpi_plan_dft_3d(N3,N2,N1, local_data, local_data, MPI_COMM_WORLD, &
                              FFTW_BACKWARD, FFTW_MEASURE)

  write(*,*) "Local_N3=,local_N3_offset=",Local_N3,Local_N3_offset,"on",myid 
  write(*,*) "N1=,N2=",N1,N2,"on",myid 
  write(*,*) "local_data(N1,N2,local_N3_offset)",local_data(32,32,2),"on",myid 

  local_nx=Local_N3
  Local_x_start=Local_N3_offset

! compute transform (as many times as desired)
!  call fftw_mpi_execute_dft(fplan, data, data)
!  call fftw_mpi_execute_dft(bplan, data, data)
!  data=data/(32*32)
! write(*,*) "data(1,1)af fftw3",data(1,1),"on",myid

 write(*,*) "size of local_datain fftw3",size(local_data),"on",myid
end subroutine fftw3_mpi_3d_dft_init

# endif /* Dim2 */

subroutine destroy_mpi_fftw3()
  use fftw3_mpi
  implicit none
  call fftw_destroy_plan(fplan)
  call fftw_destroy_plan(bplan)
  call fftw_free(cdata)
end subroutine destroy_mpi_fftw3

end module mpi_fftw3_operation
# endif /* Mem_OP */
# endif /* SLAB */



