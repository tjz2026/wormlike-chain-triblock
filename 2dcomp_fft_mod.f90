# if defined (SLAB)
# else /* SLAB */
  module decomp_fft_mod
  use decomp_2d
  use decomp_2d_fft

  integer, public :: NX
  integer, public :: NY
  integer, public :: NZ
  integer, public :: p_row
  integer, public :: p_col
  integer, dimension(3),public :: fft_start, fft_end, fft_size
  complex(mytype), allocatable, dimension(:,:,:) :: local_in, local_out
  real(mytype), allocatable, dimension(:,:,:) :: local_in_r
  complex(mytype), allocatable, dimension(:,:,:) :: local_out_c
 contains

 subroutine decomp2d_mpi_3d_dft_c2c_init()
  use mpi
  use global_para,only :SIDEx,SIDEy,SIDEz
  USE control, ONLY: myid, nprocs
  implicit none
  NX=SIDEx
  NY=SIDEy
  NZ=SIDEz

nrank=myid
nproc=nprocs
!p_row=4
!p_col=4
call decomp_2d_init(nx,ny,nz,p_row,p_col)
call decomp_2d_fft_init

allocate (local_in(xstart(1)-1:xend(1)-1,xstart(2)-1:xend(2)-1,xstart(3)-1:xend(3)-1))
allocate (local_out(zstart(1)-1:zend(1)-1,zstart(2)-1:zend(2)-1,zstart(3)-1:zend(3)-1))

end subroutine decomp2d_mpi_3d_dft_c2c_init

 subroutine decomp2d_mpi_3d_dft_r2c_init()
  use mpi
  use global_para,only :SIDEx,SIDEy,SIDEz,local_size_k_r2c
  USE control, ONLY: myid, nprocs
  implicit none
  NX=SIDEx
  NY=SIDEy
  NZ=SIDEz

nrank=myid
nproc=nprocs
!p_row=4
!p_col=4
call decomp_2d_init(nx,ny,nz,p_row,p_col)
call decomp_2d_fft_init

  allocate (local_in_r(xstart(1)-1:xend(1)-1,xstart(2)-1:xend(2)-1,xstart(3)-1:xend(3)-1))
  allocate (local_out(zstart(1)-1:zend(1)-1,zstart(2)-1:zend(2)-1,zstart(3)-1:zend(3)-1))

  call decomp_2d_fft_get_size(fft_start,fft_end,fft_size)

  allocate (local_out_c(fft_start(1)-1:fft_end(1)-1, &
       fft_start(2)-1:fft_end(2)-1, &
       fft_start(3)-1:fft_end(3)-1))

   local_size_k_r2c=fft_size(1)*fft_size(2)*fft_size(3)
  
  if(myid==0) then
  write(*,*) "decomp_2d FFT for r2c transform initialized"
  write(*,*) "xsize,zsize,fftsize=",xsize(1)*xsize(2)*xsize(3),&
             zsize(1)*zsize(2)*zsize(3), local_size_k_r2c
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
  USE control, ONLY: myid, nprocs
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
# endif /* SLAB */
