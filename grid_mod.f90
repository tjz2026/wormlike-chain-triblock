!-----------------------------------------------------------------------
!****m  src/grid_mod
! PURPOSE
!    Declare dimensions ngrid(3) of fft grid
!    Declare global data structures that are defined on an fft grid
!    Routines to allocate, deallocate and manipulate grid data
!    some variables may seem to be redundant and may be removed in future revision
! SOURCE
!-----------------------------------------------------------------------
module grid_mod
  use nrtype,only:SP,DP
   implicit none
   private

   ! Public data structures
   public :: grid_dim          ! dimensions of grid
   public :: grid_dir          ! order of grid, one of the many combinations as XYZ.YZX,ZXY,....
   public :: ngrid          ! dimensions ngrid(:)=(N1,N2,N3) of grid
   public :: local_ngrid    ! local dimensions local_ngrid(:)=(local_nx,local_ny,local_nz) of grid
   public :: ksq_grid       ! k**2 on k-grid 
   public :: global_grid              ! real space grid info
   public :: global_kspace_grid       ! k space grid info
   public :: local_grid_size
   public :: local_kgrid_size
   public :: local_r2c_kgrid_size
   public :: decomp_type_info
   public :: grid_get_info     
   public :: k_vec

   ! Public procedures
   !public :: input_grid
   !public :: generate_grid
   public :: deallocate_grid
   public :: init_local2global_index
   public :: global_grid_decomposition
   public :: k_vec_init


   integer :: grid_dim
   integer :: ngrid(3)
   integer :: local_ngrid(3)
   integer :: decomp_type_info
   integer :: local_grid_size
   integer :: local_kgrid_size
   integer :: local_r2c_kgrid_size

   character(LEN=1) :: grid_dir(3)
   integer :: RSP_start(3),KSP_start(3),RSP_end(3),KSP(3) !global index of the local real and k space grid array
   integer :: RSP_grid_local_num,KSP_grid_local_num
   integer :: RSP_grid_global_num,KSP_grid_global_num
   
type grid_type
      real(DP) :: x
      real(DP) :: y
      real(DP) :: z
end type

   type(grid_type), ALLOCATABLE,save  :: global_grid(:)
   type(grid_type), ALLOCATABLE,save  :: global_kspace_grid(:)
   type(grid_type), ALLOCATABLE,save  :: k_vec(:)
   type(grid_type), ALLOCATABLE,save  :: k_vec_r2c(:)
   real(DP), ALLOCATABLE  :: ksq_grid(:,:,:)

contains
  
  subroutine grid_get_info(dim_grid,dir_grid)
  use global_para,only: SIDEx,SIDEy,SIDEz
  implicit none
  integer :: dim_grid,i
  character(LEN=1) :: dir_grid(3)  

  grid_dim=dim_grid
  grid_dir(:)=dir_grid(:)

  ngrid(:)=1

  if(grid_dim==3) then

  do i=1,3
      if(grid_dir(i)=='X') then
       ngrid(i)=SIDEx
      else if(grid_dir(i)=='Y') then
       ngrid(i)=SIDEy
      else if(grid_dir(i)=='Z') then
       ngrid(i)=SIDEz
      else
      write(*,*) "grid direction unknown,",grid_dir(i)  
      stop
      endif  
  enddo
  else if(grid_dim==2) then

  do i=1,2
      if(grid_dir(i)=='X') then
       ngrid(i)=SIDEx
      else if(grid_dir(i)=='Y') then
       ngrid(i)=SIDEy
      else if(grid_dir(i)=='Z') then
       ngrid(i)=SIDEz
      else
      write(*,*) "grid direction unknown,",grid_dir(i)  
      stop
      endif  
  enddo
  else if(grid_dim==1) then
      if(grid_dir(1)=='X') then
       ngrid(1)=SIDEx
      else if(grid_dir(1)=='Y') then
       ngrid(1)=SIDEy
      else if(grid_dir(1)=='Z') then
       ngrid(1)=SIDEz
      else
      write(*,*) "grid direction unknown,",grid_dir(1)  
      stop
      endif  

 else
 write(*,*) "grid dimension unknown,",grid_dim
 stop
 endif


  call generate_global_grid_zerobased()

  end subroutine grid_get_info


   !==============================================================
   ! generate the grid with offset starts with 0, for example, (0:nx-1,0:ny-1,0:nz-1)
   subroutine generate_global_grid_zerobased()
   implicit none
   integer  :: i,j,k,l,K_i,K_j,K_k,k_kk,k_ii,k_jj
   integer  :: error
    


! generate global grid in real space and k space. 
   if (allocated(global_grid)) deallocate(global_grid)
   allocate(global_grid(0:ngrid(1)*ngrid(2)*ngrid(3)-1))
   if (allocated(global_kspace_grid)) deallocate(global_kspace_grid)
   allocate(global_kspace_grid(0:ngrid(1)*ngrid(2)*ngrid(3)-1))



   k=0
   do l=0,ngrid(3)-1
     do j=0,ngrid(2)-1
        do i=0,ngrid(1)-1
           global_grid(k)%x=i
           global_grid(k)%y=j
           global_grid(k)%z=l
         k=k+1
         enddo
     enddo
   enddo

   k=0
do K_k=0,ngrid(3)-1 
  if(k_k<=(ngrid(3)/2)) then
   k_kk=k_k  
  else 
   k_kk=k_k-ngrid(3)
  endif

   do K_j=0,ngrid(2)-1 
    if(k_j<=(ngrid(2)/2)) then
    k_jj=k_j  
    else 
    k_jj=k_j-ngrid(2)
    endif

     do K_i=0,ngrid(1)-1 
       if(k_i<=(ngrid(1)/2)) then
       k_ii=k_i  
       else 
       k_ii=k_i-ngrid(1)
       endif
       global_kspace_grid(k)%x=K_ii   
       global_kspace_grid(k)%y=K_jj   
       global_kspace_grid(k)%z=K_kk
       k=k+1

     enddo
    enddo
 enddo    

   end subroutine generate_global_grid_zerobased

    ! decomposite the gloabl grid into local grid according to different
    ! FFT methods.( currently, FFTW SLAB, and Pencil decomposition)   
     subroutine global_grid_decomposition(decomp_type,ndir)
     use global_para
# if defined (SLAB)
     use mpi_fftw3_operation
# else /* SLAB */ 
     use decomp_fft_mod
# endif /* SLAB */
     implicit none
     integer, intent(in) :: decomp_type
     integer,optional :: ndir
    ! decomp_type_info ==1 for FFTW slab decomposition and only array orders such as ZYX and YX are supported
    ! decomp_type_info ==2 for pencil decomposition( pencil decomposition works for 3d only) 
    ! for pencil decomposition, array order has to be XYZ!!!
    ! currently the r2c version for fftw slab decomposition hasn't been done yet.

      decomp_type_info=decomp_type
     ! check whether decomp_type_info is legal or not

     if(decomp_type_info==1) then
        if(grid_dim<2) then
        write(*,*) " slab decomposition only works with 2d and 3d model"
        stop
        endif

     local_grid_size=local_size
     local_kgrid_size=local_size
!     local_r2c_grid_size =??
     if(grid_dim==3) then
      local_ngrid(1)=SIDEz
      local_ngrid(2)=SIDEy
      local_ngrid(3)=local_nx_size
     else
      local_ngrid(1)=1
      local_ngrid(2)=SIDEy
      local_ngrid(3)=local_nx_size
     endif

!      call init_k_index()

     else if (decomp_type_info==2) then
        if(grid_dim<3) then
        write(*,*) " pencil decomposition only works with 3d model"
        stop
        endif
!ndir==1, kD is refered to x pencil
!kDz is z pencil decomposition
!ndir==3, kD is refered to z pencil
!kDz is x pencil decomposition
                  

 if(ndir .eqv. 1) then
        local_ngrid(1)=xsize(1)
        local_ngrid(2)=xsize(2)
        local_ngrid(3)=xsize(3)
 else if (ndir .eqv. 2) then
        local_ngrid(1)=ysize(1)
        local_ngrid(2)=ysize(2)
        local_ngrid(3)=ysize(3)
 else if (ndir .eqv. 3) then
        local_ngrid(1)=zsize(1)
        local_ngrid(2)=zsize(2)
        local_ngrid(3)=zsize(3)
 else
 endif


endif  ! endif decompo_type
!# if defined (R2C) 
! call init_local2global_index(ndir,.true.) 
!# else /* R2C */
! call init_local2global_index(ndir) 
!# endif /* R2C */     


     end subroutine global_grid_decomposition



   !--------------------------------------------------------------   
   subroutine deallocate_grid()
   USE global_para
   implicit none
   integer     :: error      ! deallocation error index

   DEALLOCATE(global_grid,STAT=error)
   if (error /= 0) stop "global_grid deallocation error!"

   DEALLOCATE(global_kspace_grid,STAT=error)
   if (error /= 0) stop "k space_grid deallocation error!"

   DEALLOCATE(kD,STAT=error)
   if (error /= 0) stop "kD deallocation error!"

   DEALLOCATE(kDk,STAT=error)
   if (error /= 0) stop "kDk deallocation error!"

   DEALLOCATE(kDz,STAT=error)
   if (error /= 0) stop "kDz deallocation error!"

# if defined (R2C)
   DEALLOCATE(k_vec_r2c,STAT=error)
   if (error /= 0) stop "k_vec_r2c deallocation error!"
# else /* R2C */
   DEALLOCATE(k_vec,STAT=error)
   if (error /= 0) stop "k_vec deallocation error!"
# endif /* R2C */
     
   end subroutine deallocate_grid
   !==============================================================

   
 subroutine init_local2global_index(ndir,R2C_FFT)
 USE nrtype,only :DP,PI
 USE global_para
 USE mpi
 USE control
 USE constants
 USE utility
 USE mmpi
 use decomp_fft_mod
 implicit none
 logical,optional :: R2C_FFT
 integer :: istat,k,k_i,k_j,k_k
 integer,intent(in) :: ndir
 integer :: xs1,xe1,ys1,ye1,zs1,ze1
 integer :: xs2,xe2,ys2,ye2,zs2,ze2

!ndir==1, kD is refered to x pencil
!kDz is z pencil decomposition
!ndir==3, kD is refered to z pencil
!kDz is x pencil decomposition

 allocate(kD(0:LOCAL_SIZE-1),stat=istat)
 allocate(kDz(0:LOCAL_SIZE_K-1),stat=istat)
 allocate(kD_loc(0:LOCAL_SIZE_K-1),stat=istat)
if(ndir .eqv. 1) then
      xs1=xstart(1)
      ys1=xstart(2)
      zs1=xstart(3)
      xe1=xend(1)
      ye1=xend(2)
      ze1=xend(3)

      xs2=zstart(1)
      ys2=zstart(2)
      zs2=zstart(3)
      xe2=zend(1)
      ye2=zend(2)
      ze2=zend(3)
else if (ndir .eqv. 3) then
      xs1=zstart(1)
      ys1=zstart(2)
      zs1=zstart(3)
      xe1=zend(1)
      ye1=zend(2)
      ze1=zend(3)

      xs2=xstart(1)
      ys2=xstart(2)
      zs2=xstart(3)
      xe2=xend(1)
      ye2=xend(2)
      ze2=xend(3)
endif
      
 k=0
 do k_i=zs1,ze1
    do k_j=ys1,ye1
      do k_k=xs1,xe1
        kD(k)%x=k_k-1
        kD(k)%y=k_j-1
        kD(k)%z=k_i-1

        kD_loc(k)%x=k_k-1-(xs1-1)
        kD_loc(k)%y=k_j-1-(ys1-1)
        kD_loc(k)%z=k_i-1-(zs1-1)
        if(kD_loc(k)%x <0 .or. kD_loc(k)%x>xsize(1)-1) then
        write(*,*) "kd_loc",kD_loc(k)%x
        endif            
        if(kD_loc(k)%y <0 .or. kD_loc(k)%y>ysize(2)-1) then
        write(*,*) "kd_loc",kD_loc(k)%y
        endif            
        if(kD_loc(k)%z <0 .or. kD_loc(k)%z>xsize(3)-1) then
        write(*,*) "kd_loc",kD_loc(k)%z
        endif            
        k=k+1
      enddo
    enddo
 enddo
 k=0
 do k_i=zs2,ze2
    do k_j=ys2,ye2
      do k_k=xs2,xe2
        kDz(k)%x=k_k-1
        kDz(k)%y=k_j-1
        kDz(k)%z=k_i-1
        k=k+1
      enddo
    enddo
 enddo

if( present(R2C_FFT) ) then
 allocate(kDk(0:LOCAL_SIZE_K_R2C-1),stat=istat)
  
! note that in future, fft_start(end) will be replaced by the
! variable KSP_start(end) defined in grid mod to avoid conflict with SLAB and Pencil Decomposition

 k=0
 do k_i=fft_start(3),fft_end(3)
    do k_j=fft_start(2),fft_end(2)
      do k_k=fft_start(1),fft_end(1)
        kDk(k)%x=k_k-1
        kDk(k)%y=k_j-1
        kDk(k)%z=k_i-1
        k=k+1
      enddo
    enddo
 enddo

endif
end subroutine init_local2global_index


subroutine k_vec_init(R2C_FFT)
 USE nrtype,only :DP,PI
 USE global_para
 USE mpi
 USE control
 USE constants
 USE utility
 USE mmpi
 use decomp_fft_mod
 implicit none
 logical,optional :: R2C_FFT
 real(DP) :: k_x_factor,k_y_factor,k_z_factor
 integer :: LOCAL_SIZEK,z_s,z_e,y_s,y_e,x_s,x_e
 integer :: k,K_i,K_ii,k_j,k_jj,k_k,k_kk

# if defined (R2C)
LOCAL_SIZEK=LOCAL_SIZE_K_R2C
 allocate(k_vec_r2c(0:LOCAL_SIZEK-1))
# else /* R2C */
LOCAL_SIZEK=LOCAL_SIZE_K
 allocate(k_vec(0:LOCAL_SIZEK-1))
# endif /* R2C */

! when using k_vec, it needs to be multiplied by dx_inv,dy_inv,dz_inv
k_x_factor=2.0d0*PI/(SIDEx*1.0d0)
k_y_factor=2.0d0*PI/(SIDEy*1.0d0)
k_z_factor=2.0d0*PI/(SIDEz*1.0d0)

 if(present(R2C_FFT)) then
  z_s=fft_start(3)
  z_e=fft_end(3)
  y_s=fft_start(2)
  y_e=fft_end(2)
  x_s=fft_start(1)
  x_e=fft_end(1)
 else
  z_s=zstart(3)
  z_e=zend(3)
  y_s=zstart(2)
  y_e=zend(2)
  x_s=zstart(1)
  x_e=zend(1)
endif
   
!k counter
k=0
do K_i=z_s,z_e
  if(k_i-1<=(SiDEz/2)) then
   k_ii=k_i-1
  else 
   k_ii=(k_i-1)-SIDEz 
  endif
    do k_j=y_s,y_e
        if(K_j-1<=(SIDEy/2)) then 
          k_jj=K_j-1
        else
          K_jj=K_j-1-SIDEy
        endif
           do k_k=x_s,x_e
             if(K_k-1<=(SIDEx/2)) then
             k_kk=K_k-1
             else
             k_kk=K_k-1-SIDEx
             endif
# if defined (R2C)
               k_vec_r2c(k)%x=k_kk*k_x_factor
               k_vec_r2c(k)%y=k_jj*k_y_factor
               k_vec_r2c(k)%z=k_ii*k_z_factor
# else /* R2C */
               k_vec(k)%x=k_kk*k_x_factor
               k_vec(k)%y=k_jj*k_y_factor
               k_vec(k)%z=k_ii*k_z_factor
# endif /* R2C */
               k=k+1    
          enddo
         enddo
        enddo

 end subroutine k_vec_init

end module grid_mod
