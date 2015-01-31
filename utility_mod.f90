   module utility
   implicit none
   integer,private :: isize
   integer,public :: seed  !needs to be initialize


   public ::  random_init_matrix
   private :: random_init_matrix_int
   private :: random_init_matrix_dp
   private :: random_init_matrix_cmplx


   public ::  random_init_array3D
   private :: random_init_array3D_int
   private :: random_init_array3D_dp
   private :: random_init_array3D_cmplx

 !  public ::  array2D_to_1D
 !  private :: array2D_to_1D_dp
 !  private :: array2D_to_1D_cmplx

   interface random_init_matrix
           module procedure random_init_matrix_int
           module procedure random_init_matrix_dp
           module procedure random_init_matrix_cmplx
   end interface random_init_matrix

   interface random_init_array3D
           module procedure random_init_array3D_int
           module procedure random_init_array3D_dp
           module procedure random_init_array3D_cmplx
   end interface random_init_array3D
   
   contains

       subroutine random_init_matrix_int(a,nx,ny)
       implicit none
       integer i,j,nx,ny
       integer,intent(inout) :: a(nx,ny)
       
       isize=size(a)
       if(isize/=nx*ny) then
       write(*,*)" matrix random init error,size mismatched!"
       stop
       endif

       do i=1,nx
         do j=1,ny
       a(i,j)=nint(ran2(seed)*10.0)
          enddo
       enddo
       end subroutine random_init_matrix_int
           
       subroutine random_init_matrix_dp(a,nx,ny)
       implicit none
       integer i,j,nx,ny
       real*8,intent(inout) :: a(nx,ny)
       
       isize=size(a)
       if(isize/=nx*ny) then
       write(*,*)" matrix random init error,size mismatched!"
       stop
       endif

       do i=1,nx
         do j=1,ny
       a(i,j)=ran2(seed)*1.0d0
          enddo
       enddo
       end subroutine random_init_matrix_dp
       subroutine random_init_matrix_cmplx(a,nx,ny)
       implicit none
       integer i,j,nx,ny
       double complex,intent(inout) :: a(nx,ny)
       
       isize=size(a)
       if(isize/=nx*ny) then
       write(*,*)" matrix random init error,size mismatched!"
       stop
       endif

       do i=1,nx
         do j=1,ny
        a(i,j)=cmplx(ran2(seed)*1.0,0.0)
          enddo
       enddo
       end subroutine random_init_matrix_cmplx
           
       subroutine random_init_array3D_int(a,nx,ny,nz)
       implicit none
       integer i,j,k,nx,ny,nz
       integer,intent(inout) :: a(nx,ny,nz)
       
       isize=size(a)
       if(isize/=nx*ny*nz) then
       write(*,*)" 3D array random init error,size mismatched!"
       stop
       endif

       do i=1,nx
         do j=1,ny
           do k=1,nz
       a(i,j,k)=nint(ran2(seed)*10.0)
          enddo
       enddo
      enddo
       end subroutine random_init_array3D_int

       subroutine random_init_array3D_dp(a,nx,ny,nz)
       implicit none
       integer i,j,k,nx,ny,nz
       real*8,intent(inout) :: a(nx,ny,nz)
       
       isize=size(a)
       if(isize/=nx*ny*nz) then
       write(*,*)" 3D array random init error,size mismatched!"
       stop
       endif

       do i=1,nx
         do j=1,ny
           do k=1,nz
       a(i,j,k)=ran2(seed)*1.0
          enddo
       enddo
      enddo
       end subroutine random_init_array3D_dp

       subroutine random_init_array3D_cmplx(a,nx,ny,nz)
       implicit none
       integer i,j,k,nx,ny,nz
       double complex,intent(inout) :: a(nx,ny,nz)
       
       isize=size(a)
       write(*,*) "isize=",isize,"nx*ny*nz=",nx*ny*nz
       if(isize/=nx*ny*nz) then
       write(*,*)" 3D array random init error,size mismatched!"
       stop
       endif

       do i=1,nx
         do j=1,ny
           do k=1,nz
       a(i,j,k)=cmplx(ran2(seed)*1.0,0.0)
          enddo
       enddo
      enddo
       end subroutine random_init_array3D_cmplx



           
           
       FUNCTION ran2(idum)
      INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      REAL ran2,AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1, &
        IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791,	 &
        NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum2
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/
      if (idum.le.0) then
        idum=max(-idum,1)
        idum2=idum
        do j=NTAB+8,1,-1
          k=idum/IQ1
          idum=IA1*(idum-k*IQ1)-k*IR1
          if (idum.lt.0) idum=idum+IM1
          if (j.le.NTAB) iv(j)=idum
        end do
        iy=iv(1)
      endif
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum.lt.0) idum=idum+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      if (idum2.lt.0) idum2=idum2+IM2
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1)iy=iy+IMM1
      ran2=AM*iy
      if(ran2.gt.RNMX) ran2=RNMX
      return
      END function  



!!!!Simpson integration function for PBC boundary and MPI SLAB decomposition
! it assumes the periodic boundary condition here.
! calculate the integration through transforming 3D matrix to 1D matrix

    FUNCTION simposon_3D_1D_mpi(nz,ny,nx,dz,dy,dx,a31)
    implicit none
    integer :: nz,ny,nx,i,j,k
    REAL*8  :: dx,dy,dz
    REAL*8  :: simposon_3D_1D_mpi
    REAL*8  :: a31(0:nz*ny*nx-1)
    REAL*8,allocatable :: a3(:,:,:)
    
    allocate(a3(0:nz-1,0:ny-1,0:nx-1))
    
     do k=0,nx-1
        do j=0,ny-1
          do i=0,nz-1
            a3(i,j,k)=a31(k*ny*nz+j*nz+i)
          enddo
        enddo
     enddo
     simposon_3D_1D_mpi=simposon_3D_mpi(nz,ny,nx,dz,dy,dx,a3)
    deallocate(a3)
    return
    end function simposon_3D_1D_mpi


    FUNCTION simposon_2D_1D_mpi(ny,nx,dy,dx,a21)
    implicit none
    integer :: ny,nx,i,j,k
    REAL*8  :: dx,dy
    REAL*8  :: simposon_2D_1D_mpi
    REAL*8  :: a21(0:ny*nx-1)
    REAL*8,allocatable :: a2(:,:)
    
    allocate(a2(0:ny-1,0:nx-1))
    
     do k=0,nx-1
        do j=0,ny-1
            a2(j,k)=a21(k*ny+j)
        enddo
     enddo
     simposon_2D_1D_mpi=simposon_2D_mpi(ny,nx,dy,dx,a2)
    deallocate(a2)
    return
    end function simposon_2D_1D_mpi

   FUNCTION simposon_3D_mpi(nz,ny,nx,dz,dy,dx,a3)
   implicit none
   integer :: s1,flag,i,j
   REAL*8 :: simposon_3D_mpi
   REAL*8 :: FUN_A,FUN_B,FUN2,FUN4
   INTEGER :: nz,ny,nx
   REAL*8 :: dx,dy,dz
   REAL*8 :: a3(0:nz-1,0:ny-1,0:nx-1)
   REAL*8,allocatable :: a1(:),a2(:,:)
   
!! n1, n2, n3 should be even numbers i.e. 2*N (N is an arbitrary integer larger than zero);
!!      no volume average.  
!!    double a2[n1][n2][n3] --- 0 ~ n1-1 ,  0 ~ n2-1, 0 ~ n3-1  ;  it assumes the periodic boundary condition here.


   allocate(a1(0:nx-1))
   allocate(a2(0:ny-1,0:nx-1))
   do i=0,nx-1
       do j=0,ny-1
        FUN_A=a3(0,j,i)
        FUN_B=a3(0,j,i)  !PBC boundary condition !!!!
        FUN2=0.0
        FUN4=a3(1,j,i)
          s1=1
          !note that if nz==2,the 3D SCFT turns into 2D ,and the simposon integral remains
         ! correct.(nz/2-1=0,do inner loop !)
          do flag=1,nz/2-1,1
           s1=s1+1
           FUN2=FUN2+a3(s1,j,i)
           s1=s1+1
           FUN4=FUN4+a3(s1,j,i)
           enddo
          a2(j,i)=(dz/3.0)*(FUN_A+FUN_B+4.0*FUN4+2.0*FUN2)
        enddo
    enddo
     
     do i=0,nx-1
         FUN_A=a2(0,i)
         FUN_B=a2(0,i)
         FUN2=0.0
         FUN4=a2(1,i)
         s1=1
         do flag=1,ny/2-1,1
             s1=s1+1
             FUN2=FUN2+a2(s1,i)
             s1=s1+1
             FUN4=FUN4+a2(s1,i)
          enddo
          a1(i)=(dy/3.0)*(FUN_A+FUN_B+4.0*FUN4+2.0*FUN2)
     enddo
         FUN2=0.0
         FUN4=0.0
        do s1=0,nx-1,2
          FUN2=FUN2+a1(s1)
          FUN4=FUN4+a1(s1+1)
         enddo
     simposon_3D_mpi=(dx/3.0)*(4.0*FUN4+2.0*FUN2)
     deallocate(a2)
     deallocate(a1)
     return
     end function simposon_3D_mpi
          

   FUNCTION simposon_2D_mpi(ny,nx,dy,dx,a2)
   implicit none
   integer :: s1,flag,i,j
   REAL*8 :: simposon_2D_mpi
   REAL*8 :: FUN_A,FUN_B,FUN2,FUN4
   INTEGER :: ny,nx
   REAL*8 :: dx,dy
   REAL*8 :: a2(0:ny-1,0:nx-1)
   REAL*8,allocatable :: a1(:)
   
!! n1, n2, n3 should be even numbers i.e. 2*N (N is an arbitrary integer larger than zero);
!!      no volume average.  
!!    double a2[n1][n2][n3] --- 0 ~ n1-1 ,  0 ~ n2-1, 0 ~ n3-1  ;  it assumes the periodic boundary condition here.


   allocate(a1(0:nx-1))
   do i=0,nx-1
        FUN_A=a2(0,i)
        FUN_B=a2(0,i)  !PBC boundary condition !!!!
        FUN2=0.0
        FUN4=a2(1,i)
          s1=1
          !note that if nz==2,the 3D SCFT turns into 2D ,and the simposon integral remains
         ! correct.(nz/2-1=0,do inner loop !)
          do flag=1,ny/2-1,1
           s1=s1+1
           FUN2=FUN2+a2(s1,i)
           s1=s1+1
           FUN4=FUN4+a2(s1,i)
           enddo
          a1(i)=(dy/3.0)*(FUN_A+FUN_B+4.0*FUN4+2.0*FUN2)
        enddo
     
         FUN2=0.0
         FUN4=0.0
        do s1=0,nx-1,2
          FUN2=FUN2+a1(s1)
          FUN4=FUN4+a1(s1+1)
         enddo
     simposon_2D_mpi=(dx/3.0)*(4.0*FUN4+2.0*FUN2)
     deallocate(a1)
     return
     end function simposon_2D_mpi
          
        
    FUNCTION simposon_1D_NR(n0,n1,delta1,a1)
    implicit none
    integer :: s1,flag
    real*8 :: FUN_A,FUN_B,FUN2,FUN4
    real*8 :: simposon_1D_NR
    integer :: n0,n1
    REAL*8 :: a1(n0:n1)
    REAL*8 :: delta1

   !write(*,*) "n0,n1=",n0,n1
   !write(*,*) "delta1=",delta1
   !write(*,*) "a1(n0),a1(n1)=",a1(n0),a1(n1)

    FUN_A=a1(n0)
    FUN_B=a1(n1)
    FUN2=0.0
    FUN4=a1(n0+1)
      s1=n0+1
    do flag=1,(n1-n0)/2-1,1
        s1=s1+1
        FUN2=FUN2+a1(s1)
        s1=s1+1
        FUN4=FUN4+a1(s1)
    enddo
    !    write(*,*) "FUN2=",FUN2
    !    write(*,*) "FUN4=",FUN4
     simposon_1D_NR=(delta1/3.0)*(FUN_A+FUN_B+4.0*FUN4+2.0*FUN2)
     return
    end function simposon_1D_NR

    FUNCTION Sum_sparse_2D(ndim_i,Jij_1D,Jij_1D_size,ai,aj)
    USE global_para,only :node3
    !!!untested!
    implicit none
    integer,intent(in) :: Jij_1D_size
    integer,intent(in) :: ndim_i
    type(node3) :: Jij_1D(0:Jij_1D_size-1)
    real*8,intent(in) :: ai(0:ndim_i-1),aj(0:ndim_i-1)
    integer :: i,j,index_nonzero
    real*8 :: tempx
    real*8 :: Sum_sparse_2D
    real*8,allocatable :: temp_i(:)

    allocate(temp_i(0:ndim_i-1))
    temp_i=0.0

    do index_nonzero=0,Jij_1D_size-1
          i=Jij_1D(index_nonzero)%i
          j=Jij_1D(index_nonzero)%j
          tempx=Jij_1D(index_nonzero)%value
          temp_i(i)=temp_i(i)+tempx*aj(j)
    enddo
    Sum_sparse_2D=0.0
     do i=0,ndim_i-1
        Sum_sparse_2D=Sum_sparse_2D+temp_i(i)*ai(i)
     enddo

    deallocate(temp_i)

    return 

   end function Sum_sparse_2D   



    FUNCTION double_dot_multi(a,b,n)
    implicit none
    integer :: i,j,n
    real*8 double_dot_multi
    real*8,intent(in) :: a(0:n-1,0:n-1),b(0:n-1,0:n-1)  
     double_dot_multi=0.0
    
     do j=0,n-1
       do i=0,n-1 
        double_dot_multi=double_dot_multi+a(i,j)*b(j,i)
       enddo
     enddo

      return
    end function double_dot_multi
















   
     



           
end module utility

           



















                       
     
