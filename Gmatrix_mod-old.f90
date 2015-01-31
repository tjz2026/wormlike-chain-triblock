Module G_matrix_mod
implicit none
      real*8,allocatable :: ksquare(:)
      public :: ksquare

contains

    subroutine init_k_square()
      USE nrtype,only :DP,PI
      USE global_para
      USE control
      USE constants
      USE mmpi
      USE mpi_fftw3_operation
      implicit none
      integer :: i,j,k
      integer :: K_i,K_ii
      integer :: K_j,K_jj
      integer :: K_k,K_kk
      REAL(DP) ::k_x_factor,k_y_factor,k_z_factor
     
      if (allocated(ksquare)) deallocate(ksquare)
      allocate(ksquare(0:Local_size-1)) 


     k=0 
# if defined (Dim2)
     k_x_factor=2.0*PI/(SIDEx*dx)
     k_y_factor=2.0*PI/(SIDEy*dy)
# else  /* Dim2 */
     k_x_factor=2.0*PI/(SIDEx*dx)
     k_y_factor=2.0*PI/(SIDEy*dy)
     k_z_factor=2.0*PI/(SIDEz*dz)
# endif
# if defined (Dim2)
  do K_i=0,local_nx-1 
    if((k_i+local_x_start)<=(SiDEx/2)) then
     k_ii=k_i+local_x_start  
    else 
     k_ii=(k_i+local_x_start)-SIDEx 
    endif
      do k_j=0,SIDEy-1
          if(K_j<=(SIDEy/2)) then 
            k_jj=K_j
          else
            K_jj=K_j-SIDEy
          endif
# else  /* Dim2 */
  do K_i=0,local_nx-1 
    if((k_i+local_x_start)<=(SiDEx/2)) then
     k_ii=k_i+local_x_start  
    else 
     k_ii=(k_i+local_x_start)-SIDEx 
    endif
      do k_j=0,SIDEy-1
          if(K_j<=(SIDEy/2)) then 
            k_jj=K_j
          else
            K_jj=K_j-SIDEy
          endif
  
             do k_k=0,SIDEz-1
               if(K_k<=(SIDEz/2)) then
               k_kk=K_k
               else
               k_kk=K_k-SIDEz
               endif
# endif  /* Dim2 */

# if defined (Dim2)
                ksquare(k)=dsqrt((k_ii*k_x_factor)**2 + &
                                  (k_jj*k_y_factor)**2)
# else  /* Dim2 */
                ksquare(k)=dsqrt((k_ii*k_x_factor)**2 + &
                                 (k_jj*k_y_factor)**2 + &
                                 (k_kk*k_z_factor)**2)
# endif  /* Dim2 */

                k=k+1   
# if defined (Dim2)               
          enddo   !enddo k_y
        enddo    !enddo k_x
# else  /* Dim2 */
          enddo  !enddo k_z
        enddo   !enddo k_y
      enddo    !enddo k_x
# endif  /* Dim2 */
   end subroutine init_k_square


!test the basic size imformation of matrix G
# if defined  (Test)
subroutine test_inverse_sparse_GAB()
 USE nrtype,only :DP,PI
 USE global_para
 USE mpi
 USE control
 USE constants
 USE utility
 USE mmpi
 USE mpi_fftw3_operation
 USE matrix_inverse
 USE sparse_MKL
 implicit none
 integer :: i,j,k
 integer :: K_i,K_ii
 integer :: K_j,K_jj
 integer :: K_k,K_kk
 integer :: G_index
  REAL(DP) :: kapa_temp
 REAL(DP),allocatable :: G_R(:,:)
 REAL(DP),allocatable :: G_I(:,:)
 REAL(DP) ::k_x_factor,k_y_factor,k_z_factor
 integer,parameter :: namax=1800
real(DP) :: coeff_deltaij
real(DP) :: coeff_Rij
integer :: istat2
integer :: N_low1,N_low2
integer, allocatable :: N_low_k1(:),N_low_k2(:)
integer :: loc(1)
integer job(8)
integer info

! we do not need to have set up the G matrix for A and B both,just one will do.
max_M_Bar=max(M_Bar_B,M_Bar_A,M_Bar_C)

allocate(G_R(0:max_M_Bar-1,0:max_M_Bar-1))
allocate(G_I(0:max_M_Bar-1,0:max_M_Bar-1))
allocate(N_low_k1(1:LOCAL_SIZE))
allocate(N_low_k2(1:LOCAL_SIZE))

!note that the array is storaged in cloumn before row in fortran
!so ,it is different from the c++ version code.
!always remember,otherwise it would be terribly wrong!!!!!!!!!!!!!!!!!!

G_R=0.0
G_I=0.0

# if defined (Dim2)
k_x_factor=2.0*PI/(SIDEx*dx)
k_y_factor=2.0*PI/(SIDEy*dy)
# else  /* Dim2 */
k_x_factor=2.0*PI/(SIDEx*dx)
k_y_factor=2.0*PI/(SIDEy*dy)
k_z_factor=2.0*PI/(SIDEz*dz)
# endif



do G_index=1,8

if(G_index==1) then
!GA1  
coeff_deltaij=1.0
coeff_Rij=1.0
kapa_temp=kapa_A
else if(G_index==2) then
!GA2
coeff_deltaij=1.5
coeff_Rij=1.0
kapa_temp=kapa_A
else if(G_index==3) then
!GA3
coeff_deltaij=11.0/6.0
coeff_Rij=1.0
kapa_temp=kapa_A
else if(G_index==4) then
!GA3_star
coeff_deltaij=11.0/6.0
coeff_Rij=-1.0
kapa_temp=kapa_A
else if(G_index==5) then
!GB3
coeff_deltaij=11.0/6.0
coeff_Rij=1.0
kapa_temp=kapa_B

else if(G_index==6) then
!GB3_star
coeff_deltaij=11.0/6.0
coeff_Rij=-1.0
kapa_temp=kapa_B

else if(G_index==7) then
!GB2_star
coeff_deltaij=1.5
coeff_Rij=-1.0
kapa_temp=kapa_B
else if(G_index==8) then
!GB1_star
coeff_deltaij=1.0
coeff_Rij=-1.0
kapa_temp=kapa_B

endif



G_R=0.0
G_I=0.0


!k counter
k=0
!!the first dimension must be z 

# if defined (Dim2)
do K_i=0,local_nx-1 
  if((k_i+local_x_start)<=(SiDEx/2)) then
   k_ii=k_i+local_x_start  
  else 
   k_ii=(k_i+local_x_start)-SIDEx 
  endif
    do k_j=0,SIDEy-1
        if(K_j<=(SIDEy/2)) then 
          k_jj=K_j
        else
          K_jj=K_j-SIDEy
        endif
# else  /* Dim2 */
do K_i=0,local_nx-1 
  if((k_i+local_x_start)<=(SiDEx/2)) then
   k_ii=k_i+local_x_start  
  else 
   k_ii=(k_i+local_x_start)-SIDEx 
  endif
    do k_j=0,SIDEy-1
        if(K_j<=(SIDEy/2)) then 
          k_jj=K_j
        else
          K_jj=K_j-SIDEy
        endif

           do k_k=0,SIDEz-1
             if(K_k<=(SIDEz/2)) then
             k_kk=K_k
             else
             k_kk=K_k-SIDEz
             endif
# endif  /* Dim2 */
              !enter the i,j of spehrical harmonic basis index loop
              G_R=0.0
              G_I=0.0
              do j=0,max_M_Bar-1
              !G_R has only the digonal element
              G_R(j,j)=coeff_deltaij + ds*Loa*basis_SPH(j)%x*(basis_SPH(j)%x+1)/(2.0*kapa_temp)
                do i=0,max_M_Bar-1
# if defined (Dim2)
              G_I(i,j)=coeff_Rij*ds*(matrix_Rx(i,j)*(k_ii*k_x_factor) + &
                       matrix_Ry(i,j)*(k_jj*k_y_factor))
# else  /* Dim2 */
              G_I(i,j)=coeff_Rij*ds*(matrix_Rx(i,j)*(k_ii*k_x_factor) + &
                       matrix_Ry(i,j)*(k_jj*k_y_factor) + &
                       matrix_Rz(i,j)*(k_kk*k_z_factor))
# endif  /* Dim2 */

                enddo
               enddo

              !inverse the complex matrix using the lapack lib directly
              call complex_matrix_inverse(G_R,G_I,max_M_Bar)
              
                  k=k+1 
                 if((k_I*(SIDEy)+k_j+1)>local_size .or. & 
                  (k_I*(SIDEy)+k_j+1)<1) then
                  write(*,*) "k index exceeds",k_I*(SIDEy)+k_j+1
                  stop
                 endif
                
                  N_low1=0
                  N_low2=0
                  call get_symmetric_matrix_size_dp(G_R,max_M_Bar,max_M_Bar,N_low1)
                  call get_symmetric_matrix_size_dp(G_I,max_M_Bar,max_M_Bar,N_low2)
                 ! call get_symmetric_matrix_size_dp(G_I,max_M_Bar,max_M_Bar,N_low)
                 !call get_sparse_size_one_based_dp(G_R,max_M_Bar,max_M_Bar,N_low)
                   N_low_k1(k)=N_low1
                   N_low_k2(k)=N_low2

                 !N_low is precalculated so that namax is determined 

        job(1)=0
        job(2)=1
        job(3)=1
        job(4)=0
        job(5)=namax
        job(6)=1


               do j=0,max_M_Bar-1
                  do i=0,max_M_Bar-1    
                   if(j>i) then
                   G_R(i,j)=0.0
                   G_I(i,j)=0.0
                   endif
                   enddo
               enddo
   
 !       call convert_symmetric_matrix_to_array_cmplx(G_z,max_M_Bar,max_M_Bar,N_low,sy_values_z,  &
 !       sy_columns_z,sy_rowIndex_z) !this routine has an issue when N_low doesn't match the actual nonzero
 !       number of the matrix.Needs to be fixed.
!                 
!       
# if defined (Dim2)               
       enddo   !enddo k_y
      enddo    !enddo k_x
# else  /* Dim2 */
        enddo  !enddo k_z
       enddo   !enddo k_y
      enddo    !enddo k_x
# endif  /* Dim2 */

      if(k/=local_size) then
      write(*,*) "error,total k doesn't match local_size!!"
      write(*,*) "total k=",k
      stop
      endif

      write(*,*) " maxium number of N_lowk1 is",maxval(N_low_k1),"on",myid
      write(*,*) " loc maxium number of N_lowk1 is",maxloc(N_low_k1),"on",myid
      write(*,*) " minum number of N_lowk1 is",minval(N_low_k1),"on",myid
      write(*,*) " loc minum number of N_lowk1 is",minloc(N_low_k1),"on",myid
      write(*,*) " total number of N_lowk1 is",sum(N_low_k1),"on",myid

      write(*,*) " maxium number of N_lowk2 is",maxval(N_low_k2),"on",myid
      write(*,*) " loc maxium number of N_lowk2 is",maxloc(N_low_k2),"on",myid
      write(*,*) " minum number of N_lowk2 is",minval(N_low_k2),"on",myid
      write(*,*) " loc minum number of N_lowk2 is",minloc(N_low_k2),"on",myid
      write(*,*) " total number of N_lowk2 is",sum(N_low_k2),"on",myid

enddo  !enddo G_index
   
     deallocate(G_R)
     deallocate(G_I)
     deallocate(N_low_k1)
     deallocate(N_low_k2)
     call mp_barrier() 

    end subroutine test_inverse_sparse_GAB
# endif /* Test */    

subroutine init_inverse_sparse_GAB()
!!be careful,this subroutine is going to be called ever time 
!the box size is changed
 USE nrtype,only :DP,PI
 USE global_para
 USE mpi
 USE control
 USE constants
 USE utility
 USE mmpi
 USE mpi_fftw3_operation
 USE matrix_inverse
 USE sparse_MKL
 implicit none
 integer :: i,j,k
 integer :: K_i,K_ii
 integer :: K_j,K_jj
 integer :: K_k,K_kk
 integer :: G_index
  REAL(DP) :: kapa_temp
 REAL(DP),allocatable :: G_R(:,:)
 REAL(DP),allocatable :: G_I(:,:)
 DOUBLE PRECISION,allocatable :: sy_values_R(:)
 DOUBLE PRECISION,allocatable :: sy_values_I(:)
 integer,allocatable :: sy_columns_R(:),sy_rowIndex_R(:)
 integer,allocatable :: sy_columns_I(:),sy_rowIndex_I(:)
 REAL(DP) ::k_x_factor,k_y_factor,k_z_factor
# if defined (Dim2)
  !integer,parameter :: namaxA=900
  !integer,parameter :: namaxA=1950
  integer,parameter :: namaxA=900
 !integer,parameter :: namaxA=256
 !integer,parameter :: namaxB=900
  !integer,parameter :: namaxB=625
  integer,parameter :: namaxB=900
# else  /* Dim2 */
 !integer,parameter :: namaxA=3800
 integer,parameter :: namaxA=1800
 !integer,parameter :: namaxB=1800
 !integer,parameter :: namaxB=625
 integer,parameter :: namaxB=1800
# endif /* Dim2 */
! L=6,M=49, namaxA=350
! L=7,M=64, namaxA=560
! L=8,M=81,namaxA=881(my suggest is set as 900)
! L=9,M=100,namaxA=1325
! L=10,M=121,namaxA=1921
! for 3D case:
! L=6,M=49,namaxA=640
! L=7,M=64,namaxA=1100
! L=8,M=81,namaxA=1800
! L=9,M=100,namaxA=2575
! L=10,M=121,namaxA=3751
 integer,allocatable :: ija(:)
 REAL(DP),allocatable :: sa(:)
 integer,allocatable :: namax_final_R1(:)
 integer,allocatable :: namax_final_I1(:)
real(DP) :: coeff_deltaij
real(DP) :: coeff_Rij
integer :: istat2
integer :: loc(1)
integer job(8)
integer info

! we do not need to have set up the G matrix for A and B both,just one will do.
max_M_Bar=max(M_Bar_B,M_Bar_A,M_Bar_C)

  allocate(G_R(0:max_M_Bar-1,0:max_M_Bar-1))
  allocate(G_I(0:max_M_Bar-1,0:max_M_Bar-1))
          allocate(sy_values_R(1:namaxA))
          allocate(sy_columns_R(1:namaxA))
          allocate(sy_rowIndex_R(1:max_M_Bar+1))

          allocate(sy_values_I(1:namaxA))
          allocate(sy_columns_I(1:namaxA))
          allocate(sy_rowIndex_I(1:max_M_Bar+1))

  if(allocated(R_G3_inv_val_blk)) then
  else
  allocate(R_G3_inv_val_blk(1:namaxA,1:LOCAL_SIZE,1:3))
  allocate(R_G3_inv_col_blk(1:namaxA,1:LOCAL_SIZE,1:3))
  allocate(R_G3_inv_row_blk(1:max_M_Bar+1,1:LOCAL_SIZE,1:3))

  allocate(I_G3_inv_val_blk(1:namaxA,1:LOCAL_SIZE,1:3))
  allocate(I_G3_inv_col_blk(1:namaxA,1:LOCAL_SIZE,1:3))
  allocate(I_G3_inv_row_blk(1:max_M_Bar+1,1:LOCAL_SIZE,1:3))

  allocate(R_G3star_inv_val_blk(1:namaxA,1:LOCAL_SIZE,1:3))
  allocate(R_G3star_inv_col_blk(1:namaxA,1:LOCAL_SIZE,1:3))
  allocate(R_G3star_inv_row_blk(1:max_M_Bar+1,1:LOCAL_SIZE,1:3))

  allocate(I_G3star_inv_val_blk(1:namaxA,1:LOCAL_SIZE,1:3))
  allocate(I_G3star_inv_col_blk(1:namaxA,1:LOCAL_SIZE,1:3))
  allocate(I_G3star_inv_row_blk(1:max_M_Bar+1,1:LOCAL_SIZE,1:3))
  endif

# if defined (Dim2)
k_x_factor=2.0*PI/(SIDEx*dx)
k_y_factor=2.0*PI/(SIDEy*dy)
# else  /* Dim2 */
k_x_factor=2.0*PI/(SIDEx*dx)
k_y_factor=2.0*PI/(SIDEy*dy)
k_z_factor=2.0*PI/(SIDEz*dz)
# endif

do G_index=3,8

if(G_index==1) then
!GA1  
coeff_deltaij=1.0
coeff_Rij=1.0
kapa_temp=kapa_A
else if(G_index==2) then
!GA2
coeff_deltaij=1.5
coeff_Rij=1.0
kapa_temp=kapa_A
else if(G_index==3) then
!GA3
coeff_deltaij=11.0/6.0
coeff_Rij=1.0
kapa_temp=kapa_A
else if(G_index==4) then
!GB3
coeff_deltaij=11.0/6.0
coeff_Rij=1.0
kapa_temp=kapa_B
else if(G_index==5) then
!GC3
coeff_deltaij=11.0/6.0
coeff_Rij=1.0
kapa_temp=kapa_C
else if(G_index==6) then
!GA3_star
coeff_deltaij=11.0/6.0
coeff_Rij=-1.0
kapa_temp=kapa_A
else if(G_index==7) then
!GB3_star
coeff_deltaij=11.0/6.0
coeff_Rij=-1.0
kapa_temp=kapa_B
else if(G_index==8) then
!GC3_star
coeff_deltaij=11.0/6.0
coeff_Rij=-1.0
kapa_temp=kapa_C
else if(G_index==9) then
!GC2_star
coeff_deltaij=1.5
coeff_Rij=-1.0
kapa_temp=kapa_C
else if(G_index==10) then
!GC1_star
coeff_deltaij=1.0
coeff_Rij=-1.0
kapa_temp=kapa_C

endif

G_R=0.0
G_I=0.0

!k counter
k=0
!!the first dimension must be z 
# if defined (Dim2)
do K_i=0,local_nx-1 
  if((k_i+local_x_start)<=(SiDEx/2)) then
   k_ii=k_i+local_x_start  
  else 
   k_ii=(k_i+local_x_start)-SIDEx 
  endif
    do k_j=0,SIDEy-1
        if(K_j<=(SIDEy/2)) then 
          k_jj=K_j
        else
          K_jj=K_j-SIDEy
        endif
# else  /* Dim2 */
do K_i=0,local_nx-1 
  if((k_i+local_x_start)<=(SiDEx/2)) then
   k_ii=k_i+local_x_start  
  else 
   k_ii=(k_i+local_x_start)-SIDEx 
  endif
    do k_j=0,SIDEy-1
        if(K_j<=(SIDEy/2)) then 
          k_jj=K_j
        else
          K_jj=K_j-SIDEy
        endif

           do k_k=0,SIDEz-1
             if(K_k<=(SIDEz/2)) then
             k_kk=K_k
             else
             k_kk=K_k-SIDEz
             endif
# endif  /* Dim2 */

              !enter the i,j of spehrical harmonic basis index loop
              G_R=0.0
              G_I=0.0
        do j=0,max_M_Bar-1
              !G_R has only the digonal element
         G_R(j,j)=coeff_deltaij + ds*Loa*basis_SPH(j)%x*(basis_SPH(j)%x+1)/(2.0*kapa_temp)
                do i=0,max_M_Bar-1
# if defined (Dim2)
              G_I(i,j)=coeff_Rij*ds*(matrix_Rx(i,j)*(k_ii*k_x_factor) + &
                       matrix_Ry(i,j)*(k_jj*k_y_factor))
# else  /* Dim2 */
              G_I(i,j)=coeff_Rij*ds*(matrix_Rx(i,j)*(k_ii*k_x_factor) + &
                       matrix_Ry(i,j)*(k_jj*k_y_factor) + &
                       matrix_Rz(i,j)*(k_kk*k_z_factor))
# endif  /* Dim2 */

                enddo
               enddo

              call complex_matrix_inverse(G_R,G_I,max_M_Bar)
                
                  k=k+1 
                 if((k_I*(SIDEy)+k_j+1)>local_size .or. & 
                  (k_I*(SIDEy)+k_j+1)<1) then
                  write(*,*) "k index exceeds",k_I*(SIDEy)+k_j+1
                  stop
                 endif

        job(1)=0
        job(2)=1
        job(3)=1
        job(4)=0

        job(5)=namaxA
               do j=0,max_M_Bar-1
                  do i=0,max_M_Bar-1    
                   if(j>i) then
                   G_R(i,j)=0.0
                   G_I(i,j)=0.0
                   endif
                   enddo
               enddo
        job(6)=1
   
       call mkl_ddnscsr(job,max_M_Bar,max_M_Bar,G_R,max_M_Bar,sy_values_R,sy_columns_R, &
        sy_rowIndex_R,info)
       call mkl_ddnscsr(job,max_M_Bar,max_M_Bar,G_I,max_M_Bar,sy_values_I,sy_columns_I, &
        sy_rowIndex_I,info)
           
 !       call convert_symmetric_matrix_to_array_cmplx(G_z,max_M_Bar,max_M_Bar,N_low,sy_values_z,  &
 !       sy_columns_z,sy_rowIndex_z) !this routine has an issue when N_low doesn't match the actual nonzero
 !       number of the matrix.Needs to be fixed.

                 if(G_index<=5) then
                   i=mod(G_index,3)+1
                 R_G3_inv_val(:,k,i)=sy_values_R(:)
                 R_G3_inv_col(:,k,i)=sy_columns_R(:)
                 R_G3_inv_row(:,k,i)=sy_rowIndex_R(:)

                 I_G3_inv_val(:,k,i)=sy_values_I(:)
                 I_G3_inv_col(:,k,i)=sy_columns_I(:)
                 I_G3_inv_row(:,k,i)=sy_rowIndex_I(:)
                 else if(G_index>=6) then
                   i=mod(G_index,3)+1
                 R_G3star_inv_val(:,k,i)=sy_values_R(:)
                 R_G3star_inv_col(:,k,i)=sy_columns_R(:)
                 R_G3star_inv_row(:,k,i)=sy_rowIndex_R(:)

                 I_G3star_inv_val(:,k,i)=sy_values_I(:)
                 I_G3star_inv_col(:,k,i)=sy_columns_I(:)
                 I_G3star_inv_row(:,k,i)=sy_rowIndex_I(:)

                 endif
                 
# if defined (Dim2)               
       enddo   !enddo k_y
      enddo    !enddo k_x
# else  /* Dim2 */
        enddo  !enddo k_z
       enddo   !enddo k_y
      enddo    !enddo k_x
# endif  /* Dim2 */

      if(k/=local_size) then
      write(*,*) "error,total k doesn't match local_size!!"
      write(*,*) "total k=",k
      stop
      endif

enddo  !enddo G_index
   
     deallocate(G_R)
     deallocate(G_I)
     deallocate(sy_values_R)
     deallocate(sy_columns_R)
     deallocate(sy_rowIndex_R)
     deallocate(sy_values_I)
     deallocate(sy_columns_I)
     deallocate(sy_rowIndex_I)
     call mp_barrier() 
      if(myid==0) then
        write(*,*) "init G matrix memeory not optimized"
       endif
      call mp_barrier()
    end subroutine init_inverse_sparse_GAB

subroutine init_G_mat(G_mat,M_Bar,G_index,K_i,k_j,k_k)
!!be careful,this subroutine is going to be called ever time 
!the box size is changed
 USE nrtype,only :DP,PI
 USE global_para
 USE control
 USE constants
 USE utility
 USE mmpi
 USE mpi_fftw3_operation
 USE matrix_inverse
 USE sparse_MKL
 use m_sparse
 implicit none
 integer :: i,j,k
 integer,intent(in) :: K_i,K_j
 integer,optional,intent(in) :: k_k
 integer :: K_ii
 integer :: K_jj
 integer :: K_kk
 integer,intent(in) :: G_index
 REAL(DP) :: kapa_temp
 !real(dp),parameter :: eps=1.0e-20
 real(dp),parameter :: eps=0.0d0
 integer,intent(in) :: M_Bar
 Complex(DP),intent(inout) :: G_mat(0:M_Bar-1,0:M_Bar-1)
 REAL(DP) ::k_x_factor,k_y_factor,k_z_factor
real(DP) :: coeff_deltaij
real(DP) :: coeff_Rij,G_real,G_imag
real(DP) :: t21,t1,t2,t3,t4,t34
integer :: istat2
integer :: loc(1)
integer job(8)
integer info

! we do not need to have set up the G matrix for A and B both,just one will do.
!max_M_Bar=max(M_Bar_B,M_Bar_A)

# if defined (Dim2)
k_x_factor=2.0*PI/(SIDEx*dx)
k_y_factor=2.0*PI/(SIDEy*dy)
# else  /* Dim2 */
k_x_factor=2.0*PI/(SIDEx*dx)
k_y_factor=2.0*PI/(SIDEy*dy)
k_z_factor=2.0*PI/(SIDEz*dz)
# endif

if(G_index==1) then
!GA1  
coeff_deltaij=1.0
coeff_Rij=1.0
kapa_temp=kapa_A
else if(G_index==2) then
!GA2
coeff_deltaij=1.5
coeff_Rij=1.0
kapa_temp=kapa_A
else if(G_index==3) then
!GA3
coeff_deltaij=11.0/6.0
coeff_Rij=1.0
kapa_temp=kapa_A
else if(G_index==4) then
!GB3
coeff_deltaij=11.0/6.0
coeff_Rij=1.0
kapa_temp=kapa_B
else if(G_index==5) then
!GC3
coeff_deltaij=11.0/6.0
coeff_Rij=1.0
kapa_temp=kapa_C
else if(G_index==6) then
!GA3_star
coeff_deltaij=11.0/6.0
coeff_Rij=-1.0
kapa_temp=kapa_A
else if(G_index==7) then
!GB3_star
coeff_deltaij=11.0/6.0
coeff_Rij=-1.0
kapa_temp=kapa_B
else if(G_index==8) then
!GC3_star
coeff_deltaij=11.0/6.0
coeff_Rij=-1.0
kapa_temp=kapa_C
else if(G_index==9) then
!GC2_star
coeff_deltaij=1.5
coeff_Rij=-1.0
kapa_temp=kapa_C
else if(G_index==10) then
!GC1_star
coeff_deltaij=1.0
coeff_Rij=-1.0
kapa_temp=kapa_C
endif

k=0
!!the first dimension must be z 
# if defined (Dim2)
         if((k_i+local_x_start)<=(SiDEx/2)) then
          k_ii=k_i+local_x_start  
         else 
          k_ii=(k_i+local_x_start)-SIDEx 
         endif
         if(K_j<=(SIDEy/2)) then 
          k_jj=K_j
         else
          K_jj=K_j-SIDEy
         endif
# else  /* Dim2 */
         if((k_i+local_x_start)<=(SiDEx/2)) then
          k_ii=k_i+local_x_start  
         else 
          k_ii=(k_i+local_x_start)-SIDEx 
         endif
         if(K_j<=(SIDEy/2)) then 
          k_jj=K_j
         else
           K_jj=K_j-SIDEy
         endif
         if(K_k<=(SIDEz/2)) then
         k_kk=K_k
         else
         k_kk=K_k-SIDEz
         endif
# endif  /* Dim2 */
        !enter the i,j of spehrical harmonic basis index loop
        do j=0,M_Bar-1
                do i=0,M_Bar-1
# if defined (Dim2)
              G_imag=coeff_Rij*ds*(matrix_Rx(i,j)*(k_ii*k_x_factor+eps) + &
                       matrix_Ry(i,j)*(k_jj*k_y_factor+eps))
# else  /* Dim2 */
              G_imag=coeff_Rij*ds*(matrix_Rx(i,j)*(k_ii*k_x_factor+eps) + &
                       matrix_Ry(i,j)*(k_jj*k_y_factor+eps) + &
                       matrix_Rz(i,j)*(k_kk*k_z_factor+eps))
# endif  /* Dim2 */
              G_mat(i,j)=dcmplx(0.0d0,G_imag)  
                enddo
               enddo

        do j=0,M_Bar-1
         G_real=coeff_deltaij + ds*Loa*basis_SPH(j)%x*(basis_SPH(j)%x+1)/(2.0*kapa_temp)
         G_mat(j,j)=dcmplx(G_real,aimag(G_mat(j,j)))
        enddo 
   
    end subroutine init_G_mat

end module  G_matrix_mod


