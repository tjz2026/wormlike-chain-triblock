# if defined (SLAB) 
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
!# if defined  (Test)
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
 use m_sparse
 implicit none
 integer :: i,j,k
 integer :: K_i,K_ii
 integer :: K_j,K_jj
 integer :: K_k,K_kk
 integer :: G_index
  REAL(DP) :: kapa_temp
 REAL(DP),allocatable :: G_R(:,:)
 REAL(DP),allocatable :: G_I(:,:)
 Complex(DP),allocatable :: G_Z(:,:)
 REAL(DP) ::k_x_factor,k_y_factor,k_z_factor
 integer,parameter :: namax=3330
! L_bar=4,namax=
! L_bar=5,namax=
! L_bar=6,namax=
! L_bar=7,namax=
! L_bar=8,namax=371
! L_bar=9,namax=
! L_bar=10,namax=
 real(dp),parameter :: eps=0.0d0
real(DP) :: coeff_deltaij
real(DP) :: coeff_Rij
integer :: istat2
integer :: N_low1,N_low2
integer :: N_low3
integer, allocatable :: N_low_k1(:),N_low_k2(:)
integer, allocatable :: N_low_k3(:)
integer :: loc(1)
integer job(8)
integer info

! we do not need to have set up the G matrix for A and B both,just one will do.
max_M_Bar=max(M_Bar_B,M_Bar_A)

allocate(G_R(0:max_M_Bar-1,0:max_M_Bar-1))
allocate(G_I(0:max_M_Bar-1,0:max_M_Bar-1))
allocate(G_Z(0:max_M_Bar-1,0:max_M_Bar-1))
allocate(N_low_k1(1:LOCAL_SIZE))
allocate(N_low_k2(1:LOCAL_SIZE))
allocate(N_low_k3(1:LOCAL_SIZE))



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
!k_y_factor=2.0*PI/(SIDEy*dy+0.00000001d0)
!k_z_factor=2.0*PI/(SIDEz*dz-0.00000001d0)
# endif



!do G_index=1,8
 G_index=5

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
              G_I(i,j)=coeff_Rij*ds*(matrix_Rx(i,j)*(k_ii*k_x_factor+eps) + &
                       matrix_Ry(i,j)*(k_jj*k_y_factor+eps))
# else  /* Dim2 */
              G_I(i,j)=coeff_Rij*ds*(matrix_Rx(i,j)*(k_ii*k_x_factor+eps) + &
                       matrix_Ry(i,j)*(k_jj*k_y_factor+eps) + &
                       matrix_Rz(i,j)*(k_kk*k_z_factor+eps))
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
                 !  N_low_k1(k)=N_low1
                   N_low_k2(k)=N_low2

                 !N_low is precalculated so that namax is determined 
        job(1)=0
        job(2)=1
        job(3)=1
        job(4)=0
        job(5)=namax
        job(6)=1


                  G_Z=dcmplx(0.0d0,0.0d0) 
               do j=0,max_M_Bar-1
                  do i=0,max_M_Bar-1    
                   if(j>i) then
                   G_R(i,j)=0.0
                   G_I(i,j)=0.0
                   endif
                   G_Z(i,j)=dcmplx(G_R(i,j),G_I(i,j))
                   enddo
               enddo
   
              call get_sparse_size(max_M_Bar,max_M_Bar,N_low3,G_Z)
                   N_low_k3(k)=N_low3
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

     if(myid==2) then  
     open(unit=33,file='N_low3.dat',status='replace')
        do i=1,local_size
           write(33,*) i,N_low_k3(i)
        enddo
     close(33)
     endif
!      write(*,*) " loc maxium number of N_lowk1 is",maxloc(N_low_k1),"on",myid
!      write(*,*) " minum number of N_lowk1 is",minval(N_low_k1),"on",myid
!      write(*,*) " loc minum number of N_lowk1 is",minloc(N_low_k1),"on",myid
!      write(*,*) " total number of N_lowk1 is",sum(N_low_k1),"on",myid

!      write(*,*) " maxium number of N_lowk2 is",maxval(N_low_k2),"on",myid
!      write(*,*) " loc maxium number of N_lowk2 is",maxloc(N_low_k2),"on",myid
!      write(*,*) " minum number of N_lowk2 is",minval(N_low_k2),"on",myid
!      write(*,*) " loc minum number of N_lowk2 is",minloc(N_low_k2),"on",myid
!      write(*,*) " total number of N_lowk2 is",sum(N_low_k2),"on",myid
       
      write(*,*) " maxium number of N_lowk3 is",maxval(N_low_k3),"on",myid
      write(*,*) " minum number of N_lowk3 is",minval(N_low_k3),"on",myid
      write(*,*) " ava number of N_lowk3 is",sum(N_low_k3)/local_size,"on",myid
       


!enddo  !enddo G_index
   
     deallocate(G_R)
     deallocate(G_I)
     deallocate(G_Z)
     deallocate(N_low_k1)
     deallocate(N_low_k2)
     deallocate(N_low_k3)
     call mp_barrier() 

    end subroutine test_inverse_sparse_GAB
!# endif /* Test */    


subroutine init_inverse_sparse_GAB()
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
 integer :: K_i,K_ii
 integer :: K_j,K_jj
 integer :: K_k,K_kk
 integer :: G_index
  REAL(DP) :: kapa_temp
 real(dp),parameter :: eps=0.0d0
 REAL(DP),allocatable :: G_RA(:,:)
 REAL(DP),allocatable :: G_IA(:,:)
 REAL(DP),allocatable :: G_RB(:,:)
 REAL(DP),allocatable :: G_IB(:,:)

 Complex(DP),allocatable :: G_A(:,:)
 Complex(DP),allocatable :: G_B(:,:)
 Complex(dp),allocatable :: sy_values_A(:)
 Complex(dp),allocatable :: sy_values_B(:)
 integer,allocatable :: sy_columns_A(:),sy_rowIndex_A(:)
 integer,allocatable :: sy_columns_B(:),sy_rowIndex_B(:)
 REAL(DP) ::k_x_factor,k_y_factor,k_z_factor
# if defined (Dim2)
 integer,parameter :: namaxA=900
 integer,parameter :: namaxB=900
# else  /* Dim2 */
 integer,parameter :: namaxA=3330
 integer,parameter :: namaxB=3330
! integer,parameter :: namaxA=371
# endif /* Dim2 */
 integer,allocatable :: ija(:)
 REAL(DP),allocatable :: sa(:)
 integer,allocatable :: namax_final_R1(:)
 integer,allocatable :: namax_final_I1(:)
real(DP) :: coeff_deltaij
real(DP) :: coeff_Rij
real(DP) :: t21,t1,t2,t3,t4,t34
integer :: istat2
integer :: loc(1)
integer job(8)
integer info

! we do not need to have set up the G matrix for A and B both,just one will do.
!max_M_Bar=max(M_Bar_B,M_Bar_A)

allocate(G_RA(0:M_Bar_A-1,0:M_Bar_A-1))
allocate(G_IA(0:M_Bar_A-1,0:M_Bar_A-1))
allocate(G_RB(0:M_Bar_B-1,0:M_Bar_B-1))
allocate(G_IB(0:M_Bar_B-1,0:M_Bar_B-1))

allocate(G_A(0:M_Bar_A-1,0:M_Bar_A-1))
allocate(G_B(0:M_Bar_B-1,0:M_Bar_B-1))

          allocate(sy_values_A(1:namaxA))
          allocate(sy_columns_A(1:namaxA))
          allocate(sy_rowIndex_A(1:M_Bar_A+1))
          allocate(sy_values_B(1:namaxB))
          allocate(sy_columns_B(1:namaxB))
          allocate(sy_rowIndex_B(1:M_Bar_B+1))

!complex G matrix for CSR format
 ! allocate(GA1_inverse_val(1:namaxA,1:LOCAL_SIZE))
 ! allocate(GA1_inverse_col(1:namaxA,1:LOCAL_SIZE))
 ! allocate(GA1_inverse_row(1:M_Bar_A+1,1:LOCAL_SIZE))

 ! allocate(GA2_inverse_val(1:namaxA,1:LOCAL_SIZE))
 ! allocate(GA2_inverse_col(1:namaxA,1:LOCAL_SIZE))
 ! allocate(GA2_inverse_row(1:M_Bar_A+1,1:LOCAL_SIZE))

  allocate(GA3_inverse_val(1:namaxA,1:LOCAL_SIZE))
  allocate(GA3_inverse_col(1:namaxA,1:LOCAL_SIZE))
  allocate(GA3_inverse_row(1:M_Bar_A+1,1:LOCAL_SIZE))

  allocate(GA3_star_inverse_val(1:namaxA,1:LOCAL_SIZE))
  allocate(GA3_star_inverse_col(1:namaxA,1:LOCAL_SIZE))
  allocate(GA3_star_inverse_row(1:M_Bar_A+1,1:LOCAL_SIZE))

  allocate(GB3_inverse_val(1:namaxB,1:LOCAL_SIZE))
  allocate(GB3_inverse_col(1:namaxB,1:LOCAL_SIZE))
  allocate(GB3_inverse_row(1:M_Bar_B+1,1:LOCAL_SIZE))

  allocate(GB3_star_inverse_val(1:namaxB,1:LOCAL_SIZE))
  allocate(GB3_star_inverse_col(1:namaxB,1:LOCAL_SIZE))
  allocate(GB3_star_inverse_row(1:M_Bar_B+1,1:LOCAL_SIZE))

 ! allocate(GB2_star_inverse_val(1:namaxB,1:LOCAL_SIZE))
 ! allocate(GB2_star_inverse_col(1:namaxB,1:LOCAL_SIZE))
 ! allocate(GB2_star_inverse_row(1:M_Bar_B+1,1:LOCAL_SIZE))

 ! allocate(GB1_star_inverse_val(1:namaxB,1:LOCAL_SIZE))
 ! allocate(GB1_star_inverse_col(1:namaxB,1:LOCAL_SIZE))
 ! allocate(GB1_star_inverse_row(1:M_Bar_B+1,1:LOCAL_SIZE))

!note that the array is storaged in cloumn before row in fortran
!so ,it is different from the c++ version code.
!always remember,otherwise it would be terribly wrong!!!!!!!!!!!!!!!!!!

# if defined (Dim2)
k_x_factor=2.0*PI/(SIDEx*dx)
k_y_factor=2.0*PI/(SIDEy*dy)
# else  /* Dim2 */
k_x_factor=2.0*PI/(SIDEx*dx)
k_y_factor=2.0*PI/(SIDEy*dy)
k_z_factor=2.0*PI/(SIDEz*dz)
# endif

!do G_index=1,8
do G_index=3,6
!if(G_index==1) then
!!GA1  
!coeff_deltaij=1.0
!coeff_Rij=1.0
!kapa_temp=kapa_A
!else if(G_index==2) then
!!GA2
!coeff_deltaij=1.5
!coeff_Rij=1.0
!kapa_temp=kapa_A
!else if(G_index==3) then
if(G_index==3) then
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
!else if(G_index==7) then
!!GB2_star
!coeff_deltaij=1.5
!coeff_Rij=-1.0
!kapa_temp=kapa_B
!else if(G_index==8) then
!!GB1_star
!coeff_deltaij=1.0
!coeff_Rij=-1.0
!kapa_temp=kapa_B
endif

t21=0.0d0

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
     if(G_index<=4) then
              G_RA=0.0
              G_IA=0.0
            t3=MPI_Wtime()
        do j=0,M_Bar_A-1
              !G_R has only the digonal element
         G_RA(j,j)=coeff_deltaij + ds*Loa*basis_SPH(j)%x*(basis_SPH(j)%x+1)/(2.0*kapa_temp)
                do i=0,M_Bar_A-1
# if defined (Dim2)
              G_IA(i,j)=coeff_Rij*ds*(matrix_Rx(i,j)*(k_ii*k_x_factor+eps) + &
                       matrix_Ry(i,j)*(k_jj*k_y_factor+eps))
# else  /* Dim2 */
              G_IA(i,j)=coeff_Rij*ds*(matrix_Rx(i,j)*(k_ii*k_x_factor+eps) + &
                       matrix_Ry(i,j)*(k_jj*k_y_factor+eps) + &
                       matrix_Rz(i,j)*(k_kk*k_z_factor+eps))
# endif  /* Dim2 */
              G_A(i,j)=dcmplx(G_RA(i,j),G_IA(i,j))  

                enddo
               enddo
            t4=MPI_Wtime()
      else 
              G_RB=0.0
              G_IB=0.0
        do j=0,M_Bar_B-1
        !G_R has only the digonal element
         G_RB(j,j)=coeff_deltaij + ds*Loa*basis_SPH(j)%x*(basis_SPH(j)%x+1)/(2.0*kapa_temp)
                do i=0,M_Bar_B-1
# if defined (Dim2)
              G_IB(i,j)=coeff_Rij*ds*(matrix_Rx(i,j)*(k_ii*k_x_factor+eps) + &
                       matrix_Ry(i,j)*(k_jj*k_y_factor+eps))
# else  /* Dim2 */
              G_IB(i,j)=coeff_Rij*ds*(matrix_Rx(i,j)*(k_ii*k_x_factor+eps) + &
                       matrix_Ry(i,j)*(k_jj*k_y_factor+eps) + &
                       matrix_Rz(i,j)*(k_kk*k_z_factor+eps))
# endif  /* Dim2 */
              G_B(i,j)=dcmplx(G_RB(i,j),G_IB(i,j))  
                enddo
               enddo
       endif   
              

      if(G_index<=4) then
               call complex_matrix_inverse(G_RA,G_IA,M_Bar_A)
       else
              ! t1=MPI_Wtime()
               call complex_matrix_inverse(G_RB,G_IB,M_Bar_B)
              ! call zmat_inv(M_Bar_B, G_B)
              ! t2=MPI_Wtime()
      endif

                t21=t21+t2-t1   

                  k=k+1 

                 if((k_I*(SIDEy)+k_j+1)>local_size .or. & 
                  (k_I*(SIDEy)+k_j+1)<1) then
                  write(*,*) "k index exceeds",k_I*(SIDEy)+k_j+1
                  stop
                 endif

        job(1)=0
        job(2)=1
        job(3)=1
! fill in with upper triangular part of matrix 
        job(4)=1

        if(G_index<=4) then
        job(5)=namaxA
               do j=0,M_Bar_A-1
                  do i=0,M_Bar_A-1    
                   if(i>j) then
                   G_A(i,j)=dcmplx(0.0d0,0.0d0)
                   endif
                   enddo
               enddo
        else 
        job(5)=namaxB
               do j=0,M_Bar_B-1
                  do i=0,M_Bar_B-1    
                   if(i>j) then
                   G_B(i,j)=dcmplx(0.0d0,0.0d0)
                   endif
                   enddo
               enddo
        endif
        job(6)=1

        if(G_index<=4) then
       call mkl_zdnscsr(job,M_Bar_A,M_Bar_A,G_A,M_Bar_A,sy_values_A,sy_columns_A, &
        sy_rowIndex_A,info)
       !call sparse_dns_to_csr(M_Bar_A, M_Bar_A, namaxA, G_A, sy_values_A, &
       !                          sy_columns_A, sy_rowIndex_A)
      ! call get_symmetric_matrix_size_complex_dp(G_A,M_Bar_A,M_Bar_A,j)
      ! write(*,*) "j=",j
        else
       call mkl_zdnscsr(job,M_Bar_B,M_Bar_B,G_B,M_Bar_B,sy_values_B,sy_columns_B, &
        sy_rowIndex_B,info)
        endif
           
                 if(G_index==3) then
                 GA3_inverse_val(:,k)=sy_values_A(:)
                 GA3_inverse_col(:,k)=sy_columns_A(:)
                 GA3_inverse_row(:,k)=sy_rowIndex_A(:)

                 else if(G_index==4) then
                 GA3_star_inverse_val(:,k)=sy_values_A(:)
                 GA3_star_inverse_col(:,k)=sy_columns_A(:)
                 GA3_star_inverse_row(:,k)=sy_rowIndex_A(:)

                 else if(G_index==5) then
                 GB3_inverse_val(:,k)=sy_values_B(:)
                 GB3_inverse_col(:,k)=sy_columns_B(:)
                 GB3_inverse_row(:,k)=sy_rowIndex_B(:)

                 else if(G_index==6) then
                 GB3_star_inverse_val(:,k)=sy_values_B(:)
                 GB3_star_inverse_col(:,k)=sy_columns_B(:)
                 GB3_star_inverse_row(:,k)=sy_rowIndex_B(:)

                 !else if(G_index==5) then
                 !GB3_inverse_val(:,k)=sy_values_B(:)
                 !GB3_inverse_col(:,k)=sy_columns_B(:)
                 !GB3_inverse_row(:,k)=sy_rowIndex_B(:)

                 !else if(G_index==6) then
                 !GB3_star_inverse_val(:,k)=sy_values_B(:)
                 !GB3_star_inverse_col(:,k)=sy_columns_B(:)
                 !GB3_star_inverse_row(:,k)=sy_rowIndex_B(:)

                 !else if(G_index==7) then
                 !GB2_star_inverse_val(:,k)=sy_values_B(:)
                 !GB2_star_inverse_col(:,k)=sy_columns_B(:)
                 !GB2_star_inverse_row(:,k)=sy_rowIndex_B(:)

                 !else if(G_index==8) then
                 !GB1_star_inverse_val(:,k)=sy_values_B(:)
                 !GB1_star_inverse_col(:,k)=sy_columns_B(:)
                 !GB1_star_inverse_row(:,k)=sy_rowIndex_B(:)
                 endif
                 
!                 
!                      
!            !write(*,*) "namax_finalR1",namax_final_R1(k)           
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

         write(*,*) "ava inv time",t21/local_size
         write(*,*) "init matrix time",t4-t3

enddo  !enddo G_index
   
     deallocate(G_RA)
     deallocate(G_IA)
     deallocate(G_A)
     deallocate(sy_values_A)
     deallocate(sy_columns_A)
     deallocate(sy_rowIndex_A)
     deallocate(G_RB)
     deallocate(G_IB)
     deallocate(G_B)
     deallocate(sy_values_B)
     deallocate(sy_columns_B)
     deallocate(sy_rowIndex_B)
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
 integer,intent(in) :: K_i,K_j,k_k
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

!k counter
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


subroutine init_inverse_sparse_GAB2()
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
 REAL(DP),allocatable :: G_RA(:,:)
 REAL(DP),allocatable :: G_IA(:,:)
 REAL(DP),allocatable :: G_RB(:,:)
 REAL(DP),allocatable :: G_IB(:,:)
 DOUBLE PRECISION,allocatable :: sy_values_RA(:)
 DOUBLE PRECISION,allocatable :: sy_values_IA(:)
 DOUBLE PRECISION,allocatable :: sy_values_RB(:)
 DOUBLE PRECISION,allocatable :: sy_values_IB(:)
 integer,allocatable :: sy_columns_RA(:),sy_rowIndex_RA(:)
 integer,allocatable :: sy_columns_IA(:),sy_rowIndex_IA(:)
 integer,allocatable :: sy_columns_RB(:),sy_rowIndex_RB(:)
 integer,allocatable :: sy_columns_IB(:),sy_rowIndex_IB(:)
 REAL(DP) ::k_x_factor,k_y_factor,k_z_factor
# if defined (Dim2)
 integer,parameter :: namaxA=900
 integer,parameter :: namaxB=256
# else  /* Dim2 */
 integer,parameter :: namaxA=1800
 integer,parameter :: namaxB=256
# endif /* Dim2 */
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
!max_M_Bar=max(M_Bar_B,M_Bar_A)

allocate(G_RA(0:M_Bar_A-1,0:M_Bar_A-1))
allocate(G_IA(0:M_Bar_A-1,0:M_Bar_A-1))
allocate(G_RB(0:M_Bar_B-1,0:M_Bar_B-1))
allocate(G_IB(0:M_Bar_B-1,0:M_Bar_B-1))
          allocate(sy_values_RA(1:namaxA))
          allocate(sy_columns_RA(1:namaxA))
          allocate(sy_rowIndex_RA(1:M_Bar_A+1))

          allocate(sy_values_IA(1:namaxA))
          allocate(sy_columns_IA(1:namaxA))
          allocate(sy_rowIndex_IA(1:M_Bar_A+1))

          allocate(sy_values_RB(1:namaxB))
          allocate(sy_columns_RB(1:namaxB))
          allocate(sy_rowIndex_RB(1:M_Bar_B+1))

          allocate(sy_values_IB(1:namaxB))
          allocate(sy_columns_IB(1:namaxB))
          allocate(sy_rowIndex_IB(1:M_Bar_B+1))



!REAL
  allocate(R_GA1_inverse_val(1:namaxA,1:LOCAL_SIZE))
  allocate(R_GA1_inverse_col(1:namaxA,1:LOCAL_SIZE))
  allocate(R_GA1_inverse_row(1:M_Bar_A+1,1:LOCAL_SIZE))

  allocate(R_GA2_inverse_val(1:namaxA,1:LOCAL_SIZE))
  allocate(R_GA2_inverse_col(1:namaxA,1:LOCAL_SIZE))
  allocate(R_GA2_inverse_row(1:M_Bar_A+1,1:LOCAL_SIZE))

  allocate(R_GA3_inverse_val(1:namaxA,1:LOCAL_SIZE))
  allocate(R_GA3_inverse_col(1:namaxA,1:LOCAL_SIZE))
  allocate(R_GA3_inverse_row(1:M_Bar_A+1,1:LOCAL_SIZE))

  allocate(R_GA3_star_inverse_val(1:namaxA,1:LOCAL_SIZE))
  allocate(R_GA3_star_inverse_col(1:namaxA,1:LOCAL_SIZE))
  allocate(R_GA3_star_inverse_row(1:M_Bar_A+1,1:LOCAL_SIZE))

  allocate(R_GB3_inverse_val(1:namaxB,1:LOCAL_SIZE))
  allocate(R_GB3_inverse_col(1:namaxB,1:LOCAL_SIZE))
  allocate(R_GB3_inverse_row(1:M_Bar_B+1,1:LOCAL_SIZE))

  allocate(R_GB3_star_inverse_val(1:namaxB,1:LOCAL_SIZE))
  allocate(R_GB3_star_inverse_col(1:namaxB,1:LOCAL_SIZE))
  allocate(R_GB3_star_inverse_row(1:M_Bar_B+1,1:LOCAL_SIZE))

  allocate(R_GB2_star_inverse_val(1:namaxB,1:LOCAL_SIZE))
  allocate(R_GB2_star_inverse_col(1:namaxB,1:LOCAL_SIZE))
  allocate(R_GB2_star_inverse_row(1:M_Bar_B+1,1:LOCAL_SIZE))

  allocate(R_GB1_star_inverse_val(1:namaxB,1:LOCAL_SIZE))
  allocate(R_GB1_star_inverse_col(1:namaxB,1:LOCAL_SIZE))
  allocate(R_GB1_star_inverse_row(1:M_Bar_B+1,1:LOCAL_SIZE))
!!!IMAGE
  allocate(I_GA1_inverse_val(1:namaxA,1:LOCAL_SIZE))
  allocate(I_GA1_inverse_col(1:namaxA,1:LOCAL_SIZE))
  allocate(I_GA1_inverse_row(1:max_M_Bar+1,1:LOCAL_SIZE))

  allocate(I_GA2_inverse_val(1:namaxA,1:LOCAL_SIZE))
  allocate(I_GA2_inverse_col(1:namaxA,1:LOCAL_SIZE))
  allocate(I_GA2_inverse_row(1:max_M_Bar+1,1:LOCAL_SIZE))

  allocate(I_GA3_inverse_val(1:namaxA,1:LOCAL_SIZE))
  allocate(I_GA3_inverse_col(1:namaxA,1:LOCAL_SIZE))
  allocate(I_GA3_inverse_row(1:M_Bar_A+1,1:LOCAL_SIZE))

  allocate(I_GA3_star_inverse_val(1:namaxA,1:LOCAL_SIZE))
  allocate(I_GA3_star_inverse_col(1:namaxA,1:LOCAL_SIZE))
  allocate(I_GA3_star_inverse_row(1:M_Bar_A+1,1:LOCAL_SIZE))

  allocate(I_GB3_inverse_val(1:namaxB,1:LOCAL_SIZE))
  allocate(I_GB3_inverse_col(1:namaxB,1:LOCAL_SIZE))
  allocate(I_GB3_inverse_row(1:M_Bar_B+1,1:LOCAL_SIZE))

  allocate(I_GB3_star_inverse_val(1:namaxB,1:LOCAL_SIZE))
  allocate(I_GB3_star_inverse_col(1:namaxB,1:LOCAL_SIZE))
  allocate(I_GB3_star_inverse_row(1:M_Bar_B+1,1:LOCAL_SIZE))

  allocate(I_GB2_star_inverse_val(1:namaxB,1:LOCAL_SIZE))
  allocate(I_GB2_star_inverse_col(1:namaxB,1:LOCAL_SIZE))
  allocate(I_GB2_star_inverse_row(1:M_Bar_B+1,1:LOCAL_SIZE))

  allocate(I_GB1_star_inverse_val(1:namaxB,1:LOCAL_SIZE))
  allocate(I_GB1_star_inverse_col(1:namaxB,1:LOCAL_SIZE))
  allocate(I_GB1_star_inverse_row(1:M_Bar_B+1,1:LOCAL_SIZE))


!note that the array is storaged in cloumn before row in fortran
!so ,it is different from the c++ version code.
!always remember,otherwise it would be terribly wrong!!!!!!!!!!!!!!!!!!



# if defined (Dim2)
k_x_factor=2.0*PI/(SIDEx*dx)
k_y_factor=2.0*PI/(SIDEy*dy)
# else  /* Dim2 */
k_x_factor=2.0*PI/(SIDEx*dx)
k_y_factor=2.0*PI/(SIDEy*dy)
k_z_factor=2.0*PI/(SIDEz*dz)
# endif


do G_index=3,6

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



G_RA=0.0
G_IA=0.0
G_RB=0.0
G_IB=0.0


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
     if(G_index<=4) then
              G_RA=0.0
              G_IA=0.0
        do j=0,M_Bar_A-1
              !G_R has only the digonal element
         G_RA(j,j)=coeff_deltaij + ds*Loa*basis_SPH(j)%x*(basis_SPH(j)%x+1)/(2.0*kapa_temp)
                do i=0,M_Bar_A-1
# if defined (Dim2)
              G_IA(i,j)=coeff_Rij*ds*(matrix_Rx(i,j)*(k_ii*k_x_factor) + &
                       matrix_Ry(i,j)*(k_jj*k_y_factor))
# else  /* Dim2 */
              G_IA(i,j)=coeff_Rij*ds*(matrix_Rx(i,j)*(k_ii*k_x_factor) + &
                       matrix_Ry(i,j)*(k_jj*k_y_factor) + &
                       matrix_Rz(i,j)*(k_kk*k_z_factor))
# endif  /* Dim2 */

                enddo
               enddo
      else 
              G_RB=0.0
              G_IB=0.0
        do j=0,M_Bar_B-1
              !G_R has only the digonal element
         G_RB(j,j)=coeff_deltaij + ds*Loa*basis_SPH(j)%x*(basis_SPH(j)%x+1)/(2.0*kapa_temp)
                do i=0,M_Bar_B-1
# if defined (Dim2)
              G_IB(i,j)=coeff_Rij*ds*(matrix_Rx(i,j)*(k_ii*k_x_factor) + &
                       matrix_Ry(i,j)*(k_jj*k_y_factor))
# else  /* Dim2 */
              G_IB(i,j)=coeff_Rij*ds*(matrix_Rx(i,j)*(k_ii*k_x_factor) + &
                       matrix_Ry(i,j)*(k_jj*k_y_factor) + &
                       matrix_Rz(i,j)*(k_kk*k_z_factor))
# endif  /* Dim2 */
                enddo
               enddo
       endif   


      if(G_index<=4) then
              call complex_matrix_inverse(G_RA,G_IA,M_Bar_A)
       else
              call complex_matrix_inverse(G_RB,G_IB,M_Bar_B)
      endif
                
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

        if(G_index<=4) then
        job(5)=namaxA
               do j=0,M_Bar_A-1
                  do i=0,M_Bar_A-1    
                   if(j>i) then
                   G_RA(i,j)=0.0
                   G_IA(i,j)=0.0
                   endif
                   enddo
               enddo
        else 
        job(5)=namaxB
               do j=0,M_Bar_B-1
                  do i=0,M_Bar_B-1    
                   if(j>i) then
                   G_RB(i,j)=0.0
                   G_IB(i,j)=0.0
                   endif
                   enddo
               enddo
        endif
        job(6)=1


   
        if(G_index<=4) then
       call mkl_ddnscsr(job,M_Bar_A,M_Bar_A,G_RA,M_Bar_A,sy_values_RA,sy_columns_RA, &
        sy_rowIndex_RA,info)
       call mkl_ddnscsr(job,M_Bar_A,M_Bar_A,G_IA,M_Bar_A,sy_values_IA,sy_columns_IA, &
        sy_rowIndex_IA,info)
        else
       call mkl_ddnscsr(job,M_Bar_B,M_Bar_B,G_RB,M_Bar_B,sy_values_RB,sy_columns_RB, &
        sy_rowIndex_RB,info)
       call mkl_ddnscsr(job,M_Bar_B,M_Bar_B,G_IB,M_Bar_B,sy_values_IB,sy_columns_IB, &
        sy_rowIndex_IB,info)
        endif
           
 !       call convert_symmetric_matrix_to_array_cmplx(G_z,max_M_Bar,max_M_Bar,N_low,sy_values_z,  &
 !       sy_columns_z,sy_rowIndex_z) !this routine has an issue when N_low doesn't match the actual nonzero
 !       number of the matrix.Needs to be fixed.

                 if(G_index==1) then
                 R_GA1_inverse_val(:,k)=sy_values_RA(:)
                 R_GA1_inverse_col(:,k)=sy_columns_RA(:)
                 R_GA1_inverse_row(:,k)=sy_rowIndex_RA(:)

                 I_GA1_inverse_val(:,k)=sy_values_IA(:)
                 I_GA1_inverse_col(:,k)=sy_columns_IA(:)
                 I_GA1_inverse_row(:,k)=sy_rowIndex_IA(:)
                 else if(G_index==2) then
                 R_GA2_inverse_val(:,k)=sy_values_RA(:)
                 R_GA2_inverse_col(:,k)=sy_columns_RA(:)
                 R_GA2_inverse_row(:,k)=sy_rowIndex_RA(:)

                 I_GA2_inverse_val(:,k)=sy_values_IA(:)
                 I_GA2_inverse_col(:,k)=sy_columns_IA(:)
                 I_GA2_inverse_row(:,k)=sy_rowIndex_IA(:)

                 else if(G_index==3) then
                 R_GA3_inverse_val(:,k)=sy_values_RA(:)
                 R_GA3_inverse_col(:,k)=sy_columns_RA(:)
                 R_GA3_inverse_row(:,k)=sy_rowIndex_RA(:)

                 I_GA3_inverse_val(:,k)=sy_values_IA(:)
                 I_GA3_inverse_col(:,k)=sy_columns_IA(:)
                 I_GA3_inverse_row(:,k)=sy_rowIndex_IA(:)
                 else if(G_index==4) then
                 R_GA3_star_inverse_val(:,k)=sy_values_RA(:)
                 R_GA3_star_inverse_col(:,k)=sy_columns_RA(:)
                 R_GA3_star_inverse_row(:,k)=sy_rowIndex_RA(:)

                 I_GA3_star_inverse_val(:,k)=sy_values_IA(:)
                 I_GA3_star_inverse_col(:,k)=sy_columns_IA(:)
                 I_GA3_star_inverse_row(:,k)=sy_rowIndex_IA(:)

                 else if(G_index==5) then
                 R_GB3_inverse_val(:,k)=sy_values_RB(:)
                 R_GB3_inverse_col(:,k)=sy_columns_RB(:)
                 R_GB3_inverse_row(:,k)=sy_rowIndex_RB(:)

                 I_GB3_inverse_val(:,k)=sy_values_IB(:)
                 I_GB3_inverse_col(:,k)=sy_columns_IB(:)
                 I_GB3_inverse_row(:,k)=sy_rowIndex_IB(:)
                 else if(G_index==6) then
                 R_GB3_star_inverse_val(:,k)=sy_values_RB(:)
                 R_GB3_star_inverse_col(:,k)=sy_columns_RB(:)
                 R_GB3_star_inverse_row(:,k)=sy_rowIndex_RB(:)

                 I_GB3_star_inverse_val(:,k)=sy_values_IB(:)
                 I_GB3_star_inverse_col(:,k)=sy_columns_IB(:)
                 I_GB3_star_inverse_row(:,k)=sy_rowIndex_IB(:)
                 else if(G_index==7) then
                 R_GB2_star_inverse_val(:,k)=sy_values_RB(:)
                 R_GB2_star_inverse_col(:,k)=sy_columns_RB(:)
                 R_GB2_star_inverse_row(:,k)=sy_rowIndex_RB(:)

                 I_GB2_star_inverse_val(:,k)=sy_values_IB(:)
                 I_GB2_star_inverse_col(:,k)=sy_columns_IB(:)
                 I_GB2_star_inverse_row(:,k)=sy_rowIndex_IB(:)
                 else if(G_index==8) then
                 R_GB1_star_inverse_val(:,k)=sy_values_RB(:)
                 R_GB1_star_inverse_col(:,k)=sy_columns_RB(:)
                 R_GB1_star_inverse_row(:,k)=sy_rowIndex_RB(:)

                 I_GB1_star_inverse_val(:,k)=sy_values_IB(:)
                 I_GB1_star_inverse_col(:,k)=sy_columns_IB(:)
                 I_GB1_star_inverse_row(:,k)=sy_rowIndex_IB(:)
                 endif
                 
!                 
!                      
!            !write(*,*) "namax_finalR1",namax_final_R1(k)           
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
   
     deallocate(G_RA)
     deallocate(G_IA)
     deallocate(sy_values_RA)
     deallocate(sy_columns_RA)
     deallocate(sy_rowIndex_RA)
     deallocate(sy_values_IA)
     deallocate(sy_columns_IA)
     deallocate(sy_rowIndex_IA)
     deallocate(G_RB)
     deallocate(G_IB)
     deallocate(sy_values_RB)
     deallocate(sy_columns_RB)
     deallocate(sy_rowIndex_RB)
     deallocate(sy_values_IB)
     deallocate(sy_columns_IB)
     deallocate(sy_rowIndex_IB)
     call mp_barrier() 
      if(myid==0) then
        write(*,*) "init G matrix memeory not optimized"
       endif
      call mp_barrier()
    end subroutine init_inverse_sparse_GAB2


end module  G_matrix_mod
# endif /* SLAB */

