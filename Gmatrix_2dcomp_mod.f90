# if defined (SLAB)
# else /* SLAB */
Module G_matrix_2dcomp_mod
implicit none

contains

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
 use decomp_fft_mod
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
max_M_Bar=max(M_Bar_B,M_Bar_A)
allocate(G_R(0:max_M_Bar-1,0:max_M_Bar-1))
allocate(G_I(0:max_M_Bar-1,0:max_M_Bar-1))
allocate(N_low_k1(1:LOCAL_SIZE_K))
allocate(N_low_k2(1:LOCAL_SIZE_K))

!note that the array is storaged in cloumn before row in fortran
!so ,it is different from the c++ version code.
!always remember,otherwise it would be terribly wrong!!!!!!!!!!!!!!!!!!

G_R=0.0
G_I=0.0

k_x_factor=2.0*PI/(SIDEx*dx)
k_y_factor=2.0*PI/(SIDEy*dy)
k_z_factor=2.0*PI/(SIDEz*dz)

do G_index=1,4

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


    write(*,*) "Max_M_Bar=",Max_M_bar

!k counter
k=0
do K_i=zstart(3),zend(3) 
  if(k_i-1<=(SiDEz/2)) then
   k_ii=k_i-1
  else 
   k_ii=(k_i-1)-SIDEz 
  endif
    do k_j=zstart(2),zend(2)
        if(K_j-1<=(SIDEy/2)) then 
          k_jj=K_j-1
        else
          K_jj=K_j-1-SIDEy
        endif

           do k_k=zstart(1),zend(1)
             if(K_k-1<=(SIDEx/2)) then
             k_kk=K_k-1
             else
             k_kk=K_k-1-SIDEx
             endif
              !enter the i,j of spehrical harmonic basis index loop
              G_R=0.0
              G_I=0.0
              do j=0,max_M_Bar-1
              !G_R has only the digonal element
       G_R(j,j)=coeff_deltaij + ds*Loa*basis_SPH(j)%x*(basis_SPH(j)%x+1)/(2.0*kapa_temp)
                do i=0,max_M_Bar-1
       G_I(i,j)=coeff_Rij*ds*(matrix_Rx(i,j)*(k_kk*k_x_factor) + &
                       matrix_Ry(i,j)*(k_jj*k_y_factor) + &
                       matrix_Rz(i,j)*(k_ii*k_z_factor))
                enddo
               enddo

              !inverse the complex matrix using the lapack lib directly
              call complex_matrix_inverse(G_R,G_I,max_M_Bar)
              
                  k=k+1 
                
                  N_low1=0
                  N_low2=0
                  call get_symmetric_matrix_size_dp(G_R,max_M_Bar,max_M_Bar,N_low1)
                  call get_symmetric_matrix_size_dp(G_I,max_M_Bar,max_M_Bar,N_low2)
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
        enddo  !enddo k_z
       enddo   !enddo k_y
      enddo    !enddo k_x
!
      if(k/=local_size_k) then
      write(*,*) "error,total k doesn't match local_size_k!!"
      write(*,*) "total k=",k
      stop
      endif

      if(mod(myid,4)==0) then
      write(*,*) " maxium number of N_lowk1 is",maxval(N_low_k1),"on",myid
      write(*,*) " loc maxium number of N_lowk1 is",maxloc(N_low_k1),"on",myid
      write(*,*) " minum number of N_lowk1 is",minval(N_low_k1),"on",myid
  !    write(*,*) " loc minum number of N_lowk1 is",minloc(N_low_k1),"on",myid
      write(*,*) " ava number of N_lowk1 is",sum(N_low_k1)/(1.0*Local_size_k),"on",myid

      write(*,*) " maxium number of N_lowk2 is",maxval(N_low_k2),"on",myid
 !     write(*,*) " loc maxium number of N_lowk2 is",maxloc(N_low_k2),"on",myid
      write(*,*) " minum number of N_lowk2 is",minval(N_low_k2),"on",myid
      write(*,*) " loc minum number of N_lowk2 is",minloc(N_low_k2),"on",myid
      write(*,*) " ava number of N_lowk2 is",sum(N_low_k2)/(1.0*Local_size_k),"on",myid
 
       endif
 

enddo  !enddo G_index
   
     deallocate(G_R)
     deallocate(G_I)
     deallocate(N_low_k1)
     deallocate(N_low_k2)
     call mp_barrier() 

    end subroutine test_inverse_sparse_GAB
!# endif /* Test */    


! the following subroutine will initialize the inverse of matrix G.
! when R2C FFT transform is performed, about half the k index can
! be sparsed because of the conjugate symmetry.However, one must handle
! the grid info in K space carefully.

subroutine init_inverse_sparse_GAB(R2C_FFT)
!!be careful,this subroutine is going to be called ever time 
!the box size is changed
 USE nrtype,only :DP,PI
 USE global_para
 USE mpi
 USE control
 USE constants
 USE utility
 USE mmpi
 use decomp_fft_mod
 USE matrix_inverse
 USE sparse_MKL
 implicit none
 logical,optional :: R2C_FFT
 integer :: LOCAL_SIZEK
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
# else  /* Dim2 */
 !integer,parameter :: namaxA=1800
 integer,parameter :: namaxA=1800
# endif /* Dim2 */
 integer,parameter :: namaxB=256
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
integer :: z_s,z_e,y_s,y_e,x_s,x_e

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

if(present(R2C_FFT)) then
   LOCAL_SIZEK=local_size_k_r2c
else
   LOCAL_SIZEK=local_size_k
endif


!REAL
  allocate(R_GA1_inverse_val(1:namaxA,1:LOCAL_SIZEK))
  allocate(R_GA1_inverse_col(1:namaxA,1:LOCAL_SIZEK))
  allocate(R_GA1_inverse_row(1:M_Bar_A+1,1:LOCAL_SIZEK))

  allocate(R_GA2_inverse_val(1:namaxA,1:LOCAL_SIZEK))
  allocate(R_GA2_inverse_col(1:namaxA,1:LOCAL_SIZEK))
  allocate(R_GA2_inverse_row(1:M_Bar_A+1,1:LOCAL_SIZEK))

  allocate(R_GA3_inverse_val(1:namaxA,1:LOCAL_SIZEK))
  allocate(R_GA3_inverse_col(1:namaxA,1:LOCAL_SIZEK))
  allocate(R_GA3_inverse_row(1:M_Bar_A+1,1:LOCAL_SIZEK))

  allocate(R_GA3_star_inverse_val(1:namaxA,1:LOCAL_SIZEK))
  allocate(R_GA3_star_inverse_col(1:namaxA,1:LOCAL_SIZEK))
  allocate(R_GA3_star_inverse_row(1:M_Bar_A+1,1:LOCAL_SIZEK))

  allocate(R_GB3_inverse_val(1:namaxB,1:LOCAL_SIZEK))
  allocate(R_GB3_inverse_col(1:namaxB,1:LOCAL_SIZEK))
  allocate(R_GB3_inverse_row(1:M_Bar_B+1,1:LOCAL_SIZEK))

  allocate(R_GB3_star_inverse_val(1:namaxB,1:LOCAL_SIZEK))
  allocate(R_GB3_star_inverse_col(1:namaxB,1:LOCAL_SIZEK))
  allocate(R_GB3_star_inverse_row(1:M_Bar_B+1,1:LOCAL_SIZEK))

  allocate(R_GB2_star_inverse_val(1:namaxB,1:LOCAL_SIZEK))
  allocate(R_GB2_star_inverse_col(1:namaxB,1:LOCAL_SIZEK))
  allocate(R_GB2_star_inverse_row(1:M_Bar_B+1,1:LOCAL_SIZEK))

  allocate(R_GB1_star_inverse_val(1:namaxB,1:LOCAL_SIZEK))
  allocate(R_GB1_star_inverse_col(1:namaxB,1:LOCAL_SIZEK))
  allocate(R_GB1_star_inverse_row(1:M_Bar_B+1,1:LOCAL_SIZEK))
!!!IMAGE
  allocate(I_GA1_inverse_val(1:namaxA,1:LOCAL_SIZEK))
  allocate(I_GA1_inverse_col(1:namaxA,1:LOCAL_SIZEK))
  allocate(I_GA1_inverse_row(1:max_M_Bar+1,1:LOCAL_SIZEK))

  allocate(I_GA2_inverse_val(1:namaxA,1:LOCAL_SIZEK))
  allocate(I_GA2_inverse_col(1:namaxA,1:LOCAL_SIZEK))
  allocate(I_GA2_inverse_row(1:max_M_Bar+1,1:LOCAL_SIZEK))

  allocate(I_GA3_inverse_val(1:namaxA,1:LOCAL_SIZEK))
  allocate(I_GA3_inverse_col(1:namaxA,1:LOCAL_SIZEK))
  allocate(I_GA3_inverse_row(1:M_Bar_A+1,1:LOCAL_SIZEK))

  allocate(I_GA3_star_inverse_val(1:namaxA,1:LOCAL_SIZEK))
  allocate(I_GA3_star_inverse_col(1:namaxA,1:LOCAL_SIZEK))
  allocate(I_GA3_star_inverse_row(1:M_Bar_A+1,1:LOCAL_SIZEK))

  allocate(I_GB3_inverse_val(1:namaxB,1:LOCAL_SIZEK))
  allocate(I_GB3_inverse_col(1:namaxB,1:LOCAL_SIZEK))
  allocate(I_GB3_inverse_row(1:M_Bar_B+1,1:LOCAL_SIZEK))

  allocate(I_GB3_star_inverse_val(1:namaxB,1:LOCAL_SIZEK))
  allocate(I_GB3_star_inverse_col(1:namaxB,1:LOCAL_SIZEK))
  allocate(I_GB3_star_inverse_row(1:M_Bar_B+1,1:LOCAL_SIZEK))

  allocate(I_GB2_star_inverse_val(1:namaxB,1:LOCAL_SIZEK))
  allocate(I_GB2_star_inverse_col(1:namaxB,1:LOCAL_SIZEK))
  allocate(I_GB2_star_inverse_row(1:M_Bar_B+1,1:LOCAL_SIZEK))

  allocate(I_GB1_star_inverse_val(1:namaxB,1:LOCAL_SIZEK))
  allocate(I_GB1_star_inverse_col(1:namaxB,1:LOCAL_SIZEK))
  allocate(I_GB1_star_inverse_row(1:M_Bar_B+1,1:LOCAL_SIZEK))


!note that the array is storaged in cloumn before row in fortran
!so ,it is different from the c++ version code.
!always remember,otherwise it would be terribly wrong!!!!!!!!!!!!!!!!!!



k_x_factor=2.0*PI/(SIDEx*dx)
k_y_factor=2.0*PI/(SIDEy*dy)
k_z_factor=2.0*PI/(SIDEz*dz)



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

              !enter the i,j of spehrical harmonic basis index loop
     if(G_index<=4) then
              G_RA=0.0
              G_IA=0.0
        do j=0,M_Bar_A-1
              !G_R has only the digonal element
         G_RA(j,j)=coeff_deltaij + ds*Loa*basis_SPH(j)%x*(basis_SPH(j)%x+1)/(2.0*kapa_temp)
                do i=0,M_Bar_A-1
              G_IA(i,j)=coeff_Rij*ds*(matrix_Rx(i,j)*(k_kk*k_x_factor) + &
                       matrix_Ry(i,j)*(k_jj*k_y_factor) + &
                       matrix_Rz(i,j)*(k_ii*k_z_factor))

                enddo
               enddo
      else 
              G_RB=0.0
              G_IB=0.0
        do j=0,M_Bar_B-1
              !G_R has only the digonal element
         G_RB(j,j)=coeff_deltaij + ds*Loa*basis_SPH(j)%x*(basis_SPH(j)%x+1)/(2.0*kapa_temp)
                do i=0,M_Bar_B-1
              G_IB(i,j)=coeff_Rij*ds*(matrix_Rx(i,j)*(k_kk*k_x_factor) + &
                       matrix_Ry(i,j)*(k_jj*k_y_factor) + &
                       matrix_Rz(i,j)*(k_ii*k_z_factor))
                enddo
               enddo
       endif   
      if(G_index<=4) then
              call complex_matrix_inverse(G_RA,G_IA,M_Bar_A)
       else
              call complex_matrix_inverse(G_RB,G_IB,M_Bar_B)
      endif
                  k=k+1 

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
                 
        enddo  !enddo k_z
       enddo   !enddo k_y
      enddo    !enddo k_x

      if(k/=local_sizek) then
      write(*,*) "error,total k doesn't match local_size_k!!"
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
    end subroutine init_inverse_sparse_GAB

end module  G_matrix_2dcomp_mod
# endif /* SLAB */













