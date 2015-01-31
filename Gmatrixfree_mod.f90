# if defined (SLAB)
# else /* SLAB */
Module Gmatrixfree_mod
implicit none

contains

!test the basic size imformation of matrix G
!# if defined  (Test)
subroutine init_sparse_LRij()
 USE nrtype,only :DP,PI
 USE global_para
 USE mpi
 USE control
 USE constants
 USE utility
 USE mmpi
 USE sparse_MKL
 implicit none
 integer :: i,j,k
 integer :: K_i,K_ii
 integer :: K_j,K_jj
 integer :: K_k,K_kk
 integer :: G_index
  REAL(DP) :: kapa_temp
 REAL(DP),allocatable :: G_R(:)
 REAL(DP),allocatable :: G_Ix(:,:)
 REAL(DP),allocatable :: G_Iy(:,:)
 REAL(DP),allocatable :: G_Iz(:,:)
 DOUBLE PRECISION,allocatable :: sy_values_Ix(:)
 integer,allocatable :: sy_columns_Ix(:),sy_rowIndex_Ix(:)
 DOUBLE PRECISION,allocatable :: sy_values_Iy(:)
 integer,allocatable :: sy_columns_Iy(:),sy_rowIndex_Iy(:)
 DOUBLE PRECISION,allocatable :: sy_values_Iz(:)
 integer,allocatable :: sy_columns_Iz(:),sy_rowIndex_Iz(:)
 REAL(DP) ::k_x_factor,k_y_factor,k_z_factor

 integer,parameter :: namaxA1=113
 integer,parameter :: namaxA2=113
 integer,parameter :: namaxA3=64
 !integer,parameter :: namaxB1=113
 !integer,parameter :: namaxB2=113
 !integer,parameter :: namaxB3=64
 integer,parameter :: namax=113
real(DP) :: coeff_deltaij
real(DP) :: coeff_Rij
integer :: istat2,M_Bar
integer :: N_low1,N_low2,N_low3
integer :: loc(1)
integer job(8)
integer info

! we do not need to have set up the G matrix for A and B both,just one will do.
max_M_Bar=max(M_Bar_B,M_Bar_A)
M_Bar=max(M_Bar_B,M_Bar_A)
allocate(G_R(0:max_M_Bar-1))
allocate(G_Ix(0:max_M_Bar-1,0:max_M_Bar-1))
allocate(G_Iy(0:max_M_Bar-1,0:max_M_Bar-1))
allocate(G_Iz(0:max_M_Bar-1,0:max_M_Bar-1))

allocate(sy_values_Ix(1:namaxA1))
allocate(sy_columns_Ix(1:namaxA1))
allocate(sy_rowIndex_Ix(1:1+M_Bar))

allocate(sy_values_Iy(1:namaxA2))
allocate(sy_columns_Iy(1:namaxA2))
allocate(sy_rowIndex_Iy(1:M_Bar+1))

allocate(sy_values_Iz(1:namaxA3))
allocate(sy_columns_Iz(1:namaxA3))
allocate(sy_rowIndex_Iz(1:1+M_Bar))

! allocate Matrix G for explicit method
allocate(R_G_val(0:max_M_Bar-1))

allocate(I_GA1x_val(1:namaxA1)) 
allocate(I_GA1x_col(1:namaxA1))
allocate(I_GA1x_row(1:1+M_Bar_A))
allocate(I_GA1y_val(1:namaxA2)) 
allocate(I_GA1y_col(1:namaxA2))
allocate(I_GA1y_row(1:1+M_Bar_A))
allocate(I_GA1z_val(1:namaxA3)) 
allocate(I_GA1z_col(1:namaxA3))
allocate(I_GA1z_row(1:1+M_Bar_A))

!allocate(I_GA2x_val(1:namaxA1)) 
!allocate(I_GA2x_col(1:namaxA1))
!allocate(I_GA2x_row(1:1+M_Bar_A))
!allocate(I_GA2y_val(1:namaxA2)) 
!allocate(I_GA2y_col(1:namaxA2))
!allocate(I_GA2y_row(1:1+M_Bar_A))
!allocate(I_GA2z_val(1:namaxA3)) 
!allocate(I_GA2z_col(1:namaxA3))
!allocate(I_GA2z_row(1:1+M_Bar_A))
!
!allocate(I_GA3x_val(1:namaxA1)) 
!allocate(I_GA3x_col(1:namaxA1))
!allocate(I_GA3x_row(1:1+M_Bar_A))
!allocate(I_GA3y_val(1:namaxA2)) 
!allocate(I_GA3y_col(1:namaxA2))
!allocate(I_GA3y_row(1:1+M_Bar_A))
!allocate(I_GA3z_val(1:namaxA3)) 
!allocate(I_GA3z_col(1:namaxA3))
!allocate(I_GA3z_row(1:1+M_Bar_A))
!
!allocate(I_GA3x_star_val(1:namaxA1)) 
!allocate(I_GA3x_star_col(1:namaxA1))
!allocate(I_GA3x_star_row(1:1+M_Bar_A))
!allocate(I_GA3y_star_val(1:namaxA2)) 
!allocate(I_GA3y_star_col(1:namaxA2))
!allocate(I_GA3y_star_row(1:1+M_Bar_A))
!allocate(I_GA3z_star_val(1:namaxA3)) 
!allocate(I_GA3z_star_col(1:namaxA3))
!allocate(I_GA3z_star_row(1:1+M_Bar_A))
!
!allocate(I_GB3x_star_val(1:namaxA1)) 
!allocate(I_GB3x_star_col(1:namaxA1))
!allocate(I_GB3x_star_row(1:1+M_Bar_A))
!allocate(I_GB3y_star_val(1:namaxA2)) 
!allocate(I_GB3y_star_col(1:namaxA2))
!allocate(I_GB3y_star_row(1:1+M_Bar_A))
!allocate(I_GB3z_star_val(1:namaxA3)) 
!allocate(I_GB3z_star_col(1:namaxA3))
!allocate(I_GB3z_star_row(1:1+M_Bar_A))
!
!allocate(I_GB3x_val(1:namaxA1)) 
!allocate(I_GB3x_col(1:namaxA1))
!allocate(I_GB3x_row(1:1+M_Bar_A))
!allocate(I_GB3y_val(1:namaxA2)) 
!allocate(I_GB3y_col(1:namaxA2))
!allocate(I_GB3y_row(1:1+M_Bar_A))
!allocate(I_GB3z_val(1:namaxA3)) 
!allocate(I_GB3z_col(1:namaxA3))
!allocate(I_GB3z_row(1:1+M_Bar_A))
!
!allocate(I_GB2x_star_val(1:namaxA1)) 
!allocate(I_GB2x_star_col(1:namaxA1))
!allocate(I_GB2x_star_row(1:1+M_Bar_A))
!allocate(I_GB2y_star_val(1:namaxA2)) 
!allocate(I_GB2y_star_col(1:namaxA2))
!allocate(I_GB2y_star_row(1:1+M_Bar_A))
!allocate(I_GB2z_star_val(1:namaxA3)) 
!allocate(I_GB2z_star_col(1:namaxA3))
!allocate(I_GB2z_star_row(1:1+M_Bar_A))
!
!allocate(I_GB1x_star_val(1:namaxA1)) 
!allocate(I_GB1x_star_col(1:namaxA1))
!allocate(I_GB1x_star_row(1:1+M_Bar_A))
!allocate(I_GB1y_star_val(1:namaxA2)) 
!allocate(I_GB1y_star_col(1:namaxA2))
!allocate(I_GB1y_star_row(1:1+M_Bar_A))
!allocate(I_GB1z_star_val(1:namaxA3)) 
!allocate(I_GB1z_star_col(1:namaxA3))
!allocate(I_GB1z_star_row(1:1+M_Bar_A))

!note that the array is storaged in cloumn before row in fortran
!so ,it is different from the c++ version code.
!always remember,otherwise it would be terribly wrong!!!!!!!!!!!!!!!!!!

do G_index=1,1

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
 
              !enter the i,j of spehrical harmonic basis index loop
              G_R=0.0
              do j=0,M_Bar-1
              !G_R has only the digonal element
              !when explicit method is used, the formular of matrix G is changed 
       !G_R(j)=coeff_deltaij - Loa*basis_SPH(j)%x* & 
       !         (basis_SPH(j)%x+1)/(2.0*kapa_temp)
       G_R(j)=  - Loa*basis_SPH(j)%x* & 
                (basis_SPH(j)%x+1)/(2.0d0)
                do i=0,M_Bar-1
       !G_Ix(i,j)=-1.0d0*coeff_Rij*matrix_Rx(i,j)
       !G_Iy(i,j)=-1.0d0*coeff_Rij*matrix_Ry(i,j)
       !G_Iz(i,j)=-1.0d0*coeff_Rij*matrix_Rz(i,j)
       G_Ix(i,j)=-1.0d0*matrix_Rx(i,j)
       G_Iy(i,j)=-1.0d0*matrix_Ry(i,j)
       G_Iz(i,j)=-1.0d0*matrix_Rz(i,j)
                enddo
               enddo

                  N_low1=0
    call get_symmetric_matrix_size_dp(G_Ix,M_Bar,M_Bar,N_low1)
           write(*,*) "N_low1=",N_low1,"on",myid

            job(5)=namaxA1
        job(1)=0
        job(2)=1
        job(3)=1
        job(4)=0
        job(6)=1
               do j=0,M_Bar-1
                  do i=0,M_Bar-1    
                   if(j>i) then
                   G_Ix(i,j)=0.0d0
                   endif
                   enddo
               enddo

   call mkl_ddnscsr(job,M_Bar,M_Bar,G_Ix,M_Bar,sy_values_Ix,sy_columns_Ix, &
        sy_rowIndex_Ix,info)

                  N_low1=0
    call get_symmetric_matrix_size_dp(G_Iy,M_Bar,M_Bar,N_low1)
           write(*,*) "N_low2=",N_low1
            job(5)=namaxA2
        job(1)=0
        job(2)=1
        job(3)=1
        job(4)=0
        job(6)=1
               do j=0,M_Bar-1
                  do i=0,M_Bar-1    
                   if(j>i) then
                   G_Iy(i,j)=0.0
                   endif
                   enddo
               enddo
      call mkl_ddnscsr(job,M_Bar,M_Bar,G_Iy,M_Bar,sy_values_Iy,sy_columns_Iy, &
        sy_rowIndex_Iy,info)

                  N_low1=0
    call get_symmetric_matrix_size_dp(G_Iz,M_Bar,M_Bar,N_low1)
           write(*,*) "N_low3=",N_low1
            job(5)=namaxA3
        job(1)=0
        job(2)=1
        job(3)=1
        job(4)=0
        job(6)=1
               do j=0,M_Bar-1
                  do i=0,M_Bar-1    
                   if(j>i) then
                   G_Iz(i,j)=0.0
                   endif
                   enddo
               enddo
      call mkl_ddnscsr(job,M_Bar,M_Bar,G_Iz,M_Bar,sy_values_Iz,sy_columns_Iz, &
        sy_rowIndex_Iz,info)

                 R_G_val(:)=G_R(:)
                 I_GA1x_val(:)=sy_values_Ix(:)
                 I_GA1x_col(:)=sy_columns_Ix(:)
                 I_GA1x_row(:)=sy_rowIndex_Ix(:)
                 I_GA1y_val(:)=sy_values_Iy(:)
                 I_GA1y_col(:)=sy_columns_Iy(:)
                 I_GA1y_row(:)=sy_rowIndex_Iy(:)
                 I_GA1z_val(:)=sy_values_Iz(:)
                 I_GA1z_col(:)=sy_columns_Iz(:)
                 I_GA1z_row(:)=sy_rowIndex_Iz(:)
                ! if(G_index==1) then
                ! I_GA1x_val(:)=sy_values_Ix(:)
                ! I_GA1x_col(:)=sy_columns_Ix(:)
                ! I_GA1x_row(:)=sy_rowIndex_Ix(:)
                ! I_GA1y_val(:)=sy_values_Iy(:)
                ! I_GA1y_col(:)=sy_columns_Iy(:)
                ! I_GA1y_row(:)=sy_rowIndex_Iy(:)
                ! I_GA1z_val(:)=sy_values_Iz(:)
                ! I_GA1z_col(:)=sy_columns_Iz(:)
                ! I_GA1z_row(:)=sy_rowIndex_Iz(:)
                ! else if(G_index==2) then
                ! I_GA2x_val(:)=  sy_values_Ix(:)
                ! I_GA2x_col(:)= sy_columns_Ix(:)
                ! I_GA2x_row(:)=sy_rowIndex_Ix(:)
                ! I_GA2y_val(:)=  sy_values_Iy(:)
                ! I_GA2y_col(:)= sy_columns_Iy(:)
                ! I_GA2y_row(:)=sy_rowIndex_Iy(:)
                ! I_GA2z_val(:)=  sy_values_Iz(:)
                ! I_GA2z_col(:)= sy_columns_Iz(:)
                ! I_GA2z_row(:)=sy_rowIndex_Iz(:)
                ! else if(G_index==3) then
                ! I_GA3x_val(:)=sy_values_Ix(:)
                ! I_GA3x_col(:)=sy_columns_Ix(:)
                ! I_GA3x_row(:)=sy_rowIndex_Ix(:)
                ! I_GA3y_val(:)=sy_values_Iy(:)
                ! I_GA3y_col(:)=sy_columns_Iy(:)
                ! I_GA3y_row(:)=sy_rowIndex_Iy(:)
                ! I_GA3z_val(:)=sy_values_Iz(:)
                ! I_GA3z_col(:)=sy_columns_Iz(:)
                ! I_GA3z_row(:)=sy_rowIndex_Iz(:)
                ! else if(G_index==4) then
                ! I_GA3x_star_val(:)=sy_values_Ix(:)
                ! I_GA3x_star_col(:)=sy_columns_Ix(:)
                ! I_GA3x_star_row(:)=sy_rowIndex_Ix(:)
                ! I_GA3y_star_val(:)=sy_values_Iy(:)
                ! I_GA3y_star_col(:)=sy_columns_Iy(:)
                ! I_GA3y_star_row(:)=sy_rowIndex_Iy(:)
                ! I_GA3z_star_val(:)=sy_values_Iz(:)
                ! I_GA3z_star_col(:)=sy_columns_Iz(:)
                ! I_GA3z_star_row(:)=sy_rowIndex_Iz(:)
                ! else if(G_index==5) then
                ! I_GB3x_val(:)=sy_values_Ix(:)
                ! I_GB3x_col(:)=sy_columns_Ix(:)
                ! I_GB3x_row(:)=sy_rowIndex_Ix(:)
                ! I_GB3y_val(:)=sy_values_Iy(:)
                ! I_GB3y_col(:)=sy_columns_Iy(:)
                ! I_GB3y_row(:)=sy_rowIndex_Iy(:)
                ! I_GB3z_val(:)=sy_values_Iz(:)
                ! I_GB3z_col(:)=sy_columns_Iz(:)
                ! I_GB3z_row(:)=sy_rowIndex_Iz(:)
                ! else if(G_index==6) then
                ! I_GB3x_star_val(:)=sy_values_Ix(:)
                ! I_GB3x_star_col(:)=sy_columns_Ix(:)
                ! I_GB3x_star_row(:)=sy_rowIndex_Ix(:)
                ! I_GB3y_star_val(:)=sy_values_Iy(:)
                ! I_GB3y_star_col(:)=sy_columns_Iy(:)
                ! I_GB3y_star_row(:)=sy_rowIndex_Iy(:)
                ! I_GB3z_star_val(:)=sy_values_Iz(:)
                ! I_GB3z_star_col(:)=sy_columns_Iz(:)
                ! I_GB3z_star_row(:)=sy_rowIndex_Iz(:)
                ! else if(G_index==7) then
                ! I_GB2x_star_val(:)=sy_values_Ix(:)
                ! I_GB2x_star_col(:)=sy_columns_Ix(:)
                ! I_GB2x_star_row(:)=sy_rowIndex_Ix(:)
                ! I_GB2y_star_val(:)=sy_values_Iy(:)
                ! I_GB2y_star_col(:)=sy_columns_Iy(:)
                ! I_GB2y_star_row(:)=sy_rowIndex_Iy(:)
                ! I_GB2z_star_val(:)=sy_values_Iz(:)
                ! I_GB2z_star_col(:)=sy_columns_Iz(:)
                ! I_GB2z_star_row(:)=sy_rowIndex_Iz(:)
                ! else 
                ! I_GB1x_star_val(:)=sy_values_Ix(:)
                ! I_GB1x_star_col(:)=sy_columns_Ix(:)
                ! I_GB1x_star_row(:)=sy_rowIndex_Ix(:)
                ! I_GB1y_star_val(:)=sy_values_Iy(:)
                ! I_GB1y_star_col(:)=sy_columns_Iy(:)
                ! I_GB1y_star_row(:)=sy_rowIndex_Iy(:)
                ! I_GB1z_star_val(:)=sy_values_Iz(:)
                ! I_GB1z_star_col(:)=sy_columns_Iz(:)
                ! I_GB1z_star_row(:)=sy_rowIndex_Iz(:)
                ! endif
                 
enddo  !enddo G_index
   
     deallocate(G_R)
     deallocate(G_Ix)
     deallocate(G_Iy)
     deallocate(G_Iz)
     call mp_barrier() 

end subroutine init_sparse_LRij
!# endif /* Test */    


end module  Gmatrixfree_mod
# endif /* SLAB */












