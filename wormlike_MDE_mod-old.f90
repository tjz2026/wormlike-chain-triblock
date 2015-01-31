
subroutine calc_THETAij()
use global_para
 USE control
 USE constants
 USE utility
!use fftw
 USE mpi_fftw3_operation
use mmpi
implicit none
integer :: i,j,k
REAL*8 :: M11_M22,M33,M12,M13,M23,sumM,sumM_global
integer :: nonzero_counter,nonzero_tag

!!!for A
# if defined (OPENMP)
        !$omp parallel default(shared) private(j,i,M11_M22,M33,M12,M13,M23)
        !$omp do
# endif  /* OPENMP */
 do j=0,M_Bar_A-1
   do i=0,M_Bar_A-1
       M11_M22=THETAij_M11_M22_A(i,j)
       M33=THETAij_M33_A(i,j)
       M12=THETAij_M12_A(i,j)
       M13=THETAij_M13_A(i,j)
       M23=THETAij_M23_A(i,j)
      
       do k=0,LOCAL_SIZE-1
        THETAij_A(i,j,k)=(M_OP(k,0,0)-M_OP(k,1,1))*M11_M22 &
                         + M_OP(k,2,2)*M33 &
                         + M_OP(k,0,1)*M12 &
                         + M_Op(k,0,2)*M13 &
                         + M_OP(k,1,2)*M23
    enddo
   enddo
 enddo
# if defined (OPENMP)
        !$omp end  do
        !$omp end  parallel
# endif  /* OPENMP */




nonzero_counter=0
do k=0,Local_size-1
 do j=0,M_Bar_A-1
  do i=0,M_Bar_A-1
      
     if(abs(THETAij_A(i,j,k))>=Thresh_sprase_matrix) then
      nonzero_counter=nonzero_counter+1
     endif 
    enddo
   enddo
 enddo

# if defined (Debug)
 if (myid==0) then
  write(*,*) "nonzero_size in calc thetaij=",nonzero_counter
 endif
# endif /* Debug */
  nonzero_tag=nonzero_counter
 
if(nonzero_tag>0) then
 
  if(.not. allocated(THETA_nonzero_1D_A)) then
  allocate(THETA_nonzero_1D_A(1:nonzero_counter))
  else
  deallocate(THETA_nonzero_1D_A)
  allocate(THETA_nonzero_1D_A(1:nonzero_counter))
  endif
  
nonzero_counter=0

THETA_nonzero_1D_A_indexK(-1)=0
do k=0,Local_size-1
do i=0,M_Bar_A-1
  do j=0,M_Bar_A-1
      
     if(abs(THETAij_A(i,j,k))>=Thresh_sprase_matrix) then
      nonzero_counter=nonzero_counter+1
      THETA_nonzero_1D_A(nonzero_counter)%i=i
      THETA_nonzero_1D_A(nonzero_counter)%j=j
      THETA_nonzero_1D_A(nonzero_counter)%value=THETAij_A(i,j,k)
      
     endif 
    enddo
   enddo
      THETA_nonzero_1D_A_indexK(k)=nonzero_counter
     
 enddo

else if(nonzero_tag==0) then
 if(.not. allocated(THETA_nonzero_1D_A)) then
  allocate(THETA_nonzero_1D_A(1:LOCAL_SIZE))
  else
  deallocate(THETA_nonzero_1D_A)
  allocate(THETA_nonzero_1D_A(1:LOCAL_SIZE))
  endif

nonzero_counter=0

THETA_nonzero_1D_A_indexK(-1)=0
do k=0,Local_size-1

      nonzero_counter=nonzero_counter+1
      THETA_nonzero_1D_A(nonzero_counter)%i=0
      THETA_nonzero_1D_A(nonzero_counter)%j=0
      THETA_nonzero_1D_A(nonzero_counter)%value=0.0d0
      THETA_nonzero_1D_A_indexK(k)=k+1

 enddo

endif



!!!for B

# if defined (OPENMP)
        !$omp parallel default(shared) private(j,i,M11_M22,M33,M12,M13,M23)
        !$omp do
# endif  /* OPENMP */
do j=0,M_Bar_B-1
  do i=0,M_Bar_B-1
       M11_M22=THETAij_M11_M22_B(i,j)
       M33=THETAij_M33_B(i,j)
       M12=THETAij_M12_B(i,j)
       M13=THETAij_M13_B(i,j)
       M23=THETAij_M23_B(i,j)
      
     do k=0,Local_size-1
        THETAij_B(i,j,k)=(M_OP(k,0,0)-M_OP(k,1,1))*M11_M22 &
                         + M_OP(k,2,2)*M33 &
                         + M_OP(k,0,1)*M12 &
                         + M_Op(k,0,2)*M13 &
                         + M_OP(k,1,2)*M23
    enddo
   enddo
 enddo
# if defined (OPENMP)
        !$omp end  do
        !$omp end  parallel
# endif  /* OPENMP */


nonzero_counter=0
do k=0,Local_size-1
 do j=0,M_Bar_B-1
  do i=0,M_Bar_B-1
      
     if(abs(THETAij_B(i,j,k))>=Thresh_sprase_matrix) then
      nonzero_counter=nonzero_counter+1
     endif 
    enddo
   enddo
 enddo

nonzero_tag=nonzero_counter

if(nonzero_tag>0) then


  if(.not. allocated(THETA_nonzero_1D_B)) then
  allocate(THETA_nonzero_1D_B(1:nonzero_counter))
  else
  deallocate(THETA_nonzero_1D_B)
  allocate(THETA_nonzero_1D_B(1:nonzero_counter))
  endif


  
nonzero_counter=0

THETA_nonzero_1D_B_indexK(-1)=0
do k=0,Local_size-1
do i=0,M_Bar_B-1
  do j=0,M_Bar_B-1
      
     if(abs(THETAij_B(i,j,k))>=Thresh_sprase_matrix) then
      nonzero_counter=nonzero_counter+1
      THETA_nonzero_1D_B(nonzero_counter)%i=i
      THETA_nonzero_1D_B(nonzero_counter)%j=j
      THETA_nonzero_1D_B(nonzero_counter)%value=THETAij_B(i,j,k)
      
     endif 
    enddo
   enddo
      THETA_nonzero_1D_B_indexK(k)=nonzero_counter
 enddo



else if(nonzero_tag==0) then
 if(.not. allocated(THETA_nonzero_1D_B)) then
  allocate(THETA_nonzero_1D_B(1:LOCAL_SIZE))
  else
  deallocate(THETA_nonzero_1D_B)
  allocate(THETA_nonzero_1D_B(1:LOCAL_SIZE))
  endif



nonzero_counter=0


THETA_nonzero_1D_B_indexK(-1)=0
do k=0,Local_size-1

      nonzero_counter=nonzero_counter+1
      THETA_nonzero_1D_B(nonzero_counter)%i=0
      THETA_nonzero_1D_B(nonzero_counter)%j=0
      THETA_nonzero_1D_B(nonzero_counter)%value=0.0d0
      THETA_nonzero_1D_B_indexK(k)=k+1

 enddo


endif



end subroutine calc_THETAij

Module wormlike_MDE
implicit none
contains





# if defined (Mem_OP)

subroutine MDE_q()
 USE nrtype,only :DP,PI
 USE global_para
 USE mpi
 USE control
 USE constants
 USE utility
 USE mmpi
 !USE fftw !user defined module
 USE mpi_fftw3_operation
 USE SPH
 USE matrix_inverse
 USE G_matrix_mod
 implicit none
 integer :: i,j,s
 integer :: index_nonzero
 integer :: k
 integer :: istat2
 integer,parameter :: namax=1800 !this number is obtained by a prerun of input parameters
 !for A
 real(DP) :: tempA1,tempA2,tempA3,temp,tempx,Wtemp,sumaq,sumaq_global,sumaq1
 real(DP) :: sumPA,sumPA_global
 real(DP),allocatable :: Pr_A(:,:)
 real(DP),allocatable :: Pk_Real_A(:,:)
 real(DP),allocatable :: Pk_Imag_A(:,:)
 real(DP),allocatable :: q_temp_RR_A(:,:)
 real(DP),allocatable :: q_temp_RI_A(:,:)
 real(DP),allocatable :: q_temp_IR_A(:,:)
 real(DP),allocatable :: q_temp_II_A(:,:)
 !double complex,allocatable :: q_temp_cmplx_A(:,:)
 !for B
 real(DP),allocatable :: Pr_B(:,:)
 real(DP),allocatable :: Pk_Real_B(:,:)
 real(DP),allocatable :: Pk_Imag_B(:,:)
 real(DP),allocatable :: q_temp_RR_B(:,:)
 real(DP),allocatable :: q_temp_RI_B(:,:)
 real(DP),allocatable :: q_temp_IR_B(:,:)
 real(DP),allocatable :: q_temp_II_B(:,:)
! double complex,allocatable :: q_temp_cmplx_B(:,:)
 REAL(DP) :: M_grid_inv

 REAL(DP) :: sy_values_R(1:namax)
 REAL(DP) :: sy_values_I(1:namax)
 integer ::  sy_columns_R(1:namax)
 integer ::  sy_columns_I(1:namax)
!!M_Bar_A >=M_BAr_B
 integer ::  sy_rowIndex_R(1:M_Bar_A+1)
 integer ::  sy_rowIndex_I(1:M_Bar_A+1)


 !!!! just for reference
! allocate(qA(0:NA,0:MA,0:LS),stat=istat2)

 !first ,allocate the matrixs
!for A
 allocate(Pr_A(0:M_Bar_A-1,0:LOCAL_SIZE-1),stat=istat2)
 allocate(Pk_Real_A(0:M_Bar_A-1,0:LOCAL_SIZE-1),stat=istat2)
 allocate(Pk_Imag_A(0:M_Bar_A-1,0:LOCAL_SIZE-1),stat=istat2)


 allocate(q_temp_RR_A(0:M_Bar_A-1,0:LOCAL_SIZE-1),stat=istat2)
 allocate(q_temp_RI_A(0:M_Bar_A-1,0:LOCAL_SIZE-1),stat=istat2)
 allocate(q_temp_IR_A(0:M_Bar_A-1,0:LOCAL_SIZE-1),stat=istat2)
 allocate(q_temp_II_A(0:M_Bar_A-1,0:LOCAL_SIZE-1),stat=istat2)
!for B
 allocate(Pr_B(0:M_Bar_B-1,0:LOCAL_SIZE-1),stat=istat2)
 allocate(Pk_Real_B(0:M_Bar_B-1,0:LOCAL_SIZE-1),stat=istat2)
 allocate(Pk_Imag_B(0:M_Bar_B-1,0:LOCAL_SIZE-1),stat=istat2)

 allocate(q_temp_RR_B(0:M_Bar_B-1,0:LOCAL_SIZE-1),stat=istat2)
 allocate(q_temp_RI_B(0:M_Bar_B-1,0:LOCAL_SIZE-1),stat=istat2)
 allocate(q_temp_IR_B(0:M_Bar_B-1,0:LOCAL_SIZE-1),stat=istat2)
 allocate(q_temp_II_B(0:M_Bar_B-1,0:LOCAL_SIZE-1),stat=istat2)
  

 if(istat2/=0) then
 write(*,*) "allocate matrix in mde_q failed"
 stop
 endif

 M_grid_inv=1.0d0/M_grid


do k=0,LOCAL_SIZE-1
    qA(0,0,K)=1.0
   do i=1,M_Bar_A-1
    qA(0,i,K)=0.0
   enddo
enddo

       
 !!!!!
 do s=1,Nmax

 !!!s=1
   if(s==1) then

          Pr_A=0.0

# if defined (OPENMP)
        !$omp parallel default(shared) private(k,i,j,index_nonzero,Wtemp,tempx)
        !$omp do
# endif  /* OPENMP */
     
      do k=0,LOCAL_SIZE-1
          Wtemp=WA(k)
# if defined (MaierSaupe)
          do index_nonzero=THETA_nonzero_1D_A_indexK(k-1)+1, &
            THETA_nonzero_1D_A_indexK(k)
             i=THETA_nonzero_1D_A(index_nonzero)%i
             j=THETA_nonzero_1D_A(index_nonzero)%j
             tempx=THETA_nonzero_1D_A(index_nonzero)%value
             Pr_A(i,k)=Pr_A(i,k)+ds*tempx*qA(0,j,k)
           enddo
# endif /* MaierSaupe */
            do i=0,M_Bar_A-1
            Pr_A(i,k)=Pr_A(i,k)+(1.0-ds*Wtemp)*qA(0,i,k)
            enddo
        enddo !enddo k=0,LOCAL_SIZE    
# if defined (OPENMP)
        !$omp end  do
        !$omp end  parallel
# endif  /* OPENMP */





!!!doing the FFT

      do i=0,M_Bar_A-1
        do k=0,LOCAL_SIZE-1
          local_data(k+1)=cmplx(Pr_A(i,k),0.0)
         enddo
       
       
       call fftw_mpi_execute_dft(fplan, local_data, local_data)
       !call fftwnd_f77_mpi(plancf,1,local_data,work,0,FFTW_NORMAL_ORDER)
       
        do k=0,LOCAL_SIZE-1
            Pk_Real_A(i,k)=real(local_data(k+1))
            Pk_Imag_A(i,k)=AIMAG(local_data(k+1))
            !Pk_cmplx_A(i,k)=local_data(k)
         enddo
       enddo  !enddo i=0,M_Bar_A-1  


# if defined (OPENMP)
        !$omp parallel default(shared) private(k)
        !$omp do
# endif  /* OPENMP */
         do k=0,LOCAL_SIZE-1

          call  init_GA_inverse(1,k,sy_values_R, &
                sy_values_I,sy_columns_R,sy_columns_I, &
                sy_rowIndex_R,sy_rowIndex_I )
            

         call mkl_dcsrsymv('l', M_Bar_A,  sy_values_R, sy_rowIndex_R, &
            sy_columns_R,Pk_Real_A(:,k), q_temp_RR_A(:,k))


         call mkl_dcsrsymv('l', M_Bar_A,  sy_values_R, sy_rowIndex_R, &
            sy_columns_R,Pk_Imag_A(:,k), q_temp_RI_A(:,k))

         call mkl_dcsrsymv('l', M_Bar_A,  sy_values_I, sy_rowIndex_I, &
            sy_columns_I,Pk_Real_A(:,k), q_temp_IR_A(:,k))

         call mkl_dcsrsymv('l', M_Bar_A,  sy_values_I, sy_rowIndex_I, &
            sy_columns_I,Pk_Imag_A(:,k), q_temp_II_A(:,k))

            

         enddo


# if defined (OPENMP)
        !$omp end  do
        !$omp end  parallel
# endif  /* OPENMP */

             do i=0,M_Bar_A-1
                do k=0,LOCAL_SIZE-1
                local_data(k+1)=cmplx(q_temp_RR_A(i,k)-q_temp_II_A(i,k), & 
                q_temp_RI_A(i,k)+q_temp_IR_A(i,k))
                enddo
       call fftw_mpi_execute_dft(bplan, local_data, local_data)
            
                do k=0,LOCAL_SIZE-1
                 qA(1,i,k)=real(local_data(k+1))*M_grid_inv
                enddo
         
              enddo  


!!!!s=2!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    else if(s==2) then

          Pr_A=0.0
# if defined (OPENMP)
        !$omp parallel default(shared) private(k,i,j,index_nonzero,Wtemp,tempx)
        !$omp do
# endif  /* OPENMP */
      do k=0,LOCAL_SIZE-1
          Wtemp=WA(k)
# if defined (MaierSaupe)
          do index_nonzero=THETA_nonzero_1D_A_indexK(k-1)+1, &
            THETA_nonzero_1D_A_indexK(k)
             i=THETA_nonzero_1D_A(index_nonzero)%i
             j=THETA_nonzero_1D_A(index_nonzero)%j
             tempx=THETA_nonzero_1D_A(index_nonzero)%value
             Pr_A(i,k)=Pr_A(i,k)+ds*tempx*(2.0*qA(1,j,k)-qA(0,j,k))
           enddo
# endif /* MaierSaupe */
            do i=0,M_Bar_A-1
            Pr_A(i,k)=Pr_A(i,k)+2.0*qA(1,i,k)-0.5*qA(0,i,k) - &
                      ds*Wtemp*(2.0*qA(1,i,k)-qA(0,i,k))
            enddo
        enddo !enddo k=0,LOCAL_SIZE    
# if defined (OPENMP)
        !$omp end  do
        !$omp end  parallel
# endif  /* OPENMP */

!!!doing the FFT

      do i=0,M_Bar_A-1
        do k=0,LOCAL_SIZE-1
          local_data(k+1)=cmplx(Pr_A(i,k),0.0)
         enddo
       
       call fftw_mpi_execute_dft(fplan, local_data, local_data)
       
        do k=0,LOCAL_SIZE-1
            Pk_Real_A(i,k)=real(local_data(k+1))
            Pk_Imag_A(i,k)=AIMAG(local_data(k+1))
         enddo
       enddo  !enddo i=0,M_Bar_A-1  


# if defined (OPENMP)
        !$omp parallel default(shared) private(k)
        !$omp do
# endif  /* OPENMP */
         do k=0,LOCAL_SIZE-1
          
          call  init_GA_inverse(2,k,sy_values_R, &
                sy_values_I,sy_columns_R,sy_columns_I, &
                sy_rowIndex_R,sy_rowIndex_I )
 
         call mkl_dcsrsymv('l', M_Bar_A,  sy_values_R, sy_rowIndex_R, &
            sy_columns_R,Pk_Real_A(:,k), q_temp_RR_A(:,k))

         call mkl_dcsrsymv('l', M_Bar_A,  sy_values_R, sy_rowIndex_R, &
            sy_columns_R,Pk_Imag_A(:,k), q_temp_RI_A(:,k))

         call mkl_dcsrsymv('l', M_Bar_A,  sy_values_I, sy_rowIndex_I, &
            sy_columns_I,Pk_Real_A(:,k), q_temp_IR_A(:,k))

         call mkl_dcsrsymv('l', M_Bar_A,  sy_values_I, sy_rowIndex_I, &
            sy_columns_I,Pk_Imag_A(:,k), q_temp_II_A(:,k))

         enddo
# if defined (OPENMP)
        !$omp end  do
        !$omp end  parallel
# endif  /* OPENMP */


             do i=0,M_Bar_A-1
                do k=0,LOCAL_SIZE-1
                local_data(k+1)=cmplx(q_temp_RR_A(i,k)-q_temp_II_A(i,k), & 
                q_temp_RI_A(i,k)+q_temp_IR_A(i,k))
                enddo 
       call fftw_mpi_execute_dft(bplan, local_data, local_data)
                do k=0,LOCAL_SIZE-1
                 qA(2,i,k)=real(local_data(k+1))*M_grid_inv
                enddo
              enddo  

!!!!s>=3!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    else if(s>=3) then
          temp=1.0/3.0
           
          if(s<=NA) then
            
          Pr_A=0.0
# if defined (OPENMP)
        !$omp parallel default(shared) private(k,i,j,index_nonzero,Wtemp,tempx)
        !$omp do
# endif  /* OPENMP */
      do k=0,LOCAL_SIZE-1
          Wtemp=WA(k)
# if defined (MaierSaupe)
          do index_nonzero=THETA_nonzero_1D_A_indexK(k-1)+1, &
            THETA_nonzero_1D_A_indexK(k)
             i=THETA_nonzero_1D_A(index_nonzero)%i
             j=THETA_nonzero_1D_A(index_nonzero)%j
             tempx=THETA_nonzero_1D_A(index_nonzero)%value
             Pr_A(i,k)=Pr_A(i,k)+ds*tempx*(3.0*qA(s-1,j,k)- &
                       3.0*qA(s-2,j,k)+qA(s-3,j,k))
            
           enddo
# endif /* MaierSaupe */
            do i=0,M_Bar_A-1
            Pr_A(i,k)=Pr_A(i,k)+3.0*qA(s-1,i,k)-1.5*qA(s-2,i,k)+ &
                      temp*qA(s-3,i,k) - &
                      ds*Wtemp*(3.0*qA(s-1,i,k)-3.0*qA(s-2,i,k)+qA(s-3,i,k))
            enddo
        enddo !enddo k=0,LOCAL_SIZE    
# if defined (OPENMP)
        !$omp end  do
        !$omp end  parallel
# endif  /* OPENMP */

!!!doing the FFT

      do i=0,M_Bar_A-1
        do k=0,LOCAL_SIZE-1
          local_data(k+1)=cmplx(Pr_A(i,k),0.0)
         enddo
       
       
       call fftw_mpi_execute_dft(fplan, local_data, local_data)
       
        do k=0,LOCAL_SIZE-1
            Pk_Real_A(i,k)=real(local_data(k+1))
            Pk_Imag_A(i,k)=AIMAG(local_data(k+1))
         enddo
       enddo  !enddo i=0,M_Bar_A-1  


# if defined (OPENMP)
        !$omp parallel default(shared) private(k)
        !$omp do
# endif  /* OPENMP */
         do k=0,LOCAL_SIZE-1

         call mkl_dcsrsymv('l', M_Bar_A,  R_GA3_inverse_val(:,k+1), R_GA3_inverse_row(:,k+1), &
            R_GA3_inverse_col(:,k+1),Pk_Real_A(:,k), q_temp_RR_A(:,k))

         call mkl_dcsrsymv('l', M_Bar_A,  R_GA3_inverse_val(:,k+1), R_GA3_inverse_row(:,k+1), &
            R_GA3_inverse_col(:,k+1),Pk_Imag_A(:,k), q_temp_RI_A(:,k))

         call mkl_dcsrsymv('l', M_Bar_A,  I_GA3_inverse_val(:,k+1), I_GA3_inverse_row(:,k+1), &
            I_GA3_inverse_col(:,k+1),Pk_Real_A(:,k), q_temp_IR_A(:,k))

         call mkl_dcsrsymv('l', M_Bar_A,  I_GA3_inverse_val(:,k+1), I_GA3_inverse_row(:,k+1), &
            I_GA3_inverse_col(:,k+1),Pk_Imag_A(:,k), q_temp_II_A(:,k))
         enddo
# if defined (OPENMP)
        !$omp end  do
        !$omp end  parallel
# endif  /* OPENMP */

             do i=0,M_Bar_A-1
                do k=0,LOCAL_SIZE-1
                local_data(k+1)=cmplx(q_temp_RR_A(i,k)-q_temp_II_A(i,k), & 
                q_temp_RI_A(i,k)+q_temp_IR_A(i,k))
                enddo 
       call fftw_mpi_execute_dft(bplan, local_data, local_data)
                do k=0,LOCAL_SIZE-1
                 qA(s,i,k)=real(local_data(k+1))*M_grid_inv
                enddo
              enddo  

!!!!s>=3.and s>NA!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

else  ! else s>=3:s>NA
            !!!calculate polymer block B
             Pr_B=0.0
             ! else s>=3:s>NA:s==NA+1

   if(s==NA+1) then
             
# if defined (OPENMP)
        !$omp parallel default(shared) private(k,i,j,index_nonzero,Wtemp,tempx)
        !$omp do
# endif  /* OPENMP */
      do k=0,LOCAL_SIZE-1
          Wtemp=WB(k)
# if defined (MaierSaupe)
!!!!!!!!!!question : how to joint the two blocks with different M_Bar_A,M_Bar_B???
!!!if M_bar_A<M_Bar_B,it seems OK,but if M_Bar_A>M_Bar_B,then,it remains a
!question
          do index_nonzero=THETA_nonzero_1D_B_indexK(k-1)+1, &
            THETA_nonzero_1D_B_indexK(k)
             i=THETA_nonzero_1D_B(index_nonzero)%i
             j=THETA_nonzero_1D_B(index_nonzero)%j
             tempx=THETA_nonzero_1D_B(index_nonzero)%value
             tempA1=0.0
             tempA2=0.0
             tempA3=0.0

             if(j<M_Bar_A) then
              tempA1=qA(NA,j,k)
              tempA2=qA(NA-1,j,k)
              tempA3=qA(NA-2,j,k)
              endif
              Pr_B(i,k)=Pr_B(i,k) + ds*tempx*(3.0*tempA1-3.0*tempA2+tempA3)
           enddo  !enddo index_nonzero
# endif /* MaierSaupe */

            do i=0,M_Bar_B-1
                  tempA1=0.0
                  tempA2=0.0
                  tempA3=0.0
                  if(i<M_Bar_A) then
                  tempA1=qA(NA,i,k)
                  tempA2=qA(NA-1,i,k)
                  tempA3=qA(NA-2,i,k)
                  endif
                 

            Pr_B(i,k)=Pr_B(i,k)+3.0*tempA1-1.5*tempA2+ &
                      temp*tempA3 - &
                      ds*Wtemp*(3.0*tempA1-3.0*tempA2+tempA3)
            enddo
          enddo !enddo k=0,LOCAL_SIZE    
# if defined (OPENMP)
        !$omp end  do
        !$omp end  parallel
# endif  /* OPENMP */

    else if(s==(NA+2)) then
         
# if defined (OPENMP)
        !$omp parallel default(shared) private(k,i,j,index_nonzero,Wtemp,tempx)
        !$omp do
# endif  /* OPENMP */
      do k=0,LOCAL_SIZE-1
          Wtemp=WB(k)
# if defined (MaierSaupe)
          do index_nonzero=THETA_nonzero_1D_B_indexK(k-1)+1, &
            THETA_nonzero_1D_B_indexK(k)
             i=THETA_nonzero_1D_B(index_nonzero)%i
             j=THETA_nonzero_1D_B(index_nonzero)%j
             tempx=THETA_nonzero_1D_B(index_nonzero)%value
             tempA2=0.0
             tempA3=0.0
             if(j<M_Bar_A) then
              tempA2=qA(NA,j,k)
              tempA3=qA(NA-1,j,k)
              endif
              Pr_B(i,k)=Pr_B(i,k) + ds*tempx*(3.0*qB(1,j,k)-3.0*tempA2+tempA3)
           enddo  !enddo index_nonzero
# endif /* MaierSaupe */

            do i=0,M_Bar_B-1
                  tempA2=0.0
                  tempA3=0.0
                  if(i<M_Bar_A) then
                  tempA2=qA(NA,i,k)
                  tempA3=qA(NA-1,i,k)
                  endif
                 

            Pr_B(i,k)=Pr_B(i,k)+3.0*qB(1,i,k)-1.5*tempA2+ &
                      temp*tempA3 - &
                      ds*Wtemp*(3.0*qB(1,i,k)-3.0*tempA2+tempA3)
            enddo
          enddo !enddo k=0,LOCAL_SIZE    
# if defined (OPENMP)
        !$omp end  do
        !$omp end  parallel
# endif  /* OPENMP */

     else if(s==(NA+3)) then
         
# if defined (OPENMP)
        !$omp parallel default(shared) private(k,i,j,index_nonzero,Wtemp,tempx)
        !$omp do
# endif  /* OPENMP */
      do k=0,LOCAL_SIZE-1
          Wtemp=WB(k)
# if defined (MaierSaupe)
          do index_nonzero=THETA_nonzero_1D_B_indexK(k-1)+1, &
            THETA_nonzero_1D_B_indexK(k)
             i=THETA_nonzero_1D_B(index_nonzero)%i
             j=THETA_nonzero_1D_B(index_nonzero)%j
             tempx=THETA_nonzero_1D_B(index_nonzero)%value
             tempA3=0.0
             if(j<M_Bar_A) then
              tempA3=qA(NA,j,k)
              endif
              Pr_B(i,k)=Pr_B(i,k) + ds*tempx*(3.0*qB(2,j,k)-3.0*qB(1,j,k)+tempA3)
           enddo  !enddo index_nonzero
# endif /* MaierSaupe */

            do i=0,M_Bar_B-1
                  tempA3=0.0
                  if(i<M_Bar_A) then
                  tempA3=qA(NA,i,k)
                  endif
                 

            Pr_B(i,k)=Pr_B(i,k)+3.0*qB(2,i,k)-1.5*qB(1,i,k)+ &
                      temp*tempA3 - &
                      ds*Wtemp*(3.0*qB(2,i,k)-3.0*qB(1,i,k)+tempA3)
            enddo
          enddo !enddo k=0,LOCAL_SIZE    
# if defined (OPENMP)
        !$omp end  do
        !$omp end  parallel
# endif  /* OPENMP */

     else !else s>=3,s>NA:s>NA+3
           
# if defined (OPENMP)
        !$omp parallel default(shared) private(k,i,j,index_nonzero,Wtemp,tempx)
        !$omp do
# endif  /* OPENMP */
      do k=0,LOCAL_SIZE-1
          Wtemp=WB(k)
# if defined (MaierSaupe)
          do index_nonzero=THETA_nonzero_1D_B_indexK(k-1)+1, &
            THETA_nonzero_1D_B_indexK(k)
             i=THETA_nonzero_1D_B(index_nonzero)%i
             j=THETA_nonzero_1D_B(index_nonzero)%j
             tempx=THETA_nonzero_1D_B(index_nonzero)%value
              Pr_B(i,k)=Pr_B(i,k) + ds*tempx*(3.0*qB(s-NA-1,j,k)-3.0*qB(s-NA-2,j,k)+qB(s-NA-3,j,k))
           enddo  !enddo index_nonzero
# endif /* MaierSaupe */
            do i=0,M_Bar_B-1

            Pr_B(i,k)=Pr_B(i,k)+3.0*qB(s-NA-1,i,k)-1.5*qB(s-NA-2,i,k)+ &
                      temp*qB(s-NA-3,i,k) - &
                      ds*Wtemp*(3.0*qB(s-NA-1,i,k)-3.0*qB(s-NA-2,i,k)+qB(s-NA-3,i,k))
            enddo
          enddo !enddo k=0,LOCAL_SIZE    
# if defined (OPENMP)
        !$omp end  do
        !$omp end  parallel
# endif  /* OPENMP */
     endif  !if(s==NA+1) 

!!!doing the FFT

      do i=0,M_Bar_B-1
        do k=0,LOCAL_SIZE-1
          local_data(k+1)=cmplx(Pr_B(i,k),0.0)
         enddo
       
       call fftw_mpi_execute_dft(fplan, local_data, local_data)
       
        do k=0,LOCAL_SIZE-1
            Pk_Real_B(i,k)=real(local_data(k+1))
            Pk_Imag_B(i,k)=AIMAG(local_data(k+1))
         enddo
       enddo  !enddo i=0,M_Bar_A-1  


# if defined (OPENMP)
        !$omp parallel default(shared) private(k)
        !$omp do
# endif  /* OPENMP */
         do k=0,LOCAL_SIZE-1
          
         call mkl_dcsrsymv('l', M_Bar_B,  R_GA3_inverse_val(:,k+1), R_GA3_inverse_row(:,k+1), &
            R_GA3_inverse_col(:,k+1),Pk_Real_B(:,k), q_temp_RR_B(:,k))

         call mkl_dcsrsymv('l', M_Bar_B,  R_GA3_inverse_val(:,k+1), R_GA3_inverse_row(:,k+1), &
            R_GA3_inverse_col(:,k+1),Pk_Imag_B(:,k), q_temp_RI_B(:,k))

         call mkl_dcsrsymv('l', M_Bar_B,  I_GA3_inverse_val(:,k+1), I_GA3_inverse_row(:,k+1), &
            I_GA3_inverse_col(:,k+1),Pk_Real_B(:,k), q_temp_IR_B(:,k))

         call mkl_dcsrsymv('l', M_Bar_B,  I_GA3_inverse_val(:,k+1), I_GA3_inverse_row(:,k+1), &
            I_GA3_inverse_col(:,k+1),Pk_Imag_B(:,k), q_temp_II_B(:,k))
         enddo


# if defined (OPENMP)
        !$omp end  do
        !$omp end  parallel
# endif  /* OPENMP */


             do i=0,M_Bar_B-1
                do k=0,LOCAL_SIZE-1
                local_data(k+1)=cmplx(q_temp_RR_B(i,k)-q_temp_II_B(i,k), & 
                q_temp_RI_B(i,k)+q_temp_IR_B(i,k))
                enddo 
       call fftw_mpi_execute_dft(bplan, local_data, local_data)
       !   call fftwnd_f77_mpi(plancb,1,local_data,work,0,FFTW_NORMAL_ORDER)
                do k=0,LOCAL_SIZE-1
                 qB(s-NA,i,k)=real(local_data(k+1))*M_grid_inv
                enddo
              enddo 
 
  endif   !endif s>=3:s>NA
                 
endif    !endif s==1

enddo  !enddo s=1,NMAX
!!!!!!!!!!!!!!1


!!!normalization!!!
!!warning,in case M_Bar_B /= M_Bar_A, qB(0),qA(NA) must be handled very
!carefully!

        do k=0,LOCAL_SIZE-1
         do i=0,M_Bar_B-1
         qB(0,i,k)=0.0
         qBend(i,k)=qB(NB,i,k)
         enddo
       enddo


    do k=0,LOCAL_SIZE-1
       do i=0,min(M_Bar_A-1,M_Bar_B-1)
         qB(0,i,k)=qA(NA,i,k)
        enddo
     enddo



 deallocate(Pr_A)
 deallocate(Pk_Real_A)
 deallocate(Pk_Imag_A)
 deallocate(q_temp_RR_A)
 deallocate(q_temp_RI_A)
 deallocate(q_temp_IR_A)
 deallocate(q_temp_II_A)

 deallocate(Pr_B)
 deallocate(Pk_Real_B)
 deallocate(Pk_Imag_B)
 deallocate(q_temp_RR_B)
 deallocate(q_temp_RI_B)
 deallocate(q_temp_IR_B)
 deallocate(q_temp_II_B)


end subroutine MDE_q

      


subroutine MDE_qstar()
 USE nrtype,only :DP,PI
 USE global_para
 USE mpi
 USE control
 USE constants
 USE utility
 USE mmpi
 USE mpi_fftw3_operation
 USE SPH
 USE matrix_inverse
 USE G_matrix_mod
 implicit none
 integer :: i,j,s
 integer :: index_nonzero
 integer :: k,KK
 integer :: istat2
 integer,parameter :: namax=1800 !this number is obtained by a prerun of input parameters
 !for A
 real(DP) :: tempA1,tempA2,tempA3,temp,tempx,Wtemp,sumaq,sumaq_global
 real(DP),allocatable :: Pr_A(:,:)
 real(DP),allocatable :: Pk_Real_A(:,:)
 real(DP),allocatable :: Pk_Imag_A(:,:)
 real(DP),allocatable :: qstar_temp_RR_A(:,:)
 real(DP),allocatable :: qstar_temp_RI_A(:,:)
 real(DP),allocatable :: qstar_temp_IR_A(:,:)
 real(DP),allocatable :: qstar_temp_II_A(:,:)
 real(DP),allocatable :: qstar_temp_A(:,:)
 !for B
 real(DP),allocatable :: Pr_B(:,:)
 real(DP),allocatable :: Pk_Real_B(:,:)
 real(DP),allocatable :: Pk_Imag_B(:,:)
 real(DP),allocatable :: qstar_temp_RR_B(:,:)
 real(DP),allocatable :: qstar_temp_RI_B(:,:)
 real(DP),allocatable :: qstar_temp_IR_B(:,:)
 real(DP),allocatable :: qstar_temp_II_B(:,:)
 real(DP),allocatable :: qstar_temp_B1(:,:)
 real(DP),allocatable :: qstar_temp_B2(:,:)
 real(DP),allocatable :: qstar_temp_B3(:,:)

 REAL(DP) :: M_grid_inv
 !!!! just for reference
 REAL(DP) :: sy_values_R(1:namax)
 REAL(DP) :: sy_values_I(1:namax)
 integer ::  sy_columns_R(1:namax)
 integer ::  sy_columns_I(1:namax)
!!M_Bar_A >=M_BAr_B
 integer ::  sy_rowIndex_R(1:M_Bar_A+1)
 integer ::  sy_rowIndex_I(1:M_Bar_A+1)
  



!for A
 allocate(Pr_A(0:M_Bar_A-1,0:LOCAL_SIZE-1),stat=istat2)
 allocate(Pk_Real_A(0:M_Bar_A-1,0:LOCAL_SIZE-1),stat=istat2)
 allocate(Pk_Imag_A(0:M_Bar_A-1,0:LOCAL_SIZE-1),stat=istat2)


 allocate(qstar_temp_RR_A(0:M_Bar_A-1,0:LOCAL_SIZE-1),stat=istat2)
 allocate(qstar_temp_RI_A(0:M_Bar_A-1,0:LOCAL_SIZE-1),stat=istat2)
 allocate(qstar_temp_IR_A(0:M_Bar_A-1,0:LOCAL_SIZE-1),stat=istat2)
 allocate(qstar_temp_II_A(0:M_Bar_A-1,0:LOCAL_SIZE-1),stat=istat2)
 allocate(qstar_temp_A(0:M_Bar_A-1,0:LOCAL_SIZE-1),stat=istat2)
!for B
 allocate(Pr_B(0:M_Bar_B-1,0:LOCAL_SIZE-1),stat=istat2)
 allocate(Pk_Real_B(0:M_Bar_B-1,0:LOCAL_SIZE-1),stat=istat2)
 allocate(Pk_Imag_B(0:M_Bar_B-1,0:LOCAL_SIZE-1),stat=istat2)

 allocate(qstar_temp_RR_B(0:M_Bar_B-1,0:LOCAL_SIZE-1),stat=istat2)
 allocate(qstar_temp_RI_B(0:M_Bar_B-1,0:LOCAL_SIZE-1),stat=istat2)
 allocate(qstar_temp_IR_B(0:M_Bar_B-1,0:LOCAL_SIZE-1),stat=istat2)
 allocate(qstar_temp_II_B(0:M_Bar_B-1,0:LOCAL_SIZE-1),stat=istat2)
 allocate(qstar_temp_B1(0:M_Bar_B-1,0:LOCAL_SIZE-1),stat=istat2)
 allocate(qstar_temp_B2(0:M_Bar_B-1,0:LOCAL_SIZE-1),stat=istat2)
 allocate(qstar_temp_B3(0:M_Bar_B-1,0:LOCAL_SIZE-1),stat=istat2)

 if(istat2/=0) then
 write(*,*) "allocate matrix in mde_q failed"
 stop
 endif

M_grid_inv=1.0d0/M_grid



do k=0,LOCAL_SIZE-1
    qstar_temp_B1(0,K)=1.0
    qB(NB,0,k)=qB(NB,0,k)*1.0
   do i=1,M_Bar_B-1
    qstar_temp_B1(i,K)=0.0
    qB(NB,i,K)=qB(NB,i,k)*0.0
   enddo
enddo

       
 !!!!!
 do s=Nmax-1,0,-1
 !!!s=1
   if(s==Nmax-1) then

          Pr_B=0.0
# if defined (OPENMP)
        !$omp parallel default(shared) private(k,i,j,index_nonzero,Wtemp,tempx)
        !$omp do
# endif  /* OPENMP */
      do k=0,LOCAL_SIZE-1
          Wtemp=WB(k)
# if defined (MaierSaupe)
          do index_nonzero=THETA_nonzero_1D_B_indexK(k-1)+1, &
            THETA_nonzero_1D_B_indexK(k)
             i=THETA_nonzero_1D_B(index_nonzero)%i
             j=THETA_nonzero_1D_B(index_nonzero)%j
             tempx=THETA_nonzero_1D_B(index_nonzero)%value
             Pr_B(i,k)=Pr_B(i,k)+ds*tempx*qstar_temp_B1(j,k)
           enddo
# endif /* MaierSaupe */
            do i=0,M_Bar_B-1
            !Pr_B(i,k)=Pr_B(i,k)+(1.0-ds*Wtemp)*qBstar(NB,i,k)
            Pr_B(i,k)=Pr_B(i,k)+(1.0-ds*Wtemp)*qstar_temp_B1(i,k)
            enddo
        enddo !enddo k=0,LOCAL_SIZE    
# if defined (OPENMP)
        !$omp end  do
        !$omp end  parallel
# endif  /* OPENMP */

!!!doing the FFT

      do i=0,M_Bar_B-1
        do k=0,LOCAL_SIZE-1
          local_data(k+1)=cmplx(Pr_B(i,k),0.0)
         enddo
       
       
       call fftw_mpi_execute_dft(fplan, local_data, local_data)
       
        do k=0,LOCAL_SIZE-1
            Pk_Real_B(i,k)=real(local_data(k+1))
            Pk_Imag_B(i,k)=AIMAG(local_data(k+1))
         enddo
       enddo  !enddo i=0,M_Bar_B-1  


# if defined (OPENMP)
        !$omp parallel default(shared) private(k)
        !$omp do
# endif  /* OPENMP */

         do k=0,LOCAL_SIZE-1


          call  init_GA_inverse(8,k,sy_values_R, &
                sy_values_I,sy_columns_R,sy_columns_I, &
                sy_rowIndex_R,sy_rowIndex_I )
 
         call mkl_dcsrsymv('l', M_Bar_B,  sy_values_R, sy_rowIndex_R, &
            sy_columns_R,Pk_Real_B(:,k), qstar_temp_RR_B(:,k))

         call mkl_dcsrsymv('l', M_Bar_B,  sy_values_R, sy_rowIndex_R, &
            sy_columns_R,Pk_Imag_B(:,k), qstar_temp_RI_B(:,k))

         call mkl_dcsrsymv('l', M_Bar_B,  sy_values_I, sy_rowIndex_I, &
            sy_columns_I,Pk_Real_B(:,k), qstar_temp_IR_B(:,k))

         call mkl_dcsrsymv('l', M_Bar_B,  sy_values_I, sy_rowIndex_I, &
            sy_columns_I,Pk_Imag_B(:,k), qstar_temp_II_B(:,k))



         enddo
# if defined (OPENMP)
        !$omp end  do
        !$omp end  parallel
# endif  /* OPENMP */


             do i=0,M_Bar_B-1
                do k=0,LOCAL_SIZE-1
                local_data(k+1)=cmplx(qstar_temp_RR_B(i,k)-qstar_temp_II_B(i,k), & 
                qstar_temp_RI_B(i,k)+qstar_temp_IR_B(i,k))

                enddo 
       call fftw_mpi_execute_dft(bplan, local_data, local_data)
                do k=0,LOCAL_SIZE-1
                 qstar_temp_B2(i,k)=qstar_temp_B1(i,k)
                 qstar_temp_B1(i,k)=real(local_data(k+1))*M_grid_inv
                 qB(NB-1,i,k)=qB(NB-1,i,k)*qstar_temp_B1(i,k)
                enddo
              enddo  

!!!!s=NMAX-2!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    else if(s==NMAX-2) then

          Pr_B=0.0
# if defined (OPENMP)
        !$omp parallel default(shared) private(k,i,j,index_nonzero,Wtemp,tempx)
        !$omp do
# endif  /* OPENMP */
      do k=0,LOCAL_SIZE-1
          Wtemp=WB(k)
# if defined (MaierSaupe)
          do index_nonzero=THETA_nonzero_1D_B_indexK(k-1)+1, &
            THETA_nonzero_1D_B_indexK(k)
             i=THETA_nonzero_1D_B(index_nonzero)%i
             j=THETA_nonzero_1D_B(index_nonzero)%j
             tempx=THETA_nonzero_1D_B(index_nonzero)%value
             Pr_B(i,k)=Pr_B(i,k)+ds*tempx*(2.0*qstar_temp_B1(j,k)-qstar_temp_B2(j,k))
           enddo
# endif /* MaierSaupe */          
            do i=0,M_Bar_B-1
            !Pr_B(i,k)=Pr_B(i,k)+2.0*qBstar(NB-1,i,k)-0.5*qBstar(NB,i,k) - &
            !          ds*Wtemp*(2.0*qBstar(NB-1,i,k)-qBstar(NB,i,k))
            Pr_B(i,k)=Pr_B(i,k)+2.0*qstar_temp_B1(i,k)-0.5*qstar_temp_B2(i,k) - &
                      ds*Wtemp*(2.0*qstar_temp_B1(i,k)-qstar_temp_B2(i,k))
            enddo
        enddo !enddo k=0,LOCAL_SIZE    
# if defined (OPENMP)
        !$omp end  do
        !$omp end  parallel
# endif  /* OPENMP */

!!!doing the FFT

      do i=0,M_Bar_B-1
        do k=0,LOCAL_SIZE-1
          local_data(k+1)=cmplx(Pr_B(i,k),0.0)
         enddo
       
       
       call fftw_mpi_execute_dft(fplan, local_data, local_data)
       
        do k=0,LOCAL_SIZE-1
            Pk_Real_B(i,k)=real(local_data(k+1))
            Pk_Imag_B(i,k)=AIMAG(local_data(k+1))
         enddo
       enddo  !enddo i=0,M_Bar_B-1  


# if defined (OPENMP)
        !$omp parallel default(shared) private(k)
        !$omp do
# endif  /* OPENMP */
         do k=0,LOCAL_SIZE-1

          call  init_GA_inverse(7,k,sy_values_R, &
                sy_values_I,sy_columns_R,sy_columns_I, &
                sy_rowIndex_R,sy_rowIndex_I )
 
         call mkl_dcsrsymv('l', M_Bar_B,  sy_values_R, sy_rowIndex_R, &
            sy_columns_R,Pk_Real_B(:,k), qstar_temp_RR_B(:,k))

         call mkl_dcsrsymv('l', M_Bar_B,  sy_values_R, sy_rowIndex_R, &
            sy_columns_R,Pk_Imag_B(:,k), qstar_temp_RI_B(:,k))

         call mkl_dcsrsymv('l', M_Bar_B,  sy_values_I, sy_rowIndex_I, &
            sy_columns_I,Pk_Real_B(:,k), qstar_temp_IR_B(:,k))

         call mkl_dcsrsymv('l', M_Bar_B,  sy_values_I, sy_rowIndex_I, &
            sy_columns_I,Pk_Imag_B(:,k), qstar_temp_II_B(:,k))

         enddo
# if defined (OPENMP)
        !$omp end  do
        !$omp end  parallel
# endif  /* OPENMP */

             do i=0,M_Bar_B-1
                do k=0,LOCAL_SIZE-1
                local_data(k+1)=cmplx(qstar_temp_RR_B(i,k)-qstar_temp_II_B(i,k), & 
                qstar_temp_RI_B(i,k)+qstar_temp_IR_B(i,k))
                enddo 
                
       call fftw_mpi_execute_dft(bplan, local_data, local_data)
                do k=0,LOCAL_SIZE-1
                 qstar_temp_B3(i,k)=qstar_temp_B2(i,k)
                 qstar_temp_B2(i,k)=qstar_temp_B1(i,k)
                 qstar_temp_B1(i,k)=real(local_data(k+1))*M_grid_inv
                 qB(NB-2,i,k)=qB(NB-2,i,k)*qstar_temp_B1(i,k)
                enddo
              enddo  

!!!!s<=NMAX-3!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    else if(s<=NMAX-3) then
          temp=1.0/3.0
           
          if(s>NA) then
            
          Pr_B=0.0
# if defined (OPENMP)
        !$omp parallel default(shared) private(k,i,j,index_nonzero,Wtemp,tempx)
        !$omp do
# endif  /* OPENMP */
      do k=0,LOCAL_SIZE-1
          Wtemp=WB(k)
# if defined (MaierSaupe)
          do index_nonzero=THETA_nonzero_1D_B_indexK(k-1)+1, &
            THETA_nonzero_1D_B_indexK(k)
             i=THETA_nonzero_1D_B(index_nonzero)%i
             j=THETA_nonzero_1D_B(index_nonzero)%j
             tempx=THETA_nonzero_1D_B(index_nonzero)%value
             Pr_B(i,k)=Pr_B(i,k)+ds*tempx*(3.0*qstar_temp_B1(j,k)- &
                3.0*qstar_temp_B2(j,k)+qstar_temp_B3(j,k))
            
           enddo
# endif /* MaierSaupe */
            do i=0,M_Bar_B-1
            !Pr_B(i,k)=Pr_B(i,k)+3.0*qBstar(s-NA+1,i,k)-1.5*qBstar(s-NA+2,i,k)+ &
           !           temp*qBstar(s-NA+3,i,k) - &
           !  ds*Wtemp*(3.0*qBstar(s-NA+1,i,k)-3.0*qBstar(s-NA+2,i,k)+qBstar(s-NA+3,i,k))
            Pr_B(i,k)=Pr_B(i,k)+3.0*qstar_temp_B1(i,k)-1.5*qstar_temp_B2(i,k)+ &
                      temp*qstar_temp_B3(i,k) - &
             ds*Wtemp*(3.0*qstar_temp_B1(i,k)-3.0*qstar_temp_B2(i,k)+qstar_temp_B3(i,k))
            enddo
        enddo !enddo k=0,LOCAL_SIZE    
# if defined (OPENMP)
        !$omp end  do
        !$omp end  parallel
# endif  /* OPENMP */

!!!doing the FFT

      do i=0,M_Bar_B-1
        do k=0,LOCAL_SIZE-1
          local_data(k+1)=cmplx(Pr_B(i,k),0.0)
         enddo
       
       
       call fftw_mpi_execute_dft(fplan, local_data, local_data)
       
        do k=0,LOCAL_SIZE-1
            Pk_Real_B(i,k)=real(local_data(k+1))
            Pk_Imag_B(i,k)=AIMAG(local_data(k+1))
           ! Pk_cmplx_B(i,k)=local_data(k)
         enddo
       enddo  !enddo i=0,M_Bar_B-1  


# if defined (OPENMP)
        !$omp parallel default(shared) private(k)
        !$omp do
# endif  /* OPENMP */

        do k=0,LOCAL_SIZE-1
         call mkl_dcsrsymv('l', M_Bar_B,  R_GA3_star_inverse_val(:,k+1), R_GA3_star_inverse_row(:,k+1), &
            R_GA3_star_inverse_col(:,k+1),Pk_Real_B(:,k), qstar_temp_RR_B(:,k))

         call mkl_dcsrsymv('l', M_Bar_B,  R_GA3_star_inverse_val(:,k+1), R_GA3_star_inverse_row(:,k+1), &
            R_GA3_star_inverse_col(:,k+1),Pk_Imag_B(:,k), qstar_temp_RI_B(:,k))

         call mkl_dcsrsymv('l', M_Bar_B,  I_GA3_star_inverse_val(:,k+1), I_GA3_star_inverse_row(:,k+1), &
            I_GA3_star_inverse_col(:,k+1),Pk_Real_B(:,k), qstar_temp_IR_B(:,k))

         call mkl_dcsrsymv('l', M_Bar_B,  I_GA3_star_inverse_val(:,k+1), I_GA3_star_inverse_row(:,k+1), &
            I_GA3_star_inverse_col(:,k+1),Pk_Imag_B(:,k), qstar_temp_II_B(:,k))


         enddo
# if defined (OPENMP)
        !$omp end  do
        !$omp end  parallel
# endif  /* OPENMP */


             do i=0,M_Bar_B-1
                do k=0,LOCAL_SIZE-1
                local_data(k+1)=cmplx(qstar_temp_RR_B(i,k)-qstar_temp_II_B(i,k), & 
                qstar_temp_RI_B(i,k)+qstar_temp_IR_B(i,k))
                enddo 
       call fftw_mpi_execute_dft(bplan, local_data, local_data)
       !   call fftwnd_f77_mpi(plancb,1,local_data,work,0,FFTW_NORMAL_ORDER)
                do k=0,LOCAL_SIZE-1
                 qstar_temp_B3(i,k)=qstar_temp_B2(i,k)
                 qstar_temp_B2(i,k)=qstar_temp_B1(i,k)
                 qstar_temp_B1(i,k)=real(local_data(k+1))*M_grid_inv
                 qB(s-NA,i,k)=qB(s-NA,i,k)*qstar_temp_B1(i,k)
                enddo
              enddo  



   else  ! else s<=NA :entering A block regime
            !!!calculate polymer block A
             Pr_A=0.0
            
   if(s==NA) then
           Pr_B=0.0  
# if defined (OPENMP)
        !$omp parallel default(shared) private(k,i,j,index_nonzero,Wtemp,tempx)
        !$omp do
# endif  /* OPENMP */
      do k=0,LOCAL_SIZE-1
          Wtemp=WA(k)
# if defined (MaierSaupe)
          do index_nonzero=THETA_nonzero_1D_B_indexK(k-1)+1, &
            THETA_nonzero_1D_B_indexK(k)
             i=THETA_nonzero_1D_B(index_nonzero)%i
             j=THETA_nonzero_1D_B(index_nonzero)%j
             tempx=THETA_nonzero_1D_B(index_nonzero)%value

              Pr_B(i,k)=Pr_B(i,k) + ds*tempx*(3.0*qstar_temp_B1(j,k) - & 
              3.0*qstar_temp_B2(j,k)+qstar_temp_B3(j,k))
           enddo  !enddo index_nonzero
# endif /* MaierSaupe */


            do i=0,min(M_Bar_A-1,M_Bar_B-1)
                Pr_A(i,k)=Pr_B(i,k)

            Pr_A(i,k)=Pr_A(i,k)+3.0*qstar_temp_B1(i,k)-1.5*qstar_temp_B2(i,k)+ &
                      temp*qstar_temp_B3(i,k) - &
                      ds*Wtemp*(3.0*qstar_temp_B1(i,k)-3.0*qstar_temp_B2(i,k)+qstar_temp_B3(i,k))
            enddo
          enddo !enddo k=0,LOCAL_SIZE    
# if defined (OPENMP)
        !$omp end  do
        !$omp end  parallel
# endif  /* OPENMP */

    else if(s==(NA-1)) then
         
            Pr_B=0.0

# if defined (OPENMP)
        !$omp parallel default(shared) private(k,i,j,index_nonzero,Wtemp,tempx)
        !$omp do
# endif  /* OPENMP */
      do k=0,LOCAL_SIZE-1
          Wtemp=WA(k)
# if defined (MaierSaupe)
          do index_nonzero=THETA_nonzero_1D_B_indexK(k-1)+1, &
            THETA_nonzero_1D_B_indexK(k)
             i=THETA_nonzero_1D_B(index_nonzero)%i
             j=THETA_nonzero_1D_B(index_nonzero)%j
             tempx=THETA_nonzero_1D_B(index_nonzero)%value
             tempA1=0.0
             
             if(j<M_Bar_A) then
              tempA1=qstar_temp_B1(j,k)
          
              endif
              Pr_B(i,k)=Pr_B(i,k) + ds*tempx*(3.0*tempA1-3.0*qstar_temp_B2(j,k)+ & 
               qstar_temp_B3(j,k))
           enddo  !enddo index_nonzero
# endif /* MaierSaupe */            

            do i=0,M_Bar_A-1
              Pr_A(i,k)=Pr_B(i,k)

            !Pr_A(i,k)=Pr_A(i,k)+3.0*qAstar(NA,i,k)-1.5*qBstar(1,i,k)+ &
            !          temp*qBstar(2,i,k) - &
            !          ds*Wtemp*(3.0*qAstar(NA,i,k)-3.0*qBstar(1,i,k)+qBstar(2,i,k))
            Pr_A(i,k)=Pr_A(i,k)+3.0*qstar_temp_B1(i,k)-1.5*qstar_temp_B2(i,k)+ &
                      temp*qstar_temp_B3(i,k) - &
                      ds*Wtemp*(3.0*qstar_temp_B1(i,k)-3.0*qstar_temp_B2(i,k)+qstar_temp_B3(i,k))
            enddo
          enddo !enddo k=0,LOCAL_SIZE    
# if defined (OPENMP)
        !$omp end  do
        !$omp end  parallel
# endif  /* OPENMP */

     else if(s==(NA-2)) then
         Pr_B=0.0
# if defined (OPENMP)
        !$omp parallel default(shared) private(k,i,j,index_nonzero,Wtemp,tempx)
        !$omp do
# endif  /* OPENMP */
      do k=0,LOCAL_SIZE-1
          Wtemp=WA(k)
# if defined (MaierSaupe)
          do index_nonzero=THETA_nonzero_1D_B_indexK(k-1)+1, &
            THETA_nonzero_1D_B_indexK(k)
             i=THETA_nonzero_1D_B(index_nonzero)%i
             j=THETA_nonzero_1D_B(index_nonzero)%j
             tempx=THETA_nonzero_1D_B(index_nonzero)%value
             tempA1=0.0
             tempA2=0.0
             if(j<M_Bar_A) then
              tempA1=qstar_temp_B1(j,k)
              tempA2=qstar_temp_B2(j,k)
              endif
              Pr_B(i,k)=Pr_B(i,k) + ds*tempx*(3.0*tempA1-3.0*tempA2+qstar_temp_B3(j,k))
           enddo  !enddo index_nonzero
# endif /* MaierSaupe */

            do i=0,M_Bar_A-1
                Pr_A(i,k)=Pr_B(i,k)
                 

            !Pr_A(i,k)=Pr_A(i,k)+3.0*qAstar(NA-1,i,k)-1.5*qAstar(NA,i,k)+ &
            !          temp*qBstar(1,i,k)- &
            !          ds*Wtemp*(3.0*qAstar(NA-1,i,k)-3.0*qAstar(NA,i,k)+qBstar(1,i,k))
            Pr_A(i,k)=Pr_A(i,k)+3.0*qstar_temp_B1(i,k)-1.5*qstar_temp_B2(i,k)+ &
                      temp*qstar_temp_B3(i,k) - &
                      ds*Wtemp*(3.0*qstar_temp_B1(i,k)-3.0*qstar_temp_B2(i,k)+qstar_temp_B3(i,k))
            enddo
          enddo !enddo k=0,LOCAL_SIZE    
# if defined (OPENMP)
        !$omp end  do
        !$omp end  parallel
# endif  /* OPENMP */

     else 
           
# if defined (OPENMP)
        !$omp parallel default(shared) private(k,i,j,index_nonzero,Wtemp,tempx)
        !$omp do
# endif  /* OPENMP */
      do k=0,LOCAL_SIZE-1
          Wtemp=WA(k)
# if defined (MaierSaupe)
          do index_nonzero=THETA_nonzero_1D_A_indexK(k-1)+1, &
            THETA_nonzero_1D_A_indexK(k)
             i=THETA_nonzero_1D_A(index_nonzero)%i
             j=THETA_nonzero_1D_A(index_nonzero)%j
             tempx=THETA_nonzero_1D_A(index_nonzero)%value
       Pr_A(i,k)=Pr_A(i,k) + ds*tempx*(3.0*qstar_temp_B1(j,k)- & 
           3.0*qstar_temp_B2(j,k)+qstar_temp_B3(j,k))
           enddo  !enddo index_nonzero
# endif /* MaierSaupe */

            do i=0,M_Bar_A-1

            !Pr_A(i,k)=Pr_A(i,k)+3.0*qAstar(s+1,i,k)-1.5*qAstar(s+2,i,k)+ &
            !          temp*qAstar(s+3,i,k) - &
            !          ds*Wtemp*(3.0*qAstar(s+1,i,k)-3.0*qAstar(s+2,i,k)+qAstar(s+3,i,k))
            Pr_A(i,k)=Pr_A(i,k)+3.0*qstar_temp_B1(i,k)-1.5*qstar_temp_B2(i,k)+ &
                      temp*qstar_temp_B3(i,k) - &
                      ds*Wtemp*(3.0*qstar_temp_B1(i,k)-3.0*qstar_temp_B2(i,k)+qstar_temp_B3(i,k))
            enddo
          enddo !enddo k=0,LOCAL_SIZE    
# if defined (OPENMP)
        !$omp end  do
        !$omp end  parallel
# endif  /* OPENMP */
     endif   

!!!doing the FFT

      do i=0,M_Bar_A-1
        do k=0,LOCAL_SIZE-1
          local_data(k+1)=cmplx(Pr_A(i,k),0.0)
         enddo
       
       
       call fftw_mpi_execute_dft(fplan, local_data, local_data)
       
        do k=0,LOCAL_SIZE-1
            Pk_Real_A(i,k)=real(local_data(k+1))
            Pk_Imag_A(i,k)=AIMAG(local_data(k+1))
         enddo
       enddo  !enddo i=0,M_Bar_A-1  


# if defined (OPENMP)
        !$omp parallel default(shared) private(k)
        !$omp do
# endif  /* OPENMP */
         do k=0,LOCAL_SIZE-1

         call mkl_dcsrsymv('l', M_Bar_A,  R_GA3_star_inverse_val(:,k+1), R_GA3_star_inverse_row(:,k+1), &
            R_GA3_star_inverse_col(:,k+1),Pk_Real_A(:,k), qstar_temp_RR_A(:,k))

         call mkl_dcsrsymv('l', M_Bar_A,  R_GA3_star_inverse_val(:,k+1), R_GA3_star_inverse_row(:,k+1), &
            R_GA3_star_inverse_col(:,k+1),Pk_Imag_A(:,k), qstar_temp_RI_A(:,k))

         call mkl_dcsrsymv('l', M_Bar_A,  I_GA3_star_inverse_val(:,k+1), I_GA3_star_inverse_row(:,k+1), &
            I_GA3_star_inverse_col(:,k+1),Pk_Real_A(:,k), qstar_temp_IR_A(:,k))

         call mkl_dcsrsymv('l', M_Bar_A,  I_GA3_star_inverse_val(:,k+1), I_GA3_star_inverse_row(:,k+1), &
            I_GA3_star_inverse_col(:,k+1),Pk_Imag_A(:,k), qstar_temp_II_A(:,k))

         enddo
# if defined (OPENMP)
        !$omp end  do
        !$omp end  parallel
# endif  /* OPENMP */


             do i=0,M_Bar_A-1
                do k=0,LOCAL_SIZE-1
                local_data(k+1)=cmplx(qstar_temp_RR_A(i,k)-qstar_temp_II_A(i,k), & 
                qstar_temp_RI_A(i,k)+qstar_temp_IR_A(i,k))
                enddo 
       call fftw_mpi_execute_dft(bplan, local_data, local_data)
                do k=0,LOCAL_SIZE-1
                 qstar_temp_B3(i,k)=qstar_temp_B2(i,k)
                 qstar_temp_B2(i,k)=qstar_temp_B1(i,k)
                 qstar_temp_B1(i,k)=real(local_data(k+1))*M_grid_inv
                 qA(s,i,k)=qA(s,i,k)*qstar_temp_B1(i,k)
                enddo
              enddo  
  endif   !endif 
                 
endif    !endif 

enddo  !enddo s=NMAX-1,0,-1

!
!!warning,in case M_Bar_B /= M_Bar_A, qB(0),qA(NA) must be handled very
!carefully!
      do k=0,LOCAL_SIZE-1
       do i=0,M_Bar_B-1
         qB(0,i,k)=0.0
        enddo
     enddo

    do k=0,LOCAL_SIZE-1
       do i=0,min(M_Bar_A-1,M_Bar_B-1)
         !qBstar(0,i,k)=qAstar(NA,i,k)
         qB(0,i,k)=qA(NA,i,k)
        enddo
     enddo

 deallocate(Pr_A)
 deallocate(Pk_Real_A)
 deallocate(Pk_Imag_A)
 deallocate(qstar_temp_RR_A)
 deallocate(qstar_temp_RI_A)
 deallocate(qstar_temp_IR_A)
 deallocate(qstar_temp_II_A)
 deallocate(qstar_temp_A)
!
 deallocate(Pr_B)
 deallocate(Pk_Real_B)
 deallocate(Pk_Imag_B)
 deallocate(qstar_temp_RR_B)
 deallocate(qstar_temp_RI_B)
 deallocate(qstar_temp_IR_B)
 deallocate(qstar_temp_II_B)
 deallocate(qstar_temp_B1)
 deallocate(qstar_temp_B2)
 deallocate(qstar_temp_B3)


end subroutine MDE_qstar

# else  /* Mem_OP */

 subroutine MDE_q()
 USE nrtype,only :DP,PI
 USE global_para
 USE mpi
 USE control
 USE constants
 USE utility
 USE mmpi
 USE mpi_fftw3_operation
 USE SPH
 USE matrix_inverse
 USE G_matrix_mod
 implicit none
 integer :: i,j,s
 integer :: index_nonzero
 integer :: k
 integer :: istat2
 integer,parameter :: namaxB=256 !this number is obtained by a prerun of input parameters
 integer,parameter :: namaxA=1800 !this number is obtained by a prerun of input parameters
 !for A
 real(DP) :: tempA1,tempA2,tempA3,temp,tempx,Wtemp,sumaq,sumaq_global,sumaq1
 real(DP) :: sumPA,sumPA_global
 real(DP),allocatable :: Pr_A(:,:)
 real(DP),allocatable :: Pk_Real_A(:,:)
 real(DP),allocatable :: Pk_Imag_A(:,:)
 real(DP),allocatable :: q_temp_RR_A(:,:)
 real(DP),allocatable :: q_temp_RI_A(:,:)
 real(DP),allocatable :: q_temp_IR_A(:,:)
 real(DP),allocatable :: q_temp_II_A(:,:)
 !for B
 real(DP),allocatable :: Pr_B(:,:)
 real(DP),allocatable :: Pk_Real_B(:,:)
 real(DP),allocatable :: Pk_Imag_B(:,:)
 real(DP),allocatable :: q_temp_RR_B(:,:)
 real(DP),allocatable :: q_temp_RI_B(:,:)
 real(DP),allocatable :: q_temp_IR_B(:,:)
 real(DP),allocatable :: q_temp_II_B(:,:)
 REAL(DP) :: M_grid_inv

 !first ,allocate the matrixs
!for A
 allocate(Pr_A(0:M_Bar_A-1,0:LOCAL_SIZE-1),stat=istat2)
 allocate(Pk_Real_A(0:M_Bar_A-1,0:LOCAL_SIZE-1),stat=istat2)
 allocate(Pk_Imag_A(0:M_Bar_A-1,0:LOCAL_SIZE-1),stat=istat2)

 allocate(q_temp_RR_A(0:M_Bar_A-1,0:LOCAL_SIZE-1),stat=istat2)
 allocate(q_temp_RI_A(0:M_Bar_A-1,0:LOCAL_SIZE-1),stat=istat2)
 allocate(q_temp_IR_A(0:M_Bar_A-1,0:LOCAL_SIZE-1),stat=istat2)
 allocate(q_temp_II_A(0:M_Bar_A-1,0:LOCAL_SIZE-1),stat=istat2)
!for B
 allocate(Pr_B(0:M_Bar_B-1,0:LOCAL_SIZE-1),stat=istat2)
 allocate(Pk_Real_B(0:M_Bar_B-1,0:LOCAL_SIZE-1),stat=istat2)
 allocate(Pk_Imag_B(0:M_Bar_B-1,0:LOCAL_SIZE-1),stat=istat2)

 allocate(q_temp_RR_B(0:M_Bar_B-1,0:LOCAL_SIZE-1),stat=istat2)
 allocate(q_temp_RI_B(0:M_Bar_B-1,0:LOCAL_SIZE-1),stat=istat2)
 allocate(q_temp_IR_B(0:M_Bar_B-1,0:LOCAL_SIZE-1),stat=istat2)
 allocate(q_temp_II_B(0:M_Bar_B-1,0:LOCAL_SIZE-1),stat=istat2)
  

 if(istat2/=0) then
 write(*,*) "allocate matrix in mde_q failed"
 stop
 endif

 M_grid_inv=1.0d0/M_grid

do k=0,LOCAL_SIZE-1
   do i=0,M_Bar_A-1
    qA(0,i,K)=0.0
    if(i==0) qA(0,i,K)=1.0
   enddo
enddo

# if defined (Debug)
# if defined (Dim2)
# else /* Dim2 */
if(myid==0) then
write(*,*) "0:WA(0,90,180) in MDE=",WA(0),WA(90),WA(180)
write(*,*) "0:WA(300) in MDE=",WA(300)
write(*,*) "0:WB(0,90,180)=",WB(0),WB(90),WB(180)
write(*,*) "0:WB(300) in MDE=",WB(300)
endif  
#endif /* Dim2 */
!if(myid==2) then
!write(*,*) " 2: WA(0,30,300) in MDE=",WA(0),WA(30),WA(90)
!write(*,*) "2: WB(0,30,300)=",WB(0),WB(30),WB(90)
!endif  
!if(myid==0) then
!write(*,*) "M_OP(0)(30) 00=",M_OP(0,0,0),M_OP(30,0,0)
!write(*,*) "M_OP(0)(30) 11=",M_OP(0,1,1),M_OP(30,1,1)
!endif  
# endif  /* Debug */
       
 !!!!!
 do s=1,Nmax

 !!!s=1
   if(s==1) then

          Pr_A=0.0

# if defined (OPENMP)
        !$omp parallel default(shared) private(k,i,j,index_nonzero,Wtemp,tempx)
        !$omp do
# endif  /* OPENMP */
     
      do k=0,LOCAL_SIZE-1
          Wtemp=WA(k)
# if defined (MaierSaupe)
          do index_nonzero=THETA_nonzero_1D_A_indexK(k-1)+1, &
            THETA_nonzero_1D_A_indexK(k)
             i=THETA_nonzero_1D_A(index_nonzero)%i
             j=THETA_nonzero_1D_A(index_nonzero)%j
             tempx=THETA_nonzero_1D_A(index_nonzero)%value
             Pr_A(i,k)=Pr_A(i,k)+ds*tempx*qA(0,j,k)
           enddo
# endif /* MaierSaupe */

            do i=0,M_Bar_A-1
            Pr_A(i,k)=Pr_A(i,k)+(1.0-ds*Wtemp)*qA(0,i,k)
            enddo
        enddo !enddo k=0,LOCAL_SIZE  

# if defined (Debug)
    call mp_barrier()
    if(myid==0) then
     write(*,*) "Pr_A(0,30)",Pr_A(0,30)
     endif
# endif /* Debug */




  
# if defined (OPENMP)
        !$omp end  do
        !$omp end  parallel
# endif  /* OPENMP */

!!!doing the FFT

      do i=0,M_Bar_A-1
        do k=0,LOCAL_SIZE-1
# if defined (Dim2)
          local_data(kD(k)%x,kD(k)%y)=cmplx(Pr_A(i,k),0.0)
# else  /* Dim2 */
          local_data(kD(k)%x,kD(k)%y,KD(k)%z)=cmplx(Pr_A(i,k),0.0)
# endif /* Dim2 */
         enddo

!call mp_barrier()
!if(myid==0) then
!write(*,*) "done local_data copy 1st"
!endif
!call mp_barrier()
!stop       
       
       call fftw_mpi_execute_dft(fplan, local_data, local_data)
       
        do k=0,LOCAL_SIZE-1
# if defined (Dim2)
            Pk_Real_A(i,k)=real(local_data(kD(k)%x,kD(k)%y))
            Pk_Imag_A(i,k)=AIMAG(local_data(kD(k)%x,kD(k)%y))
# else  /* Dim2 */
            Pk_Real_A(i,k)=real(local_data(kD(k)%x,kD(k)%y,KD(k)%z))
            Pk_Imag_A(i,k)=AIMAG(local_data(kD(k)%x,kD(k)%y,KD(k)%z))
# endif /* Dim2 */
         enddo
       enddo  !enddo i=0,M_Bar_A-1  

# if defined (Debug)
    call mp_barrier()
    if(myid==0) then
     write(*,*) "Pk_REAL_A(0,30)",Pk_Real_A(0,30)
     endif
# endif /* Debug */

# if defined (Debug)
    call mp_barrier()
    if(myid==0) then
     write(*,*) "bf matrix vector",R_GA1_inverse_val(1,1)
     endif
    call mp_barrier()
# endif /* Debug */


# if defined (OPENMP)
        !$omp parallel default(shared) private(k)
        !$omp do
# endif  /* OPENMP */
         do k=0,LOCAL_SIZE-1

         call mkl_dcsrsymv('l', M_Bar_A,  R_GA1_inverse_val(:,k+1), R_GA1_inverse_row(:,k+1), &
            R_GA1_inverse_col(:,k+1),Pk_Real_A(:,k), q_temp_RR_A(:,k))

         call mkl_dcsrsymv('l', M_Bar_A,  R_GA1_inverse_val(:,k+1), R_GA1_inverse_row(:,k+1), &
            R_GA1_inverse_col(:,k+1),Pk_Imag_A(:,k), q_temp_RI_A(:,k))

         call mkl_dcsrsymv('l', M_Bar_A,  I_GA1_inverse_val(:,k+1), I_GA1_inverse_row(:,k+1), &
            I_GA1_inverse_col(:,k+1),Pk_Real_A(:,k), q_temp_IR_A(:,k))

         call mkl_dcsrsymv('l', M_Bar_A,  I_GA1_inverse_val(:,k+1), I_GA1_inverse_row(:,k+1), &
            I_GA1_inverse_col(:,k+1),Pk_Imag_A(:,k), q_temp_II_A(:,k))
            

         enddo

# if defined (Debug)
    call mp_barrier()
    if(myid==0) then
     write(*,*) "done matrix vector",R_GA1_inverse_val(1,1)
     endif
    call mp_barrier()
# endif /* Debug */

# if defined (OPENMP)
        !$omp end  do
        !$omp end  parallel
# endif  /* OPENMP */

             do i=0,M_Bar_A-1
                do k=0,LOCAL_SIZE-1
# if defined (Dim2)
        local_data(kD(k)%x,kD(k)%y)=cmplx(q_temp_RR_A(i,k)-q_temp_II_A(i,k), & 
                q_temp_RI_A(i,k)+q_temp_IR_A(i,k))
# else  /* Dim2 */
        local_data(kD(k)%x,kD(k)%y,KD(k)%z)=cmplx(q_temp_RR_A(i,k)-q_temp_II_A(i,k), & 
                q_temp_RI_A(i,k)+q_temp_IR_A(i,k))
# endif  /* Dim2 */
                enddo
       call fftw_mpi_execute_dft(bplan, local_data, local_data)
            
                do k=0,LOCAL_SIZE-1
# if defined (Dim2)
                 qA(1,i,k)=real(local_data(kD(k)%x,kD(k)%y))*M_grid_inv
# else  /* Dim2 */
                 qA(1,i,k)=real(local_data(kD(k)%x,kD(k)%y,KD(k)%z))*M_grid_inv
# endif  /* Dim2 */
                enddo
              enddo  

# if defined (Debug)
    call mp_barrier()
    if(myid==0) then
     write(*,*) "qA(1,0,30)",qA(1,0,30)
     endif
# endif /* Debug */

!!!!s=2!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    else if(s==2) then

          Pr_A=0.0
# if defined (OPENMP)
        !$omp parallel default(shared) private(k,i,j,index_nonzero,Wtemp,tempx)
        !$omp do
# endif  /* OPENMP */
      do k=0,LOCAL_SIZE-1
          Wtemp=WA(k)
# if defined (MaierSaupe)
          do index_nonzero=THETA_nonzero_1D_A_indexK(k-1)+1, &
            THETA_nonzero_1D_A_indexK(k)
             i=THETA_nonzero_1D_A(index_nonzero)%i
             j=THETA_nonzero_1D_A(index_nonzero)%j
             tempx=THETA_nonzero_1D_A(index_nonzero)%value
             Pr_A(i,k)=Pr_A(i,k)+ds*tempx*(2.0*qA(1,j,k)-qA(0,j,k))
           enddo
# endif /* MaierSaupe */
            do i=0,M_Bar_A-1
            Pr_A(i,k)=Pr_A(i,k)+2.0*qA(1,i,k)-0.5*qA(0,i,k) - &
                      ds*Wtemp*(2.0*qA(1,i,k)-qA(0,i,k))
            enddo
        enddo !enddo k=0,LOCAL_SIZE    
# if defined (OPENMP)
        !$omp end  do
        !$omp end  parallel
# endif  /* OPENMP */

!!!doing the FFT

      do i=0,M_Bar_A-1
        do k=0,LOCAL_SIZE-1
# if defined (Dim2)
          local_data(kD(k)%x,kD(k)%y)=cmplx(Pr_A(i,k),0.0)
# else  /* Dim2 */
          local_data(kD(k)%x,kD(k)%y,KD(k)%z)=cmplx(Pr_A(i,k),0.0)
# endif /* Dim2 */
         enddo
       
       call fftw_mpi_execute_dft(fplan, local_data, local_data)
       
        do k=0,LOCAL_SIZE-1
# if defined (Dim2)
            Pk_Real_A(i,k)=real(local_data(kD(k)%x,kD(k)%y))
            Pk_Imag_A(i,k)=AIMAG(local_data(kD(k)%x,kD(k)%y))
# else  /* Dim2 */
            Pk_Real_A(i,k)=real(local_data(kD(k)%x,kD(k)%y,KD(k)%z))
            Pk_Imag_A(i,k)=AIMAG(local_data(kD(k)%x,kD(k)%y,KD(k)%z))
# endif /* Dim2 */
         enddo
       enddo  !enddo i=0,M_Bar_A-1  


# if defined (OPENMP)
        !$omp parallel default(shared) private(k)
        !$omp do
# endif  /* OPENMP */
         do k=0,LOCAL_SIZE-1
           
         call mkl_dcsrsymv('l', M_Bar_A,  R_GA2_inverse_val(:,k+1), R_GA2_inverse_row(:,k+1), &
            R_GA2_inverse_col(:,k+1),Pk_Real_A(:,k), q_temp_RR_A(:,k))

         call mkl_dcsrsymv('l', M_Bar_A,  R_GA2_inverse_val(:,k+1), R_GA2_inverse_row(:,k+1), &
            R_GA2_inverse_col(:,k+1),Pk_Imag_A(:,k), q_temp_RI_A(:,k))

         call mkl_dcsrsymv('l', M_Bar_A,  I_GA2_inverse_val(:,k+1), I_GA2_inverse_row(:,k+1), &
            I_GA2_inverse_col(:,k+1),Pk_Real_A(:,k), q_temp_IR_A(:,k))

         call mkl_dcsrsymv('l', M_Bar_A,  I_GA2_inverse_val(:,k+1), I_GA2_inverse_row(:,k+1), &
            I_GA2_inverse_col(:,k+1),Pk_Imag_A(:,k), q_temp_II_A(:,k))
         enddo
# if defined (OPENMP)
        !$omp end  do
        !$omp end  parallel
# endif  /* OPENMP */


             do i=0,M_Bar_A-1
                do k=0,LOCAL_SIZE-1
# if defined (Dim2)
         local_data(kD(k)%x,kD(k)%y)=cmplx(q_temp_RR_A(i,k)-q_temp_II_A(i,k), & 
                q_temp_RI_A(i,k)+q_temp_IR_A(i,k))
# else  /* Dim2 */
         local_data(kD(k)%x,kD(k)%y,KD(k)%z)=cmplx(q_temp_RR_A(i,k)-q_temp_II_A(i,k), & 
                q_temp_RI_A(i,k)+q_temp_IR_A(i,k))
# endif /* Dim2 */
                enddo 
       call fftw_mpi_execute_dft(bplan, local_data, local_data)
                do k=0,LOCAL_SIZE-1
# if defined (Dim2)
                 qA(2,i,k)=real(local_data(kD(k)%x,kD(k)%y))*M_grid_inv
# else  /* Dim2 */
                 qA(2,i,k)=real(local_data(kD(k)%x,kD(k)%y,KD(k)%z))*M_grid_inv
# endif /* Dim2 */
                enddo
              enddo  

# if defined (Debug)
    call mp_barrier()
    if(myid==0) then
     write(*,*) "qA(2,0,30)",qA(2,0,30)
     endif
# endif /* Debug */
!!!!s>=3!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    else if(s>=3) then
          temp=1.0/3.0
           
          if(s<=NA) then
            
          Pr_A=0.0
# if defined (OPENMP)
        !$omp parallel default(shared) private(k,i,j,index_nonzero,Wtemp,tempx)
        !$omp do
# endif  /* OPENMP */
      do k=0,LOCAL_SIZE-1
          Wtemp=WA(k)
# if defined (MaierSaupe)
          do index_nonzero=THETA_nonzero_1D_A_indexK(k-1)+1, &
            THETA_nonzero_1D_A_indexK(k)
             i=THETA_nonzero_1D_A(index_nonzero)%i
             j=THETA_nonzero_1D_A(index_nonzero)%j
             tempx=THETA_nonzero_1D_A(index_nonzero)%value
             Pr_A(i,k)=Pr_A(i,k)+ds*tempx*(3.0*qA(s-1,j,k)- &
                       3.0*qA(s-2,j,k)+qA(s-3,j,k))
            
           enddo
# endif /* MaierSaupe */
            do i=0,M_Bar_A-1
            Pr_A(i,k)=Pr_A(i,k)+3.0*qA(s-1,i,k)-1.5*qA(s-2,i,k)+ &
                      temp*qA(s-3,i,k) - &
                      ds*Wtemp*(3.0*qA(s-1,i,k)-3.0*qA(s-2,i,k)+qA(s-3,i,k))
            enddo
        enddo !enddo k=0,LOCAL_SIZE    
# if defined (OPENMP)
        !$omp end  do
        !$omp end  parallel
# endif  /* OPENMP */

!!!doing the FFT

      do i=0,M_Bar_A-1
        do k=0,LOCAL_SIZE-1
# if defined (Dim2)
          local_data(kD(k)%x,kD(k)%y)=cmplx(Pr_A(i,k),0.0)
# else  /* Dim2 */
          local_data(kD(k)%x,kD(k)%y,KD(k)%z)=cmplx(Pr_A(i,k),0.0)
# endif /* Dim2 */
         enddo
       
       call fftw_mpi_execute_dft(fplan, local_data, local_data)
       
        do k=0,LOCAL_SIZE-1
# if defined (Dim2)
            Pk_Real_A(i,k)=real(local_data(kD(k)%x,kD(k)%y))
            Pk_Imag_A(i,k)=AIMAG(local_data(kD(k)%x,kD(k)%y))
# else  /* Dim2 */
            Pk_Real_A(i,k)=real(local_data(kD(k)%x,kD(k)%y,KD(k)%z))
            Pk_Imag_A(i,k)=AIMAG(local_data(kD(k)%x,kD(k)%y,KD(k)%z))
# endif /* Dim2 */
         enddo
       enddo  !enddo i=0,M_Bar_A-1  


# if defined (OPENMP)
        !$omp parallel default(shared) private(k)
        !$omp do
# endif  /* OPENMP */
         do k=0,LOCAL_SIZE-1

         call mkl_dcsrsymv('l', M_Bar_A,  R_GA3_inverse_val(:,k+1), R_GA3_inverse_row(:,k+1), &
            R_GA3_inverse_col(:,k+1),Pk_Real_A(:,k), q_temp_RR_A(:,k))

         call mkl_dcsrsymv('l', M_Bar_A,  R_GA3_inverse_val(:,k+1), R_GA3_inverse_row(:,k+1), &
            R_GA3_inverse_col(:,k+1),Pk_Imag_A(:,k), q_temp_RI_A(:,k))

         call mkl_dcsrsymv('l', M_Bar_A,  I_GA3_inverse_val(:,k+1), I_GA3_inverse_row(:,k+1), &
            I_GA3_inverse_col(:,k+1),Pk_Real_A(:,k), q_temp_IR_A(:,k))

         call mkl_dcsrsymv('l', M_Bar_A,  I_GA3_inverse_val(:,k+1), I_GA3_inverse_row(:,k+1), &
            I_GA3_inverse_col(:,k+1),Pk_Imag_A(:,k), q_temp_II_A(:,k))
         enddo
# if defined (OPENMP)
        !$omp end  do
        !$omp end  parallel
# endif  /* OPENMP */




             do i=0,M_Bar_A-1
                do k=0,LOCAL_SIZE-1
# if defined (Dim2)
          local_data(kD(k)%x,kD(k)%y)=cmplx(q_temp_RR_A(i,k)-q_temp_II_A(i,k), & 
                q_temp_RI_A(i,k)+q_temp_IR_A(i,k))
# else  /* Dim2 */
          local_data(kD(k)%x,kD(k)%y,KD(k)%z)=cmplx(q_temp_RR_A(i,k)-q_temp_II_A(i,k), & 
                q_temp_RI_A(i,k)+q_temp_IR_A(i,k))
# endif /* Dim2 */
                enddo 
       call fftw_mpi_execute_dft(bplan, local_data, local_data)
                do k=0,LOCAL_SIZE-1
# if defined (Dim2)
                 qA(s,i,k)=real(local_data(kD(k)%x,kD(k)%y))*M_grid_inv
# else  /* Dim2 */
                 qA(s,i,k)=real(local_data(kD(k)%x,kD(k)%y,KD(k)%z))*M_grid_inv
# endif /* Dim2 */
                enddo
              enddo  


# if defined (Debug)
    call mp_barrier()
    if(myid==0 .and. s==NA) then
     write(*,*) "qA(NA,0,30)",qA(NA,0,30)
     endif
# endif /* Debug */
!!!!s>=3.and s>NA!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

else  ! else s>=3:s>NA
            !!!calculate polymer block B
             Pr_B=0.0
             ! else s>=3:s>NA:s==NA+1

   if(s==NA+1) then
             
# if defined (OPENMP)
        !$omp parallel default(shared) private(k,i,j,index_nonzero,Wtemp,tempx)
        !$omp do
# endif  /* OPENMP */
      do k=0,LOCAL_SIZE-1
          Wtemp=WB(k)

# if defined (MaierSaupe)
!!!!!!!!!!question : how to joint the two blocks with different M_Bar_A,M_Bar_B???
          do index_nonzero=THETA_nonzero_1D_B_indexK(k-1)+1, &
            THETA_nonzero_1D_B_indexK(k)
             i=THETA_nonzero_1D_B(index_nonzero)%i
             j=THETA_nonzero_1D_B(index_nonzero)%j
             tempx=THETA_nonzero_1D_B(index_nonzero)%value
             tempA1=0.0
             tempA2=0.0
             tempA3=0.0

             if(j<M_Bar_A) then
              tempA1=qA(NA,j,k)
              tempA2=qA(NA-1,j,k)
              tempA3=qA(NA-2,j,k)
              endif
              Pr_B(i,k)=Pr_B(i,k) + ds*tempx*(3.0*tempA1-3.0*tempA2+tempA3)
           enddo  !enddo index_nonzero
# endif /* MaierSaupe */

            do i=0,M_Bar_B-1
                  tempA1=0.0
                  tempA2=0.0
                  tempA3=0.0
                  if(i<M_Bar_A) then
                  tempA1=qA(NA,i,k)
                  tempA2=qA(NA-1,i,k)
                  tempA3=qA(NA-2,i,k)
                  endif
                 

            Pr_B(i,k)=Pr_B(i,k)+3.0*tempA1-1.5*tempA2+ &
                      temp*tempA3 - &
                      ds*Wtemp*(3.0*tempA1-3.0*tempA2+tempA3)
            enddo
          enddo !enddo k=0,LOCAL_SIZE    
# if defined (OPENMP)
        !$omp end  do
        !$omp end  parallel
# endif  /* OPENMP */

    else if(s==(NA+2)) then
         
# if defined (OPENMP)
        !$omp parallel default(shared) private(k,i,j,index_nonzero,Wtemp,tempx)
        !$omp do
# endif  /* OPENMP */
      do k=0,LOCAL_SIZE-1
          Wtemp=WB(k)
# if defined (MaierSaupe)
          do index_nonzero=THETA_nonzero_1D_B_indexK(k-1)+1, &
            THETA_nonzero_1D_B_indexK(k)
             i=THETA_nonzero_1D_B(index_nonzero)%i
             j=THETA_nonzero_1D_B(index_nonzero)%j
             tempx=THETA_nonzero_1D_B(index_nonzero)%value
             tempA2=0.0
             tempA3=0.0
             if(j<M_Bar_A) then
              tempA2=qA(NA,j,k)
              tempA3=qA(NA-1,j,k)
              endif
              Pr_B(i,k)=Pr_B(i,k) + ds*tempx*(3.0*qB(1,j,k)-3.0*tempA2+tempA3)
           enddo  !enddo index_nonzero
# endif /* MaierSaupe */

            do i=0,M_Bar_B-1
                  tempA2=0.0
                  tempA3=0.0
                  if(i<M_Bar_A) then
                  tempA2=qA(NA,i,k)
                  tempA3=qA(NA-1,i,k)
                  endif
                 

            Pr_B(i,k)=Pr_B(i,k)+3.0*qB(1,i,k)-1.5*tempA2+ &
                      temp*tempA3 - &
                      ds*Wtemp*(3.0*qB(1,i,k)-3.0*tempA2+tempA3)
            enddo
          enddo !enddo k=0,LOCAL_SIZE    
# if defined (OPENMP)
        !$omp end  do
        !$omp end  parallel
# endif  /* OPENMP */

     else if(s==(NA+3)) then
         
# if defined (OPENMP)
        !$omp parallel default(shared) private(k,i,j,index_nonzero,Wtemp,tempx)
        !$omp do
# endif  /* OPENMP */
      do k=0,LOCAL_SIZE-1
          Wtemp=WB(k)
# if defined (MaierSaupe)
          do index_nonzero=THETA_nonzero_1D_B_indexK(k-1)+1, &
            THETA_nonzero_1D_B_indexK(k)
             i=THETA_nonzero_1D_B(index_nonzero)%i
             j=THETA_nonzero_1D_B(index_nonzero)%j
             tempx=THETA_nonzero_1D_B(index_nonzero)%value
             tempA3=0.0
             if(j<M_Bar_A) then
              tempA3=qA(NA,j,k)
              endif
              Pr_B(i,k)=Pr_B(i,k) + ds*tempx*(3.0*qB(2,j,k)-3.0*qB(1,j,k)+tempA3)
           enddo  !enddo index_nonzero
# endif /* MaierSaupe */

            do i=0,M_Bar_B-1
                  tempA3=0.0
                  if(i<M_Bar_A) then
                  tempA3=qA(NA,i,k)
                  endif

            Pr_B(i,k)=Pr_B(i,k)+3.0*qB(2,i,k)-1.5*qB(1,i,k)+ &
                      temp*tempA3 - &
                      ds*Wtemp*(3.0*qB(2,i,k)-3.0*qB(1,i,k)+tempA3)
            enddo
          enddo !enddo k=0,LOCAL_SIZE    
# if defined (OPENMP)
        !$omp end  do
        !$omp end  parallel
# endif  /* OPENMP */

     else !else s>=3,s>NA:s>NA+3
           
# if defined (OPENMP)
        !$omp parallel default(shared) private(k,i,j,index_nonzero,Wtemp,tempx)
        !$omp do
# endif  /* OPENMP */
      do k=0,LOCAL_SIZE-1
          Wtemp=WB(k)
# if defined (MaierSaupe)
          do index_nonzero=THETA_nonzero_1D_B_indexK(k-1)+1, &
            THETA_nonzero_1D_B_indexK(k)
             i=THETA_nonzero_1D_B(index_nonzero)%i
             j=THETA_nonzero_1D_B(index_nonzero)%j
             tempx=THETA_nonzero_1D_B(index_nonzero)%value
              Pr_B(i,k)=Pr_B(i,k) + ds*tempx*(3.0*qB(s-NA-1,j,k)-3.0*qB(s-NA-2,j,k)+qB(s-NA-3,j,k))
           enddo  !enddo index_nonzero
# endif /* MaierSaupe */

            do i=0,M_Bar_B-1

            Pr_B(i,k)=Pr_B(i,k)+3.0*qB(s-NA-1,i,k)-1.5*qB(s-NA-2,i,k)+ &
                      temp*qB(s-NA-3,i,k) - &
                      ds*Wtemp*(3.0*qB(s-NA-1,i,k)-3.0*qB(s-NA-2,i,k)+qB(s-NA-3,i,k))
            enddo
          enddo !enddo k=0,LOCAL_SIZE    
# if defined (OPENMP)
        !$omp end  do
        !$omp end  parallel
# endif  /* OPENMP */
     endif  !if(s==NA+1) 

!!!doing the FFT

      do i=0,M_Bar_B-1
        do k=0,LOCAL_SIZE-1
# if defined (Dim2)
          local_data(kD(k)%x,kD(k)%y)=cmplx(Pr_B(i,k),0.0)
# else  /* Dim2 */
          local_data(kD(k)%x,kD(k)%y,KD(k)%z)=cmplx(Pr_B(i,k),0.0)
# endif /* Dim2 */
         enddo
       
       call fftw_mpi_execute_dft(fplan, local_data, local_data)
       
        do k=0,LOCAL_SIZE-1
# if defined (Dim2)
            Pk_Real_B(i,k)=real(local_data(kD(k)%x,kD(k)%y))
            Pk_Imag_B(i,k)=AIMAG(local_data(kD(k)%x,kD(k)%y))
# else  /* Dim2 */
            Pk_Real_B(i,k)=real(local_data(kD(k)%x,kD(k)%y,KD(k)%z))
            Pk_Imag_B(i,k)=AIMAG(local_data(kD(k)%x,kD(k)%y,KD(k)%z))
# endif /* Dim2 */
         enddo
       enddo  !enddo i=0,M_Bar_A-1  


# if defined (OPENMP)
        !$omp parallel default(shared) private(k)
        !$omp do
# endif  /* OPENMP */
         do k=0,LOCAL_SIZE-1

         call mkl_dcsrsymv('l', M_Bar_B,  R_GB3_inverse_val(:,k+1), R_GB3_inverse_row(:,k+1), &
            R_GB3_inverse_col(:,k+1),Pk_Real_B(:,k), q_temp_RR_B(:,k))

         call mkl_dcsrsymv('l', M_Bar_B,  R_GB3_inverse_val(:,k+1), R_GB3_inverse_row(:,k+1), &
            R_GB3_inverse_col(:,k+1),Pk_Imag_B(:,k), q_temp_RI_B(:,k))

         call mkl_dcsrsymv('l', M_Bar_B,  I_GB3_inverse_val(:,k+1), I_GB3_inverse_row(:,k+1), &
            I_GB3_inverse_col(:,k+1),Pk_Real_B(:,k), q_temp_IR_B(:,k))

         call mkl_dcsrsymv('l', M_Bar_B,  I_GB3_inverse_val(:,k+1), I_GB3_inverse_row(:,k+1), &
            I_GB3_inverse_col(:,k+1),Pk_Imag_B(:,k), q_temp_II_B(:,k))
         enddo
# if defined (OPENMP)
        !$omp end  do
        !$omp end  parallel
# endif  /* OPENMP */


             do i=0,M_Bar_B-1
                do k=0,LOCAL_SIZE-1
# if defined (Dim2)
        local_data(kD(k)%x,kD(k)%y)=cmplx(q_temp_RR_B(i,k)-q_temp_II_B(i,k), & 
                q_temp_RI_B(i,k)+q_temp_IR_B(i,k))
# else  /* Dim2 */
        local_data(kD(k)%x,kD(k)%y,KD(k)%z)=cmplx(q_temp_RR_B(i,k)-q_temp_II_B(i,k), & 
                q_temp_RI_B(i,k)+q_temp_IR_B(i,k))
# endif /* Dim2 */
                enddo 
       call fftw_mpi_execute_dft(bplan, local_data, local_data)
                do k=0,LOCAL_SIZE-1
# if defined (Dim2)
                 qB(s-NA,i,k)=real(local_data(kD(k)%x,kD(k)%y))*M_grid_inv
# else  /* Dim2 */
                 qB(s-NA,i,k)=real(local_data(kD(k)%x,kD(k)%y,KD(k)%z))*M_grid_inv
# endif /* Dim2 */
                enddo
              enddo 
 
  endif   !endif s>=3:s>NA
                 
endif    !endif s==1

enddo  !enddo s=1,NMAX
# if defined (Debug)
    call mp_barrier()
    if(myid==0 ) then
     write(*,*) "qB(NB,0,30)",qB(NB,0,30)
     endif
# endif /* Debug */
!!!!!!!!!!!!!!1


!!!normalization!!!
        do k=0,LOCAL_SIZE-1
         do i=0,M_Bar_B-1
         qB(0,i,k)=0.0
         enddo
       enddo

    do k=0,LOCAL_SIZE-1
       do i=0,min(M_Bar_A-1,M_Bar_B-1)
         qB(0,i,k)=qA(NA,i,k)
        enddo
     enddo

 deallocate(Pr_A)
 deallocate(Pk_Real_A)
 deallocate(Pk_Imag_A)
 deallocate(q_temp_RR_A)
 deallocate(q_temp_RI_A)
 deallocate(q_temp_IR_A)
 deallocate(q_temp_II_A)

 deallocate(Pr_B)
 deallocate(Pk_Real_B)
 deallocate(Pk_Imag_B)
 deallocate(q_temp_RR_B)
 deallocate(q_temp_RI_B)
 deallocate(q_temp_IR_B)
 deallocate(q_temp_II_B)


end subroutine MDE_q

      


subroutine MDE_qstar()
 USE nrtype,only :DP,PI
 USE global_para
 USE mpi
 USE control
 USE constants
 USE utility
 USE mmpi
 USE mpi_fftw3_operation
 USE SPH
 USE matrix_inverse
 implicit none
 integer :: i,j,s
 integer :: index_nonzero
 integer :: k,KK
 integer :: istat2
 integer,parameter :: namaxB=256 !this number is obtained by a prerun of input parameters
 integer,parameter :: namaxA=1800 !this number is obtained by a prerun of input parameters
 !for A
 real(DP) :: tempA1,tempA2,tempA3,temp,tempx,Wtemp,sumaq,sumaq_global
 real(DP),allocatable :: Pr_A(:,:)
 real(DP),allocatable :: Pk_Real_A(:,:)
 real(DP),allocatable :: Pk_Imag_A(:,:)
 real(DP),allocatable :: qstar_temp_RR_A(:,:)
 real(DP),allocatable :: qstar_temp_RI_A(:,:)
 real(DP),allocatable :: qstar_temp_IR_A(:,:)
 real(DP),allocatable :: qstar_temp_II_A(:,:)
 !for B
 real(DP),allocatable :: Pr_B(:,:)
 real(DP),allocatable :: Pk_Real_B(:,:)
 real(DP),allocatable :: Pk_Imag_B(:,:)
 real(DP),allocatable :: qstar_temp_RR_B(:,:)
 real(DP),allocatable :: qstar_temp_RI_B(:,:)
 real(DP),allocatable :: qstar_temp_IR_B(:,:)
 real(DP),allocatable :: qstar_temp_II_B(:,:)

 REAL(DP) :: M_grid_inv
 !!!! just for reference


!for A
 allocate(Pr_A(0:M_Bar_A-1,0:LOCAL_SIZE-1),stat=istat2)
 allocate(Pk_Real_A(0:M_Bar_A-1,0:LOCAL_SIZE-1),stat=istat2)
 allocate(Pk_Imag_A(0:M_Bar_A-1,0:LOCAL_SIZE-1),stat=istat2)


 allocate(qstar_temp_RR_A(0:M_Bar_A-1,0:LOCAL_SIZE-1),stat=istat2)
 allocate(qstar_temp_RI_A(0:M_Bar_A-1,0:LOCAL_SIZE-1),stat=istat2)
 allocate(qstar_temp_IR_A(0:M_Bar_A-1,0:LOCAL_SIZE-1),stat=istat2)
 allocate(qstar_temp_II_A(0:M_Bar_A-1,0:LOCAL_SIZE-1),stat=istat2)
!for B
 allocate(Pr_B(0:M_Bar_B-1,0:LOCAL_SIZE-1),stat=istat2)
 allocate(Pk_Real_B(0:M_Bar_B-1,0:LOCAL_SIZE-1),stat=istat2)
 allocate(Pk_Imag_B(0:M_Bar_B-1,0:LOCAL_SIZE-1),stat=istat2)

 allocate(qstar_temp_RR_B(0:M_Bar_B-1,0:LOCAL_SIZE-1),stat=istat2)
 allocate(qstar_temp_RI_B(0:M_Bar_B-1,0:LOCAL_SIZE-1),stat=istat2)
 allocate(qstar_temp_IR_B(0:M_Bar_B-1,0:LOCAL_SIZE-1),stat=istat2)
 allocate(qstar_temp_II_B(0:M_Bar_B-1,0:LOCAL_SIZE-1),stat=istat2)

 if(istat2/=0) then
 write(*,*) "allocate matrix in mde_q failed"
 stop
 endif

M_grid_inv=1.0d0/M_grid

do k=0,LOCAL_SIZE-1
   do i=0,M_Bar_B-1
    qBstar(NB,i,K)=0.0
    if(i==0) qBstar(NB,i,K)=1.0
   enddo
enddo

       
 !!!!!
 do s=Nmax-1,0,-1
 !!!s=1
   if(s==Nmax-1) then

          Pr_B=0.0
# if defined (OPENMP)
        !$omp parallel default(shared) private(k,i,j,index_nonzero,Wtemp,tempx)
        !$omp do
# endif  /* OPENMP */
      do k=0,LOCAL_SIZE-1
          Wtemp=WB(k)
# if defined (MaierSaupe)
          do index_nonzero=THETA_nonzero_1D_B_indexK(k-1)+1, &
            THETA_nonzero_1D_B_indexK(k)
             i=THETA_nonzero_1D_B(index_nonzero)%i
             j=THETA_nonzero_1D_B(index_nonzero)%j
             tempx=THETA_nonzero_1D_B(index_nonzero)%value
             Pr_B(i,k)=Pr_B(i,k)+ds*tempx*qBstar(NB,j,k)
           enddo
# endif /* MaierSaupe */

            do i=0,M_Bar_B-1
            Pr_B(i,k)=Pr_B(i,k)+(1.0-ds*Wtemp)*qBstar(NB,i,k)
            enddo
        enddo !enddo k=0,LOCAL_SIZE    
# if defined (OPENMP)
        !$omp end  do
        !$omp end  parallel
# endif  /* OPENMP */

!!!doing the FFT

      do i=0,M_Bar_B-1
        do k=0,LOCAL_SIZE-1
# if defined (Dim2)
          local_data(kD(k)%x,kD(k)%y)=cmplx(Pr_B(i,k),0.0)
# else  /* Dim2 */
          local_data(kD(k)%x,kD(k)%y,KD(k)%z)=cmplx(Pr_B(i,k),0.0)
# endif /* Dim2 */
         enddo
       
       
       call fftw_mpi_execute_dft(fplan, local_data, local_data)
       
        do k=0,LOCAL_SIZE-1
# if defined (Dim2)
            Pk_Real_B(i,k)=real(local_data(kD(k)%x,kD(k)%y))
            Pk_Imag_B(i,k)=AIMAG(local_data(kD(k)%x,kD(k)%y))
# else  /* Dim2 */
            Pk_Real_B(i,k)=real(local_data(kD(k)%x,kD(k)%y,KD(k)%z))
            Pk_Imag_B(i,k)=AIMAG(local_data(kD(k)%x,kD(k)%y,KD(k)%z))
# endif /* Dim2 */
         enddo
       enddo  !enddo i=0,M_Bar_B-1  


# if defined (OPENMP)
        !$omp parallel default(shared) private(k)
        !$omp do
# endif  /* OPENMP */

         do k=0,LOCAL_SIZE-1

         call mkl_dcsrsymv('l', M_Bar_B,  R_GB1_star_inverse_val(:,k+1), R_GB1_star_inverse_row(:,k+1), &
            R_GB1_star_inverse_col(:,k+1),Pk_Real_B(:,k), qstar_temp_RR_B(:,k))

         call mkl_dcsrsymv('l', M_Bar_B,  R_GB1_star_inverse_val(:,k+1), R_GB1_star_inverse_row(:,k+1), &
            R_GB1_star_inverse_col(:,k+1),Pk_Imag_B(:,k), qstar_temp_RI_B(:,k))

         call mkl_dcsrsymv('l', M_Bar_B,  I_GB1_star_inverse_val(:,k+1), I_GB1_star_inverse_row(:,k+1), &
            I_GB1_star_inverse_col(:,k+1),Pk_Real_B(:,k), qstar_temp_IR_B(:,k))

         call mkl_dcsrsymv('l', M_Bar_B,  I_GB1_star_inverse_val(:,k+1), I_GB1_star_inverse_row(:,k+1), &
            I_GB1_star_inverse_col(:,k+1),Pk_Imag_B(:,k), qstar_temp_II_B(:,k))


         enddo
# if defined (OPENMP)
        !$omp end  do
        !$omp end  parallel
# endif  /* OPENMP */


             do i=0,M_Bar_B-1
                do k=0,LOCAL_SIZE-1
# if defined (Dim2)
        local_data(kD(k)%x,kD(k)%y)=cmplx(qstar_temp_RR_B(i,k)-qstar_temp_II_B(i,k), & 
                qstar_temp_RI_B(i,k)+qstar_temp_IR_B(i,k))
# else  /* Dim2 */
        local_data(kD(k)%x,kD(k)%y,KD(k)%z)=cmplx(qstar_temp_RR_B(i,k)-qstar_temp_II_B(i,k), & 
                qstar_temp_RI_B(i,k)+qstar_temp_IR_B(i,k))
# endif /* Dim2 */
                enddo 
       call fftw_mpi_execute_dft(bplan, local_data, local_data)
                do k=0,LOCAL_SIZE-1
# if defined (Dim2)
                 qBstar(NB-1,i,k)=real(local_data(kD(k)%x,kD(k)%y))*M_grid_inv
# else  /* Dim2 */
                 qBstar(NB-1,i,k)=real(local_data(kD(k)%x,kD(k)%y,KD(k)%z))*M_grid_inv
# endif /* Dim2 */
                enddo
              enddo  

# if defined (Debug)
    call mp_barrier()
    if(myid==0 ) then
     write(*,*) "qBstar(NB-1,0,30)",qBstar(NB-1,0,30)
     endif
# endif /* Debug */
!!!!s=NMAX-2!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    else if(s==NMAX-2) then

          Pr_B=0.0
# if defined (OPENMP)
        !$omp parallel default(shared) private(k,i,j,index_nonzero,Wtemp,tempx)
        !$omp do
# endif  /* OPENMP */
      do k=0,LOCAL_SIZE-1
          Wtemp=WB(k)
# if defined (MaierSaupe)
          do index_nonzero=THETA_nonzero_1D_B_indexK(k-1)+1, &
            THETA_nonzero_1D_B_indexK(k)
             i=THETA_nonzero_1D_B(index_nonzero)%i
             j=THETA_nonzero_1D_B(index_nonzero)%j
             tempx=THETA_nonzero_1D_B(index_nonzero)%value
             Pr_B(i,k)=Pr_B(i,k)+ds*tempx*(2.0*qBstar(NB-1,j,k)-qBstar(NB,j,k))
           enddo
# endif /* MaierSaupe */

            do i=0,M_Bar_B-1
            Pr_B(i,k)=Pr_B(i,k)+2.0*qBstar(NB-1,i,k)-0.5*qBstar(NB,i,k) - &
                      ds*Wtemp*(2.0*qBstar(NB-1,i,k)-qBstar(NB,i,k))
            enddo
        enddo !enddo k=0,LOCAL_SIZE    
# if defined (OPENMP)
        !$omp end  do
        !$omp end  parallel
# endif  /* OPENMP */

!!!doing the FFT

      do i=0,M_Bar_B-1
        do k=0,LOCAL_SIZE-1
# if defined (Dim2)
          local_data(kD(k)%x,kD(k)%y)=cmplx(Pr_B(i,k),0.0)
# else  /* Dim2 */
          local_data(kD(k)%x,kD(k)%y,KD(k)%z)=cmplx(Pr_B(i,k),0.0)
# endif /* Dim2 */
         enddo
       
!fftw :
       
       call fftw_mpi_execute_dft(fplan, local_data, local_data)
       
        do k=0,LOCAL_SIZE-1
# if defined (Dim2)
            Pk_Real_B(i,k)=real(local_data(kD(k)%x,kD(k)%y))
            Pk_Imag_B(i,k)=AIMAG(local_data(kD(k)%x,kD(k)%y))
# else  /* Dim2 */
            Pk_Real_B(i,k)=real(local_data(kD(k)%x,kD(k)%y,KD(k)%z))
            Pk_Imag_B(i,k)=AIMAG(local_data(kD(k)%x,kD(k)%y,KD(k)%z))
# endif /* Dim2 */
         enddo
       enddo  !enddo i=0,M_Bar_B-1  


# if defined (OPENMP)
        !$omp parallel default(shared) private(k)
        !$omp do
# endif  /* OPENMP */
         do k=0,LOCAL_SIZE-1

         call mkl_dcsrsymv('l', M_Bar_B,  R_GB2_star_inverse_val(:,k+1), R_GB2_star_inverse_row(:,k+1), &
            R_GB2_star_inverse_col(:,k+1),Pk_Real_B(:,k), qstar_temp_RR_B(:,k))

         call mkl_dcsrsymv('l', M_Bar_B,  R_GB2_star_inverse_val(:,k+1), R_GB2_star_inverse_row(:,k+1), &
            R_GB2_star_inverse_col(:,k+1),Pk_Imag_B(:,k), qstar_temp_RI_B(:,k))

         call mkl_dcsrsymv('l', M_Bar_B,  I_GB2_star_inverse_val(:,k+1), I_GB2_star_inverse_row(:,k+1), &
            I_GB2_star_inverse_col(:,k+1),Pk_Real_B(:,k), qstar_temp_IR_B(:,k))

         call mkl_dcsrsymv('l', M_Bar_B,  I_GB2_star_inverse_val(:,k+1), I_GB2_star_inverse_row(:,k+1), &
            I_GB2_star_inverse_col(:,k+1),Pk_Imag_B(:,k), qstar_temp_II_B(:,k))
         enddo
# if defined (OPENMP)
        !$omp end  do
        !$omp end  parallel
# endif  /* OPENMP */

             do i=0,M_Bar_B-1
                do k=0,LOCAL_SIZE-1
# if defined (Dim2)
      local_data(kD(k)%x,kD(k)%y)=cmplx(qstar_temp_RR_B(i,k)-qstar_temp_II_B(i,k), & 
                qstar_temp_RI_B(i,k)+qstar_temp_IR_B(i,k))
# else /* Dim2 */
      local_data(kD(k)%x,kD(k)%y,KD(k)%z)=cmplx(qstar_temp_RR_B(i,k)-qstar_temp_II_B(i,k), & 
                qstar_temp_RI_B(i,k)+qstar_temp_IR_B(i,k))
# endif /* Dim2 */
                enddo 
                
       call fftw_mpi_execute_dft(bplan, local_data, local_data)
                do k=0,LOCAL_SIZE-1
# if defined (Dim2)
                 qBstar(NB-2,i,k)=real(local_data(kD(k)%x,kD(k)%y))*M_grid_inv
# else /* Dim2 */
                 qBstar(NB-2,i,k)=real(local_data(kD(k)%x,kD(k)%y,KD(k)%z))*M_grid_inv
# endif /* Dim2 */
                enddo
              enddo  

# if defined (Debug)
    call mp_barrier()
    if(myid==0 ) then
     write(*,*) "qBstar(NB-2,0,30)",qBstar(NB-2,0,30)
     endif
# endif /* Debug */
!!!!s<=NMAX-3!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    else if(s<=NMAX-3) then
          temp=1.0/3.0
           
          if(s>NA) then
            
          Pr_B=0.0
# if defined (OPENMP)
        !$omp parallel default(shared) private(k,i,j,index_nonzero,Wtemp,tempx)
        !$omp do
# endif  /* OPENMP */
      do k=0,LOCAL_SIZE-1
          Wtemp=WB(k)

# if defined (MaierSaupe)
          do index_nonzero=THETA_nonzero_1D_B_indexK(k-1)+1, &
            THETA_nonzero_1D_B_indexK(k)
             i=THETA_nonzero_1D_B(index_nonzero)%i
             j=THETA_nonzero_1D_B(index_nonzero)%j
             tempx=THETA_nonzero_1D_B(index_nonzero)%value
             Pr_B(i,k)=Pr_B(i,k)+ds*tempx*(3.0*qBstar(s-NA+1,j,k)- &
                       3.0*qBstar(s-NA+2,j,k)+qBstar(s-NA+3,j,k))
            
           enddo
# endif /* MaierSaupe */
            do i=0,M_Bar_B-1
            Pr_B(i,k)=Pr_B(i,k)+3.0*qBstar(s-NA+1,i,k)-1.5*qBstar(s-NA+2,i,k)+ &
                      temp*qBstar(s-NA+3,i,k) - &
                      ds*Wtemp*(3.0*qBstar(s-NA+1,i,k)-3.0*qBstar(s-NA+2,i,k)+qBstar(s-NA+3,i,k))
            enddo
        enddo !enddo k=0,LOCAL_SIZE    
# if defined (OPENMP)
        !$omp end  do
        !$omp end  parallel
# endif  /* OPENMP */

!!!doing the FFT

      do i=0,M_Bar_B-1
        do k=0,LOCAL_SIZE-1
# if defined (Dim2)
          local_data(kD(k)%x,kD(k)%y)=cmplx(Pr_B(i,k),0.0)
# else  /* Dim2 */
          local_data(kD(k)%x,kD(k)%y,KD(k)%z)=cmplx(Pr_B(i,k),0.0)
# endif /* Dim2 */
         enddo
       
!fftw :
       
       call fftw_mpi_execute_dft(fplan, local_data, local_data)
       
        do k=0,LOCAL_SIZE-1
# if defined (Dim2)
            Pk_Real_B(i,k)=real(local_data(kD(k)%x,kD(k)%y))
            Pk_Imag_B(i,k)=AIMAG(local_data(kD(k)%x,kD(k)%y))
# else  /* Dim2 */
            Pk_Real_B(i,k)=real(local_data(kD(k)%x,kD(k)%y,KD(k)%z))
            Pk_Imag_B(i,k)=AIMAG(local_data(kD(k)%x,kD(k)%y,KD(k)%z))
# endif /* Dim2 */
         enddo
       enddo  !enddo i=0,M_Bar_B-1  


# if defined (OPENMP)
        !$omp parallel default(shared) private(k)
        !$omp do
# endif  /* OPENMP */
         do k=0,LOCAL_SIZE-1

         call mkl_dcsrsymv('l', M_Bar_B,  R_GB3_star_inverse_val(:,k+1), R_GB3_star_inverse_row(:,k+1), &
            R_GB3_star_inverse_col(:,k+1),Pk_Real_B(:,k), qstar_temp_RR_B(:,k))

         call mkl_dcsrsymv('l', M_Bar_B,  R_GB3_star_inverse_val(:,k+1), R_GB3_star_inverse_row(:,k+1), &
            R_GB3_star_inverse_col(:,k+1),Pk_Imag_B(:,k), qstar_temp_RI_B(:,k))

         call mkl_dcsrsymv('l', M_Bar_B,  I_GB3_star_inverse_val(:,k+1), I_GB3_star_inverse_row(:,k+1), &
            I_GB3_star_inverse_col(:,k+1),Pk_Real_B(:,k), qstar_temp_IR_B(:,k))

         call mkl_dcsrsymv('l', M_Bar_B,  I_GB3_star_inverse_val(:,k+1), I_GB3_star_inverse_row(:,k+1), &
            I_GB3_star_inverse_col(:,k+1),Pk_Imag_B(:,k), qstar_temp_II_B(:,k))
         enddo
# if defined (OPENMP)
        !$omp end  do
        !$omp end  parallel
# endif  /* OPENMP */


             do i=0,M_Bar_B-1
                do k=0,LOCAL_SIZE-1
# if defined (Dim2)
    local_data(kD(k)%x,kD(k)%y)=cmplx(qstar_temp_RR_B(i,k)-qstar_temp_II_B(i,k), & 
                qstar_temp_RI_B(i,k)+qstar_temp_IR_B(i,k))
# else  /* Dim2 */
    local_data(kD(k)%x,kD(k)%y,KD(k)%z)=cmplx(qstar_temp_RR_B(i,k)-qstar_temp_II_B(i,k), & 
                qstar_temp_RI_B(i,k)+qstar_temp_IR_B(i,k))
# endif /* Dim2 */
                enddo 
       call fftw_mpi_execute_dft(bplan, local_data, local_data)
                do k=0,LOCAL_SIZE-1
# if defined (Dim2)
                 qBstar(s-NA,i,k)=real(local_data(kD(k)%x,kD(k)%y))*M_grid_inv
# else  /* Dim2 */
                 qBstar(s-NA,i,k)=real(local_data(kD(k)%x,kD(k)%y,KD(k)%z))*M_grid_inv
# endif /* Dim2 */
                enddo
              enddo  



   else  ! else s<=NA :entering A block regime
            !!!calculate polymer block A
             Pr_A=0.0
            
   if(s==NA) then
           Pr_B=0.0  
# if defined (OPENMP)
        !$omp parallel default(shared) private(k,i,j,index_nonzero,Wtemp,tempx)
        !$omp do
# endif  /* OPENMP */
      do k=0,LOCAL_SIZE-1
          Wtemp=WA(k)
# if defined (MaierSaupe)
          do index_nonzero=THETA_nonzero_1D_B_indexK(k-1)+1, &
            THETA_nonzero_1D_B_indexK(k)
             i=THETA_nonzero_1D_B(index_nonzero)%i
             j=THETA_nonzero_1D_B(index_nonzero)%j
             tempx=THETA_nonzero_1D_B(index_nonzero)%value

              Pr_B(i,k)=Pr_B(i,k) + ds*tempx*(3.0*qBstar(1,j,k)-3.0*qBstar(2,j,k)+qBstar(3,j,k))
           enddo  !enddo index_nonzero
# endif /* MaierSaupe */

            do i=0,min(M_Bar_A-1,M_Bar_B-1)
                Pr_A(i,k)=Pr_B(i,k)
                 

            Pr_A(i,k)=Pr_A(i,k)+3.0*qBstar(1,i,k)-1.5*qBstar(2,i,k)+ &
                      temp*qBstar(3,i,k) - &
                      ds*Wtemp*(3.0*qBstar(1,i,k)-3.0*qBstar(2,i,k)+qBstar(3,i,k))
            enddo
          enddo !enddo k=0,LOCAL_SIZE    
# if defined (OPENMP)
        !$omp end  do
        !$omp end  parallel
# endif  /* OPENMP */

    else if(s==(NA-1)) then
         
            Pr_B=0.0

# if defined (OPENMP)
        !$omp parallel default(shared) private(k,i,j,index_nonzero,Wtemp,tempx)
        !$omp do
# endif  /* OPENMP */
      do k=0,LOCAL_SIZE-1
          Wtemp=WA(k)
# if defined (MaierSaupe)
          do index_nonzero=THETA_nonzero_1D_B_indexK(k-1)+1, &
            THETA_nonzero_1D_B_indexK(k)
             i=THETA_nonzero_1D_B(index_nonzero)%i
             j=THETA_nonzero_1D_B(index_nonzero)%j
             tempx=THETA_nonzero_1D_B(index_nonzero)%value
             tempA1=0.0
             
             if(j<M_Bar_A) then
              tempA1=qAstar(NA,j,k)
          
              endif
              Pr_B(i,k)=Pr_B(i,k) + ds*tempx*(3.0*tempA1-3.0*qBstar(1,j,k)+qBstar(2,j,k))
           enddo  !enddo index_nonzero
# endif /* MaierSaupe */
             

            do i=0,min(M_Bar_A-1,M_Bar_B-1)
              Pr_A(i,k)=Pr_B(i,k)

            Pr_A(i,k)=Pr_A(i,k)+3.0*qAstar(NA,i,k)-1.5*qBstar(1,i,k)+ &
                      temp*qBstar(2,i,k) - &
                      ds*Wtemp*(3.0*qAstar(NA,i,k)-3.0*qBstar(1,i,k)+qBstar(2,i,k))
            enddo
          enddo !enddo k=0,LOCAL_SIZE    
# if defined (OPENMP)
        !$omp end  do
        !$omp end  parallel
# endif  /* OPENMP */

     else if(s==(NA-2)) then
         Pr_B=0.0
# if defined (OPENMP)
        !$omp parallel default(shared) private(k,i,j,index_nonzero,Wtemp,tempx)
        !$omp do
# endif  /* OPENMP */
      do k=0,LOCAL_SIZE-1
          Wtemp=WA(k)
# if defined (MaierSaupe)
          do index_nonzero=THETA_nonzero_1D_B_indexK(k-1)+1, &
            THETA_nonzero_1D_B_indexK(k)
             i=THETA_nonzero_1D_B(index_nonzero)%i
             j=THETA_nonzero_1D_B(index_nonzero)%j
             tempx=THETA_nonzero_1D_B(index_nonzero)%value
             tempA1=0.0
             tempA2=0.0
             if(j<M_Bar_A) then
              tempA1=qAstar(NA-1,j,k)
              tempA2=qAstar(NA,j,k)
              endif
              Pr_B(i,k)=Pr_B(i,k) + ds*tempx*(3.0*tempA1-3.0*tempA2+qBstar(1,j,k))
           enddo  !enddo index_nonzero
# endif /* MaierSaupe */

            do i=0,min(M_Bar_A-1,M_Bar_B-1)
                Pr_A(i,k)=Pr_B(i,k)
                 

            Pr_A(i,k)=Pr_A(i,k)+3.0*qAstar(NA-1,i,k)-1.5*qAstar(NA,i,k)+ &
                      temp*qBstar(1,i,k)- &
                      ds*Wtemp*(3.0*qAstar(NA-1,i,k)-3.0*qAstar(NA,i,k)+qBstar(1,i,k))
            enddo
          enddo !enddo k=0,LOCAL_SIZE    
# if defined (OPENMP)
        !$omp end  do
        !$omp end  parallel
# endif  /* OPENMP */

     else 
           
# if defined (OPENMP)
        !$omp parallel default(shared) private(k,i,j,index_nonzero,Wtemp,tempx)
        !$omp do
# endif  /* OPENMP */
      do k=0,LOCAL_SIZE-1
          Wtemp=WA(k)
# if defined (MaierSaupe)
          do index_nonzero=THETA_nonzero_1D_A_indexK(k-1)+1, &
            THETA_nonzero_1D_A_indexK(k)
             i=THETA_nonzero_1D_A(index_nonzero)%i
             j=THETA_nonzero_1D_A(index_nonzero)%j
             tempx=THETA_nonzero_1D_A(index_nonzero)%value
       Pr_A(i,k)=Pr_A(i,k) + ds*tempx*(3.0*qAstar(s+1,j,k)-3.0*qAstar(s+2,j,k)+qAstar(s+3,j,k))
           enddo  !enddo index_nonzero
# endif /* MaierSaupe */
            do i=0,M_Bar_A-1

            Pr_A(i,k)=Pr_A(i,k)+3.0*qAstar(s+1,i,k)-1.5*qAstar(s+2,i,k)+ &
                      temp*qAstar(s+3,i,k) - &
                      ds*Wtemp*(3.0*qAstar(s+1,i,k)-3.0*qAstar(s+2,i,k)+qAstar(s+3,i,k))
            enddo
          enddo !enddo k=0,LOCAL_SIZE    
# if defined (OPENMP)
        !$omp end  do
        !$omp end  parallel
# endif  /* OPENMP */
     endif   

!!!doing the FFT

      do i=0,M_Bar_A-1
        do k=0,LOCAL_SIZE-1
# if defined (Dim2)
          local_data(kD(k)%x,kD(k)%y)=cmplx(Pr_A(i,k),0.0)
# else  /* Dim2 */
          local_data(kD(k)%x,kD(k)%y,KD(k)%z)=cmplx(Pr_A(i,k),0.0)
# endif /* Dim2 */
         enddo
       
       
       call fftw_mpi_execute_dft(fplan, local_data, local_data)
       
        do k=0,LOCAL_SIZE-1
# if defined (Dim2)
            Pk_Real_A(i,k)=real(local_data(kD(k)%x,kD(k)%y))
            Pk_Imag_A(i,k)=AIMAG(local_data(kD(k)%x,kD(k)%y))
# else  /* Dim2 */
            Pk_Real_A(i,k)=real(local_data(kD(k)%x,kD(k)%y,KD(k)%z))
            Pk_Imag_A(i,k)=AIMAG(local_data(kD(k)%x,kD(k)%y,KD(k)%z))
# endif /* Dim2 */
         enddo
       enddo  !enddo i=0,M_Bar_A-1  


# if defined (OPENMP)
        !$omp parallel default(shared) private(k)
        !$omp do
# endif  /* OPENMP */
         do k=0,LOCAL_SIZE-1

         call mkl_dcsrsymv('l', M_Bar_A,  R_GA3_star_inverse_val(:,k+1), R_GA3_star_inverse_row(:,k+1), &
            R_GA3_star_inverse_col(:,k+1),Pk_Real_A(:,k), qstar_temp_RR_A(:,k))

         call mkl_dcsrsymv('l', M_Bar_A,  R_GA3_star_inverse_val(:,k+1), R_GA3_star_inverse_row(:,k+1), &
            R_GA3_star_inverse_col(:,k+1),Pk_Imag_A(:,k), qstar_temp_RI_A(:,k))

         call mkl_dcsrsymv('l', M_Bar_A,  I_GA3_star_inverse_val(:,k+1), I_GA3_star_inverse_row(:,k+1), &
            I_GA3_star_inverse_col(:,k+1),Pk_Real_A(:,k), qstar_temp_IR_A(:,k))

         call mkl_dcsrsymv('l', M_Bar_A,  I_GA3_star_inverse_val(:,k+1), I_GA3_star_inverse_row(:,k+1), &
            I_GA3_star_inverse_col(:,k+1),Pk_Imag_A(:,k), qstar_temp_II_A(:,k))
         enddo
# if defined (OPENMP)
        !$omp end  do
        !$omp end  parallel
# endif  /* OPENMP */


             do i=0,M_Bar_A-1
                do k=0,LOCAL_SIZE-1
# if defined (Dim2) 
       local_data(kD(k)%x,kD(k)%y)=cmplx(qstar_temp_RR_A(i,k)-qstar_temp_II_A(i,k), & 
                qstar_temp_RI_A(i,k)+qstar_temp_IR_A(i,k))
# else  /* Dim2 */
       local_data(kD(k)%x,kD(k)%y,KD(k)%z)=cmplx(qstar_temp_RR_A(i,k)-qstar_temp_II_A(i,k), & 
                qstar_temp_RI_A(i,k)+qstar_temp_IR_A(i,k))
# endif /* Dim2 */
                enddo 
       call fftw_mpi_execute_dft(bplan, local_data, local_data)
                do k=0,LOCAL_SIZE-1
# if defined (Dim2) 
                 qAstar(s,i,k)=real(local_data(kD(k)%x,kD(k)%y))*M_grid_inv
# else  /* Dim2 */
                 qAstar(s,i,k)=real(local_data(kD(k)%x,kD(k)%y,KD(k)%z))*M_grid_inv
# endif /* Dim2 */
                enddo
              enddo  
  endif   !endif 
                 
endif    !endif 

enddo  !enddo s=NMAX-1,0,-1

# if defined (Debug)
    call mp_barrier()
    if(myid==0 ) then
     write(*,*) "qAstar(0,0,30)",qAstar(0,0,30)
     endif
# endif /* Debug */
!

      do k=0,LOCAL_SIZE-1
       do i=0,M_Bar_B-1
         qBstar(0,i,k)=0.0
        enddo
     enddo

    do k=0,LOCAL_SIZE-1
       do i=0,min(M_Bar_A-1,M_Bar_B-1)
         qBstar(0,i,k)=qAstar(NA,i,k)
        enddo
     enddo

 deallocate(Pr_A)
 deallocate(Pk_Real_A)
 deallocate(Pk_Imag_A)
 deallocate(qstar_temp_RR_A)
 deallocate(qstar_temp_RI_A)
 deallocate(qstar_temp_IR_A)
 deallocate(qstar_temp_II_A)
!
 deallocate(Pr_B)
 deallocate(Pk_Real_B)
 deallocate(Pk_Imag_B)
 deallocate(qstar_temp_RR_B)
 deallocate(qstar_temp_RI_B)
 deallocate(qstar_temp_IR_B)
 deallocate(qstar_temp_II_B)

# if defined (Debug)
    call mp_barrier()
    if(myid==0 ) then
     write(*,*) "exit MDE_qstar"
     endif
# endif /* Debug */

end subroutine MDE_qstar

# endif  /* Mem_OP */




# if defined (Mem_OP)
subroutine MDE_q_final()
 USE nrtype,only :DP,PI
 USE global_para
 USE mpi
 USE control
 USE constants
 USE utility
 USE mmpi
 USE mpi_fftw3_operation
 USE SPH
 USE matrix_inverse
 USE G_matrix_mod
 implicit none
 integer :: i,j,s
 integer :: index_nonzero
 integer :: k
 integer :: istat2
 integer,parameter :: namax=1800 !this number is obtained by a prerun of input parameters
 !for A
 real(DP) :: tempA1,tempA2,tempA3,temp,tempx,Wtemp,sumaq,sumaq_global,sumaq1
 real(DP) :: sumPA,sumPA_global
 real(DP),allocatable :: Pr_A(:,:)
 real(DP),allocatable :: Pk_Real_A(:,:)
 real(DP),allocatable :: Pk_Imag_A(:,:)
 real(DP),allocatable :: q_temp_RR_A(:,:)
 real(DP),allocatable :: q_temp_RI_A(:,:)
 real(DP),allocatable :: q_temp_IR_A(:,:)
 real(DP),allocatable :: q_temp_II_A(:,:)
 !for B
 real(DP),allocatable :: Pr_B(:,:)
 real(DP),allocatable :: Pk_Real_B(:,:)
 real(DP),allocatable :: Pk_Imag_B(:,:)
 real(DP),allocatable :: q_temp_RR_B(:,:)
 real(DP),allocatable :: q_temp_RI_B(:,:)
 real(DP),allocatable :: q_temp_IR_B(:,:)
 real(DP),allocatable :: q_temp_II_B(:,:)
 REAL(DP) :: M_grid_inv

 REAL(DP) :: sy_values_R(1:namax)
 REAL(DP) :: sy_values_I(1:namax)
 integer ::  sy_columns_R(1:namax)
 integer ::  sy_columns_I(1:namax)
!!M_Bar_A >=M_BAr_B
 integer ::  sy_rowIndex_R(1:M_Bar_A+1)
 integer ::  sy_rowIndex_I(1:M_Bar_A+1)


 !deallocate G matrix first

if(allocated(R_GA3_inverse_val)) then
  deallocate(R_GA3_inverse_val)
endif
if(allocated(R_GA3_inverse_col)) then
  deallocate(R_GA3_inverse_col)
endif
if(allocated(R_GA3_inverse_row)) then
  deallocate(R_GA3_inverse_row)
endif
if(allocated(R_GA3_star_inverse_val)) then
  deallocate(R_GA3_star_inverse_val)
endif
if(allocated(R_GA3_star_inverse_col)) then
  deallocate(R_GA3_star_inverse_col)
endif
if(allocated(R_GA3_star_inverse_row)) then
  deallocate(R_GA3_star_inverse_row)
endif


if(allocated(I_GA3_inverse_val)) then
  deallocate(I_GA3_inverse_val)
endif
if(allocated(I_GA3_inverse_col)) then
  deallocate(I_GA3_inverse_col)
endif
if(allocated(I_GA3_inverse_row)) then
  deallocate(I_GA3_inverse_row)
endif
if(allocated(I_GA3_star_inverse_val)) then
  deallocate(I_GA3_star_inverse_val)
endif
if(allocated(I_GA3_star_inverse_col)) then
  deallocate(I_GA3_star_inverse_col)
endif
if(allocated(I_GA3_star_inverse_row)) then
  deallocate(I_GA3_star_inverse_row)
endif

  !deallocate(R_GB3_inverse_val)
  !deallocate(R_GB3_inverse_col)
  !deallocate(R_GB3_inverse_row)
  !deallocate(R_GB3_star_inverse_val)
  !deallocate(R_GB3_star_inverse_col)
  !deallocate(R_GB3_star_inverse_row)

  !deallocate(I_GB3_inverse_val)
  !deallocate(I_GB3_inverse_col)
  !deallocate(I_GB3_inverse_row)
  !deallocate(I_GB3_star_inverse_val)
  !deallocate(I_GB3_star_inverse_col)
  !deallocate(I_GB3_star_inverse_row)

 !!!! just for reference
! allocate(qA(0:NA,0:MA,0:LS),stat=istat2)

 !first ,allocate the matrixs
!for A
 allocate(Pr_A(0:M_Bar_A-1,0:LOCAL_SIZE-1),stat=istat2)
 allocate(Pk_Real_A(0:M_Bar_A-1,0:LOCAL_SIZE-1),stat=istat2)
 allocate(Pk_Imag_A(0:M_Bar_A-1,0:LOCAL_SIZE-1),stat=istat2)


 allocate(q_temp_RR_A(0:M_Bar_A-1,0:LOCAL_SIZE-1),stat=istat2)
 allocate(q_temp_RI_A(0:M_Bar_A-1,0:LOCAL_SIZE-1),stat=istat2)
 allocate(q_temp_IR_A(0:M_Bar_A-1,0:LOCAL_SIZE-1),stat=istat2)
 allocate(q_temp_II_A(0:M_Bar_A-1,0:LOCAL_SIZE-1),stat=istat2)
!for B
 allocate(Pr_B(0:M_Bar_B-1,0:LOCAL_SIZE-1),stat=istat2)
 allocate(Pk_Real_B(0:M_Bar_B-1,0:LOCAL_SIZE-1),stat=istat2)
 allocate(Pk_Imag_B(0:M_Bar_B-1,0:LOCAL_SIZE-1),stat=istat2)

 allocate(q_temp_RR_B(0:M_Bar_B-1,0:LOCAL_SIZE-1),stat=istat2)
 allocate(q_temp_RI_B(0:M_Bar_B-1,0:LOCAL_SIZE-1),stat=istat2)
 allocate(q_temp_IR_B(0:M_Bar_B-1,0:LOCAL_SIZE-1),stat=istat2)
 allocate(q_temp_II_B(0:M_Bar_B-1,0:LOCAL_SIZE-1),stat=istat2)
  

 if(istat2/=0) then
 write(*,*) "allocate matrix in mde_q failed"
 stop
 endif


 M_grid_inv=1.0d0/M_grid


do k=0,LOCAL_SIZE-1
    qA(0,0,K)=1.0
   do i=1,M_Bar_A-1
    qA(0,i,K)=0.0
   enddo
enddo

!if(myid==0) then
!write(*,*) "WA(0,30,300) in MDE=",WA(0),WA(30),WA(300)
!write(*,*) "WB(0,30,300)=",WB(0),WB(30),WB(300)
!endif  
!if(myid==0) then
!write(*,*) "M_OP(0)(30) 00=",M_OP(0,0,0),M_OP(30,0,0)
!write(*,*) "M_OP(0)(30) 11=",M_OP(0,1,1),M_OP(30,1,1)
!endif  
       
 !!!!!
 do s=1,Nmax

 !!!s=1
   if(s==1) then

          Pr_A=0.0

# if defined (OPENMP)
        !$omp parallel default(shared) private(k,i,j,index_nonzero,Wtemp,tempx)
        !$omp do
# endif  /* OPENMP */
     
      do k=0,LOCAL_SIZE-1
          Wtemp=WA(k)
# if defined (MaierSaupe)
          do index_nonzero=THETA_nonzero_1D_A_indexK(k-1)+1, &
            THETA_nonzero_1D_A_indexK(k)
             i=THETA_nonzero_1D_A(index_nonzero)%i
             j=THETA_nonzero_1D_A(index_nonzero)%j
             tempx=THETA_nonzero_1D_A(index_nonzero)%value
             Pr_A(i,k)=Pr_A(i,k)+ds*tempx*qA(0,j,k)
           enddo
# endif /* MaierSaupe */
            do i=0,M_Bar_A-1
            Pr_A(i,k)=Pr_A(i,k)+(1.0-ds*Wtemp)*qA(0,i,k)
            enddo
        enddo !enddo k=0,LOCAL_SIZE    
# if defined (OPENMP)
        !$omp end  do
        !$omp end  parallel
# endif  /* OPENMP */





!!!doing the FFT

      do i=0,M_Bar_A-1
        do k=0,LOCAL_SIZE-1
          local_data(k+1)=cmplx(Pr_A(i,k),0.0)
         enddo
       
       
       call fftw_mpi_execute_dft(fplan, local_data, local_data)
       
        do k=0,LOCAL_SIZE-1
            Pk_Real_A(i,k)=real(local_data(k+1))
            Pk_Imag_A(i,k)=AIMAG(local_data(k+1))
            !Pk_cmplx_A(i,k)=local_data(k)
         enddo
       enddo  !enddo i=0,M_Bar_A-1  



# if defined (OPENMP)
        !$omp parallel default(shared) private(k)
        !$omp do
# endif  /* OPENMP */
         do k=0,LOCAL_SIZE-1

          call  init_GA_inverse(1,k,sy_values_R, &
                sy_values_I,sy_columns_R,sy_columns_I, &
                sy_rowIndex_R,sy_rowIndex_I )
            

         call mkl_dcsrsymv('l', M_Bar_A,  sy_values_R, sy_rowIndex_R, &
            sy_columns_R,Pk_Real_A(:,k), q_temp_RR_A(:,k))

         call mkl_dcsrsymv('l', M_Bar_A,  sy_values_R, sy_rowIndex_R, &
            sy_columns_R,Pk_Imag_A(:,k), q_temp_RI_A(:,k))

         call mkl_dcsrsymv('l', M_Bar_A,  sy_values_I, sy_rowIndex_I, &
            sy_columns_I,Pk_Real_A(:,k), q_temp_IR_A(:,k))

         call mkl_dcsrsymv('l', M_Bar_A,  sy_values_I, sy_rowIndex_I, &
            sy_columns_I,Pk_Imag_A(:,k), q_temp_II_A(:,k))

            

         enddo

# if defined (OPENMP)
        !$omp end  do
        !$omp end  parallel
# endif  /* OPENMP */

             do i=0,M_Bar_A-1
                do k=0,LOCAL_SIZE-1
                local_data(k+1)=cmplx(q_temp_RR_A(i,k)-q_temp_II_A(i,k), & 
                q_temp_RI_A(i,k)+q_temp_IR_A(i,k))
                enddo
       call fftw_mpi_execute_dft(bplan, local_data, local_data)

            
                do k=0,LOCAL_SIZE-1
                 qA(1,i,k)=real(local_data(k+1))*M_grid_inv
                enddo
              enddo  


!!!!s=2!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    else if(s==2) then

          Pr_A=0.0
# if defined (OPENMP)
        !$omp parallel default(shared) private(k,i,j,index_nonzero,Wtemp,tempx)
        !$omp do
# endif  /* OPENMP */
      do k=0,LOCAL_SIZE-1
          Wtemp=WA(k)
# if defined (MaierSaupe)
          do index_nonzero=THETA_nonzero_1D_A_indexK(k-1)+1, &
            THETA_nonzero_1D_A_indexK(k)
             i=THETA_nonzero_1D_A(index_nonzero)%i
             j=THETA_nonzero_1D_A(index_nonzero)%j
             tempx=THETA_nonzero_1D_A(index_nonzero)%value
             Pr_A(i,k)=Pr_A(i,k)+ds*tempx*(2.0*qA(1,j,k)-qA(0,j,k))
           enddo
# endif /* MaierSaupe */
            do i=0,M_Bar_A-1
            Pr_A(i,k)=Pr_A(i,k)+2.0*qA(1,i,k)-0.5*qA(0,i,k) - &
                      ds*Wtemp*(2.0*qA(1,i,k)-qA(0,i,k))
            enddo
        enddo !enddo k=0,LOCAL_SIZE    
# if defined (OPENMP)
        !$omp end  do
        !$omp end  parallel
# endif  /* OPENMP */

!!!doing the FFT

      do i=0,M_Bar_A-1
        do k=0,LOCAL_SIZE-1
          local_data(k+1)=cmplx(Pr_A(i,k),0.0)
         enddo
       
       
       call fftw_mpi_execute_dft(fplan, local_data, local_data)
       
        do k=0,LOCAL_SIZE-1
            Pk_Real_A(i,k)=real(local_data(k+1))
            Pk_Imag_A(i,k)=AIMAG(local_data(k+1))
         enddo
       enddo  !enddo i=0,M_Bar_A-1  


# if defined (OPENMP)
        !$omp parallel default(shared) private(k)
        !$omp do
# endif  /* OPENMP */
         do k=0,LOCAL_SIZE-1
          
          call  init_GA_inverse(2,k,sy_values_R, &
                sy_values_I,sy_columns_R,sy_columns_I, &
                sy_rowIndex_R,sy_rowIndex_I )
 
         call mkl_dcsrsymv('l', M_Bar_A,  sy_values_R, sy_rowIndex_R, &
            sy_columns_R,Pk_Real_A(:,k), q_temp_RR_A(:,k))

         call mkl_dcsrsymv('l', M_Bar_A,  sy_values_R, sy_rowIndex_R, &
            sy_columns_R,Pk_Imag_A(:,k), q_temp_RI_A(:,k))

         call mkl_dcsrsymv('l', M_Bar_A,  sy_values_I, sy_rowIndex_I, &
            sy_columns_I,Pk_Real_A(:,k), q_temp_IR_A(:,k))

         call mkl_dcsrsymv('l', M_Bar_A,  sy_values_I, sy_rowIndex_I, &
            sy_columns_I,Pk_Imag_A(:,k), q_temp_II_A(:,k))

         enddo
# if defined (OPENMP)
        !$omp end  do
        !$omp end  parallel
# endif  /* OPENMP */


             do i=0,M_Bar_A-1
                do k=0,LOCAL_SIZE-1
                local_data(k+1)=cmplx(q_temp_RR_A(i,k)-q_temp_II_A(i,k), & 
                q_temp_RI_A(i,k)+q_temp_IR_A(i,k))
                enddo 
       call fftw_mpi_execute_dft(bplan, local_data, local_data)
                do k=0,LOCAL_SIZE-1
                 qA(2,i,k)=real(local_data(k+1))*M_grid_inv
                enddo
              enddo  

!!!!s>=3!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    else if(s>=3) then
          temp=1.0/3.0
           
          if(s<=NA) then
            
          Pr_A=0.0
# if defined (OPENMP)
        !$omp parallel default(shared) private(k,i,j,index_nonzero,Wtemp,tempx)
        !$omp do
# endif  /* OPENMP */
      do k=0,LOCAL_SIZE-1
          Wtemp=WA(k)
# if defined (MaierSaupe)
          do index_nonzero=THETA_nonzero_1D_A_indexK(k-1)+1, &
            THETA_nonzero_1D_A_indexK(k)
             i=THETA_nonzero_1D_A(index_nonzero)%i
             j=THETA_nonzero_1D_A(index_nonzero)%j
             tempx=THETA_nonzero_1D_A(index_nonzero)%value
             Pr_A(i,k)=Pr_A(i,k)+ds*tempx*(3.0*qA(s-1,j,k)- &
                       3.0*qA(s-2,j,k)+qA(s-3,j,k))
            
           enddo
# endif /* MaierSaupe */          
            do i=0,M_Bar_A-1
            Pr_A(i,k)=Pr_A(i,k)+3.0*qA(s-1,i,k)-1.5*qA(s-2,i,k)+ &
                      temp*qA(s-3,i,k) - &
                      ds*Wtemp*(3.0*qA(s-1,i,k)-3.0*qA(s-2,i,k)+qA(s-3,i,k))
            enddo
        enddo !enddo k=0,LOCAL_SIZE    
# if defined (OPENMP)
        !$omp end  do
        !$omp end  parallel
# endif  /* OPENMP */

!!!doing the FFT

      do i=0,M_Bar_A-1
        do k=0,LOCAL_SIZE-1
          local_data(k+1)=cmplx(Pr_A(i,k),0.0)
         enddo
       
       
       call fftw_mpi_execute_dft(fplan, local_data, local_data)
       
        do k=0,LOCAL_SIZE-1
            Pk_Real_A(i,k)=real(local_data(k+1))
            Pk_Imag_A(i,k)=AIMAG(local_data(k+1))
         enddo
       enddo  !enddo i=0,M_Bar_A-1  


# if defined (OPENMP)
        !$omp parallel default(shared) private(k)
        !$omp do
# endif  /* OPENMP */
         do k=0,LOCAL_SIZE-1

          call  init_GA_inverse(3,k,sy_values_R, &
                sy_values_I,sy_columns_R,sy_columns_I, &
                sy_rowIndex_R,sy_rowIndex_I )
 
         call mkl_dcsrsymv('l', M_Bar_A,  sy_values_R, sy_rowIndex_R, &
            sy_columns_R,Pk_Real_A(:,k), q_temp_RR_A(:,k))

         call mkl_dcsrsymv('l', M_Bar_A,  sy_values_R, sy_rowIndex_R, &
            sy_columns_R,Pk_Imag_A(:,k), q_temp_RI_A(:,k))

         call mkl_dcsrsymv('l', M_Bar_A,  sy_values_I, sy_rowIndex_I, &
            sy_columns_I,Pk_Real_A(:,k), q_temp_IR_A(:,k))

         call mkl_dcsrsymv('l', M_Bar_A,  sy_values_I, sy_rowIndex_I, &
            sy_columns_I,Pk_Imag_A(:,k), q_temp_II_A(:,k))
         enddo
# if defined (OPENMP)
        !$omp end  do
        !$omp end  parallel
# endif  /* OPENMP */




             do i=0,M_Bar_A-1
                do k=0,LOCAL_SIZE-1
                local_data(k+1)=cmplx(q_temp_RR_A(i,k)-q_temp_II_A(i,k), & 
                q_temp_RI_A(i,k)+q_temp_IR_A(i,k))
                enddo 
       call fftw_mpi_execute_dft(bplan, local_data, local_data)
                do k=0,LOCAL_SIZE-1
                 qA(s,i,k)=real(local_data(k+1))*M_grid_inv
                enddo
              enddo  

!!!!s>=3.and s>NA!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

else  ! else s>=3:s>NA
            !!!calculate polymer block B
             Pr_B=0.0
             ! else s>=3:s>NA:s==NA+1

   if(s==NA+1) then
             
# if defined (OPENMP)
        !$omp parallel default(shared) private(k,i,j,index_nonzero,Wtemp,tempx)
        !$omp do
# endif  /* OPENMP */
      do k=0,LOCAL_SIZE-1
          Wtemp=WB(k)
# if defined (MaierSaupe)
!!!!!!!!!!question : how to joint the two blocks with different M_Bar_A,M_Bar_B???
          do index_nonzero=THETA_nonzero_1D_B_indexK(k-1)+1, &
            THETA_nonzero_1D_B_indexK(k)
             i=THETA_nonzero_1D_B(index_nonzero)%i
             j=THETA_nonzero_1D_B(index_nonzero)%j
             tempx=THETA_nonzero_1D_B(index_nonzero)%value
             tempA1=0.0
             tempA2=0.0
             tempA3=0.0

             if(j<M_Bar_A) then
              tempA1=qA(NA,j,k)
              tempA2=qA(NA-1,j,k)
              tempA3=qA(NA-2,j,k)
              endif
              Pr_B(i,k)=Pr_B(i,k) + ds*tempx*(3.0*tempA1-3.0*tempA2+tempA3)
           enddo  !enddo index_nonzero
# endif /* MaierSaupe */


            do i=0,M_Bar_B-1
                  tempA1=0.0
                  tempA2=0.0
                  tempA3=0.0
                  if(i<M_Bar_A) then
                  tempA1=qA(NA,i,k)
                  tempA2=qA(NA-1,i,k)
                  tempA3=qA(NA-2,i,k)
                  endif
                 

            Pr_B(i,k)=Pr_B(i,k)+3.0*tempA1-1.5*tempA2+ &
                      temp*tempA3 - &
                      ds*Wtemp*(3.0*tempA1-3.0*tempA2+tempA3)
            enddo
          enddo !enddo k=0,LOCAL_SIZE    
# if defined (OPENMP)
        !$omp end  do
        !$omp end  parallel
# endif  /* OPENMP */

    else if(s==(NA+2)) then
         
# if defined (OPENMP)
        !$omp parallel default(shared) private(k,i,j,index_nonzero,Wtemp,tempx)
        !$omp do
# endif  /* OPENMP */
      do k=0,LOCAL_SIZE-1
          Wtemp=WB(k)
# if defined (MaierSaupe)
          do index_nonzero=THETA_nonzero_1D_B_indexK(k-1)+1, &
            THETA_nonzero_1D_B_indexK(k)
             i=THETA_nonzero_1D_B(index_nonzero)%i
             j=THETA_nonzero_1D_B(index_nonzero)%j
             tempx=THETA_nonzero_1D_B(index_nonzero)%value
             tempA2=0.0
             tempA3=0.0
             if(j<M_Bar_A) then
              tempA2=qA(NA,j,k)
              tempA3=qA(NA-1,j,k)
              endif
              Pr_B(i,k)=Pr_B(i,k) + ds*tempx*(3.0*qB(1,j,k)-3.0*tempA2+tempA3)
           enddo  !enddo index_nonzero
# endif /* MaierSaupe */
            do i=0,M_Bar_B-1
                  tempA2=0.0
                  tempA3=0.0
                  if(i<M_Bar_A) then
                  tempA2=qA(NA,i,k)
                  tempA3=qA(NA-1,i,k)
                  endif
                 

            Pr_B(i,k)=Pr_B(i,k)+3.0*qB(1,i,k)-1.5*tempA2+ &
                      temp*tempA3 - &
                      ds*Wtemp*(3.0*qB(1,i,k)-3.0*tempA2+tempA3)
            enddo
          enddo !enddo k=0,LOCAL_SIZE    
# if defined (OPENMP)
        !$omp end  do
        !$omp end  parallel
# endif  /* OPENMP */

     else if(s==(NA+3)) then
         
# if defined (OPENMP)
        !$omp parallel default(shared) private(k,i,j,index_nonzero,Wtemp,tempx)
        !$omp do
# endif  /* OPENMP */
      do k=0,LOCAL_SIZE-1
          Wtemp=WB(k)
# if defined (MaierSaupe)
          do index_nonzero=THETA_nonzero_1D_B_indexK(k-1)+1, &
            THETA_nonzero_1D_B_indexK(k)
             i=THETA_nonzero_1D_B(index_nonzero)%i
             j=THETA_nonzero_1D_B(index_nonzero)%j
             tempx=THETA_nonzero_1D_B(index_nonzero)%value
             tempA3=0.0
             if(j<M_Bar_A) then
              tempA3=qA(NA,j,k)
              endif
              Pr_B(i,k)=Pr_B(i,k) + ds*tempx*(3.0*qB(2,j,k)-3.0*qB(1,j,k)+tempA3)
           enddo  !enddo index_nonzero
# endif /* MaierSaupe */
            do i=0,M_Bar_B-1
                  tempA3=0.0
                  if(i<M_Bar_A) then
                  tempA3=qA(NA,i,k)
                  endif
                 

            Pr_B(i,k)=Pr_B(i,k)+3.0*qB(2,i,k)-1.5*qB(1,i,k)+ &
                      temp*tempA3 - &
                      ds*Wtemp*(3.0*qB(2,i,k)-3.0*qB(1,i,k)+tempA3)
            enddo
          enddo !enddo k=0,LOCAL_SIZE    
# if defined (OPENMP)
        !$omp end  do
        !$omp end  parallel
# endif  /* OPENMP */

     else !else s>=3,s>NA:s>NA+3
           
# if defined (OPENMP)
        !$omp parallel default(shared) private(k,i,j,index_nonzero,Wtemp,tempx)
        !$omp do
# endif  /* OPENMP */
      do k=0,LOCAL_SIZE-1
          Wtemp=WB(k)
# if defined (MaierSaupe)
          do index_nonzero=THETA_nonzero_1D_B_indexK(k-1)+1, &
            THETA_nonzero_1D_B_indexK(k)
             i=THETA_nonzero_1D_B(index_nonzero)%i
             j=THETA_nonzero_1D_B(index_nonzero)%j
             tempx=THETA_nonzero_1D_B(index_nonzero)%value
              Pr_B(i,k)=Pr_B(i,k) + ds*tempx*(3.0*qB(s-NA-1,j,k)-3.0*qB(s-NA-2,j,k)+qB(s-NA-3,j,k))
           enddo  !enddo index_nonzero
# endif /* MaierSaupe */
            do i=0,M_Bar_B-1

            Pr_B(i,k)=Pr_B(i,k)+3.0*qB(s-NA-1,i,k)-1.5*qB(s-NA-2,i,k)+ &
                      temp*qB(s-NA-3,i,k) - &
                      ds*Wtemp*(3.0*qB(s-NA-1,i,k)-3.0*qB(s-NA-2,i,k)+qB(s-NA-3,i,k))
            enddo
          enddo !enddo k=0,LOCAL_SIZE    
# if defined (OPENMP)
        !$omp end  do
        !$omp end  parallel
# endif  /* OPENMP */
     endif  !if(s==NA+1) 

!!!doing the FFT

      do i=0,M_Bar_B-1
        do k=0,LOCAL_SIZE-1
          local_data(k+1)=cmplx(Pr_B(i,k),0.0)
         enddo
       
       
       call fftw_mpi_execute_dft(fplan, local_data, local_data)
       
        do k=0,LOCAL_SIZE-1
            Pk_Real_B(i,k)=real(local_data(k+1))
            Pk_Imag_B(i,k)=AIMAG(local_data(k+1))
            !Pk_cmplx_B(i,k)=local_data(k)
         enddo
       enddo  !enddo i=0,M_Bar_A-1  


# if defined (OPENMP)
        !$omp parallel default(shared) private(k)
        !$omp do
# endif  /* OPENMP */
         do k=0,LOCAL_SIZE-1


          call  init_GA_inverse(5,k,sy_values_R, &
                sy_values_I,sy_columns_R,sy_columns_I, &
                sy_rowIndex_R,sy_rowIndex_I )
 
         call mkl_dcsrsymv('l', M_Bar_B,  sy_values_R, sy_rowIndex_R, &
            sy_columns_R,Pk_Real_B(:,k), q_temp_RR_B(:,k))

         call mkl_dcsrsymv('l', M_Bar_B,  sy_values_R, sy_rowIndex_R, &
            sy_columns_R,Pk_Imag_B(:,k), q_temp_RI_B(:,k))

         call mkl_dcsrsymv('l', M_Bar_B,  sy_values_I, sy_rowIndex_I, &
            sy_columns_I,Pk_Real_B(:,k), q_temp_IR_B(:,k))

         call mkl_dcsrsymv('l', M_Bar_B,  sy_values_I, sy_rowIndex_I, &
            sy_columns_I,Pk_Imag_B(:,k), q_temp_II_B(:,k))
         enddo
# if defined (OPENMP)
        !$omp end  do
        !$omp end  parallel
# endif  /* OPENMP */


             do i=0,M_Bar_B-1
                do k=0,LOCAL_SIZE-1
                local_data(k+1)=cmplx(q_temp_RR_B(i,k)-q_temp_II_B(i,k), & 
                q_temp_RI_B(i,k)+q_temp_IR_B(i,k))
                enddo 
       call fftw_mpi_execute_dft(bplan, local_data, local_data)
       !   call fftwnd_f77_mpi(plancb,1,local_data,work,0,FFTW_NORMAL_ORDER)
                do k=0,LOCAL_SIZE-1
                 qB(s-NA,i,k)=real(local_data(k+1))*M_grid_inv
                enddo
              enddo 
 
  endif   !endif s>=3:s>NA
                 
endif    !endif s==1

enddo  !enddo s=1,NMAX
!!!!!!!!!!!!!!1


!!!normalization!!!
!!warning,in case M_Bar_B /= M_Bar_A, qB(0),qA(NA) must be handled very
!carefully!
        do k=0,LOCAL_SIZE-1
         do i=0,M_Bar_B-1
         qB(0,i,k)=0.0
         qBend(i,k)=qB(NB,i,k)
         enddo
       enddo

    do k=0,LOCAL_SIZE-1
       do i=0,min(M_Bar_A-1,M_Bar_B-1)
         qB(0,i,k)=qA(NA,i,k)
        enddo
     enddo



 deallocate(Pr_A)
 deallocate(Pk_Real_A)
 deallocate(Pk_Imag_A)
 deallocate(q_temp_RR_A)
 deallocate(q_temp_RI_A)
 deallocate(q_temp_IR_A)
 deallocate(q_temp_II_A)

 deallocate(Pr_B)
 deallocate(Pk_Real_B)
 deallocate(Pk_Imag_B)
 deallocate(q_temp_RR_B)
 deallocate(q_temp_RI_B)
 deallocate(q_temp_IR_B)
 deallocate(q_temp_II_B)


end subroutine MDE_q_final

      


subroutine MDE_qstar_final()
 USE nrtype,only :DP,PI
 USE global_para
 USE mpi
 USE control
 USE constants
 USE utility
 USE mmpi
 USE mpi_fftw3_operation
 USE SPH
 USE G_matrix_mod
 USE matrix_inverse
 implicit none
 integer :: i,j,s
 integer :: index_nonzero
 integer :: k,KK,MA,MB,LS
 integer :: istat2
 integer,parameter :: namax=1800 !this number is obtained by a prerun of input parameters
 !for A
 real(DP) :: tempA1,tempA2,tempA3,temp,tempx,Wtemp,sumaq,sumaq_global
 real(DP),allocatable :: Pr_A(:,:)
 real(DP),allocatable :: Pk_Real_A(:,:)
 real(DP),allocatable :: Pk_Imag_A(:,:)
 real(DP),allocatable :: qstar_temp_RR_A(:,:)
 real(DP),allocatable :: qstar_temp_RI_A(:,:)
 real(DP),allocatable :: qstar_temp_IR_A(:,:)
 real(DP),allocatable :: qstar_temp_II_A(:,:)
 !for B
 real(DP),allocatable :: Pr_B(:,:)
 real(DP),allocatable :: Pk_Real_B(:,:)
 real(DP),allocatable :: Pk_Imag_B(:,:)
 real(DP),allocatable :: qstar_temp_RR_B(:,:)
 real(DP),allocatable :: qstar_temp_RI_B(:,:)
 real(DP),allocatable :: qstar_temp_IR_B(:,:)
 real(DP),allocatable :: qstar_temp_II_B(:,:)

 REAL(DP) :: M_grid_inv
 !!!! just for reference
 REAL(DP) :: sy_values_R(1:namax)
 REAL(DP) :: sy_values_I(1:namax)
 integer ::  sy_columns_R(1:namax)
 integer ::  sy_columns_I(1:namax)
!!M_Bar_A >=M_BAr_B
 integer ::  sy_rowIndex_R(1:M_Bar_A+1)
 integer ::  sy_rowIndex_I(1:M_Bar_A+1)
  


 MB=M_Bar_B-1
 MA=M_Bar_A-1
 LS=LOCAL_SIZE-1



 allocate(qAstar(0:NA,0:MA,0:LS),stat=istat2)
 allocate(qBstar(0:NB,0:MB,0:LS),stat=istat2)

!for A
 allocate(Pr_A(0:M_Bar_A-1,0:LOCAL_SIZE-1),stat=istat2)
 allocate(Pk_Real_A(0:M_Bar_A-1,0:LOCAL_SIZE-1),stat=istat2)
 allocate(Pk_Imag_A(0:M_Bar_A-1,0:LOCAL_SIZE-1),stat=istat2)


 allocate(qstar_temp_RR_A(0:M_Bar_A-1,0:LOCAL_SIZE-1),stat=istat2)
 allocate(qstar_temp_RI_A(0:M_Bar_A-1,0:LOCAL_SIZE-1),stat=istat2)
 allocate(qstar_temp_IR_A(0:M_Bar_A-1,0:LOCAL_SIZE-1),stat=istat2)
 allocate(qstar_temp_II_A(0:M_Bar_A-1,0:LOCAL_SIZE-1),stat=istat2)
!for B
 allocate(Pr_B(0:M_Bar_B-1,0:LOCAL_SIZE-1),stat=istat2)
 allocate(Pk_Real_B(0:M_Bar_B-1,0:LOCAL_SIZE-1),stat=istat2)
 allocate(Pk_Imag_B(0:M_Bar_B-1,0:LOCAL_SIZE-1),stat=istat2)

 allocate(qstar_temp_RR_B(0:M_Bar_B-1,0:LOCAL_SIZE-1),stat=istat2)
 allocate(qstar_temp_RI_B(0:M_Bar_B-1,0:LOCAL_SIZE-1),stat=istat2)
 allocate(qstar_temp_IR_B(0:M_Bar_B-1,0:LOCAL_SIZE-1),stat=istat2)
 allocate(qstar_temp_II_B(0:M_Bar_B-1,0:LOCAL_SIZE-1),stat=istat2)

 if(istat2/=0) then
 write(*,*) "allocate matrix in mde_q failed"
 stop
 endif

M_grid_inv=1.0d0/M_grid


do k=0,LOCAL_SIZE-1
   do i=0,M_Bar_B-1
    qBstar(NB,i,K)=0.0
    if(i==0) qBstar(NB,i,K)=1.0
   enddo
enddo
       
 !!!!!
 do s=Nmax-1,0,-1
 !!!s=1
   if(s==Nmax-1) then

          Pr_B=0.0
# if defined (OPENMP)
        !$omp parallel default(shared) private(k,i,j,index_nonzero,Wtemp,tempx)
        !$omp do
# endif  /* OPENMP */
      do k=0,LOCAL_SIZE-1
          Wtemp=WB(k)
# if defined (MaierSaupe)
          do index_nonzero=THETA_nonzero_1D_B_indexK(k-1)+1, &
            THETA_nonzero_1D_B_indexK(k)
             i=THETA_nonzero_1D_B(index_nonzero)%i
             j=THETA_nonzero_1D_B(index_nonzero)%j
             tempx=THETA_nonzero_1D_B(index_nonzero)%value
             Pr_B(i,k)=Pr_B(i,k)+ds*tempx*qBstar(NB,j,k)
           enddo
# endif /* MaierSaupe */
            do i=0,M_Bar_B-1
            Pr_B(i,k)=Pr_B(i,k)+(1.0-ds*Wtemp)*qBstar(NB,i,k)
            enddo
        enddo !enddo k=0,LOCAL_SIZE    
# if defined (OPENMP)
        !$omp end  do
        !$omp end  parallel
# endif  /* OPENMP */

!!!doing the FFT

      do i=0,M_Bar_B-1
        do k=0,LOCAL_SIZE-1
          local_data(k+1)=cmplx(Pr_B(i,k),0.0)
         enddo
       
       
       call fftw_mpi_execute_dft(fplan, local_data, local_data)
       
        do k=0,LOCAL_SIZE-1
            Pk_Real_B(i,k)=real(local_data(k+1))
            Pk_Imag_B(i,k)=AIMAG(local_data(k+1))
         enddo
       enddo  !enddo i=0,M_Bar_B-1  


# if defined (OPENMP)
        !$omp parallel default(shared) private(k)
        !$omp do
# endif  /* OPENMP */

         do k=0,LOCAL_SIZE-1

          call  init_GA_inverse(8,k,sy_values_R, &
                sy_values_I,sy_columns_R,sy_columns_I, &
                sy_rowIndex_R,sy_rowIndex_I )
 
         call mkl_dcsrsymv('l', M_Bar_B,  sy_values_R, sy_rowIndex_R, &
            sy_columns_R,Pk_Real_B(:,k), qstar_temp_RR_B(:,k))

         call mkl_dcsrsymv('l', M_Bar_B,  sy_values_R, sy_rowIndex_R, &
            sy_columns_R,Pk_Imag_B(:,k), qstar_temp_RI_B(:,k))

         call mkl_dcsrsymv('l', M_Bar_B,  sy_values_I, sy_rowIndex_I, &
            sy_columns_I,Pk_Real_B(:,k), qstar_temp_IR_B(:,k))

         call mkl_dcsrsymv('l', M_Bar_B,  sy_values_I, sy_rowIndex_I, &
            sy_columns_I,Pk_Imag_B(:,k), qstar_temp_II_B(:,k))



         enddo
# if defined (OPENMP)
        !$omp end  do
        !$omp end  parallel
# endif  /* OPENMP */


             do i=0,M_Bar_B-1
                do k=0,LOCAL_SIZE-1
                local_data(k+1)=cmplx(qstar_temp_RR_B(i,k)-qstar_temp_II_B(i,k), & 
                qstar_temp_RI_B(i,k)+qstar_temp_IR_B(i,k))

                enddo 
       call fftw_mpi_execute_dft(bplan, local_data, local_data)
       !   call fftwnd_f77_mpi(plancb,1,local_data,work,0,FFTW_NORMAL_ORDER)
                do k=0,LOCAL_SIZE-1
             qBstar(NB-1,i,k)=real(local_data(k+1))*M_grid_inv
                enddo
              enddo  

!!!!s=NMAX-2!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    else if(s==NMAX-2) then

          Pr_B=0.0
# if defined (OPENMP)
        !$omp parallel default(shared) private(k,i,j,index_nonzero,Wtemp,tempx)
        !$omp do
# endif  /* OPENMP */
      do k=0,LOCAL_SIZE-1
          Wtemp=WB(k)
# if defined (MaierSaupe)
          do index_nonzero=THETA_nonzero_1D_B_indexK(k-1)+1, &
            THETA_nonzero_1D_B_indexK(k)
             i=THETA_nonzero_1D_B(index_nonzero)%i
             j=THETA_nonzero_1D_B(index_nonzero)%j
             tempx=THETA_nonzero_1D_B(index_nonzero)%value
             Pr_B(i,k)=Pr_B(i,k)+ds*tempx*(2.0*qBstar(NB-1,j,k)-qBstar(NB,j,k))
           enddo
# endif /* MaierSaupe */
            do i=0,M_Bar_B-1
            Pr_B(i,k)=Pr_B(i,k)+2.0*qBstar(NB-1,i,k)-0.5*qBstar(NB,i,k) - &
                      ds*Wtemp*(2.0*qBstar(NB-1,i,k)-qBstar(NB,i,k))
            !Pr_B(i,k)=Pr_B(i,k)+2.0*qstar_temp_B1(i,k)-0.5*qstar_temp_B2(i,k) - &
            !          ds*Wtemp*(2.0*qstar_temp_B1(i,k)-qstar_temp_B2(i,k))
            enddo
        enddo !enddo k=0,LOCAL_SIZE    
# if defined (OPENMP)
        !$omp end  do
        !$omp end  parallel
# endif  /* OPENMP */

!!!doing the FFT

      do i=0,M_Bar_B-1
        do k=0,LOCAL_SIZE-1
          local_data(k+1)=cmplx(Pr_B(i,k),0.0)
         enddo
       
       
       call fftw_mpi_execute_dft(fplan, local_data, local_data)
       
        do k=0,LOCAL_SIZE-1
            Pk_Real_B(i,k)=real(local_data(k+1))
            Pk_Imag_B(i,k)=AIMAG(local_data(k+1))
         enddo
       enddo  !enddo i=0,M_Bar_B-1  


# if defined (OPENMP)
        !$omp parallel default(shared) private(k)
        !$omp do
# endif  /* OPENMP */
         do k=0,LOCAL_SIZE-1

          call  init_GA_inverse(7,k,sy_values_R, &
                sy_values_I,sy_columns_R,sy_columns_I, &
                sy_rowIndex_R,sy_rowIndex_I )
 
         call mkl_dcsrsymv('l', M_Bar_B,  sy_values_R, sy_rowIndex_R, &
            sy_columns_R,Pk_Real_B(:,k), qstar_temp_RR_B(:,k))

         call mkl_dcsrsymv('l', M_Bar_B,  sy_values_R, sy_rowIndex_R, &
            sy_columns_R,Pk_Imag_B(:,k), qstar_temp_RI_B(:,k))

         call mkl_dcsrsymv('l', M_Bar_B,  sy_values_I, sy_rowIndex_I, &
            sy_columns_I,Pk_Real_B(:,k), qstar_temp_IR_B(:,k))

         call mkl_dcsrsymv('l', M_Bar_B,  sy_values_I, sy_rowIndex_I, &
            sy_columns_I,Pk_Imag_B(:,k), qstar_temp_II_B(:,k))

         enddo
# if defined (OPENMP)
        !$omp end  do
        !$omp end  parallel
# endif  /* OPENMP */

             do i=0,M_Bar_B-1
                do k=0,LOCAL_SIZE-1
                local_data(k+1)=cmplx(qstar_temp_RR_B(i,k)-qstar_temp_II_B(i,k), & 
                qstar_temp_RI_B(i,k)+qstar_temp_IR_B(i,k))
                enddo 
                
       call fftw_mpi_execute_dft(bplan, local_data, local_data)
                do k=0,LOCAL_SIZE-1
             qBstar(NB-2,i,k)=real(local_data(k+1))*M_grid_inv
                enddo
              enddo  

!!!!s<=NMAX-3!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    else if(s<=NMAX-3) then
          temp=1.0/3.0
           
          if(s>NA) then
            
          Pr_B=0.0
# if defined (OPENMP)
        !$omp parallel default(shared) private(k,i,j,index_nonzero,Wtemp,tempx)
        !$omp do
# endif  /* OPENMP */
      do k=0,LOCAL_SIZE-1
          Wtemp=WB(k)
# if defined (MaierSaupe)
          do index_nonzero=THETA_nonzero_1D_B_indexK(k-1)+1, &
            THETA_nonzero_1D_B_indexK(k)
             i=THETA_nonzero_1D_B(index_nonzero)%i
             j=THETA_nonzero_1D_B(index_nonzero)%j
             tempx=THETA_nonzero_1D_B(index_nonzero)%value
             Pr_B(i,k)=Pr_B(i,k)+ds*tempx*(3.0*qBstar(s-NA+1,j,k)- &
                       3.0*qBstar(s-NA+2,j,k)+qBstar(s-NA+3,j,k))
            
           enddo
# endif /* MaierSaupe */
            do i=0,M_Bar_B-1
            Pr_B(i,k)=Pr_B(i,k)+3.0*qBstar(s-NA+1,i,k)-1.5*qBstar(s-NA+2,i,k)+ &
                      temp*qBstar(s-NA+3,i,k) - &
             ds*Wtemp*(3.0*qBstar(s-NA+1,i,k)-3.0*qBstar(s-NA+2,i,k)+qBstar(s-NA+3,i,k))
           ! Pr_B(i,k)=Pr_B(i,k)+3.0*qstar_temp_B1(i,k)-1.5*qstar_temp_B2(i,k)+ &
           !           temp*qstar_temp_B3(i,k) - &
           !  ds*Wtemp*(3.0*qstar_temp_B1(i,k)-3.0*qstar_temp_B2(i,k)+qstar_temp_B3(i,k))
            enddo
        enddo !enddo k=0,LOCAL_SIZE    
# if defined (OPENMP)
        !$omp end  do
        !$omp end  parallel
# endif  /* OPENMP */

!!!doing the FFT

      do i=0,M_Bar_B-1
        do k=0,LOCAL_SIZE-1
          local_data(k+1)=cmplx(Pr_B(i,k),0.0)
         enddo
       
       call fftw_mpi_execute_dft(fplan, local_data, local_data)
       
        do k=0,LOCAL_SIZE-1
            Pk_Real_B(i,k)=real(local_data(k+1))
            Pk_Imag_B(i,k)=AIMAG(local_data(k+1))
         enddo
       enddo  !enddo i=0,M_Bar_B-1  


# if defined (OPENMP)
        !$omp parallel default(shared) private(k)
        !$omp do
# endif  /* OPENMP */
         do k=0,LOCAL_SIZE-1

          call  init_GA_inverse(6,k,sy_values_R, &
                sy_values_I,sy_columns_R,sy_columns_I, &
                sy_rowIndex_R,sy_rowIndex_I )
 
         call mkl_dcsrsymv('l', M_Bar_B,  sy_values_R, sy_rowIndex_R, &
            sy_columns_R,Pk_Real_B(:,k), qstar_temp_RR_B(:,k))

         call mkl_dcsrsymv('l', M_Bar_B,  sy_values_R, sy_rowIndex_R, &
            sy_columns_R,Pk_Imag_B(:,k), qstar_temp_RI_B(:,k))

         call mkl_dcsrsymv('l', M_Bar_B,  sy_values_I, sy_rowIndex_I, &
            sy_columns_I,Pk_Real_B(:,k), qstar_temp_IR_B(:,k))

         call mkl_dcsrsymv('l', M_Bar_B,  sy_values_I, sy_rowIndex_I, &
            sy_columns_I,Pk_Imag_B(:,k), qstar_temp_II_B(:,k))

         enddo
# if defined (OPENMP)
        !$omp end  do
        !$omp end  parallel
# endif  /* OPENMP */


             do i=0,M_Bar_B-1
                do k=0,LOCAL_SIZE-1
                local_data(k+1)=cmplx(qstar_temp_RR_B(i,k)-qstar_temp_II_B(i,k), & 
                qstar_temp_RI_B(i,k)+qstar_temp_IR_B(i,k))
                enddo 
       call fftw_mpi_execute_dft(bplan, local_data, local_data)
       !   call fftwnd_f77_mpi(plancb,1,local_data,work,0,FFTW_NORMAL_ORDER)
                do k=0,LOCAL_SIZE-1
                  qBstar(s-NA,i,k)=real(local_data(k+1))*M_grid_inv
                enddo
              enddo  



   else  ! else s<=NA :entering A block regime
            !!!calculate polymer block A
             Pr_A=0.0
            
   if(s==NA) then
           Pr_B=0.0  
# if defined (OPENMP)
        !$omp parallel default(shared) private(k,i,j,index_nonzero,Wtemp,tempx)
        !$omp do
# endif  /* OPENMP */
      do k=0,LOCAL_SIZE-1
          Wtemp=WA(k)
# if defined (MaierSaupe)
          do index_nonzero=THETA_nonzero_1D_B_indexK(k-1)+1, &
            THETA_nonzero_1D_B_indexK(k)
             i=THETA_nonzero_1D_B(index_nonzero)%i
             j=THETA_nonzero_1D_B(index_nonzero)%j
             tempx=THETA_nonzero_1D_B(index_nonzero)%value

              Pr_B(i,k)=Pr_B(i,k) + ds*tempx*(3.0*qBstar(1,j,k)-3.0*qBstar(2,j,k)+qBstar(3,j,k))
           enddo  !enddo index_nonzero
# endif /* MaierSaupe */
            do i=0,min(M_Bar_A-1,M_Bar_B-1)
                Pr_A(i,k)=Pr_B(i,k)
                 

            Pr_A(i,k)=Pr_A(i,k)+3.0*qBstar(1,i,k)-1.5*qBstar(2,i,k)+ &
                      temp*qBstar(3,i,k) - &
                      ds*Wtemp*(3.0*qBstar(1,i,k)-3.0*qBstar(2,i,k)+qBstar(3,i,k))
            enddo
          enddo !enddo k=0,LOCAL_SIZE    
# if defined (OPENMP)
        !$omp end  do
        !$omp end  parallel
# endif  /* OPENMP */

    else if(s==(NA-1)) then
         
            Pr_B=0.0

# if defined (OPENMP)
        !$omp parallel default(shared) private(k,i,j,index_nonzero,Wtemp,tempx)
        !$omp do
# endif  /* OPENMP */
      do k=0,LOCAL_SIZE-1
          Wtemp=WA(k)
# if defined (MaierSaupe)
          do index_nonzero=THETA_nonzero_1D_B_indexK(k-1)+1, &
            THETA_nonzero_1D_B_indexK(k)
             i=THETA_nonzero_1D_B(index_nonzero)%i
             j=THETA_nonzero_1D_B(index_nonzero)%j
             tempx=THETA_nonzero_1D_B(index_nonzero)%value
             tempA1=0.0
             
             if(j<M_Bar_A) then
              tempA1=qAstar(NA,j,k)
          
              endif
              Pr_B(i,k)=Pr_B(i,k) + ds*tempx*(3.0*tempA1-3.0*qBstar(1,j,k)+qBstar(2,j,k))
           enddo  !enddo index_nonzero
# endif /* MaierSaupe */

            do i=0,min(M_Bar_A-1,M_Bar_B-1)
              Pr_A(i,k)=Pr_B(i,k)

            Pr_A(i,k)=Pr_A(i,k)+3.0*qAstar(NA,i,k)-1.5*qBstar(1,i,k)+ &
                      temp*qBstar(2,i,k) - &
                      ds*Wtemp*(3.0*qAstar(NA,i,k)-3.0*qBstar(1,i,k)+qBstar(2,i,k))
            !Pr_A(i,k)=Pr_A(i,k)+3.0*qstar_temp_B1(i,k)-1.5*qstar_temp_B2(i,k)+ &
            !          temp*qstar_temp_B3(i,k) - &
            !          ds*Wtemp*(3.0*qstar_temp_B1(i,k)-3.0*qstar_temp_B2(i,k)+qstar_temp_B3(i,k))
            enddo
          enddo !enddo k=0,LOCAL_SIZE    
# if defined (OPENMP)
        !$omp end  do
        !$omp end  parallel
# endif  /* OPENMP */

     else if(s==(NA-2)) then
         Pr_B=0.0
# if defined (OPENMP)
        !$omp parallel default(shared) private(k,i,j,index_nonzero,Wtemp,tempx)
        !$omp do
# endif  /* OPENMP */
      do k=0,LOCAL_SIZE-1
          Wtemp=WA(k)
# if defined (MaierSaupe)
          do index_nonzero=THETA_nonzero_1D_B_indexK(k-1)+1, &
            THETA_nonzero_1D_B_indexK(k)
             i=THETA_nonzero_1D_B(index_nonzero)%i
             j=THETA_nonzero_1D_B(index_nonzero)%j
             tempx=THETA_nonzero_1D_B(index_nonzero)%value
             tempA1=0.0
             tempA2=0.0
             if(j<M_Bar_A) then
              tempA1=qAstar(NA-1,j,k)
              tempA2=qAstar(NA,j,k)
              endif
              Pr_B(i,k)=Pr_B(i,k) + ds*tempx*(3.0*tempA1-3.0*tempA2+qBstar(1,j,k))
           enddo  !enddo index_nonzero
# endif /* MaierSaupe */

            do i=0,min(M_Bar_A-1,M_Bar_B-1)
                Pr_A(i,k)=Pr_B(i,k)
                 

            Pr_A(i,k)=Pr_A(i,k)+3.0*qAstar(NA-1,i,k)-1.5*qAstar(NA,i,k)+ &
                      temp*qBstar(1,i,k)- &
                      ds*Wtemp*(3.0*qAstar(NA-1,i,k)-3.0*qAstar(NA,i,k)+qBstar(1,i,k))
            !Pr_A(i,k)=Pr_A(i,k)+3.0*qstar_temp_B1(i,k)-1.5*qstar_temp_B2(i,k)+ &
            !          temp*qstar_temp_B3(i,k) - &
            !          ds*Wtemp*(3.0*qstar_temp_B1(i,k)-3.0*qstar_temp_B2(i,k)+qstar_temp_B3(i,k))
            enddo
          enddo !enddo k=0,LOCAL_SIZE    
# if defined (OPENMP)
        !$omp end  do
        !$omp end  parallel
# endif  /* OPENMP */

     else 
           
# if defined (OPENMP)
        !$omp parallel default(shared) private(k,i,j,index_nonzero,Wtemp,tempx)
        !$omp do
# endif  /* OPENMP */
      do k=0,LOCAL_SIZE-1
          Wtemp=WA(k)
# if defined (MaierSaupe)
          do index_nonzero=THETA_nonzero_1D_A_indexK(k-1)+1, &
            THETA_nonzero_1D_A_indexK(k)
             i=THETA_nonzero_1D_A(index_nonzero)%i
             j=THETA_nonzero_1D_A(index_nonzero)%j
             tempx=THETA_nonzero_1D_A(index_nonzero)%value
       Pr_A(i,k)=Pr_A(i,k) + ds*tempx*(3.0*qAstar(s+1,j,k)-3.0*qAstar(s+2,j,k)+qAstar(s+3,j,k))
           enddo  !enddo index_nonzero
# endif /* MaierSaupe */


            do i=0,M_Bar_A-1

            Pr_A(i,k)=Pr_A(i,k)+3.0*qAstar(s+1,i,k)-1.5*qAstar(s+2,i,k)+ &
                      temp*qAstar(s+3,i,k) - &
                      ds*Wtemp*(3.0*qAstar(s+1,i,k)-3.0*qAstar(s+2,i,k)+qAstar(s+3,i,k))
            !Pr_A(i,k)=Pr_A(i,k)+3.0*qstar_temp_B1(i,k)-1.5*qstar_temp_B2(i,k)+ &
            !          temp*qstar_temp_B3(i,k) - &
            !          ds*Wtemp*(3.0*qstar_temp_B1(i,k)-3.0*qstar_temp_B2(i,k)+qstar_temp_B3(i,k))
            enddo
          enddo !enddo k=0,LOCAL_SIZE    
# if defined (OPENMP)
        !$omp end  do
        !$omp end  parallel
# endif  /* OPENMP */
     endif   

!!!doing the FFT

      do i=0,M_Bar_A-1
        do k=0,LOCAL_SIZE-1
          local_data(k+1)=cmplx(Pr_A(i,k),0.0)
         enddo
       
       
       call fftw_mpi_execute_dft(fplan, local_data, local_data)
       
        do k=0,LOCAL_SIZE-1
            Pk_Real_A(i,k)=real(local_data(k+1))
            Pk_Imag_A(i,k)=AIMAG(local_data(k+1))
         enddo
       enddo  !enddo i=0,M_Bar_A-1  


# if defined (OPENMP)
        !$omp parallel default(shared) private(k)
        !$omp do
# endif  /* OPENMP */
         do k=0,LOCAL_SIZE-1


          call  init_GA_inverse(4,k,sy_values_R, &
                sy_values_I,sy_columns_R,sy_columns_I, &
                sy_rowIndex_R,sy_rowIndex_I )
 
         call mkl_dcsrsymv('l', M_Bar_A,  sy_values_R, sy_rowIndex_R, &
            sy_columns_R,Pk_Real_A(:,k), qstar_temp_RR_A(:,k))

         call mkl_dcsrsymv('l', M_Bar_A,  sy_values_R, sy_rowIndex_R, &
            sy_columns_R,Pk_Imag_A(:,k), qstar_temp_RI_A(:,k))

         call mkl_dcsrsymv('l', M_Bar_A,  sy_values_I, sy_rowIndex_I, &
            sy_columns_I,Pk_Real_A(:,k), qstar_temp_IR_A(:,k))

         call mkl_dcsrsymv('l', M_Bar_A,  sy_values_I, sy_rowIndex_I, &
            sy_columns_I,Pk_Imag_A(:,k), qstar_temp_II_A(:,k))

         enddo
# if defined (OPENMP)
        !$omp end  do
        !$omp end  parallel
# endif  /* OPENMP */


             do i=0,M_Bar_A-1
                do k=0,LOCAL_SIZE-1
                local_data(k+1)=cmplx(qstar_temp_RR_A(i,k)-qstar_temp_II_A(i,k), & 
                qstar_temp_RI_A(i,k)+qstar_temp_IR_A(i,k))
                enddo 
       call fftw_mpi_execute_dft(bplan, local_data, local_data)
                do k=0,LOCAL_SIZE-1
                  qAstar(s,i,k)=real(local_data(k+1))*M_grid_inv
                enddo
              enddo  
  endif   !endif 
                 
endif    !endif 

enddo  !enddo s=NMAX-1,0,-1

!

      do k=0,LOCAL_SIZE-1
       do i=0,M_Bar_B-1
         qBstar(0,i,k)=0.0
        enddo
     enddo

    do k=0,LOCAL_SIZE-1
       do i=0,min(M_Bar_A-1,M_Bar_B-1)
         qBstar(0,i,k)=qAstar(NA,i,k)
        enddo
     enddo

 deallocate(Pr_A)
 deallocate(Pk_Real_A)
 deallocate(Pk_Imag_A)
 deallocate(qstar_temp_RR_A)
 deallocate(qstar_temp_RI_A)
 deallocate(qstar_temp_IR_A)
 deallocate(qstar_temp_II_A)
!
 deallocate(Pr_B)
 deallocate(Pk_Real_B)
 deallocate(Pk_Imag_B)
 deallocate(qstar_temp_RR_B)
 deallocate(qstar_temp_RI_B)
 deallocate(qstar_temp_IR_B)
 deallocate(qstar_temp_II_B)
if(myid==0) then
write(*,*) "done MDE_final"
endif

end subroutine MDE_qstar_final

# endif /* Mem_OP */

end module  wormlike_MDE 






