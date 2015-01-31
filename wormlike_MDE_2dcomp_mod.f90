# if defined (SPH)
# else /* SPH */
# if defined (SLAB)
# else /* SLAB */
       subroutine calc_THETAij()
       use global_para
       USE control
       USE constants
       USE utility
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
# endif /* SLAB */

# if defined (SLAB)
# else /* SLAB */
Module wormlike_MDE
implicit none
contains

 subroutine MDE_q(R2C_FFT)
 USE nrtype,only :DP,PI
 USE global_para
 USE mpi
 USE control
 USE constants
 USE utility
 USE mmpi
 use decomp_fft_mod
 USE SPH
 USE matrix_inverse
 !USE G_matrix_2dcomp_mod
 implicit none
 logical, optional :: R2C_FFT
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
 integer :: Local_sizek

# if defined (R2C)
LOCAL_SIZEK=LOCAL_SIZE_K_R2C
# else /* R2C */
LOCAL_SIZEK=LOCAL_SIZE_K
# endif /* R2C */




 !first ,allocate the matrixs
!for A
 allocate(Pr_A(0:M_Bar_A-1,0:LOCAL_SIZE-1),stat=istat2)
 allocate(Pk_Real_A(0:M_Bar_A-1,0:LOCAL_SIZEK-1),stat=istat2)
 allocate(Pk_Imag_A(0:M_Bar_A-1,0:LOCAL_SIZEK-1),stat=istat2)

 allocate(q_temp_RR_A(0:M_Bar_A-1,0:LOCAL_SIZEK-1),stat=istat2)
 allocate(q_temp_RI_A(0:M_Bar_A-1,0:LOCAL_SIZEK-1),stat=istat2)
 allocate(q_temp_IR_A(0:M_Bar_A-1,0:LOCAL_SIZEK-1),stat=istat2)
 allocate(q_temp_II_A(0:M_Bar_A-1,0:LOCAL_SIZEK-1),stat=istat2)
!for B
 allocate(Pr_B(0:M_Bar_B-1,0:LOCAL_SIZE-1),stat=istat2)
 allocate(Pk_Real_B(0:M_Bar_B-1,0:LOCAL_SIZEK-1),stat=istat2)
 allocate(Pk_Imag_B(0:M_Bar_B-1,0:LOCAL_SIZEK-1),stat=istat2)

 allocate(q_temp_RR_B(0:M_Bar_B-1,0:LOCAL_SIZEK-1),stat=istat2)
 allocate(q_temp_RI_B(0:M_Bar_B-1,0:LOCAL_SIZEK-1),stat=istat2)
 allocate(q_temp_IR_B(0:M_Bar_B-1,0:LOCAL_SIZEK-1),stat=istat2)
 allocate(q_temp_II_B(0:M_Bar_B-1,0:LOCAL_SIZEK-1),stat=istat2)
  

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
if(myid==0) then
write(*,*) "0:WA(0,90,180) in MDE=",WA(0),WA(90),WA(180)
write(*,*) "0:WB(0,90,180)=",WB(0),WB(90),WB(180)
endif  
# endif  /* Debug */
 !!!!!
 do s=1,Nmax
 !!!s=1
   if(s==1) then
          Pr_A=0.0
     
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
        enddo !enddo k=0,LOCAL_SIZE-1  

# if defined (Debug)
    call mp_barrier()
    if(myid==0) then
     write(*,*) "Pr_A(0,30)",Pr_A(0,30)
     endif
# endif /* Debug */



!!!doing the FFT
!!Note that  in 2dcomp fft , local_in and local_out is formed as (0:x,0:y,0:z), which
!is different from fft slab decomposition version.
# if defined (R2C)
      do i=0,M_Bar_A-1
        do k=0,LOCAL_SIZE-1
          local_in_r(kD(k)%x,kD(k)%y,KD(k)%z)=Pr_A(i,k)
         enddo
       
        call decomp_2d_fft_3d(local_in_r, local_out_c)
 !      call fftw_mpi_execute_dft(fplan, local_in, local_out)
       
        do k=0,LOCAL_SIZE_K_R2C-1
            Pk_Real_A(i,k)=real(local_out_c(kDk(k)%x,kDk(k)%y,KDk(k)%z))
            Pk_Imag_A(i,k)=AIMAG(local_out_c(kDk(k)%x,kDk(k)%y,KDk(k)%z))
         enddo
       enddo  !enddo i=0,M_Bar_A-1  

# else /* R2C */
      do i=0,M_Bar_A-1
        do k=0,LOCAL_SIZE-1
          local_in(kD(k)%x,kD(k)%y,KD(k)%z)=cmplx(Pr_A(i,k),0.0)
         enddo
       
        call decomp_2d_fft_3d(local_in, local_out, DECOMP_2D_FFT_FORWARD)
 !      call fftw_mpi_execute_dft(fplan, local_in, local_out)
       
        do k=0,LOCAL_SIZE_K-1
            Pk_Real_A(i,k)=real(local_out(kDz(k)%x,kDz(k)%y,KDz(k)%z))
            Pk_Imag_A(i,k)=AIMAG(local_out(kDz(k)%x,kDz(k)%y,KDz(k)%z))
         enddo
       enddo  !enddo i=0,M_Bar_A-1  
# endif /* R2C */



# if defined (OPENMP)
        !$omp parallel default(shared) private(k)
        !$omp do
# endif  /* OPENMP */

         do k=0,LOCAL_SIZEK-1

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


# if defined (R2C)
             do i=0,M_Bar_A-1
                do k=0,LOCAL_SIZE_K_R2C-1
        local_out_c(kDk(k)%x,kDk(k)%y,KDk(k)%z)=cmplx(q_temp_RR_A(i,k)-q_temp_II_A(i,k), & 
                q_temp_RI_A(i,k)+q_temp_IR_A(i,k))
                enddo
       !call fftw_mpi_execute_dft(bplan, local_out, local_in)
        call decomp_2d_fft_3d(local_out_c, local_in_r)
            
                do k=0,LOCAL_SIZE-1
                 qA(1,i,k)=local_in_r(kD(k)%x,kD(k)%y,KD(k)%z)*M_grid_inv
                enddo
              enddo  
# else /* R2C */
             do i=0,M_Bar_A-1
                do k=0,LOCAL_SIZE_K-1
        local_out(kDz(k)%x,kDz(k)%y,KDz(k)%z)=cmplx(q_temp_RR_A(i,k)-q_temp_II_A(i,k), & 
                q_temp_RI_A(i,k)+q_temp_IR_A(i,k))
                enddo
                if(i==1) then 
                write(*,*) "i==1,sum local_out",sum(local_out),"on",myid
                endif
       !call fftw_mpi_execute_dft(bplan, local_out, local_in)
        call decomp_2d_fft_3d(local_out, local_in, DECOMP_2D_FFT_BACKWARD)
            
                do k=0,LOCAL_SIZE-1
                 qA(1,i,k)=real(local_in(kD(k)%x,kD(k)%y,KD(k)%z))*M_grid_inv
                enddo
              enddo  
# endif /* R2C */



# if defined (Debug)
    call mp_barrier()
    if(myid==0) then
     write(*,*) "qA(1,0,30)",qA(1,0,30)
     endif
# endif /* Debug */


!    call mp_barrier()
!    write(*,*) "q_temp_RR(1,1)",q_temp_RR_A(1,1)
!    write(*,*) "q_temp_RI(1,1)",q_temp_RI_A(1,1)
!    write(*,*) "q_temp_IR(1,1)",q_temp_IR_A(1,1)
!    write(*,*) "q_temp_II(1,1)",q_temp_II_A(1,1)
!    write(*,*) "sum of Pk_real",sum(Pk_real_A(:,1))  
!    write(*,*) "sum q_temp(1,:)",sum(q_temp_RR_A(1,:)),"on",myid    
!     write(*,*) "qA(1,0,1)",qA(1,0,1),"on",myid
!     write(*,*) "qA(1,1,1)",qA(1,1,1),"on",myid
!     write(*,*) "qA(1,2,1)",qA(1,2,1),"on",myid
!     write(*,*) "qA(1,3,1)",qA(1,3,1),"on",myid
!    write(*,*) "sum of qA(1,1,:)",sum(qA(1,:,1))
!    call mp_barrier()
!    stop

!!!!s=2!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    else if(s==2) then

          Pr_A=0.0
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

!!!doing the FFT
# if defined (R2C)
      do i=0,M_Bar_A-1
        do k=0,LOCAL_SIZE-1
          local_in_r(kD(k)%x,kD(k)%y,KD(k)%z)=Pr_A(i,k)
         enddo
       
        call decomp_2d_fft_3d(local_in_r, local_out_c)
       
        do k=0,LOCAL_SIZE_K_R2C-1
            Pk_Real_A(i,k)=real(local_out_c(kDk(k)%x,kDk(k)%y,KDk(k)%z))
            Pk_Imag_A(i,k)=AIMAG(local_out_c(kDk(k)%x,kDk(k)%y,KDk(k)%z))
         enddo
       enddo  !enddo i=0,M_Bar_A-1  
# else /* R2C */
      do i=0,M_Bar_A-1
        do k=0,LOCAL_SIZE-1
          local_in(kD(k)%x,kD(k)%y,KD(k)%z)=cmplx(Pr_A(i,k),0.0)
         enddo

        call decomp_2d_fft_3d(local_in, local_out, DECOMP_2D_FFT_FORWARD)

        do k=0,LOCAL_SIZE_K-1
            Pk_Real_A(i,k)=real(local_out(kDz(k)%x,kDz(k)%y,KDz(k)%z))
            Pk_Imag_A(i,k)=AIMAG(local_out(kDz(k)%x,kDz(k)%y,KDz(k)%z))
         enddo
       enddo  !enddo i=0,M_Bar_A-1  
# endif /* R2C */

         do k=0,LOCAL_SIZEK-1
           
         call mkl_dcsrsymv('l', M_Bar_A,  R_GA2_inverse_val(:,k+1), R_GA2_inverse_row(:,k+1), &
            R_GA2_inverse_col(:,k+1),Pk_Real_A(:,k), q_temp_RR_A(:,k))

         call mkl_dcsrsymv('l', M_Bar_A,  R_GA2_inverse_val(:,k+1), R_GA2_inverse_row(:,k+1), &
            R_GA2_inverse_col(:,k+1),Pk_Imag_A(:,k), q_temp_RI_A(:,k))

         call mkl_dcsrsymv('l', M_Bar_A,  I_GA2_inverse_val(:,k+1), I_GA2_inverse_row(:,k+1), &
            I_GA2_inverse_col(:,k+1),Pk_Real_A(:,k), q_temp_IR_A(:,k))

         call mkl_dcsrsymv('l', M_Bar_A,  I_GA2_inverse_val(:,k+1), I_GA2_inverse_row(:,k+1), &
            I_GA2_inverse_col(:,k+1),Pk_Imag_A(:,k), q_temp_II_A(:,k))
         enddo

# if defined (R2C)
             do i=0,M_Bar_A-1
                do k=0,LOCAL_SIZE_K_R2C-1
         local_out_c(kDk(k)%x,kDk(k)%y,KDk(k)%z)=cmplx(q_temp_RR_A(i,k)-q_temp_II_A(i,k), & 
                q_temp_RI_A(i,k)+q_temp_IR_A(i,k))
                enddo 
        call decomp_2d_fft_3d(local_out_c, local_in_r)

                do k=0,LOCAL_SIZE-1
                 qA(2,i,k)=local_in_r(kD(k)%x,kD(k)%y,KD(k)%z)*M_grid_inv
                enddo
              enddo  
# else /* R2C */
             do i=0,M_Bar_A-1
                do k=0,LOCAL_SIZE_K-1
         local_out(kDz(k)%x,kDz(k)%y,KDz(k)%z)=cmplx(q_temp_RR_A(i,k)-q_temp_II_A(i,k), & 
                q_temp_RI_A(i,k)+q_temp_IR_A(i,k))
                enddo 
        call decomp_2d_fft_3d(local_out, local_in, DECOMP_2D_FFT_BACKWARD)
       !call fftw_mpi_execute_dft(bplan, local_data, local_data)

                do k=0,LOCAL_SIZE-1
                 qA(2,i,k)=real(local_in(kD(k)%x,kD(k)%y,KD(k)%z))*M_grid_inv
                enddo
              enddo  
# endif /* R2C */


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

!!!doing the FFT
# if defined (R2C)
      do i=0,M_Bar_A-1
        do k=0,LOCAL_SIZE-1
          local_in_r(kD(k)%x,kD(k)%y,KD(k)%z)=Pr_A(i,k)
         enddo
       
        call decomp_2d_fft_3d(local_in_r, local_out_c)
 !      call fftw_mpi_execute_dft(fplan, local_in, local_out)
       
        do k=0,LOCAL_SIZE_K_R2C-1
            Pk_Real_A(i,k)=real(local_out_c(kDk(k)%x,kDk(k)%y,KDk(k)%z))
            Pk_Imag_A(i,k)=AIMAG(local_out_c(kDk(k)%x,kDk(k)%y,KDk(k)%z))
         enddo
       enddo  !enddo i=0,M_Bar_A-1  

# else /* R2C */
      do i=0,M_Bar_A-1
        do k=0,LOCAL_SIZE-1
          local_in(kD(k)%x,kD(k)%y,KD(k)%z)=cmplx(Pr_A(i,k),0.0)
         enddo
       
        call decomp_2d_fft_3d(local_in, local_out, DECOMP_2D_FFT_FORWARD)
 !      call fftw_mpi_execute_dft(fplan, local_in, local_out)
       
        do k=0,LOCAL_SIZE_K-1
            Pk_Real_A(i,k)=real(local_out(kDz(k)%x,kDz(k)%y,KDz(k)%z))
            Pk_Imag_A(i,k)=AIMAG(local_out(kDz(k)%x,kDz(k)%y,KDz(k)%z))
         enddo
       enddo  !enddo i=0,M_Bar_A-1  
# endif /* R2C */

         do k=0,LOCAL_SIZEK-1

         call mkl_dcsrsymv('l', M_Bar_A,  R_GA3_inverse_val(:,k+1), R_GA3_inverse_row(:,k+1), &
            R_GA3_inverse_col(:,k+1),Pk_Real_A(:,k), q_temp_RR_A(:,k))

         call mkl_dcsrsymv('l', M_Bar_A,  R_GA3_inverse_val(:,k+1), R_GA3_inverse_row(:,k+1), &
            R_GA3_inverse_col(:,k+1),Pk_Imag_A(:,k), q_temp_RI_A(:,k))

         call mkl_dcsrsymv('l', M_Bar_A,  I_GA3_inverse_val(:,k+1), I_GA3_inverse_row(:,k+1), &
            I_GA3_inverse_col(:,k+1),Pk_Real_A(:,k), q_temp_IR_A(:,k))

         call mkl_dcsrsymv('l', M_Bar_A,  I_GA3_inverse_val(:,k+1), I_GA3_inverse_row(:,k+1), &
            I_GA3_inverse_col(:,k+1),Pk_Imag_A(:,k), q_temp_II_A(:,k))
         enddo



# if defined (R2C)
             do i=0,M_Bar_A-1
                do k=0,LOCAL_SIZE_K_R2C-1
          local_out_c(kDk(k)%x,kDk(k)%y,KDk(k)%z)=cmplx(q_temp_RR_A(i,k)-q_temp_II_A(i,k), & 
                q_temp_RI_A(i,k)+q_temp_IR_A(i,k))
                enddo 
        call decomp_2d_fft_3d(local_out_c, local_in_r)
                do k=0,LOCAL_SIZE-1
                 qA(s,i,k)=local_in_r(kD(k)%x,kD(k)%y,KD(k)%z)*M_grid_inv
                enddo
              enddo  
# else /* R2C */
             do i=0,M_Bar_A-1
                do k=0,LOCAL_SIZE_K-1
          local_out(kDz(k)%x,kDz(k)%y,KDz(k)%z)=cmplx(q_temp_RR_A(i,k)-q_temp_II_A(i,k), & 
                q_temp_RI_A(i,k)+q_temp_IR_A(i,k))
                enddo 
        call decomp_2d_fft_3d(local_out, local_in, DECOMP_2D_FFT_BACKWARD)
       !call fftw_mpi_execute_dft(bplan, local_data, local_data)
                do k=0,LOCAL_SIZE-1
                 qA(s,i,k)=real(local_in(kD(k)%x,kD(k)%y,KD(k)%z))*M_grid_inv
                enddo
              enddo  
# endif /* R2C */

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

    else if(s==(NA+2)) then
         
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

     else if(s==(NA+3)) then
         
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

     else !else s>=3,s>NA:s>NA+3
           
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
     endif  !if(s==NA+1) 

!!!doing the FFT

# if defined (R2C)
      do i=0,M_Bar_B-1
        do k=0,LOCAL_SIZE-1
          local_in_r(kD(k)%x,kD(k)%y,KD(k)%z)=Pr_B(i,k)
         enddo
       
        call decomp_2d_fft_3d(local_in_r, local_out_c)
       
        do k=0,LOCAL_SIZE_K_R2C-1
            Pk_Real_B(i,k)=real(local_out_c(kDk(k)%x,kDk(k)%y,KDk(k)%z))
            Pk_Imag_B(i,k)=AIMAG(local_out_c(kDk(k)%x,kDk(k)%y,KDk(k)%z))
         enddo
       enddo  !enddo i=0,M_Bar_A-1  
# else /* R2C */
      do i=0,M_Bar_B-1
        do k=0,LOCAL_SIZE-1
          local_in(kD(k)%x,kD(k)%y,KD(k)%z)=cmplx(Pr_B(i,k),0.0)
         enddo
       
        call decomp_2d_fft_3d(local_in, local_out, DECOMP_2D_FFT_FORWARD)
       !call fftw_mpi_execute_dft(fplan, local_data, local_data)
       
        do k=0,LOCAL_SIZE_K-1
            Pk_Real_B(i,k)=real(local_out(kDz(k)%x,kDz(k)%y,KDz(k)%z))
            Pk_Imag_B(i,k)=AIMAG(local_out(kDz(k)%x,kDz(k)%y,KDz(k)%z))
         enddo
       enddo  !enddo i=0,M_Bar_A-1  
# endif /* R2C */
         do k=0,LOCAL_SIZEK-1

         call mkl_dcsrsymv('l', M_Bar_B,  R_GB3_inverse_val(:,k+1), R_GB3_inverse_row(:,k+1), &
            R_GB3_inverse_col(:,k+1),Pk_Real_B(:,k), q_temp_RR_B(:,k))

         call mkl_dcsrsymv('l', M_Bar_B,  R_GB3_inverse_val(:,k+1), R_GB3_inverse_row(:,k+1), &
            R_GB3_inverse_col(:,k+1),Pk_Imag_B(:,k), q_temp_RI_B(:,k))

         call mkl_dcsrsymv('l', M_Bar_B,  I_GB3_inverse_val(:,k+1), I_GB3_inverse_row(:,k+1), &
            I_GB3_inverse_col(:,k+1),Pk_Real_B(:,k), q_temp_IR_B(:,k))

         call mkl_dcsrsymv('l', M_Bar_B,  I_GB3_inverse_val(:,k+1), I_GB3_inverse_row(:,k+1), &
            I_GB3_inverse_col(:,k+1),Pk_Imag_B(:,k), q_temp_II_B(:,k))
         enddo


# if defined (R2C)
             do i=0,M_Bar_B-1
                do k=0,LOCAL_SIZE_K_R2C-1
        local_out_c(kDk(k)%x,kDk(k)%y,KDk(k)%z)=cmplx(q_temp_RR_B(i,k)-q_temp_II_B(i,k), & 
                q_temp_RI_B(i,k)+q_temp_IR_B(i,k))
                enddo 
        call decomp_2d_fft_3d(local_out_c, local_in_r)
                do k=0,LOCAL_SIZE-1
                 qB(s-NA,i,k)=local_in_r(kD(k)%x,kD(k)%y,KD(k)%z)*M_grid_inv
                enddo
              enddo 
# else /* R2C */
             do i=0,M_Bar_B-1
                do k=0,LOCAL_SIZE_K-1
        local_out(kDz(k)%x,kDz(k)%y,KDz(k)%z)=cmplx(q_temp_RR_B(i,k)-q_temp_II_B(i,k), & 
                q_temp_RI_B(i,k)+q_temp_IR_B(i,k))
                enddo 
        call decomp_2d_fft_3d(local_out, local_in, DECOMP_2D_FFT_BACKWARD)
       !call fftw_mpi_execute_dft(bplan, local_data, local_data)
                do k=0,LOCAL_SIZE-1
                 qB(s-NA,i,k)=real(local_in(kD(k)%x,kD(k)%y,KD(k)%z))*M_grid_inv
                enddo
              enddo 
# endif /* R2C */
 
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

    !call mp_barrier()
    ! write(*,*) "qA(5,0,1)",qA(5,0,1),"on",myid
    ! write(*,*) "qA(5,1,1)",qA(5,1,1),"on",myid
    ! write(*,*) "qA(5,2,1)",qA(5,2,1),"on",myid
    ! write(*,*) "qA(NA,0,1)",qA(NA,0,1),"on",myid
    ! write(*,*) "qA(NA,1,1)",qA(NA,1,1),"on",myid
    ! write(*,*) "qA(NA,2,1)",qA(NA,2,1),"on",myid
    !write(*,*) "sum of qA(NA,1,:)",sum(qA(NA,:,1))
    !call mp_barrier()
    !stop

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

      


subroutine MDE_qstar(R2C_FFT)
 USE nrtype,only :DP,PI
 USE global_para
 USE mpi
 USE control
 USE constants
 USE utility
 USE mmpi
 use decomp_fft_mod
 USE SPH
 USE matrix_inverse
 implicit none
 logical,optional :: R2C_FFT
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
 Integer :: LOCAL_SIZEK
 !!!! just for reference


# if defined (R2C)
LOCAL_SIZEK=LOCAL_SIZE_K_R2C
# else /* R2C */
LOCAL_SIZEK=LOCAL_SIZE_K
# endif /* R2C */
!for A
 allocate(Pr_A(0:M_Bar_A-1,0:LOCAL_SIZE-1),stat=istat2)
 allocate(Pk_Real_A(0:M_Bar_A-1,0:LOCAL_SIZEK-1),stat=istat2)
 allocate(Pk_Imag_A(0:M_Bar_A-1,0:LOCAL_SIZEK-1),stat=istat2)


 allocate(qstar_temp_RR_A(0:M_Bar_A-1,0:LOCAL_SIZEK-1),stat=istat2)
 allocate(qstar_temp_RI_A(0:M_Bar_A-1,0:LOCAL_SIZEK-1),stat=istat2)
 allocate(qstar_temp_IR_A(0:M_Bar_A-1,0:LOCAL_SIZEK-1),stat=istat2)
 allocate(qstar_temp_II_A(0:M_Bar_A-1,0:LOCAL_SIZEK-1),stat=istat2)
!for B
 allocate(Pr_B(0:M_Bar_B-1,0:LOCAL_SIZE-1),stat=istat2)
 allocate(Pk_Real_B(0:M_Bar_B-1,0:LOCAL_SIZEK-1),stat=istat2)
 allocate(Pk_Imag_B(0:M_Bar_B-1,0:LOCAL_SIZEK-1),stat=istat2)

 allocate(qstar_temp_RR_B(0:M_Bar_B-1,0:LOCAL_SIZEK-1),stat=istat2)
 allocate(qstar_temp_RI_B(0:M_Bar_B-1,0:LOCAL_SIZEK-1),stat=istat2)
 allocate(qstar_temp_IR_B(0:M_Bar_B-1,0:LOCAL_SIZEK-1),stat=istat2)
 allocate(qstar_temp_II_B(0:M_Bar_B-1,0:LOCAL_SIZEK-1),stat=istat2)

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

!!!doing the FFT
# if defined (R2C)
      do i=0,M_Bar_B-1
        do k=0,LOCAL_SIZE-1
          local_in_r(kD(k)%x,kD(k)%y,KD(k)%z)=Pr_B(i,k)
         enddo
        call decomp_2d_fft_3d(local_in_r, local_out_c)
        do k=0,LOCAL_SIZE_K_R2C-1
            Pk_Real_B(i,k)=real(local_out_c(kDk(k)%x,kDk(k)%y,KDk(k)%z))
            Pk_Imag_B(i,k)=AIMAG(local_out_c(kDk(k)%x,kDk(k)%y,KDk(k)%z))
         enddo
       enddo  !enddo i=0,M_Bar_B-1  
# else /* R2C */
      do i=0,M_Bar_B-1
        do k=0,LOCAL_SIZE-1
          local_in(kD(k)%x,kD(k)%y,KD(k)%z)=cmplx(Pr_B(i,k),0.0)
         enddo
        call decomp_2d_fft_3d(local_in, local_out, DECOMP_2D_FFT_FORWARD)
       !call fftw_mpi_execute_dft(fplan, local_data, local_data)
        do k=0,LOCAL_SIZE_K-1
            Pk_Real_B(i,k)=real(local_out(kDz(k)%x,kDz(k)%y,KDz(k)%z))
            Pk_Imag_B(i,k)=AIMAG(local_out(kDz(k)%x,kDz(k)%y,KDz(k)%z))
         enddo
       enddo  !enddo i=0,M_Bar_B-1  
# endif /* R2C */


         do k=0,LOCAL_SIZEK-1

         call mkl_dcsrsymv('l', M_Bar_B,  R_GB1_star_inverse_val(:,k+1), R_GB1_star_inverse_row(:,k+1), &
            R_GB1_star_inverse_col(:,k+1),Pk_Real_B(:,k), qstar_temp_RR_B(:,k))

         call mkl_dcsrsymv('l', M_Bar_B,  R_GB1_star_inverse_val(:,k+1), R_GB1_star_inverse_row(:,k+1), &
            R_GB1_star_inverse_col(:,k+1),Pk_Imag_B(:,k), qstar_temp_RI_B(:,k))

         call mkl_dcsrsymv('l', M_Bar_B,  I_GB1_star_inverse_val(:,k+1), I_GB1_star_inverse_row(:,k+1), &
            I_GB1_star_inverse_col(:,k+1),Pk_Real_B(:,k), qstar_temp_IR_B(:,k))

         call mkl_dcsrsymv('l', M_Bar_B,  I_GB1_star_inverse_val(:,k+1), I_GB1_star_inverse_row(:,k+1), &
            I_GB1_star_inverse_col(:,k+1),Pk_Imag_B(:,k), qstar_temp_II_B(:,k))


         enddo

# if defined (R2C)
             do i=0,M_Bar_B-1
                do k=0,LOCAL_SIZE_K_R2C-1
    local_out_c(kDk(k)%x,kDk(k)%y,KDk(k)%z)=cmplx(qstar_temp_RR_B(i,k)-qstar_temp_II_B(i,k), & 
                qstar_temp_RI_B(i,k)+qstar_temp_IR_B(i,k))
                enddo 
        call decomp_2d_fft_3d(local_out_c, local_in_r)
                do k=0,LOCAL_SIZE-1
                 qBstar(NB-1,i,k)=local_in_r(kD(k)%x,kD(k)%y,KD(k)%z)*M_grid_inv
                enddo
              enddo  
# else /* R2C */
             do i=0,M_Bar_B-1
                do k=0,LOCAL_SIZE_K-1
    local_out(kDz(k)%x,kDz(k)%y,KDz(k)%z)=cmplx(qstar_temp_RR_B(i,k)-qstar_temp_II_B(i,k), & 
                qstar_temp_RI_B(i,k)+qstar_temp_IR_B(i,k))
                enddo 
        call decomp_2d_fft_3d(local_out, local_in, DECOMP_2D_FFT_BACKWARD)
       !call fftw_mpi_execute_dft(bplan, local_data, local_data)
                do k=0,LOCAL_SIZE-1
                 qBstar(NB-1,i,k)=real(local_in(kD(k)%x,kD(k)%y,KD(k)%z))*M_grid_inv
                enddo
              enddo  
# endif /* R2C */


!!!!s=NMAX-2!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    else if(s==NMAX-2) then

          Pr_B=0.0
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

!!!doing the FFT

# if defined (R2C)
      do i=0,M_Bar_B-1
        do k=0,LOCAL_SIZE-1
          local_in_r(kD(k)%x,kD(k)%y,KD(k)%z)=Pr_B(i,k)
         enddo
        call decomp_2d_fft_3d(local_in_r, local_out_c)
        do k=0,LOCAL_SIZE_K_R2C-1
            Pk_Real_B(i,k)=real(local_out_c(kDk(k)%x,kDk(k)%y,KDk(k)%z))
            Pk_Imag_B(i,k)=AIMAG(local_out_c(kDk(k)%x,kDk(k)%y,KDk(k)%z))
         enddo
       enddo  !enddo i=0,M_Bar_B-1  
# else /* R2C */
      do i=0,M_Bar_B-1
        do k=0,LOCAL_SIZE-1
          local_in(kD(k)%x,kD(k)%y,KD(k)%z)=cmplx(Pr_B(i,k),0.0)
         enddo
        call decomp_2d_fft_3d(local_in, local_out, DECOMP_2D_FFT_FORWARD)
       !call fftw_mpi_execute_dft(fplan, local_data, local_data)
        do k=0,LOCAL_SIZE_K-1
            Pk_Real_B(i,k)=real(local_out(kDz(k)%x,kDz(k)%y,KDz(k)%z))
            Pk_Imag_B(i,k)=AIMAG(local_out(kDz(k)%x,kDz(k)%y,KDz(k)%z))
         enddo
       enddo  !enddo i=0,M_Bar_B-1  
# endif /* R2C */

         do k=0,LOCAL_SIZEK-1

         call mkl_dcsrsymv('l', M_Bar_B,  R_GB2_star_inverse_val(:,k+1), R_GB2_star_inverse_row(:,k+1), &
            R_GB2_star_inverse_col(:,k+1),Pk_Real_B(:,k), qstar_temp_RR_B(:,k))

         call mkl_dcsrsymv('l', M_Bar_B,  R_GB2_star_inverse_val(:,k+1), R_GB2_star_inverse_row(:,k+1), &
            R_GB2_star_inverse_col(:,k+1),Pk_Imag_B(:,k), qstar_temp_RI_B(:,k))

         call mkl_dcsrsymv('l', M_Bar_B,  I_GB2_star_inverse_val(:,k+1), I_GB2_star_inverse_row(:,k+1), &
            I_GB2_star_inverse_col(:,k+1),Pk_Real_B(:,k), qstar_temp_IR_B(:,k))

         call mkl_dcsrsymv('l', M_Bar_B,  I_GB2_star_inverse_val(:,k+1), I_GB2_star_inverse_row(:,k+1), &
            I_GB2_star_inverse_col(:,k+1),Pk_Imag_B(:,k), qstar_temp_II_B(:,k))
         enddo

# if defined (R2C)
             do i=0,M_Bar_B-1
                do k=0,LOCAL_SIZE_K_R2C-1
    local_out_c(kDk(k)%x,kDk(k)%y,KDk(k)%z)=cmplx(qstar_temp_RR_B(i,k)-qstar_temp_II_B(i,k), & 
                qstar_temp_RI_B(i,k)+qstar_temp_IR_B(i,k))
                enddo 
        call decomp_2d_fft_3d(local_out_c, local_in_r)
                do k=0,LOCAL_SIZE-1
                 qBstar(NB-2,i,k)=local_in_r(kD(k)%x,kD(k)%y,KD(k)%z)*M_grid_inv
                enddo
              enddo  
# else /* R2C */
             do i=0,M_Bar_B-1
                do k=0,LOCAL_SIZE_K-1
    local_out(kDz(k)%x,kDz(k)%y,KDz(k)%z)=cmplx(qstar_temp_RR_B(i,k)-qstar_temp_II_B(i,k), & 
                qstar_temp_RI_B(i,k)+qstar_temp_IR_B(i,k))
                enddo 
        call decomp_2d_fft_3d(local_out, local_in, DECOMP_2D_FFT_BACKWARD)
       !call fftw_mpi_execute_dft(bplan, local_data, local_data)
                do k=0,LOCAL_SIZE-1
                 qBstar(NB-2,i,k)=real(local_in(kD(k)%x,kD(k)%y,KD(k)%z))*M_grid_inv
                enddo
              enddo  
# endif /* R2C */

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

!!!doing the FFT
# if defined (R2C)
      do i=0,M_Bar_B-1
        do k=0,LOCAL_SIZE-1
          local_in_r(kD(k)%x,kD(k)%y,KD(k)%z)=Pr_B(i,k)
         enddo
        call decomp_2d_fft_3d(local_in_r, local_out_c)
        do k=0,LOCAL_SIZE_K_R2C-1
            Pk_Real_B(i,k)=real(local_out_c(kDk(k)%x,kDk(k)%y,KDk(k)%z))
            Pk_Imag_B(i,k)=AIMAG(local_out_c(kDk(k)%x,kDk(k)%y,KDk(k)%z))
         enddo
       enddo  !enddo i=0,M_Bar_B-1  
# else /* R2C */
      do i=0,M_Bar_B-1
        do k=0,LOCAL_SIZE-1
          local_in(kD(k)%x,kD(k)%y,KD(k)%z)=cmplx(Pr_B(i,k),0.0)
         enddo
        call decomp_2d_fft_3d(local_in, local_out, DECOMP_2D_FFT_FORWARD)
       !call fftw_mpi_execute_dft(fplan, local_data, local_data)
        do k=0,LOCAL_SIZE_K-1
            Pk_Real_B(i,k)=real(local_out(kDz(k)%x,kDz(k)%y,KDz(k)%z))
            Pk_Imag_B(i,k)=AIMAG(local_out(kDz(k)%x,kDz(k)%y,KDz(k)%z))
         enddo
       enddo  !enddo i=0,M_Bar_B-1  
# endif /* R2C */

         do k=0,LOCAL_SIZEK-1

         call mkl_dcsrsymv('l', M_Bar_B,  R_GB3_star_inverse_val(:,k+1), R_GB3_star_inverse_row(:,k+1), &
            R_GB3_star_inverse_col(:,k+1),Pk_Real_B(:,k), qstar_temp_RR_B(:,k))

         call mkl_dcsrsymv('l', M_Bar_B,  R_GB3_star_inverse_val(:,k+1), R_GB3_star_inverse_row(:,k+1), &
            R_GB3_star_inverse_col(:,k+1),Pk_Imag_B(:,k), qstar_temp_RI_B(:,k))

         call mkl_dcsrsymv('l', M_Bar_B,  I_GB3_star_inverse_val(:,k+1), I_GB3_star_inverse_row(:,k+1), &
            I_GB3_star_inverse_col(:,k+1),Pk_Real_B(:,k), qstar_temp_IR_B(:,k))

         call mkl_dcsrsymv('l', M_Bar_B,  I_GB3_star_inverse_val(:,k+1), I_GB3_star_inverse_row(:,k+1), &
            I_GB3_star_inverse_col(:,k+1),Pk_Imag_B(:,k), qstar_temp_II_B(:,k))
         enddo


# if defined (R2C)
             do i=0,M_Bar_B-1
                do k=0,LOCAL_SIZE_K_R2C-1
   local_out_c(kDk(k)%x,kDk(k)%y,KDk(k)%z)=cmplx(qstar_temp_RR_B(i,k)-qstar_temp_II_B(i,k), & 
                qstar_temp_RI_B(i,k)+qstar_temp_IR_B(i,k))
                enddo 
        call decomp_2d_fft_3d(local_out_c, local_in_r)
                do k=0,LOCAL_SIZE-1
                 qBstar(s-NA,i,k)=local_in_r(kD(k)%x,kD(k)%y,KD(k)%z)*M_grid_inv
                enddo
              enddo  
# else /* R2C */
             do i=0,M_Bar_B-1
                do k=0,LOCAL_SIZE_K-1
   local_out(kDz(k)%x,kDz(k)%y,KDz(k)%z)=cmplx(qstar_temp_RR_B(i,k)-qstar_temp_II_B(i,k), & 
                qstar_temp_RI_B(i,k)+qstar_temp_IR_B(i,k))
                enddo 
        call decomp_2d_fft_3d(local_out, local_in, DECOMP_2D_FFT_BACKWARD)
       !call fftw_mpi_execute_dft(bplan, local_data, local_data)
                do k=0,LOCAL_SIZE-1
                 qBstar(s-NA,i,k)=real(local_in(kD(k)%x,kD(k)%y,KD(k)%z))*M_grid_inv
                enddo
              enddo  
# endif /* R2C */

   else  ! else s<=NA :entering A block regime
            !!!calculate polymer block A
             Pr_A=0.0
            
   if(s==NA) then
           Pr_B=0.0  
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

    else if(s==(NA-1)) then
         
            Pr_B=0.0

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

     else if(s==(NA-2)) then
         Pr_B=0.0
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

     else 
           
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
     endif   

!!!doing the FFT
# if defined (R2C)
      do i=0,M_Bar_A-1
        do k=0,LOCAL_SIZE-1
          local_in_r(kD(k)%x,kD(k)%y,KD(k)%z)=Pr_A(i,k)
         enddo
        call decomp_2d_fft_3d(local_in_r, local_out_c)
        do k=0,LOCAL_SIZE_K_R2C-1
            Pk_Real_A(i,k)=real(local_out_c(kDk(k)%x,kDk(k)%y,KDk(k)%z))
            Pk_Imag_A(i,k)=AIMAG(local_out_c(kDk(k)%x,kDk(k)%y,KDk(k)%z))
         enddo
       enddo  !enddo i=0,M_Bar_A-1  
# else /* R2C */
      do i=0,M_Bar_A-1
        do k=0,LOCAL_SIZE-1
          local_in(kD(k)%x,kD(k)%y,KD(k)%z)=cmplx(Pr_A(i,k),0.0)
         enddo
        call decomp_2d_fft_3d(local_in, local_out, DECOMP_2D_FFT_FORWARD)
        do k=0,LOCAL_SIZE_K-1
            Pk_Real_A(i,k)=real(local_out(kDz(k)%x,kDz(k)%y,KDz(k)%z))
            Pk_Imag_A(i,k)=AIMAG(local_out(kDz(k)%x,kDz(k)%y,KDz(k)%z))
         enddo
       enddo  !enddo i=0,M_Bar_A-1  
# endif /* R2C */

         do k=0,LOCAL_SIZEK-1

         call mkl_dcsrsymv('l', M_Bar_A,  R_GA3_star_inverse_val(:,k+1), R_GA3_star_inverse_row(:,k+1), &
            R_GA3_star_inverse_col(:,k+1),Pk_Real_A(:,k), qstar_temp_RR_A(:,k))

         call mkl_dcsrsymv('l', M_Bar_A,  R_GA3_star_inverse_val(:,k+1), R_GA3_star_inverse_row(:,k+1), &
            R_GA3_star_inverse_col(:,k+1),Pk_Imag_A(:,k), qstar_temp_RI_A(:,k))

         call mkl_dcsrsymv('l', M_Bar_A,  I_GA3_star_inverse_val(:,k+1), I_GA3_star_inverse_row(:,k+1), &
            I_GA3_star_inverse_col(:,k+1),Pk_Real_A(:,k), qstar_temp_IR_A(:,k))

         call mkl_dcsrsymv('l', M_Bar_A,  I_GA3_star_inverse_val(:,k+1), I_GA3_star_inverse_row(:,k+1), &
            I_GA3_star_inverse_col(:,k+1),Pk_Imag_A(:,k), qstar_temp_II_A(:,k))
         enddo


# if defined (R2C)
             do i=0,M_Bar_A-1
                do k=0,LOCAL_SIZE_K_R2C-1
  local_out_c(kDk(k)%x,kDk(k)%y,KDk(k)%z)=cmplx(qstar_temp_RR_A(i,k)-qstar_temp_II_A(i,k), & 
                qstar_temp_RI_A(i,k)+qstar_temp_IR_A(i,k))
                enddo 
        call decomp_2d_fft_3d(local_out_c, local_in_r)
       !call fftw_mpi_execute_dft(bplan, local_data, local_data)
                do k=0,LOCAL_SIZE-1
                 qAstar(s,i,k)=local_in_r(kD(k)%x,kD(k)%y,KD(k)%z)*M_grid_inv
                enddo
              enddo
# else /* R2C */
             do i=0,M_Bar_A-1
                do k=0,LOCAL_SIZE_K-1
  local_out(kDz(k)%x,kDz(k)%y,KDz(k)%z)=cmplx(qstar_temp_RR_A(i,k)-qstar_temp_II_A(i,k), & 
                qstar_temp_RI_A(i,k)+qstar_temp_IR_A(i,k))
                enddo 
        call decomp_2d_fft_3d(local_out, local_in, DECOMP_2D_FFT_BACKWARD)
       !call fftw_mpi_execute_dft(bplan, local_data, local_data)
                do k=0,LOCAL_SIZE-1
                 qAstar(s,i,k)=real(local_in(kD(k)%x,kD(k)%y,KD(k)%z))*M_grid_inv
                enddo
              enddo
# endif /* R2C */ 
 
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
end module  wormlike_MDE 
# endif /* SLAB */


# endif /* SPH */


