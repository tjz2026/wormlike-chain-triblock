# if defined (SLAB)
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
 !for A
 real(DP) :: tempA1,tempA2,tempA3,temp,tempx,Wtemp,sumaq,sumaq_global,sumaq1
 real(DP) :: sumPA,sumPA_global
 Complex(DP),allocatable :: G_matA(:,:)
 Complex(DP),allocatable :: G_matB(:,:)
 real(DP),allocatable :: Pr_A(:,:)
 Complex(DP),allocatable :: Pk_A(:,:)

 Real(DP),allocatable :: Pk_Ar(:)
 Real(DP),allocatable :: Pk_Ai(:)
 Real(DP),allocatable :: qtmp_Arr(:)
 Real(DP),allocatable :: qtmp_Ari(:)
 Real(DP),allocatable :: qtmp_Air(:)
 Real(DP),allocatable :: qtmp_Aii(:)
 Real(DP),allocatable :: Pk_Br(:)
 Real(DP),allocatable :: Pk_Bi(:)
 Real(DP),allocatable :: qtmp_Brr(:)
 Real(DP),allocatable :: qtmp_Bri(:)
 Real(DP),allocatable :: qtmp_Bir(:)
 Real(DP),allocatable :: qtmp_Bii(:)

 Complex(DP),allocatable :: q_temp_A(:,:)
 !for B
 real(DP),allocatable :: Pr_B(:,:)
 Complex(DP),allocatable :: Pk_B(:,:)
 Complex(DP),allocatable :: q_temp_B(:,:)
 REAL(DP) :: M_grid_inv

 !!!! just for reference
 !first ,allocate the matrixs
 allocate(G_matA(0:M_Bar_A-1,0:M_Bar_A-1),stat=istat2)
 allocate(G_matB(0:M_Bar_B-1,0:M_Bar_B-1),stat=istat2)
!for A
 allocate(Pr_A(0:M_Bar_A-1,0:LOCAL_SIZE-1),stat=istat2)
 allocate(Pk_A(0:M_Bar_A-1,0:LOCAL_SIZE-1),stat=istat2)
 allocate(Pk_Ar(0:M_Bar_A-1),stat=istat2)
 allocate(qtmp_Arr(0:M_Bar_A-1),stat=istat2)
 allocate(qtmp_Ari(0:M_Bar_A-1),stat=istat2)
 allocate(Pk_Ai(0:M_Bar_A-1),stat=istat2)
 allocate(qtmp_Air(0:M_Bar_A-1),stat=istat2)
 allocate(qtmp_Aii(0:M_Bar_A-1),stat=istat2)
 !allocate(q_temp_A(0:M_Bar_A-1,0:LOCAL_SIZE-1),stat=istat2)
!for B
 allocate(Pr_B(0:M_Bar_B-1,0:LOCAL_SIZE-1),stat=istat2)
 allocate(Pk_B(0:M_Bar_B-1,0:LOCAL_SIZE-1),stat=istat2)
 allocate(Pk_Br(0:M_Bar_B-1),stat=istat2)
 allocate(qtmp_Brr(0:M_Bar_B-1),stat=istat2)
 allocate(qtmp_Bri(0:M_Bar_B-1),stat=istat2)
 allocate(Pk_Bi(0:M_Bar_B-1),stat=istat2)
 allocate(qtmp_Bir(0:M_Bar_B-1),stat=istat2)
 allocate(qtmp_Bii(0:M_Bar_B-1),stat=istat2)
 !allocate(q_temp_B(0:M_Bar_B-1,0:LOCAL_SIZE-1),stat=istat2)


 if(istat2/=0) then
 write(*,*) "allocate matrix in mde_q failed"
 stop
 endif


 M_grid_inv=1.0d0/M_grid

do k=0,LOCAL_SIZE-1
   do i=0,M_Bar_A-1
    qA(i,K,0)=0.0
    if(i==0) qA(i,K,0)=1.0
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
             Pr_A(i,k)=Pr_A(i,k)+ds*tempx*qA(j,k,0)
           enddo
# endif /* MaierSaupe */

            do i=0,M_Bar_A-1
            Pr_A(i,k)=Pr_A(i,k)+(1.0-ds*Wtemp)*qA(i,k,0)
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
            Pk_A(i,k)=local_data(kD(k)%x,kD(k)%y)
# else  /* Dim2 */
            Pk_A(i,k)=local_data(kD(k)%x,kD(k)%y,KD(k)%z)
# endif /* Dim2 */
         enddo
       enddo  !enddo i=0,M_Bar_A-1  

# if defined (OPENMP)
        !$omp parallel default(shared) private(k)
        !$omp do
# endif  /* OPENMP */
         do k=0,LOCAL_SIZE-1
         !call init_G_mat(G_matA,M_Bar_A,1,kD(k)%z,kD(k)%y,kD(k)%x)
         !call solver_z(G_matA,Pk_A(:,k),M_Bar_A)
           Pk_Ar(:)=real(Pk_A(:,k))
           Pk_Ai(:)=aimag(Pk_A(:,k))

        call mkl_dcsrsymv('l', M_Bar_A,  R_GA1_inverse_val(:,k+1), R_GA1_inverse_row(:,k+1), &
           R_GA1_inverse_col(:,k+1),Pk_Ar, qtmp_Arr)

         call mkl_dcsrsymv('l', M_Bar_A,  R_GA1_inverse_val(:,k+1), R_GA1_inverse_row(:,k+1), &
            R_GA1_inverse_col(:,k+1),Pk_Ai, qtmp_Ari)

         call mkl_dcsrsymv('l', M_Bar_A,  I_GA1_inverse_val(:,k+1), I_GA1_inverse_row(:,k+1), &
            I_GA1_inverse_col(:,k+1),Pk_Ar, qtmp_Air)

         call mkl_dcsrsymv('l', M_Bar_A,  I_GA1_inverse_val(:,k+1), I_GA1_inverse_row(:,k+1), &
            I_GA1_inverse_col(:,k+1),Pk_Ai, qtmp_Aii)

            pk_A(:,k)=dcmplx(qtmp_Arr(:)-qtmp_Aii(:),qtmp_Ari(:)+qtmp_Air(:)) 

        enddo

# if defined (Debug)
    call mp_barrier()
    if(myid==0) then
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
        local_data(kD(k)%x,kD(k)%y)=Pk_A(i,k)
# else  /* Dim2 */
        local_data(kD(k)%x,kD(k)%y,KD(k)%z)=pk_A(i,k)
# endif  /* Dim2 */
                enddo
       call fftw_mpi_execute_dft(bplan, local_data, local_data)
            
                do k=0,LOCAL_SIZE-1
# if defined (Dim2)
                 qA(i,k,1)=real(local_data(kD(k)%x,kD(k)%y))*M_grid_inv
# else  /* Dim2 */
                 qA(i,k,1)=real(local_data(kD(k)%x,kD(k)%y,KD(k)%z))*M_grid_inv
# endif  /* Dim2 */
                enddo
              enddo  

!# if defined (Debug)
    call mp_barrier()
    if(myid==0) then
     write(*,*) "qA(0,0,1)",qA(0,0,1)
     write(*,*) "qA(0,1,1)",qA(0,1,1)
     endif
!# endif /* Debug */
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
             Pr_A(i,k)=Pr_A(i,k)+ds*tempx*(2.0*qA(j,k,1)-qA(j,k,0))
           enddo
# endif /* MaierSaupe */
            do i=0,M_Bar_A-1
            Pr_A(i,k)=Pr_A(i,k)+2.0*qA(i,k,1)-0.5*qA(i,k,0) - &
                      ds*Wtemp*(2.0*qA(i,k,1)-qA(i,k,0))
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
            Pk_A(i,k)=local_data(kD(k)%x,kD(k)%y)
# else  /* Dim2 */
            Pk_A(i,k)=local_data(kD(k)%x,kD(k)%y,KD(k)%z)
# endif /* Dim2 */
         enddo
       enddo  !enddo i=0,M_Bar_A-1  


# if defined (OPENMP)
        !$omp parallel default(shared) private(k)
        !$omp do
# endif  /* OPENMP */
         do k=0,LOCAL_SIZE-1
           
       !  call init_G_mat(G_matA,M_Bar_A,2,kD(k)%z,kD(k)%y,kD(k)%x)
       !  call solver_z(G_matA,Pk_A(:,k),M_Bar_A)
       !   call mkl_dcsrsymv('l', M_Bar_A,  R_GA2_inverse_val(:,k+1), R_GA2_inverse_row(:,k+1), &
       !      R_GA2_inverse_col(:,k+1),Pk_Real_A(:,k), q_temp_RR_A(:,k))

       !  call mkl_dcsrsymv('l', M_Bar_A,  R_GA2_inverse_val(:,k+1), R_GA2_inverse_row(:,k+1), &
       !     R_GA2_inverse_col(:,k+1),Pk_Imag_A(:,k), q_temp_RI_A(:,k))

       !  call mkl_dcsrsymv('l', M_Bar_A,  I_GA2_inverse_val(:,k+1), I_GA2_inverse_row(:,k+1), &
       !     I_GA2_inverse_col(:,k+1),Pk_Real_A(:,k), q_temp_IR_A(:,k))

       !  call mkl_dcsrsymv('l', M_Bar_A,  I_GA2_inverse_val(:,k+1), I_GA2_inverse_row(:,k+1), &
       !     I_GA2_inverse_col(:,k+1),Pk_Imag_A(:,k), q_temp_II_A(:,k))
           Pk_Ar(:)=real(Pk_A(:,k))
           Pk_Ai(:)=aimag(Pk_A(:,k))

        call mkl_dcsrsymv('l', M_Bar_A,  R_GA2_inverse_val(:,k+1), R_GA2_inverse_row(:,k+1), &
           R_GA2_inverse_col(:,k+1),Pk_Ar, qtmp_Arr)

         call mkl_dcsrsymv('l', M_Bar_A,  R_GA2_inverse_val(:,k+1), R_GA2_inverse_row(:,k+1), &
            R_GA2_inverse_col(:,k+1),Pk_Ai, qtmp_Ari)

         call mkl_dcsrsymv('l', M_Bar_A,  I_GA2_inverse_val(:,k+1), I_GA2_inverse_row(:,k+1), &
            I_GA2_inverse_col(:,k+1),Pk_Ar, qtmp_Air)

         call mkl_dcsrsymv('l', M_Bar_A,  I_GA2_inverse_val(:,k+1), I_GA2_inverse_row(:,k+1), &
            I_GA2_inverse_col(:,k+1),Pk_Ai, qtmp_Aii)

            pk_A(:,k)=dcmplx(qtmp_Arr(:)-qtmp_Aii(:),qtmp_Ari(:)+qtmp_Air(:)) 
         enddo
# if defined (OPENMP)
        !$omp end  do
        !$omp end  parallel
# endif  /* OPENMP */

             do i=0,M_Bar_A-1
                do k=0,LOCAL_SIZE-1
# if defined (Dim2)
            local_data(kD(k)%x,kD(k)%y)=Pk_A(i,k)
# else  /* Dim2 */
            local_data(kD(k)%x,kD(k)%y,KD(k)%z)=Pk_A(i,k)
# endif /* Dim2 */
                enddo 
       call fftw_mpi_execute_dft(bplan, local_data, local_data)
                do k=0,LOCAL_SIZE-1
# if defined (Dim2)
                 qA(i,k,2)=real(local_data(kD(k)%x,kD(k)%y))*M_grid_inv
# else  /* Dim2 */
                 qA(i,k,2)=real(local_data(kD(k)%x,kD(k)%y,KD(k)%z))*M_grid_inv
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
             Pr_A(i,k)=Pr_A(i,k)+ds*tempx*(3.0*qA(j,k,s-1)- &
                       3.0*qA(j,k,s-2)+qA(j,k,s-3))
            
           enddo
# endif /* MaierSaupe */
            do i=0,M_Bar_A-1
            Pr_A(i,k)=Pr_A(i,k)+3.0*qA(i,k,s-1)-1.5*qA(i,k,s-2)+ &
                      temp*qA(i,k,s-3) - &
                      ds*Wtemp*(3.0*qA(i,k,s-1)-3.0*qA(i,k,s-2)+qA(i,k,s-3))
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
            Pk_A(i,k)=local_data(kD(k)%x,kD(k)%y)
# else  /* Dim2 */
            Pk_A(i,k)=local_data(kD(k)%x,kD(k)%y,KD(k)%z)
# endif /* Dim2 */
         enddo
       enddo  !enddo i=0,M_Bar_A-1  


# if defined (OPENMP)
        !$omp parallel default(shared) private(k)
        !$omp do
# endif  /* OPENMP */
         do k=0,LOCAL_SIZE-1

      !   call init_G_mat(G_matA,M_Bar_A,3,kD(k)%z,kD(k)%y,kD(k)%x)
      !   call solver_z(G_matA,Pk_A(:,k),M_Bar_A)
      !   call mkl_dcsrsymv('l', M_Bar_A,  R_GA3_inverse_val(:,k+1), R_GA3_inverse_row(:,k+1), &
      !      R_GA3_inverse_col(:,k+1),Pk_Real_A(:,k), q_temp_RR_A(:,k))

      !   call mkl_dcsrsymv('l', M_Bar_A,  R_GA3_inverse_val(:,k+1), R_GA3_inverse_row(:,k+1), &
      !      R_GA3_inverse_col(:,k+1),Pk_Imag_A(:,k), q_temp_RI_A(:,k))

      !   call mkl_dcsrsymv('l', M_Bar_A,  I_GA3_inverse_val(:,k+1), I_GA3_inverse_row(:,k+1), &
      !      I_GA3_inverse_col(:,k+1),Pk_Real_A(:,k), q_temp_IR_A(:,k))

      !   call mkl_dcsrsymv('l', M_Bar_A,  I_GA3_inverse_val(:,k+1), I_GA3_inverse_row(:,k+1), &
      !      I_GA3_inverse_col(:,k+1),Pk_Imag_A(:,k), q_temp_II_A(:,k))
           Pk_Ar(:)=real(Pk_A(:,k))
           Pk_Ai(:)=aimag(Pk_A(:,k))

        call mkl_dcsrsymv('l', M_Bar_A,  R_GA3_inverse_val(:,k+1), R_GA3_inverse_row(:,k+1), &
           R_GA3_inverse_col(:,k+1),Pk_Ar, qtmp_Arr)

         call mkl_dcsrsymv('l', M_Bar_A,  R_GA3_inverse_val(:,k+1), R_GA3_inverse_row(:,k+1), &
            R_GA3_inverse_col(:,k+1),Pk_Ai, qtmp_Ari)

         call mkl_dcsrsymv('l', M_Bar_A,  I_GA3_inverse_val(:,k+1), I_GA3_inverse_row(:,k+1), &
            I_GA3_inverse_col(:,k+1),Pk_Ar, qtmp_Air)

         call mkl_dcsrsymv('l', M_Bar_A,  I_GA3_inverse_val(:,k+1), I_GA3_inverse_row(:,k+1), &
            I_GA3_inverse_col(:,k+1),Pk_Ai, qtmp_Aii)

            pk_A(:,k)=dcmplx(qtmp_Arr(:)-qtmp_Aii(:),qtmp_Ari(:)+qtmp_Air(:)) 
         enddo
# if defined (OPENMP)
        !$omp end  do
        !$omp end  parallel
# endif  /* OPENMP */

             do i=0,M_Bar_A-1
                do k=0,LOCAL_SIZE-1
# if defined (Dim2)
          local_data(kD(k)%x,kD(k)%y)=Pk_A(i,k)
# else  /* Dim2 */
          local_data(kD(k)%x,kD(k)%y,KD(k)%z)=Pk_A(i,k)
# endif /* Dim2 */
                enddo 
       call fftw_mpi_execute_dft(bplan, local_data, local_data)
                do k=0,LOCAL_SIZE-1
# if defined (Dim2)
                 qA(i,k,s)=real(local_data(kD(k)%x,kD(k)%y))*M_grid_inv
# else  /* Dim2 */
                 qA(i,k,s)=real(local_data(kD(k)%x,kD(k)%y,KD(k)%z))*M_grid_inv
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
              tempA1=qA(j,k,NA)
              tempA2=qA(j,k,NA-1)
              tempA3=qA(j,k,NA-2)
              endif
              Pr_B(i,k)=Pr_B(i,k) + ds*tempx*(3.0*tempA1-3.0*tempA2+tempA3)
           enddo  !enddo index_nonzero
# endif /* MaierSaupe */

            do i=0,M_Bar_B-1
                  tempA1=0.0
                  tempA2=0.0
                  tempA3=0.0
                  if(i<M_Bar_A) then
                  tempA1=qA(i,k,NA)
                  tempA2=qA(i,k,NA-1)
                  tempA3=qA(i,k,NA-2)
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
              tempA2=qA(j,k,NA)
              tempA3=qA(j,k,NA-1)
              endif
              Pr_B(i,k)=Pr_B(i,k) + ds*tempx*(3.0*qB(j,k,1)-3.0*tempA2+tempA3)
           enddo  !enddo index_nonzero
# endif /* MaierSaupe */

            do i=0,M_Bar_B-1
                  tempA2=0.0
                  tempA3=0.0
                  if(i<M_Bar_A) then
                  tempA2=qA(i,k,NA)
                  tempA3=qA(i,k,NA-1)
                  endif
                 

            Pr_B(i,k)=Pr_B(i,k)+3.0*qB(i,k,1)-1.5*tempA2+ &
                      temp*tempA3 - &
                      ds*Wtemp*(3.0*qB(i,k,1)-3.0*tempA2+tempA3)
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
              tempA3=qA(j,k,NA)
              endif
              Pr_B(i,k)=Pr_B(i,k) + ds*tempx*(3.0*qB(j,k,2)-3.0*qB(j,k,1)+tempA3)
           enddo  !enddo index_nonzero
# endif /* MaierSaupe */

            do i=0,M_Bar_B-1
                  tempA3=0.0
                  if(i<M_Bar_A) then
                  tempA3=qA(i,k,NA)
                  endif

            Pr_B(i,k)=Pr_B(i,k)+3.0*qB(i,k,2)-1.5*qB(i,k,1)+ &
                      temp*tempA3 - &
                      ds*Wtemp*(3.0*qB(i,k,2)-3.0*qB(i,k,1)+tempA3)
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
              Pr_B(i,k)=Pr_B(i,k) + ds*tempx*(3.0*qB(j,k,s-NA-1)- & 
              3.0*qB(j,k,s-NA-2)+qB(j,k,s-NA-3))
           enddo  !enddo index_nonzero
# endif /* MaierSaupe */

            do i=0,M_Bar_B-1

            Pr_B(i,k)=Pr_B(i,k)+3.0*qB(i,k,s-NA-1)-1.5*qB(i,k,s-NA-2)+ &
                      temp*qB(i,k,s-NA-3) - &
                      ds*Wtemp*(3.0*qB(i,k,s-NA-1)-3.0*qB(i,k,s-NA-2)+qB(i,k,s-NA-3))
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
            Pk_B(i,k)=local_data(kD(k)%x,kD(k)%y)
# else  /* Dim2 */
            Pk_B(i,k)=local_data(kD(k)%x,kD(k)%y,KD(k)%z)
# endif /* Dim2 */
         enddo
       enddo  !enddo i=0,M_Bar_A-1  


# if defined (OPENMP)
        !$omp parallel default(shared) private(k)
        !$omp do
# endif  /* OPENMP */
         do k=0,LOCAL_SIZE-1
    !     call init_G_mat(G_matB,M_Bar_B,5,kD(k)%z,kD(k)%y,kD(k)%x)
    !     call solver_z(G_matB,Pk_B(:,k),M_Bar_B)

    !     call mkl_dcsrsymv('l', M_Bar_B,  R_GB3_inverse_val(:,k+1), R_GB3_inverse_row(:,k+1), &
    !        R_GB3_inverse_col(:,k+1),Pk_Real_B(:,k), q_temp_RR_B(:,k))

    !     call mkl_dcsrsymv('l', M_Bar_B,  R_GB3_inverse_val(:,k+1), R_GB3_inverse_row(:,k+1), &
    !        R_GB3_inverse_col(:,k+1),Pk_Imag_B(:,k), q_temp_RI_B(:,k))

    !     call mkl_dcsrsymv('l', M_Bar_B,  I_GB3_inverse_val(:,k+1), I_GB3_inverse_row(:,k+1), &
    !        I_GB3_inverse_col(:,k+1),Pk_Real_B(:,k), q_temp_IR_B(:,k))

    !     call mkl_dcsrsymv('l', M_Bar_B,  I_GB3_inverse_val(:,k+1), I_GB3_inverse_row(:,k+1), &
    !        I_GB3_inverse_col(:,k+1),Pk_Imag_B(:,k), q_temp_II_B(:,k))
           Pk_Br(:)=real(Pk_B(:,k))
           Pk_Bi(:)=aimag(Pk_B(:,k))

        call mkl_dcsrsymv('l', M_Bar_B,  R_GB3_inverse_val(:,k+1), R_GB3_inverse_row(:,k+1), &
           R_GB3_inverse_col(:,k+1),Pk_Br, qtmp_Brr)

         call mkl_dcsrsymv('l', M_Bar_B,  R_GB3_inverse_val(:,k+1), R_GB3_inverse_row(:,k+1), &
            R_GB3_inverse_col(:,k+1),Pk_Bi, qtmp_Bri)

         call mkl_dcsrsymv('l', M_Bar_B,  I_GB3_inverse_val(:,k+1), I_GB3_inverse_row(:,k+1), &
            I_GB3_inverse_col(:,k+1),Pk_Br, qtmp_Bir)

         call mkl_dcsrsymv('l', M_Bar_B,  I_GB3_inverse_val(:,k+1), I_GB3_inverse_row(:,k+1), &
            I_GB3_inverse_col(:,k+1),Pk_Bi, qtmp_Bii)

            pk_B(:,k)=dcmplx(qtmp_Brr(:)-qtmp_Bii(:),qtmp_Bri(:)+qtmp_Bir(:)) 
         enddo
# if defined (OPENMP)
        !$omp end  do
        !$omp end  parallel
# endif  /* OPENMP */


             do i=0,M_Bar_B-1
                do k=0,LOCAL_SIZE-1
# if defined (Dim2)
        local_data(kD(k)%x,kD(k)%y)= Pk_B(i,k)
# else  /* Dim2 */
        local_data(kD(k)%x,kD(k)%y,KD(k)%z)= Pk_B(i,k)
# endif /* Dim2 */
                enddo 
       call fftw_mpi_execute_dft(bplan, local_data, local_data)
                do k=0,LOCAL_SIZE-1
# if defined (Dim2)
                 qB(i,k,s-NA)=real(local_data(kD(k)%x,kD(k)%y))*M_grid_inv
# else  /* Dim2 */
                 qB(i,k,s-NA)=real(local_data(kD(k)%x,kD(k)%y,KD(k)%z))*M_grid_inv
# endif /* Dim2 */
                enddo
              enddo 
 
  endif   !endif s>=3:s>NA
                 
endif    !endif s==1

enddo  !enddo s=1,NMAX
!# if defined (Debug)
    call mp_barrier()
    if(myid==0 ) then
     write(*,*) "qB(NB,0,30)",qB(0,30,NB)
     write(*,*) "qB(NB,0,10)",qB(0,10,NB)
     write(*,*) "qB(NB,0,40)",qB(0,40,NB)
     endif
!# endif /* Debug */
!!!!!!!!!!!!!!1


!!!normalization!!!
        do k=0,LOCAL_SIZE-1
         do i=0,M_Bar_B-1
         qB(i,k,0)=0.0
         enddo
       enddo

    do k=0,LOCAL_SIZE-1
       do i=0,min(M_Bar_A-1,M_Bar_B-1)
         qB(i,k,0)=qA(i,k,NA)
        enddo
     enddo

 deallocate(G_matA)
 deallocate(G_matB)

 deallocate(Pr_A)
 deallocate(Pk_A)

 deallocate(Pr_B)
 deallocate(Pk_B)

 deallocate(Pk_Ar)
 deallocate(Pk_Ai)
 deallocate(qtmp_Arr)
 deallocate(qtmp_Ari)
 deallocate(qtmp_Air)
 deallocate(qtmp_Aii)
 deallocate(Pk_Br)
 deallocate(Pk_Bi)
 deallocate(qtmp_BRr)
 deallocate(qtmp_BRi)
 deallocate(qtmp_BIr)
 deallocate(qtmp_BIi)

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
 use G_matrix_mod
 implicit none
 integer :: i,j,s
 integer :: index_nonzero
 integer :: k,KK
 integer :: istat2
 !for A
 real(DP) :: tempA1,tempA2,tempA3,temp,tempx,Wtemp,sumaq,sumaq_global
 real(DP),allocatable :: Pr_A(:,:)
 Complex(DP),allocatable :: G_matA(:,:)
 Complex(DP),allocatable :: G_matB(:,:)
 Complex(DP),allocatable :: Pk_A(:,:)
 !for B
 real(DP),allocatable :: Pr_B(:,:)
 Complex(DP),allocatable :: Pk_B(:,:)
 REAL(DP) :: M_grid_inv

 Real(DP),allocatable :: Pk_Ar(:)
 Real(DP),allocatable :: Pk_Ai(:)
 Real(DP),allocatable :: qtmp_Arr(:)
 Real(DP),allocatable :: qtmp_Ari(:)
 Real(DP),allocatable :: qtmp_Air(:)
 Real(DP),allocatable :: qtmp_Aii(:)
 Real(DP),allocatable :: Pk_Br(:)
 Real(DP),allocatable :: Pk_Bi(:)
 Real(DP),allocatable :: qtmp_Brr(:)
 Real(DP),allocatable :: qtmp_Bri(:)
 Real(DP),allocatable :: qtmp_Bir(:)
 Real(DP),allocatable :: qtmp_Bii(:)
 !!!! just for reference


 allocate(G_matA(0:M_Bar_A-1,0:M_Bar_A-1),stat=istat2)
 allocate(G_matB(0:M_Bar_B-1,0:M_Bar_B-1),stat=istat2)
!for A
 allocate(Pr_A(0:M_Bar_A-1,0:LOCAL_SIZE-1),stat=istat2)
 allocate(Pk_A(0:M_Bar_A-1,0:LOCAL_SIZE-1),stat=istat2)
 allocate(Pk_Ar(0:M_Bar_A-1),stat=istat2)
 allocate(qtmp_Arr(0:M_Bar_A-1),stat=istat2)
 allocate(qtmp_Ari(0:M_Bar_A-1),stat=istat2)
 allocate(Pk_Ai(0:M_Bar_A-1),stat=istat2)
 allocate(qtmp_Air(0:M_Bar_A-1),stat=istat2)
 allocate(qtmp_Aii(0:M_Bar_A-1),stat=istat2)
!for B
 allocate(Pr_B(0:M_Bar_B-1,0:LOCAL_SIZE-1),stat=istat2)
 allocate(Pk_B(0:M_Bar_B-1,0:LOCAL_SIZE-1),stat=istat2)
 allocate(Pk_Br(0:M_Bar_B-1),stat=istat2)
 allocate(qtmp_Brr(0:M_Bar_B-1),stat=istat2)
 allocate(qtmp_Bri(0:M_Bar_B-1),stat=istat2)
 allocate(Pk_Bi(0:M_Bar_B-1),stat=istat2)
 allocate(qtmp_Bir(0:M_Bar_B-1),stat=istat2)
 allocate(qtmp_Bii(0:M_Bar_B-1),stat=istat2)
 if(istat2/=0) then
 write(*,*) "allocate matrix in mde_q failed"
 stop
 endif

M_grid_inv=1.0d0/M_grid

do k=0,LOCAL_SIZE-1
   do i=0,M_Bar_B-1
    qBstar(i,K,NB)=0.0
    if(i==0) qBstar(i,K,NB)=1.0
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
             Pr_B(i,k)=Pr_B(i,k)+ds*tempx*qBstar(j,k,NB)
           enddo
# endif /* MaierSaupe */

            do i=0,M_Bar_B-1
            Pr_B(i,k)=Pr_B(i,k)+(1.0-ds*Wtemp)*qBstar(i,k,NB)
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
            Pk_B(i,k)=local_data(kD(k)%x,kD(k)%y)
# else  /* Dim2 */
            Pk_B(i,k)=local_data(kD(k)%x,kD(k)%y,KD(k)%z)
# endif /* Dim2 */
         enddo
       enddo  !enddo i=0,M_Bar_B-1  


# if defined (OPENMP)
        !$omp parallel default(shared) private(k)
        !$omp do
# endif  /* OPENMP */

         do k=0,LOCAL_SIZE-1
    !     call init_G_mat(G_matB,M_Bar_B,8,kD(k)%z,kD(k)%y,kD(k)%x)
    !     call solver_z(G_matB,Pk_B(:,k),M_Bar_B)

     !    call mkl_dcsrsymv('l', M_Bar_B,  R_GB1_star_inverse_val(:,k+1), R_GB1_star_inverse_row(:,k+1), &
     !       R_GB1_star_inverse_col(:,k+1),Pk_Real_B(:,k), qstar_temp_RR_B(:,k))

     !    call mkl_dcsrsymv('l', M_Bar_B,  R_GB1_star_inverse_val(:,k+1), R_GB1_star_inverse_row(:,k+1), &
     !       R_GB1_star_inverse_col(:,k+1),Pk_Imag_B(:,k), qstar_temp_RI_B(:,k))

     !    call mkl_dcsrsymv('l', M_Bar_B,  I_GB1_star_inverse_val(:,k+1), I_GB1_star_inverse_row(:,k+1), &
     !       I_GB1_star_inverse_col(:,k+1),Pk_Real_B(:,k), qstar_temp_IR_B(:,k))

     !    call mkl_dcsrsymv('l', M_Bar_B,  I_GB1_star_inverse_val(:,k+1), I_GB1_star_inverse_row(:,k+1), &
     !       I_GB1_star_inverse_col(:,k+1),Pk_Imag_B(:,k), qstar_temp_II_B(:,k))

           Pk_Br(:)=real(Pk_B(:,k))
           Pk_Bi(:)=aimag(Pk_B(:,k))

  call mkl_dcsrsymv('l', M_Bar_B,R_GB1_star_inverse_val(:,k+1),R_GB1_star_inverse_row(:,k+1), &
            R_GB1_star_inverse_col(:,k+1),Pk_Br, qtmp_Brr)

  call mkl_dcsrsymv('l', M_Bar_B,R_GB1_star_inverse_val(:,k+1),R_GB1_star_inverse_row(:,k+1), &
            R_GB1_star_inverse_col(:,k+1),Pk_Bi, qtmp_Bri)

  call mkl_dcsrsymv('l', M_Bar_B,I_GB1_star_inverse_val(:,k+1),I_GB1_star_inverse_row(:,k+1), &
            I_GB1_star_inverse_col(:,k+1),Pk_Br, qtmp_Bir)

  call mkl_dcsrsymv('l', M_Bar_B,I_GB1_star_inverse_val(:,k+1),I_GB1_star_inverse_row(:,k+1), &
            I_GB1_star_inverse_col(:,k+1),Pk_Bi, qtmp_Bii)

            pk_B(:,k)=dcmplx(qtmp_Brr(:)-qtmp_Bii(:),qtmp_Bri(:)+qtmp_Bir(:)) 

         enddo
# if defined (OPENMP)
        !$omp end  do
        !$omp end  parallel
# endif  /* OPENMP */


             do i=0,M_Bar_B-1
                do k=0,LOCAL_SIZE-1
# if defined (Dim2)
        local_data(kD(k)%x,kD(k)%y)=Pk_B(i,k)
# else  /* Dim2 */
        local_data(kD(k)%x,kD(k)%y,KD(k)%z)=Pk_B(i,k)
# endif /* Dim2 */
                enddo 
       call fftw_mpi_execute_dft(bplan, local_data, local_data)
                do k=0,LOCAL_SIZE-1
# if defined (Dim2)
                 qBstar(i,k,NB-1)=real(local_data(kD(k)%x,kD(k)%y))*M_grid_inv
# else  /* Dim2 */
                 qBstar(i,k,NB-1)=real(local_data(kD(k)%x,kD(k)%y,KD(k)%z))*M_grid_inv
# endif /* Dim2 */
                enddo
              enddo  

!# if defined (Debug)
    call mp_barrier()
    if(myid==0 ) then
     write(*,*) "qBstar(NB,0,0)",qBstar(0,0,NB)
     write(*,*) "qBstar(NB-1,0,0)",qBstar(0,0,NB-1)
     endif
!# endif /* Debug */
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
             Pr_B(i,k)=Pr_B(i,k)+ds*tempx*(2.0*qBstar(j,k,NB-1)-qBstar(j,k,NB))
           enddo
# endif /* MaierSaupe */

            do i=0,M_Bar_B-1
            Pr_B(i,k)=Pr_B(i,k)+2.0*qBstar(i,k,NB-1)-0.5*qBstar(i,k,NB) - &
                      ds*Wtemp*(2.0*qBstar(i,k,NB-1)-qBstar(i,k,NB))
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
            Pk_B(i,k)=local_data(kD(k)%x,kD(k)%y)
# else  /* Dim2 */
            Pk_B(i,k)=local_data(kD(k)%x,kD(k)%y,KD(k)%z)
# endif /* Dim2 */
         enddo
       enddo  !enddo i=0,M_Bar_B-1  


# if defined (OPENMP)
        !$omp parallel default(shared) private(k)
        !$omp do
# endif  /* OPENMP */
         do k=0,LOCAL_SIZE-1
      !   call init_G_mat(G_matB,M_Bar_B,7,kD(k)%z,kD(k)%y,kD(k)%x)
      !   call solver_z(G_matB,Pk_B(:,k),M_Bar_B)

      !   call mkl_dcsrsymv('l', M_Bar_B,  R_GB2_star_inverse_val(:,k+1), R_GB2_star_inverse_row(:,k+1), &
      !      R_GB2_star_inverse_col(:,k+1),Pk_Real_B(:,k), qstar_temp_RR_B(:,k))

      !   call mkl_dcsrsymv('l', M_Bar_B,  R_GB2_star_inverse_val(:,k+1), R_GB2_star_inverse_row(:,k+1), &
      !      R_GB2_star_inverse_col(:,k+1),Pk_Imag_B(:,k), qstar_temp_RI_B(:,k))

      !   call mkl_dcsrsymv('l', M_Bar_B,  I_GB2_star_inverse_val(:,k+1), I_GB2_star_inverse_row(:,k+1), &
      !      I_GB2_star_inverse_col(:,k+1),Pk_Real_B(:,k), qstar_temp_IR_B(:,k))

      !   call mkl_dcsrsymv('l', M_Bar_B,  I_GB2_star_inverse_val(:,k+1), I_GB2_star_inverse_row(:,k+1), &
      !      I_GB2_star_inverse_col(:,k+1),Pk_Imag_B(:,k), qstar_temp_II_B(:,k))
           Pk_Br(:)=real(Pk_B(:,k))
           Pk_Bi(:)=aimag(Pk_B(:,k))

  call mkl_dcsrsymv('l', M_Bar_B,R_GB2_star_inverse_val(:,k+1),R_GB2_star_inverse_row(:,k+1), &
            R_GB2_star_inverse_col(:,k+1),Pk_Br, qtmp_Brr)

  call mkl_dcsrsymv('l', M_Bar_B,R_GB2_star_inverse_val(:,k+1),R_GB2_star_inverse_row(:,k+1), &
            R_GB2_star_inverse_col(:,k+1),Pk_Bi, qtmp_Bri)

  call mkl_dcsrsymv('l', M_Bar_B,I_GB2_star_inverse_val(:,k+1),I_GB2_star_inverse_row(:,k+1), &
            I_GB2_star_inverse_col(:,k+1),Pk_Br, qtmp_Bir)

  call mkl_dcsrsymv('l', M_Bar_B,I_GB2_star_inverse_val(:,k+1),I_GB2_star_inverse_row(:,k+1), &
            I_GB2_star_inverse_col(:,k+1),Pk_Bi, qtmp_Bii)
            pk_B(:,k)=dcmplx(qtmp_Brr(:)-qtmp_Bii(:),qtmp_Bri(:)+qtmp_Bir(:)) 
         enddo
# if defined (OPENMP)
        !$omp end  do
        !$omp end  parallel
# endif  /* OPENMP */

             do i=0,M_Bar_B-1
                do k=0,LOCAL_SIZE-1
# if defined (Dim2)
      local_data(kD(k)%x,kD(k)%y)=Pk_B(i,k)
# else /* Dim2 */
      local_data(kD(k)%x,kD(k)%y,KD(k)%z)=Pk_B(i,k)
# endif /* Dim2 */
                enddo 
                
       call fftw_mpi_execute_dft(bplan, local_data, local_data)
                do k=0,LOCAL_SIZE-1
# if defined (Dim2)
                 qBstar(i,k,NB-2)=real(local_data(kD(k)%x,kD(k)%y))*M_grid_inv
# else /* Dim2 */
                 qBstar(i,k,NB-2)=real(local_data(kD(k)%x,kD(k)%y,KD(k)%z))*M_grid_inv
# endif /* Dim2 */
                enddo
              enddo  

!# if defined (Debug)
    call mp_barrier()
    if(myid==0 ) then
     write(*,*) "qBstar(NB-2,0,30)",qBstar(0,30,NB-2)
     endif
!# endif /* Debug */
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
             Pr_B(i,k)=Pr_B(i,k)+ds*tempx*(3.0*qBstar(j,k,s-NA+1)- &
                       3.0*qBstar(j,k,s-NA+2)+qBstar(j,k,s-NA+3))
            
           enddo
# endif /* MaierSaupe */
            do i=0,M_Bar_B-1
            Pr_B(i,k)=Pr_B(i,k)+3.0*qBstar(i,k,s-NA+1)-1.5*qBstar(i,k,s-NA+2)+ &
                      temp*qBstar(i,k,s-NA+3) - &
                      ds*Wtemp*(3.0*qBstar(i,k,s-NA+1)-3.0*qBstar(i,k,s-NA+2)+qBstar(i,k,s-NA+3))
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
            Pk_B(i,k)=local_data(kD(k)%x,kD(k)%y)
# else  /* Dim2 */
            Pk_B(i,k)=local_data(kD(k)%x,kD(k)%y,KD(k)%z)
# endif /* Dim2 */
         enddo
       enddo  !enddo i=0,M_Bar_B-1  


# if defined (OPENMP)
        !$omp parallel default(shared) private(k)
        !$omp do
# endif  /* OPENMP */
         do k=0,LOCAL_SIZE-1
       !   call init_G_mat(G_matB,M_Bar_B,6,kD(k)%z,kD(k)%y,kD(k)%x)
       !  call solver_z(G_matB,Pk_B(:,k),M_Bar_B)

       !  call mkl_dcsrsymv('l', M_Bar_B,  R_GB3_star_inverse_val(:,k+1), R_GB3_star_inverse_row(:,k+1), &
       !     R_GB3_star_inverse_col(:,k+1),Pk_Real_B(:,k), qstar_temp_RR_B(:,k))

       !  call mkl_dcsrsymv('l', M_Bar_B,  R_GB3_star_inverse_val(:,k+1), R_GB3_star_inverse_row(:,k+1), &
       !     R_GB3_star_inverse_col(:,k+1),Pk_Imag_B(:,k), qstar_temp_RI_B(:,k))

       !  call mkl_dcsrsymv('l', M_Bar_B,  I_GB3_star_inverse_val(:,k+1), I_GB3_star_inverse_row(:,k+1), &
       !     I_GB3_star_inverse_col(:,k+1),Pk_Real_B(:,k), qstar_temp_IR_B(:,k))

       !  call mkl_dcsrsymv('l', M_Bar_B,  I_GB3_star_inverse_val(:,k+1), I_GB3_star_inverse_row(:,k+1), &
       !     I_GB3_star_inverse_col(:,k+1),Pk_Imag_B(:,k), qstar_temp_II_B(:,k))
           Pk_Br(:)=real(Pk_B(:,k))
           Pk_Bi(:)=aimag(Pk_B(:,k))

  call mkl_dcsrsymv('l', M_Bar_B,R_GB3_star_inverse_val(:,k+1),R_GB3_star_inverse_row(:,k+1), &
            R_GB3_star_inverse_col(:,k+1),Pk_Br, qtmp_Brr)

  call mkl_dcsrsymv('l', M_Bar_B,R_GB3_star_inverse_val(:,k+1),R_GB3_star_inverse_row(:,k+1), &
            R_GB3_star_inverse_col(:,k+1),Pk_Bi, qtmp_Bri)

  call mkl_dcsrsymv('l', M_Bar_B,I_GB3_star_inverse_val(:,k+1),I_GB3_star_inverse_row(:,k+1), &
            I_GB3_star_inverse_col(:,k+1),Pk_Br, qtmp_Bir)

  call mkl_dcsrsymv('l', M_Bar_B,I_GB3_star_inverse_val(:,k+1),I_GB3_star_inverse_row(:,k+1), &
            I_GB3_star_inverse_col(:,k+1),Pk_Bi, qtmp_Bii)
            pk_B(:,k)=dcmplx(qtmp_Brr(:)-qtmp_Bii(:),qtmp_Bri(:)+qtmp_Bir(:)) 
         enddo

# if defined (OPENMP)
        !$omp end  do
        !$omp end  parallel
# endif  /* OPENMP */


             do i=0,M_Bar_B-1
                do k=0,LOCAL_SIZE-1
# if defined (Dim2)
    local_data(kD(k)%x,kD(k)%y)=Pk_B(i,k)
# else  /* Dim2 */
    local_data(kD(k)%x,kD(k)%y,KD(k)%z)=Pk_B(i,k)
# endif /* Dim2 */
                enddo 
       call fftw_mpi_execute_dft(bplan, local_data, local_data)
                do k=0,LOCAL_SIZE-1
# if defined (Dim2)
                 qBstar(i,k,s-NA)=real(local_data(kD(k)%x,kD(k)%y))*M_grid_inv
# else  /* Dim2 */
                 qBstar(i,k,s-NA)=real(local_data(kD(k)%x,kD(k)%y,KD(k)%z))*M_grid_inv
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
             Pr_B(i,k)=Pr_B(i,k) + ds*tempx*(3.0*qBstar(j,k,1)-3.0*qBstar(j,k,2)+qBstar(j,k,3))
           enddo  !enddo index_nonzero
# endif /* MaierSaupe */
            do i=0,min(M_Bar_A-1,M_Bar_B-1)
                Pr_A(i,k)=Pr_B(i,k)
            Pr_A(i,k)=Pr_A(i,k)+3.0*qBstar(i,k,1)-1.5*qBstar(i,k,2)+ &
                      temp*qBstar(i,k,3) - &
                      ds*Wtemp*(3.0*qBstar(i,k,1)-3.0*qBstar(i,k,2)+qBstar(i,k,3))
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
              tempA1=qAstar(j,k,NA)
          
              endif
              Pr_B(i,k)=Pr_B(i,k) + ds*tempx*(3.0*tempA1-3.0*qBstar(j,k,1)+qBstar(j,k,2))
           enddo  !enddo index_nonzero
# endif /* MaierSaupe */
             

            do i=0,min(M_Bar_A-1,M_Bar_B-1)
              Pr_A(i,k)=Pr_B(i,k)

            Pr_A(i,k)=Pr_A(i,k)+3.0*qAstar(i,k,NA)-1.5*qBstar(i,k,1)+ &
                      temp*qBstar(i,k,2) - &
                      ds*Wtemp*(3.0*qAstar(i,k,NA)-3.0*qBstar(i,k,1)+qBstar(i,k,2))
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
              tempA1=qAstar(j,k,NA-1)
              tempA2=qAstar(j,k,NA)
              endif
              Pr_B(i,k)=Pr_B(i,k) + ds*tempx*(3.0*tempA1-3.0*tempA2+qBstar(j,k,1))
           enddo  !enddo index_nonzero
# endif /* MaierSaupe */

            do i=0,min(M_Bar_A-1,M_Bar_B-1)
                Pr_A(i,k)=Pr_B(i,k)

            Pr_A(i,k)=Pr_A(i,k)+3.0*qAstar(i,k,NA-1)-1.5*qAstar(i,k,NA)+ &
                      temp*qBstar(i,k,1)- &
                      ds*Wtemp*(3.0*qAstar(i,k,NA-1)-3.0*qAstar(i,k,NA)+qBstar(i,k,1))
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
       Pr_A(i,k)=Pr_A(i,k) + ds*tempx*(3.0*qAstar(j,k,s+1)-3.0*qAstar(j,k,s+2)+qAstar(j,k,s+3))
           enddo  !enddo index_nonzero
# endif /* MaierSaupe */
            do i=0,M_Bar_A-1
            Pr_A(i,k)=Pr_A(i,k)+3.0*qAstar(i,k,s+1)-1.5*qAstar(i,k,s+2)+ &
                      temp*qAstar(i,k,s+3) - &
                      ds*Wtemp*(3.0*qAstar(i,k,s+1)-3.0*qAstar(i,k,s+2)+qAstar(i,k,s+3))
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
            Pk_A(i,k)=local_data(kD(k)%x,kD(k)%y)
# else  /* Dim2 */
            Pk_A(i,k)=local_data(kD(k)%x,kD(k)%y,KD(k)%z)
# endif /* Dim2 */
         enddo
       enddo  !enddo i=0,M_Bar_A-1  


# if defined (OPENMP)
        !$omp parallel default(shared) private(k)
        !$omp do
# endif  /* OPENMP */
         do k=0,LOCAL_SIZE-1
  !       call init_G_mat(G_matA,M_Bar_A,4,kD(k)%z,kD(k)%y,kD(k)%x)
  !       call solver_z(G_matA,Pk_A(:,k),M_Bar_A)

  !       call mkl_dcsrsymv('l', M_Bar_A,  R_GA3_star_inverse_val(:,k+1), R_GA3_star_inverse_row(:,k+1), &
  !          R_GA3_star_inverse_col(:,k+1),Pk_Real_A(:,k), qstar_temp_RR_A(:,k))

  !       call mkl_dcsrsymv('l', M_Bar_A,  R_GA3_star_inverse_val(:,k+1), R_GA3_star_inverse_row(:,k+1), &
  !          R_GA3_star_inverse_col(:,k+1),Pk_Imag_A(:,k), qstar_temp_RI_A(:,k))

  !       call mkl_dcsrsymv('l', M_Bar_A,  I_GA3_star_inverse_val(:,k+1), I_GA3_star_inverse_row(:,k+1), &
  !          I_GA3_star_inverse_col(:,k+1),Pk_Real_A(:,k), qstar_temp_IR_A(:,k))

  !       call mkl_dcsrsymv('l', M_Bar_A,  I_GA3_star_inverse_val(:,k+1), I_GA3_star_inverse_row(:,k+1), &
  !          I_GA3_star_inverse_col(:,k+1),Pk_Imag_A(:,k), qstar_temp_II_A(:,k))
           Pk_Ar(:)=real(Pk_A(:,k))
           Pk_Ai(:)=aimag(Pk_A(:,k))

  call mkl_dcsrsymv('l', M_Bar_A,R_GA3_star_inverse_val(:,k+1),R_GA3_star_inverse_row(:,k+1), &
            R_GA3_star_inverse_col(:,k+1),Pk_Ar, qtmp_Arr)

  call mkl_dcsrsymv('l', M_Bar_A,R_GA3_star_inverse_val(:,k+1),R_GA3_star_inverse_row(:,k+1), &
            R_GA3_star_inverse_col(:,k+1),Pk_Ai, qtmp_Ari)

  call mkl_dcsrsymv('l', M_Bar_A,I_GA3_star_inverse_val(:,k+1),I_GA3_star_inverse_row(:,k+1), &
            I_GA3_star_inverse_col(:,k+1),Pk_Ar, qtmp_Air)

  call mkl_dcsrsymv('l', M_Bar_A,I_GA3_star_inverse_val(:,k+1),I_GA3_star_inverse_row(:,k+1), &
            I_GA3_star_inverse_col(:,k+1),Pk_Ai, qtmp_Aii)
            pk_A(:,k)=dcmplx(qtmp_Arr(:)-qtmp_Aii(:),qtmp_Ari(:)+qtmp_Air(:)) 
         enddo
# if defined (OPENMP)
        !$omp end  do
        !$omp end  parallel
# endif  /* OPENMP */


             do i=0,M_Bar_A-1
                do k=0,LOCAL_SIZE-1
# if defined (Dim2) 
       local_data(kD(k)%x,kD(k)%y)=Pk_A(i,k)
# else  /* Dim2 */
       local_data(kD(k)%x,kD(k)%y,KD(k)%z)=Pk_A(i,k)
# endif /* Dim2 */
                enddo 
       call fftw_mpi_execute_dft(bplan, local_data, local_data)
                do k=0,LOCAL_SIZE-1
# if defined (Dim2) 
                 qAstar(i,k,s)=real(local_data(kD(k)%x,kD(k)%y))*M_grid_inv
# else  /* Dim2 */
                 qAstar(i,k,s)=real(local_data(kD(k)%x,kD(k)%y,KD(k)%z))*M_grid_inv
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
         qBstar(i,k,0)=0.0
        enddo
     enddo

    do k=0,LOCAL_SIZE-1
       do i=0,min(M_Bar_A-1,M_Bar_B-1)
         qBstar(i,k,0)=qAstar(i,k,NA)
        enddo
     enddo

 deallocate(G_matA)
 deallocate(G_matB)
 deallocate(Pr_A)
 deallocate(Pk_A)
!
 deallocate(Pr_B)
 deallocate(Pk_B)

 deallocate(Pk_Ar)
 deallocate(Pk_Ai)
 deallocate(qtmp_Arr)
 deallocate(qtmp_Ari)
 deallocate(qtmp_Air)
 deallocate(qtmp_Aii)
 deallocate(Pk_Br)
 deallocate(Pk_Bi)
 deallocate(qtmp_BRr)
 deallocate(qtmp_BRi)
 deallocate(qtmp_BIr)
 deallocate(qtmp_BIi)
# if defined (Debug)
    call mp_barrier()
    if(myid==0 ) then
     write(*,*) "exit MDE_qstar"
     endif
# endif /* Debug */

end subroutine MDE_qstar
end module  wormlike_MDE 

# else /* SLAB */
# endif /* SLAB */
