# if defined (SLAB)
Module wormlike_MDE_ABC
    USE nrtype,only :DP,PI
    REAL(DP) :: M_grid_inv
    real(dp) :: t1,t2,t3,t4,t32,t41
contains
   
   subroutine ABC_melt_driver()
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
    real(DP) :: tempA1,tempA2,tempA3,temp,tempx,Wtemp,sumaq,sumaq_global,sumaq1
    real(DP) :: sumPA,sumPA_global
    real(DP),allocatable :: Pr(:,:)
    integer :: G_index,dir,blk_indx
    max_M_Bar =max(M_Bar_A,M_Bar_B,M_Bar_C)
 
    allocate(Pr(0:max_M_Bar-1,0:LOCAL_SIZE-1),stat=istat2)
    if(istat2/=0) then
     write(*,*) "allocate matrix in mde_q failed"
    stop
    endif

    M_grid_inv=1.0d0/M_grid

   do k=0,LOCAL_SIZE-1
     do i=0,max_M_Bar-1
      q_blk(i,K,0,1)=0.0
      if(i==0) q_blk(i,K,0,1)=1.0
     enddo
   enddo
 write(*,*) "init q forward done"
 !!!!!
do blk_indx=1,3
 do s=1+Nstart_blk(blk_indx),Ns_blk(blk_indx)+Nstart_blk(blk_indx)
     call init_pr(s,blk_indx,1,max_M_Bar,Pr)
     if(s<=Ns_blk(1))
     G_index=min(s,3)
     else if (s>Ns_blk(1) .and. s<=Ns_blk(2)+Nstart_blk(2))
     G_index=4
     else
     G_index=5
     endif
    call mde_one_step(max_M_Bar,q_f(:,:,s),Pr,G_index,blk_indx,1)
 enddo
enddo  ! enddo blk_indx
 write(*,*) "solve  qf forward done"
  call mp_barrier()
! solve the backward propagator 
do k=0,LOCAL_SIZE-1
   do i=0,max_M_Bar-1
   q_star(i,K,Nmax)=0.0
   if(i==0) q_star(i,K,Nmax)=1.0
   enddo
enddo

do blk_indx=3,1,-1
 do s=Nstart_blk(blk_indx)+Ns_blk(blk_indx)-1,Nstart_blk(blk_indx),-1! [Nmax-1,NA+NB],[NA+NB-1,NA],[NA-1,0]
     call init_pr(s,blk_indx,-1,max_M_Bar,Pr)
     if(s>=Nstart_blk(3))
     G_index=11-min(Nmax-s,3)
     else if (s>=Ns_blk(1) .and. s<=Ns_blk(2)+Nstart(2)-1)
     G_index=7
     else
     G_index=6
     endif
    call mde_one_step(max_M_Bar,q_star(:,:,s),Pr,G_index,blk_indx,-1)
 enddo
enddo  ! enddo blk_indx

 write(*,*) "solve  qBstar and qAstar done"
  call mp_barrier()
  deallocate(Pr)

   end subroutine ABC_melt_driver

 subroutine mde_one_step(M_Bar,q_out,Pr,G_index,blk_indx,dir)
 USE nrtype,only :DP,PI
 USE global_para
 USE mpi
 USE control
 USE constants
 USE utility
 USE mmpi
 USE mpi_fftw3_operation
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
 integer,intent(in) :: G_index,dir,M_Bar,blk_indx
 Complex(DP)              :: G_mat(0:M_Bar-1,0:M_Bar-1)
 Complex(DP)              :: G_mat1(0:M_Bar-1,0:M_Bar-1)
 real(DP),intent(in)    :: Pr(0:M_Bar-1,0:LOCAL_SIZE-1)
 Complex(DP)               :: Pk(0:M_Bar-1,0:LOCAL_SIZE-1)
 real(DP),intent(inout)    :: q_out(0:M_Bar-1,0:LOCAL_SIZE-1)
 real(dp) :: error
 real(DP)               :: pk_r(0:M_Bar-1)
 real(DP)               :: pk_i(0:M_Bar-1)
 real(DP)               :: qtmp_rr(0:M_Bar-1)
 real(DP)               :: qtmp_ri(0:M_Bar-1)
 real(DP)               :: qtmp_ir(0:M_Bar-1)
 real(DP)               :: qtmp_ii(0:M_Bar-1)

        t1=MPI_Wtime()
!!!doing the FFT
      do i=0,M_Bar-1
        do k=0,LOCAL_SIZE-1
# if defined (Dim2)
          local_data(kD(k)%x,kD(k)%y)=cmplx(Pr(i,k),0.0)
# else  /* Dim2 */
          local_data(kD(k)%x,kD(k)%y,KD(k)%z)=cmplx(Pr(i,k),0.0)
# endif /* Dim2 */
         enddo
       
       call fftw_mpi_execute_dft(fplan, local_data, local_data)
       
        do k=0,LOCAL_SIZE-1
# if defined (Dim2)
            Pk(i,k)=local_data(kD(k)%x,kD(k)%y)
# else  /* Dim2 */
            Pk(i,k)=local_data(kD(k)%x,kD(k)%y,KD(k)%z)
# endif /* Dim2 */
         enddo
       enddo  !enddo i=0,M_Bar_A-1  
        t2=MPI_Wtime()

   if(G_index <=2 .or. G_index >=9) then 
     do k=0,LOCAL_SIZE-1
# if defined (Dim2)
       call init_G_mat(G_mat,M_Bar,G_index,kD(k)%y,kD(k)%x)
# else  /* Dim2 */
       call init_G_mat(G_mat,M_Bar,G_index,kD(k)%z,kD(k)%y,kD(k)%x)
# endif /* Dim2 */
       call solver_z(G_mat,Pk(:,k),M_Bar)
     enddo

   else if(G_index>=3 .and. G_index<=5) then
     do k=0,LOCAL_SIZE-1
        do i=0,M_Bar-1
            Pk_r(i)=real(Pk(i,k)) 
            Pk_i(i)=aimag(Pk(i,k)) 
        enddo
        call mkl_dcsrsymv('l', M_Bar,  R_G3_inv_val(:,k+1,blk_indx), R_G3_inv_row(:,k+1,blk_indx), &
           R_G3_inv_col(:,k+1,blk_indx),Pk_r, qtmp_rr)

         call mkl_dcsrsymv('l', M_Bar, R_G3_inv_val(:,k+1,blk_indx), R_G3_inv_row(:,k+1,blk_indx), &
            R_G3_inv_col(:,k+1,blk_indx),Pk_i, qtmp_ri)

         call mkl_dcsrsymv('l', M_Bar, I_G3_inv_val(:,k+1,blk_indx), I_G3_inv_row(:,k+1,blk_indx), &
            I_G3_inv_col(:,k+1,blk_indx),Pk_r, qtmp_ir)

         call mkl_dcsrsymv('l', M_Bar, I_G3_inv_val(:,k+1,blk_indx), I_G3_inv_row(:,k+1,blk_indx), &
            I_G3_inv_col(:,k+1,blk_indx),Pk_i, qtmp_ii)

            pk(:,k)=dcmplx(qtmp_rr(:)-qtmp_ii(:),qtmp_ri(:)+qtmp_ir(:)) 
         enddo

   else if(G_index>=6 .and. G_index<=8) then
     do k=0,LOCAL_SIZE-1
        do i=0,M_Bar-1
            Pk_r(i)=real(Pk(i,k)) 
            Pk_i(i)=aimag(Pk(i,k)) 
        enddo
        call mkl_dcsrsymv('l', M_Bar,  R_G3star_inv_val(:,k+1,blk_indx), R_G3star_inv_row(:,k+1,blk_indx), &
           R_G3star_inv_col(:,k+1,blk_indx),Pk_r, qtmp_rr)

         call mkl_dcsrsymv('l', M_Bar, R_G3star_inv_val(:,k+1,blk_indx), R_G3star_inv_row(:,k+1,blk_indx), &
            R_G3star_inv_col(:,k+1,blk_indx),Pk_i, qtmp_ri)

         call mkl_dcsrsymv('l', M_Bar, I_G3star_inv_val(:,k+1,blk_indx), I_G3star_inv_row(:,k+1,blk_indx), &
            I_G3star_inv_col(:,k+1,blk_indx),Pk_r, qtmp_ir)

         call mkl_dcsrsymv('l', M_Bar, I_G3star_inv_val(:,k+1,blk_indx), I_G3star_inv_row(:,k+1,blk_indx), &
            I_G3star_inv_col(:,k+1,blk_indx),Pk_i, qtmp_ii)

            pk(:,k)=dcmplx(qtmp_rr(:)-qtmp_ii(:),qtmp_ri(:)+qtmp_ir(:)) 
         enddo

   endif 

        t3=MPI_Wtime()
     do i=0,M_Bar-1
                do k=0,LOCAL_SIZE-1
# if defined (Dim2)
        local_data(kD(k)%x,kD(k)%y)=Pk(i,k)
# else  /* Dim2 */
        local_data(kD(k)%x,kD(k)%y,KD(k)%z)=Pk(i,k)
# endif  /* Dim2 */
                enddo
       call fftw_mpi_execute_dft(bplan, local_data, local_data)
            
                do k=0,LOCAL_SIZE-1
# if defined (Dim2)
                 q_out(i,k)=real(local_data(kD(k)%x,kD(k)%y))*M_grid_inv
# else  /* Dim2 */
                 q_out(i,k)=real(local_data(kD(k)%x,kD(k)%y,KD(k)%z))*M_grid_inv
# endif  /* Dim2 */
                enddo
              enddo  
        t4=MPI_Wtime()

 end subroutine mde_one_step

 subroutine init_pr(s,blk_indx, dir,M_Bar,Pr)
 USE nrtype,only :DP,PI
 USE global_para
 USE mpi
 USE control
 USE constants
 USE utility
 USE mmpi
 USE mpi_fftw3_operation
 USE matrix_inverse
 USE G_matrix_mod
 implicit none
 integer :: i,j
 integer :: index_nonzero
 integer :: k
 integer :: istat2
 integer,intent(in) :: s,dir,M_Bar
 real(DP) :: tempA1,tempA2,tempA3,temp,tempx,Wtemp,sumaq,sumaq_global,sumaq1
 real(DP) :: sumPA,sumPA_global
 real(DP),intent(inout)    :: Pr(0:M_Bar-1,0:LOCAL_SIZE-1)
 integer :: s_block, s_BDF

if(dir==1) then

   if(s==1) then
          Pr=0.0
      do k=0,LOCAL_SIZE-1
          Wtemp=W_blk(k,blk_indx)
# if defined (MaierSaupe)
          do index_nonzero=THETA_nonzero_1D_indexK_blk(k-1,1)+1, &
             THETA_nonzero_1D_indexK_blk(k,1)
             i=THETA_nonzero_1D_blk(index_nonzero,1)%i
             j=THETA_nonzero_1D_blk(index_nonzero,1)%j
             tempx=THETA_nonzero_1D_blk(index_nonzero,1)%value
             Pr(i,k)=Pr(i,k)+ds*tempx*q_f(j,k,s-1)
           enddo
# endif /* MaierSaupe */
            do i=0,M_Bar-1
            Pr(i,k)=Pr(i,k)+(1.0-ds*Wtemp)*q_f(i,k,s-1)
            enddo
        enddo !enddo k=0,LOCAL_SIZE  
   endif

   if(s==2) then
          Pr=0.0
      do k=0,LOCAL_SIZE-1
          Wtemp=W_blk(k,blk_indx)
# if defined (MaierSaupe)
          do index_nonzero=THETA_nonzero_1D_indexK_blk(k-1,1)+1, &
             THETA_nonzero_1D_indexK_blk(k,1)
             i=THETA_nonzero_1D_blk(index_nonzero,1)%i
             j=THETA_nonzero_1D_blk(index_nonzero,1)%j
             tempx=THETA_nonzero_1D_blk(index_nonzero,1)%value
             Pr(i,k)=Pr(i,k)+ds*tempx*(2.0*q_f(j,k,s-1)-q_f(j,k,s-2))
           enddo
# endif /* MaierSaupe */
            do i=0,M_Bar-1
            Pr(i,k)=Pr(i,k)+2.0*q_f(i,k,s-1)-0.5*q_f(i,k,s-2) - &
                      ds*Wtemp*(2.0*q_f(i,k,s-1)-q_f(i,k,s-2))
            enddo
        enddo !enddo k=0,LOCAL_SIZE    
   endif

   else if(s>=3) then
          temp=1.0/3.0
          Pr=0.0
      do k=0,LOCAL_SIZE-1
          Wtemp=W_blk(k,blk_indx)
# if defined (MaierSaupe)
          do index_nonzero=THETA_nonzero_1D_indexK_blk(k-1,1)+1, &
             THETA_nonzero_1D_indexK_blk(k,1)
             i=THETA_nonzero_1D_blk(index_nonzero,1)%i
             j=THETA_nonzero_1D_blk(index_nonzero,1)%j
             tempx=THETA_nonzero_1D_blk(index_nonzero,1)%value
             Pr(i,k)=Pr(i,k)+ds*tempx*(3.0*q_f(j,k,s-1)- &
                       3.0*q_f(j,k,s-2)+q_f(j,k,s-3))
            
           enddo
# endif /* MaierSaupe */
            do i=0,M_Bar-1
            Pr(i,k)=Pr(i,k)+3.0*q_f(i,k,s-1)-1.5*q_f(i,k,s-2)+ &
                      temp*q_f(i,k,s-3) - &
           ds*Wtemp*(3.0*q_f(i,k,s-1)-3.0*q_f(i,k,s-2)+q_f(i,k,s-3))
            enddo
        enddo !enddo k=0,LOCAL_SIZE    

 endif    !endif s>=3

else if(dir==-1) then

   if(s==Nmax-1) then
          Pr=0.0
      do k=0,LOCAL_SIZE-1
          Wtemp=W_blk(k,blk_indx)
# if defined (MaierSaupe)
          do index_nonzero=THETA_nonzero_1D_indexK_blk(k-1,1)+1, &
             THETA_nonzero_1D_indexK_blk(k,1)
             i=THETA_nonzero_1D_blk(index_nonzero,1)%i
             j=THETA_nonzero_1D_blk(index_nonzero,1)%j
             tempx=THETA_nonzero_1D_blk(index_nonzero,1)%value
             Pr(i,k)=Pr(i,k)+ds*tempx*q_star(j,k,Nmax)
           enddo
# endif /* MaierSaupe */
            do i=0,M_Bar-1
            Pr(i,k)=Pr(i,k)+(1.0-ds*Wtemp)*q_star(i,k,Nmax)
            enddo
        enddo !enddo k=0,LOCAL_SIZE  
   endif

   if(s==Nmax-2) then
          Pr=0.0
      do k=0,LOCAL_SIZE-1
          Wtemp=W_blk(k,blk_indx)
# if defined (MaierSaupe)
          do index_nonzero=THETA_nonzero_1D_indexK_blk(k-1,1)+1, &
             THETA_nonzero_1D_indexK_blk(k,1)
             i=THETA_nonzero_1D_blk(index_nonzero,1)%i
             j=THETA_nonzero_1D_blk(index_nonzero,1)%j
             tempx=THETA_nonzero_1D_blk(index_nonzero,1)%value
             Pr(i,k)=Pr(i,k)+ds*tempx*(2.0*q_star(j,k,Nmax-1)-q_star(j,k,Nmax))
           enddo
# endif /* MaierSaupe */
            do i=0,M_Bar-1
            Pr(i,k)=Pr(i,k)+2.0*q_star(i,k,Nmax-1)-0.5*q_star(i,k,Nmax) - &
                      ds*Wtemp*(2.0*q_star(i,k,Nmax-1)-q_star(i,k,Nmax))
            enddo
        enddo !enddo k=0,LOCAL_SIZE    
   endif

   else if(s<=Nmax-3) then
          temp=1.0/3.0
          Pr=0.0
      do k=0,LOCAL_SIZE-1
          Wtemp=W_blk(k,blk_indx)
# if defined (MaierSaupe)
          do index_nonzero=THETA_nonzero_1D_indexK_blk(k-1,1)+1, &
             THETA_nonzero_1D_indexK_blk(k,1)
             i=THETA_nonzero_1D_blk(index_nonzero,1)%i
             j=THETA_nonzero_1D_blk(index_nonzero,1)%j
             tempx=THETA_nonzero_1D_blk(index_nonzero,1)%value
             Pr(i,k)=Pr(i,k)+ds*tempx*(3.0*q_star(j,k,s+1)- &
                       3.0*q_star(j,k,s+2)+q_star(j,k,s+3))
            
           enddo
# endif /* MaierSaupe */
            do i=0,M_Bar-1
            Pr(i,k)=Pr(i,k)+3.0*q_star(i,k,s+1)-1.5*q_star(i,k,s+2)+ &
                      temp*q_star(i,k,s+3) - &
           ds*Wtemp*(3.0*q_star(i,k,s+1)-3.0*q_star(i,k,s+2)+q_star(i,k,s+3))
            enddo
        enddo !enddo k=0,LOCAL_SIZE    

 endif    !endif s

else 
write(*,*) "illegal input of dir",dir
stop

endif  ! endif dir

 end subroutine init_pr

end module  wormlike_MDE_ABC 

# else /* SLAB */
# endif /* SLAB */
