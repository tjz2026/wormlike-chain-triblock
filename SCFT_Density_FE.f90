# if defined (SPH)
# else /* SPH */
       subroutine density
       USE nrtype,only :DP
       USE global_para
       USE mpi
       USE control
       USE constants
       USE utility
       USE mmpi
# if defined (SLAB)
       USE mpi_fftw3_operation
# else /* SLAB */
 use decomp_fft_mod
# endif /* SLAB */
       implicit none
       integer :: i,j,s,K
       real(DP) :: pff_temp_global,pff_temp,totden_global
       real(DP) :: RAtot,RAtot_global
       REAL(DP),allocatable :: ar(:)
       REAL(DP),allocatable :: br(:)
       REAL(DP),allocatable :: c_A(:)
       REAL(DP),allocatable :: c_B(:)
       REAL(DP),allocatable :: denz(:)
       REAL(DP) :: sum_iA,sum_iB,sum_half,sum_end_B
       REAL(DP) :: nor_coeff
       REAL(DP) :: sumbr,sumbr_global
       REAL(DP),allocatable :: c11_A(:)
       REAL(DP),allocatable :: c12_A(:)
       REAL(DP),allocatable :: c22_A(:)
       REAL(DP),allocatable :: c13_A(:)
       REAL(DP),allocatable :: c23_A(:)

       REAL(DP),allocatable :: c11_B(:)
       REAL(DP),allocatable :: c12_B(:)
       REAL(DP),allocatable :: c22_B(:)
       REAL(DP),allocatable :: c13_B(:)
       REAL(DP),allocatable :: c23_B(:)

!# if defined (Debug)
    call mp_barrier()
    if(myid==0 ) then
        write(*,*) "enter density"
    endif
    call mp_barrier()
!# endif /* Debug */

       allocate(ar(0:LOCAL_SIZE-1))
       do K=0,LOCAL_SIZE-1
         ar(K)=qB(0,K,NB)
       enddo

!# if defined (Debug)
    call mp_barrier()
       if(myid==0 ) then
       write(*,*) "M_v",M_v,"on",myid
       endif
       RAtot=sum(ar)
       RAtot_global=0.0 
    call mp_barrier()
!# endif  /* Debug */

# if defined (SLAB)       
# if defined (Dim2)       
       pff_temp=simposon_2D_1D_mpi(SIDEy,local_nx,dy,dx,ar)
# else /* Dim2 */
       pff_temp=simposon_3D_1D_mpi(SIDEz,SIDEy,local_nx,dz,dy,dx,ar)
# endif  /* Dim2 */
# else /* SLAB */
       pff_temp=simposon_3D_1D_mpi(xsize(1),xsize(2),xsize(3),dx,dy,dz,ar)
# endif /* SLAB */

              pff_temp_global=0.0
              call mp_barrier()
              call mp_allreduce(pff_temp,pff_temp_global)
# if defined (Debug) 
              call mp_allreduce(RAtot,RAtot_global)
# endif /* Debug */
              call mp_barrier()

        pff_global=pff_temp_global/M_v

!# if defined (Debug) 
        if(myid==0) then
         write(*,*) "arsum",RAtot_global,"pff before log",pff_global
         write(*,*) "pff_global ==",dlog(pff_global)
         endif    
              call mp_barrier()
!# endif /* Debug */
         
        pff_global=dlog(pff_global)

        deallocate(ar)
        allocate(c_A(0:NA))          
        allocate(c_B(0:NB)) 
        allocate(denz(0:LOCAL_SIZE-1))
         call mp_barrier()

       do K=0,LOCAL_SIZE-1
          sum_iA=0.0
          sum_iB=0.0
         do i=0,M_Bar_A-1
            do s=0,NA
              c_A(s)=qA(i,K,s)*qAstar(i,K,s)
            enddo
          RHOA(i,K)=simposon_1D_NR(0,NA,ds,c_A)
          sum_iA=sum_iA+RHOA(i,K)
         enddo

         do i=0,M_Bar_B-1
            do s=0,NB
              c_B(s)=qB(i,K,s)*qBstar(i,K,s)
            enddo
          RHOB(i,K)=simposon_1D_NR(0,NB,ds,c_B)
          sum_iB=sum_iB+RHOB(i,K)
         enddo
         denz(K)=sum_iA+sum_iB
       enddo
     
# if defined (SLAB)       
# if defined (Dim2)       
     totden=simposon_2D_1D_mpi(SIDEy,local_nx,dy,dx,denz)
# else  /* Dim2 */
     totden=simposon_3D_1D_mpi(SIDEz,SIDEy,local_nx,dz,dy,dx,denz)
# endif /* Dim2 */
# else /* SLAB */
     totden=simposon_3D_1D_mpi(xsize(1),xsize(2),xsize(3),dx,dy,dz,denz)
# endif /* SLAB */

              totden_global=0.0
              call mp_barrier()
              call mp_allreduce(totden,totden_global)
              call mp_barrier()
        totden=totden_global
         if(myid==0) then
         write(*,*) "totden ==",totden
         endif
     do K=0,LOCAL_SIZE-1
         sum_iA=0.0
         sum_iB=0.0
         sum_half=0.0
         sum_end_B=0.0
      do i=0,M_Bar_A-1
          RHOA(i,K)=M_v*RHOA(i,K)/totden
          sum_iA=sum_iA+RHOA(i,K)
          sum_half=sum_half+qA(i,K,NA)*qAstar(i,K,NA)
      enddo

      do i=0,M_Bar_B-1
          RHOB(i,K)=M_v*RHOB(i,K)/totden
          sum_iB=sum_iB+RHOB(i,K)
          sum_end_B=sum_end_B+qB(i,K,NB)*qBstar(i,K,NB)
      enddo

        RA(k)=sum_iA
        RB(K)=sum_iB
        R_half(K)=M_v*sum_half/totden
        R_end(K)=M_v*sum_end_B/totden
      enddo

!!! for S 
        nor_coeff=M_v/totden
        allocate(c11_A(0:NA))
        allocate(c12_A(0:NA))
        allocate(c13_A(0:NA))
        allocate(c23_A(0:NA))
        allocate(c22_A(0:NA))

        allocate(c11_B(0:NB))
        allocate(c12_B(0:NB))
        allocate(c13_B(0:NB))
        allocate(c23_B(0:NB))
        allocate(c22_B(0:NB))

    do k=0,LOCAL_SIZE-1
       do s=0,NA
  C11_A(s)=Sum_sparse_2D(M_Bar_A,J11ij_nonzero_1D_A,J11_A,qA(:,K,s),qAstar(:,K,s))
  C22_A(s)=Sum_sparse_2D(M_Bar_A,J22ij_nonzero_1D_A,J22_A,qA(:,K,s),qAstar(:,K,s))
  C12_A(s)=Sum_sparse_2D(M_Bar_A,J12ij_nonzero_1D_A,J12_A,qA(:,K,s),qAstar(:,K,s))
  C23_A(s)=Sum_sparse_2D(M_Bar_A,J23ij_nonzero_1D_A,J23_A,qA(:,K,s),qAstar(:,K,s))
  C13_A(s)=Sum_sparse_2D(M_Bar_A,J13ij_nonzero_1D_A,J13_A,qA(:,K,s),qAstar(:,K,s))
      enddo

       do s=0,NB
  C11_B(s)=Sum_sparse_2D(M_Bar_B,J11ij_nonzero_1D_B,J11_B,qB(:,K,s),qBstar(:,K,s))
  C12_B(s)=Sum_sparse_2D(M_Bar_B,J12ij_nonzero_1D_B,J12_B,qB(:,K,s),qBstar(:,K,s))
  C13_B(s)=Sum_sparse_2D(M_Bar_B,J13ij_nonzero_1D_B,J13_B,qB(:,K,s),qBstar(:,K,s))
  C22_B(s)=Sum_sparse_2D(M_Bar_B,J22ij_nonzero_1D_B,J22_B,qB(:,K,s),qBstar(:,K,s))
  C23_B(s)=Sum_sparse_2D(M_Bar_B,J23ij_nonzero_1D_B,J23_B,qB(:,K,s),qBstar(:,K,s))
       enddo
         
      SA_OP(K,0,0)=simposon_1D_NR(0,NA,ds,c11_A)*nor_coeff
      SA_OP(K,1,1)=simposon_1D_NR(0,NA,ds,c22_A)*nor_coeff
      SA_OP(K,2,2)=-(SA_OP(K,0,0)+SA_OP(K,1,1))
      SA_OP(K,0,1)=simposon_1D_NR(0,NA,ds,c12_A)*nor_coeff
      SA_OP(K,0,2)=simposon_1D_NR(0,NA,ds,c13_A)*nor_coeff
      SA_OP(K,1,2)=simposon_1D_NR(0,NA,ds,c23_A)*nor_coeff
      SA_OP(K,1,0)=SA_OP(K,0,1)
      SA_OP(K,2,0)=SA_OP(K,0,2)
      SA_OP(K,2,1)=SA_OP(K,1,2)

      SB_OP(K,0,0)=simposon_1D_NR(0,NB,ds,c11_B)*nor_coeff
      SB_OP(K,1,1)=simposon_1D_NR(0,NB,ds,c22_B)*nor_coeff
      SB_OP(K,2,2)=-(SB_OP(K,0,0)+SB_OP(K,1,1))
      SB_OP(K,0,1)=simposon_1D_NR(0,NB,ds,c12_B)*nor_coeff
      SB_OP(K,0,2)=simposon_1D_NR(0,NB,ds,c13_B)*nor_coeff
      SB_OP(K,1,2)=simposon_1D_NR(0,NB,ds,c23_B)*nor_coeff
      SB_OP(K,1,0)=SB_OP(K,0,1)
      SB_OP(K,2,0)=SB_OP(K,0,2)
      SB_OP(K,2,1)=SB_OP(K,1,2)
    
       S_OP(K,0,0)=SA_OP(K,0,0)+SB_OP(K,0,0)
       S_OP(K,1,1)=SA_OP(K,1,1)+SB_OP(K,1,1)
       S_OP(K,2,2)=SA_OP(K,2,2)+SB_OP(K,2,2)
       S_OP(K,0,1)=SA_OP(K,0,1)+SB_OP(K,0,1)
       S_OP(K,0,2)=SA_OP(K,0,2)+SB_OP(K,0,2)
       S_OP(K,1,2)=SA_OP(K,1,2)+SB_OP(K,1,2)
       S_OP(K,1,0)=S_OP(K,0,1)  
       S_OP(K,2,0)=S_OP(K,0,2)  
       S_OP(K,2,1)=S_OP(K,1,2)  

    !!!warning: the follwing code requires that M_Bar_A must equal to M_Bar_B ,needs to be modified later.
    !   do s=0,NA
    !   C11_A(s)=Sum_sparse_2D(M_Bar_A,Rx_nonzero_1D,Rx_num,qA(s,:,K),qAstar(s,:,K))
    !   C12_A(s)=Sum_sparse_2D(M_Bar_A,Ry_nonzero_1D,Ry_num,qA(s,:,K),qAstar(s,:,K))
    !   C13_A(s)=Sum_sparse_2D(M_Bar_A,Rz_nonzero_1D,Rz_num,qA(s,:,K),qAstar(s,:,K))
    !   enddo

    !  Pu_A(K,1)=simposon_1D_NR(0,NA,ds,c11_A)*nor_coeff
    !  Pu_A(K,2)=simposon_1D_NR(0,NA,ds,c12_A)*nor_coeff
    !  Pu_A(K,3)=simposon_1D_NR(0,NA,ds,c13_A)*nor_coeff

    !   do s=0,NB
    !   C11_B(s)=Sum_sparse_2D(M_Bar_B,Rx_nonzero_1D,Rx_num,qB(s,:,K),qBstar(s,:,K))
    !   C12_B(s)=Sum_sparse_2D(M_Bar_B,Ry_nonzero_1D,Ry_num,qB(s,:,K),qBstar(s,:,K))
    !   C13_B(s)=Sum_sparse_2D(M_Bar_B,Rz_nonzero_1D,Rz_num,qB(s,:,K),qBstar(s,:,K))
    !   enddo

    !  Pu_B(K,1)=simposon_1D_NR(0,NB,ds,c11_B)*nor_coeff
    !  Pu_B(K,2)=simposon_1D_NR(0,NB,ds,c12_B)*nor_coeff
    !  Pu_B(K,3)=simposon_1D_NR(0,NB,ds,c13_B)*nor_coeff

    !  Pu(K,1)=Pu_A(K,1)+Pu_B(K,1) 
    !  Pu(K,2)=Pu_A(K,2)+Pu_B(K,2) 
    !  Pu(K,3)=Pu_A(K,3)+Pu_B(K,3) 
   enddo

  deallocate(c_A)
  deallocate(c_B)
  deallocate(denz)
  deallocate(c11_A)
  deallocate(c22_A)
  deallocate(c12_A)
  deallocate(c13_A)
  deallocate(c23_A)
  deallocate(c11_B)
  deallocate(c22_B)
  deallocate(c12_B)
  deallocate(c13_B)
  deallocate(c23_B)

# if defined (Debug)
   call mp_barrier()
   if(myid==0) then
   write(*,*) "RA(0),RA(30),RA(90)=",RA(0),RA(30),RA(90)
   write(*,*) "RB(0),RB(30),RB(90)=",RB(0),RB(30),RB(90)
   endif
# endif /* Debug */

   call mp_barrier()    
   end subroutine density

   subroutine free_energy()
   USE nrtype,only :DP
       USE global_para
       USE mpi
       USE control
       USE constants
       USE utility
       USE mmpi
# if defined (SLAB)
       USE mpi_fftw3_operation
# else /* SLAB */
 use decomp_fft_mod
# endif /* SLAB */
    implicit none
    integer :: j,i,K
    real*8 :: tempEE
    real*8,allocatable :: ar(:)
    real*8,allocatable :: b(:,:) 


   allocate(ar(0:LOCAL_SIZE-1))
   allocate( b(0:N_dim_ddm-1,0:N_dim_ddm-1)) 
    
    do K=0,LOCAL_SIZE-1
       ar(K)=NXab*RA(K)*RB(K)-WA(K)*RA(K)-WB(K)*RB(K)+0.5* &
             (WA(K)+WB(K))*(RA(K)+RB(K)-1.0) 


      if(abs(NMu-0.0)<1.0e-5) then
         ar(K)=ar(K)
       else
           do j=0,N_dim_ddm-1
             do i=0,N_dim_ddm-1
               b(i,j)=M_OP(K,i,j)
              enddo
           enddo
       ar(K)=ar(K)+0.5*(double_dot_multi(b,b,N_dim_ddm)/NMu) 


      endif  
  enddo

# if defined (SLAB)       
# if defined (Dim2)   
   tempEE=simposon_2D_1D_mpi(SIDEy,local_nx,dy,dx,ar)
# else /* Dim2 */
   tempEE=simposon_3D_1D_mpi(SIDEz,SIDEy,local_nx,dz,dy,dx,ar)
# endif /* Dim2 */
# else /* SLAB */
   tempEE=simposon_3D_1D_mpi(xsize(1),xsize(2),xsize(3),dx,dy,dz,ar)
# endif /* SLAB */


   tempEE_global=0.0
              call mp_barrier()
              call mp_allreduce(tempEE,tempEE_global)
              call mp_barrier()
        tempEE_global=tempEE_global/M_v
        tempEE=tempEE_global

 FE_global=-pff_global+tempEE_global

 if(myid==0) then
 write(*,*) "FREE energy=",FE_global,pff_global,tempEE_global
 endif

     deallocate(ar)
     deallocate(b)
   call mp_barrier()    
end subroutine free_energy    
     
   subroutine free_energy_finite()
   USE nrtype,only :DP
       USE global_para
       USE mpi
       USE control
       USE constants
       USE utility
       USE mmpi
# if defined (SLAB)
       USE mpi_fftw3_operation
# else /* SLAB */
 use decomp_fft_mod
# endif /* SLAB */
    implicit none
    integer :: j,i,K
    real*8 :: tempEE
    real*8,allocatable :: ar(:)
    real*8,allocatable :: b(:,:) 


   allocate(ar(0:LOCAL_SIZE-1))
   allocate( b(0:N_dim_ddm-1,0:N_dim_ddm-1)) 
    
        call ar_from_finite(ar)
    do K=0,LOCAL_SIZE-1
        ar(K)=ar(k)-WA(K)*RA(K)-WB(K)*RB(K)+0.5* &
              (WA(K)+WB(K))*(RA(K)+RB(K)-1.0) 
       !ar(K)=NXab*RA(K)*RB(K)-WA(K)*RA(K)-WB(K)*RB(K)+0.5* &
       !      (WA(K)+WB(K))*(RA(K)+RB(K)-1.0) 
      if(abs(NMu-0.0)<1.0e-4) then
         ar(K)=ar(K)
       else
           do j=0,N_dim_ddm-1
             do i=0,N_dim_ddm-1
               b(i,j)=M_OP(K,i,j)
              enddo
           enddo
       ar(K)=ar(K)+0.5*(double_dot_multi(b,b,N_dim_ddm)/NMu) 
      endif  
  enddo

# if defined (SLAB)       
# if defined (Dim2)   
   tempEE=simposon_2D_1D_mpi(SIDEy,local_nx,dy,dx,ar)
# else /* Dim2 */
   tempEE=simposon_3D_1D_mpi(SIDEz,SIDEy,local_nx,dz,dy,dx,ar)
# endif /* Dim2 */
# else /* SLAB */
   tempEE=simposon_3D_1D_mpi(xsize(1),xsize(2),xsize(3),dx,dy,dz,ar)
# endif /* SLAB */


   tempEE_global=0.0
              call mp_barrier()
              call mp_allreduce(tempEE,tempEE_global)
              call mp_barrier()
        tempEE_global=tempEE_global/M_v
        tempEE=tempEE_global

 FE_global=-pff_global+tempEE_global

 if(myid==0) then
 write(*,*) "FREE energy=",FE_global,pff_global,tempEE_global
 endif

     deallocate(ar)
     deallocate(b)
   call mp_barrier()    
end subroutine free_energy_finite    
     
   subroutine ar_from_finite(ar)
    USE nrtype,only :DP
    USE global_para
    USE control
    USE constants
    USE utility
    USE mmpi
    USE mpi_fftw3_operation
    USE matrix_inverse
    USE G_matrix_mod
    implicit none
    REAL(DP) :: temp1,pressure_coeff,M_grid_inv,M_sq
    integer :: k,error1
    real(dp),intent(inout) :: ar(0:Local_size-1)  
    Complex(dp),DIMENSION(:),ALLOCATABLE ::  yita,phi_k_A,phi_k_B
    Complex(dp),DIMENSION(:),ALLOCATABLE ::  w_k_A,w_k_B
    real(dp) :: error    

     M_grid_inv=1.0d0/M_grid
     M_sq=dsqrt(1.0d0/M_grid)

     allocate(phi_k_B(0:LOCAL_SIZE-1),stat=error1)
     allocate(w_k_A(0:LOCAL_SIZE-1),stat=error1)
     if(error1/=0) then
     write(mystd,*) "allocate dw failed! stop"
     stop
     endif
     
        do k=0,LOCAL_SIZE-1
# if defined (Dim2)
          local_data(kD(k)%x,kD(k)%y)=cmplx(RB(k),0.0)
# else  /* Dim2 */
          local_data(kD(k)%x,kD(k)%y,KD(k)%z)=cmplx(RB(k),0.0)
# endif /* Dim2 */
         enddo
       call fftw_mpi_execute_dft(fplan, local_data, local_data)
        do k=0,LOCAL_SIZE-1
# if defined (Dim2)
          phi_k_B(k)=local_data(kD(k)%x,kD(k)%y)*M_sq
# else  /* Dim2 */
          phi_k_B(k)=local_data(kD(k)%x,kD(k)%y,KD(k)%z)*M_sq
# endif /* Dim2 */
         enddo

        do k=0,LOCAL_SIZE-1
         temp1=dexp(-0.50d0*(sigma**2)*ksquare(k))
         w_k_A(k)=NXab*temp1*phi_k_B(k)            
        enddo

        do k=0,LOCAL_SIZE-1
# if defined (Dim2)
          local_data(kD(k)%x,kD(k)%y)=w_k_A(k)
# else  /* Dim2 */
          local_data(kD(k)%x,kD(k)%y,KD(k)%z)=w_k_A(k)
# endif /* Dim2 */
         enddo

       call fftw_mpi_execute_dft(bplan, local_data, local_data)

        do k=0,LOCAL_SIZE-1
# if defined (Dim2)
          ar(k)=real(local_data(kD(k)%x,kD(k)%y))*M_sq
# else  /* Dim2 */
          ar(k)=real(local_data(kD(k)%x,kD(k)%y,KD(k)%z))*M_sq
# endif /* Dim2 */
        enddo

        do k=0,LOCAL_SIZE-1
          ar(k)=ar(k)*RA(k)
        enddo

   deallocate(phi_k_B)
   deallocate(w_k_A)


   end subroutine ar_from_finite

# endif /* SPH */    
    
     




