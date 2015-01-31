# if defined (SPH) 
      subroutine density_sph()
       USE nrtype,only :DP,PI
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
       use spherepack
       use split_operator_MDE
       
       implicit none
       integer :: i,j,s,s1,K
       real(DP) :: pff_temp_global,pff_temp,totden_global,dcs
       real(DP) :: RAtot,RAtot_global
       REAL(DP),allocatable :: ar(:)
       REAL(DP),allocatable :: br(:)
       REAL(DP),allocatable :: c_A(:)
       REAL(DP),allocatable :: c_B(:)
       REAL(DP),allocatable :: denz(:)
       REAL(DP) :: sum_iA,sum_iB,sum_half,sum_end_B
       REAL(DP) :: nor_coeff
       REAL(DP) :: sumbr,sumbr_global
      ! REAL(DP),allocatable :: c11_A(:)
      ! REAL(DP),allocatable :: c12_A(:)
      ! REAL(DP),allocatable :: c22_A(:)
      ! REAL(DP),allocatable :: c13_A(:)
      ! REAL(DP),allocatable :: c23_A(:)

      ! REAL(DP),allocatable :: c11_B(:)
      ! REAL(DP),allocatable :: c12_B(:)
      ! REAL(DP),allocatable :: c22_B(:)
      ! REAL(DP),allocatable :: c13_B(:)
      ! REAL(DP),allocatable :: c23_B(:)
      real(DP) :: ar_u(1:Ntheta+1,1:Nphi)
      real(DP) :: den_nor(0:Nmax)
      !~ warning,do not statically define large array,stack overflow may occur
 !    real(DP) :: RHOB_SPH(1:Ntheta+1,1:Nphi,0:Local_size-1)
 !    real(DP) :: RHOB_SPH(1:Ntheta+1,1:Nphi,0:local_size-1)
      
      real(DP),allocatable :: RHOA_SPH(:,:,:)
      real(DP),allocatable :: RHOB_SPH(:,:,:)

      allocate(RHOA_SPH(1:Ntheta+1,1:Nphi,0:Local_size-1))
      allocate(RHOB_SPH(1:Ntheta+1,1:Nphi,0:Local_size-1))

!# if defined (Debug)
    call mp_barrier()
    if(myid==0 ) then
        write(*,*) "enter density"
    endif
    call mp_barrier()
!# endif /* Debug */

!        do s=0,Nmax*100
!           do s1=0,Nmax*100
!             do i=0,Nmax*100
!               do j=0,Nmax*100
!                  dcs=exp(cos(j*1.0d0)**2)
!               enddo
!              enddo
!            enddo
!         enddo    
!
!        write(*,*) "done cal density"
!    call mp_barrier()



!   do s=0,Nmax
!        dcs=1.0d0
!      do s1=s+1,Nmax
!        dcs=dcs*c_nor_star(s1)/c_nor(s1)
!      enddo
!        den_nor(s)=dcs*c_nor_star(s)
!   enddo

       allocate(ar(0:LOCAL_SIZE-1))
       do K=0,LOCAL_SIZE-1
          do j=1,Nphi
           do i=1,Ntheta+1
         ar_u(i,j)=q_f(k,i,j,Nmax)   
           enddo
          enddo
         ar(K)=simposon_sph_integral(Ntheta,Nphi,dtheta,dphi,ar_u)
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

        pff_global=pff_temp_global/(M_v*4*PI)

!# if defined (Debug) 
        if(myid==0) then
         write(*,*) "arsum",RAtot_global,"pff before log",pff_global
         write(*,*) "pff_global ==",dlog(pff_global)
         endif    
              call mp_barrier()
!# endif /* Debug */
         
        pff_global=dlog(pff_global)

!         do s=0,Nmax
!           pff_global=pff_global+dlog(c_nor(s))
!         enddo 

        if(myid==0) then
         write(*,*) "af norm pff_global ==",pff_global
         write(*,*) "sum of c_nor(s)",sum(c_nor)
         endif    

        deallocate(ar)
        allocate(c_A(0:Nmax))          
        !allocate(c_B(0:NB)) 
        allocate(denz(0:LOCAL_SIZE-1))
         call mp_barrier()

       do K=0,LOCAL_SIZE-1
         do j=1,Nphi
            do i=1,Ntheta+1
              do s=0,Nmax
 !             c_A(s)=q_f(K,i,j,s)*q_b(K,i,j,s)*den_nor(s) 
               c_A(s)=q_f(K,i,j,s)*q_b(K,i,j,s)
              enddo
          RHOA_SPH(i,j,K)=simposon_1D_NR(0,NA,ds,c_A(0:NA))
          RHOB_SPH(i,j,K)=simposon_1D_NR(NA,Nmax,ds,c_A(NA:Nmax))
            enddo
         enddo
           ar_u(:,:)=RHOA_SPH(:,:,K)
           RA(K)=simposon_sph_integral(Ntheta,Nphi,dtheta,dphi,ar_u)
           ar_u(:,:)=RHOB_SPH(:,:,K)
           RB(K)=simposon_sph_integral(Ntheta,Nphi,dtheta,dphi,ar_u)
           denz(K)=RA(K)+RB(K)
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
          RA(K)=M_v*RA(K)/totden
          RB(K)=M_v*RB(K)/totden
      enddo


!!!! for S 
!        nor_coeff=M_v/totden
!        allocate(c11_A(0:NA))
!        allocate(c12_A(0:NA))
!        allocate(c13_A(0:NA))
!        allocate(c23_A(0:NA))
!        allocate(c22_A(0:NA))
!
!        allocate(c11_B(0:NB))
!        allocate(c12_B(0:NB))
!        allocate(c13_B(0:NB))
!        allocate(c23_B(0:NB))
!        allocate(c22_B(0:NB))
!
!    do k=0,LOCAL_SIZE-1
!       do s=0,NA
!# if defined (EXP)
!  C11_A(s)=Sum_sparse_2D(M_Bar_A,J11ij_nonzero_1D_A,J11_A,qA(:,K,s),qAstar(:,K,s))
!  C22_A(s)=Sum_sparse_2D(M_Bar_A,J22ij_nonzero_1D_A,J22_A,qA(:,K,s),qAstar(:,K,s))
!  C12_A(s)=Sum_sparse_2D(M_Bar_A,J12ij_nonzero_1D_A,J12_A,qA(:,K,s),qAstar(:,K,s))
!  C23_A(s)=Sum_sparse_2D(M_Bar_A,J23ij_nonzero_1D_A,J23_A,qA(:,K,s),qAstar(:,K,s))
!  C13_A(s)=Sum_sparse_2D(M_Bar_A,J13ij_nonzero_1D_A,J13_A,qA(:,K,s),qAstar(:,K,s))
!# else /* EXP */
!  C11_A(s)=Sum_sparse_2D(M_Bar_A,J11ij_nonzero_1D_A,J11_A,qA(s,:,K),qAstar(s,:,K))
!  C22_A(s)=Sum_sparse_2D(M_Bar_A,J22ij_nonzero_1D_A,J22_A,qA(s,:,K),qAstar(s,:,K))
!  C12_A(s)=Sum_sparse_2D(M_Bar_A,J12ij_nonzero_1D_A,J12_A,qA(s,:,K),qAstar(s,:,K))
!  C23_A(s)=Sum_sparse_2D(M_Bar_A,J23ij_nonzero_1D_A,J23_A,qA(s,:,K),qAstar(s,:,K))
!  C13_A(s)=Sum_sparse_2D(M_Bar_A,J13ij_nonzero_1D_A,J13_A,qA(s,:,K),qAstar(s,:,K))
!# endif /* EXP */
!      enddo
!
!       do s=0,NB
!# if defined (EXP)
!  C11_B(s)=Sum_sparse_2D(M_Bar_B,J11ij_nonzero_1D_B,J11_B,qB(:,K,s),qBstar(:,K,s))
!  C12_B(s)=Sum_sparse_2D(M_Bar_B,J12ij_nonzero_1D_B,J12_B,qB(:,K,s),qBstar(:,K,s))
!  C13_B(s)=Sum_sparse_2D(M_Bar_B,J13ij_nonzero_1D_B,J13_B,qB(:,K,s),qBstar(:,K,s))
!  C22_B(s)=Sum_sparse_2D(M_Bar_B,J22ij_nonzero_1D_B,J22_B,qB(:,K,s),qBstar(:,K,s))
!  C23_B(s)=Sum_sparse_2D(M_Bar_B,J23ij_nonzero_1D_B,J23_B,qB(:,K,s),qBstar(:,K,s))
!# else /* EXP */
!  C11_B(s)=Sum_sparse_2D(M_Bar_B,J11ij_nonzero_1D_B,J11_B,qB(s,:,K),qBstar(s,:,K))
!  C12_B(s)=Sum_sparse_2D(M_Bar_B,J12ij_nonzero_1D_B,J12_B,qB(s,:,K),qBstar(s,:,K))
!  C13_B(s)=Sum_sparse_2D(M_Bar_B,J13ij_nonzero_1D_B,J13_B,qB(s,:,K),qBstar(s,:,K))
!  C22_B(s)=Sum_sparse_2D(M_Bar_B,J22ij_nonzero_1D_B,J22_B,qB(s,:,K),qBstar(s,:,K))
!  C23_B(s)=Sum_sparse_2D(M_Bar_B,J23ij_nonzero_1D_B,J23_B,qB(s,:,K),qBstar(s,:,K))
!# endif /* EXP */
!       enddo
!         
!      SA_OP(K,0,0)=simposon_1D_NR(0,NA,ds,c11_A)*nor_coeff
!      SA_OP(K,1,1)=simposon_1D_NR(0,NA,ds,c22_A)*nor_coeff
!      SA_OP(K,2,2)=-(SA_OP(K,0,0)+SA_OP(K,1,1))
!      SA_OP(K,0,1)=simposon_1D_NR(0,NA,ds,c12_A)*nor_coeff
!      SA_OP(K,0,2)=simposon_1D_NR(0,NA,ds,c13_A)*nor_coeff
!      SA_OP(K,1,2)=simposon_1D_NR(0,NA,ds,c23_A)*nor_coeff
!      SA_OP(K,1,0)=SA_OP(K,0,1)
!      SA_OP(K,2,0)=SA_OP(K,0,2)
!      SA_OP(K,2,1)=SA_OP(K,1,2)
!
!      SB_OP(K,0,0)=simposon_1D_NR(0,NB,ds,c11_B)*nor_coeff
!      SB_OP(K,1,1)=simposon_1D_NR(0,NB,ds,c22_B)*nor_coeff
!      SB_OP(K,2,2)=-(SB_OP(K,0,0)+SB_OP(K,1,1))
!      SB_OP(K,0,1)=simposon_1D_NR(0,NB,ds,c12_B)*nor_coeff
!      SB_OP(K,0,2)=simposon_1D_NR(0,NB,ds,c13_B)*nor_coeff
!      SB_OP(K,1,2)=simposon_1D_NR(0,NB,ds,c23_B)*nor_coeff
!      SB_OP(K,1,0)=SB_OP(K,0,1)
!      SB_OP(K,2,0)=SB_OP(K,0,2)
!      SB_OP(K,2,1)=SB_OP(K,1,2)
!    
!       S_OP(K,0,0)=SA_OP(K,0,0)+SB_OP(K,0,0)
!       S_OP(K,1,1)=SA_OP(K,1,1)+SB_OP(K,1,1)
!       S_OP(K,2,2)=SA_OP(K,2,2)+SB_OP(K,2,2)
!       S_OP(K,0,1)=SA_OP(K,0,1)+SB_OP(K,0,1)
!       S_OP(K,0,2)=SA_OP(K,0,2)+SB_OP(K,0,2)
!       S_OP(K,1,2)=SA_OP(K,1,2)+SB_OP(K,1,2)
!       S_OP(K,1,0)=S_OP(K,0,1)  
!       S_OP(K,2,0)=S_OP(K,0,2)  
!       S_OP(K,2,1)=S_OP(K,1,2)  
!
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
   !enddo


  deallocate(c_A)
!  deallocate(c_B)
  deallocate(denz)
  if (allocated(RHOA_SPH)) deallocate(RHOA_SPH)
  if (allocated(RHOB_SPH)) deallocate(RHOB_SPH)
 
!  deallocate(c11_A)
!  deallocate(c22_A)
!  deallocate(c12_A)
!  deallocate(c13_A)
!  deallocate(c23_A)
!  deallocate(c11_B)
!  deallocate(c22_B)
!  deallocate(c12_B)
!  deallocate(c13_B)
!  deallocate(c23_B)

# if defined (Debug)
   call mp_barrier()
   if(myid==0) then
   write(*,*) "RA(0),RA(30),RA(60)=",RA(0),RA(30),RA(60)
   write(*,*) "RB(0),RB(30),RB(60)=",RB(0),RB(30),RB(60)
   endif
# endif /* Debug */

   call mp_barrier()    
   end subroutine density_sph

   subroutine free_energy_sph()
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


!      if(abs(NMu-0.0)<1.0e-5) then
!         ar(K)=ar(K)
!       else
!           do j=0,N_dim_ddm-1
!             do i=0,N_dim_ddm-1
!               b(i,j)=M_OP(K,i,j)
!              enddo
!           enddo
!       ar(K)=ar(K)+0.5*(double_dot_multi(b,b,N_dim_ddm)/NMu) 
!
!
!      endif  
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
end subroutine free_energy_sph    
     
# endif /* SPH */    
    
     


