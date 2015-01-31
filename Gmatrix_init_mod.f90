Module G_matrix_init
implicit none

contains

!generate the G matrix and its inverse

! if G_index >0, then G is for forward propagator q, if G_index<0, backward propagator q_star
! abs(G_index)=1,2,3 : BDF3 algorithm. 
 subroutine G_matrix_generator()
 !subroutine G_matrix_generator(G_R,G_I,ds,Loa,kapa_temp,M_bar,k_x,k_y,k_z,G_index,N_low)
 USE nrtype,only :DP,PI
 USE global_para,only : matrix_Rx,matrix_Ry,matrix_Rz,ds,Loa,basis_SPH
 USE mpi
 USE mmpi
 USE control
 USE constants
 USE utility
 USE matrix_inverse
 USE sparse_MKL
 implicit none
 integer :: i,j,k,M_bar
 integer :: G_index
 REAL(DP) :: kapa_temp,k_x,k_y,k_z
 !REAL(DP),dimension(:,:),intent(out) ::    G_R(:,:)
 !REAL(DP),dimension(:,:),intent(out) ::    G_I(:,:)
 !complex(DP),dimension(:,:),intent(out) :: G_C(:,:)
 REAL(DP),dimension(:,:),allocatable ::    G_R,G_R0,G_R1
 REAL(DP),dimension(:,:),allocatable ::    G_I
 complex(DP),dimension(:,:),allocatable :: G_C,G_C0,G_C1
real(DP) :: coeff_deltaij
real(DP) :: coeff_Rij
integer,parameter :: namaxA=3600
integer :: istat2
integer :: N_low(3)
integer :: loc(1)
integer job(8)
integer info
 DOUBLE PRECISION,allocatable :: sy_values_R(:)
 DOUBLE PRECISION,allocatable :: sy_values_I(:)
 Complex(DP),allocatable :: sy_values_C(:)
 integer,allocatable :: sy_columns_R(:),sy_rowIndex_R(:)
 integer,allocatable :: sy_columns_I(:),sy_rowIndex_I(:)
 integer,allocatable :: sy_columns_C(:),sy_rowIndex_C(:)

 DOUBLE PRECISION,dimension(:), allocatable :: Pk_R,Pk_I,qRR,qRI,qIR,qII
 Complex(DP),dimension(:), allocatable :: Pk_C,qC,qC1,qC2
 double precision :: t1,t2,t3,t4,t5,t10,t20,t30,t40,t50,sum_err 
 Complex(DP) :: aa
 
! we do not need to have set up the G matrix for A and B both,just one will do.

!note that the array is storaged in cloumn before row in fortran
!so ,it is different from the C/C++ version code.
!always remember,otherwise it would be terribly wrong!!!!!!!!!!!!!!!!!!


!k_x_factor=2.0*PI/(SIDEx*dx)
!k_y_factor=2.0*PI/(SIDEy*dy)
!k_z_factor=2.0*PI/(SIDEz*dz)
    if(myid==0) then
     write(*,*) "test begin"
    endif
     call mp_barrier() 
          M_Bar=81
          G_index=1
          k_x=0.0d0
          k_y=0.0d0
          k_z=-0.0d0
          kapa_temp=0.5
         
          allocate(sy_values_R(1:namaxA))
          allocate(sy_columns_R(1:namaxA))
          allocate(sy_rowIndex_R(1:M_Bar+1))

          allocate(sy_values_I(1:namaxA))
          allocate(sy_columns_I(1:namaxA))
          allocate(sy_rowIndex_I(1:M_Bar+1))

          allocate(sy_values_C(1:namaxA))
          allocate(sy_columns_C(1:namaxA))
          allocate(sy_rowIndex_C(1:M_Bar+1))

          allocate(Pk_R(1:M_bar))
          allocate(Pk_I(1:M_bar))
          allocate(Pk_C(1:M_bar))
          allocate(qRR(1:M_bar))
          allocate(qRI(1:M_bar))
          allocate(qIR(1:M_bar))
          allocate(qII(1:M_bar))
          allocate(qC(1:M_bar))
          allocate(qC1(1:M_bar))
          allocate(qC2(1:M_bar))

           allocate(G_R(1:M_bar,1:M_bar))
          allocate(G_R0(1:M_bar,1:M_bar))
          allocate(G_C0(1:M_bar,1:M_bar))
          allocate(G_C1(1:M_bar,1:M_bar))
          allocate(G_R1(1:M_bar,1:M_bar))
           allocate(G_I(1:M_bar,1:M_bar))
           allocate(G_C(1:M_bar,1:M_bar))


    if(myid==0) then
     write(*,*) "allocated done"
    endif
     call mp_barrier() 
do i=1,M_bar
         Pk_R(i)=ran2(seed)
         !Pk_R(i)=1.0d0
         Pk_I(i)=ran2(seed)
         Pk_C(i)=dcmplx(Pk_R(i),Pk_I(i))
enddo

          qRR=0.0d0
          qRI=0.0d0
          qIR=0.0d0
          qII=0.0d0 

          qC=dcmplx(0.0d0,0.0d0)


if(G_index>0) then
coeff_Rij=1.0
else
coeff_Rij=-1.0
endif

if(abs(G_index)==1) then
coeff_deltaij=1.0
else if (abs(G_index==2)) then
coeff_deltaij=1.5
else 
coeff_deltaij=11.0/6.0
endif


              !enter the i,j of spehrical harmonic basis index loop
              G_R=0.0d0
              G_I=0.0d0
              do j=1,M_Bar
              !G_R has only the digonal element
       G_R(j,j)=coeff_deltaij + ds*Loa*basis_SPH(j-1)%x*(basis_SPH(j-1)%x+1)/(2.0*kapa_temp)
                do i=1,M_Bar
       G_I(i,j)=coeff_Rij*ds*(matrix_Rx(i-1,j-1)*(k_x) + &
                       matrix_Ry(i-1,j-1)*(k_y) + &
                       matrix_Rz(i-1,j-1)*(k_z))
       G_C(i,j)=dcmplx(G_R(i,j),G_I(i,j))
                enddo
               enddo

    call mp_barrier()
         if(myid==0) then
           do j=1,M_Bar
            write(*,*) " G_R",G_R(j,j),"j=",j
          enddo
         endif
    call mp_barrier()
         if(myid==0) then
          do i=1,M_Bar
           do j=1,M_Bar
            write(*,*) " G_I",G_I(i,j),"i,",i,"j=",j
           enddo
          enddo
         endif
    call mp_barrier()
                
                
              G_R0=G_R
              G_C0=G_C

    
    if(myid==0) then
     write(*,*) " matrix init done"
    endif
     call mp_barrier() 
             ! before inverse, first calculate the nonzero element number of G_R and G_I
                  call get_symmetric_matrix_size_dp(G_R,M_Bar,M_Bar,N_low(1))
 !   if(myid==0) then
 !    write(*,*) "get size of real matrix done"
 !   endif
                  call get_symmetric_matrix_size_dp(G_I,M_Bar,M_Bar,N_low(2))
 !   if(myid==0) then
 !    write(*,*) "get size of imag matrix done"
 !   endif
     call mp_barrier() 
                  call get_symmetric_matrix_size_complex_dp(G_C,M_Bar,M_Bar,N_low(3))
 
            write(*,*) " bf inverse N_lowk1 is",N_low(1),"on",myid
            write(*,*) " bf inverse N_lowk2 is",N_low(2),"on",myid
            write(*,*) " bf inverse N_lowk3 is",N_low(3),"on",myid

              !inverse the complex matrix using the lapack lib directly
               
           t10=MPI_Wtime() 
        !      call complex_matrix_inverse(G_R,G_I,M_Bar)
        !   t20=MPI_Wtime() 
        !      call zmat_inv(M_bar, G_C) 
        !   t30=MPI_Wtime() 
        !          N_low=0
                  call get_symmetric_matrix_size_dp(G_R,M_Bar,M_Bar,N_low(1))
                  call get_symmetric_matrix_size_dp(G_I,M_Bar,M_Bar,N_low(2))
                  call get_symmetric_matrix_size_complex_dp(G_C,M_Bar,M_Bar,N_low(3))

    if(myid==0) then
     write(*,*) "inverse done"
    endif
     call mp_barrier() 
                 !N_low is precalculated so that namax is determined 
   
               do j=1,M_Bar
                  do i=1,M_Bar    
                  if(j>i) then
 !                 G_R(i,j)=0.0d0
                  G_I(i,j)=0.0d0
                  G_C(i,j)=dcmplx(0.0d0,0.0d0) 
                  endif
                   enddo
                enddo       
 !       call convert_symmetric_matrix_to_array_cmplx(G_z,max_M_Bar,max_M_Bar,N_low,sy_values_z,  &
 !       sy_columns_z,sy_rowIndex_z) !this routine has an issue when N_low doesn't match the actual nonzero
 !       number of the matrix.Needs to be fixed.

       write(*,*) " af  N_lowk1 is",N_low(1),"on",myid
       write(*,*) " af N_lowk2 is",N_low(2),"on",myid
       write(*,*) " af  N_lowk3 is",N_low(3),"on",myid

   
        job(1)=0
        job(2)=1
        job(3)=1
        job(4)=0
        !job(5)=namaxA
        job(5)=400
        job(6)=1


           t40=MPI_Wtime() 

    call mp_barrier()
         if(myid==0) then
          do i=1,M_Bar
           do j=1,M_Bar
            write(*,*) "inv G_R",G_R(i,j),"i,",i,"j=",j
           enddo
          enddo

         write(*,*) "inv G_R*G_R0=",sum(matmul(G_R,G_R0))  
         endif
          
    call mp_barrier()
         if(myid==0) then
          do i=1,M_Bar
           do j=1,M_Bar
            write(*,*) "inv G_I",G_I(i,j),"i,",i,"j=",j
           enddo
          enddo
         endif
          
    call mp_barrier()

       t1=MPI_Wtime()
     call convert_symmetric_matrix_to_array_dp(G_R,M_Bar,M_Bar,namaxA, &
         sy_values_R,sy_columns_R,sy_rowIndex_R)
       t2=MPI_Wtime()

               do j=1,M_Bar
                  do i=1,M_Bar 
                  if(j>i) then
                  G_R(i,j)=0.0d0
                  else
                  endif
                   enddo
                enddo       

      call mkl_ddnscsr(job,M_Bar,M_Bar,G_R,M_Bar,sy_values_R,sy_columns_R, &
        sy_rowIndex_R,info)
       t3=MPI_Wtime()

      call mkl_ddnscsr(job,M_Bar,M_Bar,G_I,M_Bar,sy_values_I,sy_columns_I, &
        sy_rowIndex_I,info)
           t50=MPI_Wtime() 
        call mp_barrier()
       call mkl_zdnscsr(job,M_Bar,M_Bar,G_C,M_Bar,sy_values_C,sy_columns_C, &
        sy_rowIndex_C,info)

         write(*,*) "compressing myine,t2-t1=",t2-t1 
         write(*,*) "compressing intel",t3-t2
 !      call convert_symmetric_matrix_to_array_cmplx(a,m,n,N_low, &
 !        sy_values,sy_col,sy_rowIndex)


         call mkl_dcsrsymv('l', M_Bar,  sy_values_R, sy_rowIndex_R, &
            sy_columns_R,Pk_R, qRR)
          ! t2=MPI_Wtime() 

           t1=MPI_Wtime() 
         call mkl_dcsrsymv('l', M_Bar,  sy_values_R, sy_rowIndex_R, &
            sy_columns_R,Pk_I, qRI)
           t2=MPI_Wtime() 

         call mkl_dcsrsymv('l', M_Bar,  sy_values_I, sy_rowIndex_I, &
            sy_columns_I,Pk_R, qIR)

         call mkl_dcsrsymv('l', M_Bar,  sy_values_I, sy_rowIndex_I, &
            sy_columns_I,Pk_I, qII)


           t3=MPI_Wtime() 
         call mkl_zcsrsymv('l', M_Bar,  sy_values_C, sy_rowIndex_C, &
            sy_columns_C,Pk_C, qC)


         call mp_barrier()  
 !       call symmetric_CSR_matvec(Pk_R,qRR,M_Bar,M_Bar,namaxA, & 
 !                      sy_values_R,sy_columns_R,sy_rowIndex_R)

           t1=MPI_Wtime() 
         call mkl_dcsrsymv('l', M_Bar,  sy_values_R, sy_rowIndex_R, &
            sy_columns_R,Pk_R, qRR)

         call mkl_dcsrsymv('l', M_Bar,  sy_values_R, sy_rowIndex_R, &
            sy_columns_R,Pk_I, qRI)

         call mkl_dcsrsymv('l', M_Bar,  sy_values_I, sy_rowIndex_I, &
            sy_columns_I,Pk_R, qIR)

         call mkl_dcsrsymv('l', M_Bar,  sy_values_I, sy_rowIndex_I, &
            sy_columns_I,Pk_I, qII)
           t2=MPI_Wtime() 


           t3=MPI_Wtime() 
         call mkl_zcsrsymv('l', M_Bar,  sy_values_C, sy_rowIndex_C, &
            sy_columns_C,Pk_C, qC)
           t4=MPI_Wtime() 

         write(*,*) "t2-t1=",t2-t1
         write(*,*) "t3-t2=",t4-t3 

        if(myid==0) then
         write(*,*) " qRR(1)",qRR(1)
         write(*,*) " qRR(2)",qRR(2)
         write(*,*) " qRR(3)",qRR(3)
         write(*,*) " qRR(4)",qRR(4)
         write(*,*) " qRR(9)",qRR(9)
         write(*,*) " qII(1)",qII(1)
         write(*,*) " qC(1)",qC(1)
         do i=1,M_Bar
         write(*,*) "Pk_R(i)",Pk_R(i),"i=",i
         enddo
        endif


!        job(1)=1
!        job(2)=1
!        job(3)=1
!        job(4)=2
!
!        job(5)=namaxA
!        job(6)=1
!
!      call mkl_zdnscsr(job,M_Bar,M_Bar,G_C1,M_Bar,sy_values_C,sy_columns_C, &
!        sy_rowIndex_C,info)
!      call mkl_ddnscsr(job,M_Bar,M_Bar,G_R,M_Bar,sy_values_R,sy_columns_R, &
!        sy_rowIndex_R,info)
!      call mkl_ddnscsr(job,M_Bar,M_Bar,G_I,M_Bar,sy_values_I,sy_columns_I, &
!        sy_rowIndex_I,info)
!
!             do j=1,M_Bar
!                do i=1,M_Bar    
!                if(j>i) then
!                G_C1(i,j)=G_C1(j,i)
!                G_C(i,j)=G_C(j,i)
!                G_R(i,j)=G_R(j,i)
!                G_I(i,j)=G_I(j,i)
!                  endif
!                   enddo
!                enddo       
!    call mp_barrier()
!         if(myid==0) then
!          do i=1,M_Bar
!           do j=1,M_Bar
!           ! write(*,*) "inv G_C",G_C1(i,j),"i,",i,"j=",j
!           ! write(*,*) "inv G_I",G_I(i,j),"i,",i,"j=",j
!           if(abs(real(G_C1(i,j))-G_R(i,j))>1.0E-10 .or. &
!           abs(aimag(G_C1(i,j))-G_I(i,j))>1.0E-10) then
!         write(*,*) "error ,G_C/=G_R+iG_I",G_R(i,j),G_I(i,j),G_C1(i,j)
!           endif
!            if(i==15) then
!            write(*,*) "inv G_C",G_C1(i,j),"i,",i,"j=",j
!            write(*,*) "inv G_R",G_R(i,j),"i,",i,"j=",j
!            write(*,*) "inv G_I",G_I(i,j),"i,",i,"j=",j
!             endif
!           enddo
!          enddo
!         endif
!    call mp_barrier()
           !qRR=matmul(G_R,Pk_R) 
           !G_I=matmul(G_C,G_C0)
           !Pk_C=cmplx(1.0d0,0.0d0)

  ! call dzgemv('N', M_Bar, M_Bar, 1.0d0, G_C, M_Bar, Pk_C, 1, 0, qC1, 1)
                !qC2=qC1
  ! call dzgemv('N', M_Bar, M_Bar, 1.0d0, G_C0, M_Bar, qC2, 1, 0, qC1, 1)
           !qC1=matmul(G_C1,Pk_C)
           !qRR=matmul(G_R,Pk_R)

                

!               qRR=0.0d0
!               qII=0.0d0
!               qC1=dcmplx(0.0d0,0.0d0)
!               do i=1,M_Bar
!                  do j=1,M_Bar
!                   pK_C(j)=dcmplx(Pk_R(j),Pk_I(j))
!                   aa=dcmplx(G_R(i,j),G_I(i,j))
!                   if(abs(aa-G_C1(i,j))>1.0E-10) then
!                       write(*,*) "error of G_C1/=G_R+iG_I"
!                      write(*,*) "G_C1(i,j)",G_C1(i,j),i,j
!                      write(*,*) "G_R/I(i,j)",G_R(i,j),G_I(i,j)
!                   endif
!                   qRR(i)=qRR(i)+G_R(i,j)*pk_R(j)
!                   qII(i)=qII(i)+G_I(i,j)*Pk_I(j)
!                   qC1(i)=qC1(i)+G_C1(i,j)*Pk_C(j)
!                   enddo
!                enddo 
!          
!         call mp_barrier()
!        if(myid==0) then
!         write(*,*) " 2qRR(1)",qRR(1)
!         write(*,*) " 2qRR(2)",qRR(2)
!         write(*,*) " 2qII(1)",qII(1)
!         write(*,*) " 2qC(1)",qC1(1)
!         write(*,*) "Pk_R(1)",Pk_R(1)
!         write(*,*) "G_R(1,1,),G_R(1,2),G_R(1,3)",G_R(1,1),G_R(1,2),G_R(1,3)  
!         write(*,*) "G_R(1,4,),G_R(1,4),G_R(1,5)",G_R(1,4),G_R(1,5),G_R(1,6)  
!         write(*,*) "G_R(1,1,),G_R(2,1),G_R(3,1)",G_R(1,1),G_R(2,1),G_R(3,1)  
!         write(*,*) "G_R(4,1,),G_R(5,1),G_R(6,1)",G_R(4,1),G_R(5,1),G_R(6,1)  
!        endif
!           call solver(G_R0,Pk_I,M_Bar)  
!         call mp_barrier()
           
       !  write(*,*) "sum GR*GR_inv",sum(G_I(:,:))     

  !        write(*,*) "t21=,",t2-t1,"t32=",t3-t2   

         write(*,*) "sum qC_R",sum(real(qC(:)))     
         write(*,*) "sum qRR-qII",sum(qRR(:)-qII(:))     
         write(*,*) "sum qC_I",sum(aimag(qC(:)))     
         write(*,*) "sum qRI+qIR",sum(qRI(:)+qIR(:))     

         write(*,*) "sum qRR",sum(qRR(:))     
 
           t4=MPI_Wtime() 
    deallocate(sy_values_R)
    deallocate(sy_columns_R)
    deallocate(sy_rowIndex_R)
    deallocate(sy_values_I)
    deallocate(sy_columns_I)
    deallocate(sy_rowIndex_I)
    deallocate(sy_values_C)
    deallocate(sy_columns_C)
    deallocate(sy_rowIndex_C)

   deallocate(Pk_R)
   deallocate(Pk_I)
   deallocate(Pk_C)
   deallocate(qRR)
   deallocate(qRI)
   deallocate(qIR)
   deallocate(qII)
   deallocate(qC)
   deallocate(qC1)
   deallocate(qC2)

   deallocate(G_R)
   deallocate(G_I)
   deallocate(G_C)
   deallocate(G_R0)
   deallocate(G_R1)
   deallocate(G_C0)
   deallocate(G_C1)

           t5=MPI_Wtime() 
 !         write(*,*) "t54=,",t5-t4 


     call mp_barrier() 

    end subroutine G_matrix_generator


subroutine dns2csr()
 
  complex(8),      allocatable :: Aall(:,:)
  integer,         allocatable :: ia  (:)
  integer,         allocatable :: ja  (:)
  complex(8),      allocatable :: zval(:)
  real(8),         allocatable :: dval(:)
  ! Variables needed by MKL
  integer                      :: job(1:8)
  integer                      :: info
  integer,           parameter :: nrow = 3
  integer,           parameter :: ncol = 3
   
  ! Allocate workspace
  allocate(Aall(1:3,1:3))
  allocate(ia (1:nrow+1))
  ! Initialize dns matrix
  Aall      = cmplx(0.0d0,0.0d0)
  Aall(1,1) = cmplx(1.0d0,0.0d0)
  Aall(3,1) = cmplx(5.0d0,0.0d0)
  Aall(2,2) = cmplx(4.0d0,0.0d0)
  Aall(3,2) = cmplx(0.0d0,4.0d0)
  Aall(1,3) = cmplx(3.0d0,0.0d0)
 
  ! Set job
  job(1) = 0 ! Convert from dense to CSR format
  job(2) = 1 ! One-based indexing for dense
  job(3) = 1 ! One-based indexing for CSR
  job(4) = 2 ! Pass the whole dense matrix to MKL routines 
  ! Compute ia only
  job(5) = 9 ! Max number of nnz
  job(6) = 1 ! Construct ia only
  allocate(ja(9),zval(9),dval(9))
  call mkl_zdnscsr(job,nrow,ncol,Aall,nrow,zval,ja,ia,info)
  print *,zval
 ! print *,ja
 ! print *,ia
 ! print *,info
  call mkl_ddnscsr(job,nrow,ncol,real(Aall),nrow,dval,ja,ia,info)
 ! print *,ia
 ! print *,info
 
  deallocate(ia,ja,dval,zval)
  deallocate(Aall)
 
end SUBROUTINE dns2csr

end module G_matrix_init
!*******************************************************************************
!   Copyright(C) 2005-2013 Intel Corporation. All Rights Reserved.
!   
!   The source code, information  and  material ("Material") contained herein is
!   owned  by Intel Corporation or its suppliers or licensors, and title to such
!   Material remains  with Intel Corporation  or its suppliers or licensors. The
!   Material  contains proprietary information  of  Intel or  its  suppliers and
!   licensors. The  Material is protected by worldwide copyright laws and treaty
!   provisions. No  part  of  the  Material  may  be  used,  copied, reproduced,
!   modified, published, uploaded, posted, transmitted, distributed or disclosed
!   in any way  without Intel's  prior  express written  permission. No  license
!   under  any patent, copyright  or  other intellectual property rights  in the
!   Material  is  granted  to  or  conferred  upon  you,  either  expressly,  by
!   implication, inducement,  estoppel or  otherwise.  Any  license  under  such
!   intellectual  property  rights must  be express  and  approved  by  Intel in
!   writing.
!   
!   *Third Party trademarks are the property of their respective owners.
!   
!   Unless otherwise  agreed  by Intel  in writing, you may not remove  or alter
!   this  notice or  any other notice embedded  in Materials by Intel or Intel's
!   suppliers or licensors in any way.
!
!*******************************************************************************
!   Content : MKL Sparse BLAS Fortran-77 example
!
!*******************************************************************************
!
! Example program for using MKL Sparse BLAS Level 2 and 3
! for matrices represented in the compressed sparse row storage scheme.
! The following Sparse  Blas routines are used in the example:
!          MKL_DCSRSM  MKL_DCSRSV  MKL_DCSRMM  MKL_DCSRMV
!          MKL_DCSRGEMV    MKL_DCSRSYMV  MKL_DCSRTRSV.
!
! Consider the matrix A (see Appendix 'Sparse Storage Formats for Sparse Blas
! level 2-3')
!
!                 |   1       -1      0   -3     0   |
!                 |  -2        5      0    0     0   |
!   A    =        |   0        0      4    6     4   |,
!                 |  -4        0      2    7     0   |
!                 |   0        8      0    0    -5   |
!
!
! decomposed as
!
!                      A = L + D + U,
!
!  where L is the strict  lower triangle of A, U is the strictly  upper triangle
!  of A, D is the main diagonal. Namely
!
!        |   0    0   0    0     0   |       |  0   -1    0   -3   0   |
!        |  -2    0   0    0     0   |       |  0    0    0    0   0   |
!   L  = |   0    0   0    0     0   |,  U=  |  0    0    0    6   4   |
!        |  -4    0   2    0     0   |       |  0    0    0    0   0   |
!        |   0    8   0    0     0   |       |  0    0    0    0   0   |
!
!
!           |   1  0  0   0   0   |
!           |   0  5  0   0   0   |
!   D    =  |   0  0  4   0   0   |.
!           |   0  0  0   7   0   |
!           |   0  0  0   0  -5   |
!
!  The matrix A is represented in the compressed sparse row storage scheme with the help of three
!  arrays  (see Appendix 'Sparse Matrix Storage') as follows:
!
!         values = (1 -1 -3 -2 5 4 6 4 -4 2 7 8 -5)
!         columns = (1 2 4 1 2 3 4 5 1 3 4 2 5)
!         rowIndex = (1  4  6  9  12 14)
!
!  It should be noted that two variations of the compressed sparse row storage scheme are supported by
!  Intel MKL Sparse Blas (see 'Sparse Storage Formats for Sparse Blas level 2-3') :
!
!        1. variation accepted in the NIST Sparse Blas
!        2. variation accepted for DSS/PARDISO, CXML and many other libraries.
!
!  The representation of the matrix A  given above is the PARDISO's variation. Two integer arrays
!  pointerB and pointerE instead of the array rowIndex are used in the NIST variation of variation
!  of the compressed sparse row format. Thus the arrays values and columns are the same for the both
!  variations. The arrays pointerB and pointerE for the matrix A are defined as follows:
!                          pointerB = (1 4  6  9 12)
!                          pointerE = (4 6  9 12 14)
!  It's easy to see that
!                    pointerB(i)= rowIndex(i) for i=1, ..5;
!                    pointerE(i)= rowIndex(i+1) for i=1, ..5.
!
!
!  The purpose of the given example is to show
!
!             1. how to call routines having interfaces suitable for the NIST's variation of the
!                compressed sparse row format
!             2. how to form the arrays pointerB and pointerE for the NIST's variation of the
!                compressed sparse row format using the  array rowIndex
!             3. how to use minors of the matrix A by redefining the arrays pointerB and pointerE
!                but the arrays values and columns are the same.
!
!  In what follows the symbol ' means taking of transposed.
!
!  The test performs the following operations :
!
!       1. The code computes (L+D)'*S = F using MKL_DCSRMM where S is a known 5 by 2
!          matrix and then the code solves the system (L+D)'*X = F with the help of MKL_DCSRSM.
!          It's evident that X should be equal to S.
!
!       2. The code computes (U+I)'*S = F using MKL_DCSRMV where S is a vector
!          and then the code calls MKL_DCSRTRSV solves the system (U+I)'*X = F with the single right
!          hand side. It's evident that X should be equal to S.
!
!       3. The code computes D*S = F using MKL_DCSRMV where S is a vector
!          and then the code solves the system D*X = F with the single right hand side.
!          It's evident that X should be equal to S.
!
!       4. The next step is the computation (U-U') S = F using MKL_DCSRMV where S is
!          a vector. It is easy to see that U-U' is a skew-symmetric matrix.
!
!       5. The next step is the computation (L+D+L') S = F using MKL_DCSRSYMV where S is
!          a vector. The vector is computed two times. At first, the sparse representation
!          of the whole matrix A is used. Then the vector is computed with the help of
!          sparse representation of L+D. These two calls must give the same vector.
!
!       6. The next step is the computation A'* S = F using MKL_DCSRGEMV where S is
!          a vector.
!
!       7. Let's T be the upper 3 by 3 minor of the matrix A. Namely, T is the following matrix
!
!                        |   1       -1      0   |
!          T    =        |  -2        5      0   |.
!                        |   0        0      4   |
!          The test performs the matrix-vector multiply T*S=F with the same arrays values, columns
!          and pointerB used before for the whole matrix A. It is enough to change two values of
!	   array pointerE in order to use the minor under consideration. Then the test solves the system
!          T*X =F using MKL_DCSRSV. The routine MKL_DCSRMV is used for getting matrix-vector multiply.
!
! The code given below uses only one sparse representation for the all operations.
!
!*******************************************************************************
!    Declaration of arrays for sparse representation of  the matrix A
!    and the lower triangle of A in the compressed sparse row format:
!*******************************************************************************
          subroutine testMKL() 
           USE sparse_MKL
          implicit none
          integer  m,  nnz, mnew, nnz1
          parameter( m = 5,  nnz=13, mnew=3, nnz1=9)
          real*8  values(nnz), values1(nnz1)
          integer columns(nnz), rowIndex(m+1), columns1(nnz1), &
            rowIndex1(m+1)
          integer pointerB(m) , pointerE(m)
!*******************************************************************************
!    Sparse representation of the matrix A
!*******************************************************************************

          data values/1.d0, -1.d0, -3.d0, -2.d0, 5.d0, &
             4.d0, 6.d0, 4.d0, -4.d0, 2.d0, 7.d0, 8.d0, -5.d0/
          data rowIndex/1, 4,  6,  9,  12, 14/
          data columns/1, 2, 4, 1, 2, 3, 4, 5, 1, 3, 4, 2, 5/
!*******************************************************************************
!    Sparse representation of the lower triangle L+D
!*******************************************************************************
          data values1/1.d0, -2.d0, 5.d0, 4.d0,  -4.d0, 2.d0, &
            7.d0,  8.d0, -5.d0/
          data columns1/1,  1, 2, 3, 1, 3, 4, 2, 5/
          data rowIndex1/1, 2, 4, 5, 8, 10/
!*******************************************************************************
!    Declaration of local variables :
!*******************************************************************************
          integer n
          parameter (n=2)
          real*8 rhs(m, n), sol(m, n), temp(m, n)
          real*8 rhs1(m), sol1(m), temp1(m)
          data sol/1.D0, 1.D0, 1.D0, 1.D0, 1.D0, &
         5.D0, 4.D0, 3.D0, 2.D0, 1.D0/
          real*8 alpha, beta
          data alpha/1.d0/, beta/0.d0/
          integer i, j, is
          print*
          print*, ' EXAMPLE PROGRAM FOR COMPRESSED SPARSE ROW'

!*******************************************************************************
! Task 2.    Obtain matrix-vector multiply (U+I)' *sol --> rhs
!    and solve triangular system   (U+I)' *temp = rhs with single right hand sides.
!    Array temp must be equal to the array sol.
!
!    Let us form the arrays pointerB and pointerE for the NIST's variation of the
!    compressed sparse row format using the  array rowIndex.
!
!*******************************************************************************
          do i=1, m
            pointerB(i)=rowIndex(i)
            pointerE(i)=rowIndex(i+1)
          enddo
          print*
          print*, '     INPUT DATA FOR MKL_DCSRMV '
          print*, '     WITH TRIANGULAR MATRIX  '
          print 102, alpha, beta
          print 103, 't'
          print*, ' Input vector '
          print 105, (sol(i,1),i=1,m)

           call mkl_dcsrmv('t', m, m, alpha, 'tuu', &
               values, columns, pointerB, pointerE, sol, beta, rhs)
          print*
          print*, '     OUTPUT DATA FOR MKL_DCSRMV '
          print*, '     WITH TRIANGULAR MATRIX  '
          print 105, (rhs(i,1),i=1,m)
          print 100
          print*, ' Solve triangular system with obtained '
          print*, ' right hand side  '
          call mkl_dcsrtrsv('u', 't', 'u', m, &
                values, rowIndex, columns, rhs, temp)
          print*
          print*, '     OUTPUT DATA FOR MKL_DCSRTRSV '
          print*, '     WITH TRIANGULAR MATRIX  '
          print 105, (temp(i,1),i=1,m)
          print 100
!*******************************************************************************
! Task 5.    Obtain matrix-vector multiply (L+D+L')*sol --> rhs with the help of
!    MKL_DCSRSYMV
! NOTE: The routine mkl_dcsrsymv as well as the similar dense Level 2
! routine dsymv has a possibilty to extract the required triangle from the input
! sparse matrix and perform symmetric matrix-vector multiply with the help of
! the triangle. Let the arrays values, rowIndex, columns be the sparse
! representation of A
!
!                 |   1       -1      0   -3     0   |
!                 |  -2        5      0    0     0   |
!   A    =        |   0        0      4    6     4   |,
!                 |  -4        0      2    7     0   |
!                 |   0        8      0    0    -5   |
! Let the arrays values1, rowIndex1, columns1 which are the following
!   data values1/1.d0, -2.d0, 5.d0, 4.d0,  -4.d0, 2.d0, 7.d0,  8.d0, -5.d0/
!   data column1/1,  1, 2, 3, 1, 3, 4, 2, 5/
!   data rowIndex/1, 2, 4, 5, 8, 10/
!
! be the sparse representation of the lower triangle of A
!
!           |   1    0   0    0     0   |
!           |  -2    5   0    0     0   |
!   L +D  = |   0    0   4    0     0   |,
!           |  -4    0   2    7     0   |
!           |   0    8   0    0    -5   |
!
!  The feature described above means that  the following two calls must give the
!  same output vector
!  call mkl_dcsrsymv(uplo, m, values, rowIndex, columns, sol_vec, rhs_vec);
!  call mkl_dcsrsymv(uplo, m, values1, rowIndex1, columns1, sol_vec, temp);
!
!  The test checks whether these two calls give the same output vector.
!*******************************************************************************
          print*
          print*, '     SYMMETRIC MATRIX-VECTOR MULTIPLY ROUTINE '
          print*, '     INPUT DATA FOR MKL_DCSRSYMV '
          print*, ' Input vector '
          print 105, (sol(i, 1),i=1,m)

           call mkl_dcsrsymv('l', m,  values, rowIndex, columns, &
                sol, rhs)
          call mkl_dcsrsymv('l', m,  values1, rowIndex1, columns1, &
                sol, temp)

             sol1(:)=sol(:,1)
          call symmetric_CSR_matvec(sol(:,1),temp(:,1),5,5,9, &
                       values1,columns1,rowIndex1)
          print*
          print*, '     SYMMETRIC MATRIX-VECTOR MULTIPLY '
          print*, '     MKL_DCSRSYMV CALLED TWO TIMES '
          do i=1, m
              print 104, rhs(i,1), temp(i,1),temp1(i)
          enddo
          print 100

 100      format('------------------------------------------------')
 101      format(7x,'M=',i1,'  N=',i1)
 102      format(7x,'ALPHA= ',f4.1,' BETA= ', f4.1)
 103      format(7x,'TRANS=',a1)
 104      format(2(f7.1, 3x))
 105      format(f4.1)
          end subroutine testMKL
