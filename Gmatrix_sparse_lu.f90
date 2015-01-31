   module  Gmatrix_lu_mod       
     use m_sparse
# if defined (SLAB) 
     use G_matrix_mod
# else /* SLAB */
     use G_matrix_2dcomp_mod
# endif /* SLAB */
     implicit none
        private 
        TYPE(pardiso_plan),public :: spa_pA1
        TYPE(pardiso_plan),public :: spa_pA2
        TYPE(pardiso_plan),public :: spa_pB1
        TYPE(pardiso_plan),public :: spa_pB2
        public :: G_matrix_lu_init 
        public :: test_sparse
        !public :: G_matrix_lu_solve
        !public :: G_matrix_lu_clean

        contains

     subroutine G_matrix_lu_init()
       USE nrtype,only :DP,PI
       USE mmpi
       USE global_para
       USE control
       USE constants
       USE mpi_fftw3_operation
       USE matrix_inverse
       USE m_sparse
       implicit none
       integer :: i,j,k
       integer :: K_i,K_ii
       integer :: n_matrixs,nrow,ncol,n_max,nn
       complex(dp) :: pk(1:M_Bar_A,1)      
       complex(dp) :: q_tmp(1:M_Bar_A,1)      
       complex(dp) :: G_C(1:M_Bar_A,1:M_Bar_A)      
       real(dp) :: aa,t1,t2
       integer      job(8)
       !integer      lda, info
       CHARACTER	    UPLO
       INTEGER	    INFO, LDA, LDB, LWORK, NRHS
       INTEGER	    IPIV(M_Bar_A )
       COMPLEX*16   A( M_Bar_A, M_Bar_A ), B(M_Bar_A, 1 ), WORK(M_Bar_A*64 )

      NRHS=1
      LDA=M_Bar_A  
      LDB=M_Bar_A  
      LWORK=M_Bar_A*64
      !LWORK=-1

       nrow=M_Bar_A
       ncol=M_Bar_A
       n_matrixs=local_size*4 ! four G matrixs for A block
       n_matrixs=1
       n_matrixs=local_size
       n_matrixs=1
! chose one matrix to do symbolic factorization
       nn= 7
       n_max=size(GA1_inverse_val,1)
       call mkl_pardiso_init(n_matrixs,nrow,ncol,n_max,GA1_inverse_row(:,nn), &
                               GA1_inverse_col(:,nn),GA1_inverse_val(:,nn) ,spa_pA1)

       call mkl_pardiso_factrz(n_matrixs,nrow,ncol,n_max,GA1_inverse_row(:,nn), &
                               GA1_inverse_col(:,nn),GA1_inverse_val,spa_pA1)

        do i=1,M_Bar_A
           pk(i,1)=dcmplx(i*0.1,i*0.2)
           q_tmp(i,1)=dcmplx(i*0.12,i*0.25)
        enddo 
 
         !t1=MPI_Wtime()

        !call mkl_zcsrsymv('u', M_Bar_A, GA1_inverse_val(:,nn), GA1_inverse_row(:,nn), &
        !     GA1_inverse_col(:,nn),Pk(:,1), q_tmp(:,1))

         !t2=MPI_Wtime()
        !write(*,*) "q_tmp(1,1)",q_tmp(1,1),q_tmp(2,1),"t21=",t2-t1
       !call mkl_set_num_threads( 1 )

         t1=MPI_Wtime()
      call mkl_pardiso_solve(n_matrixs,M_Bar_A,M_Bar_A,nmax,GA1_inverse_row(:,nn), &
                         GA1_inverse_col(:,nn), &
                         GA1_inverse_val,q_tmp,pk,spa_pA1)


         t2=MPI_Wtime()
       ! do i=1,M_Bar_A
       !   aa=aa+abs( pk(i,1)-dcmplx(i*0.1,i*0.2))
       ! enddo 
       ! write(*,*) "aa=",aa,"t21 in pardiso",t2-t1
         write(*,*) "pk(1,1)",pk(1,1),pk(2,1)

       job(1)=1 
       job(2)=1
       job(3)=1
       job(4)=2
       job(5)=n_max
       job(6)=1
       lda=nrow


         t1=MPI_Wtime()
       ! call sparse_mkl_dnscsr(M_Bar_A, M_Bar_A, nmax, G_C, 1, &
       !      GA1_inverse_val(:,nn), GA1_inverse_col(:,nn), GA1_inverse_row(:,nn))

        call mkl_zdnscsr(job, M_Bar_A, M_Bar_A, G_C, lda,GA1_inverse_val(:,nn), &
                            GA1_inverse_col(:,nn), GA1_inverse_row(:,nn), info)
        t1=MPI_Wtime()
        call zmat_inv(M_Bar_A, G_C)
         t2=MPI_Wtime()
        write(*,*) "t21 in direct inve",t2-t1,"err",info
         t1=MPI_Wtime()
         pk(:,1)=matmul(G_C,q_tmp(:,1))
         t2=MPI_Wtime()
        write(*,*) "t21 in direc  matmul",t2-t1,"err",info

        write(*,*) "pk(1,1) using inv ",pk(1,1),pk(2,1)

         t1=MPI_Wtime()
         call solver_z(G_C,q_tmp(:,1),M_Bar_A)
         t2=MPI_Wtime()
        write(*,*) "t21 in lapack solve",t2-t1

       !  t1=MPI_Wtime()
       !  call ZSYSV( 'U', M_Bar_A, 1, G_C, M_Bar_A, IPIV, q_tmp(:,1), LDB, WORK, LWORK, INFO ) 
       !  t2=MPI_Wtime()
       !  write(*,*) "t21 in lapack ZSYSV solve",t2-t1,"lwork=",work(1)
         write(*,*) "sum of err",sum(abs(q_tmp(:,1)-pk(:,1)))
          write(*,*) "G_C(3,3),G_(3,4)",G_C(3,3),G_C(3,4)  
          call zgetrf(M_Bar_A,M_Bar_A,G_C,M_Bar_A,ipiv,info)
         write(*,*) "info =",info
          write(*,*) "ad LU G_C(3,3),G_(3,4)",G_C(3,3),G_C(3,4)  
          call get_sparse_size(M_Bar_A,M_Bar_A,n_max,G_C)         
         write(*,*) "num of namx=",n_max





   !    nrow=M_Bar_B
   !    ncol=M_Bar_B
   !    n_matrixs=local_size*4
   !    nn= local_size/2+1
   !    nmax=size(GB1_inverse_val,1)
   !    call mkl_pardiso_init(n_matrixs,nrow,ncol,nmax,GB3_inverse_row(:,nn), &
   !                            GB3_inverse_col(:,nn),GB3_inverse_val(:,nn) ,spa_pB1)
   !  
   !    call mkl_pardiso_factrz(n_matrixs,nrow,ncol,nmax,GB3_inverse_row(:,nn), &
   !                            GB3_inverse_col(:,nn),GB3_inverse_val,spa_pB1)

     end subroutine G_matrix_lu_init

!     subroutine G_matrix_lu_solve(G_index)
!       USE nrtype,only :DP,PI
!       USE mmpi
!       USE global_para
!       USE control
!       USE constants
!       USE mpi_fftw3_operation
!       implicit none
!       integer :: i,j,k
!       integer :: K_i,K_ii
!       integer,intent(in) :: G_index  
!
!
!
!
!
!
!
!
!
!
!
!   call mkl_pardiso_solve_z(n_matrixs,nrow,ncol,nmax,ia,ja,a_all,rhs_all,sol_all,spa_p)
!
!     end subroutine G_matrix_lu_solve()

     subroutine test_sparse()
       USE nrtype,only :DP,PI
       USE mmpi
       USE global_para
       USE control
       USE constants
       USE mpi_fftw3_operation
       USE matrix_inverse
       USE m_sparse
       implicit none
       integer :: i,j,k
       integer :: K_i,K_ii
       integer :: n_matrixs,nrow,ncol,n_max,nn
       complex(dp) :: pk(1:100,1)      
       complex(dp) :: q_tmp(1:100,1)      
       complex(dp) :: p_tmp(1:100,1)      
       complex(dp) :: G_C(1:100,1:100)      
       complex(dp) :: G_D(1:100,1:100)      
       real(dp) :: aa,t1,t2
       integer      job(8)
       !integer      lda, info
       CHARACTER	    UPLO
       INTEGER	    INFO, LDA, LDB, LWORK, NRHS
       INTEGER	    IPIV(100 )
       COMPLEX*16   A(100,100 ), B(100, 1 ), WORK(100*64 )
       integer,parameter :: nmaxA=550
       Complex(dp) :: GA1_val(nmaxA)      
       integer     :: GA1_col(nmaxA)      
       integer     :: GA1_row(101)      
       Complex(dp) :: GA2_val(nmaxA)      
       integer     :: GA2_col(nmaxA)      
       integer     :: GA2_row(101)      

      NRHS=1
      LDA=100  
      LDB=100  
      LWORK=100*64
      !LWORK=-1

       nrow=100
       ncol=100
       n_matrixs=1
! chose one matrix to do symbolic factorization
       nn= 1


       G_C(:,:)=dcmplx(0.0d0,0.0d0)
       do i=1,100
         do j=1,100
            if(mod(abs(i-j),10)==0) then 
            G_C(i,j)=dcmplx(i*1.0,j*1.0)
            endif
        enddo
       enddo

           G_D=G_C

          call get_sparse_size(100,100,n_max,G_C)         
       do i=1,100
         Pk(i,1)=dcmplx(1.0*i,0.0d0)  
         do j=1,100
            if(j<i) then
            G_C(i,j)=dcmplx(0.0d0,0.0d0)
            endif
        enddo
       enddo
             
           n_max=nmaxA            
       call sparse_dns_to_csr(100,100, n_max,G_C,GA1_val, GA1_col, GA1_row)
       job(1)=0
       job(2)=1
       job(3)=1
       job(4)=2
       job(5)=n_max
       job(6)=1
       lda=100
        call mkl_zdnscsr(job,100,100, G_C, lda,GA2_val, &
                            GA2_col, GA2_row, info)


        GA2_val=GA2_val-GA1_val 
        GA2_row=GA2_row-GA1_row
        GA2_col=GA2_col-GA1_col

       write(*,*) "abs of err G val=",sum(abs(GA2_val))  
       write(*,*) "abs of err G row=",sum(abs(GA2_row))  
       write(*,*) "abs of err G col=",sum(abs(GA2_col))  
        
        call mkl_zcsrsymv('u',100, GA1_val, GA1_row, &
             GA1_col,Pk(:,1), q_tmp(:,1))
      ! call mkl_zcsrgemv('N', 100, GA1_val, GA1_row, GA1_col, pk(:,1), q_tmp(:,1))

       ! q_tmp(:,1)=matmul(G_D,pk(:,1))

        do i=1,100
         do j=1,100
            p_tmp(i,1)=p_tmp(i,1)+G_D(i,j)*pk(j,1)
         enddo
        enddo  

       call sparse_csr_mv_vec(100, 100,1000, GA1_val, GA1_col, GA1_row,pk(:,1), p_tmp(:,1))

        p_tmp(:,1)=p_tmp(:,1)-q_tmp(:,1)
        
       write(*,*) "abs of err=",sum(abs(p_tmp(:,1)))  

     !  call mkl_pardiso_init(n_matrixs,nrow,ncol,nmax,GA1_inverse_row(:,nn), &
     !                          GA1_inverse_col(:,nn),GA1_inverse_val(:,nn) ,spa_pA1)

     !  call mkl_pardiso_factrz(n_matrixs,nrow,ncol,nmax,GA1_inverse_row(:,nn), &
     !                          GA1_inverse_col(:,nn),GA1_inverse_val,spa_pA1)

     !   do i=1,M_Bar_A
     !      pk(i,1)=dcmplx(i*0.1,i*0.2)
     !      q_tmp(i,1)=dcmplx(i*0.12,i*0.25)
     !   enddo 
 
         !t1=MPI_Wtime()

        !call mkl_zcsrsymv('u', M_Bar_A, GA1_inverse_val(:,nn), GA1_inverse_row(:,nn), &
        !     GA1_inverse_col(:,nn),Pk(:,1), q_tmp(:,1))

         !t2=MPI_Wtime()
        !write(*,*) "q_tmp(1,1)",q_tmp(1,1),q_tmp(2,1),"t21=",t2-t1
       !call mkl_set_num_threads( 1 )

     !    t1=MPI_Wtime()
     ! call mkl_pardiso_solve(n_matrixs,M_Bar_A,M_Bar_A,nmax,GA1_inverse_row(:,nn), &
     !                    GA1_inverse_col(:,nn), &
     !                    GA1_inverse_val,q_tmp,pk,spa_pA1)


      !   t2=MPI_Wtime()
       ! do i=1,M_Bar_A
       !   aa=aa+abs( pk(i,1)-dcmplx(i*0.1,i*0.2))
       ! enddo 
       ! write(*,*) "aa=",aa,"t21 in pardiso",t2-t1
       !  write(*,*) "pk(1,1)",pk(1,1),pk(2,1)

    !   job(1)=1 
    !   job(2)=1
    !   job(3)=1
    !   job(4)=2
    !   job(5)=nmax
    !   job(6)=1
    !   lda=nrow


    !     t1=MPI_Wtime()
    !   ! call sparse_mkl_dnscsr(M_Bar_A, M_Bar_A, nmax, G_C, 1, &
    !   !      GA1_inverse_val(:,nn), GA1_inverse_col(:,nn), GA1_inverse_row(:,nn))

    !    call mkl_zdnscsr(job, M_Bar_A, M_Bar_A, G_C, lda,GA1_inverse_val(:,nn), &
    !                        GA1_inverse_col(:,nn), GA1_inverse_row(:,nn), info)
    !    t1=MPI_Wtime()
    !    call zmat_inv(M_Bar_A, G_C)
    !     t2=MPI_Wtime()
    !    write(*,*) "t21 in direct inve",t2-t1,"err",info
    !     t1=MPI_Wtime()
    !     pk(:,1)=matmul(G_C,q_tmp(:,1))
    !     t2=MPI_Wtime()
    !    write(*,*) "t21 in direc  matmul",t2-t1,"err",info

    !    write(*,*) "pk(1,1) using inv ",pk(1,1),pk(2,1)

    !     t1=MPI_Wtime()
    !     call solver_z(G_C,q_tmp(:,1),M_Bar_A)

    end  subroutine test_sparse

   end module  Gmatrix_lu_mod       
