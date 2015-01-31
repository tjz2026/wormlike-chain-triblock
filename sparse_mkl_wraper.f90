!!! source  : sparse_mkl_wraper.f90
!!! type    : subroutine
!!! author  : Jiuzhou Tang 
!!! purpose : wraper routines to use MKL spars blas and pardiso interface.
!!!           Needs to be included in sparse_solver.f90 file
!!! status  : unstable
!!! comment : only support real(dp) and complex(dp) data types
!!!-----------------------------------------------------------------------
!! Introduction
!! ============
!! Usage
!! =====
!! -----------------------------------------------------------------------


!!>>> sparse_mkl_dnscsr_dp: converts a densely stored matrix into a row
!!>>> orientied compactly sparse matrix using intel MKL mkl_?dnscsr interface
  subroutine sparse_mkl_dnscsr_dp(nrow, ncol, nmax, dns,part, a, ja, ia)
     implicit none

! external arguments
! row dimension of dense matrix
     integer, intent(in)   :: nrow

! column dimension of dense matrix
     integer, intent(in)   :: ncol

! maximum number of nonzero elements allowed
! this should be set to be the lengths of the arrays a and ja
     integer, intent(in)   :: nmax
! which part of matrix is being used
     integer, intent(in)   :: part  ! part=0,1,2, see job(4) parameter below
! input densely stored matrix
     real(dp), intent(in)  :: dns(nrow,ncol)

! a, ja, ia, output matrix in compressed sparse row format
     integer, intent(out)  :: ia(nrow+1)
     integer, intent(out)  :: ja(nmax)
     real(dp), intent(out) :: a(nmax)

! MKL variables
     integer      job(8)
     integer      lda, info

!If job(1)=0, the rectangular matrix A is converted to the CSR format;       
    job(1)=0 
!If job(2)=1,  one-based indexing for the rectangular matrix A is used.
    job(2)=1
!if job(3)=1, one-based indexing for the matrix in CSR format is used.
    job(3)=1
!If job(4)=0, adns is a lower triangular part of matrix A;
!If job(4)=1, adns is an upper triangular part of matrix A;
!If job(4)=2, adns is a whole matrix A.
    job(4)=part
!If job(5)=nzmax - maximum number of the non-zero elements allowed if job(1)=0.
    job(5)=nmax
!If job(6)=0, only array ia is generated for the output storage.
!If job(6)>0, arrays acsr, ia, ja are generated for the output storage.
    job(6)=1

! Specifies the leading dimension of adns as declared in the calling (sub)program.
    lda=nrow

    call mkl_ddnscsr(job, nrow, ncol, dns, lda, a, ja, ia, info)
    if(info/=0) then
    write(*,*) "unsucessful in convert dns to csr using mkl,error in",info,"line"
    stop
    endif

  end subroutine sparse_mkl_dnscsr_dp

!!>>> sparse_mkl_dnscsr_z: converts a densely stored complex matrix into a row
!!>>> orientied compactly sparse matrix using intel MKL mkl_?dnscsr interface
  subroutine sparse_mkl_dnscsr_z(nrow, ncol, nmax, dns, part, a, ja, ia)
     implicit none
     integer, intent(in)   :: nrow
     integer, intent(in)   :: ncol
     integer, intent(in)   :: nmax
     integer, intent(in)   :: part  ! part=0,1,2, see job(4) parameter below
     complex(dp), intent(in)  :: dns(nrow,ncol)
     integer, intent(out)  :: ia(nrow+1)
     integer, intent(out)  :: ja(nmax)
     complex(dp), intent(out) :: a(nmax)
     integer      job(8)
     integer      lda, info

    job(1)=0 
    job(2)=1
    job(3)=1
    job(4)=part
    job(5)=nmax
    job(6)=1

    lda=nrow

    call mkl_zdnscsr(job, nrow, ncol, dns, lda, a, ja, ia, info)
    if(info/=0) then
    write(*,*) "unsucessful in convert dns to csr using mkl,error in",info,"line"
    stop
    endif

  end subroutine sparse_mkl_dnscsr_z

!!! Tips : wraper for symmetrical matrix vector multiplation is not necessary to be used becasue
!!! the original mkl interface is simple and clear enough to be called directly. This wraper is just
!!! for the benefit of showing off the fortran features of function overloading. 
!!!>>> sparse_mkl_symv_dp: Computes matrix - vector product of a sparse symmetrical matrix stored in 
!!! the CSR format (3-array variation) with one-based indexing.
  subroutine sparse_mkl_symv_dp(nrow,ncol,nmax,ia,ja,a,uplo,x,y)
     implicit none
! external arguments
! row dimension of dense matrix
     integer, intent(in)   :: nrow

! column dimension of dense matrix
     integer, intent(in)   :: ncol

! maximum number of nonzero elements allowed
! this should be set to be the lengths of the arrays a and ja
     integer, intent(in)   :: nmax

! a, ja, ia, input matrix in compressed sparse row format
     integer, intent(in)   :: ia(nrow+1)
     integer, intent(in)   :: ja(nmax)
     real(dp), intent(in)  :: a(nmax)

! vector, length equal to the column dimension of the dense matrix
     real(dp), intent(in)  :: x(ncol)
! vector, real array of length nrow, containing the product y = A . x
     real(dp), intent(out) :: y(nrow)
! character*1,If uplo = 'U' or 'u', then the upper triangle of the matrix A is used.
!If uplo = 'L' or 'l', then the low triangle of the matrix A is used.
     CHARACTER*1,intent(in) :: uplo 
     
    call mkl_dcsrsymv(uplo, nrow, a, ia, ja, x, y)

  end subroutine sparse_mkl_symv_dp

!!!>>> sparse_mkl_symv_z: Computes matrix - vector product of a sparse symmetrical complex matrix
!!! stored in the CSR format (3-array variation) with one-based indexing. 
!!! Warning: it is symmetric as a_ij=a_ji, not hermite as a_ij=a_ji*.

  subroutine sparse_mkl_symv_z(nrow,ncol,nmax,ia,ja,a,uplo,x,y)
     implicit none
     integer, intent(in)   :: nrow
     integer, intent(in)   :: ncol
     integer, intent(in)   :: nmax
     integer, intent(in)   :: ia(nrow+1)
     integer, intent(in)   :: ja(nmax)
     complex(dp), intent(in)  :: a(nmax)
     complex(dp), intent(in)  :: x(ncol)
     complex(dp), intent(out) :: y(nrow)
     CHARACTER*1,intent(in) :: uplo 
     
    call mkl_zcsrsymv(uplo, nrow, a, ia, ja, x, y)

  end subroutine sparse_mkl_symv_z


!!!>>> sparse_mkl_solve_dp: Solve the CSR linear equation for double precision
  subroutine sparse_mkl_solve_dp(nrow,ncol,nmax,ia,ja,a,sym_info,x,y)
    implicit none
! external arguments
! row dimension of dense matrix
    integer, intent(in)   :: nrow
! column dimension of dense matrix
    integer, intent(in)   :: ncol
! maximum number of nonzero elements allowed
! this should be set to be the lengths of the arrays a and ja
    integer, intent(in)   :: nmax
! a, ja, ia, input matrix in compressed sparse row format
    integer, intent(in)   :: ia(nrow+1)
    integer, intent(in)   :: ja(nmax)
    real(dp), intent(in)  :: a(nmax)
! whether matrix is symmetric or not
    logical,intent(in),optional :: sym_info
! solution vector, length equal to the column dimension of the dense matrix
    real(dp), intent(out)  :: x(ncol)
! rho vector, real array of length nrow, containing the product y = A*x
    real(dp), intent(in) :: y(nrow)
! whether normal or transpose matrix used in solver,ie. transa='N' or 'T'
    character*1 :: transa
    real(dp) :: alpha
! matrix structure of a ,see MKL document for detailed infomation
! https://software.intel.com/en-us/node/520801#TBL2-6
     character  :: matdescra(6)           
     transa='N' 
     alpha=1.0d0
     
     matdescra(1)='G' !general matrix
     matdescra(2)='L' !ignored
     matdescra(3)='N' !ignored
     matdescra(4)='F' !one based indexing 
     if(present(sym_info)) then
     matdescra(1)='S' !symmetric matrix
     matdescra(2)='U' !upper trangle part 
     matdescra(3)='N' !non unit
     matdescra(4)='F' !one based indexing 
     endif    
 
!!!MKL note
!!!Comparing the array rowIndex from the Table "Storage Arrays for a Non-Symmetric Example Matrix" 
! with the arrays pointerB and pointerE from the Table "Storage Arrays for an Example Matrix in CSR Format"
! it is easy to see that
! pointerB(i) = rowIndex(i) for i=1, ..5;
! pointerE(i) = rowIndex(i+1) for i=1, ..5.
!This enables calling a routine that has values, columns, pointerB and pointerE as input parameters for a 
! sparse matrix stored in the format accepted for the direct sparse solvers. For example, a routine with the interface:
! subroutine name_routine(.... ,  values, columns, pointerB, pointerE, ...)
! can be called with parameters values, columns, rowIndex as follows:
! call name_routine(.... ,  values, columns, rowIndex, rowindex(2), ...).

! Warning,according to MKL online document,alpha is scalar, x and y are vectors, 
! A is a sparse upper or lower triangular matrix with unit or non-unit main diagonal, A' is the transpose of A.
! ??? Is matrix A has to be a upper or lower triangle matrix ? Check before using the routine!!!!! 
! X=inv(A)*Y

   call mkl_dcsrsv(transa, ncol, alpha, matdescra, a, ja, ia, ia(2), y, x)


  end subroutine sparse_mkl_solve_dp


!!!>>> sparse_mkl_solve_z: Solve the CSR linear equation for double precision complex
  subroutine sparse_mkl_solve_z(nrow,ncol,nmax,ia,ja,a,sym_info,x,y)
    implicit none
    integer, intent(in)   :: nrow
    integer, intent(in)   :: ncol
    integer, intent(in)   :: nmax
    integer, intent(in)   :: ia(nrow+1)
    integer, intent(in)   :: ja(nmax)
    complex(dp), intent(in)  :: a(nmax)
    logical,intent(in),optional :: sym_info
    complex(dp), intent(out)  :: x(ncol)
    complex(dp), intent(in) :: y(nrow)
    character*1 :: transa
    complex(dp) :: alpha
    character  :: matdescra(6)           
     transa='N' 
     alpha=cmplx(1.0d0,0.0d0)
     
     matdescra(1)='G' !general matrix
     matdescra(2)='L' !ignored
     matdescra(3)='N' !ignored
     matdescra(4)='F' !one based indexing 
     if(present(sym_info)) then
     matdescra(1)='S' !symmetric matrix
     matdescra(2)='U' !upper trangle part 
     matdescra(3)='N' !non unit
     matdescra(4)='F' !one based indexing 
     endif    

   call mkl_zcsrsv(transa, ncol, alpha, matdescra, a, ja, ia, ia(2), y, x)


  end subroutine sparse_mkl_solve_z



!!!>>> mkl_pardiso_init_dp : init the mkl pardiso solver for a series of double real matrixs
! that have a identical sparity structure. Assuming matrixs are symmetric. 
  subroutine mkl_pardiso_init_dp(n_matrixs,nrow,ncol,nmax,ia,ja,a,spa_p)
    USE mkl_pardiso
    IMPLICIT NONE
    INTEGER, PARAMETER :: dp = KIND(1.0D0)
    ! derived type of pardiso
    TYPE(pardiso_plan),intent(inout) :: spa_p
    ! input variables
    ! number of structurally identical sparse matrixs
    integer,intent(in) :: n_matrixs
    ! row number and column number, number of nonzero entries.
    integer,intent(in) :: nrow,ncol,nmax
    ! a, ja, ia, output matrix in compressed sparse row format
    Integer, intent(in)  :: ia(1:nrow+1)
    Integer, intent(in)  :: ja(1:nmax)
    Real(dp), intent(in) :: a(1:nmax)
    integer :: i
    ! set up basic parameters for pardiso.

     spa_p%nrhs=1
     spa_p%maxfct=n_matrixs
     do i=1,64
     spa_p%iparm(i)=0
     enddo
    ! the following parameters are set according to the sample code from mkl
    ! to use it on an expert level, more detailed infos can be found on intel's mkl documents.
    !Warning: do not change the following parameters unless you know what you are doing.
     spa_p%iparm(1) = 1 ! no solver default
     spa_p%iparm(2) = 2 ! fill-in reordering from METIS
     spa_p%iparm(4) = 0 ! no iterative-direct algorithm
     spa_p%iparm(5) = 0 ! no user fill-in reducing permutation
     spa_p%iparm(6) = 0 ! =0 solution on the first n compoments of x
     spa_p%iparm(8) = 9 ! numbers of iterative refinement steps
     spa_p%iparm(10) = 13 ! perturbe the pivot elements with 1E-13
     spa_p%iparm(11) = 1 ! use nonsymmetric permutation and scaling MPS
     spa_p%iparm(13) = 0 ! maximum weighted matching algorithm is switched-off (default for symmetric)
     spa_p%iparm(14) = 0 ! Output: number of perturbed pivots
     spa_p%iparm(18) = -1 ! Output: number of nonzeros in the factor LU
     spa_p%iparm(19) = -1 ! Output: Mflops for LU factorization
     spa_p%iparm(20) = 0 ! Output: Numbers of CG Iterations

     spa_p%error  = 0 ! initialize error flag
     spa_p%msglvl = 1 ! print statistical information
     spa_p%mtype  = -2 ! symmetric, indefinite
 

!.. Initiliaze the internal solver memory pointer. This is only
! necessary for the FIRST call of the PARDISO solver.
     do i = 1, 64
     spa_p%pt( i )%DUMMY =  0 
     end do

!.. Reordering and Symbolic Factorization, This step also allocates
! all memory that is necessary for the factorization
     spa_p%phase = 11 ! only reordering and symbolic factorization
! Because all input matrixs are assumed to be sparsity identical, only 
! the first matrix is used to do the symbolic factorization.
     spa_p%mnum=1
    Call pardiso (spa_p%pt, spa_p%maxfct, spa_p%mnum, spa_p%mtype, spa_p%phase, &
                  ncol, a, ia, ja, spa_p%idum,spa_p%nrhs, spa_p%iparm, &
                  spa_p%msglvl, spa_p%ddum_r, spa_p%ddum_r,spa_p%error)
    WRITE(*,*) 'Reordering completed ... '
    IF (spa_p%error /= 0) THEN
       WRITE(*,*) 'The following ERROR was detected: ', spa_p%error
       Call mkl_pardiso_clean(ncol,spa_p)
       Stop
    END IF
    WRITE(*,*) 'Number of nonzeros in factors = ',spa_p%iparm(18)
    WRITE(*,*) 'Number of factorization MFLOPS = ',spa_p%iparm(19)

  end subroutine mkl_pardiso_init_dp  


!!!>>> mkl_pardiso_factrz_dp : doing numerical factorization for all sparsity identical matrixs
  subroutine mkl_pardiso_factrz_dp(n_matrixs,nrow,ncol,nmax,ia,ja,a_all,spa_p)
    USE mkl_pardiso
    IMPLICIT NONE
    INTEGER, PARAMETER :: dp = KIND(1.0D0)
    TYPE(pardiso_plan),intent(inout) :: spa_p
    integer,intent(in) :: n_matrixs
    integer,intent(in) :: nrow,ncol,nmax
    Integer, intent(in)  :: ia(1:nrow+1)
    Integer, intent(in)  :: ja(1:nmax)
!All the sparsity identical matrixs stored in CSR format     
    Real(dp), intent(in) :: a_all(1:nmax,1:n_matrixs)
    Integer :: i
    
!.. Factorization phase.
    spa_p%phase = 22 ! only numerical factorization
    if(spa_p%maxfct/=n_matrixs) then
    write(*,*) "Error, sparse plan doesn't match the input matrixs number"
       Call mkl_pardiso_clean(ncol,spa_p)
    stop
    endif     

     do i=1,spa_p%maxfct
        spa_p%mnum=i
     call pardiso (spa_p%pt, spa_p%maxfct, spa_p%mnum, spa_p%mtype, & 
                   spa_p%phase, ncol, a_all(:,i), ia, ja,  & 
                   spa_p%idum,spa_p%nrhs, spa_p%iparm, spa_p%msglvl,& 
                   spa_p%ddum_r, spa_p%ddum_r,spa_p%error)
     !WRITE(*,*) 'Factorization completed ... '
      IF (spa_p%error /= 0) THEN
        WRITE(*,*) 'The following ERROR was detected: ', spa_p%error
       Call mkl_pardiso_clean(ncol,spa_p)
        Stop
      ENDIF
     enddo    

  end subroutine mkl_pardiso_factrz_dp


!!!>>> mkl_pardiso_solve_dp : solving the sparse linear equations.
  subroutine mkl_pardiso_solve_dp(n_matrixs,nrow,ncol,nmax,ia,ja,a_all,rhs_all,sol_all,spa_p)
    USE mkl_pardiso
    IMPLICIT NONE
    INTEGER, PARAMETER :: dp = KIND(1.0D0)
    TYPE(pardiso_plan),intent(inout) :: spa_p
    Integer,intent(in) :: n_matrixs
    Integer,intent(in) :: nrow,ncol,nmax
    Integer,intent(in)  :: ia(1:nrow+1)
    Integer,intent(in)  :: ja(1:nmax)
    Real(dp),intent(in) :: a_all(1:nmax,1:n_matrixs)
! all the right hand vectors 
    Real(dp),intent(inout) :: rhs_all(1:ncol,1:n_matrixs) 
! return all the solution vectors
    Real(dp),intent(inout) :: sol_all(1:ncol,1:n_matrixs) 
    Integer :: i

    ! Back substitution and iterative refinement
    ! Why iparm(8) is set to be 2 in the mkl sample code, find out later, right now just go with it!
    spa_p%iparm(8) = 2 ! max numbers of iterative refinement steps
    spa_p%phase = 33 ! Solve, iterative refinement

    do i=1,spa_p%maxfct
     spa_p%mnum=i 
     call pardiso (spa_p%pt, spa_p%maxfct, spa_p%mnum, spa_p%mtype,spa_p%phase, &
                   ncol, a_all(:,i), ia, ja,spa_p%idum, spa_p%nrhs, spa_p%iparm, &
                   spa_p%msglvl, rhs_all(:,i), sol_all(:,i),spa_p%error)
   ! IF (error /= 0) THEN
   ! ENDIF
   ! WRITE(*,*) 'The solution of the system is '
   
   enddo       
  
  end subroutine mkl_pardiso_solve_dp





!!!>>> mkl_pardiso_clean : Termination and release of memory
  subroutine mkl_pardiso_clean(ncol,spa_p)
    USE mkl_pardiso
    IMPLICIT NONE
    INTEGER, PARAMETER :: dp = KIND(1.0D0)
    TYPE(pardiso_plan),intent(inout) :: spa_p
    Integer,intent(in) :: ncol
    integer :: i

    do i=1,spa_p%maxfct
    spa_p%mnum=i    
    spa_p%phase = -1 ! release internal memory
    call pardiso (spa_p%pt, spa_p%maxfct, spa_p%mnum, spa_p%mtype, spa_p%phase, &
                  ncol, spa_p%ddum_z, spa_p%idum, spa_p%idum, &
    spa_p%idum, spa_p%nrhs, spa_p%iparm, spa_p%msglvl, spa_p%ddum_z, spa_p%ddum_z,spa_p%error)
    enddo  

  end subroutine mkl_pardiso_clean


!!!>>> mkl_pardiso_init_z : init the mkl pardiso solver for a series of double complex matrixs
! that have a identical sparity structure. Assuming matrixs are symmetric. 
  subroutine mkl_pardiso_init_z(n_matrixs,nrow,ncol,nmax,ia,ja,a,spa_p)
    USE mkl_pardiso
    IMPLICIT NONE
    INTEGER, PARAMETER :: dp = KIND(1.0D0)
    ! derived type of pardiso
    TYPE(pardiso_plan),intent(inout) :: spa_p
    ! input variables
    ! number of structurally identical sparse matrixs
    integer,intent(in) :: n_matrixs
    ! row number and column number, number of nonzero entries.
    integer,intent(in) :: nrow,ncol,nmax
    ! a, ja, ia, output matrix in compressed sparse row format
    Integer, intent(in)  :: ia(1:nrow+1)
    Integer, intent(in)  :: ja(1:nmax)
    Complex(dp), intent(in) :: a(1:nmax)
    integer :: i 
    ! set up basic parameters for pardiso.

     spa_p%nrhs=1
     spa_p%maxfct=n_matrixs
     do i=1,64
     spa_p%iparm(i)=0
     enddo
    ! the following parameters are set according to the sample code from mkl
    ! to use it on an expert level, more detailed infos can be found on intel's mkl documents.
    !Warning: do not change the following parameters unless you know what you are doing.
     spa_p%iparm(1) = 1 ! no solver default
     spa_p%iparm(2) = 2 ! fill-in reordering from METIS
     spa_p%iparm(4) = 0 ! no iterative-direct algorithm
     spa_p%iparm(5) = 0 ! no user fill-in reducing permutation
     spa_p%iparm(6) = 0 ! =0 solution on the first n compoments of x
     spa_p%iparm(8) = 9 ! numbers of iterative refinement steps
     spa_p%iparm(10) = 13 ! perturbe the pivot elements with 1E-13
     spa_p%iparm(11) = 1 ! use nonsymmetric permutation and scaling MPS
     spa_p%iparm(13) = 0 ! maximum weighted matching algorithm is switched-off (default for symmetric)
     spa_p%iparm(14) = 0 ! Output: number of perturbed pivots
     spa_p%iparm(18) = -1 ! Output: number of nonzeros in the factor LU
     spa_p%iparm(19) = -1 ! Output: Mflops for LU factorization
     spa_p%iparm(20) = 0 ! Output: Numbers of CG Iterations
     spa_p%iparm(27)=1 ! check out the  data structure of input matrix
     spa_p%error  = 0 ! initialize error flag
     spa_p%msglvl = 1 ! print statistical information/ 0 no info print out
     spa_p%mtype  = 6 ! complex and symmetric
 

!.. Initiliaze the internal solver memory pointer. This is only
! necessary for the FIRST call of the PARDISO solver.
     do i = 1, 64
     spa_p%pt( i )%DUMMY =  0 
     end do

!.. Reordering and Symbolic Factorization, This step also allocates
! all memory that is necessary for the factorization
     spa_p%phase = 11 ! only reordering and symbolic factorization
! Because all input matrixs are assumed to be sparsity identical, only 
! the first matrix is used to do the symbolic factorization.
     spa_p%mnum=1

    call pardiso (spa_p%pt, spa_p%maxfct, spa_p%mnum, spa_p%mtype, &
                spa_p%phase,ncol, a, ia, ja, spa_p%idum,spa_p%nrhs,&
     spa_p%iparm,spa_p%msglvl, spa_p%ddum_z, spa_p%ddum_z, spa_p%error)


    WRITE(*,*) 'Reordering completed ... '
    IF (spa_p%error /= 0) THEN
       WRITE(*,*) 'The following ERROR was detected: ', spa_p%error
       Call mkl_pardiso_clean(ncol,spa_p)
       Stop
    END IF
    WRITE(*,*) 'Number of nonzeros in factors = ',spa_p%iparm(18)
    WRITE(*,*) 'Number of factorization MFLOPS = ',spa_p%iparm(19)

  end subroutine mkl_pardiso_init_z


!!!>>> mkl_pardiso_factrz_z : doing numerical factorization for all sparsity identical matrixs
  subroutine mkl_pardiso_factrz_z(n_matrixs,nrow,ncol,nmax,ia,ja,a_all,spa_p)
    USE mkl_pardiso
    IMPLICIT NONE
    INTEGER, PARAMETER :: dp = KIND(1.0D0)
    TYPE(pardiso_plan),intent(inout) :: spa_p
    integer,intent(in) :: n_matrixs
    integer,intent(in) :: nrow,ncol,nmax
    Integer, intent(in)  :: ia(1:nrow+1)
    Integer, intent(in)  :: ja(1:nmax)
!All the sparsity identical matrixs stored in CSR format     
    Complex(dp), intent(in) :: a_all(1:nmax,1:n_matrixs)
    Integer :: i
!.. Factorization phase.
    spa_p%phase = 22 ! only numerical factorization
    if(spa_p%maxfct/=n_matrixs) then
    write(*,*) "Error, sparse plan doesn't match the input matrixs number"
       Call mkl_pardiso_clean(ncol,spa_p)
    stop
    endif     

     do i=1,spa_p%maxfct
        spa_p%mnum=i
     call pardiso (spa_p%pt, spa_p%maxfct, spa_p%mnum, spa_p%mtype, & 
                   spa_p%phase, ncol, a_all(:,i), ia, ja,  & 
                   spa_p%idum,spa_p%nrhs, spa_p%iparm, spa_p%msglvl,& 
                  spa_p%ddum_z, spa_p%ddum_z,spa_p%error)
     !WRITE(*,*) 'Factorization completed ... '
      IF (spa_p%error /= 0) THEN
        WRITE(*,*) 'The following ERROR was detected: ', spa_p%error
       Call mkl_pardiso_clean(ncol,spa_p)
        Stop
      ENDIF
     enddo    

  end subroutine mkl_pardiso_factrz_z


!!!>>> mkl_pardiso_solve_z : solving the sparse linear equations.
  subroutine mkl_pardiso_solve_z(n_matrixs,nrow,ncol,nmax,ia,ja,a_all,rhs_all,sol_all,spa_p)
    USE mkl_pardiso
    IMPLICIT NONE
    INTEGER, PARAMETER :: dp = KIND(1.0D0)
    TYPE(pardiso_plan),intent(inout) :: spa_p
    Integer,intent(in) :: n_matrixs
    Integer,intent(in) :: nrow,ncol,nmax
    Integer,intent(in)  :: ia(1:nrow+1)
    Integer,intent(in)  :: ja(1:nmax)
    Complex(dp),intent(in) :: a_all(1:nmax,1:n_matrixs)
! all the right hand vectors 
    Complex(dp),intent(inout) :: rhs_all(1:ncol,1:n_matrixs) 
! return all the solution vectors
    Complex(dp),intent(inout) :: sol_all(1:ncol,1:n_matrixs) 
    Integer :: i

    ! Back substitution and iterative refinement
    ! Why iparm(8) is set to be 2 in the mkl sample code, find out later, right now just go with it!
    spa_p%iparm(8) = 3 ! max numbers of iterative refinement steps
    spa_p%phase = 33 ! Solve, iterative refinement

    do i=1,spa_p%maxfct
     spa_p%mnum=i 
     call pardiso (spa_p%pt, spa_p%maxfct, spa_p%mnum, spa_p%mtype,spa_p%phase,  &
                   ncol, a_all(:,i), ia, ja,spa_p%idum, spa_p%nrhs, spa_p%iparm, &
                   spa_p%msglvl, rhs_all(:,i), sol_all(:,i),spa_p%error)
   ! IF (error /= 0) THEN
   ! ENDIF
   ! WRITE(*,*) 'The solution of the system is '
   
   enddo       
  
  end subroutine mkl_pardiso_solve_z


