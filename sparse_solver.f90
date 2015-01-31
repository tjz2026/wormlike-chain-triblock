!!!-----------------------------------------------------------------------
!!! source  : m_sparse.f90
!!! type    : module
!!! author  : Jiuzhou Tang 
!!! purpose : the purpose of this module is to implement important sparse
!!!           matrix  operations, including sparse linear solver, matrix 
!!!           multiplication, LU factorization,format conversion, etc. 
!!            the internal format of sparse matrix
!!!           used in this module is CSR (compressed sparse row) format.
!!!           It currently relies on the external intel MKL sparse routines.
!!!           ( sparse blas and pardiso ) 
!!! status  : unstable
!!! comment : only support real(dp) and complex(dp) data types
!!!-----------------------------------------------------------------------
!! Introduction
!! ============
!! In this module, we implement some basic sparse matrix algebra. Now it
!! supports double precision real and complex numbers.
!!
!! Usage
!! =====
!!
!! 1. import sparse support
!! ------------------------
!!
!! use sparse
!!
!! 2. convert normal matrix to sparse matrix
!! -----------------------------------------
!!
!! call sparse_csr_to_dns(...)
!!
!!
!! 3. convert sparse matrix to normal matrix
!! -----------------------------------------
!!
!! call sparse_dns_to_csr(...)
!!
!! 4. perform sparse matrix - vector multiplication
!! ------------------------------------------------
!!
!! call sparse_csr_mv_vec(...)
!!
!! 5. perform sparse linear equation solver
!! ------------------------------------------------
!!
!! call sparse_csr_solve(...)
!!

!# if defined (MKL_lib)
  INCLUDE 'mkl_pardiso.f90'
!# endif /* MKL_lib */

  module m_sparse
!# if defined (MKL_lib)
  USE mkl_pardiso
!# endif /* MKL_lib */
     implicit none

!!========================================================================
!!>>> declare global parameters                                        <<<
!!========================================================================

! dp: number precision, double precision for real and complex number
     integer, private, parameter :: dp    = kind(1.0d0)

! mystd: device descriptor, console output
     integer, private, parameter :: mystd = 6

!!========================================================================
!!>>> declare accessibility for module routines                        <<<
!!========================================================================
! get nonzero number of a dense matrix
     private :: get_sparse_size_dp     !real(dp) version
     private :: get_sparse_size_z      !complex(dp) version


! CSR -> DNS
     private :: sparse_format_csrdns   ! real(dp) version
     private :: sparse_format_csrdns_z ! complex(dp) version

! DNS -> CSR
     private :: sparse_format_dnscsr   ! real(dp) version
     private :: sparse_format_dnscsr_z ! complex(dp) version
! CSR X VEC
     private :: sparse_matmul_amuvec   ! real(dp) version
     private :: sparse_matmul_amuvec_z ! complex(dp) version
!CSR general linear solver
     !private :: sparse_gen_solve   ! real(dp) version
     !private :: sparse_gen_solve_z ! complex(dp) version
 
!CSR symmetric linear solver(symmetric for real(dp) and hermite for complex)
     !private :: sparse_sym_solve   ! real(dp) version
     !private :: sparse_he_solve ! complex(dp) version

!Wrapers for using MKL sparse blas routines 
! convert dns and csr format using MKL's standard routines
! 
!# if defined (MKL_lib)
! CSR <-> DNS
     private :: sparse_mkl_dnscsr_dp !real(dp) version
     private :: sparse_mkl_dnscsr_z  !complex(dp) version
! Symmetric CSR matrix vector multiply
     private :: sparse_mkl_symv_dp
     private :: sparse_mkl_symv_z
! CSR linear solver from MKL sparse blas level2 & level3
     private :: sparse_mkl_solve_dp
     private :: sparse_mkl_solve_z

! PARDISO Interfaces:
! In wormlike SCFT application, only symmetric 
! sparse matrix are repeatedly used during computation.
! thus currently, we only foucus on symmetric dp and complex matrix.

! Init all the variables need for pardiso 
! and doing the reordering and symbolic factorization
     private :: mkl_pardiso_init_dp
     private :: mkl_pardiso_init_z
! Doing the numerical factorization and store all needed variables
     private :: mkl_pardiso_factrz_dp
     private :: mkl_pardiso_factrz_z
! Solve linear equation for a given rhs vector using the factorized matrixs
     private :: mkl_pardiso_solve_dp
     private :: mkl_pardiso_solve_z
! clean up and release the internal memory.
     private :: mkl_pardiso_clean


      type pardiso_plan
!Internal solver memory pointer, integer*8  
      TYPE(MKL_PARDISO_HANDLE)  :: pt(64)
!Maximum number of factors with identical sparsity structure that must be kept in memory at the same time.      
      integer :: maxfct
!Number of right-hand sides that need to be solved for.
      integer :: nrhs
! Array, dimension (64). This array is used to pass various parameters to
! Intel MKL PARDISO and to return some useful information after execution of the solver.
      integer :: iparm(64)
!Indicates the actual matrix for the solution phase.
      integer :: mnum
!Message level information.
      integer :: msglvl
!Defines the matrix type, which influences the pivoting method.
!The Intel MKL PARDISO solver supports the following matrices:
! 1 real and structurally symmetric
! 2 real and symmetric positive definite
!-2 real and symmetric indefinite
! 3 complex and structurally symmetric
! 4 complex and Hermitian positive definite
!-4 complex and Hermitian indefinite
! 6 complex and symmetric
!11 real and nonsymmetric
!13 complex and nonsymmetric     
      integer :: mtype
!Controls the execution of the solver. 
      integer :: phase           
! error infomation
      integer :: error
      integer :: idum(1)
      real(dp) :: ddum_r(1) 
      complex(dp) :: ddum_z(1) 
    end type pardiso_plan
    ! pardiso plan for double real matrixs
    type(pardiso_plan) :: spa_p_dp
    ! pardiso plan for double complex matrixs
    type(pardiso_plan) :: spa_p_z

!# endif /* MKL_lib */


!!========================================================================
!!>>> declare public interface and module procedure                           <<<
!!========================================================================

     public :: get_sparse_size
     interface get_sparse_size
         module procedure get_sparse_size_dp
         module procedure get_sparse_size_z
     end interface get_sparse_size

     public :: sparse_csr_to_dns
     interface sparse_csr_to_dns
         module procedure sparse_format_csrdns
         module procedure sparse_format_csrdns_z
     end interface sparse_csr_to_dns

     public :: sparse_dns_to_csr
     interface sparse_dns_to_csr
         module procedure sparse_format_dnscsr
         module procedure sparse_format_dnscsr_z
     end interface sparse_dns_to_csr

     public :: sparse_csr_mv_vec
     interface sparse_csr_mv_vec
         module procedure sparse_matmul_amuvec
         module procedure sparse_matmul_amuvec_z
     end interface sparse_csr_mv_vec

   ! public :: sparse_gen_solve_csr
   ! interface sparse_gen_solve_csr
   !     module procedure sparse_gen_solve
   !     module procedure sparse_gen_solve_z
   ! end interface sparse_gen_solve_csr

   ! public :: sparse_sym_solve_csr
   ! interface sparse_sym_solve_csr
   !     module procedure sparse_sym_solve
   !     module procedure sparse_he_solve
   ! end interface sparse_sym_solve_csr


!# if defined (MKL_lib)
     public :: sparse_mkl_dnscsr
     interface sparse_mkl_dnscsr
         module procedure sparse_mkl_dnscsr_dp
         module procedure sparse_mkl_dnscsr_z
     end interface sparse_mkl_dnscsr
!     public :: sparse_mkl_dnscsr
!     interface sparse_mkl_dnscsr
!         module procedure sparse_mkl_dnscsr_dp
!         module procedure sparse_mkl_dnscsr_z
!     end interface sparse_mkl_dnscsr
     public :: sparse_mkl_solve
     interface sparse_mkl_solve
         module procedure sparse_mkl_solve_dp
         module procedure sparse_mkl_solve_z
     end interface sparse_mkl_solve

     public :: mkl_pardiso_init
     interface mkl_pardiso_init
         module procedure mkl_pardiso_init_dp
         module procedure mkl_pardiso_init_z
     end interface mkl_pardiso_init

     public :: mkl_pardiso_factrz
     interface mkl_pardiso_factrz
         module procedure mkl_pardiso_factrz_dp
         module procedure mkl_pardiso_factrz_z
     end interface mkl_pardiso_factrz

     public :: mkl_pardiso_solve
     interface mkl_pardiso_solve
         module procedure mkl_pardiso_solve_dp
         module procedure mkl_pardiso_solve_z
     end interface mkl_pardiso_solve

!# endif /* MKL_lib */


  contains ! encapsulated functionality

!!>>> get_sparse_size: get the number of nonzero entries of a dp precision
!!dense matrix a
  subroutine get_sparse_size_dp(nrow,ncol,nmax,a,sym_info,diag_info)
     implicit none

! external arguments
! row dimension of dense matrix
     integer, intent(in)   :: nrow
! column dimension of dense matrix
     integer, intent(in)   :: ncol
! maximum number of nonzero elements allowed
     integer, intent(out)   :: nmax
     real(dp),intent(in) :: a(nrow,ncol)
! check whether the matrix is symmetric or not
     logical,intent(out),optional :: sym_info
! check whether there is any diagonal element of the value zero
     logical,intent(out),optional :: diag_info
     real(dp),parameter :: THRESH_matrix=1.0e-15,THRESH_sym=1.0e-16
     integer :: i,j,k
     real(dp) :: sum_err
  
        if(present(sym_info)) sym_info=.false.
        if(present(diag_info)) diag_info=.false.

         k=0
         do j=1,ncol
           do i=1,nrow
             if(abs(a(i,j))>=THRESH_matrix) then
             k=k+1
             endif
           enddo
         enddo
         nmax=k

         if(present(sym_info)) then
         sum_err=0.0d0
          do j=1,ncol
           do i=j,nrow
              sum_err=sum_err+abs(a(i,j)-a(j,i))
           enddo
          enddo
           if(sum_err<=THRESH_sym) then
             sym_info=.true.
           endif
         endif

         if(present(diag_info)) then
          do j=1,ncol
            if(abs(a(j,j))<=THRESH_matrix) then
             diag_info=.true.
             exit
            endif
          enddo
         endif 

         return
  end subroutine get_sparse_size_dp

!!>>> get_sparse_size_z: get the number of nonzero entries of a dp complex 
!!dense matrix a
  subroutine get_sparse_size_z(nrow,ncol,nmax,a,sym_info,diag_info)
     implicit none

! external arguments
! row dimension of dense matrix
     integer, intent(in)   :: nrow
! column dimension of dense matrix
     integer, intent(in)   :: ncol
! maximum number of nonzero elements allowed
     integer,intent(out) :: nmax
     complex(dp),intent(in) :: a(nrow,ncol)
! check whether the matrix is symmetric or not
     logical,intent(out),optional :: sym_info
! check whether there is any diagonal element of the value zero
     logical,intent(out),optional :: diag_info
     real(dp),parameter :: THRESH_matrix=1.0e-18,THRESH_sym=1.0e-18
     real(dp) :: sum_err
     integer :: i,j,k

          k=0
          do j=1,ncol
            do i=1,nrow
             if(abs(a(i,j))>=THRESH_matrix) then
              k=k+1
             endif
            enddo
          enddo
          nmax=k

         if(present(sym_info)) then
         sum_err=0.0d0
          do j=1,ncol
           do i=j,nrow
              sum_err=sum_err+abs(a(i,j)-a(j,i))
           enddo
          enddo
           if(sum_err<=THRESH_sym) then
             sym_info=.true.
           endif
         endif

         if(present(diag_info)) then
          do j=1,ncol
            if(abs(a(j,j))<=THRESH_matrix) then
             diag_info=.true.
             exit
            endif
          enddo
         endif 

          return
  end subroutine get_sparse_size_z

!!>>> sparse_format_csrdns: converts a row-stored sparse matrix into a
!!>>> densely stored one
  subroutine sparse_format_csrdns(nrow, ncol, nmax, a, ja, ia, dns)
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

! array where to store dense matrix
     real(dp), intent(out) :: dns(nrow,ncol)

! local variables
! loop index
     integer :: i
     integer :: j
     integer :: k

! init dns matrix
     dns = 0.0_dp

! convert sparse matrix to dense matrix
     do i=1,nrow
         do k=ia(i),ia(i+1)-1
             j = ja(k)
             if ( j > ncol ) then
                 write(mystd,'(a)') 'sparse: error in sparse_format_csrdns'
                 STOP
             endif ! back if ( j > ncol ) block
             dns(i,j) = a(k)
         enddo ! over k={ia(i),ia(i+1)-1} loop
     enddo ! over i={1,nrow} loop

     return
  end subroutine sparse_format_csrdns

!!>>> sparse_format_csrdns_z: converts a row-stored sparse matrix into a
!!>>> densely stored one
  subroutine sparse_format_csrdns_z(nrow, ncol, nmax, sa, ja, ia, dns)
     implicit none

! external arguments
! row dimension of dense matrix
     integer, intent(in)      :: nrow

! column dimension of dense matrix
     integer, intent(in)      :: ncol

! maximum number of nonzero elements allowed
! this should be set to be the lengths of the arrays sa and ja
     integer, intent(in)      :: nmax

! sa, ja, ia, input matrix in compressed sparse row format
     integer, intent(in)      :: ia(nrow+1)
     integer, intent(in)      :: ja(nmax)
     complex(dp), intent(in)  :: sa(nmax)

! array where to store dense matrix
     complex(dp), intent(out) :: dns(nrow,ncol)

! local variables
! loop index
     integer :: i
     integer :: j
     integer :: k

! init dns matrix
     dns = dcmplx(0.0_dp, 0.0_dp)

! convert sparse matrix to dense matrix
     do i=1,nrow
         do k=ia(i),ia(i+1)-1
             j = ja(k)
             if ( j > ncol ) then
                 write(mystd,'(a)') 'sparse: error in sparse_format_csrdns_z'
                 STOP
             endif ! back if ( j > ncol ) block
             dns(i,j) = sa(k)
         enddo ! over k={ia(i),ia(i+1)-1} loop
     enddo ! over i={1,nrow} loop

     return
  end subroutine sparse_format_csrdns_z

!!>>> sparse_format_dnscsr: converts a densely stored matrix into a row
!!>>> orientied compactly sparse matrix
  subroutine sparse_format_dnscsr(nrow, ncol, nmax, dns, a, ja, ia)
     implicit none

! external arguments
! row dimension of dense matrix
     integer, intent(in)   :: nrow

! column dimension of dense matrix
     integer, intent(in)   :: ncol

! maximum number of nonzero elements allowed
! this should be set to be the lengths of the arrays a and ja
     integer, intent(in)   :: nmax

! input densely stored matrix
     real(dp), intent(in)  :: dns(nrow,ncol)

! a, ja, ia, output matrix in compressed sparse row format
     integer, intent(out)  :: ia(nrow+1)
     integer, intent(out)  :: ja(nmax)
     real(dp), intent(out) :: a(nmax)

! local variables
! loop index
     integer :: i
     integer :: j
     integer :: k

! init sparse matrix
     a = 0.0_dp
     ia = 0
     ja = 0

     k = 1
     ia(1) = 1
     do i=1,nrow
         do j=1,ncol
             if ( dns(i,j) == 0.0_dp ) CYCLE
             ja(k) = j
             a(k) = dns(i,j)
             k = k + 1
             if ( k > nmax+1 ) then
                 write(mystd,'(a)') 'sparse: error in sparse_format_dnscsr',"k,nmax",k,nmax
                 STOP
             endif ! back if ( k > nmax ) block
         enddo ! over j={1,ncol} loop
         ia(i+1) = k
     enddo ! over i={1,nrow} loop

     return
  end subroutine sparse_format_dnscsr

!!>>> sparse_format_dnscsr_z: converts a densely stored matrix into a row
!!>>> orientied compactly sparse matrix
  subroutine sparse_format_dnscsr_z(nrow, ncol, nmax, dns, sa, ja, ia)
     implicit none

! external arguments
! row dimension of dense matrix
     integer, intent(in)      :: nrow

! column dimension of dense matrix
     integer, intent(in)      :: ncol

! maximum number of nonzero elements allowed
! this should be set to be the lengths of the arrays sa and ja
     integer, intent(in)      :: nmax

! input densely stored matrix
     complex(dp), intent(in)  :: dns(nrow,ncol)

! sa, ja, ia, output matrix in compressed sparse row format
     integer, intent(out)     :: ia(nrow+1)
     integer, intent(out)     :: ja(nmax)
     complex(dp), intent(out) :: sa(nmax)

! local variables
! loop index
     integer :: i
     integer :: j
     integer :: k

! init sparse matrix
     sa = dcmplx(0.0_dp, 0.0_dp)
     ia = 0
     ja = 0

     k = 1
     ia(1) = 1
     do i=1,nrow
         do j=1,ncol
             if ( real( dns(i,j) ) == 0.0_dp .and. aimag( dns(i,j) ) == 0.0_dp ) CYCLE
             ja(k) = j
             sa(k) = dns(i,j)
             k = k + 1
             if ( k > nmax +1) then
                 write(mystd,'(a)') 'sparse: error in sparse_format_dnscsr_z',"k,nmax"
                 write(*,*) "k,nmax=",k,nmax
                 write(*,*) "dns(i,j)=",dns(i,j),i,j
                 STOP
             endif ! back if ( k > nmax ) block
         enddo ! over j={1,ncol} loop
         ia(i+1) = k
     enddo ! over i={1,nrow} loop

     return
  end subroutine sparse_format_dnscsr_z

!!>>> sparse_matrix_getter: this function returns the element a(i,j) of
!!>>> matrix a
  real(dp) &
  function sparse_matrix_getter(i, j, nrow, nmax, a, ja, ia) result(elm)
     implicit none

! external arguments
! the row index of the element sought
     integer, intent(in)  :: i

! the column index of the element sought
     integer, intent(in)  :: j

! row dimension of dense matrix
     integer, intent(in)  :: nrow

! maximum number of nonzero elements allowed
! this should be set to be the lengths of the arrays a and ja
     integer, intent(in)  :: nmax

! a, ja, ia, input matrix in compressed sparse row format
     integer, intent(in)  :: ia(nrow+1)
     integer, intent(in)  :: ja(nmax)
     real(dp), intent(in) :: a(nmax)

! local variables
! loop index
     integer :: k

! memory address of a(i,j)
     integer :: addr

! initialization
     addr = 0
     elm = 0.0_dp

! scan the row - exit as soon as a(i,j) is found
     do k=ia(i),ia(i+1)-1
         if ( ja(k) == j ) then
             addr = k
             EXIT
         endif ! back if ( ja(k) == j ) block
     enddo ! over k={ia(i),ia(i+1)-1} loop

! the required element is contained in sparse matrix
     if ( addr /= 0 ) then
         elm = a(addr)
     endif ! back if ( addr /= 0 ) block

     return
  end function sparse_matrix_getter

!!>>> sparse_matrix_getter_z: this function returns the element sa(i,j) of
!!>>> matrix sa
  complex(dp) &
  function sparse_matrix_getter_z(i, j, nrow, nmax, sa, ja, ia) result(elm)
     implicit none

! external arguments
! the row index of the element sought
     integer, intent(in)     :: i

! the column index of the element sought
     integer, intent(in)     :: j

! row dimension of dense matrix
     integer, intent(in)     :: nrow

! maximum number of nonzero elements allowed
! this should be set to be the lengths of the arrays sa and ja
     integer, intent(in)     :: nmax

! sa, ja, ia, input matrix in compressed sparse row format
     integer, intent(in)     :: ia(nrow+1)
     integer, intent(in)     :: ja(nmax)
     complex(dp), intent(in) :: sa(nmax)

! local variables
! loop index
     integer :: k

! memory address of sa(i,j)
     integer :: addr

! initialization
     addr = 0
     elm = dcmplx(0.0_dp, 0.0_dp)

! scan the row - exit as soon as sa(i,j) is found
     do k=ia(i),ia(i+1)-1
         if ( ja(k) == j ) then
             addr = k
             EXIT
         endif ! back if ( ja(k) == j ) block
     enddo ! over k={ia(i),ia(i+1)-1} loop

! the required element is contained in sparse matrix
     if ( addr /= 0 ) then
         elm = sa(addr)
     endif ! back if ( addr /= 0 ) block

     return
  end function sparse_matrix_getter_z

!!>>> sparse_matmul_amuvec: multiplies a matrix by a vector using the dot
!!>>> product form
  subroutine sparse_matmul_amuvec(nrow, ncol, nmax, a, ja, ia, x, y)
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

! local variables
! loop index
     integer :: i
     integer :: k

! zero out output vector
     y = 0.0_dp

! compute the inner product of row i with vector x
     do i=1,nrow
         do k=ia(i),ia(i+1)-1
             y(i) = y(i) + a(k) * x( ja(k) )
         enddo ! over k={ia(i),ia(i+1)-1} loop
     enddo ! over i={1,nrow} loop

     return
  end subroutine sparse_matmul_amuvec

!!>>> sparse_matmul_amuvec_z: multiplies a matrix by a vector using the
!!>>> dot product form
  subroutine sparse_matmul_amuvec_z(nrow, ncol, nmax, sa, ja, ia, sx, sy)
     implicit none

! external arguments
! row dimension of dense matrix
     integer, intent(in)      :: nrow

! column dimension of dense matrix
     integer, intent(in)      :: ncol

! maximum number of nonzero elements allowed
! this should be set to be the lengths of the arrays sa and ja
     integer, intent(in)      :: nmax

! sa, ja, ia, input matrix in compressed sparse row format
     integer, intent(in)      :: ia(nrow+1)
     integer, intent(in)      :: ja(nmax)
     complex(dp), intent(in)  :: sa(nmax)

! vector, length equal to the column dimension of the dense matrix
     complex(dp), intent(in)  :: sx(ncol)

! vector, complex(dp) array of length nrow, containing the product y = A . x
     complex(dp), intent(out) :: sy(nrow)

! local variables
! loop index
     integer :: i
     integer :: k

! zero out output vector
     sy = dcmplx(0.0_dp, 0.0_dp)

! compute the inner product of row i with vector sx
     do i=1,nrow
         do k=ia(i),ia(i+1)-1
             sy(i) = sy(i) + sa(k) * sx( ja(k) )
         enddo ! over k={ia(i),ia(i+1)-1} loop
     enddo ! over i={1,nrow} loop

     return
  end subroutine sparse_matmul_amuvec_z

!!>>> sparse_matmul_amumat: performs the matrix by matrix product C = A * B
  subroutine sparse_matmul_amumat(nrow, ndim, ncol, nmax, a, ja, ia, b, jb, ib, c, jc, ic)
     implicit none

! external arguments
! the row dimension of matrix A = row dimension of matrix C
     integer, intent(in)   :: nrow

! the column dimension of matrix A = row dimension of matrix B
     integer, intent(in)   :: ndim

! the column dimension of matrix B = column dimension of matrix C
     integer, intent(in)   :: ncol

! the length of the arrays c and jc
! sparse_matmul_amumat() will stop if the result matrix C has a number of
! elements that exceeds nmax
     integer, intent(in)   :: nmax

! a, ja, ia, matrix A in compressed sparse row format
     integer, intent(in)   :: ia(nrow+1)
     integer, intent(in)   :: ja(nmax)
     real(dp), intent(in)  :: a(nmax)

! b, jb, ib, matrix B in compressed sparse row format
     integer, intent(in)   :: ib(ndim+1)
     integer, intent(in)   :: jb(nmax)
     real(dp), intent(in)  :: b(nmax)

! c, jc, ic, resulting matrix C in compressed sparse row format
     integer, intent(out)  :: ic(nrow+1)
     integer, intent(out)  :: jc(nmax)
     real(dp), intent(out) :: c(nmax)

! local variables
! loop index
     integer :: i, j, k

! loop index
     integer :: ka, kb

! dummy integer variables
     integer :: p, q

! integer work array of length equal to the number of columns in matrix B,
! which is an array that has nonzero value if the column index already
! exist, in which case the value is the index of that column
     integer :: iw(ncol)

! dummy real(dp) variables, used to improve the ratio of floating point
! operations to memory accesses
     real(dp) :: atmp, btmp

! init work array
     iw = 0

! init C sparse matrix
     ic(1) = 1

     q = 0
     do i=1,nrow
         do ka=ia(i),ia(i+1)-1
             j = ja(ka)
             atmp = a(ka)
             do kb=ib(j),ib(j+1)-1
                 k = jb(kb)
                 btmp = b(kb)

                 p = iw(k)
                 if ( p == 0 ) then
                     q = q + 1
                     iw(k) = q
                     jc(q) = k
                     c(q) = atmp * btmp
                 else
                     c(p) = c(p) + atmp * btmp
                 endif ! back if ( p == 0 ) block
             enddo ! over kb={ib(j),ib(j+1)-1} loop
         enddo ! over ka={ia(i),ia(i+1)-1} loop

! done this row i, so set work array to zero again
         do k=ic(i),q
             iw( jc( k ) ) = 0
         enddo ! over k={ic(i),q} loop
         ic(i+1) = q + 1
     enddo ! over i={1,nrow} loop

! check the number of nonzero elements
     if ( q > nmax ) then
         write(mystd,'(a)') 'sparse: error in sparse_format_amumat'
         STOP
     endif ! back if ( q > nmax ) block

     return
  end subroutine sparse_matmul_amumat

!!>>> sparse_matmul_amumat_z: performs the matrix by matrix product C = A * B
  subroutine sparse_matmul_amumat_z(nrow, ndim, ncol, nmax, sa, ja, ia, sb, jb, ib, sc, jc, ic)
     implicit none

! external arguments
! the row dimension of matrix A = row dimension of matrix C
     integer, intent(in)      :: nrow

! the column dimension of matrix A = row dimension of matrix B
     integer, intent(in)      :: ndim

! the column dimension of matrix B = column dimension of matrix C
     integer, intent(in)      :: ncol

! the length of the arrays sc and jc
! sparse_matmul_amumat_z() will stop if the result matrix C has a number of
! elements that exceeds nmax
     integer, intent(in)      :: nmax

! sa, ja, ia, matrix A in compressed sparse row format
     integer, intent(in)      :: ia(nrow+1)
     integer, intent(in)      :: ja(nmax)
     complex(dp), intent(in)  :: sa(nmax)

! sb, jb, ib, matrix B in compressed sparse row format
     integer, intent(in)      :: ib(ndim+1)
     integer, intent(in)      :: jb(nmax)
     complex(dp), intent(in)  :: sb(nmax)

! sc, jc, ic, resulting matrix C in compressed sparse row format
     integer, intent(out)     :: ic(nrow+1)
     integer, intent(out)     :: jc(nmax)
     complex(dp), intent(out) :: sc(nmax)

! local variables
! loop index
     integer :: i, j, k

! loop index
     integer :: ka, kb

! dummy integer variables
     integer :: p, q

! integer work array of length equal to the number of columns in matrix B,
! which is an array that has nonzero value if the column index already
! exist, in which case the value is the index of that column
     integer :: iw(ncol)

! dummy complex(dp) variables, used to improve the ratio of floating point
! operations to memory accesses
     complex(dp) :: atmp, btmp

! init work array
     iw = 0

! init C sparse matrix
     ic(1) = 1

     q = 0
     do i=1,nrow
         do ka=ia(i),ia(i+1)-1
             j = ja(ka)
             atmp = sa(ka)
             do kb=ib(j),ib(j+1)-1
                 k = jb(kb)
                 btmp = sb(kb)

                 p = iw(k)
                 if ( p == 0 ) then
                     q = q + 1
                     iw(k) = q
                     jc(q) = k
                     sc(q) = atmp * btmp
                 else
                     sc(p) = sc(p) + atmp * btmp
                 endif ! back if ( p == 0 ) block
             enddo ! over kb={ib(j),ib(j+1)-1} loop
         enddo ! over ka={ia(i),ia(i+1)-1} loop

! done this row i, so set work array to zero again
         do k=ic(i),q
             iw( jc( k ) ) = 0
         enddo ! over k={ic(i),q} loop
         ic(i+1) = q + 1
     enddo ! over i={1,nrow} loop

! check the number of nonzero elements
     if ( q > nmax ) then
         write(mystd,'(a)') 'sparse: error in sparse_format_amumat_z'
         STOP
     endif ! back if ( q > nmax ) block

     return
  end subroutine sparse_matmul_amumat_z

!!>>> sparse_matmul_amudia: performs the matrix by matrix product B = A * Diag
  subroutine sparse_matmul_amudia(nrow, nmax, a, ja, ia, diag, b, jb, ib)
     implicit none

! external arguments
! the row dimension of dense matrix
     integer, intent(in)   :: nrow

! maximum number of nonzero elements allowed
! this should be set to be the lengths of the arrays b and jb
     integer, intent(in)   :: nmax

! a, ja, ia, matrix A in compressed sparse row format
     integer, intent(in)   :: ia(nrow+1)
     integer, intent(in)   :: ja(nmax)
     real(dp), intent(in)  :: a(nmax)

! diagonal matrix stored as a vector diag
     real(dp), intent(in)  :: diag(nrow)

! b, jb, ib, resulting matrix B in compressed sparse row format
     integer, intent(out)  :: ib(nrow+1)
     integer, intent(out)  :: jb(nmax)
     real(dp), intent(out) :: b(nmax)

! local variables
! loop index
     integer :: i
     integer :: k

! loop index
     integer :: k1
     integer :: k2

! init B sparse matrix
     b = 0.0_dp
     ib = 0
     jb = 0

! scale each element
     do i=1,nrow
         k1 = ia(i)
         k2 = ia(i+1) - 1
         do k=k1,k2
             b(k) = a(k) * diag( ja(k) )
         enddo ! over k={k1,k2} loop
     enddo ! over i={1,nrow} loop

     do i=1,nrow+1
         ib(i) = ia(i)
     enddo ! over i={1,nrow+1} loop

     do k=ia(1),ia(nrow+1)-1
         jb(k) = ja(k)
     enddo ! over k={ia(1),ia(nrow+1)-1} loop

     return
  end subroutine sparse_matmul_amudia

!!>>> sparse_matmul_amudia_z: performs the matrix by matrix product B = A * Diag
  subroutine sparse_matmul_amudia_z(nrow, nmax, sa, ja, ia, diag, sb, jb, ib)
     implicit none

! external arguments
! the row dimension of dense matrix
     integer, intent(in)      :: nrow

! maximum number of nonzero elements allowed
! this should be set to be the lengths of the arrays sb and jb
     integer, intent(in)      :: nmax

! sa, ja, ia, matrix A in compressed sparse row format
     integer, intent(in)      :: ia(nrow+1)
     integer, intent(in)      :: ja(nmax)
     complex(dp), intent(in)  :: sa(nmax)

! diagonal matrix stored as a vector diag
     complex(dp), intent(in)  :: diag(nrow)

! sb, jb, ib, resulting matrix B in compressed sparse row format
     integer, intent(out)     :: ib(nrow+1)
     integer, intent(out)     :: jb(nmax)
     complex(dp), intent(out) :: sb(nmax)

! local variables
! loop index
     integer :: i
     integer :: k

! loop index
     integer :: k1
     integer :: k2

! init B sparse matrix
     sb = dcmplx(0.0_dp, 0.0_dp)
     ib = 0
     jb = 0

! scale each element
     do i=1,nrow
         k1 = ia(i)
         k2 = ia(i+1) - 1
         do k=k1,k2
             sb(k) = sa(k) * diag( ja(k) )
         enddo ! over k={k1,k2} loop
     enddo ! over i={1,nrow} loop

     do i=1,nrow+1
         ib(i) = ia(i)
     enddo ! over i={1,nrow+1} loop

     do k=ia(1),ia(nrow+1)-1
         jb(k) = ja(k)
     enddo ! over k={ia(1),ia(nrow+1)-1} loop

     return
  end subroutine sparse_matmul_amudia_z

!!>>> sparse_matmul_diamua: performs the matrix by matrix product B = Diag * A
  subroutine sparse_matmul_diamua(nrow, nmax, diag, a, ja, ia, b, jb, ib)
     implicit none

! external arguments
! the row dimension of dense matrix
     integer, intent(in)   :: nrow

! maximum number of nonzero elements allowed
! this should be set to be the lengths of the arrays b and jb
     integer, intent(in)   :: nmax

! diagonal matrix stored as a vector diag
     real(dp), intent(in)  :: diag(nrow)

! a, ja, ia, matrix A in compressed sparse row format
     integer, intent(in)   :: ia(nrow+1)
     integer, intent(in)   :: ja(nmax)
     real(dp), intent(in)  :: a(nmax)

! b, jb, ib, resulting matrix B in compressed sparse row format
     integer, intent(out)  :: ib(nrow+1)
     integer, intent(out)  :: jb(nmax)
     real(dp), intent(out) :: b(nmax)

! local variables
! loop index
     integer :: i
     integer :: k

! loop index
     integer :: k1
     integer :: k2

! init B sparse matrix
     b = 0.0_dp
     ib = 0
     jb = 0

! normalize each row
     do i=1,nrow
         k1 = ia(i)
         k2 = ia(i+1) - 1
         do k=k1,k2
             b(k) = a(k) * diag(i)
         enddo ! over k={k1,k2} loop
     enddo ! over i={1,nrow} loop

     do i=1,nrow+1
         ib(i) = ia(i)
     enddo ! over i={1,nrow+1} loop

     do k=ia(1),ia(nrow+1)-1
         jb(k) = ja(k)
     enddo ! over k={ia(1),ia(nrow+1)-1} loop

     return
  end subroutine sparse_matmul_diamua

!!>>> sparse_matmul_diamua_z: performs the matrix by matrix product B = Diag * A
  subroutine sparse_matmul_diamua_z(nrow, nmax, diag, sa, ja, ia, sb, jb, ib)
     implicit none

! external arguments
! the row dimension of dense matrix
     integer, intent(in)      :: nrow

! maximum number of nonzero elements allowed
! this should be set to be the lengths of the arrays sb and jb
     integer, intent(in)      :: nmax

! diagonal matrix stored as a vector diag
     complex(dp), intent(in)  :: diag(nrow)

! sa, ja, ia, matrix A in compressed sparse row format
     integer, intent(in)      :: ia(nrow+1)
     integer, intent(in)      :: ja(nmax)
     complex(dp), intent(in)  :: sa(nmax)

! sb, jb, ib, resulting matrix B in compressed sparse row format
     integer, intent(out)     :: ib(nrow+1)
     integer, intent(out)     :: jb(nmax)
     complex(dp), intent(out) :: sb(nmax)

! local variables
! loop index
     integer :: i
     integer :: k

! loop index
     integer :: k1
     integer :: k2

! init B sparse matrix
     sb = dcmplx(0.0_dp, 0.0_dp)
     ib = 0
     jb = 0

! normalize each row
     do i=1,nrow
         k1 = ia(i)
         k2 = ia(i+1) - 1
         do k=k1,k2
             sb(k) = sa(k) * diag(i)
         enddo ! over k={k1,k2} loop
     enddo ! over i={1,nrow} loop

     do i=1,nrow+1
         ib(i) = ia(i)
     enddo ! over i={1,nrow+1} loop

     do k=ia(1),ia(nrow+1)-1
         jb(k) = ja(k)
     enddo ! over k={ia(1),ia(nrow+1)-1} loop

     return
  end subroutine sparse_matmul_diamua_z


!# if defined (MKL_lib)
  include "sparse_mkl_wraper.f90"
!# endif /* MKL_lib */



  end module m_sparse

