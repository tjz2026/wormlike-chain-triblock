       ! this is a wraper program  specifically written fortran module for matrix inverse 
      ! by calling the lapack lib ,developed by Jiuzhou Tang

         module matrix_inverse

         
         contains
           

         subroutine complex_matrix_inverse(a_R,a_I,m)
          integer :: m
          REAL*8 :: a_R(m,m)
          REAL*8 :: a_I(m,m)
          double complex :: zmat(m,m)
          integer :: i,j
           do j=1,m
             do i=1,m
               zmat(i,j)=dcmplx(a_R(i,j),a_I(i,j))
             enddo
            enddo

            
           call zmat_inv(m, zmat)

           
           do j=1,m
             do i=1,m
               a_R(i,j)=real(zmat(i,j))
               a_I(i,j)=aimag(zmat(i,j))
             enddo
            enddo

          end subroutine complex_matrix_inverse 
         



            !>>> invert real(dp) matrix using lapack subroutines
  subroutine dmat_inv(ndim, dmat)

     implicit none

! external arguments
! dimension of dmat matrix
     integer, intent(in) :: ndim

! object matrix, on entry, it contains the original matrix, on exit,
! it is destroyed and replaced with the inversed matrix
     real*8, intent(inout) :: dmat(ndim,ndim)

! local variables
! error flag
     integer  :: ierror

! working arrays for lapack subroutines
     integer  :: ipiv(ndim)
     real*8 :: work(ndim)

! computes the LU factorization of a general m-by-n matrix, need lapack
! package, dgetrf subroutine


     write(*,*) "ndim=",ndim 

     call dgetrf(ndim, ndim, dmat, ndim, ipiv, ierror)
     if ( ierror /= 0 ) then
         !call ctqmc_print_error('ctqmc_dmat_inv','error in lapack subroutine dgetrf')
         write(*,*)"dmat_inv','error in lapack subroutine dgetrf=",ierror
         stop
     endif

! computes the inverse of an LU-factored general matrix, need lapack
! package, dgetri subroutine
     call dgetri(ndim, dmat, ndim, ipiv, work, ndim, ierror)
     if ( ierror /= 0 ) then
         !call ctqmc_print_error('ctqmc_dmat_inv','error in lapack subroutine dgetri')
         write(*,*)"dmat_inv','error in lapack subroutine dgetri",ierror
         stop
     endif

     return
  end subroutine dmat_inv

!>>> invert complex(dp) matrix using lapack subroutines
  subroutine zmat_inv(ndim, zmat)

     implicit none

! external arguments
! dimension of zmat matrix
     integer, intent(in) :: ndim

! object matrix, on entry, it contains the original matrix, on exit,
! it is destroyed and replaced with the inversed matrix
     double complex, intent(inout) :: zmat(ndim,ndim)

! local variables
! error flag
     integer     :: ierror

! working arrays for lapack subroutines
     integer     :: ipiv(ndim)
     double complex :: work(ndim)

! computes the LU factorization of a general m-by-n matrix, need lapack
! package, zgetrf subroutine
     call zgetrf(ndim, ndim, zmat, ndim, ipiv, ierror)
     if ( ierror /= 0 ) then
         !call ctqmc_print_error('ctqmc_zmat_inv','error in lapack subroutine zgetrf')
         write(*,*)"zmat_inv','error in lapack subroutine zgetrf",ierror
         stop
     endif

! computes the inverse of an LU-factored general matrix, need lapack
! package, zgetri subroutine
     call zgetri(ndim, zmat, ndim, ipiv, work, ndim, ierror)
     if ( ierror /= 0 ) then
         !call ctqmc_print_error('ctqmc_zmat_inv','error in lapack subroutine zgetri')
         write(*,*)"zmat_inv','error in lapack subroutine zgetri",ierror
         stop
     endif

     return
  end subroutine zmat_inv




          subroutine lapacksolver(A,n,np,b)
            implicit none

! the working precision of real variables, define to give double precision
!integer, parameter :: wp = kind(1.d0)

! system matrix A and inhomgenous term b
!real(kind=8), allocatable, dimension(:,:) :: A
!real(kind=8), allocatable, dimension(:) :: b

          integer :: sys_order,n,np

          real*8 A(n,n),b(n)

!---------------------------------------------------------------------------
 
! read input and set number of processors
!call get_parameters(sys_order)
             sys_order=n
 

! initialize and solve system, return timing results
          !call init_system(A,b,sys_order)
          call solver(A,b,sys_order)

! output timing results
!call put_results(sys_order,wall_time)
    
! deallocate space
!deallocate(A)
!deallocate(b)
       
            return
             end subroutine lapacksolver


!===========================================================================

            subroutine solver(A,b,sys_order)

! argument declarations
           integer :: sys_order, nrhs, info,i
             real(kind=8) A(sys_order,sys_order),b(sys_order)
!real(kind=8), dimension(:), intent(inout) :: b
           integer, dimension(:), allocatable :: ipvt
           nrhs=1

! allocate space for pivots
            allocate(ipvt(sys_order))

!  Call LU factorization routine
!wtime(1) = rtc()
          call dgetrf(sys_order,sys_order,A,sys_order,ipvt,info)
!wtime(2) = rtc()

           if (info .ne. 0) print*,'dgetrf info: ',info

!  Call LU solver routine
       call dgetrs('N',sys_order,nrhs,A,sys_order,ipvt,b,sys_order,info)

        if (info .ne. 0) print*,'dgetrs info: ',info
           ! write(*,*) "x="
        ! do i=1,sys_order
         !  write(*,*) b(i)
         ! enddo

!wall_time = wtime(2) - wtime(1)

        deallocate(ipvt)

            end subroutine solver

            subroutine solver_z(A,b,sys_order)
            use mpi
! argument declarations
           integer :: sys_order, nrhs, info,i
           !complex(kind=16) A(sys_order,sys_order),b(sys_order)
           double complex :: A(sys_order,sys_order),b(sys_order)
!real(kind=8), dimension(:), intent(inout) :: b
           integer, dimension(:), allocatable :: ipvt
           real*8 :: t2,t1   
           nrhs=1

! allocate space for pivots
            allocate(ipvt(sys_order))

!  Call LU factorization routine
!wtime(1) = rtc()
       !   t1=MPI_Wtime()
          call zgetrf(sys_order,sys_order,A,sys_order,ipvt,info)
       !   t2=MPI_Wtime()
       !   write(*,*) "total LU time",t2-t1
!wtime(2) = rtc()

           if (info .ne. 0) print*,'dgetrf info: ',info

!  Call LU solver routine
        !  t1=MPI_Wtime()
       call zgetrs('N',sys_order,nrhs,A,sys_order,ipvt,b,sys_order,info)
        !  t2=MPI_Wtime()
        !  write(*,*) "total solve time",t2-t1
        if (info .ne. 0) print*,'dgetrs info: ',info
           ! write(*,*) "x="
        ! do i=1,sys_order
         !  write(*,*) b(i)
         ! enddo

!wall_time = wtime(2) - wtime(1)

        deallocate(ipvt)

            end subroutine solver_z





          
 end module matrix_inverse




