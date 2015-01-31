        ! Consider the matrix A (see Appendix 'Sparse Storage Formats for Sparse Blas
! level 2-3')
!
!                 |   1       -1      0   -3     0   |
!                 |  -2        5      0    0     0   |
!   A    =        |   0        0      4    6     4   |,
!                 |  -4        0      2    7     0   |
!                 |   0        8      0    0    -5   |
!         values = (1 -1 -3 -2 5 4 6 4 -4 2 7 8 -5)
!         columns = (1 2 4 1 2 3 4 5 1 3 4 2 5)
!         rowIndex = (1  4  6  9  12 14)
!let the arrays values1, rowIndex1, columns1 which are the following
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
!**************************************************************************
!**************for complex matrix******************************************
!
!                 |   1,1     -1,1    0   -3,1   0   |
!                 |  -2,1      5,1    0    0     0   |
!   A    =        |   0        0      4,1  6,1   4,1 |,
!                 |  -4,1      0      2,1  7,1   0   |
!                 |   0        8,1    0    0    -5,1 |
!  
! decomposed as 
!                      A = L + D + U,      
!
!  where L is the strict  lower triangle of A, U is the strictly  upper triangle
!  of A, D is the main diagonal. Namely 
!        |   0    0    0    0     0   |       |  0   -1,1  0   -3,1  0   | 
!        |  -2,1  0    0    0     0   |       |  0    0    0    0    0   |
!   L  = |   0    0    0    0     0   |,  U=  |  0    0    0    6,1  4,1 |
!        |  -4,1  0    2,1  0     0   |       |  0    0    0    0    0   |
!        |   0    8,1  0    0     0   |       |  0    0    0    0    0   |
!
!           |   1,1  0    0    0    0   |
!           |   0    5,1  0    0    0   |
!   D    =  |   0    0    4,1  0    0   |.
!           |   0    0    0    7,1  0   |
!           |   0    0    0    0   -5   |
!  The matrix A is represented in the compressed sparse row storage scheme with the help of three 
!  arrays  (see Appendix 'Sparse Matrix Storage') as follows:
!         values = (1,1 -1,1 -3,1 -2,1 5,1 4,1 6,1 4,1 -4,1 2,1 7,1 8,1 -5,1) 
!         columns = (1 2 4 1 2 3 4 5 1 3 4 2 5) 
!         rowIndex = (1  4  6  9  12 14)
!*****************************************************************************************************
        program main
        use sparse_MKL 
        implicit none
        !real*8,allocatable :: A(:,:)
        REAL*8 :: A(5,5) 
        integer :: m,n,i,j,k
        integer :: nonzero_size,N_low
        real*8,allocatable :: values(:),sy_values(:)
        integer,allocatable :: columns(:),rowIndex(:),sy_columns(:),sy_rowIndex(:)
        real*8 rhs(5, 2), sol(5, 2), temp(5, 2)
        data sol/1.D0, 1.D0, 1.D0, 1.D0, 1.D0, &
          5.D0, 4.D0, 3.D0, 2.D0, 1.D0/
        DATA A /0.0,-2,0,-4,0,-1,5,0,0,8,0,0,4,2,0,-3,0,6,7,0,0,0,4,0,-5 /
        double complex :: A_z(5,5)
        DATA A_z /(1.0,1.0),(-2.0,1.0),(0.0,0.0),(-4.0,1.0),(0.0,0.0), &
        (-1.0,1.0),(5.0,1.0),(0.0,0.0),(0.0,0.0),(8.0,1.0),(0.0,0.0), &
        (0.0,0.0),(4.0,1.0),(2.0,1.0),(0.0,0.0),(-3.0,1.0),(0.0,0.0), &
        (6.0,1.0),(7.0,1.0),(0.0,0.0),(0.0,0.0),(0.0,0.0),(4.0,1.0), &
        (0.0,0.0),(-5.0,1.0) /  
        double complex,allocatable :: sy_values_z(:),values_z(:)
        integer,allocatable :: columns_z(:),rowIndex_z(:),sy_columns_z(:), &
        sy_rowIndex_z(:)
        double complex :: rhs_z(5, 2), sol_z(5, 2), temp_z(5, 2) 
        data sol_z /(1.0d0, 1.0d0), (1.0d0, 1.0d0), (1.0d0, 1.0d0), &
           (1.0d0, 1.0d0), (1.0d0, 1.0d0), (5.0d0, 1.0d0), &
           (4.0d0, 1.0d0), (3.0d0, 1.0d0), (2.0d0, 1.0d0), &
           (1.0d0, 1.0d0)/
       integer job(8)
      integer info

        m=5
        n=5
       
        call get_sparse_size_one_based_dp(a,m,n,nonzero_size)
        write(*,*) "nonzero_size=",nonzero_size
        allocate(values(1:nonzero_size))
        allocate(columns(1:nonzero_size))
        allocate(rowIndex(1:m+1))
        call convert_sparse_one_based_csr3_dp(a,m,n,nonzero_size,values,columns,rowIndex) 
        do i=1,nonzero_size
        write(*,*) "values",values(i),"columns",columns(i)
        enddo
        
        do i=1,5
         do j=1,5
             if(j>i) then
             A(i,j)=A(j,i)
             endif
         enddo
        enddo  

        call get_symmetric_matrix_size_dp(a,m,n,N_low)

        allocate(sy_values(1:N_low))
        allocate(sy_columns(1:N_low))
        allocate(sy_rowIndex(1:m+1))
        

        call convert_symmetric_matrix_to_array_dp(a,m,n,N_low,sy_values,  &
        sy_columns,sy_rowIndex)
        do i=1,N_low
        write(*,*) "sy_values",sy_values(i),"sy_columns",sy_columns(i)
        enddo
        do i=1,m+1
        write(*,*) "sy_rowIndex",sy_rowIndex(i)
        enddo

        do i=1,5
         do j=1,5
             if(j>i) then
             A(i,j)=0.0
             endif
         enddo
        enddo  
        job(1)=0
        job(2)=1
        job(3)=1
        job(4)=0
        job(5)=N_low
        job(6)=1


       call mkl_ddnscsr(job,m,m,A,m,sy_values,sy_columns, &
        sy_rowIndex,info)
       write(*,*) "SYMMETRIC MATRIX-VECTOR CSR converted from MKL "
        do i=1,N_low
        write(*,*) "sy_values",sy_values(i),"sy_columns",sy_columns(i)
        enddo
        do i=1,m+1
        write(*,*) "sy_rowIndex",sy_rowIndex(i)
        enddo


          print*, '     SYMMETRIC MATRIX-VECTOR MULTIPLY ROUTINE '
          print*, '     INPUT DATA FOR MKL_DCSRSYMV '

           call mkl_dcsrsymv('l', m,  values, rowIndex, columns, &
                sol, rhs)
           call mkl_dcsrsymv('l', m,  sy_values, sy_rowIndex, sy_columns, &
                sol, temp)

          print*
          print*, '     SYMMETRIC MATRIX-VECTOR MULTIPLY '
          print*, '     MKL_DCSRSYMV CALLED TWO TIMES '
          do i=1, m
              write(*,*)  rhs(i,1), temp(i,1)
          enddo
     
!**************now test the complex matrix sparse blas routine***************

        do i=1,5
         do j=1,5
             if(j>i) then
             A_z(i,j)=A_z(j,i)
             endif
         enddo
        enddo  

        call get_symmetric_matrix_size_cmplx(A_z,m,n,N_low)

        allocate(sy_values_z(1:N_low))
        allocate(sy_columns_z(1:N_low))
        allocate(sy_rowIndex_z(1:m+1))

        call convert_symmetric_matrix_to_array_cmplx(A_z,m,n,N_low,sy_values_z,  &
        sy_columns_z,sy_rowIndex_z)
        do i=1,N_low
        write(*,*) "sy_values",sy_values_z(i),"sy_columns",sy_columns_z(i)
        enddo
        do i=1,m+1
        write(*,*) "sy_rowIndex",sy_rowIndex_z(i)
        enddo

        do i=1,5
         do j=1,5
             if(j>i) then
             A_z(i,j)=dcmplx(0.0,0.0)
             endif
         enddo
        enddo  
        job(1)=0
        job(2)=1
        job(3)=1
        job(4)=0
        job(5)=N_low
        job(6)=1


       call mkl_zdnscsr(job,m,m,A_z,m,sy_values_z,sy_columns_z, &
        sy_rowIndex_z,info)
      write(*,*) "COMPLEX SYMMETRIC MATRIX CSR converted from MKL "
        do i=1,N_low
      write(*,*) "sy_values",sy_values_z(i),"sy_columns",sy_columns_z(i)
        enddo
        do i=1,m+1
        write(*,*) "sy_rowIndex",sy_rowIndex_z(i)
        enddo

      call mkl_zcsrsymv('l', m,  sy_values_z, sy_rowIndex_z, sy_columns_z, &
                sol_z, temp_z)

          print*, '     SYMMETRIC COMPLEX MATRIX-VECTOR MULTIPLY '
          do i=1, m
              write(*,*)   temp_z(i,1)
          enddo
        stop
        end
        
        
        
