!-------------------------------------------------------------------------
!project : wormlike chain SCFT
!subroutines :
! source  : sparse_MKL_mod.f90
! type    : subroutines
! author  : Jiu-zhou Tang (tangjiuzhou@iccas.ac.cn)
! history : 08/20/2013 by Jiu-zhou Tang
! purpose : to provide utility subroutines for MKL sparse blas routines
! input   :
! output  :
! status  : unstable
! comment : 
         !this module is specifically written for calling MKL sparse blas lib.
         !Developed by Jiuzhou Tang,this may not be a good wraper,but it's better than nothing,right?
         !So,just deal with it or you can write your own damn code,like I did.Enjoy!
         !PS: Well,actually,MKL does have a routine for format converting,(see its example),but it is
        ! not suitable for symmetric matrix storage format,so,I wrote my own routine.You are welcome!

        module sparse_MKL

       !convert a sparse matrix to csr3 format for double and double complex
       !note that csr3 is the three array format for CSR storage type,what?You don't 
       !know the CSR format shit? Pick up the MKL's manual and dig in!
       

        contains

       !convert a general sparse matrix to the compressed csr format.
       !note that normally,you need to predefine a integer for nmax which 
       !is the maximum size of the arrays,thanks to my effort,we have another
       !subroutine to get the size of our storaged array before we call the following subroutine.
       !Awesome,isn't it?
         

      !question mark:if a certain row has no nonzero element,how to deal with it?
      !think about it,and come with a solution.
      !in this subroutine ,we will ignore it ,simply because it won't happen in our case.
       subroutine get_sparse_size_one_based_dp(a,m,n,nonzero_size)
        implicit none
        integer,intent(in) :: m,n
        integer,intent(inout) :: nonzero_size
        real*8,intent(in) :: a(m,n)
        real*8,parameter :: THRESH_matrix=1.0e-18
        integer :: i,j,k

         k=0
         do j=1,n
           do i=1,m
             if(abs(a(i,j))>=THRESH_matrix) then
             k=k+1
             endif
           enddo
          enddo
          nonzero_size=k

         end subroutine get_sparse_size_one_based_dp

       subroutine get_sparse_size_one_based_cmplx(a,m,n,nonzero_size)
        implicit none
        integer,intent(in) :: m,n
        integer,intent(inout) :: nonzero_size
        double complex,intent(in) :: a(m,n)
        real*8,parameter :: THRESH_matrix=1.0e-18
        integer :: i,j,k

         k=0
         do j=1,n
           do i=1,m
             if(abs(a(i,j))>=THRESH_matrix) then
             k=k+1
             endif
           enddo
          enddo
          nonzero_size=k

         end subroutine get_sparse_size_one_based_cmplx

        subroutine convert_sparse_one_based_csr3_dp(a,m,n,nonzero_size,val,col_ind,row_ptr)
        implicit none         
        integer,intent(in) :: m,n,nonzero_size
        real*8,intent(in) :: a(m,n)
        real*8,intent(inout) :: val(1:nonzero_size)
        integer,intent(inout) :: col_ind(1:nonzero_size)
        integer,intent(inout) :: row_ptr(1:m+1)
        integer :: nonzero_row(1:m+1)
        integer :: i,j,k,ii
        real*8,parameter :: THRESH_matrix=1.0e-18
        k=0
        nonzero_row=0
        do i=1,m
             ii=0
          do j=1,n
             if(abs(a(i,j))>=THRESH_matrix) then
             k=k+1

               if(k>nonzero_size) then
               write(*,*) "k exceed the nonzero_size,check out!"
               stop
               endif

             val(k)=a(i,j)
             col_ind(k)=j
             ii=ii+1
             endif
             nonzero_row(i)=ii
          enddo
       enddo

         

         row_ptr(1)=1

       do i=2,m
          row_ptr(i)=nonzero_row(i-1)+row_ptr(i-1)
       enddo
          row_ptr(m+1)=nonzero_size+1

       end subroutine convert_sparse_one_based_csr3_dp
        
        subroutine convert_one_based_csr3_dp_to_csrNIST(m,nonzero_size, & 
                           row_ptr,pointerB,pointerE)
         implicit none
         integer :: m,nonzero_size,i
         integer,intent(in) :: row_ptr(1:m+1)
         integer,intent(inout) :: pointerB(1:m)
         integer,intent(inout) :: pointerE(1:m)
         do i=1,m
           pointerB(i)=row_ptr(i)
           pointerE(i)=row_ptr(i+1)
         enddo
         end subroutine convert_one_based_csr3_dp_to_csrNIST

        !this subroutine convert a symmetric matrix "A" to a lower triangle matrix "lower",note that A is
        !a symmetric matirix and we store only the lower triangle part of the matirx
         subroutine convert_symmetric_matrix_dp(a,m,n,lower)
         implicit none
        integer,intent(in) :: m,n
        real*8,intent(in) :: a(m,n)
        real*8,intent(inout) :: lower(m,n)
        integer :: i,j,k,ii
        real*8,parameter :: THRESH_matrix=1.0e-18
         !first,check out whether A is symmetric or not
          if(m/=n) then
          write(*,*) "m /=n,not a symmetrix matrix !!!"
          stop
          endif

          lower=0.0d0
          do i=1,m
           do j=1,n
             if(abs(a(i,j)-a(j,i))>1.0e-10) then
             write(*,*) "a(i,j)/=a(j,i)","i,j:",i,j
             stop
             endif
             if(i>=j) then
             lower(i,j)=a(i,j) 
             endif
           enddo
           enddo

          end subroutine convert_symmetric_matrix_dp

         subroutine get_symmetric_matrix_size_complex_dp(a,m,n,N_low)
         implicit none
         integer,intent(in) :: m,n
         integer,intent(inout) :: N_low
         double complex,intent(in) :: a(m,n)
         integer :: i,j,k,ii
         real*8,parameter :: THRESH_matrix=1.0e-18
         real*8 :: aa
         !first,check out whether A is symmetric or not
          if(m/=n) then
          write(*,*) "m /=n,not a symmetrix matrix !!!"
          stop
          endif


          k=0
          do i=1,m
           do j=1,n
            aa=max(abs(real(a(i,j)-a(j,i))),abs(aimag(a(i,j)-a(j,i))))
           if(aa >1.0e-10) then
         write(*,*) "err:a(i,j),a(j,i),",a(i,j),a(j,i),"i,j:",i,j
             stop
             endif
             aa=max(abs(real(a(i,j))) , abs(aimag(a(i,j))))
             if(i>=j .and. aa >=THRESH_matrix) then
             k=k+1
             endif
           enddo
           enddo
           N_low=k

          end subroutine get_symmetric_matrix_size_complex_dp

         subroutine get_symmetric_matrix_size_dp(a,m,n,N_low)
         implicit none
        integer,intent(in) :: m,n
        integer,intent(inout) :: N_low
        real*8,intent(in) :: a(m,n)
        integer :: i,j,k,ii
        real*8,parameter :: THRESH_matrix=1.0e-18
         !first,check out whether A is symmetric or not
          if(m/=n) then
          write(*,*) "m /=n,not a symmetrix matrix !!!"
          stop
          endif


          k=0
          do i=1,m
           do j=1,n
             if(abs(a(i,j)-a(j,i))>1.0e-10) then
         write(*,*) "err:a(i,j),a(j,i),",a(i,j),a(j,i),"i,j:",i,j
             stop
             endif
             if(i>=j .and. abs(a(i,j))>=THRESH_matrix) then
             k=k+1
             endif
           enddo
           enddo
           N_low=k

          end subroutine get_symmetric_matrix_size_dp

   
        !to save the memeory,we can store the matrix to an array
         subroutine convert_symmetric_matrix_to_array_dp(a,m,n,N_low, &
         sy_values,sy_col,sy_rowIndex)
         implicit none
        integer,intent(in) :: m,n,N_low
        real*8,intent(in) :: a(m,n)
        real*8,intent(inout) :: sy_values(N_low)
        integer,intent(inout) :: sy_col(N_low)
        integer,intent(inout) :: sy_rowIndex(m+1)
        integer :: i,j,k,ii
        integer :: nonzero_row(1:m+1)
        real*8,parameter :: THRESH_matrix=1.0e-18
         !first,check out whether A is symmetric or not
          if(m/=n) then
          write(*,*) "m /=n,not a symmetrix matrix !!!"
          stop
          endif


          sy_values=0.0d0
          nonzero_row=0
          k=0
          do i=1,m
             ii=0
           do j=1,n
             if(abs(a(i,j)-a(j,i))>1.0e-10) then
         write(*,*) "err:a(i,j),a(j,i),",a(i,j),a(j,i),"i,j:",i,j
             stop
             endif
             if(i>=j .and. abs(a(i,j))>=THRESH_matrix) then
             k=k+1
             sy_values(k)=a(i,j) 
             sy_col(k)=j
             ii=ii+1
             endif
             
           enddo
            nonzero_row(i)=ii
           enddo
         if(k/=N_low) then
     write(*,*) "warning,index wrong of lower_array","k=",k,"N_low=",N_low
     !     stop
          endif

          sy_rowIndex(1)=1

       do i=2,m
          sy_rowIndex(i)=nonzero_row(i-1)+sy_rowIndex(i-1)
       enddo
          sy_rowIndex(m+1)=N_low+1



          end subroutine convert_symmetric_matrix_to_array_dp
            
         subroutine symmetric_CSR_matvec(vec,q,m,n,N_low, &
         sy_values,sy_col,sy_rowIndex)
         implicit none
        integer,intent(in) :: m,n,N_low
        real*8,intent(in) :: vec(n)
        real*8,intent(inout) :: sy_values(N_low),q(n)
        integer,intent(inout) :: sy_col(N_low)
        integer,intent(inout) :: sy_rowIndex(m+1)
        integer :: i,j,k,ii,N_low_true
        integer :: nonzero_row(1:m+1)
        real*8,parameter :: THRESH_matrix=1.0e-18
        real*8 :: aa
         !first,check out whether A is symmetric or not
          WRITE(*,*) "m=,n=",m,n
          if(m/=n) then
          write(*,*) "m /=n,not a symmetrix matrix !!!"
          stop
          endif

             N_low_true=N_low  
         do i=1,N_low
             !aa=abs(sy_values(i))
             if(abs(sy_values(i))<1.0E-18) then
             !if(abs(sy_values(i))<1.0E-18) then
             N_low_true=i-1
             exit
             endif
         enddo
         write(*,*) "true N_low=",N_low_true,sy_values(N_low_true)

          q=0.0d0
          do k=1,N_low_true
                j=sy_col(k)
                
              do ii=1,m
            if(sy_rowIndex(ii)<=k .and.  &
             sy_rowIndex(ii+1)>k) then
                 i=ii
                 exit
             endif
              enddo
              if(i<j) then
              write(*,*) "error detected",i,j
                else if(i>j) then
                q(i)=q(i)+sy_values(k)*vec(j)
                q(j)=q(j)+sy_values(k)*vec(i)
                else
                q(i)=q(i)+sy_values(k)*vec(j)
              endif
                
          enddo
          end subroutine symmetric_CSR_matvec
          
         subroutine get_symmetric_matrix_size_cmplx(a,m,n,N_low)
         implicit none
        integer,intent(in) :: m,n
        integer,intent(inout) :: N_low
        double complex,intent(in) :: a(m,n)
        integer :: i,j,k,ii
        real*8,parameter :: THRESH_matrix=1.0e-18
         !first,check out whether A is symmetric or not
          if(m/=n) then
          write(*,*) "m /=n,not a symmetrix matrix !!!"
          stop
          endif


          k=0
          do i=1,m
           do j=1,n
            if(max( abs(REAL(a(i,j))-REAL(a(j,i))), &
            abs(AIMAG(a(i,j))-AIMAG(a(j,i))) )>1.0e-10) then
             write(*,*) "a(i,j)/=a(j,i)","i,j:",i,j
             stop
             endif
             if(i>=j .and. abs(a(i,j))>=THRESH_matrix) then
             k=k+1
             endif
           enddo
           enddo
           N_low=k

          end subroutine get_symmetric_matrix_size_cmplx

         subroutine convert_symmetric_matrix_to_array_cmplx(a,m,n,N_low, &
         sy_values,sy_col,sy_rowIndex)
         implicit none
        integer,intent(in) :: m,n,N_low
        double complex,intent(in) :: a(m,n)
        double complex,intent(inout) :: sy_values(N_low)
        integer,intent(inout) :: sy_col(N_low)
        integer,intent(inout) :: sy_rowIndex(m+1)
        integer :: i,j,k,ii
        integer :: nonzero_row(1:m+1)
        real*8,parameter :: THRESH_matrix=1.0e-18
         !first,check out whether A is symmetric or not
          if(m/=n) then
          write(*,*) "m /=n,not a symmetrix matrix !!!"
          stop
          endif


          sy_values=dcmplx(0.0d0,0.0d0)

          nonzero_row=0
          k=0
          do i=1,m
             ii=0
           do j=1,n
            if(max( abs(REAL(a(i,j))-REAL(a(j,i))), &
            abs(AIMAG(a(i,j))-AIMAG(a(j,i))) )>1.0e-10) then
             write(*,*) "a(i,j)/=a(j,i)","i,j:",i,j
             stop
             endif
             if(i>=j .and. abs(a(i,j))>=THRESH_matrix) then
             k=k+1
             sy_values(k)=a(i,j) 
             sy_col(k)=j
             ii=ii+1
             endif
             
           enddo
            nonzero_row(i)=ii
           enddo
         if(k/=N_low) then
     write(*,*) "warning!,index wrong of lower_array","k=",k,"N_low=",N_low
          !stop
          endif

          sy_rowIndex(1)=1

       do i=2,m
          sy_rowIndex(i)=nonzero_row(i-1)+sy_rowIndex(i-1)
       enddo
          sy_rowIndex(m+1)=N_low+1



          end subroutine convert_symmetric_matrix_to_array_cmplx
            

   
      subroutine convert_sparse_zero_based_csr3_dp(a,m,n,nonzero_size, &
                                 val,col_ind,row_ptr)
        implicit none         
        integer,intent(in) :: m,n,nonzero_size
        real*8,intent(in) :: a(0:m-1,0:n-1)
        real*8,intent(inout) :: val(0:nonzero_size-1)
        integer,intent(inout) :: col_ind(0:nonzero_size-1)
        integer,intent(inout) :: row_ptr(0:m)
        integer :: nonzero_row(0:m)
        integer :: i,j,k,ii
        real*8,parameter :: THRESH_matrix=1.0e-18
        k=-1
        nonzero_row=0
        do i=0,m-1
             ii=0
          do j=0,n-1
             if(abs(a(i,j))>=THRESH_matrix) then
             k=k+1

               if(k>nonzero_size-1) then
               write(*,*) "k exceed the nonzero_size,check out!"
               stop
               endif

             val(k)=a(i,j)
             col_ind(k)=j
             ii=ii+1
             endif
             nonzero_row(i)=ii
          enddo
       enddo

         

         row_ptr(0)=0

       do i=1,m-1
          row_ptr(i)=nonzero_row(i-1)+row_ptr(i-1)
       enddo
          row_ptr(m)=nonzero_size

       end subroutine convert_sparse_zero_based_csr3_dp
        
        subroutine convert_zero_based_csr3_dp_to_csrNIST(m,nonzero_size, & 
                         row_ptr,pointerB,pointerE)
         implicit none
         integer :: m,nonzero_size,i
         integer,intent(in) :: row_ptr(0:m)
         integer,intent(inout) :: pointerB(0:m-1)
         integer,intent(inout) :: pointerE(0:m-1)
         do i=0,m-1
           pointerB(i)=row_ptr(i)
           pointerE(i)=row_ptr(i+1)
         enddo
         end subroutine convert_zero_based_csr3_dp_to_csrNIST
 
        subroutine get_sparse_size_zero_based(a,m,n,nonzero_size)
        implicit none         
        integer,intent(in) :: m,n
        integer,intent(inout) :: nonzero_size
        real*8,intent(in) :: a(0:m-1,0:n-1)
        real*8,parameter :: THRESH_matrix=1.0e-18
        integer :: i,j,k

         k=0
         do j=0,n-1
           do i=0,m-1
             if(abs(a(i,j))>=THRESH_matrix) then
             k=k+1
             endif
           enddo
          enddo
          nonzero_size=k  
                      
         end subroutine get_sparse_size_zero_based









end module sparse_MKL

















