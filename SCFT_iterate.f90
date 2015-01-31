      Function Anderson_index(index_i,n_r)
      implicit none
      integer,intent(in) :: index_i,n_r
      integer  :: Anderson_index

!     // this function is designed for locating the index of the preceding steps e.g. k-m  or k-n in the notes
!        // if n_r=n_r_WAB=5
!        // 0  -- -6   0  6   12
!        // 1  -- -5   1  7   13
!        // 2  -- -4   2  8   14
!        // 3  -- -3   3  9   15
!        // 4  -- -2   4  10  16
!        // 5  -- -1   5  11  17

        if(index_i>=0) then
        Anderson_index=mod(index_i,n_r+1)
        else 
        Anderson_index=mod(index_i,n_r+1)
        Anderson_index=n_r+1+Anderson_index
           if(Anderson_index==n_r+1) then
             Anderson_index=0
           endif
        endif
        return 
       end function Anderson_index
         
      subroutine AB_iterate(converge )
       USE nrtype,only :DP
       USE global_para
       USE mpi
       USE control
       USE constants
       USE utility
       USE mmpi
       !USE mpi_fftw3_operation
       USE matrix_inverse
       implicit none
       logical :: converge
       REAL(DP) :: temp1,temp2,pressure_coeff
       REAL(DP) :: ta_diff_tot,tb_diff_tot,tm_diff_tot
       integer :: k,error1,i,j
       REAL(DP), DIMENSION(:), ALLOCATABLE :: w1A,w2A,w1B,w2B
       REAL(DP), DIMENSION(:), ALLOCATABLE :: WA_iter,DWA_iter
       REAL(DP), DIMENSION(:), ALLOCATABLE :: WB_iter,DWB_iter
       REAL(DP), DIMENSION(:,:), ALLOCATABLE :: U_WAB_anderson
       REAL(DP), DIMENSION(:,:), ALLOCATABLE :: U_M_anderson
       REAL(DP), DIMENSION(:), ALLOCATABLE :: V_WAB_anderson
       REAL(DP), DIMENSION(:), ALLOCATABLE :: V_M_anderson
       REAL(DP), DIMENSION(:,:), ALLOCATABLE :: temp3m(:,:)
       REAL(DP) :: error_anderson_AB,error_anderson_A,error_anderson_B
       REAL(DP) :: error_anderson_M,tm11,tm33,tm12,tm23,tm13,sumaW,sumaWtot
       integer :: n_r_WAB_temp,n_r_M_temp,k_anderson,m_anderson,n_anderson,k_m,k_n
       integer :: m_Manderson,n_Manderson,k_Manderson,Anderson_nim_MM
       REAL(DP),external :: dot_product_mpi
       integer ,external :: Anderson_index
  
error1=0
allocate(w1A(0:LOCAL_SIZE-1),stat=error1)
if(error1/=0) then
write(mystd,*) "allocate dw failed! stop"
stop
endif

allocate(w1B(0:LOCAL_SIZE-1),stat=error1)
if(error1/=0) then
write(mystd,*) "allocate dw failed! stop"
stop
endif

allocate(w2A(0:LOCAL_SIZE-1),stat=error1)
if(error1/=0) then
write(mystd,*) "allocate dw failed! stop"
stop
endif

allocate(w2B(0:LOCAL_SIZE-1),stat=error1)
if(error1/=0) then
write(mystd,*) "allocate dw failed! stop"
stop
endif

allocate(WA_iter(0:LOCAL_SIZE-1),stat=error1)
if(error1/=0) then
write(mystd,*) "allocate dw failed! stop"
stop
endif

allocate(DWA_iter(0:LOCAL_SIZE-1),stat=error1)
if(error1/=0) then
write(mystd,*) "allocate dw failed! stop"
stop
endif


allocate(WB_iter(0:LOCAL_SIZE-1),stat=error1)
if(error1/=0) then
write(mystd,*) "allocate dw failed! stop"
stop
endif

allocate(DWB_iter(0:LOCAL_SIZE-1),stat=error1)
if(error1/=0) then
write(mystd,*) "allocate dw failed! stop"
stop
endif

if(n_iter>Num_simple_mixing_AB) then
n_r_WAB_temp=min(n_iter-Num_simple_mixing_AB,Anderson_nim_AB)
else
n_r_WAB_temp=1
endif

      

allocate(U_wAB_anderson(1:n_r_WAB_temp,1:n_r_WAB_temp),stat=error1)
if(error1/=0) then
write(mystd,*) "allocate dw failed! stop"
stop
endif

allocate(V_wAB_anderson(1:n_r_WAB_temp),stat=error1)
if(error1/=0) then
write(mystd,*) "allocate dw failed! stop"
stop
endif


allocate(temp3m(0:N_dim_ddm-1,0:N_dim_ddm-1),stat=error1)
if(error1/=0) then
write(mystd,*) "allocate dw failed! stop"
stop
endif


       
         ta_diff=0.0
         tb_diff=0.0
         sumaW=0.0
         sumaWtot=0.0

       if(n_iter<Num_simple_mixing_AB) then

          !!!doing simple mixing for A,B

           do k=0,LOCAL_SIZE-1
               pressure_coeff=0.5*(WA(k)+WB(k))
               temp1=NXab*(RB(k)-fB)+pressure_coeff-WA(k)
               temp2=NXab*(RA(k)-fA)+pressure_coeff-WB(k)
               WA(k)=WA(k)+lanbtWA*temp1
               WB(k)=WB(k)+lanbtWB*temp2
		if(abs(temp1)>ta_diff) ta_diff=abs(temp1)
                if(abs(temp2)>tb_diff) tb_diff=abs(temp2)
           enddo


              ta_diff_tot=0.0
              tb_diff_tot=0.0
              call mp_barrier()
              call mp_allreduce(ta_diff,ta_diff_tot)
              call mp_allreduce(tb_diff,tb_diff_tot)
              call mp_barrier()

              ta_diff=ta_diff_tot/nprocs
              tb_diff=tb_diff_tot/nprocs
              

          else if(n_iter==Num_simple_mixing_AB) then
          !after done simple mixing ,dW should be prepaired for Anderson mixing  
               
                WA_out(:,0)=NXab*RB(:)+0.5*(WA(:)+WB(:)-NXab)
                dW_A(:,0)=WA_out(:,0)-WA(:)
                WB_out(:,0)=NXab*RA(:)+0.5*(WA(:)+WB(:)-NXab)
                dW_B(:,0)=WB_out(:,0)-WB(:)
!!!still simple mixing  
           do k=0,LOCAL_SIZE-1
               pressure_coeff=0.5*(WA(k)+WB(k))
               temp1=NXab*(RB(k)-fB)+pressure_coeff-WA(k)
               temp2=NXab*(RA(k)-fA)+pressure_coeff-WB(k)
               WA(k)=WA(k)+lanbtWA*temp1
               WB(k)=WB(k)+lanbtWB*temp2
		if(abs(temp1)>ta_diff) ta_diff=abs(temp1)
                if(abs(temp2)>tb_diff) tb_diff=abs(temp2)
           enddo


              ta_diff_tot=0.0
              tb_diff_tot=0.0
              call mp_barrier()
              call mp_allreduce(ta_diff,ta_diff_tot)
              call mp_allreduce(tb_diff,tb_diff_tot)
              call mp_barrier()

              ta_diff=ta_diff_tot/nprocs
              tb_diff=tb_diff_tot/nprocs
              
             if(myid==0) then
              write(*,*) "ta_diff=",ta_diff
             endif 

       else  !else  if(n_iter>Num_simple_mixing_AB)
             !doing the anderson mixing step

                k_anderson=mod(n_iter-Num_simple_mixing_AB,Anderson_nim_AB+1)


                WA_out(:,k_anderson)=NXab*RB(:)+0.5*(WA(:)+WB(:)-NXab)
                dW_A(:,k_anderson)=WA_out(:,k_anderson)-WA(:)
                WB_out(:,k_anderson)=NXab*RA(:)+0.5*(WA(:)+WB(:)-NXab)
                dW_B(:,k_anderson)=WB_out(:,k_anderson)-WB(:)

                error_anderson_A=(dot_product_mpi(dW_A(:,k_anderson), &
                                       dW_A(:,k_anderson)))
                error_anderson_A=error_anderson_A/(dot_product_mpi(wA_out(:,k_anderson), &
                                       wA_out(:,k_anderson)))
                error_anderson_A=sqrt(abs(error_anderson_A))
                  
                error_anderson_B=(dot_product_mpi(dW_B(:,k_anderson), &
                                       dW_B(:,k_anderson)))
                error_anderson_B=error_anderson_B/(dot_product_mpi(wB_out(:,k_anderson), &
                                       wB_out(:,k_anderson)))
                error_anderson_B=sqrt(abs(error_anderson_B))


               ta_diff=error_anderson_A
               tb_diff=error_anderson_B

              ta_diff_tot=0.0
              tb_diff_tot=0.0
              call mp_barrier()
              call mp_allreduce(ta_diff,ta_diff_tot)
              call mp_allreduce(tb_diff,tb_diff_tot)
              call mp_barrier()

              ta_diff=ta_diff_tot/nprocs
              tb_diff=tb_diff_tot/nprocs

               do m_anderson=1,n_r_WAB_temp

                    k_m=Anderson_index(k_anderson-m_anderson,Anderson_nim_AB)

                    if(myid==0) then
                      if(k_m==k_anderson .or. k_m>Anderson_nim_AB .or. k_m<0) then
                      write(*,*) "terribly wrong of k_m,",k_m
                      endif
                    endif

                    w1A(:)=dW_A(:,k_anderson)-dW_A(:,k_m)
                    w1B(:)=dW_B(:,k_anderson)-dW_B(:,k_m)
                 do n_anderson=m_anderson,n_r_WAB_temp

                    k_n=Anderson_index(k_anderson-n_anderson,Anderson_nim_AB)
                    w2A(:)=dw_A(:,k_anderson)-dW_A(:,k_n)
                    w2B(:)=dw_B(:,k_anderson)-dW_B(:,k_n)
               U_WAB_Anderson(m_anderson,n_anderson)=dot_product_mpi(w1A,w2A)+dot_product_mpi(w1B,w2B)
                 enddo
               enddo
               
              do m_anderson=1,n_r_WAB_temp !!this is different from jiang's,watch out!
                  do n_anderson=1,m_anderson-1
                   U_WAB_Anderson(m_anderson,n_anderson)=U_WAB_Anderson(n_anderson,m_anderson)
                  enddo
               enddo

               w1A(:)=dW_A(:,k_anderson)
               w1B(:)=dW_B(:,k_anderson)
               do m_anderson=1,n_r_WAB_temp
                    k_m=Anderson_index(k_anderson-m_anderson,Anderson_nim_AB)
               w2A(:)=dW_A(:,k_anderson)-dW_A(:,k_m)
               w2B(:)=dW_B(:,k_anderson)-dW_B(:,k_m)
               V_WAB_Anderson(m_anderson)=dot_product_mpi(w2A,w1A)+ &
                                 dot_product_mpi(w2B,w1B)
               enddo

   !       if(myid==0) then
   !       write(*,*) "U_WAB()",U_WAB_Anderson(1,1),U_WAB_Anderson(1,10),U_WAB_Anderson(10,10)
   !       endif
   !      
   !       if(myid==0) then
   !       open(unit=34,file='U_WAB.dat',status='replace')
   !        do i=1,Anderson_nim_AB
   !         do j=1,Anderson_nim_AB
   !           write(34,*) i,j,U_WAB_Anderson(i,j)
   !         enddo
   !        enddo
   !       close(34)
   !       endif
# if defined (Debug)
           if(myid==0) then

              do m_anderson=1,n_r_WAB_temp
                  write(*,*) "m",m_anderson,"V",V_WAB_Anderson(m_anderson)
                  do n_anderson=1,n_r_WAB_temp
                  write(*,*) "m,n",m_anderson,n_anderson,"U",U_WAB_Anderson(m_anderson,n_anderson)
                 enddo
              enddo
            endif
# endif /* Debug */

          call lapacksolver(U_WAB_Anderson,n_r_WAB_temp,n_r_WAB_temp,V_WAB_Anderson) 

# if defined (Debug)
           if(myid==0) then
    write(*,*) "n_iter=",n_iter,"V(1),V(temp)=",V_WAB_Anderson(1),V_WAB_Anderson(n_r_WAB_temp)   
            endif
# endif /* Debug */

             wA_iter=0.0
             DWA_iter=0.0
             wB_iter=0.0
             DWB_iter=0.0

             do m_anderson=1,n_r_WAB_temp
                    k_m=Anderson_index(k_anderson-m_anderson,Anderson_nim_AB)
               wA_iter(:)=WA_iter(:)+V_WAB_Anderson(m_anderson)*(wA_out(:,k_m)-wA_out(:,k_anderson)) 
               wB_iter(:)=WB_iter(:)+V_WAB_Anderson(m_anderson)*(wB_out(:,k_m)-wB_out(:,k_anderson))
             enddo
               wA_iter(:)=WA_iter(:)+wA_out(:,k_anderson)
               wB_iter(:)=WB_iter(:)+wB_out(:,k_anderson)

              
             do m_anderson=1,n_r_WAB_temp
                    k_m=Anderson_index(k_anderson-m_anderson,Anderson_nim_AB)
               DwA_iter(:)=DWA_iter(:)+V_WAB_Anderson(m_anderson)*(DW_A(:,k_m)-DW_A(:,k_anderson)) 
               DwB_iter(:)=DWB_iter(:)+V_WAB_Anderson(m_anderson)*(DW_B(:,k_m)-DW_B(:,k_anderson)) 
             enddo
               DwA_iter(:)=DWA_iter(:)+dW_A(:,k_anderson)
               DwB_iter(:)=DWB_iter(:)+dW_B(:,k_anderson)
               !!!!the final mixing of the new and old w field

                wA=WA_iter+lambda_WAB_anderson*DWA_iter

                wB=WB_iter+lambda_WAB_anderson*DWB_iter

   endif  

   !!!!!for M iteration
 ! the M iteration only be conducted when (mod(n_iter,n_WABiter_per_M)==0,otherwise,do nothing
 !IF NMu is way too small,M iteration is skiped.
      
if(abs(NMu-0.0)<1.0e-5) then
                M_OP=0.0000010d0
# if defined (Debug)
            if(myid==0) then
          write(*,*) "warning!,MOP==0"
            endif
# endif /* Debug */
 else !else if Nmu >1.0e-5
      !only when n_iter/n_WABiter_per_M=0, M is iterated 
    ! /// if(mod(n_iter,n_WABiter_per_M)==0 .or. n_iter_M<=Num_simple_mixing_M+Anderson_nim_M) then
          
       if(n_iter<Num_simple_mixing_M) then

            !!doing the simple mixing for M
            tm_diff=0.0
              do i=0,N_dim_ddm-1
                do j=0,N_dim_ddm-1
                  do k=0,LOCAL_SIZE-1
                      temp3m(j,i)=(NMu*S_OP(k,j,i)-M_OP(k,j,i))
                      M_OP(k,j,i)=M_OP(k,j,i)+lanbtM*temp3m(j,i)
                     if(abs(temp3m(j,i))>tm_diff) tm_diff=abs(temp3m(j,i))
                  enddo
                enddo
               enddo
              tm_diff_tot=0.0
              call mp_barrier()
              call mp_allreduce(tm_diff,tm_diff_tot)
              call mp_barrier()
              tm_diff=tm_diff_tot/nprocs
             if(myid==0) then
             write(*,*) "tm_diff= in simple",tm_diff
             endif 
       else if(n_iter==Num_simple_mixing_M) then

                n_iter_M=n_iter_M+1
                if(method=='A') then
                !uniaxial phase: Smectic A
                M33_out(:,0)=S_OP(:,2,2)*NMu
                d33_anderson(:,0)= M33_out(:,0)-M_OP(:,2,2)
                else if(method=='B') then
                !!uniaxial phase :Smetic C
                M33_out(:,0)=S_OP(:,2,2)*NMu
                M12_out(:,0)=S_OP(:,0,1)*NMu
                M13_out(:,0)=S_OP(:,0,2)*NMu
                M23_out(:,0)=S_OP(:,1,2)*NMu
               
                d33_anderson(:,0)= & 
                M33_out(:,0)-M_OP(:,2,2)

                d12_anderson(:,0)= & 
                M12_out(:,0)-M_OP(:,0,1)

                d13_anderson(:,0)= & 
                M13_out(:,0)-M_OP(:,0,2)

                d23_anderson(:,0)= & 
                M23_out(:,0)-M_OP(:,1,2)
               else 
               !biaxial phase 
                 
                M11_out(:,0)=S_OP(:,0,0)*NMu
               
                M33_out(:,0)=S_OP(:,2,2)*NMu
                M12_out(:,0)=S_OP(:,0,1)*NMu
                M13_out(:,0)=S_OP(:,0,2)*NMu
                M23_out(:,0)=S_OP(:,1,2)*NMu
               
                d11_anderson(:,0)= & 
                M11_out(:,0)-M_OP(:,0,0)

                d33_anderson(:,0)= & 
                M33_out(:,0)-M_OP(:,2,2)

                d12_anderson(:,0)= & 
                M12_out(:,0)-M_OP(:,0,1)

                d13_anderson(:,0)= & 
                M13_out(:,0)-M_OP(:,0,2)

                d23_anderson(:,0)= & 
                M23_out(:,0)-M_OP(:,1,2)


            endif  !endif method=='A'
         


            !!doing the simple mixing for M
            tm_diff=0.0
              do i=0,N_dim_ddm-1
                do j=0,N_dim_ddm-1
                  do k=0,LOCAL_SIZE-1
                      temp3m(j,i)=(NMu*S_OP(k,j,i)-M_OP(k,j,i))
                      M_OP(k,j,i)=M_OP(k,j,i)+lanbtM*temp3m(j,i)
                     if(abs(temp3m(j,i))>tm_diff) tm_diff=abs(temp3m(j,i))
                  enddo
                enddo
               enddo
              tm_diff_tot=0.0
              call mp_barrier()
              call mp_allreduce(tm_diff,tm_diff_tot)
              call mp_barrier()
              tm_diff=tm_diff_tot/nprocs

             if(myid==0) then
             write(*,*) "tm_diff= in simple",tm_diff
             endif 




   else !else that (n_iter_M>Num_simple_mixing_M)
  !! start to do anderson mixing for M,ie.:
   !   WA_out(:,Anderson_nim_AB+1)=NXab*RB(:)+0.5*(WA(:)+WB(:)-NXab)
   !   dW_A(:,Anderson_nim_AB+1)=WA_out(:,Anderson_nim+1)-WA(:)

        if(mod(n_iter-Num_simple_mixing_M,n_WABiter_per_M)==0) then
              n_iter_M=n_iter_M+1
              n_r_M_temp=min(n_iter_M-1,Anderson_nim_M)

allocate(U_M_anderson(1:n_r_M_temp,1:n_r_M_temp),stat=error1)
if(error1/=0) then
write(mystd,*) "allocate dw failed! stop"
stop
endif

allocate(V_M_anderson(1:n_r_M_temp),stat=error1)
if(error1/=0) then
write(mystd,*) "allocate dw failed! stop"
stop
endif


                k_Manderson=mod(n_iter_M-1,Anderson_nim_M+1)


                if(method=='A') then
                !uniaxial phase: Smectic A
                M33_out(:,k_Manderson)=S_OP(:,2,2)*NMu
                d33_anderson(:,k_Manderson)= & 
                M33_out(:,k_Manderson)-M_OP(:,2,2)
                else if(method=='B') then
                !!uniaxial phase :Smetic C
                M33_out(:,k_Manderson)=S_OP(:,2,2)*NMu
                M12_out(:,k_Manderson)=S_OP(:,0,1)*NMu
                M13_out(:,k_Manderson)=S_OP(:,0,2)*NMu
                M23_out(:,k_Manderson)=S_OP(:,1,2)*NMu
               
                d33_anderson(:,k_Manderson)= & 
                M33_out(:,k_Manderson)-M_OP(:,2,2)

                d12_anderson(:,k_Manderson)= & 
                M12_out(:,k_Manderson)-M_OP(:,0,1)

                d13_anderson(:,k_Manderson)= & 
                M13_out(:,k_Manderson)-M_OP(:,0,2)

                d23_anderson(:,k_Manderson)= & 
                M23_out(:,k_Manderson)-M_OP(:,1,2)
               else 
               !biaxial phase 
                 
                M11_out(:,k_Manderson)=S_OP(:,0,0)*NMu
                M33_out(:,k_Manderson)=S_OP(:,2,2)*NMu
                M12_out(:,k_Manderson)=S_OP(:,0,1)*NMu
                M13_out(:,k_Manderson)=S_OP(:,0,2)*NMu
                M23_out(:,k_Manderson)=S_OP(:,1,2)*NMu
               
                d11_anderson(:,k_Manderson)= & 
                M11_out(:,k_Manderson)-M_OP(:,0,0)

                d33_anderson(:,k_Manderson)= & 
                M33_out(:,k_Manderson)-M_OP(:,2,2)

                d12_anderson(:,k_Manderson)= & 
                M12_out(:,k_Manderson)-M_OP(:,0,1)

                d13_anderson(:,k_Manderson)= & 
                M13_out(:,k_Manderson)-M_OP(:,0,2)

                d23_anderson(:,k_Manderson)= & 
                M23_out(:,k_Manderson)-M_OP(:,1,2)
                

            endif  !endif method=='A'
!!!calculate the error 

               Anderson_nim_MM=k_Manderson !don't ask why,I was stupid

             if(method=='A') then
               error_anderson_M=(dot_product_mpi(d33_anderson(:,Anderson_nim_MM), &
                                       d33_anderson(:,Anderson_nim_MM)))
                error_anderson_M=error_anderson_M/(dot_product_mpi(M33_out(:,Anderson_nim_MM), &
                                       M33_out(:,Anderson_nim_MM)))
                error_anderson_M=sqrt(abs(error_anderson_M))
               tm_diff=error_anderson_M

             else if (method=='B') then !!uniaxial phase: Smectic C
             
               tm33=(dot_product_mpi(d33_anderson(:,Anderson_nim_MM), &
                                       d33_anderson(:,Anderson_nim_MM)))
               tm33=tm33/(dot_product_mpi(M33_out(:,Anderson_nim_MM), &
                                       M33_out(:,Anderson_nim_MM)))
               tm33=sqrt(abs(tm33))

               tm12=(dot_product_mpi(d12_anderson(:,Anderson_nim_MM), &
                                       d12_anderson(:,Anderson_nim_MM)))
               tm12=tm12/(dot_product_mpi(M12_out(:,Anderson_nim_MM), &
                                       M12_out(:,Anderson_nim_MM)))
               tm12=sqrt(abs(tm12))

               tm13=(dot_product_mpi(d13_anderson(:,Anderson_nim_MM), &
                                       d13_anderson(:,Anderson_nim_MM)))
               tm13=tm13/(dot_product_mpi(M13_out(:,Anderson_nim_MM), &
                                       M13_out(:,Anderson_nim_MM)))
               tm13=sqrt(abs(tm13))

               
               tm23=(dot_product_mpi(d23_anderson(:,Anderson_nim_MM), &
                                       d23_anderson(:,Anderson_nim_MM)))
               tm23=tm23/(dot_product_mpi(M23_out(:,Anderson_nim_MM), &
                                       M23_out(:,Anderson_nim_MM)))
               tm23=sqrt(abs(tm23))

               tm_diff=max(tm33,tm12,tm13,tm23)

             else !! biaxial phase
              
               tm11=(dot_product_mpi(d11_anderson(:,Anderson_nim_MM), &
                                       d11_anderson(:,Anderson_nim_MM)))
               tm11=tm11/(dot_product_mpi(M11_out(:,Anderson_nim_MM), &
                                       M11_out(:,Anderson_nim_MM)))
               tm11=sqrt(abs(tm11))
                          
              
               tm33=(dot_product_mpi(d33_anderson(:,Anderson_nim_MM), &
                                       d33_anderson(:,Anderson_nim_MM)))
               tm33=tm33/(dot_product_mpi(M33_out(:,Anderson_nim_MM), &
                                       M33_out(:,Anderson_nim_MM)))
               tm33=sqrt(abs(tm33))

               tm12=(dot_product_mpi(d12_anderson(:,Anderson_nim_MM), &
                                       d12_anderson(:,Anderson_nim_MM)))
               tm12=tm12/(dot_product_mpi(M12_out(:,Anderson_nim_MM), &
                                       M12_out(:,Anderson_nim_MM)))
               tm12=sqrt(abs(tm12))

               tm13=(dot_product_mpi(d13_anderson(:,Anderson_nim_MM), &
                                       d13_anderson(:,Anderson_nim_MM)))
               tm13=tm13/(dot_product_mpi(M13_out(:,Anderson_nim_MM), &
                                       M13_out(:,Anderson_nim_MM)))
               tm13=sqrt(abs(tm13))

               
               tm23=(dot_product_mpi(d23_anderson(:,Anderson_nim_MM), &
                                       d23_anderson(:,Anderson_nim_MM)))
               tm23=tm23/(dot_product_mpi(M23_out(:,Anderson_nim_MM), &
                                       M23_out(:,Anderson_nim_MM)))
               tm23=sqrt(abs(tm23))

               tm_diff=max(tm33,tm12,tm13,tm23,tm11)  

               

          endif   

              
              tm_diff_tot=0.0
              call mp_barrier()
              call mp_allreduce(tm_diff,tm_diff_tot)
              call mp_barrier()
              tm_diff=tm_diff_tot/nprocs
             if(myid==0) then
             write(*,*) "tm_diff= anderson",tm_diff
             endif 

            
!!! now generate the U matrix and V vector for M
    U_M_Anderson=0.0
    V_M_Anderson=0.0

     if(method=='A') then
         do m_Manderson=1,n_r_M_temp
                   k_m=Anderson_index(k_Manderson-m_Manderson,Anderson_nim_M)

                    w1A(:)=d33_anderson(:,k_Manderson)-d33_anderson(:,k_m)
                 do n_Manderson=m_Manderson,n_r_M_temp
                   k_n=Anderson_index(k_Manderson-n_Manderson,Anderson_nim_M)

                    w2A(:)=d33_anderson(:,k_Manderson)-d33_anderson(:,k_n)
               U_M_Anderson(m_Manderson,n_Manderson)=dot_product_mpi(w1A,w2A)
                 enddo
               enddo

              do m_Manderson=1,n_r_M_temp
                  do n_Manderson=1,m_Manderson-1
                   U_M_Anderson(m_Manderson,n_Manderson)=U_M_Anderson(n_Manderson,m_Manderson)
                  enddo
               enddo

       else if(method=='B') then
            
         do m_Manderson=1,n_r_M_temp
                   k_m=Anderson_index(k_Manderson-m_Manderson,Anderson_nim_M)
                    w1A(:)=d33_anderson(:,k_Manderson)-d33_anderson(:,k_m)
                 do n_Manderson=m_Manderson,n_r_M_temp
                   k_n=Anderson_index(k_Manderson-n_Manderson,Anderson_nim_M)

                    w2A(:)=d33_anderson(:,k_Manderson)-d33_anderson(:,k_n)
               U_M_Anderson(m_Manderson,n_Manderson)=dot_product_mpi(w1A,w2A)
                 enddo
               enddo


         do m_Manderson=1,n_r_M_temp
                   k_m=Anderson_index(k_Manderson-m_Manderson,Anderson_nim_M)
                    w1A(:)=d12_anderson(:,k_Manderson)-d12_anderson(:,k_m)
                 do n_Manderson=m_Manderson,n_r_M_temp
                   k_n=Anderson_index(k_Manderson-n_Manderson,Anderson_nim_M)
                    w2A(:)=d12_anderson(:,k_Manderson)-d12_anderson(:,k_n)
  U_M_Anderson(m_Manderson,n_Manderson)=U_M_Anderson(m_Manderson,n_Manderson) + &
             dot_product_mpi(w1A,w2A)
                 enddo
               enddo
          

         do m_Manderson=1,n_r_M_temp
                   k_m=Anderson_index(k_Manderson-m_Manderson,Anderson_nim_M)
                    w1A(:)=d13_anderson(:,k_Manderson)-d13_anderson(:,k_m)
                 do n_Manderson=m_Manderson,n_r_M_temp
                   k_n=Anderson_index(k_Manderson-n_Manderson,Anderson_nim_M)
                    w2A(:)=d13_anderson(:,k_Manderson)-d13_anderson(:,k_n)
               U_M_Anderson(m_Manderson,n_Manderson)=U_M_Anderson(m_Manderson,n_Manderson)+ &
              dot_product_mpi(w1A,w2A)
                 enddo
               enddo

         do m_Manderson=1,n_r_M_temp
                   k_m=Anderson_index(k_Manderson-m_Manderson,Anderson_nim_M)
                    w1A(:)=d23_anderson(:,k_Manderson)-d23_anderson(:,k_m)
                 do n_Manderson=m_Manderson,n_r_M_temp
                   k_n=Anderson_index(k_Manderson-n_Manderson,Anderson_nim_M)
                    w2A(:)=d23_anderson(:,k_Manderson)-d23_anderson(:,k_n)
               U_M_Anderson(m_Manderson,n_Manderson)=U_M_Anderson(m_Manderson,n_Manderson)+ &
             dot_product_mpi(w1A,w2A)
                 enddo
               enddo


              do m_Manderson=1,n_r_M_temp
                  do n_Manderson=1,m_Manderson-1
                   U_M_Anderson(m_Manderson,n_Manderson)=U_M_Anderson(n_Manderson,m_Manderson)
                  enddo
               enddo





        else !//biaxial phase
       

         do m_Manderson=1,n_r_M_temp
                   k_m=Anderson_index(k_Manderson-m_Manderson,Anderson_nim_M)
                    w1A(:)=d33_anderson(:,k_Manderson)-d33_anderson(:,k_m)
                 do n_Manderson=m_Manderson,n_r_M_temp
                   k_n=Anderson_index(k_Manderson-n_Manderson,Anderson_nim_M)

                    w2A(:)=d33_anderson(:,k_Manderson)-d33_anderson(:,k_n)
               U_M_Anderson(m_Manderson,n_Manderson)=dot_product_mpi(w1A,w2A)
                 enddo
           enddo
             
         do m_Manderson=1,n_r_M_temp
                   k_m=Anderson_index(k_Manderson-m_Manderson,Anderson_nim_M)
                    w1A(:)=d12_anderson(:,k_Manderson)-d12_anderson(:,k_m)
                 do n_Manderson=m_Manderson,n_r_M_temp
                   k_n=Anderson_index(k_Manderson-n_Manderson,Anderson_nim_M)
                    w2A(:)=d12_anderson(:,k_Manderson)-d12_anderson(:,k_n)
    U_M_Anderson(m_Manderson,n_Manderson)=U_M_Anderson(m_Manderson,n_Manderson)+ &
      dot_product_mpi(w1A,w2A)
                 enddo
               enddo
          
         do m_Manderson=1,n_r_M_temp
                   k_m=Anderson_index(k_Manderson-m_Manderson,Anderson_nim_M)
                    w1A(:)=d13_anderson(:,k_Manderson)-d13_anderson(:,k_m)
                 do n_Manderson=m_Manderson,n_r_M_temp
                   k_n=Anderson_index(k_Manderson-n_Manderson,Anderson_nim_M)
                    w2A(:)=d13_anderson(:,k_Manderson)-d13_anderson(:,k_n)
    U_M_Anderson(m_Manderson,n_Manderson)=U_M_Anderson(m_Manderson,n_Manderson)+ & 
       dot_product_mpi(w1A,w2A)
                 enddo
               enddo
          

         do m_Manderson=1,n_r_M_temp
                   k_m=Anderson_index(k_Manderson-m_Manderson,Anderson_nim_M)
                    w1A(:)=d23_anderson(:,k_Manderson)-d23_anderson(:,k_m)
                 do n_Manderson=m_Manderson,n_r_M_temp
                   k_n=Anderson_index(k_Manderson-n_Manderson,Anderson_nim_M)
                    w2A(:)=d23_anderson(:,k_Manderson)-d23_anderson(:,k_n)
    U_M_Anderson(m_Manderson,n_Manderson)=U_M_Anderson(m_Manderson,n_Manderson)+ &
    dot_product_mpi(w1A,w2A)
                 enddo
               enddo
          
         do m_Manderson=1,n_r_M_temp
                   k_m=Anderson_index(k_Manderson-m_Manderson,Anderson_nim_M)
                    w1A(:)=d11_anderson(:,k_Manderson)-d11_anderson(:,k_m)
                 do n_Manderson=m_Manderson,n_r_M_temp
                   k_n=Anderson_index(k_Manderson-n_Manderson,Anderson_nim_M)
                    w2A(:)=d11_anderson(:,k_Manderson)-d11_anderson(:,k_n)
    U_M_Anderson(m_Manderson,n_Manderson)=U_M_Anderson(m_Manderson,n_Manderson)+ & 
    dot_product_mpi(w1A,w2A)
                 enddo
               enddo
          
          
              do m_Manderson=1,n_r_M_temp
                  do n_Manderson=1,m_Manderson-1
                   U_M_Anderson(m_Manderson,n_Manderson)=U_M_Anderson(n_Manderson,m_Manderson)
                  enddo
               enddo
             




    endif  !if(method=='A')
!!! now to generate the V_M_Anderson vector
        if(method=='A') then

               w1A(:)=d33_anderson(:,k_Manderson)
           do m_Manderson=1,n_r_M_temp
                   k_m=Anderson_index(k_Manderson-m_Manderson,Anderson_nim_M)
               w2A(:)=d33_anderson(:,k_Manderson)-d33_anderson(:,k_m)
               V_M_Anderson(m_Manderson)=dot_product_mpi(w2A,w1A)
               enddo

         else if(method=='B') then
              
               w1A(:)=d33_anderson(:,k_Manderson)
           do m_Manderson=1,n_r_M_temp
                   k_m=Anderson_index(k_Manderson-m_Manderson,Anderson_nim_M)
               w2A(:)=d33_anderson(:,k_Manderson)-d33_anderson(:,k_m)
               V_M_Anderson(m_Manderson)=dot_product_mpi(w2A,w1A)
               enddo


               w1A(:)=d13_anderson(:,k_Manderson)
           do m_Manderson=1,n_r_M_temp
                   k_m=Anderson_index(k_Manderson-m_Manderson,Anderson_nim_M)
               w2A(:)=d13_anderson(:,k_Manderson)-d13_anderson(:,k_m)
               V_M_Anderson(m_Manderson)=V_M_anderson(m_Manderson)+dot_product_mpi(w2A,w1A)
               enddo

               w1A(:)=d12_anderson(:,k_Manderson)
           do m_Manderson=1,n_r_M_temp
                   k_m=Anderson_index(k_Manderson-m_Manderson,Anderson_nim_M)
               w2A(:)=d12_anderson(:,k_Manderson)-d12_anderson(:,k_m)
               V_M_Anderson(m_Manderson)=V_M_anderson(m_Manderson)+dot_product_mpi(w2A,w1A)
               enddo

               w1A(:)=d23_anderson(:,k_Manderson)
           do m_Manderson=1,n_r_M_temp
                   k_m=Anderson_index(k_Manderson-m_Manderson,Anderson_nim_M)
               w2A(:)=d23_anderson(:,k_Manderson)-d23_anderson(:,k_m)
               V_M_Anderson(m_Manderson)=V_M_anderson(m_Manderson)+dot_product_mpi(w2A,w1A)
               enddo


        else !biaxial

               w1A(:)=d33_anderson(:,k_Manderson)
           do m_Manderson=1,n_r_M_temp
                   k_m=Anderson_index(k_Manderson-m_Manderson,Anderson_nim_M)
               w2A(:)=d33_anderson(:,k_Manderson)-d33_anderson(:,k_m)
               V_M_Anderson(m_Manderson)=dot_product_mpi(w2A,w1A)
               enddo


               w1A(:)=d13_anderson(:,k_Manderson)
           do m_Manderson=1,n_r_M_temp
                   k_m=Anderson_index(k_Manderson-m_Manderson,Anderson_nim_M)
               w2A(:)=d13_anderson(:,k_Manderson)-d13_anderson(:,k_m)
               V_M_Anderson(m_Manderson)=V_M_anderson(m_Manderson)+dot_product_mpi(w2A,w1A)
               enddo

               w1A(:)=d12_anderson(:,k_Manderson)
           do m_Manderson=1,n_r_M_temp
                   k_m=Anderson_index(k_Manderson-m_Manderson,Anderson_nim_M)
               w2A(:)=d12_anderson(:,k_Manderson)-d12_anderson(:,k_m)
               V_M_Anderson(m_Manderson)=V_M_anderson(m_Manderson)+dot_product_mpi(w2A,w1A)
               enddo

               w1A(:)=d23_anderson(:,k_Manderson)
           do m_Manderson=1,n_r_M_temp
                   k_m=Anderson_index(k_Manderson-m_Manderson,Anderson_nim_M)
               w2A(:)=d23_anderson(:,k_Manderson)-d23_anderson(:,k_m)
               V_M_Anderson(m_Manderson)=V_M_anderson(m_Manderson)+dot_product_mpi(w2A,w1A)
               enddo

               w1A(:)=d11_anderson(:,k_Manderson)
           do m_Manderson=1,n_r_M_temp
                   k_m=Anderson_index(k_Manderson-m_Manderson,Anderson_nim_M)
               w2A(:)=d11_anderson(:,k_Manderson)-d11_anderson(:,k_m)
               V_M_Anderson(m_Manderson)=V_M_anderson(m_Manderson)+dot_product_mpi(w2A,w1A)
               enddo
              
         endif  


          call lapacksolver(U_M_Anderson,n_r_M_temp,n_r_M_temp,V_M_Anderson)
         

          if(method=='A') then
          
             wA_iter=0.0
             DWA_iter=0.0

           do m_Manderson=1,n_r_M_temp
                   k_m=Anderson_index(k_Manderson-m_Manderson,Anderson_nim_M)
               wA_iter(:)=WA_iter(:)+V_M_Anderson(m_Manderson)*(M33_out(:,k_m)-M33_out(:,k_Manderson)) 
             enddo
               wA_iter(:)=WA_iter(:)+M33_out(:,k_Manderson)

              
           do m_Manderson=1,n_r_M_temp
                   k_m=Anderson_index(k_Manderson-m_Manderson,Anderson_nim_M)
               DwA_iter(:)=DWA_iter(:)+V_M_Anderson(m_Manderson)*(d33_anderson(:,k_m)- &
                            d33_anderson(:,k_Manderson)) 
             enddo
               DwA_iter(:)=DWA_iter(:)+d33_anderson(:,k_Manderson)

               !!!!the final mixing of the new and old w field
            do k=0,LOCAL_SIZE-1
                M_OP(k,2,2) =WA_iter(k)+lambda_M_anderson*DWA_iter(k)
                M_OP(k,0,0)=-0.5*M_OP(k,2,2)
                M_OP(k,1,1)=M_OP(k,0,0)
                M_OP(k,1,0)=0.0
                M_OP(k,0,1)=0.0
                M_OP(k,2,0)=0.0
                M_OP(k,0,2)=0.0
                M_OP(k,2,1)=0.0
                M_OP(k,1,2)=0.0
            enddo

       else if (method=='B') then
                         
             wA_iter=0.0
             DWA_iter=0.0

           do m_Manderson=1,n_r_M_temp
                   k_m=Anderson_index(k_Manderson-m_Manderson,Anderson_nim_M)
               wA_iter(:)=WA_iter(:)+V_M_Anderson(m_Manderson)*(M33_out(:,k_m)-M33_out(:,k_Manderson)) 
             enddo
               wA_iter(:)=WA_iter(:)+M33_out(:,k_Manderson)

           do m_Manderson=1,n_r_M_temp
                   k_m=Anderson_index(k_Manderson-m_Manderson,Anderson_nim_M)
               DwA_iter(:)=DWA_iter(:)+V_M_Anderson(m_Manderson)*(d33_anderson(:,k_m)- &
                            d33_anderson(:,k_Manderson)) 
             enddo
               DwA_iter(:)=DWA_iter(:)+d33_anderson(:,k_Manderson)
              

                 
            do k=0,LOCAL_SIZE-1
                M_OP(k,2,2) =WA_iter(k)+lambda_M_anderson*DWA_iter(k)
                M_OP(k,0,0)=-0.5*M_OP(k,2,2)
                M_OP(k,1,1)=M_OP(k,0,0)
            enddo
            
             
             wA_iter=0.0
             DWA_iter=0.0

           do m_Manderson=1,n_r_M_temp
                   k_m=Anderson_index(k_Manderson-m_Manderson,Anderson_nim_M)
               wA_iter(:)=WA_iter(:)+V_M_Anderson(m_Manderson)*(M12_out(:,k_m)-M12_out(:,k_Manderson)) 
             enddo
               wA_iter(:)=WA_iter(:)+M12_out(:,k_Manderson)

           do m_Manderson=1,n_r_M_temp
                   k_m=Anderson_index(k_Manderson-m_Manderson,Anderson_nim_M)
               DwA_iter(:)=DWA_iter(:)+V_M_Anderson(m_Manderson)*(d12_anderson(:,k_m)- &
                            d12_anderson(:,k_Manderson)) 
             enddo
               DwA_iter(:)=DWA_iter(:)+d12_anderson(:,k_Manderson)
              
                 
            do k=0,LOCAL_SIZE-1
                M_OP(k,0,1) =WA_iter(k)+lambda_M_anderson*DWA_iter(k)
                M_OP(k,1,0)=M_OP(k,0,1)
            enddo


             wA_iter=0.0
             DWA_iter=0.0

           do m_Manderson=1,n_r_M_temp
                   k_m=Anderson_index(k_Manderson-m_Manderson,Anderson_nim_M)
               wA_iter(:)=WA_iter(:)+V_M_Anderson(m_Manderson)*(M13_out(:,k_m)-M13_out(:,k_Manderson)) 
             enddo
               wA_iter(:)=WA_iter(:)+M13_out(:,k_Manderson)

           do m_Manderson=1,n_r_M_temp
                   k_m=Anderson_index(k_Manderson-m_Manderson,Anderson_nim_M)
               DwA_iter(:)=DWA_iter(:)+V_M_Anderson(m_Manderson)*(d13_anderson(:,k_m)- &
                            d13_anderson(:,k_Manderson)) 
             enddo
               DwA_iter(:)=DWA_iter(:)+d13_anderson(:,k_Manderson)
              
               !!!!the final mixing of the new and old w field
                 
                 
            do k=0,LOCAL_SIZE-1
                M_OP(k,0,2) =WA_iter(k)+lambda_M_anderson*DWA_iter(k)
                M_OP(k,2,0)=M_OP(k,0,2)
            enddo

             wA_iter=0.0
             DWA_iter=0.0

           do m_Manderson=1,n_r_M_temp
                   k_m=Anderson_index(k_Manderson-m_Manderson,Anderson_nim_M)
               wA_iter(:)=WA_iter(:)+V_M_Anderson(m_Manderson)*(M23_out(:,k_m)-M23_out(:,k_Manderson)) 
             enddo
               wA_iter(:)=WA_iter(:)+M23_out(:,k_Manderson)

           do m_Manderson=1,n_r_M_temp
                   k_m=Anderson_index(k_Manderson-m_Manderson,Anderson_nim_M)
               DwA_iter(:)=DWA_iter(:)+V_M_Anderson(m_Manderson)*(d23_anderson(:,k_m)- &
                            d23_anderson(:,k_Manderson)) 
             enddo
               DwA_iter(:)=DWA_iter(:)+d23_anderson(:,k_Manderson)
              
                 
            do k=0,LOCAL_SIZE-1
                M_OP(k,1,2) =WA_iter(k)+lambda_M_anderson*DWA_iter(k)
                M_OP(k,2,1)=M_OP(k,1,2)
            enddo


     else  !biaxial phase

             wA_iter=0.0
             DWA_iter=0.0


           do m_Manderson=1,n_r_M_temp
                   k_m=Anderson_index(k_Manderson-m_Manderson,Anderson_nim_M)
               wA_iter(:)=WA_iter(:)+V_M_Anderson(m_Manderson)*(M11_out(:,k_m)-M11_out(:,k_Manderson)) 
             enddo
               wA_iter(:)=WA_iter(:)+M11_out(:,k_Manderson)

           do m_Manderson=1,n_r_M_temp
                   k_m=Anderson_index(k_Manderson-m_Manderson,Anderson_nim_M)
               DwA_iter(:)=DWA_iter(:)+V_M_Anderson(m_Manderson)*(d11_anderson(:,k_m)- &
                            d11_anderson(:,k_Manderson)) 
             enddo
               DwA_iter(:)=DWA_iter(:)+d11_anderson(:,k_Manderson)
              

                 
               !!!!the final mixing of the new and old w field
                 
            do k=0,LOCAL_SIZE-1
                M_OP(k,0,0) =WA_iter(k)+lambda_M_anderson*DWA_iter(k)
            enddo
             
             wA_iter=0.0
             DWA_iter=0.0

           do m_Manderson=1,n_r_M_temp
                   k_m=Anderson_index(k_Manderson-m_Manderson,Anderson_nim_M)
               wA_iter(:)=WA_iter(:)+V_M_Anderson(m_Manderson)*(M33_out(:,k_m)-M33_out(:,k_Manderson)) 
             enddo
               wA_iter(:)=WA_iter(:)+M33_out(:,k_Manderson)

           do m_Manderson=1,n_r_M_temp
                   k_m=Anderson_index(k_Manderson-m_Manderson,Anderson_nim_M)
               DwA_iter(:)=DWA_iter(:)+V_M_Anderson(m_Manderson)*(d33_anderson(:,k_m)- &
                            d33_anderson(:,k_Manderson)) 
             enddo
               DwA_iter(:)=DWA_iter(:)+d33_anderson(:,k_Manderson)
              
               !!!!the final mixing of the new and old w field
                 
            do k=0,LOCAL_SIZE-1
                M_OP(k,2,2) =WA_iter(k)+lambda_M_anderson*DWA_iter(k)
                M_OP(k,1,1) =-M_OP(k,0,0)-M_OP(k,2,2)
            enddo
            
             
             wA_iter=0.0
             DWA_iter=0.0

           do m_Manderson=1,n_r_M_temp
                   k_m=Anderson_index(k_Manderson-m_Manderson,Anderson_nim_M)
               wA_iter(:)=WA_iter(:)+V_M_Anderson(m_Manderson)*(M12_out(:,k_m)-M12_out(:,k_Manderson)) 
             enddo
               wA_iter(:)=WA_iter(:)+M12_out(:,k_Manderson)

           do m_Manderson=1,n_r_M_temp
                   k_m=Anderson_index(k_Manderson-m_Manderson,Anderson_nim_M)
               DwA_iter(:)=DWA_iter(:)+V_M_Anderson(m_Manderson)*(d12_anderson(:,k_m)- &
                            d12_anderson(:,k_Manderson)) 
             enddo
               DwA_iter(:)=DWA_iter(:)+d12_anderson(:,k_Manderson)
              
               !!!!the final mixing of the new and old w field
                 
                 
            do k=0,LOCAL_SIZE-1
                M_OP(k,0,1) =WA_iter(k)+lambda_M_anderson*DWA_iter(k)
                M_OP(k,1,0)=M_OP(k,0,1)
            enddo

             wA_iter=0.0
             DWA_iter=0.0

           do m_Manderson=1,n_r_M_temp
                   k_m=Anderson_index(k_Manderson-m_Manderson,Anderson_nim_M)
               wA_iter(:)=WA_iter(:)+V_M_Anderson(m_Manderson)*(M13_out(:,k_m)-M13_out(:,k_Manderson)) 
             enddo
               wA_iter(:)=WA_iter(:)+M13_out(:,k_Manderson)

           do m_Manderson=1,n_r_M_temp
                   k_m=Anderson_index(k_Manderson-m_Manderson,Anderson_nim_M)
               DwA_iter(:)=DWA_iter(:)+V_M_Anderson(m_Manderson)*(d13_anderson(:,k_m)- &
                            d13_anderson(:,k_Manderson)) 
             enddo
               DwA_iter(:)=DWA_iter(:)+d13_anderson(:,k_Manderson)
              
               !!!!the final mixing of the new and old w field
                 
                 
                 
            do k=0,LOCAL_SIZE-1
                M_OP(k,0,2) =WA_iter(k)+lambda_M_anderson*DWA_iter(k)
                M_OP(k,2,0)=M_OP(k,0,2)
            enddo

             wA_iter=0.0
             DWA_iter=0.0

           do m_Manderson=1,n_r_M_temp
                   k_m=Anderson_index(k_Manderson-m_Manderson,Anderson_nim_M)
               wA_iter(:)=WA_iter(:)+V_M_Anderson(m_Manderson)*(M23_out(:,k_m)-M23_out(:,k_Manderson)) 
             enddo
               wA_iter(:)=WA_iter(:)+M23_out(:,k_Manderson)

           do m_Manderson=1,n_r_M_temp
                   k_m=Anderson_index(k_Manderson-m_Manderson,Anderson_nim_M)
               DwA_iter(:)=DWA_iter(:)+V_M_Anderson(m_Manderson)*(d23_anderson(:,k_m)- &
                            d23_anderson(:,k_Manderson)) 
             enddo
               DwA_iter(:)=DWA_iter(:)+d23_anderson(:,k_Manderson)
              
               !!!!the final mixing of the new and old w field
                 
            do k=0,LOCAL_SIZE-1
                M_OP(k,1,2) =WA_iter(k)+lambda_M_anderson*DWA_iter(k)
                M_OP(k,2,1)=M_OP(k,1,2)
            enddo

      endif  

     deallocate(U_M_anderson)
     deallocate(V_M_anderson)

   else !if(mod(n_iter-Num_simple_mixing_M,n_WABiter_per_M)==0)

     tm_diff=0.5000d0

    endif !if(mod(n_iter-Num_simple_mixing_M,n_WABiter_per_M)==0)

    endif  ! if(n_iter_M<=Num_simple_mixing_M) 

endif !end if abs(NMu)<1.0e-4


!!! check whether converge is true  or not

       if(abs(NMu-0.0)>1.0e-5) then

           if(myid==0) then
           write(*,*) "n_iter=",n_iter,"n_iter_M=",n_iter_M
           write(*,*) "ta,b,m_diff=",ta_diff,tb_diff,tm_diff
            endif

            if((ta_diff<CC) .and. (tb_diff<CC) .and. (tm_diff<5*CC)) then
            converge=.true.
               if(myid==0) then
               write(*,*) "converge be true!",ta_diff,tb_diff,tm_diff
               endif 
            else
            converge=.false.
            endif
        else
            if(myid==0) then
           write(*,*) "n_iter=",n_iter
           write(*,*) "ta,b_diff,CC=",ta_diff,tb_diff,CC
            endif
            if((ta_diff<CC) .and. (tb_diff<CC)) then
            converge=.true.
               if(myid==0) then
               write(*,*) "converge be true!",ta_diff,tb_diff
               endif 
             else if (n_iter>=150 .and. (ta_diff<5.0e-4) .and. & 
                       (tb_diff<5.0e-4)) then

             converge=.true.

             
            else
            converge=.false.
            endif
        endif


!!! free the allocated arrays



call mp_barrier()
# if defined (Debug)
if(myid==0) then
write(*,*) "WA(0,30,90) in iter=",WA(0),WA(30),WA(90)
write(*,*) "WB(0,30,90) in iter =",WB(0),WB(30),WB(90)
endif
# endif /* Debug */


 deallocate(w1A)
 deallocate(w1B)
 deallocate(w2A)
 deallocate(w2B)
 deallocate(wA_iter)
 deallocate(wB_iter)
 deallocate(DwA_iter)
 deallocate(DwB_iter)
 deallocate(U_WAB_anderson)
 deallocate(V_WAB_anderson)
! deallocate(U_M_anderson)
! deallocate(V_M_anderson)
deallocate(temp3m)
call mp_barrier()

end subroutine AB_iterate

        FUNCTION dot_product_mpi(w1,w2)
        USE mpi
        USE mmpi
        USE control
        USE constants
        USE global_para,only :dx,dy,dz,SIDEx,SIDEy,SIDEz,M_v,LOCAL_SIZE
# if defined (SLAB)
       USE mpi_fftw3_operation
# else /* SLAB */
       use decomp_fft_mod
# endif /* SLAB */
        USE utility
        implicit none
        integer :: i
        REAL*8 :: w1(0:LOCAL_SIZE-1),w2(0:LOCAL_SIZE-1)
        REAL*8,allocatable ::a31(:)
        REAL*8 :: dot_product_mpi,dot_tot

         allocate(a31(0:LOCAL_SIZE-1))
         do i=0,LOCAL_SIZE-1
          a31(i)=w1(i)*w2(i)
         enddo
         !be careful,M_v is not a constant until dx,dy,dz are fixed

          !M_v=dx*dy*dz*SIDEx*SIDEy*SIDEz

# if defined (SLAB)
# if defined (Dim2)
        dot_product_mpi=simposon_2D_1D_mpi(SIDEy,local_nx,dy,dx,a31)
# else /* Dim2 */
        dot_product_mpi=simposon_3D_1D_mpi(SIDEz,SIDEy,local_nx,dz,dy,dx,a31)
# endif /* Dim2 */
# else /* SLAB */
        dot_product_mpi=simposon_3D_1D_mpi(xsize(1),xsize(2),xsize(3),dx,dy,dz,a31)
# endif /* SLAB */

              dot_tot=0.0
              call mp_barrier()
              call mp_allreduce(dot_product_mpi,dot_tot)
              call mp_barrier()
             dot_product_mpi=dot_tot/M_v

         deallocate(a31)
         return
         end function dot_product_mpi
       

 subroutine mem_cal_scft_destroy()       
 USE nrtype,only :DP,PI
 USE global_para
 USE mpi
 USE control
 USE constants
 USE utility
 USE mmpi
 implicit none
 integer :: allostat
!for R
! deallocate(GA1_inverse_val,stat=allostat)
! deallocate(GA1_inverse_col,stat=allostat)
! deallocate(GA1_inverse_row,stat=allostat)
!
! deallocate(GA2_inverse_val,stat=allostat)
! deallocate(GA2_inverse_col,stat=allostat)
! deallocate(GA2_inverse_row,stat=allostat)
!
 deallocate(GA3_inverse_val,stat=allostat)
 deallocate(GA3_inverse_col,stat=allostat)
 deallocate(GA3_inverse_row,stat=allostat)

 deallocate(GA3_star_inverse_val,stat=allostat)
 deallocate(GA3_star_inverse_col,stat=allostat)
 deallocate(GA3_star_inverse_row,stat=allostat)

 deallocate(GB3_star_inverse_val,stat=allostat)
 deallocate(GB3_star_inverse_col,stat=allostat)
 deallocate(GB3_star_inverse_row,stat=allostat)

! deallocate(GB2_star_inverse_val,stat=allostat)
! deallocate(GB2_star_inverse_col,stat=allostat)
! deallocate(GB2_star_inverse_row,stat=allostat)
!
! deallocate(GB1_star_inverse_val,stat=allostat)
! deallocate(GB1_star_inverse_col,stat=allostat)
! deallocate(GB1_star_inverse_row,stat=allostat)

 deallocate(GB3_inverse_val,stat=allostat)
 deallocate(GB3_inverse_col,stat=allostat)
 deallocate(GB3_inverse_row,stat=allostat)

 if(allostat/=0) then
 write(*,*) "failed to deallocate G_inverse array,check !"
 stop
 endif


 end subroutine mem_cal_scft_destroy

       

