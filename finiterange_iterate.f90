      subroutine finite_range_iterate(converge )
       USE nrtype,only :DP
       USE global_para
       USE control
       USE constants
       USE utility
       USE mmpi
       USE mpi_fftw3_operation
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
       REAL(DP), DIMENSION(:), ALLOCATABLE :: V_WAB_anderson
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

       
         ta_diff=0.0
         tb_diff=0.0
         sumaW=0.0
         sumaWtot=0.0

       if(n_iter<Num_simple_mixing_AB) then
          !!!doing simple mixing for A,B

            call fin_simp_iter(WA_out(:,0),WB_out(:,0))
           do k=0,LOCAL_SIZE-1
               !pressure_coeff=0.5*(WA(k)+WB(k))-Nxab
               !temp1=NXab*(RB(k)-fB)+pressure_coeff-WA(k)
               !temp2=NXab*(RA(k)-fA)+pressure_coeff-WB(k)
               !WA(k)=WA(k)+lanbtWA*temp1
               !WB(k)=WB(k)+lanbtWB*temp2
               !temp1=WA_out(k,0)-NXab*(RB(k))-pressure_coeff
               !temp2=WB_out(k,0)-NXab*(RA(k))-pressure_coeff
               temp1=WA_out(k,0)-WA(k)
               temp2=WB_out(k,0)-WB(k)
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

              write(*,*) "total error",ta_diff,tb_diff

          else if(n_iter==Num_simple_mixing_AB) then
          !after done simple mixing ,dW should be prepaired for Anderson mixing  
               
                call fin_simp_iter(WA_out(:,0),WB_out(:,0))
                !WA_out(:,0)=NXab*RB(:)+0.5*(WA(:)+WB(:)-NXab)
                dW_A(:,0)=WA_out(:,0)-WA(:)
                !WB_out(:,0)=NXab*RA(:)+0.5*(WA(:)+WB(:)-NXab)
                dW_B(:,0)=WB_out(:,0)-WB(:)
!!!still simple mixing  
           do k=0,LOCAL_SIZE-1
               temp1=WA_out(k,0)-WA(k)
               temp2=WB_out(k,0)-WB(k)
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

                call fin_simp_iter(WA_out(:,k_anderson),WB_out(:,k_anderson))

                !WA_out(:,k_anderson)=NXab*RB(:)+0.5*(WA(:)+WB(:)-NXab)
                dW_A(:,k_anderson)=WA_out(:,k_anderson)-WA(:)
                !WB_out(:,k_anderson)=NXab*RA(:)+0.5*(WA(:)+WB(:)-NXab)
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

!!! check whether converge is true  or not

       if(abs(NMu-0.0)>1.0e-4) then

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
             else if (n_iter>=200 .and. (ta_diff<5.0e-3) .and. & 
                       (tb_diff<5.0e-3)) then

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
 call mp_barrier()

end subroutine finite_range_iterate


   subroutine fin_simp_iter(WA_new,WB_new)
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
    REAL(DP) :: temp1,pressure_coeff,M_grid_inv
    integer :: k,error1
    real(dp),intent(inout) :: WA_new(0:Local_size-1)  
    real(dp),intent(inout) :: WB_new(0:Local_size-1)  
    Complex(dp),DIMENSION(:),ALLOCATABLE ::  yita,phi_k_A,phi_k_B
    real(dp) :: error
 
     M_grid_inv=1.0d0/M_grid
     allocate(yita(0:LOCAL_SIZE-1),stat=error1)
     allocate(phi_k_A(0:LOCAL_SIZE-1),stat=error1)
     allocate(phi_k_B(0:LOCAL_SIZE-1),stat=error1)
     if(error1/=0) then
     write(mystd,*) "allocate dw failed! stop"
     stop
     endif
     
        do k=0,LOCAL_SIZE-1
# if defined (Dim2)
          local_data(kD(k)%x,kD(k)%y)=cmplx(RA(k)-fA,0.0)
# else  /* Dim2 */
          local_data(kD(k)%x,kD(k)%y,KD(k)%z)=cmplx(RA(k)-fA,0.0)
# endif /* Dim2 */
         enddo
       call fftw_mpi_execute_dft(fplan, local_data, local_data)
        do k=0,LOCAL_SIZE-1
# if defined (Dim2)
          phi_k_A(k)=local_data(kD(k)%x,kD(k)%y)
# else  /* Dim2 */
          phi_k_A(k)=local_data(kD(k)%x,kD(k)%y,KD(k)%z)
# endif /* Dim2 */
         enddo

        do k=0,LOCAL_SIZE-1
# if defined (Dim2)
          local_data(kD(k)%x,kD(k)%y)=cmplx(RB(k)-fB,0.0)
# else  /* Dim2 */
          local_data(kD(k)%x,kD(k)%y,KD(k)%z)=cmplx(RB(k)-fB,0.0)
# endif /* Dim2 */
         enddo
       call fftw_mpi_execute_dft(fplan, local_data, local_data)
        do k=0,LOCAL_SIZE-1
# if defined (Dim2)
          phi_k_B(k)=local_data(kD(k)%x,kD(k)%y)
# else  /* Dim2 */
          phi_k_B(k)=local_data(kD(k)%x,kD(k)%y,KD(k)%z)
# endif /* Dim2 */
        enddo

        do k=0,LOCAL_SIZE-1
               pressure_coeff=0.5*(WA(k)+WB(k))
# if defined (Dim2)
          local_data(kD(k)%x,kD(k)%y)=cmplx(pressure_coeff,0.0)
# else  /* Dim2 */
          local_data(kD(k)%x,kD(k)%y,KD(k)%z)=cmplx(pressure_coeff,0.0)
# endif /* Dim2 */
         enddo

       call fftw_mpi_execute_dft(fplan, local_data, local_data)
        do k=0,LOCAL_SIZE-1
# if defined (Dim2)
          yita(k)=local_data(kD(k)%x,kD(k)%y)
# else  /* Dim2 */
          yita(k)=local_data(kD(k)%x,kD(k)%y,KD(k)%z)
# endif /* Dim2 */
         enddo
   
           do k=0,LOCAL_SIZE-1
               temp1=dexp(-0.50d0*(sigma**2)*ksquare(k))
               phi_k_B(k)=NXab*temp1*phi_k_B(k)+yita(k)
               phi_k_A(k)=NXab*temp1*phi_k_A(k)+yita(k)
            enddo

        do k=0,LOCAL_SIZE-1
# if defined (Dim2)
          local_data(kD(k)%x,kD(k)%y)=phi_k_B(k)
# else  /* Dim2 */
          local_data(kD(k)%x,kD(k)%y,KD(k)%z)=phi_k_B(k)
# endif /* Dim2 */
         enddo
       call fftw_mpi_execute_dft(bplan, local_data, local_data)
        do k=0,LOCAL_SIZE-1
# if defined (Dim2)
          WA_new(k)=real(local_data(kD(k)%x,kD(k)%y))*M_grid_inv
# else  /* Dim2 */
          WA_new(k)=real(local_data(kD(k)%x,kD(k)%y,KD(k)%z))*M_grid_inv
# endif /* Dim2 */
        enddo

        do k=0,LOCAL_SIZE-1
# if defined (Dim2)
          local_data(kD(k)%x,kD(k)%y)=phi_k_A(k)
# else  /* Dim2 */
          local_data(kD(k)%x,kD(k)%y,KD(k)%z)=phi_k_A(k)
# endif /* Dim2 */
         enddo
       call fftw_mpi_execute_dft(bplan, local_data, local_data)
        do k=0,LOCAL_SIZE-1
# if defined (Dim2)
          WB_new(k)=real(local_data(kD(k)%x,kD(k)%y))*M_grid_inv
# else  /* Dim2 */
          WB_new(k)=real(local_data(kD(k)%x,kD(k)%y,KD(k)%z))*M_grid_inv
# endif /* Dim2 */
        enddo

           error=0.0d0
        do k=0,LOCAL_SIZE-1
          error=error+(WA_new(k)-NXab*(RB(k)-fB)-(0.5*(WA(k)+WB(k))))
          error=error+(WB_new(k)-NXab*(RA(k)-fA)-(0.5*(WA(k)+WB(k))))
        enddo

      write(*,*) "error in iterate",error
      write(*,*) "sum of WA_new",sum(WA_new),"on",myid
      

      deallocate(phi_k_A)
      deallocate(phi_k_B)
      deallocate(yita)


   end subroutine fin_simp_iter


