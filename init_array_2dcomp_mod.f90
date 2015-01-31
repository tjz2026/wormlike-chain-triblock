# if defined (SLAB)
# else /* SLAB */

# if defined (SPH)

# else /* SPH */
subroutine init_arrays_2decomp()
 USE nrtype,only :DP
 USE global_para
 USE mpi
 USE control
 USE constants
 USE utility
 USE mmpi
 USE SPH
 use decomp_fft_mod
 implicit none
 integer :: istat2,i,j,k
 REAL(DP) :: temp,temp1,temp2,sumTHETA
 integer :: LS,MB,MA
 integer :: n_anderson
 integer :: nonzero_counter

 LS=LOCAL_SIZE-1
 max_M_Bar=max(M_Bar_A,M_Bar_B)

 MB=M_Bar_B-1
 MA=M_Bar_A-1
 if(myid==0) then
 write(*,*) "MA=,MB=",MA,MB
 endif
 call mp_barrier()

 !*****************note that N_dim_ddm==3!
  
 allocate(matrix_Rx(0:max_M_Bar-1,0:max_M_Bar-1),stat=istat2)
 allocate(matrix_Ry(0:max_M_Bar-1,0:max_M_Bar-1),stat=istat2)
 allocate(matrix_Rz(0:max_M_Bar-1,0:max_M_Bar-1),stat=istat2)

 allocate(GAMA(0:max_M_Bar-1,0:max_M_Bar-1,0:max_M_Bar-1),stat=istat2)
! for explicit method only
# if defined (EXP)
 allocate(qA(0:MA,0:LS,0:NA),stat=istat2)
 allocate(qB(0:MB,0:LS,0:NB),stat=istat2)
 allocate(qAstar(0:MA,0:LS,0:NA),stat=istat2)
 allocate(qBstar(0:MB,0:LS,0:NB),stat=istat2)
# else /* EXP */
 allocate(qA(0:NA,0:MA,0:LS),stat=istat2)
 allocate(qB(0:NB,0:MB,0:LS),stat=istat2)
 allocate(qAstar(0:NA,0:MA,0:LS),stat=istat2)
 allocate(qBstar(0:NB,0:MB,0:LS),stat=istat2)
# endif /* EXP */

 allocate(wA(0:LS),stat=istat2)
 allocate(wB(0:LS),stat=istat2)

!modified for SCFTread mode
 allocate(wA_save(0:LS,1:20),stat=istat2)
 allocate(wB_save(0:LS,1:20),stat=istat2)
 wA_save=0.0
 wB_save=0.0
!

 allocate(wA_out(0:LS,0:n_r_WAB),stat=istat2)
 allocate(wB_out(0:LS,0:n_r_WAB),stat=istat2)
 allocate(dw_A(0:LS,0:n_r_WAB),stat=istat2)
 allocate(dw_B(0:LS,0:n_r_WAB),stat=istat2)

 allocate(M11_out(0:LS,0:n_r_M),stat=istat2)
 allocate(M33_out(0:LS,0:n_r_M),stat=istat2)
 allocate(M12_out(0:LS,0:n_r_M),stat=istat2)
 allocate(M13_out(0:LS,0:n_r_M),stat=istat2)
 allocate(M23_out(0:LS,0:n_r_M),stat=istat2)
 allocate(d11_anderson(0:LS,0:n_r_M),stat=istat2)
 allocate(d12_anderson(0:LS,0:n_r_M),stat=istat2)
 allocate(d13_anderson(0:LS,0:n_r_M),stat=istat2)
 allocate(d23_anderson(0:LS,0:n_r_M),stat=istat2)
 allocate(d33_anderson(0:LS,0:n_r_M),stat=istat2)


 allocate(RA(0:LS),stat=istat2)
 allocate(RB(0:LS),stat=istat2)
 allocate(R_half(0:LS),stat=istat2)
 allocate(R_end(0:LS),stat=istat2)

 allocate(RHOA(0:MA,0:LS),stat=istat2)
 allocate(RHOB(0:MB,0:LS),stat=istat2)

 allocate(M_OP(0:LS,0:N_dim_ddm-1,0:N_dim_ddm-1),stat=istat2)
 allocate(S_OP(0:LS,0:N_dim_ddm-1,0:N_dim_ddm-1),stat=istat2)
 allocate(SA_OP(0:LS,0:N_dim_ddm-1,0:N_dim_ddm-1),stat=istat2)
 allocate(SB_OP(0:LS,0:N_dim_ddm-1,0:N_dim_ddm-1),stat=istat2)


 allocate(Pu(0:LS,1:N_dim_ddm),stat=istat2)
 allocate(Pu_A(0:LS,1:N_dim_ddm),stat=istat2)
 allocate(Pu_B(0:LS,1:N_dim_ddm),stat=istat2)

 allocate(THETAij_A(0:MA,0:MA,0:LS),stat=istat2)
 allocate(THETAij_M11_M22_A(0:MA,0:MA),stat=istat2)
 allocate(THETAij_M33_A(0:MA,0:MA),stat=istat2)
 allocate(THETAij_M12_A(0:MA,0:MA),stat=istat2)
 allocate(THETAij_M13_A(0:MA,0:MA),stat=istat2)
 allocate(THETAij_M23_A(0:MA,0:MA),stat=istat2)

 allocate(J11ij_A(0:MA,0:MA),stat=istat2)
 allocate(J22ij_A(0:MA,0:MA),stat=istat2)
 allocate(J12ij_A(0:MA,0:MA),stat=istat2)
 allocate(J13ij_A(0:MA,0:MA),stat=istat2)
 allocate(J23ij_A(0:MA,0:MA),stat=istat2)

!! note that THETAij_B could be optimized for memeory ,this is soon to de done
 allocate(THETAij_B(0:MB,0:MB,0:LS),stat=istat2)
 allocate(THETAij_M11_M22_B(0:MB,0:MB),stat=istat2)
 allocate(THETAij_M33_B(0:MB,0:MB),stat=istat2)
 allocate(THETAij_M12_B(0:MB,0:MB),stat=istat2)
 allocate(THETAij_M13_B(0:MB,0:MB),stat=istat2)
 allocate(THETAij_M23_B(0:MB,0:MB),stat=istat2)

 allocate(J11ij_B(0:MB,0:MB),stat=istat2)
 allocate(J22ij_B(0:MB,0:MB),stat=istat2)
 allocate(J12ij_B(0:MB,0:MB),stat=istat2)
 allocate(J13ij_B(0:MB,0:MB),stat=istat2)
 allocate(J23ij_B(0:MB,0:MB),stat=istat2)
 


allocate(THETA_nonzero_1D_A_indexK(-1:LOCAL_SIZE-1),stat=istat2)
allocate(THETA_nonzero_1D_B_indexK(-1:LOCAL_SIZE-1),stat=istat2)

call mp_barrier()
if(istat2/=0) then
write(*,*) "allocate failed,exit!","on",myid
stop
else
write(*,*) "allocate arrays succeeded","on",myid

endif 
call mp_barrier()
 

call basis_SPH_create()
call mp_barrier()
if(myid==0) then
write(*,*) "done basis SPH create"
endif
call mp_barrier()

!!! intialize for GAMA
do k=0,max_M_Bar-1
  do j=0,max_M_Bar-1
     do i=0,max_M_Bar-1
          GAMA(i,j,k)=Triple_product(basis_SPH(i)%x,basis_SPH(j)%x,basis_SPH(k)%x, &
                      basis_SPH(i)%y,basis_SPH(j)%y,basis_SPH(k)%y)
          enddo
      enddo
enddo

call mp_barrier()
if(myid==0) then
write(*,*) "done GAMA create"
endif
call mp_barrier()


temp=-1.0/sqrt(3.0)
  do j=0,max_M_Bar-1
     do i=0,max_M_Bar-1
     matrix_Rx(i,j)=temp*Triple_product(basis_SPH(i)%x,basis_SPH(j)%x,1,basis_SPH(i)%y, &
                    basis_SPH(j)%y,1)
     matrix_Ry(i,j)=temp*Triple_product(basis_SPH(i)%x,basis_SPH(j)%x,1,basis_SPH(i)%y, &
                    basis_SPH(j)%y,-1)
    enddo
  enddo

temp=1.0/sqrt(3.0)
  do j=0,max_M_Bar-1
     do i=0,max_M_Bar-1
     matrix_Rz(i,j)=temp*Triple_product(basis_SPH(i)%x,basis_SPH(j)%x,1,basis_SPH(i)%y, &
                    basis_SPH(j)%y,0)
    enddo
enddo

!!! THETAij_A 

temp=1.0/sqrt(15.0)
  do j=0,M_Bar_A-1
     do i=0,M_BAr_A-1
        THETAij_M11_M22_A(i,j)=temp*Triple_product(basis_SPH(i)%x,basis_SPH(j)%x, &
                                2,basis_SPH(i)%y,basis_SPH(j)%y,2)
     enddo
  enddo

temp=1.0/sqrt(5.0)
  do j=0,M_Bar_A-1
     do i=0,M_BAr_A-1
        THETAij_M33_A(i,j)=temp*Triple_product(basis_SPH(i)%x,basis_SPH(j)%x, &
                                2,basis_SPH(i)%y,basis_SPH(j)%y,0)
     enddo
  enddo

temp=2.0/sqrt(15.0)
  do j=0,M_Bar_A-1
     do i=0,M_BAr_A-1
        THETAij_M12_A(i,j)=temp*Triple_product(basis_SPH(i)%x,basis_SPH(j)%x, &
                                2,basis_SPH(i)%y,basis_SPH(j)%y,-2)
     enddo
  enddo
 

temp=-2.0/sqrt(15.0)
  do j=0,M_Bar_A-1
     do i=0,M_BAr_A-1
        THETAij_M13_A(i,j)=temp*Triple_product(basis_SPH(i)%x,basis_SPH(j)%x, &
                                2,basis_SPH(i)%y,basis_SPH(j)%y,1)
     enddo
  enddo

temp=-2.0/sqrt(15.0)
  do j=0,M_Bar_A-1
     do i=0,M_BAr_A-1
        THETAij_M23_A(i,j)=temp*Triple_product(basis_SPH(i)%x,basis_SPH(j)%x, &
                                2,basis_SPH(i)%y,basis_SPH(j)%y,-1)
     enddo
  enddo


!!!for THETAij_B

temp=1.0/sqrt(15.0)
  do j=0,M_Bar_B-1
     do i=0,M_BAr_B-1
        THETAij_M11_M22_B(i,j)=temp*Triple_product(basis_SPH(i)%x,basis_SPH(j)%x, &
                                2,basis_SPH(i)%y,basis_SPH(j)%y,2)
     enddo
  enddo

temp=1.0/sqrt(5.0)
  do j=0,M_Bar_B-1
     do i=0,M_BAr_B-1
        THETAij_M33_B(i,j)=temp*Triple_product(basis_SPH(i)%x,basis_SPH(j)%x, &
                                2,basis_SPH(i)%y,basis_SPH(j)%y,0)
     enddo
  enddo
                               

temp=2.0/sqrt(15.0)
  do j=0,M_Bar_B-1
     do i=0,M_BAr_B-1
        THETAij_M12_B(i,j)=temp*Triple_product(basis_SPH(i)%x,basis_SPH(j)%x, &
                                2,basis_SPH(i)%y,basis_SPH(j)%y,-2)
     enddo
  enddo
 

temp=-2.0/sqrt(15.0)
  do j=0,M_Bar_B-1
     do i=0,M_BAr_B-1
        THETAij_M13_B(i,j)=temp*Triple_product(basis_SPH(i)%x,basis_SPH(j)%x, &
                                2,basis_SPH(i)%y,basis_SPH(j)%y,1)
     enddo
  enddo

temp=-2.0/sqrt(15.0)
  do j=0,M_Bar_B-1
     do i=0,M_BAr_B-1
        THETAij_M23_B(i,j)=temp*Triple_product(basis_SPH(i)%x,basis_SPH(j)%x, &
                                2,basis_SPH(i)%y,basis_SPH(j)%y,-1)
     enddo
  enddo




!!!J for A
temp1=1.0/sqrt(15.0)
temp2=-1.0/(3.0*sqrt(5.0))
   do j=0,M_Bar_A-1
      do i=0,M_Bar_A-1
         J11ij_A(i,j)=temp1*triple_product(basis_SPH(i)%x,basis_SPH(j)%x,2, &
                           basis_SPH(i)%y,basis_SPH(j)%y,2) + &
                      temp2*triple_product(basis_SPH(i)%x,basis_SPH(j)%x,2, &
                           basis_SPH(i)%y,basis_SPH(j)%y,0)
      enddo
   enddo


temp1=-1.0/sqrt(15.0)
temp2=-1.0/(3.0*sqrt(5.0))
   do j=0,M_Bar_A-1
      do i=0,M_Bar_A-1
         J22ij_A(i,j)=temp1*triple_product(basis_SPH(i)%x,basis_SPH(j)%x,2, &
                           basis_SPH(i)%y,basis_SPH(j)%y,2) + &
                      temp2*triple_product(basis_SPH(i)%x,basis_SPH(j)%x,2, &
                           basis_SPH(i)%y,basis_SPH(j)%y,0)
      enddo
   enddo

temp1=1.0/sqrt(15.0)
   do j=0,M_Bar_A-1
      do i=0,M_Bar_A-1
         J12ij_A(i,j)=temp1*triple_product(basis_SPH(i)%x,basis_SPH(j)%x,2, &
                           basis_SPH(i)%y,basis_SPH(j)%y,-2)
      enddo
   enddo


temp1=-1.0/sqrt(15.0)
   do j=0,M_Bar_A-1
      do i=0,M_Bar_A-1
         J13ij_A(i,j)=temp1*triple_product(basis_SPH(i)%x,basis_SPH(j)%x,2, &
                           basis_SPH(i)%y,basis_SPH(j)%y,1)
      enddo
   enddo


temp1=-1.0/sqrt(15.0)
   do j=0,M_Bar_A-1
      do i=0,M_Bar_A-1
         J23ij_A(i,j)=temp1*triple_product(basis_SPH(i)%x,basis_SPH(j)%x,2, &
                           basis_SPH(i)%y,basis_SPH(j)%y,-1)
      enddo
   enddo



!!!J for B
temp1=1.0/sqrt(15.0)
temp2=-1.0/(3.0*sqrt(5.0))
   do j=0,M_Bar_B-1
      do i=0,M_Bar_B-1
         J11ij_B(i,j)=temp1*triple_product(basis_SPH(i)%x,basis_SPH(j)%x,2, &
                           basis_SPH(i)%y,basis_SPH(j)%y,2) + &
                      temp2*triple_product(basis_SPH(i)%x,basis_SPH(j)%x,2, &
                           basis_SPH(i)%y,basis_SPH(j)%y,0)
      enddo
   enddo


temp1=-1.0/sqrt(15.0)
temp2=-1.0/(3.0*sqrt(5.0))
   do j=0,M_Bar_B-1
      do i=0,M_Bar_B-1
         J22ij_B(i,j)=temp1*triple_product(basis_SPH(i)%x,basis_SPH(j)%x,2, &
                           basis_SPH(i)%y,basis_SPH(j)%y,2) + &
                      temp2*triple_product(basis_SPH(i)%x,basis_SPH(j)%x,2, &
                           basis_SPH(i)%y,basis_SPH(j)%y,0)
      enddo
   enddo

temp1=1.0/sqrt(15.0)
   do j=0,M_Bar_B-1
      do i=0,M_Bar_B-1
         J12ij_B(i,j)=temp1*triple_product(basis_SPH(i)%x,basis_SPH(j)%x,2, &
                           basis_SPH(i)%y,basis_SPH(j)%y,-2)
      enddo
   enddo


temp1=-1.0/sqrt(15.0)
   do j=0,M_Bar_B-1
      do i=0,M_Bar_B-1
         J13ij_B(i,j)=temp1*triple_product(basis_SPH(i)%x,basis_SPH(j)%x,2, &
                           basis_SPH(i)%y,basis_SPH(j)%y,1)
      enddo
   enddo


temp1=-1.0/sqrt(15.0)
   do j=0,M_Bar_B-1
      do i=0,M_Bar_B-1
         J23ij_B(i,j)=temp1*triple_product(basis_SPH(i)%x,basis_SPH(j)%x,2, &
                           basis_SPH(i)%y,basis_SPH(j)%y,-1)
      enddo
   enddo


!note that we will skip the intial of THETA_nonzero_2D_A because
!it will not be allocated before we call the cal_theta() subroutine 
!also,the cal_theta must be called before every iteration of M

!note that the following variable stays constant once initialized
!init Rx_nonzero_1D

 nonzero_counter=0
 do j=0,max_M_Bar-1
    do i=0,max_M_Bar-1
          if(abs(matrix_Rx(i,j))>=Thresh_sprase_matrix) then
              nonzero_counter=nonzero_counter+1
          endif
     enddo
  enddo

# if defined (Debug)
 if(myid==0) then
 write(*,*) "nonzero elements of Rx=",nonzero_counter
 endif
# endif  /* Debug */
       
  if(.not. allocated(Rx_nonzero_1D)) then
  allocate(Rx_nonzero_1D(0:nonzero_counter-1))
  endif
  
 nonzero_counter=0
 do j=0,max_M_Bar-1
    do i=0,max_M_Bar-1
          if(abs(matrix_Rx(i,j))>=Thresh_sprase_matrix) then
              Rx_nonzero_1D(nonzero_counter)%i=i
              Rx_nonzero_1D(nonzero_counter)%j=j
              Rx_nonzero_1D(nonzero_counter)%value=matrix_Rx(i,j)
              nonzero_counter=nonzero_counter+1
           endif
      enddo
  enddo

        if(size(Rx_nonzero_1D)/=nonzero_counter) then
        write(*,*) "counter wrong in Rx_nonzero_1D"
        stop
        endif
        
 Rx_num=nonzero_counter
          

 nonzero_counter=0
 do j=0,max_M_Bar-1
    do i=0,max_M_Bar-1
          if(abs(matrix_Ry(i,j))>=Thresh_sprase_matrix) then
              nonzero_counter=nonzero_counter+1
          endif
     enddo
  enddo
       
# if defined (Debug)
 if(myid==0) then
 write(*,*) "nonzero elements of Ry=",nonzero_counter
 endif
# endif  /* Debug */

  if(.not. allocated(Ry_nonzero_1D)) then
  allocate(Ry_nonzero_1D(0:nonzero_counter-1))
  endif
  
 nonzero_counter=0
 do j=0,max_M_Bar-1
    do i=0,max_M_Bar-1
          if(abs(matrix_Ry(i,j))>=Thresh_sprase_matrix) then
              Ry_nonzero_1D(nonzero_counter)%i=i
              Ry_nonzero_1D(nonzero_counter)%j=j
              Ry_nonzero_1D(nonzero_counter)%value=matrix_Ry(i,j)
              nonzero_counter=nonzero_counter+1
           endif
      enddo
  enddo

        if(size(Ry_nonzero_1D)/=nonzero_counter) then
        write(*,*) "counter wrong in Ry_nonzero_1D"
        stop
        endif
        
 Ry_num=nonzero_counter
          

 nonzero_counter=0
 do j=0,max_M_Bar-1
    do i=0,max_M_Bar-1
          if(abs(matrix_Rz(i,j))>=Thresh_sprase_matrix) then
              nonzero_counter=nonzero_counter+1
          endif
     enddo
  enddo
       
# if defined (Debug)
 if(myid==0) then
 write(*,*) "nonzero elements of Rz=",nonzero_counter
 endif
# endif  /* Debug */

  if(.not. allocated(Rz_nonzero_1D)) then
  allocate(Rz_nonzero_1D(0:nonzero_counter-1))
  endif
  
 nonzero_counter=0
 do j=0,max_M_Bar-1
    do i=0,max_M_Bar-1
          if(abs(matrix_Rz(i,j))>=Thresh_sprase_matrix) then
              Rz_nonzero_1D(nonzero_counter)%i=i
              Rz_nonzero_1D(nonzero_counter)%j=j
              Rz_nonzero_1D(nonzero_counter)%value=matrix_Rz(i,j)
              nonzero_counter=nonzero_counter+1
           endif
      enddo
  enddo

        if(size(Rz_nonzero_1D)/=nonzero_counter) then
        write(*,*) "counter wrong in Rx_nonzero_1D"
        stop
        endif
        
 Rz_num=nonzero_counter


          
!!initialize J??ij_nonzero_1D_A
 nonzero_counter=0
 do j=0,M_Bar_A-1
    do i=0,M_bar_A-1
          if(abs(J11ij_A(i,j))>=Thresh_sprase_matrix) then
              nonzero_counter=nonzero_counter+1
          endif
     enddo
  enddo
       
# if defined (Debug)
 if(myid==0) then
 write(*,*) "nonzero elements of J11ij_A=",nonzero_counter
 endif
# endif  /* Debug */

  if(.not. allocated(J11ij_nonzero_1D_A)) then
  allocate(J11ij_nonzero_1D_A(0:nonzero_counter-1))
  endif
  
 nonzero_counter=0
 do j=0,M_Bar_A-1
    do i=0,M_bar_A-1
          if(abs(J11ij_A(i,j))>=Thresh_sprase_matrix) then
              J11ij_nonzero_1D_A(nonzero_counter)%i=i
              J11ij_nonzero_1D_A(nonzero_counter)%j=j
              J11ij_nonzero_1D_A(nonzero_counter)%value=J11ij_A(i,j)
              nonzero_counter=nonzero_counter+1
           endif
      enddo
  enddo

        if(size(J11ij_nonzero_1D_A)/=nonzero_counter) then
        write(*,*) "counter wrong in J11ij_nonzero_1D_A"
        stop
        endif
        
 J11_A=nonzero_counter
          
 nonzero_counter=0
 do j=0,M_Bar_A-1
    do i=0,M_bar_A-1
          if(abs(J22ij_A(i,j))>=Thresh_sprase_matrix) then
              nonzero_counter=nonzero_counter+1
          endif
     enddo
  enddo
       
  if(.not. allocated(J22ij_nonzero_1D_A)) then
  allocate(J22ij_nonzero_1D_A(0:nonzero_counter-1))
  endif
  
# if defined (Debug)
 if(myid==0) then
 write(*,*) "nonzero elements of J22ij_A=",nonzero_counter
 endif
# endif  /* Debug */

 nonzero_counter=0
 do j=0,M_Bar_A-1
    do i=0,M_bar_A-1
          if(abs(J22ij_A(i,j))>=Thresh_sprase_matrix) then
              J22ij_nonzero_1D_A(nonzero_counter)%i=i
              J22ij_nonzero_1D_A(nonzero_counter)%j=j
              J22ij_nonzero_1D_A(nonzero_counter)%value=J22ij_A(i,j)
              nonzero_counter=nonzero_counter+1
           endif
      enddo
  enddo

 J22_A=nonzero_counter
              
 nonzero_counter=0
 do j=0,M_Bar_A-1
    do i=0,M_bar_A-1
          if(abs(J12ij_A(i,j))>=Thresh_sprase_matrix) then
              nonzero_counter=nonzero_counter+1
          endif
     enddo
  enddo
       
  if(.not. allocated(J12ij_nonzero_1D_A)) then
  allocate(J12ij_nonzero_1D_A(0:nonzero_counter-1))
  endif
  
 J12_A=nonzero_counter
 nonzero_counter=0
 do j=0,M_Bar_A-1
    do i=0,M_bar_A-1
          if(abs(J12ij_A(i,j))>=Thresh_sprase_matrix) then
              J12ij_nonzero_1D_A(nonzero_counter)%i=i
              J12ij_nonzero_1D_A(nonzero_counter)%j=j
              J12ij_nonzero_1D_A(nonzero_counter)%value=J12ij_A(i,j)
              nonzero_counter=nonzero_counter+1
           endif
      enddo
  enddo


 nonzero_counter=0
 do j=0,M_Bar_A-1
    do i=0,M_bar_A-1
          if(abs(J13ij_A(i,j))>=Thresh_sprase_matrix) then
              nonzero_counter=nonzero_counter+1
          endif
     enddo
  enddo
       
  if(.not. allocated(J13ij_nonzero_1D_A)) then
  allocate(J13ij_nonzero_1D_A(0:nonzero_counter-1))
  endif
  
 nonzero_counter=0
 do j=0,M_Bar_A-1
    do i=0,M_bar_A-1
          if(abs(J13ij_A(i,j))>=Thresh_sprase_matrix) then
              J13ij_nonzero_1D_A(nonzero_counter)%i=i
              J13ij_nonzero_1D_A(nonzero_counter)%j=j
              J13ij_nonzero_1D_A(nonzero_counter)%value=J13ij_A(i,j)
              nonzero_counter=nonzero_counter+1
           endif
      enddo
  enddo


 J13_A=nonzero_counter
 nonzero_counter=0
 do j=0,M_Bar_A-1
    do i=0,M_bar_A-1
          if(abs(J23ij_A(i,j))>=Thresh_sprase_matrix) then
              nonzero_counter=nonzero_counter+1
          endif
     enddo
  enddo
       
  if(.not. allocated(J23ij_nonzero_1D_A)) then
  allocate(J23ij_nonzero_1D_A(0:nonzero_counter-1))
  endif
  
 nonzero_counter=0
 do j=0,M_Bar_A-1
    do i=0,M_bar_A-1
          if(abs(J23ij_A(i,j))>=Thresh_sprase_matrix) then
              J23ij_nonzero_1D_A(nonzero_counter)%i=i
              J23ij_nonzero_1D_A(nonzero_counter)%j=j
              J23ij_nonzero_1D_A(nonzero_counter)%value=J23ij_A(i,j)
              nonzero_counter=nonzero_counter+1
           endif
      enddo
  enddo


 J23_A=nonzero_counter

!!initialize J??ij_nonzero_1D_B
 nonzero_counter=0
 do j=0,M_Bar_B-1
    do i=0,M_bar_B-1
          if(abs(J11ij_B(i,j))>=Thresh_sprase_matrix) then
              nonzero_counter=nonzero_counter+1
          endif
     enddo
  enddo
       
  if(.not. allocated(J11ij_nonzero_1D_B)) then
  allocate(J11ij_nonzero_1D_B(0:nonzero_counter-1))
  endif
  
 nonzero_counter=0
 do j=0,M_Bar_B-1
    do i=0,M_bar_B-1
          if(abs(J11ij_B(i,j))>=Thresh_sprase_matrix) then
              J11ij_nonzero_1D_B(nonzero_counter)%i=i
              J11ij_nonzero_1D_B(nonzero_counter)%j=j
              J11ij_nonzero_1D_B(nonzero_counter)%value=J11ij_B(i,j)
              nonzero_counter=nonzero_counter+1
           endif
      enddo
  enddo

        if(size(J11ij_nonzero_1D_B)/=nonzero_counter) then
        write(*,*) "counter wrong in J11ij_nonzero_1D_B"
        stop
        endif
        
           
 J11_B=nonzero_counter

 nonzero_counter=0
 do j=0,M_Bar_B-1
    do i=0,M_bar_B-1
          if(abs(J22ij_B(i,j))>=Thresh_sprase_matrix) then
              nonzero_counter=nonzero_counter+1
          endif
     enddo
  enddo
       
  if(.not. allocated(J22ij_nonzero_1D_B)) then
  allocate(J22ij_nonzero_1D_B(0:nonzero_counter-1))
  endif
  
 nonzero_counter=0
 do j=0,M_Bar_B-1
    do i=0,M_bar_B-1
          if(abs(J22ij_B(i,j))>=Thresh_sprase_matrix) then
              J22ij_nonzero_1D_B(nonzero_counter)%i=i
              J22ij_nonzero_1D_B(nonzero_counter)%j=j
              J22ij_nonzero_1D_B(nonzero_counter)%value=J22ij_B(i,j)
              nonzero_counter=nonzero_counter+1
           endif
      enddo
  enddo

              
 J22_B=nonzero_counter
 nonzero_counter=0
 do j=0,M_Bar_B-1
    do i=0,M_bar_B-1
          if(abs(J12ij_B(i,j))>=Thresh_sprase_matrix) then
              nonzero_counter=nonzero_counter+1
          endif
     enddo
  enddo
       
  if(.not. allocated(J12ij_nonzero_1D_B)) then
  allocate(J12ij_nonzero_1D_B(0:nonzero_counter-1))
  endif
  
 nonzero_counter=0
 do j=0,M_Bar_B-1
    do i=0,M_bar_B-1
          if(abs(J12ij_B(i,j))>=Thresh_sprase_matrix) then
              J12ij_nonzero_1D_B(nonzero_counter)%i=i
              J12ij_nonzero_1D_B(nonzero_counter)%j=j
              J12ij_nonzero_1D_B(nonzero_counter)%value=J12ij_B(i,j)
              nonzero_counter=nonzero_counter+1
           endif
      enddo
  enddo

 J12_B=nonzero_counter

 nonzero_counter=0
 do j=0,M_Bar_B-1
    do i=0,M_bar_B-1
          if(abs(J13ij_B(i,j))>=Thresh_sprase_matrix) then
              nonzero_counter=nonzero_counter+1
          endif
     enddo
  enddo
       
  if(.not. allocated(J13ij_nonzero_1D_B)) then
  allocate(J13ij_nonzero_1D_B(0:nonzero_counter-1))
  endif
  
 nonzero_counter=0
 do j=0,M_Bar_B-1
    do i=0,M_bar_B-1
          if(abs(J13ij_B(i,j))>=Thresh_sprase_matrix) then
              J13ij_nonzero_1D_B(nonzero_counter)%i=i
              J13ij_nonzero_1D_B(nonzero_counter)%j=j
              J13ij_nonzero_1D_B(nonzero_counter)%value=J13ij_B(i,j)
              nonzero_counter=nonzero_counter+1
           endif
      enddo
  enddo

 J13_B=nonzero_counter

 nonzero_counter=0
 do j=0,M_Bar_B-1
    do i=0,M_bar_B-1
          if(abs(J23ij_B(i,j))>=Thresh_sprase_matrix) then
              nonzero_counter=nonzero_counter+1
          endif
     enddo
  enddo
       
  if(.not. allocated(J23ij_nonzero_1D_B)) then
  allocate(J23ij_nonzero_1D_B(0:nonzero_counter-1))
  endif
  
 nonzero_counter=0
 do j=0,M_Bar_B-1
    do i=0,M_bar_B-1
          if(abs(J23ij_B(i,j))>=Thresh_sprase_matrix) then
              J23ij_nonzero_1D_B(nonzero_counter)%i=i
              J23ij_nonzero_1D_B(nonzero_counter)%j=j
              J23ij_nonzero_1D_B(nonzero_counter)%value=J23ij_B(i,j)
              nonzero_counter=nonzero_counter+1
           endif
      enddo
  enddo

 J23_B=nonzero_counter

 nonzero_counter=0
do k=0,max_M_Bar-1
  do j=0,max_M_Bar-1
     do i=0,max_M_Bar-1
          if(abs(GAMA(i,j,k))>=Thresh_sprase_matrix) then
              nonzero_counter=nonzero_counter+1
          endif
     enddo
  enddo
enddo
       
# if defined (Debug)
 if(myid==0) then
 write(*,*) "nonzero elements of GAMA=",nonzero_counter
 endif
# endif  /* Debug */

  if(.not. allocated(GAMA_nonzero_1D)) then
  allocate(GAMA_nonzero_1D(0:nonzero_counter-1))
  endif
  
 nonzero_counter=0
do k=0,max_M_Bar-1
  do j=0,max_M_Bar-1
     do i=0,max_M_Bar-1
          if(abs(GAMA(i,j,k))>=Thresh_sprase_matrix) then
              GAMA_nonzero_1D(nonzero_counter)%i=i
              GAMA_nonzero_1D(nonzero_counter)%j=j
              GAMA_nonzero_1D(nonzero_counter)%k=k
              GAMA_nonzero_1D(nonzero_counter)%value=GAMA(i,j,k)
              nonzero_counter=nonzero_counter+1
           endif
      enddo
  enddo
enddo



end subroutine init_arrays_2decomp


          
!!!!!!!!!!!!!!!!!!!initialize the input field!!!!!!!!!!!!!!!!!!

         subroutine init_w()
         USE nrtype,only :PI
         USE global_para
         USE mpi
         USE control
         USE constants
         USE utility
         USE mmpi
         use decomp_fft_mod
         implicit none
         integer :: k_i,k_j,k_k
         integer :: KK,j,k,i
         integer :: initseed
         integer :: dt(8)
         
          call date_and_time(values=dt)
          initseed=(dt(8)-500)*54321*(myid+1)+88888
   if(wAB_init_type==0) then
          !random init
          call date_and_time(values=dt)
          initseed=(dt(8)-500)*54321*(myid+1)+88888
          
            KK=0  
                do KK=0,LOCAL_SIZE-1
                   WA(KK)=1.0*ran2(initseed)
                   WB(KK)=1.0*ran2(initseed)
               enddo
        
          KK=0   
            do k=0,N_dim_ddm-1
              do j=0,N_dim_ddm-1
                do KK=0,LOCAL_SIZE-1
                  S_OP(KK,j,k)=0.0
                enddo
              enddo
             enddo  
          
       if(method=='A') then
           !Smectic A
              do k=0,LOCAL_SIZE-1
                S_OP(k,0,0)=1.0*ran2(initseed)
                S_OP(k,1,1)=1.0*ran2(initseed)
                S_OP(k,2,2)=1.0*ran2(initseed)
              enddo
        else if(method=='B') then
            !Smectic C
              do k=0,LOCAL_SIZE-1
                S_OP(k,2,2)=1.0*ran2(initseed)
                S_OP(k,0,1)=1.0*ran2(initseed)
                S_OP(k,0,2)=1.0*ran2(initseed)
                S_OP(k,1,2)=1.0*ran2(initseed)
                S_OP(k,0,0)=-0.5*S_OP(k,2,2)
                S_OP(k,1,1)=S_OP(k,0,0)
                S_OP(k,1,0)=S_OP(k,0,1)
                S_OP(k,2,0)=S_OP(k,0,2)
                S_OP(k,2,1)=S_OP(k,1,2)
               enddo
        else 
         !biaxial phase
            do k=0,LOCAL_SIZE-1
                S_OP(k,0,0)=1.0*ran2(initseed)
                S_OP(k,2,2)=1.0*ran2(initseed)
                S_OP(k,0,1)=1.0*ran2(initseed)
                S_OP(k,0,2)=1.0*ran2(initseed)
                S_OP(k,1,2)=1.0*ran2(initseed)

                S_OP(k,1,1)=-S_OP(k,0,0)-S_OP(k,2,2)
                S_OP(k,1,0)=S_OP(k,0,1)
                S_OP(k,2,0)=S_OP(k,0,2)
                S_OP(k,2,1)=S_OP(k,1,2)
             enddo
         endif

        do i=0,N_dim_ddm-1  
          do j=0,N_dim_ddm-1
            do k=0,LOCAL_SIZE-1
               M_OP(k,j,i)=M_initial*S_OP(k,j,i)
            enddo
           enddo
         enddo

# if defined (Debug)
if(myid==0) then
    write(*,*) "random init done!","on",myid
    write(*,*) "S_OP()",S_OP(0,0,0),S_OP(1,1,2),"on",myid
endif
# endif  /* Debug */


   else if(wAB_init_type==1) then   
      !!!FCC initial 3D
         do KK=0,local_size-1
           WA(KK)=NXab*(1-fA*(1+0.7*(cos(2.0*PI*(kD(KK)%x+1)/SIDEx)* &
           cos(2.0*PI*(kD(KK)%y+1)/SIDEy)*cos(2.0*PI*(kD(KK)%z+1)/SIDEz))))
                   
           WB(KK)=NXab*fA*(1+0.7*(cos(2.0*PI*(kD(KK)%x+1)/SIDEx)* &
                         cos(2.0*PI*(kD(KK)%y+1)/SIDEy)*cos(2.0*PI*(kD(KK)%z+1)/SIDEz)))
          enddo
   
                
    !!! BCC initial 3D
   else if(wAB_init_type==2) then
                 do KK=0,local_size-1
        WA(KK)=NXab*(1-fA*(1+0.7*(cos(2.0*PI*(kD(KK)%x+1)/SIDEx)* &
        cos(2.0*PI*(kD(KK)%y+1)/SIDEy) + & 
        cos(2.0*PI*(kD(KK)%y+1)/SIDEy)*cos(2.0*PI*(kD(KK)%z+1)/SIDEz) + &
        cos(2.0*PI*(kD(KK)%x+1)/SIDEx)*cos(2.0*PI*(kD(KK)%z+1)/SIDEz))))

         WB(KK)=NXab*fA*(1+0.7*(cos(2.0*PI*(kD(KK)%x+1)/SIDEx)* &
         cos(2.0*PI*(kD(KK)%y+1)/SIDEy) + & 
         cos(2.0*PI*(kD(KK)%y+1)/SIDEy)*cos(2.0*PI*(kD(KK)%z+1)/SIDEz) + &
         cos(2.0*PI*(kD(KK)%x+1)/SIDEx)*cos(2.0*PI*(kD(KK)%z+1)/SIDEz)))
                  enddo
   
                
      !Gyroid initial 3D
   else if(wAB_init_type==3) then
      
              do KK=0,local_size-1
                WA(KK)=NXab*(1-fA*(1+0.7*(cos(2.0*PI*(kD(KK)%x+1)/SIDEx)* &
                   sin(2.0*PI*(kD(KK)%y+1)/SIDEy)*sin(4.0*PI*(kD(KK)%z+1)/SIDEz) &
                   + cos(2.0*PI*(kD(KK)%y+1)/SIDEy)*sin(2.0*PI*(kD(KK)%z+1)/SIDEz)* &
                     sin(4.0*Pi*(kD(KK)%x+1)/SIDEx) &
                   + cos(2.0*PI*(kD(KK)%z+1)/SIDEz)*sin(2.0*Pi*(kD(KK)%x+1)/SIDEx)* &
                     sin(4.0*PI*(kD(KK)%y+1)/SIDEy))))

                WB(KK)=NXab*fA*(1+0.7*(cos(2.0*PI*(kD(KK)%x+1)/SIDEx)* &
                   sin(2.0*PI*(kD(KK)%y+1)/SIDEy)*sin(4.0*PI*(kD(KK)%z+1)/SIDEz) &
                   + cos(2.0*PI*(kD(KK)%y+1)/SIDEy)*sin(2.0*PI*(kD(KK)%z+1)/SIDEz)* &
                     sin(4.0*Pi*(kD(KK)%x+1)/SIDEx) &
                   + cos(2.0*PI*(kD(KK)%z+1)/SIDEz)*sin(2.0*Pi*(kD(KK)%x+1)/SIDEx)* &
                     sin(4.0*PI*(kD(KK)%y+1)/SIDEy)))

                  enddo

   else if(wAB_init_type==4) then
      !!!Cylinder initial 3D
                 do KK=0,local_size-1
                   WA(KK)=NXab*(1-fA*(1+0.8*cos(2.0*PI*(kD(KK)%x+1)/SIDEx)* &
                         cos(2.0*PI*(kD(KK)%y+1)/SIDEy)))

                   WB(KK)=NXab*fA*(1+0.8*cos(2.0*PI*(kD(KK)%x+1)/SIDEx)* &
                         cos(2.0*PI*(kD(KK)%y+1)/SIDEy))
                  enddo


   else if(wAB_init_type==5) then
      !!!Lamellar initial 3D
                 do KK=0,local_size-1
                   WA(KK)=NXab*(1-fA*(1+0.8*sin(2.0*PI*(kD(KK)%z+1)/SIDEz)))
                   WB(KK)=NXab*fA*(1+0.8*sin(2.0*PI*(kD(KK)%z+1)/SIDEz))
                  enddo


   else if(wAB_init_type==6) then
      !!!P4 initial 3D
                 do k_i=0,SIDEz-1
             WA(KK)=NXab*(1-fA*(cos(2.0*PI*(kD(KK)%y+1)/(1.0*SIDEy))+ & 
               cos(2.0*PI*(kD(KK)%x+1)/(1.0*SIDEx))))
             WB(KK)=NXab*fA*(cos(2.0*PI*(kD(KK)%y+1)/(1.0*SIDEy))+ & 
                 cos(2.0*PI*(kD(KK)%x+1)/(1.0*SIDEx)))
                  enddo

endif


!!!!!!!!!!!!!!!!!!!!!for M ,not randomly initialized


          do k=0,N_dim_ddm-1
             do j=0,N_dim_ddm-1
               do KK=0,LOCAL_SIZE-1
                  S_OP(KK,j,k)=0.05
                enddo
              enddo
            enddo  
          
          if(method=='A') then
           !Smectic A
              do k=0,LOCAL_SIZE-1
                S_OP(k,0,0)=-0.5
                S_OP(k,1,1)=-0.5
                S_OP(k,2,2)=1.0
              enddo
           else if(method=='B') then
            !Smectic C
              do k=0,LOCAL_SIZE-1
                S_OP(k,2,2)=1.0
                S_OP(k,0,1)=0.5
                S_OP(k,0,2)=0.5
                S_OP(k,1,2)=0.5
                S_OP(k,0,0)=-0.5*S_OP(k,2,2)
                S_OP(k,1,1)=S_OP(k,0,0)
                S_OP(k,1,0)=S_OP(k,0,1)
                S_OP(k,2,0)=S_OP(k,0,2)
                S_OP(k,2,1)=S_OP(k,1,2)
               enddo
          else 
         !biaxial phase
            do k=0,LOCAL_SIZE-1
             
                S_OP(k,0,0)=-0.5
                S_OP(k,2,2)=1.0
                S_OP(k,0,1)=0.5
                S_OP(k,0,2)=0.5
                S_OP(k,1,2)=0.5

                S_OP(k,1,1)=-S_OP(k,0,0)-S_OP(k,2,2)
                S_OP(k,1,0)=S_OP(k,0,1)
                S_OP(k,2,0)=S_OP(k,0,2)
                S_OP(k,2,1)=S_OP(k,1,2)
             enddo
            endif

        do i=0,N_dim_ddm-1  
          do j=0,N_dim_ddm-1
            do k=0,LOCAL_SIZE-1
               M_OP(k,j,i)=M_initial*S_OP(k,j,i)*(1.0+0.15*(ran2(seed)-0.5))
            enddo
           enddo
         enddo
     
   end subroutine init_w




  subroutine readW()
   USE nrtype,only :DP
 USE global_para
 USE mpi
 USE control
 USE constants
 USE utility
 USE mmpi
  implicit none
  integer :: i,j,k,aaa
  integer :: ii,iitot
  REAL*8 :: x,y,z,sumaW,sumaWtot
  REAL*8 :: sumM00,sumM01,sumM02,sumM10,sumM11,sumM12,sumM20,sumM21,sumM22,sumMtot
  character(len=30)::aa

  aaa=myid
  ii=myid
  write(aa,*) aaa
  open(unit=23,file= 'W' //trim(adjustl(aa)) // '.dat',status='old')
  do K=0,LOCAL_SIZE-1
  read(23,*) i,WA(K),WB(K)
  
  enddo
  close(23)

  open(unit=23,file= 'M_OP' //trim(adjustl(aa)) // '.dat',status='old')
  do K=0,LOCAL_SIZE-1
  read(23,*) x,y,z,M_OP(K,0,0),M_OP(K,0,1),M_OP(K,0,2)
  read(23,*) x,y,z,M_OP(K,1,0),M_OP(K,1,1),M_OP(K,1,2)
  read(23,*) x,y,z,M_OP(K,2,0),M_OP(K,2,1),M_OP(K,2,2)
  
  enddo
  close(23)



 end subroutine readW


 subroutine global_arrays_clean()
 USE nrtype,only :DP
 USE global_para
 USE mpi
 USE control
 USE constants
 USE utility
 USE mmpi
 implicit none
 deallocate(matrix_Rx)
 deallocate(matrix_Ry)
 deallocate(matrix_Rz)
 deallocate(GAMA)
 deallocate(qA)
 deallocate(qB)
 deallocate(qAstar)
 deallocate(qBstar)
 deallocate(wA)
 deallocate(wB)
 deallocate(RA)
 deallocate(RB)
 deallocate(R_half)
 deallocate(R_end)
 deallocate(RHOA)
 deallocate(RHOB)

 deallocate(M_OP)
 deallocate(S_OP)
 deallocate(SA_OP)
 deallocate(SB_OP)

 deallocate(THETAij_A)
 deallocate(THETAij_M11_M22_A)
 deallocate(THETAij_M33_A)
 deallocate(THETAij_M12_A)
 deallocate(THETAij_M13_A)
 deallocate(THETAij_M23_A)
 deallocate(J11ij_A)
 deallocate(J22ij_A)
 deallocate(J12ij_A)
 deallocate(J13ij_A)
 deallocate(J23ij_A)


 deallocate(THETAij_B)
 deallocate(THETAij_M11_M22_B)
 deallocate(THETAij_M33_B)
 deallocate(THETAij_M12_B)
 deallocate(THETAij_M13_B)
 deallocate(THETAij_M23_B)
 deallocate(J11ij_B)
 deallocate(J22ij_B)
 deallocate(J12ij_B)
 deallocate(J13ij_B)
 deallocate(J23ij_B)

 deallocate(dw_B)
 deallocate(dw_A)
 deallocate(wA_out)
 deallocate(wB_out)

deallocate(M33_out)
deallocate(M13_out)
deallocate(M23_out)
deallocate(M12_out)
deallocate(M11_out)

deallocate(d33_anderson)
deallocate(d13_anderson)
deallocate(d23_anderson)
deallocate(d12_anderson)
deallocate(d11_anderson)

end subroutine global_arrays_clean
# endif /* SPH */
# endif /* SLAB */
