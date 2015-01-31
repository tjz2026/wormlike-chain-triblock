# if defined (EXP)
Module wormlike_MDE_single_step
 USE nrtype,only :DP,PI
 implicit none
 public

 public:: M_Bar
 public:: namaxA1
 public:: namaxA2
 public:: namaxA3

 public :: Pr
 public :: q_k
 public :: Pk_Real
 public :: Pk_Imag
 public :: q_tmp
 public :: q_temp_RR
 public :: q_temp_RI
 public :: q_temp_IR
 public :: q_temp_II

 public :: q_temp_IRx
 public :: q_temp_IIx
 public :: q_temp_IRy
 public :: q_temp_IIy
 public :: q_temp_IRz
 type node
  real(DP) :: x
  real(DP) :: y
  real(DP) :: z
 end type 
  integer :: M_Bar
  integer,parameter :: namaxA1=113
  integer,parameter :: namaxA2=113
  integer,parameter :: namaxA3=64

 real(DP),allocatable :: Pr(:,:)
 real(DP),allocatable :: q_k(:,:)
 real(DP),allocatable :: Pk_Real(:,:)
 real(DP),allocatable :: Pk_Imag(:,:)
 real(DP),allocatable :: q_tmp(:,:)
 real(DP),allocatable :: q_temp_RR(:,:)
 real(DP),allocatable :: q_temp_RI(:,:)
 real(DP),allocatable :: q_temp_IR(:,:)
 real(DP),allocatable :: q_temp_II(:,:)

 real(DP),allocatable :: q_temp_IRx(:)
 real(DP),allocatable :: q_temp_IIx(:)
 real(DP),allocatable :: q_temp_IRy(:)
 real(DP),allocatable :: q_temp_IIy(:)
 real(DP),allocatable :: q_temp_IRz(:)
 real(DP),allocatable :: q_temp_IIz(:)
 type(node),allocatable :: k_vec(:)
 real(DP),allocatable :: c_nor_A(:)
 real(DP),allocatable :: c_nor_A_star(:)
 real(DP),allocatable :: c_nor_B(:)
 real(DP),allocatable :: c_nor_B_star(:)

contains

 subroutine MDE_onestep_init_array(R2C_FFT)
 USE nrtype,only :DP,PI
 USE global_para
 USE mpi
 USE control
 USE constants
 USE utility
 USE mmpi
 use decomp_fft_mod
 implicit none
 logical, optional :: R2C_FFT
 integer :: LOCAL_SIZEK,z_s,z_e,y_s,y_e,x_s,x_e
 integer :: k,K_i,K_ii,k_j,k_jj,k_k,k_kk
 integer :: istat2
 real(DP) :: k_x_factor,k_y_factor,k_z_factor

 call mp_barrier()
 write(*,*) "enter MDE onestep array init subroutine"
 call mp_barrier()
 M_Bar=max(M_Bar_A,M_Bar_B)
 
# if defined (R2C)
LOCAL_SIZEK=LOCAL_SIZE_K_R2C
# else /* R2C */
LOCAL_SIZEK=LOCAL_SIZE_K
# endif /* R2C */
 write(*,*) "LOCAL_SIZEK=",LOCAL_SIZEK

 allocate(k_vec(0:LOCAL_SIZEK-1)) 
 allocate(Pr(0:M_Bar-1,0:LOCAL_SIZE-1),stat=istat2)
 allocate(q_k(0:M_Bar-1,0:LOCAL_SIZE-1),stat=istat2)
 allocate(Pk_Real(0:M_Bar-1,0:LOCAL_SIZEK-1),stat=istat2)
 allocate(Pk_Imag(0:M_Bar-1,0:LOCAL_SIZEK-1),stat=istat2)

 allocate(q_tmp(0:M_Bar-1,0:LOCAL_SIZEK-1),stat=istat2)
 allocate(q_temp_RR(0:M_Bar-1,0:LOCAL_SIZEK-1),stat=istat2)
 allocate(q_temp_RI(0:M_Bar-1,0:LOCAL_SIZEK-1),stat=istat2)
 allocate(q_temp_IR(0:M_Bar-1,0:LOCAL_SIZEK-1),stat=istat2)
 allocate(q_temp_II(0:M_Bar-1,0:LOCAL_SIZEK-1),stat=istat2)

 allocate(q_temp_IRx(0:M_Bar-1),stat=istat2)
 allocate(q_temp_IIx(0:M_Bar-1),stat=istat2)
 allocate(q_temp_IRy(0:M_Bar-1),stat=istat2)
 allocate(q_temp_IIy(0:M_Bar-1),stat=istat2)
 allocate(q_temp_IRz(0:M_Bar-1),stat=istat2)
 allocate(q_temp_IIz(0:M_Bar-1),stat=istat2)

 allocate(c_nor_A(0:NA),stat=istat2)
 allocate(c_nor_A_star(0:NA),stat=istat2)

 allocate(c_nor_B(0:NB),stat=istat2)
 allocate(c_nor_B_star(0:NB),stat=istat2)
 if(istat2/=0) then
 write(*,*) "allocate matrix in mde_q failed"
 stop
 endif

 
k_x_factor=2.0d0*PI/(SIDEx*1.0d0)
k_y_factor=2.0d0*PI/(SIDEy*1.0d0)
k_z_factor=2.0d0*PI/(SIDEz*1.0d0)

 if(present(R2C_FFT)) then
  z_s=fft_start(3)
  z_e=fft_end(3)
  y_s=fft_start(2)
  y_e=fft_end(2)
  x_s=fft_start(1)
  x_e=fft_end(1)
 else
  z_s=zstart(3)
  z_e=zend(3)
  y_s=zstart(2)
  y_e=zend(2)
  x_s=zstart(1)
  x_e=zend(1)
endif
   
!k counter
k=0
do K_i=z_s,z_e
  if(k_i-1<=(SiDEz/2)) then
   k_ii=k_i-1
  else 
   k_ii=(k_i-1)-SIDEz 
  endif
    do k_j=y_s,y_e
        if(K_j-1<=(SIDEy/2)) then 
          k_jj=K_j-1
        else
          K_jj=K_j-1-SIDEy
        endif
           do k_k=x_s,x_e
             if(K_k-1<=(SIDEx/2)) then
             k_kk=K_k-1
             else
             k_kk=K_k-1-SIDEx
             endif

               k_vec(k)%x=k_kk*k_x_factor
               k_vec(k)%y=k_jj*k_y_factor
               k_vec(k)%z=k_ii*k_z_factor
               k=k+1    
          enddo
         enddo
        enddo

 end subroutine MDE_onestep_init_array

 subroutine MDE_q_norm(ar,s,dir)
 USE nrtype,only :DP,PI
 USE global_para
 USE mpi
 USE control
 USE constants
 USE utility
 USE mmpi
 use decomp_fft_mod
 implicit none
 real(DP) :: ar(0:Local_size-1)
 integer,intent(in) :: s,dir
 integer :: i,j,k
 real(DP) :: temp,temp_global

       temp=simposon_3D_1D_mpi(xsize(1),xsize(2),xsize(3),dx,dy,dz,ar)
              temp_global=0.0d0
              call mp_barrier()
              call mp_allreduce(temp,temp_global)
              call mp_barrier()
   temp=temp_global*4*PI
  if(dir==1) then
   c_nor_A(s)=temp
  else if(dir==-1) then
   c_nor_A_star(s)=temp
  
  else if(dir==2) then
   c_nor_B(s)=temp

  else if(dir==-2) then
   c_nor_B_star(s)=temp
  endif

 end subroutine MDE_q_norm

 subroutine MDE_onestep_AB_diblock_solver(R2C_FFT)
 USE nrtype,only :DP,PI
 USE global_para
 USE mpi
 USE control
 USE constants
 USE utility
 USE mmpi
 use decomp_fft_mod
 implicit none
 logical, optional :: R2C_FFT
 integer :: i,j,s
 integer :: index_nonzero
 integer :: k
 real(DP) :: ds_temp
 character :: MDEmethod
 real(DP) :: q_tmp(0:M_Bar-1,0:Local_size-1)
 real(DP) :: q_err(0:M_Bar-1,0:Local_size-1)
 real(DP) :: ar(0:Local_size-1)


ds_temp=ds

do k=0,LOCAL_SIZE-1
   do i=0,M_Bar-1
    qA(i,K,0)=0.0
    if(i==0) qA(i,K,0)=1.0
   enddo
enddo

  MDEmethod='E'

  do k=0,Local_size-1
      ar(k)=qA(0,k,0)
  enddo
  call MDE_q_norm(ar,0,1)
  qA(:,:,0)=qA(:,:,0)/c_nor_A(0)

!  c_nor_A=1.0
!  c_nor_A_star=1.0 
!  c_nor_B=1.0
!  c_nor_B_star=1.0 


 write(*,*) "calling AB slover driver","method=",MDEmethod

do s=1,NA

! maybe qA should be allocated as qA(i,k,s) not qA(s,i,k) ?? needs to be tested later

  if(MDEmethod=='I') then
   call MDE_onestep_q_im(qA(:,:,s-1), qA(:,:,s),WA,1,kapa_A_inv, & 
                      ds_temp,MDEmethod, &
                     R_G_val,I_GA1x_val,I_GA1x_col,I_GA1x_row, &
                     I_GA1y_val,I_GA1y_col,I_GA1y_row, &
                     I_GA1z_val,I_GA1z_col,I_GA1z_row, &
                      R2C_FFT)
   else if (MDEmethod=='M') then
  call mp_barrier()
   call MDE_onestep_q_im(qA(:,:,s-1), qA(:,:,s),WA,1,kapa_A_inv, & 
                      ds_temp,MDEmethod, &
                     R_G_val,I_GA1x_val,I_GA1x_col,I_GA1x_row, &
                     I_GA1y_val,I_GA1y_col,I_GA1y_row, &
                     I_GA1z_val,I_GA1z_col,I_GA1z_row, &
                      R2C_FFT)
    q_tmp(:,:)=qA(:,:,s)

  call mp_barrier()
               if(myid==0) then
                open(unit=44,file='q1-im.dat',status='replace')
                do j=0,M_Bar-1
                     write(44,*) q_tmp(j,0),j
                enddo
                close(44)
                endif 
  call mp_barrier()

   call MDE_onestep_q_ex2(qA(:,:,s-1), qA(:,:,s),q_tmp,WA,1,kapa_A_inv, & 
                      ds_temp,MDEmethod, &
                     R_G_val,I_GA1x_val,I_GA1x_col,I_GA1x_row, &
                     I_GA1y_val,I_GA1y_col,I_GA1y_row, &
                     I_GA1z_val,I_GA1z_col,I_GA1z_row, &
                      R2C_FFT)

  q_err(:,:)=qA(:,:,s)-q_tmp(:,:)
  write(*,*) "qA1 err",sqrt(sum(q_err(:,:)**2))
  call mp_barrier()
               if(myid==0) then
                open(unit=44,file='q1-exp.dat',status='replace')
                do j=0,M_Bar-1
                     write(44,*) qA(j,0,s),j
                enddo
                close(44)
                endif 
  call mp_barrier()

  else
!    q_tmp(:,:)=qA(:,:,s)
   call MDE_onestep_q_ex2(qA(:,:,s-1), qA(:,:,s),q_tmp,WA,1,kapa_A_inv, & 
                      ds_temp,MDEmethod, &
                     R_G_val,I_GA1x_val,I_GA1x_col,I_GA1x_row, &
                     I_GA1y_val,I_GA1y_col,I_GA1y_row, &
                     I_GA1z_val,I_GA1z_col,I_GA1z_row, &
                      R2C_FFT)
   endif

  do k=0,Local_size-1
      ar(k)=qA(0,k,s)
  enddo
  call MDE_q_norm(ar,s,1)
  qA(:,:,s)=qA(:,:,s)/c_nor_A(s)
enddo


!  call mp_barrier()
!   write(*,*) "qA(0,1,1",qA(0,1,1),"on",myid
!   write(*,*) "qA(0,2,1",qA(0,2,1),"on",myid
!   write(*,*) "qA(0,3,1",qA(0,3,1),"on",myid
!   write(*,*) "qA(1,1,1",qA(1,1,1)
!   write(*,*) "qA(2,1,1",qA(2,1,1)
!   write(*,*) "sum of qA(0,:,2)",sum(qA(0,:,1))
!   !stop
!  call mp_barrier()
!   write(*,*) "qA(0,1,NA",qA(0,1,NA),"on",myid
!   write(*,*) "qA(1,1,NA",qA(1,1,NA),"on",myid
!   write(*,*) "sum of qA(0,:,NA)",sum(qA(0,:,NA)),"on",myid
!  call mp_barrier()
!  stop


    qB(:,:,0)=qA(:,:,NA)
!do s=1,1
!   call MDE_onestep_q_im(qB(:,:,s-1), qB(:,:,s),WB,1,kapa_B_inv, & 
!                      ds_temp,MDEmethod, &
!                     R_G_val,I_GA1x_val,I_GA1x_col,I_GA1x_row, &
!                     I_GA1y_val,I_GA1y_col,I_GA1y_row, &
!                     I_GA1z_val,I_GA1z_col,I_GA1z_row, &
!                      R2C_FFT)
!enddo  
!  call mp_barrier()
!   write(*,*) "qB(0,1,1",qB(0,1,1)
!   write(*,*) "qB(0,2,1",qB(0,2,1)
!   write(*,*) "qB(1,1,1",qB(1,1,1)
!  call mp_barrier()

  do k=0,Local_size-1
      ar(k)=qB(0,k,0)
  enddo
  call MDE_q_norm(ar,0,2)

  qB(:,:,0)=qB(:,:,0)/c_nor_B(0)

do s=1,NB
 ! if(mod(s,2)/=0) then
  if(MDEmethod=='I') then
   call MDE_onestep_q_im(qB(:,:,s-1), qB(:,:,s),WB,1,kapa_B_inv, & 
                      ds_temp,MDEmethod, &
                     R_G_val,I_GA1x_val,I_GA1x_col,I_GA1x_row, &
                     I_GA1y_val,I_GA1y_col,I_GA1y_row, &
                     I_GA1z_val,I_GA1z_col,I_GA1z_row, &
                      R2C_FFT)
   else if (MDEmethod=='M') then
   call MDE_onestep_q_im(qB(:,:,s-1), qB(:,:,s),WB,1,kapa_B_inv, & 
                      ds_temp,MDEmethod, &
                     R_G_val,I_GA1x_val,I_GA1x_col,I_GA1x_row, &
                     I_GA1y_val,I_GA1y_col,I_GA1y_row, &
                     I_GA1z_val,I_GA1z_col,I_GA1z_row, &
                      R2C_FFT)
    q_tmp(:,:)=qB(:,:,s)

   call MDE_onestep_q_ex2(qB(:,:,s-1), qB(:,:,s),q_tmp,WB,1,kapa_B_inv, & 
                      ds_temp,MDEmethod, &
                     R_G_val,I_GA1x_val,I_GA1x_col,I_GA1x_row, &
                     I_GA1y_val,I_GA1y_col,I_GA1y_row, &
                     I_GA1z_val,I_GA1z_col,I_GA1z_row, &
                      R2C_FFT)
   
   else
   call MDE_onestep_q_ex2(qB(:,:,s-1), qB(:,:,s),q_tmp,WB,1,kapa_B_inv, & 
                      ds_temp,MDEmethod, &
                     R_G_val,I_GA1x_val,I_GA1x_col,I_GA1x_row, &
                     I_GA1y_val,I_GA1y_col,I_GA1y_row, &
                     I_GA1z_val,I_GA1z_col,I_GA1z_row, &
                      R2C_FFT)
!   else
!   call MDE_onestep_q_im(qB(:,:,s-1), qB(:,:,s),WB,1,kapa_B_inv, & 
!                      ds_temp,MDEmethod, &
!                     R_G_val,I_GA1x_val,I_GA1x_col,I_GA1x_row, &
!                     I_GA1y_val,I_GA1y_col,I_GA1y_row, &
!                     I_GA1z_val,I_GA1z_col,I_GA1z_row, &
!                      R2C_FFT)
   endif
  do k=0,Local_size-1
      ar(k)=qB(0,k,s)
  enddo
  call MDE_q_norm(ar,s,2)
  qB(:,:,s)=qB(:,:,s)/c_nor_B(s)
enddo  

  call mp_barrier()
   write(*,*) "qB(0,1,NB",qB(0,1,NB)
   write(*,*) "qB(0,2,NB",qB(0,2,NB)
   write(*,*) "qB(1,1,NB",qB(1,1,NB)
  call mp_barrier()

do k=0,LOCAL_SIZE-1
   do i=0,M_Bar-1
    qBstar(i,K,NB)=0.0
    if(i==0) qBstar(i,K,NB)=1.0
   enddo
enddo

  do k=0,Local_size-1
      ar(k)=qBstar(0,k,NB)
  enddo
  call MDE_q_norm(ar,NB,-2)

  qBstar(:,:,NB)=qBstar(:,:,NB)/c_nor_B_star(NB)
do s=NB,1,-1
!   if(mod(s,2)==0) then
  if(MDEmethod=='I') then
   call MDE_onestep_q_im(qBstar(:,:,s), qBstar(:,:,s-1),WB,-1, kapa_B_inv,& 
                      ds_temp,MDEmethod, &
                     R_G_val,I_GA1x_val,I_GA1x_col,I_GA1x_row, &
                     I_GA1y_val,I_GA1y_col,I_GA1y_row, &
                     I_GA1z_val,I_GA1z_col,I_GA1z_row, &
                     R2C_FFT)
   else if (MDEmethod=='M') then
   call MDE_onestep_q_im(qBstar(:,:,s), qBstar(:,:,s-1),WB,-1, kapa_B_inv,& 
                      ds_temp,MDEmethod, &
                     R_G_val,I_GA1x_val,I_GA1x_col,I_GA1x_row, &
                     I_GA1y_val,I_GA1y_col,I_GA1y_row, &
                     I_GA1z_val,I_GA1z_col,I_GA1z_row, &
                     R2C_FFT)
     q_tmp(:,:)=qBstar(:,:,s-1)
   call MDE_onestep_q_ex2(qBstar(:,:,s), qBstar(:,:,s-1),q_tmp,WB,-1, kapa_B_inv,& 
                      ds_temp,MDEmethod, &
                     R_G_val,I_GA1x_val,I_GA1x_col,I_GA1x_row, &
                     I_GA1y_val,I_GA1y_col,I_GA1y_row, &
                     I_GA1z_val,I_GA1z_col,I_GA1z_row, &
                     R2C_FFT)
   else
   call MDE_onestep_q_ex2(qBstar(:,:,s), qBstar(:,:,s-1),q_tmp,WB,-1, kapa_B_inv,& 
                      ds_temp,MDEmethod, &
                     R_G_val,I_GA1x_val,I_GA1x_col,I_GA1x_row, &
                     I_GA1y_val,I_GA1y_col,I_GA1y_row, &
                     I_GA1z_val,I_GA1z_col,I_GA1z_row, &
                     R2C_FFT)

 endif
  do k=0,Local_size-1
      ar(k)=qBstar(0,k,s)
  enddo
  call MDE_q_norm(ar,s,-2)
  qBstar(:,:,s)=qBstar(:,:,s)/c_nor_B_star(s)

enddo
  call mp_barrier()
   write(*,*) "qBstar(0,1,0)",qBstar(0,1,0)
   write(*,*) "qBstar(1,1,0",qB(1,1,0)
  call mp_barrier()

    qAstar(:,:,NA)=qBstar(:,:,0)

  do k=0,Local_size-1
      ar(k)=qAstar(0,k,NA)
  enddo
  call MDE_q_norm(ar,NA,-1)
  qAstar(:,:,NA)=qAstar(:,:,NA)/c_nor_A_star(NA)
do s=NA,1,-1
  if(MDEmethod=='I') then
   call MDE_onestep_q_im(qAstar(:,:,s), qAstar(:,:,s-1),WA,-1,kapa_A_inv, & 
                      ds_temp,MDEmethod, &
                     R_G_val,I_GA1x_val,I_GA1x_col,I_GA1x_row, &
                     I_GA1y_val,I_GA1y_col,I_GA1y_row, &
                     I_GA1z_val,I_GA1z_col,I_GA1z_row, &
                     R2C_FFT)
   else if (MDEmethod=='M') then
   call MDE_onestep_q_im(qAstar(:,:,s), qAstar(:,:,s-1),WA,-1,kapa_A_inv, & 
                      ds_temp,MDEmethod, &
                     R_G_val,I_GA1x_val,I_GA1x_col,I_GA1x_row, &
                     I_GA1y_val,I_GA1y_col,I_GA1y_row, &
                     I_GA1z_val,I_GA1z_col,I_GA1z_row, &
                     R2C_FFT)
    q_tmp(:,:)=qAstar(:,:,s-1)
   call MDE_onestep_q_ex2(qAstar(:,:,s), qAstar(:,:,s-1),q_tmp,WA,-1,kapa_A_inv, & 
                      ds_temp,MDEmethod, &
                     R_G_val,I_GA1x_val,I_GA1x_col,I_GA1x_row, &
                     I_GA1y_val,I_GA1y_col,I_GA1y_row, &
                     I_GA1z_val,I_GA1z_col,I_GA1z_row, &
                     R2C_FFT)
else
   call MDE_onestep_q_ex2(qAstar(:,:,s), qAstar(:,:,s-1),q_tmp,WA,-1,kapa_A_inv, & 
                      ds_temp,MDEmethod, &
                     R_G_val,I_GA1x_val,I_GA1x_col,I_GA1x_row, &
                     I_GA1y_val,I_GA1y_col,I_GA1y_row, &
                     I_GA1z_val,I_GA1z_col,I_GA1z_row, &
                     R2C_FFT)
endif
  do k=0,Local_size-1
      ar(k)=qAstar(0,k,s)
  enddo
  call MDE_q_norm(ar,s,-1)
  qAstar(:,:,s)=qAstar(:,:,s)/c_nor_A_star(s)

enddo

!  call mp_barrier()
!   write(*,*) "qAstar(0,1,0)",qAstar(0,1,0)
!   write(*,*) "qAstar(1,1,0",qA(1,1,0)
!  call mp_barrier()
! 
! call mp_barrier()
! write(*,*) "sum of qB(0,:,NB)",sum(qB(0,:,NB))
! write(*,*) "done MDE onestep solver"
! call mp_barrier()
!
 end subroutine MDE_onestep_AB_diblock_solver



 subroutine MDE_onestep_q(q_in, q_out,W_temp,orient,kapa_inv, & 
                             ds_temp,MDEmethod, & 
            R_GG_val,I_Gx_val,I_Gx_col,I_Gx_row, &
                    I_Gy_val,I_Gy_col,I_Gy_row, &
                    I_Gz_val,I_Gz_col,I_Gz_row, &
                          R2C_FFT)
 USE nrtype,only :DP,PI
 USE global_para
 USE mpi
 USE control
 USE constants
 USE utility
 USE mmpi
 use decomp_fft_mod
 !USE G_matrix_2dcomp_mod
 USE matrix_inverse
 USE sparse_MKL
 implicit none
 logical, optional :: R2C_FFT
 integer :: i,j,k,s
 integer :: index_nonzero
 integer :: istat2
 real(DP) :: Wtemp,error
 REAL(DP) :: M_grid_inv,ds_tmp
 integer :: Local_sizek
 real(DP),intent(in) :: ds_temp,kapa_inv

 real(DP),intent(in) :: R_GG_val(0:M_Bar-1)
 real(DP),intent(in) :: I_Gx_val(1:namaxA1)
 integer,intent(in) ::  I_Gx_col(1:namaxA1)
 integer,intent(in) ::  I_Gx_row(1:M_Bar+1)
 real(DP),intent(in) :: I_Gy_val(1:namaxA2)
 integer,intent(in) ::  I_Gy_col(1:namaxA2)
 integer,intent(in) ::  I_Gy_row(1:M_Bar+1)
 real(DP),intent(in) :: I_Gz_val(1:namaxA3)
 integer,intent(in) ::  I_Gz_col(1:namaxA3)
 integer,intent(in) ::  I_Gz_row(1:M_Bar+1)
 
 real(DP),intent(in) :: q_in(0:M_Bar-1,0:Local_size-1)
 real(DP),intent(in) :: W_temp(0:LOCAL_SIZE-1)
 real(DP),intent(inout) :: q_out(0:M_Bar-1,0:Local_size-1)
 character,intent(in) :: MDEmethod
 integer, intent(in) :: orient
 real(DP) :: G_R(0:M_Bar-1,0:M_Bar-1)
 real(DP) :: G_Ix(0:M_Bar-1,0:M_Bar-1)
 real(DP) :: G_Iy(0:M_Bar-1,0:M_Bar-1)
 real(DP) :: G_Iz(0:M_Bar-1,0:M_Bar-1)
integer job(8)
integer info

 if(orient==1) then
 ds_tmp=ds_temp
 else if(orient==-1) then
 ds_tmp=-1.0d0*ds_temp
 else
 write(*,*) "error of orient input"
 stop
 endif


# if defined (R2C)
LOCAL_SIZEK=LOCAL_SIZE_K_R2C
# else /* R2C */
LOCAL_SIZEK=LOCAL_SIZE_K
# endif /* R2C */

 M_grid_inv=1.0d0/M_grid

          Pr=0.0
          q_k(:,:)=q_in(:,:)
      do k=0,LOCAL_SIZE-1
          Wtemp=W_temp(k)
            do i=0,M_Bar-1
            Pr(i,k)=Pr(i,k)+(-ds_temp*Wtemp)*q_in(i,k)
            enddo
        enddo !enddo k=0,LOCAL_SIZE-1  

!!!doing the FFT
!!Note that  in 2dcomp fft , local_in and local_out is formed as (0:x,0:y,0:z), which
!is different from fft slab decomposition version.
# if defined (R2C)
      do i=0,M_Bar-1
        do k=0,LOCAL_SIZE-1
          local_in_r(kD(k)%x,kD(k)%y,KD(k)%z)=q_in(i,k)
          !local_in_r(kD(k)%x,kD(k)%y,KD(k)%z)=q_k(i,k)
         enddo
       
        call decomp_2d_fft_3d(local_in_r, local_out_c)
 !      call fftw_mpi_execute_dft(fplan, local_in, local_out)
       
        do k=0,LOCAL_SIZE_K_R2C-1
            Pk_Real(i,k)=real(local_out_c(kDk(k)%x,kDk(k)%y,KDk(k)%z))
            Pk_Imag(i,k)=AIMAG(local_out_c(kDk(k)%x,kDk(k)%y,KDk(k)%z))
         enddo
       enddo  !enddo i=0,M_Bar_A-1  

# else /* R2C */
      do i=0,M_Bar-1
        do k=0,LOCAL_SIZE-1
         !  local_in(kD(k)%x,kD(k)%y,KD(k)%z)=dcmplx(Pr(i,k),0.0d0)
         ! local_in(kD(k)%x,kD(k)%y,KD(k)%z)=dcmplx(q_in(i,k),0.0d0)
          local_in(kD(k)%x,kD(k)%y,KD(k)%z)=dcmplx(q_tmp(i,k),0.0d0)
         enddo
       
        call decomp_2d_fft_3d(local_in, local_out, DECOMP_2D_FFT_FORWARD)
 !      call fftw_mpi_execute_dft(fplan, local_in, local_out)
       
        do k=0,LOCAL_SIZE_K-1
            Pk_Real(i,k)=real(local_out(kDz(k)%x,kDz(k)%y,KDz(k)%z))
            Pk_Imag(i,k)=AIMAG(local_out(kDz(k)%x,kDz(k)%y,KDz(k)%z))
         enddo
       enddo  !enddo i=0,M_Bar_A-1  
# endif /* R2C */

         do k=0,LOCAL_SIZEK-1
          do  j=0,M_Bar-1
          q_temp_RR(j,k)=ds_temp*kapa_inv*R_GG_val(j)*Pk_Real(j,k)
          q_temp_RI(j,k)=ds_temp*kapa_inv*R_GG_val(j)*Pk_Imag(j,k)
          enddo 
         call mkl_dcsrsymv('l', M_Bar,  I_Gx_val, I_Gx_row, &
            I_Gx_col,Pk_Real(:,k), q_temp_IRx)
         call mkl_dcsrsymv('l', M_Bar,  I_Gy_val, I_Gy_row, &
            I_Gy_col,Pk_Real(:,k), q_temp_IRy)
         call mkl_dcsrsymv('l', M_Bar,  I_Gz_val, I_Gz_row, &
            I_Gz_col,Pk_Real(:,k), q_temp_IRz)

         call mkl_dcsrsymv('l', M_Bar,  I_Gx_val, I_Gx_row, &
            I_Gx_col,Pk_Imag(:,k), q_temp_IIx)
         call mkl_dcsrsymv('l', M_Bar,  I_Gy_val, I_Gy_row, &
            I_Gy_col,Pk_Imag(:,k), q_temp_IIy)
         call mkl_dcsrsymv('l', M_Bar,  I_Gz_val, I_Gz_row, &
            I_Gz_col,Pk_Imag(:,k), q_temp_IIz)

            
          do j=0,M_Bar-1
             q_temp_IR(j,k)=ds_tmp*dx_inv*k_vec(k)%x*q_temp_IRx(j) + &
                            ds_tmp*dy_inv*k_vec(k)%y*q_temp_IRy(j) + &
                            ds_tmp*dz_inv*k_vec(k)%z*q_temp_IRz(j)    
             q_temp_II(j,k)=ds_tmp*dx_inv*k_vec(k)%x*q_temp_IIx(j) + &
                            ds_tmp*dy_inv*k_vec(k)%y*q_temp_IIy(j) + &
                            ds_tmp*dz_inv*k_vec(k)%z*q_temp_IIz(j)    
          enddo


         enddo


# if defined (R2C)
             do i=0,M_Bar-1
                do k=0,LOCAL_SIZE_K_R2C-1
     local_out_c(kDk(k)%x,kDk(k)%y,KDk(k)%z)=dcmplx(q_temp_RR(i,k)-q_temp_II(i,k), & 
                q_temp_RI(i,k)+q_temp_IR(i,k))
                enddo
       !call fftw_mpi_execute_dft(bplan, local_out, local_in)
        call decomp_2d_fft_3d(local_out_c, local_in_r)
            
                do k=0,LOCAL_SIZE-1
                 q_out(i,k)=local_in_r(kD(k)%x,kD(k)%y,KD(k)%z)*M_grid_inv + &
                            Pr(i,k) + q_in(i,k)  
                enddo
              enddo  
# else /* R2C */
             do i=0,M_Bar-1
                do k=0,LOCAL_SIZE_K-1
        local_out(kDz(k)%x,kDz(k)%y,KDz(k)%z)=dcmplx(q_temp_RR(i,k)-q_temp_II(i,k), & 
                q_temp_RI(i,k)+q_temp_IR(i,k))
                enddo
        !        if(i==1) then 
        !        write(*,*) "i==1,sum local_out",sum(local_out),"on",myid
        !        endif
       !call fftw_mpi_execute_dft(bplan, local_out, local_in)
        call decomp_2d_fft_3d(local_out, local_in, DECOMP_2D_FFT_BACKWARD)
 !           error=0.0d0
                do k=0,LOCAL_SIZE-1
              !  q_out(i,k)=real(local_in(kD(k)%x,kD(k)%y,KD(k)%z))*M_grid_inv 
               q_out(i,k)=real(local_in(kD(k)%x,kD(k)%y,KD(k)%z))*M_grid_inv + &
                           Pr(i,k) + q_in(i,k)
  !              error=abs(q_out(i,k)-q_k(i,k)) 
                enddo
              enddo  
# endif /* R2C */

write(*,*) "sum(q_out(i,k)",sum(q_out(:,:))
 end subroutine MDE_onestep_q


 subroutine MDE_onestep_q_im(q_in, q_out,W_temp,orient,kapa_inv, & 
                             ds_temp,MDEmethod, & 
            R_GG_val,I_Gx_val,I_Gx_col,I_Gx_row, &
                    I_Gy_val,I_Gy_col,I_Gy_row, &
                    I_Gz_val,I_Gz_col,I_Gz_row, &
                          R2C_FFT)
 USE nrtype,only :DP,PI
 USE global_para
 USE mpi
 USE control
 USE constants
 USE utility
 USE mmpi
 use decomp_fft_mod
 !USE G_matrix_2dcomp_mod
 USE matrix_inverse
 USE sparse_MKL
 implicit none
 logical, optional :: R2C_FFT
 integer :: i,j,k,s
 integer :: index_nonzero
 integer :: istat2
 real(DP) :: Wtemp,error
 REAL(DP) :: M_grid_inv,ds_tmp
 integer :: Local_sizek
 real(DP),intent(in) :: ds_temp,kapa_inv
 real(DP),intent(inout) :: R_GG_val(0:M_Bar-1)
 real(DP),intent(in) :: I_Gx_val(1:namaxA1)
 integer,intent(in) ::  I_Gx_col(1:namaxA1)
 integer,intent(in) ::  I_Gx_row(1:M_Bar+1)
 real(DP),intent(in) :: I_Gy_val(1:namaxA2)
 integer,intent(in) ::  I_Gy_col(1:namaxA2)
 integer,intent(in) ::  I_Gy_row(1:M_Bar+1)
 real(DP),intent(in) :: I_Gz_val(1:namaxA3)
 integer,intent(in) ::  I_Gz_col(1:namaxA3)
 integer,intent(in) ::  I_Gz_row(1:M_Bar+1)
 
 real(DP),intent(in) :: q_in(0:M_Bar-1,0:Local_size-1)
 real(DP),intent(in) :: W_temp(0:LOCAL_SIZE-1)
 real(DP),intent(inout) :: q_out(0:M_Bar-1,0:Local_size-1)
 character,intent(in) :: MDEmethod
 integer, intent(in) :: orient
 real(DP) :: G_R(0:M_Bar-1,0:M_Bar-1)
 real(DP) :: G_Ix(0:M_Bar-1,0:M_Bar-1)
 real(DP) :: G_Iy(0:M_Bar-1,0:M_Bar-1)
 real(DP) :: G_Iz(0:M_Bar-1,0:M_Bar-1)
integer job(8)
integer info

 if(orient==1) then
 ds_tmp=ds_temp
 else if(orient==-1) then
 ds_tmp=-1.0d0*ds_temp
 else
 write(*,*) "error of orient input"
 stop
 endif


# if defined (R2C)
LOCAL_SIZEK=LOCAL_SIZE_K_R2C
# else /* R2C */
LOCAL_SIZEK=LOCAL_SIZE_K
# endif /* R2C */

 M_grid_inv=1.0d0/M_grid

          Pr=0.0

      do k=0,LOCAL_SIZE-1
          Wtemp=W_temp(k)
            do i=0,M_Bar-1
            Pr(i,k)=Pr(i,k)+(1.0-ds_temp*Wtemp)*q_in(i,k)
            enddo
        enddo !enddo k=0,LOCAL_SIZE-1  

         write(*,*) "sum Pr(0,:)",sum(Pr(0,:)),"on",myid      
         write(*,*) " Pr_(0:1,localsize-1)",Pr(0,0),Pr(0,LOCAL_SIZE-1),"on",myid

!!!doing the FFT
!!Note that  in 2dcomp fft , local_in and local_out is formed as (0:x,0:y,0:z), which
!is different from fft slab decomposition version.
# if defined (R2C)
      do i=0,M_Bar-1
        do k=0,LOCAL_SIZE-1
          local_in_r(kD(k)%x,kD(k)%y,KD(k)%z)=Pr(i,k)
          !local_in_r(kD(k)%x,kD(k)%y,KD(k)%z)=q_k(i,k)
         enddo
       
        call decomp_2d_fft_3d(local_in_r, local_out_c)
 !      call fftw_mpi_execute_dft(fplan, local_in, local_out)
       
        do k=0,LOCAL_SIZE_K_R2C-1
            Pk_Real(i,k)=real(local_out_c(kDk(k)%x,kDk(k)%y,KDk(k)%z))
            Pk_Imag(i,k)=AIMAG(local_out_c(kDk(k)%x,kDk(k)%y,KDk(k)%z))
         enddo
       enddo  !enddo i=0,M_Bar_A-1  

# else /* R2C */
      do i=0,M_Bar-1
        do k=0,LOCAL_SIZE-1
          local_in(kD(k)%x,kD(k)%y,KD(k)%z)=dcmplx(Pr(i,k),0.0d0)
         ! local_in(kD(k)%x,kD(k)%y,KD(k)%z)=dcmplx(q_k(i,k),0.0d0)
         enddo
       
        call decomp_2d_fft_3d(local_in, local_out, DECOMP_2D_FFT_FORWARD)
 !      call fftw_mpi_execute_dft(fplan, local_in, local_out)
       
        do k=0,LOCAL_SIZE_K-1
            Pk_Real(i,k)=real(local_out(kDz(k)%x,kDz(k)%y,KDz(k)%z))
            Pk_Imag(i,k)=AIMAG(local_out(kDz(k)%x,kDz(k)%y,KDz(k)%z))
         enddo
       enddo  !enddo i=0,M_Bar_A-1  
# endif /* R2C */

              write(*,*) "sum of Pk",sum(Pk_Real(0,:)),sum(Pk_Imag(0,:))
              write(*,*) " Pk_(0,0)",Pk_Real(0,0),Pk_Imag(0,0),"on",myid

                 q_temp_RR=0.0d0     
                 q_temp_RI=0.0d0     
                 q_temp_IR=0.0d0     
                 q_temp_II=0.0d0     

         do k=0,LOCAL_SIZEK-1

        job(1)=1
        job(2)=1
        job(3)=1
        job(4)=2
        job(5)=113
        job(6)=1
      call mkl_ddnscsr(job,M_Bar,M_Bar,G_Ix,M_Bar,I_Gx_val,I_Gx_col, &
        I_Gx_row,info)
      call mkl_ddnscsr(job,M_Bar,M_Bar,G_Iy,M_Bar,I_Gy_val,I_Gy_col, &
        I_Gy_row,info)
        job(5)=64
      call mkl_ddnscsr(job,M_Bar,M_Bar,G_Iz,M_Bar,I_Gz_val,I_Gz_col, &
        I_Gz_row,info)
          
             G_R=0.0d0

              G_Ix(:,:)=-ds_tmp*dx_inv*k_vec(k)%x*G_Ix(:,:)
              G_Iy(:,:)=-ds_tmp*dy_inv*k_vec(k)%y*G_Iy(:,:)
              G_Iz(:,:)=-ds_tmp*dz_inv*k_vec(k)%z*G_Iz(:,:)

             do j=0,M_Bar-1
                 G_R(j,j)=1.0d0-R_GG_val(j)*ds_temp*kapa_inv
                do i=0,M_Bar-1    
                if(j>i) then
                G_Ix(i,j)=G_Ix(j,i)
                G_Iy(i,j)=G_Iy(j,i)
                G_Iz(i,j)=G_Iz(j,i)
                endif
                enddo
              enddo

                 G_Ix=G_Ix+G_Iy+G_Iz
              call complex_matrix_inverse(G_R,G_Ix,M_Bar)

               if(myid==10) then
                open(unit=44,file='G.dat',status='replace')
               do i=0,M_Bar-1
                do j=0,M_Bar-1
                     write(44,*) G_R(i,j),G_Ix(i,j),i,j
                enddo
                enddo
                close(44)
                endif 
           


              
               do i=0,M_Bar-1
                do j=0,M_Bar-1
             q_temp_RR(i,k)=q_temp_RR(i,k)+G_R(i,j)*Pk_real(j,k)
             q_temp_RI(i,k)=q_temp_RI(i,k)+G_R(i,j)*Pk_Imag(j,k)
             q_temp_IR(i,k)=q_temp_IR(i,k)+G_Ix(i,j)*Pk_real(j,k)
             q_temp_II(i,k)=q_temp_II(i,k)+G_Ix(i,j)*Pk_Imag(j,k)
                  enddo
                enddo                    

         enddo

              write(*,*) "sum of q_temp(:,0)",sum(q_temp_RR(:,0)),sum(q_temp_RI(:,0)),"on",myid
# if defined (R2C)
             do i=0,M_Bar-1
                do k=0,LOCAL_SIZE_K_R2C-1
     local_out_c(kDk(k)%x,kDk(k)%y,KDk(k)%z)=dcmplx(q_temp_RR(i,k)-q_temp_II(i,k), & 
                q_temp_RI(i,k)+q_temp_IR(i,k))
                enddo
       !call fftw_mpi_execute_dft(bplan, local_out, local_in)
        call decomp_2d_fft_3d(local_out_c, local_in_r)
            
                do k=0,LOCAL_SIZE-1
                 q_out(i,k)=local_in_r(kD(k)%x,kD(k)%y,KD(k)%z)*M_grid_inv 
                            
                enddo
              enddo  
# else /* R2C */
             do i=0,M_Bar-1
                do k=0,LOCAL_SIZE_K-1
    local_out(kDz(k)%x,kDz(k)%y,KDz(k)%z)=dcmplx(q_temp_RR(i,k)-q_temp_II(i,k), & 
               q_temp_RI(i,k)+q_temp_IR(i,k))
                enddo
       !call fftw_mpi_execute_dft(bplan, local_out, local_in)
        call decomp_2d_fft_3d(local_out, local_in, DECOMP_2D_FFT_BACKWARD)
            
                do k=0,LOCAL_SIZE-1
                 q_out(i,k)=real(local_in(kD(k)%x,kD(k)%y,KD(k)%z))*M_grid_inv 
                enddo
              enddo  
# endif /* R2C */

write(*,*) "sum(q_out(i,k)",sum(q_out(:,:)),"on",myid
 end subroutine MDE_onestep_q_im

 subroutine MDE_onestep_q_ex2(q_in, q_out,q_tmp,W_temp,orient,kapa_inv, & 
                             ds_temp,MDEmethod, & 
            R_GG_val,I_Gx_val,I_Gx_col,I_Gx_row, &
                    I_Gy_val,I_Gy_col,I_Gy_row, &
                    I_Gz_val,I_Gz_col,I_Gz_row, &
                          R2C_FFT)
 USE nrtype,only :DP,PI
 USE global_para
 USE mpi
 USE control
 USE constants
 USE utility
 USE mmpi
 use decomp_fft_mod
 !USE G_matrix_2dcomp_mod
 USE matrix_inverse
 USE sparse_MKL
 implicit none
 logical, optional :: R2C_FFT
 integer :: i,j,k,s
 integer :: index_nonzero
 integer :: istat2
 real(DP) :: Wtemp,error
 REAL(DP) :: M_grid_inv,ds_tmp
 integer :: Local_sizek
 real(DP),intent(in) :: ds_temp,kapa_inv
 real(DP),intent(inout) :: R_GG_val(0:M_Bar-1)
 real(DP),intent(in) :: I_Gx_val(1:namaxA1)
 integer,intent(in) ::  I_Gx_col(1:namaxA1)
 integer,intent(in) ::  I_Gx_row(1:M_Bar+1)
 real(DP),intent(in) :: I_Gy_val(1:namaxA2)
 integer,intent(in) ::  I_Gy_col(1:namaxA2)
 integer,intent(in) ::  I_Gy_row(1:M_Bar+1)
 real(DP),intent(in) :: I_Gz_val(1:namaxA3)
 integer,intent(in) ::  I_Gz_col(1:namaxA3)
 integer,intent(in) ::  I_Gz_row(1:M_Bar+1)
 
 real(DP),intent(in) :: q_in(0:M_Bar-1,0:Local_size-1)
 real(DP),intent(in) :: q_tmp(0:M_Bar-1,0:Local_size-1)
 real(DP),intent(in) :: W_temp(0:LOCAL_SIZE-1)
 real(DP),intent(inout) :: q_out(0:M_Bar-1,0:Local_size-1)
 character,intent(in) :: MDEmethod
 integer, intent(in) :: orient
 real(DP) :: ar(0:Local_size-1)
 real(DP) :: G_R(0:M_Bar-1,0:M_Bar-1)
 real(DP) :: G_Ix(0:M_Bar-1,0:M_Bar-1)
 real(DP) :: G_Iy(0:M_Bar-1,0:M_Bar-1)
 real(DP) :: G_Iz(0:M_Bar-1,0:M_Bar-1)
integer job(8)
integer info

 if(orient==1) then
 ds_tmp=ds_temp
 else if(orient==-1) then
 ds_tmp=-1.0d0*ds_temp
 else
 write(*,*) "error of orient input"
 stop
 endif


# if defined (R2C)
LOCAL_SIZEK=LOCAL_SIZE_K_R2C
# else /* R2C */
LOCAL_SIZEK=LOCAL_SIZE_K
# endif /* R2C */

 M_grid_inv=1.0d0/M_grid

          Pr=0.0

      do k=0,LOCAL_SIZE-1
          Wtemp=W_temp(k)
            do i=0,M_Bar-1
            !Pr(i,k)=Pr(i,k)+( -ds_temp*Wtemp)*q_in(i,k)
            Pr(i,k)= -ds_temp*Wtemp*q_in(i,k)
            enddo
        enddo !enddo k=0,LOCAL_SIZE-1  

!!!doing the FFT
!!Note that  in 2dcomp fft , local_in and local_out is formed as (0:x,0:y,0:z), which
!is different from fft slab decomposition version.
# if defined (R2C)
      do i=0,M_Bar-1
        do k=0,LOCAL_SIZE-1
          local_in_r(kD(k)%x,kD(k)%y,KD(k)%z)=Pr(i,k)
          !local_in_r(kD(k)%x,kD(k)%y,KD(k)%z)=q_k(i,k)
         enddo
       
        call decomp_2d_fft_3d(local_in_r, local_out_c)
 !      call fftw_mpi_execute_dft(fplan, local_in, local_out)
       
        do k=0,LOCAL_SIZE_K_R2C-1
            Pk_Real(i,k)=real(local_out_c(kDk(k)%x,kDk(k)%y,KDk(k)%z))
            Pk_Imag(i,k)=AIMAG(local_out_c(kDk(k)%x,kDk(k)%y,KDk(k)%z))
         enddo
       enddo  !enddo i=0,M_Bar_A-1  

# else /* R2C */
      do i=0,M_Bar-1
        do k=0,LOCAL_SIZE-1
         ! local_in(kD(k)%x,kD(k)%y,KD(k)%z)=dcmplx(Pr(i,k),0.0d0)
           if(MDEmethod=='M') then
          local_in(kD(k)%x,kD(k)%y,KD(k)%z)=dcmplx(q_in(i,k),0.0d0)
          !local_in(kD(k)%x,kD(k)%y,KD(k)%z)=dcmplx(q_tmp(i,k),0.0d0)
           else
          local_in(kD(k)%x,kD(k)%y,KD(k)%z)=dcmplx(q_in(i,k),0.0d0)
           endif
         enddo
       
        call decomp_2d_fft_3d(local_in, local_out, DECOMP_2D_FFT_FORWARD)
 !      call fftw_mpi_execute_dft(fplan, local_in, local_out)
       
        do k=0,LOCAL_SIZE_K-1
            Pk_Real(i,k)=real(local_out(kDz(k)%x,kDz(k)%y,KDz(k)%z))
            Pk_Imag(i,k)=AIMAG(local_out(kDz(k)%x,kDz(k)%y,KDz(k)%z))
         enddo
       enddo  !enddo i=0,M_Bar_A-1  
# endif /* R2C */

                if(myid==0) then
                open(unit=44,file='Pk_temp.dat',status='replace')
                do j=0,M_Bar-1
                     write(44,*) Pk_real(j,0),Pk_Imag(j,0),q_in(j,0)
                     write(44,*) j
                enddo
                close(44)
                endif

                 q_temp_RR=0.0d0     
                 q_temp_RI=0.0d0     
                 q_temp_IR=0.0d0     
                 q_temp_II=0.0d0     

         do k=0,LOCAL_SIZEK-1

         G_R=0.0d0
         do j=0,M_Bar-1        
        G_R(j,j)=  - Loa*basis_SPH(j)%x* &
                (basis_SPH(j)%x+1)/(2.0d0)*ds_temp*kapa_inv
                do i=0,M_Bar-1
       G_Ix(i,j)=-1.0d0*matrix_Rx(i,j)
       G_Iy(i,j)=-1.0d0*matrix_Ry(i,j)
       G_Iz(i,j)=-1.0d0*matrix_Rz(i,j)
                enddo
         enddo




!        job(1)=1
!        job(2)=1
!        job(3)=1
!        job(4)=2
!        job(5)=113
!        job(6)=1
!      call mkl_ddnscsr(job,M_Bar,M_Bar,G_Ix,M_Bar,I_Gx_val,I_Gx_col, &
!        I_Gx_row,info)
!      call mkl_ddnscsr(job,M_Bar,M_Bar,G_Iy,M_Bar,I_Gy_val,I_Gy_col, &
!        I_Gy_row,info)
!        job(5)=64
!      call mkl_ddnscsr(job,M_Bar,M_Bar,G_Iz,M_Bar,I_Gz_val,I_Gz_col, &
!        I_Gz_row,info)
!          
!             G_R=0.0d0
!
              G_Ix(:,:)=ds_tmp*dx_inv*k_vec(k)%x*G_Ix(:,:)
              G_Iy(:,:)=ds_tmp*dy_inv*k_vec(k)%y*G_Iy(:,:)
              G_Iz(:,:)=ds_tmp*dz_inv*k_vec(k)%z*G_Iz(:,:)
!
!             do j=0,M_Bar-1
!                ! G_R(j,j)=1.0d0-R_GG_val(j)*ds_temp*kapa_inv
!                  G_R(j,j)=R_GG_val(j)*ds_temp*kapa_inv
!                do i=0,M_Bar-1    
!                if(j>i) then
!                G_Ix(i,j)=G_Ix(j,i)
!                G_Iy(i,j)=G_Iy(j,i)
!                G_Iz(i,j)=G_Iz(j,i)
!                endif
!                enddo
!              enddo

                 G_Ix=G_Ix+G_Iy+G_Iz
           !   call complex_matrix_inverse(G_R,G_Ix,M_Bar)
!               if(myid==0) then
!                open(unit=44,file='H.dat',status='replace')
!               do i=0,M_Bar-1
!                do j=0,M_Bar-1
!                     write(44,*) G_R(i,j),G_Ix(i,j),i,j
!                enddo
!                enddo
!                close(44)
!                endif 

               do i=0,M_Bar-1
                do j=0,M_Bar-1
             q_temp_RR(i,k)=q_temp_RR(i,k)+G_R(i,j)*Pk_real(j,k)
             q_temp_RI(i,k)=q_temp_RI(i,k)+G_R(i,j)*Pk_Imag(j,k)
             q_temp_IR(i,k)=q_temp_IR(i,k)+G_Ix(i,j)*Pk_real(j,k)
             q_temp_II(i,k)=q_temp_II(i,k)+G_Ix(i,j)*Pk_Imag(j,k)
                  enddo
                enddo                    

!             do j=0,M_Bar-1
                ! G_R(j,j)=1.0d0-R_GG_val(j)*ds_temp*kapa_inv
!                  G_R(j,j)=1.0-G_R(j,j)
!                do i=0,M_Bar-1    
!                  G_Ix(i,j)=-G_Ix(i,j)
!                enddo
!              enddo
!              call complex_matrix_inverse(G_R,G_Ix,M_Bar)
!               if(myid==0) then
!                open(unit=44,file='H_inv.dat',status='replace')
!               do i=0,M_Bar-1
!                do j=0,M_Bar-1
!                     write(44,*) G_R(i,j),G_Ix(i,j),i,j
!                enddo
!                enddo
!                close(44)
!                open(unit=44,file='q_temp.dat',status='replace')
!                do j=0,M_Bar-1
!                     write(44,*) q_temp_RR(j,0),q_temp_RI(j,0)
!                     write(44,*) q_temp_IR(j,0),q_temp_II(j,0)
!                     write(44,*) j
!                enddo
!                close(44)
!                endif 


         enddo

# if defined (R2C)
             do i=0,M_Bar-1
                do k=0,LOCAL_SIZE_K_R2C-1
     local_out_c(kDk(k)%x,kDk(k)%y,KDk(k)%z)=dcmplx(q_temp_RR(i,k)-q_temp_II(i,k), & 
                q_temp_RI(i,k)+q_temp_IR(i,k))
                enddo
       !call fftw_mpi_execute_dft(bplan, local_out, local_in)
        call decomp_2d_fft_3d(local_out_c, local_in_r)
            
                do k=0,LOCAL_SIZE-1
                 q_out(i,k)=local_in_r(kD(k)%x,kD(k)%y,KD(k)%z)*M_grid_inv 
                            
                enddo
              enddo  
# else /* R2C */
             do i=0,M_Bar-1
                do k=0,LOCAL_SIZE_K-1
    local_out(kDz(k)%x,kDz(k)%y,KDz(k)%z)=dcmplx(q_temp_RR(i,k)-q_temp_II(i,k), & 
               q_temp_RI(i,k)+q_temp_IR(i,k))
                enddo
       !call fftw_mpi_execute_dft(bplan, local_out, local_in)
        call decomp_2d_fft_3d(local_out, local_in, DECOMP_2D_FFT_BACKWARD)
            
                do k=0,LOCAL_SIZE-1
                 q_out(i,k)=real(local_in(kD(k)%x,kD(k)%y,KD(k)%z))*M_grid_inv + &
                      Pr(i,k) +q_in(i,k)
                enddo
              enddo  
# endif /* R2C */

!write(*,*) "sum(q_out(i,k)",sum(q_out(:,:))

 end subroutine MDE_onestep_q_ex2

 subroutine MDE_onestep_AB_diblock_solver_test(R2C_FFT)
 USE nrtype,only :DP,PI
 USE global_para
 USE mpi
 USE control
 USE constants
 USE utility
 USE mmpi
 use decomp_fft_mod
 implicit none
 logical, optional :: R2C_FFT
 integer :: i,j,s
 integer :: index_nonzero
 integer :: k
 real(DP) :: ds_temp
 character :: MDEmethod
 real(DP) :: q_tmp(0:M_Bar-1,0:Local_size-1)
 real(DP) :: q_err(0:M_Bar-1,0:Local_size-1)
 real(DP) :: q_expf(0:M_Bar-1,0:Local_size-1)
 real(DP) :: q_imf(0:M_Bar-1,0:Local_size-1)
 real(DP) :: q_exps(0:M_Bar-1,0:Local_size-1)
 real(DP) :: q_ims(0:M_Bar-1,0:Local_size-1)
 real*8 :: t1,t2,t3,t4

 ds_temp=ds


do k=0,LOCAL_SIZE-1
   do i=0,M_Bar-1
    qA(i,K,0)=0.0
    if(i==0) qA(i,K,0)=1.0
   enddo
enddo
   q_exps=qA(:,:,0)
   q_ims=qA(:,:,0)


  MDEmethod='M'
           t1=MPI_Wtime() 
do s=1,NA
! maybe qA should be allocated as qA(i,k,s) not qA(s,i,k) ?? needs to be tested later
   call MDE_onestep_q_im(q_ims, q_imf,WA,1,kapa_A_inv, & 
                      ds_temp,MDEmethod, &
                     R_G_val,I_GA1x_val,I_GA1x_col,I_GA1x_row, &
                     I_GA1y_val,I_GA1y_col,I_GA1y_row, &
                     I_GA1z_val,I_GA1z_col,I_GA1z_row, &
                      R2C_FFT)
         q_ims=q_imf
enddo 
           t2=MPI_Wtime() 
  call mp_barrier()
               if(myid==0) then
                open(unit=44,file='q1-im.dat',status='replace')
                do j=0,M_Bar-1
                     write(44,*) q_imf(j,0),j
                enddo
                close(44)
                endif 
  call mp_barrier()

 ds_temp=ds*0.01

           t3=MPI_Wtime() 
do s=1,NA*100
   call MDE_onestep_q_ex2(q_exps, q_expf,q_tmp,WA,1,kapa_A_inv, & 
                      ds_temp,MDEmethod, &
                     R_G_val,I_GA1x_val,I_GA1x_col,I_GA1x_row, &
                     I_GA1y_val,I_GA1y_col,I_GA1y_row, &
                     I_GA1z_val,I_GA1z_col,I_GA1z_row, &
                      R2C_FFT)
         q_exps=q_expf
enddo
           t4=MPI_Wtime() 

  q_err(:,:)=q_expf-q_imf
  write(*,*) "qA1 err",sqrt(sum(q_err(:,:)**2))
  write(*,*) "time",t2-t1,t3-t2
  call mp_barrier()
               if(myid==0) then
                open(unit=44,file='q1-exp.dat',status='replace')
                do j=0,M_Bar-1
                     write(44,*) q_expf(j,0),j
                enddo
                close(44)
                endif 
  call mp_barrier()

  stop


 end subroutine MDE_onestep_AB_diblock_solver_test
end Module wormlike_MDE_single_step
# endif /* EXP */
