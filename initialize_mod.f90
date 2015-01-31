 subroutine initialize()
 USE nrtype,only :DP,PI
 USE global_para
 USE mpi
 USE control
 USE constants
 USE utility
 USE mmpi
 !USE fftw !user defined module 1D-slab or 2d pencil
# if defined (SLAB)
 use fftw3_mpi
 use mpi_fftw3_operation
# else /* SLAB */
 use decomp_fft_mod
# endif /* SLAB */
 implicit none

 logical  :: exists
 integer :: istat1,i,j,k
 integer :: Mem_optimized
 Integer*8 :: planff     		! forward plan for use of fft => xyz
 character(LEN=1) :: dir_grid(3) 
 !logical ::  exists ! why this statement causes error ,find out later!!!

! initialize mpi envirnoment
# if defined (MPI)

! initialize the mpi execution environment
     call mp_init()

! determines the rank of the calling process in the communicator
     call mp_comm_rank(myid)

! determines the size of the group associated with a communicator
     call mp_comm_size(nprocs)

# endif  /* MPI */

# if defined (OPENMP)
call OMP_SET_NUM_THREADS(1)
# endif  /* OPENMP */


call OMP_SET_NUM_THREADS(1)
call mkl_set_num_threads( 1 )


if(myid == master) then
call SCMFT_print_header()
endif

  if(myid==master) then
         exists = .false.

         inquire (file = 'input.txt', exist = exists)

! read in parameters, default setting should be overrided
         if ( exists .eqv. .true. ) then
      open(unit=mytmp,file='input.txt')
read(mytmp,*) Loa             ! L/a
read(mytmp,*) kapa_B          ! kapa_B
read(mytmp,*) kapa_C          ! kapa_B
read(mytmp,*) NMu             ! Maier-Saupe interaction para
read(mytmp,*) Nmax            ! number of discretized bins along the chain for variable s
read(mytmp,*) wABC_init_type   ! wAB initialized type,group symmetry adapted.
read(mytmp,*) iter_method     ! iteration scheme,simple,anderson,or SIS etc.
read(mytmp,*) Anderson_nim_AB     ! iteration scheme,simple,anderson,or SIS etc.
read(mytmp,*) lambda_simple   ! simple mixing parameter if simple mixing is used.
read(mytmp,*) SIDEx           !spatial grid number for x axis
read(mytmp,*) SIDEy    
read(mytmp,*) SIDEz    
read(mytmp,*) Lx             ! side-length of box on x axis 
read(mytmp,*) Ly    
read(mytmp,*) Lz    
read(mytmp,*) x1             ! X1,X2,X3 for 1D brent 
read(mytmp,*) x2    
read(mytmp,*) x3    
read(mytmp,*) NxAB           !XAB*N  flory-Huggins parameter times N kuhn number
read(mytmp,*) NxAC           !XAC*N  flory-Huggins parameter times N kuhn number
read(mytmp,*) NxBC           !XBC*N  flory-Huggins parameter times N kuhn number
read(mytmp,*) fA             ! composition fraction of polymer A
read(mytmp,*) fB             ! composition fraction of polymer A
!read(mytmp,*) p_row         ! P_row for pencil decomposition 
!read(mytmp,*) p_col         ! P_col for pencil decomposition
close(mytmp)
      else 
      write(mystd,*) "input file missing ,program stop!"
      stop
      endif
endif  !endif myid==master


! since input parameters may be updated in master node, it is important
! to broadcast input parameters from root to all children processes

!-----------------------------------------------------------------------
!using the so called function overload scheme generally used in the sofware engineering,
!we have the uniform subroutine mp_bcast for every type of data we want to
! bcast.
# if defined (MPI)
      call mp_barrier()
      call mp_bcast( Loa, master )                                    
      call mp_bcast( kapa_B, master )                                    
      call mp_bcast( kapa_C, master )                                    
      call mp_bcast( NMu, master )                                     
      call mp_bcast( NxAB, master )                                    
      call mp_bcast( NxAC, master )                                    
      call mp_bcast( NxBC, master )                                    
      call mp_bcast( fA, master )                                      
      call mp_bcast( fB, master )                                      
      call mp_bcast( Nmax, master )                                   
      call mp_bcast( iter_method, master )                                
      call mp_bcast( Anderson_nim_AB, master )                                
      call mp_bcast( wABC_init_type, master )                                
      call mp_bcast( lambda_simple, master )                          
      call mp_bcast( SIDEx, master )                                  
      call mp_bcast( SIDEy, master )         
      call mp_bcast( SIDEz, master )         
      call mp_bcast( lx , master )                                     
      call mp_bcast( ly , master )                                    
      call mp_bcast( lz , master )                                  
      call mp_bcast( x1 , master )                                  
      call mp_bcast( x2 , master )                                  
      call mp_bcast( x3 , master )                                  
 !     call mp_bcast( p_row , master )                                  
 !     call mp_bcast( p_col , master )                                  
      call mp_barrier() 
# endif  /* MPI */

# if defined (Debug)
write(*,*) "Loa=",Loa,"kapa_B",kapa_B,"on",myid
write(*,*) "fA=",fA,"on",myid
write(*,*) "SIDEx,SIDEy,SIDEz=",SIDEx,SIDEy,SIDEz,"on",myid
write(*,*) "labntWA,B,M",lanbtwA,lanbtwB,lanbtM,"on",myid
write(*,*) "CC=",CC,"on",myid
# endif  /* Debug */


!set up some parameters: check out !
 ds=1.0d0/Nmax
 fC=1.0d0-fA-fB
 NA=NINT(fA*Nmax)
 NB=NINT(fB*Nmax)
 NC=Nmax-NA-NB
 ! kapa_A is a parameter 
 kapa_B_inv=1.0d0/kapa_B
 kapa_C_inv=1.0d0/kapa_C

     chi_blk(1,1)=0.0d0
     chi_blk(1,2)=NxAB
     chi_blk(1,3)=NxAC
     chi_blk(2,1)=NxAB
     chi_blk(2,2)=0.0d0
     chi_blk(2,3)=NxBC
     chi_blk(3,1)=NxAC
     chi_blk(3,2)=NxBC
     chi_blk(3,3)=0.0d0

      f_blk(1)=fA
      f_blk(2)=fB
      f_blk(3)=fC

      Ns_blk(1)=NA 
      Ns_blk(2)=NB 
      Ns_blk(3)=NC 
      Nstart_blk(1)=0 
      Nstart_blk(2)=NA 
      Nstart_blk(3)=NA+NB 

      kapa_blk(1)=kapa_A 
      kapa_blk(2)=kapa_B
      kapa_blk(3)=kapa_C 

      L_bar_blk(1)=L_bar_A
      L_bar_blk(2)=L_bar_B
      L_bar_blk(3)=L_bar_C
      M_Bar_A=(L_Bar_A+1)**2
      M_Bar_B=(L_Bar_B+1)**2
      M_Bar_C=(L_Bar_C+1)**2
      M_bar_blk(1)=M_bar_A
      M_bar_blk(2)=M_bar_B
      M_bar_blk(3)=M_bar_C

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1      
 x1=x1/SIDEx
 x2=x2/SIDEx
 x3=x3/SIDEx

if(myid==0) then
 if(mod(NA,2)==0 .and. mod(NB,2)==0 .and. mod(NC,2)==0) then
  write(*,*) "NA,NB,NC checked!",NA,NB,NC
  else
  write(*,*) "NA,NB,NC=",NA,NB,NC
  stop
  endif
endif

 M_initial=NMu*0.1
 dx_erro=1.0

!! 3D calcualtion by default,2D if -DDim2 option added for fftw slab
!! for pencil fft, 3D calculation only.
! set grid infomation using module grid

# if defined (SLAB)
# if defined (Dim2)  
 M_grid=SIDEx*SIDEy
!call grid_get_info(2,dir_grid)
# else  /* Dim2 */
 M_grid=SIDEx*SIDEy*SIDEz
!call grid_get_info(3,dir_grid)
# endif /* Dim2 */
# else /* SLAB */
! this means if pencil FFT is used, the dimension has be be 3D
 M_grid=SIDEx*SIDEy*SIDEz
!call grid_get_info(3,dir_grid)
# endif /* SLAB */


 if(abs(NMu)<1.0e-5) then
  NMu_add=.false.
 else
  NMu_add=.true.
 endif


# if defined (SLAB)
# if defined (Dim2)  
      call fftw3_mpi_2d_dft_init()
      local_size=local_nx*SiDEy
      write(*,*) "local_size=",local_size,"on",myid
      !call global_grid_decomposition(1) 
      call init_k_index()
# else  /* Dim2 */
      call fftw3_mpi_3d_dft_init()
      local_size=local_nx*SiDEy*SIDEz
      write(*,*) "local_size=",local_size,"on",myid
      !call global_grid_decomposition(1) 
      call init_k_index()
# endif  /* Dim2 */

# else /* SLAB */

# if defined (R2C)
call decomp2d_mpi_3d_dft_r2c_init()
!if initialized as X pencil array
      Local_size=xsize(1)*xsize(2)*xsize(3)
      Local_size_k=zsize(1)*zsize(2)*zsize(3)
      write(*,*) "local_size=",local_size,"on",myid
      write(*,*) "local_size_k=",local_size_k,"on",myid
      write(*,*) "local_size_k_r2c=",local_size_k_r2c,"on",myid
! be careful, when multidimensional R2C FFT is performed, 
! the array size of the complex output is only about half the size
! of the real input array, more detail can be found on http://www.2decomp.org/fft_api.html 

     ! call global_grid_decomposition(2,1) 
      call init_local2global_index(1,.true.)
# else /* R2C */
      call decomp2d_mpi_3d_dft_c2c_init()
!if initialized as X pencil array
      Local_size=xsize(1)*xsize(2)*xsize(3)
      Local_size_k=zsize(1)*zsize(2)*zsize(3)
      write(*,*) "local_size=",local_size,"on",myid
      write(*,*) "local_size_k=",local_size_k,"on",myid
    !  call global_grid_decomposition(2,1) 
      call init_local2global_index(1)
# endif /* R2C */

# endif /* SLAB */

      call mp_barrier()
!!!!***************************begin to intialize w and scft variables**
!INIT needed arrays for different methods

# if defined (SLAB)
call init_arrays_BDF()
call k_index()
# else /* SLAB */
!call init_arrays_2decomp_BDF()
call init_arrays_2decomp()
# endif /* SLAB */
! init w field form input or generate itself
if(wAB_init_type<0) then
call readW()
else
call init_w()
endif

call mp_barrier()
if(myid==0) then
write(*,*) "done calling initw"
endif
call mp_barrier()
!call calc_THETAij()

!!! init the dump file 
!***************************************************************
call mp_barrier()
if(myid==0) then
open(unit=48,file='FreeEG.dat',status='new')
close(48)
open(unit=43,file='error.dat',status='new')
close(43)
open(unit=42,file='result.txt',status='new')
close(42)
endif

call mp_barrier()
!***********************************************************
!checkout the input paras
# if defined (Debug)
if(myid==0) then
write(*,*) "wABC_init_type=",wABC_init_type,"on",myid
write(*,*) "M_initial=",M_initial,"on",myid
write(*,*) "Method=",Method,"on",myid
write(*,*) "NMu=",NMu,"on",myid
write(*,*) "SIDEy,SIZEz",SIDEy,SIDEz,"on",myid
write(*,*) "local_nx=",local_nx,"on",myid
write(*,*) "local_size=",local_size,"on",myid
write(*,*) "x1,x2,x3",x1,x2,x3,"on",myid
endif
# endif  /* Debug */ 

call mp_barrier()

end subroutine initialize


!!!!test region ,remove after it is done
# if defined (SLAB)
subroutine init_arrays_BDF()
 USE nrtype,only :DP
 USE global_para
 USE mpi
 USE control
 USE constants
 USE utility
 USE mmpi
 USE mpi_fftw3_operation
 !USE fftw !user defined module
 USE SPH
 implicit none
 integer :: istat2,i,j,k
 REAL(DP) :: temp,temp1,temp2,sumTHETA
 integer :: LS,MB,MA
 integer :: n_anderson
 integer :: nonzero_counter

 LS=LOCAL_SIZE-1
 max_M_Bar=max(M_Bar_blk(:))

 if(myid==0) then
 write(*,*) "MA=,MB=,MC=",M_Bar_blk
 endif
 call mp_barrier()
 !*****************note that N_dim_ddm==3!
  
# if defined (Dim2)
 allocate(matrix_Rx(0:max_M_Bar-1,0:max_M_Bar-1),stat=istat2)
 allocate(matrix_Ry(0:max_M_Bar-1,0:max_M_Bar-1),stat=istat2)
# else /* Dim2 */
 allocate(matrix_Rx(0:max_M_Bar-1,0:max_M_Bar-1),stat=istat2)
 allocate(matrix_Ry(0:max_M_Bar-1,0:max_M_Bar-1),stat=istat2)
 allocate(matrix_Rz(0:max_M_Bar-1,0:max_M_Bar-1),stat=istat2)
# endif  /* Dim2 */

 allocate(GAMA(0:max_M_Bar-1,0:max_M_Bar-1,0:max_M_Bar-1),stat=istat2)

 allocate(q_f(0:max_M_Bar-1,0:LS,0:Nmax),stat=istat2)
 allocate(q_star(0:max_M_Bar-1,0:LS,0:Nmax),stat=istat2)

! allocate(qA(0:MA,0:LS,0:NA),stat=istat2)
! allocate(qB(0:MB,0:LS,0:NB),stat=istat2)
 !allocate(qA(0:NA,0:MA,0:LS),stat=istat2)
 !allocate(qB(0:NB,0:MB,0:LS),stat=istat2)

 !allocate(qAstar(0:MA,0:LS,0:NA),stat=istat2)
 !allocate(qBstar(0:MB,0:LS,0:NB),stat=istat2)
 !allocate(qAstar(0:NA,0:MA,0:LS),stat=istat2)
 !allocate(qBstar(0:NB,0:MB,0:LS),stat=istat2)
 !allocate(wA(0:LS),stat=istat2)
 !allocate(wB(0:LS),stat=istat2)
  allocate(w_blk(0:LS,1:3),stat=istat2)
  allocate(w_out_blk(0:LS,0:n_r_WAB,1:3),stat=istat2)
  allocate(dw_blk(0:LS,0:n_r_WAB,1:3),stat=istat2)
  
 !allocate(wA_out(0:LS,0:n_r_WAB),stat=istat2)
 !allocate(wB_out(0:LS,0:n_r_WAB),stat=istat2)
 !allocate(dw_A(0:LS,0:n_r_WAB),stat=istat2)
 !allocate(dw_B(0:LS,0:n_r_WAB),stat=istat2)

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


 !allocate(RA(0:LS),stat=istat2)
 !allocate(RB(0:LS),stat=istat2)
 allocate(R_blk(0:LS,1:3),stat=istat2)
 allocate(R_joint(0:LS,1:3),stat=istat2)
 allocate(RHO_blk(0:MA,0:LS,1:3),stat=istat2)
 !allocate(R_half(0:LS),stat=istat2)
 !allocate(R_end(0:LS),stat=istat2)

 !allocate(RHOA(0:MA,0:LS),stat=istat2)
 !allocate(RHOB(0:MB,0:LS),stat=istat2)

 allocate(M_OP(0:LS,0:N_dim_ddm-1,0:N_dim_ddm-1),stat=istat2)
 allocate(S_OP(0:LS,0:N_dim_ddm-1,0:N_dim_ddm-1),stat=istat2)
 allocate(S_OP_blk(0:LS,0:N_dim_ddm-1,0:N_dim_ddm-1,1:3),stat=istat2)
 !allocate(SA_OP(0:LS,0:N_dim_ddm-1,0:N_dim_ddm-1),stat=istat2)
 !allocate(SB_OP(0:LS,0:N_dim_ddm-1,0:N_dim_ddm-1),stat=istat2)

 allocate(Pu(0:LS,1:N_dim_ddm),stat=istat2)
 allocate(Pu_blk(0:LS,1:N_dim_ddm,1:3),stat=istat2)
 !allocate(Pu_A(0:LS,1:N_dim_ddm),stat=istat2)
 !allocate(Pu_B(0:LS,1:N_dim_ddm),stat=istat2)

 allocate(THETAij_blk(0:max_M_Bar-1,0:max_M_Bar-1,0:LS,1:3),stat=istat2)
 allocate(THETAij_M11_M22_blk(0:max_M_Bar-1,0:max_M_Bar-1,1:3),stat=istat2)
 allocate(THETAij_M33_blk(0:max_M_Bar-1,0:max_M_Bar-1,1:3),stat=istat2)
 allocate(THETAij_M12_blk(0:max_M_Bar-1,0:max_M_Bar-1,1:3),stat=istat2)
 allocate(THETAij_M13_blk(0:max_M_Bar-1,0:max_M_Bar-1,1:3),stat=istat2)
 allocate(THETAij_M23_blk(0:max_M_Bar-1,0:max_M_Bar-1,1:3),stat=istat2)

 allocate(J11ij_blk(0:max_M_Bar-1,0:max_M_Bar-1,1:3),stat=istat2)
 allocate(J22ij_blk(0:max_M_Bar-1,0:max_M_Bar-1,1:3),stat=istat2)
 allocate(J12ij_blk(0:max_M_Bar-1,0:max_M_Bar-1,1:3),stat=istat2)
 allocate(J13ij_blk(0:max_M_Bar-1,0:max_M_Bar-1,1:3),stat=istat2)
 allocate(J23ij_blk(0:max_M_Bar-1,0:max_M_Bar-1,1:3),stat=istat2)
! allocate(THETAij_A(0:MA,0:MA,0:LS),stat=istat2)
! allocate(THETAij_M11_M22_A(0:MA,0:MA),stat=istat2)
! allocate(THETAij_M33_A(0:MA,0:MA),stat=istat2)
! allocate(THETAij_M12_A(0:MA,0:MA),stat=istat2)
! allocate(THETAij_M13_A(0:MA,0:MA),stat=istat2)
! allocate(THETAij_M23_A(0:MA,0:MA),stat=istat2)

! allocate(J11ij_A(0:MA,0:MA),stat=istat2)
! allocate(J22ij_A(0:MA,0:MA),stat=istat2)
! allocate(J12ij_A(0:MA,0:MA),stat=istat2)
! allocate(J13ij_A(0:MA,0:MA),stat=istat2)
! allocate(J23ij_A(0:MA,0:MA),stat=istat2)

!! note that THETAij_B could be optimized for memeory ,this is soon to de done
! allocate(THETAij_B(0:MB,0:MB,0:LS),stat=istat2)
! allocate(THETAij_M11_M22_B(0:MB,0:MB),stat=istat2)
! allocate(THETAij_M33_B(0:MB,0:MB),stat=istat2)
! allocate(THETAij_M12_B(0:MB,0:MB),stat=istat2)
! allocate(THETAij_M13_B(0:MB,0:MB),stat=istat2)
! allocate(THETAij_M23_B(0:MB,0:MB),stat=istat2)
!
! allocate(J11ij_B(0:MB,0:MB),stat=istat2)
! allocate(J22ij_B(0:MB,0:MB),stat=istat2)
! allocate(J12ij_B(0:MB,0:MB),stat=istat2)
! allocate(J13ij_B(0:MB,0:MB),stat=istat2)
! allocate(J23ij_B(0:MB,0:MB),stat=istat2)
 

allocate(THETA_nonzero_1D_indexK_blk(-1:LOCAL_SIZE-1,1:3),stat=istat2)

!allocate(THETA_nonzero_1D_A_indexK(-1:LOCAL_SIZE-1),stat=istat2)
!allocate(THETA_nonzero_1D_B_indexK(-1:LOCAL_SIZE-1),stat=istat2)

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

!!!!for matrix RX,RY,RZ
# if defined (Dim2)
temp=-1.0/sqrt(3.0)
  do j=0,max_M_Bar-1
     do i=0,max_M_Bar-1
     matrix_Rx(i,j)=temp*Triple_product(basis_SPH(i)%x,basis_SPH(j)%x,1,basis_SPH(i)%y, &
                    basis_SPH(j)%y,1)
     matrix_Ry(i,j)=temp*Triple_product(basis_SPH(i)%x,basis_SPH(j)%x,1,basis_SPH(i)%y, &
                    basis_SPH(j)%y,-1)
 !  if(myid==0) then
 !    if(matrix_Rx(i,j)>1.0e-10) then
 !     write(*,*) "matrix_Rx",matrix_Rx(i,j),"i",i,"j",j
 !    endif
 ! endif
    enddo
enddo

!write(*,*) "basis_SPH(1)%x,y=",basis_SPH(1)%x,basis_SPH(1)%y
# else  /* Dim2 */

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
# endif  /* Dim2 */

if(myid==0) then
write(*,*) "marix_Rx(1,7)=",matrix_Rx(1,7)
write(*,*) "marix_Rx(7,1)=",matrix_Rx(7,1)
endif 
!!! THETAij_A 

temp=1.0/sqrt(15.0)
  do j=0,max_M_Bar-1
     do i=0,max_M_Bar-1
        THETAij_M11_M22_blk(i,j,1)=temp*Triple_product(basis_SPH(i)%x,basis_SPH(j)%x, &
                                2,basis_SPH(i)%y,basis_SPH(j)%y,2)
     enddo
  enddo

temp=1.0/sqrt(5.0)
  do j=0,max_M_Bar-1
     do i=0,max_M_Bar-1
        THETAij_M33_blk(i,j,1)=temp*Triple_product(basis_SPH(i)%x,basis_SPH(j)%x, &
                                2,basis_SPH(i)%y,basis_SPH(j)%y,0)
     enddo
  enddo

temp=2.0/sqrt(15.0)
  do j=0,max_M_Bar-1
     do i=0,max_M_Bar-1
        THETAij_M12_blk(i,j,1)=temp*Triple_product(basis_SPH(i)%x,basis_SPH(j)%x, &
                                2,basis_SPH(i)%y,basis_SPH(j)%y,-2)
     enddo
  enddo
 

temp=-2.0/sqrt(15.0)
  do j=0,max_M_Bar-1
     do i=0,max_M_Bar-1
        THETAij_M13_blk(i,j,1)=temp*Triple_product(basis_SPH(i)%x,basis_SPH(j)%x, &
                                2,basis_SPH(i)%y,basis_SPH(j)%y,1)
     enddo
  enddo

temp=-2.0/sqrt(15.0)
  do j=0,max_M_Bar-1
     do i=0,max_M_Bar-1
        THETAij_M23_blk(i,j,1)=temp*Triple_product(basis_SPH(i)%x,basis_SPH(j)%x, &
                                2,basis_SPH(i)%y,basis_SPH(j)%y,-1)
     enddo
  enddo


temp1=1.0/sqrt(15.0)
temp2=-1.0/(3.0*sqrt(5.0))
  do j=0,max_M_Bar-1
     do i=0,max_M_Bar-1
         J11ij_blk(i,j,1)=temp1*triple_product(basis_SPH(i)%x,basis_SPH(j)%x,2, &
                           basis_SPH(i)%y,basis_SPH(j)%y,2) + &
                      temp2*triple_product(basis_SPH(i)%x,basis_SPH(j)%x,2, &
                           basis_SPH(i)%y,basis_SPH(j)%y,0)
      enddo
   enddo


temp1=-1.0/sqrt(15.0)
temp2=-1.0/(3.0*sqrt(5.0))
  do j=0,max_M_Bar-1
     do i=0,max_M_Bar-1
         J22ij_blk(i,j,1)=temp1*triple_product(basis_SPH(i)%x,basis_SPH(j)%x,2, &
                           basis_SPH(i)%y,basis_SPH(j)%y,2) + &
                      temp2*triple_product(basis_SPH(i)%x,basis_SPH(j)%x,2, &
                           basis_SPH(i)%y,basis_SPH(j)%y,0)
      enddo
   enddo

temp1=1.0/sqrt(15.0)
  do j=0,max_M_Bar-1
     do i=0,max_M_Bar-1
         J12ij_blk(i,j,1)=temp1*triple_product(basis_SPH(i)%x,basis_SPH(j)%x,2, &
                           basis_SPH(i)%y,basis_SPH(j)%y,-2)
      enddo
   enddo


temp1=-1.0/sqrt(15.0)
  do j=0,max_M_Bar-1
     do i=0,max_M_Bar-1
         J13ij_blk(i,j,1)=temp1*triple_product(basis_SPH(i)%x,basis_SPH(j)%x,2, &
                           basis_SPH(i)%y,basis_SPH(j)%y,1)
      enddo
   enddo


temp1=-1.0/sqrt(15.0)
  do j=0,max_M_Bar-1
     do i=0,max_M_Bar-1
         J23ij_blk(i,j,1)=temp1*triple_product(basis_SPH(i)%x,basis_SPH(j)%x,2, &
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
          
# if defined (Dim2)

# else  /* Dim2 */
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
# endif  /* Dim2 */


          
!!initialize J??ij_nonzero_1D_A
 nonzero_counter=0
 do j=0,max_M_Bar-1
    do i=0,max_M_bar-1
          if(abs(J11ij_blk(i,j,1))>=Thresh_sprase_matrix) then
              nonzero_counter=nonzero_counter+1
          endif
     enddo
  enddo
       
# if defined (Debug)
 if(myid==0) then
 write(*,*) "nonzero elements of J11ij_A=",nonzero_counter
 endif
# endif  /* Debug */

  if(.not. allocated(J11ij_nonzero_1D_blk)) then
  allocate(J11ij_nonzero_1D_blk(0:nonzero_counter-1,1:1))
  endif
  
 nonzero_counter=0
 do j=0,max_M_Bar-1
    do i=0,max_M_bar-1
          if(abs(J11ij_blk(i,j,1))>=Thresh_sprase_matrix) then
              J11ij_nonzero_1D_blk(nonzero_counter,1)%i=i
              J11ij_nonzero_1D_blk(nonzero_counter,1)%j=j
              J11ij_nonzero_1D_blk(nonzero_counter,1)%value=J11ij_blk(i,j,1)
              nonzero_counter=nonzero_counter+1
           endif
      enddo
  enddo
        
 J11_A=nonzero_counter
          
 nonzero_counter=0
 do j=0,max_M_Bar-1
    do i=0,max_M_bar-1
          if(abs(J22ij_blk(i,j,1))>=Thresh_sprase_matrix) then
              nonzero_counter=nonzero_counter+1
          endif
     enddo
  enddo
       
  if(.not. allocated(J22ij_nonzero_1D_blk)) then
  allocate(J22ij_nonzero_1D_blk(0:nonzero_counter-1,1:1))
  endif
  
# if defined (Debug)
 if(myid==0) then
 write(*,*) "nonzero elements of J22ij_A=",nonzero_counter
 endif
# endif  /* Debug */

 nonzero_counter=0
 do j=0,max_M_Bar-1
    do i=0,max_M_bar-1
          if(abs(J22ij_blk(i,j,1))>=Thresh_sprase_matrix) then
              J22ij_nonzero_1D_blk(nonzero_counter,1)%i=i
              J22ij_nonzero_1D_blk(nonzero_counter,1)%j=j
              J22ij_nonzero_1D_blk(nonzero_counter,1)%value=J22ij_blk(i,j,1)
              nonzero_counter=nonzero_counter+1
           endif
      enddo
  enddo

 J22_A=nonzero_counter
              
 nonzero_counter=0
 do j=0,max_M_Bar-1
    do i=0,max_M_bar-1
          if(abs(J12ij_blk(i,j,1))>=Thresh_sprase_matrix) then
              nonzero_counter=nonzero_counter+1
          endif
     enddo
  enddo
       
  if(.not. allocated(J12ij_nonzero_1D_blk)) then
  allocate(J12ij_nonzero_1D_blk(0:nonzero_counter-1,1:1))
  endif
  
 J12_A=nonzero_counter
 nonzero_counter=0
 do j=0,max_M_Bar-1
    do i=0,max_M_bar-1
          if(abs(J12ij_blk(i,j,1))>=Thresh_sprase_matrix) then
              J12ij_nonzero_1D_blk(nonzero_counter,1)%i=i
              J12ij_nonzero_1D_blk(nonzero_counter,1)%j=j
              J12ij_nonzero_1D_blk(nonzero_counter,1)%value=J12ij_blk(i,j,1)
              nonzero_counter=nonzero_counter+1
           endif
      enddo
  enddo


 nonzero_counter=0
 do j=0,max_M_Bar-1
    do i=0,max_M_bar-1
          if(abs(J13ij_blk(i,j,1))>=Thresh_sprase_matrix) then
              nonzero_counter=nonzero_counter+1
          endif
     enddo
  enddo
       
  if(.not. allocated(J13ij_nonzero_1D_blk)) then
  allocate(J13ij_nonzero_1D_blk(0:nonzero_counter-1,1:1))
  endif
  
 nonzero_counter=0
 do j=0,max_M_Bar-1
    do i=0,max_M_bar-1
          if(abs(J13ij_blk(i,j,1))>=Thresh_sprase_matrix) then
              J13ij_nonzero_1D_blk(nonzero_counter,1)%i=i
              J13ij_nonzero_1D_blk(nonzero_counter,1)%j=j
              J13ij_nonzero_1D_blk(nonzero_counter,1)%value=J13ij_blk(i,j,1)
              nonzero_counter=nonzero_counter+1
           endif
      enddo
  enddo

 J13_A=nonzero_counter
 nonzero_counter=0
 do j=0,max_M_Bar-1
    do i=0,max_M_bar-1
          if(abs(J23ij_blk(i,j,1))>=Thresh_sprase_matrix) then
              nonzero_counter=nonzero_counter+1
          endif
     enddo
  enddo
       
  if(.not. allocated(J23ij_nonzero_1D_blk)) then
  allocate(J23ij_nonzero_1D_blk(0:nonzero_counter-1,1:1))
  endif
  
 nonzero_counter=0
 do j=0,max_M_Bar-1
    do i=0,max_M_bar-1
          if(abs(J23ij_blk(i,j,1))>=Thresh_sprase_matrix) then
              J23ij_nonzero_1D_blk(nonzero_counter,1)%i=i
              J23ij_nonzero_1D_blk(nonzero_counter,1)%j=j
              J23ij_nonzero_1D_blk(nonzero_counter,1)%value=J23ij_blk(i,j,1)
              nonzero_counter=nonzero_counter+1
           endif
      enddo
  enddo


 J23_A=nonzero_counter

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

end subroutine init_arrays_BDF
          
!!!!!!!!!!!!!!!!!!!initialize the input field!!!!!!!!!!!!!!!!!!

       subroutine init_w()
          USE nrtype,only :PI
         USE global_para
         USE mpi
         USE control
         USE constants
         USE utility
         USE mmpi
         USE mpi_fftw3_operation
         implicit none
         integer :: k_i,k_j,k_k
         integer :: KK,j,k,i
         integer :: initseed
         integer :: dt(8)
         
   if(wAB_init_type==0) then
          !random init
          call date_and_time(values=dt)
          initseed=(dt(8)-500)*12345*(myid+1)+43210
          
          KK=0  
# if defined (Dim2)     
          do K_k=0,Local_nx-1
             do k_j=0,SIDEy-1
                   WA(KK)=1.0*ran2(initseed)
                   WB(KK)=1.0*ran2(initseed)
                   KK=KK+1
             enddo
          enddo
# else  /* Dim2 */
         do K_k=0,Local_nx-1
             do k_j=0,SIDEy-1
               do k_i=0,SIDEz-1
                   WA(KK)=1.0*ran2(initseed)
                   WB(KK)=1.0*ran2(initseed)
                   KK=KK+1
             enddo
          enddo
        enddo
# endif  /* Dim2 */
        

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


# if defined (Dim2)
   else if(wAB_init_type==4) then
      !!!Cylinder initial 3D
     KK=0
          do K_k=0,Local_nx-1
             do k_j=0,SIDEy-1
                   WA(KK)=NXab*(1-fA*(1+0.8*cos(2.0*PI*(local_x_start+K_k+1)/SIDEx)* &
                         cos(2.0*PI*(K_j+1)/SIDEy)))

                   WB(KK)=NXab*fA*(1+0.8*cos(2.0*PI*(local_x_start+K_k+1)/SIDEx)* &
                         cos(2.0*PI*(K_j+1)/SIDEy))
                   KK=KK+1
             enddo
          enddo


   else if(wAB_init_type==5) then
      !!!Lamellar initial 3D
     KK=0
          do K_k=0,Local_nx-1
             do k_j=0,SIDEy-1
                   WA(KK)=NXab*(1-fA*(1+0.8*sin(2.0*PI*(K_j+1)/SIDEy)))

                   WB(KK)=NXab*fA*(1+0.8*sin(2.0*PI*(K_j+1)/SIDEy))
                   KK=KK+1
             enddo
          enddo


   else if(wAB_init_type==6) then
      !!!P4 initial 3D
     KK=0
          do K_k=0,Local_nx-1
             do k_j=0,SIDEy-1
             WA(KK)=NXab*(1-fA*(cos(2.0*PI*(K_j+1)/(1.0*SIDEy))+ &
              cos(2.0*PI*(local_x_start+K_k+1)/(1.0*SIDEx))))
              WB(KK)=NXab*fA*(cos(2.0*PI*(K_j+1)/(1.0*SIDEy))+ & 
              cos(2.0*PI*(local_x_start+K_k+1)/(1.0*SIDEx)))
                   KK=KK+1
             enddo
          enddo

# else  /* Dim2 */

   else if(wAB_init_type==1) then   
      !!!FCC initial 3D
     KK=0
          do K_k=0,Local_nx-1
             do k_j=0,SIDEy-1
                 do k_i=0,SIDEz-1
                   WA(KK)=NXab*(1-fA*(1+0.7*(cos(2.0*PI*(local_x_start+K_k+1)/SIDEx)* &
                         cos(2.0*PI*(K_j+1)/SIDEy)*cos(2.0*PI*(K_i+1)/SIDEz))))
                   
                   WB(KK)=NXab*fA*(1+0.7*(cos(2.0*PI*(local_x_start+K_k+1)/SIDEx)* &
                         cos(2.0*PI*(K_j+1)/SIDEy)*cos(2.0*PI*(K_i+1)/SIDEz)))
                   KK=KK+1
                  enddo
             enddo
          enddo
   
                
    !!! BCC initial 3D
   else if(wAB_init_type==2) then
      
     KK=0
          do K_k=0,Local_nx-1
             do k_j=0,SIDEy-1
                 do k_i=0,SIDEz-1
         WA(KK)=NXab*(1-fA*(1+0.7*(cos(2.0*PI*(local_x_start+K_k+1)/SIDEx)* &
         cos(2.0*PI*(K_j+1)/SIDEy) + cos(2.0*PI*(K_j+1)/SIDEy)*cos(2.0*PI*(k_i+1)/SIDEz) &
          + cos(2.0*PI*(local_x_start+k_k+1)/SIDEx)*cos(2.0*PI*(K_i+1)/SIDEz))))

         WB(KK)=NXab*fA*(1+0.7*(cos(2.0*PI*(local_x_start+K_k+1)/SIDEx)* &
         cos(2.0*PI*(K_j+1)/SIDEy) + cos(2.0*PI*(K_j+1)/SIDEy)*cos(2.0*PI*(k_i+1)/SIDEz) &
         + cos(2.0*PI*(local_x_start+k_k+1)/SIDEx)*cos(2.0*PI*(K_i+1)/SIDEz)))

                   KK=KK+1
                  enddo
             enddo
          enddo
   
                
      !Gyroid initial 3D
   else if(wAB_init_type==3) then
      
     KK=0
          do K_k=0,Local_nx-1
             do k_j=0,SIDEy-1
                 do k_i=0,SIDEz-1
                   WA(KK)=NXab*(1-fA*(1+0.7*(cos(2.0*PI*(local_x_start+K_k+1)/SIDEx)* &
                   sin(2.0*PI*(K_j+1)/SIDEy)*sin(4.0*PI*(k_i+1)/SIDEz) &
                   + cos(2.0*PI*(K_j+1)/SIDEy)*sin(2.0*PI*(k_i+1)/SIDEz)* &
                     sin(4.0*Pi*(local_x_start+k_k+1)/SIDEx) &
                   + cos(2.0*PI*(K_i+1)/SIDEz)*sin(2.0*Pi*(local_x_start+k_k+1)/SIDEx)* &
                     sin(4.0*PI*(k_j+1)/SIDEy))))

                 WB(KK)=NXab*fA*(1+0.7*(cos(2.0*PI*(local_x_start+K_k+1)/SIDEx)* &
                   sin(2.0*PI*(K_j+1)/SIDEy)*sin(4.0*PI*(k_i+1)/SIDEz) &
                   + cos(2.0*PI*(K_j+1)/SIDEy)*sin(2.0*PI*(k_i+1)/SIDEz)* &
                     sin(4.0*Pi*(local_x_start+k_k+1)/SIDEx) &
                   + cos(2.0*PI*(K_i+1)/SIDEz)*sin(2.0*Pi*(local_x_start+k_k+1)/SIDEx)* &
                     sin(4.0*PI*(k_j+1)/SIDEy)))

                   KK=KK+1
                  enddo
             enddo
          enddo

   else if(wAB_init_type==4) then
      !!!Cylinder initial 3D
     KK=0
          do K_k=0,Local_nx-1
             do k_j=0,SIDEy-1
                 do k_i=0,SIDEz-1
                   WA(KK)=NXab*(1-fA*(1+0.8*cos(2.0*PI*(local_x_start+K_k+1)/SIDEx)* &
                         cos(2.0*PI*(K_j+1)/SIDEy)))

                   WB(KK)=NXab*fA*(1+0.8*cos(2.0*PI*(local_x_start+K_k+1)/SIDEx)* &
                         cos(2.0*PI*(K_j+1)/SIDEy))
                   KK=KK+1
                  enddo
             enddo
          enddo


   else if(wAB_init_type==5) then
      !!!Lamellar initial 3D
     KK=0
          do K_k=0,Local_nx-1
             do k_j=0,SIDEy-1
                 do k_i=0,SIDEz-1
                   WA(KK)=NXab*(1-fA*(1+0.8*sin(2.0*PI*(K_i+1)/SIDEz)))

                   WB(KK)=NXab*fA*(1+0.8*sin(2.0*PI*(K_i+1)/SIDEz))
                   KK=KK+1
                  enddo
             enddo
          enddo


   else if(wAB_init_type==6) then
      !!!P4 initial 3D
     KK=0
          do K_k=0,Local_nx-1
             do k_j=0,SIDEy-1
                 do k_i=0,SIDEz-1
             WA(KK)=NXab*(1-fA*(cos(2.0*PI*(K_j+1)/(1.0*SIDEy))+ & 
               cos(2.0*PI*(local_x_start+K_k+1)/(1.0*SIDEx))))
             WB(KK)=NXab*fA*(cos(2.0*PI*(K_j+1)/(1.0*SIDEy))+ & 
                 cos(2.0*PI*(local_x_start+K_k+1)/(1.0*SIDEx)))
                   KK=KK+1
                  enddo
             enddo
          enddo

# endif /* Dim2 */

endif


!!!!!!!!!!!!!!!!!!!!!for M ,not randomly initialized


          KK=0   
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
               M_OP(k,j,i)=M_initial*S_OP(k,j,i)
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
# if defined (SLAB)
 use fftw3_mpi
 use mpi_fftw3_operation
# else /* SLAB */
 use decomp_fft_mod
# endif /* SLAB */

  implicit none
  integer :: i,j,k,aaa
  integer :: k_k,k_j,k_i
  REAL*8 :: x,y,z,sumaW,sumaWtot
  REAL*8 :: sumM00,sumM01,sumM02,sumM10,sumM11,sumM12,sumM20,sumM21,sumM22,sumMtot
  character(len=30)::aa

  double precision, allocatable :: WA_global(:,:,:) 
  double precision, allocatable :: WB_global(:,:,:) 

  double precision, allocatable :: M00(:,:,:),M01(:,:,:),M02(:,:,:)
  double precision, allocatable :: M10(:,:,:),M11(:,:,:),M12(:,:,:)
  double precision, allocatable :: M20(:,:,:),M21(:,:,:),M22(:,:,:)
  



      allocate(WA_global(1:SIDEx,1:SIDEy,1:SIDEz))
      allocate(WB_global(1:SIDEx,1:SIDEy,1:SIDEz))

      allocate(M00(1:SIDEx,1:SIDEy,1:SIDEz))
      allocate(M01(1:SIDEx,1:SIDEy,1:SIDEz))
      allocate(M02(1:SIDEx,1:SIDEy,1:SIDEz))
      allocate(M10(1:SIDEx,1:SIDEy,1:SIDEz))
      allocate(M11(1:SIDEx,1:SIDEy,1:SIDEz))
      allocate(M12(1:SIDEx,1:SIDEy,1:SIDEz))
      allocate(M20(1:SIDEx,1:SIDEy,1:SIDEz))
      allocate(M21(1:SIDEx,1:SIDEy,1:SIDEz))
      allocate(M22(1:SIDEx,1:SIDEy,1:SIDEz))



!read input field file : the old way.
!  aaa=myid
!  ii=myid
!  write(aa,*) aaa
!  open(unit=23,file= 'W' //trim(adjustl(aa)) // '.dat',status='old')
!  do K=0,LOCAL_SIZE-1
!  read(23,*) i,WA(K),WB(K)
!  
!  enddo
!  close(23)
!
!  open(unit=23,file= 'M_OP' //trim(adjustl(aa)) // '.dat',status='old')
!  do K=0,LOCAL_SIZE-1
!  read(23,*) x,y,z,M_OP(K,0,0),M_OP(K,0,1),M_OP(K,0,2)
!  read(23,*) x,y,z,M_OP(K,1,0),M_OP(K,1,1),M_OP(K,1,2)
!  read(23,*) x,y,z,M_OP(K,2,0),M_OP(K,2,1),M_OP(K,2,2)
!  
!  enddo
!  close(23)

if(myid==master) then
  open(unit=23,file= 'W_global.dat',status='old')
        k=0
        do K_i=1,SIDEz
         do k_j=1,SIDEy
          do K_k=1,SIDEx
          read(23,*) x,y,z,WA_global(k_k,k_j,k_i), &
                           WB_global(k_k,k_j,k_i)
          k=k+1
          enddo
         enddo
        enddo
  close(23)

  open(unit=23,file= 'M0_global.dat',status='old')
        k=0
        do K_i=1,SIDEz
         do k_j=1,SIDEy
          do K_k=1,SIDEx
          read(23,*) x,y,z,M00(k_k,k_j,k_i),M01(k_k,k_j,k_i),M02(k_k,k_j,k_i) 
          k=k+1
          enddo
         enddo
        enddo
  close(23)
  open(unit=23,file= 'M1_global.dat',status='old')
        k=0
        do K_i=1,SIDEz
         do k_j=1,SIDEy
          do K_k=1,SIDEx
          read(23,*) x,y,z,M10(k_k,k_j,k_i),M11(k_k,k_j,k_i),M12(k_k,k_j,k_i) 
          k=k+1
          enddo
         enddo
        enddo
  close(23)
  open(unit=23,file= 'M2_global.dat',status='old')
        k=0
        do K_i=1,SIDEz
         do k_j=1,SIDEy
          do K_k=1,SIDEx
          read(23,*) x,y,z,M20(k_k,k_j,k_i),M21(k_k,k_j,k_i),M22(k_k,k_j,k_i) 
          k=k+1
          enddo
         enddo
        enddo
  close(23)
endif

      call mp_bcast( WA_global , master )                                    
      call mp_bcast( WB_global , master )         

      call mp_bcast(M00,master)                           
      call mp_bcast(M01,master)                           
      call mp_bcast(M02,master)                           
      call mp_bcast(M10,master)                           
      call mp_bcast(M11,master)                           
      call mp_bcast(M12,master)                           
      call mp_bcast(M20,master)                           
      call mp_bcast(M21,master)                           
      call mp_bcast(M22,master)                           
 
# if defined (SLAB)
# else /* SLAB */
         k=0
        do K_i=xstart(3),xend(3)
          do k_j=xstart(2),xend(2)
            do K_k=xstart(1),xend(1)
             wA(k)=wA_global(k_k,k_j,k_i)
             wB(k)=wA_global(k_k,k_j,k_i)
             M_OP(k,0,0)=M00(k_k,k_j,k_i)
             M_OP(k,0,1)=M01(k_k,k_j,k_i)
             M_OP(k,0,2)=M02(k_k,k_j,k_i)
             M_OP(k,1,0)=M10(k_k,k_j,k_i)
             M_OP(k,1,1)=M11(k_k,k_j,k_i)
             M_OP(k,1,2)=M12(k_k,k_j,k_i)
             M_OP(k,2,0)=M20(k_k,k_j,k_i)
             M_OP(k,2,1)=M21(k_k,k_j,k_i)
             M_OP(k,2,2)=M22(k_k,k_j,k_i)
             k=k+1
            enddo
          enddo
        enddo
# endif /* SLAB */

     deallocate(WA_global)
     deallocate(WB_global)
      deallocate(M00)
      deallocate(M01)
      deallocate(M02)
      deallocate(M10)
      deallocate(M11)
      deallocate(M12)
      deallocate(M20)
      deallocate(M21)
      deallocate(M22)

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
# if defined (Dim2)
# else /* Dim2 */
 deallocate(matrix_Rz)
# endif /* Dim2 */
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

# if defined (Dim2)
subroutine init_k_index()
 USE nrtype,only :DP,PI
 USE global_para
 USE mpi
 USE control
 USE constants
 USE utility
 USE mmpi
 USE mpi_fftw3_operation
 implicit none
 integer :: istat,k,k_i,k_j
 allocate(kD(0:LOCAL_SIZE-1),stat=istat)
 if(istat/=0) then
  if(myid==0) then
  write(*,*) "error alocate k_index"
  call mp_barrier()
  stop
  endif
 endif

 k=0
 do k_i=0,local_nx-1
    do k_j=0,SIDEy-1
        kD(k)%x=k_j+1
        kD(k)%y=k_i+1
        k=k+1
    enddo
 enddo
end subroutine init_k_index
# else  /* Dim2 */
subroutine init_k_index()
 USE nrtype,only :DP,PI
 USE global_para
 USE mpi
 USE control
 USE constants
 USE utility
 USE mmpi
 USE mpi_fftw3_operation
 implicit none
 integer :: istat,k,k_i,k_j,k_k
 allocate(kD(0:LOCAL_SIZE-1),stat=istat)
 if(istat/=0) then
  if(myid==0) then
  write(*,*) "error alocate k_index"
  call mp_barrier()
  stop
  endif
 endif
!!note that local_data is formed as (1:SIDEz,1:SIDEy,1:local_nx) not
!as (0:SIDEz-1,0:SIDEy-1,0:local_nx-1)
 k=0
 do k_i=0,local_nx-1
    do k_j=0,SIDEy-1
      do k_k=0,SIDEz-1
        kD(k)%x=k_k+1
        kD(k)%y=k_j+1
        kD(k)%z=k_i+1
        k=k+1
      enddo
    enddo
 enddo
end subroutine init_k_index
# endif  /* Dim2 */

# endif /* SLAB */

# if defined (SLAB)                       
subroutine K_index()
 USE global_para
 USE mpi
 USE control
 USE constants
 USE utility
 USE mmpi
 USE mpi_fftw3_operation
 implicit none
 integer :: i,j,k
 integer :: K_i,K_ii
 integer :: K_j,K_jj
 integer :: K_k,K_kk
integer :: istat2
!do K_i=0,local_nx 
   
allocate(k_index_1D(0:local_size-1))
k=0
do K_i=0,local_nx-1 
   
  if((k_i+local_x_start)<=(SiDEx/2)) then
   k_ii=k_i+local_x_start  
  else 
   k_ii=(k_i+local_x_start)-SIDEx 
  endif

    do k_j=0,SIDEy-1
    ! do k_j=1,1
        if(K_j<=(SIDEy/2)) then 
          k_jj=K_j
        else
          K_jj=K_j-SIDEy
        endif
# if defined (Dim2)
# else  /* Dim2 */
           do k_k=0,SIDEz-1
           !do k_k=0,0
             if(K_k<=(SIDEz/2)) then
             k_kk=K_k
             else
             k_kk=K_k-SIDEz
             endif
# endif /* Dim2 */
             
            k_index_1D(k)%i=k_ii
            k_index_1D(k)%j=k_jj
# if defined (Dim2)
# else  /* Dim2 */
            k_index_1D(k)%k=k_kk
# endif /* Dim2 */
            k=k+1
            enddo
      enddo
# if defined (Dim2)
# else  /* Dim2 */
 enddo
# endif /* Dim2 */

# if defined (Dim2)
     write(*,*) "done calling k_index 2D ","on",myid
# else  /* Dim2 */
     write(*,*) "done calling k_index 3D ","on",myid
# endif /* Dim2 */
     end subroutine K_index
# endif /* SLAB */
