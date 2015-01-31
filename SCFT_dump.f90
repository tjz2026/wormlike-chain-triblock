       subroutine SCFT_dump(converge)
       USE nrtype,only :DP
       USE global_para
       USE mpi
       USE control
       USE constants
       USE utility
       USE mmpi
# if defined (SLAB)
       USE mpi_fftw3_operation
# else /* SLAB */
       use decomp_fft_mod
# endif /* SLAB */
       !USE fftw !user defined module
       implicit none
       logical,intent(in) :: converge
       integer :: k,aaa,K_i,K_j,K_k
       double precision, allocatable :: RA_global(:,:,:) 
       double precision, allocatable :: RB_global(:,:,:) 

       character(len=30)::aa

       
       



        aaa=myid
       write(aa,*) aaa

      if(converge) then
       
         if(myid==0) then

         open(unit=41,file='error.dat',status='old',position='append')
         write(41,*) n_iter,ta_diff,tb_diff,tm_diff
         close(41)

         open(unit=44,file='FreeEG.dat',status='old',position='append')
         write(44,*) n_iter,FE_global,pff_global,tempEE_global,FE_global-NXAB*fa*fb
         write(44,*) n_iter,dx*SIDEx,"FE-FH",FE_global-NXAB*fa*fb
         close(44)

        open(unit=45,file='result.txt',status='old',position='append')  
        write(45,*) "FREE_ENERGY",FE_global,pff_global,tempEE_global,FE_global-NXAB*fa*fb
        write(45,*) "FREE_ENERGY-FH",FE_global-NXAB*fa*fb
        write(45,*) "box size",SIDEx*dx,SIDEy*dy
         write(45,*) "converge",converge 
        close(45)

        open(unit=45,file='result.dat',status='replace')  
        write(45,*) "FREE_ENERGY",FE_global,pff_global,tempEE_global,FE_global-NXAB*fa*fb
        write(45,*) "FREE_ENERGY-FH",FE_global-NXAB*fa*fb
        write(45,*) "box size",SIDEx*dx,SIDEy*dy
        write(45,*) "converge",converge 
        close(45)
         endif
# if defined (SLAB)
# if defined (Dim2)
  open(unit=23,file= 'RHO_total' //trim(adjustl(aa)) // '.dat',status='replace')
  ! write(23,*) "local_x_start=",local_x_start
         k=0
        do K_i=0,local_nx-1
         do k_j=0,SIDEy-1
            
          write(23,'(2f10.6,3f12.8)') (K_i+local_x_start)*dx,K_j*dy,RA(k),RB(k),1.0-RA(K)-RB(k)
          k=k+1
         enddo
        enddo
close(23)


  open(unit=23,file= 'RHO_half_end' //trim(adjustl(aa)) // '.dat',status='replace')
  ! write(23,*) "local_x_start=",local_x_start
         k=0
        do K_i=0,local_nx-1
        do k_j=0,SIDEy-1
          write(23,'(2f10.6,2f12.8)') (K_i+local_x_start)*dx,K_j*dy,R_half(k),R_end(k)
          k=k+1
         enddo
        enddo
close(23)
# else  /* Dim2 */

  open(unit=23,file= 'RHO_total' //trim(adjustl(aa)) // '.dat',status='replace')
  ! write(23,*) "local_x_start=",local_x_start
         k=0
        do K_i=0,local_nx-1
        do k_j=0,SIDEy-1
         do K_k=0,SIDEz-1
     write(23,'(3f10.6,3f12.8)') (K_i+local_x_start)*dx,K_j*dy,K_k*dz,RA(k),RB(k),1.0-RA(K)-RB(k)
          k=k+1
         enddo
         enddo
        enddo
close(23)

  open(unit=23,file= 'RHO_half_end' //trim(adjustl(aa)) // '.dat',status='replace')
  ! write(23,*) "local_x_start=",local_x_start
         k=0
        do K_i=0,local_nx-1
        do k_j=0,SIDEy-1
         do K_k=0,SIDEz-1
          write(23,'(3f10.6,2f12.8)') (K_i+local_x_start)*dx,K_j*dy,K_k*dz,R_half(k),R_end(k)
          k=k+1
         enddo
         enddo
        enddo
close(23)
# endif /* Dim2 */
# else /* SLAB */

      allocate(RA_global(1:SIDEx,1:SIDEy,1:SIDEz))
      allocate(RB_global(1:SIDEx,1:SIDEy,1:SIDEz))
      RA_global=0.0
      RB_global=0.0
!      call assembly_local_1d_to_global_3d(1,RA,RA_global,SIDEx,SIDEy,SIDEz)
!      call assembly_local_1d_to_global_3d(1,RB,RB_global,SIDEx,SIDEy,SIDEz)
if(myid==0) then
  open(unit=23,file= 'RHO_total.dat',status='replace')
         k=0
        do K_i=1,SIDEz
         do k_j=1,SIDEy
          do K_k=1,SIDEx
     write(23,'(3f10.6,3f12.8)') (K_k-1)*dx,(K_j-1)*dy,(K_i-1)*dz,RA_global(k_k,k_j,k_i), &
                 RB_global(k_k,k_j,k_i),1.0-RA_global(k_k,k_j,k_i)-RB_global(k_k,k_j,k_i) 
          k=k+1
          enddo
         enddo
        enddo
  close(23)
endif

      RA_global=0.0
      RB_global=0.0
!      call assembly_local_1d_to_global_3d(1,R_half,RA_global,SIDEx,SIDEy,SIDEz)
!      call assembly_local_1d_to_global_3d(1,R_end,RB_global,SIDEx,SIDEy,SIDEz)

if(myid==0) then
  open(unit=23,file= 'RHO_half_end.dat',status='replace')
         k=0
        do K_i=1,SIDEz
         do k_j=1,SIDEy
          do K_k=1,SIDEx
     write(23,'(3f10.6,2f12.8)') (K_k-1)*dx,(K_j-1)*dy,(K_i-1)*dz,RA_global(k_k,k_j,k_i), &
                                 RB_global(k_k,k_j,k_i)
          k=k+1
          enddo
         enddo
        enddo
close(23)
endif

deallocate(RA_global)
deallocate(RB_global)

# endif /* SLAB */
  open(unit=23,file= 'S_OP' //trim(adjustl(aa)) // '.dat',status='replace')
  do K=0,LOCAL_SIZE-1
  write(23,*) SIDEx*dx,SIDEy*dy,SIDEz*dz,S_OP(K,0,0),S_OP(K,0,1),S_OP(K,0,2)
  write(23,*) SIDEx*dx,SIDEy*dy,SIDEz*dz,S_OP(K,1,0),S_OP(K,1,1),S_OP(K,1,2)
  write(23,*) SIDEx*dx,SIDEy*dy,SIDEz*dz,S_OP(K,2,0),S_OP(K,2,1),S_OP(K,2,2)
  enddo
  close(23)

  open(unit=23,file= 'SA_OP' //trim(adjustl(aa)) // '.dat',status='replace')
  do K=0,LOCAL_SIZE-1
  write(23,*) SIDEx*dx,SIDEy*dy,SIDEz*dz,SA_OP(K,0,0),SA_OP(K,0,1),SA_OP(K,0,2)
  write(23,*) SIDEx*dx,SIDEy*dy,SIDEz*dz,SA_OP(K,1,0),SA_OP(K,1,1),SA_OP(K,1,2)
  write(23,*) SIDEx*dx,SIDEy*dy,SIDEz*dz,SA_OP(K,2,0),SA_OP(K,2,1),SA_OP(K,2,2)
  enddo
  close(23)

  open(unit=23,file= 'SB_OP' //trim(adjustl(aa)) // '.dat',status='replace')
  do K=0,LOCAL_SIZE-1
  write(23,*) SIDEx*dx,SIDEy*dy,SIDEz*dz,SB_OP(K,0,0),SB_OP(K,0,1),SB_OP(K,0,2)
  write(23,*) SIDEx*dx,SIDEy*dy,SIDEz*dz,SB_OP(K,1,0),SB_OP(K,1,1),SB_OP(K,1,2)
  write(23,*) SIDEx*dx,SIDEy*dy,SIDEz*dz,SB_OP(K,2,0),SB_OP(K,2,1),SB_OP(K,2,2)
  enddo
  close(23)

     else
        if(myid==0) then
        
         open(unit=41,file='error.dat',status='old',position='append')
         write(41,*) n_iter,ta_diff,tb_diff,tm_diff
         close(41)
        endif

         if(myid==0) then
         open(unit=44,file='FreeEG.dat',status='old',position='append')
         write(44,*) n_iter,FE_global,pff_global,tempEE_global
         write(44,*) n_iter,dx*SIDEx,"FE-FH",FE_global-NXAB*fa*fb
         write(44,*) "dx,dy",dx,dy
        write(44,*) "box size",SIDEx*dx,SIDEy*dy
         close(44)
         endif
# if defined (SLAB)
# if defined (Dim2)
  open(unit=23,file= 'RHO_total' //trim(adjustl(aa)) // '.dat',status='replace')
  ! write(23,*) "local_x_start=",local_x_start
         k=0
        do K_i=0,local_nx-1
        do k_j=0,SIDEy-1
            
          write(23,'(2f10.6,3f12.8)') (K_i+local_x_start)*dx,K_j*dy,RA(k),RB(k),1.0-RA(K)-RB(k)
          k=k+1
         enddo
        enddo
close(23)
# else /* Dim2 */
  open(unit=23,file= 'RHO_total' //trim(adjustl(aa)) // '.dat',status='replace')
  ! write(23,*) "local_x_start=",local_x_start
         k=0
        do K_i=0,local_nx-1
        do k_j=0,SIDEy-1
         do K_k=0,SIDEz-1
    write(23,'(3f10.6,3f12.8)') (K_i+local_x_start)*dx,K_j*dy,K_k*dz,RA(k),RB(k),1.0-RA(K)-RB(k)
          k=k+1
         enddo
         enddo
        enddo
close(23)
# endif  /* Dim2 */
# else /* SLAB */
      allocate(RA_global(1:SIDEx,1:SIDEy,1:SIDEz))
      allocate(RB_global(1:SIDEx,1:SIDEy,1:SIDEz))
      RA_global=0.0
      RB_global=0.0
!      call assembly_local_1d_to_global_3d(1,RA,RA_global,SIDEx,SIDEy,SIDEz)
!      call assembly_local_1d_to_global_3d(1,RB,RB_global,SIDEx,SIDEy,SIDEz)
if(myid==0) then
  open(unit=23,file= 'RHO_total.dat',status='replace')
         k=0
        do K_i=1,SIDEz
         do k_j=1,SIDEy
          do K_k=1,SIDEx
     write(23,'(3f10.6,3f12.8)') (K_k-1)*dx,(K_j-1)*dy,(K_i-1)*dz,RA_global(k_k,k_j,k_i), &
                 RB_global(k_k,k_j,k_i),1.0-RA_global(k_k,k_j,k_i)-RB_global(k_k,k_j,k_i) 
          k=k+1
          enddo
         enddo
        enddo
  close(23)
endif
      RA_global=0.0
      RB_global=0.0
!      call assembly_local_1d_to_global_3d(1,R_half,RA_global,SIDEx,SIDEy,SIDEz)
!      call assembly_local_1d_to_global_3d(1,R_end,RB_global,SIDEx,SIDEy,SIDEz)

if(myid==0) then
  open(unit=23,file= 'RHO_half_end.dat',status='replace')
         k=0
        do K_i=1,SIDEz
         do k_j=1,SIDEy
          do K_k=1,SIDEx
     write(23,'(3f10.6,2f12.8)') (K_k-1)*dx,(K_j-1)*dy,(K_i-1)*dz,RA_global(k_k,k_j,k_i), &
                                 RB_global(k_k,k_j,k_i)
          k=k+1
          enddo
         enddo
        enddo
close(23)
endif

 deallocate(RA_global)
 deallocate(RB_global)

# endif /* SLAB */
endif !if converge
end subroutine SCFT_dump

  subroutine W_dump()
   USE nrtype,only :DP
 USE global_para
 USE mpi
 USE control
 USE constants
 USE utility
 USE mmpi
# if defined (SLAB)
       USE mpi_fftw3_operation
# else /* SLAB */
       use decomp_fft_mod
# endif /* SLAB */
  implicit none
  integer :: i,j,k,aaa,k_k,k_j,k_i
  REAL*8 :: x,y,z
  character(len=30)::aa
  double precision, allocatable :: WA_global(:,:,:) 
  double precision, allocatable :: WB_global(:,:,:) 

  aaa=myid
  write(aa,*) aaa
  open(unit=23,file= 'W' //trim(adjustl(aa)) // '.dat',status='replace')
  do K=0,LOCAL_SIZE-1
  write(23,*) K,WA(K),WB(K)
  enddo
  close(23)

  
  open(unit=23,file= 'M_OP' //trim(adjustl(aa)) // '.dat',status='replace')
  do K=0,LOCAL_SIZE-1
  write(23,*) SIDEx*dx,SIDEy*dy,SIDEz*dz,M_OP(K,0,0),M_OP(K,0,1),M_OP(K,0,2)
  write(23,*) SIDEx*dx,SIDEy*dy,SIDEz*dz,M_OP(K,1,0),M_OP(K,1,1),M_OP(K,1,2)
  write(23,*) SIDEx*dx,SIDEy*dy,SIDEz*dz,M_OP(K,2,0),M_OP(K,2,1),M_OP(K,2,2)
  enddo
  close(23)


      allocate(WA_global(1:SIDEx,1:SIDEy,1:SIDEz))
      allocate(WB_global(1:SIDEx,1:SIDEy,1:SIDEz))

      !allocate(M00_local(0:LOCAL_SIZE-1))

      WB_global=0.0
      !call assembly_local_1d_to_global_3d(1,WA,WA_global,SIDEx,SIDEy,SIDEz)
      !call assembly_local_1d_to_global_3d(1,WB,WB_global,SIDEx,SIDEy,SIDEz)

if(myid==0) then
  open(unit=23,file= 'W_global.dat',status='replace')
         k=0
        do K_i=1,SIDEz
         do k_j=1,SIDEy
          do K_k=1,SIDEx
     write(23,'(3f10.6,2f12.8)') (K_k-1)*dx,(K_j-1)*dy,(K_i-1)*dz,WA_global(k_k,k_j,k_i), &
                                  WB_global(k_k,k_j,k_i)
          k=k+1
          enddo
         enddo
        enddo
  close(23)
 endif
     deallocate(WA_global)
     deallocate(WB_global)

 end subroutine W_dump




