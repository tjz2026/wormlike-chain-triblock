# if defined (SLAB)
       FUNCTION cal_scft(tz_scal)
       USE nrtype,only :DP
       USE global_para
       USE mpi
       USE control
       USE constants
       USE utility
       USE mmpi
       USE mpi_fftw3_operation
       USE G_matrix_mod
       USE wormlike_MDE
       USE wormlike_MDE_AB
       USE  Gmatrix_lu_mod       
       implicit none
       REAL(DP),intent(inout) :: tz_scal
       logical :: converge
       REAL(DP) :: cal_scft
       Real(dp) :: t11,t22
       if(abs(WAB_init_type)==4) then
       !4 and -4 for HEX ,-6 means input from file
# if defined (Dim2)
       dx=tz_scal
       dy=dx*SIDEx/(sqrt(3.0)*SIDEy)
       M_v=dx*dy*SIDEx*SIDEy
# else  /* Dim2 */
       dx=tz_scal
       dy=dx*SIDEx/(sqrt(3.0)*SIDEy)
       dz=tz_scal
       M_v=dz*dx*dy*SIDEx*SIDEy*SIDEz
# endif  /* Dim2 */
       else if(abs(WAB_init_type)==5) then
# if defined (Dim2)
       dx=tz_scal*0.1
       dy=tz_scal
       M_v=dx*dy*SIDEx*SIDEy
# else  /* Dim2 */
       dx=tz_scal*0.1
       dy=tz_scal*0.1
       dz=tz_scal
       M_v=dz*dx*dy*SIDEx*SIDEy*SIDEz
# endif  /* Dim2 */
       else 
         
       !cubic cell
# if defined (Dim2)
       dx=tz_scal
       dy=tz_scal
       M_v=dx*dy*SIDEx*SIDEy
# else  /* Dim2 */
       dx=tz_scal
       dy=tz_scal
       dz=tz_scal
       M_v=dz*dx*dy*SIDEx*SIDEy*SIDEz
# endif  /* Dim2 */
       endif

       dx_inv=1.0d0/dx
       dy_inv=1.0d0/dy
       dz_inv=1.0d0/dz

       if(myid==0) then 
       write(*,*) "M_v=",M_v
       write(*,*) "dx,dy,dz",dx,dy,dz
       write(*,*) "SIDEx,SIDEy,SIDEz",dx,dy,dz
       endif
       converge=.false.
       n_iter=1
       n_iter_M=0

# if defined (Test)
      call test_inverse_sparse_GAB()
      call mp_barrier()
      call mp_finalize()
      stop 
# endif /* Test */
       
!!Remember to init GAB and W every time the box size changes.!!!    
       if(wAB_init_type<0) then
       call readW()
       else
       call init_w()
       endif

      call init_inverse_sparse_GAB()

      do while (converge/=.true. .and. n_iter<=MAXITS)
         if(myid==0) then
         write(*,*) "the",n_iter,"th iteration","on",myid
         write(*,*) "the",n_iter_M,"th iteration","on",myid
         endif

         call calc_THETAij()
         call mp_barrier()
         !call MDE_q()
         !call MDE_qstar()
         call AB_melt_driver()  
         call density()
         if(range_or_not==1)then   
          call free_energy_finite()
         call finite_range_iterate(converge )
         !call free_energy()
         !call AB_iterate(converge)
         else
         call free_energy()
         call AB_iterate(converge)
         endif

         call SCFT_dump(converge)
         
         n_iter=n_iter+1
         if(n_iter_M<=Num_simple_mixing_M) then
         n_iter_M=n_iter_M+1
         else
            if(mod(n_iter,n_WABiter_per_M)==0) then
             n_iter_M=n_iter_M+1
            endif
         endif 

        if(n_iter<=20) then
         converge=.false.
        endif
          
      enddo

        cal_scft=FE_global 

        if(myid==0) then
        write(*,*) "exit the calscft","dx=",dx,"converge=",converge
        write(*,*) "FreeEG:",FE_global,"the",n_iter,"iteration"
        write(*,*) "ta_diff:",ta_diff,"tb_diff",tb_diff,"tm_diff",tm_diff
        write(*,*) "now is the ",n_SCFT,"th scft"
        endif

        n_SCFT=n_SCFT+1
        dx_save(n_SCFT)=dx 

        if(n_SCFT<=20) then
        wA_save(:,n_SCFT)=wA(:)
        wB_save(:,n_SCFT)=wB(:)
        endif
        
       call SCFT_dump(converge)
       call w_dump()
      ! call mem_cal_scft_destroy()

      end FUNCTION cal_scft

# else /* SLAB */
       FUNCTION cal_scft(tz_scal)
       USE nrtype,only :DP
       USE global_para
       USE mpi
       USE control
       USE constants
       USE utility
       USE mmpi
       use decomp_fft_mod
       use Gmatrixfree_mod 
       USE G_matrix_init
       USE G_matrix_2dcomp_mod
       USE wormlike_MDE
       implicit none
       REAL(DP),intent(inout) :: tz_scal
       logical :: converge
       REAL(DP) :: cal_scft
       if(abs(WAB_init_type)==4) then
       !4 and -4 for HEX ,- means input from file
       dx=tz_scal
       dy=dx*SIDEx/(sqrt(3.0)*SIDEy)
       dz=tz_scal
       M_v=dz*dx*dy*SIDEx*SIDEy*SIDEz
       else
       !cubic cell
       dx=tz_scal
       dy=tz_scal
       dz=tz_scal
       M_v=dz*dx*dy*SIDEx*SIDEy*SIDEz
       endif

       dx_inv=1.0d0/dx
       dy_inv=1.0d0/dy
       dz_inv=1.0d0/dz


       if(myid==0) then 
       write(*,*) "M_v=",M_v
       write(*,*) "dx,dy,dz",dx,dy,dz
       write(*,*) "SIDEx,SIDEy,SIDEz",dx,dy,dz
       endif
       converge=.false.
       n_iter=1
       n_iter_M=0

# if defined (Test)
      call test_inverse_sparse_GAB()
      call mp_barrier()
      call mp_finalize()
      stop
# endif /* Test */
 
!       call test_inverse_sparse_GAB()
!            call mp_barrier()
       
!       if(myid==0) then
!        call dns2csr()
!        call testMKL() 
!       endif
!       call mp_barrier()
!       call init_sparse_LRij()
!       call G_matrix_generator()
!       call mp_barrier()
!       call mp_finalize()
!       stop
      
     call init_sparse_LRij()

      
!!Remember to init GAB and W every time the box size changes.!!!    
      call w_init_or_not()

      write(*,*) "done init w_init"
      call mp_barrier()
      !call init_inverse_sparse_GAB() ! Complex to Complex FFT
      !call init_inverse_sparse_GAB(.true.) ! Real to Complex FFT

# if defined (R2C)
       call init_inverse_sparse_GAB(.true.) ! Complex to Complex FFT
# else /* R2C */
       call init_inverse_sparse_GAB() ! Real to Complex FFT
# endif /* R2C */

      write(*,*) "done init GAB"
      call mp_barrier()


      do while (converge/=.true. .and. n_iter<=MAXITS)
         if(myid==0) then
         write(*,*) "the",n_iter,"th iteration","on",myid
         write(*,*) "the",n_iter_M,"th iteration","on",myid
         endif

         call calc_THETAij()
         write(*,*) "done calc_THETAij"
         call MDE_q(.true.)
         call MDE_qstar(.true.)

         call density()
      write(*,*) "done density"
      call mp_barrier()
         !stop

         call free_energy()
      write(*,*) "done free energy"
      call mp_barrier()
         call iterate(converge)
         call SCFT_dump(converge)
      write(*,*) "done SCFT dump"
      call mp_barrier()
         
         n_iter=n_iter+1
         if(n_iter_M<=Num_simple_mixing_M) then
         n_iter_M=n_iter_M+1
         else
            if(mod(n_iter,n_WABiter_per_M)==0) then
             n_iter_M=n_iter_M+1
            endif
         endif 

        if(n_iter<=20) then
         converge=.false.
        endif
      enddo
      
      ! call free_energy()
        cal_scft=FE_global 
        if(myid==0) then
        write(*,*) "exit the calscft","dx=",dx,"converge=",converge
        write(*,*) "FreeEG:",FE_global,"the",n_iter,"iteration"
        write(*,*) "ta_diff:",ta_diff,"tb_diff",tb_diff,"tm_diff",tm_diff
        write(*,*) "now is the ",n_SCFT,"th scft"
        endif
          
        n_SCFT=n_SCFT+1
        dx_save(n_SCFT)=dx 
        if(n_SCFT<=20) then
        wA_save(:,n_SCFT)=wA(:)
        wB_save(:,n_SCFT)=wB(:)
        endif
       call SCFT_dump(converge)
       call w_dump()

       call mem_cal_scft_destroy()

      end FUNCTION cal_scft
# endif /* SLAB */

