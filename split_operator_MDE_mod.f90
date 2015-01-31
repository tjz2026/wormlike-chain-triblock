# if defined (SPH)
module split_operator_MDE
implicit none
 
      double precision, allocatable :: q_f(:,:,:,:)
      double precision, allocatable :: q_b(:,:,:,:)
      double precision, dimension(:), allocatable :: c_nor,c_nor_star
      real*8,allocatable,save :: exp_L_A(:),exp_L_B(:)
      real*8,allocatable,save :: exp_u_dot_k_r(:,:,:)  
      real*8,allocatable,save :: exp_u_dot_k_i(:,:,:)  

       public  split_operator_MDE_init  
       public  AB_diblock_sph_driver
       public  c_nor
       public  c_nor_star
       public  q_f
       public  q_b



contains
 
        ! needs to be called every time the cell size ( dx,dy,dz) is changed
        subroutine split_operator_MDE_init()
        use global_para
        use spherepack
        use grid_mod
        use mmpi
        implicit none
        integer :: LS,LOCAL_SIZEK,istat,i,j,k,nlat
        double precision :: temp
        REAL(DP) :: a(Ntheta+1,Ntheta+1),b(Ntheta+1,Ntheta+1)
        REAL(DP) :: r(Ntheta+1,Nphi)
  
        LS=LOCAL_SIZE-1
# if defined (R2C)
LOCAL_SIZEK=LOCAL_SIZE_K_R2C
# else /* R2C */
LOCAL_SIZEK=LOCAL_SIZE_K
# endif /* R2C */
        nlat=Ntheta+1
       if (allocated(q_f)) deallocate(q_f)
       allocate(q_f(0:LS,1:Ntheta+1,1:Nphi,0:Nmax),stat=istat)
       if (allocated(q_b)) deallocate(q_b)
       allocate(q_b(0:LS,1:Ntheta+1,1:Nphi,0:Nmax),stat=istat)
       if(istat/=0) then
       write(*,*) "unsucessful in allocate qf and qb"
       stop
       else 
       write(*,*) "sucessful in allocate qf and qb"
       endif  

         
       if (allocated(c_nor)) deallocate(c_nor)
       allocate(c_nor(0:Nmax))
       if (allocated(c_nor_star)) deallocate(c_nor_star)
       allocate(c_nor_star(0:Nmax))
         c_nor=0.0d0
         c_nor_star=0.0d0

   if(allocated(exp_u_dot_k_i)) deallocate(exp_u_dot_k_i)
   allocate(exp_u_dot_k_i(0:LOCAL_SIZEK-1,1:Ntheta+1,1:Nphi),stat=istat)
   if(allocated(exp_u_dot_k_r)) deallocate(exp_u_dot_k_r)
   allocate(exp_u_dot_k_r(0:LOCAL_SIZEK-1,1:Ntheta+1,1:Nphi),stat=istat)

   if(allocated(exp_L_A)) deallocate(exp_L_A)
   allocate(exp_L_A(1:nlat),stat=istat)
   if(allocated(exp_L_B)) deallocate(exp_L_B)
   allocate(exp_L_B(1:nlat),stat=istat)


    do i=1,Nphi
      do j=1,Ntheta+1
        do k=0,LOCAL_SIZEK-1
# if defined (R2C)
         temp=unit_vec(j,i)%x*k_vec_r2c(k)%x*dx_inv + &
              unit_vec(j,i)%y*k_vec_r2c(k)%y*dy_inv + &
              unit_vec(j,i)%z*k_vec_r2c(k)%z*dz_inv
# else /* R2C */
         temp=unit_vec(j,i)%x*k_vec(k)%x*dx_inv + &
              unit_vec(j,i)%y*k_vec(k)%y*dy_inv + &
              unit_vec(j,i)%z*k_vec(k)%z*dz_inv
# endif /* R2C */
        temp=-0.5d0*ds*temp
        exp_u_dot_k_r(k,j,i)=cos(temp)
        exp_u_dot_k_i(k,j,i)=sin(temp)
        enddo 
      enddo
    enddo

       

       do i=1,nlat
          temp=ds*Loa*(i-1)*i
          exp_L_A(i)= dexp(-temp/(2.0d0*kapa_A)) 
          exp_L_B(i)= dexp(-temp/(2.0d0*kapa_B)) 
       enddo   

       call sphscalar_init(r,a,b)


        end subroutine split_operator_MDE_init

        subroutine split_operator_MDE_clean()
        use global_para
        implicit none
       if (allocated(q_f)) deallocate(q_f)
       if (allocated(q_b)) deallocate(q_b)
       if (allocated(exp_u_dot_k_r)) deallocate(exp_u_dot_k_r)
       if (allocated(exp_u_dot_k_i)) deallocate(exp_u_dot_k_i)
       if (allocated(exp_L_A)) deallocate(exp_L_A)
       if (allocated(exp_L_B)) deallocate(exp_L_B)
       if (allocated(c_nor)) deallocate(c_nor)
       if (allocated(c_nor_star)) deallocate(c_nor_star)
          

        end subroutine split_operator_MDE_clean
       

      subroutine AB_diblock_sph_driver()
       USE nrtype,only :DP,PI
       USE global_para
       USE mpi
       USE mmpi
       USE control
       USE constants
       USE utility
       USE spherepack
       use decomp_fft_mod,only:xsize
       implicit none
       integer :: k,i,j,s,nlat,istat2,n1,n2,n3
       real(DP),allocatable :: q_temp1(:,:,:)
       real(DP) :: q_temp2(0:Local_size-1)
        !for A
       real(DP),allocatable :: exp_WA(:)
        !for B
       real(DP),allocatable :: exp_WB(:)
       REAL(DP) :: M_grid_inv,temp,c_nor_inv,c_nor_local,c_nor_global
       REAL(DP) :: a(Ntheta+1,Ntheta+1),b(Ntheta+1,Ntheta+1)
       REAL(DP) :: as(Ntheta+1,Ntheta+1),bs(Ntheta+1,Ntheta+1)
       REAL(DP) :: r_u(Ntheta+1,Nphi)
      ! remember nz is refered to as the local grid number always the first dimension
      ! beleow is just set temporarily
      n1=xsize(1)
      n2=xsize(2)
      n3=xsize(3)

       M_grid_inv=1.0d0/M_grid 
       nlat=Ntheta+1 

!for A
 allocate(exp_WA(0:LOCAL_SIZE-1),stat=istat2)
!for B
 allocate(exp_WB(0:LOCAL_SIZE-1),stat=istat2)
 
 allocate(q_temp1(0:LOCAL_SIZE-1,1:Ntheta+1,1:Nphi),stat=istat2)

 if(istat2/=0) then
 write(*,*) "allocate matrix in mde_q failed"
 stop
 endif

        do k=0,LOCAL_SIZE-1
            exp_WA(k)=dexp(-0.5d0*ds*wA(k))
            exp_WB(k)=dexp(-0.5d0*ds*wB(k))
        enddo
             
        ! init the initial condition for qf 
    do i=1,Nphi
      do j=1,Ntheta+1
        do k=0,LOCAL_SIZE-1
           q_f(k,j,i,0)=1.0d0
           q_temp1(k,j,i)=q_f(k,j,i,0)
       enddo
      enddo
     enddo   
       call mp_barrier()
   !     c_nor(0)=simposon_3D_1D_mpi(xsize(1),xsize(2),xsize(3),dx,dy,dz,q_temp1(:,1,1))

      c_nor(0)= simposon_sph_spa3d_mpi(Ntheta,Nphi,n1,n2,n3,dtheta,dphi, &
                                       dx,dy,dz,q_temp1)

       c_nor_local=c_nor(0)
              c_nor_global=0.0d0
              call mp_barrier()
              call mp_allreduce(c_nor_local,c_nor_global)
              call mp_barrier()
      c_nor(0)=c_nor_global
      c_nor_inv=1.0d0/c_nor(0)
      !q_f(:,:,:,0)=q_f(:,:,:,0)*c_nor_inv

         
       do s=1,NA
       !call sph_onestep_MDE(q_f(:,:,:,s-1),q_f(:,:,:,s),s,exp_WA,exp_L_A,1)
       call sph_onestep_MDE(s,exp_WA,exp_L_A,1)
       enddo ! s=1,Nmax

       do s=NA+1,Nmax
       call sph_onestep_MDE(s,exp_WB,exp_L_B,1)
       enddo ! s=1,Nmax

        write(*,*) "done solver qB in one step split operator MDE"
       call mp_barrier()
      ! note that q_star(:,:,:,Nmax)=1.0,so c_nor_star(Nmax)=c_nor(0)
      ! init the initial condition for qf 
    do i=1,Nphi
      do j=1,Ntheta+1
        do k=0,LOCAL_SIZE-1
           q_b(k,j,i,Nmax)=1.0d0
           q_temp1(k,j,i)=q_b(k,j,i,Nmax)
       enddo
      enddo
     enddo   
                 
      c_nor_star(Nmax)=c_nor(0)
      c_nor_inv=1.0d0/c_nor_star(Nmax)

      !q_b(:,:,:,Nmax)=q_b(:,:,:,Nmax)*c_nor_inv

       
       do s=Nmax-1,NA+1,-1
       call sph_onestep_MDE(s,exp_WB,exp_L_B,-1)
       enddo 

       do s=NA,0,-1
       call sph_onestep_MDE(s,exp_WA,exp_L_A,-1)
       enddo 

        write(*,*) "done solver qBstar,qAstar in one step split operator MDE"
       call mp_barrier()

      deallocate(exp_WA)
      deallocate(exp_WB)
      deallocate(q_temp1)

        write(*,*) "done driver one step split operator MDE"
       call mp_barrier()
      end subroutine AB_diblock_sph_driver


       !subroutine sph_onestep_MDE(q_temp_in,q_temp_out,s,exp_W,exp_L_temp,orientation)
       subroutine sph_onestep_MDE(s,exp_W,exp_L_temp,orientation)
       USE nrtype,only :DP,PI
       USE global_para
       USE mpi
       USE mmpi
       USE control
       USE constants
       USE utility
       USE decomp_fft_mod
       USE spherepack
       implicit none
       integer :: k,i,j,nlat,LOCAL_SIZEK
       integer,intent(in) :: orientation
       !real(DP),intent(in) ::  q_temp_in(0:LOCAL_SIZE-1,1:Ntheta+1,1:Nphi)
       !real(DP),intent(out) :: q_temp_out(0:LOCAL_SIZE-1,1:Ntheta+1,1:Nphi)
       real(DP),allocatable :: q_temp_out(:,:,:)  
       real(DP) :: q_temp2(0:Local_size-1)
       integer,intent(in) ::  s
       real(DP),intent(in) :: exp_W(0:local_size-1)
       real(DP),intent(in) :: exp_L_temp(1:Ntheta+1)
       REAL(DP) :: M_grid_inv,temp,c_nor_inv,c_nor_global,c_nor_local
       REAL(DP) :: a(Ntheta+1,Ntheta+1),b(Ntheta+1,Ntheta+1)
       REAL(DP) :: as(Ntheta+1,Ntheta+1),bs(Ntheta+1,Ntheta+1)
       REAL(DP) :: r(Ntheta+1,Nphi)

       

     allocate(q_temp_out(0:local_size-1,1:Ntheta+1,1:Nphi))

     if(orientation==1) then
     q_temp_out(:,:,:)=q_f(:,:,:,s-1)
     else
     q_temp_out(:,:,:)=q_b(:,:,:,s+1)
     endif
     




# if defined (R2C)
LOCAL_SIZEK=LOCAL_SIZE_K_R2C
# else /* R2C */
LOCAL_SIZEK=LOCAL_SIZE_K
# endif /* R2C */

       M_grid_inv=1.0d0/M_grid 
       nlat=Ntheta+1 
 
        ! step 1 and 2, q multiply exp_wA and then a pair of fft 
           do i=1,Nphi
             do j=1,Ntheta+1

                do k=0,LOCAL_SIZE-1
# if defined (R2C)
                 !local_in_r(kD(k)%x,kD(k)%y,KD(k)%z)=q_temp_in(k,j,i)*exp_W(k)
                 local_in_r(kD(k)%x,kD(k)%y,KD(k)%z)=q_temp_out(k,j,i)*exp_W(k)
# else /* R2C */
                 !local_in(kD(k)%x,kD(k)%y,KD(k)%z)=dcmplx(q_temp_in(k,j,i)*exp_W(k),0.0d0)
                 local_in(kD(k)%x,kD(k)%y,KD(k)%z)=dcmplx(q_temp_out(k,j,i)*exp_W(k),0.0d0)
# endif /* R2C */
                enddo

# if defined (R2C)
        call decomp_2d_fft_3d(local_in_r, local_out_c)
# else /* R2C */
        call decomp_2d_fft_3d(local_in, local_out, DECOMP_2D_FFT_FORWARD)
# endif /* R2C */

                
                 do k=0,LOCAL_SIZEK-1
# if defined (R2C)
         local_out_c(kDk(k)%x,kDk(k)%y,kDk(k)%z)=dcmplx(exp_u_dot_k_r(k,j,i), & 
                                               orientation*exp_u_dot_k_i(k,j,i))* & 
                                            local_out_c(kDk(k)%x,kDk(k)%y,kDk(k)%z) 
# else /* R2C */
         local_out(kDz(k)%x,kDz(k)%y,kDz(k)%z)=dcmplx(exp_u_dot_k_r(k,j,i), & 
                                               orientation*exp_u_dot_k_i(k,j,i))* & 
                                            local_out(kDz(k)%x,kDz(k)%y,kDz(k)%z) 
# endif /* R2C */
                 enddo
                

# if defined (R2C)
        call decomp_2d_fft_3d(local_out_c, local_in_r)
# else /* R2C */
        call decomp_2d_fft_3d(local_out, local_in, DECOMP_2D_FFT_BACKWARD)
# endif /* R2C */

                do k=0,LOCAL_SIZE-1
# if defined (R2C)
                 q_temp_out(k,j,i)=local_in_r(kD(k)%x,kD(k)%y,kD(k)%z)*M_grid_inv
# else /* R2C */
                 q_temp_out(k,j,i)=real(local_in(kD(k)%x,kD(k)%y,kD(k)%z))*M_grid_inv
# endif /* R2C */
                enddo 

              enddo !Ntheta+1
            enddo !Nphi
               

            do k=0,LOCAL_SIZE-1

                do i=1,Nphi
                  do j=1,Ntheta+1
                    r(j,i)=q_temp_out(k,j,i) 
                  enddo
                 enddo
 
                !step 3 doing  sphere transformation 
                call sph_scalar_analysis(r,a,b)

                      as=0.0d0
                      bs=0.0d0
                      !r=0.0d0
                    do i=1,Ntheta+1
                      do j=i,Ntheta+1
                        as(i,j)=exp_L_temp(j)*a(i,j)  
                        bs(i,j)=exp_L_temp(j)*b(i,j)
                      enddo
                     enddo 
                call sph_scalar_synthesis(r,as,bs)
                 
                do i=1,Nphi
                  do j=1,Ntheta+1
                    q_temp_out(k,j,i)=r(j,i) 
                  enddo
                 enddo

              enddo ! enddo local_size-1

              ! step 4,5: doing the other half part like step 1 and step 2

            do i=1,Nphi
             do j=1,Ntheta+1

                do k=0,LOCAL_SIZE-1
# if defined (R2C)
                 local_in_r(kD(k)%x,kD(k)%y,KD(k)%z)=q_temp_out(k,j,i)
# else /* R2C */
                 local_in(kD(k)%x,kD(k)%y,KD(k)%z)=dcmplx(q_temp_out(k,j,i),0.0d0)
# endif /* R2C */
                enddo

# if defined (R2C)
        call decomp_2d_fft_3d(local_in_r, local_out_c)
# else /* R2C */
        call decomp_2d_fft_3d(local_in, local_out, DECOMP_2D_FFT_FORWARD)
# endif /* R2C */


                 do k=0,LOCAL_SIZEK-1
# if defined (R2C)
         local_out_c(kDk(k)%x,kDk(k)%y,kDk(k)%z)=dcmplx(exp_u_dot_k_r(k,j,i), & 
                                               orientation*exp_u_dot_k_i(k,j,i))* & 
                                            local_out_c(kDk(k)%x,kDk(k)%y,kDk(k)%z) 
# else /* R2C */
         local_out(kDz(k)%x,kDz(k)%y,kDz(k)%z)=dcmplx(exp_u_dot_k_r(k,j,i), & 
                                               orientation*exp_u_dot_k_i(k,j,i))* & 
                                            local_out(kDz(k)%x,kDz(k)%y,kDz(k)%z) 
# endif /* R2C */
                 enddo


# if defined (R2C)
        call decomp_2d_fft_3d(local_out_c, local_in_r)
# else /* R2C */
        call decomp_2d_fft_3d(local_out, local_in, DECOMP_2D_FFT_BACKWARD)
# endif /* R2C */


                do k=0,LOCAL_SIZE-1
# if defined (R2C)
                 q_temp_out(k,j,i)=local_in_r(kD(k)%x,kD(k)%y,kD(k)%z)*M_grid_inv
# else /* R2C */
                 q_temp_out(k,j,i)=real(local_in(kD(k)%x,kD(k)%y,kD(k)%z))*M_grid_inv
# endif /* R2C */
                 q_temp_out(k,j,i)=q_temp_out(k,j,i)*exp_W(k)
                enddo

              enddo !Ntheta+1
            enddo !Nphi

           ! doing normalization of qf

!      if(orientation==1) then
!      c_nor(s)= simposon_sph_spa3d_mpi(Ntheta,Nphi,xsize(1),xsize(2), & 
!                                       xsize(3),dtheta,dphi, & 
!                                       dx,dy,dz,q_temp_out)
!      c_nor_local=c_nor(s)
!      else
!      c_nor_star(s)= simposon_sph_spa3d_mpi(Ntheta,Nphi,xsize(1),xsize(2), &
!                                           xsize(3),dtheta, & 
!                                            dphi,dx,dy,dz,q_temp_out)
!      c_nor_local=c_nor_star(s)
!      endif 
!
!              c_nor_global=0.0d0
!              call mp_allreduce(c_nor_local,c_nor_global)
!              call mp_barrier()
!      if(orientation==1) then
!              c_nor(s)=c_nor_global
!       else 
!              c_nor_star(s)=c_nor_global
!      endif
!
!      c_nor_inv=1.0d0/c_nor_global
!      !write(*,*) "bg norm q_temp_out(0,2,2)",q_temp_out(0,2,2)
!
!      q_temp_out(:,:,:)=q_temp_out(:,:,:)*c_nor_inv
!   !    if(s>91 .and. orientation==-1) then
!   !   write(*,*) "s=",s 
!   !   write(*,*) "c_nor_inv=",c_nor_inv
!   !   write(*,*) "af norm q_temp_out(0,2,2)",q_temp_out(0,2,2)
!   !   endif
     if(orientation==1) then
     q_f(:,:,:,s)=q_temp_out(:,:,:)
     else
     q_b(:,:,:,s)=q_temp_out(:,:,:)
     endif
 
     deallocate(q_temp_out)      
      end subroutine sph_onestep_MDE


end module split_operator_MDE
# endif /* SPH */
