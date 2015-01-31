# if defined (SLAB)
# else /* SLAB */

# if defined (SPH) 
subroutine init_arrays_2decomp_SPH()
 USE nrtype,only :DP
 USE global_para
 USE mpi
 USE control
 USE constants
 USE utility
 USE mmpi
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

 allocate(M_OP(0:LS,0:N_dim_ddm-1,0:N_dim_ddm-1),stat=istat2)
 allocate(S_OP(0:LS,0:N_dim_ddm-1,0:N_dim_ddm-1),stat=istat2)
 allocate(SA_OP(0:LS,0:N_dim_ddm-1,0:N_dim_ddm-1),stat=istat2)
 allocate(SB_OP(0:LS,0:N_dim_ddm-1,0:N_dim_ddm-1),stat=istat2)


 allocate(Pu(0:LS,1:N_dim_ddm),stat=istat2)
 allocate(Pu_A(0:LS,1:N_dim_ddm),stat=istat2)
 allocate(Pu_B(0:LS,1:N_dim_ddm),stat=istat2)

call mp_barrier()
if(istat2/=0) then
write(*,*) "allocate failed,exit!","on",myid
stop
else
write(*,*) "allocate arrays succeeded","on",myid

endif 

call mp_barrier()

end subroutine init_arrays_2decomp_SPH


          
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
 deallocate(wA)
 deallocate(wB)
 deallocate(RA)
 deallocate(RB)
 deallocate(R_half)
 deallocate(R_end)

 deallocate(M_OP)
 deallocate(S_OP)
 deallocate(SA_OP)
 deallocate(SB_OP)

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
