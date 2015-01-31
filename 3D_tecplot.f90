!!!! this program is used for gathering the density profile from each processor
!and write the profile in the format of tecplot .

         program main
         implicit none
         integer :: ii,K_i,K_j,K_k,K,aaa
         real*8 :: dx,dy,dz
         real*8,allocatable :: RHO_A(:),RHO_B(:),RHO_C(:)
         real*8,allocatable :: R_X(:),R_Y(:),R_Z(:)
         integer,parameter :: SIDEx=32,SIDEy=32,SIDEz=32,local_nx=2,filenumber=16
          character(len=30)::aa

         
        allocate(RHO_A(0:local_nx*SIDEy*SIDEz-1))
        allocate(RHO_B(0:local_nx*SIDEy*SIDEz-1))
        allocate(RHO_C(0:local_nx*SIDEy*SIDEz-1))
        allocate(R_X(0:local_nx-1))
        allocate(R_Y(0:SIDEy-1))
        allocate(R_Z(0:SIDEz-1))





         open(unit=45,file='dfpde_total.dat',status='new')
          WRITE (45,2501) 
          WRITE (45,2502) SIDEz,SIDEy,SIDEx
         close(45)

         open(unit=45,file='Pu_total.dat',status='new')
          WRITE (45,2501) 
          WRITE (45,2502) SIDEz,SIDEy,SIDEx
         close(45)

       do ii=0,filenumber-1
          aaa=ii
          write(aa,*) aaa
          open(unit=33,file='RHO_total' //trim(adjustl(aa)) // '.dat',status='old') 
          k=0
          do k_i=0,Local_nx-1
             do k_j=0,SIDEy-1
              do K_k=0,SIDEz-1
               read(33,*) R_X(K_i),R_Y(K_j),R_Z(K_k),RHO_A(k),RHO_B(k),RHO_C(k)
                  k=k+1
               enddo
              enddo
           enddo
          close(33)

     open(unit=45,file='dfpde_total.dat',status='old',position='append')
         
          k=0
          do k_i=0,Local_nx-1
             do k_j=0,SIDEy-1
              do K_k=0,SIDEz-1
               write(45,2504) R_X(K_i),R_Y(K_j),R_Z(K_k),RHO_A(k),RHO_B(k),1.0-RHO_C(k)
                  k=k+1
              enddo
             enddo
          enddo
      close(45)
       enddo

       do ii=0,filenumber-1
          aaa=ii
          write(aa,*) aaa
          open(unit=33,file='Pu' //trim(adjustl(aa)) // '.dat',status='old') 
          k=0
          do k_i=0,Local_nx-1
             do k_j=0,SIDEy-1
              do K_k=0,SIDEz-1
               read(33,*) R_X(K_i),R_Y(K_j),R_Z(K_k),RHO_A(k),RHO_B(k),RHO_C(k)
                  k=k+1
               enddo
              enddo
           enddo
          close(33)

     open(unit=45,file='Pu_total.dat',status='old',position='append')
         
          k=0
          do k_i=0,Local_nx-1
             do k_j=0,SIDEy-1
              do K_k=0,SIDEz-1
               write(45,2504) R_X(K_i),R_Y(K_j),R_Z(K_k),RHO_A(k),RHO_B(k),RHO_C(k)
                  k=k+1
              enddo
             enddo
          enddo
      close(45)







         enddo        
           
       

         







2501 FORMAT ('VARIABLES =',1X,'"X"',1X,'"Y"',1X,'"Z"',1X,'"A"',1X,'"B"',1X,' "C"')
2502 FORMAT ('ZONE',1x,'I=',I3,1X,'J=',I3,1X,'K=',I3,1X,'F=POINT')
2503 FORMAT (I4,2X,I4,2X,I4,2X,F9.5,2X,F9.5)
2504 FORMAT (F9.5,2X,F9.5,2X,F9.5,2X,F9.5,2X,F9.5,2X,F9.5)

deallocate(RHO_A)
deallocate(RHO_B)
deallocate(RHO_C)
deallocate(R_X)
deallocate(R_Y)
deallocate(R_Z)

      end 



