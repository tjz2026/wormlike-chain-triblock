
 subroutine init_pr_backward(s, dir,M_Bar,Pr)
 USE nrtype,only :DP,PI
 USE global_para
 USE mpi
 USE control
 USE constants
 USE utility
 USE mmpi
 USE mpi_fftw3_operation
 USE matrix_inverse
 USE G_matrix_mod
 implicit none
 integer :: i,j,s
 integer :: index_nonzero
 integer :: k
 integer :: istat2
 integer,intent(in) :: s,dir,M_Bar
 real(DP) :: tempA1,tempA2,tempA3,temp,tempx,Wtemp,sumaq,sumaq_global,sumaq1
 real(DP) :: sumPA,sumPA_global
 real(DP),intent(inout)    :: Pr(0:M_Bar-1,0:LOCAL_SIZE-1)
 integer :: s_block, s_BDF

   if(s==NMAX-1) then
          Pr=0.0
      do k=0,LOCAL_SIZE-1
          Wtemp=WB(k)
# if defined (MaierSaupe)
          do index_nonzero=THETA_nonzero_1D_B_indexK(k-1)+1, &
             THETA_nonzero_1D_B_indexK(k)
             i=THETA_nonzero_1D_B(index_nonzero)%i
             j=THETA_nonzero_1D_B(index_nonzero)%j
             tempx=THETA_nonzero_1D_B(index_nonzero)%value
             Pr(i,k)=Pr(i,k)+ds*tempx*qBstar(j,k,NB)
           enddo
# endif /* MaierSaupe */
            do i=0,M_Bar_B-1
            Pr(i,k)=Pr(i,k)+(1.0-ds*Wtemp)*qBstar(i,k,NB)
            enddo
        enddo !enddo k=0,LOCAL_SIZE  
   endif
   if(s==NMAX-2) then
          Pr=0.0
      do k=0,LOCAL_SIZE-1
          Wtemp=WB(k)
# if defined (MaierSaupe)
          do index_nonzero=THETA_nonzero_1D_B_indexK(k-1)+1, &
            THETA_nonzero_1D_B_indexK(k)
             i=THETA_nonzero_1D_B(index_nonzero)%i
             j=THETA_nonzero_1D_B(index_nonzero)%j
             tempx=THETA_nonzero_1D_B(index_nonzero)%value
             Pr(i,k)=Pr(i,k)+ds*tempx*(2.0*qBstar(j,k,NB-1)-qBstar(j,k,NB))
           enddo
# endif /* MaierSaupe */
            do i=0,M_Bar_B-1
            Pr(i,k)=Pr(i,k)+2.0*qBstar(i,k,NB-1)-0.5*qBstar(i,k,NB) - &
                      ds*Wtemp*(2.0*qBstar(i,k,NB-1)-qBstar(i,k,NB))
            enddo
        enddo !enddo k=0,LOCAL_SIZE    
   endif

   if(s<=NMAX-3) then
          temp=1.0/3.0
     if(s>NA) then
          Pr=0.0
      do k=0,LOCAL_SIZE-1
          Wtemp=WB(k)
# if defined (MaierSaupe)
          do index_nonzero=THETA_nonzero_1D_B_indexK(k-1)+1, &
            THETA_nonzero_1D_B_indexK(k)
             i=THETA_nonzero_1D_B(index_nonzero)%i
             j=THETA_nonzero_1D_B(index_nonzero)%j
             tempx=THETA_nonzero_1D_B(index_nonzero)%value
             Pr(i,k)=Pr(i,k)+ds*tempx*(3.0*qBstar(j,k,s-NA+1)- &
                       3.0*qBstar(j,k,s-NA+2)+qBstar(j,k,s-NA+3))
            
           enddo
# endif /* MaierSaupe */
           do i=0,M_Bar_B-1
          Pr(i,k)=Pr(i,k)+3.0*qBstar(i,k,s-NA+1)-1.5*qBstar(i,k,s-NA+2)+ &
                    temp*qBstar(i,k,s-NA+3) - &
          ds*Wtemp*(3.0*qBstar(i,k,s-NA+1)-3.0*qBstar(i,k,s-NA+2)+qBstar(i,k,s-NA+3))
           enddo
        enddo !enddo k=0,LOCAL_SIZE    
!!!!s<=Nmax-3.and s<=NA!
    else  ! else s<=NA
    !!!calculate polymer block A
             Pr=0.0
             ! else s>=3:s>NA:s==NA+1
     if(s==NA) then
       do k=0,LOCAL_SIZE-1
          Wtemp=WA(k)
# if defined (MaierSaupe)
          do index_nonzero=THETA_nonzero_1D_A_indexK(k-1)+1, &
            THETA_nonzero_1D_A_indexK(k)
             i=THETA_nonzero_1D_A(index_nonzero)%i
             j=THETA_nonzero_1D_A(index_nonzero)%j
             tempx=THETA_nonzero_1D_A(index_nonzero)%value
             tempA1=0.0
             tempA2=0.0
             tempA3=0.0

             if(j<M_Bar_B) then
              tempA1=qBstar(j,k,1)
              tempA2=qBstar(j,k,2)
              tempA3=qBstar(j,k,3)
              endif
              Pr(i,k)=Pr(i,k) + ds*tempx*(3.0*tempA1-3.0*tempA2+tempA3)
           enddo  !enddo index_nonzero
# endif /* MaierSaupe */
            do i=0,M_Bar_A-1
                  tempA1=0.0
                  tempA2=0.0
                  tempA3=0.0
                  if(i<M_Bar_B) then
                  tempA1=qBstar(i,k,1)
                  tempA2=qBstar(i,k,2)
                  tempA3=qBstar(i,k,3)
                  endif
            Pr(i,k)=Pr(i,k)+3.0*tempA1-1.5*tempA2+ &
                      temp*tempA3 - &
                      ds*Wtemp*(3.0*tempA1-3.0*tempA2+tempA3)
            enddo
        enddo !enddo k=0,LOCAL_SIZE    
     else if(s==(NA-1)) then
       do k=0,LOCAL_SIZE-1
          Wtemp=WA(k)
# if defined (MaierSaupe)
          do index_nonzero=THETA_nonzero_1D_A_indexK(k-1)+1, &
            THETA_nonzero_1D_A_indexK(k)
             i=THETA_nonzero_1D_A(index_nonzero)%i
             j=THETA_nonzero_1D_A(index_nonzero)%j
             tempx=THETA_nonzero_1D_A(index_nonzero)%value
             tempA2=0.0
             tempA3=0.0
             if(j<M_Bar_B) then
              tempA2=qBstar(j,k,1)
              tempA3=qBstar(j,k,2)
              endif
              Pr(i,k)=Pr(i,k) + ds*tempx*(3.0*qAstar(j,k,NA)-3.0*tempA2+tempA3)
           enddo  !enddo index_nonzero
# endif /* MaierSaupe */
            do i=0,M_Bar_A-1
                  tempA2=0.0
                  tempA3=0.0
                  if(i<M_Bar_B) then
                  tempA2=qBstar(i,k,1)
                  tempA3=qBstar(i,k,2)
                  endif
            Pr(i,k)=Pr(i,k)+3.0*qAstar(i,k,NA)-1.5*tempA2+ &
                      temp*tempA3 - &
                      ds*Wtemp*(3.0*qAstar(i,k,NA)-3.0*tempA2+tempA3)
            enddo
          enddo !enddo k=0,LOCAL_SIZE    
     else if(s==(NA-2)) then
      do k=0,LOCAL_SIZE-1
          Wtemp=WA(k)
# if defined (MaierSaupe)
          do index_nonzero=THETA_nonzero_1D_A_indexK(k-1)+1, &
            THETA_nonzero_1D_A_indexK(k)
             i=THETA_nonzero_1D_A(index_nonzero)%i
             j=THETA_nonzero_1D_A(index_nonzero)%j
             tempx=THETA_nonzero_1D_A(index_nonzero)%value
             tempA3=0.0
             if(j<M_Bar_B) then
              tempA3=qBstar(j,k,1)
              endif
              Pr(i,k)=Pr(i,k) + ds*tempx*(3.0*qAstar(j,k,NA-1)-3.0*qAstar(j,k,NA)+tempA3)
           enddo  !enddo index_nonzero
# endif /* MaierSaupe */
            do i=0,M_Bar_A-1
                  tempA3=0.0
                  if(i<M_Bar_B) then
                  tempA3=qBstar(i,k,1)
                  endif
            Pr(i,k)=Pr(i,k)+3.0*qAstar(i,k,NA-1)-1.5*qAstar(i,k,NA)+ &
                      temp*tempA3 - &
                      ds*Wtemp*(3.0*qAstar(i,k,NA-1)-3.0*qAstar(i,k,NA)+tempA3)
            enddo
          enddo !enddo k=0,LOCAL_SIZE    
     else !else s>=3,s>NA:s>NA+3
      do k=0,LOCAL_SIZE-1
          Wtemp=WA(k)
# if defined (MaierSaupe)
          do index_nonzero=THETA_nonzero_1D_A_indexK(k-1)+1, &
            THETA_nonzero_1D_A_indexK(k)
             i=THETA_nonzero_1D_A(index_nonzero)%i
             j=THETA_nonzero_1D_A(index_nonzero)%j
             tempx=THETA_nonzero_1D_A(index_nonzero)%value
             Pr(i,k)=Pr(i,k) + ds*tempx*(3.0*qAstar(j,k,s+1)-3.0*qAstar(j,k,s+2)+ &
                             qAstar(j,k,s+3))
           enddo  !enddo index_nonzero
# endif /* MaierSaupe */
            do i=0,M_Bar_A-1
            Pr(i,k)=Pr(i,k)+3.0*qAstar(i,k,s+1)-1.5*qAstar(i,k,s+2)+ &
                      temp*qAstar(i,k,s+3) - &
            ds*Wtemp*(3.0*qAstar(i,k,s+1)-3.0*qAstar(i,k,s+2)+qAstar(i,k,s+3))
            enddo
          enddo !enddo k=0,LOCAL_SIZE    
     endif  !if(s==NA+1) 
  endif   !endif s<=NA
endif    !endif s>=3


 end subroutine init_pr_backward
