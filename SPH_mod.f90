! a compiling bug exsits when using -DSPH option, needs to be fixed some day,currentlt,
! when use SPH module, don't add -DSPH compiling flag,OK?
# if defined (SPH)
# else /* SPH */
module SPH
implicit none  


contains

       FUNCTION triple_product(l1,l2,l3,m1,m2,m3)
       !this subroutine should better use the case statement than
       !the if statement,but,what the hell,I'll change it later.
       implicit none
       integer :: l1,l2,l3,m1,m2,m3
       real*8 trpd,temp,suma
       real*8 sqrt1,sqrt2
       real*8 triple_product
        
           !write(*,*) "start trpd"
           !write(*,*) "l1,l2,l3,",l1,l2,l3
           !write(*,*) "m1,m2,m3",m1,m2,m3


        sqrt1=sqrt(1.0*(2*l1+1)*(2*l2+1)*(2*l3+1))
        sqrt2=sqrt1/sqrt(2.0)
 
        trpd=0.0
       if(sign_ov(l1+l2+l3)==-1) then
            trpd=0.0      
       else
              if(m1==0 .and. m2==0 .and. m3==0) then  ! 1
               temp=three_j(l1,l2,l3,0,0,0)
               trpd=sqrt(1.0*(2*l1+1)*(2*l2+1)*(2*l3+1))*temp*temp

               else if (m1>0 .and. m2>0 .and. m3>0) then  !2
               temp=sqrt(1.0*(2*l1+1)*(2*l2+1)*(2*l3+1)/2.0)
               temp=temp*three_j(l1,l2,l3,0,0,0)
               suma=sign_ov(m1)*three_j(l1,l2,l3,-m1,m2,m3) +sign_ov(m2)*three_j(l1,l2,l3,m1,-m2,m3) + &
               sign_ov(m3)*three_j(l1,l2,l3,m1,m2,-m3)
               trpd=temp*suma

               else if(m1==0 .and. m2>0 .and. m3>0) then !3
               temp=sign_ov(m2)*sqrt(1.0*(2*l1+1)*(2*l2+1)*(2*l3+1))
               trpd=temp*three_j(l1,l2,l3,0,0,0)*three_j(l1,l2,l3,0,-m2,m3)

               else if (m1>0 .and. m2==0 .and. m3>0) then !4
               temp=sign_ov(m1)*sqrt(1.0*(2*l1+1)*(2*l2+1)*(2*l3+1))
               trpd=temp*three_j(l1,l2,l3,0,0,0)*three_j(l1,l2,l3,-m1,0,m3)

               else if(m1>0 .and. m2>0 .and. m3==0) then !5
               temp=sign_ov(m1)*sqrt(1.0*(2*l1+1)*(2*l2+1)*(2*l3+1))
                trpd=temp*three_j(l1,l2,l3,0,0,0)*three_j(l1,l2,l3,-m1,m2,0) 
                
                else if(m1==0 .and. m2<0 .and. m3<0) then !6
                 temp=sign_ov(m2)*sqrt1
                 trpd=temp*three_j(l1,l2,l3,0,0,0)*three_j(l1,l2,l3,0,m2,-m3)
                
                 else if(m1<0 .and. m2==0 .and. m3<0) then !7
                 temp=sign_ov(m1)*sqrt1
                 trpd=temp*three_j(l1,l2,l3,0,0,0)*three_j(l1,l2,l3,m1,0,-m3)
               
                 else if(m1<0 .and. m2<0 .and. m3==0) then !8
                  temp=sign_ov(m1)*sqrt1
                  trpd=temp*three_j(l1,l2,l3,0,0,0)*three_j(l1,l2,l3,m1,-m2,0)
                    
                 else if(m1>0 .and. m2<0 .and. m3<0) then !9
                   temp=sqrt2
                   temp=temp*three_j(l1,l2,l3,0,0,0)
                   suma=sign_ov(m2)*three_j(l1,l2,l3,m1,m2,-m3)+sign_ov(m3)* &
                   three_j(l1,l2,l3,m1,-m2,m3) -sign_ov(m1)*three_j(l1,l2,l3,m1,m2,m3)
                   trpd=temp*suma

                 else if(m1<0 .and. m2>0 .and. m3<0) then !10
                   temp=sqrt2
                   temp=temp*three_j(l1,l2,l3,0,0,0)
                   suma=sign_ov(m1)*three_j(l1,l2,l3,m1,m2,-m3)-sign_ov(m2)* &
                   three_j(l1,l2,l3,m1,m2,m3) + sign_ov(m3)*three_j(l1,l2,l3,-m1,m2,m3)
                   trpd=temp*suma
                  
                   
                 else if(m1<0 .and. m2<0 .and. m3>0) then !11
                   temp=sqrt2
                   temp=temp*three_j(l1,l2,l3,0,0,0)
                   suma=sign_ov(m1)*three_j(l1,l2,l3,m1,-m2,m3)+sign_ov(m2)* &
                   three_j(l1,l2,l3,-m1,m2,m3) - sign_ov(m3)*three_j(l1,l2,l3,m1,m2,m3)
                   trpd=temp*suma

                else
                   trpd=0.0
               endif
   
            endif

           triple_product=trpd
           !write(*,*) "trpd=",triple_product 
           return 
           end FUNCTION triple_product                   

           !function to calculate three j symbol
           FUNCTION three_j(l1,l2,l3,m1,m2,m3)
           implicit none
           integer :: l1,l2,l3,m1,m2,m3
           integer,parameter :: FACTVALUE=20
           integer :: k,k_min,k_max
           real*8 :: temp1,temp2
           real*8 :: temp,tempk
           real*8 :: suma
           real*8 :: three_j
           
         if((m1+m2+m3)/=0 .or. (l1+l2-l3)<0 .or. (l1-l2+l3)<0 .or. &
             (l2+l3-l1)<0 .or. abs(l1)<abs(m1) .or. abs(l2)<abs(m2) .or. &
             abs(l3)<abs(m3)) then
             three_j=0.0
          else 
               if((l1+l2+L3+1)>FACTVALUE) then
                  temp1=1.0/fact(FACTVALUE)
                  temp2=1.0/fact_con(FACTVALUE+1,l1+l2+l3+1)
                else
                   temp1=1.0/fact(l1+l2+l3+1)
                   temp2=1.0
                endif
          temp=fact(l1+l2-l3)*fact(l1-l2+l3)*fact(l2+l3-l1)*temp1
          temp=temp*fact(l1-m1)*fact(l1+m1)*fact(l2-m2)*fact(m2+l2)*fact(l3-m3)*fact(l3+m3)*temp2
          temp=sign_ov(l1-l2-m3)*sqrt(temp)
          k_min=max(0,l2-l3-m1,l1-l3+m2)
          k_max=min(l1+l2-l3,l1-m1,l2+m2)
          suma=0.0
          do k=k_min,k_max
            tempk=temp/fact(k)
            tempk=tempk/fact(l1+l2-l3-k)
            tempk=tempk/fact(l1-m1-k)
            tempk=tempk/fact(l2+m2-k)
            tempk=tempk/fact(l3-l2+m1+k)
            tempk=tempk/fact(l3-l1-m2+k)
            suma=suma+sign_ov(k)*tempk
         enddo
             three_j=suma
           
       endif                 
       return
      end function three_j



           FUNCTION sign_ov(i)
           implicit none
           integer i
           integer sign_ov
           if(mod(i,2)==0) then
           sign_ov=1
           else
           sign_ov=-1
           endif
           return
           end function sign_ov

            FUNCTION fact(l)
            implicit none
            integer i,l
            real*8 fact
            fact=1.0d0
            if(l==0 .or. l==1) then
                fact=1.0
            else
             do i=2,l
             fact=fact*i
             enddo
             endif
             return
             end function fact
 
             
            FUNCTION fact_con(m1,m2)
            implicit none
            integer i,m1,m2
            real*8 fact_con
            fact_con=1.0d0

            if(m2<m1) then
             write(*,*) "error :m2<m1"
             stop
             endif

             do i=m1,m2
             fact_con=fact_con*i
             enddo
             return
             end function fact_con


           subroutine basis_sph_create()
           use global_para,only: L_bar_A,L_Bar_B,basis_SPH, & 
           M_Bar_B,M_Bar_A
           !use global_para
           USE control
           USE constants
           USE utility
           USE mmpi
           implicit none
           integer :: i,l,m,istat3
           integer :: max_M_Bar,max_L_Bar
              
          max_M_Bar=max(M_Bar_A,M_Bar_B)
          max_L_Bar=max(L_Bar_A,L_Bar_B)
             
          allocate(basis_SPH(0:max_M_Bar-1),stat=istat3)
          if(istat3/=0) then
          write(*,*) "allocate basis_SPH failed!"
          stop
          endif
           
              i=0
           do l=0,max_L_Bar
              do m=-l,l
                  
                  basis_SPH(i)%x=l
                  basis_SPH(i)%y=m
                  i=i+1
                 !write(*,*) "i=",i
               enddo
            enddo
            if((i-1)/=max_M_Bar-1) then
              if(myid==0) then
             write(*,*) "basis_SPH index created wrong!!,exit!","on",myid
             write(*,*) "i-1=",i-1,"M_bar_B-1=",M_bar_B-1
             write(*,*) "basis_SPH(i)%x",basis_SPH(i)%x
              endif
              call mp_barrier()
              call mp_finalize() 
            stop
            endif
            end subroutine basis_sph_create
                    
             
           end module SPH
           
# endif /* SPH */           






