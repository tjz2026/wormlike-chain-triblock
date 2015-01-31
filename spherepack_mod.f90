

module spherepack
     implicit none
      
      private 
      integer nnlat,nnlon,nn15,llsave,llwork,lldwork
! integer variables are defined as follows for calling shaec ans shsec only:
!      parameter (nnlat=19,nnlon=36)
!     set saved and unsaved work space lengths in terms of nnlat,nnlon
!     (see documentation for shaec,shsec,islapec)
!      parameter (nn15=nnlon+15)
!      parameter (llsave=nnlat*(nnlat+1)+3*((nnlat-2)*(nnlat-1)+nn15))
!      parameter (llwork=nnlat*(2*nnlon+3*(nnlat+1)+2*nnlat+1))
!     set double precision work space length for initializations
!      parameter (lldwork = nnlat+1)
      
      public unit_vec
      type vec_xyz
      real*8 :: x 
      real*8 :: y
      real*8 :: z 
      end type 
      type(vec_xyz),dimension(:,:),allocatable,save :: unit_vec 
      double precision,dimension(:),allocatable,save :: sint,cost,sinp,cosp

      double precision, dimension(:), allocatable, save :: wshaec,wshsec
      public init_spherepack_parameter,sphscalar_init,sph_scalar_analysis, &
             sph_scalar_synthesis,sphscalar_clean,simposon_sph_integral,&
             simposon_sph_spa3d_mpi,simposon_sph_spa2d_mpi             
contains

          subroutine init_spherepack_parameter()
          use global_para,only: Ntheta,Nphi
          implicit none
          integer :: i,j,k
          double precision pi,x,y,z,dlat,dlon,theta,phi

          nnlat=Ntheta+1
          nnlon=Nphi
          nn15=nnlon+15
          llsave=nnlat*(nnlat+1)+3*((nnlat-2)*(nnlat-1)+nn15)
          llwork=nnlat*(2*nnlon+3*(nnlat+1)+2*nnlat+1)
          lldwork = nnlat+1 
          ! generate some useful arrays
          allocate(sint(1:Ntheta+1)) 
          allocate(cost(1:Ntheta+1)) 
          allocate(sinp(1:Nphi))
          allocate(cosp(1:Nphi))
          allocate(unit_vec(1:Ntheta+1,1:Nphi))
       

      pi = 4.0d0*datan(1.0d0)
      !set sine and cosine vectors
      dlat = pi/(Ntheta)
      dlon = (pi+pi)/Nphi
      do i=1,Ntheta+1
       theta = -0.5d0*pi+(i-1)*dlat
       sint(i) = dsin(theta)
       cost(i) = dcos(theta)
      end do
      do j=1,Nphi
       phi = (j-1)*dlon
       sinp(j) = dsin(phi)
       cosp(j) = dcos(phi)
      end do
!     set right hand side as helmholtz operator
!     applied to ue = (1.+x*y)*exp(z)
      do j=1,Nphi
       phi = (j-1)*dlon
       sinp(j) = dsin(phi)
       cosp(j) = dcos(phi)
      end do
!     set right hand side as helmholtz operator
!     applied to ue = (1.+x*y)*exp(z)
      do j=1,Nphi
       do i=1,Ntheta+1
	x = cost(i)*cosp(j)
	y = cost(i)*sinp(j)
	z = sint(i)
	unit_vec(i,j)%x=x
	unit_vec(i,j)%y=y
	unit_vec(i,j)%z=z
       enddo
      enddo

          end subroutine init_spherepack_parameter
          ! allocate the needed array for scalar analysis and synthesis

          subroutine sphscalar_init(r,a,b)
          implicit none
          double precision, dimension(:,:), intent(in) :: r
          double precision, dimension(:,:), intent(in) :: a,b
          double precision work(llwork)
          double precision dwork(lldwork)
          !double precision a(nnlat,nnlat),b(nnlat,nnlat)
          integer nlat,nlon,lshaec,lshsec,lwork,ierror,isym,nt
          integer ldwork
! first check whether r,a,b are of the correct dimension
        nlat = size(r,1)
        nlon = size(r,2)
        if(nlat/=nnlat .or. nlon/=nnlon) then
        write(*,*) "dimension error in input array r,exit!"
        stop
        endif
        if(nlat/=size(a,1) .or. nlat/=size(a,2)) then
        write(*,*) "dimension error in input array a,exit!"
        stop
        endif
        if(nlat/=size(b,1) .or. nlat/=size(b,2)) then
        write(*,*) "dimension error in input array b,exit!"
        stop
        endif

! allocate workspace array wshaec, wshsec

      if (allocated(wshaec)) deallocate(wshaec)
      allocate(wshaec(llsave))
      if (allocated(wshsec)) deallocate(wshsec)
      allocate(wshsec(llsave))
!     set work space length arguments
      lwork = llwork
      ldwork = lldwork
      lshaec = llsave
      lshsec = llsave
!     set grid size arguments
      nlat = nnlat
      nlon = nnlon
!     set no symmetry and one array
      isym = 0
      nt = 1
            
      CALL SHAECI(NLAT,NLON,WSHAEC,LSHAEC,DWORK,LDWORK,IERROR)
      if (ierror .gt. 0) then
      write (*,*) "error in call shaeci",ierror
      call exit(0)
      end if


      CALL SHSECI(NLAT,NLON,WSHSEC,LSHSEC,DWORK,LDWORK,IERROR)
      if (ierror .gt. 0) then
      write (*,*) "error in call shseci",ierror
      call exit(0)
      end if

          end subroutine sphscalar_init


           subroutine sphscalar_clean()
           implicit none 
           if (allocated(wshaec)) deallocate(wshaec)
           if (allocated(wshsec)) deallocate(wshsec)
           if (allocated(sint)) deallocate(sint)
           if (allocated(cost)) deallocate(cost)
           if (allocated(sinp)) deallocate(sinp)
           if (allocated(cosp)) deallocate(cosp)
           if (allocated(unit_vec)) deallocate(unit_vec)

           end subroutine sphscalar_clean
 
         subroutine sph_scalar_analysis(r,a,b) 
         implicit none
         double precision, dimension(:,:), intent(in) :: r
         double precision, dimension(:,:), intent(out) :: a,b
         double precision work(llwork)
         double precision dwork(lldwork)
         integer nlat,nlon,lshaec,lshsec,lwork,ierror,isym,nt
         integer ldwork
!     set work space length arguments
      lwork = llwork
      ldwork = lldwork
      lshaec = llsave
      lshsec = llsave
!     set grid size arguments
      nlat = nnlat
      nlon = nnlon
!     set no symmetry and one array
      isym = 0
      nt = 1

          call shaec(nlat,nlon,isym,nt,r,nlat,nlon,a,b,nlat,nlat, &
                wshaec,lshaec,work,lwork,ierror)

         end subroutine sph_scalar_analysis


         subroutine sph_scalar_synthesis(r,a,b) 
         implicit none
         double precision, dimension(:,:), intent(out) :: r
         double precision, dimension(:,:), intent(in) :: a,b
         double precision work(llwork)
         double precision dwork(lldwork)
         integer nlat,nlon,lshaec,lshsec,lwork,ierror,isym,nt
         integer ldwork
!     set work space length arguments
      lwork = llwork
      ldwork = lldwork
      lshaec = llsave
      lshsec = llsave
!     set grid size arguments
      nlat = nnlat
      nlon = nnlon
!     set no symmetry and one array
      isym = 0
      nt = 1

         call shsec(nlat,nlon,isym,nt,r,nlat,nlon,a,b,nlat,nlat, &
                         wshsec,lshsec,work,lwork,ierror)

         end subroutine sph_scalar_synthesis
             






!     this file contains a program for solving the Helmholtz
!     equation with constant 1.0 on a ten degree grid on the full sphere
!
! ... required spherepack files
!
!     islapec.f, shaec.f, shsec.f, sphcom.f, hrfft.f
! ... description
!     let theta be latitude and phi be east longitude in radians.
!     and let
!       x = cos(theta)*sin(phi)
!       y = cos(theta)*cos(phi)
!       z = sint(theta)
!
!     be the cartesian coordinates corresponding to theta and phi.
!     on the unit sphere.  The exact solution
!
!        ue(theta,phi) = (1.+x*y)*exp(z)
!
!     is used to set the right hand side and compute error.
!
!
! **********************************************************************
!
! OUTPUT FROM EXECUTING THE PROGRAM BELOW
! WITH 32 AND 64 BIT FLOATING POINT ARITHMETIC
!
! Helmholtz approximation on a ten degree grid
! nlat = 19   nlon = 36
! xlmbda =  1.00   pertrb =  0.000E+00
! maximum error =  0.715E-06 *** (32 BIT)
! maximum error =  0.114E-12 *** (64 BIT)
!
! ***********************************************
! ***********************************************
!      subroutine helmsph()
!!     set grid size with parameter statements
!      implicit none
!      integer nnlat,nnlon,nn15,llsave,llwork,lldwork
!      parameter (nnlat=19,nnlon=36)
!!     set saved and unsaved work space lengths in terms of nnlat,nnlon
!!     (see documentation for shaec,shsec,islapec)
!      parameter (nn15=nnlon+15)
!      parameter (llsave=nnlat*(nnlat+1)+3*((nnlat-2)*(nnlat-1)+nn15))
!      parameter (llwork=nnlat*(2*nnlon+3*(nnlat+1)+2*nnlat+1))
!!     set double precision work space length for initializations
!      parameter (lldwork = nnlat+1)
!!     dimension arrays
!      double precision u(nnlat,nnlon),r(nnlat,nnlon)
!      double precision sint(nnlat),cost(nnlat),sinp(nnlon),cosp(nnlon)
!     ! double precision work(llwork)
!      double precision work(llwork),wshaec(llsave),wshsec(llsave)
!      double precision dwork(lldwork)
!      double precision a(nnlat,nnlat),b(nnlat,nnlat)
!      integer nlat,nlon,i,j,lshaec,lshsec,lwork,ierror,isym,nt
!      integer ldwork
!      double precision pi,x,y,z,dlat,dlon,theta,phi,xlmbda,pertrb,ez,ue,errm
!
!      pi = 4.0d0*atan(1.0d0)
!!     set helmholtz constant
!      xlmbda = 1.0d0
!!     set work space length arguments
!      lwork = llwork
!      ldwork = lldwork
!      lshaec = llsave
!      lshsec = llsave
!!     set grid size arguments
!      nlat = nnlat
!      nlon = nnlon
!!     set sine and cosine vectors
!      dlat = pi/(nlat-1)
!      dlon = (pi+pi)/nlon
!      do i=1,nlat
!       theta = -0.5d0*pi+(i-1)*dlat
!       sint(i) = dsin(theta)
!       cost(i) = dcos(theta)
!      end do
!      do j=1,nlon
!       phi = (j-1)*dlon
!       sinp(j) = dsin(phi)
!       cosp(j) = dcos(phi)
!      end do
!!     set right hand side as helmholtz operator
!!     applied to ue = (1.+x*y)*exp(z)
!      do j=1,nlon
!       do i=1,nlat
!	x = cost(i)*cosp(j)
!	y = cost(i)*sinp(j)
!	z = sint(i)
!	r(i,j) = -(x*y*(z*z+6.*(z+1.))+z*(z+2.))*dexp(z)
!       end do
!      end do
!
!!     initialize saved work space arrays for scalar harmonic
!!     analysis and Helmholtz inversion of r
!!
!      CALL SHAECI(NLAT,NLON,WSHAEC,LSHAEC,DWORK,LDWORK,IERROR)
!      if (ierror .gt. 0) then
!      write (6,200) ierror
!  200 format(' shaeci, ierror = ',i2)
!      call exit(0)
!      end if
!      CALL SHSECI(NLAT,NLON,WSHSEC,LSHSEC,DWORK,LDWORK,IERROR)
!      if (ierror .gt. 0) then
!      write (6,201) ierror
!  201 format(' shseci, ierror = ',i2)
!      call exit(0)
!      end if
!!
!!     set no symmetry and one array
!!
!      isym = 0
!      nt = 1
!!
!!     compute coefficients of r for input to islapec
!      call shaec(nlat,nlon,isym,nt,r,nlat,nlon,a,b,nlat,nlat, &
!                wshaec,lshaec,work,lwork,ierror)
!
!      call shsec(nlat,nlon,isym,nt,r,nlat,nlon,a,b,nlat,nlat, &
!                        wshsec,lshsec,work,lwork,ierror)
!
!      call shaec(nlat,nlon,isym,nt,r,nlat,nlon,a,b,nlat,nlat, &
!                wshaec,lshaec,work,lwork,ierror)
!
!
!      if (ierror .gt. 0) then
!      write(*,202) ierror
!  202 format(' shaec , ierror = ',i2)
!      call exit(0)
!      end if
!!
!!     solve Helmholtz equation on the sphere in u
!!
!      write (6,100) nlat,nlon
!  100 format(' helmholtz approximation on a ten degree grid' &
!            /' nlat = ',i3,2x,' nlon = ', i3)
!      call islapec(nlat,nlon,isym,nt,xlmbda,u,nlat,nlon,a,b,nlat,nlat, &
!                  wshsec,lshsec,work,lwork,pertrb,ierror)
!      if (ierror .ne. 0) then
!      write (6,103) ierror
!  103 format(' islapec, ierror = ',i2)
!      if (ierror .gt. 0) call exit(0)
!      end if
!!
!!     compute and print maximum error in u
!!
!      errm = 0.0
!      do j=1,nlon
!       do i=1,nlat
!	x = cost(i)*cosp(j)
!	y = cost(i)*sinp(j)
!	z = sint(i)
!	ez = dexp(z)
!	ue = (1.+x*y)*ez
!	errm = max(errm,abs(u(i,j)-ue))
!       end do
!      end do
!      write(*,204) xlmbda,pertrb,errm
!  204 format(' xlmbda = ',f5.2,2x, ' pertrb = ' ,e10.3, &
!           /' maximum error = ',e10.3)
!
!      end subroutine helmsph


   FUNCTION simposon_sph_spa3d_mpi(N_theta,N_phi,nz,ny,nx,d_theta,d_phi,dz,dy,dx,a5)
   use global_para, only: Ntheta,Nphi,Local_size,kD_loc
   use grid_mod
   use utility
   implicit none
   integer :: s1,flag,i,j,k,l,m
   REAL*8 :: simposon_sph_spa3d_mpi
   REAL*8 :: FUN_A,FUN_B,FUN2,FUN4,temp
   INTEGER :: nz,ny,nx
   integer,intent(in) :: N_theta,N_phi
   REAL*8 :: dx,dy,dz,d_theta,d_phi
   REAL*8 :: a5(0:Local_size-1,0:N_theta,0:N_phi-1)
   REAL*8 :: a_u(0:N_theta,0:N_phi-1)
   REAL*8 :: a3(0:nz-1,0:ny-1,0:nx-1)
   REAL*8 :: a1(0:nx-1),a2(0:ny-1,0:nx-1),a1u(0:N_phi-1)
!! n1, n2, n3 should be even numbers i.e. 2*N (N is an arbitrary integer larger than zero);
!!      no volume average.  
!!    double a2[n1][n2][n3] --- 0 ~ n1-1 ,  0 ~ n2-1, 0 ~ n3-1  ;  it assumes the periodic boundary condition here.
! Note that nz,ny,nx only represent the first, second, and last dimension of the spatial space, not necessarily in the order of
! Z,Y,X, be aware of that, and if you aren't, please don't change my code in case you may think you find a bug somewhere,OK?



do k=0,N_phi-1
  do l=0,N_theta

    do m=0,Local_size-1
! kD_loc(m)%x represents the first dimension in spatial space, not necessarily in X  axis.
  ! a3(kD_loc(m)%x,kD_loc(m)%y,kD_loc(m)%z)=a5(m,l,k)*sint(l+1) !!! this is wrong, we define different spherical coordinate here
   a3(kD_loc(m)%x,kD_loc(m)%y,kD_loc(m)%z)=a5(m,l,k)*cost(l+1)
    enddo

       FUN_A=0.0d0
       FUN_B=0.0d0
       FUN2=0.0d0
       FUN4=0.0d0

   do i=0,nx-1
       do j=0,ny-1
        FUN_A=a3(0,j,i)
        FUN_B=a3(0,j,i)  !PBC boundary condition !!!!
        FUN2=0.0
        FUN4=a3(1,j,i)
          s1=1
          !note that if nz==2,the 3D SCFT turns into 2D ,and the simposon integral remains
         ! correct.(nz/2-1=0,do inner loop !)
          do flag=1,nz/2-1,1
           s1=s1+1
           FUN2=FUN2+a3(s1,j,i)
           s1=s1+1
           FUN4=FUN4+a3(s1,j,i)
           enddo
          a2(j,i)=(dz/3.0)*(FUN_A+FUN_B+4.0*FUN4+2.0*FUN2)
        enddo
    enddo
     
     do i=0,nx-1
         FUN_A=a2(0,i)
         FUN_B=a2(0,i)
         FUN2=0.0
         FUN4=a2(1,i)
         s1=1
         do flag=1,ny/2-1,1
             s1=s1+1
             FUN2=FUN2+a2(s1,i)
             s1=s1+1
             FUN4=FUN4+a2(s1,i)
          enddo
          a1(i)=(dy/3.0)*(FUN_A+FUN_B+4.0*FUN4+2.0*FUN2)
     enddo
         FUN2=0.0
         FUN4=0.0
        do s1=0,nx-1,2
          FUN2=FUN2+a1(s1)
          FUN4=FUN4+a1(s1+1)
         enddo
     !simposon_3D_mpi=(dx/3.0)*(4.0*FUN4+2.0*FUN2)
      a_u(l,k)=(dx/3.0)*(4.0*FUN4+2.0*FUN2)

 enddo !enddo k=0,N_theta

enddo ! enddo l=1,N_phi-1



! now we integrate the vector u on the unit sphere
    do i=0,N_phi-1
         FUN_A=a_u(0,i)
         FUN_B=a_u(N_theta,i)
         FUN2=0.0
         FUN4=a_u(1,i)
         s1=1
         do flag=1,N_theta/2-1,1
             s1=s1+1
             FUN2=FUN2+a_u(s1,i)
             s1=s1+1
             FUN4=FUN4+a_u(s1,i)
          enddo
          a1u(i)=(d_theta/3.0)*(FUN_A+FUN_B+4.0*FUN4+2.0*FUN2)
     enddo
         FUN2=0.0
         FUN4=0.0
        do s1=0,N_phi-1,2
          FUN2=FUN2+a1u(s1)
          FUN4=FUN4+a1u(s1+1)
         enddo
     simposon_sph_spa3d_mpi=(d_phi/3.0)*(4.0*FUN4+2.0*FUN2)
     return


   end FUNCTION simposon_sph_spa3d_mpi


   FUNCTION simposon_sph_spa2d_mpi(N_theta,N_phi,ny,nx,d_theta,d_phi,dy,dx,a4)
   use global_para, only: Ntheta,Nphi,Local_size,kD_loc
   use grid_mod
   use utility
   implicit none
   integer :: s1,flag,i,j,k,l,m
   REAL*8 :: simposon_sph_spa2d_mpi
   REAL*8 :: FUN_A,FUN_B,FUN2,FUN4,temp
   INTEGER :: ny,nx
   integer,intent(in) :: N_theta,N_phi
   REAL*8 :: dx,dy,d_theta,d_phi
   REAL*8 :: a4(0:Local_size-1,0:N_theta,0:N_phi-1)
   REAL*8 :: a_u(0:N_theta,0:N_phi-1)
   !REAL*8 :: a3(0:nz-1,0:ny-1,0:nx-1)
   REAL*8 :: a1(0:nx-1),a2(0:ny-1,0:nx-1),a1u(0:N_phi-1)
!! n1, n2, n3 should be even numbers i.e. 2*N (N is an arbitrary integer larger than zero);
!!      no volume average.  
!!    double a2[n1][n2][n3] --- 0 ~ n1-1 ,  0 ~ n2-1, 0 ~ n3-1  ;  it assumes the periodic boundary condition here.
! Note that nz,ny,nx only represent the first, second, and last dimension of the spatial space, not necessarily in the order of
! Z,Y,X, be aware of that, and if you aren't, please don't change my code in case you may think you find a bug somewhere,OK?



do k=0,N_phi-1
  do l=0,N_theta

    do m=0,Local_size-1
! kD_loc(m)%x represents the first dimension in spatial space, not necessarily in X  axis.
   !a2(kD_loc(m)%x,kD_loc(m)%y)=a4(m,l,k)*sint(l+1)
   a2(kD_loc(m)%x,kD_loc(m)%y)=a4(m,l,k)*cost(l+1)
    enddo

       FUN_A=0.0d0
       FUN_B=0.0d0
       FUN2=0.0d0
       FUN4=0.0d0
     
     do i=0,nx-1
         FUN_A=a2(0,i)
         FUN_B=a2(0,i)
         FUN2=0.0
         FUN4=a2(1,i)
         s1=1
         do flag=1,ny/2-1,1
             s1=s1+1
             FUN2=FUN2+a2(s1,i)
             s1=s1+1
             FUN4=FUN4+a2(s1,i)
          enddo
          a1(i)=(dy/3.0)*(FUN_A+FUN_B+4.0*FUN4+2.0*FUN2)
     enddo
         FUN2=0.0
         FUN4=0.0
        do s1=0,nx-1,2
          FUN2=FUN2+a1(s1)
          FUN4=FUN4+a1(s1+1)
         enddo
     !simposon_3D_mpi=(dx/3.0)*(4.0*FUN4+2.0*FUN2)
      a_u(l,k)=(dx/3.0)*(4.0*FUN4+2.0*FUN2)

 enddo !enddo k=0,N_theta

enddo ! enddo l=1,N_phi-1

! now we integrate the vector u on the unit sphere
    do i=0,N_phi-1
         FUN_A=a_u(0,i)
         FUN_B=a_u(N_theta,i)
         FUN2=0.0
         FUN4=a_u(1,i)
         s1=1
         do flag=1,N_theta/2-1,1
             s1=s1+1
             FUN2=FUN2+a_u(s1,i)
             s1=s1+1
             FUN4=FUN4+a_u(s1,i)
          enddo
          a1u(i)=(d_theta/3.0)*(FUN_A+FUN_B+4.0*FUN4+2.0*FUN2)
     enddo
         FUN2=0.0
         FUN4=0.0
        do s1=0,N_phi-1,2
          FUN2=FUN2+a1u(s1)
          FUN4=FUN4+a1u(s1+1)
         enddo
     simposon_sph_spa2d_mpi=(d_phi/3.0)*(4.0*FUN4+2.0*FUN2)
     return


   end FUNCTION simposon_sph_spa2d_mpi


   FUNCTION simposon_sph_integral(N_theta,N_phi,d_theta,d_phi,a_u)
   !use global_para, only: Ntheta,Nphi
   use utility
   implicit none
   integer :: s1,flag,i,j,k,l,m
   integer :: N_theta,N_phi
   REAL*8 :: d_theta,d_phi
   REAL*8 :: simposon_sph_integral
   REAL*8 :: FUN_A,FUN_B,FUN2,FUN4,temp
   REAL*8 :: a_u(0:N_theta,0:N_phi-1)
   REAL*8 :: a1u(0:N_phi-1)

! now we integrate the vector u on the unit sphere
  do i=0,N_phi-1
   do j=0,N_theta
        !a_u(j,i)=a_u(j,i)*sint(j+1)
        a_u(j,i)=a_u(j,i)*cost(j+1)
   enddo
  enddo       

  
    do i=0,N_phi-1
         FUN_A=a_u(0,i)
         FUN_B=a_u(N_theta,i)
         FUN2=0.0
         FUN4=a_u(1,i)
         s1=1
         do flag=1,N_theta/2-1,1
             s1=s1+1
             FUN2=FUN2+a_u(s1,i)
             s1=s1+1
             FUN4=FUN4+a_u(s1,i)
          enddo
          a1u(i)=(d_theta/3.0)*(FUN_A+FUN_B+4.0*FUN4+2.0*FUN2)
     enddo
         FUN2=0.0
         FUN4=0.0
        do s1=0,N_phi-1,2
          FUN2=FUN2+a1u(s1)
          FUN4=FUN4+a1u(s1+1)
         enddo
     simposon_sph_integral=(d_phi/3.0)*(4.0*FUN4+2.0*FUN2)
     return



   end FUNCTION simposon_sph_integral

end module spherepack
