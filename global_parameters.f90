
! project : wormlike-chain-SCFT
! program : global_para,mpi
! source  : global_parameters.f90
! type    : module
! author  : Jiuzhou Tang (tangjiuzhou@iccas.ac.cn)
! history : 21/05/2013 by Jz Tang
! purpose : define and declare all the global parameters and the global variables 
!           in the global_para module
! input   : none
! output  : none
! status  : unstable
! comment :
!-------------------------------------------------------------------------


module global_para
use nrtype,only:SP,DP
implicit none
!***********parameters for AB-diblock wormlike chain*****
INTEGER  :: N_blk  ! N block copolymer
REAL(DP) :: chiN(3,3)   !flory-hugins para
REAL(DP) :: NxAB,NxAC,NxBC   !flory-hugins para
REAL(DP) :: fa,fb,fc 
REAL(DP) :: f_blk(3)
Integer(DP) :: Ns_blk(3)
Integer(DP) :: Nstart_blk(3)
INTEGER  :: Nmax,NA,NB,NC !discretize grid number for polymer contour length L
REAL(DP) :: ds
REAL(DP) :: Loa  !L over persistent length
REAL(DP) :: nu   ! Maier-Saupe interaction parameter
REAL(DP),parameter :: kapa_A=0.5d0,kapa_A_inv=2.0d0  ! Lp/a for A
REAL(DP) :: kapa_B,kapa_B_inv
REAL(DP) :: kapa_C,kapa_C_inv
REAL(DP) :: kapa_blk(3)
REAL(DP) :: NMu_NXab     !NMu/NXab
REAL(DP) :: NMu     !NMu=NXab*NMu_NXab
REAL(DP) :: M_initial  !M_initial=NXAb*NMu_NXab*0.1
!*************parameters for spatial grid***************************
INTEGER,parameter ::  Ndim=3 !system dimention,1 for 1D,2 for 2D,3 for 3D
INTEGER ::  SIDEx,SIDEy,SIDEz !spatial grid number,make sure they are equal to 2^N for fftw 
INTEGER :: M_grid !M_grid=SIDEx*SIDEy*SIDEz for 3D case
INTEGER :: LOCAL_SIZE !local size of subdomain
INTEGER :: LOCAL_SIZE_K !local size of subdomian in k space for pencil fft only.
INTEGER :: LOCAL_SIZE_K_R2C !local size of subdomian in k space when R2C FFT is used.
INTEGER :: LOCAL_NX_SIZE !it euqals to local_nx in slab decomposition fftw
REAL(DP) :: Lx,Ly,Lz ! unit cell length for x,y,z
REAL(DP) :: dx,dy,dz,dx_inv,dy_inv,dz_inv
REAL(DP) :: x1,x2,x3 !box size range for brent
REAL(DP) :: M_v  !total volume
!************parameters for spherical harmonic basis******
integer,parameter :: L_bar_A=3
integer ::  M_bar_A,Dim_ordering_A
integer,parameter :: L_bar_B=3
integer ::  M_bar_B,Dim_ordering_B
integer,parameter :: L_bar_C=3
integer ::  M_bar_C,Dim_ordering_C
integer ::  L_bar_blk(3),M_bar_blk(3),Dim_ordering_blk(3)
integer ::  max_M_Bar
REAL(DP),parameter :: Thresh_sprase_matrix=1.0e-10
!!!para for free Energy
REAL(DP):: totden,pff_global,tempEE_global,FE_global,bigQ

!************parameters for spherical harmonic decomposition scheme ******
INTEGER,parameter  :: Ntheta=16,Nphi=32 !discretize grid number for theta and phi angle
REAL(DP),parameter :: dtheta=3.141592653589793d0/16,dphi=2*3.14159265d0/32

!*************parameters for SCFT *********************************
INTEGER,parameter :: MAXITS=5000 !maxmum iteration step number
INTEGER :: iter_method !iteration scheme,1 for simple mixing,2 for anderson_mixing
REAL(DP) :: lambda_simple 
REAL(DP),PARAMETER :: TOL=1.0e-2
!!converge criterion
!REAL(DP),PARAMETER :: CC=2.0e-5
REAL(DP) :: CC,CC_final
logical :: NMu_add,final_run

INTEGER :: n_iter,n_iter_M,n_SCFT
!**************variable for SCFT**************************
!*****for tensor like M,S
! declare the 3*3 tensor for M(r) and S(r) in Maier-Saupe
integer ,parameter :: N_dim_ddm=3

type node
		DOUBLE PRECISION :: xx
		DOUBLE PRECISION :: xy
		DOUBLE PRECISION :: xz
		!DOUBLE PRECISION :: yx
		!DOUBLE PRECISION :: yy
		DOUBLE PRECISION :: yz
		!DOUBLE PRECISION :: zx
		!DOUBLE PRECISION :: zy
		DOUBLE PRECISION :: zz
	end type

type node2
       INTEGER :: x
       INTEGER :: y
end type

type(node2),allocatable :: basis_SPH(:)


REAL(DP),allocatable :: wA(:),wB(:),WC
REAL(DP),allocatable :: wA_save(:,:),wB_save(:,:),wC_save(:,:)
REAL(DP) :: dx_save(20),dx_erro(20) 
REAL(DP),allocatable :: RA(:),RB(:),RC(:)
REAL(DP),allocatable :: R_AB_joint(:),R_BC_joint(:),R_end(:)
REAL(DP),allocatable :: RHOA(:,:),RHOB(:,:),RHOC(:,:)
REAL(DP),allocatable :: M_OP(:,:,:),S_OP(:,:,:)
REAL(DP),allocatable :: SA_OP(:,:,:),SB_OP(:,:,:),SC_OP(:,:)
REAL(DP),allocatable :: Pu_A(:,:),Pu_B(:,:),Pu_C(:,:),Pu(:,:)

REAL(DP),allocatable :: w_blk(:,:)
REAL(DP),allocatable :: R_blk(:,:)
REAL(DP),allocatable :: R_joint(:,:)
REAL(DP),allocatable :: RHO_blk(:,:,:)
REAL(DP),allocatable :: M_OP(:,:,:),S_OP(:,:,:)
REAL(DP),allocatable :: S_OP_blk(:,:,:,:)
REAL(DP),allocatable :: Pu(:,:),Pu_blk(:,:,:)



!*********for Anderson mixing 
integer,parameter :: n_r_WAB=5
!INTEGER,parameter ::  Anderson_nim_AB=5
INTEGER ::  Anderson_nim_AB
integer,parameter :: Num_step_simple_mixing_for_WAB=20
INTEGER,parameter ::  Anderson_nim_M=4
integer,parameter :: n_r_M=4
INTEGER,parameter ::  N_WABITER_PER_M=1
INTEGER,parameter ::  Num_simple_mixing_AB=20
INTEGER,parameter ::  Num_simple_mixing_M=5000

REAL(DP),parameter :: lambda_WAB_anderson=0.1d0
REAL(DP),parameter :: lambda_M_anderson=0.10d0
!!!
REAL(DP),parameter :: lanbtwA=0.10d0
REAL(DP),parameter :: lanbtwB=0.10d0
REAL(DP),parameter :: lanbtM=0.080d0


REAL(DP) :: ta_diff,tb_diff,tm_diff
REAL(DP) :: lambda_WAB_anderson_const
integer ::Num_step_simple_mixing_for_M
REAL(DP),allocatable :: wA_out(:,:),wB_out(:,:)
REAL(DP),allocatable :: dw_A(:,:),dw_B(:,:)
REAL(DP),allocatable :: dA_anderson(:,:),dB_anderson(:,:)

REAL(DP),allocatable :: w_out_blk(:,:,:)
REAL(DP),allocatable :: dw_blk(:,:,:)
REAL(DP),allocatable :: d_anderson_blk(:,:,:)

REAL(DP),allocatable :: d11_anderson(:,:)
REAL(DP),allocatable :: d12_anderson(:,:)
REAL(DP),allocatable :: d13_anderson(:,:)
REAL(DP),allocatable :: d23_anderson(:,:)
REAL(DP),allocatable :: d33_anderson(:,:)
REAL(DP),allocatable :: M11_out(:,:)
REAL(DP),allocatable :: M12_out(:,:)
REAL(DP),allocatable :: M13_out(:,:)
REAL(DP),allocatable :: M23_out(:,:)
REAL(DP),allocatable :: M33_out(:,:)
!*****************************
!***********matrix and tensors related to i,j index of spherical harmonic baisi
REAL(DP),allocatable :: matrix_Rx(:,:)
REAL(DP),allocatable :: matrix_Ry(:,:)
REAL(DP),allocatable :: matrix_Rz(:,:)
REAL(DP),allocatable :: GAMA(:,:,:)
REAL(DP),allocatable :: qA(:,:,:) !balock 0~NA
REAL(DP),allocatable :: qB(:,:,:)
REAL(DP),allocatable :: qC(:,:,:)
REAL(DP),allocatable :: qAstar(:,:,:)
REAL(DP),allocatable :: qBstar(:,:,:)
REAL(DP),allocatable :: qCstar(:,:,:)


REAL(DP),allocatable :: q_f(:,:,:) 
REAL(DP),allocatable :: q_star(:,:,:) 


REAL(DP),allocatable :: THETAij_A(:,:,:)
REAL(DP),allocatable :: THETAij_M11_M22_A(:,:)
REAL(DP),allocatable :: THETAij_M33_A(:,:)
REAL(DP),allocatable :: THETAij_M12_A(:,:)
REAL(DP),allocatable :: THETAij_M13_A(:,:)
REAL(DP),allocatable :: THETAij_M23_A(:,:)
REAL(DP),allocatable :: J11ij_A(:,:)
REAL(DP),allocatable :: J22ij_A(:,:)
REAL(DP),allocatable :: J12ij_A(:,:)
REAL(DP),allocatable :: J13ij_A(:,:)
REAL(DP),allocatable :: J23ij_A(:,:)

REAL(DP),allocatable :: THETAij_blk(:,:,:,:)
REAL(DP),allocatable :: THETAij_M11_M22_blk(:,:,:)
REAL(DP),allocatable :: THETAij_M33_blk(:,:,:)
REAL(DP),allocatable :: THETAij_M12_blk(:,:,:)
REAL(DP),allocatable :: THETAij_M13_blk(:,:,:)
REAL(DP),allocatable :: THETAij_M23_blk(:,:,:)
REAL(DP),allocatable :: J11ij_blk(:,:,:)
REAL(DP),allocatable :: J22ij_blk(:,:,:)
REAL(DP),allocatable :: J12ij_blk(:,:,:)
REAL(DP),allocatable :: J13ij_blk(:,:,:)
REAL(DP),allocatable :: J23ij_blk(:,:,:)

! define the real and image part of G matrix separately
REAL(DP) ,allocatable :: R_G3_inv_val_blk(:,:,:)
INTEGER,allocatable :: R_G3_inv_col_blk(:,:,:)
INTEGER,allocatable :: R_G3_inv_row_blk(:,:,:)

REAL(DP) ,allocatable :: I_G3_inv_val_blk(:,:,:)
INTEGER,allocatable :: I_G3_inv_col_blk(:,:,:)
INTEGER,allocatable :: I_G3_inv_row_blk(:,:,:)

REAL(DP) ,allocatable :: R_G3star_inv_val_blk(:,:,:)
INTEGER,allocatable :: R_G3star_inv_col_blk(:,:,:)
INTEGER,allocatable :: R_G3star_inv_row_blk(:,:,:)

REAL(DP) ,allocatable :: I_G3start_inv_val_blk(:,:,:)
INTEGER,allocatable :: I_G3star_inv_col_blk(:,:,:)
INTEGER,allocatable :: I_G3star_inv_row_blk(:,:,:)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
REAL(DP) ,allocatable :: R_GA1_inverse_val(:,:)
INTEGER,allocatable :: R_GA1_inverse_col(:,:)
INTEGER,allocatable :: R_GA1_inverse_row(:,:)

REAL(DP) ,allocatable :: R_GA1_star_inverse_val(:,:)
INTEGER,allocatable :: R_GA1_star_inverse_col(:,:)
INTEGER,allocatable :: R_GA1_star_inverse_row(:,:)
! image part
REAL(DP) ,allocatable :: I_GA1_inverse_val(:,:)
INTEGER,allocatable :: I_GA1_inverse_col(:,:)
INTEGER,allocatable :: I_GA1_inverse_row(:,:)

REAL(DP) ,allocatable :: I_GA2_inverse_val(:,:)
INTEGER,allocatable :: I_GA2_inverse_col(:,:)
INTEGER,allocatable :: I_GA2_inverse_row(:,:)

REAL(DP) ,allocatable :: I_GA3_inverse_val(:,:)
INTEGER,allocatable :: I_GA3_inverse_col(:,:)
INTEGER,allocatable :: I_GA3_inverse_row(:,:)

REAL(DP) ,allocatable :: I_GA3_star_inverse_val(:,:)
INTEGER,allocatable :: I_GA3_star_inverse_col(:,:)
INTEGER,allocatable :: I_GA3_star_inverse_row(:,:)

REAL(DP) ,allocatable :: I_GB3_inverse_val(:,:)
INTEGER,allocatable :: I_GB3_inverse_col(:,:)
INTEGER,allocatable :: I_GB3_inverse_row(:,:)

REAL(DP) ,allocatable :: I_GB3_star_inverse_val(:,:)
INTEGER,allocatable :: I_GB3_star_inverse_col(:,:)
INTEGER,allocatable :: I_GB3_star_inverse_row(:,:)

REAL(DP) ,allocatable :: I_GB2_star_inverse_val(:,:)
INTEGER,allocatable :: I_GB2_star_inverse_col(:,:)
INTEGER,allocatable :: I_GB2_star_inverse_row(:,:)

REAL(DP) ,allocatable :: I_GB1_star_inverse_val(:,:)
INTEGER,allocatable :: I_GB1_star_inverse_col(:,:)
INTEGER,allocatable :: I_GB1_star_inverse_row(:,:)

! define the sparse array for G Matrix
REAL(DP) ,allocatable :: R_G_val(:)
REAL(DP) ,allocatable :: I_GA1x_val(:)
INTEGER,allocatable ::   I_GA1x_col(:)
INTEGER,allocatable ::   I_GA1x_row(:)
REAL(DP) ,allocatable :: I_GA1y_val(:)
INTEGER,allocatable ::   I_GA1y_col(:)
INTEGER,allocatable ::   I_GA1y_row(:)
REAL(DP) ,allocatable :: I_GA1z_val(:)
INTEGER,allocatable ::   I_GA1z_col(:)
INTEGER,allocatable ::   I_GA1z_row(:)

type node3
      integer :: i
      integer :: j
      real(DP) :: value
end type
type node4
      integer :: i
      integer :: j
      integer :: k
      real(DP) :: value
end type

type node5
     integer :: i
     integer :: j
end type

# if defined (Dim2)
type node6
     integer :: x
     integer :: y
end type
# else  /* Dim2 */
type node6
     integer :: x
     integer :: y
     integer :: z
end type
# endif /* Dim2 */


type node7
     integer :: i
     integer :: j
     integer :: k
end type
!type(node3),allocatable :: THETA_nonzero_2D_A(:,:)
!type(node3),allocatable :: THETA_nonzero_2D_B(:,:)
type(node6),allocatable :: kD(:)
type(node6),allocatable :: kD_loc(:)
type(node6),allocatable :: kDz(:)
type(node6),allocatable :: kDk(:)
type(node7),allocatable :: k_index_1D(:)

type(node3),allocatable :: THETA_nonzero_1D_A(:)
type(node3),allocatable :: THETA_nonzero_1D_B(:)
INTEGER,allocatable :: THETA_nonzero_1D_A_indexK(:)
INTEGER,allocatable :: THETA_nonzero_1D_B_indexK(:)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
type(node3),allocatable :: THETA_nonzero_1D_blk(:,:)
INTEGER,allocatable :: THETA_nonzero_1D_indexK_blk(:,:)
type(node3),allocatable :: J11ij_nonzero_1D_blk(:,:)
type(node3),allocatable :: J22ij_nonzero_1D_blk(:,:)
type(node3),allocatable :: J12ij_nonzero_1D_blk(:,:)
type(node3),allocatable :: J13ij_nonzero_1D_blk(:,:)
type(node3),allocatable :: J23ij_nonzero_1D_blk(:,:)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
type(node3),allocatable :: J11ij_nonzero_1D_A(:)
type(node3),allocatable :: J22ij_nonzero_1D_A(:)
type(node3),allocatable :: J12ij_nonzero_1D_A(:)
type(node3),allocatable :: J13ij_nonzero_1D_A(:)
type(node3),allocatable :: J23ij_nonzero_1D_A(:)

type(node3),allocatable :: J11ij_nonzero_1D_B(:)
type(node3),allocatable :: J22ij_nonzero_1D_B(:)
type(node3),allocatable :: J12ij_nonzero_1D_B(:)
type(node3),allocatable :: J13ij_nonzero_1D_B(:)
type(node3),allocatable :: J23ij_nonzero_1D_B(:)


type(node3),allocatable :: Rx_nonzero_1D(:)
type(node3),allocatable :: Ry_nonzero_1D(:)
type(node3),allocatable :: Rz_nonzero_1D(:)
integer :: Rx_num,Ry_num,Rz_num

integer :: J11_A,J12_A,J13_A,J22_A,J23_A
integer :: J11_B,J12_B,J13_B,J22_B,J23_B
integer :: J11_B,J12_blk(3),J13_blk(3),J22_blk(3),J23_blk(3)
type(node4),allocatable :: GAMA_nonzero_1D(:)

!type(node5),allocatable :: counter(:)
!different types to init the w field 
integer :: wAB_init_type  
integer :: wABC_init_type  
 !6 read from input file named wfield.txt
 !0 for random
 !1 for FCC initial 3D
 !2 for  BCC initial 3D
 !3 for Gyroid initial 3D
 !4 for Cylinder initial 3D 
 !5 for Lamellar initial 3D 
 !6 for P4 initial 3D 
character,parameter :: method='C' 
!!! for brent minimalization
!!????????is tol=1.0e-2 too big ?
!! only to store the sparse matirxi G_R_diag,G_I
REAL(DP),allocatable :: G_RR(:,:)
REAL(DP),allocatable :: G_II(:,:)

end module global_para
   module mpi
  include 'mpif.h'
  end module mpi
