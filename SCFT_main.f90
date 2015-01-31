!-------------------------------------------------------------------------
! project : wormlike-chain-SCFT
! program : main
! source  : SCFT_main.f90
! type    : main program
! author  : Jiuzhou Tang (tangjiuzhou@iccas.ac.cn)
! history : 21/05/2013 by Jz Tang
! purpose : the main program of wormlike-chain-SCFT
! input   : none
! output  : none
! status  : unstable
! comment :
!-------------------------------------------------------------------------
program main
use nrtype,only : DP
use global_para
use constants
use control
use mpi
use mmpi
implicit none
REAL*8 :: minFE,minx
REAL(DP),external :: cal_scft
REAL(DP),external :: brent
REAL(DP) :: t1,t2


t1=MPI_Wtime()
call initialize()
# if defined (Mem_OP)
final_run=.false.
# endif /* Mem_OP */

n_SCFT=0
CC=2.0e-4
!minFE=brent(x1,x2,x3,cal_scft,tol,minx)

# if defined (Mem_OP)
final_run=.true.
# endif /* Mem_OP */
CC=1.0e-4

!minFE=cal_scft(minx)
minFE=cal_scft(x2)

call global_arrays_clean()

if(myid == master) then
call SCMFT_print_footer()
endif
t2=MPI_Wtime()

if(myid == master) then
write(*,*) "total running time",t2-t1
endif

call mp_barrier()
call mp_finalize()
stop
end

















