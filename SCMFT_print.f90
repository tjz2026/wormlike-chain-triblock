!-------------------------------------------------------------------------
! source  : scmft_print.f90
! type    : subroutines
! author  : Jiuzhou Tang (email:tangjiuzhou@iccas.ac.cn)
! history : 04/06/2013 by Jiuzhou Tang
! purpose : provide printing infrastructure for hybridization expansion
!           version of SCMFT(single chain in mean field theory)

! input   :
! output  :
! status  : very unstable
! comment :
!-------------------------------------------------------------------------

!>>> print the startup information for SCMFT
  subroutine SCMFT_print_header()
     use nrtype
     use constants
     use control 
     use mpi

     implicit none

! string for current date and time
     character (len = 20) :: date_time_string

! obtain current date and time
     call SCMFT_time_builder(date_time_string)

     write(mystd,'(2X,a)') 'wormlike chain SCFT'
     write(mystd,'(2X,a)') '>>> A SCFT program for wormlike chain'
     write(mystd,*)

     write(mystd,'(2X,a)') 'fortran version: 2013.05.21T '//'(built at '//__TIME__//" "//__DATE__//')'
     write(mystd,'(2X,a)') 'develop: by JZ Tang, ICCAS & BNU'
     write(mystd,'(2X,a)') 'support: tangjiuzhou@iccas.ac.cn'
     write(mystd,'(2X,a)') 'license: GPL2 and later versions'
     write(mystd,*)

     write(mystd,'(2X,a)') 'SCFT >>> start running at '//date_time_string

# if defined (MPI)

     write(mystd,'(2X,a,i4)') 'SCFT >>> parallelism: Yes >>> processors:', nprocs

# else   /* MPI */

     write(mystd,'(2X,a,i4)') 'SCFT >>> parallelism: No  >>> processors:', 1

# endif  /* MPI */

     write(mystd,*)

     return
  end subroutine SCMFT_print_header

!>>> print the ending information for SCMFT
  subroutine SCMFT_print_footer()
     use nrtype
     use constants

     implicit none

! string for current date and time
     character (len = 20) :: date_time_string

! used to record the time usage information
     real(dp) :: tot_time

! obtain time usage information
     call cpu_time(tot_time)

! obtain current date and time
     call SCMFT_time_builder(date_time_string)

     write(mystd,'(2X,a,f10.2,a)') 'SCMFT >>> total time spent:', tot_time, 's'
     write(mystd,*)

     write(mystd,'(2X,a)') 'SCMFT >>> Master,I am tired and want to go to bed. Bye!'
     write(mystd,'(2X,a)') 'SCMFT >>> happy ending at '//date_time_string

     return
  end subroutine SCMFT_print_footer

subroutine SCMFT_time_builder(date_time_string)
     implicit none

! external arguments
! output date and time
     character (len = 20), intent(out) :: date_time_string

! local variables
! used to extract data from a standard fortran call: date_and_time()
     integer :: date_time(8)

! string for current date
     character (len = 12) :: cdate

! string for current time
     character (len = 08) :: ctime

! month array
     character (len = 03) :: months(12)

! init the month array
     months( 1) = 'Jan'; months( 2) = 'Feb'; months( 3) = 'Mar'
     months( 4) = 'Apr'; months( 5) = 'May'; months( 6) = 'Jun'
     months( 7) = 'Jul'; months( 8) = 'Aug'; months( 9) = 'Sep'
     months(10) = 'Oct'; months(11) = 'Nov'; months(12) = 'Dec'

! obtain current date and time
     call date_and_time(values = date_time)

! convert date and time from integer to string
     write (cdate,'(1X,a3,1X,i2,1X,i4)') months(date_time(2)), date_time(3), date_time(1)
     write (ctime, '(i2,":",i2,":",i2)') date_time(5), date_time(6), date_time(7)

! build final output string by concating them
     date_time_string = ctime // cdate

     return
  end subroutine SCMFT_time_builder
