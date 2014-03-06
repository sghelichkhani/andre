!     This include file contains data pertaining to the
!     parallel (MPI and PVM) implementation.

      INCLUDE           'mpif.h'
 
      INTEGER            MAX_PROCS
      INTEGER            MPI_DATATYPE
      PARAMETER        ( MAX_PROCS = 1024 )
      PARAMETER        ( MPI_DATATYPE = MPI_DOUBLE_PRECISION )
!     PARAMETER        ( MPI_DATATYPE = MPI_REAL )
 
      INTEGER            MSG_REQUEST, MSG_STATUS( MPI_STATUS_SIZE )
      INTEGER            ITIDS, ierr, myid2
      INTEGER            MYID, NPROC, MSGNUM, LVPROC, MPROC,&
     &                   SEND_REQUESTS( 2 ), RECV_REQUESTS( 2 ),&
     &                   SEND_STATUS( MPI_STATUS_SIZE,2 ),&
     &                   RECV_STATUS( MPI_STATUS_SIZE,2 )
      COMMON             /pcom/ MYID, NPROC, MSGNUM, LVPROC,&
     &                          MPROC, MSG_REQUEST, MSG_STATUS,&
     &                          SEND_REQUESTS, RECV_REQUESTS,&
     &                          SEND_STATUS, RECV_STATUS, ierr, myid2
      COMMON             /pvm1/ itids(0:1023)
!    ========================== End pcom.h ==========================
