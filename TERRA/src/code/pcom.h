*     This include file contains data pertaining to the
*     parallel (MPI and PVM) implementation.
*
      INCLUDE           'mpif.h'
 
      INTEGER            MAX_PROCS
      INTEGER            MPI_DATATYPE
      PARAMETER        ( MAX_PROCS = 1024 )
      PARAMETER        ( MPI_DATATYPE = MPI_DOUBLE_PRECISION )
c     PARAMETER        ( MPI_DATATYPE = MPI_REAL )
 
      INTEGER            MSG_REQUEST, MSG_STATUS( MPI_STATUS_SIZE )
      INTEGER            MYNUM, NPROC, MSGNUM, LVPROC, MPROC, IERROR,
     &                   SEND_REQUESTS( 2 ), RECV_REQUESTS( 2 ),
     &                   SEND_STATUS( MPI_STATUS_SIZE,2 ),
     &                   RECV_STATUS( MPI_STATUS_SIZE,2 )
      COMMON             /pcom/ MYNUM, NPROC, MSGNUM, LVPROC,
     &                          MPROC, MSG_REQUEST, MSG_STATUS,
     &                          SEND_REQUESTS, RECV_REQUESTS,
     &                          SEND_STATUS, RECV_STATUS

c     ========================== End pcom.h ==========================
