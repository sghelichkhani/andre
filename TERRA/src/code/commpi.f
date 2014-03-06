c -----------------------------------------------------------------------------
*dk pcom1
      subroutine pcom1(coml,comr,crnr,buff,nj,kr,kt,ighost)

c...  This routine broadcasts and receives boundary information
c...  to and from the neighboring subdomains.

      include 'size.h'
      include 'pcom.h'
      real  coml((kt+1)*nd*nj*(kr+1),2)
      real  comr((kt+1)*nd*nj*(kr+1),2)
      real  crnr(nd*nj*(kr + 1))
      real  buff((kt+1)*nd*nj*(kr+1),2)
      integer sender, ibegin, iend, ipending
      common /nbr1/ ileft(4,0:6), iright(4,0:6)

      ibegin   = 1
      iend     = 2
      ipending = 0

c...  Increment the message number.
      msgnum = msgnum + 1

c...  Compute message length along subdomain edges.
      mm = (kt + 1)*nd*nj*(kr + 1)

c...  Broadcast row/column data in coml to two subdomains to the left,
c...  and receive similar information from two subdomains to the right.
c...  In order to avoid self-communication check, if neighoring processor 
c...  number equals mynum. In this case don't send the array coml, but 
c...  rather copy it.

      if(ileft(1,lvproc) .eq. mynum)then
      ibegin=iend
      call scopy(mm,coml(1,1),1,comr(1,1),1)
      endif
      if(ileft(2,lvproc) .eq. mynum)then
      iend=ibegin
      call scopy(mm,coml(1,2),1,comr(1,2),1)
      endif

      do ib=ibegin,iend
         CALL MPI_IRECV( BUFF(1,IB),  MM, MPI_DATATYPE, MPI_ANY_SOURCE,
     &                   MSGNUM, MPI_COMM_WORLD, RECV_REQUESTS( 1 ),
     &                   IERROR )

c...     Count number of pending send operations
         ipending = ipending + 1

         CALL MPI_ISEND( COML( 1, IB ), MM, MPI_DATATYPE,
     &                   ILEFT( IB, LVPROC ), MSGNUM, MPI_COMM_WORLD,
     &                   SEND_REQUESTS( ipending ), IERROR )

c...     Complete the nonblocking receives before processing.
         CALL MPI_WAITANY( 1, RECV_REQUESTS, SENDER, RECV_STATUS,
     &                     IERROR )

c...     Save this information in array comr.
c...     Determine the "iright" value of the isrc_proc subdomain.
         isrc_proc = RECV_STATUS( MPI_SOURCE, 1 )

         do i=1,2
            if(iright(i,lvproc) .eq. isrc_proc)
     &         call scopy(mm, buff(1,ib), 1, comr(1,i), 1)
         end do
      end do

c.... Complete the non-blocking sends before proceeding.
      CALL MPI_WAITALL( iend - ibegin + 1, SEND_REQUESTS, SEND_STATUS,
     &                  IERROR )

c...  Increment the message number.
      msgnum = msgnum + 1

c...  Compute message length for subdomain corner.
      mm = nd*nj*(kr + 1)

c...  Send left corner data to appropriate left subdomain.
      iproc = ileft(3,lvproc)
      if(ighost .eq. 1) iproc = ileft(4,lvproc)

c...  Avoid self-communication
      if(iproc.ne.mynum)then

c...  Receive right corner data to appropriate right subdomain.
      CALL MPI_IRECV( COML(1,1), MM, MPI_DATATYPE, MPI_ANY_SOURCE,
     &                MSGNUM, MPI_COMM_WORLD, RECV_REQUESTS( 1 ),
     &                IERROR )

      CALL MPI_ISEND( CRNR, MM, MPI_DATATYPE, IPROC,
     &                MSGNUM, MPI_COMM_WORLD,RECV_REQUESTS(2),
     &                IERROR )

c.... Complete the non-blocking send and receive before proceeding.
      CALL MPI_WAITALL( 2, RECV_REQUESTS, RECV_STATUS,IERROR)

c...  Copy into correct location.
      CALL SCOPY( MM, COML(1,1), 1, CRNR, 1 )

      endif

      end

c -----------------------------------------------------------------------------
*dk pcom1to4
      subroutine pcom1to4(a4,a1,lva1,kr,nj)
 
c...  This routine copies array a1 that contains but one point in the
c...  lateral grid per active processor to array a4 that contains four
c...  points in the lateral grid per active processor (but on only
c...  one-fourth the number of active processors). The parameter lva1
c...  specifies the grid level of the array a1.
 
      include 'size.h'
      include 'pcom.h'
      real a1(0:1,2,nd*nj*(kr+1)), a4(0:2,3,nd*nj*(kr+1)),
     &     a(0:1,2,189*nd*9*4)
      integer sender
 
c...  Compute the message length.
      mm = nd*nj*(kr+1)
 
c...  Increment the message number.
      msgnum = msgnum + 1
 
c...  Compute the number of processors npedg1
c...  along a diamond edge for array a1.
      npedg1 = 2**lva1
 
c...  Compute the destination processor number idest.
      idest  = (mynum/(2*npedg1))*npedg1/2 + mod(mynum, npedg1)/2
 
c...  If idest is the local processor, we need to invoke special
c...  logic in the receiving loop to avoid sending messages locally. 
      if(idest.eq.mynum)then
         iend=3
      else 
         iend=4
      endif
     
      call nulvec(a4, 9*mm)

c...  If idest is not the local processor, then send array a1 to idest. 
      IF(idest.ne.mynum) CALL MPI_ISEND( A1, 4*MM, MPI_DATATYPE, 
     &                   IDEST, MSGNUM, MPI_COMM_WORLD, 
     &                   RECV_REQUESTS(2), IERROR )

      IF(mynum .lt. npedg1**2*10/(4*nd)) THEN
 
c...     Processor mynum is an active processor for array a4.  Receive
c...     three contributions from the three non-local processors if you
c...     you are processor zero. Otherwise, receive four contributions
c...     from the four non-local processors, if you are processor one. 
c...     Proceed and load array a4.
 
         do iproc=1,iend
 
            CALL MPI_IRECV( A, 4*MM, MPI_DATATYPE, MPI_ANY_SOURCE,
     &                      MSGNUM, MPI_COMM_WORLD, RECV_REQUESTS( 1 ),
     &                      IERROR )
 
c...        Complete the nonblocking receives before processing.
            CALL MPI_WAITANY( 1, RECV_REQUESTS (1), SENDER,
     &                        RECV_STATUS, IERROR )
 
            isrc_proc = RECV_STATUS( MPI_SOURCE, SENDER )

c...        Compute the (i1,i2) coordinates for this contribution.
            i1 = mod(isrc_proc, 2) + 1
            i2 = mod(isrc_proc/npedg1, 2) + 1
 
c...        Load array a4.
            if(nj .ne. 189) then
               do ii=1,mm
                  a4(i1-1,i2  ,ii) = a(0,1,ii)
                  a4(i1  ,i2  ,ii) = a(1,1,ii)
                  a4(i1-1,i2+1,ii) = a(0,2,ii)
                  a4(i1  ,i2+1,ii) = a(1,2,ii)
               end do
            else
               do ii=1,mm
                  a4(i1-1,i2  ,ii) = a4(i1-1,i2  ,ii) + a(0,1,ii)
                  a4(i1  ,i2  ,ii) = a4(i1  ,i2  ,ii) + a(1,1,ii)
                  a4(i1-1,i2+1,ii) = a4(i1-1,i2+1,ii) + a(0,2,ii)
                  a4(i1  ,i2+1,ii) = a4(i1  ,i2+1,ii) + a(1,2,ii)
               end do
            endif
 
         end do
 
c...     On processor zero, the remaining local portion of a1 can be copied 
c...     locally into the buffer array "a", and is then loaded into array a4

         if(idest.eq.mynum)then
            call scopy(4*mm, a1, 1, a, 1)

c...        Compute the (i1,i2) coordinates for this local contribution.
            i1 = mod(mynum, 2) + 1
            i2 = mod(mynum/npedg1, 2) + 1

c...        Load array a4.
            if(nj .ne. 189) then
               do ii=1,mm
                  a4(i1-1,i2  ,ii) = a(0,1,ii)
                  a4(i1  ,i2  ,ii) = a(1,1,ii)
                  a4(i1-1,i2+1,ii) = a(0,2,ii)
                  a4(i1  ,i2+1,ii) = a(1,2,ii)
               end do
            else
               do ii=1,mm
                  a4(i1-1,i2  ,ii) = a4(i1-1,i2  ,ii) + a(0,1,ii)
                  a4(i1  ,i2  ,ii) = a4(i1  ,i2  ,ii) + a(1,1,ii)
                  a4(i1-1,i2+1,ii) = a4(i1-1,i2+1,ii) + a(0,2,ii)
                  a4(i1  ,i2+1,ii) = a4(i1  ,i2+1,ii) + a(1,2,ii)
               end do
            endif
         endif

      ENDIF

c...  Complete non-blocking sends before exiting.
      CALL MPI_WAITANY( 1, RECV_REQUESTS(2), SENDER,
     &                        RECV_STATUS, IERROR )

      end

c -----------------------------------------------------------------------------
*dk pcom2
      subroutine pcom2(coml,comr,crnr,buff,nj,kr,kt,ighost)
 
c...  This routine broadcasts and receives boundary information
c...  to and from the neighboring subdomains.
 
      include 'size.h'
      include 'pcom.h'
      real  coml((kt+1)*nd*nj*(kr+1),2)
      real  comr((kt+1)*nd*nj*(kr+1),2)
      real  crnr(nd*nj*(kr + 1))
      real  buff((kt+1)*nd*nj*(kr+1),2)
      integer sender, ibegin, iend, ipending
      common /nbr1/ ileft(4,0:6), iright(4,0:6)

      ibegin   = 1
      iend     = 2
      ipending = 0

c...  Increment the message number.
      msgnum = msgnum + 1
 
c...  Compute message length along subdomain edges.
      mm = (kt + 1)*nd*nj*(kr + 1)
 
c...  Broadcast row/column data in coml to two subdomains to the left,
c...  and receive similar information from two subdomains to the right.
c...  To avoid selfcommunication check, if neighoring processor number
c...  equals mynum. In this case don't send array comr, but rather copy it.

      if( iright(1,lvproc) .eq. mynum) then
      ibegin=iend
      call scopy(mm,comr(1,1),1,coml(1,1),1)
      endif
      if(iright(2,lvproc).eq.mynum) then
      iend=ibegin
      call scopy(mm,comr(1,2),1,coml(1,2),1)
      endif

      do ib=ibegin,iend
         CALL MPI_IRECV( BUFF(1,IB), MM, MPI_DATATYPE,
     &                   MPI_ANY_SOURCE, MSGNUM, MPI_COMM_WORLD,
     &                   RECV_REQUESTS( 1 ), IERROR )

c...     Count number of pending send operations
         ipending = ipending + 1

         CALL MPI_ISEND( COMR( 1, IB ), MM, MPI_DATATYPE,
     &                   IRIGHT( IB, LVPROC ), MSGNUM, MPI_COMM_WORLD,
     &                   SEND_REQUESTS( ipending ), IERROR )

c...     Complete the nonblocking receives before processing.
         CALL MPI_WAITANY( 1, RECV_REQUESTS, SENDER,
     &                     RECV_STATUS, IERROR )

c...     Save this information in array comr.
c...     Determine the "ileft" value of the isrc_proc subdomain.
         isrc_proc = RECV_STATUS( MPI_SOURCE, SENDER )
 
         do i=1,2
            if(ileft(i,lvproc) .eq. isrc_proc)
     &         call scopy(mm, buff(1,ib), 1, coml(1,i), 1)
         end do
      end do

c.... Complete the non-blocking sends before proceeding.
      CALL MPI_WAITALL( iend - ibegin + 1, SEND_REQUESTS,
     &     SEND_STATUS, IERROR )

c...  Increment the message number.
      msgnum = msgnum + 1

c...  Compute message length for subdomain corner.
      mm = nd*nj*(kr + 1)
 
c...  Send right corner data to appropriate right subdomain.

      iproc = iright(3,lvproc)
      if(ighost .eq. 1) iproc = iright(4,lvproc)
 
c...  Avoid self-communication
      if(iproc.ne.mynum)then

c...  Receive left corner data to appropriate left subdomain.
      CALL MPI_IRECV( COMR(1,1), MM, MPI_DATATYPE, MPI_ANY_SOURCE,
     &                MSGNUM, MPI_COMM_WORLD, RECV_REQUESTS( 1 ),
     &                IERROR )

      CALL MPI_ISEND( CRNR, MM, MPI_DATATYPE, IPROC, MSGNUM,
     &                MPI_COMM_WORLD, RECV_REQUESTS(2), IERROR )

c.... Complete the non-blocking send and receive before proceeding.
      CALL MPI_WAITALL( 2, RECV_REQUESTS, RECV_STATUS, IERROR )

c...  Copy into correct location.
      CALL SCOPY( MM, COMR(1,1), 1, CRNR, 1 )
      endif
 
      end

c -----------------------------------------------------------------------------
*dk pcom4to1
      subroutine pcom4to1(a1,a4,lva1,kr,nj)
 
c...  This routine copies array a4 that contains four points in the
c...  lateral grid per active processor to array a1 that contains but
c...  one point in the lateral grid per active processor (but on four
c...  times the number of active processors).  The parameter lva1
c...  specifies the grid level of the array a1. A special case arises
c...  for processor zero, where one needs to avoid selfcommunication.
c...  Special logic is invoked for proc-zero to avoid this.

      include 'size.h'
      include 'pcom.h'
      real a1(0:1,2,nd*nj*(kr+1)), a4(0:2,3,nd*nj*(kr+1))
      integer sender, iflag, j1, j2
 
c...  Compute the message length.
      mm = nd*nj*(kr+1)
 
c...  Increment the message number.
      msgnum = msgnum + 1
 
c...  Compute the number of processors npedg1
c...  along a diamond edge for array a1.
 
      npedg1 = 2**lva1

      if(mynum .lt. npedg1**2*10/(4*nd)) then

c...  The receiving flag is set to zero for all sending processors, 
c...  in order to avoid communication to themselves.
 
c...     Processor mynum is an active processor for array a4.  Send
c...     contributions from array a4 to four other processors.
 
         do i2=1,2
            do i1=1,2
 
c...           Copy the appropriate elements of a4 into array a1.
 
               do ii=1,mm
                  a1(0,1,ii) = a4(i1-1,i2  ,ii)
                  a1(1,1,ii) = a4(i1  ,i2  ,ii)
                  a1(0,2,ii) = a4(i1-1,i2+1,ii)
                  a1(1,2,ii) = a4(i1  ,i2+1,ii)
               end do
 
c...           Compute the destination processor number idest.
 
               idest = (2*mynum/npedg1)*2*npedg1 + (i2 - 1)*npedg1
     &               + mod(mynum, npedg1/2)*2 + i1 - 1
              
c...           Send to idest only, If idest is not the local processor.
               if(idest.eq.mynum)then
                  j1=i1
                  j2=i2
               else
                   CALL MPI_SEND(A1, 4*MM, MPI_DATATYPE, IDEST,
     &                           MSGNUM, MPI_COMM_WORLD, IERROR )
               endif
 
            end do
         end do
 
      endif

c...  Receive array a1 only for non-local processors. Complete the 
c...  loading of array 'a1' on processor zero through a local copy.

      if(mynum.ne.0) then
         CALL MPI_IRECV( A1, 4*MM, MPI_DATATYPE, MPI_ANY_SOURCE,
     &                   MSGNUM, MPI_COMM_WORLD, RECV_REQUESTS(1),
     &                   IERROR )

      else

c...     Copy the appropriate elements of a4 into array a1.
         do ii=1,mm
            a1(0,1,ii) = a4(j1-1,j2  ,ii)
            a1(1,1,ii) = a4(j1  ,j2  ,ii)
            a1(0,2,ii) = a4(j1-1,j2+1,ii)
            a1(1,2,ii) = a4(j1  ,j2+1,ii)
         end do

      endif
 
c...  Complete the nonblocking receive before processing.
      CALL MPI_WAITANY( 1, RECV_REQUESTS(1), SENDER,
     &                  RECV_STATUS, IERROR )
 
      end

c -----------------------------------------------------------------------------
*dk pcomfromzero
      subroutine pcomfromzero(a,ap,nc)
 
c...  At the lowest (lv = 0) grid level, this routine takes a global
c...  array a on processor zero and distributes the appropriate portions
c...  to all the active processors as array ap.
 
      include 'size.h'
      include 'pcom.h'
      real a(4,10,nc), ap(4,nd,nc)
      integer sender
 
c...  Increment the message number.
      msgnum = msgnum + 1
 
c...  Compute the message length for the subarray.
      mm = 4*nd*nc
 
      if(mynum .eq. 0) then
 
c...     Send the appropriate portion of the global array to
c...     the other active processors.
         do iproc=1,10/nd-1
 
            do id=1,nd
               jd = id + iproc*nd
               do ii=1,4
                  do j=1,nc
                     ap(ii,id,j) = a(ii,jd,j)
                  end do
               end do
            end do
 
            CALL MPI_SEND( AP, MM, MPI_DATATYPE, IPROC,
     &                     MSGNUM, MPI_COMM_WORLD, IERROR )
 
         end do
 
c...     Now load the correct portion of the global arrary
c...     into the local array.
         do id=1,nd
            do ii=1,4
               do j=1,nc
                  ap(ii,id,j) = a(ii,id,j)
               end do
            end do
         end do
 
      else
 
c...     For all processes other than zero, receive the local
c...     array from process zero.
         CALL MPI_IRECV( AP, MM, MPI_DATATYPE, MPI_ANY_SOURCE,
     &                   MSGNUM, MPI_COMM_WORLD, RECV_REQUESTS( 1 ),
     &                   IERROR )
 
c...     Complete the nonblocking receive before processing.
         CALL MPI_WAITANY( 1, RECV_REQUESTS(1), SENDER,
     &                     RECV_STATUS, IERROR )
 
      end if
 
      end

c -----------------------------------------------------------------------------
*dk pcomptcl
      subroutine pcomptcl(ibndr,w,nsnd,nrcv)
 
c...  This routine performs interprocess communication to move
c...  particles across subdomain boundaries.
 
      include 'size.h'
      include 'pcom.h'
      parameter (np=(nt+1)**2*nd*5)
      real w(5,(nt+1)**2*5,4)
      common /prt1/ npart, p(np), q(np), idp(np), kpl(np), ap(np)
      common /nbr1/ ileft(4,0:6), iright(4,0:6)
      integer ibndr((nt+1)**2*5,4), nsnd(4), nrcv(4)
      integer sender
 
c...  Increment the message number.
      msgnum = msgnum + 1
 
      do ib=1,4
 
c...     Send particle data to the four neighboring subdomains.
 
         if(ib.eq.1) iproc = iright(1,lvproc)
         if(ib.eq.2) iproc =  ileft(2,lvproc)
         if(ib.eq.3) iproc =  ileft(1,lvproc)
         if(ib.eq.4) iproc = iright(2,lvproc)
 
c...     Load particle attributes into array w to be sent to a
c...     neighboring subdomain.
 
         do i=1,nsnd(ib)
            k = ibndr(i,ib)
            w(1,i,1) = p(k)
            w(2,i,1) = q(k)
            w(3,i,1) = idp(k)
            w(4,i,1) = kpl(k)
            w(5,i,1) = ap(k)
         end do
 
         CALL MPI_SEND( W( 1,1,1 ), 5*NSND( IB ),
     &                  MPI_DATATYPE, IPROC, MSGNUM,
     &                  MPI_COMM_WORLD, IERROR )
 
      end do
 
      do ib=1,4
 
c...     Receive particle data from the four neighboring subdomains.
 
         CALL MPI_PROBE( MPI_ANY_SOURCE, MSGNUM, MPI_COMM_WORLD,
     &                   MSG_STATUS, IERROR )
 
         CALL MPI_GET_COUNT( MSG_STATUS, MPI_DATATYPE, ICOUNT ,IERROR )
 
         iproc = MSG_STATUS( MPI_SOURCE )
 
         CALL MPI_IRECV( W( 1,1,IB ), ICOUNT, MPI_DATATYPE,
     &                   IPROC, MSGNUM, MPI_COMM_WORLD,
     &                   RECV_REQUESTS( 1 ), IERROR )
 
c...     Complete the nonblocking receives before processing.
         CALL MPI_WAITANY( 1, RECV_REQUESTS, SENDER, RECV_STATUS,
     &                     IERROR )
 
         nrcv(ib) = ICOUNT/5
 
      end do
 
      end

c -----------------------------------------------------------------------------
*dk pcomtozero
      subroutine pcomtozero(a,ap,nc)
 
c...  At the lowest (lv = 0) grid level, this routine loads a global
c...  array a on processor zero from contributions ap from all the
c...  other active processors.
 
      include 'size.h'
      include 'pcom.h'
      real a(4,10,nc), ap(4,nd,nc)
      integer sender
 
c...  Increment the message number.
      msgnum = msgnum + 1
 
c...  Compute the message length for the subarray.
      mm = 4*nd*nc
 
      if(mynum .ne. 0) then
 
c...     For all processes other than zero, broadcast the local
c...     array to process zero.
 
         CALL MPI_SEND( AP, MM, MPI_DATATYPE, 0, MSGNUM,
     &                  MPI_COMM_WORLD, IERROR )
 
      else
 
c...     On processor zero, first load the local array into
c...     the global array.
 
         do id=1,nd
            do ii=1,4
               do j=1,nc
                  a(ii,id,j) = ap(ii,id,j)
               end do
            end do
         end do
 
c...     Next receive the arrays from all the other processes
c...     and copy them into the global array.
 
         do iproc=1,10/nd-1
 
            CALL MPI_IRECV( AP , MM, MPI_DATATYPE, MPI_ANY_SOURCE,
     &                      MSGNUM, MPI_COMM_WORLD, RECV_REQUESTS( 1 ),
     &                      IERROR )
 
c...        Complete the nonblocking receives before processing.
            CALL MPI_WAITANY( 1, RECV_REQUESTS(1), SENDER, RECV_STATUS,
     &                        IERROR )
 
c...        Save this information in array comr.
c...        Determine the "iright" value of the isrc_proc subdomain.
            jproc = RECV_STATUS( MPI_SOURCE, SENDER )
 
            do id=1,nd
               jd = id + jproc*nd
               do ii=1,4
                  do j=1,nc
                     a(ii,jd,j) = ap(ii,id,j)
                  end do
               end do
            end do
 
         end do
 
      end if
 
      end

c -----------------------------------------------------------------------------
*dk pcomtozeronew
      subroutine pcomtozeronew(a,ap,kt,lt,kr)

c...  At the lowest (lv = 0) grid level, this routine loads a global
c...  array a on processor zero from contributions ap from all the
c...  other active processors.

      include 'size.h'
      include 'pcom.h'
      real a((kt+1)**2,10,kr+1), ap((lt+1)**2,nd,kr+1)
      integer sender

c...  Increment the message number.
      msgnum = msgnum + 1

c...  Compute the message length for the subarray.
      mm = (lt+1)**2*nd*(kr+1)

      if(mynum .ne. 0) then

c...     For all processes other than zero, broadcast the local
c...     array to process zero.

         CALL MPI_SEND( AP, MM, MPI_DATATYPE, 0, MSGNUM,
     &                  MPI_COMM_WORLD, IERROR )

      else

c...     On processor zero, first load the local array into
c...     the global array.

         call collectnew(a,ap,0,kt,lt,nd,kr+1)

c...     Next receive the arrays from all the other processes
c...     and copy them into the global array.

         do iproc=1,nproc-1

            CALL MPI_IRECV( AP , MM, MPI_DATATYPE, MPI_ANY_SOURCE,
     &                      MSGNUM, MPI_COMM_WORLD, RECV_REQUESTS( 1 ),
     &                      IERROR )

c...        Complete the nonblocking receives before processing.
            CALL MPI_WAITANY( 1, RECV_REQUESTS, SENDER, RECV_STATUS,
     &                        IERROR )

c...        Save this information in array comr.
c...        Determine the "iright" value of the isrc_proc subdomain.
            jproc = RECV_STATUS( MPI_SOURCE, SENDER )

            call collectnew(a,ap,jproc,kt,lt,nd,kr+1)

         end do

      end if

      end

c -----------------------------------------------------------------------------
*dk pexit
      subroutine pexit
 
      CALL MPI_FINALIZE( IERROR )
 
      end

c -----------------------------------------------------------------------------
*dk pinit
      subroutine pinit
 
c...  This routine initializes the MPI processes.

      implicit none

      include 'size.h'
      include 'pcom.h'
      integer ileft, iright, isize, lv, lvg
      common /nbr1/ ileft(4,0:6), iright(4,0:6)
 
c...  Set the message number.
 
      msgnum = 1
 
c...  Enroll this process in MPI.
 
      CALL MPI_INIT( IERROR )
      CALL MPI_COMM_RANK( MPI_COMM_WORLD, MYNUM, IERROR )
      CALL MPI_COMM_SIZE( MPI_COMM_WORLD, NPROC, IERROR )

c...  Check for consistency in number of processes
      isize = (mt/nt)**2
      IF ( nd .EQ. 5 ) isize = isize * 2
      IF ( nproc .NE. isize ) THEN
         IF ( mynum .EQ. 0 ) THEN
            WRITE(*,'(''TERRA: Sorry, got '',I5,'' processes, '''//
     &           '''but need '',I5)') nproc, isize
         ENDIF
         CALL MPI_BARRIER( MPI_COMM_WORLD, IERROR )
         CALL MPI_ABORT( MPI_COMM_WORLD, -1, IERROR )
      ENDIF

c...  Determine the process numbers of neigboring subdomains.
 
      lvg = 1.45*log(real(mt/nt))
 
      do lv=0,lvg
         if(lv.ge.1 .and. mynum.lt.2**(2*lv)*10/nd) then
            call neighborid(2**lv,ileft(1,lv),iright(1,lv))
         else
            ileft(1,lv)  = mynum
            ileft(2,lv)  = mynum
            ileft(3,lv)  = mynum
            ileft(4,lv)  = mynum
            iright(1,lv) = mynum
            iright(2,lv) = mynum
            iright(3,lv) = mynum
            iright(4,lv) = mynum
         end if
         if(lv.eq.0 .and. mynum.le.1) then
            ileft(2,lv)  = 1 - mynum
            ileft(3,lv)  = 1 - mynum
            ileft(4,lv)  = 1 - mynum
            iright(2,lv) = 1 - mynum
            iright(3,lv) = 1 - mynum
            iright(4,lv) = 1 - mynum
         end if
      end do
 
      end

c -----------------------------------------------------------------------------
*dk pmax
      subroutine pmax(smax)
 
c...  This routine obtains the global maximum of the scalar
c...  quantity smax from all of the nproc processes.
 
      include 'pcom.h'
 
      qmax = smax
 
      CALL MPI_ALLREDUCE( QMAX, SMAX, 1, MPI_DATATYPE,
     &                    MPI_MAX, MPI_COMM_WORLD, IERROR )
 
      end
*dk pmin
      subroutine pmin(smin)
 
c...  This routine obtains the global minimum of the scalar
c...  quantity smin from all of the nproc processes.
 
      include 'pcom.h'
 
      qmin = smin
 
      CALL MPI_ALLREDUCE( QMIN, SMIN, 1, MPI_DATATYPE,
     &                    MPI_MIN, MPI_COMM_WORLD, IERROR )
 
      end


c -----------------------------------------------------------------------------
*dk psum
      subroutine psum(u,n)
 
c...  This routine obtains the global sum for each of the n
c...  components of array u from all of the nproc processes.
 
      include 'size.h'
      include 'pcom.h'
      real u(n), v(300)
 
      kproc = mproc*10/nd
 
      if(kproc .eq. 1) return
 
      if(mynum .lt. kproc) then
         do i=1,n
            v(i) = u(i)
         end do
      else
         do i=1,n
            v(i) = 0.
         end do
      end if
 
      CALL MPI_ALLREDUCE( V, U, N, MPI_DATATYPE,
     &                    MPI_SUM, MPI_COMM_WORLD, IERROR )
 
      end

c -----------------------------------------------------------------------------
*dk psumlong
      subroutine psumlong(u,wk,n)
 
c...  This routine obtains the global sum for each of the n
c     components of array u from all of the nproc processes.
c     The difference to routine psum lies in the explicit
c     passing in of work array wk. All computations are done
c     on processor zero. Results are passed back to the rest.
 
      include 'size.h'
      include 'pcom.h'
      real u(n), wk(n)
 
      kproc = mproc*10/nd
 
      if(kproc .eq. 1) return
 
      if(mynum .lt. kproc) then
         do i=1,n
            wk(i) = u(i)
         end do
      else
         do i=1,n
            wk(i) = 0.
         end do
      end if
 
      CALL MPI_ALLREDUCE( WK, U, N, MPI_DATATYPE,
     &                    MPI_SUM, MPI_COMM_WORLD, IERROR )
 
      end

c -----------------------------------------------------------------------------
c Auxiliary routine for debugging
c -----------------------------------------------------------------------------

*dk SYNCPROC
      SUBROUTINE SYNCPROC( msg )
      include 'pcom.h'
      CHARACTER*(*) msg
      CALL MPI_BARRIER( MPI_COMM_WORLD, ierror )
      IF ( mynum .EQ. 0 ) THEN
         WRITE (6,'(A)'), msg
         CALL FLUSH(6)
      ENDIF
      END SUBROUTINE

c -----------------------------------------------------------------------------
