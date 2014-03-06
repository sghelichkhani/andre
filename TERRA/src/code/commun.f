c -----------------------------------------------------------------------------
*DK collect

c> This routine reads the subarray ap for process kproc into
c> the appropriate portion of the global array a.  Array a has
c> dimensions (mt+1,mt+1,10,nc).
      subroutine collect(a,ap,kproc,nc)

      include 'size.h'
      real a(mt+1,mt+1,10,nc), ap(nt+1,nt+1,nd,nc)

c...  Determine subdomain limits.

      mproc = (mt/nt)**2
      iproc = mod(kproc, mproc)
      id0   = (kproc/mproc)*nd
      rnm   = real(nt)/real(mt)
      j0    = iproc*rnm
      i0    = (real(iproc)*rnm - j0)/rnm
      ibeg  = i0*nt + 1
      jbeg  = j0*nt + 1
      iend  = ibeg + nt - 1
      jend  = jbeg + nt - 1
      if(iend .eq. mt) iend = iend + 1
      if(jend .eq. mt) jend = jend + 1

c...  Load process subarray ap into global array a.

      do ic=1,nc

         do id=1,nd

            jd = id + id0

            j = 0
            do jj=jbeg,jend
               j = j + 1

               i = 0
               do ii=ibeg,iend
                  i = i + 1

                  a(ii,jj,jd,ic) = ap(i,j,id,ic)

               end do

            end do

         end do

      end do

      end

c -----------------------------------------------------------------------------
*dk collectnew
 
c> This routine reads the subarray ap for process kproc into
c> the appropriate portion of the global array a.  Array a has
c> dimensions (mt+1,mt+1,10,nc).
      subroutine collectnew(a,ap,kproc,mt,nt,nd,nc)
 
      real a(mt+1,mt+1,10,nc), ap(nt+1,nt+1,nd,nc)
 
c...  Determine subdomain limits.
 
      mproc = (mt/nt)**2
      iproc = mod(kproc, mproc)
      id0   = (kproc/mproc)*nd
      rnm   = real(nt)/real(mt)
      j0    = iproc*rnm
      i0    = (real(iproc)*rnm - j0)/rnm
      ibeg  = i0*nt + 1
      jbeg  = j0*nt + 1
      iend  = ibeg + nt - 1
      jend  = jbeg + nt - 1
      if(iend .eq. mt) iend = iend + 1
      if(jend .eq. mt) jend = jend + 1
 
c...  Load process subarray ap into global array a.
 
      do ic=1,nc
 
         do id=1,nd
 
            jd = id + id0
 
            j = 0
            do jj=jbeg,jend
               j = j + 1
 
               i = 0
               do ii=ibeg,iend
                  i = i + 1
 
                  a(ii,jj,jd,ic) = ap(i,j,id,ic)
 
               end do
 
            end do
 
         end do
 
      end do
 
      end

c -----------------------------------------------------------------------------
*dk comm3s
c> This routine serves as a driver routine to perform boundary
c> communication among mesh subdomains. It sends and receives
c> boundary information to and from neighboring subdomains,
c> sums together the boundary values and sends and receives
c> appropriately the resulting boundary data.
      subroutine comm3s(u,kr,kt,nj)
 
      include 'size.h'
      include 'pcom.h'
      common /bufr/ edgr(2*(nt+1)*nd*3*(nr+1)),
     &              edgl(2*(nt+1)*nd*3*(nr+1)), crnr(nd*3*(nr+1)),
     &              buff(2*(nt+1)*nd*3*(nr+1))
      common /clck/ itmng, sec(50)
      if(itmng.eq.1) call mytime(tin)
 
c...  If calculation is invoking but a single process, call the
c...  serial communication routine.
 
      if(nproc.eq.1 .or. (lvproc.eq.0 .and. nd.eq.10)) then
 
         call edgadd0(u,nj,kr,kt,1)
 
         if(itmng.eq.1) call mytime(tout)
         if(itmng.eq.1) sec(14) = sec(14) + (tout - tin)
         return
 
      endif
 
      if(mynum .ge. mproc*10/nd) return
 
c...  Copy left subdomain boundary into the arrays edgl and c.
 
      call edgload(u,edgl,crnr,kr,kt,nj)
 
c...  Broadcast boundary information to the three processors
c...  associated with the three subdomains to the left and
c...  receive boundary information from processors associated
c...  with the three subdomains to the right.
 
      call pcom1(edgl,edgr,crnr,buff,nj,kr,kt,0)
 
c...  Add appropriate boundary values together.
 
      call edgsum(u,edgr,crnr,kr,kt,nj)
 
c...  Broadcast and receive results to and from the three
c...  neighboring subdomains.
 
      call pcom2(edgl,edgr,crnr,buff,nj,kr,kt,0)
 
c...  Load this information back into the full array u.
 
      call edgrepl(u,edgl,crnr,kr,kt,nj)
 
c...  Sum the values at the polar nodes in the five contiguous diamonds,
 
      if(mod(mynum, mproc) .eq. 0) call sumpole(u,nj,nd,kr,kt,1)
 
      if(itmng.eq.1) call mytime(tout)
      if(itmng.eq.1) sec(14) = sec(14) + (tout - tin)
      end

c -----------------------------------------------------------------------------
*dk edgadd0

c> This routine adds the values at the nodes in row one of each
c> diamond to the values at the corresponding nodes in column
c> one of the diamond to the right and loads the result back into
c> both positions.  It does a similar operation for column kt+1
c> of each diamond and row kt+1 of the diamond to the right.  It
c> also sums the values at the poles from each diamond and loads
c> the results back into the pole locations.
      subroutine edgadd0(u,nj,kr,kt,ip)
 
      real u(0:kt,kt+1,10,nj*(kr+1))
 
      k = kt + 1
 
      do id=1,10
 
         idlr = id + 5
         idur = mod(id, 5) + 1
         if(id .ge. 6) then
            idlr = idur
            idur = idur + 5
         endif
 
         do j=1,nj*(kr+1)
 
            do i=1,kt-1
               u(0,i+1,id  ,j) = u(0,i+1,id,j) + u( i,  1,idur,j)
               u(i,  k,id  ,j) = u(i,  k,id,j) + u(kt,k-i,idlr,j)
            end do
 
            do i=1,kt-1
               u( i,  1,idur,j) = u(0,i+1,id,j)
               u(kt,k-i,idlr,j) = u(i,  k,id,j)
            end do
 
            u( 0,k,id  ,j) = u(0,k,id,j) + u(kt,1,idur,j)
     &                                   + u(kt,k,idlr,j)
            u(kt,1,idur,j) = u(0,k,id,j)
            u(kt,k,idlr,j) = u(0,k,id,j)
 
         end do
 
      end do
 
      call sumpole(u,nj,10,kr,kt,ip)
 
      end

c -----------------------------------------------------------------------------
*dk edgload

c> This subroutine copies the left boundary values of the array u
c> into the boundary communication arrays edgl and crnr.
      subroutine edgload(u,edgl,crnr,kr,kt,nj)
 
      include 'size.h'
      include 'pcom.h'
      real u(kt+1,kt+1,nd*nj*(kr+1)), edgl(kt+1,nd*nj*(kr+1),2), crnr(*)
 
      k  = kt + 1
      npedg = 2**lvproc
 
      do kk=1,nd*nj*(kr+1)
         do ii=1,kt+1
            edgl(ii,kk,1) = u(ii,1,kk)
            edgl(ii,kk,2) = u(k,ii,kk)
         end do
      end do
 
      if(mod(mynum, npedg**2) .lt. npedg) then
 
         do kk=1,nd*nj*(kr+1)
            crnr(kk) = u(1, 1,kk)
         end do
 
      else
 
         do kk=1,nd*nj*(kr+1)
            crnr(kk) = u(k, 1,kk)
         end do
 
      endif
 
      end

c -----------------------------------------------------------------------------
*dk edgrepl

c> This routine stores the updated boundary information from the
c> the left subdomains back into the left boundary of the full
c> array u of the local subdomain.
      subroutine edgrepl(u,edgl,crnr,kr,kt,nj)
 
      include 'size.h'
      include 'pcom.h'
      real u(kt+1,kt+1,nd*nj*(kr+1)), edgl(kt+1,nd*nj*(kr+1),2), crnr(*)
 
      k     = kt + 1
      npedg = 2**lvproc
      iproc = mod(mynum, npedg**2)
 
      do kk=1,nd*nj*(kr+1)
         do ii=2,kt+1
            u(ii,1,kk) = edgl(ii,kk,1)
            u(k,ii,kk) = edgl(ii,kk,2)
         end do
      end do
 
      if(iproc.lt.npedg .and. iproc.gt.0) then
 
         do kk=1,nd*nj*(kr+1)
            u(1,1,kk) = crnr(kk)
         end do
 
      elseif(iproc .ne. 0) then
 
         do kk=1,nd*nj*(kr+1)
            u(1,1,kk) = edgl(1,kk,1)
         end do
 
      endif
 
      if(iproc .ge. npedg) then
 
         do kk=1,nd*nj*(kr+1)
            u(k,1,kk) = crnr(kk)
         end do
 
      endif
 
      end

c -----------------------------------------------------------------------------
*dk edgsum

c> This routine adds together boundary values of neighboring
c> subdomains to update boundary information on local subdomain.
      subroutine edgsum(u,edgr,crnr,kr,kt,nj)
 
      include 'size.h'
      include 'pcom.h'
      real u(kt+1,kt+1,nd,nj*(kr+1))
      real edgr(kt+1,nd,nj*(kr+1),2), crnr(nd,*)
 
      k     = kt + 1
      m     = kt + 2
      npedg = 2**lvproc
      iproc = mod(mynum, mproc)
 
c...  Loop over diamonds.
 
      do id=1,nd
 
      idur = mod(id, 5) + 1
 
      if(nd .eq. 10) then
         idlr = id + 5
         if(id .ge. 6) then
            idlr = idur
            idur = idur + 5
         end if
      elseif(nd .eq. 5) then
         idlr = id
         if(mynum .ge. mproc) idlr = idur
      end if
 
c...  Check to determine whether right subdomain boundaries lie
c...  within the interior of diamond id.
 
      inur = 0
      inlr = 0
 
      if(mod(mynum, npedg) .ne. 0) then
         inur = 1
         idur = id
      endif
 
      if(iproc .lt. mproc-npedg) then
         inlr = 1
         idlr = id
      endif
 
c...  FIRST:  Add values on upper and lower boundaries.
 
      do i=2,kt
         do jr=1,nj*(kr+1)
            u(1,i,id,jr)      = u(1,i,id,jr) + edgr(i,idur,jr,1)
            edgr(i,idur,jr,1) = u(1,i,id,jr)
         end do
      end do
 
      if(inlr .eq. 1) then
 
c...     If inside diamond, don't change numbering for lower subdomain.
 
         do i=2,kt
            do jr=1,nj*(kr+1)
               u(i,k,id,jr)      = u(i,k,id,jr) + edgr(i,idlr,jr,2)
               edgr(i,idlr,jr,2) = u(i,k,id,jr)
            end do
         end do
 
      else
 
         do i=2,kt
            do jr=1,nj*(kr+1)
               u(i,k,id,jr)        = u(i,k,id,jr) + edgr(m-i,idlr,jr,2)
               edgr(m-i,idlr,jr,2) = u(i,k,id,jr)
            end do
         end do
 
      endif
 
c...  SECOND:  Treat corners where three or four subdomains join.
 
      if(iproc .eq. mproc-npedg) then
 
c...     Treat right corner of diamond (only three subdomains join).
 
         do jr=1,nj*(kr+1)
            u(1,k,id,jr)    = u(1,k,id,jr) + edgr(k,idur,jr,1)
     &                                     + edgr(k,idlr,jr,2)
            edgr(k,idur,jr,1) = u(1,k,id,jr)
            edgr(k,idlr,jr,2) = u(1,k,id,jr)
         end do
 
      elseif(inur.eq.0 .and. inlr.eq.1) then
 
c...     Corner lies on upper right boundary of diamond.
 
         do jr=1,nj*(kr+1)
            u(1,k,id,jr)      = u(1,k,id,jr) + edgr(k,idur,jr,1)
     &                                       + edgr(1,idlr,jr,2)
     &                                       + crnr(idur,jr)
            edgr(k,idur,jr,1) = u(1,k,id,jr)
            edgr(1,idlr,jr,2) = u(1,k,id,jr)
            crnr(idur,jr)     = u(1,k,id,jr)
         end do
 
      elseif(inur.eq.1 .and. inlr.eq.1) then
 
c...     Corner lies inside diamond.
 
         do jr=1,nj*(kr+1)
            u(1,k,id,jr)    = u(1,k,id,jr) + edgr(k,id,jr,1)
     &                                     + edgr(1,id,jr,2)
     &                                     + crnr(id,jr)
            edgr(k,id,jr,1) = u(1,k,id,jr)
            edgr(1,id,jr,2) = u(1,k,id,jr)
            crnr(id,jr)     = u(1,k,id,jr)
         end do
 
      elseif(inur.eq.1 .and. inlr.eq.0) then
 
c...     Corner lies on lower right boundary of diamond.
 
         do jr=1,nj*(kr+1)
            u(1,k,id,jr)      = u(1,k,id,jr) + edgr(k,idur,jr,1)
     &                                       + edgr(k,idlr,jr,2)
     &                                       + crnr(idlr,jr)
            edgr(k,idur,jr,1) = u(1,k,id,jr)
            edgr(k,idlr,jr,2) = u(1,k,id,jr)
            crnr(idlr,jr)     = u(1,k,id,jr)
         end do
 
      endif
 
      end do
 
      end

c -----------------------------------------------------------------------------
*dk ghostadd

c> This routine adds the values stored in ghost node locations of
c> array s to those at corresponding nodes in subdomain interiors.
      subroutine ghostadd(s,sp,nc)
 
      include 'size.h'
      include 'pcom.h'
      real s(0:nt+1,0:nt+1,nd,nc), sp(0:2,0:2,2,nc)
      common /work/ wk((nt+1)**2*(nr+1)*2)
      common /spok/ ispoke(10)
 
      if(nt*mod(mynum, mproc) .eq. mt-nt) then
         do id=1,nd
            do ic=1,nc
               s(nt+1,1,id,ic) = s(nt+1,1,id,ic) + s(nt+1,0,id,ic)
               s(nt+1,0,id,ic) = 0.
            end do
         end do
      end if
 
      if(mod(mynum, mproc) .eq. 0) then
 
         do id=1,nd
            jp = 1 + id/6
            j1 = 1 + ispoke(mod(id-1,5)+1)
            j2 = 1 + ispoke(mod(id-1,5)+5)
            do ic=1,nc
               s(1,1,id,ic) = s(1,1,id,ic) + sp(j1,j2,jp,ic)
            end do
         end do
 
         do id=1,nd,5
            jd = 1 + (id/5)*5
            jp = 1 +  id/6
            do ic=1,nc
               s(0,1,jd  ,ic) = s(0,1,jd  ,ic) + s(0,1,jd+1,ic)
     &                        + s(0,1,jd+2,ic) + s(0,1,jd+3,ic)
     &                        + s(0,1,jd+4,ic) + sp(1,1,jp,ic)
               s(0,1,jd+1,ic) = 0.
               s(0,1,jd+2,ic) = 0.
               s(0,1,jd+3,ic) = 0.
               s(0,1,jd+4,ic) = 0.
            end do
         end do
 
      end if
 
      if(nproc .eq. 1) then
 
c...     The following is for calculation on a single processor.
 
         do id=1,10
 
            idur = mod(id,   5) + 1
            idul = mod(id+3, 5) + 1
            idlr = id   + 5
            idll = idul + 5
            if(id .ge. 6) then
               idlr = idur
               idur = idur + 5
               idul = idll
               idll = id   - 5
            end if
 
            do ic=1,nc
               do i=1,nt
                  s( 1,     i,idul,ic) = s( 1,     i,idul,ic)
     &                                 + s( i,     0,id  ,ic)
                  s( i,     0,id  ,ic) = 0.
               end do
               do i=1,nt
                  s( i,     1,idur,ic) = s( i,     1,idur,ic)
     &                                 + s( 0,   i+1,id  ,ic)
                  s( 0,   i+1,id  ,ic) = 0.
               end do
               do i=1,nt
                  s(nt+1-i,nt,idll,ic) = s(nt+1-i,nt,idll,ic)
     &                                 + s(nt+1,   i,id  ,ic)
                  s(nt+1,   i,id  ,ic) = 0.
               end do
               do i=1,nt
                  s(nt,nt+1-i,idlr,ic) = s(nt,nt+1-i,idlr,ic)
     &                                 + s( i,  nt+1,id  ,ic)
                  s( i,  nt+1,id  ,ic) = 0.
               end do
            end do
 
         end do
 
      else
 
c...     For calculation on multiple processors, call parallel
c...     communication routine.
 
         call pghostadd(s,nc,0)
 
      end if
 
      end

c -----------------------------------------------------------------------------
*dk pghostadd

c> This routine adds the values stored in ghost node locations of
c> array s to those at corresponding nodes in subdomain interiors.
      subroutine pghostadd(s,nc,np)
 
      include 'size.h'
      include 'pcom.h'
      real s(0:nt+1,0:nt+1,nd,nc)
      common /bufr/ edgr(2*(nt+1)*nd*3*(nr+1)),
     &              edgl(2*(nt+1)*nd*3*(nr+1)), crnr(nd*3*(nr+1)),
     &              buff(2*(nt+1)*nd*3*(nr+1))
 
c...  Load data from the ghost points on the right edge of array s
c...  into buffer arrays edgr and crnr.
 
      call gstloadright(s,edgr,crnr,nc)
 
c...  Broadcast this data to the three processors associated with the
c...  three subdomains to the right and receive similar data from
c...  processors associated with the three subdomains to the left.
 
      call pcom2(edgl,edgr,crnr,buff,nc,0,nt,1)
 
c...  Add this information at appropriate locations in array s.
 
      if(np .eq. 0) then
         call gstaddleft(s,edgl,crnr,nc)
      else
         call ptclgstaddleft(s,edgl,crnr,np)
      end if
 
c...  Load data from the ghost points on the left edge of array s
c...  into buffer arrays edgl and crnr.
 
      call gstloadleft(s,edgl,crnr,nc)
 
c...  Broadcast this data to the three processors associated with the
c...  three subdomains to the left and receive similar data from
c...  processors associated with the three subdomains to the right.
 
      call pcom1(edgl,edgr,crnr,buff,nc,0,nt,1)
 
c...  Add this information at appropriate locations in array s.
 
      if(np .eq. 0) then
         call gstaddright(s,edgr,crnr,nc)
      else
         call ptclgstaddright(s,edgr,crnr,np)
      end if
 
      end

c -----------------------------------------------------------------------------
*dk gstloadright

c> This subroutine loads data from the ghost points on the right
c> edge of array s into buffer arrays edgr and crnr.
      subroutine gstloadright(s,edgr,crnr,nc)
 
      include 'size.h'
      include 'pcom.h'
      real s(0:nt+1,0:nt+1,nd*nc)
      real edgr(nt+1,nd*nc,2), crnr(*)
 
      npedg = 2**lvproc
 
      do kk=1,nd*nc
         do i=1,nt+1
            edgr(i,kk,1) = s(0,   i,kk)
            edgr(i,kk,2) = s(i,nt+1,kk)
         end do
      end do
 
      if(mod(mynum, npedg) .eq. 0) then
 
c...     Subdomain lies along upper right diamond boundary.
 
         do kk=1,nd*nc
            crnr(kk) = s(0,1,kk)
         end do
 
      else
 
         do kk=1,nd*nc
            crnr(kk) = s(0,nt+1,kk)
         end do
 
      end if
 
      end

c -----------------------------------------------------------------------------
*dk gstloadleft

c> This subroutine loads data from the ghost points on the left
c> edge of array s into buffer arrays edgl and crnr.
      subroutine gstloadleft(s,edgl,crnr,nc)
 
      include 'size.h'
      real s(0:nt+1,0:nt+1,nd*nc), edgl(nt+1,nd*nc,2), crnr(*)
 
      do kk=1,nd*nc
         do i=1,nt
            edgl(i,kk,1) = s(   i,0,kk)
            edgl(i,kk,2) = s(nt+1,i,kk)
         end do
      end do
 
      do kk=1,nd*nc
         crnr(kk) = s(nt+1,0,kk)
      end do
 
      end

c -----------------------------------------------------------------------------
*dk gstaddleft

c> This routine adds the ghost information received from the left
c> subdomains to the leftmost non-ghost row and column of array s.
      subroutine gstaddleft(s,edgl,crnr,nc)
 
      include 'size.h'
      include 'pcom.h'
      real s(0:nt+1,0:nt+1,nd,nc)
      real edgl(nt+1,nd,nc,2), crnr(nd,*)
 
      npedg = 2**lvproc
      iproc = mod(mynum, npedg**2)
 
c...  Loop over diamonds.
 
      do id=1,nd
 
         idul = mod(id+3, 5) + 1
 
         if(nd .eq. 10) then
            idll = idul + 5
            if(id .ge. 6) then
               idll = id - 5
               idul = idul + 5
            end if
         elseif(nd .eq. 5) then
            idll = idul
            if(mynum .ge. mproc) idll = id
         end if
 
c...     Add ghost contributions to column 1 from upper left subdomain.
 
         if(iproc .lt. npedg) then
 
c...        Subdomain lies along upper left diamond boundary.
 
            do ic=1,nc
               do i=1,nt
                  s(i,1,id,ic) = s(i,1,id,ic) + edgl(i+1,idul,ic,1)
               end do
            end do
 
         else
 
            do ic=1,nc
               do i=1,nt
                  s(i,1,id,ic) = s(i,1,id,ic) + edgl(i,id,ic,1)
               end do
            end do
 
         end if
 
c...     Add ghost contribution to row nt from lower left subdomain.
 
         if(mod(iproc+1, npedg) .eq. 0) then
 
c...        Subdomain lies along lower left diamond boundary.
 
            do ic=1,nc
               do i=1,nt
                  s(nt,nt+1-i,id,ic) = s(nt,nt+1-i,id,ic)
     &                               + edgl(i,idll,ic,2)
               end do
            end do
 
         else
 
            do ic=1,nc
               do i=1,nt
                  s(nt,i,id,ic) = s(nt,i,id,ic) + edgl(i,id,ic,2)
               end do
            end do
 
         end if
 
c...     Add ghost contributions at corner points.
 
         if(iproc .lt. npedg-1) then
 
c...        Subdomain lies along upper left diamond boundary.
 
            do ic=1,nc
               s(nt,1,id,ic) = s(nt,1,id,ic) + crnr(idul,ic)
            end do
 
         elseif(mod(iproc+1, npedg).eq.0 .and. iproc.ge.npedg) then
 
c...        Subdomain lies along lower left diamond boundary.
 
            do ic=1,nc
               s(nt,1,id,ic) = s(nt,1,id,ic) + crnr(idll,ic)
            end do
 
         elseif(iproc .ne. npedg-1) then
 
            do ic=1,nc
               s(nt,1,id,ic) = s(nt,1,id,ic) + crnr(id,ic)
            end do
 
         end if
 
      end do
 
      end

c -----------------------------------------------------------------------------
*dk gstaddright

c> This routine adds the ghost information received from the right
c> subdomains to the rightmost non-ghost row and column of array s.
      subroutine gstaddright(s,edgr,crnr,nc)
 
      include 'size.h'
      include 'pcom.h'
      real s(0:nt+1,0:nt+1,nd,nc)
      real edgr(nt+1,nd,nc,2), crnr(nd,*)
 
      npedg = 2**lvproc
      iproc = mod(mynum, mproc)
 
c...  Loop over diamonds.
 
      do id=1,nd
 
         idur = mod(id, 5) + 1
 
         if(nd .eq. 10) then
            idlr = id + 5
            if(id .ge. 6) then
               idlr = idur
               idur = idur + 5
            end if
         elseif(nd .eq. 5) then
            idlr = id
            if(mynum .ge. mproc) idlr = idur
         end if
 
         if(mod(mynum, npedg) .ne. 0) idur = id
 
c...     Add values on upper right side of subdomain.
 
         do i=1,nt
            do ic=1,nc
               s(1,i,id,ic) = s(1,i,id,ic) + edgr(i,idur,ic,1)
            end do
         end do
 
c...     Add values on lower right side of subdomain.
 
         if(iproc .lt. mproc-npedg) then
 
c...        Subdomain edge lies inside diamond.
 
            do i=1,nt
               do ic=1,nc
                  s(i,nt,id,ic) = s(i,nt,id,ic) + edgr(i,id,ic,2)
               end do
            end do
 
         else
 
            do i=1,nt
               do ic=1,nc
                  s(i,nt,id,ic) = s(i,nt,id,ic) + edgr(nt+1-i,idlr,ic,2)
               end do
            end do
 
         endif
 
c...     Treat corners where four subdomains join.
 
         if(mod(iproc, npedg).eq.0 .and. iproc.ne.0) then
 
c...        Corner lies on upper right boundary of diamond.
 
            do ic=1,nc
               s(1,1,id,ic) = s(1,1,id,ic) + crnr(idur,ic)
            end do
 
         elseif(mod(iproc, npedg).ne.0 .and. iproc.lt.mproc-npedg) then
 
c...        Corner lies inside diamond.
 
            do ic=1,nc
               s(1,nt,id,ic) = s(1,nt,id,ic) + crnr(id,ic)
            end do
 
         elseif(iproc.gt.mproc-npedg) then
 
c...        Corner lies on lower right boundary of diamond.
 
            do ic=1,nc
               s(1,nt,id,ic) = s(1,nt,id,ic) + crnr(idlr,ic)
            end do
 
         endif
 
      end do
 
      end

c -----------------------------------------------------------------------------
*dk ptclgstaddleft

c> This routine adds the ghost information received from the left
c> subdomains to the leftmost non-ghost row and column of array a.
      subroutine ptclgstaddleft(a,edgl,crnr,np)
 
      include 'size.h'
      include 'pcom.h'
      real a(0:nt+1,0:nt+1,nd,4,np)
      real edgl(nt+1,nd,4,np,2), crnr(nd,4,np)
 
      npedg = 2**lvproc
      iproc = mod(mynum, npedg**2)
 
c...  Loop over diamonds.
 
      do id=1,nd
 
         idul = mod(id+3, 5) + 1
 
         if(nd .eq. 10) then
            idll = idul + 5
            if(id .ge. 6) then
               idll = id - 5
               idul = idul + 5
            end if
         elseif(nd .eq. 5) then
            idll = idul
            if(mynum .ge. mproc) idll = id
         end if
 
c...     Add ghost contributions to column 1 from upper left subdomain.
 
         if(iproc .lt. npedg) then
 
c...        Subdomain lies along upper left diamond boundary.
 
            do ip=1,np
               do ia=1,4
                  do i=1,nt
                     a(i,1,id,ia,ip) = edgl(i+1,idul,ia,ip,1)
                  end do
               end do
            end do
 
         else
 
            do ip=1,np
               do ia=1,4
                  do i=1,nt
                     a(i,1,id,ia,ip) = edgl(i,id,ia,ip,1)
                  end do
               end do
            end do
 
         end if
 
c...     Add ghost contribution to row nt from lower left subdomain.
 
         na = 0
         do jp=1,4
            if(a(nt,1,id,1,jp) .lt. 1.) na = na + 1
         end do
 
         if(mod(iproc+1, npedg) .eq. 0) then
 
c...        Subdomain lies along lower left diamond boundary.
 
            do ip=1,np
               do ia=1,4
                  do i=1,nt-1
                     a(nt,nt+1-i,id,ia,ip) = edgl(i,idll,ia,ip,2)
                  end do
               end do
               if(edgl(nt,idll,1,ip,2) .lt. 1.) then
                  na = min(np, na + 1)
                  do ia=1,4
                     a(nt,1,id,ia,na) = edgl(nt,idll,ia,ip,2)
                  end do
               end if
            end do
 
         else
 
            do ip=1,np
               do ia=1,4
                  do i=2,nt
                     a(nt,i,id,ia,ip) = edgl(i,id,ia,ip,2)
                  end do
               end do
               if(edgl(1,id,1,ip,2) .lt. 1.) then
                  na = min(np, na + 1)
                  do ia=1,4
                     a(nt,1,id,ia,na) = edgl(1,id,ia,ip,2)
                  end do
               end if
            end do
 
         end if
 
c...     Add ghost contributions at corner points.
 
         idcr = id
         if(iproc .lt. npedg-1) idcr = idul
         if(mod(iproc+1, npedg)  .eq. 0) idcr = idll
 
         if(iproc .ne. npedg-1) then
 
            do ip=1,np
               if(crnr(idcr,1,ip) .lt. 1.) then
                  na = min(np, na + 1)
                  do ia=1,4
                     a(nt,1,id,ia,na) = crnr(idcr,ia,ip)
                  end do
               end if
            end do
 
         end if
 
      end do
 
      end

c -----------------------------------------------------------------------------
*dk ptclgstaddright

c> This routine adds the ghost information received from the right
c> subdomains to the rightmost non-ghost row and column of array a.
      subroutine ptclgstaddright(a,edgr,crnr,np)
 
      include 'size.h'
      include 'pcom.h'
      real a(0:nt+1,0:nt+1,nd,4,np)
      real edgr(nt+1,nd,4,np,2), crnr(nd,4,np)
 
      npedg = 2**lvproc
      iproc = mod(mynum, mproc)
 
c...  Loop over diamonds.
 
      do id=1,nd
 
         idur = mod(id, 5) + 1
 
         if(nd .eq. 10) then
            idlr = id + 5
            if(id .ge. 6) then
               idlr = idur
               idur = idur + 5
            end if
         elseif(nd .eq. 5) then
            idlr = id
            if(mynum .ge. mproc) idlr = idur
         end if
 
         if(mod(mynum, npedg) .ne. 0) idur = id
 
         n1 = 0
         n2 = 0
         n3 = 0
         do jp=1,4
            if(a( 1, 1,id,1,jp) .lt. 1.) n1 = n1 + 1
            if(a( 1,nt,id,1,jp) .lt. 1.) n2 = n2 + 1
            if(a(nt,nt,id,1,jp) .lt. 1.) n3 = n3 + 1
         end do
 
c...     Add values on upper right side of subdomain.
 
         do ip=1,np
            do ia=1,4
               do i=2,nt
                  a(1,i,id,ia,ip) = edgr(i,idur,ia,ip,1)
               end do
            end do
            if(edgr(1,idur,1,ip,1) .lt. 1.) then
               n1 = min(np, n1 + 1)
               do ia=1,4
                  a(1,1,id,ia,n1) = edgr(1,idur,ia,ip,1)
               end do
            end if
         end do
 
c...     Add values on lower right side of subdomain.
 
         if(iproc .lt. mproc-npedg) then
 
c...        Subdomain edge lies inside diamond.
 
            do ip=1,np
               do ia=1,4
                  do i=2,nt-1
                     a(i,nt,id,ia,ip) = edgr(i,id,ia,ip,2)
                  end do
               end do
               if(edgr(1,id,1,ip,2) .lt. 1.) then
                  n2 = min(np, n2 + 1)
                  do ia=1,4
                     a(1,nt,id,ia,n2) = edgr(1,id,ia,ip,2)
                  end do
               end if
               if(edgr(nt,id,1,ip,2) .lt. 1.) then
                  n3 = min(np, n3 + 1)
                  do ia=1,4
                     a(nt,nt,id,ia,n3) = edgr(nt,id,ia,ip,2)
                  end do
               end if
            end do
 
         else
 
            do ip=1,np
               do ia=1,4
                  do i=2,nt-1
                     a(i,nt,id,ia,ip) = edgr(nt+1-i,idlr,ia,ip,2)
                  end do
               end do
               if(edgr(nt,idlr,1,ip,2) .lt. 1.) then
                  n2 = min(np, n2 + 1)
                  do ia=1,4
                     a(1,nt,id,ia,n2) = edgr(nt,idlr,ia,ip,2)
                  end do
               end if
               if(edgr(1,idlr,1,ip,2) .lt. 1.) then
                  n3 = min(np, n3 + 1)
                  do ia=1,4
                     a(nt,nt,id,ia,n3) = edgr(1,idlr,ia,ip,2)
                  end do
               end if
            end do
 
         endif
 
c...     Treat corners where four subdomains join.
 
         if(mod(iproc, npedg).eq.0 .and. iproc.ne.0) then
 
c...        Corner lies on upper right boundary of diamond.
 
            do ip=1,np
               if(crnr(idur,1,ip) .lt. 1.) then
                  n1 = min(np, n1 + 1)
                  do ia=1,4
                     a(1,1,id,ia,n1) = crnr(idur,ia,ip)
                  end do
               end if
            end do
 
         elseif(mod(iproc, npedg).ne.0 .and. iproc.lt.mproc-npedg) then
 
c...        Corner lies inside diamond.
 
            do ip=1,np
               if(crnr(id,1,ip) .lt. 1.) then
                  n2 = min(np, n2 + 1)
                  do ia=1,4
                     a(1,nt,id,ia,n2) = crnr(id,ia,ip)
                  end do
               end if
            end do
 
         elseif(iproc.gt.mproc-npedg) then
 
c...        Corner lies on lower right boundary of diamond.
 
            do ip=1,np
               if(crnr(idlr,1,ip) .lt. 1.) then
                  n2 = min(np, n2 + 1)
                  do ia=1,4
                     a(1,nt,id,ia,n2) = crnr(idlr,ia,ip)
                  end do
               end if
            end do
 
         endif
 
      end do
 
      end

c -----------------------------------------------------------------------------
*dk neighborid

c> This routine computes the process number of the three neighboring
c> subdomains to the left and the three neighboring subdomains to the
c> right of the subdomain associated with process mynum.

c> \param npedg is the number of subdomains along a diamond edge
      subroutine neighborid(npedg,ileft,iright)
 
      include 'size.h'
      include 'pcom.h'
      integer ileft(4), iright(4)
 
c...  kproc is the number of processors onto which each set of nd
c           diamonds is mapped.
c...  iproc is the local process number relative to the set of kproc
c...        processes.
c...  kown  is the increment in process number between mynum and iproc.
c...  koth  is the increment if process number for diamond located
c...        in the other hemisphere.
 
      kproc = npedg**2
      iproc = mod(mynum, kproc)
 
      if(nd .eq. 10) then
         kown = 0
         koth = 0
      elseif(nd .eq. 5) then
         kown = mynum - iproc
         koth = kproc - kown
      endif
 
c...  Treat case of subdomain boundaries within the diamond interior.
 
      iright(1) = mynum - 1
      iright(2) = mynum + npedg
      iright(3) = mynum + npedg - 1
      iright(4) = iright(3)
      ileft(1)  = mynum - npedg
      ileft(2)  = mynum + 1
      ileft(3)  = mynum - npedg + 1
      ileft(4)  = ileft(3)
 
c...  1) Special case of upper right boundary of diamond:
 
      if(mod(iproc,npedg) .eq. 0) then
         iright(1) = iproc/npedg + kown
         iright(2) = mynum     + npedg
         iright(3) = iright(1) + 1
         iright(4) = iright(1) - 1
      endif
 
c...  2) Special case of lower right boundary of diamond:
 
      if(iproc .gt. kproc-npedg-1) then
         iright(1) =  mynum - 1
         iright(2) = (kproc - iproc)*npedg - 1 + koth
         iright(3) = iright(2) + npedg
         iright(4) = iright(3)
      endif
 
c...  3) Special case of upper left boundary of diamond:
 
      if(iproc .lt. npedg) then
         ileft(1) = iproc*npedg + kown
         ileft(2) = mynum + 1
         ileft(3) = ileft(1) - npedg
         ileft(4) = ileft(1) + npedg
      endif
 
c...  4) Special case of lower left boundary of diamond:
 
      if(mod(iproc+1,npedg) .eq. 0) then
         ileft(1) = mynum - npedg
         ileft(2) = kproc - (iproc + 1)/npedg + koth
         ileft(3) = ileft(2) + 1
         ileft(4) = ileft(3)
      endif
 
c...  5) Special case of right corner node of diamond:
 
      if(iproc .eq. kproc-npedg) then
         iright(1) = iproc/npedg  + kown
         iright(2) = kproc - 1 + koth
         iright(3) = kown
         iright(4) = iproc/npedg - 1 + kown
      endif
 
c...  6) Special case of pole:
 
      if(iproc .eq. 0) then
         ileft(1)  = kown
         ileft(2)  = kown + 1
         ileft(3)  = kown + (npedg - 1)*npedg
         iright(4) = npedg - 1 + kown
      endif
 
c...  7) Special case of left corner node of diamond:
 
      if(iproc .eq. npedg-1) then
         ileft(1) = iproc*npedg + kown
         ileft(2) = kproc - (iproc + 1)/npedg + koth
         ileft(3) = ileft(1) - npedg
         ileft(4) = kown
      endif
 
      end

c -----------------------------------------------------------------------------
*dk subarray
 
c> This routine reads the appropriate portion of the global array
c> a into the sub-array ap needed for process mynum.  Array a has
c> dimensions (kt+1,kt+1,kd,nc) and ibctype specifies the type of
c> boundary conditions to be applied along subdomain edges.  For
c> ibctype = 0, none of the edge elements are set to zero; for
c> ibctype = 1, all of the edge elements are set to zero; and for
c> ibctype = 2, array ap is treated as a 7-point stencil and edge
c> node components outside the subdomain are set to zero.
      subroutine subarray(a,ap,ibctype,kd,kdp,kt,ktp,nc)
 
      include 'size.h'
      include 'pcom.h'
      real a(kt+1,kt+1,kd,nc), ap(ktp+1,ktp+1,kdp,nc)
 
c...  Determine subdomain limits.
 
      npedg = 2**lvproc
      iproc = mod(mynum, mproc)
      rnm   = 1./real(npedg)
      j0    = iproc*rnm
      i0    = (real(iproc)*rnm - j0)/rnm
 
      ibeg  = i0*ktp + 1
      jbeg  = j0*ktp + 1
      iend  = ibeg + ktp
      jend  = jbeg + ktp
 
c...  Load sub-array.
 
      do ic=1,nc
         do id=1,kdp
 
            jd = id
            if(kdp.eq.5 .and. mynum.ge.mproc) jd = id + 5
 
            j = 0
            do jj=jbeg,jend
               j = j + 1
 
               i = 0
               do ii=ibeg,iend
                  i = i + 1
 
                  ap(i,j,id,ic) = a(ii,jj,jd,ic)
 
               end do
 
            end do
 
         end do
      end do
 
c...  Apply boundary conditions.
 
      if(ibctype .eq. 0) return
 
      if(ibctype .eq. 1) then
 
         k  = ktp + 1
         i1 = 1
         if(iproc .eq. 0) i1 = 2
 
         do ic=1,nc
            do id=1,kdp
               do i=i1,ktp+1
                  ap(1,i,id,ic) = 0.
                  ap(i,k,id,ic) = 0.
               end do
            end do
         end do
 
      elseif(ibctype .eq. 2) then
 
         call subarraybc(ap,nc/7,ktp)
 
      endif
 
      end

c -----------------------------------------------------------------------------
*dk subarraybc

c> This routine sets to zero the appropriate stencil components
c> along the sub-domain boundaries for the operator ap.  The flag
c> idiamond, when set to one, causes redundant components on the
c> diamond edges to be set to zero.
      subroutine subarraybc(ap,nc,kt)
 
      include 'size.h'
      include 'pcom.h'
      real ap(kt+1,kt+1,7,nc)
 
      k     = kt + 1
      npedg = 2**lvproc
      iproc = mod(mynum, mproc)
      i0    = 1
      if(iproc .eq. 0) i0 = 2
 
      do j=1,nc
 
         do i=i0,kt+1
 
c...        Treat upper right edge.
 
            ap(1,i,4,j) = 0.
            ap(1,i,5,j) = 0.
 
            if(mod(iproc, npedg) .ne. 0) then
               ap(1,i,1,j) = 0.
               ap(1,i,3,j) = 0.
               ap(1,i,6,j) = 0.
            endif
 
c...        Treat upper left edge.
 
            ap(i,1,6,j) = 0.
            ap(i,1,7,j) = 0.
 
         end do
 
         do i=1,kt+1
 
c...        Treat lower left edge.
 
            ap(k,i,2,j) = 0.
            ap(k,i,7,j) = 0.
 
c...        Treat lower right edge.
 
            ap(i,k,3,j) = 0.
            ap(i,k,4,j) = 0.
 
            if(iproc .le. mproc-1-npedg) then
               ap(i,k,1,j) = 0.
               ap(i,k,2,j) = 0.
               ap(i,k,5,j) = 0.
            endif
 
         end do
 
      end do
 
      end

c -----------------------------------------------------------------------------
*dk sumpole

c>    This routine sums the values in array u at the polar nodes
c>    in the five diamonds and stores these sums back into the polar
c>    node locations when the flag ip = 1 and weights the sum by
c>    0.20 when ip has some other value.
      subroutine sumpole(u,nj,kd,kr,kt,ip)
 
      real u( 0:kt, kt+1, kd, nj*(kr+1) )
 
      aa = 0.2
      if(ip .eq. 1) aa = 1.
 
      do id=1,kd,5
         do ii=1,nj*(kr+1)
            u(0,1,id  ,ii) = (u(0,1,id  ,ii) + u(0,1,id+1,ii)
     &                     +  u(0,1,id+2,ii) + u(0,1,id+3,ii)
     &                     +  u(0,1,id+4,ii))*aa
            u(0,1,id+1,ii) =  u(0,1,id  ,ii)
            u(0,1,id+2,ii) =  u(0,1,id  ,ii)
            u(0,1,id+3,ii) =  u(0,1,id  ,ii)
            u(0,1,id+4,ii) =  u(0,1,id  ,ii)
         end do
      end do
 
      end

c -----------------------------------------------------------------------------
