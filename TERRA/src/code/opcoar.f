*dk coarseopr
      subroutine coarseopr(b,a,d,rw,hw,w,kd,nrf,ntf)
 
c...  This routine assembles coarse grid forward operators.
 
      include 'size.h'
      real a(ntf+1  ,ntf+1  ,kd,9,7,3, nrf+1)
      real b(ntf/2+1,ntf/2+1,ndo,9,7,3,(nrf+1)/2+1)
      real w(ntf/2+1,ntf/2+1,ndo,4,7)
      common /ofst/ j1n(7), j2n(7), md(7)
 
      if(nrf.eq.1) then
         call zerolvopr(a,b,hw)
         return
      endif
 
      nrc = (nrf+1)/2
      ntc =  ntf/2
 
      call nulvec(b,189*ndo*(ntc+1)**2*(nrc+1))
 
      if(ntf .eq. nt) call fineoprlayer(a,1)
 
      do 60 irc=1,nrc
 
      jrc = irc + 1
      irf = irc + irc - 1
      jrf = irf + 1
      krf = jrf + 1
 
      call weights(w,rw,hw,irc,ntc)
 
      if(ntf .eq. nt) then
         call fineoprlayer(a,jrf)
         irf = mod(irf-1, 3) + 1
         jrf = mod(jrf-1, 3) + 1
         krf = mod(krf-1, 3) + 1
      end if
 
      do 60  id=1,ndo
 
      if(ntf .eq. nt) call fineoprdiag(a,d,id,irc+irc-1)
 
      jd = id
      if(ntf .eq. nt) jd = 1
 
      do 60   j=1,9
 
c     Coarse from coarse.
 
      do 10 i2c=1,ntc+1
      i2f = i2c + i2c - 1
      do 10 i1c=1,ntc+1
      i1f = i1c + i1c - 1
 
      amt =                   a(i1f,i2f,jd,j,1,1,jrf)
     &    + w(i1c,i2c,id,2,1)*a(i1f,i2f,jd,j,1,2,jrf)
      amb = w(i1c,i2c,id,3,1)*a(i1f,i2f,jd,j,1,2,jrf)
     &    +                   a(i1f,i2f,jd,j,1,3,jrf)
 
      b(i1c,i2c,id,j,1,2,irc) = ((b(i1c,i2c,id,j,1,2,irc)
     &                          + a(i1f,i2f,jd,j,1,2,irf))
     &       + w(i1c,i2c,id,2,1)*(a(i1f,i2f,jd,j,1,3,irf) + amt))
      b(i1c,i2c,id,j,1,3,irc) = ((b(i1c,i2c,id,j,1,3,irc)
     &       + w(i1c,i2c,id,3,1)* a(i1f,i2f,jd,j,1,3,irf))
     &       + w(i1c,i2c,id,2,1)* amb)
      b(i1c,i2c,id,j,1,1,jrc) = ((b(i1c,i2c,id,j,1,1,jrc)
     &       + w(i1c,i2c,id,2,1)* a(i1f,i2f,jd,j,1,1,krf))
     &       + w(i1c,i2c,id,3,1)* amt)
      b(i1c,i2c,id,j,1,2,jrc) =  (b(i1c,i2c,id,j,1,2,jrc)
     &       + w(i1c,i2c,id,3,1)*(a(i1f,i2f,jd,j,1,1,krf) + amb))
 
 10   continue
 
      do 20   m=2,7
      do 20 i2c=1,ntc+1
      j2c = i2c + j2n(m)
      i2f = i2c + i2c - 1
      do 20 i1c=1,ntc+1
      j1c = i1c + j1n(m)
      i1f = i1c + i1c - 1
 
      amt = ((w(i1c,i2c,id,1,   m) *a(i1f,i2f,jd,j,m,1,jrf))
     &      + w(i1c,i2c,id,2,   m) *a(i1f,i2f,jd,j,m,2,jrf))
      amb = ((w(i1c,i2c,id,3,   m) *a(i1f,i2f,jd,j,m,2,jrf))
     &      + w(i1c,i2c,id,4,   m) *a(i1f,i2f,jd,j,m,3,jrf))
      bmt = ((w(j1c,j2c,id,1,md(m))*a(i1f,i2f,jd,j,m,1,jrf))
     &      + w(j1c,j2c,id,2,md(m))*a(i1f,i2f,jd,j,m,2,jrf))
      bmb = ((w(j1c,j2c,id,3,md(m))*a(i1f,i2f,jd,j,m,2,jrf))
     &      + w(j1c,j2c,id,4,md(m))*a(i1f,i2f,jd,j,m,3,jrf))
 
      b(i1c,i2c,id,j,1,2,irc) =  (((b(i1c,i2c,id,j,1,2,irc)
     &      + w(i1c,i2c,id,1,   m) *a(i1f,i2f,jd,j,m,2,irf))
     &      + w(i1c,i2c,id,2,   m) *a(i1f,i2f,jd,j,m,3,irf))
     &      + w(i1c,i2c,id,2,   1) *amt)
      b(i1c,i2c,id,j,1,3,irc) =   ((b(i1c,i2c,id,j,1,3,irc)
     &      + w(i1c,i2c,id,3,   m) *a(i1f,i2f,jd,j,m,3,irf))
     &      + w(i1c,i2c,id,2,   1) *amb)
      b(i1c,i2c,id,j,1,1,jrc) =   ((b(i1c,i2c,id,j,1,1,jrc)
     &      + w(i1c,i2c,id,2,   m) *a(i1f,i2f,jd,j,m,1,krf))
     &      + w(i1c,i2c,id,3,   1) *amt)
      b(i1c,i2c,id,j,1,2,jrc) =   ((b(i1c,i2c,id,j,1,2,jrc)
     &      + w(i1c,i2c,id,3,   m) *a(i1f,i2f,jd,j,m,1,krf))
     &      + w(i1c,i2c,id,3,   1) *amb)
      b(i1c,i2c,id,j,m,2,irc) =  (((b(i1c,i2c,id,j,m,2,irc)
     &      + w(j1c,j2c,id,1,md(m))*a(i1f,i2f,jd,j,m,2,irf))
     &      + w(j1c,j2c,id,2,md(m))*a(i1f,i2f,jd,j,m,3,irf))
     &      + w(i1c,i2c,id,2,   1) *bmt)
      b(i1c,i2c,id,j,m,3,irc) =   ((b(i1c,i2c,id,j,m,3,irc)
     &      + w(j1c,j2c,id,3,md(m))*a(i1f,i2f,jd,j,m,3,irf))
     &      + w(i1c,i2c,id,2,   1) *bmb)
      b(i1c,i2c,id,j,m,1,jrc) =   ((b(i1c,i2c,id,j,m,1,jrc)
     &      + w(j1c,j2c,id,2,md(m))*a(i1f,i2f,jd,j,m,1,krf))
     &      + w(i1c,i2c,id,3,   1) *bmt)
      b(i1c,i2c,id,j,m,2,jrc) =   ((b(i1c,i2c,id,j,m,2,jrc)
     &      + w(j1c,j2c,id,3,md(m))*a(i1f,i2f,jd,j,m,1,krf))
     &      + w(i1c,i2c,id,3,   1) *bmb)
 
 20   continue
 
c     Project fine node values to coarse node values.
 
      do 30 i2c=1,ntc+1
      j2c = i2c - 1
      k2c = i2c + 1
      i2f = i2c + i2c - 1
      do 30 i1c=2,ntc+1
      j1c = i1c - 1
      j1f = i1c + i1c - 2
 
      att = ((((((((w(i1c,i2c,id,1,5)*a(j1f,i2f,jd,j,1,2,irf))
     &            +                   a(j1f,i2f,jd,j,2,2,irf))
     &            + w(i1c,i2c,id,1,4)*a(j1f,i2f,jd,j,3,2,irf))
     &            + w(i1c,i2c,id,1,6)*a(j1f,i2f,jd,j,7,2,irf))
     &            + w(i1c,i2c,id,2,5)*a(j1f,i2f,jd,j,1,3,irf))
     &            + w(i1c,i2c,id,2,1)*a(j1f,i2f,jd,j,2,3,irf))
     &            + w(i1c,i2c,id,2,4)*a(j1f,i2f,jd,j,3,3,irf))
     &            + w(i1c,i2c,id,2,6)*a(j1f,i2f,jd,j,7,3,irf))
      amt = ((((((((w(i1c,i2c,id,1,5)*a(j1f,i2f,jd,j,1,1,jrf))
     &            +                   a(j1f,i2f,jd,j,2,1,jrf))
     &            + w(i1c,i2c,id,1,4)*a(j1f,i2f,jd,j,3,1,jrf))
     &            + w(i1c,i2c,id,1,6)*a(j1f,i2f,jd,j,7,1,jrf))
     &            + w(i1c,i2c,id,2,5)*a(j1f,i2f,jd,j,1,2,jrf))
     &            + w(i1c,i2c,id,2,1)*a(j1f,i2f,jd,j,2,2,jrf))
     &            + w(i1c,i2c,id,2,4)*a(j1f,i2f,jd,j,3,2,jrf))
     &            + w(i1c,i2c,id,2,6)*a(j1f,i2f,jd,j,7,2,jrf))
      abt =     ((((w(i1c,i2c,id,2,5)*a(j1f,i2f,jd,j,1,1,krf))
     &            + w(i1c,i2c,id,2,1)*a(j1f,i2f,jd,j,2,1,krf))
     &            + w(i1c,i2c,id,2,4)*a(j1f,i2f,jd,j,3,1,krf))
     &            + w(i1c,i2c,id,2,6)*a(j1f,i2f,jd,j,7,1,krf))
      atb =     ((((w(i1c,i2c,id,3,5)*a(j1f,i2f,jd,j,1,3,irf))
     &            + w(i1c,i2c,id,3,1)*a(j1f,i2f,jd,j,2,3,irf))
     &            + w(i1c,i2c,id,3,4)*a(j1f,i2f,jd,j,3,3,irf))
     &            + w(i1c,i2c,id,3,6)*a(j1f,i2f,jd,j,7,3,irf))
      amb = ((((((((w(i1c,i2c,id,3,5)*a(j1f,i2f,jd,j,1,2,jrf))
     &            + w(i1c,i2c,id,3,1)*a(j1f,i2f,jd,j,2,2,jrf))
     &            + w(i1c,i2c,id,3,4)*a(j1f,i2f,jd,j,3,2,jrf))
     &            + w(i1c,i2c,id,3,6)*a(j1f,i2f,jd,j,7,2,jrf))
     &            + w(i1c,i2c,id,4,5)*a(j1f,i2f,jd,j,1,3,jrf))
     &            +                   a(j1f,i2f,jd,j,2,3,jrf))
     &            + w(i1c,i2c,id,4,4)*a(j1f,i2f,jd,j,3,3,jrf))
     &            + w(i1c,i2c,id,4,6)*a(j1f,i2f,jd,j,7,3,jrf))
      abb =     ((((w(i1c,i2c,id,3,5)*a(j1f,i2f,jd,j,1,1,krf))
     &            + w(i1c,i2c,id,3,1)*a(j1f,i2f,jd,j,2,1,krf))
     &            + w(i1c,i2c,id,3,4)*a(j1f,i2f,jd,j,3,1,krf))
     &            + w(i1c,i2c,id,3,6)*a(j1f,i2f,jd,j,7,1,krf))
 
      b(i1c,i2c,id,j,1,2,irc) =     ((b(i1c,i2c,id,j,1,2,irc)
     &     + w(i1c,i2c,id,1,5)*att) + w(i1c,i2c,id,2,5)*amt)
      b(i1c,i2c,id,j,1,3,irc) =     ((b(i1c,i2c,id,j,1,3,irc)
     &     + w(i1c,i2c,id,1,5)*atb) + w(i1c,i2c,id,2,5)*amb)
      b(i1c,i2c,id,j,1,1,jrc) =     ((b(i1c,i2c,id,j,1,1,jrc)
     &     + w(i1c,i2c,id,4,5)*abt) + w(i1c,i2c,id,3,5)*amt)
      b(i1c,i2c,id,j,1,2,jrc) =     ((b(i1c,i2c,id,j,1,2,jrc)
     &     + w(i1c,i2c,id,4,5)*abb) + w(i1c,i2c,id,3,5)*amb)
      b(j1c,i2c,id,j,2,2,irc) =     ((b(j1c,i2c,id,j,2,2,irc)
     &     + w(j1c,i2c,id,1,2)*att) + w(j1c,i2c,id,2,2)*amt)
      b(j1c,i2c,id,j,2,3,irc) =     ((b(j1c,i2c,id,j,2,3,irc)
     &     + w(j1c,i2c,id,1,2)*atb) + w(j1c,i2c,id,2,2)*amb)
      b(j1c,i2c,id,j,2,1,jrc) =     ((b(j1c,i2c,id,j,2,1,jrc)
     &     + w(j1c,i2c,id,4,2)*abt) + w(j1c,i2c,id,3,2)*amt)
      b(j1c,i2c,id,j,2,2,jrc) =     ((b(j1c,i2c,id,j,2,2,jrc)
     &     + w(j1c,i2c,id,4,2)*abb) + w(j1c,i2c,id,3,2)*amb)
 
      btt = ((((w(j1c,k2c,id,1,7)*a(j1f,i2f,jd,j,3,2,irf))
     &        + w(j1c,k2c,id,1,6)*a(j1f,i2f,jd,j,4,2,irf))
     &        + w(j1c,k2c,id,2,7)*a(j1f,i2f,jd,j,3,3,irf))
     &        + w(j1c,k2c,id,2,6)*a(j1f,i2f,jd,j,4,3,irf))
      bmt = ((((w(j1c,k2c,id,1,7)*a(j1f,i2f,jd,j,3,1,jrf))
     &        + w(j1c,k2c,id,1,6)*a(j1f,i2f,jd,j,4,1,jrf))
     &        + w(j1c,k2c,id,2,7)*a(j1f,i2f,jd,j,3,2,jrf))
     &        + w(j1c,k2c,id,2,6)*a(j1f,i2f,jd,j,4,2,jrf))
      bbt =   ((w(j1c,k2c,id,2,7)*a(j1f,i2f,jd,j,3,1,krf))
     &        + w(j1c,k2c,id,2,6)*a(j1f,i2f,jd,j,4,1,krf))
      btb =   ((w(j1c,k2c,id,3,7)*a(j1f,i2f,jd,j,3,3,irf))
     &        + w(j1c,k2c,id,3,6)*a(j1f,i2f,jd,j,4,3,irf))
      bmb = ((((w(j1c,k2c,id,3,7)*a(j1f,i2f,jd,j,3,2,jrf))
     &        + w(j1c,k2c,id,3,6)*a(j1f,i2f,jd,j,4,2,jrf))
     &        + w(j1c,k2c,id,4,7)*a(j1f,i2f,jd,j,3,3,jrf))
     &        + w(j1c,k2c,id,4,6)*a(j1f,i2f,jd,j,4,3,jrf))
      bbb =   ((w(j1c,k2c,id,3,7)*a(j1f,i2f,jd,j,3,1,krf))
     &        + w(j1c,k2c,id,3,6)*a(j1f,i2f,jd,j,4,1,krf))
 
      b(i1c,i2c,id,j,4,2,irc) =     ((b(i1c,i2c,id,j,4,2,irc)
     &     + w(i1c,i2c,id,1,5)*btt) + w(i1c,i2c,id,2,5)*bmt)
      b(i1c,i2c,id,j,4,3,irc) =     ((b(i1c,i2c,id,j,4,3,irc)
     &     + w(i1c,i2c,id,1,5)*btb) + w(i1c,i2c,id,2,5)*bmb)
      b(i1c,i2c,id,j,4,1,jrc) =     ((b(i1c,i2c,id,j,4,1,jrc)
     &     + w(i1c,i2c,id,4,5)*bbt) + w(i1c,i2c,id,3,5)*bmt)
      b(i1c,i2c,id,j,4,2,jrc) =     ((b(i1c,i2c,id,j,4,2,jrc)
     &     + w(i1c,i2c,id,4,5)*bbb) + w(i1c,i2c,id,3,5)*bmb)
      b(j1c,i2c,id,j,3,2,irc) =     ((b(j1c,i2c,id,j,3,2,irc)
     &     + w(j1c,i2c,id,1,2)*btt) + w(j1c,i2c,id,2,2)*bmt)
      b(j1c,i2c,id,j,3,3,irc) =     ((b(j1c,i2c,id,j,3,3,irc)
     &     + w(j1c,i2c,id,1,2)*btb) + w(j1c,i2c,id,2,2)*bmb)
      b(j1c,i2c,id,j,3,1,jrc) =     ((b(j1c,i2c,id,j,3,1,jrc)
     &     + w(j1c,i2c,id,4,2)*bbt) + w(j1c,i2c,id,3,2)*bmt)
      b(j1c,i2c,id,j,3,2,jrc) =     ((b(j1c,i2c,id,j,3,2,jrc)
     &     + w(j1c,i2c,id,4,2)*bbb) + w(j1c,i2c,id,3,2)*bmb)
 
      ctt = ((((((((w(j1c,i2c,id,1,2)*a(j1f,i2f,jd,j,1,2,irf))
     &            + w(j1c,i2c,id,1,3)*a(j1f,i2f,jd,j,4,2,irf))
     &            +                   a(j1f,i2f,jd,j,5,2,irf))
     &            + w(j1c,i2c,id,1,7)*a(j1f,i2f,jd,j,6,2,irf))
     &            + w(j1c,i2c,id,2,2)*a(j1f,i2f,jd,j,1,3,irf))
     &            + w(j1c,i2c,id,2,3)*a(j1f,i2f,jd,j,4,3,irf))
     &            + w(j1c,i2c,id,2,1)*a(j1f,i2f,jd,j,5,3,irf))
     &            + w(j1c,i2c,id,2,7)*a(j1f,i2f,jd,j,6,3,irf))
      cmt = ((((((((w(j1c,i2c,id,1,2)*a(j1f,i2f,jd,j,1,1,jrf))
     &            + w(j1c,i2c,id,1,3)*a(j1f,i2f,jd,j,4,1,jrf))
     &            +                   a(j1f,i2f,jd,j,5,1,jrf))
     &            + w(j1c,i2c,id,1,7)*a(j1f,i2f,jd,j,6,1,jrf))
     &            + w(j1c,i2c,id,2,2)*a(j1f,i2f,jd,j,1,2,jrf))
     &            + w(j1c,i2c,id,2,3)*a(j1f,i2f,jd,j,4,2,jrf))
     &            + w(j1c,i2c,id,2,1)*a(j1f,i2f,jd,j,5,2,jrf))
     &            + w(j1c,i2c,id,2,7)*a(j1f,i2f,jd,j,6,2,jrf))
      cbt =     ((((w(j1c,i2c,id,2,2)*a(j1f,i2f,jd,j,1,1,krf))
     &            + w(j1c,i2c,id,2,3)*a(j1f,i2f,jd,j,4,1,krf))
     &            + w(j1c,i2c,id,2,1)*a(j1f,i2f,jd,j,5,1,krf))
     &            + w(j1c,i2c,id,2,7)*a(j1f,i2f,jd,j,6,1,krf))
      ctb =     ((((w(j1c,i2c,id,3,2)*a(j1f,i2f,jd,j,1,3,irf))
     &            + w(j1c,i2c,id,3,3)*a(j1f,i2f,jd,j,4,3,irf))
     &            + w(j1c,i2c,id,3,1)*a(j1f,i2f,jd,j,5,3,irf))
     &            + w(j1c,i2c,id,3,7)*a(j1f,i2f,jd,j,6,3,irf))
      cmb = ((((((((w(j1c,i2c,id,3,2)*a(j1f,i2f,jd,j,1,2,jrf))
     &            + w(j1c,i2c,id,3,3)*a(j1f,i2f,jd,j,4,2,jrf))
     &            + w(j1c,i2c,id,3,1)*a(j1f,i2f,jd,j,5,2,jrf))
     &            + w(j1c,i2c,id,3,7)*a(j1f,i2f,jd,j,6,2,jrf))
     &            + w(j1c,i2c,id,4,2)*a(j1f,i2f,jd,j,1,3,jrf))
     &            + w(j1c,i2c,id,4,3)*a(j1f,i2f,jd,j,4,3,jrf))
     &            +                   a(j1f,i2f,jd,j,5,3,jrf))
     &            + w(j1c,i2c,id,4,7)*a(j1f,i2f,jd,j,6,3,jrf))
      cbb =     ((((w(j1c,i2c,id,3,2)*a(j1f,i2f,jd,j,1,1,krf))
     &            + w(j1c,i2c,id,3,3)*a(j1f,i2f,jd,j,4,1,krf))
     &            + w(j1c,i2c,id,3,1)*a(j1f,i2f,jd,j,5,1,krf))
     &            + w(j1c,i2c,id,3,7)*a(j1f,i2f,jd,j,6,1,krf))
 
      b(i1c,i2c,id,j,5,2,irc) =     ((b(i1c,i2c,id,j,5,2,irc)
     &     + w(i1c,i2c,id,1,5)*ctt) + w(i1c,i2c,id,2,5)*cmt)
      b(i1c,i2c,id,j,5,3,irc) =     ((b(i1c,i2c,id,j,5,3,irc)
     &     + w(i1c,i2c,id,1,5)*ctb) + w(i1c,i2c,id,2,5)*cmb)
      b(i1c,i2c,id,j,5,1,jrc) =     ((b(i1c,i2c,id,j,5,1,jrc)
     &     + w(i1c,i2c,id,4,5)*cbt) + w(i1c,i2c,id,3,5)*cmt)
      b(i1c,i2c,id,j,5,2,jrc) =     ((b(i1c,i2c,id,j,5,2,jrc)
     &     + w(i1c,i2c,id,4,5)*cbb) + w(i1c,i2c,id,3,5)*cmb)
      b(j1c,i2c,id,j,1,2,irc) =     ((b(j1c,i2c,id,j,1,2,irc)
     &     + w(j1c,i2c,id,1,2)*ctt) + w(j1c,i2c,id,2,2)*cmt)
      b(j1c,i2c,id,j,1,3,irc) =     ((b(j1c,i2c,id,j,1,3,irc)
     &     + w(j1c,i2c,id,1,2)*ctb) + w(j1c,i2c,id,2,2)*cmb)
      b(j1c,i2c,id,j,1,1,jrc) =     ((b(j1c,i2c,id,j,1,1,jrc)
     &     + w(j1c,i2c,id,4,2)*cbt) + w(j1c,i2c,id,3,2)*cmt)
      b(j1c,i2c,id,j,1,2,jrc) =     ((b(j1c,i2c,id,j,1,2,jrc)
     &     + w(j1c,i2c,id,4,2)*cbb) + w(j1c,i2c,id,3,2)*cmb)
 
      dtt = ((((w(i1c,j2c,id,1,4)*a(j1f,i2f,jd,j,6,2,irf))
     &        + w(i1c,j2c,id,1,3)*a(j1f,i2f,jd,j,7,2,irf))
     &        + w(i1c,j2c,id,2,4)*a(j1f,i2f,jd,j,6,3,irf))
     &        + w(i1c,j2c,id,2,3)*a(j1f,i2f,jd,j,7,3,irf))
      dmt = ((((w(i1c,j2c,id,1,4)*a(j1f,i2f,jd,j,6,1,jrf))
     &        + w(i1c,j2c,id,1,3)*a(j1f,i2f,jd,j,7,1,jrf))
     &        + w(i1c,j2c,id,2,4)*a(j1f,i2f,jd,j,6,2,jrf))
     &        + w(i1c,j2c,id,2,3)*a(j1f,i2f,jd,j,7,2,jrf))
      dbt =   ((w(i1c,j2c,id,2,4)*a(j1f,i2f,jd,j,6,1,krf))
     &        + w(i1c,j2c,id,2,3)*a(j1f,i2f,jd,j,7,1,krf))
      dtb =   ((w(i1c,j2c,id,3,4)*a(j1f,i2f,jd,j,6,3,irf))
     &        + w(i1c,j2c,id,3,3)*a(j1f,i2f,jd,j,7,3,irf))
      dmb = ((((w(i1c,j2c,id,3,4)*a(j1f,i2f,jd,j,6,2,jrf))
     &        + w(i1c,j2c,id,3,3)*a(j1f,i2f,jd,j,7,2,jrf))
     &        + w(i1c,j2c,id,4,4)*a(j1f,i2f,jd,j,6,3,jrf))
     &        + w(i1c,j2c,id,4,3)*a(j1f,i2f,jd,j,7,3,jrf))
      dbb =   ((w(i1c,j2c,id,3,4)*a(j1f,i2f,jd,j,6,1,krf))
     &        + w(i1c,j2c,id,3,3)*a(j1f,i2f,jd,j,7,1,krf))
 
      b(i1c,i2c,id,j,6,2,irc) =     ((b(i1c,i2c,id,j,6,2,irc)
     &     + w(i1c,i2c,id,1,5)*dtt) + w(i1c,i2c,id,2,5)*dmt)
      b(i1c,i2c,id,j,6,3,irc) =     ((b(i1c,i2c,id,j,6,3,irc)
     &     + w(i1c,i2c,id,1,5)*dtb) + w(i1c,i2c,id,2,5)*dmb)
      b(i1c,i2c,id,j,6,1,jrc) =     ((b(i1c,i2c,id,j,6,1,jrc)
     &     + w(i1c,i2c,id,4,5)*dbt) + w(i1c,i2c,id,3,5)*dmt)
      b(i1c,i2c,id,j,6,2,jrc) =     ((b(i1c,i2c,id,j,6,2,jrc)
     &     + w(i1c,i2c,id,4,5)*dbb) + w(i1c,i2c,id,3,5)*dmb)
      b(j1c,i2c,id,j,7,2,irc) =     ((b(j1c,i2c,id,j,7,2,irc)
     &     + w(j1c,i2c,id,1,2)*dtt) + w(j1c,i2c,id,2,2)*dmt)
      b(j1c,i2c,id,j,7,3,irc) =     ((b(j1c,i2c,id,j,7,3,irc)
     &     + w(j1c,i2c,id,1,2)*dtb) + w(j1c,i2c,id,2,2)*dmb)
      b(j1c,i2c,id,j,7,1,jrc) =     ((b(j1c,i2c,id,j,7,1,jrc)
     &     + w(j1c,i2c,id,4,2)*dbt) + w(j1c,i2c,id,3,2)*dmt)
      b(j1c,i2c,id,j,7,2,jrc) =     ((b(j1c,i2c,id,j,7,2,jrc)
     &     + w(j1c,i2c,id,4,2)*dbb) + w(j1c,i2c,id,3,2)*dmb)
 
 30   continue
 
      do 40 i2c=1,ntc
      k2c = i2c + 1
      k2f = i2c + i2c
      do 40 i1c=2,ntc+1
      j1c = i1c - 1
      j1f = i1c + i1c - 2
      ett = ((((((((w(i1c,i2c,id,1,4)*a(j1f,k2f,jd,j,1,2,irf))
     &            + w(i1c,i2c,id,1,3)*a(j1f,k2f,jd,j,2,2,irf))
     &            + w(i1c,i2c,id,1,5)*a(j1f,k2f,jd,j,6,2,irf))
     &            +                   a(j1f,k2f,jd,j,7,2,irf))
     &            + w(i1c,i2c,id,2,4)*a(j1f,k2f,jd,j,1,3,irf))
     &            + w(i1c,i2c,id,2,3)*a(j1f,k2f,jd,j,2,3,irf))
     &            + w(i1c,i2c,id,2,5)*a(j1f,k2f,jd,j,6,3,irf))
     &            + w(i1c,i2c,id,2,1)*a(j1f,k2f,jd,j,7,3,irf))
      emt = ((((((((w(i1c,i2c,id,1,4)*a(j1f,k2f,jd,j,1,1,jrf))
     &            + w(i1c,i2c,id,1,3)*a(j1f,k2f,jd,j,2,1,jrf))
     &            + w(i1c,i2c,id,1,5)*a(j1f,k2f,jd,j,6,1,jrf))
     &            +                   a(j1f,k2f,jd,j,7,1,jrf))
     &            + w(i1c,i2c,id,2,4)*a(j1f,k2f,jd,j,1,2,jrf))
     &            + w(i1c,i2c,id,2,3)*a(j1f,k2f,jd,j,2,2,jrf))
     &            + w(i1c,i2c,id,2,5)*a(j1f,k2f,jd,j,6,2,jrf))
     &            + w(i1c,i2c,id,2,1)*a(j1f,k2f,jd,j,7,2,jrf))
      ebt =     ((((w(i1c,i2c,id,2,4)*a(j1f,k2f,jd,j,1,1,krf))
     &            + w(i1c,i2c,id,2,3)*a(j1f,k2f,jd,j,2,1,krf))
     &            + w(i1c,i2c,id,2,5)*a(j1f,k2f,jd,j,6,1,krf))
     &            + w(i1c,i2c,id,2,1)*a(j1f,k2f,jd,j,7,1,krf))
      etb =     ((((w(i1c,i2c,id,3,4)*a(j1f,k2f,jd,j,1,3,irf))
     &            + w(i1c,i2c,id,3,3)*a(j1f,k2f,jd,j,2,3,irf))
     &            + w(i1c,i2c,id,3,5)*a(j1f,k2f,jd,j,6,3,irf))
     &            + w(i1c,i2c,id,3,1)*a(j1f,k2f,jd,j,7,3,irf))
      emb = ((((((((w(i1c,i2c,id,3,4)*a(j1f,k2f,jd,j,1,2,jrf))
     &            + w(i1c,i2c,id,3,3)*a(j1f,k2f,jd,j,2,2,jrf))
     &            + w(i1c,i2c,id,3,5)*a(j1f,k2f,jd,j,6,2,jrf))
     &            + w(i1c,i2c,id,3,1)*a(j1f,k2f,jd,j,7,2,jrf))
     &            + w(i1c,i2c,id,4,4)*a(j1f,k2f,jd,j,1,3,jrf))
     &            + w(i1c,i2c,id,4,3)*a(j1f,k2f,jd,j,2,3,jrf))
     &            + w(i1c,i2c,id,4,5)*a(j1f,k2f,jd,j,6,3,jrf))
     &            +                   a(j1f,k2f,jd,j,7,3,jrf))
      ebb =     ((((w(i1c,i2c,id,3,4)*a(j1f,k2f,jd,j,1,1,krf))
     &            + w(i1c,i2c,id,3,3)*a(j1f,k2f,jd,j,2,1,krf))
     &            + w(i1c,i2c,id,3,5)*a(j1f,k2f,jd,j,6,1,krf))
     &            + w(i1c,i2c,id,3,1)*a(j1f,k2f,jd,j,7,1,krf))
 
      b(i1c,i2c,id,j,1,2,irc) =     ((b(i1c,i2c,id,j,1,2,irc)
     &     + w(i1c,i2c,id,1,4)*ett) + w(i1c,i2c,id,2,4)*emt)
      b(i1c,i2c,id,j,1,3,irc) =     ((b(i1c,i2c,id,j,1,3,irc)
     &     + w(i1c,i2c,id,1,4)*etb) + w(i1c,i2c,id,2,4)*emb)
      b(i1c,i2c,id,j,1,1,jrc) =     ((b(i1c,i2c,id,j,1,1,jrc)
     &     + w(i1c,i2c,id,4,4)*ebt) + w(i1c,i2c,id,3,4)*emt)
      b(i1c,i2c,id,j,1,2,jrc) =     ((b(i1c,i2c,id,j,1,2,jrc)
     &     + w(i1c,i2c,id,4,4)*ebb) + w(i1c,i2c,id,3,4)*emb)
      b(j1c,k2c,id,j,7,2,irc) =     ((b(j1c,k2c,id,j,7,2,irc)
     &     + w(j1c,k2c,id,1,7)*ett) + w(j1c,k2c,id,2,7)*emt)
      b(j1c,k2c,id,j,7,3,irc) =     ((b(j1c,k2c,id,j,7,3,irc)
     &     + w(j1c,k2c,id,1,7)*etb) + w(j1c,k2c,id,2,7)*emb)
      b(j1c,k2c,id,j,7,1,jrc) =     ((b(j1c,k2c,id,j,7,1,jrc)
     &     + w(j1c,k2c,id,4,7)*ebt) + w(j1c,k2c,id,3,7)*emt)
      b(j1c,k2c,id,j,7,2,jrc) =     ((b(j1c,k2c,id,j,7,2,jrc)
     &     + w(j1c,k2c,id,4,7)*ebb) + w(j1c,k2c,id,3,7)*emb)
 
      ftt = ((((w(i1c,k2c,id,1,6)*a(j1f,k2f,jd,j,2,2,irf))
     &        + w(i1c,k2c,id,1,5)*a(j1f,k2f,jd,j,3,2,irf))
     &        + w(i1c,k2c,id,2,6)*a(j1f,k2f,jd,j,2,3,irf))
     &        + w(i1c,k2c,id,2,5)*a(j1f,k2f,jd,j,3,3,irf))
      fmt = ((((w(i1c,k2c,id,1,6)*a(j1f,k2f,jd,j,2,1,jrf))
     &        + w(i1c,k2c,id,1,5)*a(j1f,k2f,jd,j,3,1,jrf))
     &        + w(i1c,k2c,id,2,6)*a(j1f,k2f,jd,j,2,2,jrf))
     &        + w(i1c,k2c,id,2,5)*a(j1f,k2f,jd,j,3,2,jrf))
      fbt =   ((w(i1c,k2c,id,2,6)*a(j1f,k2f,jd,j,2,1,krf))
     &        + w(i1c,k2c,id,2,5)*a(j1f,k2f,jd,j,3,1,krf))
      ftb =   ((w(i1c,k2c,id,3,6)*a(j1f,k2f,jd,j,2,3,irf))
     &        + w(i1c,k2c,id,3,5)*a(j1f,k2f,jd,j,3,3,irf))
      fmb = ((((w(i1c,k2c,id,3,6)*a(j1f,k2f,jd,j,2,2,jrf))
     &        + w(i1c,k2c,id,3,5)*a(j1f,k2f,jd,j,3,2,jrf))
     &        + w(i1c,k2c,id,4,6)*a(j1f,k2f,jd,j,2,3,jrf))
     &        + w(i1c,k2c,id,4,5)*a(j1f,k2f,jd,j,3,3,jrf))
      fbb =   ((w(i1c,k2c,id,3,6)*a(j1f,k2f,jd,j,2,1,krf))
     &        + w(i1c,k2c,id,3,5)*a(j1f,k2f,jd,j,3,1,krf))
 
      b(i1c,i2c,id,j,3,2,irc) =     ((b(i1c,i2c,id,j,3,2,irc)
     &     + w(i1c,i2c,id,1,4)*ftt) + w(i1c,i2c,id,2,4)*fmt)
      b(i1c,i2c,id,j,3,3,irc) =     ((b(i1c,i2c,id,j,3,3,irc)
     &     + w(i1c,i2c,id,1,4)*ftb) + w(i1c,i2c,id,2,4)*fmb)
      b(i1c,i2c,id,j,3,1,jrc) =     ((b(i1c,i2c,id,j,3,1,jrc)
     &     + w(i1c,i2c,id,4,4)*fbt) + w(i1c,i2c,id,3,4)*fmt)
      b(i1c,i2c,id,j,3,2,jrc) =     ((b(i1c,i2c,id,j,3,2,jrc)
     &     + w(i1c,i2c,id,4,4)*fbb) + w(i1c,i2c,id,3,4)*fmb)
      b(j1c,k2c,id,j,2,2,irc) =     ((b(j1c,k2c,id,j,2,2,irc)
     &     + w(j1c,k2c,id,1,7)*ftt) + w(j1c,k2c,id,2,7)*fmt)
      b(j1c,k2c,id,j,2,3,irc) =     ((b(j1c,k2c,id,j,2,3,irc)
     &     + w(j1c,k2c,id,1,7)*ftb) + w(j1c,k2c,id,2,7)*fmb)
      b(j1c,k2c,id,j,2,1,jrc) =     ((b(j1c,k2c,id,j,2,1,jrc)
     &     + w(j1c,k2c,id,4,7)*fbt) + w(j1c,k2c,id,3,7)*fmt)
      b(j1c,k2c,id,j,2,2,jrc) =     ((b(j1c,k2c,id,j,2,2,jrc)
     &     + w(j1c,k2c,id,4,7)*fbb) + w(j1c,k2c,id,3,7)*fmb)
 
      gtt = ((((((((w(j1c,k2c,id,1,7)*a(j1f,k2f,jd,j,1,2,irf))
     &            + w(j1c,k2c,id,1,2)*a(j1f,k2f,jd,j,3,2,irf))
     &            +                   a(j1f,k2f,jd,j,4,2,irf))
     &            + w(j1c,k2c,id,1,6)*a(j1f,k2f,jd,j,5,2,irf))
     &            + w(j1c,k2c,id,2,7)*a(j1f,k2f,jd,j,1,3,irf))
     &            + w(j1c,k2c,id,2,2)*a(j1f,k2f,jd,j,3,3,irf))
     &            + w(j1c,k2c,id,2,1)*a(j1f,k2f,jd,j,4,3,irf))
     &            + w(j1c,k2c,id,2,6)*a(j1f,k2f,jd,j,5,3,irf))
      gmt = ((((((((w(j1c,k2c,id,1,7)*a(j1f,k2f,jd,j,1,1,jrf))
     &            + w(j1c,k2c,id,1,2)*a(j1f,k2f,jd,j,3,1,jrf))
     &            +                   a(j1f,k2f,jd,j,4,1,jrf))
     &            + w(j1c,k2c,id,1,6)*a(j1f,k2f,jd,j,5,1,jrf))
     &            + w(j1c,k2c,id,2,7)*a(j1f,k2f,jd,j,1,2,jrf))
     &            + w(j1c,k2c,id,2,2)*a(j1f,k2f,jd,j,3,2,jrf))
     &            + w(j1c,k2c,id,2,1)*a(j1f,k2f,jd,j,4,2,jrf))
     &            + w(j1c,k2c,id,2,6)*a(j1f,k2f,jd,j,5,2,jrf))
      gbt =     ((((w(j1c,k2c,id,2,7)*a(j1f,k2f,jd,j,1,1,krf))
     &            + w(j1c,k2c,id,2,2)*a(j1f,k2f,jd,j,3,1,krf))
     &            + w(j1c,k2c,id,2,1)*a(j1f,k2f,jd,j,4,1,krf))
     &            + w(j1c,k2c,id,2,6)*a(j1f,k2f,jd,j,5,1,krf))
      gtb =     ((((w(j1c,k2c,id,3,7)*a(j1f,k2f,jd,j,1,3,irf))
     &            + w(j1c,k2c,id,3,2)*a(j1f,k2f,jd,j,3,3,irf))
     &            + w(j1c,k2c,id,3,1)*a(j1f,k2f,jd,j,4,3,irf))
     &            + w(j1c,k2c,id,3,6)*a(j1f,k2f,jd,j,5,3,irf))
      gmb = ((((((((w(j1c,k2c,id,3,7)*a(j1f,k2f,jd,j,1,2,jrf))
     &            + w(j1c,k2c,id,3,2)*a(j1f,k2f,jd,j,3,2,jrf))
     &            + w(j1c,k2c,id,3,1)*a(j1f,k2f,jd,j,4,2,jrf))
     &            + w(j1c,k2c,id,3,6)*a(j1f,k2f,jd,j,5,2,jrf))
     &            + w(j1c,k2c,id,4,7)*a(j1f,k2f,jd,j,1,3,jrf))
     &            + w(j1c,k2c,id,4,2)*a(j1f,k2f,jd,j,3,3,jrf))
     &            +                   a(j1f,k2f,jd,j,4,3,jrf))
     &            + w(j1c,k2c,id,4,6)*a(j1f,k2f,jd,j,5,3,jrf))
      gbb =     ((((w(j1c,k2c,id,3,7)*a(j1f,k2f,jd,j,1,1,krf))
     &            + w(j1c,k2c,id,3,2)*a(j1f,k2f,jd,j,3,1,krf))
     &            + w(j1c,k2c,id,3,1)*a(j1f,k2f,jd,j,4,1,krf))
     &            + w(j1c,k2c,id,3,6)*a(j1f,k2f,jd,j,5,1,krf))
 
      b(i1c,i2c,id,j,4,2,irc) =     ((b(i1c,i2c,id,j,4,2,irc)
     &     + w(i1c,i2c,id,1,4)*gtt) + w(i1c,i2c,id,2,4)*gmt)
      b(i1c,i2c,id,j,4,3,irc) =     ((b(i1c,i2c,id,j,4,3,irc)
     &     + w(i1c,i2c,id,1,4)*gtb) + w(i1c,i2c,id,2,4)*gmb)
      b(i1c,i2c,id,j,4,1,jrc) =     ((b(i1c,i2c,id,j,4,1,jrc)
     &     + w(i1c,i2c,id,4,4)*gbt) + w(i1c,i2c,id,3,4)*gmt)
      b(i1c,i2c,id,j,4,2,jrc) =     ((b(i1c,i2c,id,j,4,2,jrc)
     &     + w(i1c,i2c,id,4,4)*gbb) + w(i1c,i2c,id,3,4)*gmb)
      b(j1c,k2c,id,j,1,2,irc) =     ((b(j1c,k2c,id,j,1,2,irc)
     &     + w(j1c,k2c,id,1,7)*gtt) + w(j1c,k2c,id,2,7)*gmt)
      b(j1c,k2c,id,j,1,3,irc) =     ((b(j1c,k2c,id,j,1,3,irc)
     &     + w(j1c,k2c,id,1,7)*gtb) + w(j1c,k2c,id,2,7)*gmb)
      b(j1c,k2c,id,j,1,1,jrc) =     ((b(j1c,k2c,id,j,1,1,jrc)
     &     + w(j1c,k2c,id,4,7)*gbt) + w(j1c,k2c,id,3,7)*gmt)
      b(j1c,k2c,id,j,1,2,jrc) =     ((b(j1c,k2c,id,j,1,2,jrc)
     &     + w(j1c,k2c,id,4,7)*gbb) + w(j1c,k2c,id,3,7)*gmb)
 
      htt = ((((w(j1c,i2c,id,1,3)*a(j1f,k2f,jd,j,5,2,irf))
     &        + w(j1c,i2c,id,1,2)*a(j1f,k2f,jd,j,6,2,irf))
     &        + w(j1c,i2c,id,2,3)*a(j1f,k2f,jd,j,5,3,irf))
     &        + w(j1c,i2c,id,2,2)*a(j1f,k2f,jd,j,6,3,irf))
      hmt = ((((w(j1c,i2c,id,1,3)*a(j1f,k2f,jd,j,5,1,jrf))
     &        + w(j1c,i2c,id,1,2)*a(j1f,k2f,jd,j,6,1,jrf))
     &        + w(j1c,i2c,id,2,3)*a(j1f,k2f,jd,j,5,2,jrf))
     &        + w(j1c,i2c,id,2,2)*a(j1f,k2f,jd,j,6,2,jrf))
      hbt =   ((w(j1c,i2c,id,2,3)*a(j1f,k2f,jd,j,5,1,krf))
     &        + w(j1c,i2c,id,2,2)*a(j1f,k2f,jd,j,6,1,krf))
      htb =   ((w(j1c,i2c,id,3,3)*a(j1f,k2f,jd,j,5,3,irf))
     &        + w(j1c,i2c,id,3,2)*a(j1f,k2f,jd,j,6,3,irf))
      hmb = ((((w(j1c,i2c,id,3,3)*a(j1f,k2f,jd,j,5,2,jrf))
     &        + w(j1c,i2c,id,3,2)*a(j1f,k2f,jd,j,6,2,jrf))
     &        + w(j1c,i2c,id,4,3)*a(j1f,k2f,jd,j,5,3,jrf))
     &        + w(j1c,i2c,id,4,2)*a(j1f,k2f,jd,j,6,3,jrf))
      hbb =   ((w(j1c,i2c,id,3,3)*a(j1f,k2f,jd,j,5,1,krf))
     &        + w(j1c,i2c,id,3,2)*a(j1f,k2f,jd,j,6,1,krf))
 
      b(i1c,i2c,id,j,5,2,irc) =     ((b(i1c,i2c,id,j,5,2,irc)
     &     + w(i1c,i2c,id,1,4)*htt) + w(i1c,i2c,id,2,4)*hmt)
      b(i1c,i2c,id,j,5,3,irc) =     ((b(i1c,i2c,id,j,5,3,irc)
     &     + w(i1c,i2c,id,1,4)*htb) + w(i1c,i2c,id,2,4)*hmb)
      b(i1c,i2c,id,j,5,1,jrc) =     ((b(i1c,i2c,id,j,5,1,jrc)
     &     + w(i1c,i2c,id,4,4)*hbt) + w(i1c,i2c,id,3,4)*hmt)
      b(i1c,i2c,id,j,5,2,jrc) =     ((b(i1c,i2c,id,j,5,2,jrc)
     &     + w(i1c,i2c,id,4,4)*hbb) + w(i1c,i2c,id,3,4)*hmb)
      b(j1c,k2c,id,j,6,2,irc) =     ((b(j1c,k2c,id,j,6,2,irc)
     &     + w(j1c,k2c,id,1,7)*htt) + w(j1c,k2c,id,2,7)*hmt)
      b(j1c,k2c,id,j,6,3,irc) =     ((b(j1c,k2c,id,j,6,3,irc)
     &     + w(j1c,k2c,id,1,7)*htb) + w(j1c,k2c,id,2,7)*hmb)
      b(j1c,k2c,id,j,6,1,jrc) =     ((b(j1c,k2c,id,j,6,1,jrc)
     &     + w(j1c,k2c,id,4,7)*hbt) + w(j1c,k2c,id,3,7)*hmt)
      b(j1c,k2c,id,j,6,2,jrc) =     ((b(j1c,k2c,id,j,6,2,jrc)
     &     + w(j1c,k2c,id,4,7)*hbb) + w(j1c,k2c,id,3,7)*hmb)
 
 40   continue
 
      do 50 i2c=1,ntc
      k2c = i2c + 1
      k2f = i2c + i2c
      do 50 i1c=1,ntc+1
      j1c = i1c - 1
      k1c = i1c + 1
      i1f = i1c + i1c - 1
 
      ptt = ((((((((w(i1c,i2c,id,1,3)*a(i1f,k2f,jd,j,1,2,irf))
     &            + w(i1c,i2c,id,1,4)*a(i1f,k2f,jd,j,5,2,irf))
     &            +                   a(i1f,k2f,jd,j,6,2,irf))
     &            + w(i1c,i2c,id,1,2)*a(i1f,k2f,jd,j,7,2,irf))
     &            + w(i1c,i2c,id,2,3)*a(i1f,k2f,jd,j,1,3,irf))
     &            + w(i1c,i2c,id,2,4)*a(i1f,k2f,jd,j,5,3,irf))
     &            + w(i1c,i2c,id,2,1)*a(i1f,k2f,jd,j,6,3,irf))
     &            + w(i1c,i2c,id,2,2)*a(i1f,k2f,jd,j,7,3,irf))
      pmt = ((((((((w(i1c,i2c,id,1,3)*a(i1f,k2f,jd,j,1,1,jrf))
     &            + w(i1c,i2c,id,1,4)*a(i1f,k2f,jd,j,5,1,jrf))
     &            +                   a(i1f,k2f,jd,j,6,1,jrf))
     &            + w(i1c,i2c,id,1,2)*a(i1f,k2f,jd,j,7,1,jrf))
     &            + w(i1c,i2c,id,2,3)*a(i1f,k2f,jd,j,1,2,jrf))
     &            + w(i1c,i2c,id,2,4)*a(i1f,k2f,jd,j,5,2,jrf))
     &            + w(i1c,i2c,id,2,1)*a(i1f,k2f,jd,j,6,2,jrf))
     &            + w(i1c,i2c,id,2,2)*a(i1f,k2f,jd,j,7,2,jrf))
      pbt =     ((((w(i1c,i2c,id,2,3)*a(i1f,k2f,jd,j,1,1,krf))
     &            + w(i1c,i2c,id,2,4)*a(i1f,k2f,jd,j,5,1,krf))
     &            + w(i1c,i2c,id,2,1)*a(i1f,k2f,jd,j,6,1,krf))
     &            + w(i1c,i2c,id,2,2)*a(i1f,k2f,jd,j,7,1,krf))
      ptb =     ((((w(i1c,i2c,id,3,3)*a(i1f,k2f,jd,j,1,3,irf))
     &            + w(i1c,i2c,id,3,4)*a(i1f,k2f,jd,j,5,3,irf))
     &            + w(i1c,i2c,id,3,1)*a(i1f,k2f,jd,j,6,3,irf))
     &            + w(i1c,i2c,id,3,2)*a(i1f,k2f,jd,j,7,3,irf))
      pmb = ((((((((w(i1c,i2c,id,3,3)*a(i1f,k2f,jd,j,1,2,jrf))
     &            + w(i1c,i2c,id,3,4)*a(i1f,k2f,jd,j,5,2,jrf))
     &            + w(i1c,i2c,id,3,1)*a(i1f,k2f,jd,j,6,2,jrf))
     &            + w(i1c,i2c,id,3,2)*a(i1f,k2f,jd,j,7,2,jrf))
     &            + w(i1c,i2c,id,4,3)*a(i1f,k2f,jd,j,1,3,jrf))
     &            + w(i1c,i2c,id,4,4)*a(i1f,k2f,jd,j,5,3,jrf))
     &            +                   a(i1f,k2f,jd,j,6,3,jrf))
     &            + w(i1c,i2c,id,4,2)*a(i1f,k2f,jd,j,7,3,jrf))
      pbb =     ((((w(i1c,i2c,id,3,3)*a(i1f,k2f,jd,j,1,1,krf))
     &            + w(i1c,i2c,id,3,4)*a(i1f,k2f,jd,j,5,1,krf))
     &            + w(i1c,i2c,id,3,1)*a(i1f,k2f,jd,j,6,1,krf))
     &            + w(i1c,i2c,id,3,2)*a(i1f,k2f,jd,j,7,1,krf))
 
      b(i1c,i2c,id,j,1,2,irc) =     ((b(i1c,i2c,id,j,1,2,irc)
     &     + w(i1c,i2c,id,1,3)*ptt) + w(i1c,i2c,id,2,3)*pmt)
      b(i1c,i2c,id,j,1,3,irc) =     ((b(i1c,i2c,id,j,1,3,irc)
     &     + w(i1c,i2c,id,1,3)*ptb) + w(i1c,i2c,id,2,3)*pmb)
      b(i1c,i2c,id,j,1,1,jrc) =     ((b(i1c,i2c,id,j,1,1,jrc)
     &     + w(i1c,i2c,id,4,3)*pbt) + w(i1c,i2c,id,3,3)*pmt)
      b(i1c,i2c,id,j,1,2,jrc) =     ((b(i1c,i2c,id,j,1,2,jrc)
     &     + w(i1c,i2c,id,4,3)*pbb) + w(i1c,i2c,id,3,3)*pmb)
      b(i1c,k2c,id,j,6,2,irc) =     ((b(i1c,k2c,id,j,6,2,irc)
     &     + w(i1c,k2c,id,1,6)*ptt) + w(i1c,k2c,id,2,6)*pmt)
      b(i1c,k2c,id,j,6,3,irc) =     ((b(i1c,k2c,id,j,6,3,irc)
     &     + w(i1c,k2c,id,1,6)*ptb) + w(i1c,k2c,id,2,6)*pmb)
      b(i1c,k2c,id,j,6,1,jrc) =     ((b(i1c,k2c,id,j,6,1,jrc)
     &     + w(i1c,k2c,id,4,6)*pbt) + w(i1c,k2c,id,3,6)*pmt)
      b(i1c,k2c,id,j,6,2,jrc) =     ((b(i1c,k2c,id,j,6,2,jrc)
     &     + w(i1c,k2c,id,4,6)*pbb) + w(i1c,k2c,id,3,6)*pmb)
 
      qtt = ((((w(k1c,i2c,id,1,4)*a(i1f,k2f,jd,j,2,2,irf))
     &        + w(k1c,i2c,id,1,5)*a(i1f,k2f,jd,j,7,2,irf))
     &        + w(k1c,i2c,id,2,4)*a(i1f,k2f,jd,j,2,3,irf))
     &        + w(k1c,i2c,id,2,5)*a(i1f,k2f,jd,j,7,3,irf))
      qmt = ((((w(k1c,i2c,id,1,4)*a(i1f,k2f,jd,j,2,1,jrf))
     &        + w(k1c,i2c,id,1,5)*a(i1f,k2f,jd,j,7,1,jrf))
     &        + w(k1c,i2c,id,2,4)*a(i1f,k2f,jd,j,2,2,jrf))
     &        + w(k1c,i2c,id,2,5)*a(i1f,k2f,jd,j,7,2,jrf))
      qbt =   ((w(k1c,i2c,id,2,4)*a(i1f,k2f,jd,j,2,1,krf))
     &        + w(k1c,i2c,id,2,5)*a(i1f,k2f,jd,j,7,1,krf))
      qtb =   ((w(k1c,i2c,id,3,4)*a(i1f,k2f,jd,j,2,3,irf))
     &        + w(k1c,i2c,id,3,5)*a(i1f,k2f,jd,j,7,3,irf))
      qmb = ((((w(k1c,i2c,id,3,4)*a(i1f,k2f,jd,j,2,2,jrf))
     &        + w(k1c,i2c,id,3,5)*a(i1f,k2f,jd,j,7,2,jrf))
     &        + w(k1c,i2c,id,4,4)*a(i1f,k2f,jd,j,2,3,jrf))
     &        + w(k1c,i2c,id,4,5)*a(i1f,k2f,jd,j,7,3,jrf))
      qbb =   ((w(k1c,i2c,id,3,4)*a(i1f,k2f,jd,j,2,1,krf))
     &        + w(k1c,i2c,id,3,5)*a(i1f,k2f,jd,j,7,1,krf))
 
      b(i1c,i2c,id,j,2,2,irc) =     ((b(i1c,i2c,id,j,2,2,irc)
     &     + w(i1c,i2c,id,1,3)*qtt) + w(i1c,i2c,id,2,3)*qmt)
      b(i1c,i2c,id,j,2,3,irc) =     ((b(i1c,i2c,id,j,2,3,irc)
     &     + w(i1c,i2c,id,1,3)*qtb) + w(i1c,i2c,id,2,3)*qmb)
      b(i1c,i2c,id,j,2,1,jrc) =     ((b(i1c,i2c,id,j,2,1,jrc)
     &     + w(i1c,i2c,id,4,3)*qbt) + w(i1c,i2c,id,3,3)*qmt)
      b(i1c,i2c,id,j,2,2,jrc) =     ((b(i1c,i2c,id,j,2,2,jrc)
     &     + w(i1c,i2c,id,4,3)*qbb) + w(i1c,i2c,id,3,3)*qmb)
      b(i1c,k2c,id,j,7,2,irc) =     ((b(i1c,k2c,id,j,7,2,irc)
     &     + w(i1c,k2c,id,1,6)*qtt) + w(i1c,k2c,id,2,6)*qmt)
      b(i1c,k2c,id,j,7,3,irc) =     ((b(i1c,k2c,id,j,7,3,irc)
     &     + w(i1c,k2c,id,1,6)*qtb) + w(i1c,k2c,id,2,6)*qmb)
      b(i1c,k2c,id,j,7,1,jrc) =     ((b(i1c,k2c,id,j,7,1,jrc)
     &     + w(i1c,k2c,id,4,6)*qbt) + w(i1c,k2c,id,3,6)*qmt)
      b(i1c,k2c,id,j,7,2,jrc) =     ((b(i1c,k2c,id,j,7,2,jrc)
     &     + w(i1c,k2c,id,4,6)*qbb) + w(i1c,k2c,id,3,6)*qmb)
 
      rtt = ((((((((w(i1c,k2c,id,1,6)*a(i1f,k2f,jd,j,1,2,irf))
     &            + w(i1c,k2c,id,1,7)*a(i1f,k2f,jd,j,2,2,irf))
     &            +                   a(i1f,k2f,jd,j,3,2,irf))
     &            + w(i1c,k2c,id,1,5)*a(i1f,k2f,jd,j,4,2,irf))
     &            + w(i1c,k2c,id,2,6)*a(i1f,k2f,jd,j,1,3,irf))
     &            + w(i1c,k2c,id,2,7)*a(i1f,k2f,jd,j,2,3,irf))
     &            + w(i1c,k2c,id,2,1)*a(i1f,k2f,jd,j,3,3,irf))
     &            + w(i1c,k2c,id,2,5)*a(i1f,k2f,jd,j,4,3,irf))
      rmt = ((((((((w(i1c,k2c,id,1,6)*a(i1f,k2f,jd,j,1,1,jrf))
     &            + w(i1c,k2c,id,1,7)*a(i1f,k2f,jd,j,2,1,jrf))
     &            +                   a(i1f,k2f,jd,j,3,1,jrf))
     &            + w(i1c,k2c,id,1,5)*a(i1f,k2f,jd,j,4,1,jrf))
     &            + w(i1c,k2c,id,2,6)*a(i1f,k2f,jd,j,1,2,jrf))
     &            + w(i1c,k2c,id,2,7)*a(i1f,k2f,jd,j,2,2,jrf))
     &            + w(i1c,k2c,id,2,1)*a(i1f,k2f,jd,j,3,2,jrf))
     &            + w(i1c,k2c,id,2,5)*a(i1f,k2f,jd,j,4,2,jrf))
      rbt =     ((((w(i1c,k2c,id,2,6)*a(i1f,k2f,jd,j,1,1,krf))
     &            + w(i1c,k2c,id,2,7)*a(i1f,k2f,jd,j,2,1,krf))
     &            + w(i1c,k2c,id,2,1)*a(i1f,k2f,jd,j,3,1,krf))
     &            + w(i1c,k2c,id,2,5)*a(i1f,k2f,jd,j,4,1,krf))
      rtb =     ((((w(i1c,k2c,id,3,6)*a(i1f,k2f,jd,j,1,3,irf))
     &            + w(i1c,k2c,id,3,7)*a(i1f,k2f,jd,j,2,3,irf))
     &            + w(i1c,k2c,id,3,1)*a(i1f,k2f,jd,j,3,3,irf))
     &            + w(i1c,k2c,id,3,5)*a(i1f,k2f,jd,j,4,3,irf))
      rmb = ((((((((w(i1c,k2c,id,3,6)*a(i1f,k2f,jd,j,1,2,jrf))
     &            + w(i1c,k2c,id,3,7)*a(i1f,k2f,jd,j,2,2,jrf))
     &            + w(i1c,k2c,id,3,1)*a(i1f,k2f,jd,j,3,2,jrf))
     &            + w(i1c,k2c,id,3,5)*a(i1f,k2f,jd,j,4,2,jrf))
     &            + w(i1c,k2c,id,4,6)*a(i1f,k2f,jd,j,1,3,jrf))
     &            + w(i1c,k2c,id,4,7)*a(i1f,k2f,jd,j,2,3,jrf))
     &            +                   a(i1f,k2f,jd,j,3,3,jrf))
     &            + w(i1c,k2c,id,4,5)*a(i1f,k2f,jd,j,4,3,jrf))
      rbb =     ((((w(i1c,k2c,id,3,6)*a(i1f,k2f,jd,j,1,1,krf))
     &            + w(i1c,k2c,id,3,7)*a(i1f,k2f,jd,j,2,1,krf))
     &            + w(i1c,k2c,id,3,1)*a(i1f,k2f,jd,j,3,1,krf))
     &            + w(i1c,k2c,id,3,5)*a(i1f,k2f,jd,j,4,1,krf))
 
      b(i1c,i2c,id,j,3,2,irc) =     ((b(i1c,i2c,id,j,3,2,irc)
     &     + w(i1c,i2c,id,1,3)*rtt) + w(i1c,i2c,id,2,3)*rmt)
      b(i1c,i2c,id,j,3,3,irc) =     ((b(i1c,i2c,id,j,3,3,irc)
     &     + w(i1c,i2c,id,1,3)*rtb) + w(i1c,i2c,id,2,3)*rmb)
      b(i1c,i2c,id,j,3,1,jrc) =     ((b(i1c,i2c,id,j,3,1,jrc)
     &     + w(i1c,i2c,id,4,3)*rbt) + w(i1c,i2c,id,3,3)*rmt)
      b(i1c,i2c,id,j,3,2,jrc) =     ((b(i1c,i2c,id,j,3,2,jrc)
     &     + w(i1c,i2c,id,4,3)*rbb) + w(i1c,i2c,id,3,3)*rmb)
      b(i1c,k2c,id,j,1,2,irc) =     ((b(i1c,k2c,id,j,1,2,irc)
     &     + w(i1c,k2c,id,1,6)*rtt) + w(i1c,k2c,id,2,6)*rmt)
      b(i1c,k2c,id,j,1,3,irc) =     ((b(i1c,k2c,id,j,1,3,irc)
     &     + w(i1c,k2c,id,1,6)*rtb) + w(i1c,k2c,id,2,6)*rmb)
      b(i1c,k2c,id,j,1,1,jrc) =     ((b(i1c,k2c,id,j,1,1,jrc)
     &     + w(i1c,k2c,id,4,6)*rbt) + w(i1c,k2c,id,3,6)*rmt)
      b(i1c,k2c,id,j,1,2,jrc) =     ((b(i1c,k2c,id,j,1,2,jrc)
     &     + w(i1c,k2c,id,4,6)*rbb) + w(i1c,k2c,id,3,6)*rmb)
 
      stt = ((((w(j1c,k2c,id,1,2)*a(i1f,k2f,jd,j,4,2,irf))
     &        + w(j1c,k2c,id,1,7)*a(i1f,k2f,jd,j,5,2,irf))
     &        + w(j1c,k2c,id,2,2)*a(i1f,k2f,jd,j,4,3,irf))
     &        + w(j1c,k2c,id,2,7)*a(i1f,k2f,jd,j,5,3,irf))
      smt = ((((w(j1c,k2c,id,1,2)*a(i1f,k2f,jd,j,4,1,jrf))
     &        + w(j1c,k2c,id,1,7)*a(i1f,k2f,jd,j,5,1,jrf))
     &        + w(j1c,k2c,id,2,2)*a(i1f,k2f,jd,j,4,2,jrf))
     &        + w(j1c,k2c,id,2,7)*a(i1f,k2f,jd,j,5,2,jrf))
      sbt =   ((w(j1c,k2c,id,2,2)*a(i1f,k2f,jd,j,4,1,krf))
     &        + w(j1c,k2c,id,2,7)*a(i1f,k2f,jd,j,5,1,krf))
      stb =   ((w(j1c,k2c,id,3,2)*a(i1f,k2f,jd,j,4,3,irf))
     &        + w(j1c,k2c,id,3,7)*a(i1f,k2f,jd,j,5,3,irf))
      smb = ((((w(j1c,k2c,id,3,2)*a(i1f,k2f,jd,j,4,2,jrf))
     &        + w(j1c,k2c,id,3,7)*a(i1f,k2f,jd,j,5,2,jrf))
     &        + w(j1c,k2c,id,4,2)*a(i1f,k2f,jd,j,4,3,jrf))
     &        + w(j1c,k2c,id,4,7)*a(i1f,k2f,jd,j,5,3,jrf))
      sbb =   ((w(j1c,k2c,id,3,2)*a(i1f,k2f,jd,j,4,1,krf))
     &        + w(j1c,k2c,id,3,7)*a(i1f,k2f,jd,j,5,1,krf))
 
      b(i1c,i2c,id,j,4,2,irc) =     ((b(i1c,i2c,id,j,4,2,irc)
     &     + w(i1c,i2c,id,1,3)*stt) + w(i1c,i2c,id,2,3)*smt)
      b(i1c,i2c,id,j,4,3,irc) =     ((b(i1c,i2c,id,j,4,3,irc)
     &     + w(i1c,i2c,id,1,3)*stb) + w(i1c,i2c,id,2,3)*smb)
      b(i1c,i2c,id,j,4,1,jrc) =     ((b(i1c,i2c,id,j,4,1,jrc)
     &     + w(i1c,i2c,id,4,3)*sbt) + w(i1c,i2c,id,3,3)*smt)
      b(i1c,i2c,id,j,4,2,jrc) =     ((b(i1c,i2c,id,j,4,2,jrc)
     &     + w(i1c,i2c,id,4,3)*sbb) + w(i1c,i2c,id,3,3)*smb)
      b(i1c,k2c,id,j,5,2,irc) =     ((b(i1c,k2c,id,j,5,2,irc)
     &     + w(i1c,k2c,id,1,6)*stt) + w(i1c,k2c,id,2,6)*smt)
      b(i1c,k2c,id,j,5,3,irc) =     ((b(i1c,k2c,id,j,5,3,irc)
     &     + w(i1c,k2c,id,1,6)*stb) + w(i1c,k2c,id,2,6)*smb)
      b(i1c,k2c,id,j,5,1,jrc) =     ((b(i1c,k2c,id,j,5,1,jrc)
     &     + w(i1c,k2c,id,4,6)*sbt) + w(i1c,k2c,id,3,6)*smt)
      b(i1c,k2c,id,j,5,2,jrc) =     ((b(i1c,k2c,id,j,5,2,jrc)
     &     + w(i1c,k2c,id,4,6)*sbb) + w(i1c,k2c,id,3,6)*smb)
 
 50   continue
 
 60   continue
 
c     Add the unaccounted bottom layer contribution.
 
      irc = nrc+1
      irf = nrf+1
      if(ntf .eq. nt) irf = mod(irf-1, 3) + 1
 
      do 130 id=1,ndo
 
      if(ntf .eq. nt) call fineoprdiag(a,d,id,nrf-1)
 
      jd = id
      if(ntf .eq. nt) jd = 1
 
      do 130  j=1,9
 
      do 90 i2c=1,ntc+1
      i2f = i2c + i2c - 1
      do 90 i1c=1,ntc+1
      i1f = i1c + i1c - 1
 
      b(i1c,i2c,id,j,1,2,irc) =   b(i1c,i2c,id,j,1,2,irc)
     &                          + a(i1f,i2f,jd,j,1,2,irf)
 
      do 80   m=2,7
      j1c = i1c + j1n(m)
      j2c = i2c + j2n(m)
      b(i1c,i2c,id,j,1,2,irc) =   b(i1c,i2c,id,j,1,2,irc)
     &    + w(i1c,i2c,id,4,   m )*a(i1f,i2f,jd,j,m,2,irf)
      b(i1c,i2c,id,j,m,2,irc) =   b(i1c,i2c,id,j,m,2,irc)
     &    + w(j1c,j2c,id,4,md(m))*a(i1f,i2f,jd,j,m,2,irf)
 80   continue
 
 90   continue
 
      do 100 i2c=1,ntc+1
      j2c = i2c - 1
      k2c = i2c + 1
      i2f = i2c + i2c - 1
      do 100 i1c=2,ntc+1
      j1c = i1c - 1
      j1f = i1c + i1c - 2
      abb = ((((w(i1c,i2c,id,4,5)*a(j1f,i2f,jd,j,1,2,irf))
     &        +                   a(j1f,i2f,jd,j,2,2,irf))
     &        + w(i1c,i2c,id,4,4)*a(j1f,i2f,jd,j,3,2,irf))
     &        + w(i1c,i2c,id,4,6)*a(j1f,i2f,jd,j,7,2,irf))
      bbb =   ((w(j1c,k2c,id,4,7)*a(j1f,i2f,jd,j,3,2,irf))
     &        + w(j1c,k2c,id,4,6)*a(j1f,i2f,jd,j,4,2,irf))
      cbb = ((((w(j1c,i2c,id,4,2)*a(j1f,i2f,jd,j,1,2,irf))
     &        + w(j1c,i2c,id,4,3)*a(j1f,i2f,jd,j,4,2,irf))
     &        +                   a(j1f,i2f,jd,j,5,2,irf))
     &        + w(j1c,i2c,id,4,7)*a(j1f,i2f,jd,j,6,2,irf))
      dbb =   ((w(i1c,j2c,id,4,4)*a(j1f,i2f,jd,j,6,2,irf))
     &        + w(i1c,j2c,id,4,3)*a(j1f,i2f,jd,j,7,2,irf))
      b(i1c,i2c,id,j,1,2,irc) =   b(i1c,i2c,id,j,1,2,irc)
     &                          + w(i1c,i2c,id,4,5)*abb
      b(i1c,i2c,id,j,4,2,irc) =   b(i1c,i2c,id,j,4,2,irc)
     &                          + w(i1c,i2c,id,4,5)*bbb
      b(i1c,i2c,id,j,5,2,irc) =   b(i1c,i2c,id,j,5,2,irc)
     &                          + w(i1c,i2c,id,4,5)*cbb
      b(i1c,i2c,id,j,6,2,irc) =   b(i1c,i2c,id,j,6,2,irc)
     &                          + w(i1c,i2c,id,4,5)*dbb
      b(j1c,i2c,id,j,2,2,irc) =   b(j1c,i2c,id,j,2,2,irc)
     &                          + w(j1c,i2c,id,4,2)*abb
      b(j1c,i2c,id,j,3,2,irc) =   b(j1c,i2c,id,j,3,2,irc)
     &                          + w(j1c,i2c,id,4,2)*bbb
      b(j1c,i2c,id,j,1,2,irc) =   b(j1c,i2c,id,j,1,2,irc)
     &                          + w(j1c,i2c,id,4,2)*cbb
      b(j1c,i2c,id,j,7,2,irc) =   b(j1c,i2c,id,j,7,2,irc)
     &                          + w(j1c,i2c,id,4,2)*dbb
 100  continue
 
      do 110 i2c=1,ntc
      k2c = i2c + 1
      k2f = i2c + i2c
      do 110 i1c=2,ntc+1
      j1c = i1c - 1
      j1f = i1c + i1c - 2
      ebb = ((((w(i1c,i2c,id,4,4)*a(j1f,k2f,jd,j,1,2,irf))
     &        + w(i1c,i2c,id,4,3)*a(j1f,k2f,jd,j,2,2,irf))
     &        + w(i1c,i2c,id,4,5)*a(j1f,k2f,jd,j,6,2,irf))
     &        +                   a(j1f,k2f,jd,j,7,2,irf))
      fbb =   ((w(i1c,k2c,id,4,6)*a(j1f,k2f,jd,j,2,2,irf))
     &        + w(i1c,k2c,id,4,5)*a(j1f,k2f,jd,j,3,2,irf))
      gbb = ((((w(j1c,k2c,id,4,7)*a(j1f,k2f,jd,j,1,2,irf))
     &        + w(j1c,k2c,id,4,2)*a(j1f,k2f,jd,j,3,2,irf))
     &        +                   a(j1f,k2f,jd,j,4,2,irf))
     &        + w(j1c,k2c,id,4,6)*a(j1f,k2f,jd,j,5,2,irf))
      hbb =   ((w(j1c,i2c,id,4,3)*a(j1f,k2f,jd,j,5,2,irf))
     &        + w(j1c,i2c,id,4,2)*a(j1f,k2f,jd,j,6,2,irf))
      b(i1c,i2c,id,j,1,2,irc) =   b(i1c,i2c,id,j,1,2,irc)
     &                          + w(i1c,i2c,id,4,4)*ebb
      b(i1c,i2c,id,j,3,2,irc) =   b(i1c,i2c,id,j,3,2,irc)
     &                          + w(i1c,i2c,id,4,4)*fbb
      b(i1c,i2c,id,j,4,2,irc) =   b(i1c,i2c,id,j,4,2,irc)
     &                          + w(i1c,i2c,id,4,4)*gbb
      b(i1c,i2c,id,j,5,2,irc) =   b(i1c,i2c,id,j,5,2,irc)
     &                          + w(i1c,i2c,id,4,4)*hbb
      b(j1c,k2c,id,j,7,2,irc) =   b(j1c,k2c,id,j,7,2,irc)
     &                          + w(j1c,k2c,id,4,7)*ebb
      b(j1c,k2c,id,j,2,2,irc) =   b(j1c,k2c,id,j,2,2,irc)
     &                          + w(j1c,k2c,id,4,7)*fbb
      b(j1c,k2c,id,j,1,2,irc) =   b(j1c,k2c,id,j,1,2,irc)
     &                          + w(j1c,k2c,id,4,7)*gbb
      b(j1c,k2c,id,j,6,2,irc) =   b(j1c,k2c,id,j,6,2,irc)
     &                          + w(j1c,k2c,id,4,7)*hbb
 110  continue
 
      do 120 i2c=1,ntc
      k2c = i2c + 1
      k2f = i2c + i2c
      do 120 i1c=1,ntc+1
      j1c = i1c - 1
      k1c = i1c + 1
      i1f = i1c + i1c - 1
      pbb = ((((w(i1c,i2c,id,4,3)*a(i1f,k2f,jd,j,1,2,irf))
     &        + w(i1c,i2c,id,4,4)*a(i1f,k2f,jd,j,5,2,irf))
     &        +                   a(i1f,k2f,jd,j,6,2,irf))
     &        + w(i1c,i2c,id,4,2)*a(i1f,k2f,jd,j,7,2,irf))
      qbb =   ((w(k1c,i2c,id,4,4)*a(i1f,k2f,jd,j,2,2,irf))
     &        + w(k1c,i2c,id,4,5)*a(i1f,k2f,jd,j,7,2,irf))
      rbb = ((((w(i1c,k2c,id,4,6)*a(i1f,k2f,jd,j,1,2,irf))
     &        + w(i1c,k2c,id,4,7)*a(i1f,k2f,jd,j,2,2,irf))
     &        +                   a(i1f,k2f,jd,j,3,2,irf))
     &        + w(i1c,k2c,id,4,5)*a(i1f,k2f,jd,j,4,2,irf))
      sbb =   ((w(j1c,k2c,id,4,2)*a(i1f,k2f,jd,j,4,2,irf))
     &        + w(j1c,k2c,id,4,7)*a(i1f,k2f,jd,j,5,2,irf))
      b(i1c,i2c,id,j,1,2,irc) =   b(i1c,i2c,id,j,1,2,irc)
     &                          + w(i1c,i2c,id,4,3)*pbb
      b(i1c,i2c,id,j,2,2,irc) =   b(i1c,i2c,id,j,2,2,irc)
     &                          + w(i1c,i2c,id,4,3)*qbb
      b(i1c,i2c,id,j,3,2,irc) =   b(i1c,i2c,id,j,3,2,irc)
     &                          + w(i1c,i2c,id,4,3)*rbb
      b(i1c,i2c,id,j,4,2,irc) =   b(i1c,i2c,id,j,4,2,irc)
     &                          + w(i1c,i2c,id,4,3)*sbb
      b(i1c,k2c,id,j,6,2,irc) =   b(i1c,k2c,id,j,6,2,irc)
     &                          + w(i1c,k2c,id,4,6)*pbb
      b(i1c,k2c,id,j,7,2,irc) =   b(i1c,k2c,id,j,7,2,irc)
     &                          + w(i1c,k2c,id,4,6)*qbb
      b(i1c,k2c,id,j,1,2,irc) =   b(i1c,k2c,id,j,1,2,irc)
     &                          + w(i1c,k2c,id,4,6)*rbb
      b(i1c,k2c,id,j,5,2,irc) =   b(i1c,k2c,id,j,5,2,irc)
     &                          + w(i1c,k2c,id,4,6)*sbb
 120  continue
 
 130  continue
 
      end
*dk fineoprdiag
      subroutine fineoprdiag(a,d,id,ir)
 
c...  This routine loads the diagonal tensor elements from d into a for
c...  diamond id and the three radial layers beginning with layer ir.
 
      include 'size.h'
      real a((nt+1)**2,9,7,3,3), d((nt+1)**2,ndo,9,nr+1)
 
      do jr=0,2
         kr = mod(ir+jr-1, 3) + 1
         do j=1,9
            do ii=1,(nt+1)**2
               a(ii,j,1,2,kr) = d(ii,id,j,ir+jr)
            end do
         end do
      end do
 
      end
*dk fineoprlayer
      subroutine fineoprlayer(a,ir)
 
c...  This routine assembles the tensor operator a (without its
c...  diagonal elements) at the finest grid level for radial layers
c...  ir and ir+1.
 
      include 'size.h'
      real a((nt+1)**2,9,7,3,3)
      common /oper/ ra(3,8,nr+1), alp(7,(nt+1)**2,2),
     &              atn(4,7,(nt+1)**2,9)
 
      do kr=ir,ir+1
 
         kb = 1
         ke = 3
         if(kr .eq. 1)    kb = 2
         if(kr .eq. nr+1) ke = 2
         jr = mod(kr-1, 3) + 1
 
         do k=kb,ke
 
            m1 = 1
            if(k .eq. 2) m1 = 2
 
            do m=m1,7
               do j=1,9
                  do ii=1,(nt+1)**2
                     a(ii,j,m,k,jr) = ((((ra(k,1,kr)*atn(1,m,ii,j))
     &                                  + ra(k,2,kr)*atn(2,m,ii,j))
     &                                  + ra(k,3,kr)*atn(3,m,ii,j))
     &                                  + ra(k,4,kr)*atn(4,m,ii,j))
                  end do
               end do
            end do
 
         end do
 
      end do
 
      end
*dk weights
      subroutine weights(w,rw,hw,ir,kt)
 
      include 'size.h'
      real hw(2:7,(kt+1)**2*ndo,*), rw(2,(2*kt+1)**2*ndo,*)
      real w((kt+1)**2*ndo,4,7)
      common /ofst/ j1n(7), j2n(7), md(7)
 
      jr   = ir + 1
      ntcp = (  kt+1)
      ntfp = (2*kt+1)
      nnc  = (  kt+1)**2
      nnf  = (2*kt+1)**2
 
      m2 = j1n(2) + (2*kt+1)*j2n(2)
      m3 = j1n(3) + (2*kt+1)*j2n(3)
      m4 = j1n(4) + (2*kt+1)*j2n(4)
      m5 = j1n(5) + (2*kt+1)*j2n(5)
      m6 = j1n(6) + (2*kt+1)*j2n(6)
      m7 = j1n(7) + (2*kt+1)*j2n(7)
 
      do ii=1,nnc*ndo
 
         ij = mod(ii-1,nnc)
         jj = 2*((ij/ntcp)*ntfp + mod(ij,ntcp))
     &      + ((ii-1)/nnc)*nnf  + 1
 
         w(ii,1,1) = 1.0
         w(ii,1,2) = hw(2,ii,ir)
         w(ii,1,3) = hw(3,ii,ir)
         w(ii,1,4) = hw(4,ii,ir)
         w(ii,1,5) = hw(5,ii,ir)
         w(ii,1,6) = hw(6,ii,ir)
         w(ii,1,7) = hw(7,ii,ir)
         w(ii,2,1) =             rw(2,jj   ,ir)
         w(ii,2,2) = hw(2,ii,ir)*rw(2,jj+m2,ir)
         w(ii,2,3) = hw(3,ii,ir)*rw(2,jj+m3,ir)
         w(ii,2,4) = hw(4,ii,ir)*rw(2,jj+m4,ir)
         w(ii,2,5) = hw(5,ii,ir)*rw(2,jj+m5,ir)
         w(ii,2,6) = hw(6,ii,ir)*rw(2,jj+m6,ir)
         w(ii,2,7) = hw(7,ii,ir)*rw(2,jj+m7,ir)
         w(ii,3,1) =             rw(1,jj   ,jr)
         w(ii,3,2) = hw(2,ii,jr)*rw(1,jj+m2,jr)
         w(ii,3,3) = hw(3,ii,jr)*rw(1,jj+m3,jr)
         w(ii,3,4) = hw(4,ii,jr)*rw(1,jj+m4,jr)
         w(ii,3,5) = hw(5,ii,jr)*rw(1,jj+m5,jr)
         w(ii,3,6) = hw(6,ii,jr)*rw(1,jj+m6,jr)
         w(ii,3,7) = hw(7,ii,jr)*rw(1,jj+m7,jr)
         w(ii,4,1) = 1.0
         w(ii,4,2) = hw(2,ii,jr)
         w(ii,4,3) = hw(3,ii,jr)
         w(ii,4,4) = hw(4,ii,jr)
         w(ii,4,5) = hw(5,ii,jr)
         w(ii,4,6) = hw(6,ii,jr)
         w(ii,4,7) = hw(7,ii,jr)
 
      end do
 
      end
*dk zerolvopr
      subroutine zerolvopr(a,b,hw)
 
c     This routine assembles coarse grid forward operators.
 
      include 'size.h'
      real a(3,3,ndo*9,7,3,2), b(2,2,ndo*9,7,3,2), hw(2:7,2,2,ndo,2)
      common /ofst/ j1n(7), j2n(7), md(7)
 
      call nulvec(b,1512*ndo)
 
c     Coarse from coarse.
 
      do 10 i2c=1,2
      i2f = i2c + i2c - 1
      do 10 i1c=1,2
      i1f = i1c + i1c - 1
 
      do 10  id=1,ndo*9
      b(i1c,i2c,id,1,2,1) = b(i1c,i2c,id,1,2,1) + a(i1f,i2f,id,1,2,1)
      b(i1c,i2c,id,1,3,1) = b(i1c,i2c,id,1,3,1) + a(i1f,i2f,id,1,3,1)
      b(i1c,i2c,id,1,1,2) = b(i1c,i2c,id,1,1,2) + a(i1f,i2f,id,1,1,2)
      b(i1c,i2c,id,1,2,2) = b(i1c,i2c,id,1,2,2) + a(i1f,i2f,id,1,2,2)
 10   continue
 
      do 20   m=2,7
      do 20 i2c=1,2
      i2f = i2c + i2c - 1
      j2c = i2c + j2n(m)
      do 20 i1c=1,2
      i1f = i1c + i1c - 1
      j1c = i1c + j1n(m)
 
      do 20  id=1,ndo*9
      jd = min(ndo, mod(id-1,nd) + 1)
 
      b(i1c,i2c,id,1,2,1) = b(i1c,i2c,id,1,2,1)
     &                    + a(i1f,i2f,id,m,2,1)*hw(   m ,i1c,i2c,jd,1)
      b(i1c,i2c,id,1,3,1) = b(i1c,i2c,id,1,3,1)
     &                    + a(i1f,i2f,id,m,3,1)*hw(   m ,i1c,i2c,jd,2)
      b(i1c,i2c,id,1,1,2) = b(i1c,i2c,id,1,1,2)
     &                    + a(i1f,i2f,id,m,1,2)*hw(   m ,i1c,i2c,jd,1)
      b(i1c,i2c,id,1,2,2) = b(i1c,i2c,id,1,2,2)
     &                    + a(i1f,i2f,id,m,2,2)*hw(   m ,i1c,i2c,jd,2)
      b(i1c,i2c,id,m,2,1) = b(i1c,i2c,id,m,2,1)
     &                    + a(i1f,i2f,id,m,2,1)*hw(md(m),j1c,j2c,jd,1)
      b(i1c,i2c,id,m,3,1) = b(i1c,i2c,id,m,3,1)
     &                    + a(i1f,i2f,id,m,3,1)*hw(md(m),j1c,j2c,jd,2)
      b(i1c,i2c,id,m,1,2) = b(i1c,i2c,id,m,1,2)
     &                    + a(i1f,i2f,id,m,1,2)*hw(md(m),j1c,j2c,jd,1)
      b(i1c,i2c,id,m,2,2) = b(i1c,i2c,id,m,2,2)
     &                    + a(i1f,i2f,id,m,2,2)*hw(md(m),j1c,j2c,jd,2)
 
 20   continue
 
c     Project fine node values to coarse node values.
 
      i1c = 2
      j1c = 1
      i1f = 3
      j1f = 2
 
      do 30 i2c=1,2
      j2c = i2c - 1
      k2c = i2c + 1
      i2f = i2c + i2c - 1
      do 30 id=1,ndo*9
      jd = min(ndo, mod(id-1,nd) + 1)
 
      att = ((((hw(5,i1c,i2c,jd,1)*a(j1f,i2f,id,1,2,1))
     &        +                    a(j1f,i2f,id,2,2,1))
     &        + hw(4,i1c,i2c,jd,1)*a(j1f,i2f,id,3,2,1))
     &        + hw(6,i1c,i2c,jd,1)*a(j1f,i2f,id,7,2,1))
      atb = ((((hw(5,i1c,i2c,jd,2)*a(j1f,i2f,id,1,3,1))
     &        +                    a(j1f,i2f,id,2,3,1))
     &        + hw(4,i1c,i2c,jd,2)*a(j1f,i2f,id,3,3,1))
     &        + hw(6,i1c,i2c,jd,2)*a(j1f,i2f,id,7,3,1))
      abt = ((((hw(5,i1c,i2c,jd,1)*a(j1f,i2f,id,1,1,2))
     &        +                    a(j1f,i2f,id,2,1,2))
     &        + hw(4,i1c,i2c,jd,1)*a(j1f,i2f,id,3,1,2))
     &        + hw(6,i1c,i2c,jd,1)*a(j1f,i2f,id,7,1,2))
      abb = ((((hw(5,i1c,i2c,jd,2)*a(j1f,i2f,id,1,2,2))
     &        +                    a(j1f,i2f,id,2,2,2))
     &        + hw(4,i1c,i2c,jd,2)*a(j1f,i2f,id,3,2,2))
     &        + hw(6,i1c,i2c,jd,2)*a(j1f,i2f,id,7,2,2))
 
      b(i1c,i2c,id,1,2,1) = b(i1c,i2c,id,1,2,1) + att*hw(5,i1c,i2c,jd,1)
      b(i1c,i2c,id,1,3,1) = b(i1c,i2c,id,1,3,1) + atb*hw(5,i1c,i2c,jd,1)
      b(i1c,i2c,id,1,1,2) = b(i1c,i2c,id,1,1,2) + abt*hw(5,i1c,i2c,jd,2)
      b(i1c,i2c,id,1,2,2) = b(i1c,i2c,id,1,2,2) + abb*hw(5,i1c,i2c,jd,2)
      b(j1c,i2c,id,2,2,1) = b(j1c,i2c,id,2,2,1) + att*hw(2,j1c,i2c,jd,1)
      b(j1c,i2c,id,2,3,1) = b(j1c,i2c,id,2,3,1) + atb*hw(2,j1c,i2c,jd,1)
      b(j1c,i2c,id,2,1,2) = b(j1c,i2c,id,2,1,2) + abt*hw(2,j1c,i2c,jd,2)
      b(j1c,i2c,id,2,2,2) = b(j1c,i2c,id,2,2,2) + abb*hw(2,j1c,i2c,jd,2)
 
      btt = ((hw(7,j1c,k2c,jd,1)*a(j1f,i2f,id,3,2,1))
     &      + hw(6,j1c,k2c,jd,1)*a(j1f,i2f,id,4,2,1))
      btb = ((hw(7,j1c,k2c,jd,2)*a(j1f,i2f,id,3,3,1))
     &      + hw(6,j1c,k2c,jd,2)*a(j1f,i2f,id,4,3,1))
      bbt = ((hw(7,j1c,k2c,jd,1)*a(j1f,i2f,id,3,1,2))
     &      + hw(6,j1c,k2c,jd,1)*a(j1f,i2f,id,4,1,2))
      bbb = ((hw(7,j1c,k2c,jd,2)*a(j1f,i2f,id,3,2,2))
     &      + hw(6,j1c,k2c,jd,2)*a(j1f,i2f,id,4,2,2))
 
      b(i1c,i2c,id,4,2,1) = b(i1c,i2c,id,4,2,1) + btt*hw(5,i1c,i2c,jd,1)
      b(i1c,i2c,id,4,3,1) = b(i1c,i2c,id,4,3,1) + btb*hw(5,i1c,i2c,jd,1)
      b(i1c,i2c,id,4,1,2) = b(i1c,i2c,id,4,1,2) + bbt*hw(5,i1c,i2c,jd,2)
      b(i1c,i2c,id,4,2,2) = b(i1c,i2c,id,4,2,2) + bbb*hw(5,i1c,i2c,jd,2)
      b(j1c,i2c,id,3,2,1) = b(j1c,i2c,id,3,2,1) + btt*hw(2,j1c,i2c,jd,1)
      b(j1c,i2c,id,3,3,1) = b(j1c,i2c,id,3,3,1) + btb*hw(2,j1c,i2c,jd,1)
      b(j1c,i2c,id,3,1,2) = b(j1c,i2c,id,3,1,2) + bbt*hw(2,j1c,i2c,jd,2)
      b(j1c,i2c,id,3,2,2) = b(j1c,i2c,id,3,2,2) + bbb*hw(2,j1c,i2c,jd,2)
 
      ctt = ((((hw(2,j1c,i2c,jd,1)*a(j1f,i2f,id,1,2,1))
     &        + hw(3,j1c,i2c,jd,1)*a(j1f,i2f,id,4,2,1))
     &        +                    a(j1f,i2f,id,5,2,1))
     &        + hw(7,j1c,i2c,jd,1)*a(j1f,i2f,id,6,2,1))
      ctb = ((((hw(2,j1c,i2c,jd,2)*a(j1f,i2f,id,1,3,1))
     &        + hw(3,j1c,i2c,jd,2)*a(j1f,i2f,id,4,3,1))
     &        +                    a(j1f,i2f,id,5,3,1))
     &        + hw(7,j1c,i2c,jd,2)*a(j1f,i2f,id,6,3,1))
      cbt = ((((hw(2,j1c,i2c,jd,1)*a(j1f,i2f,id,1,1,2))
     &        + hw(3,j1c,i2c,jd,1)*a(j1f,i2f,id,4,1,2))
     &        +                    a(j1f,i2f,id,5,1,2))
     &        + hw(7,j1c,i2c,jd,1)*a(j1f,i2f,id,6,1,2))
      cbb = ((((hw(2,j1c,i2c,jd,2)*a(j1f,i2f,id,1,2,2))
     &        + hw(3,j1c,i2c,jd,2)*a(j1f,i2f,id,4,2,2))
     &        +                    a(j1f,i2f,id,5,2,2))
     &        + hw(7,j1c,i2c,jd,2)*a(j1f,i2f,id,6,2,2))
 
      b(i1c,i2c,id,5,2,1) = b(i1c,i2c,id,5,2,1) + ctt*hw(5,i1c,i2c,jd,1)
      b(i1c,i2c,id,5,3,1) = b(i1c,i2c,id,5,3,1) + ctb*hw(5,i1c,i2c,jd,1)
      b(i1c,i2c,id,5,1,2) = b(i1c,i2c,id,5,1,2) + cbt*hw(5,i1c,i2c,jd,2)
      b(i1c,i2c,id,5,2,2) = b(i1c,i2c,id,5,2,2) + cbb*hw(5,i1c,i2c,jd,2)
      b(j1c,i2c,id,1,2,1) = b(j1c,i2c,id,1,2,1) + ctt*hw(2,j1c,i2c,jd,1)
      b(j1c,i2c,id,1,3,1) = b(j1c,i2c,id,1,3,1) + ctb*hw(2,j1c,i2c,jd,1)
      b(j1c,i2c,id,1,1,2) = b(j1c,i2c,id,1,1,2) + cbt*hw(2,j1c,i2c,jd,2)
      b(j1c,i2c,id,1,2,2) = b(j1c,i2c,id,1,2,2) + cbb*hw(2,j1c,i2c,jd,2)
 
      dtt = ((hw(4,i1c,j2c,jd,1)*a(j1f,i2f,id,6,2,1))
     &      + hw(3,i1c,j2c,jd,1)*a(j1f,i2f,id,7,2,1))
      dtb = ((hw(4,i1c,j2c,jd,2)*a(j1f,i2f,id,6,3,1))
     &      + hw(3,i1c,j2c,jd,2)*a(j1f,i2f,id,7,3,1))
      dbt = ((hw(4,i1c,j2c,jd,1)*a(j1f,i2f,id,6,1,2))
     &      + hw(3,i1c,j2c,jd,1)*a(j1f,i2f,id,7,1,2))
      dbb = ((hw(4,i1c,j2c,jd,2)*a(j1f,i2f,id,6,2,2))
     &      + hw(3,i1c,j2c,jd,2)*a(j1f,i2f,id,7,2,2))
 
      b(i1c,i2c,id,6,2,1) = b(i1c,i2c,id,6,2,1) + dtt*hw(5,i1c,i2c,jd,1)
      b(i1c,i2c,id,6,3,1) = b(i1c,i2c,id,6,3,1) + dtb*hw(5,i1c,i2c,jd,1)
      b(i1c,i2c,id,6,1,2) = b(i1c,i2c,id,6,1,2) + dbt*hw(5,i1c,i2c,jd,2)
      b(i1c,i2c,id,6,2,2) = b(i1c,i2c,id,6,2,2) + dbb*hw(5,i1c,i2c,jd,2)
      b(j1c,i2c,id,7,2,1) = b(j1c,i2c,id,7,2,1) + dtt*hw(2,j1c,i2c,jd,1)
      b(j1c,i2c,id,7,3,1) = b(j1c,i2c,id,7,3,1) + dtb*hw(2,j1c,i2c,jd,1)
      b(j1c,i2c,id,7,1,2) = b(j1c,i2c,id,7,1,2) + dbt*hw(2,j1c,i2c,jd,2)
      b(j1c,i2c,id,7,2,2) = b(j1c,i2c,id,7,2,2) + dbb*hw(2,j1c,i2c,jd,2)
 
 30   continue
 
      i1c = 2
      j1c = 1
      i1f = 3
      j1f = 2
 
      i2c = 1
      k2c = 2
      i2f = 1
      k2f = 2
 
      do 40 id=1,ndo*9
      jd = min(ndo, mod(id-1,nd) + 1)
 
      ett = ((((hw(4,i1c,i2c,jd,1)*a(j1f,k2f,id,1,2,1))
     &        + hw(3,i1c,i2c,jd,1)*a(j1f,k2f,id,2,2,1))
     &        + hw(5,i1c,i2c,jd,1)*a(j1f,k2f,id,6,2,1))
     &        +                    a(j1f,k2f,id,7,2,1))
      etb = ((((hw(4,i1c,i2c,jd,2)*a(j1f,k2f,id,1,3,1))
     &        + hw(3,i1c,i2c,jd,2)*a(j1f,k2f,id,2,3,1))
     &        + hw(5,i1c,i2c,jd,2)*a(j1f,k2f,id,6,3,1))
     &        +                    a(j1f,k2f,id,7,3,1))
      ebt = ((((hw(4,i1c,i2c,jd,1)*a(j1f,k2f,id,1,1,2))
     &        + hw(3,i1c,i2c,jd,1)*a(j1f,k2f,id,2,1,2))
     &        + hw(5,i1c,i2c,jd,1)*a(j1f,k2f,id,6,1,2))
     &        +                    a(j1f,k2f,id,7,1,2))
      ebb = ((((hw(4,i1c,i2c,jd,2)*a(j1f,k2f,id,1,2,2))
     &        + hw(3,i1c,i2c,jd,2)*a(j1f,k2f,id,2,2,2))
     &        + hw(5,i1c,i2c,jd,2)*a(j1f,k2f,id,6,2,2))
     &        +                    a(j1f,k2f,id,7,2,2))
 
      b(i1c,i2c,id,1,2,1) = b(i1c,i2c,id,1,2,1) + ett*hw(4,i1c,i2c,jd,1)
      b(i1c,i2c,id,1,3,1) = b(i1c,i2c,id,1,3,1) + etb*hw(4,i1c,i2c,jd,1)
      b(i1c,i2c,id,1,1,2) = b(i1c,i2c,id,1,1,2) + ebt*hw(4,i1c,i2c,jd,2)
      b(i1c,i2c,id,1,2,2) = b(i1c,i2c,id,1,2,2) + ebb*hw(4,i1c,i2c,jd,2)
      b(j1c,k2c,id,7,2,1) = b(j1c,k2c,id,7,2,1) + ett*hw(7,j1c,k2c,jd,1)
      b(j1c,k2c,id,7,3,1) = b(j1c,k2c,id,7,3,1) + etb*hw(7,j1c,k2c,jd,1)
      b(j1c,k2c,id,7,1,2) = b(j1c,k2c,id,7,1,2) + ebt*hw(7,j1c,k2c,jd,2)
      b(j1c,k2c,id,7,2,2) = b(j1c,k2c,id,7,2,2) + ebb*hw(7,j1c,k2c,jd,2)
 
      ftt = ((hw(6,i1c,k2c,jd,1)*a(j1f,k2f,id,2,2,1))
     &      + hw(5,i1c,k2c,jd,1)*a(j1f,k2f,id,3,2,1))
      ftb = ((hw(6,i1c,k2c,jd,2)*a(j1f,k2f,id,2,3,1))
     &      + hw(5,i1c,k2c,jd,2)*a(j1f,k2f,id,3,3,1))
      fbt = ((hw(6,i1c,k2c,jd,1)*a(j1f,k2f,id,2,1,2))
     &      + hw(5,i1c,k2c,jd,1)*a(j1f,k2f,id,3,1,2))
      fbb = ((hw(6,i1c,k2c,jd,2)*a(j1f,k2f,id,2,2,2))
     &      + hw(5,i1c,k2c,jd,2)*a(j1f,k2f,id,3,2,2))
 
      b(i1c,i2c,id,3,2,1) = b(i1c,i2c,id,3,2,1) + ftt*hw(4,i1c,i2c,jd,1)
      b(i1c,i2c,id,3,3,1) = b(i1c,i2c,id,3,3,1) + ftb*hw(4,i1c,i2c,jd,1)
      b(i1c,i2c,id,3,1,2) = b(i1c,i2c,id,3,1,2) + fbt*hw(4,i1c,i2c,jd,2)
      b(i1c,i2c,id,3,2,2) = b(i1c,i2c,id,3,2,2) + fbb*hw(4,i1c,i2c,jd,2)
      b(j1c,k2c,id,2,2,1) = b(j1c,k2c,id,2,2,1) + ftt*hw(7,j1c,k2c,jd,1)
      b(j1c,k2c,id,2,3,1) = b(j1c,k2c,id,2,3,1) + ftb*hw(7,j1c,k2c,jd,1)
      b(j1c,k2c,id,2,1,2) = b(j1c,k2c,id,2,1,2) + fbt*hw(7,j1c,k2c,jd,2)
      b(j1c,k2c,id,2,2,2) = b(j1c,k2c,id,2,2,2) + fbb*hw(7,j1c,k2c,jd,2)
 
      gtt = ((((hw(7,j1c,k2c,jd,1)*a(j1f,k2f,id,1,2,1))
     &        + hw(2,j1c,k2c,jd,1)*a(j1f,k2f,id,3,2,1))
     &        +                    a(j1f,k2f,id,4,2,1))
     &        + hw(6,j1c,k2c,jd,1)*a(j1f,k2f,id,5,2,1))
      gtb = ((((hw(7,j1c,k2c,jd,2)*a(j1f,k2f,id,1,3,1))
     &        + hw(2,j1c,k2c,jd,2)*a(j1f,k2f,id,3,3,1))
     &        +                    a(j1f,k2f,id,4,3,1))
     &        + hw(6,j1c,k2c,jd,2)*a(j1f,k2f,id,5,3,1))
      gbt = ((((hw(7,j1c,k2c,jd,1)*a(j1f,k2f,id,1,1,2))
     &        + hw(2,j1c,k2c,jd,1)*a(j1f,k2f,id,3,1,2))
     &        +                    a(j1f,k2f,id,4,1,2))
     &        + hw(6,j1c,k2c,jd,1)*a(j1f,k2f,id,5,1,2))
      gbb = ((((hw(7,j1c,k2c,jd,2)*a(j1f,k2f,id,1,2,2))
     &        + hw(2,j1c,k2c,jd,2)*a(j1f,k2f,id,3,2,2))
     &        +                    a(j1f,k2f,id,4,2,2))
     &        + hw(6,j1c,k2c,jd,2)*a(j1f,k2f,id,5,2,2))
 
      b(i1c,i2c,id,4,2,1) = b(i1c,i2c,id,4,2,1) + gtt*hw(4,i1c,i2c,jd,1)
      b(i1c,i2c,id,4,3,1) = b(i1c,i2c,id,4,3,1) + gtb*hw(4,i1c,i2c,jd,1)
      b(i1c,i2c,id,4,1,2) = b(i1c,i2c,id,4,1,2) + gbt*hw(4,i1c,i2c,jd,2)
      b(i1c,i2c,id,4,2,2) = b(i1c,i2c,id,4,2,2) + gbb*hw(4,i1c,i2c,jd,2)
      b(j1c,k2c,id,1,2,1) = b(j1c,k2c,id,1,2,1) + gtt*hw(7,j1c,k2c,jd,1)
      b(j1c,k2c,id,1,3,1) = b(j1c,k2c,id,1,3,1) + gtb*hw(7,j1c,k2c,jd,1)
      b(j1c,k2c,id,1,1,2) = b(j1c,k2c,id,1,1,2) + gbt*hw(7,j1c,k2c,jd,2)
      b(j1c,k2c,id,1,2,2) = b(j1c,k2c,id,1,2,2) + gbb*hw(7,j1c,k2c,jd,2)
 
      htt = ((hw(3,j1c,i2c,jd,1)*a(j1f,k2f,id,5,2,1))
     &      + hw(2,j1c,i2c,jd,1)*a(j1f,k2f,id,6,2,1))
      htb = ((hw(3,j1c,i2c,jd,2)*a(j1f,k2f,id,5,3,1))
     &      + hw(2,j1c,i2c,jd,2)*a(j1f,k2f,id,6,3,1))
      hbt = ((hw(3,j1c,i2c,jd,1)*a(j1f,k2f,id,5,1,2))
     &      + hw(2,j1c,i2c,jd,1)*a(j1f,k2f,id,6,1,2))
      hbb = ((hw(3,j1c,i2c,jd,2)*a(j1f,k2f,id,5,2,2))
     &      + hw(2,j1c,i2c,jd,2)*a(j1f,k2f,id,6,2,2))
 
      b(i1c,i2c,id,5,2,1) = b(i1c,i2c,id,5,2,1) + htt*hw(4,i1c,i2c,jd,1)
      b(i1c,i2c,id,5,3,1) = b(i1c,i2c,id,5,3,1) + htb*hw(4,i1c,i2c,jd,1)
      b(i1c,i2c,id,5,1,2) = b(i1c,i2c,id,5,1,2) + hbt*hw(4,i1c,i2c,jd,2)
      b(i1c,i2c,id,5,2,2) = b(i1c,i2c,id,5,2,2) + hbb*hw(4,i1c,i2c,jd,2)
      b(j1c,k2c,id,6,2,1) = b(j1c,k2c,id,6,2,1) + htt*hw(7,j1c,k2c,jd,1)
      b(j1c,k2c,id,6,3,1) = b(j1c,k2c,id,6,3,1) + htb*hw(7,j1c,k2c,jd,1)
      b(j1c,k2c,id,6,1,2) = b(j1c,k2c,id,6,1,2) + hbt*hw(7,j1c,k2c,jd,2)
      b(j1c,k2c,id,6,2,2) = b(j1c,k2c,id,6,2,2) + hbb*hw(7,j1c,k2c,jd,2)
 
 40   continue
 
      i2c = 1
      k2c = 2
      i2f = 1
      k2f = 2
 
      do 50 i1c=1,2
      j1c = i1c - 1
      k1c = i1c + 1
      i1f = i1c + i1c - 1
      do 50 id=1,ndo*9
      jd = min(ndo, mod(id-1,nd) + 1)
 
      ptt = ((((hw(3,i1c,i2c,jd,1)*a(i1f,k2f,id,1,2,1))
     &        + hw(4,i1c,i2c,jd,1)*a(i1f,k2f,id,5,2,1))
     &        +                    a(i1f,k2f,id,6,2,1))
     &        + hw(2,i1c,i2c,jd,1)*a(i1f,k2f,id,7,2,1))
      ptb = ((((hw(3,i1c,i2c,jd,2)*a(i1f,k2f,id,1,3,1))
     &        + hw(4,i1c,i2c,jd,2)*a(i1f,k2f,id,5,3,1))
     &        +                    a(i1f,k2f,id,6,3,1))
     &        + hw(2,i1c,i2c,jd,2)*a(i1f,k2f,id,7,3,1))
      pbt = ((((hw(3,i1c,i2c,jd,1)*a(i1f,k2f,id,1,1,2))
     &        + hw(4,i1c,i2c,jd,1)*a(i1f,k2f,id,5,1,2))
     &        +                    a(i1f,k2f,id,6,1,2))
     &        + hw(2,i1c,i2c,jd,1)*a(i1f,k2f,id,7,1,2))
      pbb = ((((hw(3,i1c,i2c,jd,2)*a(i1f,k2f,id,1,2,2))
     &        + hw(4,i1c,i2c,jd,2)*a(i1f,k2f,id,5,2,2))
     &        +                    a(i1f,k2f,id,6,2,2))
     &        + hw(2,i1c,i2c,jd,2)*a(i1f,k2f,id,7,2,2))
 
      b(i1c,i2c,id,1,2,1) = b(i1c,i2c,id,1,2,1) + ptt*hw(3,i1c,i2c,jd,1)
      b(i1c,i2c,id,1,3,1) = b(i1c,i2c,id,1,3,1) + ptb*hw(3,i1c,i2c,jd,1)
      b(i1c,i2c,id,1,1,2) = b(i1c,i2c,id,1,1,2) + pbt*hw(3,i1c,i2c,jd,2)
      b(i1c,i2c,id,1,2,2) = b(i1c,i2c,id,1,2,2) + pbb*hw(3,i1c,i2c,jd,2)
      b(i1c,k2c,id,6,2,1) = b(i1c,k2c,id,6,2,1) + ptt*hw(6,i1c,k2c,jd,1)
      b(i1c,k2c,id,6,3,1) = b(i1c,k2c,id,6,3,1) + ptb*hw(6,i1c,k2c,jd,1)
      b(i1c,k2c,id,6,1,2) = b(i1c,k2c,id,6,1,2) + pbt*hw(6,i1c,k2c,jd,2)
      b(i1c,k2c,id,6,2,2) = b(i1c,k2c,id,6,2,2) + pbb*hw(6,i1c,k2c,jd,2)
 
      qtt = ((hw(4,k1c,i2c,jd,1)*a(i1f,k2f,id,2,2,1))
     &      + hw(5,k1c,i2c,jd,1)*a(i1f,k2f,id,7,2,1))
      qtb = ((hw(4,k1c,i2c,jd,2)*a(i1f,k2f,id,2,3,1))
     &      + hw(5,k1c,i2c,jd,2)*a(i1f,k2f,id,7,3,1))
      qbt = ((hw(4,k1c,i2c,jd,1)*a(i1f,k2f,id,2,1,2))
     &      + hw(5,k1c,i2c,jd,1)*a(i1f,k2f,id,7,1,2))
      qbb = ((hw(4,k1c,i2c,jd,2)*a(i1f,k2f,id,2,2,2))
     &      + hw(5,k1c,i2c,jd,2)*a(i1f,k2f,id,7,2,2))
 
      b(i1c,i2c,id,2,2,1) = b(i1c,i2c,id,2,2,1) + qtt*hw(3,i1c,i2c,jd,1)
      b(i1c,i2c,id,2,3,1) = b(i1c,i2c,id,2,3,1) + qtb*hw(3,i1c,i2c,jd,1)
      b(i1c,i2c,id,2,1,2) = b(i1c,i2c,id,2,1,2) + qbt*hw(3,i1c,i2c,jd,2)
      b(i1c,i2c,id,2,2,2) = b(i1c,i2c,id,2,2,2) + qbb*hw(3,i1c,i2c,jd,2)
      b(i1c,k2c,id,7,2,1) = b(i1c,k2c,id,7,2,1) + qtt*hw(6,i1c,k2c,jd,1)
      b(i1c,k2c,id,7,3,1) = b(i1c,k2c,id,7,3,1) + qtb*hw(6,i1c,k2c,jd,1)
      b(i1c,k2c,id,7,1,2) = b(i1c,k2c,id,7,1,2) + qbt*hw(6,i1c,k2c,jd,2)
      b(i1c,k2c,id,7,2,2) = b(i1c,k2c,id,7,2,2) + qbb*hw(6,i1c,k2c,jd,2)
 
      rtt = ((((hw(6,i1c,k2c,jd,1)*a(i1f,k2f,id,1,2,1))
     &        + hw(7,i1c,k2c,jd,1)*a(i1f,k2f,id,2,2,1))
     &        +                    a(i1f,k2f,id,3,2,1))
     &        + hw(5,i1c,k2c,jd,1)*a(i1f,k2f,id,4,2,1))
      rtb = ((((hw(6,i1c,k2c,jd,2)*a(i1f,k2f,id,1,3,1))
     &        + hw(7,i1c,k2c,jd,2)*a(i1f,k2f,id,2,3,1))
     &        +                    a(i1f,k2f,id,3,3,1))
     &        + hw(5,i1c,k2c,jd,2)*a(i1f,k2f,id,4,3,1))
      rbt = ((((hw(6,i1c,k2c,jd,1)*a(i1f,k2f,id,1,1,2))
     &        + hw(7,i1c,k2c,jd,1)*a(i1f,k2f,id,2,1,2))
     &        +                    a(i1f,k2f,id,3,1,2))
     &        + hw(5,i1c,k2c,jd,1)*a(i1f,k2f,id,4,1,2))
      rbb = ((((hw(6,i1c,k2c,jd,2)*a(i1f,k2f,id,1,2,2))
     &        + hw(7,i1c,k2c,jd,2)*a(i1f,k2f,id,2,2,2))
     &        +                    a(i1f,k2f,id,3,2,2))
     &        + hw(5,i1c,k2c,jd,2)*a(i1f,k2f,id,4,2,2))
 
      b(i1c,i2c,id,3,2,1) = b(i1c,i2c,id,3,2,1) + rtt*hw(3,i1c,i2c,jd,1)
      b(i1c,i2c,id,3,3,1) = b(i1c,i2c,id,3,3,1) + rtb*hw(3,i1c,i2c,jd,1)
      b(i1c,i2c,id,3,1,2) = b(i1c,i2c,id,3,1,2) + rbt*hw(3,i1c,i2c,jd,2)
      b(i1c,i2c,id,3,2,2) = b(i1c,i2c,id,3,2,2) + rbb*hw(3,i1c,i2c,jd,2)
      b(i1c,k2c,id,1,2,1) = b(i1c,k2c,id,1,2,1) + rtt*hw(6,i1c,k2c,jd,1)
      b(i1c,k2c,id,1,3,1) = b(i1c,k2c,id,1,3,1) + rtb*hw(6,i1c,k2c,jd,1)
      b(i1c,k2c,id,1,1,2) = b(i1c,k2c,id,1,1,2) + rbt*hw(6,i1c,k2c,jd,2)
      b(i1c,k2c,id,1,2,2) = b(i1c,k2c,id,1,2,2) + rbb*hw(6,i1c,k2c,jd,2)
 
      stt = ((hw(2,j1c,k2c,jd,1)*a(i1f,k2f,id,4,2,1))
     &      + hw(7,j1c,k2c,jd,1)*a(i1f,k2f,id,5,2,1))
      stb = ((hw(2,j1c,k2c,jd,2)*a(i1f,k2f,id,4,3,1))
     &      + hw(7,j1c,k2c,jd,2)*a(i1f,k2f,id,5,3,1))
      sbt = ((hw(2,j1c,k2c,jd,1)*a(i1f,k2f,id,4,1,2))
     &      + hw(7,j1c,k2c,jd,1)*a(i1f,k2f,id,5,1,2))
      sbb = ((hw(2,j1c,k2c,jd,2)*a(i1f,k2f,id,4,2,2))
     &      + hw(7,j1c,k2c,jd,2)*a(i1f,k2f,id,5,2,2))
 
      b(i1c,i2c,id,4,2,1) = b(i1c,i2c,id,4,2,1) + stt*hw(3,i1c,i2c,jd,1)
      b(i1c,i2c,id,4,3,1) = b(i1c,i2c,id,4,3,1) + stb*hw(3,i1c,i2c,jd,1)
      b(i1c,i2c,id,4,1,2) = b(i1c,i2c,id,4,1,2) + sbt*hw(3,i1c,i2c,jd,2)
      b(i1c,i2c,id,4,2,2) = b(i1c,i2c,id,4,2,2) + sbb*hw(3,i1c,i2c,jd,2)
      b(i1c,k2c,id,5,2,1) = b(i1c,k2c,id,5,2,1) + stt*hw(6,i1c,k2c,jd,1)
      b(i1c,k2c,id,5,3,1) = b(i1c,k2c,id,5,3,1) + stb*hw(6,i1c,k2c,jd,1)
      b(i1c,k2c,id,5,1,2) = b(i1c,k2c,id,5,1,2) + sbt*hw(6,i1c,k2c,jd,2)
      b(i1c,k2c,id,5,2,2) = b(i1c,k2c,id,5,2,2) + sbb*hw(6,i1c,k2c,jd,2)
 
 50   continue
 
      end
