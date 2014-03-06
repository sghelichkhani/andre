*DK addalp
      subroutine addalp(a,alp,nt)
 
c...  This routine assembles the node-to-node operator a by summing
c...  the Laplacian and grad div operators in appropriate proportion.
 
      real a(4,7*(nt+1)**2,3,3), alp(7*(nt+1)**2,2)
 
c...  Note:  A nonzero bulk viscosity is applied here.  This has the
c...  effect of suppressing divergence in the velocity field similar
c...  to what is done when the 'penalty method' is applied.
 
      bulkvisc      = 1./3.
      divfac        = bulkvisc - 2./3.
      oneplusdivfac = 1. + divfac
 
      do j=1,3
         do ii=1,(nt+1)**2*7
            a(1,ii,j,j) = oneplusdivfac*a(1,ii,j,j) + alp(ii,1)
            a(2,ii,j,j) = oneplusdivfac*a(2,ii,j,j)
            a(3,ii,j,j) = oneplusdivfac*a(3,ii,j,j)
            a(4,ii,j,j) = oneplusdivfac*a(4,ii,j,j) + alp(ii,2)
         end do
      end do
 
      do j=1,4
         do ii=1,(nt+1)**2*7
            a12 = a(j,ii,1,2)
            a21 = a(j,ii,2,1)
            a13 = a(j,ii,1,3)
            a31 = a(j,ii,3,1)
            a23 = a(j,ii,2,3)
            a32 = a(j,ii,3,2)
            a(j,ii,1,2) = a12 + divfac*a21
            a(j,ii,2,1) = a21 + divfac*a12
            a(j,ii,1,3) = a13 + divfac*a31
            a(j,ii,3,1) = a31 + divfac*a13
            a(j,ii,2,3) = a23 + divfac*a32
            a(j,ii,3,2) = a32 + divfac*a23
         end do
      end do
 
      end
*dk areacalc
      subroutine areacalc(area,xn,nt)
 
c...  This routine computes the areas of the spherical triangles in
c...  a diamond on the unit sphere.
 
      real area(nt+1,nt+1,2), xn(nt+1,nt+1,10,3), xv(3,3)
 
      call nulvec(area, 2*(nt+1)**2)
 
      do i2=1,nt
         do i1=2,nt+1
 
            xv(1,1) = xn(i1  ,i2  ,1,1)
            xv(1,2) = xn(i1  ,i2  ,1,2)
            xv(1,3) = xn(i1  ,i2  ,1,3)
            xv(2,1) = xn(i1-1,i2+1,1,1)
            xv(2,2) = xn(i1-1,i2+1,1,2)
            xv(2,3) = xn(i1-1,i2+1,1,3)
            xv(3,1) = xn(i1-1,i2  ,1,1)
            xv(3,2) = xn(i1-1,i2  ,1,2)
            xv(3,3) = xn(i1-1,i2  ,1,3)
            t1 = xv(1,1)*xv(2,1) + xv(1,2)*xv(2,2) + xv(1,3)*xv(2,3)
            t2 = xv(2,1)*xv(3,1) + xv(2,2)*xv(3,2) + xv(2,3)*xv(3,3)
            t3 = xv(3,1)*xv(1,1) + xv(3,2)*xv(1,2) + xv(3,3)*xv(1,3)
            t1 = 0.5*acos(t1)
            t2 = 0.5*acos(t2)
            t3 = 0.5*acos(t3)
            s  = 0.5*(t1 + t2 + t3)
            a  = tan(s)*tan(s-t1)*tan(s-t2)*tan(s-t3)
            area(i1,i2,1) = 4.*atan(sqrt(a))

            xv(1,1) = xn(i1  ,i2  ,1,1)
            xv(1,2) = xn(i1  ,i2  ,1,2)
            xv(1,3) = xn(i1  ,i2  ,1,3)
            xv(2,1) = xn(i1  ,i2+1,1,1)
            xv(2,2) = xn(i1  ,i2+1,1,2)
            xv(2,3) = xn(i1  ,i2+1,1,3)
            xv(3,1) = xn(i1-1,i2+1,1,1)
            xv(3,2) = xn(i1-1,i2+1,1,2)
            xv(3,3) = xn(i1-1,i2+1,1,3)
            t1 = xv(1,1)*xv(2,1) + xv(1,2)*xv(2,2) + xv(1,3)*xv(2,3)
            t2 = xv(2,1)*xv(3,1) + xv(2,2)*xv(3,2) + xv(2,3)*xv(3,3)
            t3 = xv(3,1)*xv(1,1) + xv(3,2)*xv(1,2) + xv(3,3)*xv(1,3)
            t1 = 0.5*acos(t1)
            t2 = 0.5*acos(t2)
            t3 = 0.5*acos(t3)
            s  = 0.5*(t1 + t2 + t3)
            a  = tan(s)*tan(s-t1)*tan(s-t2)*tan(s-t3)
            area(i1,i2,2) = 4.*atan(sqrt(a))

         end do
      end do

      end
*dk feoper
      subroutine feoper(xne,xnbr)
 
c...  This routine generates the Laplacian, grad div, and gradient
c...  operators.
 
      include 'size.h'
      include 'pcom.h'
      real xne(-1:nt+1,0:nt+2,nd,3), xnbr(0:nt,nt+1,6,2)
      real xnv(3,3), t(3), xns(3,2), fn(3,3), adl(3,3), ur(3,3)
      real xg(-1:nt+1,0:nt+2,3,3)
      common /oper/ ra(3,8,nr+1), alp(0:6,0:nt,nt+1,2),
     &              atn(4,0:6,0:nt,nt+1,3,3)
      common /grad/ rg(3,2,nr+1), grd(0:6,0:nt,nt+1,3,2)
      common /ofst/ j1n(0:6), j2n(0:6), md(0:6)
 
      third  = 1./3.
      sixth  = 1./6.
      twlfth = 1./12.
 
      call localcoord(xnbr,xg,xne)
 
      call nulvec(alp, (nt+1)**2*14)
      call nulvec(atn, (nt+1)**2*252)
      call nulvec(grd, (nt+1)**2*42)
 
      do i2=1,nt+1
         do i1=0,nt
 
            do m=1,6
 
               n = mod(m, 6) + 1
 
               if(i1+j1n(m).ge.0    .and. i1+j1n(n).ge.0    .and.
     &            i1+j1n(m).le.nt   .and. i1+j1n(n).le.nt   .and.
     &            i2+j2n(m).ge.1    .and. i2+j2n(n).ge.1    .and.
     &            i2+j2n(m).le.nt+1 .and. i2+j2n(n).le.nt+1) then
 
c...              Load the Cartesian coordinates of the triangle
c...              vertices.
 
                  do j=1,3
                     xnv(1,j) = xne(i1       ,i2       ,1,j)
                     xnv(2,j) = xne(i1+j1n(m),i2+j2n(m),1,j)
                     xnv(3,j) = xne(i1+j1n(n),i2+j2n(n),1,j)
                  end do
 
c...              Compute the triangle area on the unit sphere.
 
                  s = 0.
 
                  do n1=1,3
                     n2 = mod(n1,3) + 1
                     t(n1) = acos(xnv(n1,1)*xnv(n2,1)
     &                          + xnv(n1,2)*xnv(n2,2)
     &                          + xnv(n1,3)*xnv(n2,3))
                     s  = s + 0.5*t(n1)
                  end do
 
                  w = tan(0.5*s)*tan(0.5*(s - t(1)))
     &                          *tan(0.5*(s - t(2)))
     &                          *tan(0.5*(s - t(3)))
                  area  = 4.*atan(sqrt(w))
                  rarea = 1./area
                  ar6   = area*sixth
                  ar12  = area*twlfth
 
c...              Load local spherical coordinates of triangle
c...              vertices.
 
                  do j=1,2
                     xns(1,j) = 0.
                     xns(2,j) = xnbr(i1,i2,m,j)
                     xns(3,j) = xnbr(i1,i2,n,j)
                  end do
 
c...              Compute triangle face normals.
 
                  do n1=1,3
                     n2 = mod(n1,3) + 1
                     n3 = mod(n2,3) + 1
                     f1 = xns(n2,2)   - xns(n3,2)
                     f2 = xns(n3,1)*cos(xns(n3,2))
     &                  - xns(n2,1)*cos(xns(n2,2))
                     fn(n1,1) = f1*xg(i1,i2,1,2)
     &                        + f2*xg(i1,i2,1,3)
                     fn(n1,2) = f1*xg(i1,i2,2,2)
     &                        + f2*xg(i1,i2,2,3)
                     fn(n1,3) = f1*xg(i1,i2,3,2)
     &                        + f2*xg(i1,i2,3,3)
                  end do
 
c...              Compute triangle normals.
 
                  ur1 = xnv(1,1) + xnv(2,1) + xnv(3,1)
                  ur2 = xnv(1,2) + xnv(2,2) + xnv(3,2)
                  ur3 = xnv(1,3) + xnv(2,3) + xnv(3,3)
                  urn = 1./sqrt(ur1**2 + ur2**2 + ur3**2)
                  ur1 = urn*ur1
                  ur2 = urn*ur2
                  ur3 = urn*ur3
 
                  do n1=1,3
                     ur(n1,1) = ur1
                     ur(n1,2) = ur2
                     ur(n1,3) = ur3
                  end do
 
c...              Construct the triangle gradient operator.
 
                  do n1=1,3
                     n2 = mod(n1,3) + 1
                     n3 = mod(n2,3) + 1
                     adl(n1,1) = (2.*fn(n1,1) - fn(n2,1)
     &                             - fn(n3,1))*sixth
                     adl(n1,2) = (2.*fn(n1,2) - fn(n2,2)
     &                             - fn(n3,2))*sixth
                     adl(n1,3) = (2.*fn(n1,3) - fn(n2,3)
     &                             - fn(n3,3))*sixth
                  end do
 
c...              Construct the nodal Laplacian and mass matrix
c...              operators.
 
                  dot = adl(1,1)*adl(1,1) + adl(1,2)*adl(1,2)
     &                + adl(1,3)*adl(1,3)
                  alp(0,i1,i2,1) = alp(0,i1,i2,1) + dot*rarea
                  alp(0,i1,i2,2) = alp(0,i1,i2,2) + ar6
 
                  dot = adl(1,1)*adl(2,1) + adl(1,2)*adl(2,2)
     &                + adl(1,3)*adl(2,3)
                  alp(m,i1,i2,1) = alp(m,i1,i2,1) + dot*rarea
                  alp(m,i1,i2,2) = alp(m,i1,i2,2) + ar12
 
                  dot = adl(1,1)*adl(3,1) + adl(1,2)*adl(3,2)
     &                + adl(1,3)*adl(3,3)
                  alp(n,i1,i2,1) = alp(n,i1,i2,1) + dot*rarea
                  alp(n,i1,i2,2) = alp(n,i1,i2,2) + ar12
 
c...              Construct the nodal gradient and grad div
c...              operators.
 
                  do k1=1,3
 
                     grd(0,i1,i2,k1,1) = grd(0,i1,i2,k1,1)
     &                                 + third*adl(1,k1)
                     grd(m,i1,i2,k1,1) = grd(m,i1,i2,k1,1)
     &                                 + third*adl(2,k1)
                     grd(n,i1,i2,k1,1) = grd(n,i1,i2,k1,1)
     &                                 + third*adl(3,k1)
 
                     grd(0,i1,i2,k1,2) = grd(0,i1,i2,k1,2)
     &                                 +   ar6* ur(1,k1)
                     grd(m,i1,i2,k1,2) = grd(m,i1,i2,k1,2)
     &                                 +  ar12* ur(1,k1)
                     grd(n,i1,i2,k1,2) = grd(n,i1,i2,k1,2)
     &                                 +  ar12* ur(1,k1)
 
                     do k2=1,3
 
                        atn(1,0,i1,i2,k2,k1) = atn(1,0,i1,i2,k2,k1)
     &                                  + rarea*adl(1,k2)*adl(1,k1)
                        atn(1,m,i1,i2,k2,k1) = atn(1,m,i1,i2,k2,k1)
     &                                  + rarea*adl(1,k2)*adl(2,k1)
                        atn(1,n,i1,i2,k2,k1) = atn(1,n,i1,i2,k2,k1)
     &                                  + rarea*adl(1,k2)*adl(3,k1)
 
                        atn(2,0,i1,i2,k2,k1) = atn(2,0,i1,i2,k2,k1)
     &                                  + third*adl(1,k2)* ur(1,k1)
                        atn(2,m,i1,i2,k2,k1) = atn(2,m,i1,i2,k2,k1)
     &                                  + third*adl(1,k2)* ur(2,k1)
                        atn(2,n,i1,i2,k2,k1) = atn(2,n,i1,i2,k2,k1)
     &                                  + third*adl(1,k2)* ur(3,k1)
 
                        atn(3,0,i1,i2,k2,k1) = atn(3,0,i1,i2,k2,k1)
     &                                  + third* ur(1,k2)*adl(1,k1)
                        atn(3,m,i1,i2,k2,k1) = atn(3,m,i1,i2,k2,k1)
     &                                  + third* ur(1,k2)*adl(2,k1)
                        atn(3,n,i1,i2,k2,k1) = atn(3,n,i1,i2,k2,k1)
     &                                  + third* ur(1,k2)*adl(3,k1)
 
                        atn(4,0,i1,i2,k2,k1) = atn(4,0,i1,i2,k2,k1)
     &                                  + ar6  * ur(1,k2)* ur(1,k1)
                        atn(4,m,i1,i2,k2,k1) = atn(4,m,i1,i2,k2,k1)
     &                                  + ar12 * ur(1,k2)* ur(2,k1)
                        atn(4,n,i1,i2,k2,k1) = atn(4,n,i1,i2,k2,k1)
     &                                  + ar12 * ur(1,k2)* ur(3,k1)
 
                     end do
                  end do
 
               end if
 
            end do
         end do
      end do
 
      call addalp(atn,alp,nt)
 
      end
*dk grdgen
      subroutine grdgen(xn,nt)
 
c...  This routine generates the nodal coordinates xn for an
c...  icosahedral grid on the unit sphere.  The grid resolution
c...  corresponds to a subdivision of the edges of the original
c...  icosahedral triangles into nt equal parts.
 
      real xn(nt+1,nt+1,10,3)
 
      fifthpi = 0.4*asin(1.)
      w       = 2.0*acos(1./(2.*sin(fifthpi)))
      cosw    = cos(w)
      sinw    = sin(w)
      lvt     = 1.45*log(real(nt))
      nn      = (nt+1)**2*10
 
      do id=1,10
 
         sgn = 1.
         if(id .ge. 6) sgn = -1.
         phi = (2*mod(id, 5) - 3 + (id - 1)/5)*fifthpi
 
         xn(   1,   1,id,1) =  0.
         xn(   1,   1,id,2) =  0.
         xn(   1,   1,id,3) =  sgn
         xn(nt+1,   1,id,1) =  sinw*cos(phi)
         xn(nt+1,   1,id,2) =  sinw*sin(phi)
         xn(nt+1,   1,id,3) =  cosw*sgn
         xn(   1,nt+1,id,1) =  sinw*cos(phi + fifthpi + fifthpi)
         xn(   1,nt+1,id,2) =  sinw*sin(phi + fifthpi + fifthpi)
         xn(   1,nt+1,id,3) =  cosw*sgn
         xn(nt+1,nt+1,id,1) =  sinw*cos(phi + fifthpi)
         xn(nt+1,nt+1,id,2) =  sinw*sin(phi + fifthpi)
         xn(nt+1,nt+1,id,3) = -cosw*sgn
 
         do k=0,lvt-1
 
            m  = 2**k
            l  = nt/m
            l2 = l/2
 
c...        rows of diamond--
 
            do j1=1,m+1
               do j2=1,m
                     i1 = (j1-1)*l + 1
                     i2 = (j2-1)*l + l2 + 1
                     call midpt(xn(i1,i2,id,1),xn(i1,i2-l2,id,1),
     &                          xn(i1,i2+l2,id,1),nn)
               end do
            end do
 
c...        columns of diamond--
 
            do j1=1,m+1
               do j2=1,m
                     i1 = (j2-1)*l + l2 + 1
                     i2 = (j1-1)*l + 1
                     call midpt(xn(i1,i2,id,1),xn(i1-l2,i2,id,1),
     &                          xn(i1+l2,i2,id,1),nn)
               end do
            end do
 
c...        diagonals of diamond--
 
            do j1=1,m
               do j2=1,m
                     i1 = (j1-1)*l + l2 + 1
                     i2 = (j2-1)*l + l2 + 1
                     call midpt(xn(i1,i2,id,1),xn(i1-l2,i2+l2,id,1),
     &                          xn(i1+l2,i2-l2,id,1),nn)
               end do
            end do
 
         end do
 
      end do
 
      end
*dk gridinit
      subroutine gridinit
 
c...  This routine initializes all arrays in common blocks /radl/,
c...  /mesh/, /ndar/, /volm/, and /xnex/.
 
      include 'size.h'
      include 'pcom.h'
      parameter (nxm=4000+(nt+1)**2*41)
      parameter (nopr=(nt/2+1)**2*ndo*189*(nr/2+1)*7/5+8000)
      common /fopr/ a(nopr+nv*ndo/nd*9+nv*ndo/nd*18*5/4),
     &              mopr(0:10), mb(0:10)
      common /grid/ mxm(0:10), xm(nxm)
      common /mesh/ xn((nt+1)**2,nd,3)
      common /ndar/  arn((nt+1)**2),  arne((nt+1)**2),
     &              rarn((nt+1)**2), rarne((nt+1)**2)
      common /radl/ rshl(nr+1), ird
      common /volm/ vol((nt+1)**2*(nr+1),2)
      common /xnex/ xne((nt+3)**2*nd*3)
 
      lvf = 1.45*log(real(mt))
      lvg = 1.45*log(real(mt/nt))
 
c...  Generate the nodal coordinate array xm for all grid levels.
 
      do lv=0,lvf
 
         kt = 2**lv
 
         call grdgen(a,kt)
 
         if(nproc.eq.1 .or. lv.eq.0) then
            call scopy((kt+1)**2*30, a, 1, xm(mxm(lv)), 1)
         else
            jt     = 2**max(0, lv - lvg)
            lvproc = min(lvg, lv)
            mproc  = 2**(2*lvproc)
            call subarray(a, xm(mxm(lv)), 0, 10, nd, kt, jt, 3)
         endif
 
      end do
 
      lvproc = lvg
      mproc  = 2**(2*lvproc)
 
c...  Generate nodal coordinate array xn for finest level. (This
c...  duplicates part of array xm, but is done for convenience
c...  because xn does not require use of a pointer.)
 
      mm = (mt + 1)**2
      m0 = mm*30 + 1
 
      call grdgen(a,mt)
 
      if(nproc .eq. 1) then
         call scopy((mt+1)**2*30, a, 1, xn, 1)
      else
         call subarray(a,xn,0,10,nd,mt,nt,3)
      endif
 
c...  Generate the nodal area arrays.
 
      call ndarea(a(m0),a(m0+mm),a(m0+2*mm),a(m0+3*mm),a(m0+4*mm),a,mt)
 
      if(nproc .eq. 1) then
         call scopy((mt+1)**2, a(m0),      1, arn,   1)
         call scopy((mt+1)**2, a(m0+mm),   1, arne,  1)
         call scopy((mt+1)**2, a(m0+2*mm), 1, rarn,  1)
         call scopy((mt+1)**2, a(m0+3*mm), 1, rarne, 1)
      else
         call subarray(a(m0),     arn,  1,1,1,mt,nt,1)
         call subarray(a(m0+mm),  arne, 0,1,1,mt,nt,1)
         call subarray(a(m0+2*mm),rarn, 1,1,1,mt,nt,1)
         call subarray(a(m0+3*mm),rarne,0,1,1,mt,nt,1)
      endif
 
c...  Generate array xne that contains an extra ring of ghost nodes
c...  surrounding the subdomain.
 
      call xnextend(a(m0),a)
 
      call xnesplit(xne,a(m0))
 
c...  Generate the tangential finite element operator arrays needed
c...  to compute Laplacian, grad div, gradient, and mass matrix.
 
      if(mynum .ge. mproc) then
 
c...     The following forces the operators for southern hemisphere
c...     diamonds to be built using diamond 1 coordinates so that
c...     rotations to diamond 1 work correctly.
 
         mproc = nproc
         call xnesplit(a,a(m0))
         call feoper(a,a(m0))
         mproc = 2**(2*lvproc)
 
      else
 
         call feoper(xne,a)
 
      endif
 
c...  Generate the radial grid array rshl and the radial operators.
 
      call rshlgen
 
      call rgen
 
c...  Generate nodal volumes and reciprocal volumes.
 
      call volume(vol(1,1),0,0)
      call volume(vol(1,2),1,1)
 
      end

*dk localcoord
      subroutine localcoord(xnbr,xg,xne)
 
c...  This routine calculates the unit vectors that define the
c...  local spherical coordinate system at each node.  It also
c...  computes the coordinates of the neighboring nodes.
 
      include 'size.h'
      include 'pcom.h'
      real  xg(-1:nt+1,0:nt+2,3,3), xnbr(0:nt,nt+1,6,2)
      real xne(-1:nt+1,0:nt+2,nd,3)
      common /ofst/ j1n(0:6), j2n(0:6), md(0:6)
 
      chi   = 0.8*asin(1.)
      iproc = mod(mynum, 2**(2*lvproc))
 
c...  Obtain unit longitude and unit latitude vectors for extended
c...  array xg.
 
      do i2=0,nt+2
         do i1=-1,nt+1
            t1 = -xne(i1,i2,1,2)
            t2 =  xne(i1,i2,1,1)
            tn =  1./sqrt((1.e-100 + t1**2) + t2**2)
            t1 =  tn*t1
            t2 =  tn*t2
            xg(i1,i2,1,2)  =  t1
            xg(i1,i2,2,2)  =  t2
            xg(i1,i2,3,2)  =  0.
            xg(i1,i2,1,3)  = -xne(i1,i2,1,3)*t2
            xg(i1,i2,2,3)  =  xne(i1,i2,1,3)*t1
            xg(i1,i2,3,3)  =  xne(i1,i2,1,1)*t2 - xne(i1,i2,1,2)*t1
            xg(i1,i2,1,1)  =  xne(i1,i2,1,1)
            xg(i1,i2,2,1)  =  xne(i1,i2,1,2)
            xg(i1,i2,3,1)  =  xne(i1,i2,1,3)
         end do
      end do
 
      if(iproc .eq. 0) then
         xg( 0, 1,1,2) = -sin(-0.5*chi)
         xg( 0, 1,2,2) =  cos(-0.5*chi)
         xg( 0, 1,1,3) = -cos(-0.5*chi)
         xg( 0, 1,2,3) = -sin(-0.5*chi)
      endif
 
c...  Compute the local coordinates of the neighboring nodes in the
c...  coordinate system of node (i1,i2).
 
      do  m=1,6
         do i2=1,nt+1
            j2  = i2 + j2n(m)
 
            do i1=0,nt
               j1  = i1 + j1n(m)
               d1  = xg(i1,i2,1,1)*xg(j1,j2,1,1)
     &             + xg(i1,i2,2,1)*xg(j1,j2,2,1)
     &             + xg(i1,i2,3,1)*xg(j1,j2,3,1)
               d2  = xg(i1,i2,1,2)*xg(j1,j2,1,1)
     &             + xg(i1,i2,2,2)*xg(j1,j2,2,1)
     &             + xg(i1,i2,3,2)*xg(j1,j2,3,1)
               d3  = xg(i1,i2,1,3)*xg(j1,j2,1,1)
     &             + xg(i1,i2,2,3)*xg(j1,j2,2,1)
     &             + xg(i1,i2,3,3)*xg(j1,j2,3,1)
               d3  = sign(min(1., abs(d3)), d3)
               tht = asin(d3)
               phi = atan2(d2, d1 + 1.e-100)
               if(iproc.eq.0 .and. i1.eq.0 .and. i2.eq.1) phi = asin(d2)
               xnbr(i1,i2,m,1) = phi
               xnbr(i1,i2,m,2) = tht
            end do
 
         end do
      end do
 
      end
*dk opinit
      subroutine opinit
 
c...  This routine initializes arrays in common blocks /grid/, /fopr/
c...  /wght/, and /rttn/ and calls routine gridinit which initializes
c...  grid and operator arrays.
 
      include 'size.h'
      include 'pcom.h'
      parameter (nxm=4000+(nt+1)**2*41)
      parameter (nopr=(nt/2+1)**2*ndo*189*(nr/2+1)*7/5+8000)
 
      common /grid/ mxm(0:10), xm(nxm)
      common /fopr/ a(nopr+nv*ndo/nd*9+nv*ndo/nd*18*5/4),
     &              mopr(0:10), mb(0:10)
      common /wght/ mrw(0:10), rw(nv*ndo/nd*7/5),
     &              mhw(0:10), hw(nv*ndo/nd*7/5)
      common /rttn/ phi(10), trtn(3,3,10)
      common /clck/ itmng, sec(50)
      call mytime(tin)
 
      lvf = 1.45*log(real(mt))
      lvr = 1.45*log(real(nr))
      lvg = 1.45*log(real(mt/nt))
 
      mxm(0)  = 1
      mb(0)   = 1
      mb(1)   = 1
      mopr(0) = 1
      mrw(0)  = 1
      mrw(1)  = 11
      mhw(0)  = 13
 
      do lv=0,lvf
         kr = 2**max(0, lv - lvf + lvr)
         kt = 2**max(0, lv - lvg)
         ic = 189*ndo
         if(lv .eq. lvf) ic = 0
         if(lv .eq.   0) ic = 1890
         mopr(lv+1) =  mopr(lv) + (kt+1)**2*ic*(kr+1)
         mhw(lv+1)  =  mhw(lv)  + (kt+1)**2*ndo*6*(kr+1)
         if(lv .ge. 1) then
            mxm(lv+1) = mxm(lv) + (kt+1)**2*nd*3
            mb(lv+1)  = mb(lv)  + (kt+1)**2*ndo*18*(kr+1)
            mrw(lv+1) = mrw(lv) + (2*kt+1)**2*ndo*2*(kr+1)
         else
            mxm(lv+1) = mxm(lv) + (kt+1)**2*10*3
         end if
      end do
      mrw(lvf+1) = mrw(lvf)
      mhw(lvf+1) = mhw(lvf)
 
      if(mynum .eq. 0) then
         nb   = nv*ndo/nd*18*5/4
         nhw  = nv*ndo/nd*7/5
         nrw  = nv*ndo/nd*7/5
!         write(6,10) nopr, mopr(lvf+1), nb,   mb(lvf+1),
!     &               nhw,   mhw(lvf+1), nrw, mrw(lvf+1)
 10      format(/'          nopr, mopr(lvf+1) = ', 2i10/
     &           '          nb,     mb(lvf+1) = ', 2i10/
     &           '          nhw,   mhw(lvf+1) = ', 2i10/
     &           '          nrw,   mrw(lvf+1) = ', 2i10/)
         if(nopr.lt.mopr(lvf+1) .or. nhw.lt.mhw(lvf+1)
     &     .or. nb.lt.mb(lvf+1) .or. nrw.lt.mrw(lvf+1)) then
!            write(6,20)
 20         format(/' STOP--not enough memory allocated for operators')
            stop
         end if
      end if
 
c     Generate rotation matrix.
 
      call nulvec(trtn,90)
 
      do id=1,10
 
         trtn(1,1,id) =  cos(phi(id))
         trtn(1,2,id) = -sin(phi(id))
         trtn(3,3,id) =  1.0
 
         if(id .ge. 6) then
            trtn(1,1,id) = -trtn(1,1,id)
            trtn(1,2,id) = -trtn(1,2,id)
            trtn(3,3,id) = -1.0
         endif
 
         trtn(2,2,id) =  trtn(1,1,id)
         trtn(2,1,id) = -trtn(1,2,id)
 
      end do
 
c...  Compute all other grid arrays.
 
      call gridinit
 
      call mytime(tout)
      sec(1) = sec(1) + tout - tin
      end
*dk oprot
      subroutine oprot(a,t,nskp,np)
 
c     This routine performs on the operator a the rotation defined
c     by the matrix t.
 
      real a(nskp,np,3,3), t(3,3), c(3,3)
 
      do ii=1,np
 
         do i=1,3
            do j=1,3
               c(i,j) = 0.
               do k=1,3
                  c(i,j) = c(i,j) + a(1,ii,i,k)*t(j,k)
               end do
            end do
         end do
 
         do i=1,3
            do j=1,3
               a(1,ii,i,j) = 0.
               do k=1,3
                  a(1,ii,i,j) = a(1,ii,i,j) + t(i,k)*c(k,j)
               end do
            end do
         end do
 
      end do
 
      end
*dk ndarea
      subroutine ndarea(arn,arne,rarn,rarne,area,xn,nt)
 
c...  This routine computes the areas, arn, associated with the nodes
c...  on the unit sphere as well as the reciprocal areas, rarn.  These
c...  arrays have zero values along the upper right and lower right
c...  diamond edges.  Arrays arne and rarne are identical except they
c...  contain the actual areas and reciprocal areas, respectively, in
c...  these edge locations.
 
      real  arn(nt+1,nt+1),  rarn(nt+1,nt+1), area(nt+1,nt+1,2)
      real arne(nt+1,nt+1), rarne(nt+1,nt+1), xn(*)
 
      call areacalc(area,xn,nt)
      call nulvec(arn,  (nt+1)**2)
      call nulvec(rarn, (nt+1)**2)
 
c...  Treat interior nodes.
 
      do i2=2,nt
         do i1=2,nt
            arn(  i1,i2) = (area(i1  ,i2  ,1) + area(i1  ,i2  ,2)
     &                   +  area(i1+1,i2  ,1) + area(i1  ,i2-1,2)
     &                   +  area(i1+1,i2-1,1) + area(i1+1,i2-1,2))/3.
            rarn( i1,i2) =  1./arn(i1,i2)
            arne( i1,i2) =     arn(i1,i2)
            rarne(i1,i2) =    rarn(i1,i2)
         end do
      end do
 
c...  Treat edge nodes.
 
      do i=2,nt
         arn(     i,1) = (area(   i,1  ,1) + area(i   ,1,2)
     &                 +  area( i+1,1  ,1))/1.5
         arn(  nt+1,i) = (area(nt+1,i  ,1) + area(nt+1,i,2)
     &                 +  area(nt+1,i-1,2))/1.5
         rarn(    i,1) =  1./arn(   i,1)
         rarn( nt+1,i) =  1./arn(nt+1,i)
         arne(    i,1) =     arn(   i,1)
         arne( nt+1,i) =     arn(nt+1,i)
         rarne(   i,1) =    rarn(   i,1)
         rarne(nt+1,i) =    rarn(nt+1,i)
         arne(    1,i) =     arn(   i,1)
         arne( i,nt+1) =     arn(nt+1,i)
         rarne(   1,i) =    rarn(   i,1)
         rarne(i,nt+1) =    rarn(nt+1,i)
      end do
 
c...  Treat pentagonal nodes.
 
      arn(     1,   1) = 5.*area(2,1,1)/3.
      arn(  nt+1,   1) =    arn(1,1)
      rarn(    1,   1) = 1./arn(1,1)
      rarn( nt+1,   1) =   rarn(1,1)
      arne(    1,   1) =    arn(1,1)
      arne( nt+1,   1) =    arn(1,1)
      arne(    1,nt+1) =    arn(1,1)
      arne( nt+1,nt+1) =    arn(1,1)
      rarne(   1,   1) =   rarn(1,1)
      rarne(nt+1,   1) =   rarn(1,1)
      rarne(   1,nt+1) =   rarn(1,1)
      rarne(nt+1,nt+1) =   rarn(1,1)
 
	end subroutine


*dk pseugen
	subroutine pseugen(nj,ibc)
!	This routine generates a Moore-Penrose pseudo-inverse
!	of the stiffness operator at the lowest, or k=0, level.
 
	include 'size.h'
	include 'pcom.h'

	common /pseu/ v(2,2,10,3,2), u(2,2,10,3,2), w(10)
	common /opwt/ wd, wm
	common /pseu/ ainvd1(576), ainvm1(576), ainvd3(5184), ainvm3(5184)
	common /radl/ rshl(nr+1), ird

	integer map(0:nt,nt+1,nd), ibc
	real aa(5184), bb(5184), det(2)
	real urot(3,pl_size)

	map=0
	urot=0.0

      nn   = 24*nj
      epsq = 1.e03*wd*rshl(1)**2
 
      do id=0,11
 
         if(id.eq.0.or.id.eq.11) then
            ll = iabs(id-1)
            i1 = 1
         else
            ll = id
            i1 = 2
         endif
 
         do ir=1,2
            do i=1,nj
               call vdelgen(u,i1,1,ll,i,ir,nj)
               call nulvec(v,80*nj)
				! map and urot don't matter because kt=1 and so, the if
				! condition in axu3s is never satisfied
               call axu3s(v,u,nj,0,10,1,1,map,urot,ibc)
               icol = i + (ir-1)*nj + id*2*nj
               call vnn(aa(1+nn*(icol-1)),v,nj)
            end do
         end do
 
      end do
 
      call matprd3(bb,aa,aa,nn,nn)
 
      bmax = 1.
      do j=1,nn
         jj = j + nn*(j-1)
         bmax = max(bmax, bb(jj))
      end do
 
      do j=1,nn
         jj = j + nn*(j-1)
         bb(jj) = bb(jj) + 1.e-06*bmax
c        bb(jj) = bb(jj)*(1. + 1.e-06)
c        if(bb(jj) .eq. 0.) bb(jj) = 1.e-06*bmax
      end do
 
      call spoco(bb,nn,nn,rcond,u,info)
 
      if(mynum.eq.0 .and. (rcond.lt.1.e-10.or.info.ne.0))
     &   write(6,10) info,rcond
 10   format(/" in pseugen   info =",i4,",   rcond =",1pe9.2)
 
      call spodi(bb,nn,nn,det,1)
 
      call ltrfill(bb,nn)
 
      if(nj.eq.1) then
         if(wd.ne.0.) then
            call matprd3(ainvd1,bb,aa,nn,nn)
         else
            call matprd3(ainvm1,bb,aa,nn,nn)
         endif
      else
         if(wd.ne.0.) then
            call matprd3(ainvd3,bb,aa,nn,nn)
         else
            call matprd3(ainvm3,bb,aa,nn,nn)
         endif
      endif
 
      end
*dk rflayr
      subroutine rflayr(rf,rt,rb,nr)
 
c...  This routine generates the radial factors needed to assemble the
c...  3-D node-to-node operators for the shell between radii rb and rt.
 
      real rf(3,9)
 
      d2 = 1./(rt - rb)**2
      r1 = d2*(rt    - rb)
      r2 = d2*(rt**2 - rb**2)/2.
      r3 = d2*(rt**3 - rb**3)/3.
      r4 = d2*(rt**4 - rb**4)/4.
      r5 = d2*(rt**5 - rb**5)/5.
 
      call nulvec(rf,24)
 
      rf(1,1) = -r3 + (rb+rb)*r2 - rb*rb*r1
      rf(2,1) =  r3 - (rt+rb)*r2 + rt*rb*r1
      rf(3,1) = -r3 + (rt+rt)*r2 - rt*rt*r1
      rf(1,2) = -r3 + rb*r2
      rf(2,2) =  r3 - rb*r2
      rf(3,2) = -r3 + rt*r2
      rf(1,3) = -r3 + rb*r2
      rf(2,3) =  r3 - rt*r2
      rf(3,3) = -r3 + rt*r2
      rf(1,4) = -r3
      rf(2,4) =  r3
      rf(3,4) = -r3
      rf(1,5) = -r3 + (rb+rb)*r2 - rb*rb*r1
      rf(2,5) =  r3 - (rt+rb)*r2 + rt*rb*r1
      rf(3,5) = -r3 + (rt+rt)*r2 - rt*rt*r1
      rf(1,6) = -r3
      rf(2,6) =  r3
      rf(3,6) = -r3
      rf(1,7) =  r5 - (rb+rb)*r4 + rb*rb*r3
      rf(2,7) = -r5 + (rt+rb)*r4 - rt*rb*r3
      rf(3,7) =  r5 - (rt+rt)*r4 + rt*rt*r3
      rf(1,8) =  r4 - (rb+rb)*r3 + rb*rb*r2
      rf(2,8) = -r4 + (rt+rb)*r3 - rt*rb*r2
      rf(3,8) =  r4 - (rt+rt)*r3 + rt*rt*r2
      rf(1,9) =  r4 - rb*r3
      rf(2,9) = -rf(1,9)
      rf(3,9) =  r4 - rt*r3
 
      end
*dk rgen
      subroutine rgen
 
c...  This routine generates the radial factors needed
c...  to assemble the finite-element operators.
 
      include 'size.h'
      common /radl/ rshl(nr+1), ird
      common /oper/ ra(3,8,nr+1), alp(7,(nt+1)**2,2),
     &              atn(4,7,(nt+1)**2,9)
      common /grad/ rg(3,2,nr+1), grd(7,(nt+1)**2,3,2)
      real rfu(3,9), rfl(3,9)
 
      call nulvec(rfu, 27)
 
      do ir=1,nr+1
 
         if(ir .le. nr) then
            rtop = rshl(ir)
            rbot = rshl(ir+1)
            call rflayr(rfl,rtop,rbot,nr)
         endif
 
         qu = 1.0
         ql = 1.0
         if(ir .eq. 1)    qu = 0.0
         if(ir .eq. nr+1) ql = 0.0
 
         do m=1,7
            ra(1,m,ir) = qu*rfu(2,m)
            ra(2,m,ir) = qu*rfu(3,m) + ql*rfl(1,m)
            ra(3,m,ir) = ql*rfl(2,m)
            if((ir.eq.1 .or. ir.eq.nr+1)
     &         .and. (m.ge.2 .and. m.le.3)) ra(2,m,ir) = 0.
         end do
 
         rg(1,1,ir) = qu*rfu(2,8)
         rg(2,1,ir) = qu*rfu(3,8) + ql*rfl(1,8)
         rg(3,1,ir) = ql*rfl(2,8)
         rg(1,2,ir) = qu*rfu(2,9)
         rg(2,2,ir) = qu*rfu(3,9) + ql*rfl(1,9)
         rg(3,2,ir) = ql*rfl(2,9)
 
         do j=1,3
            do i=1,9
            rfu(j,i) =  rfl(j,i)
            end do
         end do
 
         rfu(2,2) = -rfu(3,2)
         rfu(2,3) = -rfu(1,3)
         rfu(2,9) = -rfu(3,9)
 
      end do
 
      end
*dk rshlgen
      subroutine rshlgen
 
c...  This routine generates the nodal radial locations for the grid
c...  type specified by parameter ird for all grids with N .le. nr.
 
      include 'size.h'
      common /radl/ rshl(nr+1), ird
      common /dpth/ dep(128,4)
 
      pi   = 3.141592653589793
      lvr  = 1.45*log(real(nr))
      rmax = rshl(1)
      rmin = rshl(nr+1)
 
      do ir=2,nr
 
         fr = real(nr+1-ir)/real(nr)
 
         if(ird .eq. 1) then
            rshl(ir) = rmin + (rmax - rmin)*fr
         elseif(ird .le. 4) then
            aa    = 0.3 + 0.20*(ird - 2)
            rshl(ir) = (1. - aa)*(rmin + (rmax - rmin)*fr)
     &            + aa*(rmax + rmin - (rmax - rmin)*cos(pi*fr))*0.5
         else
         	! A.H.: common block 'dep' not defined for resolution finer than nr=128
         	! so for nr=256 (mt=512) and higher, ird=5 is equal to ird=1 (equal spacing)
            if(nr>128) then
            	rshl(ir)=rmin + (rmax - rmin)*fr
            else
            	rshl(ir) = rshl(ir-1) - dep(ir-1,lvr-3)*1.e3
            endif
         endif
 
      end do
 
      end
*dk tricntr
      subroutine tricntr(xc,xne,id)
 
c...  This routine computes the coordinates of the centers of
c...  the triangles of diamond id from vertex coordinates xne.
 
      include 'size.h'
      real xne(-1:nt+1,0:nt+2,nd,3), xc(nt+1,0:nt,2,3)
 
      do i1=1,nt+1
         do i2=0,nt
 
            x1 = xne(i1  ,i2,id,1) + xne(i1-1,i2+1,id,1)
     &         + xne(i1-1,i2,id,1)
            x2 = xne(i1  ,i2,id,2) + xne(i1-1,i2+1,id,2)
     &         + xne(i1-1,i2,id,2)
            x3 = xne(i1  ,i2,id,3) + xne(i1-1,i2+1,id,3)
     &         + xne(i1-1,i2,id,3)
            xx = 1./sqrt((x1**2 + x2**2) + x3**2)
            xc(i1,i2,1,1) = xx*x1
            xc(i1,i2,1,2) = xx*x2
            xc(i1,i2,1,3) = xx*x3
 
            x1 = xne(i1  ,i2,id,1) + xne(i1-1,i2+1,id,1)
     &         + xne(i1,i2+1,id,1)
            x2 = xne(i1  ,i2,id,2) + xne(i1-1,i2+1,id,2)
     &         + xne(i1,i2+1,id,2)
            x3 = xne(i1  ,i2,id,3) + xne(i1-1,i2+1,id,3)
     &         + xne(i1,i2+1,id,3)
            xx = 1./sqrt((x1**2 + x2**2) + x3**2)
            xc(i1,i2,2,1) = xx*x1
            xc(i1,i2,2,2) = xx*x2
            xc(i1,i2,2,3) = xx*x3
 
         end do
      end do
 
      end
*dk vdelgen
      subroutine vdelgen(v,i1,i2,id,j,ir,nj)
 
c...  This routine generates a field v with a (negative)
c...  delta function for component j at node (i1,i2,id,ir).
c...  This routine is required by routine pseugen which
c...  creates the pseudoinverse for the k=0 level of the multigrid
c...  solver.
 
      real v(2,2,10,nj,2)
 
      call nulvec(v,80*nj)
 
      v(i1,i2,id,j,ir) = -1.0
 
      call edgadd0(v,nj,1,1,1)
 
      end
*dk vnn
      subroutine vnn(vv,v,nj)
 
c...  This routine loads the array vv with the field values from
c...  array v for the 24 nodes at the lowest grid level such that
c...  all redundant nodes are removed.
 
      real vv(24*nj), v(2,2,10,nj,2)
 
      i = 1
 
      do ir=1,2
         do j=1,nj
            vv(i) = v(1,1,1,j,ir)
            i = i + 1
         end do
      end do
 
      do id=1,10
         do ir=1,2
            do j=1,nj
               vv(i) = v(2,1,id,j,ir)
               i = i + 1
            end do
         end do
      end do
 
      do ir=1,2
         do j=1,nj
            vv(i) = v(1,1,10,j,ir)
            i = i + 1
         end do
      end do
 
      end
*dk xnextend
      subroutine xnextend(xne,xn)
 
c...  This routine loads the extended nodal coordinate array xne
c...  which includes an extra row or column of nodes beyond each
c...  diamond boundary.
 
      include 'size.h'
      real xne(-1:mt+1,0:mt+2,10,3), xn(0:mt,mt+1,10,3)
 
      do id=1,10
 
         idlr = id + 5
         idur = mod(id  ,5) + 1
         idul = mod(id+3,5) + 1
         idll = idul + 5
         id3r = mod(id+2,5) + 1
 
         if(id .ge. 6) then
            idlr = idur
            idur = idur + 5
            idul = idul + 5
            idll = id   - 5
            id3r = id3r + 5
         endif
 
         do j=1,3
 
            do i2=1,mt+1
               do i1=0,mt
                  xne(i1,i2,id,j) = xn(i1,i2,id,j)
               end do
            end do
 
            do i=1,mt+1
               xne(   i,   0,id,j) = xn(   1,     i,idul,j)
               xne(  -1, i+1,id,j) = xn( i-1,     2,idur,j)
               xne(mt+1,   i,id,j) = xn(mt+1-i,  mt,idll,j)
               xne( i-1,mt+2,id,j) = xn(mt-1,mt+2-i,idlr,j)
            end do
 
            xne(   0,   0,id,j) =  xn( 1,   1,id3r,j)
            xne(  -1,   0,id,j) =  xn( 1,   1,id3r,j)
            xne(  -1,   1,id,j) =  xn( 1,   1,id3r,j)
            xne(mt+1,mt+1,id,j) = xne(mt,mt+2,id  ,j)
            xne(mt+1,mt+2,id,j) = xne(mt,mt+2,id  ,j)
 
         end do
 
      end do
 
      end
*dk xnesplit
      subroutine xnesplit(xnep,xne)
 
c...  This routine reads the appropriate portion of the global array
c...  xne into the sub-array xnep needed for process mynum.
 
      include 'size.h'
      include 'pcom.h'
      real xnep(0:nt+2,0:nt+2,nd,3), xne(0:mt+2,0:mt+2,10,3)
 
c...  Determine subdomain limits.
 
      iproc = mod(mynum, 2**(2*lvproc))
      rnm   = real(nt)/real(mt)
      j0    = iproc*rnm
      i0    = (real(iproc)*rnm - j0)/rnm

      ibeg  = i0*nt
      jbeg  = j0*nt
      iend  = ibeg + nt + 2
      jend  = jbeg + nt + 2
 
c...  Load sub-array.
 
      do id=1,nd
 
         jd = id
         if(nd.eq.5 .and. mynum.ge.mproc) jd = id + 5
         i  = -1
 
         do ii=ibeg,iend
 
            i =  i + 1
            j = -1
 
            do jj=jbeg,jend
 
               j = j + 1
 
               xnep(i,j,id,1) = xne(ii,jj,jd,1)
               xnep(i,j,id,2) = xne(ii,jj,jd,2)
               xnep(i,j,id,3) = xne(ii,jj,jd,3)
 
            end do
 
         end do
 
      end do
 
      end
*dk bdopinit
      blockdata bdopinit
 
      common /dpth/ dep(128,4)
      data dep/2*100., 150., 100., 150., 3*100., 200., 250., 300., 350.,
     &     300., 250., 200., 140.,112*0.,
     &	6*50., 75., 2*50., 2*75., 50., 75., 100., 9*150., 2*125., 100.,
     &	2*75., 3*50., 40., 96*0.,
     &     15*30., 10*40., 3*50., 22*60., 4*50., 7*40., 3*30., 64*0.,
     &     30*30., 20*40., 6*50., 44*60., 8*50.,14*40., 6*30./
 
      common /ofst/ j1n(7), j2n(7), md(7)
      data j1n/0,1,0,-1,-1,0,1/, j2n/0,0,1,1,0,-1,-1/, md/1,5,6,7,2,3,4/
 
      end blockdata
