*dk fldsout
	subroutine fldsout(idump,iter)
!	This routine writes out the temperature, velocity, pressure,
!	and plate map fields.

	include 'size.h'
	include 'pcom.h'

	common /temp/ tpb(nt+1), temp((nt+1)**2,nd,nr+1), tpe(nt+1)
	common /velo/ upb(nt+1), u((nt+1)**2,nd,3,nr+1),  upe(nt+1)
	common /pres/ ppb(nt+1), pres((nt+1)**2,nd,nr+1), ppe(nt+1)
	common /mesh/ xn((nt+1)**2,nd,3)
	common /radl/ rshl(nr+1), ird                                      
	common /io01/ casenum, gpath,  lpath, opath, tpath
	common /io02/ idump0, ngpath, nlpath, vtkw_output, cfile_output

	character gpath*80, lpath*80, cname*12, casenum*3, vname*7,cname2*2
	character opath*80
	logical vtkw_output, cfile_output
	integer gen_pvtu, iter
      
!	Output using VTKW functionality
	if(vtkw_output) then

#ifdef USE_VTKW
		call proprty(iter)

     		write(cname2,'(I2.2)') idump
     		vname='v'//casenum//'.'//cname2

		if(nproc>1) then
			gen_pvtu = 1
		else
			gen_pvtu = 0
		endif

		call vtkw_prepare( nproc, mynum, nt, nr, nd,
     &        xn, rshl, gen_pvtu, 0, 1 )

		call vtkw_fields( vname//char(0), 2, 1,
     &           temp, 'temperature'//char(0),
     &           pres, 'pressure'//char(0),
     &           u, 'velocity'//char(0) )

#else 
		if(mynum==0) write(*,2020)
		if(nproc>1) call pexit
		stop

 2020    format(/,' Error:'//' Invalid attempt to use VTKW'
     &        ' functionality! You need to recompile after'/' setting '
     &        'use_vtkw = yes inside the Makefile to enable this '
     &        ' feature!'/)

#endif
	endif


!	Output using 'c-file' format
	if(cfile_output) then

		write(cname, '(''c'',A3,''.'',I4.4,''.'',I2.2)' )
     &        casenum, mynum, idump
        
		iunit = 20
		open (iunit, file=trim(lpath)//cname, status='unknown')
         
		call proprty(iter)
         
		call vecout(temp,rshl,1,nr,nt,nr,iunit,1)
		call vecout(u   ,rshl,3,nr,nt,nr,iunit,0)
		call vecout(pres,rshl,1,nr,nt,nr,iunit,0)

		close (iunit)
	endif

	if((.not.cfile_output).and.(.not.vtkw_output)) then
		print *, 'WARNING: No output done in fldsout!!!'
	endif

	end subroutine



*dk fldsout2
	subroutine fldsout2(idump,iadj)
!	This routine writes out the temperature, velocity, pressure,
!	and plate map fields.

	include 'size.h'
	include 'pcom.h'

	common /temp/ tpb(nt+1), temp((nt+1)**2,nd,nr+1), tpe(nt+1)
	common /velo/ upb(nt+1), u((nt+1)**2,nd,3,nr+1),  upe(nt+1)
	common /pres/ ppb(nt+1), pres((nt+1)**2,nd,nr+1), ppe(nt+1)
	common /mesh/ xn((nt+1)**2,nd,3)
	common /radl/ rshl(nr+1), ird                                      
	common /io01/ casenum, gpath,  lpath, opath, tpath
	common /io02/ idump0, ngpath, nlpath, vtkw_output, cfile_output

	character gpath*80, lpath*80, cname*15, casenum*3, vname*7,cname2*2
	character opath*80
	logical vtkw_output, cfile_output
	integer gen_pvtu
	integer iadj
      
!	Output using VTKW functionality
	if(vtkw_output) then

#ifdef USE_VTKW
		call proprty(iadj)

     		write(cname2,'(I2.2)') idump
     		vname='v'//casenum//'.'//cname2

		if(nproc>1) then
			gen_pvtu = 1
		else
			gen_pvtu = 0
		endif

		call vtkw_prepare( nproc, mynum, nt, nr, nd,
     &        xn, rshl, gen_pvtu, 0, 1 )

		call vtkw_fields( vname//char(0), 2, 1,
     &           temp, 'temperature'//char(0),
     &           pres, 'pressure'//char(0),
     &           u, 'velocity'//char(0) )

#else 
		if(mynum==0) write(*,2020)
		if(nproc>1) call pexit
		stop

 2020    format(/,' Error:'//' Invalid attempt to use VTKW'
     &        ' functionality! You need to recompile after'/' setting '
     &        'use_vtkw = yes inside the Makefile to enable this '
     &        ' feature!'/)

#endif
	endif


!	Output using 'c-file' format
	if(cfile_output) then

		write(cname, '(''c'',A3,''.'',I4.4,''.'',I2.2,"_",I2.2)' )
     &        casenum, mynum, idump, iadj-1
        
		iunit = 20
		open (iunit, file=trim(opath)//'c-files/'//cname,
     &			 status='unknown')
         
		call proprty(iadj)
         
		call vecout(temp,rshl,1,nr,nt,nr,iunit,1)
		call vecout(u   ,rshl,3,nr,nt,nr,iunit,0)
		call vecout(pres,rshl,1,nr,nt,nr,iunit,0)

		close (iunit)
	endif

	if((.not.cfile_output).and.(.not.vtkw_output)) then
		print *, 'WARNING: No output done in fldsout!!!'
	endif

	end subroutine
	
	
	
*dk check_vtkw_usage
C> This function reads the value of the exptype parameter from the inputfile
C> and returns a logical value that specifies whether we will use vtkw output
C> functionality or not.
C> For backward compatibility it returns .false. if the parameter is not
C> present. However, this works only, if the interra file ends directly
C> after the tmpmlt line.
C> \params unit unit descriptor from which to read data
      subroutine get_output_types( unit, vtkw_output, cfile_output )
      implicit none

      include 'pcom.h'

C     formal parameters
      integer unit
      logical vtkw_output, cfile_output

C     local variables
      integer sinit, ios
      logical valid_value_found
      character string*12, blanks*12, line*80

      blanks = '            '

      valid_value_found = .false.

      read(5,'(A)', iostat=ios) line

      if ( ios .lt. 0 .or. index(line,'exptype') .eq. 0 ) then

        vtkw_output = .false.
        cfile_output = .true.

      else if ( ios .eq. 0 ) then

         string = line(1:12)
         sinit = index( string, 'vtkw' )

         if ( sinit .gt. 0 ) then
            if ( string(sinit:sinit+3) .eq. 'vtkw'
     &           .and. string(1:sinit-1) .eq. blanks(1:sinit-1)
     &           .and. string(sinit+4:) .eq. blanks(sinit+4:) ) then
               vtkw_output = .true.
               cfile_output = .false.
               valid_value_found = .true.
            endif
         endif

         sinit = index( string, 'c-file' )
         if ( sinit .gt. 0 ) then
            if ( string(sinit:sinit+5) .eq. 'c-file'
     &           .and. string(1:sinit-1) .eq. blanks(1:sinit-1)
     &           .and. string(sinit+6:) .eq. blanks(sinit+6:) ) then
               vtkw_output = .false.
               cfile_output = .true.
               valid_value_found = .true.
            endif
         endif

         sinit = index( string, 'both' )
         if ( sinit .gt. 0 ) then
            if ( string(sinit:sinit+3) .eq. 'both'
     &           .and. string(1:sinit-1) .eq. blanks(1:sinit-1)
     &           .and. string(sinit+4:) .eq. blanks(sinit+4:) ) then
               vtkw_output = .true.
               cfile_output = .true.
               valid_value_found = .true.
            endif
         endif

         if ( valid_value_found .neqv. .true. ) then
            write(*,2000) string
 2000       format(/,' Error:'//' Invalid value ''',A,
     &           ''' for parameter ''exptype''!'/
     &           ' Only ''vtkw'', ''c-file'' and ''both'' allowed!'/ )
            if(nproc .gt. 1) call pexit
            stop
         endif

      else
         print *,'Read from unit', unit, 'failed with IOSTAT=', ios
         stop
      endif

#ifndef USE_VTKW
      if ( vtkw_output ) then
         if(mynum .eq. 0) write(*,2010)
         if(nproc .gt. 1) call pexit
         stop
      endif

 2010       format(/,' Error:'//' Found option ''vtkw'' for parameter '
     &     '''exptype''. You need to recompile after'/' setting '
     &     'use_vtkw = yes inside the Makfile to enable this feature!'/)

#endif
       
      end

C -----------------------------------------------------------------------------
*dk vecout
      subroutine vecout(u,rshl,nj,kr,kt,mr,nf,ifmt)
 
c...  This routine writes the nodal field u to logical unit nf
c...  using 1pe10.3 format when ifmt = 0 and f10.3 when ifmt = 1.
 
      include 'size.h'
      real u(*), rshl(*)
      common /name/ titl(4,4)
      common /prty/ propr(20)
      character*8 titl

      write(nf,10) kr, kt
 10   format(2i5)
 
      write(nf,20) titl
 20   format(4a8)
 
      write(nf,30) (rshl(i),i=1,kr+1)
      write(nf,30) propr
 30   format(1p10e15.8)
 
      if(ifmt .eq. 0) then
 
         write(nf,40) (u(ii),ii=1,(kt+1)**2*nd*nj*(mr+1))
 40      format(1p15e10.3)
 
      elseif(ifmt .eq. 1) then
 
         write(nf,50) (u(ii),ii=1,(kt+1)**2*nd*nj*(mr+1))
 50      format(15f10.3)
 
      endif
 
      end

C -----------------------------------------------------------------------------
*dk vecoutgmt
      subroutine vecoutgmt(s,f,rshl,kr,kt,nf,ifmt)
 
c...  This routine writes the nodal field s to logical unit nf
c...  using 1pe10.3 format when ifmt = 0 and f10.3 when ifmt = 1.
c...  Output is written as an x/y/data array for 2-D gmt-plots.
 
      include 'size.h'
      real s((kt+1)**2*nd,kr+1), f((kt+1)**2*nd,2)
      real rshl(kr+1)
      common /name/ titl(4,4)
      common /prty/ propr(20)
      character*8 titl
 
      if(ifmt .eq. 0) then
 
         do ii=1,(kt+1)**2*nd
            write(nf,40) f(ii,2),f(ii,1),s(ii,5)
         end do
 
 40      format(3(1pe10.3,3x))
 
      elseif(ifmt .eq. 1) then
 
         do ii=1,(kt+1)**2*nd
            write(nf,50) f(ii,2),f(ii,1),s(ii,5)
         end do
 
 50      format(3(f10.3,3x))
 
      endif
 
      end

C -----------------------------------------------------------------------------
*dk gmtout
      subroutine gmtout(idump,iter)
 
c     This routine is a driver to write out temperature and velocity
c     fields in a format that can be used by the graphics tool GMT.
 
      include 'size.h'
      include 'pcom.h'

      real rsh(nr/2+1)
c     common /flds/ f((nt+1)**2*nd,3,(nr+1),2)
      common /flds/ f((nt+1)**2*nd,3,(nr+1),2), fdummy(nv*2)
      common /velo/ upb(nt+1), u((nt+1)**2*nd,3,(nr+1)),  upe(nt+1)
      common /temp/ tpb(nt+1), temp((nt+1)**2*nd*(nr+1)), tpe(nt+1)
      common /pres/ ppb(nt+1), pres((nt+1)**2*nd*(nr+1)), ppe(nt+1)
      common /shht/ shb(nt+1),  shr((nt+1)**2*nd*(nr+1)), she(nt+1)
      common /strn/ srate((nt+1)**2*nd*(nr+1))
      common /radl/ rshl(nr+1), ird
      common /mesh/ xn((nt+1)**2,nd,3)
      common /prop/ ieos,  rho0, visc,  grav,  texpn, tcond,
     &              sheat, hgen, tb(2), cl410, cl660
      common /vis0/ vb(nt+1), vv((nt+1)**2*nd*(nr+1),2), ve(nt+1)
      common /vis1/ vscmax, rvscscl, tvscscl, pwrlawn, pwrlawsr, yldstrs
      common /vis2/ rdvsc(nr+1), tactv(nr+1), vscl(nr+1)
      common /io01/ casenum, gpath,  lpath, opath, tpath
      common /io02/ idump0, ngpath, nlpath, vtkw_output, cfile_output

	integer iter
      logical vtkw_output, cfile_output
      character gpath*80, lpath*80, casenum*3
      character char1*1, char2*1, char3*1, char4*1, char5*1, char6*1
      character gmtpath*6, name1*19, name2*19, name3*19, opath*80
 
c...  place gmt output into the local GMT directory
      gmtpath  = './GMT/'
 
      char1 = char(48 + mynum/1000)
      char2 = char(48 + mod(mynum/100,10))
      char3 = char(48 + mod(mynum/10,10))
      char4 = char(48 + mod(mynum,10))
      char5 = char(48 + idump/10)
      char6 = char(48 + mod(idump,10))
      name1 = 'gmt'//casenum//'.'//char1//char2//char3//char4
     &             //'.'//char5//char6//'.temp'
      name2 = 'gmt'//casenum//'.'//char1//char2//char3//char4
     &             //'.'//char5//char6//'.vec1'
      name3 = 'gmt'//casenum//'.'//char1//char2//char3//char4
     &             //'.'//char5//char6//'.vec2'
      iunit = 20
 
      call proprty(iter)
      call horizontvel(f(1,1,1,1),u)
      call uthetaphi(f(1,1,1,2),f(1,1,1,1),xn)
      call xntolatlong(f(1,1,1,1),xn)
 
      open (iunit, file=gmtpath//name1, status='unknown')
 
c     if(mt .ge. 128) then
c        do ir=1,nr/2+1
c           rsh(ir) = rshl(ir+ir-1)
c        end do
c        call sample(f,temp,1,nr,nt)
c        call vecoutgmt(f(1,4),f(1,1),rsh,nr/2,nt/2,nr/2,iunit,1)
c     else
         call vecoutgmt(temp,f(1,1,1,1),rshl,nr,nt,iunit,1)
c     end if
 
      close (iunit)
      open  (iunit, file=gmtpath//name2, status='unknown')
 
c     if(mt .ge. 128) then
c        call sample(f,u,3,nr,nt)
c        call vecoutgmt(f,rsh,nr/2,nt/2,nr/2,iunit,0)
c     else
         call vecoutgmt(f(1,1,1,2),f(1,1,1,1),rshl,nr,nt,iunit,0)
c     end if
 
      close (iunit)
      open  (iunit, file=gmtpath//name3, status='unknown')
 
c     if(mt .ge. 128) then
c        call sample(f,u,3,nr,nt)
c        call vecoutgmt(f,rsh,nr/2,nt/2,nr/2,iunit,0)
c     else
         call vecoutgmt(f(1,2,1,2),f(1,1,1,1),rshl,nr,nt,iunit,0)
c     end if
 
      close (iunit)
 
      end


C -----------------------------------------------------------------------------
*dk graphout
      subroutine graphout(idump,iter)
 
c     This routine writes out the temperature, pressure, velocity,
c     shear heating, viscosity, and plate map fields.
 
      include 'size.h'
      include 'pcom.h'

      real rsh(nr/2+1)
      common /flds/ f((nt+1)**2*nd*(nr+1),8)
      common /velo/ upb(nt+1), u((nt+1)**2*nd*3*(nr+1)),  upe(nt+1)
      common /temp/ tpb(nt+1), temp((nt+1)**2*nd*(nr+1)), tpe(nt+1)
      common /pres/ ppb(nt+1), pres((nt+1)**2*nd,(nr+1)), ppe(nt+1)
      common /shht/ shb(nt+1),  shr((nt+1)**2*nd*(nr+1)), she(nt+1)
      common /strn/ srate((nt+1)**2*nd*(nr+1))
      common /topo/ h((nt+1)**2*nd,3), g((nt+1)**2*nd,2)
      common /radl/ rshl(nr+1), ird
      common /prop/ ieos,  rho0, visc,  grav,  texpn, tcond,
     &              sheat, hgen, tb(2), cl410, cl660
      common /vis0/ vb(nt+1), vv((nt+1)**2*nd*(nr+1),2), ve(nt+1)
      common /vis1/ vscmax, rvscscl, tvscscl, pwrlawn, pwrlawsr, yldstrs
      common /vis2/ rdvsc(nr+1), tactv(nr+1), vscl(nr+1)
      common /io01/ casenum, gpath,  lpath, opath, tpath
      common /io02/ idump0, ngpath, nlpath, vtkw_output, cfile_output

	integer iter
      logical vtkw_output, cfile_output
      character gpath*80, lpath*80, cname*12, casenum*3, opath*80
 
      write( cname, '(''g'',A3,''.'',I4.4,''.'',I2.2)' )
     &     casenum, mynum, idump

      iunit = 20
 
      call proprty(iter)
 
      open (iunit, file=trim(lpath)//cname, status='unknown')
 
      do ir=1,nr/2+1
         rsh(ir) = rshl(ir+ir-1)
      end do
 
c     call topogrph
 
c...  Write geoid field to graphics file.
c
c     if(mt .ge. 128) then
c        call sample(f,g,1,0,nt)
c        call vecout(f,rsh,1,nr/2,nt/2,0,iunit,0)
c     else
c        call vecout(g,rshl,1,nr,nt,0,iunit,0)
c     end if
c
c...  Write topography fields to graphics file.
c
c     if(mt .ge. 128) then
c        call sample(f,h,2,0,nt)
c        call vecout(f,rsh,1,nr/2,nt/2,1,iunit,0)
c     else
c        call vecout(h,rshl,1,nr,nt,1,iunit,0)
c     end if
c
c...  Write plate map field to graphics file.
c
c     if(nplate .gt. 0) then
c        if(mt .ge. 128) then
c           call sample(f,pmap,2,0,nt)
c           call vecout(f,rsh,1,nr/2,nt/2,1,iunit,1)
c        else
c           call vecout(pmap,rshl,1,nr,nt,1,iunit,1)
c        end if
c     end if
 
c...  Write temperature field to graphics file.
 
      if(mt .ge. 128) then
         call sample(f,temp,1,nr,nt)
         call vecout(f,rsh,1,nr/2,nt/2,nr/2,iunit,1)
      else
         call vecout(temp,rshl,1,nr,nt,nr,iunit,1)
      end if
 
c...  Write velocity field to graphics file.
 
      if(mt .ge. 128) then
         call sample(f,u,3,nr,nt)
         call vecout(f,rsh,3,nr/2,nt/2,nr/2,iunit,0)
      else
         call vecout(u,rshl,3,nr,nt,nr,iunit,0)
      end if
 
c...  Write viscosity field to graphics file.
 
      if(vscmax.gt.0. .and. tvscscl.gt.0.) then
 
         lvf = 1.45*log(real(mt))
         do ir=1,nr+1
            jj = (ir-1)*(nt+1)**2*nd
            do ii=1,(nt+1)**2*nd
               if(vv(ii+jj,1) .gt. 0.)
     &            f(ii+jj,1) = log10(visc*vv(ii+jj,1)**2)
            end do
         end do
 
         if(mt .ge. 128) then
            call sample(f(1,2),f,1,nr,nt)
            call vecout(f(1,2),rsh,1,nr/2,nt/2,nr/2,iunit,0)
         else
            call vecout(f,rshl,1,nr,nt,nr,iunit,0)
         end if
 
      end if
 
c...  Write shear heating rate field to graphics file.
 
      if(ieos .ge. 10) then
         if(mt .ge. 128) then
            call sample(f,shr,1,nr,nt)
            call vecout(f,rsh,1,nr/2,nt/2,nr/2,iunit,0)
         else
            call vecout(shr,rshl,1,nr,nt,nr,iunit,0)
         end if
      end if
 
c...  Write strain rate field to graphics file.
 
      if(mt .ge. 128) then
         call sample(f,srate,1,nr,nt)
         call vecout(f,rsh,1,nr/2,nt/2,nr/2,iunit,0)
      else
         call vecout(srate,rshl,1,nr,nt,nr,iunit,0)
      end if
 
      close (iunit)
 
      end

C -----------------------------------------------------------------------------
*dk radialout
      subroutine radialout
 
c...  This routine writes out the radial profiles and other information
c...  needed to plot these profiles.
 
      include 'size.h'
      common /name/ titl(4,4)
      common /prop/ ieos,  rho0, visc,  grav,  texpn, tcond,
     &              sheat, hgen, tb(2), cl410, cl660
      common /heat/ htop, hbot, hrad, heat, hshr, hvol, hnet, tnrm,
     &              tav(nr+1), qc(nr)
      common /radl/ rshl(nr+1), ird
      common /eos1/ rhorf(nr+1), tmprf(nr+1), drhdp(nr+1), grv(nr+1)
      common /eos2/ rhocv(nr+1), gamma(nr+1), alpha(nr+1), prf(nr+1),
     &               cond(nr+1)
      common /vis2/ rdvsc(nr+1), tactv(nr+1), vscl(nr+1)
      common /shrh/ shrlayr(nr+1)
      common /urut/ ur(nr+1), ut(nr+1)
      common /io01/ casenum, gpath,  lpath, opath, tpath
      common /io02/ idump0, ngpath, nlpath, vtkw_output, cfile_output
      
      logical vtkw_output, cfile_output
      character casenum*3, gpath*80, lpath*80, cname*9
      character titl*8, opath*80
 
      cname = 'radial'//casenum
 
      open (33, file=cname, status='unknown')
 
      write(33, '(4a8)') titl
 
      write(33, '(/)')
 
      write(33, '(i5,1p5e12.5)') ieos, rho0, visc, grav, texpn, tcond
      write(33, '(5x,1p5e12.5)') sheat, hgen, tb, 0
 
      write(33, '(/"     rshl       rhorf       tmprf   ",
     &             "     tav         prf         grv "/)')
 
      do ir=1,nr+1
         write(33, '(1p6e12.5)')  rshl(ir), rhorf(ir), tmprf(ir),
     &                             tav(ir),   prf(ir),   grv(ir)
      end do
 
      write(33, '(/"     rshl         ur          ut    ",
     &             "      qc       shrlayr      rdvsc"/)')
 
      do ir=1,nr+1
         qcnd = hbot
         if(ir .le. nr) qcnd = qc(ir)
         write(33, '(1p6e12.5)')  rshl(ir),    ur(ir),    ut(ir),
     &                            qcnd, shrlayr(ir), rdvsc(ir)
      end do
 
      write(33, '(/"     rshl        cond       alpha   ",
     &             "    gamma       rhocv       drhdp"/)')
 
      do ir=1,nr+1
         write(33, '(1p6e12.5)')  rshl(ir),  cond(ir), alpha(ir),
     &                           gamma(ir), rhocv(ir), drhdp(ir)
      end do
 
      close (33)
 
      end
C -----------------------------------------------------------------------------
