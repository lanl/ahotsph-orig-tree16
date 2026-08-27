c © 2026. Triad National Security, LLC. All rights reserved.
c This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare. derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.

      subroutine abinit(irank,netrcfile)

      implicit none

      include 'dimen'
      include 'comod'

      integer*4 irank
      character*20 netrcfile

      real*8 y(ndim),xx(ndim)

c..   solar system abundance data
      integer*4 iagdim
      parameter( iagdim = 286 )
      integer*4 nzag(iagdim),naag(iagdim)
      real*8 xxag(iagdim)
      character*2 chag(iagdim)
      integer*4 k, i,j,n

      include 'crate'
      include 'caeps'
      include 'cgen'
      include 'cburn'

      real*8 summ, zhe

      integer*4 izbu(ndim),inbu(ndim),ibu,idummy
      character*5 cbu

      integer*4 nsp
      parameter( nsp=20 )
c..Z and N for special nuclei
      integer*4 nspz(nsp), nspn(nsp)
c..   nscr=scratch array for index reordering
c cie: nscr should have N=nucpg elements?
      integer*4 nscr(nucpg),iscr,itno,jz,ja
 
      integer*4 lun

      logical tobe

      integer*4 arrshape1(1), arrshape2(2)

c..default values for new network; used to reset abundances
      data zpop/0.01542d0/,zhyd/0.7095d0/

c..special nuclei to be used in diagnostics and i/o
      data nspz/ 1, 2, 6, 7, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 26, 
     1     26, 28, 28, 28, 28/
      data nspn/ 0, 2, 6, 7, 8, 10, 12, 14, 16, 18, 20, 26, 28, 28, 30,
     1     32, 30, 32, 34, 36/

c-------------------------------------------------------------------

c      write(*,*) netrcfile

      if(irank.eq.0) write(*,*)'ENTERING ABINIT, netsize ',netsize
ccccccccccccccccccccccccccccc

c..   resource file for analysis programs
      inquire(file=trim(netrcfile),exist=tobe)
      if( .not. tobe )then
         if(irank.eq.0) write(*,*)'abinit: no ',trim(netrcfile)
         newnet = 1
      endif

      arrshape1 = shape (cnuc)
      if (arrshape1(1) .lt. ndim) then
          write(*,*)"Error: cnuc has wrong dimensions. Expected ",
     1        ndim, " got ", arrshape1
          call exit(22)
      endif
      arrshape1 = shape (lz)
      if (arrshape1(1) .lt. ndim) then
          write(*,*)"Error: lz has wrong dimensions. Expected ",
     1        ndim, " got ", arrshape1
          call exit(22)
      endif
      arrshape1 = shape (ln)
      if (arrshape1(1) .lt. ndim) then
          write(*,*)"Error: ln has wrong dimensions. Expected ",
     1        ndim, " got ", arrshape1
          call exit(22)
      endif
      if( newnet .ne. 0 )then

         if(irank.eq.0) write(*,*)
     1      "abinit: newnet =",newnet,', redefining abundances'

         do j = 1, netsize
            cnuc(j) = xid(j)
            lz(j) = nz(j)
            ln(j) = nn(j)
         enddo
       
         do j = netsize +2, ndim
            cnuc(j) = '     '
            lz(j) = 0
            ln(j) = 0
         enddo
         cnuc(netsize+1) = '   Ye'

c..   solar  "mass" (nucleon) fractions
         open(10,file='ssystem.dat')
c..naag is the "mass number" A=Z+N, not the mass in amu
         do i = 1, iagdim
            read(10,'(i5,i3,a3,i4,1pe11.3)')idummy,nzag(i),chag(i),
     1           naag(i),xxag(i)
         enddo
         do n = 1, netsize+1
            xx(n) = 0.0d0
         enddo
         do n = 1, netsize
            do i = 1, iagdim
               if( lz(n) .eq. nzag(i) .and.
     1              ln(n)+lz(n) .eq. naag(i))then
                  xx(n) = xxag(i)
               endif
            enddo
         enddo
         arrshape1 = shape (solarx)
         if (arrshape1(1) .lt. netsize) then
             write(*,*)"Error: solarx has wrong dimensions. Expected ",
     1         netsize, " got ", arrshape1
             call exit (22)
         endif
         do n = 1, netsize
            solarx(n) = xx(n)
         enddo
c..   determine which nuclei are stable
c..   chose them for pgplot.f and for dout3.f
c         do j = 1, nucpg
c            nucp(j) = 0
c         enddo
c..   always choose H1 and He4
c         nucp(1) = netsize-1
c         nucp(2) = netsize
c         j = 2
c         do n = 1, netsize
c            if( xx(n) .gt. 0.0d0 .and. j+1 .le. nucpg )then
c               j = j+1
c               nucp(j) = n
c               if(irank.eq.0) write(*,'(5i5,a5,1p8e12.3)')
c     1              j,n,nucp(j),lz(n),ln(n),
c     1              cnuc(n),xx(n)
c            endif
c         enddo

c..   number of nonradioactive nuclei
c         mnucpg = j


c..force values if new network (ag88m)
         zhyd = 0.0d0
         zpop = 0.0d0
         zhe  = 0.0d0
         do n = 1, nnuc
            if( lz(n) .eq. 1 )then
               zhyd = zhyd + solarx(n)
            elseif( lz(n) .eq. 2 )then
               zhe  = zhe  + solarx(n)
            else
               zpop = zpop + solarx(n)
            endif
         enddo

c         write(*,'(a30,2(a5,1pe12.3))')'forcing values:','zpop',zpop,
c     1        'zhyd',zhyd

c..   adjust elemental helium to fit
         do n = 1, netsize
            if( lz(n) .eq. 2 )then
               xx(n) = xx(n)*(1.0d0-zhyd-zpop)/zhe
            endif
         enddo

         summ = 0.0d0
         do n = 1, netsize
            summ = summ + xx(n)
         enddo
c         write(*,'(a30,1pe12.3)')'error after He adjustment:',
c     1        summ-1.0d0

c..   mole fractions
c.. cie: I think this block is unnecessary?!
c         do n = 1, netsize
c            xa(n) = nz(n) + nn(n) + qex(n)/931.487d0
c            y(n) = xx(n)/xa(n)
c         enddo
c         y(netsize+1) = 0.0d0
c         do n = 1, netsize
c            y(netsize+1) =  y(netsize+1) + y(n)*dble(lz(n))
c         enddo
c         if(irank.eq.0) write(*,*)'Ye =',y(netsize+1)

c..   spead over spatial grid 
c         do k = 2, kk+1
c            do n = 1, netsize+1
c               x(n,k) = y(n)
c            enddo
c         enddo

c..   redefine net.rc values
         lun = 30+irank
         open(lun,file=trim(netrcfile))
         do j = 1, netsize
            write(lun,'(3i5,a5,0pf10.4,1pe12.4)')
     1           j,nz(j),nn(j),xid(j),qex(j),solarx(j)
         enddo
c..   write whole array, including zeros, for ease in reading
c..   by many routines
         do j = 1, nucpg, 10
            write(lun,'(10i5)')(nucp(i),i=j,j+9)
         enddo


         close(lun)

         if(irank.eq.0) write(*,*)
     1      'abinit: new network and abundances defined'
         newnet = 0

      else

         if(irank.eq.0) write(*,*)'abinit: checking netrc'

         open(lun,file=trim(netrcfile),status='old')
         ibu = 1
c..initializing qex and solarx
 100     read(lun,'(3i5,a5,0pf10.4,1pe12.4)',end=101)
     1        idummy,izbu(ibu),inbu(ibu),cbu,qex(ibu),solarx(ibu)
         if( izbu(ibu) .ne. nz(ibu) .or.
     1        inbu(ibu) .ne. nn(ibu) )then
c..   different net.rc
            if(irank.eq.0) write(*,*)trim(netrcfile)
            if(irank.eq.0) write(*,'(3i5,a5)')
     1         ibu,izbu(ibu),inbu(ibu),cbu
            if(irank.eq.0) write(*,*)'abinit'
            if(irank.eq.0) write(*,'(3i5,a5)')
     1         ibu,nz(ibu),nn(ibu),xid(ibu)
            if(irank.eq.0) write(*,*)
     1         'abinit: different ',trim(netrcfile)
            if(irank.eq.0) write(*,*)
     1         'rm ',trim(netrcfile),' or change its name, and rerun'
            stop'abinit: reading net.rc'
         endif
         if( izbu(ibu) .eq. 2 .and. inbu(ibu) .eq. 2 )goto 101
         ibu = ibu + 1
         goto 100
 101     continue
c..   net.rc is consistent
         do j = 1, netsize
            cnuc(j) = xid(j)
            lz(j) = nz(j)
            ln(j) = nn(j)
         enddo
       
         do j = netsize +2, ndim
            cnuc(j) = '     '
            lz(j) = 0
            ln(j) = 0
         enddo
         cnuc(netsize+1) = '   Ye'

         read(lun,'(10i5)')nucp

c..   determine nonzero entries
c         mnucpg = 0
c         do j = 1, nucpg
c            if( nucp(j) .ne. 0 )then
c               mnucpg = mnucpg + 1
c            endif
c         enddo

c..define solar metallicity from solar tables for consistency
         zsol = 0.0d0
         do j = 1, ibu
            if( lz(j) .ne. 1 .and. lz(j) .ne. 2 )then
               zsol = zsol + solarx(j)
            endif
         enddo

         if(irank.eq.0) write(*,*)
     1      'abinit: ',trim(netrcfile),' exists and is left unchanged'
         close(lun)
      endif

c..   define Ye
c      do k = 2, kk
         xa(netsize+1) = 0.0d0
         do n = 1, netsize
            xa(netsize+1) = xa(netsize+1) + xa(n)*dble(lz(n))
         enddo
         if( xa(netsize+1) .lt. 0.0d0 
     1        .or. xa(netsize+1) .gt. 1.0d0)then
            if(irank.eq.0) write(*,*)
     1         ' abinit: Ye error, k ', xa(netsize+1)
            stop'abinit: Ye error'
         endif
c      enddo

c..adjust nucp array to control which nuclei are monitored.............

c..find special nuclei
c      j = 0
c      do i = 1, nsp
c         do n = 1, netsize
c            if( nspz(i) .eq. lz(n) .and. nspn(i) .eq. ln(n) )then
c               nucp(i) = n
c               j = j + 1
c               go to 110
c            endif
c         enddo
c         if(irank.eq.0) write(*,'(a8,i5,a3,i5,a3,i5,a15)')'nucleus',i,
c     1        'Z',nspz(i),'N',nspn(i),'not found'
c 110     continue
c      enddo
c      if(irank.eq.0) write(*,*)j
c      if( j .eq. nsp )then
c         do i = 1, nsp
c               if(irank.eq.0) write(*,'(3i5,a5,1p8e12.3)')
c     1            i,nspz(i),nspn(i),
c     1              cnuc(nucp(i))
c         enddo
c         if(irank.eq.0) write(*,*)'ABINIT: All special nuclei found'
c      else
c         if(irank.eq.0) write(*,*)
c     1      'ABINIT: ',nsp-mnucpg,'  special nuclei NOT found'
c      endif

c      j = mnucpg
c..   j is number of special nuclei which were found
c      do n = 1, netsize
c         if( solarx(n) .gt. 0.0d0 .and. j+1 .le. nucpg )then
c            do i = 1, j
c               if( nucp(i) .eq. n )then
c..   avoid duplication
c                  go to 120
c               endif
c            enddo
c            j = j+1
c            nucp(j) = n
c         endif
c 120     continue
c      enddo
c..j has been increased to include some stable isotopes, up to nucpg=60
c      if(irank.eq.0) write(*,*)j,' nuclei chosen'
c      if(irank.eq.0) write(*,'(20a6)')(cnuc(nucp(i)),i=1, j)

c..reorder by charge and atomic number
c      if(irank.eq.0) write(*,*)'begin reordering'
c      do n = 1, j
c         nscr(n) = nucp(n)
c      enddo

c      do j = 1, nucpg !cie: why is this nucpg?? lz/ln are declared in
cburn to have ndim elements
c      do j = 1, ndim
c         itno = 0
c         do n = 1, nucpg-1
c            if( lz(nscr(n)) .gt. lz(nscr(n+1)) )then
c..switch n and n+1 in nscr (index array)
c               iscr = nscr(n)
c               nscr(n) = nscr(n+1)
c               nscr(n+1)   = iscr
c               itno = itno + 1
c            endif
c         enddo
c         if( itno .le. 0 )then
c            if(irank.eq.0) write(*,'(20a6)')(cnuc(nscr(i)),i=1, nucpg)
c            goto 130
c         endif
c      enddo
c 130  continue
c      jz = j-1
c      if(irank.eq.0) write(*,*)'reordered in Z in ',jz,' steps'

c      do j = 1, nucpg
c         itno = 0
c         do n = 1, nucpg-1
c            if(  lz(nscr(n))   + ln(nscr(n)) .gt. 
c     1           lz(nscr(n+1)) + ln(nscr(n+1)) .and.
c     2           lz(nscr(n)) .eq. lz(nscr(n+1)))then
c..switch n and n+1 in nscr (index array)
c               iscr = nscr(n)
c               nscr(n) = nscr(n+1)
c               nscr(n+1)   = iscr
c               itno = itno + 1
c            endif
c         enddo

c         if( itno .le. 0 )then
c            if(irank.eq.0) write(*,'(20a6)')(cnuc(nscr(i)),i=1, nucpg)
c            goto 140
c         endif
c      enddo
c 140  continue
c      ja = j-1
c      if(irank.eq.0) write(*,*)'reordered in A in ',ja,' steps'

c      do n = 1, nucpg
c         nucp(n) = nscr(n)
c      enddo
c      if(irank.eq.0) write(*,*)'replaced original by reordered index'
c..   reordering finished

c      if( jz .ne. 0 .or. ja .ne. 0 )then
c..   redefine net.rc values
c         open(lun,file=trim(netrcfile))
c         do j = 1, netsize
c            read(lun,'(3i5,a5,0pf10.4,1pe12.4)')
c     1           i,nz(j),nn(j),xid(j),qex(j),solarx(j)
c         enddo
c..   overwrite array, including zeros, for ease in reading
c..   by many routines
c         do j = 1, nucpg, 10
c            write(lun,'(10i5)')(nucp(i),i=j,j+9)
c         enddo
c         close(lun)
c         if(irank.eq.0) write(*,*) trim(netrcfile),
c     1      ' adjusted for new nucp index array'
c      endif

      if(irank.eq.0) write(*,*)'LEAVING ABINIT'
      return
      end



