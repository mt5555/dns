C
C $Id: fft99.f,v 1.2 2001-06-28 22:58:57 mt Exp $
C
      subroutine fft99(a,work,trigs,ifax,inc,jump,n,lot,isign)
      implicit real*8 (a-h), real*8 (o-z)
c
c purpose      performs multiple fast fourier transforms.  this package
c              will perform a number of simultaneous real/half-complex
c              periodic fourier transforms or corresponding inverse
c              transforms, i.e.  given a set of real data vectors, the
c              package returns a set of 'half-complex' fourier
c              coefficient vectors, or vice versa.  the length of the
c              transforms must be an even number greater than 4 that has
c              no other factors except possibly powers of 2, 3, and 5.
c              this is an all-fortran version of a optimized routine
c              fft99 written for xmp/ymps by dr. clive temperton of
c              ecmwf.
c
c              the package fft99f contains several user-level routines:
c
c            subroutine set99
c                an initialization routine that must be called once
c                before a sequence of calls to the fft routines
c                (provided that n is not changed).
c
c            subroutines fft99 and fft991
c                two fft routines that return slightly different
c                arrangements of the data in gridpoint space.
c
c usage        let n be of the form 2**p * 3**q * 5**r, where p .ge. 1,
c              q .ge. 0, and r .ge. 0.  then a typical sequence of
c              calls to transform a given set of real vectors of length
c              n to a set of 'half-complex' fourier coefficient vectors
c              of length n is
c
c                   dimension ifax(13),trigs(3*n/2+1),a(m*(n+2)),
c                  +          work(m*(n+1))
c
c                   call set99 (trigs, ifax, n)
c                   call fft99 (a,work,trigs,ifax,inc,jump,n,m,isign)
c
c              see the individual write-ups for set99, fft99, and
c              fft991 below, for a detailed description of the
c              arguments.
c
c history      the package was written by clive temperton at ecmwf in
c              november, 1978.  it was modified, documented, and tested
c              for ncar by russ rew in september, 1980.
c
c-----------------------------------------------------------------------
c
c subroutine set99 (trigs, ifax, n)
c
c purpose      a set-up routine for fft99 and fft991.  it need only be
c              called once before a sequence of calls to the fft
c              routines (provided that n is not changed).
c
c argument     ifax(13),trigs(3*n/2+1)
c dimensions
c
c arguments
c
c on input     trigs
c               a floating point array of dimension 3*n/2 if n/2 is
c               even, or 3*n/2+1 if n/2 is odd.
c
c              ifax
c               an integer array.  the number of elements actually used
c               will depend on the factorization of n.  dimensioning
c               ifax for 13 suffices for all n less than a million.
c
c              n
c               an even number greater than 4 that has no prime factor
c               greater than 5.  n is the length of the transforms (see
c               the documentation for fft99 and fft991 for the
c               definitions of the transforms).
c
c on output    ifax
c               contains the factorization of n/2.  ifax(1) is the
c               number of factors, and the factors themselves are stored
c               in ifax(2),ifax(3),...  if set99 is called with n odd,
c               or if n has any prime factors greater than 5, ifax(1)
c               is set to -99.
c
c              trigs
c               an array of trigonometric function values subsequently
c               used by the fft routines.
c
c-----------------------------------------------------------------------
c
c subroutine fft991 (a,work,trigs,ifax,inc,jump,n,m,isign)
c                       and
c subroutine fft99 (a,work,trigs,ifax,inc,jump,n,m,isign)
c
c purpose      perform a number of simultaneous real/half-complex
c              periodic fourier transforms or corresponding inverse
c              transforms, using ordinary spatial order of gridpoint
c              values (fft991) or explicit cyclic continuity in the
c              gridpoint values (fft99).  given a set
c              of real data vectors, the package returns a set of
c              'half-complex' fourier coefficient vectors, or vice
c              versa.  the length of the transforms must be an even
c              number that has no other factors except possibly powers
c              of 2, 3, and 5.  this is an all-fortran version of 
c              optimized routine fft991 written for xmp/ymps by
c              dr. clive temperton of ecmwf.
c
c argument     a(m*(n+2)), work(m*(n+1)), trigs(3*n/2+1), ifax(13)
c dimensions
c
c arguments
c
c on input     a
c               an array of length m*(n+2) containing the input data
c               or coefficient vectors.  this array is overwritten by
c               the results.
c
c              work
c               a work array of dimension m*(n+1)
c
c              trigs
c               an array set up by set99, which must be called first.
c
c              ifax
c               an array set up by set99, which must be called first.
c
c              inc
c               the increment (in words) between successive elements of
c               each data or coefficient vector (e.g.  inc=1 for
c               consecutively stored data).
c
c              jump
c               the increment (in words) between the first elements of
c               successive data or coefficient vectors.  on crays, 
c               try to arrange data so that jump is not a multiple of 8
c               (to avoid memory bank conflicts).  for clarification of
c               inc and jump, see the examples below.
c
c              n
c               the length of each transform (see definition of
c               transforms, below).
c
c              m
c               the number of transforms to be done simultaneously.
c
c              isign
c               = +1 for a transform from fourier coefficients to
c                    gridpoint values.
c               = -1 for a transform from gridpoint values to fourier
c                    coefficients.
c
c on output    a
c               if isign = +1, and m coefficient vectors are supplied
c               each containing the sequence:
c
c               a(0),b(0),a(1),b(1),...,a(n/2),b(n/2)  (n+2 values)
c
c               then the result consists of m data vectors each
c               containing the corresponding n+2 gridpoint values:
c
c               for fft991, x(0), x(1), x(2),...,x(n-1),0,0.
c               for fft99, x(n-1),x(0),x(1),x(2),...,x(n-1),x(0).
c                   (explicit cyclic continuity)
c
c               when isign = +1, the transform is defined by:
c                 x(j)=sum(k=0,...,n-1)(c(k)*exp(2*i*j*k*pi/n))
c                 where c(k)=a(k)+i*b(k) and c(n-k)=a(k)-i*b(k)
c                 and i=sqrt (-1)
c
c               if isign = -1, and m data vectors are supplied each
c               containing a sequence of gridpoint values x(j) as
c               defined above, then the result consists of m vectors
c               each containing the corresponding fourier cofficients
c               a(k), b(k), 0 .le. k .le n/2.
c
c               when isign = -1, the inverse transform is defined by:
c                 c(k)=(1/n)*sum(j=0,...,n-1)(x(j)*exp(-2*i*j*k*pi/n))
c                 where c(k)=a(k)+i*b(k) and i=sqrt(-1)
c
c               a call with isign=+1 followed by a call with isign=-1
c               (or vice versa) returns the original data.
c
c               note: the fact that the gridpoint values x(j) are real
c               implies that b(0)=b(n/2)=0.  for a call with isign=+1,
c               it is not actually necessary to supply these zeros.
c
c examples      given 19 data vectors each of length 64 (+2 for explicit
c               cyclic continuity), compute the corresponding vectors of
c               fourier coefficients.  the data may, for example, be
c               arranged like this:
c
c first data   a(1)=    . . .                a(66)=             a(70)
c vector       x(63) x(0) x(1) x(2) ... x(63) x(0)  (4 empty locations)
c
c second data  a(71)=   . . .                                  a(140)
c vector       x(63) x(0) x(1) x(2) ... x(63) x(0)  (4 empty locations)
c
c               and so on.  here inc=1, jump=70, n=64, m=19, isign=-1,
c               and fft99 should be used (because of the explicit cyclic
c               continuity).
c
c               alternatively the data may be arranged like this:
c
c                first         second                          last
c                data          data                            data
c                vector        vector                          vector
c
c                 a(1)=         a(2)=                           a(19)=
c
c                 x(63)         x(63)       . . .               x(63)
c        a(20)=   x(0)          x(0)        . . .               x(0)
c        a(39)=   x(1)          x(1)        . . .               x(1)
c                  .             .                               .
c                  .             .                               .
c                  .             .                               .
c
c               in which case we have inc=19, jump=1, and the remaining
c               parameters are the same as before.  in either case, each
c               coefficient vector overwrites the corresponding input
c               data vector.
c
c-----------------------------------------------------------------------
      dimension a(n),work(n),trigs(n),ifax(1)
c
c     subroutine "fft99" - multiple fast real periodic transform
c     corresponding to old scalar routine fft9
c     procedure used to convert to half-length complex transform
c     is given by cooley, lewis and welch (j. sound vib., vol. 12
c     (1970), 315-337)
c
c     a is the array containing input and output data
c     work is an area of size (n+1)*lot
c     trigs is a previously prepared list of trig function values
c     ifax is a previously prepared list of factors of n/2
c     inc is the increment within each data 'vector'
c         (e.g. inc=1 for consecutively stored data)
c     jump is the increment between the start of each data vector
c     n is the length of the data vectors
c     lot is the number of data vectors
c     isign = +1 for transform from spectral to gridpoint
c           = -1 for transform from gridpoint to spectral
c
c     ordering of coefficients:
c         a(0),b(0),a(1),b(1),a(2),b(2),...,a(n/2),b(n/2)
c         where b(0)=b(n/2)=0; (n+2) locations required
c
c     ordering of data:
c         x(n-1),x(0),x(1),x(2),...,x(n),x(0)
c         i.e. explicit cyclic continuity; (n+2) locations required
c
c     vectorization is achieved on cray by doing the transforms in
c     parallel
c
c     *** n.b. n is assumed to be an even number
c
c     definition of transforms:
c     -------------------------
c
c     isign=+1: x(j)=sum(k=0,...,n-1)(c(k)*exp(2*i*j*k*pi/n))
c         where c(k)=a(k)+i*b(k) and c(n-k)=a(k)-i*b(k)
c
c     isign=-1: a(k)=(1/n)*sum(j=0,...,n-1)(x(j)*cos(2*j*k*pi/n))
c               b(k)=-(1/n)*sum(j=0,...,n-1)(x(j)*sin(2*j*k*pi/n))
c
c
c
c
      nfax=ifax(1)
      nx=n+1
      nh=n/2
      ink=inc+inc
      if (isign.eq.+1) go to 30
c
c     if necessary, transfer data to work area
      igo=50
      if (mod(nfax,2).eq.1) goto 40
      ibase=inc+1
      jbase=1
      do 20 l=1,lot
      i=ibase
      j=jbase
cdir$ ivdep
      do 10 m=1,n
      work(j)=a(i)
      i=i+inc
      j=j+1
   10 continue
      ibase=ibase+jump
      jbase=jbase+nx
   20 continue
c
      igo=60
      go to 40
c
c     preprocessing (isign=+1)
c     ------------------------
c
   30 continue
      call fft99a(a,work,trigs,inc,jump,n,lot)
      igo=60
c
c     complex transform
c     -----------------
c
   40 continue
      ia=inc+1
      la=1
      do 80 k=1,nfax
      if (igo.eq.60) go to 60
   50 continue
      call vpassm(a(ia),a(ia+inc),work(1),work(2),trigs,
     *   ink,2,jump,nx,lot,nh,ifax(k+1),la)
      igo=60
      go to 70
   60 continue
      call vpassm(work(1),work(2),a(ia),a(ia+inc),trigs,
     *    2,ink,nx,jump,lot,nh,ifax(k+1),la)
      igo=50
   70 continue
      la=la*ifax(k+1)
   80 continue
c
      if (isign.eq.-1) go to 130
c
c     if necessary, transfer data from work area
      if (mod(nfax,2).eq.1) go to 110
      ibase=1
      jbase=ia
      do 100 l=1,lot
      i=ibase
      j=jbase
cdir$ ivdep
      do 90 m=1,n
      a(j)=work(i)
      i=i+1
      j=j+inc
   90 continue
      ibase=ibase+nx
      jbase=jbase+jump
  100 continue
c
c     fill in cyclic boundary points
  110 continue
      ia=1
      ib=n*inc+1
cdir$ ivdep
      do 120 l=1,lot
      a(ia)=a(ib)
      a(ib+inc)=a(ia+inc)
      ia=ia+jump
      ib=ib+jump
  120 continue
      go to 140
c
c     postprocessing (isign=-1):
c     --------------------------
c
  130 continue
      call fft99b(work,a,trigs,inc,jump,n,lot)
c
  140 continue
      return
      end
      subroutine fft99a(a,work,trigs,inc,jump,n,lot)
      implicit real*8 (a-h), real*8 (o-z)
      dimension a(n),work(n),trigs(n)
c
c     subroutine fft99a - preprocessing step for fft99, isign=+1
c     (spectral to gridpoint transform)
c
      nh=n/2
      nx=n+1
      ink=inc+inc
c
c     a(0) and a(n/2)
      ia=1
      ib=n*inc+1
      ja=1
      jb=2
cdir$ ivdep
      do 10 l=1,lot
      work(ja)=a(ia)+a(ib)
      work(jb)=a(ia)-a(ib)
      ia=ia+jump
      ib=ib+jump
      ja=ja+nx
      jb=jb+nx
   10 continue
c
c     remaining wavenumbers
      iabase=2*inc+1
      ibbase=(n-2)*inc+1
      jabase=3
      jbbase=n-1
c
      do 30 k=3,nh,2
      ia=iabase
      ib=ibbase
      ja=jabase
      jb=jbbase
      c=trigs(n+k)
      s=trigs(n+k+1)
cdir$ ivdep
      do 20 l=1,lot
      work(ja)=(a(ia)+a(ib))-
     *    (s*(a(ia)-a(ib))+c*(a(ia+inc)+a(ib+inc)))
      work(jb)=(a(ia)+a(ib))+
     *    (s*(a(ia)-a(ib))+c*(a(ia+inc)+a(ib+inc)))
      work(ja+1)=(c*(a(ia)-a(ib))-s*(a(ia+inc)+a(ib+inc)))+
     *    (a(ia+inc)-a(ib+inc))
      work(jb+1)=(c*(a(ia)-a(ib))-s*(a(ia+inc)+a(ib+inc)))-
     *    (a(ia+inc)-a(ib+inc))
      ia=ia+jump
      ib=ib+jump
      ja=ja+nx
      jb=jb+nx
   20 continue
      iabase=iabase+ink
      ibbase=ibbase-ink
      jabase=jabase+2
      jbbase=jbbase-2
   30 continue
c
      if (iabase.ne.ibbase) go to 50
c     wavenumber n/4 (if it exists)
      ia=iabase
      ja=jabase
cdir$ ivdep
      do 40 l=1,lot
      work(ja)=2*a(ia)
      work(ja+1)=-2*a(ia+inc)
      ia=ia+jump
      ja=ja+nx
   40 continue
c
   50 continue
      return
      end
      subroutine fft99b(work,a,trigs,inc,jump,n,lot)
      implicit real*8 (a-h), real*8 (o-z)
      dimension work(n),a(n),trigs(n)
c
c     subroutine fft99b - postprocessing step for fft99, isign=-1
c     (gridpoint to spectral transform)
c
      nh=n/2
      nx=n+1
      ink=inc+inc
c
c     a(0) and a(n/2)
      scale=1d0/n
      ia=1
      ib=2
      ja=1
      jb=n*inc+1
cdir$ ivdep
      do 10 l=1,lot
      a(ja)=scale*(work(ia)+work(ib))
      a(jb)=scale*(work(ia)-work(ib))
      a(ja+inc)=0
      a(jb+inc)=0
      ia=ia+nx
      ib=ib+nx
      ja=ja+jump
      jb=jb+jump
   10 continue
c
c     remaining wavenumbers
      scale=0.5*scale
      iabase=3
      ibbase=n-1
      jabase=2*inc+1
      jbbase=(n-2)*inc+1
c
      do 30 k=3,nh,2
      ia=iabase
      ib=ibbase
      ja=jabase
      jb=jbbase
      c=trigs(n+k)
      s=trigs(n+k+1)
cdir$ ivdep
      do 20 l=1,lot
      a(ja)=scale*((work(ia)+work(ib))
     *   +(c*(work(ia+1)+work(ib+1))+s*(work(ia)-work(ib))))
      a(jb)=scale*((work(ia)+work(ib))
     *   -(c*(work(ia+1)+work(ib+1))+s*(work(ia)-work(ib))))
      a(ja+inc)=scale*((c*(work(ia)-work(ib))-s*(work(ia+1)+work(ib+1)))
     *    +(work(ib+1)-work(ia+1)))
      a(jb+inc)=scale*((c*(work(ia)-work(ib))-s*(work(ia+1)+work(ib+1)))
     *    -(work(ib+1)-work(ia+1)))
      ia=ia+nx
      ib=ib+nx
      ja=ja+jump
      jb=jb+jump
   20 continue
      iabase=iabase+2
      ibbase=ibbase-2
      jabase=jabase+ink
      jbbase=jbbase-ink
   30 continue
c
      if (iabase.ne.ibbase) go to 50
c     wavenumber n/4 (if it exists)
      ia=iabase
      ja=jabase
      scale=2*scale
cdir$ ivdep
      do 40 l=1,lot
      a(ja)=scale*work(ia)
      a(ja+inc)=-scale*work(ia+1)
      ia=ia+nx
      ja=ja+jump
   40 continue
c
   50 continue
      return
      end
      subroutine fft991(a,work,trigs,ifax,inc,jump,n,lot,isign)
      implicit real*8 (a-h), real*8 (o-z)
      dimension a(n),work(n),trigs(n),ifax(1)
c
c     subroutine "fft991" - multiple real/half-complex periodic
c     fast fourier transform
c
c     same as fft99 except that ordering of data corresponds to
c     that in mrfft2
c
c     procedure used to convert to half-length complex transform
c     is given by cooley, lewis and welch (j. sound vib., vol. 12
c     (1970), 315-337)
c
c     a is the array containing input and output data
c     work is an area of size (n+1)*lot
c     trigs is a previously prepared list of trig function values
c     ifax is a previously prepared list of factors of n/2
c     inc is the increment within each data 'vector'
c         (e.g. inc=1 for consecutively stored data)
c     jump is the increment between the start of each data vector
c     n is the length of the data vectors
c     lot is the number of data vectors
c     isign = +1 for transform from spectral to gridpoint
c           = -1 for transform from gridpoint to spectral
c
c     ordering of coefficients:
c         a(0),b(0),a(1),b(1),a(2),b(2),...,a(n/2),b(n/2)
c         where b(0)=b(n/2)=0; (n+2) locations required
c
c     ordering of data:
c         x(0),x(1),x(2),...,x(n-1)
c
c     vectorization is achieved on cray by doing the transforms in
c     parallel
c
c     *** n.b. n is assumed to be an even number
c
c     definition of transforms:
c     -------------------------
c
c     isign=+1: x(j)=sum(k=0,...,n-1)(c(k)*exp(2*i*j*k*pi/n))
c         where c(k)=a(k)+i*b(k) and c(n-k)=a(k)-i*b(k)
c
c     isign=-1: a(k)=(1/n)*sum(j=0,...,n-1)(x(j)*cos(2*j*k*pi/n))
c               b(k)=-(1/n)*sum(j=0,...,n-1)(x(j)*sin(2*j*k*pi/n))
c
c
c
      nfax=ifax(1)
      nx=n+1
      nh=n/2
      ink=inc+inc
      if (isign.eq.+1) go to 30
c
c     if necessary, transfer data to work area
      igo=50
      if (mod(nfax,2).eq.1) goto 40
      ibase=1
      jbase=1
      do 20 l=1,lot
      i=ibase
      j=jbase
cdir$ ivdep
      do 10 m=1,n
      work(j)=a(i)
      i=i+inc
      j=j+1
   10 continue
      ibase=ibase+jump
      jbase=jbase+nx
   20 continue
c
      igo=60
      go to 40
c
c     preprocessing (isign=+1)
c     ------------------------
c
   30 continue
      call fft99a(a,work,trigs,inc,jump,n,lot)
      igo=60
c
c     complex transform
c     -----------------
c
   40 continue
      ia=1
      la=1
      do 80 k=1,nfax
      if (igo.eq.60) go to 60
   50 continue
      call vpassm(a(ia),a(ia+inc),work(1),work(2),trigs,
     *   ink,2,jump,nx,lot,nh,ifax(k+1),la)
      igo=60
      go to 70
   60 continue
      call vpassm(work(1),work(2),a(ia),a(ia+inc),trigs,
     *    2,ink,nx,jump,lot,nh,ifax(k+1),la)
      igo=50
   70 continue
      la=la*ifax(k+1)
   80 continue
c

      if (isign.eq.-1) go to 130
c
c     if necessary, transfer data from work area
      if (mod(nfax,2).eq.1) go to 110
      ibase=1
      jbase=1
      do 100 l=1,lot
      i=ibase
      j=jbase
cdir$ ivdep
      do 90 m=1,n
      a(j)=work(i)
      i=i+1
      j=j+inc
   90 continue
      ibase=ibase+nx
      jbase=jbase+jump
  100 continue
c
c     fill in zeros at end
  110 continue
      ib=n*inc+1
cdir$ ivdep
      do 120 l=1,lot
      a(ib)=0
      a(ib+inc)=0
      ib=ib+jump
  120 continue
      go to 140
c
c     postprocessing (isign=-1):
c     --------------------------
c
  130 continue
      call fft99b(work,a,trigs,inc,jump,n,lot)
c
  140 continue
      return
      end
      subroutine set99 (trigs, ifax, n)
      implicit real*8 (a-h), real*8 (o-z)
      dimension ifax(13),trigs(1)
c
c mode 3 is used for real/half-complex transforms.  it is possible
c to do complex/complex transforms with other values of mode, but
c documentation of the details were not available when this routine
c was written.
c
      data mode /3/
      call fax (ifax, n, mode)
      i = ifax(1)
      if (ifax(i+1) .gt. 5 .or. n .le. 4) ifax(1) = -99
      if (ifax(1) .le. 0 ) then 
        write(6,*) ' set99 -- invalid n'
        stop'set99'
      endif
      call fftrig (trigs, n, mode)
      return
      end
      subroutine fax(ifax,n,mode)
      implicit real*8 (a-h), real*8 (o-z)
      dimension ifax(10)
      nn=n
      if (iabs(mode).eq.1) go to 10
      if (iabs(mode).eq.8) go to 10
      nn=n/2
      if ((nn+nn).eq.n) go to 10
      ifax(1)=-99
      return
   10 k=1
c     test for factors of 4
   20 if (mod(nn,4).ne.0) go to 30
      k=k+1
      ifax(k)=4
      nn=nn/4
      if (nn.eq.1) go to 80
      go to 20
c     test for extra factor of 2
   30 if (mod(nn,2).ne.0) go to 40
      k=k+1
      ifax(k)=2
      nn=nn/2
      if (nn.eq.1) go to 80
c     test for factors of 3
   40 if (mod(nn,3).ne.0) go to 50
      k=k+1
      ifax(k)=3
      nn=nn/3
      if (nn.eq.1) go to 80
      go to 40
c     now find remaining factors
   50 l=5
      inc=2
c     inc alternately takes on values 2 and 4
   60 if (mod(nn,l).ne.0) go to 70
      k=k+1
      ifax(k)=l
      nn=nn/l
      if (nn.eq.1) go to 80
      go to 60
   70 l=l+inc
      inc=6-inc
      go to 60
   80 ifax(1)=k-1
c     ifax(1) contains number of factors
      nfax=ifax(1)
c     sort factors into ascending order
      if (nfax.eq.1) go to 110
      do 100 ii=2,nfax
      istop=nfax+2-ii
      do 90 i=2,istop
      if (ifax(i+1).ge.ifax(i)) go to 90
      item=ifax(i)
      ifax(i)=ifax(i+1)
      ifax(i+1)=item
   90 continue
  100 continue
  110 continue
      return
      end
      subroutine fftrig(trigs,n,mode)
      implicit real*8 (a-h), real*8 (o-z)
      dimension trigs(1)
      pi=2*asin(1.0d0)
      imode=iabs(mode)
      nn=n
      if (imode.gt.1.and.imode.lt.6) nn=n/2
      del=(pi+pi)/(nn)
      l=nn+nn
      do 10 i=1,l,2
      angle=0.5*(i-1)*del
      trigs(i)=cos(angle)
      trigs(i+1)=sin(angle)
   10 continue
      if (imode.eq.1) return
      if (imode.eq.8) return
      del=0.5*del
      nh=(nn+1)/2
      l=nh+nh
      la=nn+nn
      do 20 i=1,l,2
      angle=0.5*(i-1)*del
      trigs(la+i)=cos(angle)
      trigs(la+i+1)=sin(angle)
   20 continue
      if (imode.le.3) return
      del=0.5*del
      la=la+nn
      if (mode.eq.5) go to 40
      do 30 i=2,nn
      angle=(i-1)*del
      trigs(la+i)=2*sin(angle)
   30 continue
      return
   40 continue
      del=0.5*del
      do 50 i=2,n
      angle=(i-1)*del
      trigs(la+i)=sin(angle)
   50 continue
      return
      end
      subroutine vpassm(a,b,c,d,trigs,inc1,inc2,inc3,inc4,lot,n,ifac,la)
      implicit real*8 (a-h), real*8 (o-z)
      dimension a(n),b(n),c(n),d(n),trigs(n)
c
c     subroutine "vpassm" - multiple version of "vpassa"
c     performs one pass through data
c     as part of multiple complex fft routine
c     a is first real input vector
c     b is first imaginary input vector
c     c is first real output vector
c     d is first imaginary output vector
c     trigs is precalculated table of sines " cosines
c     inc1 is addressing increment for a and b
c     inc2 is addressing increment for c and d
c     inc3 is addressing increment between a"s & b"s
c     inc4 is addressing increment between c"s & d"s
c     lot is the number of vectors
c     n is length of vectors
c     ifac is current factor of n
c     la is product of previous factors
c
      data sin36/0.587785252292473d0/,cos36/0.809016994374947d0/,
     *     sin72/0.951056516295154d0/,cos72/0.309016994374947d0/,
     *     sin60/0.866025403784437d0/
c
      m=n/ifac
      iink=m*inc1
      jink=la*inc2
      jump=(ifac-1)*jink
      ibase=0
      jbase=0
      igo=ifac-1
      if (igo.gt.4) return
      go to (10,50,90,130),igo
c
c     coding for factor 2
c
   10 ia=1
      ja=1
      ib=ia+iink
      jb=ja+jink
      do 20 l=1,la
      i=ibase
      j=jbase
cdir$ ivdep
      do 15 ijk=1,lot
      c(ja+j)=a(ia+i)+a(ib+i)
      d(ja+j)=b(ia+i)+b(ib+i)
      c(jb+j)=a(ia+i)-a(ib+i)
      d(jb+j)=b(ia+i)-b(ib+i)
      i=i+inc3
      j=j+inc4
   15 continue
      ibase=ibase+inc1
      jbase=jbase+inc2
   20 continue
      if (la.eq.m) return
      la1=la+1
      jbase=jbase+jump
      do 40 k=la1,m,la
      kb=k+k-2
      c1=trigs(kb+1)
      s1=trigs(kb+2)
      do 30 l=1,la
      i=ibase
      j=jbase
cdir$ ivdep
      do 25 ijk=1,lot
      c(ja+j)=a(ia+i)+a(ib+i)
      d(ja+j)=b(ia+i)+b(ib+i)
      c(jb+j)=c1*(a(ia+i)-a(ib+i))-s1*(b(ia+i)-b(ib+i))
      d(jb+j)=s1*(a(ia+i)-a(ib+i))+c1*(b(ia+i)-b(ib+i))
      i=i+inc3
      j=j+inc4
   25 continue
      ibase=ibase+inc1
      jbase=jbase+inc2
   30 continue
      jbase=jbase+jump
   40 continue
      return
c
c     coding for factor 3
c
   50 ia=1
      ja=1
      ib=ia+iink
      jb=ja+jink
      ic=ib+iink
      jc=jb+jink
      do 60 l=1,la
      i=ibase
      j=jbase
cdir$ ivdep
      do 55 ijk=1,lot
      c(ja+j)=a(ia+i)+(a(ib+i)+a(ic+i))
      d(ja+j)=b(ia+i)+(b(ib+i)+b(ic+i))
      c(jb+j)=(a(ia+i)-0.5*(a(ib+i)+a(ic+i)))
     $     -(sin60*(b(ib+i)-b(ic+i)))
      c(jc+j)=(a(ia+i)-0.5*(a(ib+i)+a(ic+i)))
     $     +(sin60*(b(ib+i)-b(ic+i)))
      d(jb+j)=(b(ia+i)-0.5*(b(ib+i)+b(ic+i)))
     $     +(sin60*(a(ib+i)-a(ic+i)))
      d(jc+j)=(b(ia+i)-0.5*(b(ib+i)+b(ic+i)))
     $     -(sin60*(a(ib+i)-a(ic+i)))
      i=i+inc3
      j=j+inc4
   55 continue
      ibase=ibase+inc1
      jbase=jbase+inc2
   60 continue
      if (la.eq.m) return
      la1=la+1
      jbase=jbase+jump
      do 80 k=la1,m,la
      kb=k+k-2
      kc=kb+kb
      c1=trigs(kb+1)
      s1=trigs(kb+2)
      c2=trigs(kc+1)
      s2=trigs(kc+2)
      do 70 l=1,la
      i=ibase
      j=jbase
cdir$ ivdep
      do 65 ijk=1,lot
      c(ja+j)=a(ia+i)+(a(ib+i)+a(ic+i))
      d(ja+j)=b(ia+i)+(b(ib+i)+b(ic+i))
      c(jb+j)=
     *    c1*((a(ia+i)-0.5*(a(ib+i)+a(ic+i)))
     $     -(sin60*(b(ib+i)-b(ic+i))))
     *   -s1*((b(ia+i)-0.5*(b(ib+i)+b(ic+i)))
     $     +(sin60*(a(ib+i)-a(ic+i))))
      d(jb+j)=
     *    s1*((a(ia+i)-0.5*(a(ib+i)+a(ic+i)))
     $     -(sin60*(b(ib+i)-b(ic+i))))
     *   +c1*((b(ia+i)-0.5*(b(ib+i)+b(ic+i)))
     $     +(sin60*(a(ib+i)-a(ic+i))))
      c(jc+j)=
     *    c2*((a(ia+i)-0.5*(a(ib+i)+a(ic+i)))
     $     +(sin60*(b(ib+i)-b(ic+i))))
     *   -s2*((b(ia+i)-0.5*(b(ib+i)+b(ic+i)))
     $     -(sin60*(a(ib+i)-a(ic+i))))
      d(jc+j)=
     *    s2*((a(ia+i)-0.5*(a(ib+i)+a(ic+i)))
     $     +(sin60*(b(ib+i)-b(ic+i))))
     *   +c2*((b(ia+i)-0.5*(b(ib+i)+b(ic+i)))
     $     -(sin60*(a(ib+i)-a(ic+i))))
      i=i+inc3
      j=j+inc4
   65 continue
      ibase=ibase+inc1
      jbase=jbase+inc2
   70 continue
      jbase=jbase+jump
   80 continue
      return
c
c     coding for factor 4
c
   90 ia=1
      ja=1
      ib=ia+iink
      jb=ja+jink
      ic=ib+iink
      jc=jb+jink
      id=ic+iink
      jd=jc+jink
      do 100 l=1,la
      i=ibase
      j=jbase
cdir$ ivdep
      do 95 ijk=1,lot
      c(ja+j)=(a(ia+i)+a(ic+i))+(a(ib+i)+a(id+i))
      c(jc+j)=(a(ia+i)+a(ic+i))-(a(ib+i)+a(id+i))
      d(ja+j)=(b(ia+i)+b(ic+i))+(b(ib+i)+b(id+i))
      d(jc+j)=(b(ia+i)+b(ic+i))-(b(ib+i)+b(id+i))
      c(jb+j)=(a(ia+i)-a(ic+i))-(b(ib+i)-b(id+i))
      c(jd+j)=(a(ia+i)-a(ic+i))+(b(ib+i)-b(id+i))
      d(jb+j)=(b(ia+i)-b(ic+i))+(a(ib+i)-a(id+i))
      d(jd+j)=(b(ia+i)-b(ic+i))-(a(ib+i)-a(id+i))
      i=i+inc3
      j=j+inc4
   95 continue
      ibase=ibase+inc1
      jbase=jbase+inc2
  100 continue
      if (la.eq.m) return
      la1=la+1
      jbase=jbase+jump
      do 120 k=la1,m,la
      kb=k+k-2
      kc=kb+kb
      kd=kc+kb
      c1=trigs(kb+1)
      s1=trigs(kb+2)
      c2=trigs(kc+1)
      s2=trigs(kc+2)
      c3=trigs(kd+1)
      s3=trigs(kd+2)
      do 110 l=1,la
      i=ibase
      j=jbase
cdir$ ivdep
      do 105 ijk=1,lot
      c(ja+j)=(a(ia+i)+a(ic+i))+(a(ib+i)+a(id+i))
      d(ja+j)=(b(ia+i)+b(ic+i))+(b(ib+i)+b(id+i))
      c(jc+j)=
     *    c2*((a(ia+i)+a(ic+i))-(a(ib+i)+a(id+i)))
     *   -s2*((b(ia+i)+b(ic+i))-(b(ib+i)+b(id+i)))
      d(jc+j)=
     *    s2*((a(ia+i)+a(ic+i))-(a(ib+i)+a(id+i)))
     *   +c2*((b(ia+i)+b(ic+i))-(b(ib+i)+b(id+i)))
      c(jb+j)=
     *    c1*((a(ia+i)-a(ic+i))-(b(ib+i)-b(id+i)))
     *   -s1*((b(ia+i)-b(ic+i))+(a(ib+i)-a(id+i)))
      d(jb+j)=
     *    s1*((a(ia+i)-a(ic+i))-(b(ib+i)-b(id+i)))
     *   +c1*((b(ia+i)-b(ic+i))+(a(ib+i)-a(id+i)))
      c(jd+j)=
     *    c3*((a(ia+i)-a(ic+i))+(b(ib+i)-b(id+i)))
     *   -s3*((b(ia+i)-b(ic+i))-(a(ib+i)-a(id+i)))
      d(jd+j)=
     *    s3*((a(ia+i)-a(ic+i))+(b(ib+i)-b(id+i)))
     *   +c3*((b(ia+i)-b(ic+i))-(a(ib+i)-a(id+i)))
      i=i+inc3
      j=j+inc4
  105 continue
      ibase=ibase+inc1
      jbase=jbase+inc2
  110 continue
      jbase=jbase+jump
  120 continue
      return
c
c     coding for factor 5
c
  130 ia=1
      ja=1
      ib=ia+iink
      jb=ja+jink
      ic=ib+iink
      jc=jb+jink
      id=ic+iink
      jd=jc+jink
      ie=id+iink
      je=jd+jink
      do 140 l=1,la
      i=ibase
      j=jbase
cdir$ ivdep
      do 135 ijk=1,lot
      c(ja+j)=a(ia+i)+(a(ib+i)+a(ie+i))+(a(ic+i)+a(id+i))
      d(ja+j)=b(ia+i)+(b(ib+i)+b(ie+i))+(b(ic+i)+b(id+i))
      c(jb+j)=(a(ia+i)+cos72*(a(ib+i)+a(ie+i))-cos36*(a(ic+i)+a(id+i)))
     *  -(sin72*(b(ib+i)-b(ie+i))+sin36*(b(ic+i)-b(id+i)))
      c(je+j)=(a(ia+i)+cos72*(a(ib+i)+a(ie+i))-cos36*(a(ic+i)+a(id+i)))
     *  +(sin72*(b(ib+i)-b(ie+i))+sin36*(b(ic+i)-b(id+i)))
      d(jb+j)=(b(ia+i)+cos72*(b(ib+i)+b(ie+i))-cos36*(b(ic+i)+b(id+i)))
     *  +(sin72*(a(ib+i)-a(ie+i))+sin36*(a(ic+i)-a(id+i)))
      d(je+j)=(b(ia+i)+cos72*(b(ib+i)+b(ie+i))-cos36*(b(ic+i)+b(id+i)))
     *  -(sin72*(a(ib+i)-a(ie+i))+sin36*(a(ic+i)-a(id+i)))
      c(jc+j)=(a(ia+i)-cos36*(a(ib+i)+a(ie+i))+cos72*(a(ic+i)+a(id+i)))
     *  -(sin36*(b(ib+i)-b(ie+i))-sin72*(b(ic+i)-b(id+i)))
      c(jd+j)=(a(ia+i)-cos36*(a(ib+i)+a(ie+i))+cos72*(a(ic+i)+a(id+i)))
     *  +(sin36*(b(ib+i)-b(ie+i))-sin72*(b(ic+i)-b(id+i)))
      d(jc+j)=(b(ia+i)-cos36*(b(ib+i)+b(ie+i))+cos72*(b(ic+i)+b(id+i)))
     *  +(sin36*(a(ib+i)-a(ie+i))-sin72*(a(ic+i)-a(id+i)))
      d(jd+j)=(b(ia+i)-cos36*(b(ib+i)+b(ie+i))+cos72*(b(ic+i)+b(id+i)))
     *  -(sin36*(a(ib+i)-a(ie+i))-sin72*(a(ic+i)-a(id+i)))
      i=i+inc3
      j=j+inc4
  135 continue
      ibase=ibase+inc1
      jbase=jbase+inc2
  140 continue
      if (la.eq.m) return
      la1=la+1
      jbase=jbase+jump
      do 160 k=la1,m,la
      kb=k+k-2
      kc=kb+kb
      kd=kc+kb
      ke=kd+kb
      c1=trigs(kb+1)
      s1=trigs(kb+2)
      c2=trigs(kc+1)
      s2=trigs(kc+2)
      c3=trigs(kd+1)
      s3=trigs(kd+2)
      c4=trigs(ke+1)
      s4=trigs(ke+2)
      do 150 l=1,la
      i=ibase
      j=jbase
cdir$ ivdep
      do 145 ijk=1,lot
      c(ja+j)=a(ia+i)+(a(ib+i)+a(ie+i))+(a(ic+i)+a(id+i))
      d(ja+j)=b(ia+i)+(b(ib+i)+b(ie+i))+(b(ic+i)+b(id+i))
      c(jb+j)=
     *    c1*((a(ia+i)+cos72*(a(ib+i)+a(ie+i))-cos36*(a(ic+i)+a(id+i)))
     *      -(sin72*(b(ib+i)-b(ie+i))+sin36*(b(ic+i)-b(id+i))))
     *   -s1*((b(ia+i)+cos72*(b(ib+i)+b(ie+i))-cos36*(b(ic+i)+b(id+i)))
     *      +(sin72*(a(ib+i)-a(ie+i))+sin36*(a(ic+i)-a(id+i))))
      d(jb+j)=
     *    s1*((a(ia+i)+cos72*(a(ib+i)+a(ie+i))-cos36*(a(ic+i)+a(id+i)))
     *      -(sin72*(b(ib+i)-b(ie+i))+sin36*(b(ic+i)-b(id+i))))
     *   +c1*((b(ia+i)+cos72*(b(ib+i)+b(ie+i))-cos36*(b(ic+i)+b(id+i)))
     *      +(sin72*(a(ib+i)-a(ie+i))+sin36*(a(ic+i)-a(id+i))))
      c(je+j)=
     *    c4*((a(ia+i)+cos72*(a(ib+i)+a(ie+i))-cos36*(a(ic+i)+a(id+i)))
     *      +(sin72*(b(ib+i)-b(ie+i))+sin36*(b(ic+i)-b(id+i))))
     *   -s4*((b(ia+i)+cos72*(b(ib+i)+b(ie+i))-cos36*(b(ic+i)+b(id+i)))
     *      -(sin72*(a(ib+i)-a(ie+i))+sin36*(a(ic+i)-a(id+i))))
      d(je+j)=
     *    s4*((a(ia+i)+cos72*(a(ib+i)+a(ie+i))-cos36*(a(ic+i)+a(id+i)))
     *      +(sin72*(b(ib+i)-b(ie+i))+sin36*(b(ic+i)-b(id+i))))
     *   +c4*((b(ia+i)+cos72*(b(ib+i)+b(ie+i))-cos36*(b(ic+i)+b(id+i)))
     *      -(sin72*(a(ib+i)-a(ie+i))+sin36*(a(ic+i)-a(id+i))))
      c(jc+j)=
     *    c2*((a(ia+i)-cos36*(a(ib+i)+a(ie+i))+cos72*(a(ic+i)+a(id+i)))
     *      -(sin36*(b(ib+i)-b(ie+i))-sin72*(b(ic+i)-b(id+i))))
     *   -s2*((b(ia+i)-cos36*(b(ib+i)+b(ie+i))+cos72*(b(ic+i)+b(id+i)))
     *      +(sin36*(a(ib+i)-a(ie+i))-sin72*(a(ic+i)-a(id+i))))
      d(jc+j)=
     *    s2*((a(ia+i)-cos36*(a(ib+i)+a(ie+i))+cos72*(a(ic+i)+a(id+i)))
     *      -(sin36*(b(ib+i)-b(ie+i))-sin72*(b(ic+i)-b(id+i))))
     *   +c2*((b(ia+i)-cos36*(b(ib+i)+b(ie+i))+cos72*(b(ic+i)+b(id+i)))
     *      +(sin36*(a(ib+i)-a(ie+i))-sin72*(a(ic+i)-a(id+i))))
      c(jd+j)=
     *    c3*((a(ia+i)-cos36*(a(ib+i)+a(ie+i))+cos72*(a(ic+i)+a(id+i)))
     *      +(sin36*(b(ib+i)-b(ie+i))-sin72*(b(ic+i)-b(id+i))))
     *   -s3*((b(ia+i)-cos36*(b(ib+i)+b(ie+i))+cos72*(b(ic+i)+b(id+i)))
     *      -(sin36*(a(ib+i)-a(ie+i))-sin72*(a(ic+i)-a(id+i))))
      d(jd+j)=
     *    s3*((a(ia+i)-cos36*(a(ib+i)+a(ie+i))+cos72*(a(ic+i)+a(id+i)))
     *      +(sin36*(b(ib+i)-b(ie+i))-sin72*(b(ic+i)-b(id+i))))
     *   +c3*((b(ia+i)-cos36*(b(ib+i)+b(ie+i))+cos72*(b(ic+i)+b(id+i)))
     *      -(sin36*(a(ib+i)-a(ie+i))-sin72*(a(ic+i)-a(id+i))))
      i=i+inc3
      j=j+inc4
  145 continue
      ibase=ibase+inc1
      jbase=jbase+inc2
  150 continue
      jbase=jbase+jump
  160 continue
      return
      end



      subroutine fftfax (n,ifax,trigs)
      implicit none
      integer i
      integer ifax(13)
      integer mode
      integer n
      real*8 trigs(1)

c----
c 
c mode 3 is used for real/half-complex transforms.  it is possible
c to do complex/complex transforms with other values of mode, but
c documentation of the details were not available when this routine
c was written.
c 
      data mode/3/
      call fax (ifax,n,mode)
      i=ifax(1)
      if (ifax(i+1).gt.5.or.n.le.4) ifax(1)=-99
      if (ifax(1).le.0) print *,' fftfax -- invalid n'
      call fftrig (trigs,n,mode)
      return
      end

