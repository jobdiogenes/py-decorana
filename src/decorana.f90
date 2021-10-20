!     Decorana: the classic goodie by Mark O. Hill
!     Ported to R by Jari Oksanen, September 2001:
!     * Kept only the numerical engine: subroutines eigy, cutup and
!     yxmult, and the subroutines eigy calls.
!     * Changed to lowercase so that it looks prettier.
!     * Changed to double precision throughout, with
!     'implicit double precision (a-h,o-z)' for automatic variables.
!     * Changed 'ifix' to 'int', since g77 doesn't know 'ifix' for double.
!     * Removed all write statements (they don't show in GUI term --
!     else it would make sense to add a `trace' parameter.)

!     ORIGINAL INTRODUCTORY COMMENTS:

!     Cornell Ecology Program Decorana - written by M.O. Hill, July 1979
!     Source code and accompanying documentation are available from
!     Hugh G. Gauch, Jr., Ecology and Systematics, Cornell University,
!     Ithaca, New York 14850.

!     Performs detrended correspondence analysis;  also will do reciprocal
!     averaging as a special case.

!     **** further modified by Dr Peter R. Minchin May-June 1997
!     - changed tolerance to 0.000005 and iteration limit to 999 in
!     eigy - strict settings of oksanen & minchin 1997
!     - corrected the order-dependent bug in smooth
!     - added parameter statements for maximum dimensions
!     - changed max. dimensions to 5000, 5000, 330000
!     - updated handling of file names retrieved from command line
!     - moved character data to character variables
!     - now accepts relaxed ccf format, with maximum number of pairs
!     per record on line 3
!     see: Oksanen, J. & Minchin, P.R. 1997. Instability of ordination
!     results under changes in input data order: explanations and
!     remedies. Journal of Vegetation Science 8: 447-454.


    subroutine cutup(x,ix,mi,mk)
!     takes a vector x and cuts up into (mk-4) segments, putting a
!     segmented version of the vector in ix.
    implicit double precision (a-h,o-z)
    double precision :: x(mi)
    integer :: ix(mi)
    mmk=mk-4
    maxk=mk-2
    call xmaxmi(x,axmax,axmin,mi)
    axbit=(axmax-axmin)/float(mmk)
    do 10 i=1,mi
    !---  iax=ifix((x(i)-axmin)/axbit)+3
        iax=int((x(i)-axmin)/axbit)+3
        if(iax < 3) iax=3
        if(iax > maxk) iax=maxk
        ix(i)=iax
    10 END DO
    return
    end subroutine cutup

    subroutine trans(y,yy, &
    x,neig,ira,aidot,xeig1,xeig2,xeig3,ix1,ix2,ix3, &
    mi,mk,n,nid,ibegin,iend,idat,qidat)
!     this subroutine is the crux of the whole program, in that it
!     takes a set of species scores y and iterates to find a new set
!     of scores yy.  repeated iteration of this subroutine would lead
!     eventually to the correct solution (except that the scores need
!     to be divided by the y-totals adotj at each iteration).  the
!     calling program eigy is made lengthy by some fancy algebra put
!     there to speed up the calculation.  essentially trans is the
!     standard reciprocal averaging iteration with either detrending
!     with respect to previously derived axes (in the case of detrended
!     correspondence analysis) or orthogonalization with respect to
!     them (in the case of reciprocal averaging).
    implicit double precision (a-h,o-z)
    double precision :: x(mi),xeig1(mi),xeig2(mi),xeig3(mi)
    double precision :: y(n),yy(n),aidot(mi),qidat(nid)
    integer :: ix1(mi),ix2(mi),ix3(mi),idat(nid),ibegin(mi),iend(mi)
    call yxmult(y,x,mi,n,nid,ibegin,iend,idat,qidat)
    do 10 i=1,mi
        x(i)=x(i)/aidot(i)
    10 END DO
    if(neig == 0) goto 200
    if(ira == 1) goto 100
    call detrnd(x,aidot,ix1,mi,mk)
    if(neig == 1) goto 200
    call detrnd(x,aidot,ix2,mi,mk)
    if(neig == 2) goto 90
    call detrnd(x,aidot,ix3,mi,mk)
    call detrnd(x,aidot,ix2,mi,mk)
    90 call detrnd(x,aidot,ix1,mi,mk)
    goto 200
    100 a1=0.0
    do 110 i=1,mi
        a1=a1+aidot(i)*x(i)*xeig1(i)
    110 END DO
    do 120 i=1,mi
        x(i)=x(i)-a1*xeig1(i)
    120 END DO
    if(neig == 1) goto 200
    a2=0.0
    do 130 i=1,mi
        a2=a2+aidot(i)*x(i)*xeig2(i)
    130 END DO
    do 140 i=1,mi
        x(i)=x(i)-a2*xeig2(i)
    140 END DO
    if(neig == 2) goto 200
    a3=0.0
    do 150 i=1,mi
        a3=a3+aidot(i)*x(i)*xeig3(i)
    150 END DO
    do 160 i=1,mi
        x(i)=x(i)-a3*xeig3(i)
    160 END DO
    200 call xymult(x,yy,mi,n,nid,ibegin,iend,idat,qidat)
    return
    end subroutine trans

    subroutine detrnd(x,aidot,ix,mi,mk)
    implicit double precision (a-h,o-z)
    double precision :: x(mi),z(50),zn(50),zbar(50),aidot(mi)
    integer :: ix(mi)
!     starts with a vector x and detrends with respect to groups defined
!     by ix.  detrending is in blocks of 3 units at a time, and the
!     result calculated is the average of the 3 possible results that
!     can be obtained, corresponding to 3 possible starting positions
!     for the blocks of 3.
    do 10 k=1,mk
        z(k)=0.0
        zn(k)=0.0
    10 END DO
    do 20 i=1,mi
        k=ix(i)
        z(k)=z(k)+x(i)*aidot(i)
        zn(k)=zn(k)+aidot(i)
    20 END DO
    mmk=mk-1
    do 30 k=2,mmk
        zbar(k)=(z(k-1)+z(k)+z(k+1))/(zn(k-1)+zn(k)+zn(k+1)+1.0e-12)
    30 END DO
    mmk=mmk-1
    do 35 k=3,mmk
        z(k)=(zbar(k-1)+zbar(k)+zbar(k+1))/3.0
    35 END DO
    do 40 i=1,mi
        k=ix(i)
        x(i)=x(i)-z(k)
    40 END DO
    return
    end subroutine detrnd

! Function segfit is identical to detrnd, but it also returns the fitted
! values z. x is respone in input, and residuals in output. --Added by
! J. Oksanen 1 Oct, 2010.
    subroutine segfit(x,aidot,ix,mi,mk,fit)
    implicit double precision (a-h,o-z)
    double precision :: x(mi),z(50),zn(50),zbar(50),aidot(mi)
    double precision :: fit(mi)
    integer :: ix(mi)
    do 10 k=1,mk
        z(k)=0.0
        zn(k)=0.0
    10 END DO
    do 20 i=1,mi
        k=ix(i)
        z(k)=z(k)+x(i)*aidot(i)
        zn(k)=zn(k)+aidot(i)
    20 END DO
    mmk=mk-1
    do 30 k=2,mmk
        zbar(k)=(z(k-1)+z(k)+z(k+1))/(zn(k-1)+zn(k)+zn(k+1)+1.0e-12)
    30 END DO
    mmk=mmk-1
    do 35 k=3,mmk
        z(k)=(zbar(k-1)+zbar(k)+zbar(k+1))/3.0
    35 END DO
    do 40 i=1,mi
        k=ix(i)
        fit(i) = z(k)
        x(i)=x(i)-z(k)
    40 END DO
    return
    end subroutine segfit

    subroutine eigy(x,y,eig,neig,ira,iresc,short, &
    mi,mk,n,nid,ibegin,iend,idat,qidat,y2,y3,y4,y5, &
    xeig1,xeig2,xeig3,ix1,ix2,ix3,aidot,adotj)
!     extracts an eigenvector y with eigenvalue eig.  the algebra
!     is a little complicated, but consists essentially of repre-
!     senting the transformation (subroutine trans) approximately
!     by a tridiagonal 4x4 matrix.  the eigenproblem for the
!     tridiagonal matrix is solved and this solution is plugged
!     back in to obtain a new trial vector.
!     after getting the eigenvector, the scores may be rescaled
!     (subroutine strtch).
    implicit double precision (a-h,o-z)
    double precision :: x(mi),y(n),y2(n),y3(n),y4(n),y5(n)
    double precision :: xeig1(mi),xeig2(mi),xeig3(mi),aidot(mi),adotj(n)
    double precision :: qidat(nid)
    integer :: ibegin(mi),iend(mi),idat(nid),ix1(mi),ix2(mi),ix3(mi)
!     strings to print R warnings
    character (len=64) warning
    character (len=2) axnam
    tot=0.0
    do 10 j=1,n
        tot=tot+adotj(j)
        y(j)=float(j)
    10 END DO
    y(1)=1.1
!---  tolerance reduced by p.minchin jan 1997
!     tol=0.0001
    tol=0.000005
    call trans(y,y, &
    x,neig,ira,aidot,xeig1,xeig2,xeig3,ix1,ix2,ix3, &
    mi,mk,n,nid,ibegin,iend,idat,qidat)
    icount=0
    20 a=0.0
    do 30 j=1,n
        a=a+y(j)*adotj(j)
    30 END DO
    a=a/tot
    ex=0.0
    do 40 j=1,n
        ay=y(j)-a
        ex=ex+ay*ay*adotj(j)
        y(j)=ay
    40 END DO
    ex=sqrt(ex)
    do 50 j=1,n
        y(j)=y(j)/ex
    50 END DO
    call trans(y,y2, &
    x,neig,ira,aidot,xeig1,xeig2,xeig3,ix1,ix2,ix3, &
    mi,mk,n,nid,ibegin,iend,idat,qidat)
    a=0.0
    a11=0.0
    a12=0.0
    a22=0.0
    a23=0.0
    a33=0.0
    a34=0.0
    a44=0.0
    do 60 j=1,n
        ay=y2(j)
        y2(j)=ay/adotj(j)
        a=a+ay
        a11=a11+ay*y(j)
    60 END DO
    a=a/tot
    do 70 j=1,n
        ay=y2(j)-(a+a11*y(j))
        a12=a12+ay*ay*adotj(j)
        y2(j)=ay
    70 END DO
    a12=sqrt(a12)
    do 80 j=1,n
        y2(j)=y2(j)/a12
    80 END DO
!---------removed write statements--------------------------
!     if(icount.eq.0) write(*,1000)
!     1000 format(/1x)
!     write(*,1011) a12,icount
!     1011 format(1x,'residual',f10.6,'       at iteration',i3)
    if(a12 < tol) goto 200
!---  maximum iteration limit increased by P.Minchin jan 1997
!     if(icount.gt.9) goto 200
    if(icount > 999) goto 200
    icount=icount+1
    call trans(y2,y3, &
    x,neig,ira,aidot,xeig1,xeig2,xeig3,ix1,ix2,ix3, &
    mi,mk,n,nid,ibegin,iend,idat,qidat)
    a=0.0
    b13=0.0
    do 90 j=1,n
        ay=y3(j)
        y3(j)=ay/adotj(j)
        a=a+ay
        a22=a22+ay*y2(j)
        b13=b13+ay*y(j)
    90 END DO
    a=a/tot
    do 100 j=1,n
        ay=y3(j)-(a+a22*y2(j)+b13*y(j))
        a23=a23+ay*ay*adotj(j)
        y3(j)=ay
    100 END DO
    a23=sqrt(a23)
    if(a23 > tol) goto 105
    a23=0.0
    goto 160
    105 continue
    do 110 j=1,n
        y3(j)=y3(j)/a23
    110 END DO
    call trans(y3,y4, &
    x,neig,ira,aidot,xeig1,xeig2,xeig3,ix1,ix2,ix3, &
    mi,mk,n,nid,ibegin,iend,idat,qidat)
    a=0.0
    b14=0.0
    b24=0.0
    do 120 j=1,n
        ay=y4(j)
        y4(j)=y4(j)/adotj(j)
        a=a+ay
        a33=a33+ay*y3(j)
        b14=b14+ay*y(j)
        b24=b24+ay*y2(j)
    120 END DO
    a=a/tot
    do 130 j=1,n
        ay=y4(j)-(a+a33*y3(j)+b14*y(j)+b24*y2(j))
        a34=a34+ay*ay*adotj(j)
        y4(j)=ay
    130 END DO
    a34=sqrt(a34)
    if(a34 > tol) goto 135
    a34=0.0
    goto 160
    135 continue
    do 140 j=1,n
        y4(j)=y4(j)/a34
    140 END DO
    call trans(y4,y5, &
    x,neig,ira,aidot,xeig1,xeig2,xeig3,ix1,ix2,ix3, &
    mi,mk,n,nid,ibegin,iend,idat,qidat)
    do 150 j=1,n
        a44=a44+y4(j)*y5(j)
    150 END DO
!     we now have the tridiagonal representation of trans.  solve
!     eigenproblem for tridiagonal matrix.
    160 ax1=1.0
    ax2=0.1
    ax3=0.01
    ax4=0.001
    do 170 itimes=1,100
        axx1=a11*ax1+a12*ax2
        axx2=a12*ax1+a22*ax2+a23*ax3
        axx3=a23*ax2+a33*ax3+a34*ax4
        axx4=a34*ax3+a44*ax4
        ax1=a11*axx1+a12*axx2
        ax2=a12*axx1+a22*axx2+a23*axx3
        ax3=a23*axx2+a33*axx3+a34*axx4
        ax4=a34*axx3+a44*axx4
        ex=sqrt(ax1**2+ax2**2+ax3**2+ax4**2)
        ax1=ax1/ex
        ax2=ax2/ex
        ax3=ax3/ex
        ax4=ax4/ex
        if(itimes /= (itimes/5)*5) goto 170
        exx=sqrt(ex)
        resi=sqrt((ax1-axx1/exx)**2+(ax2-axx2/exx)**2+ &
        (ax3-axx3/exx)**2+(ax4-axx4/exx)**2)
        if(resi < tol*0.05) goto 180
    170 END DO
    180 continue
    do 190 j=1,n
        y(j)=ax1*y(j)+ax2*y2(j)+ax3*y3(j)+ax4*y4(j)
    190 END DO
    goto 20
!-----------Removed write statements, added 200 continue--------
    200 continue
!     200 write(*,1010) a11
!     1010 format(1x,'eigenvalue',f10.5)
!     if(a12.gt.tol) write(*,1012) tol
!     1012 format(1x,'*** beware ***     residual bigger than tolerance',
!     1', which is',f10.6)
!     R version of the above warning. You must change the length of
!     character*n warning definition if you change the warning text
    if (a12 > tol) then
        if (neig == 0) axnam = "1"
        if (neig == 1) axnam = "2"
        if (neig == 2) axnam = "3"
        if (neig == 3) axnam = "4"
        warning = "residual bigger than tolerance on axis "//axnam
        call rwarn(warning)
    end if
!     we calculate x from y, and set x to unit length if reciprocal
!     averaging option is in force (ira=1)
    call xmaxmi(y,aymax,aymin,n)
    sign=1.0
    if(-aymin > aymax) sign=-1.0
    do 210 j=1,n
        y(j)=y(j)*sign
    210 END DO
    call yxmult(y,x,mi,n,nid,ibegin,iend,idat,qidat)
    do 220 i=1,mi
        x(i)=x(i)/aidot(i)
    220 END DO
    if(iresc == 0) goto 225
    if(a11 > 0.999) goto 225
    do 223 i=1,iresc
        monit=0
        if(i == 1 .OR. i == iresc) monit=1
        call strtch(x,y,short,monit, &
        mi,n,nid,aidot,ibegin,iend,idat,qidat)
    223 END DO
    eig=a11
    return
    225 axlong=0.0
    do 230 i=1,mi
        axlong=axlong+aidot(i)*x(i)**2
    230 END DO
    axlong=sqrt(axlong)
    do 240 i=1,mi
        x(i)=x(i)/axlong
    240 END DO
    do 250 j=1,n
        y(j)=y(j)/axlong
    250 END DO
!     it remains to scale y to unit within-sample standard deviation
    sumsq=0.0
    do 260 i=1,mi
        id1=ibegin(i)
        id2=iend(i)
        ax=x(i)
        do 255 id=id1,id2
            j=idat(id)
            sumsq=sumsq+qidat(id)*(ax-y(j))**2
        255 END DO
    260 END DO
    sd=sqrt(sumsq/tot)
    if(a11 <= 0.999) goto 265
    sd=aymax/axlong
    sd1=-aymin/axlong
    if(sd1 > sd) sd=sd1
    265 continue
    do 270 j=1,n
        y(j)=y(j)/sd
    270 END DO
    eig=a11
    return
    end subroutine eigy

    subroutine xymult(x,y,mi,n,nid,ibegin,iend,idat,qidat)
!     starts with vector x and forms matrix product y=ax
    implicit double precision (a-h,o-z)
    double precision :: x(mi),y(n),qidat(nid)
    integer :: ibegin(mi),iend(mi),idat(nid)
    do 10 j=1,n
        y(j)=0.0
    10 END DO
    do 30 i=1,mi
        id1=ibegin(i)
        id2=iend(i)
        ax=x(i)
        do 20 id=id1,id2
            j=idat(id)
            y(j)=y(j)+ax*qidat(id)
        20 END DO
    30 END DO
    return
    end subroutine xymult

    subroutine yxmult(y,x,mi,n,nid,ibegin,iend,idat,qidat)
!     starts with vector y and forms matrix product x=ay
    implicit double precision (a-h,o-z)
    double precision :: x(mi),y(n),qidat(nid)
    integer :: ibegin(mi),iend(mi),idat(nid)
    do 20 i=1,mi
        id1=ibegin(i)
        id2=iend(i)
        ax=0.0
        do 10 id=id1,id2
            j=idat(id)
            ax=ax+y(j)*qidat(id)
        10 END DO
        x(i)=ax
    20 END DO
    return
    end subroutine yxmult

    subroutine xmaxmi(x,axmax,axmin,m)
!     forms maximum and minimum of x(m)
    implicit double precision (a-h,o-z)
    double precision :: x(m)
    axmax=-1.0e10
    axmin=-axmax
    do 10 i=1,m
        ax=x(i)
        if(ax > axmax) axmax=ax
        if(ax < axmin) axmin=ax
    10 END DO
    return
    end subroutine xmaxmi


    subroutine strtch(x,y,short,monit, &
    mi,n,nid,aidot,ibegin,iend,idat,qidat)
!     takes an axis (x,y) and scales to unit mean square dev of species
!     scores per sample.  an attempt is made for longer axes (l > short)
!     to squeeze them in and out so that they have the right mean square
!     deviation all the way along the axis and not only on average.
    implicit double precision (a-h,o-z)
    double precision :: x(mi),y(n),aidot(mi),qidat(nid)
    double precision :: zn(50),zv(50)
    integer :: ibegin(mi),iend(mi),idat(nid)
!     Use monit so that gfortran -Wall does not warn for it being unused
    monit = 0
!---  common block added by p.minchin feb 1988
!     common /lunits/ iuinp1,iuinp2,iuout1,iuout2,iuout3
    do 200 icount=1,2
        mk=20
        call segmnt(x,y,zn,zv,mi,mk,n,nid, &
        aidot,ibegin,iend,idat,qidat)
        call smooth(zv,mk)
        call smooth(zn,mk)
        zvsum=0.0
        do 50 k=1,mk
            zvsum=zvsum+zv(k)/zn(k)
        50 END DO
        sd=sqrt(zvsum/float(mk))
    !     we want mean within-sample square deviation to be 1.0, so we divide
    !     everything by sd
        along=0.0
        do 60 i=1,mi
            ax=x(i)/sd
            x(i)=ax
            if(along < ax) along=ax
        60 END DO
    !--------Removed write statements----------------------------
    !     if(icount.eq.1.and.monit.eq.1) write(*,1000)
    !     1000 format(/1x)
    !     if(monit.eq.1) write(*,1001) along
    !     1001 format(1x,'length of gradient',f10.3)
        do 70 j=1,n
            y(j)=y(j)/sd
        70 END DO
        if(along < short) return
        if(icount == 2) return
    !     mk=ifix(along*5.0)+1
        mk=int(along*5.0)+1
        if(mk < 10) mk=10
        if(mk > 45) mk=45
        call segmnt(x,y,zn,zv,mi,mk,n,nid,aidot,ibegin,iend,idat,qidat)
        call smooth(zv,mk)
        call smooth(zn,mk)
        zvsum=0.0
        do 100 k=1,mk
            azv=1.0/sqrt(0.2/along+zv(k)/zn(k))
            zvsum=zvsum+azv
            zv(k)=azv
        100 END DO
        do 110 k=1,mk
            zv(k)=zv(k)*along/zvsum
        110 END DO
    !----------Removed write statements-------------------------
    !     if(monit.eq.1) write(*,1002) (zv(k),k=1,mk)
    !     1002 format(1x,'length of segments',10f6.2)
        az=0.0
        zn(1)=0.0
        do 120 k=1,mk
            az=az+zv(k)
            zn(k+1)=az
        120 END DO
        axbit=along/float(mk)
        do 130 j=1,n
        !     iay=ifix(y(j)/axbit)+1
            iay=int(y(j)/axbit)+1
            if(iay < 1) iay=1
            if(iay > mk) iay=mk
            y(j)=zn(iay)+zv(iay)*(y(j)/axbit-float(iay-1))
        130 END DO
        call yxmult(y,x,mi,n,nid,ibegin,iend,idat,qidat)
        do 140 i=1,mi
            x(i)=x(i)/aidot(i)
        140 END DO
    200 END DO
    return
    end subroutine strtch

    subroutine smooth(z,mk)
    implicit double precision (a-h,o-z)
    double precision :: z(mk)
!     takes a vector z and does (1,2,1)-smoothing until no blanks left
!     and then 2 more iterations of (1,2,1)-smoothing.  if no blanks to
!     begin with, then does 3 smoothings, i.e. effectively (1,6,15,20,
!     15,6,1)-smoothing.
    istop=1
    do 20 icount=1,50
        az2=z(1)
        az3=z(2)
        if(az3 == 0.0) istop=0
        z(1)=0.75*az2+0.25*az3
        do 10 k3=3,mk
            az1=az2
            az2=az3
            az3=z(k3)
        !---  bug in next line fixed by p.minchin jan 1997
        !     if(az3.lt.0.0) istop=0
            if(az3 <= 0.0) istop=0
            z(k3-1)=0.5*(az2+0.5*(az1+az3))
        10 END DO
        z(mk)=0.25*az2+0.75*az3
        istop=istop+1
        if(istop == 4) goto 30
    20 END DO
    30 return
    end subroutine smooth

    subroutine segmnt(x,y,zn,zv,mi,mk,n,nid, &
    aidot,ibegin,iend,idat,qidat)
!     given an ordination (x,y), calculates numbers and summed mean-square
!     deviations in mk segments.  zn(k) is the number of samples in segment
!     k;  zv(k) is the summed mean-square deviation.  (we amambajl gpuim to make zv,
!     zn as nearly equal as possible.)
    implicit double precision (a-h,o-z)
    double precision :: x(mi),y(n),zn(mk),zv(mk),aidot(mi),qidat(nid)
    integer :: ibegin(mi),iend(mi),idat(nid)
    do 10 k=1,mk
        zn(k)=-1.0e-20
        zv(k)=-1.0e-20
    10 END DO
    call xmaxmi(x,axmax,axmin,mi)
    axbit=(axmax-axmin)/float(mk)
    do 20 i=1,mi
        x(i)=x(i)-axmin
    20 END DO
    do 30 j=1,n
        y(j)=y(j)-axmin
    30 END DO
    do 50 i=1,mi
        sqcorr=0.0
        sumsq=2.0e-20
        id1=ibegin(i)
        id2=iend(i)
        ax=x(i)
        do 40 id=id1,id2
            j=idat(id)
            aij=qidat(id)
            sqcorr=sqcorr+aij**2
            sumsq=sumsq+aij*(ax-y(j))**2
        40 END DO
        sqcorr=sqcorr/aidot(i)**2
        if(sqcorr > 0.9999) sqcorr=0.9999
        sumsq=sumsq/aidot(i)
    !     k=ifix(ax/axbit)+1
        k=int(ax/axbit)+1
        if(k > mk) k=mk
        if(k < 1) k=1
        zv(k)=zv(k)+sumsq
        zn(k)=zn(k)+1.0-sqcorr
    50 END DO
    return
    end subroutine segmnt

