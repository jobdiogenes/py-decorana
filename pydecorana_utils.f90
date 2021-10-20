module decorana_interface
use iso_c_binding, only: c_double, c_int
use vetorana_module, only: cutup,yxmult,eigy
implicit none
contains
subroutine c_cutup(x,ix,mi,mk)
    real(c_double), dimension(mi), intent(in out) :: x
    integer(c_int), dimension(mi), intent(out) :: ix
    integer(c_int), intent(in) :: mi
    integer(c_int), intent(in) :: mk 
    call cutup(x,ix,mi,m) 
end subroutine

subroutine c_yxmult(y,x,mi,n,nid,ibegin,iend,idat,qidat)
    real(c_double), dimension(mi), intent(in) :: x
    real(c_double), dimension(n), intent(out) :: y
    integer(c_int), intent(in) :: mi
    integer(c_int), intent(in) :: n
    integer(c_int), intent(in out) :: nid
    integer(c_int), dimension(mi), intent(in) :: ibegin
    integer(c_int), dimension(mi), intent(in) :: iend
    integer(c_int), dimension(nid), intent(in) :: idat
    double precision, dimension(nid), intent(in) :: qidat
    call yxmult(y,x,mi,n,nid,ibegin,iend,idat,qidat)
end subroutine

subroutine c_eigy(x,y,eig,neig,ira,iresc,short,  &
    mi,mk,n,nid,ibegin,iend,idat,qidat,y2,y3,y4,y5,  &
    xeig1,xeig2,xeig3,ix1,ix2,ix3,aidot,adotj)
    real(c_double), intent(out)    :: x(mi)
    real(c_double), intent(out)    :: y(n)
    real(c_double), intent(out)    :: eig
    integer(c_int), intent(in out) :: neig
    integer(c_int), intent(in out) :: ira
    integer(c_int), intent(in)     :: iresc
    real(c_double), intent(in out) :: short
    integer(c_int), intent(in)     :: mi
    integer(c_int), intent(in out) :: mk
    integer(c_int), intent(in)     :: n
    integer(c_int), intent(in out) :: nid
    integer(c_int), intent(in)     :: ibegin(mi)
    integer(c_int), intent(in)     :: iend(mi)
    integer(c_int), intent(in)     :: idat(nid)
    real(c_double), intent(in)     :: qidat(nid)
    real(c_double), intent(in out) :: y2(n)
    real(c_double), intent(in out) :: y3(n)
    real(c_double), intent(in out) :: y4(n)
    real(c_double), intent(in)     :: y5(n)
    real(c_double), intent(in out) :: xeig1(mi)
    real(c_double), intent(in out) :: xeig2(mi)
    real(c_double), intent(in out) :: xeig3(mi)
    integer(c_int), intent(in out) :: ix1(mi)
    integer(c_int), intent(in out) :: ix2(mi)
    integer(c_int), intent(in out) :: ix3(mi)
    real(c_double), intent(in)     :: aidot(mi)
    real(c_double), intent(in)     :: adotj(n)
    call eigy(x,y,eig,neig,ira,iresc,short,mi,mk,n,nid,ibegin,iend,idat,qidat,y2,y3,y4,y5,xeig1,xeig2,xeig3,ix1,ix2,ix3,aidot,adotj)
end subroutine
end module