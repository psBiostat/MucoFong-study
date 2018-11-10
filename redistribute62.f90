! File: redistribute6.f
! Version: 1.0
! Date: 2016-06-17
! Author: Stephen Rush <srush01@uoguelph.ca>
! Maintainer: Stephen Rush <srush01@uoguelph.ca>
! Description: Subroutine to redistribute mass among the group and
! individual coefficients according to mass equilbrium result

subroutine redistribute62(taxo, itax, dol, aol, ane, ntl, ni, nta)
implicit none
integer, parameter :: dp = kind(1.d0)
integer :: ntl, ni, nta, i, j, k, jerr
integer :: taxo(ni,ntl), itax(2,ntl-1)
real :: att, dtt
real(kind=dp) :: dol(nta), aol(ni), ane(ni)
integer, dimension(:), allocatable :: cind
allocate(cind(1:ni), stat = jerr)
aol = sqrt(abs(aol) * abs(ane))
! Determine new alpha
do i = (ntl-1),1,-1
do j = itax(1,i),itax(2,i)
att = 0
do k = 1,ni
if (taxo(k,i) /= j) cycle
att = att + aol(k)
enddo
att = max(att ** (1.0/(ntl-i+2)), 1e-20)
dtt = dol(j) ** (1.0/(ntl-i+2))
do k =1,ni
if (taxo(k,i) /= j) cycle
aol(k) = aol(k) * dtt / att
enddo
enddo
enddo
ane = sign(aol, ane)
! Build new group coefficients
do i = 1,(ntl-1)
do j = itax(1,i),itax(2,i)
dtt = 0
do k = 1,ni
if (taxo(k,i) /= j) cycle
dtt = dtt + aol(k)
enddo
dol(j) = dtt
enddo
enddo
deallocate(cind)
return
end subroutine redistribute62