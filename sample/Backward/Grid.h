!!! Grid configuration !!!
!!! For SPARC reference climatology


integer,parameter::nlon=64,nlat=41,nlev=33

real(4),parameter::axlatl=-80.0,axlatr=80.0
real(4),parameter::axlonl=0.0,axlonr=360.0
real(4),parameter::zbot=0.0,ztop=32.0*7.0e3*log(10.0)/6.0

real(4),parameter::axlatSN(nlat) &
  = (/(axlatl+real(ilat-1)*(axlatr-axlatl)/real(nlat-1),ilat=1,nlat)/)
real(4),parameter::axlon(nlon)   &
  = (/(axlonl+(real(ilon-1))*(axlonr-axlonl)/real(nlon),ilon=1,nlon)/)
real(4),parameter::zlev(nlev)    &
  = (/(real(ilev-1)*7.0e3*log(10.0)/6.0,ilev=1,nlev)/)                
