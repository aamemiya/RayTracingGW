!!! Grid configuration !!!


integer,parameter::nlon=128,nlat=64,nlev=80

real(4),parameter::axlatl=-90.0,axlatr=90.0
real(4),parameter::axlonl=0.0,axlonr=360.0
real(4),parameter::zbot=0.0e3,ztop=90.0e3


real(4),parameter::axlatSN(nlat) &
  = (/(axlatl+(real(ilat)-0.5)*(axlatr-axlatl)/real(nlat),ilat=1,nlat)/)
real(4),parameter::axlon(nlon)   &
  = (/(axlonl+(real(ilon)-1.0)*(axlonr-axlonl)/real(nlon),ilon=1,nlon)/)
real(4),parameter::zlev(nlev)    &
                 = (/(zbot+(real(ilev)-0.5)*(ztop-zbot)/real(nlev),ilev=1,nlev)/)
