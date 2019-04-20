program generate


integer,parameter::nk_tide=2
integer,parameter::nw_tide=2
real(4),parameter::amp_tide=20.0

real(4)::prof_amp(46)

integer,parameter::nday=10
integer,parameter::ntime=241

real(4),parameter::axtime(ntime)=(/(0.0 + 86400.0 * real(nday)/real(ntime-1) * real(itime-1) , itime=1,ntime )/)

character*20::ufile_t='./wind_tide.ascii'
real(4)::u_sparc_tide_xyzt(64,41,46,ntime)


pi = 4.0 * atan(1.0)

prof_amp( 1:20) = 0.0
prof_amp(21:46) = 1.0

do itime=1,ntime
do ix = 1,64
 phase = 2.0 * pi * ( real(nk_tide)*real(ix-1)/64.0 + real(nw_tide)*axtime(itime)/86400.0)
 do ilev=1,46
  u_sparc_tide_xyzt(ix,:,ilev,itime) =  amp_tide * prof_amp(ilev) * sin(phase) 
 end do
end do
end do

  open(14,file=trim(ufile_t),form='formatted')
  write(14,'(8f10.3)')u_sparc_tide_xyzt
  close(14)

stop
end program generate
