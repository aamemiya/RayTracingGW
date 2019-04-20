
!============================================!

subroutine quicklook_zm
use setup

real(4)::varyz(nlat,nlev)
real(4)::varxy(nlon,nlat)
character*3::cvar
character*20::title1,title2
namelist /quicklook_nml/ iout, cvar, title1, title2


open(11,file='calc.conf')
 read(11,nml=quicklook_nml)
close(11)

select case(trim(cvar))
case ('u') 
 varyz=sum(ubg,1)/real(nlon)
 varmx=100.0
 varmn=-100.0
 dvar=10.0
 write(*,*) 'BG U max, min', maxval(varyz), minval(varyz)

case('t')
 varyz=sum(tbg,1)/real(nlon)
 varmx=300.0
 varmn=160.0
 dvar=10.0
 write(*,*) 'BG T max, min', maxval(varyz), minval(varyz)

case('bn2')
 varyz=sum(bn2,1)/real(nlon)
 varmx=1.0e3
 varmn=0.0
 dvar=1.0e4
 write(*,*) 'BG N2 max, min', maxval(varyz), minval(varyz)
 
end select


! output
iouts=1
ioute=2
if(iout.eq.1) ioute=1
if(iout.eq.2) iouts=2
do ioutl=iouts,ioute
! *** general settings ***
      call sgiset ('IFONT',1)
      call swistx ('ICLRMAP',14) ! colormap blue-white-red 
      call swcmll
      call swcset ('FNAME','figure')
      call swiset ('IWIDTH',1000) ! output window size
      call swiset ('IHEIGHT',800)     
      call gropn(ioutl) 
      call sglset ('LFULL',.TRUE.) ! using fullsize
      call swlset ('LSEP',.TRUE.) ! psfile numbering
      call slmgn (0.0,0.0,0.0,0.0) ! margin 
      call grfrm 
      call grswnd (axlatl,axlatr,1.0e-3*zbot,1.0e-3*ztop) ! set window
      call grsvpt (0.15,0.85,0.13,0.73) ! set viewport
      call grstrn (1) ! linear or log
      call grstrf


! **** contour ****

      call udlset ('LMSG',.FALSE.) ! 'contour interval'
      call udrset ('RSIZEL',0.015) ! label
      call uwsgxa (axlatSN,nlat) ! yaxis value
      call uwsgya (1.0e-3*zlev,nlev) ! zaxis value (km)
      call udiset ('INDXMJ',4) ! contour line index
      call udiset ('INDXMN',2) ! contour line index (minor)
      call udsfmt ('B') ! contour label format
      call udgcla (varmn,varmx,dvar) ! contour level
      call udcntr (varyz,nlat,nlat,nlev) ! draw contour 

! *** lines on xy plane ***

  call uuslnt(1) ! line type
  call uuslni(45) ! line index (color*10 + index*1)
  do iray=1,nray
   do itime=1,ntmax
    if (rayi(itime,iray).le.0)goto 199
   end do    
199 ncount=itime-1
   do itime=1,ncount
   end do
   call uulin (ncount,rayy(1:ncount,iray),0.001*rayz(1:ncount,iray))
  end do

! **** x ,y axis ****



      call uziset ('INDEXT2',3)
      call uziset ('INDEXL1',5)
      call uzrset ('RSIZEL1',0.018)
      call uzrset ('RSIZEC1',0.018)
      call uzrset ('RSIZET2',0.010)
      call uzrset ('RSIZET1',0.004)
!      call uzlset ('LABELXB',.FALSE.)
      call uxaxdv ('B',10.0,30.0)
      call uxaxdv ('T',10.0,30.0)
      call uxsttl ('B','Latitude',0.0)
      call uysfmt ('(I2)')
      call uzlset ('LABELYR',.FALSE.)
!     call uzlset ('LABELYL',.FALSE.)
      call uyaxdv ('L',10.0,20.0)
      call uyaxdv ('R',10.0,20.0)
      call uziset ('IROTCYL',1)
      call uysttl ('L','Height(km)',0.0)

      call sgtxzv (0.50,0.76,trim(title1),0.024,0,0,5) ! title
      call sgtxzv (0.17,0.76,trim(title2),0.020,0,-1,5) ! title
      call grcls
end do

return
end subroutine quicklook_zm

!============================================!

subroutine quicklook_pl
use setup

real(4)::varxy(nlon,nlat)
character*3::cvar
character*20::title1,title2
namelist /quicklook_nml/ iout, cvar, title1, title2


open(11,file='calc.conf')
 read(11,nml=quicklook_nml)
close(11)

! output
iouts=1
ioute=2
if(iout.eq.1) ioute=1
if(iout.eq.2) iouts=2
do ioutl=iouts,ioute
! *** general settings ***
      call sgiset ('IFONT',1)
      call swistx ('ICLRMAP',14) ! colormap blue-white-red 
      call swcmll
      call swcset ('FNAME','figure')
      call swiset ('IWIDTH',1000) ! output window size
      call swiset ('IHEIGHT',800)     
      call gropn(ioutl) 
      call sglset ('LFULL',.TRUE.) ! using fullsize
      call swlset ('LSEP',.TRUE.) ! psfile numbering
      call slmgn (0.0,0.0,0.0,0.0) ! margin 
      call grfrm 

      call grstrn (30) ! linear or log

      call umiset ('INDEXOUT',32) !!! HOSOI MIDORI

      call umscwd (rayx(1,1),rayy(1,1),40.0)

      call grswnd (0.0,360.0,-90.0,90.0) ! set window
      call grsvpt (0.20,0.60,0.20,0.60) ! set viewport
      call umpfit
      call grstrf
      call umpglb
      call umpmap('coast_world')



! *** lines on xy plane ***

  call uuslnt(1) ! line type
  call uuslni(45) ! line index (color*10 + index*1)
  do iray=1,nray
   do itime=1,ntmax
    if (rayi(itime,iray).le.0)goto 199
   end do    
199 ncount=itime-1
   do itime=1,ncount
   end do
   call uulin (ncount,rayx(1:ncount,iray),rayy(1:ncount,iray))
  end do

! **** x ,y axis ****

      call sgtxzv (0.50,0.76,trim(title1),0.024,0,0,5) ! title
      call sgtxzv (0.17,0.76,trim(title2),0.020,0,-1,5) ! title
      call grcls
end do


return
end subroutine quicklook_pl
