!============================================!
subroutine read_BG
use setup
use getBG

character*20::tfile='SPARC/temp.ascii'
character*20::ufile='SPARC/wind.ascii'
character*50::ufile_tide='SPARC/plus_tide/wind_tide.ascii'
real(4)::t_sparc(12,41,33)
real(4)::u_sparc(12,41,46)

real(4)::u_sparc_tide_xyzt(nlon,41,46,ntime)

select case(trim(cbgmode))
 case('test')
  !!! constant wind and stability 

   u0   = 0.0
   bn20 = 0.018**2
   tbot = 300.0

   ubg_xyzt = u0
   vbg_xyzt = 0.0

   tbg_xyzt(:,:,1,:)=tbot

   do ilev=2,nlev !!! forward Euler integration
    tbg_xyzt(:,:,ilev,:) = tbg_xyzt(:,:,ilev-1,:)+(zlev(ilev)-zlev(ilev-1))*(bn20*tbg_xyzt(:,:,ilev-1,:)/grav-grav/cp)   
   end do

 case('SPARC_July')

  if (nlat.ne.41.or.nlev.ne.33)then
   write(*,*) 'ERROR:: reading SPARC clim -- dimension must be nlat=41 & nlev=33'
   stop
  end if

  imonth=7 !!! July

  open(13,file=trim(tfile),form='formatted')
  read(13,'(10f8.3)')t_sparc
  close(13)
  open(13,file=trim(ufile),form='formatted')
  read(13,'(8f10.3)')u_sparc
  close(13)

  do ilon=1,nlon
  do itime=1,ntime
   tbg_xyzt(ilon,1:41,1:33,itime)=t_sparc(imonth,1:41,1:33)
   ubg_xyzt(ilon,1:41,1:33,itime)=u_sparc(imonth,1:41,1:33)
  end do
  end do
  vbg_xyzt=0.0

  case('SPARC_July_tide')

  if (nlat.ne.41.or.nlev.ne.33)then
   write(*,*) 'ERROR:: reading SPARC clim -- dimension must be nlat=41 & nlev=33'
   stop
  end if

  imonth=7 !!! July

  open(13,file=trim(tfile),form='formatted')
  read(13,'(10f8.3)')t_sparc
  close(13)
  open(13,file=trim(ufile),form='formatted')
  read(13,'(8f10.3)')u_sparc
  close(13)

  do ilon=1,nlon
  do itime=1,ntime
   tbg_xyzt(ilon,1:41,1:33,itime)=t_sparc(imonth,1:41,1:33)
   ubg_xyzt(ilon,1:41,1:33,itime)=u_sparc(imonth,1:41,1:33)
  end do
  end do
  vbg_xyzt=0.0

  open(13,file=trim(ufile_tide),form='formatted')
   read(13,'(8f10.3)')u_sparc_tide_xyzt
  close(13)

  ubg_xyzt = ubg_xyzt + u_sparc_tide_xyzt(:,:,1:33,:)

 case('file')
  
  !!! read Ubg_xyzt, Vbg_xyzt and Tbg_xyzt from an external file --- EDIT ME 

end select

!!!call diagnose_BG


 write(*,*) 'U max, min',maxval(ubg),minval(ubg)
 write(*,*) 'V max, min',maxval(vbg),minval(vbg)
 write(*,*) 'T max, min',maxval(tbg),minval(tbg)
 write(*,*) 'N2 max, min',maxval(bn2),minval(bn2)


!!! Initialize BG spline

!!!! call getBG_init


return
end subroutine read_BG
!============================================!
 subroutine diagnose_BG
 use setup
  
 real(4)::tbgdz(nlon,nlat,nlev)
 real(4)::tbgdzt(nlon,nlat,nlev)

  do ilon=1,nlon
  do ilat=1,nlat
   call grad1D(tbgdz(ilon,ilat,:),tbg(ilon,ilat,:),zlev,nlev)
   call grad1D(tbgdzt(ilon,ilat,:),Ttbg(ilon,ilat,:),zlev,nlev)
  end do
  end do

  bn2 = grav * (tbgdz+grav/cp) /tbg
  alp = 0.5 * grav / rd / tbg  * ishe
  bn2t = bn2 * Ttbg / tbg
  alpt = - alp * Ttbg / tbg

  alp0 = 0.5 * grav / rd / 250.0 * ishe !!! TORI AEZU 
  if (ialp.eq.0) alp = alp0 !!! KANTAN 


 return
 end subroutine diagnose_BG

!============================================!
subroutine read_init
use setup


 rayx = rmiss
 rayy = rmiss
 rayz = rmiss
 rayk = rmiss
 rayl = rmiss
 raym = rmiss
 rayw = rmiss
 rayh = rmiss
 rayb = rmiss
 rayt = rmiss
 rayi = imiss


open(12,file='init.dat',form='formatted')
 read(12,*)  !!! 1st raw is for a header

nray=0
 do iray=1,nraymax
  read(12,*,end=999) rayx(1,iray),rayy(1,iray),rayz(1,iray),rayt(1,iray),rayk(1,iray),rayl(1,iray),raym(1,iray),rayw(1,iray)
  nray=iray !!! number of rays to calculate (nray =< nraymax)
 end do
999 close(12)

 write(*,*) '# of rays =',nray


 call diagnose_init

 call amplitude_init 

return
end subroutine read_init

!============================================!

subroutine diagnose_init
use setup
use getBG

!!! diagnose omega(hat) 
!!! and diagnose m or omega if needed

 if (iupwd.ne.1) iupwd = -1 !!! downward

 do iray=1,nray
   call getBG_val(vu,vv,vn2,valp,rayx(1,iray),rayy(1,iray),rayz(1,iray),rayt(1,iray))
   if (trim(cdiag).eq.'m')then
    rayh(1,iray)=rayw(1,iray)-rayk(1,iray)*vu-rayl(1,iray)*vv
    vm2 = (rayk(1,iray)**2+rayl(1,iray)**2)*(vn2-real(inhd)*rayh(1,iray)**2)/(rayh(1,iray)**2-vf(rayy(1,iray))**2) - valp**2
     if (vm2.le.0.0)then
      write(*,*) 'wave # ',iray,' is invalid :: removed'
      rayi(1,iray)=-1
     else
      raym(1,iray)=sign(sqrt(vm2),real(iupwd)*(-1)*rayh(1,iray)) !!! dependent of Cgz sign
      rayi(1,iray) = 1
     end if
   elseif (trim(cdiag).eq.'omega')then
    wh2 = (vn2*(rayk(1,iray)**2+rayl(1,iray)**2)+vf(rayy(1,iray))**2*(raym(1,iray)**2+valp**2)) &
        / (real(inhd)*(rayk(1,iray)**2+rayl(1,iray)**2)+raym(1,iray)**2+valp**2)
     if (wh2.ge.vn2.or.wh2.le.vf(rayy(1,iray)))then
      write(*,*) 'wave # ',iray,' is invalid :: removed'
      rayi(1,iray)=-1
     else
      rayh(1,iray)=sign(sqrt(wh2),real(iupwd)*(-1)*raym(1,iray)) !!! dependent of Cgz sign
      rayw(1,iray)=rayh(1,iray)+rayk(1,iray)*vu+rayl(1,iray)*vv
      rayi(1,iray) = 1
     end if
   else !!! no need to diagnose
      rayi(1,iray) = 1
   end if

 end do

return
end subroutine diagnose_init

!============================================!

subroutine amplitude_init
use setup

 if (iamp.eq.1) then

  !!! amplitude diagnose -- TORI AEZU IRANAI  

 else 
  !!! amplitudes are assumed to be sufficiently small 
  !!! and nonlinear-breaking criterion is not accounted (except for critical level)
  rayb(1,1:nray)=1.0e-10
 end if

return
end subroutine amplitude_init

!============================================!

subroutine grad1D(arrayout,arrayin,arrayax,ndim)
use EZspline_obj
use EZspline


!!! arrayax must be ascending

integer::ndim
real(4)::arrayout(ndim)
real(4)::arrayin(ndim)
real(4)::arrayax(ndim)
integer::ibcs1(2)=(/0,0/)
type(EZspline1_r4)::spl

call EZspline_init(spl,ndim,ibcs1,ier)
call EZspline_error(ier)
spl%x1=arrayax
call EZspline_setup(spl, arrayin, ier)
call EZspline_error(ier)
call EZspline_derivative(spl,1,ndim,arrayax,arrayout, ier)
call EZspline_error(ier)

call EZspline_free(spl,ier)


return
end subroutine grad1D

!============================================!

subroutine output
use setup

!!! EDIT ME !!! 

return
end subroutine output
!============================================!
