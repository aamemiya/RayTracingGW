!===============================================!
module getBG
use EZspline_obj
use EZspline

! matrices for EZspline
type(EZspline3_r4)::u_spl,v_spl,n2_spl,alp_spl
type(EZspline3_r4)::ut_spl,vt_spl,n2t_spl,alpt_spl
integer,parameter::ibcs1(2)=(/-1, -1/) ! periodic
integer,parameter::ibcs2(2)=(/0 , 0 /)
integer,parameter::ibcs3(2)=(/0 , 0 /)


integer::ilonc,ilatc,ilevc
real(4)::fact_xint,fact_yint,fact_zint

contains!---------------------------------------!

subroutine getBG_snap(time)
use setup

real(4),save::time_b=-99999.9

if (time.eq.time_b) return !!! no need to update

time_b = time

if (ntime.eq.1)then

 ubg = ubg_xyzt(:,:,:,1)
 vbg = vbg_xyzt(:,:,:,1)
 Tbg = Tbg_xyzt(:,:,:,1)

 utbg = 0.0
 vtbg = 0.0
 Ttbg = 0.0

else

itimeb = iblkge(axtime,ntime,time)
if(itimeb.le.0)then
 write(*,*) 'getBG_snap : time',time, '< axtime(1)',axtime(1)  
 stop
! itimeb=1
! fact_tint=0.0
elseif(itimeb.ge.ntime)then
 write(*,*) 'getBG_snap : time',time, '> axtime(ntime)',axtime(ntime)  
 stop
! itimeb=ntime-1
! fact_tint=1.0
else

fact_tint = (time - axtime(itimeb)) / (axtime(itimeb+1)-axtime(itimeb))
end if

ubg = (1.0-fact_tint) * ubg_xyzt(:,:,:,itimeb) + fact_tint * ubg_xyzt(:,:,:,itimeb+1)
vbg = (1.0-fact_tint) * vbg_xyzt(:,:,:,itimeb) + fact_tint * vbg_xyzt(:,:,:,itimeb+1)
Tbg = (1.0-fact_tint) * Tbg_xyzt(:,:,:,itimeb) + fact_tint * Tbg_xyzt(:,:,:,itimeb+1)

!!! simple linear differetiation
 Utbg = (ubg_xyzt(:,:,:,itimeb+1) - ubg_xyzt(:,:,:,itimeb)) / (axtime(itimeb+1)-axtime(itimeb))
 Vtbg = (vbg_xyzt(:,:,:,itimeb+1) - vbg_xyzt(:,:,:,itimeb)) / (axtime(itimeb+1)-axtime(itimeb))
 Ttbg = (tbg_xyzt(:,:,:,itimeb+1) - tbg_xyzt(:,:,:,itimeb)) / (axtime(itimeb+1)-axtime(itimeb))

end if

!!! prepare N2 and alp
call diagnose_BG

return
end subroutine getBG_snap
!---------------------------------------!

subroutine getBG_prep(vx,vy,vz,time)
use setup

!!! For cyclic condition
real(4)::axlon_ext(nlon+1) 
real(4)::ubg_ext(nlon+1,nlat,nlev) 
real(4)::vbg_ext(nlon+1,nlat,nlev) 
real(4)::bn2_ext(nlon+1,nlat,nlev) 
real(4)::alp_ext(nlon+1,nlat,nlev) 
integer::ibcsx(2)

 axlon_ext(1:nlon)=axlon
 axlon_ext(nlon+1)=360.0

 ilonc=iblkge(axlon_ext,nlon+1,vx)
 ilatc=iblkge(axlatSN,nlat,vy)
 ilevc=iblkge(zlev,nlev,vz)
 fact_xint = (vx-axlon(ilonc))   / (axlon_ext(ilonc+1)-axlon(ilonc)) 
 fact_yint = (vy-axlatSN(ilatc)) / (axlatSN(ilatc+1)-axlatSN(ilatc)) 
 fact_zint = (vz-zlev(ilevc))    / (zlev(ilevc+1)-zlev(ilevc)) 


call getbg_snap(time)

select case(trim(cintpmode))
case('Linear','linear')

case('Spline','spline')
! nearest grid point
if (fact_xint.gt.0.5) ilonc=ilonc+1
if (fact_yint.gt.0.5) ilatc=ilatc+1
if (fact_zint.gt.0.5) ilevc=ilevc+1

if (ilatsp.le.0) then !!! merid. global
 ilatl=1
 ilatr=nlat
 nlatsp=nlat
else
 ilatl=max(ilatc-ilatsp,1)
 ilatr=min(ilatc+ilatsp,nlat)
 nlatsp=ilatr-ilatl+1
end if

if (ilatsp.le.0) then !!! vert. global
 ilevb=1
 ilevt=nlev
 nlevsp=nlev
else
 ilevb=max(ilevc-ilevsp,1)
 ilevt=min(ilevc+ilevsp,nlev)
 nlevsp=ilevt-ilevb+1
end if

if (ilonsp.le.0) then !!! zonal global (periodic)
 ibcsx=ibcs1
 nlonsp=nlon+1
 ubg_ext(1:nlon,1:nlatsp,1:nlevsp)=ubg
 vbg_ext(1:nlon,1:nlatsp,1:nlevsp)=vbg
 bn2_ext(1:nlon,1:nlatsp,1:nlevsp)=bn2
 alp_ext(1:nlon,1:nlatsp,1:nlevsp)=alp
 ubg_ext(nlon+1,1:nlatsp,1:nlevsp)=ubg(1,ilatl:ilatr,ilevb:ilevt)
 vbg_ext(nlon+1,1:nlatsp,1:nlevsp)=vbg(1,ilatl:ilatr,ilevb:ilevt)
 bn2_ext(nlon+1,1:nlatsp,1:nlevsp)=bn2(1,ilatl:ilatr,ilevb:ilevt)
 alp_ext(nlon+1,1:nlatsp,1:nlevsp)=alp(1,ilatl:ilatr,ilevb:ilevt)

else
 ibcsx=ibcs2 !!! non-periodic
 ilonl=ilonc-ilonsp
 ilonr=ilonc+ilonsp
 nlonsp=ilonr-ilonl+1
  if (ilonl.le.0)then
   axlon_ext(1:1-ilonl)=axlon(nlon+ilonl:nlon)-360.0
   axlon_ext(2-ilonl:nlonsp)=axlon(1:nlonsp+ilonl+1)
   ubg_ext(1:1-ilonl,1:nlatsp,1:nlevsp)=ubg(nlon+ilonl:nlon,ilatl:ilatr,ilevb:ilevt)
   vbg_ext(1:1-ilonl,1:nlatsp,1:nlevsp)=vbg(nlon+ilonl:nlon,ilatl:ilatr,ilevb:ilevt)
   bn2_ext(1:1-ilonl,1:nlatsp,1:nlevsp)=bn2(nlon+ilonl:nlon,ilatl:ilatr,ilevb:ilevt)
   alp_ext(1:1-ilonl,1:nlatsp,1:nlevsp)=alp(nlon+ilonl:nlon,ilatl:ilatr,ilevb:ilevt)
   ubg_ext(2-ilonl:nlonsp,1:nlatsp,1:nlevsp)=ubg(1:nlonsp+ilonl+1,ilatl:ilatr,ilevb:ilevt)
   vbg_ext(2-ilonl:nlonsp,1:nlatsp,1:nlevsp)=vbg(1:nlonsp+ilonl+1,ilatl:ilatr,ilevb:ilevt)
   bn2_ext(2-ilonl:nlonsp,1:nlatsp,1:nlevsp)=bn2(1:nlonsp+ilonl+1,ilatl:ilatr,ilevb:ilevt)
   alp_ext(2-ilonl:nlonsp,1:nlatsp,1:nlevsp)=alp(1:nlonsp+ilonl+1,ilatl:ilatr,ilevb:ilevt)
  elseif (ilonr.ge.nlon+1)then
   axlon_ext(1:nlon-ilonl+2) = axlon(ilonl-1:nlon)
   axlon_ext(nlon-ilonl+3:nlonsp)=axlon(1:ilonr-nlon)+360.0
   ubg_ext(1:nlon-ilonl+2,1:nlatsp,1:nlevsp)=ubg(ilonl-1:nlon,ilatl:ilatr,ilevb:ilevt)
   vbg_ext(1:nlon-ilonl+2,1:nlatsp,1:nlevsp)=vbg(ilonl-1:nlon,ilatl:ilatr,ilevb:ilevt)
   bn2_ext(1:nlon-ilonl+2,1:nlatsp,1:nlevsp)=bn2(ilonl-1:nlon,ilatl:ilatr,ilevb:ilevt)
   alp_ext(1:nlon-ilonl+2,1:nlatsp,1:nlevsp)=alp(ilonl-1:nlon,ilatl:ilatr,ilevb:ilevt)
   ubg_ext(nlon-ilonl+3:nlonsp,1:nlatsp,1:nlevsp)=ubg(1:ilonr-nlon,ilatl:ilatr,ilevb:ilevt)
   vbg_ext(nlon-ilonl+3:nlonsp,1:nlatsp,1:nlevsp)=vbg(1:ilonr-nlon,ilatl:ilatr,ilevb:ilevt)
   bn2_ext(nlon-ilonl+3:nlonsp,1:nlatsp,1:nlevsp)=bn2(1:ilonr-nlon,ilatl:ilatr,ilevb:ilevt)
   alp_ext(nlon-ilonl+3:nlonsp,1:nlatsp,1:nlevsp)=alp(1:ilonr-nlon,ilatl:ilatr,ilevb:ilevt)
  else
   axlon_ext(1:nlonsp)=axlon(ilonl:ilonr)
   ubg_ext(1:nlonsp,1:nlatsp,1:nlevsp)=ubg(ilonl:ilonr,ilatl:ilatr,ilevb:ilevt)
   vbg_ext(1:nlonsp,1:nlatsp,1:nlevsp)=vbg(ilonl:ilonr,ilatl:ilatr,ilevb:ilevt)
   bn2_ext(1:nlonsp,1:nlatsp,1:nlevsp)=bn2(ilonl:ilonr,ilatl:ilatr,ilevb:ilevt)
   alp_ext(1:nlonsp,1:nlatsp,1:nlevsp)=alp(ilonl:ilonr,ilatl:ilatr,ilevb:ilevt)
  end if


  if (any(axlon_ext(2:nlonsp)-axlon_ext(1:nlonsp-1).ne.axlon(2)-axlon(1)))then
   write(*,*) 'spline prep error: invalid axlon_ext',axlon_ext(1:nlonsp) 
   stop
  end if

end if

!write(*,*)nlonsp,nlatsp,nlevsp
!stop

 call EZspline_init(u_spl,nlonsp,nlatsp,nlevsp,ibcsx,ibcs2,ibcs3,ier)
 call EZspline_init(v_spl,nlonsp,nlatsp,nlevsp,ibcsx,ibcs2,ibcs3,ier)
 call EZspline_init(n2_spl,nlonsp,nlatsp,nlevsp,ibcsx,ibcs2,ibcs3,ier)
 call EZspline_init(alp_spl,nlonsp,nlatsp,nlevsp,ibcsx,ibcs2,ibcs3,ier)

 ! set grid location
 u_spl%x1=axlon_ext(1:nlonsp)
 v_spl%x1=axlon_ext(1:nlonsp)
 n2_spl%x1=axlon_ext(1:nlonsp)
 alp_spl%x1=axlon_ext(1:nlonsp)
 u_spl%x2=axlatSN(ilatl:ilatr)
 v_spl%x2=axlatSN(ilatl:ilatr)
 n2_spl%x2=axlatSN(ilatl:ilatr)
 alp_spl%x2=axlatSN(ilatl:ilatr)
 u_spl%x3=zlev(ilevb:ilevt)
 v_spl%x3=zlev(ilevb:ilevt)
 n2_spl%x3=zlev(ilevb:ilevt)
 alp_spl%x3=zlev (ilevb:ilevt)

 call EZspline_setup(u_spl,   ubg_ext, ier)
 call EZspline_setup(v_spl,   vbg_ext, ier)
 call EZspline_setup(n2_spl,  bn2_ext, ier)
 call EZspline_setup(alp_spl, alp_ext, ier)

end select

return
end subroutine getBG_prep
!---------------------------------------!

subroutine getBG_prep_tend(vx,vy,vz,time)
use setup

!!! For cyclic condition
real(4)::axlon_ext(nlon+1) 
real(4)::utbg_ext(nlon+1,nlat,nlev) 
real(4)::vtbg_ext(nlon+1,nlat,nlev) 
real(4)::bn2t_ext(nlon+1,nlat,nlev) 
real(4)::alpt_ext(nlon+1,nlat,nlev) 
integer::ibcsx(2)

 axlon_ext(1:nlon)=axlon
 axlon_ext(nlon+1)=360.0

 ilonc=iblkge(axlon_ext,nlon+1,vx)
 ilatc=iblkge(axlatSN,nlat,vy)
 ilevc=iblkge(zlev,nlev,vz)
 fact_xint = (vx-axlon(ilonc))   / (axlon_ext(ilonc+1)-axlon(ilonc)) 
 fact_yint = (vy-axlatSN(ilatc)) / (axlatSN(ilatc+1)-axlatSN(ilatc)) 
 fact_zint = (vz-zlev(ilevc))    / (zlev(ilevc+1)-zlev(ilevc)) 


call getbg_snap(time)

select case(trim(cintpmode))
case('Linear','linear')

case('Spline','spline')
! nearest grid point
if (fact_xint.gt.0.5) ilonc=ilonc+1
if (fact_yint.gt.0.5) ilatc=ilatc+1
if (fact_zint.gt.0.5) ilevc=ilevc+1

if (ilatsp.le.0) then !!! merid. global
 ilatl=1
 ilatr=nlat
 nlatsp=nlat
else
 ilatl=max(ilatc-ilatsp,1)
 ilatr=min(ilatc+ilatsp,nlat)
 nlatsp=ilatr-ilatl+1
end if

if (ilatsp.le.0) then !!! vert. global
 ilevb=1
 ilevt=nlev
 nlevsp=nlev
else
 ilevb=max(ilevc-ilevsp,1)
 ilevt=min(ilevc+ilevsp,nlev)
 nlevsp=ilevt-ilevb+1
end if

if (ilonsp.le.0) then !!! zonal global (periodic)
 ibcsx=ibcs1
 nlonsp=nlon+1
 utbg_ext(1:nlon,1:nlatsp,1:nlevsp)=utbg
 vtbg_ext(1:nlon,1:nlatsp,1:nlevsp)=vtbg
 bn2t_ext(1:nlon,1:nlatsp,1:nlevsp)=bn2t
 alpt_ext(1:nlon,1:nlatsp,1:nlevsp)=alpt
 utbg_ext(nlon+1,1:nlatsp,1:nlevsp)=utbg(1,ilatl:ilatr,ilevb:ilevt)
 vtbg_ext(nlon+1,1:nlatsp,1:nlevsp)=vtbg(1,ilatl:ilatr,ilevb:ilevt)
 bn2t_ext(nlon+1,1:nlatsp,1:nlevsp)=bn2t(1,ilatl:ilatr,ilevb:ilevt)
 alpt_ext(nlon+1,1:nlatsp,1:nlevsp)=alpt(1,ilatl:ilatr,ilevb:ilevt)

else
 ibcsx=ibcs2 !!! non-periodic
 ilonl=ilonc-ilonsp
 ilonr=ilonc+ilonsp
 nlonsp=ilonr-ilonl+1
  if (ilonl.le.0)then
   axlon_ext(1:1-ilonl)=axlon(nlon+ilonl:nlon)-360.0
   axlon_ext(2-ilonl:nlonsp)=axlon(1:nlonsp+ilonl+1)
   utbg_ext(1:1-ilonl,1:nlatsp,1:nlevsp)=utbg(nlon+ilonl:nlon,ilatl:ilatr,ilevb:ilevt)
   vtbg_ext(1:1-ilonl,1:nlatsp,1:nlevsp)=vtbg(nlon+ilonl:nlon,ilatl:ilatr,ilevb:ilevt)
   bn2t_ext(1:1-ilonl,1:nlatsp,1:nlevsp)=bn2t(nlon+ilonl:nlon,ilatl:ilatr,ilevb:ilevt)
   alpt_ext(1:1-ilonl,1:nlatsp,1:nlevsp)=alpt(nlon+ilonl:nlon,ilatl:ilatr,ilevb:ilevt)
   utbg_ext(2-ilonl:nlonsp,1:nlatsp,1:nlevsp)=utbg(1:nlonsp+ilonl+1,ilatl:ilatr,ilevb:ilevt)
   vtbg_ext(2-ilonl:nlonsp,1:nlatsp,1:nlevsp)=vtbg(1:nlonsp+ilonl+1,ilatl:ilatr,ilevb:ilevt)
   bn2t_ext(2-ilonl:nlonsp,1:nlatsp,1:nlevsp)=bn2t(1:nlonsp+ilonl+1,ilatl:ilatr,ilevb:ilevt)
   alpt_ext(2-ilonl:nlonsp,1:nlatsp,1:nlevsp)=alpt(1:nlonsp+ilonl+1,ilatl:ilatr,ilevb:ilevt)
  elseif (ilonr.ge.nlon+1)then
   axlon_ext(1:nlon-ilonl+2) = axlon(ilonl-1:nlon)
   axlon_ext(nlon-ilonl+3:nlonsp)=axlon(1:ilonr-nlon)+360.0
   utbg_ext(1:nlon-ilonl+2,1:nlatsp,1:nlevsp)=utbg(ilonl-1:nlon,ilatl:ilatr,ilevb:ilevt)
   vtbg_ext(1:nlon-ilonl+2,1:nlatsp,1:nlevsp)=vtbg(ilonl-1:nlon,ilatl:ilatr,ilevb:ilevt)
   bn2t_ext(1:nlon-ilonl+2,1:nlatsp,1:nlevsp)=bn2t(ilonl-1:nlon,ilatl:ilatr,ilevb:ilevt)
   alpt_ext(1:nlon-ilonl+2,1:nlatsp,1:nlevsp)=alpt(ilonl-1:nlon,ilatl:ilatr,ilevb:ilevt)
   utbg_ext(nlon-ilonl+3:nlonsp,1:nlatsp,1:nlevsp)=utbg(1:ilonr-nlon,ilatl:ilatr,ilevb:ilevt)
   vtbg_ext(nlon-ilonl+3:nlonsp,1:nlatsp,1:nlevsp)=vtbg(1:ilonr-nlon,ilatl:ilatr,ilevb:ilevt)
   bn2t_ext(nlon-ilonl+3:nlonsp,1:nlatsp,1:nlevsp)=bn2t(1:ilonr-nlon,ilatl:ilatr,ilevb:ilevt)
   alpt_ext(nlon-ilonl+3:nlonsp,1:nlatsp,1:nlevsp)=alpt(1:ilonr-nlon,ilatl:ilatr,ilevb:ilevt)
  else
   axlon_ext(1:nlonsp)=axlon(ilonl:ilonr)
   utbg_ext(1:nlonsp,1:nlatsp,1:nlevsp)=utbg(ilonl:ilonr,ilatl:ilatr,ilevb:ilevt)
   vtbg_ext(1:nlonsp,1:nlatsp,1:nlevsp)=vtbg(ilonl:ilonr,ilatl:ilatr,ilevb:ilevt)
   bn2t_ext(1:nlonsp,1:nlatsp,1:nlevsp)=bn2t(ilonl:ilonr,ilatl:ilatr,ilevb:ilevt)
   alpt_ext(1:nlonsp,1:nlatsp,1:nlevsp)=alpt(ilonl:ilonr,ilatl:ilatr,ilevb:ilevt)
  end if


  if (any(axlon_ext(2:nlonsp)-axlon_ext(1:nlonsp-1).ne.axlon(2)-axlon(1)))then
   write(*,*) 'spline prep error: invalid axlon_ext',axlon_ext(1:nlonsp) 
   stop
  end if

end if

!write(*,*)nlonsp,nlatsp,nlevsp
!stop

 call EZspline_init(ut_spl,nlonsp,nlatsp,nlevsp,ibcsx,ibcs2,ibcs3,ier)
 call EZspline_init(vt_spl,nlonsp,nlatsp,nlevsp,ibcsx,ibcs2,ibcs3,ier)
 call EZspline_init(n2t_spl,nlonsp,nlatsp,nlevsp,ibcsx,ibcs2,ibcs3,ier)
 call EZspline_init(alpt_spl,nlonsp,nlatsp,nlevsp,ibcsx,ibcs2,ibcs3,ier)

 ! set grid location
 ut_spl%x1=axlon_ext(1:nlonsp)
 vt_spl%x1=axlon_ext(1:nlonsp)
 n2t_spl%x1=axlon_ext(1:nlonsp)
 alpt_spl%x1=axlon_ext(1:nlonsp)
 ut_spl%x2=axlatSN(ilatl:ilatr)
 vt_spl%x2=axlatSN(ilatl:ilatr)
 n2t_spl%x2=axlatSN(ilatl:ilatr)
 alpt_spl%x2=axlatSN(ilatl:ilatr)
 ut_spl%x3=zlev(ilevb:ilevt)
 vt_spl%x3=zlev(ilevb:ilevt)
 n2t_spl%x3=zlev(ilevb:ilevt)
 alpt_spl%x3=zlev(ilevb:ilevt)

 call EZspline_setup(ut_spl,   utbg_ext, ier)
 call EZspline_setup(vt_spl,   vtbg_ext, ier)
 call EZspline_setup(n2t_spl,  bn2t_ext, ier)
 call EZspline_setup(alpt_spl, alpt_ext, ier)

end select

return
end subroutine getBG_prep_tend
!---------------------------------------!
subroutine getBG_val(vu,vv,vn2,va,vx,vy,vz,time)
use setup

real(4)::vu,vv,vn2,va
real(4)::vx,vy,vz,time

ier=0

select case(trim(cintpmode))
case('Linear','linear')
   vu =  (1.0-fact_xint) * (1.0-fact_yint) * (1.0-fact_zint) * ubg(ilonc,ilatc,ilevc)        &
       + (1.0-fact_xint) * (1.0-fact_yint) *      fact_zint  * ubg(ilonc,ilatc,ilevc+1)      &
       + (1.0-fact_xint) *      fact_yint  * (1.0-fact_zint) * ubg(ilonc,ilatc+1,ilevc)      &
       + (1.0-fact_xint) *      fact_yint  *      fact_zint  * ubg(ilonc,ilatc+1,ilevc+1)    &
       +      fact_xint  * (1.0-fact_yint) * (1.0-fact_zint) * ubg(mod(ilonc,nlon)+1,ilatc,ilevc)     &
       +      fact_xint  * (1.0-fact_yint) *      fact_zint  * ubg(mod(ilonc,nlon)+1,ilatc,ilevc+1)   &
       +      fact_xint  *      fact_yint  * (1.0-fact_zint) * ubg(mod(ilonc,nlon)+1,ilatc+1,ilevc)   &
       +      fact_xint  *      fact_yint  *      fact_zint  * ubg(mod(ilonc,nlon)+1,ilatc+1,ilevc+1) 
   vv =  (1.0-fact_xint) * (1.0-fact_yint) * (1.0-fact_zint) * vbg(ilonc,ilatc,ilevc)        &
       + (1.0-fact_xint) * (1.0-fact_yint) *      fact_zint  * vbg(ilonc,ilatc,ilevc+1)      &
       + (1.0-fact_xint) *      fact_yint  * (1.0-fact_zint) * vbg(ilonc,ilatc+1,ilevc)      &
       + (1.0-fact_xint) *      fact_yint  *      fact_zint  * vbg(ilonc,ilatc+1,ilevc+1)    &
       +      fact_xint  * (1.0-fact_yint) * (1.0-fact_zint) * vbg(mod(ilonc,nlon)+1,ilatc,ilevc)     &
       +      fact_xint  * (1.0-fact_yint) *      fact_zint  * vbg(mod(ilonc,nlon)+1,ilatc,ilevc+1)   &
       +      fact_xint  *      fact_yint  * (1.0-fact_zint) * vbg(mod(ilonc,nlon)+1,ilatc+1,ilevc)   &
       +      fact_xint  *      fact_yint  *      fact_zint  * vbg(mod(ilonc,nlon)+1,ilatc+1,ilevc+1) 
   vn2 = (1.0-fact_xint) * (1.0-fact_yint) * (1.0-fact_zint) * bn2(ilonc,ilatc,ilevc)        &
       + (1.0-fact_xint) * (1.0-fact_yint) *      fact_zint  * bn2(ilonc,ilatc,ilevc+1)      &
       + (1.0-fact_xint) *      fact_yint  * (1.0-fact_zint) * bn2(ilonc,ilatc+1,ilevc)      &
       + (1.0-fact_xint) *      fact_yint  *      fact_zint  * bn2(ilonc,ilatc+1,ilevc+1)    &
       +      fact_xint  * (1.0-fact_yint) * (1.0-fact_zint) * bn2(mod(ilonc,nlon)+1,ilatc,ilevc)     &
       +      fact_xint  * (1.0-fact_yint) *      fact_zint  * bn2(mod(ilonc,nlon)+1,ilatc,ilevc+1)   &
       +      fact_xint  *      fact_yint  * (1.0-fact_zint) * bn2(mod(ilonc,nlon)+1,ilatc+1,ilevc)   &
       +      fact_xint  *      fact_yint  *      fact_zint  * bn2(mod(ilonc,nlon)+1,ilatc+1,ilevc+1) 
   va =  (1.0-fact_xint) * (1.0-fact_yint) * (1.0-fact_zint) * alp(ilonc,ilatc,ilevc)        &
       + (1.0-fact_xint) * (1.0-fact_yint) *      fact_zint  * alp(ilonc,ilatc,ilevc+1)      &
       + (1.0-fact_xint) *      fact_yint  * (1.0-fact_zint) * alp(ilonc,ilatc+1,ilevc)      &
       + (1.0-fact_xint) *      fact_yint  *      fact_zint  * alp(ilonc,ilatc+1,ilevc+1)    &
       +      fact_xint  * (1.0-fact_yint) * (1.0-fact_zint) * alp(mod(ilonc,nlon)+1,ilatc,ilevc)     &
       +      fact_xint  * (1.0-fact_yint) *      fact_zint  * alp(mod(ilonc,nlon)+1,ilatc,ilevc+1)   &
       +      fact_xint  *      fact_yint  * (1.0-fact_zint) * alp(mod(ilonc,nlon)+1,ilatc+1,ilevc)   &
       +      fact_xint  *      fact_yint  *      fact_zint  * alp(mod(ilonc,nlon)+1,ilatc+1,ilevc+1) 

case('Spline','spline')
   call getBG_prep(vx,vy,vz,time)
   call EZspline_interp(u_spl,vx,vy,vz,vu,ier)
   call EZspline_interp(v_spl,vx,vy,vz,vv,ier)
   call EZspline_interp(n2_spl,vx,vy,vz,vn2,ier)
   call EZspline_interp(alp_spl,vx,vy,vz,valp,ier)
   call getBG_free
end select

if (ier.ne.0)then
 write(*,*)'ERROR in getBG_val ::ier=',ier
 stop
end if



return
end subroutine getBG_val
!---------------------------------------!
subroutine getBG_grad(vdu,vdv,vdn2,vda,idx,idy,idz,vx,vy,vz,time)

use setup, only: cintpmode
include 'Pcon.h'

real(4)::vdu,vdv,vdn2,vda
real(4)::vx,vy,vz,time
integer::idx,idy,idz

select case(trim(cintpmode))
case('Linear','linear')
 write(*,*) 'sorry, linear differentiation is not yet prepared. use spline. '
 stop
case('Spline','spline')

call getBG_prep(vx,vy,vz,time)

call EZspline_derivative(u_spl,idx,idy,idz,vx,vy,vz,vdu,ier)
if (ier.ne.0)then
 write(*,*)'ERROR in getBG_val u_spl ::ier=',ier
 stop
end if


call EZspline_derivative(v_spl,idx,idy,idz,vx,vy,vz,vdv,ier)
call EZspline_derivative(n2_spl,idx,idy,idz,vx,vy,vz,vdn2,ier)
call EZspline_derivative(alp_spl,idx,idy,idz,vx,vy,vz,vda,ier)

call getBG_free

if (ier.ne.0)then
 write(*,*)'ERROR in getBG_val u_spl ::ier=',ier
 stop
end if


!!! derivative factor correction ( to 1/m )  !!!
yfac=1.0/er/drad
zfac=1.0
xfac=yfac/cos(vy*drad)

vfac = xfac**idx * yfac**idy * zfac**idz

vdu  = vdu  * vfac
vdv  = vdu  * vfac
vdn2 = vdn2 * vfac
vda  = vda  * vfac

if (ier.ne.0)then
 write(*,*)'ERROR in getBG_val ::ier=',ier
 stop
end if

end select

return
end subroutine getBG_grad
!---------------------------------------!
subroutine getBG_tend(vdut,vdvt,vdn2t,vdat,vx,vy,vz,time)
use setup 

real(4)::vdut,vdvt,vdn2t,vdat
real(4)::vx,vy,vz,time


ier=0

select case(trim(cintpmode))
case('Linear','linear')
   vdut =  (1.0-fact_xint) * (1.0-fact_yint) * (1.0-fact_zint) * utbg(ilonc,ilatc,ilevc)        &
       + (1.0-fact_xint) * (1.0-fact_yint) *      fact_zint  * utbg(ilonc,ilatc,ilevc+1)      &
       + (1.0-fact_xint) *      fact_yint  * (1.0-fact_zint) * utbg(ilonc,ilatc+1,ilevc)      &
       + (1.0-fact_xint) *      fact_yint  *      fact_zint  * utbg(ilonc,ilatc+1,ilevc+1)    &
       +      fact_xint  * (1.0-fact_yint) * (1.0-fact_zint) * utbg(mod(ilonc,nlon)+1,ilatc,ilevc)     &
       +      fact_xint  * (1.0-fact_yint) *      fact_zint  * utbg(mod(ilonc,nlon)+1,ilatc,ilevc+1)   &
       +      fact_xint  *      fact_yint  * (1.0-fact_zint) * utbg(mod(ilonc,nlon)+1,ilatc+1,ilevc)   &
       +      fact_xint  *      fact_yint  *      fact_zint  * utbg(mod(ilonc,nlon)+1,ilatc+1,ilevc+1) 
   vdvt =  (1.0-fact_xint) * (1.0-fact_yint) * (1.0-fact_zint) * vtbg(ilonc,ilatc,ilevc)        &
       + (1.0-fact_xint) * (1.0-fact_yint) *      fact_zint  * vtbg(ilonc,ilatc,ilevc+1)      &
       + (1.0-fact_xint) *      fact_yint  * (1.0-fact_zint) * vtbg(ilonc,ilatc+1,ilevc)      &
       + (1.0-fact_xint) *      fact_yint  *      fact_zint  * vtbg(ilonc,ilatc+1,ilevc+1)    &
       +      fact_xint  * (1.0-fact_yint) * (1.0-fact_zint) * vtbg(mod(ilonc,nlon)+1,ilatc,ilevc)     &
       +      fact_xint  * (1.0-fact_yint) *      fact_zint  * vtbg(mod(ilonc,nlon)+1,ilatc,ilevc+1)   &
       +      fact_xint  *      fact_yint  * (1.0-fact_zint) * vtbg(mod(ilonc,nlon)+1,ilatc+1,ilevc)   &
       +      fact_xint  *      fact_yint  *      fact_zint  * vtbg(mod(ilonc,nlon)+1,ilatc+1,ilevc+1) 
   vdn2t = (1.0-fact_xint) * (1.0-fact_yint) * (1.0-fact_zint) * bn2t(ilonc,ilatc,ilevc)        &
       + (1.0-fact_xint) * (1.0-fact_yint) *      fact_zint  * bn2t(ilonc,ilatc,ilevc+1)      &
       + (1.0-fact_xint) *      fact_yint  * (1.0-fact_zint) * bn2t(ilonc,ilatc+1,ilevc)      &
       + (1.0-fact_xint) *      fact_yint  *      fact_zint  * bn2t(ilonc,ilatc+1,ilevc+1)    &
       +      fact_xint  * (1.0-fact_yint) * (1.0-fact_zint) * bn2t(mod(ilonc,nlon)+1,ilatc,ilevc)     &
       +      fact_xint  * (1.0-fact_yint) *      fact_zint  * bn2t(mod(ilonc,nlon)+1,ilatc,ilevc+1)   &
       +      fact_xint  *      fact_yint  * (1.0-fact_zint) * bn2t(mod(ilonc,nlon)+1,ilatc+1,ilevc)   &
       +      fact_xint  *      fact_yint  *      fact_zint  * bn2t(mod(ilonc,nlon)+1,ilatc+1,ilevc+1) 
   vdat =  (1.0-fact_xint) * (1.0-fact_yint) * (1.0-fact_zint) * alpt(ilonc,ilatc,ilevc)        &
       + (1.0-fact_xint) * (1.0-fact_yint) *      fact_zint  * alpt(ilonc,ilatc,ilevc+1)      &
       + (1.0-fact_xint) *      fact_yint  * (1.0-fact_zint) * alpt(ilonc,ilatc+1,ilevc)      &
       + (1.0-fact_xint) *      fact_yint  *      fact_zint  * alpt(ilonc,ilatc+1,ilevc+1)    &
       +      fact_xint  * (1.0-fact_yint) * (1.0-fact_zint) * alpt(mod(ilonc,nlon)+1,ilatc,ilevc)     &
       +      fact_xint  * (1.0-fact_yint) *      fact_zint  * alpt(mod(ilonc,nlon)+1,ilatc,ilevc+1)   &
       +      fact_xint  *      fact_yint  * (1.0-fact_zint) * alpt(mod(ilonc,nlon)+1,ilatc+1,ilevc)   &
       +      fact_xint  *      fact_yint  *      fact_zint  * alpt(mod(ilonc,nlon)+1,ilatc+1,ilevc+1) 

case('Spline','spline')
   call getBG_prep_tend(vx,vy,vz,time)
    call EZspline_interp(ut_spl,vx,vy,vz,vdut,ier)
    call EZspline_interp(vt_spl,vx,vy,vz,vdvt,ier)
    call EZspline_interp(n2t_spl,vx,vy,vz,vdn2t,ier)
    call EZspline_interp(alpt_spl,vx,vy,vz,vdat,ier)
   call getBG_free_tend
end select

if (ier.ne.0)then
 write(*,*)'ERROR in getBG_tend ::ier=',ier
 stop
end if


return
end subroutine getBG_tend
!---------------------------------------!
subroutine getBG_free

call EZspline_free(u_spl,ier)
call EZspline_free(v_spl,ier)
call EZspline_free(n2_spl,ier)
call EZspline_free(alp_spl,ier)

return
end subroutine getBG_free
!---------------------------------------!
subroutine getBG_free_tend

call EZspline_free(ut_spl,ier)
call EZspline_free(vt_spl,ier)
call EZspline_free(n2t_spl,ier)
call EZspline_free(alpt_spl,ier)

return
end subroutine getBG_free_tend
!-----------------------------------------------!
end module getBG
!===============================================!
