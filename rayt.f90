!============================================!
subroutine raytracing
use setup
use getBG

do itime=1,ntmax-1

do iray=1,nray
! write(*,*)'ray #', iray,'/',nray

  if (rayi(1,iray).ge.1.and.rayi(itime,iray).ge.1)then

!      write(*,*)'time ', rayt(itime,iray) / 86400.0 ,' (day)'
    
    if (dr0.ne.0.0)then
     call RT_cg(rayx(itime,iray),rayy(itime,iray),rayz(itime,iray),rayt(itime,iray),rayk(itime,iray),rayl(itime,iray),raym(itime,iray),rayw(itime,iray),cgx,cgy,cgz)
     cgabs=sqrt(cgx**2+cgy**2+cgz**2)
     dt= dr0 / cgabs
    else
     dt = dt0
    end if


!!!! test
!!!   write(*,*) rayt(itime,iray),rayx(itime,iray),rayy(itime,iray),rayz(itime,iray),rayk(itime,iray),rayl(itime,iray),raym(itime,iray)

    call integrate (rayx(itime+1,iray),rayy(itime+1,iray),rayz(itime+1,iray), &
                    rayt(itime+1,iray), &
                    rayk(itime+1,iray),rayl(itime+1,iray),raym(itime+1,iray), &
                    rayw(itime+1,iray),rayh(itime+1,iray),rayb(itime+1,iray), &
                    rayi(itime+1,iray), &
                    rayx(itime,iray),rayy(itime,iray),rayz(itime,iray), &
                    rayt(itime,iray), &
                    rayk(itime,iray),rayl(itime,iray),raym(itime,iray), &
                    rayw(itime,iray),rayh(itime,iray),rayb(itime,iray), &
                    rayi(itime,iray),  dt)
!    if (rayi(itime+1,iray).le.0)  goto 9999
!   end do 
!9999 continue
  end if

end do

end do


!!! Test

!do iray=1,nray
! if (rayi(1,iray).ge.1)then
! write(*,*)'ray #', iray,'/',nray
!  do itime=1,ntime
!   if (rayi(itime,iray).ge.1)then
!    write(*,*) rayt(itime,iray)/86400.0, rayz(itime,iray)/1000.0 , rayw(itime,iray)*1000.0 , rayh(itime,iray)*1000.0 
!   end if
!  end do
! end if
!end do



return
end subroutine raytracing

!============================================!

subroutine integrate (vx1,vy1,vz1,time1,vk1,vl1,vm1,vw1,vh1,vb1,istat1, &
                      vx0,vy0,vz0,time0,vk0,vl0,vm0,vw0,vh0,vb0,istat0, dt)
use setup
use getBG

real(4)::xyzklmw0(7),xyzklmw1(7)
real(4)::xyzklmw_1(7),xyzklmw_2(7)
real(4)::dxyzklmw_1(7),dxyzklmw_2(7),dxyzklmw_3(7)

real(4)::rwork(7,3)
real(4)::rwork1(7,1)

external RT_full

if (istat0.le.0) goto 999

xyzklmw0=(/vx0,vy0,vz0,vk0,vl0,vm0,vw0/)

!write(*,*)xyzklmw0

dth=0.5*dt

call RT_full(7,time,xyzklmw0,dxyzklmw_1)

select case(trim(cintgmode))
 case('RK4')
 call odrk4(7,RT_full,time0,dt,xyzklmw0,dxyzklmw_1,xyzklmw1,rwork)
 case('Euler')
 call odrk1(7,RT_full,time0,dt,xyzklmw0,dxyzklmw_1,xyzklmw1,rwork1)
 case ('RK2')
 call odrk2(7,RT_full,time0,dt,xyzklmw0,dxyzklmw_1,xyzklmw1,rwork1)
end select


vx1=xyzklmw1(1); vy1=xyzklmw1(2); vz1=xyzklmw1(3)
vk1=xyzklmw1(4); vl1=xyzklmw1(5); vm1=xyzklmw1(6)
vw1=xyzklmw1(7)

time1 = time0 + dt

call range(vx1,vy1,vz1,istat1)
if (istat1.le.0)goto 999

!!! diagnose omega(hat)

!!! vw1=vw0

call getBG_val(vu,vv,vn2,valp,vx1,vy1,vz1,time1)

vh1=vw1-vk1*vu-vl1*vv

if (vh1**2.ge.vn2.or.vh1**2.le.vf(vy1)**2) goto 999
if (vh1*vh0.le.0.0) goto 999 !!! cross a critical level

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! correct m by dispersion relation (optional)
vm2 = (vk1**2+vl1**2)*(vn2-real(inhd)*vh1**2)/(vh1**2-vf(vy1)**2) - valp**2
if (vm2.le.0.0) goto 999
vm1_diag=sign(sqrt(vm2),real(iupwd)*(-1)*vh1)!!! dependent of Cgz sign

!write(*,*) vm1,vm1_diag
vm1=vm1_diag

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

istat1=1

if (iamp.gt.0)then !!! consider amplitude 
!!! call saturation
else
 vb1=vb0
end if

return

999 continue !!
istat1=-1
vx1=rmiss
vy1=rmiss
vz1=rmiss
vk1=rmiss
vl1=rmiss
vm1=rmiss
vw1=rmiss
vh1=rmiss
vb1=rmiss

return
end subroutine integrate

!============================================!
subroutine RT_full(n,time,xyzklmw,dxyzklmw)
use setup
use getBG

integer,intent(in)::n
real(4),intent(in)::time

real(4),intent(in)::xyzklmw(n)
real(4),intent(out)::dxyzklmw(n)

vx=xyzklmw(1); vy=xyzklmw(2); vz=xyzklmw(3)
vk=xyzklmw(4); vl=xyzklmw(5); vm=xyzklmw(6)
vw=xyzklmw(7)

istat=1
call range(vx,vy,vz,istat)
if (istat.le.0) then
 write(*,*) 'out.'
 dxyzklmw=0.0
 return
end if

call getBG_val(vu,vv,vn2,valp,vx,vy,vz,time)
call getBG_tend(dut,dvt,dvn2t,dvalpt,vx,vy,vz,time)
call getBG_grad(dux,dvx,dvn2x,dvalpx,1,0,0,vx,vy,vz,time)
call getBG_grad(duy,dvy,dvn2y,dvalpy,0,1,0,vx,vy,vz,time)
call getBG_grad(duz,dvz,dvn2z,dvalpz,0,0,1,vx,vy,vz,time)


vh=vw-vk*vu-vl*vv

x1=vk**2+vl**2
x2=vn2-real(inhd)*vh**2
x3=vh**2-vf(vy)**2
vodl=(real(inhd)*x1+vm**2+valp**2)*vh


dxt=vu+vk*x2/vodl
dyt=vv+vl*x2/vodl
dzt=-vm*x3/vodl
dkt=-vk*dux-vl*dvx-(x1*dvn2x-x3*dvalpx*2.0*valp)*0.5/vodl
dlt=-vk*duy-vl*dvy-(x1*dvn2y-x3*dvalpy*2.0*valp)*0.5/vodl-vf(vy)*(vm**2+valp**2)*beta(vy)/vodl
dmt=-vk*duz-vl*dvz-(x1*dvn2z-x3*dvalpz*2.0*valp)*0.5/vodl

!!! frequency change by (slowly) tenporally-varying background
dwt = vk * dut + vl * dvt + (x1*dvn2t-x3*dvalpt*2.0*valp)*0.5/vodl

if (icur.eq.1) then !!! curvature terms 
 dkt=dkt+(vk*tan(vy*drad)*vv+(tan(vy*drad)*vk*vl*x2-vk*vm*x3)/vodl)/er
 dlt=dlt+(-vk*tan(vy*drad)*vu+(-tan(vy*drad)*vk*vk*x2-vl*vm*x3)/vodl)/er
 dmt=dmt+(vk*vu+vl*vv+((vk**2+vl**2)*x2)/vodl)/er
end if

!!! m/s ==> deg/s
yfac=1.0/er/drad
zfac=1.0
xfac=yfac/cos(vy*drad)

dxt = dxt * xfac
dyt = dyt * yfac
dzt = dzt * zfac

dxyzklmw(1:7)=(/dxt,dyt,dzt,dkt,dlt,dmt,dwt/)

return
end subroutine RT_full
!============================================!
subroutine RT_cg(vx,vy,vz,time,vk,vl,vm,vw,dxt,dyt,dzt)
use setup
use getBG


call getBG_val(vu,vv,vn2,valp,vx,vy,vz,time)
call getBG_grad(dux,dvx,dvn2x,dvalpx,1,0,0,vx,vy,vz,time)
call getBG_grad(duy,dvy,dvn2y,dvalpy,0,1,0,vx,vy,vz,time)
call getBG_grad(duz,dvz,dvn2z,dvalpxz,0,0,1,vx,vy,vz,time)

vh=vw-vk*vu-vl*vv

x1=vk**2+vl**2
x2=vn2-real(inhd)*vh**2
x3=vh**2-vf(vy)**2
vodl=(real(inhd)*x1+vm**2+valp**2)*vh


dxt=vu+vk*x2/vodl
dyt=vv+vl*x2/vodl
dzt=-vm*x3/vodl


return
end subroutine RT_cg
!============================================!

subroutine range (vx,vy,vz,istat)
use setup

if (axlonl.eq.0.0.and.axlonr.eq.360.0)then !!! cyclic
 if (vx.lt.0.0) vx=vx+360.0 
 if (vx.ge.360.0) vx=vx-360.0 
else
 if (vx.lt.axlonl) istat=-1
 if (vx.gt.axlonr) istat=-1
end if

if (vy.ge.axlatr.or.vy.le.axlatl) istat=-1
if (vz.ge.ztop.or.vz.le.zbot) istat=-1

return
end subroutine range

!============================================!
