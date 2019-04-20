!===============================================!
module getBG
use EZspline_obj
use EZspline

! matrices for EZspline
type(EZspline3_r4)::u_spl,v_spl,n2_spl,alp_spl
integer,parameter::ibcs1(2)=(/-1, -1/) ! periodic
integer,parameter::ibcs2(2)=(/0 , 0 /)
integer,parameter::ibcs3(2)=(/0 , 0 /)


contains!---------------------------------------!

subroutine getBG_init
use setup

!!! For cyclic condition
real(4)::axlon_ext(nlon+1) 
real(4)::ubg_ext(nlon+1,nlat,nlev) 
real(4)::vbg_ext(nlon+1,nlat,nlev) 
real(4)::bn2_ext(nlon+1,nlat,nlev) 
real(4)::alp_ext(nlon+1,nlat,nlev) 



if (axlonl.eq.0.0.and.axlonr.eq.360.0)then!!! cyclic in x direction

 axlon_ext(1:nlon)=axlon
 axlon_ext(nlon+1)=360.0
 ubg_ext(1:nlon,:,:)=ubg
 vbg_ext(1:nlon,:,:)=vbg
 bn2_ext(1:nlon,:,:)=bn2
 alp_ext(1:nlon,:,:)=alp
 ubg_ext(nlon+1,:,:)=ubg(1,:,:)
 vbg_ext(nlon+1,:,:)=vbg(1,:,:)
 bn2_ext(nlon+1,:,:)=bn2(1,:,:)
 alp_ext(nlon+1,:,:)=alp(1,:,:)

 call EZspline_init(u_spl,nlon+1,nlat,nlev,ibcs1,ibcs2,ibcs3,ier)
 call EZspline_init(v_spl,nlon+1,nlat,nlev,ibcs1,ibcs2,ibcs3,ier)
 call EZspline_init(n2_spl,nlon+1,nlat,nlev,ibcs1,ibcs2,ibcs3,ier)
 call EZspline_init(alp_spl,nlon+1,nlat,nlev,ibcs1,ibcs2,ibcs3,ier)


 ! set grid location
 u_spl%x1=axlon_ext
 v_spl%x1=axlon_ext
 n2_spl%x1=axlon_ext
 alp_spl%x1=axlon_ext
 u_spl%x2=axlatSN
 v_spl%x2=axlatSN
 n2_spl%x2=axlatSN
 alp_spl%x2=axlatSN
 u_spl%x3=zlev
 v_spl%x3=zlev
 n2_spl%x3=zlev
 alp_spl%x3=zlev 

 call EZspline_setup(u_spl,   ubg_ext, ier)
 call EZspline_setup(v_spl,   vbg_ext, ier)
 call EZspline_setup(n2_spl,  bn2_ext, ier)
 call EZspline_setup(alp_spl, alp_ext, ier)

else !!! not cyclic

 call EZspline_init(u_spl,nlon,nlat,nlev,ibcs2,ibcs2,ibcs3,ier)
 call EZspline_init(v_spl,nlon,nlat,nlev,ibcs2,ibcs2,ibcs3,ier)
 call EZspline_init(n2_spl,nlon,nlat,nlev,ibcs2,ibcs2,ibcs3,ier)
 call EZspline_init(alp_spl,nlon,nlat,nlev,ibcs2,ibcs2,ibcs3,ier)


 ! set grid location
 u_spl%x1=axlon
 v_spl%x1=axlon
 n2_spl%x1=axlon
 alp_spl%x1=axlon
 u_spl%x2=axlatSN
 v_spl%x2=axlatSN
 n2_spl%x2=axlatSN
 alp_spl%x2=axlatSN
 u_spl%x3=zlev
 v_spl%x3=zlev
 n2_spl%x3=zlev
 alp_spl%x3=zlev

 call EZspline_setup(u_spl, ubg, ier)
 call EZspline_setup(v_spl, vbg, ier)
 call EZspline_setup(n2_spl, bn2, ier)
 call EZspline_setup(alp_spl, alp, ier)

end if
return
end subroutine getBG_init
!---------------------------------------!
subroutine getBG_val(vu,vv,vn2,va,vx,vy,vz)

real(4)::vu,vv,vn2,va
real(4)::vx,vy,vz

   call EZspline_interp(u_spl,vx,vy,vz,vu,ier)
   call EZspline_interp(v_spl,vx,vy,vz,vv,ier)
   call EZspline_interp(n2_spl,vx,vy,vz,vn2,ier)
   call EZspline_interp(alp_spl,vx,vy,vz,valp,ier)

if (ier.ne.0)then
 write(*,*)'ERROR in getBG_val ::ier=',ier
 stop
end if

return
end subroutine getBG_val
!---------------------------------------!
subroutine getBG_grad(vdu,vdv,vdn2,vda,idx,idy,idz,vx,vy,vz)
include 'Pcon.h'

real(4)::vu,vv,vn2,va
real(4)::vx,vy,vz
integer::idx,idy,idz

call EZspline_derivative(u_spl,idx,idy,idz,vx,vy,vz,vdu,ier)
call EZspline_derivative(v_spl,idx,idy,idz,vx,vy,vz,vdv,ier)
call EZspline_derivative(n2_spl,idx,idy,idz,vx,vy,vz,vdn2,ier)
call EZspline_derivative(alp_spl,idx,idy,idz,vx,vy,vz,vda,ier)

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

return
end subroutine getBG_grad
!---------------------------------------!
subroutine getBG_free

call EZspline_free(u_spl,ier)
call EZspline_free(v_spl,ier)
call EZspline_free(n2_spl,ier)
call EZspline_free(alp_spl,ier)

return
end subroutine getBG_free
!-----------------------------------------------!
end module getBG
!===============================================!
