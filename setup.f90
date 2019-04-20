!============================================!

module setup

implicit real(a-h,o-z)

include 'Grid.h'
include 'Pcon.h'
include 'RayT.h'

!! Time-dependent BG field variables

real(4)::Ubg_xyzt(nlon,nlat,nlev,ntime) !!! zonal  wind
real(4)::Vbg_xyzt(nlon,nlat,nlev,ntime) !!! merid. wind
real(4)::Tbg_xyzt(nlon,nlat,nlev,ntime) !!! temperature 

!! BG field variables

real(4)::ubg(nlon,nlat,nlev) !!! zonal  wind
real(4)::vbg(nlon,nlat,nlev) !!! merid. wind
real(4)::tbg(nlon,nlat,nlev) !!! temperature 

!! diagnosed variables

real(4)::bn2(nlon,nlat,nlev) !!! buoyancy freq.
real(4)::alp(nlon,nlat,nlev) !!! 1/(2H) -- H: density scale height
real(4)::alp0 !!! constant alpha

!! time derivative

real(4)::Utbg(nlon,nlat,nlev) !!! zonal  wind tendency
real(4)::Vtbg(nlon,nlat,nlev) !!! merid. wind tendency
real(4)::Ttbg(nlon,nlat,nlev) !!! temperature tendency
real(4)::bn2t(nlon,nlat,nlev) !!! buoyancy freq. tendency
real(4)::alpt(nlon,nlat,nlev) !!! alpha tendency

!! ray tracers

integer::nray !! determined by file length

real(4)::rayx(ntmax,nraymax) !!! x location (deg)
real(4)::rayy(ntmax,nraymax) !!! y (deg)
real(4)::rayz(ntmax,nraymax) !!! z (m)
real(4)::rayk(ntmax,nraymax) !!! x wavenumber (rad/m)
real(4)::rayl(ntmax,nraymax) !!! y (rad/m)
real(4)::raym(ntmax,nraymax) !!! z (rad/m)
real(4)::rayw(ntmax,nraymax) !!! omega (rad/s)
real(4)::rayh(ntmax,nraymax) !!! omega(Hat)
real(4)::rayb(ntmax,nraymax) !!! amplitude |CgA|
real(4)::rayt(ntmax,nraymax) !!! time (s)
integer::rayi(ntmax,nraymax) !!! state of the wave 
                             !!! (1 : propagate, -1 : vanished, imiss: undefined)

character*15::cbgmode

character*5::cdiag

character*5::cintgmode
character*15::cintpmode

namelist /init_nml/ cbgmode,iamp,iupwd,cdiag

namelist /RayT_nml/ cintgmode,cintpmode, ilonsp,ilatsp,ilevsp,dr0,dt0,icur,irot,inhd,ishe

  contains 
   real function vf(vlat_deg)
    real(4)::vlat_deg !!! latitude
    vf=2.0*eomg*sin(drad*vlat_deg) * irot
   end function vf
   real function beta(vlat_deg)
    real(4)::vlat_deg !!! latitude
    beta=2.0*eomg/er*cos(drad*vlat_deg) * irot  
   end function beta

end module setup
!============================================!
