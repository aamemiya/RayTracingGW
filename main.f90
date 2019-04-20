
!============================================!

program main
use setup


open (11,file='calc.conf')
 read(11,nml=init_nml)
 read(11,nml=RayT_nml)
close(11)

call read_BG

call read_init

call RayTracing

call output

call QuickLook_zm
call QuickLook_pl

stop
end program main
!============================================!
