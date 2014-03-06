#!/bin/bash

path=/import/netapp-terra/horbach/adjoint/runs/20.1_BACKWARD_LRZ/plotdata
path2=/import/netapp-terra/horbach/adjoint/Endfeld/TOMO_NORESTRICT/plotdata

size=12
res=100
ydist=-6.5c

psize=57c
dsize=54.5c	    # psize-2.5
ydist3=52c      # psize-5

gmtset PS_MEDIA 67cx25c
gmtset FONT_ANNOT_PRIMARY 16p
gmtset FONT_LABEL 20p

#makecpt -Cpolar2.cpt -T-500/500/10 -Z -D > test.cpt
makecpt -Cpolar2.cpt -T-20/20/0.5 -Z -D > test.cpt

xyz2grd $path/T_new_0090 -Gtest1_1.grd -I1 -Rd
xyz2grd $path/T_new_0338 -Gtest1_2.grd -I1 -Rd
xyz2grd $path/T_new_0699 -Gtest1_3.grd -I1 -Rd
xyz2grd $path/T_new_1038 -Gtest1_4.grd -I1 -Rd
xyz2grd $path/T_new_1399 -Gtest1_5.grd -I1 -Rd
xyz2grd $path/T_new_1761 -Gtest1_6.grd -I1 -Rd
xyz2grd $path/T_new_2099 -Gtest1_7.grd -I1 -Rd
xyz2grd $path/T_new_2461 -Gtest1_8.grd -I1 -Rd
xyz2grd $path/T_new_2799 -Gtest1_9.grd -I1 -Rd
xyz2grd $path2/T_old_0090 -Gtest2_1.grd -I1 -Rd
xyz2grd $path2/T_old_0338 -Gtest2_2.grd -I1 -Rd
xyz2grd $path2/T_old_0699 -Gtest2_3.grd -I1 -Rd
xyz2grd $path2/T_old_1038 -Gtest2_4.grd -I1 -Rd
xyz2grd $path2/T_old_1399 -Gtest2_5.grd -I1 -Rd
xyz2grd $path2/T_old_1761 -Gtest2_6.grd -I1 -Rd
xyz2grd $path2/T_old_2099 -Gtest2_7.grd -I1 -Rd
xyz2grd $path2/T_old_2461 -Gtest2_8.grd -I1 -Rd
xyz2grd $path2/T_old_2799 -Gtest2_9.grd -I1 -Rd

#-180/+180/-89/+89
#grd2xyz test.grd -R-180/+180/-89/+89 > test.xyz

# IMAGES
grdimage test1_1.grd -Y$psize -X0.2c -Ctest.cpt -Rd -JW$size -K -E$res > test.ps
pscoast -JW$size -Rd -W0.9p -Dc -K -O -A0/0/1 >> test.ps

grdimage test1_2.grd -Y$ydist  -Ctest.cpt -Rd -JW$size -K -E$res -O >> test.ps
pscoast -JW$size -Rd -W0.9p -Dc -O -A0/0/1 -K >> test.ps

grdimage test1_3.grd -Y$ydist -Ctest.cpt -Rd -JW$size -K -E$res -O >> test.ps
pscoast -JW$size -Rd -W0.9p -Dc -O -A0/0/1 -K >> test.ps

grdimage test1_4.grd -Y$ydist -Ctest.cpt -Rd -JW$size -K -E$res -O >> test.ps
pscoast -JW$size -Rd -W0.9p -Dc -O -A0/0/1 -K >> test.ps

grdimage test1_5.grd -Y$ydist -Ctest.cpt -Rd -JW$size -K -E$res -O >> test.ps
pscoast -JW$size -Rd -W0.9p -Dc -O -A0/0/1 -K >> test.ps

grdimage test1_6.grd -Y$ydist  -Ctest.cpt -Rd -JW$size -K -E$res -O >> test.ps
pscoast -JW$size -Rd -W0.9p -Dc -O -A0/0/1 -K >> test.ps

grdimage test1_7.grd -Y$ydist -Ctest.cpt -Rd -JW$size -K -E$res -O >> test.ps
pscoast -JW$size -Rd -W0.9p -Dc -O -A0/0/1 -K >> test.ps

grdimage test1_8.grd -Y$ydist -Ctest.cpt -Rd -JW$size -K -E$res -O >> test.ps
pscoast -JW$size -Rd -W0.9p -Dc -O -A0/0/1 -K >> test.ps

grdimage test1_9.grd -Y$ydist -Ctest.cpt -Rd -JW$size -K -E$res -O >> test.ps
pscoast -JW$size -Rd -W0.9p -Dc -O -A0/0/1 -K >> test.ps

grdimage test2_1.grd -X12.5c -Y$ydist3 -Ctest.cpt -Rd -JW$size -K -E$res -O >> test.ps
pscoast -JW$size -Rd -W0.9p -Dc -K -O -A0/0/1 >> test.ps

grdimage test2_2.grd -Y$ydist -Ctest.cpt -Rd -JW$size -K -E$res -O >> test.ps
pscoast -JW$size -Rd -W0.9p -Dc -O -A0/0/1 -K >> test.ps

grdimage test2_3.grd -Y$ydist -Ctest.cpt -Rd -JW$size -K -E$res -O >> test.ps
pscoast -JW$size -Rd -W0.9p -Dc -O -A0/0/1 -K >> test.ps

grdimage test2_4.grd -Y$ydist -Ctest.cpt -Rd -JW$size -K -E$res -O >> test.ps
pscoast -JW$size -Rd -W0.9p -Dc -O -A0/0/1 -K >> test.ps

grdimage test2_5.grd -Y$ydist -Ctest.cpt -Rd -JW$size -K -E$res -O >> test.ps
pscoast -JW$size -Rd -W0.9p -Dc -O -A0/0/1 -K >> test.ps

grdimage test2_6.grd -Y$ydist -Ctest.cpt -Rd -JW$size -K -E$res -O >> test.ps
pscoast -JW$size -Rd -W0.9p -Dc -O -A0/0/1 -K >> test.ps

grdimage test2_7.grd -Y$ydist -Ctest.cpt -Rd -JW$size -K -E$res -O >> test.ps
pscoast -JW$size -Rd -W0.9p -Dc -O -A0/0/1 -K >> test.ps

grdimage test2_8.grd -Y$ydist -Ctest.cpt -Rd -JW$size -K -E$res -O >> test.ps
pscoast -JW$size -Rd -W0.9p -Dc -O -A0/0/1 -K >> test.ps

grdimage test2_9.grd -Y$ydist -Ctest.cpt -Rd -JW$size -K -E$res -O >> test.ps
pscoast -JW$size -Rd -W0.9p -Dc -O -A0/0/1 -K >> test.ps

# SCALE
psscale -D0c/-1c/14c/1ch  -Ctest.cpt -O -S -B5:"temperature variations [%]":/:: -K -E >> test.ps

lsize=18p
# LAYERS
pstext -Rd -JW$size -F+f$lsize -O -K -N -Y$dsize -X-6.25c >> test.ps << END
0 0 90 km
END
pstext -Rd -JW$size -F+f$lsize -O -K -N -Y$ydist >> test.ps << END
0 0 350 km
END
pstext -Rd -JW$size -F+f$lsize -O -K -N -Y$ydist >> test.ps << END
0 0 700 km
END
pstext -Rd -JW$size -F+f$lsize -O -K -N -Y$ydist >> test.ps << END
0 0 1050 km
END
pstext -Rd -JW$size -F+f$lsize -O -K -N -Y$ydist >> test.ps << END
0 0 1400 km
END
pstext -Rd -JW$size -F+f$lsize -O -K -N -Y$ydist >> test.ps << END
0 0 1750 km
END
pstext -Rd -JW$size -F+f$lsize -O -K -N -Y$ydist >> test.ps << END
0 0 2100 km
END
pstext -Rd -JW$size -F+f$lsize -O -K -N -Y$ydist >> test.ps << END
0 0 2450 km
END
pstext -Rd -JW$size -F+f$lsize -O -K -N -Y$ydist >> test.ps << END
0 0 2800 km
END


# HEADLINE
pstext -Rd -JW$size -F+f22p -O -K -N -Y$dsize -X-6c >> test.ps << END
0 0 Model (40 Ma)
END
pstext -Rd -JW$size -F+f22p -O -K -N -X12.5c >> test.ps << END
0 0 Model (0 Ma)
END


ps2pdf test.ps
mv test.pdf slices.pdf
rm test.ps
evince slices.pdf

