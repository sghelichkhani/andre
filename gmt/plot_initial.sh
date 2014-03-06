#!/bin/bash

path1=/import/netapp-m-02-terra/horbach/adjoint/Endfeld/TOMO_NORESTRICT/plotdata
path2=/import/netapp-m-02-terra/horbach/adjoint/Endfeld/TOMO_SHIFT_NORESTRICT/plotdata
path3=/import/netapp-m-02-terra/horbach/adjoint/Endfeld/Blank0512_mt256/plotdata
path4=/import/netapp-m-02-terra/horbach/adjoint/runs/20.1_40Ma_BACKWARD_LRZ/plotdata
#path4=/import/netapp-m-02-terra/horbach/adjoint/runs/20.2_40Ma_BACKWARD_tethys_81out/gmt


size=12
res=200
ydist=-6.5c

lsize=18p

psize_tmp=19

psize=$psize_tmp"c"
dsize=`echo $psize_tmp - 3.5 | bc`"c"
ydist3=`echo $psize_tmp - 6 | bc`"c"   
ydist4=`echo $psize_tmp - 8.5 | bc`"c"

gmtset PS_MEDIA 29cx50c
gmtset FONT_ANNOT_PRIMARY 16p
gmtset FONT_LABEL 20p

#makecpt -Cpolar2.cpt -T-500/500/10 -Z -D > test.cpt
makecpt -Cpolar2.cpt -T-15/15/0.5 -Z -D > test.cpt

for i in {1..4}; do
if [ "$i" -eq "1" ]; then
xyz2grd $path1/T_old_0316 -Gtest01_$i.grd -I1 -Rd
xyz2grd $path1/T_old_1399 -Gtest02_$i.grd -I1 -Rd
xyz2grd $path1/T_old_2800 -Gtest03_$i.grd -I1 -Rd
else if [ "$i" -eq "2" ]; then
xyz2grd $path2/T_0320_00 -Gtest01_$i.grd -I1 -Rd
xyz2grd $path2/T_1400_00 -Gtest02_$i.grd -I1 -Rd
xyz2grd $path2/T_2800_00 -Gtest03_$i.grd -I1 -Rd
else if [ "$i" -eq "3" ]; then
xyz2grd $path3/T_0316_01 -Gtest01_$i.grd -I1 -Rd
xyz2grd $path3/T_1399_01 -Gtest02_$i.grd -I1 -Rd
xyz2grd $path3/T_2799_01 -Gtest03_$i.grd -I1 -Rd
else if [ "$i" -eq "4" ]; then
xyz2grd $path4/T_0270_01 -Gtest01_$i.grd -I1 -Rd
xyz2grd $path4/T_1399_01 -Gtest02_$i.grd -I1 -Rd
xyz2grd $path4/T_2799_01 -Gtest03_$i.grd -I1 -Rd
#xyz2grd $path4/t127.0360.80_00 -Gtest01_$i.grd -I1 -Rd
#xyz2grd $path4/t127.1400.80_00 -Gtest02_$i.grd -I1 -Rd
#xyz2grd $path4/t127.2800.80_00 -Gtest03_$i.grd -I1 -Rd
fi
fi
fi
fi;
done

# IMAGES
for i in {1..4}; do
if [ "$i" -eq "1" ]; then
grdimage test01_$i".grd" -Y$psize -X0.2c -Ctest.cpt -Rd -JW$size -K -E$res > test.ps
pscoast -JW$size -Rd -W0.9p -Dc -K -O -A0/0/1 >> test.ps
else if [ "$i" -eq "2" ]; then
grdimage test01_$i".grd" -X12.5c -Y$ydist3 -Ctest.cpt -Rd -JW$size -K -E$res -O >> test.ps
pscoast -JW$size -Rd -W0.9p -Dc -K -O -A0/0/1 >> test.ps
else
grdimage test01_$i".grd" -X18.75c -Y$ydist4 -Ctest.cpt -Rd -JW$size -K -E$res -O >> test.ps
pscoast -JW$size -Rd -W0.9p -Dc -K -O -A0/0/1 >> test.ps
fi
fi

for j in {2..3}; do
grdimage test0$j"_"$i".grd" -Y$ydist  -Ctest.cpt -Rd -JW$size -K -E$res -O >> test.ps
pscoast -JW$size -Rd -W0.9p -Dc -O -A0/0/1 -K >> test.ps;
done;

if [ "$i" -ge "2" ]; then
# LAYERS
pstext -Rd -JW$size -F+f$lsize -O -K -N -Y$dsize -X-6.25c >> test.ps << END
0 0 300 km
END
pstext -Rd -JW$size -F+f$lsize -O -K -N -Y$ydist >> test.ps << END
0 0 1450 km
END
pstext -Rd -JW$size -F+f$lsize -O -K -N -Y$ydist >> test.ps << END
0 0 2800 km
END
fi
done


# SCALE
psscale -D-6c/-4c/20c/1ch  -Ctest.cpt -O -S -B5:"temperature variations [%]":/:: -K -E >> test.ps

# HEADLINE
pstext -Rd -JW$size -F+f22p -O -K -N -Y$dsize -X-31.5c >> test.ps << END
0 0 a) tomography
END
pstext -Rd -JW$size -F+f22p -O -K -N -X12.5c >> test.ps << END
0 0 b) rotated tomography
END
pstext -Rd -JW$size -F+f22p -O -K -N -X12.5c >> test.ps << END
0 0 c) homogeneous mantle
END
pstext -Rd -JW$size -F+f22p -O -K -N -X12.5c >> test.ps << END
0 0 d) backward advection
END


ps2pdf test.ps
mv test.pdf initial.pdf
rm test.ps
evince initial.pdf

