#!/bin/bash

path=/import/netapp-m-02-terra/horbach/adjoint/runs/18.4_40Ma_mt512_10pres_1021_F121/gmt_It2
path2=/import/netapp-m-02-terra/horbach/adjoint/runs/16.1_40Ma_BLANKMANTLE_LRZ/Iteration_07/plotdata

size=12
res=200
ilay=2

ydist_tmp=-8
psize_tmp=57

ydist=$ydist_tmp"c"
dist_help=`echo -1*$ydist_tmp*$ilay + 9*$ydist_tmp | bc`
psize=`echo $psize_tmp + $dist_help | bc`"c"
dsize=`echo $psize_tmp + $dist_help - 2.5 | bc`"c"
ydist3=`echo $psize_tmp + $dist_help -5 | bc`"c"   

ydist=-9c
ydist3=9c
psize=13.5c
dsize=13.5c

layer[1]=2030
layer[2]=2080
layer[3]=2032
layer[4]=2077


plot1=t363.
time1=.00_03
plot2=T_
time2=_40

gmtset PS_MEDIA 22cx25c
gmtset FONT_ANNOT_PRIMARY 16p
gmtset FONT_LABEL 20p

range11=-15
range12=15
range21=-15
range22=15
#range21=-2
#range22=2
d1=`echo $range12/100 - $range11/100 | bc`
d2=`echo $range22/100 - $range21/100 | bc`

makecpt -Cpolar2.cpt -T$range11/$range12/$d1 -Z -D > test1.cpt
makecpt -Cpolar2.cpt -T$range21/$range22/$d2 -Z -D > test2.cpt

# IMAGES
for ((i=1; i <= $ilay; i++)); do

j=`echo $i + 2 | bc`
xyz2grd $path/$plot1''${layer[$i]}''$time1 -Gtest1_$i.grd -I+5 -Rd
#surface $path/$plot1''${layer[$i]}''$time1 -Gtest1_$i.grd -I0.5 -Rd
xyz2grd $path2/$plot2''${layer[$j]}''$time2 -Gtest2_$i.grd -I+5 -Rd
#surface $path2/$plot2''${layer[$i]}''$time2 -Gtest2_$i.grd -I0.5 -Rd

#grdmath test1_$i.grd DIV = test1_$i.grd
#grdmath test2_$i.grd ABS = test2_$i.grd
#grdmath test1_$i.grd test2_$i.grd SUB = test1_$i.grd

if [ "$i" -eq "1" ]; then
grdimage test1_$i.grd -Y$psize -X0.2c -Ctest1.cpt -Rd -JW$size -K -E$res > test.ps
pscoast -JW$size -Rd -W0.9p -Dc -K -O -A0/0/1 >> test.ps
else
grdimage test1_$i.grd -Y$ydist  -Ctest1.cpt -Rd -JW$size -K -E$res -O >> test.ps
pscoast -JW$size -Rd -W0.9p -Dc -O -A0/0/1 -K >> test.ps
fi;
done

for((i=1;i<=$ilay;i++)); do
if [ "$i" -eq "1" ]; then
grdimage test2_$i.grd -X12.5c -Y$ydist3 -Ctest2.cpt -Rd -JW$size -K -E$res -O >> test.ps
pscoast -JW$size -Rd -W0.9p -Dc -K -O -A0/0/1 >> test.ps
else
grdimage test2_$i.grd -Y$ydist -Ctest2.cpt -Rd -JW$size -K -E$res -O >> test.ps
pscoast -JW$size -Rd -W0.9p -Dc -O -A0/0/1 -K >> test.ps
fi;
done


# SCALE
#psscale -D0c/-1c/14c/1ch  -Ctest.cpt -O -S -B5:"temperature variations [%]":/:: -K -E >> test.ps
#psscale -D-6c/-1c/8c/1ch  -Ctest1.cpt -O -S -B2e-8:"temperature variations [%]":/:: -K -E >> test.ps

#psscale -D-6c/-1c/8c/1ch  -Ctest1.cpt -O -S -B`echo $d1/20 | bc`:"velocity [cm/yr]":/:: -K -E >> test.ps
#psscale -D6c/-1c/8c/1ch  -Ctest2.cpt -O -S -B`echo $d2/20 | bc`:"velocity [cm/yr]":/:: -K -E >> test.ps

#psscale -D-6c/-1c/8c/1ch  -Ctest1.cpt -O -S -B`echo $d1/20 | bc`:"density var. [%]":/:: -K -E >> test.ps
#psscale -D6c/-1c/8c/1ch  -Ctest2.cpt -O -S -B`echo $d2/20 | bc`:"temperature var. [%]":/:: -K -E >> test.ps

#psscale -D-6c/-1c/8c/1ch  -Ctest1.cpt -O -S -B5e-10:"temp. gradient.":/:: -K -E >> test.ps
#psscale -D6c/-1c/8c/1ch  -Ctest2.cpt -O -S -B5e-10:"adj. buoyancy":/:: -K -E >> test.ps
#psscale -D0c/-1c/14c/1ch  -Ctest.cpt -O -S -B250:"adjoint temperature [K]":/:: -K -E >> test.ps

#psscale -D-1c/-1c/16c/1ch  -Ctest1.cpt -O -S -B`echo $d1/20 | bc`:"temperature anomalies [%]":/:: -K -E >> test.ps
psscale -D0c/-1c/14c/1ch  -Ctest1.cpt -O -S -B5:"temperature anomalies [%]":/:: -K -E >> test.ps

lsize=18p

# HEADLINE
pstext -Rd -JW$size -F+f22p -O -K -N -Y$dsize -X-12.5c >> test.ps << END
0 0 a) tomography
END
pstext -Rd -JW$size -F+f22p -O -K -N -X12.5c >> test.ps << END
0 0 b) rotated tomography
END
pstext -Rd -JW$size -F+f22p -O -K -N -Y$ydist -X-12.5c >> test.ps << END
0 0 c) homogeneous mantle
END
pstext -Rd -JW$size -F+f22p -O -K -N -X12.5c >> test.ps << END
0 0 d) backward advection
END


ps2pdf test.ps
mv test.pdf slices.pdf
rm test.ps
evince slices.pdf &
#mv slices.pdf ../../TERRA_Benchmark/src/code/

