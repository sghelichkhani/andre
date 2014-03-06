#!/bin/bash

path=./backw/

size=25
res=200
ilay=1

ydist_tmp=-6.5
psize_tmp=57

ydist=$ydist_tmp"c"
dist_help=`echo -1*$ydist_tmp*$ilay + 9*$ydist_tmp | bc`
psize=`echo $psize_tmp + $dist_help | bc`"c"
dsize=`echo $psize_tmp + $dist_help - 2.5 | bc`"c"
ydist3=`echo $psize_tmp + $dist_help -5 | bc`"c"   

layer[1]=0180

gmtset PS_MEDIA 18cx25.5c
gmtset FONT_ANNOT_PRIMARY 16p
gmtset FONT_LABEL 20p

range11=-10
range12=10
range21=-500
range22=500
#range21=-2
#range22=2
d1=`echo $range12/100 - $range11/100 | bc`
d2=`echo $range22/100 - $range21/100 | bc`

makecpt -Cpolar2.cpt -T$range11/$range12/$d1 -Z -D > test1.cpt
makecpt -Cpolar2.cpt -T$range21/$range22/$d2 -Z -D > test2.cpt

for ((j=0; j<=0;j++)); do

plot1=t127

# IMAGES
for ((i=1; i <= $ilay; i++)); do

if [ "$j" -le "9" ]; then
xyz2grd $path/$plot1''.${layer[$i]}'.0'$j'_00' -Gtest1_$i.grd -I+5 -Rd
else
xyz2grd $path/$plot1''.${layer[$i]}'.'$j'_00' -Gtest1_$i.grd -I+5 -Rd
fi;

#surface $path/$plot1''.${layer[$i]}''$time1 -Gtest1_$i.grd -I0.5 -Rd
#xyz2grd $path2/$plot2''.${layer[$i]}''$time2 -Gtest2_$i.grd -I+5 -Rd
#surface $path2/$plot2''.${layer[$i]}''$time2 -Gtest2_$i.grd -I0.5 -Rd

#grdmath test1_$i.grd ABS = test1_$i.grd
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


# SCALE
psscale -D13c/-1c/14c/1ch  -Ctest1.cpt -O -S -B`echo $d1/20 | bc`:"temperature variations [%]":/:: -K -E >> test.ps
#psscale -D-6c/-1c/8c/1ch  -Ctest1.cpt -O -S -B2e-8:"temperature variations [%]":/:: -K -E >> test.ps

#psscale -D-6c/-1c/8c/1ch  -Ctest1.cpt -O -S -B`echo $d1/20 | bc`:"velocity [cm/yr]":/:: -K -E >> test.ps

#psscale -D-6c/-1c/8c/1ch  -Ctest1.cpt -O -S -B5e-10:"temp. gradient.":/:: -K -E >> test.ps
#psscale -D6c/-1c/8c/1ch  -Ctest2.cpt -O -S -B5e-10:"adj. buoyancy":/:: -K -E >> test.ps
#psscale -D0c/-1c/14c/1ch  -Ctest.cpt -O -S -B250:"adjoint temperature [K]":/:: -K -E >> test.ps

pstext -Rd -JW$size -F+f22p -O -K -N -X11c -Y5.5c >> test.ps << END
0 0 `echo $j/2 | bc` Ma
END

ps2pdf test.ps
mv test.pdf backw/slices.pdf
rm test.ps
# evince slices.pdf
convert backw/slices.pdf backw/slices$j.png;

done
