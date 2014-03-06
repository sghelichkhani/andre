#!/bin/bash

path=/import/netapp-m-02-terra/horbach/adjoint/runs/18.1_20Ma_TETHYS_mt512/Iteration_00/plotdata
#path2=/import/netapp-m-02-terra/horbach/adjoint/runs/18.2_20Ma_mt256/Iteration_05/plotdata
path2=/import/netapp-m-02-terra/horbach/adjoint/runs/plotdata

# Fig 1
#path=/import/netapp-m-02-terra/horbach/adjoint/runs/12.2_40Ma_Africa_SCALED_b-files/Iteration_01/plotdata
path=../../TERRA512/
#path2=/import/netapp-m-02-terra/horbach/adjoint/runs/12.2_40Ma_Africa_SCALED_b-files/Iteration_07/plotdata
path2=../../TERRA512/

# Fig 2
#path=/import/netapp-m-02-terra/horbach/adjoint/runs/16.2_40Ma_TOMO512_LRZ/Iteration_07/plotdata
#path2=/import/netapp-m-02-terra/horbach/adjoint/Endfeld/TOMO_NORESTRICT/plotdata

size=12
res=200
ilay=4

ydist_tmp=-6.5
psize_tmp=57

ydist=$ydist_tmp"c"
dist_help=`echo -1*$ydist_tmp*$ilay + 9*$ydist_tmp | bc`
psize=`echo $psize_tmp + $dist_help | bc`"c"
dsize=`echo $psize_tmp + $dist_help - 2.5 | bc`"c"
ydist3=`echo $psize_tmp + $dist_help -5 | bc`"c"   

#layer[1]=0090
#layer[2]=0950
#layer[3]=1760
#layer[4]=2800

layer[1]=0180
layer[2]=1083
layer[3]=1806
layer[4]=2890

#layer[1]=0090
#layer[2]=0135
#layer[3]=0225
#layer[4]=0316
#layer[5]=0406
#layer[6]=0541
#layer[7]=0632
#layer[8]=2483
#layer[9]=2799

plot1=h602
time1='' #_00
plot2=d602	#T_old
time2='' #_00

gmtset PS_MEDIA `echo 67 + $dist_help | bc`"c"x25c
gmtset FONT_ANNOT_PRIMARY 16p
gmtset FONT_LABEL 20p

makecpt -Cpolar2.cpt -T-500/500/10 -Z -D > test.cpt
makecpt -Cpolar2.cpt -T300/4200/50 -Z -D > test.cpt
makecpt -Cpolar2.cpt -T1e-10/1e-09/0.1e-10 -Z -D > test1.cpt
makecpt -Cpolar2.cpt -T1e-10/1e-09/0.1e-10 -Z -D > test2.cpt
#makecpt -Cpolar2.cpt -T-20/20/0.5 -Z -D > test.cpt
makecpt -Cpolar2.cpt -T0/5/0.2 -Z -D > vel.cpt
makecpt -Cpolar2.cpt -T-2000/2000/200 -Z -D > shrh.cpt


# IMAGES
for ((i=1; i <= $ilay; i++)); do

xyz2grd $path/$plot1''.${layer[$i]}''$time1 -Gtest1_$i.grd -I+4 -Rd
#surface $path/$plot1''.${layer[$i]}''$time1 -Gtest1_$i.grd -I0.5 -Rd
xyz2grd $path/$plot2''.${layer[$i]}''$time2 -Gtest2_$i.grd -I+4 -Rd
#surface $path2/$plot2''.${layer[$i]}''$time2 -Gtest2_$i.grd -I0.5 -Rd

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
psscale -D-6c/-1c/8c/1ch  -Ctest1.cpt -O -S -B10:"temperature variations [%]":/:: -K -E >> test.ps
psscale -D6c/-1c/8c/1ch  -Ctest2.cpt -O -S -B500:"shear heating":/:: -K -E >> test.ps
#psscale -D0c/-1c/14c/1ch  -Ctest.cpt -O -S -B250:"adjoint temperature [K]":/:: -K -E >> test.ps


lsize=18p
# LAYERS
for((i=1;i<=$ilay;i++)); do
if [ "$i" -eq "1" ]; then
pstext -Rd -JW$size -F+f$lsize -O -K -N -Y$dsize -X-6.25c >> test.ps << END
0 0 `echo ${layer[$i]} + 0 | bc` km
END
else
pstext -Rd -JW$size -F+f$lsize -O -K -N -Y$ydist >> test.ps << END
0 0 `echo ${layer[$i]} + 0 | bc` km
END
fi;
done


# HEADLINE
pstext -Rd -JW$size -F+f22p -O -K -N -Y$dsize -X-6c >> test.ps << END
0 0 Iteration 1
END
pstext -Rd -JW$size -F+f22p -O -K -N -X12.5c >> test.ps << END
0 0 Iteration 7
END


ps2pdf test.ps
mv test.pdf slices.pdf
rm test.ps
evince slices.pdf

