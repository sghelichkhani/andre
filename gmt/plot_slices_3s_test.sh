#!/bin/bash

path=/import/netapp-m-02-terra/horbach/adjoint/runs/18.4_40Ma_mt512_10pres_1021_F121/gmt_It2
path2=$path
path3=/import/netapp-m-02-terra/horbach/adjoint/runs/16.1_40Ma_BLANKMANTLE_LRZ/Iteration_07/plotdata
path2=$path3

#path3=/import/netapp-m-02-terra/horbach/adjoint/Endfeld/Endfeld0512_mt256/plotdata
#path2=$path
#path3=$path

# Fig 1
#path=/import/netapp-m-02-terra/horbach/adjoint/runs/12.2_40Ma_Africa_SCALED_b-files/Iteration_01/plotdata
#path=../../TERRA_Benchmark/802/gmt
#path=./
#path=/import/netapp-m-02-terra/horbach/adjoint/Endfeld/Endfeld0064_mt128/plotdata
#path=../../TERRA512/357/gmt
#path=/import/buserror-data/horbach/TERRA_BM/802/gmt
#path2=/import/netapp-m-02-terra/horbach/adjoint/runs/12.2_40Ma_Africa_SCALED_b-files/Iteration_07/plotdata
#path2=../../TERRA_Benchmark/802/gmt
#path3=../../TERRA_Benchmark/802/gmt
#path2=./
#path2=/import/buserror-data/horbach/TERRA_BM/802/gmt

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

layer[1]=0630
layer[2]=1450
layer[3]=2030
layer[4]=2800

#layer[1]=0000
#layer[2]=0380
#layer[3]=0750
#layer[4]=1000

#layer[1]=0090
#layer[2]=0135
#layer[3]=0225
#layer[4]=0316
#layer[5]=0406
#layer[6]=0541
#layer[7]=0632
#layer[8]=2483
#layer[9]=2799

plot1=t363.
time1=.00_02
plot2=T_	#T_old
time2=_40
plot3=T_	# Tomo: T_xxxx_xx
time3=_40


gmtset PS_MEDIA `echo 67 + $dist_help | bc`"c"x37.5c
gmtset FONT_ANNOT_PRIMARY 16p
gmtset FONT_LABEL 20p

range11=-15
range12=15
range21=-15
range22=15
range31=-15
range32=15

#range11=-300
#range12=300
#range21=-300
#range22=300
#range31=-300
#range32=300
d1=`echo $range12/100 - $range11/100 | bc`
d2=`echo $range22/100 - $range21/100 | bc`
d3=`echo $range32/100 - $range31/100 | bc`

makecpt -Cpolar2.cpt -T$range11/$range12/$d1 -Z -D > test1.cpt
makecpt -Cpolar2.cpt -T$range21/$range22/$d2 -Z -D > test2.cpt
makecpt -Cpolar2.cpt -T$range31/$range32/$d3 -Z -D > test3.cpt

# IMAGES
for ((i=1; i <= $ilay; i++)); do

xyz2grd $path/$plot1''${layer[$i]}''$time1 -Gtest1_$i.grd -I+5 -Rd
#surface $path/$plot1''${layer[$i]}''$time1 -Gtest1_$i.grd -I0.5 -Rd
xyz2grd $path2/$plot2''${layer[$i]}''$time2 -Gtest2_$i.grd -I+5 -Rd
#surface $path2/$plot2''${layer[$i]}''$time2 -Gtest2_$i.grd -I0.5 -Rd
xyz2grd $path3/$plot3''${layer[$i]}''$time3 -Gtest3_$i.grd -I+5 -Rd

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

for((i=1;i<=$ilay;i++)); do
if [ "$i" -eq "1" ]; then
grdimage test3_$i.grd -X12.5c -Y$ydist3 -Ctest3.cpt -Rd -JW$size -K -E$res -O >> test.ps
pscoast -JW$size -Rd -W0.9p -Dc -K -O -A0/0/1 >> test.ps
else
grdimage test3_$i.grd -Y$ydist -Ctest3.cpt -Rd -JW$size -K -E$res -O >> test.ps
pscoast -JW$size -Rd -W0.9p -Dc -O -A0/0/1 -K >> test.ps
fi;
done


# SCALE
#psscale -D0c/-1c/14c/1ch  -Ctest.cpt -O -S -B5:"temperature variations [%]":/:: -K -E >> test.ps
#psscale -D-6c/-1c/8c/1ch  -Ctest1.cpt -O -S -B2e-8:"temperature variations [%]":/:: -K -E >> test.ps

#psscale -D-19c/-1c/8c/1ch  -Ctest1.cpt -O -S -B`echo $d1/20 | bc`:"temp. anomalies [%]":/:: -K -E >> test.ps
#psscale -D-6.5c/-1c/8c/1ch  -Ctest2.cpt -O -S -B`echo $d2/20 | bc`:"temp. anomalies [%]":/:: -K -E >> test.ps
#psscale -D6c/-1c/8c/1ch  -Ctest3.cpt -O -S -B`echo $d3/20 | bc`:"temp. anomalies [%]":/:: -K -E >> test.ps

#psscale -D-6c/-1c/8c/1ch  -Ctest1.cpt -O -S -B5e-10:"temp. gradient.":/:: -K -E >> test.ps
#psscale -D6c/-1c/8c/1ch  -Ctest2.cpt -O -S -B5e-10:"adj. buoyancy":/:: -K -E >> test.ps

# one large scale
#psscale -D-6c/-1c/18c/1ch  -Ctest1.cpt -O -S -B150:"adjoint temperature [K]":/:: -K -E >> test.ps
psscale -D-6c/-1c/18c/1ch  -Ctest1.cpt -O -S -B`echo $d1/20 | bc`:"temperature anomalies [%]":/:: -K -E >> test.ps


lsize=18p
# LAYERS
for((i=1;i<=$ilay;i++)); do
if [ "$i" -eq "1" ]; then
pstext -Rd -JW$size -F+f$lsize -O -K -N -Y$dsize -X-18.75c >> test.ps << END
0 0 `echo ${layer[$i]} + 0 | bc` km
END
else
pstext -Rd -JW$size -F+f$lsize -O -K -N -Y$ydist >> test.ps << END
0 0 `echo ${layer[$i]} + 0 | bc` km
END
fi;
done
for((i=1;i<=$ilay;i++)); do
if [ "$i" -eq "1" ]; then
pstext -Rd -JW$size -F+f$lsize -O -K -N -Y19.5 -X12.5c >> test.ps << END
0 0 `echo ${layer[$i]} + 0 | bc` km
END
else
pstext -Rd -JW$size -F+f$lsize -O -K -N -Y$ydist >> test.ps << END
0 0 `echo ${layer[$i]} + 0 | bc` km
END
fi;
done


# HEADLINE
pstext -Rd -JW$size -F+f22p -O -K -N -Y$dsize -X-18.5c >> test.ps << END
0 0 iteration 1
END
pstext -Rd -JW$size -F+f22p -O -K -N -X12.5c >> test.ps << END
0 0 iteration 2
END
pstext -Rd -JW$size -F+f22p -O -K -N -X12.5c >> test.ps << END
0 0 iteration 7
END


ps2pdf test.ps
mv test.pdf slices.pdf
rm test.ps
evince slices.pdf &
#mv slices.pdf ../../TERRA_Benchmark/src/code/

