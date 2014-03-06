#!/bin/bash

#10^21: 18.4
path=/import/netapp-m-02-terra/horbach/adjoint/runs/18.4_40Ma_mt512_10pres_1021_F121/gmt_It2
plot1=t363
it1=_02
#10^21 PHASE4: 15.1
path=/import/netapp-m-02-terra/horbach/adjoint/runs/15.1_40Ma_mt256_1021_PHASE4/Iteration_07/plotdata
plot1=T
it1=''
#10^20, pres10: 18.2
#path=/import/netapp-m-02-terra/horbach/adjoint/runs/18.2_40Ma_mt512_10pres_1020/Iteration_03/plotdata
#plot1=T
#it1=''
#10^20, pres5: 18.1
path=/import/netapp-m-02-terra/horbach/adjoint/runs/18.1_40Ma_mt512_5pres_1020/Iteration_03/plotdata
plot1=T
it1=''
#10^20 PHASE2: 15.2
#path=/import/netapp-m-02-terra/horbach/adjoint/runs/15.2_40Ma_mt512_1020_PHASE2/gmt_It1
#plot1=t381
#it1=_02
#10^19: 18.3
path=/import/netapp-m-02-terra/horbach/adjoint/runs/18.3_40Ma_mt512_15pres_1019/gmt_It1
plot1=t373
it1=_02

path2=$path

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

layer[1]=0270
layer[2]=0540
layer[3]=1310
layer[4]=2840

#time1=.40
#time2=.32
#time3=.20
#time4=.08
#time5=.00

time1=.00
time2=.10
time3=.20
time4=.30
time5=.40

plot2=$plot1
plot3=$plot1
plot4=$plot1
plot5=$plot1
it2=$it1
it3=$it1
it4=$it1
it5=$it1

gmtset PS_MEDIA `echo 67 + $dist_help | bc`"c"x72c
gmtset FONT_ANNOT_PRIMARY 16p
gmtset FONT_LABEL 20p

range11=-15
range12=15
range21=-15
range22=15

d1=`echo $range12/100 - $range11/100 | bc`
d2=`echo $range22/100 - $range21/100 | bc`

makecpt -Cpolar2.cpt -T$range11/$range12/$d1 -Z -D > test1.cpt
makecpt -Cpolar2.cpt -T$range21/$range22/$d2 -Z -D > test2.cpt

# IMAGES
for ((i=1; i <= $ilay; i++)); do

#surface $path/$plot1''.${layer[$i]}''$time1 -Gtest1_$i.grd -I0.5 -Rd
#grdmath test1_$i.grd DIV = test1_$i.grd
#grdmath test2_$i.grd ABS = test2_$i.grd
#grdmath test1_$i.grd test2_$i.grd SUB = test1_$i.grd

xyz2grd $path/$plot1''.${layer[$i]}''$time1''$it1 -Gtest1_$i.grd -I+5 -Rd
if [ "$i" -eq "1" ]; then
grdimage test1_$i.grd -Y$psize -X0.2c -Ctest1.cpt -Rd -JW$size -K -E$res > test.ps
pscoast -JW$size -Rd -W0.9p -Dc -K -O -A0/0/1 >> test.ps
else
grdimage test1_$i.grd -Y$ydist  -Ctest1.cpt -Rd -JW$size -K -E$res -O >> test.ps
pscoast -JW$size -Rd -W0.9p -Dc -O -A0/0/1 -K >> test.ps
fi;
done


for((i=1;i<=$ilay;i++)); do
xyz2grd $path2/$plot2''.${layer[$i]}''$time2''$it2 -Gtest1_$i.grd -I+5 -Rd
if [ "$i" -eq "1" ]; then
grdimage test1_$i.grd -X12.5c -Y$ydist3 -Ctest2.cpt -Rd -JW$size -K -E$res -O >> test.ps
pscoast -JW$size -Rd -W0.9p -Dc -K -O -A0/0/1 >> test.ps
else
grdimage test1_$i.grd -Y$ydist -Ctest2.cpt -Rd -JW$size -K -E$res -O >> test.ps
pscoast -JW$size -Rd -W0.9p -Dc -O -A0/0/1 -K >> test.ps
fi;
done

for((i=1;i<=$ilay;i++)); do
xyz2grd $path2/$plot3''.${layer[$i]}''$time3''$it3 -Gtest1_$i.grd -I+5 -Rd
if [ "$i" -eq "1" ]; then
grdimage test1_$i.grd -X12.5c -Y$ydist3 -Ctest2.cpt -Rd -JW$size -K -E$res -O >> test.ps
pscoast -JW$size -Rd -W0.9p -Dc -K -O -A0/0/1 >> test.ps
else
grdimage test1_$i.grd -Y$ydist -Ctest2.cpt -Rd -JW$size -K -E$res -O >> test.ps
pscoast -JW$size -Rd -W0.9p -Dc -O -A0/0/1 -K >> test.ps
fi;
done

for((i=1;i<=$ilay;i++)); do
xyz2grd $path2/$plot4''.${layer[$i]}''$time4''$it4 -Gtest1_$i.grd -I+5 -Rd
if [ "$i" -eq "1" ]; then
grdimage test1_$i.grd -X12.5c -Y$ydist3 -Ctest2.cpt -Rd -JW$size -K -E$res -O >> test.ps
pscoast -JW$size -Rd -W0.9p -Dc -K -O -A0/0/1 >> test.ps
else
grdimage test1_$i.grd -Y$ydist -Ctest2.cpt -Rd -JW$size -K -E$res -O >> test.ps
pscoast -JW$size -Rd -W0.9p -Dc -O -A0/0/1 -K >> test.ps
fi;
done

for((i=1;i<=$ilay;i++)); do
xyz2grd $path2/$plot5''.${layer[$i]}''$time5''$it5 -Gtest1_$i.grd -I+5 -Rd
if [ "$i" -eq "1" ]; then
grdimage test1_$i.grd -X12.5c -Y$ydist3 -Ctest2.cpt -Rd -JW$size -K -E$res -O >> test.ps
pscoast -JW$size -Rd -W0.9p -Dc -K -O -A0/0/1 >> test.ps
else
grdimage test1_$i.grd -Y$ydist -Ctest2.cpt -Rd -JW$size -K -E$res -O >> test.ps
pscoast -JW$size -Rd -W0.9p -Dc -O -A0/0/1 -K >> test.ps
fi;
done


# SCALE
#psscale -D0c/-1c/14c/1ch  -Ctest.cpt -O -S -B5:"temperature variations [%]":/:: -K -E >> test.ps
#psscale -D-6c/-1c/8c/1ch  -Ctest1.cpt -O -S -B2e-8:"temperature variations [%]":/:: -K -E >> test.ps

#psscale -D-6c/-1c/8c/1ch  -Ctest1.cpt -O -S -B`echo $d1/20 | bc`:"velocity [cm/yr]":/:: -K -E >> test.ps
#psscale -D6c/-1c/8c/1ch  -Ctest2.cpt -O -S -B`echo $d2/20 | bc`:"velocity [cm/yr]":/:: -K -E >> test.ps

#psscale -D-6c/-1c/8c/1ch  -Ctest1.cpt -O -S -B`echo $d1/20 | bc`:"density var. [%]":/:: -K -E >> test.ps
psscale -D6c/-1c/8c/1ch  -Ctest2.cpt -O -S -B`echo $d2/20 | bc`:"temperature var. [%]":/:: -K -E >> test.ps

#psscale -D-6c/-1c/8c/1ch  -Ctest1.cpt -O -S -B5e-10:"temp. gradient.":/:: -K -E >> test.ps
#psscale -D6c/-1c/8c/1ch  -Ctest2.cpt -O -S -B5e-10:"adj. buoyancy":/:: -K -E >> test.ps
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
pstext -Rd -JW$size -F+f22p -O -K -N -Y$dsize -X-43.5c >> test.ps << END
0 0 40 Ma
END
pstext -Rd -JW$size -F+f22p -O -K -N -X12.5c >> test.ps << END
0 0 30 Ma
END
pstext -Rd -JW$size -F+f22p -O -K -N -X12.5c >> test.ps << END
0 0 20 Ma
END
pstext -Rd -JW$size -F+f22p -O -K -N -X12.5c >> test.ps << END
0 0 10 Ma
END
pstext -Rd -JW$size -F+f22p -O -K -N -X12.5c >> test.ps << END
0 0 0 Ma
END


ps2pdf test.ps
mv test.pdf slices.pdf
rm test.ps
evince slices.pdf &
#mv slices.pdf ../../TERRA_Benchmark/src/code/

