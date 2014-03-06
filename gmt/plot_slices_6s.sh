#!/bin/bash

#path=/import/netapp-terra/horbach/adjoint/runs/20.1_40Ma_BACKWARD_LRZ/plotdata
path=/import/netapp-m-02-terra/horbach/adjoint/runs/15.2_40Ma_mt512_1020_PHASE2/gmt_It1
#path=/import/netapp-m-02-terra/horbach/adjoint/runs/18.2_20Ma_mt256/Iteration_05/plotdata
#path2=/import/netapp-terra/horbach/adjoint/Endfeld/TOMO_NORESTRICT/plotdata

size=12
res=100
ydist=-6.5c

lsize=18p

psize_tmp=57

layer[1]=0270
layer[2]=0630
layer[3]=1720
layer[4]=2800

plot1=t381
time1=_00

psize=$psize_tmp"c"
dsize=`echo $psize_tmp - 2.5 | bc`"c"
ydist3=`echo $psize_tmp - 5 | bc`"c"   
ydist4=`echo $psize_tmp - 7.5 | bc`"c"

gmtset PS_MEDIA 67cx75c
gmtset FONT_ANNOT_PRIMARY 16p
gmtset FONT_LABEL 20p

#makecpt -Cpolar2.cpt -T-500/500/10 -Z -D > test.cpt
makecpt -Cpolar2.cpt -T-20/20/0.5 -Z -D > test.cpt

for i in {1..6}; do
j=`echo 4*$i - 4| bc`
if [ "$j" -ge "10" ]; then
xyz2grd $path/$plot1''.${layer[$i]}''.$j''$time1 -Gtest1_$i.grd -I1 -Rd
else
xyz2grd $path/$plot1''.${layer[$i]}''.0$j''$time1 -Gtest1_$i.grd -I1 -Rd

fi;
done

# IMAGES
for i in {1..6}; do
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

for j in {2..9}; do
grdimage test0$j"_"$i".grd" -Y$ydist  -Ctest.cpt -Rd -JW$size -K -E$res -O >> test.ps
pscoast -JW$size -Rd -W0.9p -Dc -O -A0/0/1 -K >> test.ps;
done;


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


# SCALE
psscale -D-19c/-4c/20c/1ch  -Ctest.cpt -O -S -B5:"temperature variations [%]":/:: -K -E >> test.ps

# HEADLINE
pstext -Rd -JW$size -F+f22p -O -K -N -Y$dsize -X-56c >> test.ps << END
0 0 Model (20 Ma)
END
pstext -Rd -JW$size -F+f22p -O -K -N -X12.5c >> test.ps << END
0 0 Model (16 Ma)
END
pstext -Rd -JW$size -F+f22p -O -K -N -X12.5c >> test.ps << END
0 0 Model (12 Ma)
END
pstext -Rd -JW$size -F+f22p -O -K -N -X12.5c >> test.ps << END
0 0 Model (8 Ma)
END
pstext -Rd -JW$size -F+f22p -O -K -N -X12.5c >> test.ps << END
0 0 Model (4 Ma)
END
pstext -Rd -JW$size -F+f22p -O -K -N -X12.5c >> test.ps << END
0 0 Model (0 Ma)
END


ps2pdf test.ps
mv test.pdf slicesx2.pdf
rm test.ps
evince slicesx2.pdf

