#!/bin/bash

for i in ../pyro.*
	#|grep -v 1401 | sed 's/.*\.\./\.\./' |grep -v total`
do
	nl=`wc -l $i | sed 's/ .*//'`

	if [ $nl -ne 1401 ]
	then

	nla=`echo 1401-$nl | bc`

	suffix=`echo $i |sed 's/.*\.//'`

	sed -n 1,${nla}p my_pyro.0000 | sed 's/AAAA/'$suffix'/' > mytest_`echo $i |sed 's/\.\.\///'`
	cat $i >> mytest_`echo $i |sed 's/\.\.\///'`

	else
		cp $i mytest_`echo $i |sed 's/\.\.\///'`
	fi

done
