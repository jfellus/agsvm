#!/bin/bash

mkdir -p results

for i in `ls -d thread*`; do
	echo "Doing $i";
	if [ -d $i/alldata ]; then
		cp s.sh $i/alldata
		cp $i/test.cfg $i/alldata
		cd $i/alldata
		./s.sh test.cfg
		cd ../..
		ls -d $i/alldata/*
		mv $i/alldata/* results
	else
		echo "no data in $i"
	fi
done
