#!/bin/bash

#cp -f config.properties.mnist config.properties

#mkdir -p /local/jfellus/agsvm
#cp -f config.properties /local/jfellus/agsvm
#cp -f gossip_svm /local/jfellus/agsvm

#cd /local/jfellus/agsvm

#rm -rf data/*

sed -i "s|ALGO = .*|ALGO = STAG|g" config.properties

#for sbs in 1 2 3 4 5 6 7 8 9 10 20 30 40 50 100 200 500 1000 2000 5000 10000 30000 50000; do
#sed -i "s|STAG_BUFFER_SIZE = .*|STAG_BUFFER_SIZE = $sbs|g" config.properties


for N in 1000; do
for l in 1 0.5 0.1 0.01 0.001 0.0001 0.05 0.005 0.0005 0.00001 0.000001; do 
for i in 30; do
	sed -i "s|N = .*|N = $N|g" config.properties
	sed -i "s|LEARNING_RATE = .*|LEARNING_RATE = $l|g" config.properties
#	sed -i "s|STAG_BUFFER_SIZE = .*|STAG_BUFFER_SIZE = $i|g" config.properties
	sed -i "s|NB_MESSAGES = .*|NB_MESSAGES = $i|g" config.properties
	sed -i "s|STAG_BUFFER_SIZE = .*|STAG_BUFFER_SIZE = 100|g" config.properties

	./agsvm & 
	sleep 0.1
done
done
done

read zob
killall agsvm
