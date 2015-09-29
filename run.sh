#!/bin/bash

#cp -f config.properties.mnist config.properties

#mkdir -p /local/jfellus/agsvm
#cp -f config.properties /local/jfellus/agsvm
#cp -f gossip_svm /local/jfellus/agsvm

#cd /local/jfellus/agsvm

#rm -rf data/*

sed -i "s|ALGO = .*|ALGO = SGD|g" config.properties

#for sbs in 1 2 3 4 5 6 7 8 9 10 20 30 40 50 100 200 500 1000 2000 5000 10000 30000 50000; do


#for u in 1 2 3 4 5 6 7 8 9 10; do
for N in 100; do
for l in 1 0.5 0.05 0.01 0.06 0.07 0.08 0.09 0.02 0.03 0.04; do
for i in 100; do
	sed -i "s|N = .*|N = $N|g" config.properties
	sed -i "s|LEARNING_RATE = .*|LEARNING_RATE = $l|g" config.properties
#	sed -i "s|STAG_BUFFER_SIZE = .*|STAG_BUFFER_SIZE = $i|g" config.properties
	sed -i "s|MESSAGES = .*|MESSAGES = $i|g" config.properties

	./agsvm & 
	sleep 0.1
done
#done
done
done

read zob
killall agsvm
