#!/bin/bash

#cp -f config.properties.mnist config.properties

#mkdir -p /local/jfellus/agsvm
#cp -f config.properties /local/jfellus/agsvm
#cp -f gossip_svm /local/jfellus/agsvm

#cd /local/jfellus/agsvm

#rm -rf data/*

sed -i "s|ALGO = .*|ALGO = SAG|g" config.properties

#for sbs in 1 2 3 4 5 6 7 8 9 10 20 30 40 50 100 200 500 1000 2000 5000 10000 30000 50000; do
#sed -i "s|STAG_BUFFER_SIZE = .*|STAG_BUFFER_SIZE = $sbs|g" config.properties

for i in 0.1 0.01 0.001 0.0001 0.05 0.005 0.0005; do
	sed -i "s|LAMBDA = .*|LAMBDA = $i|g" config.properties
	./agsvm & 
	sleep 0.1
done
#done

read zob
killall agsvm
