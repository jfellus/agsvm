#!/bin/bash

cp -f config.properties.mnist config.properties

mkdir -p /local/jfellus/agsvm
cp -f config.properties /local/jfellus/agsvm
cp -f gossip_svm /local/jfellus/agsvm

cd /local/jfellus/agsvm



sed -i "s|ALGO = .*|ALGO = SAG|g" config.properties

#for sbs in 1 2 3 4 5 6 7 8 9 10 20 30 40 50 100 200 500 1000 2000 5000 10000 30000 50000; do
#sed -i "s|STAG_BUFFER_SIZE = .*|STAG_BUFFER_SIZE = $sbs|g" config.properties

for i in 0.0001 0.01 0.1 0.00001 0.000001; do
	sed -i "s|LEARNING_RATE = .*|LEARNING_RATE = $i|g" config.properties
	./gossip_svm & 
	sleep 0.1
done
#done
