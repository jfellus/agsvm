#!/bin/bash

cp -f config.properties.mnist config.properties

sed -i "s|ALGO = .*|ALGO = STAG|g" config.properties

for sbs in 1 2 3 4 5 6 7 8 9 10 20 30 40 50 100 200 500 1000 2000 5000 10000 30000 50000; do
sed -i "s|STAG_BUFFER_SIZE = .*|STAG_BUFFER_SIZE = $sbs|g" config.properties

for i in 10 100 50 200 40 30 20; do
	sed -i "s|LEARNING_RATE = .*|LEARNING_RATE = $i|g" config.properties
	./gossip_svm & 
	sleep 0.1
done
done
