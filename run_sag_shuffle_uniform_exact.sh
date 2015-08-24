#!/bin/bash

#cp -f config.properties.mnist config.properties

#mkdir -p /local/jfellus/agsvm
#cp -f config.properties /local/jfellus/agsvm
#cp -f gossip_svm /local/jfellus/agsvm

#cd /local/jfellus/agsvm

#rm -rf data/*

lr=0.01


sed -i "s|ALGO = .*|ALGO = SAG|g" config.properties
sed -i "s|LEARNING_RATE = .*|LEARNING_RATE = ${lr}|g" config.properties

sed -i "s|SHUFFLE_DATASET = .*|SHUFFLE_DATASET = 0|g" config.properties
sed -i "s|EXACT_REGUL = .*|EXACT_REGUL = 0|g" config.properties
./agsvm uniform_reg & 
sleep 0.1

sed -i "s|SHUFFLE_DATASET = .*|SHUFFLE_DATASET = 1|g" config.properties
sed -i "s|EXACT_REGUL = .*|EXACT_REGUL = 0|g" config.properties
./agsvm shuffle_reg & 
sleep 0.1

sed -i "s|SHUFFLE_DATASET = .*|SHUFFLE_DATASET = 0|g" config.properties
sed -i "s|EXACT_REGUL = .*|EXACT_REGUL = 1|g" config.properties
./agsvm uniform_exact & 
sleep 0.1

sed -i "s|SHUFFLE_DATASET = .*|SHUFFLE_DATASET = 1|g" config.properties
sed -i "s|EXACT_REGUL = .*|EXACT_REGUL = 1|g" config.properties
./agsvm shuffle_exact & 
sleep 0.1

read zob
killall agsvm
