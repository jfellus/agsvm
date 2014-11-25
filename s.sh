#!/bin/bash

a=""
i=0
while read line; do
	if [ "a$line" = "a-----" ]; then
		if [ -d data_$i ]; then mv "data_$i" "${a// /}"; fi 
		i=$(($i+1))
		a=""
	else
		if [ -z "$a" ]; then a=$line
		else a="$a
$line"
		fi
	fi
done < $1
