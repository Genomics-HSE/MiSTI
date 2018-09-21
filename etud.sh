#!/bin/sh
CLEAN=0
echo $1
if [ "$#" -eq 1 ]; then
	if [ "$1" = "-c" ]; then
		CLEAN=1
	fi
fi

if [ $CLEAN -eq 1 ]; then
	echo True
else
	echo False
fi