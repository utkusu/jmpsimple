#!/bin/bash
if [ $1 -eq  1 ]
then
	echo "Switching from serial to MPI makefile"
	cp makefile makefileserial
	rm makefile
	cp makefilempi makefile
else
	echo "Switching from MPI to Serial makefile"
	cp makefile makefilempi
	rm makefile
	cp makefileserial makefile
fi

