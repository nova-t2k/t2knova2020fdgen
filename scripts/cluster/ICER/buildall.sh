#!/bin/bash

for i in *.C; do

	HDRSUPDATED=0

	for j in ./*.h ../../../*.hxx; do
		if [ "${j}" -nt "${i%%.C}.exe" ]; then
			HDRSUPDATED=1
		fi
	done

	if [ "${i}" -nt "${i%%.C}.exe" ] || [ "${HDRSUPDATED}" == "1" ]; then
		CMD="g++ -I . -I ../../../ $(root-config --glibs --cflags) -std=c++11 -O2 -g ${i} -o ${i%%.C}.exe"
		echo $CMD
		${CMD}
		if [[ ! "$?" == "0" ]]; then
			exit;
		fi
	fi
done