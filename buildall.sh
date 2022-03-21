#!/bin/bash

#I know I should have written a make file

if [ ! -e "Prob3-Wrapper" ]; then
    git clone https://github.com/luketpickering/Prob3-Wrapper.git
fi
if [ ! -e "Prob3-Wrapper/WrappedProb3plusplus.so" ]; then
    cd Prob3-Wrapper
    ./build.sh
fi

mkdir -p bin

LOF=$(ls src/*/*.C)
if [ ! -z ${1} ] && [ "${1}" != "all" ]; then
    LOF=${1}
fi

for i in ${LOF}; do

    HDRSUPDATED=0

    if [ "${1}" == "all" ]; then
        HDRSUPDATED=1
    fi

    o=${i##*/}

    for j in include/*.h T2KNOvA/*.hxx; do
        if [ "${j}" -nt "bin/${o%%.C}.exe" ]; then
            HDRSUPDATED=1
        fi
    done


    if [ ! -e bin/${o%%.C}.exe ] \
        || [ "${i}" -nt "bin/${o%%.C}.exe" ] \
        || [ "${HDRSUPDATED}" == "1" ]; then

        PROB3PP=""
        if [ "${o}" == "fakedatavalid.C" ]; then
            PROB3PP="-I Prob3-Wrapper Prob3-Wrapper/WrappedProb3plusplus.so"
        fi

        CMD="g++ -I include/ -I ./ -O2 -g ${i} -o bin/${o%%.C}.exe $(root-config --glibs --cflags) $(root-config --glibs --cflags) ${PROB3PP}"
        echo $CMD
        ${CMD}
        if [[ ! "$?" == "0" ]]; then
            exit;
        fi
    fi
done