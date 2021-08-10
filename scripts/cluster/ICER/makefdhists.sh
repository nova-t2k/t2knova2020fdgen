#!/bin/bash

rm FakeDataHists.root
g++ fakedatahists.C -I ../../../ $(root-config --cflags --glibs) -o fakedatahists.exe
if [[ ! "$?" == "0" ]]; then
	exit;
fi

g++ fakedatarwgen.C -I ../../../ -o fakedatarwgen.exe $(root-config --cflags) $(root-config --libs)
if [[ ! "$?" == "0" ]]; then
        exit;
fi

declare -A TGTEL
TGTEL["NOvAND_CH"]="C H"
TGTEL["ND280_CH"]="C"
TGTEL["ND280_H2O"]="H O"
#TGTEL["ND280_H2O"]="O"

GENERATORS=( GENIE NEUT )
SPECIES=( numu numub nue nueb )
SPECIES=( numu )
DETECTORS=( NOvAND ND280 )
DETECTORS=( ND280 )

declare -A DET_MATS
DET_MATS["NOvAND"]="CH"
DET_MATS["ND280"]="H2O CH"

for gen in ${GENERATORS[@]}; do
  for spec in ${SPECIES[@]}; do
    for det in ${DETECTORS[@]}; do
      for mat in ${DET_MATS[${det}]}; do

        if [ ! -e "t2knova.flattree.${gen}.${det}.${mat}.${spec}.root" ]; then
          continue
        fi

        for tgtel in ${TGTEL["${det}_${mat}"]}; do

          CMD="./fakedatahists.exe t2knova.flattree.${gen}.${det}.${mat}.${spec}.root \
                                   ${det} \
                                   FakeDataHists.root \
                                   ${tgtel} \
                                   ${gen}/${det}/${tgtel}/${spec}"
          echo $CMD
          ${CMD}
        done
      done
    done
  done
done

echo ./fakedatarwgen.exe FakeDataHists.root FakeDataInputs.root
./fakedatarwgen.exe FakeDataHists.root FakeDataInputs.root
