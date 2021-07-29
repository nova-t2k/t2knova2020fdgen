#!/bin/bash

rm FakeDataHists.root
g++ fakedatahists.C $(root-config --cflags --glibs) -o fakedatahists.exe
if [[ ! "$?" == "0" ]]; then
	exit;
fi


declare -A TGTEL
TGTEL["NOvAND_CH"]="C H"
TGTEL["ND280_CH"]="C"
TGTEL["ND280_H2O"]="H O"

GENERATORS=( GENIE NEUT )
SPECIES=( numu numub nue nueb )
DETECTORS=( NOvAND ND280 )

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

          echo ./fakedatahists.exe t2knova.flattree.${gen}.${det}.${mat}.${spec}.root \
                                   ${det} \
                                   FakeDataHists.root \
                                   ${tgtel} \
                                   ${gen}/${det}/${tgtel}/${spec}

        done
      done
    done
  done
done

echo root -l -b -q "fakedatarwgen.C(\"FakeDataHists.root\",\"FakeDataInputs.root\")"
