#!/bin/bash

set -x
set -e

rm -rf FDSInputs
mkdir -p FDSInputs

declare -A TGTEL
TGTEL["NOvAND_CH"]="C H"

TGTEL["ND280_CH"]="C"
TGTEL["ND280_H2O"]="H O"

GENERATORS=( GENIE NEUT )

SPECIES=( numu numub nue nueb )
# SPECIES=( numu numub )

DETECTORS=( NOvAND ND280 )
# DETECTORS=( ND280 )

declare -A DET_MATS
DET_MATS["NOvAND"]="CH"
DET_MATS["ND280"]="H2O CH"
# DET_MATS["ND280"]="CH"

declare -A TUNES
TUNES["NEUT"]="BANFF_PRE BANFF_POST"
TUNES["GENIE"]="2020"

declare -A FakeDataSets
FakeDataSets["BANFF_PRE"]="Generated NDTuned Mnv1Pi NonQE"
FakeDataSets["BANFF_POST"]="NDTuned"
FakeDataSets["2020"]="Generated NDTuned"

for gen in ${GENERATORS[@]}; do
  for det in ${DETECTORS[@]}; do
    for mat in ${DET_MATS[${det}]}; do
      for tgtel in ${TGTEL["${det}_${mat}"]}; do
        for tune in ${TUNES["${gen}"]}; do
          for fds in ${FakeDataSets["${tune}"]}; do
            for spec in ${SPECIES[@]}; do

              if [ ! -e "flattrees/t2knova.flattree.${gen}.${det}.${mat}.${spec}.root" ]; then
                echo "Failed to find: flattrees/t2knova.flattree.${gen}.${det}.${mat}.${spec}.root"
                continue
              fi

              LASTDIRNAME="${fds}"
              if [ "${LASTDIRNAME}" = "NDTuned" ]; then
                LASTDIRNAME=${tune}
              fi

              CMD="bin/fakedatahists.exe -i flattrees/t2knova.flattree.${gen}.${det}.${mat}.${spec}.root \
                                         -H config/FakeDataConfig.toml \
                                         -e ${det} \
                                         --FDS ${fds} \
                                         -o FDSInputs/FakeDataHists_${spec}.root \
                                         -M \
                                         -a ${tgtel} \
                                         -d ${gen}/${det}/${tgtel}/${spec}/${LASTDIRNAME}"
              echo $CMD
              ${CMD} &

            done
            wait
          done
        done
      done
    done
  done
done

HADD_CMD="hadd -j 4 FDSInputs/FakeDataHists.root"
for spec in ${SPECIES[@]}; do
  HADD_CMD="${HADD_CMD} FDSInputs/FakeDataHists_${spec}.root"
done
echo $HADD_CMD
${HADD_CMD}

if [ ! -e FDSInputs/FakeDataHists.root ]; then
  echo "Didn't produce FakeDataHists.root. Failed"
  exit 1
fi

echo bin/fakedatarwgen.exe FDSInputs/FakeDataHists.root FDSInputs/FakeDataInputs.root
bin/fakedatarwgen.exe FDSInputs/FakeDataHists.root FDSInputs/FakeDataInputs.root
 