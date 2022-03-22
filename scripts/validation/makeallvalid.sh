#!/bin/bash

rm -rf FDSValid
mkdir -p FDSValid

set -e
set -x

GENERATORS=( GENIE NEUT )
SPECIES=( numu numub nue nueb )
# SPECIES=( numu numub )
# SPECIES=( numu )
DETECTORS=( NOvAND ND280 )
DETECTORS=( ND280 )
# DETECTORS=( NOvAND )

declare -A DET_MATS
DET_MATS["NOvAND"]="CH"
DET_MATS["ND280"]="H2O CH"
# DET_MATS["ND280"]="CH"
# DET_MATS["ND280"]="H2O"

for DET in ${DETECTORS}; do
    for TGT in ${DET_MATS[${DET}]}; do
        for SPC in ${SPECIES[@]}; do

            if [ "${DET}" == "ND280" ]; then

                bin/fakedatavalid.exe -i flattrees/t2knova.flattree.NEUT.ND280.${TGT}.${SPC}.root \
                                    -F FDSInputs/FakeDataInputs.root \
                                    -M -H config/FakeDataValidConfig_${DET}.toml \
                                    -W T2KND_to_NOvA \
                                    -a any \
                                    -o FDSValid/FakeDataValid_T2KND_to_NOvA_${SPC}.root \
                                    -d ND280/T2KND_to_NOvA/${TGT}/${SPC} &

                bin/fakedatavalid.exe -i flattrees/t2knova.flattree.NEUT.ND280.${TGT}.${SPC}.root \
                                    -F FDSInputs/FakeDataInputs.root \
                                    -M -H config/FakeDataValidConfig_${DET}.toml \
                                    -a any \
                                    -o FDSValid/FakeDataValid_NEUT.root \
                                    -d ND280/T2KNDTune/${TGT}/${SPC} &

                bin/fakedatavalid.exe -i flattrees/t2knova.flattree.GENIE.ND280.${TGT}.${SPC}.root \
                                    -F FDSInputs/FakeDataInputs.root \
                                    -M -H config/FakeDataValidConfig_${DET}.toml \
                                    -a any \
                                    -o FDSValid/FakeDataValid_GENIE.root \
                                    -d ND280/NOvATune/${TGT}/${SPC} &

                wait

                if [ "${TGT}" == "H2O" ]; then
                    bin/fakedatavalid.exe -i flattrees/t2knova.flattree.NEUT.ND280.${TGT}.${SPC}.root \
                                        -F FDSInputs/FakeDataInputs.root \
                                        -M -H config/FakeDataValidConfig_${DET}.toml \
                                        -W T2KND_to_NOvA \
                                        -a any --oscillate \
                                        -o FDSValid/FakeDataValid_T2KND_to_NOvA_${SPC}.root \
                                        -d ND280/T2KND_to_NOvA_osc/${TGT}/${SPC} &

                    bin/fakedatavalid.exe -i flattrees/t2knova.flattree.NEUT.ND280.${TGT}.${SPC}.root \
                                        -F FDSInputs/FakeDataInputs.root \
                                        -M -H config/FakeDataValidConfig_${DET}.toml \
                                        -a any --oscillate \
                                        -o FDSValid/FakeDataValid_NEUT.root \
                                        -d ND280/T2KNDTune_osc/${TGT}/${SPC} &

                    bin/fakedatavalid.exe -i flattrees/t2knova.flattree.GENIE.ND280.${TGT}.${SPC}.root \
                                        -F FDSInputs/FakeDataInputs.root \
                                        -M -H config/FakeDataValidConfig_${DET}.toml \
                                        -a any --oscillate \
                                        -o FDSValid/FakeDataValid_GENIE.root \
                                        -d ND280/NOvATune_osc/${TGT}/${SPC} &
                    wait
                fi

            else

                bin/fakedatavalid.exe -i flattrees/t2knova.flattree.GENIE.NOvAND.${TGT}.${SPC}.root \
                                    -F FDSInputs/FakeDataInputs.root \
                                    -M -H config/FakeDataValidConfig_${DET}.toml \
                                    -W NOvA_to_T2KND_ptlep \
                                    -a any \
                                    -o FDSValid/FakeDataValid_NOvA_to_T2KND_ptlep.root \
                                    -d NOvAND/NOvA_to_T2KND_ptlep/${TGT}/${SPC} &

                bin/fakedatavalid.exe -i flattrees/t2knova.flattree.NEUT.NOvAND.${TGT}.${SPC}.root \
                                    -F FDSInputs/FakeDataInputs.root \
                                    -M -H config/FakeDataValidConfig_${DET}.toml \
                                    -a any \
                                    -o FDSValid/FakeDataValid_NEUT.root \
                                    -d NOvAND/T2KNDTune/${TGT}/${SPC} &

                bin/fakedatavalid.exe -i flattrees/t2knova.flattree.GENIE.NOvAND.${TGT}.${SPC}.root \
                                    -F FDSInputs/FakeDataInputs.root \
                                    -M -H config/FakeDataValidConfig_${DET}.toml \
                                    -a any \
                                    -o FDSValid/FakeDataValid_GENIE.root \
                                    -d NOvAND/NOvATune/${TGT}/${SPC} &

                wait

            fi
        done
        wait
    done
done
 
hadd -j 2 FDSValid/FakeDataValid.root \
          FDSValid/FakeDataValid_*.root


if [ ! -e FDSValid/FakeDataValid.root ]; then
  echo "Didn't produce FakeDataValid.root. Failed"
  exit 1
fi