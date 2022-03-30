#!/bin/bash

if [ "${1}" = "force" ]; then
    rm -rf FDSValid
    mkdir -p FDSValid
fi

set -e
set -x

SPECIES=( numu numub nue nueb )
SPECIES=( numu numub )
# SPECIES=( numu )
DETECTORS=( NOvAND ND280 )
# DETECTORS=( ND280 )
# DETECTORS=( NOvAND )

declare -A DET_MATS
DET_MATS["NOvAND"]="CH"
DET_MATS["ND280"]="H2O CH"
# DET_MATS["ND280"]="CH"
# DET_MATS["ND280"]="H2O"

declare -A MAT_ELEMENTS
MAT_ELEMENTS["CH"]="C H"
# MAT_ELEMENTS["CH"]="C"
MAT_ELEMENTS["H2O"]="H O"
# MAT_ELEMENTS["CH"]=""
# MAT_ELEMENTS["H2O"]=""

DO_MAIN=0
DO_MAT_ELEMENTS=1
DO_OSC=0
DO_UNTUNED=1
DO_ND280=1
DO_NOvAND=1

for DET in ${DETECTORS[@]}; do
    for TGT in ${DET_MATS[${DET}]}; do
        for SPC in ${SPECIES[@]}; do

            if [ "${DET}" == "ND280" ] && [ "${DO_ND280}" == "1" ]; then

                if [ "${DO_MAIN}" == "1" ]; then
                    bin/fakedatavalid.exe -i flattrees/t2knova.flattree.NEUT.ND280.${TGT}.${SPC}.root \
                                        -F FDSInputs/FakeDataInputs.root \
                                        -M -H config/FakeDataValidConfig_${DET}.toml \
                                        -W T2KND_to_NOvA \
                                        -a any \
                                        -o FDSValid/FakeDataValid_T2KND_to_NOvA.root \
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
                fi

                if [ "${DO_UNTUNED}" == "1" ]; then
                    bin/fakedatavalid.exe -i flattrees/t2knova.flattree.NEUT.ND280.${TGT}.${SPC}.root \
                                        -F FDSInputs/FakeDataInputs.root \
                                        -M -H config/FakeDataValidConfig_${DET}.toml \
                                        -a any --No-Tune \
                                        -o FDSValid/FakeDataValid_NEUT.root \
                                        -d ND280/NEUT/${TGT}/${SPC} &

                    bin/fakedatavalid.exe -i flattrees/t2knova.flattree.GENIE.ND280.${TGT}.${SPC}.root \
                                        -F FDSInputs/FakeDataInputs.root \
                                        -M -H config/FakeDataValidConfig_${DET}.toml \
                                        -a any --No-Tune \
                                        -o FDSValid/FakeDataValid_GENIE.root \
                                        -d ND280/GENIE/${TGT}/${SPC} &

                    wait
                fi

                if [ "${DO_OSC}" == "1" ] && [ "${TGT}" == "H2O" ]; then
                    bin/fakedatavalid.exe -i flattrees/t2knova.flattree.NEUT.ND280.${TGT}.${SPC}.root \
                                        -F FDSInputs/FakeDataInputs.root \
                                        -M -H config/FakeDataValidConfig_${DET}.toml \
                                        -W T2KND_to_NOvA \
                                        -a any --oscillate \
                                        -o FDSValid/FakeDataValid_T2KND_to_NOvA.root \
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

                if [ "${DO_MAT_ELEMENTS}" == "1" ]; then 
                    for ELE in ${MAT_ELEMENTS["${TGT}"]}; do

                        if [ "${DO_MAIN}" == "1" ]; then
                            bin/fakedatavalid.exe -i flattrees/t2knova.flattree.NEUT.ND280.${TGT}.${SPC}.root \
                                                -F FDSInputs/FakeDataInputs.root \
                                                -M -H config/FakeDataValidConfig_${DET}.toml \
                                                -W T2KND_to_NOvA \
                                                -a ${ELE} \
                                                -o FDSValid/FakeDataValid_T2KND_to_NOvA.root \
                                                -d ND280/T2KND_to_NOvA/${ELE}/${SPC} &

                            bin/fakedatavalid.exe -i flattrees/t2knova.flattree.NEUT.ND280.${TGT}.${SPC}.root \
                                                -F FDSInputs/FakeDataInputs.root \
                                                -M -H config/FakeDataValidConfig_${DET}.toml \
                                                -a ${ELE} \
                                                -o FDSValid/FakeDataValid_NEUT.root \
                                                -d ND280/T2KNDTune/${ELE}/${SPC} &

                            bin/fakedatavalid.exe -i flattrees/t2knova.flattree.GENIE.ND280.${TGT}.${SPC}.root \
                                                -F FDSInputs/FakeDataInputs.root \
                                                -M -H config/FakeDataValidConfig_${DET}.toml \
                                                -a ${ELE} \
                                                -o FDSValid/FakeDataValid_GENIE.root \
                                                -d ND280/NOvATune/${ELE}/${SPC} &

                            wait
                        fi

                        if [ "${DO_UNTUNED}" == "1" ]; then
                            bin/fakedatavalid.exe -i flattrees/t2knova.flattree.NEUT.ND280.${TGT}.${SPC}.root \
                                                -F FDSInputs/FakeDataInputs.root \
                                                -M -H config/FakeDataValidConfig_${DET}.toml \
                                                -a ${ELE} --No-Tune \
                                                -o FDSValid/FakeDataValid_NEUT.root \
                                                -d ND280/NEUT/${ELE}/${SPC} &

                            bin/fakedatavalid.exe -i flattrees/t2knova.flattree.GENIE.ND280.${TGT}.${SPC}.root \
                                                -F FDSInputs/FakeDataInputs.root \
                                                -M -H config/FakeDataValidConfig_${DET}.toml \
                                                -a ${ELE} --No-Tune \
                                                -o FDSValid/FakeDataValid_GENIE.root \
                                                -d ND280/GENIE/${ELE}/${SPC} &

                            wait
                        fi
                    done
                fi

            elif [ "${DET}" == "NOvAND" ] && [ "${DO_NOvAND}" == "1" ]; then

                if [ "${DO_MAIN}" == "1" ]; then
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

                if [ "${DO_UNTUNED}" == "1" ]; then
                    bin/fakedatavalid.exe -i flattrees/t2knova.flattree.NEUT.NOvAND.${TGT}.${SPC}.root \
                                        -F FDSInputs/FakeDataInputs.root \
                                        -M -H config/FakeDataValidConfig_${DET}.toml \
                                        -a any --No-Tune \
                                        -o FDSValid/FakeDataValid_NEUT.root \
                                        -d NOvAND/NEUT/${TGT}/${SPC} &

                    bin/fakedatavalid.exe -i flattrees/t2knova.flattree.GENIE.NOvAND.${TGT}.${SPC}.root \
                                        -F FDSInputs/FakeDataInputs.root \
                                        -M -H config/FakeDataValidConfig_${DET}.toml \
                                        -a any --No-Tune \
                                        -o FDSValid/FakeDataValid_GENIE.root \
                                        -d NOvAND/GENIE/${TGT}/${SPC} &

                    wait
                fi

                if [ "${DO_MAT_ELEMENTS}" == "1" ]; then

                    for ELE in ${MAT_ELEMENTS["${TGT}"]}; do
                        if [ "${DO_MAIN}" == "1" ]; then
                            bin/fakedatavalid.exe -i flattrees/t2knova.flattree.GENIE.NOvAND.${TGT}.${SPC}.root \
                                                -F FDSInputs/FakeDataInputs.root \
                                                -M -H config/FakeDataValidConfig_${DET}.toml \
                                                -W NOvA_to_T2KND_ptlep \
                                                -a ${ELE} \
                                                -o FDSValid/FakeDataValid_NOvA_to_T2KND_ptlep.root \
                                                -d NOvAND/NOvA_to_T2KND_ptlep/${ELE}/${SPC} &

                            bin/fakedatavalid.exe -i flattrees/t2knova.flattree.NEUT.NOvAND.${TGT}.${SPC}.root \
                                                -F FDSInputs/FakeDataInputs.root \
                                                -M -H config/FakeDataValidConfig_${DET}.toml \
                                                -a ${ELE} \
                                                -o FDSValid/FakeDataValid_NEUT.root \
                                                -d NOvAND/T2KNDTune/${ELE}/${SPC} &

                            bin/fakedatavalid.exe -i flattrees/t2knova.flattree.GENIE.NOvAND.${TGT}.${SPC}.root \
                                                -F FDSInputs/FakeDataInputs.root \
                                                -M -H config/FakeDataValidConfig_${DET}.toml \
                                                -a ${ELE} \
                                                -o FDSValid/FakeDataValid_GENIE.root \
                                                -d NOvAND/NOvATune/${ELE}/${SPC} &

                            wait
                        fi

                        if [ "${DO_UNTUNED}" == "1" ]; then
                            bin/fakedatavalid.exe -i flattrees/t2knova.flattree.NEUT.NOvAND.${TGT}.${SPC}.root \
                                                -F FDSInputs/FakeDataInputs.root \
                                                -M -H config/FakeDataValidConfig_${DET}.toml \
                                                -a ${ELE} --No-Tune \
                                                -o FDSValid/FakeDataValid_NEUT.root \
                                                -d NOvAND/NEUT/${ELE}/${SPC} &

                            bin/fakedatavalid.exe -i flattrees/t2knova.flattree.GENIE.NOvAND.${TGT}.${SPC}.root \
                                                -F FDSInputs/FakeDataInputs.root \
                                                -M -H config/FakeDataValidConfig_${DET}.toml \
                                                -a ${ELE} --No-Tune \
                                                -o FDSValid/FakeDataValid_GENIE.root \
                                                -d NOvAND/GENIE/${ELE}/${SPC} &

                            wait
                        fi

                    done
                fi
            fi
        done
    done
done

rm FDSValid/FakeDataValid.root
hadd -j 2 FDSValid/FakeDataValid.root \
          FDSValid/FakeDataValid_*.root

if [ ! -e FDSValid/FakeDataValid.root ]; then
  echo "Didn't produce FakeDataValid.root. Failed"
  exit 1
fi