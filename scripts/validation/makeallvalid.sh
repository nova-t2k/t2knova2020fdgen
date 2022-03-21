#!/bin/bash

rm -rf FDSValid
mkdir -p FDSValid

set -e
set -x

GENERATORS=( GENIE NEUT )
SPECIES=( numu numub nue nueb )
# SPECIES=( numu numub )
DETECTORS=( NOvAND ND280 )
DETECTORS=( ND280 )
# DETECTORS=( NOvAND )

declare -A DET_MATS
DET_MATS["NOvAND"]="CH"
DET_MATS["ND280"]="H2O CH"
# DET_MATS["ND280"]="CH"

declare -A MAT_ELEMENTS
MAT_ELEMENTS["CH"]="C H"
# MAT_ELEMENTS["CH"]="C"
MAT_ELEMENTS["H2O"]="H O"
# MAT_ELEMENTS["CH"]=""
# MAT_ELEMENTS["H2O"]=""

for DET in ${DETECTORS}; do
    for TGT in ${DET_MATS[${DET}]}; do
        for SPC in ${SPECIES[@]}; do

            if [ "${DET}" == "ND280" ]; then

                OSCARG=""
                if [ "${TGT}" == "H2O" ]; then
                    OSCARG="--oscillate"
                fi

                bin/fakedatavalid.exe -i flattrees/t2knova.flattree.NEUT.ND280.${TGT}.${SPC}.root \
                                    -F FDSInputs/FakeDataInputs.root \
                                    -M -H config/FakeDataValidConfig.toml \
                                    -W T2KND_to_NOvA \
                                    -a any ${OSCARG} \
                                    -o FDSValid/FakeDataValid_T2KND_to_NOvA_${SPC}.root \
                                    -d ND280/T2KND_to_NOvA/${TGT}/${SPC} &


                bin/fakedatavalid.exe -i flattrees/t2knova.flattree.NEUT.ND280.${TGT}.${SPC}.root \
                                    -F FDSInputs/FakeDataInputs.root \
                                    -M -H config/FakeDataValidConfig.toml \
                                    -a any ${OSCARG} \
                                    -o FDSValid/FakeDataValid_NEUT.root \
                                    -d ND280/T2KNDTune/${TGT}/${SPC} &

                bin/fakedatavalid.exe -i flattrees/t2knova.flattree.GENIE.ND280.${TGT}.${SPC}.root \
                                    -F FDSInputs/FakeDataInputs.root \
                                    -M -H config/FakeDataValidConfig.toml \
                                    -a any ${OSCARG} \
                                    -o FDSValid/FakeDataValid_GENIE.root \
                                    -d ND280/NOvATune/${TGT}/${SPC} &

                wait

                # bin/fakedatavalid.exe -i flattrees/t2knova.flattree.NEUT.ND280.${TGT}.${SPC}.root \
                #                     -F FDSInputs/FakeDataInputs.root \
                #                     --No-Tune -M -H config/FakeDataValidConfig.toml \
                #                     -a any \
                #                     -o FDSValid/FakeDataValid_NEUT_notune.root \
                #                     -d ND280/NEUT/${TGT}/${SPC} &

                # bin/fakedatavalid.exe -i flattrees/t2knova.flattree.GENIE.ND280.${TGT}.${SPC}.root \
                #                     -F FDSInputs/FakeDataInputs.root \
                #                     --No-Tune -M -H config/FakeDataValidConfig.toml \
                #                     -a any \
                #                     -o FDSValid/FakeDataValid_GENIE_notune.root \
                #                     -d ND280/GENIE/${TGT}/${SPC} &
              
                # wait

                # for ELE in ${MAT_ELEMENTS["${TGT}"]}; do

                #     bin/fakedatavalid.exe -i flattrees/t2knova.flattree.NEUT.ND280.${TGT}.${SPC}.root \
                #                         -F FDSInputs/FakeDataInputs.root \
                #                         -M -H config/FakeDataValidConfig.toml \
                #                         -W T2KND_to_NOvA \
                #                         -a ${ELE} \
                #                         -o FDSValid/FakeDataValid_T2KND_to_NOvA.root \
                #                         -d ND280/T2KND_to_NOvA/${ELE}/${SPC} &

                #     bin/fakedatavalid.exe -i flattrees/t2knova.flattree.NEUT.ND280.${TGT}.${SPC}.root \
                #                         -F FDSInputs/FakeDataInputs.root \
                #                         -M -H config/FakeDataValidConfig.toml \
                #                         -a ${ELE} \
                #                         -o FDSValid/FakeDataValid_NEUT.root \
                #                         -d ND280/T2KNDTune/${ELE}/${SPC} &

                #     bin/fakedatavalid.exe -i flattrees/t2knova.flattree.GENIE.ND280.${TGT}.${SPC}.root \
                #                         -F FDSInputs/FakeDataInputs.root \
                #                         -M -H config/FakeDataValidConfig.toml \
                #                         -a ${ELE} \
                #                         -o FDSValid/FakeDataValid_GENIE.root \
                #                         -d ND280/NOvATune/${ELE}/${SPC} &

                #     wait

                #     bin/fakedatavalid.exe -i flattrees/t2knova.flattree.NEUT.ND280.${TGT}.${SPC}.root \
                #                         -F FDSInputs/FakeDataInputs.root \
                #                         --No-Tune -M -H config/FakeDataValidConfig.toml \
                #                         -a ${ELE} \
                #                         -o FDSValid/FakeDataValid_NEUT_notune.root \
                #                         -d ND280/NEUT/${ELE}/${SPC} &

                #     bin/fakedatavalid.exe -i flattrees/t2knova.flattree.GENIE.ND280.${TGT}.${SPC}.root \
                #                         -F FDSInputs/FakeDataInputs.root \
                #                         --No-Tune -M -H config/FakeDataValidConfig.toml \
                #                         -a ${ELE} \
                #                         -o FDSValid/FakeDataValid_GENIE_notune.root \
                #                         -d ND280/GENIE/${ELE}/${SPC} &
                #     wait

                # done
            else

                bin/fakedatavalid.exe -i flattrees/t2knova.flattree.GENIE.NOvAND.${TGT}.${SPC}.root \
                                    -F FDSInputs/FakeDataInputs.root \
                                    -M -H config/FakeDataValidConfig.toml \
                                    -W NOvA_to_T2KND_ptlep \
                                    -a any \
                                    -o FDSValid/FakeDataValid_NOvA_to_T2KND_ptlep.root \
                                    -d NOvAND/NOvA_to_T2KND_ptlep/${TGT}/${SPC} &

                bin/fakedatavalid.exe -i flattrees/t2knova.flattree.NEUT.NOvAND.${TGT}.${SPC}.root \
                                    -F FDSInputs/FakeDataInputs.root \
                                    -M -H config/FakeDataValidConfig.toml \
                                    -a any \
                                    -o FDSValid/FakeDataValid_NEUT.root \
                                    -d NOvAND/T2KNDTune/${TGT}/${SPC} &

                bin/fakedatavalid.exe -i flattrees/t2knova.flattree.GENIE.NOvAND.${TGT}.${SPC}.root \
                                    -F FDSInputs/FakeDataInputs.root \
                                    -M -H config/FakeDataValidConfig.toml \
                                    -a any \
                                    -o FDSValid/FakeDataValid_GENIE.root \
                                    -d NOvAND/NOvATune/${TGT}/${SPC} &

                wait

                bin/fakedatavalid.exe -i flattrees/t2knova.flattree.NEUT.NOvAND.${TGT}.${SPC}.root \
                                    -F FDSInputs/FakeDataInputs.root \
                                    --No-Tune -M -H config/FakeDataValidConfig.toml \
                                    -a any \
                                    -o FDSValid/FakeDataValid_NEUT_notune.root \
                                    -d NOvAND/NEUT/${TGT}/${SPC} &

                bin/fakedatavalid.exe -i flattrees/t2knova.flattree.GENIE.NOvAND.${TGT}.${SPC}.root \
                                    -F FDSInputs/FakeDataInputs.root \
                                    --No-Tune -M -H config/FakeDataValidConfig.toml \
                                    -a any \
                                    -o FDSValid/FakeDataValid_GENIE_notune.root \
                                    -d NOvAND/GENIE/${TGT}/${SPC} &
              
                wait

                for ELE in ${MAT_ELEMENTS["${TGT}"]}; do

                    bin/fakedatavalid.exe -i flattrees/t2knova.flattree.GENIE.NOvAND.${TGT}.${SPC}.root \
                                        -F FDSInputs/FakeDataInputs.root \
                                        -M -H config/FakeDataValidConfig.toml \
                                        -W NOvA_to_T2KND_ptlep \
                                        -a ${ELE} \
                                        -o FDSValid/FakeDataValid_NOvA_to_T2KND_ptlep.root \
                                        -d NOvAND/NOvA_to_T2KND_ptlep/${ELE}/${SPC} &

                    bin/fakedatavalid.exe -i flattrees/t2knova.flattree.NEUT.NOvAND.${TGT}.${SPC}.root \
                                        -F FDSInputs/FakeDataInputs.root \
                                        -M -H config/FakeDataValidConfig.toml \
                                        -a ${ELE} \
                                        -o FDSValid/FakeDataValid_NEUT.root \
                                        -d NOvAND/T2KNDTune/${ELE}/${SPC} &

                    bin/fakedatavalid.exe -i flattrees/t2knova.flattree.GENIE.NOvAND.${TGT}.${SPC}.root \
                                        -F FDSInputs/FakeDataInputs.root \
                                        -M -H config/FakeDataValidConfig.toml \
                                        -a ${ELE} \
                                        -o FDSValid/FakeDataValid_GENIE.root \
                                        -d NOvAND/NOvATune/${ELE}/${SPC} &

                    wait

                    bin/fakedatavalid.exe -i flattrees/t2knova.flattree.NEUT.NOvAND.${TGT}.${SPC}.root \
                                        -F FDSInputs/FakeDataInputs.root \
                                        --No-Tune -M -H config/FakeDataValidConfig.toml \
                                        -a ${ELE} \
                                        -o FDSValid/FakeDataValid_NEUT_notune.root \
                                        -d NOvAND/NEUT/${ELE}/${SPC} &

                    bin/fakedatavalid.exe -i flattrees/t2knova.flattree.GENIE.NOvAND.${TGT}.${SPC}.root \
                                        -F FDSInputs/FakeDataInputs.root \
                                        --No-Tune -M -H config/FakeDataValidConfig.toml \
                                        -a ${ELE} \
                                        -o FDSValid/FakeDataValid_GENIE_notune.root \
                                        -d NOvAND/GENIE/${ELE}/${SPC} &
                    wait

                done

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