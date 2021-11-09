#!/bin/bash

rm -rf FDSValid
mkdir -p FDSValid

set -e
set -x

GENERATORS=( GENIE NEUT )
SPECIES=( numu numub nue nueb )
SPECIES=( numu numub )
DETECTORS=( NOvAND ND280 )
# DETECTORS=( ND280 )

declare -A DET_MATS
DET_MATS["NOvAND"]="CH"
DET_MATS["ND280"]="H2O CH"
DET_MATS["ND280"]="CH"

declare -A MAT_ELEMENTS
MAT_ELEMENTS["CH"]="C H"
# MAT_ELEMENTS["CH"]="C"
MAT_ELEMENTS["H2O"]="H O"
# MAT_ELEMENTS["CH"]=""
# MAT_ELEMENTS["H2O"]=""

for DET in ${DETECTORS}; do
    for TGT in ${DET_MATS[${DET}]}; do
      for SPC in ${SPECIES[@]}; do

        for WGHTSCHEME in T2KND_to_NOvA \
                          T2KND_to_NOvA_Enu \
                          T2KND_to_NOvA_Q2 \
                          T2KND_to_NOvA_EnuKludge; do

            bin/fakedatavalid.exe -i flattrees/t2knova.flattree.NEUT.ND280.${TGT}.${SPC}.root \
                                -F FDSInputs/FakeDataInputs.root \
                                -H config/FakeDataValidConfig.toml \
                                -W ${WGHTSCHEME} \
                                -a any \
                                -o FDSValid/FakeDataValid_${WGHTSCHEME}.root \
                                -d ND280/${WGHTSCHEME}/${TGT}/${SPC} &
        done

        wait

        bin/fakedatavalid.exe -i flattrees/t2knova.flattree.NEUT.ND280.${TGT}.${SPC}.root \
                            -F FDSInputs/FakeDataInputs.root \
                            -H config/FakeDataValidConfig.toml \
                            -a any \
                            -o FDSValid/FakeDataValid_NEUT.root \
                            -d ND280/T2KNDTune/${TGT}/${SPC} &

        bin/fakedatavalid.exe -i flattrees/t2knova.flattree.GENIE.ND280.${TGT}.${SPC}.root \
                            -F FDSInputs/FakeDataInputs.root \
                            -H config/FakeDataValidConfig.toml \
                            -a any \
                            -o FDSValid/FakeDataValid_GENIE.root \
                            -d ND280/NOvATune/${TGT}/${SPC} &
      
        wait

        for ELE in ${MAT_ELEMENTS["${TGT}"]}; do

            for WGHTSCHEME in T2KND_to_NOvA \
                              T2KND_to_NOvA_Enu \
                              T2KND_to_NOvA_Q2 \
                              T2KND_to_NOvA_EnuKludge; do

                bin/fakedatavalid.exe -i flattrees/t2knova.flattree.NEUT.ND280.${TGT}.${SPC}.root \
                                    -F FDSInputs/FakeDataInputs.root \
                                    -H config/FakeDataValidConfig.toml \
                                    -W ${WGHTSCHEME} \
                                    -a ${ELE} \
                                    -o FDSValid/FakeDataValid_${WGHTSCHEME}.root \
                                    -d ND280/${WGHTSCHEME}/${ELE}/${SPC} &
            done

            wait

            bin/fakedatavalid.exe -i flattrees/t2knova.flattree.NEUT.ND280.${TGT}.${SPC}.root \
                                -F FDSInputs/FakeDataInputs.root \
                                -H config/FakeDataValidConfig.toml \
                                -a ${ELE} \
                                -o FDSValid/FakeDataValid_NEUT.root \
                                -d ND280/T2KNDTune/${ELE}/${SPC} &

            bin/fakedatavalid.exe -i flattrees/t2knova.flattree.GENIE.ND280.${TGT}.${SPC}.root \
                                -F FDSInputs/FakeDataInputs.root \
                                -H config/FakeDataValidConfig.toml \
                                -a ${ELE} \
                                -o FDSValid/FakeDataValid_GENIE.root \
                                -d ND280/NOvATune/${ELE}/${SPC} &
            wait

        done

        for WGHTSCHEME in NOvA_to_T2KND_plep \
                          NOvA_to_T2KND_Q2 \
                          NOvA_to_T2KND_ptlep; do

            bin/fakedatavalid.exe -i flattrees/t2knova.flattree.GENIE.NOvAND.${TGT}.${SPC}.root \
                                -F FDSInputs/FakeDataInputs.root \
                                -H config/FakeDataValidConfig.toml \
                                -W ${WGHTSCHEME} \
                                -a any \
                                -o FDSValid/FakeDataValid_${WGHTSCHEME}.root \
                                -d NOvAND/${WGHTSCHEME}/${TGT}/${SPC} &
        done

        wait

        bin/fakedatavalid.exe -i flattrees/t2knova.flattree.NEUT.NOvAND.${TGT}.${SPC}.root \
                            -F FDSInputs/FakeDataInputs.root \
                            -H config/FakeDataValidConfig.toml \
                            -a any \
                            -o FDSValid/FakeDataValid_NEUT.root \
                            -d NOvAND/T2KNDTune/${TGT}/${SPC} &

        bin/fakedatavalid.exe -i flattrees/t2knova.flattree.GENIE.NOvAND.${TGT}.${SPC}.root \
                            -F FDSInputs/FakeDataInputs.root \
                            -H config/FakeDataValidConfig.toml \
                            -a any \
                            -o FDSValid/FakeDataValid_GENIE.root \
                            -d NOvAND/NOvATune/${TGT}/${SPC} &
      
        wait

        for ELE in ${MAT_ELEMENTS["${TGT}"]}; do

            for WGHTSCHEME in NOvA_to_T2KND_plep \
                          NOvA_to_T2KND_Q2 \
                          NOvA_to_T2KND_ptlep; do

                bin/fakedatavalid.exe -i flattrees/t2knova.flattree.GENIE.NOvAND.${TGT}.${SPC}.root \
                                    -F FDSInputs/FakeDataInputs.root \
                                    -H config/FakeDataValidConfig.toml \
                                    -W ${WGHTSCHEME} \
                                    -a ${ELE} \
                                    -o FDSValid/FakeDataValid_${WGHTSCHEME}.root \
                                    -d NOvAND/${WGHTSCHEME}/${ELE}/${SPC} &
            done

            wait

            bin/fakedatavalid.exe -i flattrees/t2knova.flattree.NEUT.NOvAND.${TGT}.${SPC}.root \
                                -F FDSInputs/FakeDataInputs.root \
                                -H config/FakeDataValidConfig.toml \
                                -a ${ELE} \
                                -o FDSValid/FakeDataValid_NEUT.root \
                                -d NOvAND/T2KNDTune/${ELE}/${SPC} &

            bin/fakedatavalid.exe -i flattrees/t2knova.flattree.GENIE.NOvAND.${TGT}.${SPC}.root \
                                -F FDSInputs/FakeDataInputs.root \
                                -H config/FakeDataValidConfig.toml \
                                -a ${ELE} \
                                -o FDSValid/FakeDataValid_GENIE.root \
                                -d NOvAND/NOvATune/${ELE}/${SPC} &
            wait

        done
      done
    done
done
 
hadd -j 2 FDSValid/FakeDataValid.root \
          FDSValid/FakeDataValid_*.root


if [ ! -e FDSValid/FakeDataValid.root ]; then
  echo "Didn't produce FakeDataValid.root. Failed"
  exit 1
fi