#!/bin/bash

g++ -I ../../../ fakedatavalid.C -O3 $(root-config --cflags --glibs) -o fakedatavalid.exe
if [[ ! "$?" == "0" ]]; then
        exit;
fi

declare -A TGTEL
TGTEL["NOvAND_CH"]="C H"
TGTEL["ND280_CH"]="C"
TGTEL["ND280_H2O"]="H O"

GENERATORS=( GENIE NEUT )
SPECIES=( numu numub nue nueb )
SPECIES=( numu numub )
DETECTORS=( NOvAND ND280 )
DETECTORS=( ND280 )

declare -A DET_MATS
DET_MATS["NOvAND"]="CH"
DET_MATS["ND280"]="H2O CH"
DET_MATS["ND280"]="CH"

declare -A MAT_ELEMENTS
MAT_ELEMENTS["CH"]= "C H"
MAT_ELEMENTS["CH"]= "C"
MAT_ELEMENTS["H20"]= "H O"

rm -f FakeDataValid.root
for TGT in ${DET_MATS["ND280"]}; do
  for SPC in ${SPECIES[@]}; do

    for WGHTSCHEME in T2KND_To_NOvA \
                      T2KND_To_NOvA_Enu \
                      T2KND_To_NOvA_Q2; do

        ./fakedatavalid.exe -i flattrees/t2knova.flattree.NEUT.ND280.${TGT}.${SPC}.root \
                            -H FakeDataValidConfig.toml \
                            -W ${WGHTSCHEME} \
                            -a any \
                            -o FakeDataValid.root \
                            -d ND280/${WGHTSCHEME}/${TGT}/${SPC}
    done

    ./fakedatavalid.exe -i flattrees/t2knova.flattree.NEUT.ND280.${TGT}.${SPC}.root \
                        -H FakeDataValidConfig.toml \
                        -a any \
                        -o FakeDataValid.root \
                        -d ND280/T2KNDTune/${TGT}/${SPC}

    ./fakedatavalid.exe -i flattrees/t2knova.flattree.GENIE.ND280.${TGT}.${SPC}.root \
                        -H FakeDataValidConfig.toml \
                        -a any \
                        -o FakeDataValid.root \
                        -d ND280/NOvATune/${TGT}/${SPC}
  
    for ELE in ${MAT_ELEMENTS["${TGT}"]}; do

        for WGHTSCHEME in T2KND_To_NOvA \
                  T2KND_To_NOvA_Enu \
                  T2KND_To_NOvA_Q2; do

            ./fakedatavalid.exe -i flattrees/t2knova.flattree.NEUT.ND280.${TGT}.${SPC}.root \
                                -H FakeDataValidConfig.toml \
                                -W ${WGHTSCHEME} \
                                -a ${ELE} \
                                -o FakeDataValid.root \
                                -d ND280/${WGHTSCHEME}/${TGT}/${ELE}/${SPC}
        done

        ./fakedatavalid.exe -i flattrees/t2knova.flattree.NEUT.ND280.${TGT}.${SPC}.root \
                            -H FakeDataValidConfig.toml \
                            -a ${ELE} \
                            -o FakeDataValid.root \
                            -d ND280/T2KNDTune/${TGT}/${ELE}/${SPC}

        ./fakedatavalid.exe -i flattrees/t2knova.flattree.GENIE.ND280.${TGT}.${SPC}.root \
                            -H FakeDataValidConfig.toml \
                            -a ${ELE} \
                            -o FakeDataValid.root \
                            -d ND280/NOvATune/${TGT}/${ELE}/${SPC}

    done
  done
done
 
