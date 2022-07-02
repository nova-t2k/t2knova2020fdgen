#!/bin/bash

if [ "${1}" = "force" ]; then
    rm -rf FDSValid
fi
mkdir -p FDSValid

set -e
set -x

# SPECIES=( numu numub nue nueb )
SPECIES=( numu numub )
# SPECIES=( numu )

DETECTORS=( NOvAND ND280 )

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

DO_MAIN=1
DO_MAT_ELEMENTS=0
DO_ND280=1
DO_NOvAND=0

for DET in ${DETECTORS[@]}; do
    for TGT in ${DET_MATS[${DET}]}; do
        for SPC in ${SPECIES[@]}; do

            if [ "${DET}" == "ND280" ] && [ "${DO_ND280}" == "1" ]; then

                if [ "${DO_MAIN}" == "1" ]; then
                    bin/fakedatavalid.exe -i flattrees/t2knova.flattree.NEUT.ND280.${TGT}.BANFF_PRE.${SPC}.root \
                                        -F FDSInputs/FakeDataInputs.root \
                                        -H config/FakeDataValidConfig_ND280.toml \
                                        -a any --No-Tune \
                                        -o FDSValid/FakeDataValid_Generated.root \
                                        -d ND280/NEUT/Generated/${TGT}/${SPC} &

                    bin/fakedatavalid.exe -i flattrees/t2knova.flattree.NEUT.ND280.${TGT}.BANFF_PRE.${SPC}.root \
                                        -F FDSInputs/FakeDataInputs.root \
                                        -H config/FakeDataValidConfig_ND280.toml \
                                        -a any \
                                        -o FDSValid/FakeDataValid_BANFF_PRE.root \
                                        -d ND280/NEUT/BANFF_PRE/${TGT}/${SPC} &

                    bin/fakedatavalid.exe -i flattrees/t2knova.flattree.NEUT.ND280.${TGT}.BANFF_POST.${SPC}.root \
                                        -F FDSInputs/FakeDataInputs.root \
                                        -H config/FakeDataValidConfig_ND280.toml \
                                        -a any \
                                        -o FDSValid/FakeDataValid_BANFF_POST.root \
                                        -d ND280/NEUT/BANFF_POST/${TGT}/${SPC} &

                    bin/fakedatavalid.exe -i flattrees/t2knova.flattree.NEUT.ND280.${TGT}.BANFF_POST.${SPC}.root \
                                        -F FDSInputs/FakeDataInputs.root \
                                        -H config/FakeDataValidConfig_ND280.toml \
                                        -a any --No-Tune \
                                        -W T2KND_to_NOvA \
                                        -o FDSValid/FakeDataValid_Generated_to_2020.root \
                                        -d ND280/NEUT/Generated_to_2020/${TGT}/${SPC} &

                    bin/fakedatavalid.exe -i flattrees/t2knova.flattree.GENIE.ND280.${TGT}.2020.${SPC}.root \
                                        -F FDSInputs/FakeDataInputs.root \
                                        -H config/FakeDataValidConfig_ND280.toml \
                                        -a any \
                                        -o FDSValid/FakeDataValid_GENIE.root \
                                        -d ND280/GENIE/2020/${TGT}/${SPC} &

                    wait
                fi

                if [ "${DO_MAT_ELEMENTS}" == "1" ]; then 
                    for ELE in ${MAT_ELEMENTS["${TGT}"]}; do
                        bin/fakedatavalid.exe -i flattrees/t2knova.flattree.NEUT.ND280.${TGT}.BANFF_PRE.${SPC}.root \
                                            -F FDSInputs/FakeDataInputs.root \
                                            -H config/FakeDataValidConfig_ND280.toml \
                                            -a ${ELE} --No-Tune \
                                            -o FDSValid/FakeDataValid_Generated.root \
                                            -d ND280/NEUT/Generated/${ELE}/${SPC} &

                        bin/fakedatavalid.exe -i flattrees/t2knova.flattree.NEUT.ND280.${TGT}.BANFF_PRE.${SPC}.root \
                                            -F FDSInputs/FakeDataInputs.root \
                                            -H config/FakeDataValidConfig_ND280.toml \
                                            -a ${ELE} \
                                            -o FDSValid/FakeDataValid_BANFF_PRE.root \
                                            -d ND280/NEUT/BANFF_PRE/${ELE}/${SPC} &

                        bin/fakedatavalid.exe -i flattrees/t2knova.flattree.NEUT.ND280.${TGT}.BANFF_POST.${SPC}.root \
                                            -F FDSInputs/FakeDataInputs.root \
                                            -H config/FakeDataValidConfig_ND280.toml \
                                            -a ${ELE} \
                                            -o FDSValid/FakeDataValid_BANFF_POST.root \
                                            -d ND280/NEUT/BANFF_POST/${ELE}/${SPC} &

                        bin/fakedatavalid.exe -i flattrees/t2knova.flattree.NEUT.ND280.${TGT}.BANFF_POST.${SPC}.root \
                                            -F FDSInputs/FakeDataInputs.root \
                                            -H config/FakeDataValidConfig_ND280.toml \
                                            -a ${ELE} --No-Tune \
                                            -W T2KND_to_NOvA \
                                            -o FDSValid/FakeDataValid_Generated_to_2020.root \
                                            -d ND280/NEUT/Generated_to_2020/${ELE}/${SPC} &

                        bin/fakedatavalid.exe -i flattrees/t2knova.flattree.GENIE.ND280.${TGT}.2020.${SPC}.root \
                                            -F FDSInputs/FakeDataInputs.root \
                                            -H config/FakeDataValidConfig_ND280.toml \
                                            -a ${ELE} \
                                            -o FDSValid/FakeDataValid_GENIE.root \
                                            -d ND280/GENIE/2020/${ELE}/${SPC} &

                        wait
                    done
                fi

            elif [ "${DET}" == "NOvAND" ] && [ "${DO_NOvAND}" == "1" ]; then

                if [ "${DO_MAIN}" == "1" ]; then
                    bin/fakedatavalid.exe -i flattrees/t2knova.flattree.GENIE.NOvAND.${TGT}.2020.${SPC}.root \
                                        -F FDSInputs/FakeDataInputs.root \
                                        -H config/FakeDataValidConfig_NOvAND.toml \
                                        -a any --No-Tune \
                                        -o FDSValid/FakeDataValid_Generated.root \
                                        -d NOvAND/GENIE/Generated/${TGT}/${SPC} &

                    bin/fakedatavalid.exe -i flattrees/t2knova.flattree.GENIE.NOvAND.${TGT}.2020.${SPC}.root \
                                        -F FDSInputs/FakeDataInputs.root \
                                        -H config/FakeDataValidConfig_NOvAND.toml \
                                        -a any \
                                        -o FDSValid/FakeDataValid_2020.root \
                                        -d NOvAND/GENIE/2020/${TGT}/${SPC} &

                    bin/fakedatavalid.exe -i flattrees/t2knova.flattree.NEUT.NOvAND.${TGT}.BANFF_PRE.${SPC}.root \
                                        -F FDSInputs/FakeDataInputs.root \
                                        -H config/FakeDataValidConfig_NOvAND.toml \
                                        -a any \
                                        -o FDSValid/FakeDataValid_BANFF_PRE.root \
                                        -d NOvAND/NEUT/BANFF_PRE/${TGT}/${SPC} &

                    bin/fakedatavalid.exe -i flattrees/t2knova.flattree.NEUT.NOvAND.${TGT}.BANFF_POST.${SPC}.root \
                                        -F FDSInputs/FakeDataInputs.root \
                                        -H config/FakeDataValidConfig_NOvAND.toml \
                                        -a any \
                                        -o FDSValid/FakeDataValid_BANFF_POST.root \
                                        -d NOvAND/NEUT/BANFF_POST/${TGT}/${SPC} &

                    wait

                    bin/fakedatavalid.exe -i flattrees/t2knova.flattree.GENIE.NOvAND.${TGT}.2020.${SPC}.root \
                                        -F FDSInputs/FakeDataInputs.root \
                                        -H config/FakeDataValidConfig_NOvAND.toml \
                                        -a any --No-Tune -W NOvA_to_T2KND_ptlep \
                                        -o FDSValid/FakeDataValid_Generated_to_BANFF_POST.root \
                                        -d NOvAND/GENIE/Generated_to_BANFF_POST/${TGT}/${SPC} &

                    bin/fakedatavalid.exe -i flattrees/t2knova.flattree.GENIE.NOvAND.${TGT}.2020.${SPC}.root \
                                        -F FDSInputs/FakeDataInputs.root \
                                        -H config/FakeDataValidConfig_NOvAND.toml \
                                        -a any --No-Tune -W NOvA_to_T2KPre_ptlep \
                                        -o FDSValid/FakeDataValid_Generated_to_BANFF_PRE.root \
                                        -d NOvAND/GENIE/Generated_to_BANFF_PRE/${TGT}/${SPC} &

                    bin/fakedatavalid.exe -i flattrees/t2knova.flattree.GENIE.NOvAND.${TGT}.2020.${SPC}.root \
                                        -F FDSInputs/FakeDataInputs.root \
                                        -H config/FakeDataValidConfig_NOvAND.toml \
                                        -a any --No-Tune -W NOvA_to_T2KMnv1Pi_ptlep \
                                        -o FDSValid/FakeDataValid_Generated_to_Mnv1Pi.root \
                                        -d NOvAND/GENIE/Generated_to_Mnv1Pi/${TGT}/${SPC} &


                    bin/fakedatavalid.exe -i flattrees/t2knova.flattree.GENIE.NOvAND.${TGT}.2020.${SPC}.root \
                                        -F FDSInputs/FakeDataInputs.root \
                                        -H config/FakeDataValidConfig_NOvAND.toml \
                                        -a any --No-Tune -W NOvA_to_T2KNonQE_ptlep \
                                        -o FDSValid/FakeDataValid_Generated_to_NonQE.root \
                                        -d NOvAND/GENIE/Generated_to_NonQE/${TGT}/${SPC} &

                    wait
                fi


                if [ "${DO_MAT_ELEMENTS}" == "1" ]; then

                    for ELE in ${MAT_ELEMENTS["${TGT}"]}; do

                        bin/fakedatavalid.exe -i flattrees/t2knova.flattree.GENIE.NOvAND.${TGT}.2020.${SPC}.root \
                                            -F FDSInputs/FakeDataInputs.root \
                                            -H config/FakeDataValidConfig_NOvAND.toml \
                                            -a ${ELE} --No-Tune \
                                            -o FDSValid/FakeDataValid_Generated.root \
                                            -d NOvAND/GENIE/Generated/${ELE}/${SPC} &

                        bin/fakedatavalid.exe -i flattrees/t2knova.flattree.GENIE.NOvAND.${TGT}.2020.${SPC}.root \
                                            -F FDSInputs/FakeDataInputs.root \
                                            -H config/FakeDataValidConfig_NOvAND.toml \
                                            -a ${ELE} \
                                            -o FDSValid/FakeDataValid_2020.root \
                                            -d NOvAND/GENIE/2020/${ELE}/${SPC} &

                        bin/fakedatavalid.exe -i flattrees/t2knova.flattree.NEUT.NOvAND.${TGT}.BANFF_PRE.${SPC}.root \
                                            -F FDSInputs/FakeDataInputs.root \
                                            -H config/FakeDataValidConfig_NOvAND.toml \
                                            -a ${ELE} \
                                            -o FDSValid/FakeDataValid_BANFF_PRE.root \
                                            -d NOvAND/NEUT/BANFF_PRE/${ELE}/${SPC} &

                        bin/fakedatavalid.exe -i flattrees/t2knova.flattree.NEUT.NOvAND.${TGT}.BANFF_POST.${SPC}.root \
                                            -F FDSInputs/FakeDataInputs.root \
                                            -H config/FakeDataValidConfig_NOvAND.toml \
                                            -a ${ELE} \
                                            -o FDSValid/FakeDataValid_BANFF_POST.root \
                                            -d NOvAND/NEUT/BANFF_POST/${ELE}/${SPC} &

                        wait

                        bin/fakedatavalid.exe -i flattrees/t2knova.flattree.GENIE.NOvAND.${TGT}.2020.${SPC}.root \
                                            -F FDSInputs/FakeDataInputs.root \
                                            -H config/FakeDataValidConfig_NOvAND.toml \
                                            -a ${ELE} --No-Tune -W NOvA_to_T2KND_ptlep \
                                            -o FDSValid/FakeDataValid_Generated_to_BANFF_POST.root \
                                            -d NOvAND/GENIE/Generated_to_BANFF_POST/${ELE}/${SPC} &

                        bin/fakedatavalid.exe -i flattrees/t2knova.flattree.GENIE.NOvAND.${TGT}.2020.${SPC}.root \
                                            -F FDSInputs/FakeDataInputs.root \
                                            -H config/FakeDataValidConfig_NOvAND.toml \
                                            -a ${ELE} --No-Tune -W NOvA_to_T2KPre_ptlep \
                                            -o FDSValid/FakeDataValid_Generated_to_BANFF_PRE.root \
                                            -d NOvAND/GENIE/Generated_to_BANFF_PRE/${ELE}/${SPC} &

                        bin/fakedatavalid.exe -i flattrees/t2knova.flattree.GENIE.NOvAND.${TGT}.2020.${SPC}.root \
                                            -F FDSInputs/FakeDataInputs.root \
                                            -H config/FakeDataValidConfig_NOvAND.toml \
                                            -a ${ELE} --No-Tune -W NOvA_to_T2KMnv1Pi_ptlep \
                                            -o FDSValid/FakeDataValid_Generated_to_Mnv1Pi.root \
                                            -d NOvAND/GENIE/Generated_to_Mnv1Pi/${ELE}/${SPC} &


                        bin/fakedatavalid.exe -i flattrees/t2knova.flattree.GENIE.NOvAND.${TGT}.2020.${SPC}.root \
                                            -F FDSInputs/FakeDataInputs.root \
                                            -H config/FakeDataValidConfig_NOvAND.toml \
                                            -a ${ELE} --No-Tune -W NOvA_to_T2KNonQE_ptlep \
                                            -o FDSValid/FakeDataValid_Generated_to_NonQE.root \
                                            -d NOvAND/GENIE/Generated_to_NonQE/${ELE}/${SPC} &

                        wait

                    done
                fi
            fi
        done
    done
done

rm -f FDSValid/FakeDataValid.root
hadd -j 2 FDSValid/FakeDataValid.root \
          FDSValid/FakeDataValid_*.root

if [ ! -e FDSValid/FakeDataValid.root ]; then
  echo "Didn't produce FakeDataValid.root. Failed"
  exit 1
fi