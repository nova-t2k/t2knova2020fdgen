#!/bin/bash

OUTDIR=FDSValid

if [ "${1}" = "force" ]; then
    rm -rf ${OUTDIR}
fi
mkdir -p ${OUTDIR}

set -e
set -x

export T2KNOVA_INPUTS=$(readlink -f inputs)

SPECIES=( numu numub nue nueb )
# SPECIES=( numu numub )
# SPECIES=( numu )

DETECTORS=( NOvAND ND280 )
# DETECTORS=( ND280 )

declare -A DET_MATS
DET_MATS["NOvAND"]="CH"
DET_MATS["ND280"]="H2O CH"
DET_MATS["ND280"]="CH"
# DET_MATS["ND280"]="H2O"

declare -A MAT_ELEMENTS
MAT_ELEMENTS["CH"]="C H"
# MAT_ELEMENTS["CH"]="C"
MAT_ELEMENTS["H2O"]="H O"
# MAT_ELEMENTS["CH"]=""
# MAT_ELEMENTS["H2O"]=""

DO_MAIN=1
DO_MAT_ELEMENTS=1
DO_ND280=1
DO_NOvAND=1

declare -a NUCS

for DET in ${DETECTORS[@]}; do
    for TGT in ${DET_MATS[${DET}]}; do
        for SPC in ${SPECIES[@]}; do

            if [ "${DET}" == "ND280" ] && [ "${DO_ND280}" == "1" ]; then

              for NUC in any ${MAT_ELEMENTS[${TGT}]}; do

                AFLAG="any"
                TGTDIR=${TGT}

                if [ "${NUC}" == "any" ] && [ "${DO_MAIN}" != "1" ]; then
                  continue
                fi

                if [ "${NUC}" != "any" ]; then
                  if [ "${DO_MAT_ELEMENTS}" != "1" ]; then
                    continue
                  fi
                  AFLAG=${NUC}
                  TGTDIR=${NUC}
                fi

                  #########Generated
                  bin/fakedatavalid.exe -i flattrees/t2knova.flattree.NEUT.ND280.${TGT}.BANFF_PRE.${SPC}.root \
                                      -H config/FakeDataValidConfig_ND280.toml \
                                      -a ${AFLAG} -T Generated \
                                      -o ${OUTDIR}/FakeDataValid_Generated.root \
                                      -d ND280/NEUT/Generated/${TGTDIR}/${SPC} &

                  #########FDS Targets
                  #BANFF_PRE
                  bin/fakedatavalid.exe -i flattrees/t2knova.flattree.NEUT.ND280.${TGT}.BANFF_PRE.${SPC}.root \
                                      -H config/FakeDataValidConfig_ND280.toml \
                                      -a ${AFLAG} -T NDTuned \
                                      -o ${OUTDIR}/FakeDataValid_BANFF_PRE.root \
                                      -d ND280/NEUT/BANFF_PRE/${TGTDIR}/${SPC} &

                  #BANFF_POST
                  bin/fakedatavalid.exe -i flattrees/t2knova.flattree.NEUT.ND280.${TGT}.BANFF_POST.${SPC}.root \
                                      -H config/FakeDataValidConfig_ND280.toml \
                                      -a ${AFLAG} -T NDTuned \
                                      -o ${OUTDIR}/FakeDataValid_BANFF_POST.root \
                                      -d ND280/NEUT/BANFF_POST/${TGTDIR}/${SPC} &
                  #BANFF_POST + NonQE
                  bin/fakedatavalid.exe -i flattrees/t2knova.flattree.NEUT.ND280.${TGT}.BANFF_POST.${SPC}.root \
                                      -H config/FakeDataValidConfig_ND280.toml \
                                      -a ${AFLAG} -T NonQE \
                                      -o ${OUTDIR}/FakeDataValid_NonQE.root \
                                      -d ND280/NEUT/NonQE/${TGTDIR}/${SPC} &
                  #BANFF_PRE + Mnv1Pi
                  bin/fakedatavalid.exe -i flattrees/t2knova.flattree.NEUT.ND280.${TGT}.BANFF_PRE.${SPC}.root \
                                      -H config/FakeDataValidConfig_ND280.toml \
                                      -a ${AFLAG} -T Mnv1Pi \
                                      -o ${OUTDIR}/FakeDataValid_Mnv1Pi.root \
                                      -d ND280/NEUT/Mnv1Pi/${TGTDIR}/${SPC} &

                  #GENIE 2020
                  bin/fakedatavalid.exe -i flattrees/t2knova.flattree.GENIE.ND280.${TGT}.2020.${SPC}.root \
                                      -H config/FakeDataValidConfig_ND280.toml \
                                      -a ${AFLAG} -T NDTuned \
                                      -o ${OUTDIR}/FakeDataValid_GENIE.root \
                                      -d ND280/GENIE/2020/${TGTDIR}/${SPC} &

                  #########FDS ReWeights
                  #Base Tune + RW to NOvA2020
                  bin/fakedatavalid.exe -i flattrees/t2knova.flattree.NEUT.ND280.${TGT}.BANFF_PRE.${SPC}.root \
                                      --From BANFFPre -F FDSInputs/FakeDataInputs_FromTuned_BANFFPre.root \
                                      -H config/FakeDataValidConfig_ND280.toml \
                                      -a ${AFLAG} -T NDTuned \
                                      -W T2KND_to_NOvA \
                                      -o ${OUTDIR}/FakeDataValid_ReWeighted_to_2020.root \
                                      -d ND280/NEUT/ReWeighted_to_2020/${TGTDIR}/${SPC} &

                  #BANFFPre to T2KNonQE
                  bin/fakedatavalid.exe -i flattrees/t2knova.flattree.NEUT.ND280.${TGT}.BANFF_PRE.${SPC}.root \
                                      --From BANFFPre -F FDSInputs/FakeDataInputs_FromTuned_BANFFPre.root \
                                      -H config/FakeDataValidConfig_ND280.toml \
                                      -a ${AFLAG} -T NDTuned \
                                      -W kT2KND_to_T2KNonQE \
                                      -o ${OUTDIR}/FakeDataValid_ReWeighted_to_NonQE.root \
                                      -d ND280/NEUT/ReWeighted_to_NonQE/${TGTDIR}/${SPC} &

                  #BANFFPre to T2KMnv1Pi
                  bin/fakedatavalid.exe -i flattrees/t2knova.flattree.NEUT.ND280.${TGT}.BANFF_PRE.${SPC}.root \
                                      --From BANFFPre -F FDSInputs/FakeDataInputs_FromTuned_BANFFPre.root \
                                      -H config/FakeDataValidConfig_ND280.toml \
                                      -a ${AFLAG} -T NDTuned \
                                      -W kT2KND_to_T2KMnv1Pi \
                                      -o ${OUTDIR}/FakeDataValid_ReWeighted_to_Mnv1Pi.root \
                                      -d ND280/NEUT/ReWeighted_to_Mnv1Pi/${TGTDIR}/${SPC} &

                  wait
              done

            elif [ "${DET}" == "NOvAND" ] && [ "${DO_NOvAND}" == "1" ]; then

                for NUC in any ${MAT_ELEMENTS[${TGT}]}; do

                  AFLAG="any"
                  TGTDIR=${TGT}

                  if [ "${NUC}" == "any" ] && [ "${DO_MAIN}" != "1" ]; then
                    continue
                  fi

                  if [ "${NUC}" != "any" ]
                    if [ "${DO_MAT_ELEMENTS}" != "1" ]; then
                      continue
                    fi
                    AFLAG=${NUC}
                    TGTDIR=${NUC}
                  fi
                  bin/fakedatavalid.exe -i flattrees/t2knova.flattree.GENIE.NOvAND.${TGT}.2020.${SPC}.root \
                                      -H config/FakeDataValidConfig_NOvAND.toml \
                                      -a ${AFLAG} -T Generated \
                                      -o ${OUTDIR}/FakeDataValid_Generated.root \
                                      -d NOvAND/GENIE/Generated/${TGTDIR}/${SPC} &

                  bin/fakedatavalid.exe -i flattrees/t2knova.flattree.GENIE.NOvAND.${TGT}.2020.${SPC}.root \
                                      -H config/FakeDataValidConfig_NOvAND.toml \
                                      -a ${AFLAG} -T NDTuned \
                                      -o ${OUTDIR}/FakeDataValid_2020.root \
                                      -d NOvAND/GENIE/2020/${TGTDIR}/${SPC} &

                  bin/fakedatavalid.exe -i flattrees/t2knova.flattree.NEUT.NOvAND.${TGT}.BANFF_PRE.${SPC}.root \
                                      -H config/FakeDataValidConfig_NOvAND.toml \
                                      -a ${AFLAG} -T Generated \
                                      -o ${OUTDIR}/FakeDataValid_NEUT_Generated.root \
                                      -d NOvAND/NEUT/Generated/${TGTDIR}/${SPC} &

                  bin/fakedatavalid.exe -i flattrees/t2knova.flattree.NEUT.NOvAND.${TGT}.BANFF_PRE.${SPC}.root \
                                      -H config/FakeDataValidConfig_NOvAND.toml \
                                      -a ${AFLAG} -T NDTuned \
                                      -o ${OUTDIR}/FakeDataValid_BANFF_PRE.root \
                                      -d NOvAND/NEUT/BANFF_PRE/${TGTDIR}/${SPC} &

                  bin/fakedatavalid.exe -i flattrees/t2knova.flattree.NEUT.NOvAND.${TGT}.BANFF_POST.${SPC}.root \
                                      -H config/FakeDataValidConfig_NOvAND.toml \
                                      -a ${AFLAG} -T NDTuned \
                                      -o ${OUTDIR}/FakeDataValid_BANFF_POST.root \
                                      -d NOvAND/NEUT/BANFF_POST/${TGTDIR}/${SPC} &

                  bin/fakedatavalid.exe -i flattrees/t2knova.flattree.NEUT.NOvAND.${TGT}.BANFF_PRE.${SPC}.root \
                                      -H config/FakeDataValidConfig_NOvAND.toml \
                                      -a ${AFLAG} -T Mnv1Pi \
                                      -o ${OUTDIR}/FakeDataValid_Mnv1Pi.root \
                                      -d NOvAND/NEUT/Mnv1Pi/${TGTDIR}/${SPC} &

                  bin/fakedatavalid.exe -i flattrees/t2knova.flattree.NEUT.NOvAND.${TGT}.BANFF_POST.${SPC}.root \
                                      -H config/FakeDataValidConfig_NOvAND.toml \
                                      -a ${AFLAG} -T NonQE \
                                      -o ${OUTDIR}/FakeDataValid_NonQE.root \
                                      -d NOvAND/NEUT/NonQE/${TGTDIR}/${SPC} &

                  bin/fakedatavalid.exe -i flattrees/t2knova.flattree.GENIE.NOvAND.${TGT}.2020.${SPC}.root \
                                      --From Generated -F FDSInputs/FakeDataInputs_FromGenerated.root \
                                      -H config/FakeDataValidConfig_NOvAND.toml \
                                      -a ${AFLAG} -T Generated -W NOvA_to_T2KND_ptlep \
                                      -o ${OUTDIR}/FakeDataValid_ReWeighted_to_BANFF_POST.root \
                                      -d NOvAND/GENIE/ReWeighted_to_BANFF_POST/${TGTDIR}/${SPC} &

                  bin/fakedatavalid.exe -i flattrees/t2knova.flattree.GENIE.NOvAND.${TGT}.2020.${SPC}.root \
                                      --From Generated -F FDSInputs/FakeDataInputs_FromGenerated.root \
                                      -H config/FakeDataValidConfig_NOvAND.toml \
                                      -a ${AFLAG} -T Generated -W NOvA_to_T2KPre_ptlep \
                                      -o ${OUTDIR}/FakeDataValid_ReWeighted_to_BANFF_PRE.root \
                                      -d NOvAND/GENIE/ReWeighted_to_BANFF_PRE/${TGTDIR}/${SPC} &

                  bin/fakedatavalid.exe -i flattrees/t2knova.flattree.GENIE.NOvAND.${TGT}.2020.${SPC}.root \
                                      --From Generated -F FDSInputs/FakeDataInputs_FromGenerated.root \
                                      -H config/FakeDataValidConfig_NOvAND.toml \
                                      -a ${AFLAG} -T Generated -W NOvA_to_T2KMnv1Pi_ptlep \
                                      -o ${OUTDIR}/FakeDataValid_ReWeighted_to_Mnv1Pi.root \
                                      -d NOvAND/GENIE/ReWeighted_to_Mnv1Pi/${TGTDIR}/${SPC} &

                  bin/fakedatavalid.exe -i flattrees/t2knova.flattree.GENIE.NOvAND.${TGT}.2020.${SPC}.root \
                                      --From Generated -F FDSInputs/FakeDataInputs_FromGenerated.root \
                                      -H config/FakeDataValidConfig_NOvAND.toml \
                                      -a ${AFLAG} -T Generated -W NOvA_to_T2KNonQE_ptlep \
                                      -o ${OUTDIR}/FakeDataValid_ReWeighted_to_NonQE.root \
                                      -d NOvAND/GENIE/ReWeighted_to_NonQE/${TGTDIR}/${SPC} &

                  wait
              done

            fi
        done
    done
done

rm -f ${OUTDIR}/FakeDataValid.root
hadd -j 2 ${OUTDIR}/FakeDataValid.root \
          ${OUTDIR}/FakeDataValid_*.root

if [ ! -e ${OUTDIR}/FakeDataValid.root ]; then
  echo "Didn't produce FakeDataValid.root. Failed"
  exit 1
fi