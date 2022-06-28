#!/bin/bash

IDIR=/mnt/research/NuInt/generation/

GENERATORS=( NEUT GENIE )
GENERATORS=( NEUT )
SPECIES=( numu numub nue nueb )
DETECTORS=( NOvAND ND280 )
DETECTORS=( ND280 )

declare -A DET_MATS
DET_MATS["NOvAND"]="CH"
DET_MATS["ND280"]="H2O CH"

declare -A TUNES

TUNES["NEUT_CH_numu"]="BANFF_PRE BANFF_POST"
TUNES["NEUT_CH_numub"]="BANFF_PRE BANFF_POST"
TUNES["NEUT_CH_nue"]="BANFF_PRE BANFF_POST"
TUNES["NEUT_CH_nueb"]="BANFF_PRE BANFF_POST"
TUNES["NEUT_H2O_numu"]="BANFF_PRE BANFF_POST"
TUNES["NEUT_H2O_numub"]="BANFF_PRE BANFF_POST"
TUNES["NEUT_H2O_nue"]="BANFF_PRE BANFF_POST"
TUNES["NEUT_H2O_nueb"]="BANFF_PRE BANFF_POST"
TUNES["GENIE_CH_numu"]="2020"
TUNES["GENIE_CH_numub"]="2020"
TUNES["GENIE_CH_nue"]="2020"
TUNES["GENIE_CH_nueb"]="2020"
TUNES["GENIE_H2O_numu"]="2020"
TUNES["GENIE_H2O_numub"]="2020"
TUNES["GENIE_H2O_nue"]="2020"
TUNES["GENIE_H2O_nueb"]="2020"

for gen in ${GENERATORS[@]}; do
  for spec in ${SPECIES[@]}; do
    for det in ${DETECTORS[@]}; do
      for mat in ${DET_MATS[${det}]}; do
        for tune in ${TUNES["${gen}_${mat}_${spec}"]}; do


          INPAT="${IDIR}/t2knova/${gen}/${det}/${mat}/${tune}/${spec}/t2knova.flattree.${gen}.${det}.${mat}.${spec}.*.root"
          NFILES=$(ls ${INPAT} | wc -l)

          echo "Found ${NFILES} matching ${INPAT}"

          if [ "${NFILES}" -lt 5 ]; then
            continue
          fi

          nuis_flat_tree_combiner -i ${INPAT} \
                                  -o t2knova.flattree.${gen}.${det}.${mat}.${tune}.${spec}.root \
                                  -t T2KNOvATruthTree -b fScaleFactor &

        done
      done
    done
    wait
  done
done
