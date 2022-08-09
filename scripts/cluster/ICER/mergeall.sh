#!/bin/bash

IDIR=/mnt/research/NuInt/generation/

GENERATORS=( NEUT GENIE )
SPECIES=( numu numub nue nueb )
DETECTORS=( NOvAND ND280 )

declare -A DET_MATS
DET_MATS["NOvAND"]="CH"
DET_MATS["ND280"]="H2O CH"

declare -A TUNES

TUNES["NEUT"]="BANFF_PRE BANFF_POST"
TUNES["GENIE"]="2020"

for gen in ${GENERATORS[@]}; do
  for spec in ${SPECIES[@]}; do
    for det in ${DETECTORS[@]}; do
      for mat in ${DET_MATS[${det}]}; do
        for tune in ${TUNES["${gen}"]}; do

          if [ -e t2knova.flattree.${gen}.${det}.${mat}.${tune}.${spec}.root ]; then
            continue
          fi

          INPAT="${IDIR}/t2knova/${gen}/${det}/${mat}/${spec}/${tune}/t2knova.flattree.${gen}.${det}.${mat}.${spec}.*.root"
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
