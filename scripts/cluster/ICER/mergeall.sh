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

for gen in ${GENERATORS[@]}; do
  for spec in ${SPECIES[@]}; do
    for det in ${DETECTORS[@]}; do
      for mat in ${DET_MATS[${det}]}; do


        INPAT="${IDIR}/t2knova/${gen}/${det}/${mat}/${spec}/t2knova.flattree.${gen}.${det}.${mat}.${spec}.*.root"
        NFILES=$(ls ${INPAT} | wc -l)

        echo "Found ${NFILES} matching ${INPAT}"

        if [ "${NFILES}" -lt 5 ]; then
          continue
        fi

        nuis_flat_tree_combiner -i ${INPAT} \
                                -o t2knova.flattree.${gen}.${det}.${mat}.${spec}.root \
                                -t T2KNOvATruthTree -b fScaleFactor &

      done

    done
    wait

  done

done
