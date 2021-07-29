#!/bin/bash

IDIR=/mnt/research/NuInt/generation/t2knova

for gen in NEUT GENIE; do

  TUNE="banffpost"
  if [ "${gen}" == "GENIE" ]; then
    TUNE="nova2020"
  fi


  for nuspec in numu nue numub nueb; do

    HC="FHC"
    if [ "${nuspec}" == "numub" ] || [ "${nuspec}" == "nueb" ]; then
      HC="RHC"
    fi

    specswap="${nuspec}"
    if [ "${nuspec}" == "nue" ] || [ "${nuspec}" == "nueb" ]; then
      specswap="${nuspec}swap"
    fi
    

    for det in NOvAND ND280; do
      #NEUT_ND280_RHC_nuebswap.flattree.22747639_92.root

      INPAT="${IDIR}/${gen}/${det}/${nuspec}/${gen}_${det}_${HC}_${specswap}.flattree.*.root"
      NFILES=$(ls ${INPAT} | wc -l)

      echo "Found ${NFILES} matching ${INPAT}"

      if [ "${NFILES}" -lt 5 ]; then
        continue
      fi

      nuis_flat_tree_combiner -i ${INPAT} \
                              -o flat.${TUNE}.${det}.${nuspec}.CH.${HC}.${gen}.root \
                              -t T2KNOvATruthTree -b fScaleFactor

    done

  done

done
