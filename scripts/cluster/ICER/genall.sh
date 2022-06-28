#!/bin/bash

#if a job submission fails, stop trying to submit as we've hit our QOL quota
set -e

DO_SUBMIT=1

GENERATORS=( NEUT GENIE )
SPECIES=( numu numub nue nueb )
DETECTORS=( NOvAND ND280 )

declare -A FLUXES
FLUXES["NOvAND_numu"]="/var/t2k-nova/fluxes/NOvA/FHC_Flux_NOvA_ND_2017.root,flux_numu"
FLUXES["NOvAND_numub"]="/var/t2k-nova/fluxes/NOvA/RHC_Flux_NOvA_ND_2017.root,flux_numubar"
FLUXES["NOvAND_nue"]="/var/t2k-nova/fluxes/NOvA/FHC_Flux_NOvA_ND_2017.root,flux_numu"
FLUXES["NOvAND_nueb"]="/var/t2k-nova/fluxes/NOvA/RHC_Flux_NOvA_ND_2017.root,flux_numubar"
FLUXES["ND280_numu"]="/var/t2k-nova/fluxes/T2K/t2kflux_2016_plus250kA.root,enu_nd280_numu"
FLUXES["ND280_numub"]="/var/t2k-nova/fluxes/T2K/t2kflux_2016_minus250kA.root,enu_nd280_numub"
FLUXES["ND280_nue"]="/var/t2k-nova/fluxes/T2K/t2kflux_2016_plus250kA.root,enu_nd280_numu"
FLUXES["ND280_nueb"]="/var/t2k-nova/fluxes/T2K/t2kflux_2016_minus250kA.root,enu_nd280_numub"

declare -A NEVS
NEVS["GENIE_numu_NOvAND"]=200000
NEVS["GENIE_numub_NOvAND"]=200000
NEVS["GENIE_nue_NOvAND"]=200000
NEVS["GENIE_nueb_NOvAND"]=200000
NEVS["GENIE_numu_ND280"]=200000
NEVS["GENIE_numub_ND280"]=100000
NEVS["GENIE_nue_ND280"]=100000
NEVS["GENIE_nueb_ND280"]=100000

NEVS["NEUT_numu_NOvAND"]=200000
NEVS["NEUT_numub_NOvAND"]=200000
NEVS["NEUT_nue_NOvAND"]=200000
NEVS["NEUT_nueb_NOvAND"]=200000
NEVS["NEUT_numu_ND280"]=200000
NEVS["NEUT_numub_ND280"]=200000
NEVS["NEUT_nue_ND280"]=200000
NEVS["NEUT_nueb_ND280"]=200000

declare -A NJOBS
NJOBS["GENIE_numu_NOvAND"]=200
NJOBS["GENIE_numub_NOvAND"]=200
NJOBS["GENIE_nue_NOvAND"]=200
NJOBS["GENIE_nueb_NOvAND"]=200
NJOBS["GENIE_numu_ND280"]=200
NJOBS["GENIE_numub_ND280"]=400
NJOBS["GENIE_nue_ND280"]=400
NJOBS["GENIE_nueb_ND280"]=400

NJOBS["NEUT_numu_NOvAND"]=200
NJOBS["NEUT_numub_NOvAND"]=200
NJOBS["NEUT_nue_NOvAND"]=200
NJOBS["NEUT_nueb_NOvAND"]=200
NJOBS["NEUT_numu_ND280"]=200 #200
NJOBS["NEUT_numub_ND280"]=200 #200
NJOBS["NEUT_nue_ND280"]=200 #200
NJOBS["NEUT_nueb_ND280"]=200 #200

declare -A TUNES

TUNES["NEUT"]="BANFF_PRE BANFF_POST"
TUNES["GENIE"]="2020"

declare -A DET_MATS
DET_MATS["NOvAND"]="CH"
DET_MATS["ND280"]="CH H2O"

for gen in ${GENERATORS[@]}; do
  for spec in ${SPECIES[@]}; do
    for det in ${DETECTORS[@]}; do
      for mat in ${DET_MATS[${det}]}; do
        for tune in ${TUNES["${gen}"]}; do

          ./genone.sh --nfiles ${NJOBS["${gen}_${spec}_${det}"]} \
            -G ${gen} \
              -f ${FLUXES["${det}_${spec}"]} \
              -t ${mat} \
              -n ${NEVS["${gen}_${spec}_${det}"]} \
              -P ${spec} \
            -S ${DO_SUBMIT} \
            -T  ${tune} \
            --out-dir t2knova/${gen}/${det}/${mat}/${tune}/${spec} \
            --out-file-stub t2knova.flattree.${gen}.${det}.${mat}.${spec}

        done
      done
    done
  done
done
