#!/bin/bash

#if a job submission fails, stop trying to submit as we've hit our QOL quota
set -e

GENERATORS=( GENIE NEUT )
SPECIES=( numu numub nue nueb )
DETECTORS=( NOvAND ND280 )

declare -A FLUXES
FLUXES["NOvAND_numu"]="/var/t2k-nova/fluxes/NOvA/FHC_Flux_NOvA_ND_2017.root,flux_numu"
FLUXES["NOvAND_numub"]="/var/t2k-nova/fluxes/NOvA/RHC_Flux_NOvA_ND_2017.root,flux_numubar"
FLUXES["NOvAND_nue"]="/var/t2k-nova/fluxes/NOvA/FHC_Flux_NOvA_ND_2017.root,flux_numu"
FLUXES["NOvAND_nueb"]="/var/t2k-nova/fluxes/NOvA/RHC_Flux_NOvA_ND_2017.root,flux_numubar"
FLUXES["ND280_numu"]="/var/t2k-nova/fluxes/T2K/t2kflux_2016_minus250kA.root,enu_nd280_numu"
FLUXES["ND280_numub"]="/var/t2k-nova/fluxes/T2K/t2kflux_2016_plus250kA.root,enu_nd280_numub"
FLUXES["ND280_nue"]="/var/t2k-nova/fluxes/T2K/t2kflux_2016_minus250kA.root,enu_nd280_numu"
FLUXES["ND280_nueb"]="/var/t2k-nova/fluxes/T2K/t2kflux_2016_plus250kA.root,enu_nd280_numub"

declare -A NEVS
NEVS["GENIE_numu_NOvAND"]=200000
NEVS["GENIE_numub_NOvAND"]=200000
NEVS["GENIE_nue_NOvAND"]=200000
NEVS["GENIE_nueb_NOvAND"]=200000
NEVS["GENIE_numu_ND280"]=200000
NEVS["GENIE_numub_ND280"]=200000
NEVS["GENIE_nue_ND280"]=200000
NEVS["GENIE_nueb_ND280"]=200000

NEVS["NEUT_numu_NOvAND"]=200000
NEVS["NEUT_numub_NOvAND"]=200000
NEVS["NEUT_nue_NOvAND"]=200000
NEVS["NEUT_nueb_NOvAND"]=200000
NEVS["NEUT_numu_ND280"]=200000
NEVS["NEUT_numub_ND280"]=200000
NEVS["NEUT_nue_ND280"]=200000
NEVS["NEUT_nueb_ND280"]=200000

declare -A NJOBS
NJOBS["GENIE_numu_NOvAND"]=1
NJOBS["GENIE_numub_NOvAND"]=1
NJOBS["GENIE_nue_NOvAND"]=1
NJOBS["GENIE_nueb_NOvAND"]=1
NJOBS["GENIE_numu_ND280"]=1
NJOBS["GENIE_numub_ND280"]=1
NJOBS["GENIE_nue_ND280"]=1
NJOBS["GENIE_nueb_ND280"]=1

NJOBS["NEUT_numu_NOvAND"]=1
NJOBS["NEUT_numub_NOvAND"]=1
NJOBS["NEUT_nue_NOvAND"]=1
NJOBS["NEUT_nueb_NOvAND"]=1
NJOBS["NEUT_numu_ND280"]=1
NJOBS["NEUT_numub_ND280"]=1
NJOBS["NEUT_nue_ND280"]=1
NJOBS["NEUT_nueb_ND280"]=1

declare -A TUNES

TUNES["NEUT_CH_numu"]="BANFF_POST_C_NU"
TUNES["NEUT_CH_numub"]="BANFF_POST_C_NUB"
TUNES["NEUT_CH_nue"]="BANFF_POST_C_NU"
TUNES["NEUT_CH_nueb"]="BANFF_POST_C_NUB"
TUNES["NEUT_H2O_numu"]="BANFF_POST_O_NU"
TUNES["NEUT_H2O_numub"]="BANFF_POST_O_NUB"
TUNES["NEUT_H2O_nue"]="BANFF_POST_O_NU"
TUNES["NEUT_H2O_nueb"]="BANFF_POST_O_NUB"
TUNES["GENIE_CH_numu"]="2020"
TUNES["GENIE_CH_numub"]="2020"
TUNES["GENIE_CH_nue"]="2020"
TUNES["GENIE_CH_nueb"]="2020"
TUNES["GENIE_H2O_numu"]="2020"
TUNES["GENIE_H2O_numub"]="2020"
TUNES["GENIE_H2O_nue"]="2020"
TUNES["GENIE_H2O_nueb"]="2020"

declare -A DET_MATS
DET_MATS["NOvAND"]="CH"
DET_MATS["ND280"]="H2O CH"

declare -A SCRIPTS
SCRIPTS["NEUT"]=nuis_genev_neut
SCRIPTS["GENIE"]=nuis_genev_genie

declare -A SPEC_PDG
SPEC_PDG["numu"]=14
SPEC_PDG["numub"]=-14
SPEC_PDG["nue"]=12
SPEC_PDG["nueb"]=-12

for gen in ${GENERATORS[@]}; do
  for spec in ${SPECIES[@]}; do
    for det in ${DETECTORS[@]}; do
      for mat in ${DET_MATS[${det}]}; do

      ./genone.sh --nfiles ${NJOBS["${gen}_${spec}_${det}"]} \
          SCRIPTS["${gen}"] \
              -g ${gen} \
              -t ${mat} \
              -p ${SPEC_PDG["${spec}"]} \
              -n ${NEVS["${gen}_${spec}_${det}"]} \
              -f ${FLUXES["${det}_${spec}"]} \
            -T ${TUNES["${gen}_${mat}_${spec}"]} \
            --out-dir t2knova/${gen}/${det}/${mat}/${spec} \
            --out-file-stub t2knova.flattree.${gen}.${det}.${mat}.${spec}

      done
    done
  done
done