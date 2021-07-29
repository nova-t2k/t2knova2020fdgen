#!/bin/bash

NJOBS=1
NEVS=100000

DO_NEUT=1
DO_ND280=1
DO_SK=1
DO_GENIE=0
DO_NOvAND=1

if [[ "${DO_GENIE}" == "1" ]]; then
  if [[ "${DO_SK}" == "1" ]]; then

    sbatch --array=1-${NJOBS} ./t2knovagen.sh nuis_genev_genie \
     --out-dir t2knova/GENIE/SK/nueb \
     -T 2020 -g GENIE -t H2O -p -12 -n ${NEVS} \
    -f /var/t2k-nova/fluxes/T2K/t2kflux_2016_plus250kA.root,enu_sk_numub \
     --out-file-stub GENIE_SK_RHC_nuebswap

    sbatch --array=1-${NJOBS} ./t2knovagen.sh nuis_genev_genie \
     --out-dir t2knova/GENIE/SK/numub \
     -T 2020 -g GENIE -t H2O -p -14 -n ${NEVS} \
     -f /var/t2k-nova/fluxes/T2K/t2kflux_2016_plus250kA.root,enu_sk_numub \
     --out-file-stub GENIE_SK_RHC_numub

    sbatch --array=1-${NJOBS} ./t2knovagen.sh nuis_genev_genie \
     --out-dir t2knova/GENIE/SK/nue \
     -T 2020 -g GENIE -t H2O -p 12 -n ${NEVS} \
     -f /var/t2k-nova/fluxes/T2K/t2kflux_2016_minus250kA.root,enu_sk_numu \
     --out-file-stub GENIE_SK_FHC_nueswap

    sbatch --array=1-${NJOBS} ./t2knovagen.sh nuis_genev_genie \
     --out-dir t2knova/GENIE/SK/numu \
     -T 2020 -g GENIE -t H2O -p 14 -n ${NEVS} \
     -f /var/t2k-nova/fluxes/T2K/t2kflux_2016_minus250kA.root,enu_sk_numu \
     --out-file-stub GENIE_SK_FHC_numu

	fi
  if [[ "${DO_ND280}" == "1" ]]; then

    sbatch --array=1-${NJOBS} ./t2knovagen.sh nuis_genev_genie \
     --out-dir t2knova/GENIE/ND280/nueb \
     -T 2020 -g GENIE -t CH -p -12 -n ${NEVS} \
     -f /var/t2k-nova/fluxes/T2K/t2kflux_2016_plus250kA.root,enu_nd280_numub \
     --out-file-stub GENIE_ND280_RHC_nuebswap

    sbatch --array=1-${NJOBS} ./t2knovagen.sh nuis_genev_genie \
     --out-dir t2knova/GENIE/ND280/numub \
     -T 2020 -g GENIE -t CH -p -14 -n ${NEVS} \
     -f /var/t2k-nova/fluxes/T2K/t2kflux_2016_plus250kA.root,enu_nd280_numub \
     --out-file-stub GENIE_ND280_RHC_numub

    sbatch --array=1-${NJOBS} ./t2knovagen.sh nuis_genev_genie \
     --out-dir t2knova/GENIE/ND280/nue \
     -T 2020 -g GENIE -t CH -p 12 -n ${NEVS} \
     -f /var/t2k-nova/fluxes/T2K/t2kflux_2016_minus250kA.root,enu_nd280_numu \
     --out-file-stub GENIE_ND280_FHC_nueswap

    sbatch --array=1-${NJOBS} ./t2knovagen.sh nuis_genev_genie \
     --out-dir t2knova/GENIE/ND280/numu \
     -T 2020 -g GENIE -t CH -p 14 -n ${NEVS} \
     -f /var/t2k-nova/fluxes/T2K/t2kflux_2016_minus250kA.root,enu_nd280_numu \
     --out-file-stub GENIE_ND280_FHC_numu

  fi
  if [[ "${DO_NOvAND}" == "1" ]]; then

    sbatch --array=1-${NJOBS} ./t2knovagen.sh nuis_genev_genie \
     --out-dir t2knova/GENIE/NOvAND/nueb \
     -T 2020 -g GENIE -t CH -p -12 -n ${NEVS} \
     -f /var/t2k-nova/fluxes/NOvA/RHC_Flux_NOvA_ND_2017.root,flux_numubar \
     --out-file-stub GENIE_NOvAND_RHC_nuebswap

    sbatch --array=1-${NJOBS} ./t2knovagen.sh nuis_genev_genie \
     --out-dir t2knova/GENIE/NOvAND/numub \
     -T 2020 -g GENIE -t CH -p -14 -n ${NEVS} \
     -f /var/t2k-nova/fluxes/NOvA/RHC_Flux_NOvA_ND_2017.root,flux_numubar \
     --out-file-stub GENIE_NOvAND_RHC_numub

    sbatch --array=1-${NJOBS} ./t2knovagen.sh nuis_genev_genie \
     --out-dir t2knova/GENIE/NOvAND/nue \
     -T 2020 -g GENIE -t CH -p 12 -n ${NEVS} \
     -f /var/t2k-nova/fluxes/NOvA/FHC_Flux_NOvA_ND_2017.root,flux_numu \
     --out-file-stub GENIE_NOvAND_FHC_nueswap

    sbatch --array=1-${NJOBS} ./t2knovagen.sh nuis_genev_genie \
     --out-dir t2knova/GENIE/NOvAND/numu \
     -T 2020 -g GENIE -t CH -p 14 -n ${NEVS} \
     -f /var/t2k-nova/fluxes/NOvA/FHC_Flux_NOvA_ND_2017.root,flux_numu \
     --out-file-stub GENIE_NOvAND_FHC_numu

  fi
fi

if [[ "${DO_NEUT}" == "1" ]]; then
  if [[ "${DO_SK}" == "1" ]]; then

    sbatch --array=1-${NJOBS} ./t2knovagen.sh nuis_genev_neut \
     --out-dir t2knova/NEUT/SK/nueb \
     -T BANFF_POST_O_NUB -g NEUT -t H2O -p -12 -n ${NEVS} \
     -f /var/t2k-nova/fluxes/T2K/t2kflux_2016_plus250kA.root,enu_sk_numub \
     --out-file-stub NEUT_SK_RHC_nuebswap

    sbatch --array=1-${NJOBS} ./t2knovagen.sh nuis_genev_neut \
     --out-dir t2knova/NEUT/SK/numub \
     -T BANFF_POST_O_NUB -g NEUT -t H2O -p -14 -n ${NEVS} \
     -f /var/t2k-nova/fluxes/T2K/t2kflux_2016_plus250kA.root,enu_sk_numub \
     --out-file-stub NEUT_SK_RHC_numub

    sbatch --array=1-${NJOBS} ./t2knovagen.sh nuis_genev_neut \
     --out-dir t2knova/NEUT/SK/nue \
     -T BANFF_POST_O_NU -g NEUT -t H2O -p 12 -n ${NEVS} \
     -f /var/t2k-nova/fluxes/T2K/t2kflux_2016_minus250kA.root,enu_sk_numu \
     --out-file-stub NEUT_SK_FHC_nueswap

    sbatch --array=1-${NJOBS} ./t2knovagen.sh nuis_genev_neut \
     --out-dir t2knova/NEUT/SK/numu \
     -T BANFF_POST_O_NU -g NEUT -t H2O -p 14 -n ${NEVS} \
     -f /var/t2k-nova/fluxes/T2K/t2kflux_2016_minus250kA.root,enu_sk_numu \
     --out-file-stub NEUT_SK_FHC_numu

  fi
  if [[ "${DO_ND280}" == "1" ]]; then

    sbatch --array=1-${NJOBS} ./t2knovagen.sh nuis_genev_neut \
     --out-dir t2knova/NEUT/ND280/nueb \
     -T BANFF_POST_C_NUB -g NEUT -t CH -p -12 -n ${NEVS} \
     -f /var/t2k-nova/fluxes/T2K/t2kflux_2016_plus250kA.root,enu_nd280_numub \
     --out-file-stub NEUT_ND280_RHC_nuebswap

    sbatch --array=1-${NJOBS} ./t2knovagen.sh nuis_genev_neut \
     --out-dir t2knova/NEUT/ND280/numub \
     -T BANFF_POST_C_NUB -g NEUT -t CH -p -14 -n ${NEVS} \
     -f /var/t2k-nova/fluxes/T2K/t2kflux_2016_plus250kA.root,enu_nd280_numub \
     --out-file-stub NEUT_ND280_RHC_numub

    sbatch --array=1-${NJOBS} ./t2knovagen.sh nuis_genev_neut \
     --out-dir t2knova/NEUT/ND280/nue \
     -T BANFF_POST_C_NU -g NEUT -t CH -p 12 -n ${NEVS} \
     -f /var/t2k-nova/fluxes/T2K/t2kflux_2016_minus250kA.root,enu_nd280_numu \
     --out-file-stub NEUT_ND280_FHC_nueswap

    sbatch --array=1-${NJOBS} ./t2knovagen.sh nuis_genev_neut \
     --out-dir t2knova/NEUT/ND280/numu \
     -T BANFF_POST_C_NU -g NEUT -t CH -p 14 -n ${NEVS} \
     -f /var/t2k-nova/fluxes/T2K/t2kflux_2016_minus250kA.root,enu_nd280_numu \
     --out-file-stub NEUT_ND280_FHC_numu

  fi
  if [[ "${DO_NOvAND}" == "1" ]]; then

    sbatch --array=1-${NJOBS} ./t2knovagen.sh nuis_genev_neut \
     --out-dir t2knova/NEUT/NOvAND/nueb \
     -T BANFF_POST_C_NUB -g NEUT -t CH -p -12 -n ${NEVS} \
     -f /var/t2k-nova/fluxes/NOvA/RHC_Flux_NOvA_ND_2017.root,flux_numubar \
     --out-file-stub NEUT_NOvAND_RHC_nuebswap

    sbatch --array=1-${NJOBS} ./t2knovagen.sh nuis_genev_neut \
     --out-dir t2knova/NEUT/NOvAND/numub \
     -T BANFF_POST_C_NUB -g NEUT -t CH -p -14 -n ${NEVS} \
     -f /var/t2k-nova/fluxes/NOvA/RHC_Flux_NOvA_ND_2017.root,flux_numubar \
     --out-file-stub NEUT_NOvAND_RHC_numub

    sbatch --array=1-${NJOBS} ./t2knovagen.sh nuis_genev_neut \
     --out-dir t2knova/NEUT/NOvAND/nue \
     -T BANFF_POST_C_NU -g NEUT -t CH -p 12 -n ${NEVS} \
     -f /var/t2k-nova/fluxes/NOvA/FHC_Flux_NOvA_ND_2017.root,flux_numu \
     --out-file-stub NEUT_NOvAND_FHC_nueswap

    sbatch --array=1-${NJOBS} ./t2knovagen.sh nuis_genev_neut \
     --out-dir t2knova/NEUT/NOvAND/numu \
     -T BANFF_POST_C_NU -g NEUT -t CH -p 14 -n ${NEVS} \
     -f /var/t2k-nova/fluxes/NOvA/FHC_Flux_NOvA_ND_2017.root,flux_numu \
     --out-file-stub NEUT_NOvAND_FHC_numu

  fi
fi
