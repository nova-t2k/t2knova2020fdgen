#!/bin/bash

# Intercept the output directory and number of jobs
declare -a OPTARRAY

while [[ ${#} -gt 0 ]]; do

  key="$1"
  case $key in

      --nfiles)
      NTARGETFILES="$2"
        echo "[OPT]: Submitting up to ${NTARGETFILES} jobs"
      shift # past argument
    ;;

      --out-dir)
        OUTDIR="$2"
        echo "[OPT]: Writing outputs to: ${OUTDIR}"
        OPTARRAY+=("${key}")
        OPTARRAY+=("${2}")
        shift # past argument
      ;;

      --out-file-stub)
        OUTFILESTUB="$2"
        echo "[OPT]: Writing file name stub: ${OUTFILESTUB}"
        OPTARRAY+=("${key}")
        OPTARRAY+=("${2}")
        shift # past argument
      ;;

      *)
        #pass through unknown options
        OPTARRAY+=("${key}")
      ;;
  esac
  shift
done

ODIR=/mnt/research/NuInt/generation/${OUTDIR}/
mkdir -p ${ODIR}

NFILES=$(find $ODIR -name "${OUTFILESTUB}.*.root" | wc -l)

echo "Dir: $ODIR contains $NFILES ${OUTFILESTUB}.*.root files."

JTORUN=$(( NJOBS - NFILES ))

if [ $JTORUN -lt 1 ]; then
  echo "Running jobs: $JTORUN"

  CMD="--array=1-${JTORUN} $(readlink -f t2knovagen.sh) ${OPTARRAY[@]}"
  echo "sbatch ${CMD}"
  sbatch ${CMD}

  if [ "$?" == "1" ]; then
    echo "Failed to submit job"
    exit 1
  fi

fi