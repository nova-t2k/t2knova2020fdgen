#!/bin/bash

DO_SUBMIT=0

# Intercept the output directory and number of jobs
declare -a OPTARRAY

while [[ ${#} -gt 0 ]]; do

  key="$1"
  case $key in

      --nfiles)
        NTARGETFILES="$2"
        shift # past argument
      ;;

      -S)
        DO_SUBMIT=${2}
        shift
      ;;

      -G|--generator)
        GENERATOR="$2"
        OPTARRAY+=("${key}")
        OPTARRAY+=("${2}")
        shift # past argument
      ;;

      --out-dir)
        OUTDIR="$2"
        OPTARRAY+=("${key}")
        OPTARRAY+=("${2}")
        shift # past argument
      ;;

      --out-file-stub)
        OUTFILESTUB="$2"
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

declare -A TUNES

TUNES["NEUT"]="BANFF_PRE BANFF_POST"
TUNES["GENIE"]="2020"

for tune in ${TUNES[${GENERATOR}]}; do

  ODIR=/mnt/research/NuInt/generation/${OUTDIR}/${tune}/
  mkdir -p ${ODIR}

  NFILES=$(find $ODIR -name "${OUTFILESTUB}.*.root" | wc -l)


  JTORUN=$(( NTARGETFILES - NFILES ))

  if [ ${JTORUN} -lt 1 ]; then 
    echo -e "\u001b[32m[DONE]\u001b[0m: $ODIR contains ${NFILES}/${NTARGETFILES} (${OUTFILESTUB}.*.root) files."
  else
    echo -e "\u001b[31m[GENE]\u001b[0m: Will generate \u001b[31m${JTORUN}\u001b[0m new files: $ODIR contains ${NFILES}/${NTARGETFILES} (${OUTFILESTUB}.*.root) files."
  fi

  if [ $DO_SUBMIT -gt 0 ] && [ $JTORUN -gt 0 ]; then
    CMD="--array=1-${JTORUN} $(readlink -f t2knovagen.sh) ${OPTARRAY[@]}"
    echo "sbatch ${CMD}"
    sbatch ${CMD}

    if [ "$?" == "1" ]; then
      echo "Failed to submit job"
      exit 1
    fi

    #If we have submitted for one tune, don't for the other
    break;
  fi

done
