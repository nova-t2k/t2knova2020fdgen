#!/bin/bash --login
########## SBATCH Lines for Resource Request ##########
 
#SBATCH --time=02:00:00 
#SBATCH --mem-per-cpu=2G
#SBATCH --job-name t2knovagengen
 
########## Command Lines for Job Running ##########

#set -x

IMAGE=~/t2k-nova-mc-tools_oa2020_latest.sif

OUTDIR=""
TUNE=""
GENERATOR=""
PROBE=""
OUTFILESTUB=""
TARGET=""

declare -a OPTARRAY

while [[ ${#} -gt 0 ]]; do

  key="$1"
  case $key in

      -t|--target)
        TARGET="$2"
        echo "[OPT]: Generating on ${TARGET} target"
        OPTARRAY+=("${key}")
        OPTARRAY+=("${TARGET}")
        shift # past argument
      ;;

      -G|--generator)
        GENERATOR="$2"
        echo "[OPT]: Using Generator: ${GENERATOR}"
        shift # past argument
      ;;

      -P|--probe)
        if [ "${2}" == "nueb" ]; then
          PROBE=nuebar
        elif [ "${2}" == "numub" ]; then
          PROBE=numubar
        else
          PROBE="$2"
        fi
           
        echo "[OPT]: Using Probe: ${PROBE}"
        OPTARRAY+=("${key}")
        OPTARRAY+=("${PROBE}")
        shift # past argument
      ;;
      
      --out-dir)
        OUTDIR="$2"
        echo "[OPT]: Writing outputs to: ${OUTDIR}"
        shift # past argument
      ;;

      --out-file-stub)
        OUTFILESTUB="$2"
        echo "[OPT]: Writing file name stub: ${OUTFILESTUB}"
        shift # past argument
      ;;

      *)
      #pass through unknown options
      OPTARRAY+=("${key}")
      ;;
  esac
  shift
done

if [[ -z ${PROBE} ]] || [[ -z ${TARGET} ]] || [[ -z ${GENERATOR} ]]; then
  echo "[ERROR]: Not all required options passed, each of: -T, -t, -p, and -g are required."
  exit 1
fi

if [[ -z "${SLURM_ARRAY_TASK_ID}" ]]; then
  SLURM_ARRAY_TASK_ID=1
fi

if [[ -z "${SLURM_ARRAY_JOB_ID}" ]]; then
  SLURM_ARRAY_JOB_ID=${RANDOM}
fi

if [[ -z "${TMPDIR}" ]]; then
  TMPDIR=$(readlink -f .)
fi

DIR=${TMPDIR}/t2knovagen_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}
mkdir ${DIR}
cd ${DIR}
echo ${PWD}

date

if [[ -z ${OUTFILESTUB} ]]; then
  OUTFILESTUB=${GENERATOR}.${PROBE}
fi

OUTFILENAME=${GENERATOR}.vector.${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.root

OPTARRAY+=("-o")
OPTARRAY+=("${OUTFILENAME}")

if [ "${GENERATOR}" == "NEUT" ]; then #ensure that the generation options match T2K production
  OPTARRAY+=("-x")
  OPTARRAY+=("NEUT-MAQE=1.21")
  OPTARRAY+=("-x")
  OPTARRAY+=("NEUT-MVSPI=0.84")
  OPTARRAY+=("-x")
  OPTARRAY+=("NEUT-IFF=1")
  OPTARRAY+=("-x")
  OPTARRAY+=("NEUT-MASPI=0.95")
  OPTARRAY+=("-x")
  OPTARRAY+=("NEUT-BGSCL=1.3")
  OPTARRAY+=("-x")
  OPTARRAY+=("NEUT-NRTYPE=1")
  OPTARRAY+=("-x")
  OPTARRAY+=("NEUT-CA5=1.01")
  OPTARRAY+=("-x")
  OPTARRAY+=("NEUT-MDLQE=402")
fi

echo "Running: singularity run ${IMAGE} run ${IMAGE} nuis gen ${GENERATOR} ${OPTARRAY[@]}"
singularity run ${IMAGE} nuis gen ${GENERATOR} ${OPTARRAY[@]} &> job_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.log

date

T=$(( RANDOM % 10 )).$(( RANDOM % 1000 )); echo sleep $T; sleep $T

declare -A TUNES

TUNES["NEUT"]="BANFF_PRE BANFF_POST"
TUNES["GENIE"]="2020"

for tune in ${TUNES[${GENERATOR}]}; do

  FLATFILENAME=${OUTFILESTUB}.${tune}.${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.root

  echo "Running: singularity run ${IMAGE} anaev.sh -g ${GENERATOR} -i ${OUTFILENAME} -P ${PROBE} -T ${tune} -t ${TARGET} -o ${FLATFILENAME}"
  singularity run ${IMAGE} anaev.sh -g ${GENERATOR} -i ${OUTFILENAME} -P ${PROBE} -T ${tune} -t ${TARGET} -o ${FLATFILENAME}  &>> job_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.${tune}.log

  if [[ ! -e /mnt/research/NuInt/generation/${OUTDIR}/${tune} ]]; then
    mkdir -p /mnt/research/NuInt/generation/${OUTDIR}/${tune}
  fi

  date

  mv ${FLATFILENAME} *.${tune}.log /mnt/research/NuInt/generation/${OUTDIR}/${tune}/

done

cat *.card
