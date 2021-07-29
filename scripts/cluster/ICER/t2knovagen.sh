#!/bin/bash --login
########## SBATCH Lines for Resource Request ##########
 
#SBATCH --time=02:00:00 
#SBATCH --mem-per-cpu=2G
#SBATCH --job-name t2knovagengen
 
########## Command Lines for Job Running ##########

#set -x

IMAGE=~/t2k-nova-mc-tools_oa2020_debian_bullseye-slim.sif

OUTDIR=""
TUNE=""
GENERATOR=""
PROBE=""
OUTFILESTUB=""

declare -a OPTARRAY

while [[ ${#} -gt 0 ]]; do

  key="$1"
  case $key in

      -T|--tune)
        TUNE="$2"
        echo "[OPT]: Using Tune: ${TUNE}"
        shift # past argument
      ;;

      -g|--generator)
        GENERATOR="$2"
        echo "[OPT]: Using Generator: ${GENERATOR}"
        shift # past argument
      ;;

      -p|--probe)
        PROBE="$2"
        echo "[OPT]: Using Probe PDG: ${PROBE}"
        OPTARRAY+=("${key}")
        OPTARRAY+=("${2}")
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

if [[ -z ${TUNE} ]] || [[ -z ${PROBE} ]] || [[ -z ${GENERATOR} ]]; then
  echo "[ERROR]: Not all required options passed, each of: -T, -p, and -g are required."
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
  OUTFILESTUB=${GENERATOR}.nupdg_${PROBE}.${TUNE}
fi

OUTFILENAME=${OUTFILESTUB}.${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.root

OPTARRAY+=("-o")
OPTARRAY+=("${OUTFILENAME}")

echo "Running: singularity run ${IMAGE} ${OPTARRAY[@]}"
singularity run ${IMAGE} ${OPTARRAY[@]} &> job_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.log

date

FLATFILENAME=${OUTFILESTUB}.flattree.${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.root

echo "Running: singularity run ${IMAGE} anaev.sh -g ${GENERATOR} -i ${OUTFILENAME} -p ${PROBE} -T ${TUNE} -o ${FLATFILENAME}"
singularity run ${IMAGE} anaev.sh -g ${GENERATOR} -i ${OUTFILENAME} -p ${PROBE} -T ${TUNE} -o ${FLATFILENAME}  &>> job_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.log

cat *.card

if [[ ! -e /mnt/research/NuInt/generation/${OUTDIR}/ ]]; then
  mkdir -p /mnt/research/NuInt/generation/${OUTDIR}/
fi

date

mv ${FLATFILENAME} *.log /mnt/research/NuInt/generation/${OUTDIR}/
