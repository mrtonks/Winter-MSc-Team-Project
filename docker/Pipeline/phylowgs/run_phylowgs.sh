#!/bin/sh

OPTIONS=$(getopt -o s:p:c:b:m:i:d:r: -l,--sample-name,--purity,--num-chains,--burnin-samples,--mcmc-samples,--num-iter,--path-to-data,--random-seed -- "$@")

if [ $? -ne 0 ]; then
  echo "getopt error"
  exit 1
fi

eval set -- $OPTIONS

RANDOM_SEED=-1

while true; do
  case "$1" in
    -s|--sample-name) SAMPLE_NAME="$2" ; shift ;;
    -p|--purity)  PURITY="$2"; shift ;;
    -c|--num-chains) NUM_CHAINS="$2"; shift ;;
    -b|--burnin-samples)  BURNIN_SMPLS="$2"; shift ;;
    -m|--mcmc-samples)  MCMC_SMPLS="$2"; shift ;;
    -i|--num-iter) NUM_ITER="$2"; shift ;;
    -d|--path-to-data) PATH_TO_DATA="$2"; shift ;;
    -r|--random-seed) RANDOM_SEED="$2"; shift ;;
    --)        shift ; break ;;
    *)         echo "unknown option: $1" ; exit 1 ;;
  esac
  shift
done

if [ $# -ne 0 ]; then
  echo "unknown option(s): $@"
  exit 1
fi

echo "----- Sample to evolve: -----"
echo "Sample Name: $SAMPLE_NAME"
echo "Purity: $PURITY"

rm -rf "./phylowgs/phylowgs.results/$SAMPLE_NAME.results"
mkdir -p "./phylowgs/phylowgs.results/$SAMPLE_NAME.results"

python ./phylowgs/phylowgs.source.code.project18/parser/parse_cnvs.py -c "$PURITY" -f team18 "$PATH_TO_DATA/Segments/$SAMPLE_NAME.segments.txt" --cnv-output "./phylowgs/phylowgs.results/$SAMPLE_NAME.results/cnvs.txt"

python ./phylowgs/phylowgs.source.code.project18/parser/create_phylowgs_inputs.py --regions=all --cnvs "$SAMPLE_NAME=./phylowgs/phylowgs.results/$SAMPLE_NAME.results/cnvs.txt" --vcf-type "$SAMPLE_NAME=project18" "$SAMPLE_NAME=$PATH_TO_DATA/VCF/$SAMPLE_NAME.no_real_info.vcf" --output-cnvs "./phylowgs/phylowgs.results/$SAMPLE_NAME.results/cnv_data.txt" --output-variants "./phylowgs/phylowgs.results/$SAMPLE_NAME.results/ssm_data.txt"

mv ./params.json "./phylowgs/phylowgs.results/$SAMPLE_NAME.results/params.json"

mkdir -p "./phylowgs/phylowgs.results/$SAMPLE_NAME.results/temp_files"

python ./phylowgs/phylowgs.source.code.project18/multievolve.py --num-chains "$NUM_CHAINS" --ssms "./phylowgs/phylowgs.results/$SAMPLE_NAME.results/ssm_data.txt" --cnvs "./phylowgs/phylowgs.results/$SAMPLE_NAME.results/cnv_data.txt" --params "./phylowgs/phylowgs.results/$SAMPLE_NAME.results/params.json" --output-dir "./phylowgs/phylowgs.results/$SAMPLE_NAME.results/" --tmp-dir "./phylowgs/phylowgs.results/$SAMPLE_NAME.results/temp_files" --burnin-samples "$BURNIN_SMPLS" --mcmc-samples "$MCMC_SMPLS" --mh-iterations "$NUM_ITER"

rm -f "./phylowgs/phylowgs.results/$SAMPLE_NAME.results"/chain_?/"trees.zip"

mkdir -p "./phylowgs/phylowgs.results/$SAMPLE_NAME.results/test_results"

python ./phylowgs/phylowgs.source.code.project18/write_results.py "$SAMPLE_NAME" "./phylowgs/phylowgs.results/$SAMPLE_NAME.results/trees.zip" "./phylowgs/phylowgs.results/$SAMPLE_NAME.results/test_results/$SAMPLE_NAME.summ.json.gz" "./phylowgs/phylowgs.results/$SAMPLE_NAME.results/test_results/$SAMPLE_NAME.muts.json.gz" "./phylowgs/phylowgs.results/$SAMPLE_NAME.results/test_results/$SAMPLE_NAME.mutass.zip"
