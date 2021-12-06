#!/usr/bin/env bash
# EXAMPLE USAGE:
# bash tests/runRandomClients.sh -c 1 -s sampleBC -p 9123 -t .005\
#   -g tests/test_outputs/sampleBC/bc_gene_map.txt \
#   -d tests/test_outputs/sampleBC/datasets.txt

# get commandline args - process the -h help arg
args=("${@}")
for i in ${!args[@]}; do
  #echo "$i = ${args[i]}"
  if [[ ${args[i]} = "-c" ]]; then
    # Get the number of clients value and remove the arg from the list
    NUM_CLIENTS=${args[i+1]}
    unset 'args[i]'
    unset 'args[i+1]'
  elif [[ ${args[i]} = "-h" ]]; then
    echo "USAGE: $0 [-c <num_clients>] [-s <species>] [-p <port>] [-g <gene_file>] [-d <datasets_file>] [-t <hours_to_run>]"
    exit 0
  fi
done

if [ -z $NUM_CLIENTS ]; then
    echo "Error: Must specify number of clients -c <num_clients>"
    exit -1
fi

# activate genomics conda env if needed
if [ -z $CONDA_DEFAULT_ENV ] || [ $CONDA_DEFAULT_ENV != "genomics" ]; then
    source ~/.bashrc
    conda activate genomics
fi

WORKDIR=$(dirname "$0")

# Spawn off NUM_CLIENTS processes

PIDS=()
for i in $(seq $NUM_CLIENTS); do
    echo "Start client $i"
    python $WORKDIR/randomSeekRpcClient.py ${args[@]} &
    PIDS+=($!)
done

wait
