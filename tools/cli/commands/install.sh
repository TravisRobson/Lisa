#!/usr/bin/env bash

# set -o errexit
. "$ROOT_DIR/tools/cli/utils.sh"


cd $ROOT_DIR/build
NUM_CORES=$(echo -e "import multiprocessing\nprint(int(1.5 * multiprocessing.cpu_count()))" | python3)

cmd="make install -j${NUM_CORES}"
cli_log "Building and installing Lisa project: ${cmd}"
# eval $cmd
# exit $?
($cmd)
exit $?