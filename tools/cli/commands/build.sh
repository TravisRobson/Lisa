#!/usr/bin/env bash

set -o errexit
. "$ROOT_DIR/tools/cli/utils.sh"


cd $ROOT_DIR/build
NUM_CORES=$(echo -e "import multiprocessing\nprint(int(1.5 * multiprocessing.cpu_count()))" | python3)

cd $ROOT_DIR/build
cmd="make -j${NUM_CORES}" # \todo prepend, optionally, with VERBOSE=1 
cli_log "${cmd}"
eval $cmd