#!/usr/bin/env bash

set -o errexit
. "$ROOT_DIR/tools/cli/utils.sh"


cmd="rm -rf $ROOT_DIR/build $ROOT_DIR/bin $ROOT_DIR/lib $ROOT_DIR/include *.log $ROOT_DIR/instrumentscli*"
cli_log "Cleaning: ${cmd}"
eval $cmd