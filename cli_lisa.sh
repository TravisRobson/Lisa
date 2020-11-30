#!/usr/bin/env bash

# This file must be located in the Lisa project's root directory.
# The goal of this driver to provide a common interface to every standard 
#   operation you do with the raytracer.
# References
# 1) https://medium.com/@brotandgames/build-a-custom-cli-with-bash-e3ce60cfb9a4
# 2) https://github.com/SierraSoftworks/bash-cli

set -o errexit

ROOT_DIR=$(dirname "$(echo -e "import os\nprint(os.path.abspath('"$0"'))" | python3)")

export ROOT_DIR
. "$ROOT_DIR/tools/cli/utils.sh"


cli_help() {
  cli_name=${0##*/}
  echo "
$cli_name
Lisa CLI
Version: $(cat $ROOT_DIR/tools/cli/VERSION)

Usage: $cli_name [command]

commands:
  build      Build the Lisa project.
  clean      Clean up logs and build folder.
  configure  Configure the Lisa's build.
  driver     Execute the main driver
  install    Install Lisa library, executables, headers, and cmake files.
  test       Run tests.
  *          Display help message
"
  exit 1
}


cli_log "Exporting config..."
export $(cat "$ROOT_DIR/tools/cli/config" | xargs)


case "$1" in
  build)
    "$ROOT_DIR/tools/cli/commands/build.sh" | tee -ia "$ROOT_DIR/cli_build.log"
    ;;
  clean)
    "$ROOT_DIR/tools/cli/commands/clean.sh" # you don't want to tee, that's unclean.
    ;;    
  configure)
    shift # move past the "configure" command line field
    "$ROOT_DIR/tools/cli/commands/configure.sh" "$@" | tee -ia "$ROOT_DIR/cli_configure.log"
    ;;  
  install)
    "$ROOT_DIR/tools/cli/commands/install.sh" | tee -ia "$ROOT_DIR/cli_install.log"
    ;;
  *)
    cli_help
    ;;
esac

if [ ${PIPESTATUS} -ne 0 ]; then
  exit 1
fi


  # driver)
  #   shift # move past the "driver" command line field
  #   "$ROOT_DIR/tools/cli/commands/driver" "$@" | tee -ia "$ROOT_DIR/cli_driver_${@}.log"
  #   ;;
  # test)
  #   "$ROOT_DIR/tools/cli/commands/test" "$2" | tee -ia "$ROOT_DIR/cli_test_${2}.log"
  #   ;;



