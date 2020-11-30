#!/usr/bin/env bash

set -o errexit
. "$ROOT_DIR/tools/cli/utils.sh"


# Help message
cli_help_configure() {
  echo "
Command: configure

Usage:
  configure <options> <args>

options:
  --debug,   -d  Debug build [default CMAKE_BUILD_TYPE is RelWithDebInfo].
  --log,     -l  Turn on logging.
  --prefix,  -p  Installation prefix (default is project's root).
  --release, -r  Release build.
  --no-tests     Build without tests.
  --no-openmp    Build without OpenMP.
  --no-profile   Build without Lisa profiling.
"
  exit 1
}


# defaults
BUILD_TYPE=RelWithDebInfo
INSTALL_PREFIX=$ROOT_DIR
ENABLE_LISA_TESTS=ON
ENABLE_LISA_OPENMP=OFF
ENABLE_LISA_LOGGING=OFF
ENABLE_LISA_PROFILE=ON
BUILD_LEGACY_LISA=ON


# parse the command line options
# https://stackoverflow.com/a/33826763
cli_log "Parsing command options: $@"
while [[ "$#" -gt 0 ]]; do
  case "$1" in
    -d|--debug) BUILD_TYPE=Debug ;;
    -l|--log) ENABLE_LISA_LOGGING=ON ;;
    -p|--prefix) INSTALL_PREFIX="$2"; shift ;;
    -r|--release) BUILD_TYPE=Release ;;
    --no-tests) ENABLE_LISA_TESTS=OFF ;;
    --no-openmp) ENABLE_LISA_OPENMP=OFF ;;
    --no-profile) ENABLE_LISA_PROFILE=OFF ;;
    *) cli_help_configure ;;
  esac
  shift
done


blue="\033[0;94m"
yellow="\033[0;33m"
bold="\033[1m"
reset="\033[0m"

# Report the chosen configuration
cli_log "${blue}${bold}Build type${reset}................${yellow}${bold}${BUILD_TYPE}${reset}"
cli_log "${blue}${bold}Install prefix${reset}............${yellow}${bold}${INSTALL_PREFIX}${reset}"
cli_log "${blue}${bold}Enable tests${reset}..............${yellow}${bold}${ENABLE_LISA_TESTS}${reset}"
cli_log "${blue}${bold}Use OpenMP${reset}................${yellow}${bold}${ENABLE_LISA_OPENMP}${reset}"
cli_log "${blue}${bold}Enable Lisa profiling${reset}.....${yellow}${bold}${ENABLE_LISA_PROFILE}${reset}"
cli_log "${blue}${bold}Enable logging${reset}............${yellow}${bold}${ENABLE_LISA_LOGGING}${reset}"
cli_log "${blue}${bold}Build legacy code${reset}.........${yellow}${bold}${BUILD_LEGACY_LISA}${reset}"

cli_log "Creating build/"
mkdir -p $ROOT_DIR/build
cd $ROOT_DIR/build


# construct CMake command and execute
export CC=/usr/local/bin/gcc-10
export CXX=/usr/local/bin/g++-10
cmd="cmake -G \"Unix Makefiles\" .."
cmd="${cmd} -D CMAKE_BUILD_TYPE:STRING=${BUILD_TYPE}"
cmd="${cmd} -D CMAKE_INSTALL_PREFIX:PATH=${INSTALL_PREFIX}"
cmd="${cmd} -D ENABLE_LISA_TESTS=${ENABLE_LISA_TESTS}"
cmd="${cmd} -D ENABLE_LISA_OPENMP=${ENABLE_LISA_OPENMP}"
cmd="${cmd} -D ENABLE_LISA_PROFILE=${ENABLE_LISA_PROFILE}"
cmd="${cmd} -D ENABLE_LISA_PROFILE=${ENABLE_LISA_LOGGING}"
cmd="${cmd} -D BUILD_LEGACY_LISA=${BUILD_LEGACY_LISA}"
cmd="${cmd} .."
cli_log "${cmd}"
echo "------------------------------ CMake Log ------------------------------"
eval $cmd
echo "-----------------------------------------------------------------------"

