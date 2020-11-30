#!/usr/bin/env bash

cli_log() {
  script_name=${0##*/}
  timestamp=$(date -u +"%Y-%m-%dT%H:%M:%SZ")

  echo -e "== $script_name $timestamp $1"
}