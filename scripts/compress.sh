#!/bin/sh

# Copyright (c) Google LLC 2019
#
# Use of this source code is governed by an MIT-style
# license that can be found in the LICENSE file or at
# https://opensource.org/licenses/MIT.

function updateBr() {
  local ORIG="$1"
  local BR="${ORIG}.br"
  local TMP="$(mktemp)"
  local FAIL=1
  brotli -dfk "${BR}" -o "${TMP}" >/dev/null 2>/dev/null
  cmp -s "${ORIG}" "${TMP}" && FAIL=0
  rm -f "${TMP}" >/dev/null 2>/dev/null
  [[ ${FAIL} = 0 ]] && return
  ORIG_SIZE=`stat -f %z ${ORIG}`
  echo "Creating ${BR}"
  brotli -Zfk "${ORIG}"
  BR_SIZE=`stat -f %z "${BR}"`
  [[ $BR_SIZE -ge $ORIG_SIZE ]] && rm "${BR}"
}

function updateGz() {
  local ORIG="$1"
  local GZ="${ORIG}.gz"
  local TMP="$(mktemp)"
  local FAIL=1
  gzip -cdf "${GZ}" > "${TMP}" 2>/dev/null
  cmp -s "${ORIG}" "${TMP}" && FAIL=0
  rm -f "${TMP}" >/dev/null 2>/dev/null
  [[ ${FAIL} = 0 ]] && return
  ORIG_SIZE=`stat -f %z ${ORIG}`
  echo "Creating ${GZ}"
  zopfli --i500 "${ORIG}"
  GZ_SIZE=`stat -f %z "${GZ}"`
  [[ $GZ_SIZE -ge $ORIG_SIZE ]] && rm "${GZ}"
}

function updateFile() {
  FILE=$1
  case "${FILE##*.}" in
    "br") exit 0;;
    "gz") exit 0;;
    "png") exit 0;;
    "jpg") exit 0;;
    "j") exit 0;;
  esac
  updateBr "${FILE}"
  updateGz "${FILE}"
}

export -f updateBr
export -f updateGz
export -f updateFile
all_files=`ls docs/*`
echo $all_files | xargs -n 1 -P `nproc` bash -c 'updateFile "$@"' _
