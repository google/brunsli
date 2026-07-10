#!/usr/bin/env bash

# Copyright (c) Google LLC 2026
#
# Use of this source code is governed by an MIT-style
# license that can be found in the LICENSE file or at
# https://opensource.org/licenses/MIT.

set -eu

SELF=$(realpath "$0")
SELF_DIR=$(dirname "${SELF}")
MYDIR=$(realpath "${SELF_DIR}/..")
TMPDIR=$(mktemp -d)
TEXT_BOLD_PURPLE="\033[1;35m"
TEXT_RESET="\033[0m"
COLORDIFF_BIN="cat"

cleanup() {
  rm -rf "${TMPDIR}"
}
trap cleanup EXIT

if [[ -t 1 ]]; then
  COLORDIFF_BIN=$(which colordiff cat 2>/dev/null | head -n 1)
fi

# It is ok, if buildifier is not installed.
if which buildifier >/dev/null; then
  BUILDIFIER_PATCH="${TMPDIR}/buildifier.patch"
  BAZEL_FILES=`git -C "${MYDIR}" ls-files | grep -E "/BUILD$|WORKSPACE|.bzl$|.bazel$"`
  set -x
  buildifier -d ${BAZEL_FILES} >"${BUILDIFIER_PATCH}" || true
  { set +x; } 2>/dev/null
  if [ -s "${BUILDIFIER_PATCH}" ]; then
    echo 'buildifier have found some problems in Bazel build files:' >&2
    "${COLORDIFF_BIN}" <"${BUILDIFIER_PATCH}"
    echo 'To fix them run (from the base directory):' >&2
    echo '  buildifier `git ls-files | grep -E "/BUILD$|WORKSPACE|.bzl$|.bazel$"`' >&2
    exit 1
  fi
else
  echo -e "${TEXT_BOLD_PURPLE}SKIPPED:${TEXT_RESET} buildifier (not installed)"
fi

# It is ok, if spell-checker is not installed.
if which typos >/dev/null; then
  SRC_EXT="bazel|bzl|c|cc|cmake|gni|h|html|in|java|js|m|md|nix|py|rst|sh|ts|txt|yaml|yml"
  SOURCES=`git -C "${MYDIR}" ls-files | grep -E "\.(${SRC_EXT})$"`
  typos -c "${MYDIR}/scripts/typos.toml" ${SOURCES}
else
  echo -e "${TEXT_BOLD_PURPLE}SKIPPED:${TEXT_RESET} typos not installed; try: cargo install typos-cli"
fi

# It is ok, if zizmor is not installed.
if which zizmor >/dev/null; then
  zizmor "${MYDIR}/.github/workflows/"
else
  echo -e "${TEXT_BOLD_PURPLE}SKIPPED:${TEXT_RESET} zizmor (not installed)"
fi

GIT_EMPTY_TREE=`git commit-tree 4b825dc642cb6eb9a060e54bf8d69288fbee4904 -m "empty"`
# We consider version served by GitHub to be golden. As for now it is clang-format-18.
# TODO(eustas): make this (and other checks) hermetic / reproducible
CLANG_FORMAT=`git --help -a | grep -Eo "(clang-format(-[0-9.]+)?)" | tail -n 1`
if ! which "${CLANG_FORMAT}" >/dev/null; then
  echo "You must install clang-format for \"git clang-format\"" >&2
  exit 1
fi
echo "Using ${CLANG_FORMAT}" >&2
DEFAULT_CLANG_PATCH="${TMPDIR}/${CLANG_FORMAT}.patch"
CLANG_PATCH="${CLANG_PATCH:-${DEFAULT_CLANG_PATCH}}"
set -x
git -C "${MYDIR}" "${CLANG_FORMAT}" --binary "${CLANG_FORMAT}" \
  --style=file --diff "${GIT_EMPTY_TREE}" -- >"${CLANG_PATCH}" || true
{ set +x; } 2>/dev/null
if grep -E '^--- ' "${CLANG_PATCH}" >/dev/null; then
  echo "clang-format findings:" >&2
  "${COLORDIFF_BIN}" < "${CLANG_PATCH}"
  echo "clang-format found issues."
  echo "Run \`./scripts/lint.sh | patch -p1\` to apply fixes." >&2
  exit 1
else
  echo "clang-format check OK" >&2
fi
