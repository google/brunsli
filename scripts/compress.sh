#!/bin/sh

function compressOne() {
  FILE=$1
  case "${FILE##*.}" in
    "br") exit 0;;
    "gz") exit 0;;
    "png") exit 0;;
    "jpg") exit 0;;
    "j") exit 0;;
  esac
  echo "Compressing $FILE"
  ORIG_SIZE=`stat -f %z $FILE`
  brotli -Zfk $FILE
  BR_SIZE=`stat -f %z $FILE.br`
  zopfli --i500 $FILE
  GZ_SIZE=`stat -f %z $FILE.gz`
  if [ $BR_SIZE -ge $ORIG_SIZE ]; then rm $FILE.br; fi
  if [ $GZ_SIZE -ge $ORIG_SIZE ]; then rm $FILE.gz; fi
}

export -f compressOne
all_files=`ls docs/*`
echo $all_files | xargs -n 1 -P `nproc` bash -c 'compressOne "$@"' _
