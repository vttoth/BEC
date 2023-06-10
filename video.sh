#!/bin/sh

rate=100

while [ $# -gt 0 ] ; do
    case "$1" in
    -r) rate="$2"; shift;;
    -i) srcfile="$2"; shift;;
    -o) dstfile="$2"; shift;;
    -h) echo -e "Usage: $0 -i source-file -o destination-file [-r frame-rate]" 1>&2; exit 1;;
    -*) echo "$0: error - unrecognized option $1" 1>&2; exit 1;;
    esac
    shift
done

if [ "$srcfile" == "" ] ; then
    echo "$0: error - missing source file" 1>&2; exit 1
fi
if [ "$dstfile" == "" ] ; then
    echo "$0: error - missing destination file" 1>&2; exit 1
fi

ffmpeg -r:v "${rate}/1" -i $srcfile -r:v "30/1" $dstfile

