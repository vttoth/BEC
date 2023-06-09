#!/bin/sh

rate=20

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

CURDIR=`pwd`
TMPDIR=`mktemp -d`
srcfile=`readlink -f $srcfile`
dstfile=`readlink -f $dstfile`
cd $TMPDIR
mplayer -speed 100 -vo png $srcfile
cd $CURDIR
ffmpeg -sameq -r $rate -i $TMPDIR/%08d.png -y -an $dstfile
rm $TMPDIR/*.png
rmdir $TMPDIR
