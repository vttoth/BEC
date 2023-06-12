#!/bin/sh

# Copyright (c) 2023 Eniko J. M. Madarassy and Viktor T. Toth
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

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

