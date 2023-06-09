#! /bin/bash

outfile=""
plttype="dens"
datadim=3
plane=Z
norm=10001000
delay=20
scale=1e70
range=*:*
contour=0

while [ $# -gt 0 ] ; do
	case "$1" in
	-t) plttype="$2"; shift;;
	-o) outfile="$2"; shift;;
	-d) datadim="$2"; shift;;
	-p) plane="$2"; shift;;
    -n) norm="$((10000000+$2))"; shift;;
	-l) delay="$2"; shift;;
	-s) scale="$2"; shift;;
	-r) range="-$2:$2"; shift;;
	--) shift;;
	-h) echo -e "Usage: $0 -o outfile [-d dim] [-t type] [-p plane]\n                               [-n norm] [-l delay] [-s scale]\n  dim: 3 (default) or 2\n type:  dens (default), cont, grav, phas, accel, vrot or 3d\nplane: Z (default), X or Y\n norm: normalize graph at this slide (default: 1000)\ndelay: frame delay (units of 10 ms; default: 20)\nscale: value scale for 3D plots (default: 1e70)\nrange: coordinate range for 3D plots (default: 100)" 1>&2; exit 1;;
	-?) echo -e "Run '$0 -h' for help" 1>&2; exit 1;;
	-*) echo "$0: error - unrecognized option $1" 1>&2; exit 1;;
	*) ;;
	esac
	shift
done

if [ "$outfile" == "" ] ; then
	echo "$0: error - missing output file" 1>&2; exit 1
fi

if [[ "$plttype" != "dens" && "$plttype" != "phas" && "$plttype" != "grav" && $plttype != "accel" && $plttype != "vrot" && $plttype != "3d" && $plttype != "cont" ]]
then
	echo "$0: error - unknown plot type $plttype" 1>&2; exit 1
fi

if [ "$plttype" == "cont" ] ; then
	plttype="dens"
	contour=1
fi

if [[ "$datadim" != "2" && "$datadim" != "3" ]] ; then
	echo "$0: error - incorrect data file dimension $datadim" 1>&2; exit 1
fi

if [[ "$plane" == "x" ]] ; then
	plane="X"
fi

if [[ "$plane" == "y" ]] ; then
	plane="Y"
fi

if [[ "$plane" == "z" ]] ; then
	plane="Z"
fi

if [[ "$plane" != "" && "$plane" != "X" && "$plane" != "Y" && "$plane" != "Z" ]]
then
	echo "$0: error - incorrect projection plane $plane" 1>&2; exit 1
fi

using="1:2:(abs(\$3)"
mincp=`cat densZ0000000.dat | cut -d ' ' -f 1 | perl -e '$min=1e99;while(<STDIN>){$next=$_;if(($next!="")and(0.0+$next<$min)and(0.0+$next>=0)){$min=0.0+$next;}}print 1.1*$min;'`

if [[ "$plttype" == "accel" ]] ; then
(
	echo set terminal gif medium size 640,480
	echo set output \"/dev/null\"
	echo unset title
	echo set format x
	echo set format y
	echo set yrange [*:*] writeback
	echo old_x=NaN
	echo old_y=NaN
	echo "plot 'grav${plane}${norm:1:7}.dat' using ((abs(\$2)<=$mincp&&abs(\$3)<=$mincp)?(dx=\$1-old_x,old_x=\$1,\$1-dx/2):1/0):((abs(\$2)<=$mincp&&abs(\$3)<=$mincp)?(dy=\$4-old_y,old_y=\$4,-dy/dx):1/0) w l notitle"
	echo set terminal gif medium size 640,480 animate delay ${delay} loop 0 optimize
	echo set output \"$outfile\"
	echo set yrange restore
	for i in {10000000..11500000..1} ; do
		if [ -f grav${plane}${i:1:7}.dat ] ; then
			stdbuf -o 0 echo -n -e "Processing frame $((i-10000000))\r" 1>&2
			echo set label 1 \"${i:1:7}\" at screen 0.88, screen 0.88 tc rgb \"#008000\"
			echo old_x=NaN
			echo old_y=NaN
			echo "plot 'grav${plane}${i:1:7}.dat' using ((abs(\$2)<=$mincp&&abs(\$3)<=$mincp)?(dx=\$1-old_x,old_x=\$1,\$1-dx/2):1/0):((abs(\$2)<=$mincp&&abs(\$3)<=$mincp)?(dy=\$4-old_y,old_y=\$4,-dy/dx):1/0) w l notitle"
		fi
	done
	echo 1>&2
) | gnuplot
elif [[ "$plttype" == "vrot" ]] ; then
(
	echo set terminal gif medium size 640,480
	echo set output \"/dev/null\"
	echo unset title
	echo set format x
	echo set format y
	echo set yrange [*:*] writeback
	echo old_x=NaN
	echo old_y=NaN
	echo "plot 'grav${plane}${norm:1:7}.dat' using ((abs(\$2)<=$mincp&&abs(\$3)<=$mincp)?(dx=\$1-old_x,old_x=\$1,\$1-dx/2):1/0):((abs(\$2)<=$mincp&&abs(\$3)<=$mincp)?(dy=\$4-old_y,old_y=\$4,sqrt(abs(dy/dx*\$1))):1/0) w l notitle"
	echo set terminal gif medium size 640,480 animate delay ${delay} loop 0 optimize
	echo set output \"$outfile\"
	echo set yrange restore
	for i in {10000000..11500000..1} ; do
		if [ -f grav${plane}${i:1:7}.dat ] ; then
			stdbuf -o 0 echo -n -e "Processing frame $((i-10000000))\r" 1>&2
			echo set label 1 \"${i:1:7}\" at screen 0.88, screen 0.88 tc rgb \"#008000\"
			echo old_x=NaN
			echo old_y=NaN
			echo "plot 'grav${plane}${i:1:7}.dat' using ((abs(\$2)<=$mincp&&abs(\$3)<=$mincp)?(dx=\$1-old_x,old_x=\$1,\$1-dx/2):1/0):((abs(\$2)<=$mincp&&abs(\$3)<=$mincp)?(dy=\$4-old_y,old_y=\$4,sqrt(abs(dy/dx*\$1))):1/0) w l notitle"
		fi
	done
	echo 1>&2
) | gnuplot
elif [[ "$plttype" == "3d" ]] ; then
(
	echo set terminal gif medium size 480,480 animate delay ${delay} loop 0 optimize
	echo set size square
	echo set ticslevel 0
	echo set view ,,,1.125
	echo unset title
	echo set xrange [$range]
	echo set yrange [$range]
	echo set zrange [$range]
	echo set cbrange [-$scale:-3*$scale]
	echo set palette model XYZ functions gray**0.35,gray**0.5,gray**0.8
	echo unset colorbox
	echo set output \"$outfile\"
	for i in {10000000..11500000..1} ; do
		if [ -f densZ${i:1:7}.dat ] ; then
			stdbuf -o 0 echo -n -e "Processing frame $((i-10000000))\r" 1>&2
			echo set label 1 \"${i:1:7}\" at screen 0.88, screen 0.88 tc rgb \"#008000\"
			echo "splot 'densX${i:1:7}.dat' using (-\$1):(\$2):(-\$3):(\$4>$scale?-\$4:1/0) with points pt 7 ps 2 palette notitle"
		fi
	done
	echo 1>&2
) | gnuplot
else
	case "$plane" in
		(X) using="2:3:(abs(\$1)";;
		(Y) using="(-\$1):3:(abs(\$2)";;
		(Z) using="1:2:(abs(\$3)";;
	esac
	(
		echo set size square
		echo unset title
		echo set pm3d map
		echo set pm3d interpolate 10,10
		echo set format x
		echo set format y
		echo set xrange [$range]
		echo set yrange [$range]
		lc=""
		if [ "$contour" == "1" ] ; then
			echo set contour base
			echo set cntrparam bspline
			echo set view map
			echo "MAX=`cat densZ${norm:1:7}.dat | cut -d ' ' -f 4 | perl -e '$max=0;while(<STDIN>){$next=$_;if(0.0+$next>$max){$max=$next;}}print $max;'`"
			echo set cntrparam levels increment 0,MAX/12,MAX
			echo unset clabel
			echo unset pm3d
			echo unset surface
			lc=" pt 0"
		fi
		if [ "$plttype" != "phas" ] ; then
			if [ "$contour" == "1" ] ; then
				echo set terminal gif tiny size 240,240
			else
				echo set terminal gif medium size 480,480
			fi
			echo set output \"/dev/null\"
			echo "set cbrange [*:*] writeback"
			echo -n "splot \"${plttype}${plane}${norm:1:7}.dat\""
			if [ "$datadim" == "3" ] ; then
				echo -n " using ${using}<$mincp?\$4:1/0) $lc"
			fi
			echo " notitle"
		fi
		if [ "$contour" == "1" ] ; then
			echo set terminal gif tiny size 240,240 animate delay ${delay} loop 0 optimize
		else
			echo set terminal gif medium size 480,480 animate delay ${delay} loop 0 optimize
		fi
		echo set output \"$outfile\"
		if [ "$plttype" == "phas" ] ; then
			echo "set zrange [-pi:pi]"
			echo "set cbrange [-pi:pi]"
			echo "set palette rgbformulae 15,12,14"
		else
			echo set cbrange restore
		fi
		for i in {10000000..11500000..1} ; do
			if [ -f ${plttype}${plane}${i:1:7}.dat ] ; then
				stdbuf -o 0 echo -n -e "Processing frame $((i-10000000))\r" 1>&2
				if [ "$contour" == "1" ] ; then
   	 				echo set label 1 \"${i:1:7}\" at screen 0.82, screen 0.88 tc rgb \"#008000\"
				else
   	 				echo set label 1 \"${i:1:7}\" at screen 0.88, screen 0.88 tc rgb \"#008000\"
				fi
				echo -n "splot \"${plttype}${plane}${i:1:7}.dat\""
				if [ "$datadim" == "3" ] ; then
					echo -n " using ${using}<$mincp?\$4:1/0)"
				fi
				echo "$lc notitle"
			fi
		done
		echo 1>&2
	) | gnuplot
fi
