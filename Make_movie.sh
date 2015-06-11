#!/bin/bash
# This script produces standard mpeg1 movie out of a sequence of frames.
# This mpeg1 format is just the video stream (ES) without a container (PS), because MS-Windows does not like mpeg1 PS file without a audio stream. 
#
# There are four different output specifier:
#	hq: highest possible bitrate in mpeg1 (vbv full stream length). RECOMMENDED
# vcd: pure mpeg1 for VCD with greetings from the 80´th. This will run even on your toaster, but will never look good.
# dvd: DVD standard settings, should work even on hardware DVD-players, but it´s just mpeg1 ES stream so ... ?
# preview: this is the default settings for a small but fast mpeg1 stream.
#
# The possible input file formats are: jpeg, png, tga, sgi
#
# Frank Wiederschein
# 08.11.2006
#
###############################################################################
#

shopt -s -o nounset
declare -rx SCRIPT=${0##*/}

if [ $# -lt 3 ]; then
	echo "usage: $SCRIPT 'frame-identifier' frame-type file_name.mpg [quality/vcd/dvd/preview]"
	echo "e.g.: $SCRIPT 'movie*.png' png file_name.mpg [hq]"
	exit 192
fi
if [[ $# -eq 4 && $4 = hq ]]; then
	LAVOPTS="vcodec=mpeg1video:vbitrate=16000:keyint=15:mbd=2:mv0:trell:mbcmp=2:precmp=2:subcmp=2:cmp=2:dia=-60:predia=-60:cbp:vb_strategy=2:bidir_refine=4:vqmin=2"
elif [[ $# -eq 4 && $4 = vcd ]]; then
	LAVOPTS="vcodec=mpeg1video:vrc_buf_size=327:vrc_minrate=1152:vrc_maxrate=1152:vbitrate=1152:keyint=15:mbd=2:mv0:trell:mbcmp=2:precmp=2:subcmp=2:cmp=2:dia=-60:predia=-60:cbp:vb_strategy=2:bidir_refine=4:vqmin=2"
elif [[ $# -eq 4 && $4 = dvd ]]; then
	LAVOPTS="vcodec=mpeg1video:vrc_buf_size=1835:vrc_minrate=1152:vrc_maxrate=9800:vbitrate=5000:keyint=15:mbd=2:mv0:trell:mbcmp=2:precmp=2:subcmp=2:cmp=2:dia=-60:predia=-60:cbp:vb_strategy=2:bidir_refine=4:vqmin=2"
elif [[ $# -eq 4 && $4 = preview ]]; then
	LAVOPTS="vcodec=mpeg1video:vbitrate=1152:keyint=15"
else 
	LAVOPTS="vcodec=mpeg1video:vbitrate=1152:keyint=15"
fi

MENCODER="/usr/bin/mencoder"
MPEGOPTS="format=mpeg1"

$MENCODER mf://${1} -mf type=${2}:fps=25 -of rawvideo -mpegopts $MPEGOPTS -ovc lavc -lavcopts $LAVOPTS -o $3
$MENCODER  -oac copy -ovc copy -o repeat.mpg  $3 $3 $3 $3 $3 $3 $3
