#!/bin/bash -l
# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of osop and is released under the BSD 3-Clause license.
# See LICENSE in the root of the repository for full licensing details.

# Script to download common hindcast period (1993-2016) for 
# ArabCOF region. Hindcasts initialised May for JJA
set -eu

# this conda env gives an error on load, so
# can't use -u option
set +u
conda activate osop
set -u

# pick download location
downloaddir=$SCRATCH/test_downloads
mkdir -p $downloaddir

# loop over all centres of interest and get data
for centre in meteo_france dwd cmcc ncep ukmo  ecmwf jma eccc ;do 
    set +e
    python get_any_hindcast.py \
        --centre $centre \
        --month 5 \
        --leads "2,3,4" \
        --area "45,-30,-2.5,60" \
        --downloaddir $downloaddir \
        > $downloaddir/download_log_jja_${centre}.txt 2>&1
    exitcode=$?
    set -e
    if [ $exitcode -eq 0 ]; then
        echo $centre : download sucessful
    else
        echo $centre : download failed
    fi
done
