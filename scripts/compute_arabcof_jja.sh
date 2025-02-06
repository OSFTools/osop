#!/bin/bash -l
#SBATCH --mail-user=emma.dyer@metoffice.gov.uk
#SBATCH --mail-type=ALL
#SBATCH --qos=normal
#SBATCH --mem=64G
#SBATCH --ntasks=8
#SBATCH --output=product_batch.txt
#SBATCH --error=product_batch.err
#SBATCH --time=00-06:00:00
#SBATCH --export=NONE
# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of osop and is released under the BSD 3-Clause license.
# See LICENSE in the root of the repository for full licensing details.

# Script to calculate anomalies and terciles for 
# ArabCOF region. Hindcasts initialised May for JJA
set -eu

# this conda env gives an error on load, so
# can't use -u option
set +u
conda activate climada_2024 ## to do build osop env, error using conda main
set -u

# pick download location
downloaddir=$SCRATCH/seafoam/data

# loop over all centres of interest and get data
for centre in meteo_france dwd cmcc ncep ukmo  ecmwf jma eccc ;do 
    set +e
    python compute_products.py \
        --centre $centre \
        --month 5 \
        --leads "2,3,4" \
        --area "45,-30,-2.5,60" \
        --downloaddir $downloaddir \
        > $downloaddir/product_log_jja_${centre}.txt 2>&1
    exitcode=$?
    set -e
    if [ $exitcode -eq 0 ]; then
        echo $centre : products generated
    else
        echo $centre : product generation failed
    fi
done
