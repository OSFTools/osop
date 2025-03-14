#!/bin/bash -l
#SBATCH --qos=normal
#SBATCH --mem=128G
#SBATCH --ntasks=8
#SBATCH --output=master.txt
#SBATCH --error=master.err
#SBATCH --time=00-06:00:00
#SBATCH --export=NONE
#SBATCH --mail-user=emma.dyer@metoffice.gov.uk
#SBATCH --mail-type=ALL
# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of osop and is released under the BSD 3-Clause license.
# See LICENSE in the root of the repository for full licensing details.

# Script to calculate download hindcasts, calculate terciles and plot verification measures.
set -eu

# this conda env gives an error on load, so
# can't use -u option
set +u
conda activate climada_2024 ## to do build osop env, error using conda main

export PYTHONPATH=${PYTHONPATH}:${HOME}/IASAS/osop/
set -u

# pick download location
downloaddir=$SCRATCH/seafoam

# set parameters
month=11 # initialisation month
leads="2,3,4" # e.g. if month=5 and leads="2,3,4", valid months are JJA (6,7,8)
area="45,-30,-2.5,60" # sub-area in degrees for area of interest (comma separated N,W,S,E)
variable="2m_temperature" # variable of interest, typically "2m_temperature" or "total_precipitation"

# loop over all centres of interest and get data
for centre in meteo_france dwd cmcc ncep ukmo ecmwf jma eccc ;do 
    set +e
    python get_any_hindcast.py \
        --centre $centre \
        --month $month \
        --leads $leads \
        --area $area \
        --variable $variable\
        --downloaddir $downloaddir \
        > $downloaddir/download_log_${variable}_${centre}.txt 2>&1
    exitcode=$?
    set -e
    if [ $exitcode -eq 0 ]; then
        echo $centre : download sucessful
    else
        echo $centre : download failed
    fi
    # compute terciles and anomalies
    set +e
    python compute_products.py \
        --centre $centre \
        --month $month \
        --leads $leads \
        --area $area \
        --variable $variable \
        --downloaddir $downloaddir \
        > $downloaddir/product_log_${variable}_${centre}.txt 2>&1
    exitcode=$?
    set -e
    if [ $exitcode -eq 0 ]; then
        echo $centre : products generated
    else
        echo $centre : product generation failed
    fi
    # get ERA5 data
    set +e
    python get_era5.py \
        --month $month \
        --leads $leads \
        --area $area \
        --downloaddir $downloaddir \
        --variable $variable \
        > $downloaddir/era5_log_${variable}.txt 2>&1
    exitcode=$?
    set -e
    if [ $exitcode -eq 0 ]; then
        echo $centre : era5 downloaded
    else
        echo $centre : era5 download failed
    fi
    # calculate verification scores
    set +e
    python compute_scores.py \
        --centre $centre \
        --month $month \
        --leads $leads \
        --area $area \
        --downloaddir $downloaddir \
        --variable $variable \
        > $downloaddir/verification_log_${variable}_${centre}.txt 2>&1
    exitcode=$?
    set -e
    if [ $exitcode -eq 0 ]; then
        echo $centre : scores generated
    else
        echo $centre : score generation failed
    fi
    # plot scores
        set +e
    python plot_verification.py \
        --centre $centre \
        --month $month \
        --leads $leads \
        --area $area \
        --downloaddir $downloaddir \
        --variable $variable \
        > $downloaddir/plot_log_${variable}_${centre}.txt 2>&1
    exitcode=$?
    set -e
    if [ $exitcode -eq 0 ]; then
        echo $centre : plots generated
    else
        echo $centre : plot generation failed
    fi
done
