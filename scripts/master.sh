#!/bin/bash -l
#SBATCH --mail-user=emma.dyer@metoffice.gov.uk
#SBATCH --mail-type=ALL
#SBATCH --qos=normal
#SBATCH --mem=64G
#SBATCH --ntasks=8
#SBATCH --output=master_batch.txt
#SBATCH --error=master_batch.err
#SBATCH --time=00-06:00:00
#SBATCH --export=NONE
# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of osop and is released under the BSD 3-Clause license.
# See LICENSE in the root of the repository for full licensing details.

# Script to calculate verifiction  metrics for 
# ArabCOF region. Hindcasts initialised May for JJA
set -eu

# this conda env gives an error on load, so
# can't use -u option
set +u
conda activate climada_2024 ## to do build osop env, error using conda main
set -u

# pick download location
downloaddir=$SCRATCH/seafoam/data/master

# set parameters
month=5
leads="2,3,4"
lead_obs="1,2,3"
area="45,-30,-2.5,60"
variable="2m_temperature"
aggregation="3m"

# loop over all centres of interest and get data
for centre in meteo_france dwd cmcc ncep ukmo  ecmwf jma eccc ;do 
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
        --leads_obs $lead_obs \
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
        --leads_obs $lead_obs \
        --area $area \
        --aggregation $aggregation \
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
        --aggregation $aggregation \
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
