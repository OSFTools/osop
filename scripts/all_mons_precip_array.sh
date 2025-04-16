#!/bin/bash -l
#SBATCH --qos=normal
#SBATCH --mem=128G
#SBATCH --ntasks=8
#SBATCH --output=master_%A_%a.out
#SBATCH --error=master_%A_%a.err
#SBATCH --time=00-06:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=ALL
# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of osop and is released under the BSD 3-Clause license.
# See LICENSE in the root of the repository for full licensing details.

# Script to calculate download hindcasts, calculate terciles and plot verification measures.
# use sbatch job array to submit 12 months at one time
# sbatch --array=1-12 -N1 all_mons_precip_array.sh
set -eu

# this conda env gives an error on load, so
# can't use -u option
set +u
conda activate osop 
set -u

# pick download location
downloaddir=$SCRATCH/seafoam/data/master
mkdir -p $downloaddir

# set PYTHONPATH relative to this location
lib_path=$(pushd ./../lib > /dev/null && pwd && popd > /dev/null)
set +u
export PYTHONPATH=$PYTHONPATH:$lib_path
set -u

# set parameters
leads="2,3,4" # e.g. if month=5 and leads="2,3,4", valid months are JJA (6,7,8)
area="45,-30,-2.5,60" # sub-area in degrees for area of interest (comma separated N,W,S,E)
# variable of interest, typically "2m_temperature" or "total_precipitation"
variable="total_precipitation" 

# set month from SLURM_ARRAY_TASK_ID
month=$SLURM_ARRAY_TASK_ID
echo month: $month
# get ERA5 data
set +e
python get_era5.py \
    --month $month \
    --leads $leads \
    --area $area \
    --downloaddir $downloaddir \
    --variable $variable \
    > $downloaddir/era5_log_${variable}_${month}.txt 2>&1
exitcode=$?
set -e
if [ $exitcode -eq 0 ]; then
    echo era5 downloaded
else
    echo era5 download failed
fi

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
        > $downloaddir/download_log_${variable}_${centre}_${month}.txt 2>&1
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
        > $downloaddir/product_log_${variable}_${centre}_${month}.txt 2>&1
    exitcode=$?
    set -e
    if [ $exitcode -eq 0 ]; then
        echo $centre : products generated
    else
        echo $centre : product generation failed
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
        > $downloaddir/verification_log_${variable}_${centre}_${month}.txt 2>&1
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
        > $downloaddir/plot_log_${variable}_${centre}_${month}.txt 2>&1
    exitcode=$?
    set -e
    if [ $exitcode -eq 0 ]; then
        echo $centre : plots generated
    else
        echo $centre : plot generation failed
    fi
done