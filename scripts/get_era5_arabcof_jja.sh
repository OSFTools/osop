#!/bin/bash -l
#SBATCH --mail-user=emma.dyer@metoffice.gov.uk
#SBATCH --mail-type=ALL
#SBATCH --qos=normal
#SBATCH --mem=64G
#SBATCH --ntasks=8
#SBATCH --output=era5_batch.txt
#SBATCH --error=era5_batch.err
#SBATCH --time=00-06:00:00
#SBATCH --export=NONE
# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of osop and is released under the BSD 3-Clause license.
# See LICENSE in the root of the repository for full licensing details.

# Script to download ERA5 for the common hindcast period (1993-2016) for 
# ArabCOF region. For JJA
set -eu

# this conda env gives an error on load, so
# can't use -u option
set +u
conda activate climada_2024 ## to do build osop env, error using conda main
set -u

# pick download location
downloaddir=$SCRATCH/seafoam/data

# get data
set +e
python get_era5.py \
    --month 5 \
    --leads_obs "1,2,3" \
    --area "45,-30,-2.5,60" \
    --downloaddir $downloaddir \
    --variable "2m_temperature" \
    > $downloaddir/era5_log_jja.txt 2>&1
exitcode=$?
set -e
if [ $exitcode -eq 0 ]; then
    echo $centre : era5 downloaded
else
    echo $centre : era5 download failed
fi
done
