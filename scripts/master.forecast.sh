#!/bin/bash -l
#SBATCH --qos=normal
#SBATCH --mem=128G
#SBATCH --ntasks=8
#SBATCH --output=master.txt
#SBATCH --error=master.err
#SBATCH --time=00-06:00:00
#SBATCH --export=NONE
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
conda activate osop 
set -u

# pick download location
downloaddir=$SCRATCH/seafoam/data/master/forecast/downloads
productsdir=$SCRATCH/seafoam/data/master/forecast/products
scoresdir=$SCRATCH/seafoam/data/master/forecast/scores
plotdir=$SCRATCH/seafoam/data/master/forecast/plots
logdir=$SCRATCH/seafoam/data/master/forecast/logfiles
mkdir -p $downloaddir
mkdir -p $plotdir
mkdir -p $logdir
mkdir -p $productsdir
mkdir -p $scoresdir

#locations for hindcast relatives
productshcdir=$SCRATCH/seafoam/data/master/hindcast/products
downloadhcdir=$SCRATCH/seafoam/data/master/hindcast/downloads

# set PYTHONPATH relative to this location
lib_path=$(pushd ./../lib > /dev/null && pwd && popd > /dev/null)
export PYTHONPATH=${PYTHONPATH:+$PYTHONPATH:}$lib_path

#create a yml file to pass dictionary parameters (being used for centers/systems)
parseyml="$downloaddir/parseyml.yml"

# set parameters 
month=5 # initialisation month
leads="2,3,4" # e.g. if month=5 and leads="2,3,4", valid months are JJA (6,7,8)
area="45,-30,-2.5,60" # sub-area in degrees for area of interest (comma separated N,W,S,E)
variable="total_precipitation" # variable of interest, typically "2m_temperature" or "total_precipitation"
location="Morocco" #Current options include 'None' - no borders, 'UK','Morocco' and 'SAU' - Saudi Arabia
years=2025


# Services in use: 
# edit as approprite to the most up to date systems. 
cat <<EOF > "$parseyml"
Services:
    ecmwf: 51
    meteo_france: 9
    dwd: 22
    cmcc: 35
    ncep: 2
    jma: 3
    eccc_can: 4
    eccc_gem5: 5
    ukmo: 604
EOF
echo "YML file created: $parseyml"

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
        --years $years \
        > $logdir/download_log_${variable}_${centre}.txt 2>&1
    exitcode=$?
    set -e
    if [ $exitcode -eq 0 ]; then
        echo $centre : download sucessful
    else
        echo $centre : download failed
    fi
    # calculate verification scores
    set +e
    python forecast_products.py \
        --centre $centre \
        --month $month \
        --variable $variable \
        --leads $leads \
        --area $area \
        --downloaddir $downloaddir \
        --downloadhcdir $downloadhcdir \
        --productshcdir $productshcdir \
        --productsfcdir $productsdir \
        --yearsfc $years \
        > $logdir/products_log_${variable}_${centre}.txt 2>&1
    exitcode=$?
    set -e
    if [ $exitcode -eq 0 ]; then
        echo $centre : forecast products generated
    else
        echo $centre : forecast products failed
    fi
done   
