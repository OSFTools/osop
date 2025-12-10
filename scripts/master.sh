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
downloaddir=$SCRATCH/seafoam/data/master/hindcast/downloads
productsdir=$SCRATCH/seafoam/data/master/hindcast/products
scoresdir=$SCRATCH/seafoam/data/master/hindcast/scores
plotdir=$SCRATCH/seafoam/data/master/hindcast/plots
logdir=$SCRATCH/seafoam/data/master/hindcast/logfiles
mkdir -p $downloaddir
mkdir -p $plotdir
mkdir -p $logdir
mkdir -p $productsdir
mkdir -p $scoresdir

# set PYTHONPATH relative to this location
lib_path=$(pushd ./../lib > /dev/null && pwd && popd > /dev/null)
set +u
export PYTHONPATH=$PYTHONPATH:$lib_path
set -u

#create a yml file to pass dictionary parameters
parseyml="$downloaddir/parseyml.yml"

# set parameters
month=5 # initialisation month
leads="2,3,4" # e.g. if month=5 and leads="2,3,4", valid months are JJA (6,7,8)
area="45,-30,-2.5,60" # sub-area in degrees for area of interest (comma separated N,W,S,E)
variable="2m_temperature" # variable of interest, typically "2m_temperature" or "total_precipitation"
location="Morocco" #Current options include 'None' - no borders, 'UK','Morocco' and 'SAU' - Saudi Arabia
method="pmesh" #Remove for smooth plotting on correlation plots

# Services in use:
cat <<EOF > "$parseyml"
Services:
    ecmwf: [51,1]
    meteo_france: [9,1]
    dwd: [22,1]
    cmcc: [35,1]
    ncep: [2,1]
    jma: [3,1]
    eccc_can: [4,1]
    eccc_gem5: [5,1]
    ukmo: [604,1]
    mme: [1,1]
EOF
echo "YML file created: $parseyml"

# get ERA5 data
set +e
python get_era5.py \
    --month $month \
    --leads $leads \
    --area $area \
    --downloaddir $downloaddir \
    --logdir $logdir \
    --variable $variable
exitcode=$?
set -e
if [ $exitcode -eq 0 ]; then
    echo era5 downloaded
else
    echo era5 download failed
fi

# loop over all centres of interest and get data #for centre in meteo_france dwd cmcc ncep ukmo ecmwf jma eccc mme ;do 
for centre in meteo_france dwd cmcc ncep ukmo ecmwf jma eccc mme ;do 
    if [ "$centre" != "mme" ]; then
        set +e
        python get_any_hindcast.py \
            --centre $centre \
            --month $month \
            --leads $leads \
            --area $area \
            --variable $variable\
            --downloaddir $downloaddir \
            --logdir $logdir
        exitcode=$?
        set -e
        if [ $exitcode -eq 0 ]; then
            echo $centre : download sucessful
        else
            echo $centre : download failed
        fi
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
        --productsdir $productsdir \
        --logdir $logdir
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
        --scoresdir $scoresdir \
        --productsdir $productsdir \
        --variable $variable \
        --logdir $logdir
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
        --location $location \
        --centre $centre \
        --month $month \
        --leads $leads \
        --area $area \
        --downloaddir $downloaddir \
        --scoresdir $scoresdir \
        --plotdir $plotdir \
        --variable $variable \
        --method $method \
        --logdir $logdir 
    exitcode=$?
    set -e
    if [ $exitcode -eq 0 ]; then
        echo $centre : plots generated
    else
        echo $centre : plot generation failed
    fi
done
