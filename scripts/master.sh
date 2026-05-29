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

test=0
while getopts ":t" option; do
   case $option in
      t) # display Help
         test=1;;
     \?) # Invalid option
         echo "Error: Invalid keyword option, must be -t for test version or no option for full version"
         exit;;
   esac
done

# Script to calculate download hindcasts, calculate terciles and plot verification measures.
set -e

# this conda env gives an error on load, so
# can't use -u option
conda activate osop-simple 
set -u

# pick download location
downloaddir=$SCRATCH/seafoam/data/master/hindcast/downloads
productsdir=$SCRATCH/seafoam/data/master/hindcast/products
scoresdir=$SCRATCH/seafoam/data/master/hindcast/scores
plotdir=$SCRATCH/seafoam/data/master/hindcast/plots
logdir=$SCRATCH/seafoam/data/master/hindcast/logfiles
pycptdir=$SCRATCH/seafoam/data/master/hindcast/pycpt
mkdir -p $downloaddir
mkdir -p $plotdir
mkdir -p $logdir
mkdir -p $productsdir
mkdir -p $scoresdir
mkdir -p $pycptdir

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
area="55,-90,30,-60" # sub-area in degrees for area of interest (comma separated N,W,S,E)
variable="total_precipitation" # variable of interest, typically "2m_temperature" or "total_precipitation"
location="Morocco" #Current options include 'None' - no borders, 'UK','Morocco' and 'SAU' - Saudi Arabia
method="pmesh" #Remove for smooth plotting on correlation plots
pycpt="True" #True or False --> True you want pycpt, auto sets to off
predictor_area="45,-30,-2.5,60" #gcm area for predictor - if pycpt set to off, ignores (N,W,S,E)


# for the test version only run two models and get mme - ukmo
if [ $test -eq 1 ]; then
    centres="meteo_france ukmo mme"
    cat <<EOF > "$parseyml"
    Services:
        meteo_france: [9,1]
        jma: [3,0]  #need to leave this in as deleted currently and causes error if not included in ymls
        ukmo: [604,1]
        mme: [1,0]
EOF
else
    centres="meteo_france dwd cmcc ncep ukmo ecmwf jma eccc mme"
    # Services in use:
    # First column service, second column weight
    # mme weight should be set to 0, 1 on all other for equal weights
    # JMA set to 0 until regridding issue resolved

    cat <<EOF > "$parseyml"
    Services:
        ecmwf: [51,1]
        meteo_france: [9,1]
        dwd: [22,1]
        cmcc: [35,1]
        ncep: [2,1]
        jma: [3,0]  
        eccc_can: [4,1]
        eccc_gem5: [5,1]
        ukmo: [604,1]
        mme: [1,0]
EOF
fi
echo "YML file created: $parseyml"

# get ERA5 data
set +e
python get_era5.py \
    --month $month \
    --leads $leads \
    --area $area \
    --downloaddir $downloaddir \
    --logdir $logdir \
    --variable $variable \
    --pycpt $pycpt \
    --pycptdir $pycptdir
exitcode=$?
set -e
if [ $exitcode -eq 0 ]; then
    echo era5 downloaded
else
    echo era5 download failed
fi


#
# loop over all centres of interest and get data #for centre in meteo_france dwd cmcc ncep ukmo ecmwf jma eccc mme ;do 
for centre in $centres ;do  #meteo_france dwd cmcc ncep ukmo ecmwf jma eccc mme
    if [ "$centre" != "mme" ]; then
        set +e
        python get_any_hindcast.py \
            --centre $centre \
            --month $month \
            --leads $leads \
            --area $area \
            --variable $variable\
            --downloaddir $downloaddir \
            --logdir $logdir \
            --predictor_area $predictor_area \
            --pycpt $pycpt \
            --pycptdir $pycptdir
        exitcode=$?
        set -e
        if [ $exitcode -eq 0 ]; then
            echo $centre : download successful
        else
            echo $centre : download failed
            continue
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
        --logdir $logdir \
        --pycpt $pycpt \
        --pycptdir $pycptdir \
        --predictor_area $predictor_area 
    exitcode=$?
    set -e
    if [ $exitcode -eq 0 ]; then
        echo $centre : products generated
    else
        echo $centre : product generation failed
        continue
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
        continue
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
        continue
    fi
done
echo DONE