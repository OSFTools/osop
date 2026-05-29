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
## NOTE: TO SUCCESSFULLY RUN THIS SCRIPT; MASTER.SH MUST BE RUN IN FULL FIRST WITH THE SAME "SET PARAMETERS".

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

set -eu

# this conda env gives an error on load, so
# can't use -u option
set +u
conda activate osop-simple 
set -u

# pick download location
downloaddir=$SCRATCH/seafoam/data/master/forecast/downloads
productsdir=$SCRATCH/seafoam/data/master/forecast/products
scoresdir=$SCRATCH/seafoam/data/master/forecast/scores
plotdir=$SCRATCH/seafoam/data/master/forecast/plots
logdir=$SCRATCH/seafoam/data/master/forecast/logfiles
pycptdir=$SCRATCH/seafoam/data/master/forecast/pycpt
mkdir -p $downloaddir
mkdir -p $plotdir
mkdir -p $logdir
mkdir -p $productsdir
mkdir -p $scoresdir
mkdir -p $pycptdir

#locations for hindcast relatives
productshcdir=$SCRATCH/seafoam/data/master/hindcast/products
downloadhcdir=$SCRATCH/seafoam/data/master/hindcast/downloads

# set PYTHONPATH relative to this location
lib_path=$(pushd ./../lib > /dev/null && pwd && popd > /dev/null)
export PYTHONPATH=${PYTHONPATH:+$PYTHONPATH:}$lib_path

#create a yml file to pass dictionary parameters (being used for centers/systems)
parseyml="$downloaddir/parseyml.yml"


## NOTE: TO SUCCESSFULLY RUN THIS SCRIPT; MASTER.SH MUST BE RUN IN FULL FIRST WITH THE SAME "SET PARAMETERS".

# set parameters 
month=5 # initialisation month
leads="2,3,4" # e.g. if month=5 and leads="2,3,4", valid months are JJA (6,7,8)
area="55,-90,30,-60" # sub-area in degrees for area of interest (comma separated N,W,S,E) #"45,-30,-2.5,60" 
variable="total_precipitation" # variable of interest, typically "2m_temperature" or "total_precipitation"
location="None" #Current options include 'None' - no borders, 'UK','Morocco' and 'SAU' - Saudi Arabia
years=2025
pycpt="True" #True or False --> True you want pycpt, auto sets to off
predictor_area="45,-30,-2.5,60" #gcm area for predictor - if pycpt set to off, ignores (N,W,S,E)


# for the test version only run two models and get mme - ukmo
if [ $test -eq 1 ]; then
    centres="meteo_france ukmo mme"
    cat <<EOF > "$parseyml"
    Services:
        meteo_france: [9,1]
        jma: [3,0]  #need to leave this in as deleted currently and causes error if not included in yml
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

# loop over all centres of interest and get data
for centre in $centres  ;do  
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
            --years $years \
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
        fi
        # calculate products 
    fi
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
        --logdir $logdir \
        --predictor_area $predictor_area \
        --pycpt $pycpt \
        --pycptdir $pycptdir
    exitcode=$?
    set -e
    if [ $exitcode -eq 0 ]; then
        echo $centre : forecast products generated
    else
        echo $centre : forecast products failed -  check master.sh has been run correctly
    fi
    set +e
    python forecast_plots.py \
        --location $location \
        --centre $centre \
        --month $month \
        --variable $variable \
        --leads $leads \
        --area $area \
        --downloaddir $downloaddir \
        --productsfcdir $productsdir \
        --plotsdir $plotdir \
        --yearsfc $years \
        --logdir $logdir
    exitcode=$?
    set -e
    if [ $exitcode -eq 0 ]; then
        echo $centre : forecast plots generated
    else
        echo $centre : forecast plots failed 
    fi
done
echo DONE
