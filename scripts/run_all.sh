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

# set parameters
do_hc=1 # 1 to run hindcast, 0 to skip. Need to have run at least once to calculate terciles for forecast to work
do_fc=0 # 1 to run forecast, 0 to skip. 
month=5 # initialisation month
leads="2,3,4" # e.g. if month=5 and leads="2,3,4", valid months are JJA (6,7,8)
area="39,60,-11,141" # sub-area in degrees for area of interest (comma separated N,W,S,E) 
variable="total_precipitation" # variable of interest, typically "2m_temperature" or "total_precipitation"
location="None" #Current options include 'None' - no borders, 'UK','Morocco' and 'SAU' - Saudi Arabia
method="pmesh" #Remove for smooth plotting on correlation plots
pycpt="True" #True or False --> True you want pycpt, auto sets to off
predictor_area="40,0,-40,359" #gcm area for predictor - if pycpt set to off, ignores (N,W,S,E)
fc_year=2025 #year to run forecast for

exp_name=single_script
# pick download location
base_path=$SCRATCH/osop/${exp_name}
logdir=${base_path}/logfiles

downloaddir=${base_path}/hindcast/downloads
productsdir=${base_path}/hindcast/products
scoresdir=${base_path}/hindcast/scores
plotdir=${base_path}/hindcast/plots
pycptdir=${base_path}/hindcast/pycpt

fc_downloaddir=${base_path}/forecast/downloads
fc_productsdir=${base_path}/forecast/products
fc_scoresdir=${base_path}/forecast/scores
fc_plotdir=${base_path}/forecast/plots
fc_pycptdir=${base_path}/forecast/pycpt

# end of user chosen options

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

# this conda env gives an error on load, so
# can't use -u option
set +u
conda activate osop
set -u

mkdir -p "$logdir"

mkdir -p "$downloaddir"
mkdir -p "$plotdir"
mkdir -p "$productsdir"
mkdir -p "$scoresdir"
mkdir -p "$pycptdir"

mkdir -p "$fc_downloaddir"
mkdir -p "$fc_productsdir"
mkdir -p "$fc_plotdir"
mkdir -p "$fc_scoresdir"
mkdir -p "$fc_pycptdir"


# set PYTHONPATH relative to this location
lib_path=$(pushd ./../lib > /dev/null && pwd && popd > /dev/null)
set +u
export PYTHONPATH="${PYTHONPATH:+$PYTHONPATH:}$lib_path"
set -u

#create a yml file to pass dictionary parameters
parseyml="$downloaddir/parseyml.yml"

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
    centres="meteo_france dwd cmcc ncep ukmo ecmwf jma eccc bom mme"
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
        bom: [2,1]
        mme: [1,0]
EOF
fi
echo "YML file created: $parseyml"

# currently forecast expects a copy of this in fc_downloaddir
cp "$parseyml" "$fc_downloaddir/$(basename "$parseyml")"

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
for centre in $centres ;do  
  if [ "$do_hc" -eq 1 ]; then
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
  fi
  if [ "$do_fc" -eq 1 ]; then
    # run forecast
    if [ "$centre" != "mme" ]; then
        set +e
        python get_any_hindcast.py \
            --centre $centre \
            --month $month \
            --leads $leads \
            --area $area \
            --variable $variable\
            --downloaddir $fc_downloaddir \
            --logdir $logdir \
            --years $fc_year \
            --predictor_area $predictor_area \
            --pycpt $pycpt \
            --pycptdir $fc_pycptdir
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
        --downloaddir $fc_downloaddir \
        --downloadhcdir $downloaddir \
        --productshcdir $productsdir \
        --productsfcdir $fc_productsdir \
        --yearsfc $fc_year \
        --logdir $logdir \
        --predictor_area $predictor_area \
        --pycpt $pycpt \
        --pycptdir $fc_pycptdir \
        --hindcast_pycptdir $pycptdir
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
        --downloaddir $fc_downloaddir \
        --productsfcdir $fc_productsdir \
        --plotsdir $fc_plotdir \
        --yearsfc $fc_year \
        --logdir $logdir
    exitcode=$?
    set -e
    if [ $exitcode -eq 0 ]; then
        echo $centre : forecast plots generated
    else
        echo $centre : forecast plots failed 
    fi
  fi  
done
echo DONE