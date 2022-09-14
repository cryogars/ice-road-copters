#!/bin/bash

ASP_DIR='/SNOWDATA/IDALS/ASP/bin/'
MCS_SHP='/SNOWDATA/IDALS/ice-road-copters/transform_area/hwy_21/hwy_21_utm_edit_v2.shp'
DC_SHP='/SNOWDATA/IDALS/ice-road-copters/transform_areaBogusBasin/Roads_DryCreek.shp'
MCS_DEM='/SNOWDATA/IDALS/REF_DEM/MCS_REFDEM_WGS84.tif'
DC_DEM='/SNOWDATA/IDALS/REF_DEM/DC_REFDEM_WGS84.tif'

ICEROAD='/SNOWDATA/IDALS/ice-road-copters/scripts/ice-road-pipeline.py'

D2021='/SNOWDATA/IDALS/2021'
D2022='/SNOWDATA/IDALS/2022'
SCRATCH='/home/zacharykeskinen/scratch'
echo "Starting 2021 Flights..."

for FP in $(ls $D2021)
do
    $(cp -r $D2021/$FP $SCRATCH)
    if [[ $FP == *"MCS"* ]]; then
        MCS_CMD="python $ICEROAD $SCRATCH/$FP -e $MCS_DEM -s $MCS_SHP -a $ASP_DIR"
        echo $MCS_CMD
        $($MCS_CMD)
    fi

    if [[ $FP == *"DC"* ]]; then
        DC_CMD="python $ICEROAD $SCRATCH/$FP -e $DC_DEM -s $DC_SHP -a $ASP_DIR"
        echo $DC_CMD
        $($DC_CMD)
    fi
    echo $?
    $(cp -r $SCRATCH/$FP/ice-road $D2021/$FP)
    $(rm -r $SCRATCH/$FP)
done

echo $'\nStarting 2022 Flights...'

for FP in $(ls $D2022)
do
    $(cp -r $D2022/$FP $SCRATCH)
    if [[ $FP == *"MCS"* ]]; then
        MCS_CMD="python $ICEROAD $SCRATCH/$FP -e $MCS_DEM -s $MCS_SHP -a $ASP_DIR"
        echo $MCS_CMD
        $($MCS_CMD)
    fi

    if [[ $FP == *"DC"* ]]; then
        DC_CMD="python $ICEROAD $SCRATCH/$FP -e $DC_DEM -s $DC_SHP -a $ASP_DIR"
        echo $DC_CMD
        $($DC_CMD)
    fi
    echo $?
    $(cp -r $SCRATCH/$FP/ice-road $D2022/$FP/)
    $(rm -r $SCRATCH/$FP)
done