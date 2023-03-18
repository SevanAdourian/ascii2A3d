#!/bin/bash

function parse_yaml {
   local prefix=$2
   local s='[[:space:]]*' w='[a-zA-Z0-9_]*' fs=$(echo @|tr @ '\034')
   sed -ne "s|^\($s\):|\1|" \
        -e "s|^\($s\)\($w\)$s:$s[\"']\(.*\)[\"']$s\$|\1$fs\2$fs\3|p" \
        -e "s|^\($s\)\($w\)$s:$s\(.*\)$s\$|\1$fs\2$fs\3|p"  $1 |
   awk -F$fs '{
      indent = length($1)/2;
      vname[indent] = $2;
      for (i in vname) {if (i > indent) {delete vname[i]}}
      if (length($3) > 0) {
         vn=""; for (i=0; i<indent; i++) {vn=(vn)(vname[i])("_")}
         printf("%s%s%s=\"%s\"\n", "'$prefix'",vn, $2, $3);
      }
   }'
}

eval $(parse_yaml $1)

# Output_file
model_aggregated_file=samples_174_win.dat
radial_grid_file=sample_radii.dat
lateral_grid_file=sample_locations.dat
flag=1

if [ -e $model_aggregated_file ]; then
    rm $model_aggregated_file
fi

if [ -e $lateral_grid_file ]; then
    rm $lateral_grid_file
fi

if [ -e $radial_grid_file ]; then
    rm $radial_grid_file
fi

# List all file for the sample
# for f in *.xyz; do
for f in $(ls -1 $path_to_ascii_files/merged*.xyz | sort -r); do
    npts=`awk 'NR==1 {print $2}' $f`
    depth=`echo $f | awk -F. '{print $1}' | grep -o '....$'`
    radius=`echo 6371-$depth | bc -l`
    echo $f $depth $radius
    if [ $flag -eq 1 ]
    then
    	awk 'NR > 1 {print $1"\t"$2}' $f > $lateral_grid_file	
    	flag=0
    fi
    echo $radius >> $radial_grid_file
    awk 'NR > 1 {print $3}' $f >> $model_aggregated_file
    # awk 'NR > 5 {print 1.0}' $f >> $model_aggregated_file
done
