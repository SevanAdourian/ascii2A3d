#! /bin/bash

infile=$1
# range="-Rd"
range="-R160/270/30/85"
frame="-JR230/20c"
gmt makecpt -Cpolar -D -I -Z -T-0.1/0.1/0.01 > tomo.cpt
symb="-Sc0.2c -Ctomo.cpt"
pen="-W0.03c,black"
bor="-B45g45WSen"
interp="-I0.638"
grout="xyz_plot.grd"
outps="xyz_plot2.ps"

gmt surface $infile $range $interp -G$grout
gmt grdimage $range $frame $bor $grout -Ctomo.cpt -K -Y4c > $outps

echo "plotting the velocity coefficients"
# gmt psxy $range $frame $bor $infile $symb -K -Y3c > $outps
echo "plotting the coast line"
gmt pscoast $range $frame -Dl -A5000 $pen -O -K  >> $outps
echo "plotting the scale"
gmt psscale -Ctomo.cpt -D10c/-1c/5c/0.25ch  -B2:"dVs %": -O >> $outps

# gv $outps &
