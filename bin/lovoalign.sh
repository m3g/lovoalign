#!/bin/bash
#
# Script to run lovoalign using an input file and VMD style selections
# Requires VMD (http://www.ks.uiuc.edu/Research/vmd/)
#
#  Run with: ./lovoalign.sh input.dat 
#
#  L. Martinez, Institute of Chemistry, State University of Campinas
#               First version: Feb 26, 2009
#
# More information: http://www.ime.unicamp.br/~martinez/lovoalign
#
# IMPORTANT:
# Path for lovoalign program:
#
lovoalign=./lovoalign
#
# Path to temporary files
tmp=.
#
# Checking if the lovoalign executable is available

if [ ! -f $lovoalign ]; then

# Check if lovoalign is in the path

  lovoalign=`which lovoalign`

  if [ ! -e $lovoalign ]; then

    echo " ERROR: The lovoalign executable is not in the specified "
    echo "        path. Modify the path to the executable in the "
    echo "        lovoalign.sh script or add lovoalign to your path. "
    exit

  fi
fi

# Read lovoalign.sh input file and extract relevant data

file=$1

#Check if the file exists 

if [ ! -e $file ]; then
  echo " ERROR: Could not find file: $file "
  exit
fi

gap=1.e30
score=2
while read line; do
  keyword=`echo $line | awk '{print $1}'`
  case "$keyword" in
    proteinA) proteinA=`echo $line | cut -d' ' -f 2` ;;
    proteinB) proteinB=`echo $line | cut -d' ' -f 2` ;;
    selA) selA=`echo $line | cut -d' ' -f 2-` ;;
    selB) selB=`echo $line | cut -d' ' -f 2-` ;;
    output) output=`echo $line | cut -d' ' -f 2` ;;
    gap) gap=`echo $line | cut -d' ' -f 2` ;;
    score) score=`echo $line | cut -d' ' -f 2` ;;
    tmp) tmp=`echo $line | cut -d' ' -f 2` ;; 
    seqoff) seqoff=`echo $line | cut -d' ' -f 2` ;; 
  esac
done < <(cat $file) 

# Set unique tmp folder

tmp=$tmp/tmp_lovoalign_$RANDOM
mkdir $tmp

# Write VMD input file

echo "  ####################################################"
echo " "

for file in A B; do

if [ "$file" == "A" ]; then
  prot=$proteinA
  tempname=A
  sel=$selA
fi
if [ "$file" == "B" ]; then
  prot=$proteinB
  tempname=B
  sel=$selB
fi

vmdfile=lovoalignvmd.temp1
selection=lovoalignvmd.$tempname
echo "
set prot $prot
set select \"$sel\"

mol new \$prot type pdb first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all
 
# Clear beta
set notsel [atomselect top all frame first]
\$notsel set beta 0
\$notsel set occupancy 0
 
# Add 1.00 to beta field of the solute
set sel1 [ atomselect top \$select ]
\$sel1 set beta 1
 
# Write temporary pdb file with beta values
\$notsel writepdb $tmp/$selection
 
exit " > $tmp/$vmdfile  

# Run VMD to build a temporary pdb file with group definitions

echo "    Running VMD to define selections for protein $tempname " 

vmd=`which vmd`
vmdlog=`$vmd -dispdev text < $tmp/$vmdfile`
IFS=$'\n'
for line in $vmdlog; do
 error=`echo $line |grep ERROR` 
 if [ "$error" \> " " ]; then 
   echo $error
   rm -rf $tmp
   exit
 fi
done

done
echo " "
echo "  ####################################################"

# Run lovoalign

if [ $seqoff == "true" ]; then

  $lovoalign -p1 $tmp/lovoalignvmd.A -beta1 \
             -p2 $tmp/lovoalignvmd.B -beta2 \
             -m $score -g $gap \
             -all -seqoff \
             -o $output  

else

  $lovoalign -p1 $tmp/lovoalignvmd.A -beta1 \
             -p2 $tmp/lovoalignvmd.B -beta2 \
             -m $score -g $gap \
             -all \
             -o $output  

fi

# Removing temporary files

rm -rf $tmp

