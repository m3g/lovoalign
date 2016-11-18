#!/bin/bash
####################################################################################################

# Software name:

package=lovoalign

# HTML file containing download list:

downloads=~/public_html/lovoalign/versionhistory/downloads.html

# GIT URL:

giturl=https://github.com/leandromartinez98/lovoalign/archive/

# Name of file containing version number

versionfile=./title.f90

####################################################################################################

year=`date +%y`
day=`date +%j`
version="${year:0:1}${year:1:1}.$day"

file="$package-$version.tar.gz" 
version_file=$version

i=2
while grep -q $file $downloads ; do
  file=$package-$version.$i.tar.gz
  version_file=$version.$i
  i=`expr $i + 1`
done
version=$version_file
file=$package-$version.tar.gz
echo "Will create file: $file"

cat $versionfile | sed -e "s/Version.*/Version\ $version \',\/\&/" > version_title_temp.f90
\mv -f version_title_temp.f90 $versiontitle

git add -A .
git commit -m "Changed version file to $version"
git tag -a $version -m "Release $version"
git push --all --tags

newline="<tr><td width=190px valign=top><a href=$giturl/$version.tar.gz> $filename</a></td><td> Release $version </td></tr>"
htmlfile=$downloads

writeline=yes
while IFS= read -r line ; do
 
  if [ "$writeline" = "yes" ] ; then
    echo $line >> ./htmlfile_new_temp
  fi
  if [ $writeline = "no" ] ; then
    writeline="yes"
  fi

  if [[ $line == *NEW_VERSION_HERE* ]] ; then
    echo $newline >> ./htmlfile_new_temp
  fi

done < "$htmlfile"
rm $htmlfile
mv htmlfile_new_temp $htmlfile   

echo " Done. " 

