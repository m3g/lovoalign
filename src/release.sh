#!/bin/bash
####################################################################################################

# Software name:

package=lovoalign

# HTML file containing download list:

downloads=~/public_html/lovoalign/versionhistory/downloads.html

# GIT URL:

giturl=https://github.com/leandromartinez98/lovoalign

# Name of file containing version number

versionfile=./title.f90

####################################################################################################

# Usage: ./release.sh <major|minor|patch|X.Y.Z>
#   major|minor|patch : bump the corresponding field of the current version
#   X.Y.Z             : set the version explicitly
# The version follows semantic versioning (https://semver.org) and the git
# tag is the version prefixed with "v" (e.g. v22.0.3).

current=`grep -o "Version [0-9][0-9.]*" $versionfile | awk '{print $2}'`
if [[ ! $current =~ ^[0-9]+\.[0-9]+\.[0-9]+$ ]] ; then
  echo "Current version in $versionfile ($current) is not in X.Y.Z format."
  echo "Provide the new version explicitly: ./release.sh X.Y.Z"
  current=""
fi
IFS=. read -r major minor patch <<< "$current"

case "$1" in
  major) [ -z "$current" ] && exit 1 ; version="$((major+1)).0.0" ;;
  minor) [ -z "$current" ] && exit 1 ; version="$major.$((minor+1)).0" ;;
  patch) [ -z "$current" ] && exit 1 ; version="$major.$minor.$((patch+1))" ;;
  *)
    if [[ $1 =~ ^[0-9]+\.[0-9]+\.[0-9]+$ ]] ; then
      version=$1
    else
      echo "Usage: ./release.sh <major|minor|patch|X.Y.Z>   (current version: ${current:-unknown})"
      exit 1
    fi
    ;;
esac

tag="v$version"
if git rev-parse -q --verify "refs/tags/$tag" > /dev/null ; then
  echo "Tag $tag already exists."
  exit 1
fi

file=$package-$version.tar.gz
echo "Current version: ${current:-unknown}"
echo "New version: $version (tag $tag)"
echo "Will create file: $file"
read -p "Proceed? [y/N] " answer
[[ $answer =~ ^[Yy]$ ]] || exit 1

cat $versionfile | sed -e "s/Version [0-9][0-9.]* /Version $version /" > version_title_temp.f90
\mv -f version_title_temp.f90 $versionfile

git add -A .
git commit -m "Changed version file to $version"
git tag -a $tag -m "Release $version"
git push origin master tag $tag

today=`date +"%b %d, %Y"`
changelog="https://github.com/leandromartinez98/$package/releases/tag/$tag"
newline="<tr><td width=190px valign=top><a href=$giturl/archive/$tag.tar.gz> $file </a></td><td> Released on $today - <a target=_blank_ href=$changelog> [change log at github]  </a></td></tr>"
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

echo "----------------------"
echo "CHANGE LOG:"
echo "----------------------"
range=`git tag --sort=v:refname -l 'v[0-9]*.[0-9]*.[0-9]*' | tail -n 2 | xargs | sed 's! !...!'`
git log --pretty=oneline $range | awk '{$1=""; print "-"$0}'
echo "----------------------"

echo " Done. " 

