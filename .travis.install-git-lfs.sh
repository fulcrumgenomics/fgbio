#!/bin/bash

GITHUB_ACCESS_TOKEN=$1;

GIT_LFS_VERSION="2.0.2"
GIT_LFS_LINK=https://github.com/github/git-lfs/releases/download/v${GIT_LFS_VERSION}/git-lfs-linux-amd64-${GIT_LFS_VERSION}.tar.gz
GIT_LFS="git-lfs-${GIT_LFS_VERSION}/git-lfs"
echo "Downloading and untarring git-lfs binary" 
wget -qO- $GIT_LFS_LINK | tar xvz
git_dir=`dirname $(which git)`
cp $GIT_LFS $git_dir/

echo "Resetting remote URL"
git remote set-url origin "https://${GITHUB_ACCESS_TOKEN}@github.com/fulcrumgenomics/fgbio.git"

echo "Fetching LFS files."
$GIT_LFS install
$GIT_LFS pull
