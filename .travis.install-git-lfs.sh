#!/bin/bash

GIT_LFS_VERSION="1.5.4"
GIT_LFS_LINK=https://github.com/github/git-lfs/releases/download/v${GIT_LFS_VERSION}/git-lfs-linux-amd64-${GIT_LFS_VERSION}.tar.gz
GIT_LFS="git-lfs-${GIT_LFS_VERSION}/git-lfs"
echo "Downloading and untarring git-lfs binary" 
wget -qO- $GIT_LFS_LINK | tar xvz

echo "Resetting remote URL"
git remote set-url origin "https://github.com/fulcrumgenomics/fgbio.git"

echo "Fetching LFS files."
$GIT_LFS install
$GIT_LFS pull