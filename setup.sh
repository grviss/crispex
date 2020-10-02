#!/bin/bash
# Get the version information accessible to crispex

git_id=`git describe --tags --long`
last_commit=`git log -1 --date=format:"%Y/%m/%d %T" --format="%ad %an"`

echo "$git_id ($last_commit)" > resources/versioninfo.txt
echo "crispex setup: done!"
