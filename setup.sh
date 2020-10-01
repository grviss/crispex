#!/bin/bash
# Get the version information accessible to crispex

git_id=`git describe --tags --long`
last_commit_date=`git log -1 --date=format:"%Y/%m/%d %T" --format="%ad"`
last_commit_author=`git log -1 --format="%an"`

echo "$git_id ($last_commit_date $last_commit_author)" > resources/versioninfo.txt
echo "crispex setup: done!"
