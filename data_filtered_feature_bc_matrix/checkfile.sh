#!/bin/bash
find ./ -type f -print0 | xargs -0 md5sum > ./all.md5