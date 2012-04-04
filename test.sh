#!/bin/sh
ls ../data/h264/*.264 | parallel src/app/h264dec -O
for f in *.out; do echo $f; cmp ../data/h264/${f%out}md5 $f; done
