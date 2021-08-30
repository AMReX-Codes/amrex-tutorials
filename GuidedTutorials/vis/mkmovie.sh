#!/usr/bin/bash

for f in $(ls -d plt?????); do
    python3 slice.py $f -f u -ax z -min 0.0 -max 3.5 &
done

wait

python3 ffmpeg_make_mp4 plt*.png -s -ifps 5 > /dev/null 2>&1
