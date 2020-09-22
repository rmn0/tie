
# Tile Image Encoder

## Building

    gcc -O2 -lgif tie.c -o tie

## Example Usage

    ffmpeg -i video.mp4 -vf "scale=256:224,split[s0][s1];[s0]palettegen=max_colors=17:reserve_transparent=1:stats_mode=full[p];[s1][p]paletteuse=dither=bayer:new=0" -loop 0 video.gif
    ./tie -c 50 -i video.gif -o video.tie 
    ./tie -d video.tie video.decoded.gif

Input videos cannot be not larger than 256x256.

## Options

**-i** *input filename*
Set input filename (required).

**-o** *output filename*
Set output filename (required).

**-f** *frames*
Set number of frames to encode / decode (default: number of frames in input file).

**-c** *compression level*
Compression level, from 0 to 100 (default 50).

**-v**
Be more verbose.

**-d**
Convert back to gif.

