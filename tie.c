

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>

#include <gif_lib.h>

#define BITPLANES 4
#define COLORS (1 << BITPLANES)
#define MIRRORBIT 13
#define TILE_SIZE 8

#define TILE_PIXELS (TILE_SIZE * TILE_SIZE)
#define TILE_BYTES (BITPLANES * TILE_SIZE)

#define CELLS 16
#define FRAME 32
#define TILES_IN_FRAME (FRAME * FRAME)

#define VRAM_TILES 1024

#define NSTATS 6

enum update_code {
  NONE,
  UPDATE,
  CACHE,
  UPDATE_SKIP
};

struct tile {
  unsigned char byte[TILE_BYTES];
};

struct vram_info {
  struct tile tiles[4][VRAM_TILES * 2];
  unsigned char hash[4][VRAM_TILES * 2];
  int refcount[VRAM_TILES * 2];
  unsigned short tilemap[TILES_IN_FRAME * 2];
  int index;
};

struct encode_info {
  int frameno;

  struct vram_info vram, vram_prev;

  struct tile frame[TILES_IN_FRAME * 2];
  unsigned char update[TILES_IN_FRAME * 2];
  unsigned short skip[TILES_IN_FRAME * 2];
  
  int stats[NSTATS];

  unsigned char flag[TILES_IN_FRAME * 2];
};

unsigned char reverse_bits(unsigned char b)
{
  unsigned char r = 0;
  unsigned byte_len = 8;

  while (byte_len--) {
    r = (r << 1) | (b & 1);
    b >>= 1;
  }
  return r;
}

unsigned char reverse[256];

struct tie_header {
  unsigned char magic[4];
  unsigned short width;
  unsigned short height;
  unsigned short frames;
  unsigned short delay;
  unsigned char palette[COLORS * 3];
};

void swap(unsigned char *a, unsigned char *b) {
  unsigned char t = *a;
  *a = *b;
  *b = t;
}

struct tile mirror(struct tile t, const unsigned int m)
{
  if(m & 2) {
    for(int ii = 0; ii < TILE_SIZE / 2; ++ii) {
      swap(&t.byte[ii * 2 + 0], &t.byte[(TILE_SIZE - 1 - ii) * 2 + 0]);
      swap(&t.byte[ii * 2 + 1], &t.byte[(TILE_SIZE - 1 - ii) * 2 + 1]);
      swap(&t.byte[ii * 2 + TILE_SIZE * 2 + 0], &t.byte[(TILE_SIZE - 1 - ii) * 2 + TILE_SIZE * 2 + 0]);
      swap(&t.byte[ii * 2 + TILE_SIZE * 2 + 1], &t.byte[(TILE_SIZE - 1 - ii) * 2 + TILE_SIZE * 2 + 1]);
    }
  }

  if (m & 1) {
    for(int ii = 0; ii < TILE_BYTES; ++ii) t.byte[ii] = reverse[t.byte[ii]];
  }

  return t;
}

int sametile(const struct tile a, const struct tile b)
{
  for(int ii = 0; ii < TILE_BYTES; ++ii) if(a.byte[ii] != b.byte[ii]) return 0;
  return 1;
}

unsigned char tilehash(const struct tile a) {
  unsigned char hash = 0;
  for(int ii = 0; ii < TILE_BYTES; ++ii) hash += a.byte[ii];
  return hash;
}

unsigned int tilediff(const struct tile a, const struct tile b, int maxdiff)
{
  unsigned int count = 0;
  for(int ii = 0; ii < TILE_SIZE; ++ii) {
    for(int jj = 0; jj < TILE_SIZE; ++jj) {
      int bit = 1 << jj;
      int ac = (a.byte[ii * 2 + 0] & bit ? 1 : 0)
        | (a.byte[ii * 2 + 1] & bit ? 2 : 0)
        | (a.byte[ii * 2 + 16] & bit ? 4 : 0)
        | (a.byte[ii * 2 + 17] & bit ? 8 : 0);
      int bc = (b.byte[ii * 2 + 0] & bit ? 1 : 0)
        | (b.byte[ii * 2 + 1] & bit ? 2 : 0)
        | (b.byte[ii * 2 + 16] & bit ? 4 : 0)
        | (b.byte[ii * 2 + 17] & bit ? 8 : 0);

      count += ac - bc < 0 ? bc - ac : ac - bc;
    }
    if(count >= maxdiff) return maxdiff;
  }
  return count;
}

void init(struct encode_info* info)
{
  for(int ii = 0; ii < 256; ++ii) reverse[ii] = reverse_bits(ii);

  for(int ii = 0; ii < COLORS; ++ii) {
    for(int kk = 0; kk < 4; ++kk) {
      for(int jj = 0; jj < 8; ++jj) {
        info->vram.tiles[kk][ii].byte[jj * 2 + 0] = ii & 1 ? 0xff : 0x00;
        info->vram.tiles[kk][ii].byte[jj * 2 + 1] = ii & 2 ? 0xff : 0x00;
        info->vram.tiles[kk][ii].byte[jj * 2 + 16] = ii & 4 ? 0xff : 0x00;
        info->vram.tiles[kk][ii].byte[jj * 2 + 17] = ii & 8 ? 0xff : 0x00;
      }
    }
    info->vram.refcount[ii] = 1;
  }

  for(int ii = COLORS; ii < VRAM_TILES; ++ii) {
    for(int jj = 0; jj < TILE_BYTES; ++jj) {
      info->vram.tiles[0][ii].byte[jj] = 0;
    }
    info->vram.refcount[ii] = 0;
  }

  for(int kk = 0; kk < 4; ++kk) {
    for(int ii = 0; ii < VRAM_TILES; ++ii) {
      info->vram.hash[kk][ii] = tilehash(info->vram.tiles[kk][ii]);
    }
  }

  for(int ii = 0; ii < TILES_IN_FRAME; ++ii) {
    for(int jj = 0; jj < TILE_BYTES; ++jj) {
      info->frame[ii].byte[jj] = 0;
    }
    info->vram.tilemap[ii] = 0;
    info->vram.refcount[0]++;
  }

  info->vram.index = COLORS;

  for(int ii = 0; ii < NSTATS; ++ii) info->stats[ii] = 0;

}



int compress_tile(struct tile a, FILE *fo)
{
  int seg = 16;
  int s = 0;
  for(int kk = 0; kk < TILE_BYTES; kk += seg) {
    s += 3;
    unsigned char count[256];
    for(int ii = 0; ii < 256; ++ii) count[ii] = 0;
    for(int ii = 0; ii < seg; ++ii) {
      count[a.byte[ii + kk]]++;
    }
    unsigned char byte = 0;
    unsigned char bestcount = 0;
    for(int ii = 0; ii < 256; ++ii) {
      if(count[ii] > bestcount) {
        bestcount = count[ii];
        byte = ii;
      }
    }
    unsigned short flag = 0;
    for(int ii = 0; ii < seg; ++ii) {
      if(a.byte[ii + kk] != byte) {
        flag |= 1 << ii;
        s++;
      }
    }

    fwrite(&flag, 1, 2, fo);
    fwrite(&byte, 1, 1, fo);

    for(int ii = 0; ii < seg; ++ii) {
      if(a.byte[ii + kk] != byte) {
        fwrite(&a.byte[ii + kk], 1, 1, fo);
      }
    }
  }
  return s;
}

enum update_code encode_tile(struct encode_info* info, int ii)
{
  int tm = info->vram.tilemap[ii];
  unsigned char hash = tilehash(info->frame[ii]);

  if(info->vram.hash[tm >> MIRRORBIT][tm & (VRAM_TILES - 1)] == hash)
  if(sametile(info->vram.tiles[tm >> MIRRORBIT][tm & (VRAM_TILES - 1)], info->frame[ii])) {
    return NONE;
  }

  if(!info->vram.refcount[tm & (VRAM_TILES - 1)]) {
    fprintf(stderr, "refcount negative for tile %i, this should not occur\n", tm);
    exit(1);
  }

  info->vram.refcount[tm & (VRAM_TILES - 1)]--;

  for(int kk = 0; kk < 4; ++kk) {
    for(int jj = 0; jj < VRAM_TILES; ++jj) {
      if(info->vram.hash[kk][jj] == hash)
      if(sametile(info->vram.tiles[kk][jj], info->frame[ii])) {
        info->vram.refcount[jj]++;
        info->vram.tilemap[ii] = jj | (kk << MIRRORBIT);
        return CACHE;
      }
    }
  }

  int skip = 0;
  while(info->vram.refcount[(info->vram.index + skip) & (VRAM_TILES - 1)]) {
    skip++;

    if(skip >= 1024) {
      fprintf(stderr, "skipped to many tiles, this should not occur.\n");
      exit(1);
    }
  }

  info->skip[ii] = skip;

  int vi = (info->vram.index + skip) & (VRAM_TILES - 1);
  info->vram.index += skip + 1;

  info->vram.refcount[vi]++;
  info->vram.tilemap[ii] = vi;
  for(int kk = 0; kk < 4; ++kk) {
    struct tile t = mirror(info->frame[ii], kk);
    info->vram.tiles[kk][vi] = t;
    info->vram.hash[kk][vi] = tilehash(t);
  }

  return skip ? UPDATE_SKIP : UPDATE;
}

void encode(struct encode_info* info)
{
  info->vram = info->vram_prev;

  for(int ii = 0; ii < NSTATS; ++ii) info->stats[ii] = 0;

  for(int yy = 0; yy < FRAME; yy += 2) {
    for(int xx = 0; xx < FRAME; xx += 2) {
      for(int yyt = yy; yyt < yy + 2; ++yyt) {
        for(int xxt = xx; xxt < xx + 2; ++xxt) {
          int ii = yyt * FRAME + xxt;
          enum update_code c = encode_tile(info, ii);
          info->update[ii] = c;
          info->stats[c]++;
        }
      }
    }
  }
}

struct priority_info
{
  unsigned short tile;
  unsigned short match;
};

unsigned int count_set_bits(unsigned char n)
{
  unsigned int count = 0;
  while (n) {
    count += n & 1;
    n >>= 1;
  }
  return count;
}

void reduce_update(struct encode_info* info, int compression)
{
  int maxdiff = 155 + compression;

  int reduce_tiles = (info->stats[1] + info->stats[3]) * compression / 100;
  if(reduce_tiles == 0) return;

  struct priority_info priority_bucket[256][TILES_IN_FRAME];
  int priority_count[256];

  for(int ii = 0; ii < 256; ++ii) priority_count[ii] = 0;

  for(int ii = 0; ii < TILES_IN_FRAME; ++ii) {
    if(info->update[ii] == UPDATE || info->update[ii] == UPDATE_SKIP) {
      int match = info->vram_prev.tilemap[ii];
      int jj = match & (VRAM_TILES - 1);
      int kk = match >> MIRRORBIT;
      int best_diff = tilediff(info->frame[ii], info->vram_prev.tiles[kk][jj], maxdiff + 4) - 4;
      if(best_diff < 0) best_diff = 0;
      for(int kk = 0; kk < 4; ++kk) {
        for(int jj = 0; jj < VRAM_TILES; ++jj) {
          int diff = tilediff(info->frame[ii], info->vram_prev.tiles[kk][jj], maxdiff);
          if(diff < best_diff) {
            match = jj | (kk << MIRRORBIT);
            best_diff = diff;
          }
        }
        if(best_diff == 0) break;
      }
      struct priority_info p = { ii, match };
      priority_bucket[best_diff][priority_count[best_diff]] = p;
      priority_count[best_diff]++;
    }
  }

  for(int ii = 0; ii < maxdiff; ++ii) {
    if(!priority_count[ii]) continue;
    int n = reduce_tiles > priority_count[ii] ? priority_count[ii] : reduce_tiles;
    for(int jj = 0; jj < n; ++jj) {
      struct priority_info p = priority_bucket[ii][jj];
      info->frame[p.tile] = info->vram_prev.tiles[p.match >> MIRRORBIT][p.match & (VRAM_TILES - 1)];
      info->flag[p.tile] = 1;
    }
    reduce_tiles -= n;
    if(reduce_tiles == 0) {
      break;
    }
  }
}

void palette(ColorMapObject* cm, int tc, unsigned char *pal)
{
  for(int i = 0; i <= COLORS; i ++) {
    /* unsigned short p = (cm->Colors[i].Red >> 3) */
    /*   | ((cm->Colors[i].Green >> 3) << 5) */
    /*   | ((cm->Colors[i].Blue >> 3) << 10); */

    int k = i;
    if(k > tc) k--;

    pal[k * 3 + 0] = cm->Colors[i].Red;
    pal[k * 3 + 1] = cm->Colors[i].Green;
    pal[k * 3 + 2] = cm->Colors[i].Blue;
  }
}

void getframe(GifFileType *gif, SavedImage *img, unsigned char *canvas)
{
  int tc = gif->SBackGroundColor;
  int delay = 0;

  for(int ii = 0; ii < img->ExtensionBlockCount; ++ii) {
    if(img->ExtensionBlocks[ii].Function == GRAPHICS_EXT_FUNC_CODE) {
      delay = *(unsigned short*)(img->ExtensionBlocks[ii].Bytes + 1);
      break;
    }
  }

  /* ColorMapObject* cm = img->ImageDesc.ColorMap; */
  /* if(cm) palette(cm, tc); */

  for(int yy = 0; yy < img->ImageDesc.Height; ++yy) {
    for(int xx = 0; xx < img->ImageDesc.Width; ++xx) {
      int c = img->RasterBits[yy * img->ImageDesc.Width + xx];
      if(c != tc) {
        if(c > tc) c--;
        if(c >= COLORS - 1) c = COLORS - 1;
        canvas[(xx + img->ImageDesc.Left) + (yy + img->ImageDesc.Top) * gif->SWidth] = c;
      }
    }
  }
}

void bitplane(unsigned char* source, struct tile* target, int width, int height)
{
  int t = 0;

  for(int ln = 0; ln < height / TILE_SIZE; ++ln)
    for(int i = 0; i < width / TILE_SIZE; ++i)
      for(int m = 0; m < BITPLANES; m += 2)
        for(int l = 0; l < TILE_SIZE; ++l)
          for(int k = 0; k < (BITPLANES < 2 ? BITPLANES : 2); ++k) {
            unsigned int byte = 0;

            for(int j = 0; j < TILE_SIZE; ++j) {
              int index = i * TILE_SIZE + l * TILE_SIZE * (width / TILE_SIZE) + j + ln * width * TILE_SIZE;
              if(source[index] & (1 << (k + m)))
                byte |= 1 << (TILE_SIZE - j - 1);
            }

            target[t / TILE_BYTES].byte[t & (TILE_BYTES - 1)] = byte;
            t++;
          }

}

void debitplane(unsigned char* source, unsigned char* target, int width, int height)
{
  for(int yy = 0; yy < height; ++yy)
    for(int xx = 0; xx < width; ++xx) {
      int tile = (yy / TILE_SIZE) * (width / TILE_SIZE) + (xx / TILE_SIZE);
      int bit = 0x80 >> (xx & (TILE_SIZE - 1));
      int row = yy & (TILE_SIZE - 1);
      unsigned char c = 0;
      if(bit & source[tile * FRAME + row * 2 + 0]) c |= 1;
      if(bit & source[tile * FRAME + row * 2 + 1]) c |= 2;
      if(bit & source[tile * FRAME + row * 2 + 16]) c |= 4;
      if(bit & source[tile * FRAME + row * 2 + 17]) c |= 8;
      target[yy * width + xx] = c;
    }
}

int compress(const char* inputfile, const char* outputfile, int nframes, int compression, int verbose)
{
  GifFileType* gif = DGifOpenFileName(inputfile, NULL);

  if(GIF_OK != DGifSlurp(gif)) {
    fprintf(stderr, "error in DGifSlurp()\n");
    return 1;
  }

  if (verbose) fprintf(stderr, "background color = %i\n", gif->SBackGroundColor);

  FILE *fo = fopen(outputfile, "wb");

  if(!fo) {
    fprintf(stderr, "error in fopen()\n");
    return 1;
  }

  struct tie_header header;

  header.width = gif->SWidth;
  header.height = gif->SHeight;
  header.frames = nframes ? nframes : gif->ImageCount;

  if(gif->SColorMap) palette(gif->SColorMap, gif->SBackGroundColor, header.palette);

  for(int kk = 0; kk < gif->SavedImages[0].ExtensionBlockCount; ++kk) {
    if(gif->SavedImages[0].ExtensionBlocks[kk].Function == GRAPHICS_EXT_FUNC_CODE) {
      header.delay = *(unsigned short*)&gif->SavedImages[0].ExtensionBlocks[kk].Bytes[1];
      break;
    }
  }

  fwrite(&header, 1, sizeof(header), fo);

  unsigned char *canvas = malloc(gif->SWidth * gif->SHeight);
  memset(canvas, 0, gif->SWidth * gif->SHeight);

  struct encode_info info;

  init(&info);

  for(int ii = 0; ii < header.frames; ++ii) {

    info.vram_prev = info.vram;
    getframe(gif, &gif->SavedImages[ii], canvas);
    bitplane(canvas, info.frame, gif->SWidth, gif->SHeight);
    encode(&info);
    reduce_update(&info, compression);
    encode(&info);

    unsigned short cellflag[CELLS];

    for(int jj = 0; jj < CELLS; jj++) cellflag[jj] = 0xffff;

    for(int yy = 0; yy < FRAME; yy += 2) {
      for(int xx = 0; xx < FRAME; xx += 2) {
        if(info.update[yy * FRAME + xx] >= 1) continue;
        if(info.update[yy * FRAME + xx + 1] >= 1) continue;
        if(info.update[yy * FRAME + xx + FRAME] >= 1) continue;
        if(info.update[yy * FRAME + xx + FRAME + 1] >= 1) continue;
        cellflag[yy / 2] &= ~(1 << (xx / 2));
        info.stats[5]++;
      }
    }

    fwrite(cellflag, 1, sizeof(cellflag), fo);
    info.stats[4] += CELLS * CELLS / 8;

    for(int yy = 0; yy < FRAME; yy += 2) {
      for(int xx = 0; xx < FRAME; xx += 2) {
        if(cellflag[yy / 2] & (1 << (xx / 2))) {
          unsigned char byte = info.update[yy * FRAME + xx]
            | (info.update[yy * FRAME + xx + 1] << 2)
            | (info.update[yy * FRAME + xx + FRAME] << 4)
            | (info.update[yy * FRAME + xx + FRAME + 1] << 6);
          fwrite(&byte, 1, 1, fo);

          for(int yyt = yy; yyt < yy + 2; ++yyt) {
            for(int xxt = xx; xxt < xx + 2; ++xxt) {
              switch(info.update[yyt * FRAME + xxt]) {
              case CACHE:
                fwrite(&info.vram.tilemap[yyt * FRAME + xxt], 1, 2, fo);
                info.stats[4] += 2;
                break;
              case UPDATE_SKIP:
                fwrite(&info.skip[yyt * FRAME + xxt], 1, 2, fo);
                info.stats[4] += 2;
                break;
              }
            }
          }
        }
      }
    }

    for(int yy = 0; yy < FRAME; yy += 2) {
      for(int xx = 0; xx < FRAME; xx += 2) {
        for(int yyt = yy; yyt < yy + 2; ++yyt) {
          for(int xxt = xx; xxt < xx + 2; ++xxt) {
            int kk = yyt * FRAME + xxt;
            switch(info.update[kk]) {
            case UPDATE_SKIP:
              info.stats[4]++;
            case UPDATE:
              info.stats[4] += compress_tile(info.frame[kk], fo);
              break;
            }
          }
        }
      }
    }

    printf("frame %i: %i unchanged, %i update, %i cache, %i update and skip, %i compressed, %i empty cells\n",
           ii, info.stats[0], info.stats[1], info.stats[2], info.stats[3], info.stats[4], info.stats[5]);

  }

  fclose(fo);

  free(canvas);

  DGifCloseFile(gif, NULL);
}

int decompress(const char* inputfile, const char* outputfile, int nframes, int verbose)
{
  FILE *fi = fopen(inputfile, "rb");

  struct tie_header header;

  fread(&header, 1, sizeof(header), fi);

  if(verbose) {
    printf("%ix%i, %i frames, delay = %i\n", header.width, header.height, header.frames, header.delay);
  }

  if(!nframes) nframes = header.frames;

  GifColorType pal[256];

  for(int ii = 0; ii < COLORS; ++ii) {
    pal[ii].Red = header.palette[ii * 3 + 0];
    pal[ii].Green = header.palette[ii * 3 + 1];
    pal[ii].Blue = header.palette[ii * 3 + 2];
  }

  int error;
  GifFileType* gif = EGifOpenFileName(outputfile, false, &error);

  gif->SWidth = header.width;
  gif->SHeight = header.height;
  gif->SColorResolution = 256;
  gif->SBackGroundColor = 255;
  gif->SColorMap = GifMakeMapObject(256, pal);

  SavedImage image;

  unsigned char *canvas = malloc(gif->SWidth * gif->SHeight);
  memset(canvas, 0, gif->SWidth * gif->SHeight);

  image.ImageDesc.Left = 0;
  image.ImageDesc.Top = 0;
  image.ImageDesc.Width = gif->SWidth;
  image.ImageDesc.Height = gif->SHeight;
  image.ImageDesc.Interlace = false;
  image.ImageDesc.ColorMap = NULL;
  image.RasterBits = canvas;
  image.ExtensionBlockCount = 0;
  image.ExtensionBlocks = NULL;

  struct encode_info info;

  init(&info);

  for(int ii = 0; ii < nframes; ++ii) {
    unsigned short cellflag[CELLS];

    fread(cellflag, 1, sizeof(cellflag), fi);

    unsigned short vram_update_start[256];
    unsigned short vram_update_length[256];
    unsigned short vram_update_count = 0;

    vram_update_start[0] = info.vram.index & (VRAM_TILES - 1);
    vram_update_length[0] = 0;

    for(int yy = 0; yy < FRAME; yy += 2) {
      for(int xx = 0; xx < FRAME; xx += 2) {
        if(cellflag[yy / 2] & (1 << (xx / 2))) {
          unsigned char cell;
          fread(&cell, 1, 1, fi);
          unsigned char update;
          unsigned short skip;
          for(int yyt = yy; yyt < yy + 2; ++yyt) {
            for(int xxt = xx; xxt < xx + 2; ++xxt) {
              update = cell & 3;
              switch(update) {
              case CACHE:
                fread(&info.vram.tilemap[yyt * FRAME + xxt], 1, 2, fi);
                break;
              case UPDATE:
                vram_update_length[vram_update_count]++;
                info.vram.tilemap[yyt * FRAME + xxt] = info.vram.index & (VRAM_TILES - 1);
                info.vram.index++;
                break;
              case UPDATE_SKIP:
                fread(&skip, 1, 2, fi);
                vram_update_count++;
                info.vram.index += skip;
                vram_update_start[vram_update_count] = info.vram.index & (VRAM_TILES - 1);
                vram_update_length[vram_update_count] = 1;
                info.vram.tilemap[yyt * FRAME + xxt] = info.vram.index & (VRAM_TILES - 1);
                info.vram.index++;
                break;
              }
              cell >>= 2;
            }
          }
        }
      }
    }

    if(verbose) printf("frame %i - %i skips\n", ii, vram_update_count);

    for(int jj = 0; jj <= vram_update_count; ++jj) {
      if(verbose) printf("update %i - %i tiles\n", jj, vram_update_length[jj]);

      for(int kk = 0; kk < vram_update_length[jj]; ++kk) {
        struct tile t;
        for(int seg = 0; seg < 32; seg += 16) {
          unsigned short flags;
          unsigned char byte;
          fread(&flags, 1, 2, fi);
          fread(&byte, 1, 1, fi);
          for(int ll = 0; ll < 16; ++ll) {
            if(flags & (1 << ll)) {
              fread(&t.byte[ll + seg], 1, 1, fi);
            } else {
              t.byte[ll + seg] = byte;
            }
          }
        }
        info.vram.tiles[0][(vram_update_start[jj] + kk) & (VRAM_TILES - 1)] = t;
      }
    }

    for(int jj = 0; jj < TILES_IN_FRAME; ++jj) {
      unsigned short tm = info.vram.tilemap[jj];
      info.frame[jj] = mirror(info.vram.tiles[0][tm & (VRAM_TILES - 1)], tm >> MIRRORBIT);
    }

    debitplane((unsigned char*)info.frame, canvas, gif->SWidth, gif->SHeight);

    SavedImage* image2 = GifMakeSavedImage(gif, &image);

      unsigned char graphicscontrol[4] = {0, 0, 0, 0};
      *(unsigned short*)&graphicscontrol[1] = header.delay;

        GifAddExtensionBlock(&image2->ExtensionBlockCount,
                             &image2->ExtensionBlocks,
                             GRAPHICS_EXT_FUNC_CODE,
                             4,
                             graphicscontrol);

  }

  EGifSpew(gif);

  free(canvas);
}


int main(int argc, char *argv[])
{
  int opt;
  int verbose = 0;
  int d = 0;
  int frames = 0;
  int compression = 50;

  char inputfile[1024], outputfile[1024];
  inputfile[0] = 0; outputfile[0] = 0;

  while ((opt = getopt(argc, argv, "i:o:c:f:vd")) != -1) {
    switch (opt) {
    case 'c':
      compression = atoi(optarg);
      break;
    case 'f':
      frames = atoi(optarg);
      break;
    case 'v':
      verbose = 1;
      break;
    case 'd':
      d = 1;
      break;
    case 'i':
      strcpy(inputfile, optarg);
      break;
    case 'o':
      strcpy(outputfile, optarg);
      break;
    default:
      fprintf(stderr, "usage: %s [-v] [-d] [-f frames] [-i input filename] [-o output filename]\n",
              argv[0]);
      return 1;
    }
  }

  if(!inputfile[0]) {
    fprintf(stderr, "must specify input file\n");
    return 1;
  }

  if(!outputfile[0]) {
    fprintf(stderr, "must specify output file\n");
    return 1;
  }

  if(!d) {
    return compress(inputfile, outputfile, frames, compression, verbose);
  } else {
    return decompress(inputfile, outputfile, frames, verbose);
  }
}
