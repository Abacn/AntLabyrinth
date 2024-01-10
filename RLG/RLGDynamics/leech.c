/*
 * Leech lattice quantization.
 * Author: @Abacn.
 *
 * Refer to quantizer_L24 function for usage.
 *
 * It referenced the Leech decoder implementation by
 * Alex van Poppelen (avanpoppelen@gmail.com), based on
 * Vardy and Be'ery, "Maximum Likelihood Decoding of the
 * Leech Lattice", Trans. on Information Theory, Vol.39,
 * No. 4, July 1993.
 */

#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

/* bit operation builtins used */
#ifdef _MSC_VER
#  include <intrin.h>
#  define __builtin_popcount(A) _mm_popcnt_u32(A)
#  define __builtin_popcountll(A) _mm_popcnt_u64(A)
#  define __builtin_parity(A) (__builtin_popcount(A) & 1)
#  define __builtin_parityll(A) (__builtin_popcountll(A) & 1)
#endif

#ifndef __POPCNT__
#pragma warning( POPCNT instruction not enabled)
#endif

/* The representations of the 16 subsets of points
 * in D_2, in a 16-QAM constellation.
 */
const uint8_t D2_SUBSETS[2][8][2] = {
  {{0, 0}, {4, 0}, {2, 0}, {6, 0}, {4, 6}, {4, 2}, {2, 2}, {2, 6}},  // A_000 to A_111
  {{1, 1}, {5, 1}, {3, 1}, {3, 5}, {5, 7}, {1, 7}, {3, 3}, {3, 7}}}; // B_000 to B_111


/* Enumeration of the hexacode over GF(4). The four
 * characters of GF(4) are represented as {0,1,2,3}.
 */
const uint8_t HEXACODES[64][6] = {
  {0, 0, 0, 0, 0, 0}, {0, 0, 1, 1, 1, 1}, {0, 0, 2, 2, 2, 2}, {0, 0, 3, 3, 3, 3},
  {1, 1, 0, 0, 1, 1}, {2, 2, 0, 0, 2, 2}, {3, 3, 0, 0, 3, 3}, {1, 1, 1, 1, 0, 0},
  {2, 2, 2, 2, 0, 0}, {3, 3, 3, 3, 0, 0}, {1, 1, 2, 2, 3, 3}, {1, 1, 3, 3, 2, 2},
  {2, 2, 1, 1, 3, 3}, {2, 2, 3, 3, 1, 1}, {3, 3, 1, 1, 2, 2}, {3, 3, 2, 2, 1, 1},
  {2, 3, 2, 3, 2, 3}, {3, 1, 3, 1, 3, 1}, {1, 2, 1, 2, 1, 2}, {2, 3, 3, 2, 3, 2},
  {3, 2, 2, 3, 3, 2}, {3, 2, 3, 2, 2, 3}, {3, 1, 1, 3, 1, 3}, {1, 3, 3, 1, 1, 3},
  {1, 3, 1, 3, 3, 1}, {1, 2, 2, 1, 2, 1}, {2, 1, 1, 2, 2, 1}, {2, 1, 2, 1, 1, 2},
  {0, 1, 0, 1, 2, 3}, {0, 1, 1, 0, 3, 2}, {1, 0, 0, 1, 3, 2}, {1, 0, 1, 0, 2, 3},
  {0, 1, 2, 3, 0, 1}, {0, 1, 3, 2, 1, 0}, {1, 0, 2, 3, 1, 0}, {1, 0, 3, 2, 0, 1},
  {2, 3, 0, 1, 0, 1}, {2, 3, 1, 0, 1, 0}, {3, 2, 0, 1, 1, 0}, {3, 2, 1, 0, 0, 1},
  {0, 2, 0, 2, 3, 1}, {0, 2, 2, 0, 1, 3}, {2, 0, 0, 2, 1, 3}, {2, 0, 2, 0, 3, 1},
  {0, 2, 3, 1, 0, 2}, {0, 2, 1, 3, 2, 0}, {2, 0, 3, 1, 2, 0}, {2, 0, 1, 3, 0, 2},
  {3, 1, 0, 2, 0, 2}, {3, 1, 2, 0, 2, 0}, {1, 3, 0, 2, 2, 0}, {1, 3, 2, 0, 0, 2},
  {0, 3, 0, 3, 1, 2}, {0, 3, 3, 0, 2, 1}, {3, 0, 0, 3, 2, 1}, {3, 0, 3, 0, 1, 2},
  {0, 3, 1, 2, 0, 3}, {0, 3, 2, 1, 3, 0}, {3, 0, 1, 2, 3, 0}, {3, 0, 2, 1, 0, 3},
  {1, 2, 0, 3, 0, 3}, {1, 2, 3, 0, 3, 0}, {2, 1, 0, 3, 3, 0}, {2, 1, 3, 0, 0, 3}};

typedef struct _Pen {
  double val;
  uint8_t ijk;
  uint8_t l;
  uint8_t c;
} Pen;

typedef struct _Q24 {
  double d_ij [12][4];
  double delta_ij [12][4];
  double mu_x[6][4];
  double S_j[3][4][4];

  uint8_t offsets[12];
  uint8_t ijk1[6][4];
  uint8_t ijk2[6][4];
  Pen pens[3 * 24];
  Pen *pen[3][24];

  uint64_t cv;
  double d;

  uint8_t coset_h;
  uint8_t coset_q;
} Q24;

/* Returns the minimum of two numbers.
 */
static double min2(double u1, double u2)
{
  return u1 < u2 ? u1 : u2;
}

/* Same as min2() but sets a flag.
 */
static double min2_flag(double u1, double u2, uint8_t *flag)
{
  if (u1 < u2) {
    *flag = 0;
    return u1;
  } else {
    *flag = 1;
    return u2;
  }
}

/* Computes the squared Euclidean distance (SED)
 * between two points in ZZ^2_Q.
 */
static double sqrdist(double x1, double y1, double x2, double y2)
{
  double dx = fabs((double)x1 - (double)x2);
  double dy = fabs((double)y1 - (double)y2);
  dx = min2(dx, 8 - dx);
  dy = min2(dy, 8 - dy);
  return dx * dx + dy * dy;
}

/* Decodes a X_ijk subset, for the 32-QAM constellation
 * (see comments above). This function will need to be
 * modified if another constellation is used.
 */
static void decode_subset(const double *t, uint8_t coset_h, uint8_t ijk, double *d, uint8_t *o)
{
  double d1 = sqrdist(t[0], t[1], D2_SUBSETS[coset_h][ijk][0], D2_SUBSETS[coset_h][ijk][1]);
  double d2 = sqrdist(t[0], t[1], (D2_SUBSETS[coset_h][ijk][0] + 4) % 8, (D2_SUBSETS[coset_h][ijk][1] + 4) % 8);

  *d = min2(d1, d2);
  *o ^= (d1 < d2 ? 0 : 1) << ijk;
}

/* The precomputation step to calculate the d_ij's and
 * delta_ij's. This step is done for each Leech half-
 * lattice.
 */
static void precomputation_H24(const double *t, Q24* q, uint8_t coset_h)
{
  int n;
  uint8_t ij;
  double d_ij0, d_ij1;

  // d_ij and delta_ij are set explicitly, but the offsets need to be zeroed.
  memset(q->offsets, 0, sizeof(uint8_t) * 12);
  q->coset_h = coset_h;
  for (n = 0; n < 12; ++n) {
    for (ij = 0; ij < 4; ++ij) {
      decode_subset(t + 2 * n, coset_h, ij<<1, &d_ij0, q->offsets + n);
      decode_subset(t + 2 * n, coset_h, ij<<1 | 1, &d_ij1, q->offsets + n);
      q->d_ij[n][ij] = min2(d_ij0, d_ij1);
      q->delta_ij[n][ij] = d_ij1 - d_ij0;
    }
  }
}

/* Initialize Q24 struct.
 */
static Q24 *init_Q24(Q24* q, uint8_t coset_q)
{
  q->coset_q = coset_q;

  int i, l, c;
  for (i = 0; i < 3; ++i) {
    for (l = 0; l < 6; ++l) {
      for (c = 0; c < 4; ++c) {
        q->pen[i][4 * l + c] = &q->pens[24 * i + 4 * l + c];
        q->pen[i][4 * l + c]->l = l;
        q->pen[i][4 * l + c]->c = c;
      }
    }
  }
  return q;
}

/* Computes confidence values, preferable representations,
 * and penalties. c refers to the character of GF(4).
 *
 * The preferable representations are calculated from the
 * (i1,j1,i2,j2) interpretation of c, and its complement.
 */
static void compute_vals(Q24 *q, int l, int c, uint8_t ij1, uint8_t ij2)
{
  double d1 = q->d_ij[2 * l][ij1] + q->d_ij[2 * l + 1][ij2];
  double d2 = q->d_ij[2 * l][ij1 ^ 0x3] + q->d_ij[2 * l + 1][ij2 ^ 0x3];

  // preferable representations
  if (d2 < d1) {
    ij1 ^= 0x3; ij2 ^= 0x3;
  }
  q->ijk1[l][c] = ij1<<1 | (q->delta_ij[2 * l][ij1] < 0);
  q->ijk2[l][c] = ij2<<1 | (q->delta_ij[2 * l + 1][ij2] < 0);

  // confidence value
  q->mu_x[l][c] = min2(d1, d2);

  // penalties
  uint8_t bit, ck1, ck2;
  ck1 = (q->delta_ij[2 * l][ij1 ^ 0x3] < 0) ^ (q->ijk1[l][c] & 1);
  ck2 = (q->delta_ij[2 * l + 1][ij2 ^ 0x3] < 0) ^ (q->ijk2[l][c] & 1);

  q->pen[0][4 * l + c]->val = min2_flag(fabs(q->delta_ij[2 * l][ij1]), fabs(q->delta_ij[2 * l + 1][ij2]), &bit);
  q->pen[0][4 * l + c]->ijk = bit ? 010 : 001;

  double val1 = q->d_ij[2 * l][ij1 ^ 0x3] + q->d_ij[2 * l + 1][ij2 ^ 0x3] - q->d_ij[2 * l][ij1] - q->d_ij[2 * l + 1][ij2];
  uint8_t ijk1 = 066 ^ (ck1 ? 001 : 0) ^ (ck2 ? 010 : 0);
  double val2 = val1 + min2_flag(fabs(q->delta_ij[2 * l][ij1 ^ 0x3]), fabs(q->delta_ij[2 * l + 1][ij2 ^ 0x3]), &bit);
  uint8_t ijk2 = 066 ^ (ck1 ? 001 : 0) ^ (ck2 ? 010 : 0) ^ (bit ? 010 : 001);

  if (ck1 == ck2) {
    q->pen[1][4 * l + c]->val = val1;
    q->pen[1][4 * l + c]->ijk = ijk1;
    q->pen[2][4 * l + c]->val = val2;
    q->pen[2][4 * l + c]->ijk = ijk2;
  } else {
    // switch last two penalties based on intrinsic
    // change to k-parity
    q->pen[1][4 * l + c]->val = val2;
    q->pen[1][4 * l + c]->ijk = ijk2;
    q->pen[2][4 * l + c]->val = val1;
    q->pen[2][4 * l + c]->ijk = ijk1;
  }
}

/* Seek to the first penalty matching a hexacode-
 * word digit in a sorted array of penalties.
 */
static Pen **penalty_seek(Pen **pens, const uint8_t *word)
{
  while ((*pens)->c != word[(*pens)->l]) {
    ++pens;
  }
  return pens;
}

static void sort_penalties(Q24 *q)
{
  int l;

  // Step 2: Sorting the penalties
  for (l = 0; l < 3; ++l) {
    Pen **pens = q->pen[l];
    int i, j;
    for (i = 1; i < 24; ++i) {
      Pen *p = pens[i];
      j = i;
      while (j > 0 && p->val < pens[j - 1]->val) {
        pens[j] = pens[j - 1];
        --j;
      }
      pens[j] = p;
    }
  }
}

/* Decodes the Leech quarter-lattice.
 */
static void decoder_Q24(Q24 *q)
{
  int penalty_sorted = 0;

  // Step 1: Computing the confidence values,
  // preferable representations, and penalties
  int l;
  for (l = 0; l < 6; ++l) {
    if (q->coset_q == 0) {
      compute_vals(q, l, 0, 0, 0);
      compute_vals(q, l, 1, 0, 3);
      compute_vals(q, l, 2, 1, 1);
      compute_vals(q, l, 3, 1, 2);
    } else {
      compute_vals(q, l, 0, 2, 0);
      compute_vals(q, l, 1, 1, 0);
      compute_vals(q, l, 2, 0, 2);
      compute_vals(q, l, 3, 0, 1);
    }
  }

  // Step 3: Computing the confidence values
  // of the blocks
  int x1, x2;
  for (l = 0; l < 3; ++l) {
    for (x1 = 0; x1 < 4; ++x1) {
      for (x2 = 0; x2 < 4; ++x2) {
        q->S_j[l][x1][x2] = q->mu_x[2 * l][x1] + q->mu_x[2 * l + 1][x2];
      }
    }
  }

  // Step 4: Finding the images of the hexacodewords
  // and computing their metrics
  int w;
  double mind = q->d;
  uint64_t min_pts = 0;

  for (w = 0; w < 64; ++w) {
    const uint8_t *word = HEXACODES[w];
        // calculate base metric
    double met = q->S_j[0][word[0]][word[1]] + q->S_j[1][word[2]][word[3]] + q->S_j[2][word[4]][word[5]];

    if (met > mind) {
      continue;
    }
    // calculate parities
    uint8_t h_parity = 0;
    uint8_t k_parity = 0;
    for (l = 0; l < 6; ++l) {
      h_parity ^= q->ijk1[l][word[l]] >> 2;
      k_parity ^= (q->ijk1[l][word[l]] ^ q->ijk2[l][word[l]]) & 1;
    }
    // check parity relative to Q24 coset
    h_parity ^= q->coset_q;
    k_parity ^= q->coset_h;

    uint64_t pts = 0;
    // find and apply appropriate penalties
    if (h_parity || k_parity) {
      // In the array containing pointers to the penalties,
      // 0 stands for altering k_parity, 1 for h_parity, and
      // 2 for both.
      uint8_t a, b1, b2;
      if (h_parity && !k_parity) {
        a = 1; b1 = 0; b2 = 2;
      } else if (!h_parity && k_parity) {
        a = 0; b1 = 1; b2 = 2;
      } else {
        a = 2; b1 = 0; b2 = 1;
      }

      // sort penalties before first use
      if (!penalty_sorted)
      {
        sort_penalties(q);
        penalty_sorted = 1;
      }

      // check two penalty case
      Pen **pos1 = penalty_seek(q->pen[b1], word);
      Pen **pos2 = penalty_seek(q->pen[b2], word);
      Pen *p1 = *pos1;
      Pen *p2 = *pos2;

      // -> account for coordinate collision
      if (p1->l == p2->l) {
        Pen *p1_next = *penalty_seek(pos1 + 1, word);
        Pen *p2_next = *penalty_seek(pos2 + 1, word);
        if (p1->val + p2_next->val < p2->val + p1_next->val) {
          p2 = p2_next;
        } else {
          p1 = p1_next;
        }
      }

      // compare to single penalty case
      Pen *pa = *penalty_seek(q->pen[a], word);
      if (pa->val < p1->val + p2->val) {
        p1 = pa;
        p2 = NULL;
      }

      // apply penalties
      met += p1->val + (p2 == NULL ? 0 : p2->val);
      pts = ((uint64_t)p1->ijk << (6 * p1->l)) ^ (p2 == NULL ? 0 : ((uint64_t)p2->ijk << (6 * p2->l)));
    }

    // Step 5: Final minimization
    if (met < mind) {
      mind = met;
      // construct base point
      for (l = 0; l < 6; ++l) {
        pts ^= ((uint64_t)q->ijk2[l][word[l]] << 3 | (uint64_t)q->ijk1[l][word[l]]) << (l * 6);
      }
      min_pts = pts;
    }
  }
  if (q->d > mind)
  {
    q->d = mind;
    q->cv = min_pts;
    // attach offset information
    int n;
    for (n = 0; n < 12; ++n) {
      uint8_t ijk = q->cv >> (3 * n) & 0x7;
      min_pts ^= (uint64_t)(q->offsets[n] >> ijk & 1) << (36 + n);
    }
    q->cv = min_pts;
  }
}

/* Find the nearest mod 8 image of y (in [0, 8)) of (x in [0, 8)) */
static int32_t nearest_mod8(double x, int32_t y)
{
  if(x < y)
  {
    if(x-(y-8) < y-x) return y-8;
  }
  else
  {
    if((y+8)-x < x-y) return y+8;
  }
  return y;
}

/**
 * Leech lattice quantization.
 * Given a point t, return its nearest lattice point lpoint.
 * The lattice points are scale with nearest neighbor distance=4 sqrt(2),
 *
 */
void quantizer_L24(const double t[24], int32_t lpoint[24])
{
  double offset[24];
  double in[24];
  int rp;
  Q24 q;
  uint8_t coset_h, coset_q;

  // fit input into [0, 8)^24
  for(rp=0; rp<24; ++rp)
  {
    offset[rp] = floor(t[rp] / 8.0) * 8.0;
    in[rp] = (double)(t[rp] - offset[rp]);
  }

  // Decodes the Leech lattice by means of 4 Leech quarter-lattice decoders.
  q.d = DBL_MAX;
  for (coset_h = 0; coset_h < 2; ++coset_h) {
    precomputation_H24(in, &q, coset_h);

    for (coset_q = 0; coset_q < 2; ++coset_q) {
      init_Q24(&q, coset_q);
      decoder_Q24(&q);
    }
  }

  // get coordinate from q.cv
  const uint64_t mask_bits = 0777777777777ull; // 36 bits
  const uint64_t mask_overallk = 0111111111111ull;
  uint8_t k_parity = __builtin_parityll(q.cv & mask_overallk);
  uint64_t bits = q.cv & mask_bits;
  uint64_t offsetbits = q.cv >> 36;
  uint8_t ijk;
  int32_t x, y;
  for(rp=0; rp<12; ++rp)
  {
    ijk = bits & 07;
    x = D2_SUBSETS[k_parity][ijk][0];
    y = D2_SUBSETS[k_parity][ijk][1];
    if(offsetbits & 1)
    {
      x = (x + 4) % 8;
      y = (y + 4) % 8;
    }
    lpoint[2 * rp] = offset[2 * rp] + nearest_mod8(in[2 * rp], x);
    lpoint[2 * rp + 1] = offset[2 * rp + 1] + nearest_mod8(in[2 * rp + 1], y);
    offsetbits >>= 1;
    bits >>= 3;
  }
}
