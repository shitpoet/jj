//todo: copy EXIF orientation if present
//todo: mb copy colorspace? (Adobe RGB)
//todo: remove of if-checks in decoding by making reading functions tolerant
//todo: mb move counting to decoding or blend counting and encoding code
//todo: maybe merge dc/ac_count_ptrs arrays
//todo: seems like we decode huffman symbols and then convert them to coefs, mb skip that?
//      the problem is that decoded values are refined bit by bit in progressive mode
//      baseline would allow at least DC values without last_dc_val-delting for sure (tested it)
//      but i'm reluctant to drop progressive mode support bc progressive images
//      can easily be non-optimized but produced by some editing software (or even camera app?)
//      nevertheless i would give 700 bytes stright away by dropping progressive decoding
//      and in the end i think about several ks drop would be possible (mb even better...)
//todo: mb read components' info stright into C structure
//todo: 
/* Scan component selector – Selects which of the Nf image components specified in the frame parameters
shall be the jth component in the scan. Each Csj shall match one of the Ci values specified in the frame header,
and the ordering in the scan header shall follow the ordering in the frame header. If Ns > 1, the order of
interleaved components in the MCU is Cs1 first, Cs2 secon */
// https://www.w3.org/Graphics/JPEG/itu-t81.pdf (pg. 41)

// spec fact: size of DCT coefficient is 12 bit (including sign)

//#include <emscripten.h>
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <stdint.h>
#include <stdbool.h>
#include <ctype.h>
#include <limits.h>
#include <string.h>


#define R restrict
#define noinline __attribute__((noinline)) 

#define MAX(a, b) ((a) > (b) ? (a) : (b))
// noinline static int MAX(int a, int b) {
//   return a >= b ? a : b;
// }
#define MIN(a, b) ((a) < (b) ? (a) : (b))

// Define if your (broken) compiler shifts signed values as if they were unsigned.
//#define RIGHT_SHIFT_IS_UNSIGNED 1
#define USE_CLZ_INTRINSIC

#define JPEG_MAX_DIMENSION 16383L 
#define BITS_IN_JSAMPLE 8
#define MAX_COEF_BITS 10
#define NCOMPS 3
#define DCTSIZE 8                               /* The basic DCT block is 8x8 samples */
#define DCTSIZE2 64                             /* DCTSIZE squared; # of elements in a block */
#define NUM_QUANT_TBLS 4                        /* Quantization tables are numbered 0..3 */
#define NUM_HUFF_TBLS 4                         /* Huffman tables are numbered 0..3 */
#define MAX_COMPS_IN_SCAN NCOMPS                /* JPEG limit on # of components in one scan */
#define MAX_SAMP_FACTOR 4                       /* JPEG limit on sampling factors */
#define MAX_BLOCKS_IN_MCU 10 

// #define boolean int
typedef int32_t int32min_t; /* at least signed 32-bit values. */
#define boolean bool
// typedef short coef_t;
typedef int coef_t;
typedef unsigned int freq_t;

#define BIT_BUF_64_BIT 
#ifdef BIT_BUF_64_BIT 
typedef uint64_t bit_buffer_t; /* type of bit-extraction buffer */
#define BIT_BUF_SIZE 64      /* size of buffer in bits */
#define BIT_BUF_FILL ((bit_buffer_t)0xffffffffffffffff)
// #define MIN_GET_BITS (BIT_BUF_SIZE - 8)
#else
typedef uint32_t bit_buffer_t; /* type of bit-extraction buffer */
#define BIT_BUF_SIZE 32 /* size of buffer in bits */
#define BIT_BUF_FILL ((bit_buffer_t))0xffffffff)
#endif

/* Data structures for images (arrays of samples and of DCT coefficients). */
typedef coef_t block_t[DCTSIZE2];   /* one block of coefficients */

typedef enum { 
  M_SOF0 = 0xc0, M_SOF1 = 0xc1, M_SOF2 = 0xc2, M_SOF3 = 0xc3,
  M_SOF5 = 0xc5, M_SOF6 = 0xc6, M_SOF7 = 0xc7, 
  M_JPG = 0xc8, M_SOF9 = 0xc9, M_SOF10 = 0xca, M_SOF11 = 0xcb, 
  M_SOF13 = 0xcd, M_SOF14 = 0xce, M_SOF15 = 0xcf,

  M_DHT = 0xc4, M_DAC = 0xcc,

  M_RST0 = 0xd0, M_RST1 = 0xd1, M_RST2 = 0xd2, M_RST3 = 0xd3, 
  M_RST4 = 0xd4, M_RST5 = 0xd5, M_RST6 = 0xd6, M_RST7 = 0xd7,

  M_SOI = 0xd8, M_EOI = 0xd9, M_SOS = 0xda,
  M_DQT = 0xdb, 
  M_DNL = 0xdc, M_DRI = 0xdd, M_DHP = 0xde, M_EXP = 0xdf,

  M_APP0 = 0xe0, M_APP1 = 0xe1, M_APP2 = 0xe2, M_APP3 = 0xe3,
  M_APP4 = 0xe4, M_APP5 = 0xe5, M_APP6 = 0xe6, M_APP7 = 0xe7,
  M_APP8 = 0xe8, M_APP9 = 0xe9, M_APP10 = 0xea, M_APP11 = 0xeb,
  M_APP12 = 0xec, M_APP13 = 0xed, M_APP14 = 0xee, M_APP15 = 0xef,

  M_JPG0 = 0xf0, M_JPG13 = 0xfd, 
  
  M_COM = 0xfe, 

  M_TEM = 0x01,
} marker_t;

typedef struct {
  // nb: `huffval` should be right after `bits` in memory! 
  uint8_t bits[17];     /* bits[k] = # of symbols with codes of */
                        /* length k bits; bits[0] is unused */
  uint8_t huffval[256]; /* The symbols, in order of incr code length */
  boolean read;
} huff_table_t;

#define AC 0x10 // ac hufftable table index offset

typedef struct {
  int component_id;    /* identifier for this component (0..255) */
  int component_index; /* its index in SOF or comp_info[] */
  int h_samp_factor;   /* horizontal sampling factor (1..4) */
  int v_samp_factor;   /* vertical sampling factor (1..4) */
  int dc_tbl_no;       /* DC entropy table selector (0..3) */
  int ac_tbl_no;       /* AC entropy table selector (0..3) */
  int width_in_blocks;
  int mcu_width;       /* number of blocks per MCU, horizontally */
  int mcu_height;      /* number of blocks per MCU, vertically */
} component_t;

#define HUFF_LOOKAHEAD  8       /* # of bits of lookahead */

/* Derived data constructed for each Huffman table */
typedef struct {
  uint8_t*R huffval; // link to huffman table (needed only in huff_decode)

  /* Basic tables: (element [0] of each array is unused) */
  int32min_t maxcode[18]; /* largest code of length k (-1 if none) */
  /* (maxcode[17] is a sentinel to ensure huff_decode terminates) */
  int32min_t valoffset[18]; /* huffval[] offset for codes of length k */
  /* valoffset[k] = huffval[] index of 1st symbol of code length k, less
   * the smallest code of length k; so given a code of length k, the
   * corresponding symbol is huffval[code + valoffset[k]] */

  /* Lookahead table: indexed by the next HUFF_LOOKAHEAD bits of
   * the input data stream. If the next Huffman code is no more
   * than HUFF_LOOKAHEAD bits long, we can obtain its length and
   * the corresponding symbol directly from this tables.
   *
   * The lower 8 bits of each table entry contain the number of
   * bits in the corresponding Huffman code, or HUFF_LOOKAHEAD + 1
   * if too long. The next 8 bits of each entry contain the
   * symbol. */
  uint16_t lookup[1 << HUFF_LOOKAHEAD];

  // used for compression
  unsigned int ehufco[256]; /* code for each symbol */
  int ehufsi[256];         /* length of code for each symbol */
  /* If no code has been allocated for a symbol S, ehufsi[S] contains 0 */
} derived_table_t; // universal

typedef derived_table_t c_derived_table_t;
typedef derived_table_t d_derived_table_t;

/* If long is > 32 bits on your machine, and shifting/masking longs is
 * reasonably fast, making bit_buffer_t be long and setting BIT_BUF_SIZE
 * appropriately should be a win.  Unfortunately we can't define the size
 * with something like  #define BIT_BUF_SIZE (sizeof(bit_buffer_t)*8)
 * because not all machines measure sizeof in 8-bit bytes. */
typedef struct {
  union {
    uint8_t* restrict next_input_byte; /* => next byte to read from buffer */
    const uint8_t* restrict next_output_byte_2; /* => next byte to read from buffer */
  };
  uint8_t* restrict next_output_byte; /* => next byte to write in buffer */
  union {
    const uint8_t* restrict input_end; /* => next byte to write in buffer */
    uint8_t* restrict output_end; /* => next byte to write in buffer */
  };
  
  bit_buffer_t get_buffer; /* current bit-extraction buffer */
  int bits_left;           /* # of unused bits in it */

  /* ptrs to Huffman coding tables, or NULL if not defined */
  // huff_table_t dc_huff_tbls[NUM_HUFF_TBLS];
  // huff_table_t ac_huff_tbls[NUM_HUFF_TBLS];
  huff_table_t huff_tbls[AC + NUM_HUFF_TBLS];
} state_t;

typedef state_t*R d_state_ptr;
typedef state_t*R c_state_ptr;

/* Left shift macro that handles a negative operand without causing any
 * sanitizer warnings */
#define LEFT_SHIFT(a, b) ((int32min_t)((unsigned long)(a) << (b)))
/* We assume that right shift corresponds to signed division by 2 with
 * rounding towards minus infinity.  This is correct for typical "arithmetic
 * shift" instructions that shift in copies of the sign bit.  But some
 * C compilers implement >> with an unsigned shift.  For these machines you
 * must define RIGHT_SHIFT_IS_UNSIGNED.
 * RIGHT_SHIFT provides a proper signed right shift of a int32min_t quantity.
 * It is only applied with constant shift counts.  SHIFT_TEMPS must be
 * included in the variables of any routine using RIGHT_SHIFT.
 */
#ifdef RIGHT_SHIFT_IS_UNSIGNED
#define SHIFT_TEMPS JLONG shift_temp;
#define RIGHT_SHIFT(x, shft) \
  ((shift_temp = (x)) < 0 ? (shift_temp >> (shft)) | ((~((int32min_t)0)) << (32 - (shft))) : (shift_temp >> (shft)))
#else
#define SHIFT_TEMPS
#define RIGHT_SHIFT(x, shft) ((x) >> (shft))
#endif

/* jpeg_natural_order[i] is the natural-order position of the i'th element
 * of zigzag order.
 * When reading corrupted data, the Huffman decoders could attempt
 * to reference an entry beyond the end of this array (if the decoded
 * zero run length reaches past the end of the block).  To prevent
 * wild stores without adding an inner-loop test, we put some extra
 * "63"s after the real entries.  This will cause the extra coefficient
 * to be stored in location 63 of the block, not somewhere random.
 * The worst case would be a run-length of 15, which means we need 16
 * fake entries. */
static const uint8_t jpeg_natural_order[DCTSIZE2 + 16] = {
  0, 1, 8, 16, 9, 2, 3, 10,
  17, 24, 32, 25, 18, 11, 4, 5,
  12, 19, 26, 33, 40, 48, 41, 34,
  27, 20, 13, 6, 7, 14, 21, 28,
  35, 42, 49, 56, 57, 50, 43, 36,
  29, 22, 15, 23, 30, 37, 44, 51,
  58, 59, 52, 45, 38, 31, 39, 46,
  53, 60, 61, 54, 47, 55, 62, 63,
  63, 63, 63, 63, 63, 63, 63, 63, /* extra entries for safety in decoder */
  63, 63, 63, 63, 63, 63, 63, 63
};

noinline static int read_byte(d_state_ptr state) {
  uint8_t* p = state->next_input_byte;
  if (p < state->input_end) {
    const uint8_t x = *p++;
    state->next_input_byte = p;
    return x;
  } else {
    return -1;
  }
}

int read_short(d_state_ptr dinfo) {
  int hi = read_byte(dinfo);
  if (hi == -1) return -1;
  int lo = read_byte(dinfo);
  if (lo == -1) return -1;
  return (((unsigned int)hi) << 8) + (unsigned int)lo;
}

// note: get_bits should be able to return 0 if nbits == 0
noinline static int peek_bits(d_state_ptr state, int nbits) {
  bit_buffer_t get_buffer = state->get_buffer;
  int bits_left = state->bits_left;
  //todo: it is not more efficient to get bit_buf_size
  //      bc it increases still buffer by one byte every 8 bits
  //note: nbits can be up to 12 for coefficients
  // printf("bits_left %d \n", bits_left);
  if (bits_left < HUFF_LOOKAHEAD || bits_left < nbits) { 
    // printf("fill buffer - ");
    //inv: BIT_BUF_SIZE - 8 >= 12
    while (bits_left <= BIT_BUF_SIZE - 8) { 
    // while (bits_left <= nbits) { 
      int c = read_byte(state);
      if (c < 0) return -1;      
      // if it's 0xFF, check and discard stuffed zero byte 
      if (c == 0xff) {
        // c = input_non_ff_byte(dstate);
        do {
          c = read_byte(state);
          // if (c < 0) return false;
        } while (c == 0xff);
        if (c == 0) {
          c = 0xff; /* Found FF/00, which represents an FF data byte */
        } else { // restart or other marker
          // // unreachable if jpeg file is correct and not corrupted
          // // it is working with 8 bit max prereading
          // return -1;
        
          // some marker reached - rewind buffer pointer
          state->next_input_byte -= 2; // go to the 0xff marker prefix
          // printf("rewind %x %x \n", *state->next_input_byte, state->next_input_byte[1]);
          // printf("nbits %d  bits_left %d \n", nbits, bits_left);
          if (nbits <= bits_left) {
            break; // dont shift bits over the marker
          } else {          
            if (bits_left > 0 && nbits == HUFF_LOOKAHEAD) {
              // expand exisintg bits to fill in HUFF_LOOKAHEAD
              int fill = HUFF_LOOKAHEAD - bits_left;    
              get_buffer = (get_buffer << fill);
              bits_left = HUFF_LOOKAHEAD;
              break;
            } else {
              printf("not enough bits nbits: %d, bits_left: %d \n", nbits, bits_left);
              return -1; // not enough bits
            }
          }
        }
      }
      /* OK, load c into get_buffer */
      get_buffer = (get_buffer << 8) | c;
      bits_left += 8;
    } /* end while */
    // printf("now %d bits_left \n", bits_left);
  }
  state->get_buffer = get_buffer;
  state->bits_left = bits_left;
  return (((int)(get_buffer >> (bits_left - nbits))) & ((1 << nbits) - 1)); 
}

noinline static int get_bits(d_state_ptr state, int nbits) {
  int val = peek_bits(state, nbits);
  state->bits_left -= nbits;
  return val;
}

noinline static int get_bit(d_state_ptr state) {
  return get_bits(state, 1);
}

#define NEG_1 ((unsigned int)-1)
// #define extend(x, s) \
//   ((x) + ((((x) - (1 << ((s)-1))) >> 31) & (((NEG_1) << (s)) + 1)))

// note: extend should return 0 for x == 0 && s == 0
inline static int extend(int x, int s) {
  //return ((x) + ((((x) - (1 << ((s)-1))) >> 31) & (((NEG_1) << (s)) + 1)));
  return x + (
      ((x - (1 << (s - 1))) >> 31) & 
      ((NEG_1 << s) + 1)
  );
}  
 
noinline static int get_bits_extend(d_state_ptr state, int s) {
  return extend(get_bits(state, s), s);
}

// #define GET_BITS_ACTION(result, state, nbits, action) { \
//   result = get_bits(state, nbits); \
//   if (result < 0) { \
//     action; \
//   } \
// }

noinline static void fill(void*R mem, const int size, const uint8_t val) {
  uint8_t*R ptr = mem;
  for (int i = 0; i < size; i++) {
    *ptr++ = val;
  }
}

noinline static void zero(void*R mem, const int size) {
  fill(mem, size, 0);
}

#define NUM_MEMORY_BLOCKS ( \
  1 + /* src_buf */ \
  1 + /* pool */ \
  NCOMPS + /* whole_image subimages */ \
  NUM_HUFF_TBLS * 2 + /* dc/ac_derived_tbls */ \
  1 + /* state */ \
  NCOMPS + /* whole_image subimages */ \
  0 \
)

void* memory_blocks[NUM_MEMORY_BLOCKS];
int num_memory_blocks = 0;

__attribute__((malloc))
noinline static void*R alloc_mem(int sizeofobject) {
  void* ptr = calloc(sizeofobject, 1);
  memory_blocks[num_memory_blocks++] = ptr;
  return ptr;
}

static void free_mem() {
  for (int i = 0; i < NUM_MEMORY_BLOCKS; i++) {
    free(memory_blocks[i]);
  }
  num_memory_blocks = 0;
}

__attribute__((malloc))
inline static block_t* alloc_block_array(int blocksperrow, int numrows) {
  return alloc_mem(
    numrows * blocksperrow * (int)sizeof(block_t)
  );
}

/* NOTE: If USE_CLZ_INTRINSIC is defined, then clz/bsr instructions will be
 * used for bit counting rather than the lookup table.  This will reduce the
 * memory footprint by 64k, which is important for some mobile applications
 * that create many isolated instances of libjpeg-turbo (web browsers, for
 * instance.)  This may improve performance on some mobile platforms as well.
 * This feature is enabled by default only on Arm processors, because some x86
 * chips have a slow implementation of bsr, and the use of clz/bsr cannot be
 * shown to have a significant performance impact even on the x86 chips that
 * have a fast implementation of it.  When building for Armv6, you can
 * explicitly disable the use of clz/bsr by adding -mthumb to the compiler
 * locals (this defines __thumb__). */

/* NOTE: Both GCC and Clang define __GNUC__ */
#if (defined(__GNUC__) && (defined(__arm__) || defined(__aarch64__))) || \
    defined(_M_ARM) || defined(_M_ARM64)
#if !defined(__thumb__) || defined(__thumb2__)
#define USE_CLZ_INTRINSIC
#endif
#endif

// #ifdef HAVE_INTRIN_H
// #include <intrin.h>
// #ifdef _MSC_VER
// #ifdef HAVE_BITSCANFORWARD64
// #pragma intrinsic(_BitScanForward64)
// #endif
// #ifdef HAVE_BITSCANFORWARD
// #pragma intrinsic(_BitScanForward)
// #endif
// #endif
// #endif

#ifdef USE_CLZ_INTRINSIC
#if defined(_MSC_VER) && !defined(__clang__)
#define JPEG_NBITS_NONZERO(x) (32 - _CountLeadingZeros(x))
#else
#define JPEG_NBITS_NONZERO(x) (32 - __builtin_clz(x))
#endif
#define JPEG_NBITS(x) (x ? JPEG_NBITS_NONZERO(x) : 0)
#else
#define JPEG_NBITS(x) (jpeg_nbits_table[x])
#define JPEG_NBITS_NONZERO(x) JPEG_NBITS(x)
#endif

noinline static boolean make_derived_table(
  d_state_ptr state,
  boolean is_dc, int tblno,
  d_derived_table_t* restrict* restrict pdtbl
) {
  if (tblno < 0 || tblno >= NUM_HUFF_TBLS) return false; // JERR_NO_HUFF_TABLE
  huff_table_t*R htbl = is_dc ? 
    &state->huff_tbls[tblno] : 
    &state->huff_tbls[tblno + AC]; 
  if (!htbl->read) return 0; // JERR_NO_HUFF_TABLE
  
  unsigned int huffcode[257];
  char huffsize[257];
  
  /* Figure C.1: make table of Huffman code length for each symbol */
  int p = 0;
  //note: bits[0] should be skipped bc of another simplification elsewhere
  for (int l = 1; l <= 16; l++) {
    int i = (int)htbl->bits[l];
    /* protect against table overrun */
    if (i < 0 | i + p > 256) return 0; // JERR_BAD_HUFF_TABLE
    while (i--) {
       huffsize[p++] = (char)l;
    }
  }
  huffsize[p] = 0; // required
  
  const int max_symbols = p;

  /* Figure C.2: generate the codes themselves */
  /* We also validate that the counts represent a legal Huffman code tree. */

  unsigned int code = 0;
  int si = huffsize[0];
  p = 0;
  while (huffsize[p]) {
    while (((int)huffsize[p]) == si) {
      huffcode[p++] = code;
      code++;
    }
    /* code is now 1 more than the last code used for codelength si; but
     * it must still fit in si bits, since no code is allowed to be all ones. */
    if (((int32min_t)code) >= (((int32min_t)1) << si)) return 0; // ERREXIT(cstate, JERR_BAD_HUFF_TABLE);
    code <<= 1;
    si++;
  }
  
  /* Allocate a workspace if we haven't already done so. */
  if (*pdtbl == NULL) {
    *pdtbl = alloc_mem(sizeof(c_derived_table_t));
    if (*pdtbl == NULL) return false;
  }
  c_derived_table_t*R dtbl = *pdtbl;  
  zero(dtbl, sizeof(c_derived_table_t));
  
  //////////////////////////////////////
  
  {
    int p = 0;
    // for (l = 1; l <= 16; l++) {
    // ors: start from zero to have maxcode[0] == -1 to simpify huff_decode a bit
    for (int l = 0; l <= 16; l++) {
      int max_code;
      if (htbl->bits[l]) {
        /* valoffset[l] = huffval[] index of 1st symbol of code length l,
        * minus the minimum code of length l */
        dtbl->valoffset[l] = (int32min_t)p - (int32min_t)huffcode[p];
        p += htbl->bits[l];
        max_code = huffcode[p - 1]; /* maximum code of length l */
      } else {
        max_code = -1;
      }
      dtbl->maxcode[l] = max_code; // -1 if no codes of this length 
    }
    dtbl->valoffset[17] = 0;
    dtbl->maxcode[17] = 0xFFFFF; /* ensures huff_decode terminates */

    dtbl->huffval = htbl->huffval; /* fill in back link */  
  }
  
  for (int i = 0; i < (1 << HUFF_LOOKAHEAD); i++)
    dtbl->lookup[i] = (HUFF_LOOKAHEAD + 1) << HUFF_LOOKAHEAD;

  p = 0;
  for (int l = 1; l <= HUFF_LOOKAHEAD; l++) {
    for (int i = 1; i <= (int)htbl->bits[l]; i++, p++) {
      /* l = current code's length, p = its index in huffcode[] & huffval[]. */
      /* Generate left-justified code followed by all possible bit sequences */
      int lookbits = huffcode[p] << (HUFF_LOOKAHEAD - l);
      for (int ctr = 1 << (HUFF_LOOKAHEAD - l); ctr > 0; ctr--) {
        dtbl->lookup[lookbits] = (l << HUFF_LOOKAHEAD) | htbl->huffval[p];
        lookbits++;
      }
    }
  }  
  
  //////////////////////////////////////
  /* Figure C.3: generate encoding tables */

  /* These are code and size indexed by symbol value */

  /* Set all codeless symbols to have code length 0;
   * this lets us detect duplicate VAL entries here, and later
   * allows emit_bits to detect any attempt to emit such symbols. */
  
  for (int p = 0; p < max_symbols; p++) {
    int v = htbl->huffval[p];
    if (dtbl->ehufsi[v] || is_dc && v > 15) return false; // ERREXIT(cstate, JERR_BAD_HUFF_TABLE);
    dtbl->ehufco[v] = huffcode[p];
    dtbl->ehufsi[v] = huffsize[p];
  }

  return true;
}

noinline static bool write_byte(c_state_ptr state, int val) {
  uint8_t* restrict p = state->next_output_byte;
  *p++ = (uint8_t)val;
  if (p >= state->output_end) return false;
  state->next_output_byte = p;
  return true;
}

/* Emit a 2-byte integer; these are always MSB first in JPEG files */
static void write_short(c_state_ptr state, int value) {
  write_byte(state, (value >> 8) & 0xFF);
  write_byte(state, value & 0xFF);
}

static void write_marker(c_state_ptr state, marker_t mark) {
  write_byte(state, 0xFF);
  write_byte(state, (int)mark);
}

// doesnt help
// static void write_marker_len(c_state_ptr cstate, marker_t code, int length) {
//   write_marker(cstate, code);
//   write_short(cstate, length);
// }

inline
static void write_block_unsafe(c_state_ptr state, const void *ptr, int size) {
  memmove(state->next_output_byte, ptr, size);
  state->next_output_byte += size;
  // cstate->free_in_buffer -= size;
}
 
noinline static int sum_bits(const uint8_t*R bits) {
  int length = 0;
  //todo: if bits[0] == 0, we can start from 0
  for (int i = 1; i <= 16; i++)
    length += bits[i];
  return length;
} 
 
/* Emit a DHT marker */
inline static void emit_dht(c_state_ptr state, huff_table_t*R htbl, int index, boolean is_ac) {
  if (is_ac) {
    // htbl = cstate->ac_huff_tbls[index];
    index += 0x10; /* output index has AC bit set */
  } else {
    // htbl = cstate->dc_huff_tbls[index];
  }
    
  int length = sum_bits(htbl->bits);

  write_marker(state, M_DHT);
  write_short(state, length + 2 + 1 + 16);
  // write_marker_len(cstate, M_DHT, length + 2 + 1 + 16); 
  write_byte(state, index); 
  // hack: write bits and huffval by single copy
  write_block_unsafe(state, &htbl->bits[1], 16 + length);
}

noinline static int huff_decode(d_state_ptr state, const d_derived_table_t*R htbl) {


  // if ((nb = (htbl->lookup[look] >> HUFF_LOOKAHEAD)) <= HUFF_LOOKAHEAD) { 
  //   DROP_BITS(nb); 
  //   result = htbl->lookup[look] & ((1 << HUFF_LOOKAHEAD) - 1); 

  int byte = peek_bits(state, HUFF_LOOKAHEAD);
  if (byte == -1) return -1;
  uint16_t look = htbl->lookup[byte];
  int nb = look >> HUFF_LOOKAHEAD;
  if (nb <= HUFF_LOOKAHEAD) {
    state->bits_left -= nb; // drop the bits
    const int val = look & ((1 << HUFF_LOOKAHEAD) - 1);
    return val;
  }
  state->bits_left -= nb; // drop the bits

  int len = HUFF_LOOKAHEAD;
  int code = byte;

  while (code > htbl->maxcode[len]) {
    const int r = get_bit(state); 
    if (r < 0) return -1;
    code = (code << 1) | r;
    len++;
  }
  /* With garbage input we may reach the sentinel value l = 17. */
  if (len > 16) return -1; // JWRN_HUFF_BAD_CODE

  const int s = htbl->huffval[(int)(code + htbl->valoffset[len])];
  return s;
}

static int peek_be(const uint8_t* restrict p) {
  uint16_t w = *((uint16_t*)p);  
  return ((w & 0xff) << 8) | (w >> 8);
}

#define APP0_DATA_LEN 14  /* Length of interesting data in APP0 */
#define APP14_DATA_LEN 12 /* Length of interesting data in APP14 */
#define APPN_DATA_LEN 14  /* Must be the largest of the above!! */

/* Note that the result might not be a valid marker code,
 * but it will never be 0 or FF. */
// manual inlining of it increases size for some reason
noinline static int next_marker(d_state_ptr state) {
  for (;;) {
    // 1. read 0xff
    int c = read_byte(state);
    // note: libjpeg here was skipping any non-FF bytes
    // they are invalid according to specs
    if (c != 0xFF) return -1;
    // 2. skip duplicate 0xff's until non-zero
    /* This loop swallows any duplicate FF bytes.  Extra FFs are legal as
     * pad bytes, so don't count them in discarded_bytes.  We assume there
     * will not be so many consecutive FF bytes as to overflow a suspending
     * data source's input buffer. */
    // c = input_non_ff_byte(dstate);
    do {
      c = read_byte(state);
      // if (c == -1) return false; // can be skipped
    } while (c == 0xFF);
    if (c != 0) {
      // dstate->unread_marker = c;
      return c;
    }
  }
}

noinline static int div_size(int a, int b) {
  b *= DCTSIZE;
  return (a + b - 1) / b;
}

// baseline jpeg encoing

#define MAX_MCU_SIZE DCTSIZE2 * 16 // account for 0xff byte stuffing

static uint8_t*R emit_data_byte(uint8_t*R next_output_byte, uint8_t b) {
  // *((uint16_t*R)next_output_byte++) = (uint16_t)b;
  // return next_output_byte + (b == 0xff ? 1 : 0);    
  *next_output_byte++ = b;
  *next_output_byte = 0; // memory is zeroed already? - not in the current version
  return next_output_byte + (b == 0xff);
  // return next_output_byte + (b == 0xff ? 2 : 1);
  // next_output_byte[0] = b;
  // next_output_byte[1] = 0;
  // return next_output_byte - (-2 + ((uint8_t)(b) < 0xFF));
}

// #define FLUSH_START_SHIFT (BIT_BUF_SIZE - 8)

noinline static bit_buffer_t put_bits(c_state_ptr state, bit_buffer_t put_buffer, bit_buffer_t code, char size, int*R free_bits_ptr) {
  int free_bits = *free_bits_ptr;
  uint8_t* next_output_byte = state->next_output_byte;
  while (free_bits < size) { 
    // state->next_output_byte = emit_data_byte(state->next_output_byte, put_buffer >> (used - 8));
    free_bits += 8;
    next_output_byte = emit_data_byte(
      next_output_byte, 
      put_buffer >> ((BIT_BUF_SIZE - free_bits))
    );
  }
  state->next_output_byte = next_output_byte;
  put_buffer = (put_buffer << size) | code;
  free_bits -= size;
  /*free_bits -= size;
  if (free_bits < 0) {
    put_buffer = (put_buffer << (size + free_bits)) | (code >> -free_bits);
    for (int shift = FLUSH_START_SHIFT; shift >= 0; shift -= 8) {
      state->next_output_byte = emit_data_byte(state->next_output_byte, put_buffer >> shift);
    }
    free_bits += BIT_BUF_SIZE;
    put_buffer = code;
  } else {
    put_buffer = (put_buffer << size) | code;
  }*/
  *free_bits_ptr = free_bits;
  return put_buffer;
}

// noinline static int count_abs_bits(int temp, int min) {
//   if (temp == 0) { 
//     return min;
//   } else {
//     return 32 - __builtin_clz(abs(temp));
//   }
// }

noinline 
static int count_abs_bits(int temp) {
  if (temp == 0) { 
    return 0;
  } else {
    return 32 - __builtin_clz(abs(temp));
  }
}

// noinline static int count_abs_bits(int temp) {  
  // temp = abs(temp);
  // int nbits = 0;
  // while (temp) {
  //   nbits++;
  //   temp >>= 1;
  // }
  // return nbits;
// }

inline static int prepare_comp(component_t*R compptr, int ci, int*R mcu_membership, int blocks_in_mcu) {
  const int mcu_blocks = compptr->mcu_width * compptr->mcu_height;
  for (int i = mcu_blocks; i > 0; i--) {
    mcu_membership[blocks_in_mcu++] = ci;
  }
  return blocks_in_mcu;
}  

noinline static int prepare_comps(component_t*R*R cur_comp_info, int*R mcu_membership) {
  const int comps_in_scan = NCOMPS;
  int blocks_in_mcu = 0;
  for (int ci = 0; ci < comps_in_scan; ci++) {
    component_t* compptr = cur_comp_info[ci];
    if (comps_in_scan == 1) {
      compptr->mcu_width = 1;
      compptr->mcu_height = 1;
      // mcus_per_row = compptr->width_in_blocks;
    } else {
      compptr->mcu_width = compptr->h_samp_factor;
      compptr->mcu_height = compptr->v_samp_factor;
      // mcus_per_row = (int)jdiv_round_up((long)image_width, (long)(max_h_samp_factor * DCTSIZE));
    }
    blocks_in_mcu = prepare_comp(compptr, ci, mcu_membership, blocks_in_mcu);
  }
  return blocks_in_mcu;
}

noinline
static int find_next_smallest_nonzero_freq(freq_t freq[], int c_prev) {
  int c1 = -1;
  freq_t v = ~(freq_t)0;
  for (int i = 0; i <= 256; i++) {
    freq_t w = freq[i];
    if (w && w <= v && i != c_prev) {
      v = w;
      c1 = i;
    }
  }
  return c1;
}

/* Increment the codesize of everything in c1's tree branch */
__attribute__((always_inline))
inline
static int increment_branch_codesize(int*R codesize, int*R others, int c1) {
  codesize[c1]++;
  while (others[c1] >= 0) {
    c1 = others[c1];
    codesize[c1]++;
  }
  return c1;
}

noinline
static boolean jpeg_gen_optimal_table(c_state_ptr state, huff_table_t*R htbl, freq_t freq[], int ci, int k) {
  // huff_table_t *htbl = *htblptr;
  // long *freq = counts;

  int p, i, j;
  // int c1, c2;
  // long v;
  #define MAX_CLEN 32                             /* assumed maximum initial code length */
  uint8_t bits[MAX_CLEN + 1]; /* bits[k] = # of symbols with code length k */

  /* This algorithm is explained in section K.2 of the JPEG standard */

  int codesize[257];          /* codesize[k] = code length of symbol k */
  int others[257];            /* next symbol in current branch of tree */
  zero(bits, sizeof(bits));
  zero(codesize, sizeof(codesize));
  fill(others, 257 * sizeof(others[0]), -1);
  // for (int i = 0; i < 257; i++)
    // others[i] = -1; /* init links to empty */

  /* Including the pseudo-symbol 256 in the Huffman procedure guarantees
    * that no real symbol is given code-value of all ones, because 256
    * will be placed last in the largest codeword category. */
  freq[256] = 1; /* make sure 256 has a nonzero count */

  /* Huffman's basic algorithm to assign optimal code lengths to symbols */

  for (;;) {
    int c1 = find_next_smallest_nonzero_freq(freq, -1);
    int c2 = find_next_smallest_nonzero_freq(freq, c1);

    /* Done if we've merged everything into one frequency */
    if (c2 < 0) break;

    /* Else merge the two counts/trees */
    freq[c1] += freq[c2];
    freq[c2] = 0;

    c1 = increment_branch_codesize(codesize, others, c1);

    others[c1] = c2; /* chain c2 onto c1's tree branch */

    increment_branch_codesize(codesize, others, c2);
  }

  /* Now count the number of symbols of each code length */
  for (i = 0; i <= 256; i++) {
    const int cs = codesize[i];
    if (cs) {
      /* The JPEG standard seems to think that this can't happen, */
      /* but I'm paranoid... */
      if (cs > MAX_CLEN) return false; // ERREXIT(cstate, JERR_HUFF_CLEN_OVERFLOW);
      bits[cs]++;
    }
  }

  /* JPEG doesn't allow symbols with code lengths over 16 bits, so if the pure
    * Huffman procedure assigned any such lengths, we must adjust the coding.
    * Here is what Rec. ITU-T T.81 | ISO/IEC 10918-1 says about how this next
    * bit works: Since symbols are paired for the longest Huffman code, the
    * symbols are removed from this length category two at a time.  The prefix
    * for the pair (which is one bit shorter) is allocated to one of the pair;
    * then, skipping the BITS entry for that prefix length, a code word from the
    * next shortest nonzero BITS entry is converted into a prefix for two code
    * words one bit longer.
    */

  for (i = MAX_CLEN; i > 16; i--) {
    while (bits[i] > 0) {
      j = i - 2; /* find length of new prefix to be used */
      while (bits[j] == 0) j--;
      bits[i] -= 2;     /* remove two symbols */
      bits[i - 1]++;    /* one goes in this length */
      bits[j + 1] += 2; /* two new symbols in this length */
      bits[j]--;        /* symbol of this length is now a prefix */
    }
  }

  /* Remove the count for the pseudo-symbol 256 from the largest codelength */
  while (bits[i] == 0) /* find largest codelength still in use */
    i--;
  bits[i]--;

  /* Return final symbol counts (only for lengths 0..16) */
  memcpy(htbl->bits, bits, sizeof(htbl->bits));

  /* Return a list of the symbols sorted by code length */
  /* It's not real clear to me why we don't need to consider the codelength
    * changes made above, but Rec. ITU-T T.81 | ISO/IEC 10918-1 seems to think
    * this works.
    */
  p = 0;
  for (int i = 1; i <= MAX_CLEN; i++)
  {
    for (int j = 0; j <= 255; j++)
    {
      if (codesize[j] == i)
      {
        htbl->huffval[p] = (uint8_t)j;
        p++;
      }
    }
  }

  emit_dht(state, htbl, ci, k);

  return true;
}

inline static boolean check_sof(int length, int hsamp, int vsamp, int qtblno) {
  length -= 8;
  if (length != (NCOMPS * 3)) return false; // ERREXIT(dstate, JERR_BAD_LENGTH);
  // if (id == -1) return false;
  // if (c == -1) return false;
  if (hsamp <= 0 || hsamp == 3 || hsamp > MAX_SAMP_FACTOR ||
      vsamp <= 0 || vsamp == 3 || vsamp > MAX_SAMP_FACTOR)
    return false; // JERR_BAD_SAMPLING
  if (qtblno < 0 || qtblno >= NUM_QUANT_TBLS) return false;
  return true;
}

noinline 
static void input_block_unsafe(uint8_t* dest, d_state_ptr state, int size) {
  /* const */ uint8_t* restrict p = state->next_input_byte;
  memcpy(dest, p, size);
  state->next_input_byte = p + size;
}                

noinline
static boolean get_dht(d_state_ptr state, int length) {
  length -= 2;
  while (length > 16) {
    const int index = read_byte(state);
    // if (index == -1) return false; // checked by `if` below
    const int real_index = index & ~AC;
    if (real_index < 0 || real_index >= NUM_HUFF_TBLS) return false; // JERR_DHT_INDEX
    huff_table_t*R htbl;
    // if (index & 0x10) { 
    //   htbl = &state->ac_huff_tbls[real_index];
    // } else { 
    //   htbl = &state->dc_huff_tbls[real_index];
    // }
    htbl = &state->huff_tbls[index];
    htbl->read = true;
    uint8_t*R const bits = htbl->bits;
    uint8_t*R const huffval = htbl->huffval;

    bits[0] = 0;
    input_block_unsafe(&bits[1], state, 16);
    const int count = sum_bits(bits);

    length -= 1 + 16;

    /* Here we just do minimal validation of the counts to avoid walking
      * off the end of our table space.  jdhuff.c will check more carefully. */
    if (count > 256 || ((int32min_t)count) > length) return false; // ERREXIT(dstate, JERR_BAD_HUFF_TABLE);

    zero(huffval, 256);
    input_block_unsafe(huffval, state, count);

    length -= count;
  }
  if (length != 0) return false; // ERREXIT(dstate, JERR_BAD_LENGTH);
  return true;
}

inline static boolean copy_dqt(d_state_ptr state, int length) {
  state->next_input_byte -= 4; // copy marker and length too
  length += 2;
  while (length-- > 0) {
    int b = read_byte(state);
    if (b < 0) return false;
    if (!write_byte(state, b)) return false;
  }
  return true;
}

noinline static void write_sos_header(c_state_ptr state) {
  const uint8_t sos_header[14] = {
    0xff, M_SOS, 
    (2 * NCOMPS + 2 + 1 + 3) >> 8, 
    (2 * NCOMPS + 2 + 1 + 3) & 0xff,
    NCOMPS,
    0, (0 << 4) | 0,
    1, (1 << 4) | 1,
    2, (1 << 4) | 1,
    0, DCTSIZE2 - 1, 0,
  };
  write_block_unsafe(state, sos_header, 14);
}

noinline static int setup_mcu_buffer(
  coef_t*R mcu_data[10], 
  block_t*R whole_image[NCOMPS],
  component_t*R const cur_comp_info[MAX_COMPS_IN_SCAN], 
  int imcu_row, int mcu_col_num
) {
  const int comps_in_scan = NCOMPS;
  const int yoffset = 0;
  int blkn = 0; /* index of current DCT block within MCU */
  for (int ci = 0; ci < comps_in_scan; ci++) {
    const component_t*R compptr = cur_comp_info[ci];
    block_t*R subimage = whole_image[compptr->component_index];
    const int mcu_width = compptr->mcu_width;
    const int mcu_height = compptr->mcu_height;
    const int start_col = mcu_col_num * mcu_width;
    for (int yindex = yoffset; yindex < yoffset + mcu_height; yindex++) {     
      // select row
      const int blocks_per_row =  compptr->width_in_blocks;
      block_t*R block_array = (
        subimage + (
          imcu_row * compptr->v_samp_factor + yindex
        ) * blocks_per_row
      );        
      // seek to col
      block_t*R block_row = block_array + start_col;
      // fill
      for (int xindex = 0; xindex < mcu_width; xindex++) {
        mcu_data[blkn++] = *(block_row++);
      }
    }
  }  
  return blkn; // total blocks in mcu
}

// noinline static coef_t* ref_mcu_block(
//   block_t*R whole_image[NCOMPS],
//   int comps_in_scan, 
//   component_t*R const cur_comp_info[MAX_COMPS_IN_SCAN], 
//   int imcu_row, int yoffset, int mcu_col_num,
//   int blkn
// ) {
//   int i = 0; /* index of current DCT block within MCU */
//   for (int ci = 0; ci < comps_in_scan; ci++) {
//     const component_t*R compptr = cur_comp_info[ci];
//     const int mcu_width = compptr->mcu_width;
//     const int mcu_height = compptr->mcu_height;
//     const int start_col = mcu_col_num * mcu_width;
//     for (int yindex = yoffset; yindex < yoffset + mcu_height; yindex++) {     
//       block_t* subimage = whole_image[compptr->component_index];
//       // select row
//       const int blocks_per_row =  compptr->width_in_blocks;
//       block_t* block_array = (
//         subimage + (
//           imcu_row * compptr->v_samp_factor + yindex
//         ) * blocks_per_row
//       );        
//       // seek to col
//       block_t* block = block_array + start_col;
//       // fill
//       if (blkn >= i && blkn < i + mcu_width) {
//         return (coef_t*)(block + (blkn - i));
//       } else {
//         i += mcu_width;
//       }
//     }
//   }  
//   return NULL;
// }

// inline static void write_sof_header(state_t* cstate, int image_width, int image_height) {
//   uint8_t file_header[] = {
//     0xff, M_SOF0,
//     0, 3 * NCOMPS + 2 + 5 + 1, // length word
//     BITS_IN_JSAMPLE,
//     (int)image_height >> 8,
//     (int)image_height & 0xff,
//     (int)image_width >> 8,
//     (int)image_width & 0xff,
//     NCOMPS,
//   };
//   write_block_unsafe(cstate, file_header, sizeof(file_header));              
// }

/////////////////////////////////////////////

// #define JJ_FAIL 0

/* Note that jpeg_gen_optimal_table expects 257 entries in each table! */
typedef freq_t freq_tabs_t[2][257];

// returns real size of compressed size or 0 or negative number in case of any errors
// optimization takes right in the buffer provided
__attribute__((hot))
noinline int do_jj_optimize(uint8_t* restrict buf_orig, const int buf_size) {
  if ((buf_orig == NULL) | (buf_size <= 0)) return 0;
  // boolean params_check = (src_buf_orig == NULL) + (src_size == 0) + (dst_buf == NULL) + (dst_size == 0);
  // if (params_check) return -params_check;
  
  int result = 0;
      
  uint8_t* restrict const src_buf = alloc_mem(buf_size + 2048);
  if (src_buf == NULL) return 0;
  uint8_t* const dst_buf = src_buf;
  const int src_size = buf_size;
  // const int dst_size = buf_size;
  memcpy(src_buf, buf_orig, src_size);
    
  state_t*R const state = alloc_mem(sizeof(state_t));
  #define dstate state
  #define cstate state
  if (state == NULL) goto done;

  dstate->next_input_byte = (/* const */ uint8_t*)src_buf;
  dstate->input_end = src_buf + src_size;

  cstate->next_output_byte = dst_buf;
  // cstate->output_end = dst_buf + dst_size;
  // cstate->output_end = dstate->input_end;

  bit_buffer_t put_buffer; /* current bit accumulation buffer */
  int free_bits; /* # of bits available in it */ 
  // int input_scan_number = 0; /* Number of SOS markers seen so far */
  
  int image_width = 0; 
  int image_height;
  // int imcu_row; /* Number of iMCU rows completed */
  /* These parameters are never carried across datastreams, since they
   * are given in SOF/SOS markers or defined to be reset by SOI. */
  unsigned int restart_interval = 0; /* MCUs per restart interval, or 0 for no restart */

  int max_h_samp_factor; /* largest h_samp_factor */
  int max_v_samp_factor; /* largest v_samp_factor */
  int total_imcu_rows; /* # of iMCU rows in image */

  const int comps_in_scan = NCOMPS; /* # of JPEG components in this scan */
  component_t* cur_comp_info[MAX_COMPS_IN_SCAN];


  /* State variables made visible to other modules */
  boolean saw_soi = false;            /* found SOI? */
  boolean saw_sos = false;            /* true until first SOS is reached */
  boolean eoi_reached = false;        /* True when EOI has been consumed */
  // boolean saw_sof; 

  unsigned int restarts_to_go; /* MCUs left in this restart interval */

  void* restrict pool = alloc_mem(
    // src_size + 255 + // src_buf copy
    
    // sizeof(struct jpeg_decompress_struct) + // state
    
    NCOMPS * 2 * DCTSIZE2 * sizeof(int) + // coef_bits
    
    MAX_BLOCKS_IN_MCU * sizeof(int) + // mcu_membership
    // NUM_QUANT_TBLS * sizeof(JQUANT_TBL*) + // quant_tbl_ptrs
    MAX_COMPS_IN_SCAN * sizeof(int) + // last_dc_val
    // MAX_BLOCKS_IN_MCU * sizeof(int) + // 
    MAX_BLOCKS_IN_MCU * sizeof(block_t*) + // mcu_data
    NCOMPS * sizeof(block_t*) + // whole_image
    
    NCOMPS * sizeof(component_t) + // comp_info
    
    NUM_HUFF_TBLS * sizeof(d_derived_table_t*) + // dc_derived_tbls
    NUM_HUFF_TBLS * sizeof(d_derived_table_t*) + // ac_derived_tbls
    
    MAX_BLOCKS_IN_MCU * sizeof(d_derived_table_t*) + // dc_cur_tbls
    MAX_BLOCKS_IN_MCU * sizeof(d_derived_table_t*) + // ac_cur_tbls
    
    // MAX_BLOCKS_IN_MCU * sizeof(block_row_t) + // c_coef_dummy_buffer
    
    // NUM_HUFF_TBLS * sizeof(long*) + // dc_count_ptrs
    // NUM_HUFF_TBLS * sizeof(long*) + // ac_count_ptrs
    sizeof(freq_tabs_t) * 2 + // dc_count_ptrs, ac_count_ptrs
    
    // 17 + 256 + // dht_bits + dht_huffval
    
    // MAX_COMPS_IN_SCAN * sizeof(block_row_array_t) + // // consume_data_buffer
    
    // DCTSIZE2 * sizeof(int) + // newnz_pos
    
    0
  );
  if (pool == NULL) goto done;
  /* Current progression status.  coef_bits[c][i] indicates the precision
   * with which component c's DCT coefficient i (in zigzag order) is known.
   * It is -1 when no data has yet been received, otherwise it is the point
   * transform (shift) value for the most recent scan of the coefficient
   * (thus, 0 at completion of the progression).
   * This pointer is NULL when reading a non-progressive file. */
  // int coef_bits[NCOMPS * 2][DCTSIZE2]; /* -1 or current Al value for each coef */
  // fill(coef_bits, -1, sizeof(coef_bits));
  int* coef_bits = pool;
  pool += NCOMPS * 2 * DCTSIZE2 * sizeof(int);
  fill(coef_bits, NCOMPS * 2 * DCTSIZE2 * sizeof(int), -1);
  
  /* mcu_membership[i] is index in cur_comp_info of component owning i'th block in an MCU */
  // int mcu_membership[MAX_BLOCKS_IN_MCU];
  int* restrict mcu_membership = pool;
  pool += MAX_BLOCKS_IN_MCU * sizeof(int);
   
  // int last_dc_val[MAX_COMPS_IN_SCAN]; /* last DC coef for each component */
  int* restrict last_dc_val = pool; /* last DC coef for each component */
  pool += MAX_COMPS_IN_SCAN * sizeof(int);
  
  /* In single-pass modes, it's sufficient to buffer just one MCU.
   * We allocate a workspace of MAX_BLOCKS_IN_MCU coefficient blocks,
   * and let the entropy decoder write into that workspace each time.
   * In multi-pass modes, this array points to the current MCU's blocks
   * within the virtual arrays; it is used only by the input side. */
  // block_row_t d_coef_mcu_buffer[MAX_BLOCKS_IN_MCU];
  coef_t**R mcu_data = pool;
  pool += MAX_BLOCKS_IN_MCU * sizeof(block_t*);
  
  /* In multi-pass modes, we need a virtual block array for each component. */
  // block_row_array_t whole_image[NCOMPS];
  // block_t*R whole_image[NCOMPS];
  block_t**R whole_image = pool;
  pool += NCOMPS * sizeof(block_t*); 

  // component_t comp_info[NCOMPS];
  // fill(&comp_info[0], 0, NCOMPS * sizeof(comp_info[0]));
  component_t* restrict comp_info = pool;
  pool += NCOMPS * sizeof(component_t);
  
  // baseline mode
  
  /* Pointers to derived tables (these workspaces have image lifespan) */
  // d_derived_table_t* restrict dc_derived_tbls[NUM_HUFF_TBLS];
  // d_derived_table_t* restrict ac_derived_tbls[NUM_HUFF_TBLS];
  // /* Pointers to derived tables to be used for each block within an MCU */
  // d_derived_table_t* restrict dc_cur_tbls[MAX_BLOCKS_IN_MCU];
  // d_derived_table_t* restrict ac_cur_tbls[MAX_BLOCKS_IN_MCU];
  
  d_derived_table_t* restrict* restrict dc_derived_tbls = pool;
  pool += NUM_HUFF_TBLS * sizeof(d_derived_table_t*);
  d_derived_table_t* restrict* restrict ac_derived_tbls = pool;
  pool += NUM_HUFF_TBLS * sizeof(d_derived_table_t*);
  /* Pointers to derived tables to be used for each block within an MCU */
  d_derived_table_t* restrict* restrict dc_cur_tbls = pool; //alloc_mem(MAX_BLOCKS_IN_MCU * sizeof(d_derived_table_t*));
  pool += MAX_BLOCKS_IN_MCU * sizeof(d_derived_table_t*);
  d_derived_table_t* restrict* restrict ac_cur_tbls = pool;
  pool += MAX_BLOCKS_IN_MCU * sizeof(d_derived_table_t*);
  
  // fill((void*)dc_derived_tbls, 0, sizeof(dc_derived_tbls));
  // fill((void*)ac_derived_tbls, 0, sizeof(ac_derived_tbls));
  // fill((void*)dc_cur_tbls, 0, sizeof(dc_cur_tbls)); 
  // fill((void*)ac_cur_tbls, 0, sizeof(ac_cur_tbls));
  
  // fill(derived_tbls, 0, sizeof(derived_tbls));
  
  // d_derived_table_t* restrict ac_derived_tbl = NULL; /* active table during an AC scan */

  /* Statistics tables for optimization */
  // long* restrict dc_count_ptrs[NUM_HUFF_TBLS];
  // long* restrict ac_count_ptrs[NUM_HUFF_TBLS];
  // fill((void*)dc_count_ptrs, 0, NUM_HUFF_TBLS * sizeof(dc_count_ptrs[0])); 
  // fill((void*)ac_count_ptrs, 0, NUM_HUFF_TBLS * sizeof(ac_count_ptrs[0])); 
  freq_tabs_t* restrict dc_count_ptrs = pool;
  pool += sizeof(freq_tabs_t);
  freq_tabs_t* restrict ac_count_ptrs = pool;
  pool += sizeof(freq_tabs_t);
  
  // hack: dht_bits and dht_huffval should be right next to each other in memory
  // uint8_t* restrict dht_bits = pool;
  // pool += 17;
  // uint8_t* dht_huffval = pool;
  // pool += 256;
  
  
  for (;;) { // consume loop
    do { // consume markers loop
      {      
        int m = next_marker(dstate);
        /* Successfully processed marker, so reset state variable */
        if (m == M_SOI) {
          write_marker(cstate, M_SOI);
          if (saw_soi) goto done; // JERR_SOI_DUPLICATE
          saw_soi = true;
        } else if (m == M_EOI) {
          eoi_reached = true;
          break;
        } else if (((m >= M_RST0) & (m <= M_RST7)) | (m == M_TEM)) {
          /* these are all parameterless */
        } else { // marker with length or unknown marker
          int32min_t length = read_short(dstate);
          if (length < 0) goto done;
          if ((m == M_SOF0) || (m == M_SOF1)) { // no -2 
            /* baseline or extended sequential, huffman */ 
            { // get_sof
              if (image_width) goto done; // JERR_SOF_DUPLICATE
              // saw_sof = true;
              // if (length < 0) goto done;

              uint8_t* restrict p = dstate->next_input_byte;
              
              const int data_precision = *p++;
              if (data_precision != BITS_IN_JSAMPLE) goto done; 
              
              const int h = peek_be(p);
              p += 2;
              if ((h <= 0) | (h > JPEG_MAX_DIMENSION)) goto done;
              image_height = h;
              
              const int w = peek_be(p);
              p += 2;
              if ((w <= 0) | (w > JPEG_MAX_DIMENSION)) goto done;
              image_width = w;
              
              const int ncomps = *p++;
              if (ncomps != 3) goto done;

              // write_block_unsafe(cstate, p - (2 + 2 + 1 + 4 + 1), 2 + 2 + 1 + 4 + 1);

              // length -= 8;
              // if (length != (NCOMPS * 3)) goto done; // ERREXIT(dstate, JERR_BAD_LENGTH);

              int ci;
              component_t*R compptr;
              for (ci = 0, compptr = comp_info; ci < NCOMPS; ci++, compptr++) {
                compptr->component_index = ci; // used in multi-scan mode
                const int id = *p;
                // if (id == -1) goto done; // check_sof
                compptr->component_id = id;
                *p++ = ci; // encoder uses 0-based ids
                const int c = *p++;
                // if (c == -1) goto done; // check_sof
                const int hsamp = c >> 4;
                const int vsamp = c & 15;
                // if (hsamp <= 0 || hsamp > MAX_SAMP_FACTOR ||
                //     vsamp <= 0 || vsamp > MAX_SAMP_FACTOR)
                //   goto done; // JERR_BAD_SAMPLING // check_sof
                compptr->h_samp_factor = hsamp;
                compptr->v_samp_factor = vsamp;
                const int qtblno = *p++;
                // if (qtblno < 0 || qtblno >= NUM_QUANT_TBLS) goto done; // check_sof
                if (!check_sof(length, hsamp, vsamp, qtblno)) goto done;
                // compptr->quant_tbl_no = qtblno;
                // write_byte(cstate, ci); // component id
                // write_byte(cstate, (hsamp << 4) + vsamp);
                // write_byte(cstate, qtblno);                
              }
              dstate->next_input_byte = p;
              
              #define SOF_HEADER_SIZE (2 + 2 + 1 + 4 + 1 + NCOMPS * 3)
              write_block_unsafe(state, p - SOF_HEADER_SIZE, SOF_HEADER_SIZE);
              
            } // get_sof
          } else if (m == M_SOS) { // no -2
            // get_sos
            {
              if (!image_width) goto done; //  JERR_SOS_NO_SOF
              
              /* const */ uint8_t* restrict p = dstate->next_input_byte;
              
              const int n = *p++;
              
              if ((length != (n * 2 + 6)) | (n != NCOMPS))
                goto done; // JERR_BAD_LENGTH

              /* Collect the component-spec parameters */

              for (int i = MAX_COMPS_IN_SCAN - 1; i >= 0; i--)
                cur_comp_info[i] = NULL;

              for (int i = 0; i < n; i++) {
                const int id = *p++; // id

                component_t *compptr;
                int ci;
                for (ci = 0, compptr = &comp_info[0]; ci < NCOMPS; ci++, compptr++) {
                  if (id == compptr->component_id) {
                    cur_comp_info[i] = compptr;
                    const int c = *p++; // dc/ac
                    compptr->dc_tbl_no = c >> 4;
                    compptr->ac_tbl_no = c & 15;
                    // dc/ac_tbl_no are checked by jpeg_make_d_derived_tbl
                    break;
                  }
                }
                
                if (ci == NCOMPS) goto done; // JERR_BAD_COMPONENT_ID

                /* This CSi (cc) should differ from the previous CSi */
                // for (int pi = 0; pi < i; pi++) {
                for (int pi = i - 1; pi >= 0; pi--) {
                  if (cur_comp_info[pi] == compptr) goto done; // ERREXIT1(dstate, JERR_BAD_COMPONENT_ID, cc);
                }
              }

              if (*p++ != 0) goto done; // spectral start should be 0
              if (*p++ != 63) goto done;
              
              // const int Ah_Al = *p++;
              // Ah = Ah_Al >> 4;
              // Al = Ah_Al & 15;
              p++; // skip ah, al
              
              dstate->next_input_byte = p;  

              /* Prepare to scan data & restart markers */
              // next_restart_num = 0;
            }
            
            if (!saw_sos) {  /* 1st SOS */
            // if (input_scan_number++ == 0) {
              // if (input_scan_number > 255) goto done;            
              saw_sos = true; // we've now seen SOI, SOS and - hopefully - SOF

              // if (!initial_setup_d(dstate)) goto done;
              // static boolean initial_setup_d(d_state_ptr dstate)
              {
                component_t*R compptr = comp_info;

                /* Compute maximum sampling factors; check factor validity */
                max_h_samp_factor = 1;
                max_v_samp_factor = 1;
                for (int ci = 0; ci < NCOMPS; ci++, compptr++) {
                  const int hsamp = compptr->h_samp_factor;
                  max_h_samp_factor = MAX(max_h_samp_factor, hsamp);
                  const int vsamp = compptr->v_samp_factor;
                  max_v_samp_factor = MAX(max_v_samp_factor, vsamp);
                }
                
                /* Compute number of fully interleaved MCU rows. */
                total_imcu_rows = div_size(image_height, max_v_samp_factor);                

                /* Compute dimensions of components */
                compptr = comp_info;
                for (int ci = 0; ci < NCOMPS; ci++) {
                  const int height_in_blocks = div_size(image_height * compptr->v_samp_factor, max_v_samp_factor);
                  compptr->width_in_blocks = div_size(image_width * compptr->h_samp_factor, max_h_samp_factor);
                  block_t*R subimage = alloc_block_array(compptr->width_in_blocks, height_in_blocks);;
                  if (subimage == NULL) goto done;
                  whole_image[ci] = subimage;
                  compptr++;
                }
              }
            } else { /* 2nd or later SOS marker */
              goto done; // no multiple scans support
            }
            
            break;
            // end of M_SOS
            
          } else if (m == M_DHT) {
            if (!get_dht(dstate, length)) goto done;
          } else if (m == M_DRI) { // no -2
            // if (!get_dri(dstate)) goto done;
            if (length != 4) goto done; // JERR_BAD_LENGTH
            int bb = read_short(dstate);
            if (bb < 0) goto done;
            restart_interval = bb;
          } else if (m == M_DQT) {
            if (!copy_dqt(dstate, length)) goto done; 
          } else if (m == M_APP0 || m == M_DAC || m >= M_APP1 && m <= M_APP15 || m == M_COM || m == M_DNL) {
            /* Skip over an unknown or uninteresting variable-length marker */
            length -= 2;
            // if (length > dstate->bytes_in_buffer) goto done;
            // dstate->next_input_byte += length;
            // dstate->bytes_in_buffer -= length;
            /* const */ uint8_t*R p = dstate->next_input_byte + length;
            dstate->next_input_byte = p;
            if (p >= dstate->input_end) goto done;
          } else { 
            goto done;
          }
        }

      } // consume_markers

    } while (true /*retcode != JPEG_REACHED_SOS && retcode != JPEG_REACHED_EOI*/);

    //inv: saw_sos or eoi_reached
    if (!saw_sos) goto done; // headers not parsed
    // if (input_scan_number == 0) goto done; // headers not parsed
    if (eoi_reached) break; // break main decompression loop

    //inv: new sos marker reached    
    // if (!per_scan_setup_d(dstate)) goto done;
    const int blocks_in_mcu = prepare_comps(cur_comp_info, mcu_membership);
    if (blocks_in_mcu > MAX_BLOCKS_IN_MCU) goto done;

    // make derived tables
    // const bool is_DC_band = (Ss == 0);
     
    for (int ci = comps_in_scan - 1; ci >= 0; ci--) {
      component_t* restrict compptr = cur_comp_info[ci];
      d_derived_table_t* restrict * restrict pdtbl;
      int dc_tbl_no = compptr->dc_tbl_no;
      int ac_tbl_no = compptr->ac_tbl_no;
      
      pdtbl = (d_derived_table_t **)(dc_derived_tbls) + dc_tbl_no;
      if (!make_derived_table(dstate, true, dc_tbl_no, pdtbl))
        goto done;
      pdtbl = (d_derived_table_t **)(ac_derived_tbls) + ac_tbl_no;
      if (!make_derived_table(dstate, false, ac_tbl_no, pdtbl))
        goto done;
    
      last_dc_val[ci] = 0;
    }    
    /* Precalculate decoding info for each block in an MCU of this scan */
    // for (int blkn = 0; blkn < blocks_in_MCU; blkn++) {
    for (int blkn = 0; blkn < blocks_in_mcu; blkn++) {
      int ci = mcu_membership[blkn];
      component_t* restrict compptr = cur_comp_info[ci];
      /* Precalculate which table to use for each block */
      dc_cur_tbls[blkn] = dc_derived_tbls[compptr->dc_tbl_no]; //ors: needed?
      ac_cur_tbls[blkn] = ac_derived_tbls[compptr->ac_tbl_no];
    }    
    dstate->bits_left = 0;
    restarts_to_go = restart_interval;    

    {    
      {
        /* In an interleaved scan, an MCU row is the same as an iMCU row.
          * In a noninterleaved scan, an iMCU row has v_samp_factor MCU rows.
          * But at the bottom of the image, process only what's left. */
        for (int imcu_row = 0; imcu_row < total_imcu_rows; imcu_row++) {
          /* Loop to process one whole iMCU row */
          {
            int mcus_per_row; /* # of MCUs across the image */
            mcus_per_row = div_size(image_width, max_h_samp_factor);
            for (int mcu_col_num = 0; mcu_col_num < mcus_per_row; mcu_col_num++) {
              // coef_t* mcu_data[10];
              const int blocks_in_mcu = setup_mcu_buffer(mcu_data, 
                whole_image,
                cur_comp_info, 
                imcu_row, mcu_col_num);
              // if (blocks_in_mcu > MAX_BLOCKS_IN_MCU) goto done;                  
            
              // first, process restart marker if needed; may have to suspend 
              if ((restart_interval != 0) & (restarts_to_go == 0)) {
                // printf("process restart marker \n");
                // process_restart
                // if ((read_short(dstate) & 0xfff8) != 0xffd0) goto done;
                if ((read_short(dstate) & 0xfff8) != 0xffd0) goto done;
                dstate->bits_left = 0;
                zero(last_dc_val, NCOMPS * sizeof(int));
                // for (int ci = 0; ci < NCOMPS; ci++)
                  // last_dc_val[ci] = 0;
                restarts_to_go = restart_interval; 
              }
              /* Account for restart interval (no-op if not using restarts) */
              restarts_to_go--;

              /* MCU decoding */
              {
                { // decode_mcu_DC_first or decode_mcu
                  { 
                    /* Outer loop handles each block in the MCU */
                    const int blocks_in_mcu_copy = blocks_in_mcu;
                    for (int blkn = 0; blkn < blocks_in_mcu_copy; blkn++) {
                      coef_t* block = (mcu_data[blkn]);
                      const int ci = mcu_membership[blkn];
                      
                      const d_derived_table_t*R tbl = dc_cur_tbls[blkn];
                      const d_derived_table_t*R actbl = ac_cur_tbls[blkn];

                      /* Decode a single block's worth of coefficients */

                      /* Section F.2.2.1: decode the DC coefficient difference */
                      int s = huff_decode(dstate, tbl); 
                      if (s < 0) goto done;
                      // printf("dc %4d \n", s);
                      // if (s) {
                      s = get_bits_extend(dstate, s);
                      // }

                      /* Convert DC difference to actual value, update last_dc_val */
                      /* Certain malformed JPEG images produce repeated DC coefficient
                        * differences of 2047 or -2047, which causes last_dc_val[ci] to
                        * grow until it overflows or underflows a 32-bit signed integer.  This
                        * behavior is, to the best of our understanding, innocuous, and it is
                        * unclear how to work around it without potentially affecting
                        * performance. */                      
                      const int last_dc_val_val = last_dc_val[ci];
                      // if (
                      //   ((last_dc_val_val >= 0) && (s > INT_MAX - last_dc_val_val)) ||
                      //   ((last_dc_val_val < 0) && (s < INT_MIN - last_dc_val_val))
                      // )
                      //   goto done; // ERREXIT(dstate, JERR_BAD_DCT_COEF);
                      // simplified version:
                      // if (abs(last_dc_val_val) > INT_MAX / 2) goto done;
                      if (last_dc_val_val > INT_MAX / 2 || last_dc_val_val < -INT_MAX / 2) goto done;
                      // printf("ci %d decode s %d last_dc_val %d ", ci, s, last_dc_val_val);
                      s += last_dc_val_val;
                      // printf("%d\n", s);
                      last_dc_val[ci] = s;
                      
                      /* Output the DC coefficient (assumes jpeg_natural_order[0] = 0) */
                      block[0] = (coef_t)s;
                      
                      /* Section F.2.2.2: decode the AC coefficients */
                      /* Since zeroes are skipped, output area must be cleared beforehand */
                      for (int k = 1; k < DCTSIZE2; k++) {
                        int s = huff_decode(dstate, actbl); 
                        if (s < 0) goto done;
                        const int r = s >> 4;
                        s &= 15;
                        if (s) {
                          k += r;
                          s = get_bits_extend(dstate, s);
                          /* Output coefficient in natural (dezigzagged) order.
                            * Note: the extra entries in jpeg_natural_order[] will save us
                            * if k >= DCTSIZE2, which could happen if the data is corrupted. */
                          block[jpeg_natural_order[k]] = (coef_t)s;
                        } else {
                          if (r != 15) break;
                          k += 15;
                        }
                      }                      
                      
                    }
                  }
                } 
              } 
              
              // check bounds every now and then
              if (dstate->next_input_byte > dstate->input_end) goto done;
            }
          } // for yoffset
        } // consume_data for loop until scan is completed
      }
    }
  } // consume loop

  // image is decoded //

  // jpeg_finish_compress(cstate);
  {
    boolean gather_statistics = true; 

    // fill `cur_comp_info` for prepare_comps
    for (int ci = 0; ci < NCOMPS; ci++) {
      cur_comp_info[ci] = &comp_info[ci];
    }      
    // const int blocks_in_mcu = prepare_comps(cur_comp_info, NCOMPS, mcu_membership);
    // if (blocks_in_mcu > MAX_BLOCKS_IN_MCU) goto done;      
    prepare_comps(cur_comp_info, mcu_membership);

    while (true) {

      // METHODDEF(void) start_pass_huff(c_state_ptr cstate, bool gather_statistics)
      
      zero(last_dc_val, NCOMPS * sizeof(int));
      /* Initialize bit buffer to empty */
      put_buffer = 0;
      free_bits = BIT_BUF_SIZE;

      for (int imcu_row = 0; imcu_row < total_imcu_rows; imcu_row++) {
        //  METHODDEF(bool) jctrans_compress_output(c_state_ptr cstate)
        {
          const int mcus_per_row = div_size(image_width, max_h_samp_factor);          
          /* process one whole iMCU row */
          {
            for (int mcu_col_num = 0; mcu_col_num < mcus_per_row; mcu_col_num++) {
            
              /* Construct list of pointers to DCT blocks belonging to this MCU */
              // coef_t* mcu_data[10];              
              const int blocks_in_mcu = setup_mcu_buffer(mcu_data, 
                whole_image,
                cur_comp_info, 
                imcu_row, mcu_col_num);    
              if (blocks_in_mcu > MAX_BLOCKS_IN_MCU) goto done;                
              
              for (int blkn = 0; blkn < blocks_in_mcu; blkn++) {
                int ci = mcu_membership[blkn];
                
                // const coef_t* mcu_start = ref_mcu_block(
                //   whole_image,
                //   NCOMPS, cur_comp_info, 
                //   imcu_row, 0, mcu_col_num, blkn);  
                // const coef_t* block = mcu_start;
                const coef_t*R block = mcu_data[blkn];
                
                const int last_cur_dc_val = last_dc_val[ci];
                // int temp = block[0] - last_cur_dc_val;
                // if (gather_statistics) printf("ci %d encode block[0] %d last_dc_val %d ", ci, block[0], last_dc_val[ci]);
                // int temp = block[0];// - last_cur_dc_val;
                int temp = block[0] - last_cur_dc_val;
                // if (gather_statistics) printf("temp %d\n", temp);
                if (gather_statistics) { 
                  // htest_one_block(cstate, MCU_data[blkn][0], last_dc_val[ci],
                  //   dc_count_ptrs[compptr->dc_tbl_no],
                  //   ac_count_ptrs[compptr->ac_tbl_no])
                  /* Process a single block's worth of coefficients */
                  // static boolean htest_one_block(c_state_ptr cstate, JCOEFPTR block, int last_dc_val, long dc_counts[], long ac_counts[])
                  {
                    const int c_tbl_no = ci != 0;
                    freq_t*R dc_counts = (*dc_count_ptrs)[c_tbl_no];
                    freq_t*R ac_counts = (*ac_count_ptrs)[c_tbl_no];
                    int nbits = count_abs_bits(temp);
                    // int nbits = count_abs_bits(temp, 0);
                    /* Check for out-of-range coefficient values.
                      * Since we're encoding a difference, the range limit is twice as much. */
                    if (nbits > MAX_COEF_BITS + 1) return false; // JERR_BAD_DCT_COEF
                    dc_counts[nbits]++;
                    /* Encode the AC coefficients per section F.1.2.2 */
                    int r = 0; /* r = run length of zeros */
                    for (int k = 1; k < DCTSIZE2; k++) {
                      int temp = block[jpeg_natural_order[k]];
                      if (temp == 0) {
                        r++;
                      } else {
                        /* Find the number of bits needed for the magnitude of the coefficient */
                        const int nbits = MAX(count_abs_bits(temp), 1);
                        // const int nbits = count_abs_bits(temp, 1);
                        if (nbits > MAX_COEF_BITS) return false; // JERR_BAD_DCT_COEF
                        /* if run length > 15, must emit special run-length-16 codes (0xF0) */
                        while (r > 15) {
                          ac_counts[0xF0]++;
                          r -= 16;
                        }
                        /* Count Huffman symbol for run length / number of bits */
                        ac_counts[(r << 4) + nbits]++;
                        r = 0;
                      }
                    }
                    /* If the last coef(s) were zero, emit an end-of-block code */
                    if (r > 0) ac_counts[0]++;
                  }
                } else {
                  // encode_one_block(cstate, MCU_data[blkn][0], last_dc_val[ci],
                  //   dc_derived_tbls[compptr->dc_tbl_no],
                  //   ac_derived_tbls[compptr->ac_tbl_no])
                  // static boolean encode_one_block(c_state_ptr cstate, JCOEFPTR block, int last_dc_val, c_derived_table_t *dctbl, c_derived_table_t *actbl)
                  {
                    const int c_tbl_no = ci != 0;
                    c_derived_table_t *dctbl = dc_derived_tbls[c_tbl_no];
                    c_derived_table_t *actbl = ac_derived_tbls[c_tbl_no];
                    c_derived_table_t *tbl = dctbl;

                    // bit_buffer_t put_buffer = put_buffer; //L-mark
                    // int free_bits = free_bits; // L-mark

                    // uint8_t *start_output_byte = cstate->next_output_byte;

                    /* Encode the AC coefficients per section F.1.2.2 */
                    {
                      int r = 0; /* r = run length of zeros */
                      // for (int k = 1; k < DCTSIZE2; k++) {
                      for (int k = 0; k < DCTSIZE2; k++)
                      {
                        int temp0 = block[jpeg_natural_order[k]];
                        int temp;
                        int nbits;
                        int code;
                        if (k == 0) {
                          temp = temp0 - last_cur_dc_val;
                          
                          nbits = temp >> (CHAR_BIT * sizeof(int) - 1);
                          temp += nbits;
                          nbits ^= temp;
                          
                          nbits = JPEG_NBITS(nbits);
                          
                          code = nbits;
                          // PUT_CODE(dctbl->ehufco[nbits], dctbl->ehufsi[nbits])
                        } else {
                          temp = temp0;
                          if (temp == 0) {
                            r += 16;
                          } else {
                            /* Branch-less absolute value, bitwise complement, etc., same as above */
                            
                            nbits = temp >> (CHAR_BIT * sizeof(int) - 1);
                            temp += nbits;
                            nbits ^= temp;
                            
                            nbits = JPEG_NBITS_NONZERO(nbits);
                            
                            /* if run length > 15, must emit special run-length-16 codes (0xF0) */
                            while (r >= 16 * 16) {
                              r -= 16 * 16;
                              // put_bits(actbl->ehufco[0xf0], actbl->ehufsi[0xf0])
                              put_buffer = put_bits(cstate, put_buffer, actbl->ehufco[0xf0], actbl->ehufsi[0xf0], &free_bits);
                            }
                            /* Emit Huffman symbol for run length / number of bits */
                            r += nbits;
                            code = r;
                            r = 0;
                          }
                        }
                        if (k == 0 || (k != 0 && temp0 != 0)) {
                          // PUT_CODE(tbl->ehufco[code], tbl->ehufsi[code]);
                          // int code0 = code;
                          const int co = tbl->ehufco[code];
                          const int si = tbl->ehufsi[code];

                          temp &= (((int32min_t)1) << nbits) - 1;
                          temp |= co << nbits;
                          nbits += si;
                          
                          // put_bits(temp, nbits);
                          put_buffer = put_bits(cstate, put_buffer, temp, nbits, &free_bits);
                          tbl = actbl; // switch
                        }
                      }
                      if (r > 0) { /* If the last coef(s) were zero, emit an end-of-block code */
                        // put_bits(actbl->ehufco[0], actbl->ehufsi[0])
                        put_buffer = put_bits(cstate, put_buffer, actbl->ehufco[0], actbl->ehufsi[0], &free_bits);                                            
                      }
                    }
                  }
                } // block encoded
                
                // last_dc_val[ci] = mcu_data[blkn][0][0];
                last_dc_val[ci] = block[0];
                
                if (state->next_output_byte >= state->output_end) {
                  // printf("state->next_output_byte > state->output_end 3\n");
                  goto done;
                }
              }                    
            }
            /* Completed an MCU row, but perhaps not an iMCU row */
            // mcu_ctr = 0;
          }
          /* Completed the iMCU row, advance counters for next one */
        } // jctrans_compress_output

      } // for imcu row

      if (gather_statistics)
      {
        // finish_pass_gather(cstate);
        // METHODDEF(void) finish_pass_gather(c_state_ptr cstate)
        {
          // for (int ci = 0; ci < NCOMPS; ci++) {
          for (int ci = 0; ci < 2; ci++) {
            
            // huff_table_t*R htbl = &cstate->dc_huff_tbls[ci];
            huff_table_t*R htbl = &cstate->huff_tbls[ci];
            freq_t* counts = (*dc_count_ptrs)[ci];

            for (int k = 0; k < 2; k++) { // dc/ac
              if (!jpeg_gen_optimal_table(cstate, htbl, counts, ci, k)) goto done;

              // htbl = &cstate->ac_huff_tbls[ci];
              htbl = &cstate->huff_tbls[ci + AC];
              counts = (*ac_count_ptrs)[ci];
            }
            
          }
        }
        
        { 
          // write_scan_header
          // only possible on the second pass when tables are optimized
          // for (int ci = 0; ci < 2; ci++) {
          //   emit_dht(cstate, ci, false);
          //   emit_dht(cstate, ci, true);
          // }
          write_sos_header(cstate);
        }    
        
        for (int ci = 0; ci < 2; ci++) {
          const int c_tbl_no = ci != 0;
          const int dctbl = c_tbl_no;
          const int actbl = c_tbl_no;
          /* Compute derived values for Huffman tables */
          /* We may do this more than once for a table, but it's not expensive */
          if (!make_derived_table(cstate, true, dctbl, &dc_derived_tbls[dctbl])) goto done;
          if (!make_derived_table(cstate, false, actbl, &ac_derived_tbls[actbl])) goto done;            
        }
        
      } else {
        break;
      }

      gather_statistics = false;
    } // while (!cstate->is_last_pass)
    
    // finish_pass_huff(cstate);
    {
      // flush (BIT_BUF_SIZE - free_bits) bits rounded up to bytes
      //note: this fills partial byte with ones 
      //note: two calls to avoid negative shift problem
      for (int i = 1; i >= 0; i--) {
        put_buffer = put_bits(cstate, put_buffer, BIT_BUF_FILL >> (BIT_BUF_SIZE / 2), BIT_BUF_SIZE / 2, &free_bits);
      }
    }
    write_marker(cstate, M_EOI);
    // {
      // uint8_t* p = cstate->next_output_byte;
      // zero(p, cstate->output_end - p);
    // }
    
    const int result_size = cstate->next_output_byte - dst_buf;
    memcpy(buf_orig, src_buf, result_size);
    result = result_size;
    goto done;
  }

done:
  free_mem();
  return result;
}

int test(const char *src_fn, const char *dest_fn)
{
  // const char* src_fn = "test.jpg";
  // const char* src_fn = "test2.jpg";
  // const char* dest_fn = "result.jpg";
  // int retval = 0;
  uint8_t *src_buf = NULL;
  unsigned long src_size;
  const int pad = 65536;

  {
    FILE *src_file = NULL;
    if ((src_file = fopen(src_fn, "rb")) == NULL) return 1;
    if (fseek(src_file, 0, SEEK_END) < 0 || ((src_size = ftell(src_file)) < 0) ||
        fseek(src_file, 0, SEEK_SET) < 0)
      return 1;
    if (src_size == 0) return 1;
    if ((src_buf = (unsigned char *)malloc(src_size + pad)) == NULL)
      return 1;
    if (fread(src_buf, src_size, 1, src_file) < 1)
      return 1;
    fclose(src_file);
  }

  // prealloc mem ?
  // unsigned long dst_size = src_size * 2;
  // unsigned char *dst_buf = NULL;
  // if ((dst_buf = (unsigned char *)calloc(1, dst_size * 2 + 1024)) == NULL)
  //   return 1;

  int real_size = do_jj_optimize(src_buf, src_size + pad);
  // printf("real_size %d \n", real_size);
  bool ok = real_size > 0;

  if (ok) {
    uint8_t* dst_buf = src_buf;
    FILE *dst_file = NULL;
    if ((dst_file = fopen(dest_fn, "wb")) == NULL) {
      printf("opening output file");
      return 1;
    }
    if (fwrite(dst_buf, real_size, 1, dst_file) < 1) {
      printf("writing output file");
      return 1;
    }
    fclose(dst_file);
  }

  free(src_buf);
  // free(dst_buf);

  if (ok)
  {
    // printf("ok\n");
    return 0;
  }
  else
  {
    printf("error\n");
    return 1;
  }
}

int main(int __argc, char **__argv)
{
  // test("test8.jpg", "result8.jpg");
  test("test.jpg", "result.jpg");
  // test("test2.jpg", "result2.jpg"); // prog
  // test("test3.jpg", "result3.jpg"); // prog restartss
  // test("test_odd.jpg", "result_odd.jpg"); // odd size
  test("test4br.jpg", "result4.jpg"); // baseline restarts
}

/////////////////////////////////////////////
/*EMSCRIPTEN_KEEPALIVE
int jj_optimize(int a, int b) {
  return a + b;
}
*/
