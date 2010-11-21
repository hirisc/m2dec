/**VLD DCT coefficient tables.
 */

#ifdef __cplusplus
extern "C" {
#endif

enum {
	VLD_BITLEN = 7
};

/** Mapping of Q scale indicator and Q scale.
 */
static const int m2d_q_mapping[2][32] = {
	{
		2, 2, 4, 6, 8, 10, 12, 14,
		16, 18, 20, 22, 24, 26, 28, 30,
		32, 34, 36, 38, 40, 42, 44, 46,
		48, 50, 52, 54, 56, 58, 60, 62
	},
	{
		1, 1, 2, 3, 4, 5, 6, 7,
		8, 10, 12, 14, 16, 18, 20, 22,
		24, 28, 32, 36, 40, 44, 48, 52,
		56, 64, 72, 80, 88, 96, 104, 112
	}
};

/** No need to inverse zigzag-scan for looking up.
 */
static const uint8_t m2d_qmat_default[2][64] = {
#if 1
	{
		8, 16, 19, 22, 26, 27, 29, 34,
		16, 16, 22, 24, 27, 29, 34, 37,
		19, 22, 26, 27, 29, 34, 34, 38,
		22, 22, 26, 27, 29, 34, 37, 40,
		22, 26, 27, 29, 32, 35, 40, 48,
		26, 27, 29, 32, 35, 40, 48, 58,
		26, 27, 29, 34, 38, 46, 56, 69,
		27, 29, 35, 38, 46, 56, 69, 83,
	},
#else
	{
		8, 16, 16, 19, 16, 19, 22, 22,
		22, 22, 22, 22, 26, 24, 26, 27,
		27, 27, 26, 26, 26, 26, 27, 27,
		27, 29, 29, 29, 34, 34, 34, 29,
		29, 29, 27, 27, 29, 29, 32, 32,
		34, 34, 37, 38, 37, 35, 35, 34,
		35, 38, 38, 40, 40, 40, 48, 48,
		46, 46, 56, 56, 58, 69, 69, 83,
	},
#endif
	{
		16, 16, 16, 16, 16, 16, 16, 16,
		16, 16, 16, 16, 16, 16, 16, 16,
		16, 16, 16, 16, 16, 16, 16, 16,
		16, 16, 16, 16, 16, 16, 16, 16,
		16, 16, 16, 16, 16, 16, 16, 16,
		16, 16, 16, 16, 16, 16, 16, 16,
		16, 16, 16, 16, 16, 16, 16, 16,
		16, 16, 16, 16, 16, 16, 16, 16,
	}
};

static const int8_t m2d_zigzag[2][64] = {
	{
		0,  1, 8, 16, 9, 2, 3, 10,
		17, 24, 32, 25, 18, 11, 4, 5,
		12, 19, 26, 33, 40, 48, 41, 34,
		27, 20, 13, 6, 7, 14, 21, 28,
		35, 42, 49, 56, 57, 50, 43, 36,
		29, 22, 15, 23, 30, 37, 44, 51,
		58, 59, 52, 45, 38, 31, 39, 46,
		53, 60, 61, 54, 47, 55, 62, 63
	},
	{
		0, 8, 16, 24, 1, 9, 2, 10,
		17, 25, 32, 40, 48, 56, 57, 49,
		41, 33, 26, 18, 3, 11, 4, 12,
		19, 27, 34, 42, 50, 58, 35, 43,
		51, 59, 20, 28, 5, 13, 6, 14,
		21, 29, 36, 44, 52, 60, 37, 45,
		53, 61, 22, 30, 7, 15, 23, 31,
		38, 46, 54, 62, 39, 47, 55, 63,
	},
};

#if 0
#include <stdio.h>
/**Inverse zigzag scan default matrices. Just for development only.
 */
void qmat_trans() {
	int tbl;

	printf("const uint8_t m2d_qmat_default[2][64] = {\n");
	for (tbl = 0; tbl < 2; ++tbl) {
		int i = 0;
		printf("\t{\n");
		do {
			int j;
			printf("\t\t");
			for (j = 0; j < 8; ++j) {
				printf("%d,", m2d_qmat_default[tbl][m2d_zigzag[0][i]]);
				++i;
			}
			printf("\n");
		} while (i < 64);
		printf("\t},\n");
	}
	printf("};\n\n");
}
#endif

/** This is generated from standard document using makevld.rb. */
static const vlc_dct_t m2d_dct_table0_bit7[392] = {
	{128, 10, 0}, {352, 4, 0}, {-1, 0, 6}, {-1, 0, 6},
	{368, 1, 0}, {370, 1, 0}, {372, 1, 0}, {374, 1, 0},
	{7, 2, 7}, {7, 3, 7}, {6, 2, 7}, {6, 3, 7},
	{1, 4, 7}, {1, 5, 7}, {5, 2, 7}, {5, 3, 7},
	{376, 2, 0}, {380, 2, 0}, {384, 2, 0}, {388, 2, 0},
	{0, 6, 6}, {0, 6, 6}, {0, 7, 6}, {0, 7, 6},
	{4, 2, 6}, {4, 2, 6}, {4, 3, 6}, {4, 3, 6},
	{3, 2, 6}, {3, 2, 6}, {3, 3, 6}, {3, 3, 6},
	{0, 4, 5}, {0, 4, 5}, {0, 4, 5}, {0, 4, 5},
	{0, 5, 5}, {0, 5, 5}, {0, 5, 5}, {0, 5, 5},
	{2, 2, 5}, {2, 2, 5}, {2, 2, 5}, {2, 2, 5},
	{2, 3, 5}, {2, 3, 5}, {2, 3, 5}, {2, 3, 5},
	{1, 2, 4}, {1, 2, 4}, {1, 2, 4}, {1, 2, 4},
	{1, 2, 4}, {1, 2, 4}, {1, 2, 4}, {1, 2, 4},
	{1, 3, 4}, {1, 3, 4}, {1, 3, 4}, {1, 3, 4},
	{1, 3, 4}, {1, 3, 4}, {1, 3, 4}, {1, 3, 4},
	{-1, 3, 2}, {-1, 3, 2}, {-1, 3, 2}, {-1, 3, 2},
	{-1, 3, 2}, {-1, 3, 2}, {-1, 3, 2}, {-1, 3, 2},
	{-1, 3, 2}, {-1, 3, 2}, {-1, 3, 2}, {-1, 3, 2},
	{-1, 3, 2}, {-1, 3, 2}, {-1, 3, 2}, {-1, 3, 2},
	{-1, 3, 2}, {-1, 3, 2}, {-1, 3, 2}, {-1, 3, 2},
	{-1, 3, 2}, {-1, 3, 2}, {-1, 3, 2}, {-1, 3, 2},
	{-1, 3, 2}, {-1, 3, 2}, {-1, 3, 2}, {-1, 3, 2},
	{-1, 3, 2}, {-1, 3, 2}, {-1, 3, 2}, {-1, 3, 2},
	{0, 2, 3}, {0, 2, 3}, {0, 2, 3}, {0, 2, 3},
	{0, 2, 3}, {0, 2, 3}, {0, 2, 3}, {0, 2, 3},
	{0, 2, 3}, {0, 2, 3}, {0, 2, 3}, {0, 2, 3},
	{0, 2, 3}, {0, 2, 3}, {0, 2, 3}, {0, 2, 3},
	{0, 3, 3}, {0, 3, 3}, {0, 3, 3}, {0, 3, 3},
	{0, 3, 3}, {0, 3, 3}, {0, 3, 3}, {0, 3, 3},
	{0, 3, 3}, {0, 3, 3}, {0, 3, 3}, {0, 3, 3},
	{0, 3, 3}, {0, 3, 3}, {0, 3, 3}, {0, 3, 3},
	{-1, -1, -1}, {-1, -1, -1}, {-1, -1, -1}, {-1, -1, -1},
	{128, 3, 0}, {136, 3, 0}, {144, 3, 0}, {152, 3, 0},
	{160, 2, 0}, {164, 2, 0}, {168, 2, 0}, {172, 2, 0},
	{176, 2, 0}, {180, 2, 0}, {184, 2, 0}, {188, 2, 0},
	{192, 1, 0}, {194, 1, 0}, {196, 1, 0}, {198, 1, 0},
	{200, 1, 0}, {202, 1, 0}, {204, 1, 0}, {206, 1, 0},
	{208, 1, 0}, {210, 1, 0}, {212, 1, 0}, {214, 1, 0},
	{216, 1, 0}, {218, 1, 0}, {220, 1, 0}, {222, 1, 0},
	{10, 4, 7}, {10, 5, 7}, {9, 4, 7}, {9, 5, 7},
	{5, 6, 7}, {5, 7, 7}, {3, 8, 7}, {3, 9, 7},
	{2, 10, 7}, {2, 11, 7}, {1, 14, 7}, {1, 15, 7},
	{1, 12, 7}, {1, 13, 7}, {0, 30, 7}, {0, 31, 7},
	{0, 28, 7}, {0, 29, 7}, {0, 26, 7}, {0, 27, 7},
	{0, 24, 7}, {0, 25, 7}, {26, 2, 7}, {26, 3, 7},
	{25, 2, 7}, {25, 3, 7}, {24, 2, 7}, {24, 3, 7},
	{23, 2, 7}, {23, 3, 7}, {22, 2, 7}, {22, 3, 7},
	{0, 22, 6}, {0, 22, 6}, {0, 23, 6}, {0, 23, 6},
	{8, 4, 6}, {8, 4, 6}, {8, 5, 6}, {8, 5, 6},
	{4, 6, 6}, {4, 6, 6}, {4, 7, 6}, {4, 7, 6},
	{0, 20, 6}, {0, 20, 6}, {0, 21, 6}, {0, 21, 6},
	{2, 8, 6}, {2, 8, 6}, {2, 9, 6}, {2, 9, 6},
	{7, 4, 6}, {7, 4, 6}, {7, 5, 6}, {7, 5, 6},
	{21, 2, 6}, {21, 2, 6}, {21, 3, 6}, {21, 3, 6},
	{20, 2, 6}, {20, 2, 6}, {20, 3, 6}, {20, 3, 6},
	{0, 18, 6}, {0, 18, 6}, {0, 19, 6}, {0, 19, 6},
	{19, 2, 6}, {19, 2, 6}, {19, 3, 6}, {19, 3, 6},
	{18, 2, 6}, {18, 2, 6}, {18, 3, 6}, {18, 3, 6},
	{1, 10, 6}, {1, 10, 6}, {1, 11, 6}, {1, 11, 6},
	{3, 6, 6}, {3, 6, 6}, {3, 7, 6}, {3, 7, 6},
	{0, 16, 6}, {0, 16, 6}, {0, 17, 6}, {0, 17, 6},
	{6, 4, 6}, {6, 4, 6}, {6, 5, 6}, {6, 5, 6},
	{17, 2, 6}, {17, 2, 6}, {17, 3, 6}, {17, 3, 6},
	{1, 36, 3}, {1, 37, 3}, {1, 34, 3}, {1, 35, 3},
	{1, 32, 3}, {1, 33, 3}, {1, 30, 3}, {1, 31, 3},
	{6, 6, 3}, {6, 7, 3}, {16, 4, 3}, {16, 5, 3},
	{15, 4, 3}, {15, 5, 3}, {14, 4, 3}, {14, 5, 3},
	{13, 4, 3}, {13, 5, 3}, {12, 4, 3}, {12, 5, 3},
	{11, 4, 3}, {11, 5, 3}, {31, 2, 3}, {31, 3, 3},
	{30, 2, 3}, {30, 3, 3}, {29, 2, 3}, {29, 3, 3},
	{28, 2, 3}, {28, 3, 3}, {27, 2, 3}, {27, 3, 3},
	{0, 80, 2}, {0, 81, 2}, {0, 78, 2}, {0, 79, 2},
	{0, 76, 2}, {0, 77, 2}, {0, 74, 2}, {0, 75, 2},
	{0, 72, 2}, {0, 73, 2}, {0, 70, 2}, {0, 71, 2},
	{0, 68, 2}, {0, 69, 2}, {0, 66, 2}, {0, 67, 2},
	{0, 64, 2}, {0, 65, 2}, {1, 28, 2}, {1, 29, 2},
	{1, 26, 2}, {1, 27, 2}, {1, 24, 2}, {1, 25, 2},
	{1, 22, 2}, {1, 23, 2}, {1, 20, 2}, {1, 21, 2},
	{1, 18, 2}, {1, 19, 2}, {1, 16, 2}, {1, 17, 2},
	{0, 62, 1}, {0, 63, 1}, {0, 60, 1}, {0, 61, 1},
	{0, 58, 1}, {0, 59, 1}, {0, 56, 1}, {0, 57, 1},
	{0, 54, 1}, {0, 55, 1}, {0, 52, 1}, {0, 53, 1},
	{0, 50, 1}, {0, 51, 1}, {0, 48, 1}, {0, 49, 1},
	{0, 46, 1}, {0, 47, 1}, {0, 44, 1}, {0, 45, 1},
	{0, 42, 1}, {0, 43, 1}, {0, 40, 1}, {0, 41, 1},
	{0, 38, 1}, {0, 39, 1}, {0, 36, 1}, {0, 37, 1},
	{0, 34, 1}, {0, 35, 1}, {0, 32, 1}, {0, 33, 1},
	{16, 2, 4}, {16, 3, 4}, {5, 4, 4}, {5, 5, 4},
	{0, 14, 4}, {0, 15, 4}, {2, 6, 4}, {2, 7, 4},
	{1, 8, 4}, {1, 9, 4}, {15, 2, 4}, {15, 3, 4},
	{14, 2, 4}, {14, 3, 4}, {4, 4, 4}, {4, 5, 4},
	{2, 4, 1}, {2, 5, 1}, {9, 2, 1}, {9, 3, 1},
	{0, 8, 1}, {0, 9, 1}, {8, 2, 1}, {8, 3, 1},
	{13, 2, 2}, {13, 3, 2}, {0, 12, 2}, {0, 13, 2},
	{12, 2, 2}, {12, 3, 2}, {11, 2, 2}, {11, 3, 2},
	{3, 4, 2}, {3, 5, 2}, {1, 6, 2}, {1, 7, 2},
	{0, 10, 2}, {0, 11, 2}, {10, 2, 2}, {10, 3, 2},
};

/** This is generated from standard document using makevld.rb. */
static const vlc_dct_t m2d_dct_table1_bit7[414] = {
	{128, 10, 0}, {352, 4, 0}, {-1, 0, 6}, {-1, 0, 6},
	{368, 1, 0}, {370, 1, 0}, {372, 1, 0}, {374, 1, 0},
	{0, 14, 7}, {0, 15, 7}, {0, 12, 7}, {0, 13, 7},
	{4, 2, 7}, {4, 3, 7}, {5, 2, 7}, {5, 3, 7},
	{376, 2, 0}, {380, 2, 0}, {384, 2, 0}, {388, 2, 0},
	{2, 2, 6}, {2, 2, 6}, {2, 3, 6}, {2, 3, 6},
	{1, 4, 6}, {1, 4, 6}, {1, 5, 6}, {1, 5, 6},
	{3, 2, 6}, {3, 2, 6}, {3, 3, 6}, {3, 3, 6},
	{1, 2, 4}, {1, 2, 4}, {1, 2, 4}, {1, 2, 4},
	{1, 2, 4}, {1, 2, 4}, {1, 2, 4}, {1, 2, 4},
	{1, 3, 4}, {1, 3, 4}, {1, 3, 4}, {1, 3, 4},
	{1, 3, 4}, {1, 3, 4}, {1, 3, 4}, {1, 3, 4},
	{-1, 3, 4}, {-1, 3, 4}, {-1, 3, 4}, {-1, 3, 4},
	{-1, 3, 4}, {-1, 3, 4}, {-1, 3, 4}, {-1, 3, 4},
	{0, 6, 5}, {0, 6, 5}, {0, 6, 5}, {0, 6, 5},
	{0, 7, 5}, {0, 7, 5}, {0, 7, 5}, {0, 7, 5},
	{0, 2, 3}, {0, 2, 3}, {0, 2, 3}, {0, 2, 3},
	{0, 2, 3}, {0, 2, 3}, {0, 2, 3}, {0, 2, 3},
	{0, 2, 3}, {0, 2, 3}, {0, 2, 3}, {0, 2, 3},
	{0, 2, 3}, {0, 2, 3}, {0, 2, 3}, {0, 2, 3},
	{0, 3, 3}, {0, 3, 3}, {0, 3, 3}, {0, 3, 3},
	{0, 3, 3}, {0, 3, 3}, {0, 3, 3}, {0, 3, 3},
	{0, 3, 3}, {0, 3, 3}, {0, 3, 3}, {0, 3, 3},
	{0, 3, 3}, {0, 3, 3}, {0, 3, 3}, {0, 3, 3},
	{0, 4, 4}, {0, 4, 4}, {0, 4, 4}, {0, 4, 4},
	{0, 4, 4}, {0, 4, 4}, {0, 4, 4}, {0, 4, 4},
	{0, 5, 4}, {0, 5, 4}, {0, 5, 4}, {0, 5, 4},
	{0, 5, 4}, {0, 5, 4}, {0, 5, 4}, {0, 5, 4},
	{0, 8, 6}, {0, 8, 6}, {0, 9, 6}, {0, 9, 6},
	{0, 10, 6}, {0, 10, 6}, {0, 11, 6}, {0, 11, 6},
	{392, 1, 0}, {394, 1, 0}, {396, 1, 0}, {398, 1, 0},
	{400, 1, 0}, {402, 2, 0}, {406, 2, 0}, {410, 2, 0},
	{-1, -1, -1}, {-1, -1, -1}, {-1, -1, -1}, {-1, -1, -1},
	{128, 3, 0}, {136, 3, 0}, {144, 3, 0}, {152, 3, 0},
	{160, 2, 0}, {164, 2, 0}, {168, 2, 0}, {172, 2, 0},
	{176, 2, 0}, {180, 2, 0}, {184, 2, 0}, {188, 2, 0},
	{192, 1, 0}, {194, 1, 0}, {196, 1, 0}, {198, 1, 0},
	{200, 1, 0}, {202, 1, 0}, {204, 1, 0}, {206, 1, 0},
	{208, 1, 0}, {210, 1, 0}, {212, 1, 0}, {214, 1, 0},
	{216, 1, 0}, {218, 1, 0}, {220, 1, 0}, {222, 1, 0},
	{10, 4, 7}, {10, 5, 7}, {9, 4, 7}, {9, 5, 7},
	{5, 6, 7}, {5, 7, 7}, {3, 8, 7}, {3, 9, 7},
	{2, 10, 7}, {2, 11, 7}, {1, 14, 7}, {1, 15, 7},
	{1, 12, 7}, {1, 13, 7}, {-1, -1, -1}, {-1, -1, -1},
	{-1, -1, -1}, {-1, -1, -1}, {-1, -1, -1}, {-1, -1, -1},
	{-1, -1, -1}, {-1, -1, -1}, {26, 2, 7}, {26, 3, 7},
	{25, 2, 7}, {25, 3, 7}, {24, 2, 7}, {24, 3, 7},
	{23, 2, 7}, {23, 3, 7}, {22, 2, 7}, {22, 3, 7},
	{-1, -1, -1}, {-1, -1, -1}, {-1, -1, -1}, {-1, -1, -1},
	{8, 4, 6}, {8, 4, 6}, {8, 5, 6}, {8, 5, 6},
	{4, 6, 6}, {4, 6, 6}, {4, 7, 6}, {4, 7, 6},
	{-1, -1, -1}, {-1, -1, -1}, {-1, -1, -1}, {-1, -1, -1},
	{-1, -1, -1}, {-1, -1, -1}, {-1, -1, -1}, {-1, -1, -1},
	{7, 4, 6}, {7, 4, 6}, {7, 5, 6}, {7, 5, 6},
	{21, 2, 6}, {21, 2, 6}, {21, 3, 6}, {21, 3, 6},
	{20, 2, 6}, {20, 2, 6}, {20, 3, 6}, {20, 3, 6},
	{-1, -1, -1}, {-1, -1, -1}, {-1, -1, -1}, {-1, -1, -1},
	{19, 2, 6}, {19, 2, 6}, {19, 3, 6}, {19, 3, 6},
	{18, 2, 6}, {18, 2, 6}, {18, 3, 6}, {18, 3, 6},
	{-1, -1, -1}, {-1, -1, -1}, {-1, -1, -1}, {-1, -1, -1},
	{3, 6, 6}, {3, 6, 6}, {3, 7, 6}, {3, 7, 6},
	{-1, -1, -1}, {-1, -1, -1}, {-1, -1, -1}, {-1, -1, -1},
	{6, 4, 6}, {6, 4, 6}, {6, 5, 6}, {6, 5, 6},
	{17, 2, 6}, {17, 2, 6}, {17, 3, 6}, {17, 3, 6},
	{1, 36, 3}, {1, 37, 3}, {1, 34, 3}, {1, 35, 3},
	{1, 32, 3}, {1, 33, 3}, {1, 30, 3}, {1, 31, 3},
	{6, 6, 3}, {6, 7, 3}, {16, 4, 3}, {16, 5, 3},
	{15, 4, 3}, {15, 5, 3}, {14, 4, 3}, {14, 5, 3},
	{13, 4, 3}, {13, 5, 3}, {12, 4, 3}, {12, 5, 3},
	{11, 4, 3}, {11, 5, 3}, {31, 2, 3}, {31, 3, 3},
	{30, 2, 3}, {30, 3, 3}, {29, 2, 3}, {29, 3, 3},
	{28, 2, 3}, {28, 3, 3}, {27, 2, 3}, {27, 3, 3},
	{0, 80, 2}, {0, 81, 2}, {0, 78, 2}, {0, 79, 2},
	{0, 76, 2}, {0, 77, 2}, {0, 74, 2}, {0, 75, 2},
	{0, 72, 2}, {0, 73, 2}, {0, 70, 2}, {0, 71, 2},
	{0, 68, 2}, {0, 69, 2}, {0, 66, 2}, {0, 67, 2},
	{0, 64, 2}, {0, 65, 2}, {1, 28, 2}, {1, 29, 2},
	{1, 26, 2}, {1, 27, 2}, {1, 24, 2}, {1, 25, 2},
	{1, 22, 2}, {1, 23, 2}, {1, 20, 2}, {1, 21, 2},
	{1, 18, 2}, {1, 19, 2}, {1, 16, 2}, {1, 17, 2},
	{0, 62, 1}, {0, 63, 1}, {0, 60, 1}, {0, 61, 1},
	{0, 58, 1}, {0, 59, 1}, {0, 56, 1}, {0, 57, 1},
	{0, 54, 1}, {0, 55, 1}, {0, 52, 1}, {0, 53, 1},
	{0, 50, 1}, {0, 51, 1}, {0, 48, 1}, {0, 49, 1},
	{0, 46, 1}, {0, 47, 1}, {0, 44, 1}, {0, 45, 1},
	{0, 42, 1}, {0, 43, 1}, {0, 40, 1}, {0, 41, 1},
	{0, 38, 1}, {0, 39, 1}, {0, 36, 1}, {0, 37, 1},
	{0, 34, 1}, {0, 35, 1}, {0, 32, 1}, {0, 33, 1},
	{5, 4, 3}, {5, 4, 3}, {5, 5, 3}, {5, 5, 3},
	{14, 2, 3}, {14, 2, 3}, {14, 3, 3}, {14, 3, 3},
	{2, 8, 4}, {2, 9, 4}, {16, 2, 4}, {16, 3, 4},
	{15, 2, 3}, {15, 2, 3}, {15, 3, 3}, {15, 3, 3},
	{7, 2, 1}, {7, 3, 1}, {8, 2, 1}, {8, 3, 1},
	{6, 2, 1}, {6, 3, 1}, {2, 4, 1}, {2, 5, 1},
	{1, 10, 2}, {1, 11, 2}, {11, 2, 2}, {11, 3, 2},
	{0, 22, 2}, {0, 23, 2}, {0, 20, 2}, {0, 21, 2},
	{13, 2, 2}, {13, 3, 2}, {12, 2, 2}, {12, 3, 2},
	{3, 4, 2}, {3, 5, 2}, {1, 8, 2}, {1, 9, 2},
	{9, 2, 1}, {9, 3, 1}, {1, 6, 1}, {1, 7, 1},
	{10, 2, 1}, {10, 3, 1}, {0, 16, 1}, {0, 17, 1},
	{0, 18, 1}, {0, 19, 1}, {0, 24, 2}, {0, 25, 2},
	{0, 26, 2}, {0, 27, 2}, {2, 6, 2}, {2, 7, 2},
	{4, 4, 2}, {4, 5, 2}, {0, 28, 2}, {0, 29, 2},
	{0, 30, 2}, {0, 31, 2},
};

static const vlc_dct_t * const m2d_dct_tables[4] = {
	m2d_dct_table0_bit7,
	m2d_dct_table1_bit7,
	m2d_dct_table0_bit7,
	m2d_dct_table1_bit7
};


/** This is generated from standard document using vldbuild.rb. */
static const vlc_t mb_inc_bit4[62] = {
	{16, -6}, {54, -3}, {7, 4}, {6, 4},
	{5, 3}, {5, 3}, {4, 3}, {4, 3},
	{3, 2}, {3, 2}, {3, 2}, {3, 2},
	{2, 2}, {2, 2}, {2, 2}, {2, 2},
	{0, 0}, {0, 0}, {16, -2}, {0, 0},
	{0, 0}, {0, 0}, {20, -2}, {24, -2},
	{28, -2}, {32, -1}, {34, -1}, {36, -1},
	{15, 3}, {15, 3}, {14, 3}, {14, 3},
	{0, 2}, {0, 0}, {0, 0}, {0, 0},
	{33, 2}, {32, 2}, {31, 2}, {30, 2},
	{29, 2}, {28, 2}, {27, 2}, {26, 2},
	{25, 2}, {24, 2}, {23, 2}, {22, 2},
	{21, 1}, {20, 1}, {19, 1}, {18, 1},
	{17, 1}, {16, 1}, {13, 3}, {12, 3},
	{11, 3}, {10, 3}, {9, 2}, {9, 2},
	{8, 2}, {8, 2},
};

/** This is generated from standard document using vldbuild.rb. */
static const vlc_t dct_dc_size_luma_bit5[48] = {
	{1, 2}, {1, 2}, {1, 2}, {1, 2},
	{1, 2}, {1, 2}, {1, 2}, {1, 2},
	{2, 2}, {2, 2}, {2, 2}, {2, 2},
	{2, 2}, {2, 2}, {2, 2}, {2, 2},
	{0, 3}, {0, 3}, {0, 3}, {0, 3},
	{3, 3}, {3, 3}, {3, 3}, {3, 3},
	{4, 3}, {4, 3}, {4, 3}, {4, 3},
	{5, 4}, {5, 4}, {6, 5}, {32, -4},
	{7, 1}, {7, 1}, {7, 1}, {7, 1},
	{7, 1}, {7, 1}, {7, 1}, {7, 1},
	{8, 2}, {8, 2}, {8, 2}, {8, 2},
	{9, 3}, {9, 3}, {10, 4}, {11, 4},
};

/** This is generated from standard document using vldbuild.rb. */
static const vlc_t dct_dc_size_chroma_bit4[36] = {
	{0, 2}, {0, 2}, {0, 2}, {0, 2},
	{1, 2}, {1, 2}, {1, 2}, {1, 2},
	{2, 2}, {2, 2}, {2, 2}, {2, 2},
	{3, 3}, {3, 3}, {4, 4}, {16, -6},
	{5, 1}, {5, 1}, {5, 1}, {5, 1},
	{5, 1}, {5, 1}, {5, 1}, {5, 1},
	{6, 2}, {6, 2}, {6, 2}, {6, 2},
	{7, 3}, {7, 3}, {8, 4}, {16, -2},
	{9, 1}, {9, 1}, {10, 2}, {11, 2},
};


/** This is generated from standard document using vldbuild.rb. */
	
static const vlc_t mb_type_b_bit4[22] = {
        {16, -2}, {20, -1}, {1, 4}, {9, 4},
        {2, 3}, {2, 3}, {10, 3}, {10, 3},
        {3, 2}, {3, 2}, {3, 2}, {3, 2},
        {11, 2}, {11, 2}, {11, 2}, {11, 2},
        {0, 0}, {20, 2}, {26, 2}, {25, 2},
        {27, 1}, {4, 1},
};

static const vlc_t mb_type_p_bit3[16] = {
	{8, -3}, {1, 3}, {8, 2}, {8, 2},
	{9, 1}, {9, 1}, {9, 1}, {9, 1},
	{0, 0}, {20, 3}, {24, 2}, {24, 2},
	{25, 2}, {25, 2}, {4, 2}, {4, 2},
};

/** This is generated from standard document using vldbuild.rb. */
static const vlc_t motion_code_bit5[102] = {
	{32, -5}, {64, -5}, {96, -2}, {100, -1},
	{3, 4}, {3, 4}, {-3, 4}, {-3, 4},
	{2, 3}, {2, 3}, {2, 3}, {2, 3},
	{-2, 3}, {-2, 3}, {-2, 3}, {-2, 3},
	{1, 2}, {1, 2}, {1, 2}, {1, 2},
	{1, 2}, {1, 2}, {1, 2}, {1, 2},
	{-1, 2}, {-1, 2}, {-1, 2}, {-1, 2},
	{-1, 2}, {-1, 2}, {-1, 2}, {-1, 2},
	{0, 0}, {0, 0}, {0, 0}, {0, 0},
	{0, 0}, {0, 0}, {0, 0}, {0, 0},
	{0, 0}, {0, 0}, {0, 0}, {0, 0},
	{0, 0}, {0, 0}, {0, 0}, {0, 0},
	{0, 0}, {0, 0}, {0, 0}, {0, 0},
	{0, 0}, {0, 0}, {0, 0}, {0, 0},
	{16, 5}, {-16, 5}, {15, 5}, {-15, 5},
	{14, 5}, {-14, 5}, {13, 5}, {-13, 5},
	{12, 5}, {-12, 5}, {11, 5}, {-11, 5},
	{10, 4}, {10, 4}, {-10, 4}, {-10, 4},
	{9, 4}, {9, 4}, {-9, 4}, {-9, 4},
	{8, 4}, {8, 4}, {-8, 4}, {-8, 4},
	{7, 2}, {7, 2}, {7, 2}, {7, 2},
	{7, 2}, {7, 2}, {7, 2}, {7, 2},
	{-7, 2}, {-7, 2}, {-7, 2}, {-7, 2},
	{-7, 2}, {-7, 2}, {-7, 2}, {-7, 2},
	{6, 2}, {-6, 2}, {5, 2}, {-5, 2},
	{4, 1}, {-4, 1},
};

/** This is generated from standard document using vldbuild.rb. */
static const vlc_t coded_block_pattern_bit5[84] = {
	{32, -4}, {48, -3}, {56, -3}, {64, -3},
	{72, -2}, {76, -2}, {80, -1}, {82, -1},
	{62, 5}, {2, 5}, {61, 5}, {1, 5},
	{56, 5}, {52, 5}, {44, 5}, {28, 5},
	{40, 5}, {20, 5}, {48, 5}, {12, 5},
	{32, 4}, {32, 4}, {16, 4}, {16, 4},
	{8, 4}, {8, 4}, {4, 4}, {4, 4},
	{60, 3}, {60, 3}, {60, 3}, {60, 3},
	{0, 0}, {0, 4}, {39, 4}, {27, 4},
	{59, 4}, {55, 4}, {47, 4}, {31, 4},
	{58, 3}, {58, 3}, {54, 3}, {54, 3},
	{46, 3}, {46, 3}, {30, 3}, {30, 3},
	{57, 3}, {53, 3}, {45, 3}, {29, 3},
	{38, 3}, {26, 3}, {37, 3}, {25, 3},
	{43, 3}, {23, 3}, {51, 3}, {15, 3},
	{42, 3}, {22, 3}, {50, 3}, {14, 3},
	{41, 3}, {21, 3}, {49, 3}, {13, 3},
	{35, 3}, {19, 3}, {11, 3}, {7, 3},
	{34, 2}, {18, 2}, {10, 2}, {6, 2},
	{33, 2}, {17, 2}, {9, 2}, {5, 2},
	{63, 1}, {3, 1}, {36, 1}, {24, 1},
};

#ifdef __cplusplus
}
#endif
