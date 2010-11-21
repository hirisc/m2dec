#ifndef _TXT2BIN_H_
#define _TXT2BIN_H_

#ifdef __cplusplus
extern "C" {
#endif

enum {
	/* max number of character */
	MAXTXT = 512
};

/**Convert binary digits in text expression to raw data.
 *'0', '1', and blank letter are allowed.
 *\return -1 on error. Otherwise, size of raw data in byte units.
 */
int txt2bin( const char *src, unsigned char *dst );

#ifdef __cplusplus
}
#endif

#endif /* _TXT2BIN_H_ */
